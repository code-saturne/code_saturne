!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine cs_user_initialization &
!================================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Initialize variables

! This subroutine is called at beginning of the computation
! (restart or not) before the loop time step

! This subroutine enables to initialize or modify (for restart)
!     unkown variables and time step values

! rom and viscl values are equal to ro0 and viscl0 or initialize
! by reading the restart file
! viscls and cp variables (when there are defined) have no value
! excepted if they are read from a restart file

! Physical quantities are defined in the following arrays:
!  propce (physical quantities defined at cell center),
!  propfb (physical quantities defined at border face center).
!
! Examples:
!  propce(iel, ipproc(irom  )) means rom  (iel)
!  propce(iel, ipproc(iviscl)) means viscl(iel)
!  propce(iel, ipproc(icp   )) means cp   (iel)
!  propce(iel, ipproc(ivisls(iscal))) means visls(iel, iscal)
!  propfb(ifac, ipprob(irom )) means romb  (ifac)

! Modification of the behaviour law of physical quantities (rom, viscl,
! viscls, cp) is not done here. It is the purpose of the user subroutine
! usphyv

! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp(ncelet, *)   ! ra ! <-- ! computed variables at cell centers at current  !
!                  !    !     ! time steps                                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfb(nfabor,*)

! Local variables

!< [loc_var_dec]
integer          iel, mode
integer          iesp , idimve

double precision tinit, hinit, coefe(ngazem)

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================
!
!< [init1]
! For Joule heating by direct conduction, it stops
!   you have tot be sure that the enthalpy function H(T) is the right one
!
if ( ippmod(ieljou).ge.1 ) then

  write(nfecra,9010)
  call csexit (1)

! For electric arc, we continue because the value are given by defauft
!       (H(T) is given from the data file dp_ELE)
elseif (ippmod(ielarc).ge.1) then

  if (ntcabs.eq.1) then
    write(nfecra,9011)
  endif

  return

endif
!< [init1]

 9010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION : Stop in the definition of Thermal properties  ',/,&
'@    =========                                               ',/,&
'@                      for Electric module                   ',/,&
'@                                                            ',/,&
'@     The user routine uselph has to be completed            ',/,&
'@                                                            ',/,&
'@     This user routine is used to define thermal properties ',/,&
'@     It is unavoidable.                                     ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(/,                                                   &
' ELECTRIC ARC MODULE : THERMAL PROPERTIES ARE READ IN A FILE',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!---------------
! Initialization
!---------------

!< [init2]
allocate(lstelt(ncel)) ! temporary array for cells selection

! Control output

write(nfecra,9001)

!===============================================================================
! Initialization
!    (only at the beginning of the calculation)
!===============================================================================

if ( isuite.eq.0 ) then

! --> Enthalpy = H(T0) ou 0

!     For electric arc,
!     for the whole compution domain enthalpy is set to H(T0) of the 1st
!     constituant of the gas
!
!     For Joule jeating by direct conduction,
!     enthalpy is set to zero, and the user will enter his H(T) function
!     tabulation.
!
!   --  HINIT calculations

  if ( ippmod(ielarc).ge.1 ) then
    mode = -1
    tinit = t0
    coefe(1) = 1.d0
    if ( ngazg .gt. 1 ) then
      do iesp = 2, ngazg
        coefe(iesp) = 0.d0
      enddo
    endif
    call elthht(mode,ngazg,coefe,hinit,tinit)
  else
    mode = -1
    tinit = t0
    call usthht(mode,hinit,tinit)
  endif

!    -- Entahlpy value

  do iel = 1, ncel
    rtp(iel,isca(ihm)) = hinit
  enddo


! --> Mass fraction  = 1 ou 0

  if ( ngazg .gt. 1 ) then
    do iel = 1, ncel
      rtp(iel,isca(iycoel(1))) = 1.d0
    enddo
    do iesp = 2, ngazg-1
      do iel = 1, ncel
        rtp(iel,isca(iycoel(iesp))) = 0.d0
      enddo
    enddo
  endif


! --> Electric potentials = 0

!     -- Real Component
  do iel = 1, ncel
    rtp(iel,isca(ipotr)) = 0.d0
  enddo

!     -- Imaginary (for Joule heating by direct conduction)
  if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
    do iel = 1, ncel
      rtp(iel,isca(ipoti)) = 0.d0
    enddo
  endif

!     -- Vector potential (3D electric arc 3D)
  if ( ippmod(ielarc).ge.2 ) then
    do idimve = 1, ndimve
      do iel = 1, ncel
        rtp(iel,isca(ipotva(idimve))) = 0.d0
      enddo
    enddo
  endif

endif
!< [init2]

!--------
! Formats
!--------

 9001 format(                                                   /,&
'                       ELECTRIC MODULE'                       ,/,&
'  variables initialization by user'                           ,/,&
                                                                /)

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_initialization
