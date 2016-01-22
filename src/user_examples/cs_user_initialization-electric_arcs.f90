!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_initialization-electrics_arcs.f90
!> \brief Electrics arcs example
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

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
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, mode
integer          iesp , idimve

double precision tinit, hinit, coefe(ngazem)
character(len=80) :: f_name

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cvar_scalt, cvar_ycoel
double precision, dimension(:), pointer :: cvar_potr, cvar_poti, cvar_potva
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
'@     The user routine cs_user_physical_properties           ',/,&
'@     has to be completed                                    ',/,&
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
    if ( ngazge .gt. 1 ) then
      do iesp = 2, ngazge
        coefe(iesp) = 0.d0
      enddo
    endif
    call elthht(mode,coefe,hinit,tinit)
  else
    mode = -1
    tinit = t0
    call usthht(mode,hinit,tinit)
  endif

!    -- Entahlpy value

  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
  do iel = 1, ncel
    cvar_scalt(iel) = hinit
  enddo


! --> Mass fraction  = 1 ou 0

  if ( ngazge .gt. 1 ) then
    write(f_name,'(a13,i2.2)') 'esl_fraction_', 1
    call field_get_val_s_by_name(trim(f_name), cvar_ycoel)
    do iel = 1, ncel
      cvar_ycoel(iel) = 1.d0
    enddo
    do iesp = 2, ngazge-1
      write(f_name,'(a13,i2.2)') 'esl_fraction_',iesp
      call field_get_val_s_by_name(trim(f_name), cvar_ycoel)
      do iel = 1, ncel
        cvar_ycoel(iel) = 0.d0
      enddo
    enddo
  endif


! --> Electric potentials = 0

!     -- Real Component
  call field_get_val_s_by_name('elec_pot_r', cvar_potr)
  do iel = 1, ncel
    cvar_potr(iel) = 0.d0
  enddo

!     -- Imaginary (for Joule heating by direct conduction)
  if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
    call field_get_val_s_by_name('elec_pot_r', cvar_poti)
    do iel = 1, ncel
      cvar_poti(iel) = 0.d0
    enddo
  endif

!     -- Vector potential (3D electric arc 3D)
  if ( ippmod(ielarc).ge.2 ) then
    call field_get_val_s_by_name('vec_potential_01', cvar_potva)
    do iel = 1, ncel
      cvar_potva(iel) = 0.d0
    enddo
    call field_get_val_s_by_name('vec_potential_02', cvar_potva)
    do iel = 1, ncel
      cvar_potva(iel) = 0.d0
    enddo
    call field_get_val_s_by_name('vec_potential_03', cvar_potva)
    do iel = 1, ncel
      cvar_potva(iel) = 0.d0
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
