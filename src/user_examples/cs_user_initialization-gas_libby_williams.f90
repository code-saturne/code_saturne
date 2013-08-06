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
   dt     , rtp    , propce , propfa , propfb )

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
!  propfa (physical quantities defined at interior face center),
!  propfa (physical quantities defined at border face center).
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
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
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

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)

! Local variables

integer          iel, mode, igg, izone
double precision hinit, coefg(ngazgm)
double precision sommqf, sommqt, sommq, tentm, fmelm


character*80     chaine
integer          iscal, ivar, ii
double precision valmax, valmin

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!---------------
! Initialization
!---------------

allocate(lstelt(ncel)) ! temporary array for cells selection

! Control output

write(nfecra,9001)

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if ( isuite.eq.0 ) then

  ! a. Preliminary calculations

  sommqf = zero
  sommq  = zero
  sommqt = zero

  !  Deals with multiple inlets
  do izone = 1, nozapm
    sommqf = sommqf + qimp(izone)*fment(izone)
    sommqt = sommqt + qimp(izone)*tkent(izone)
    sommq  = sommq  + qimp(izone)
  enddo

  if (abs(sommq).gt.epzero) then
    fmelm = sommqf / sommq
    tentm = sommqt / sommq
  else
    fmelm = zero
    tentm = t0
  endif

  ! ----- Calculation of the Enthalpy of the gas mixture
  !       (unburned mean gas)

  if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3               &
                           .or. ippmod(icolwc).eq.5 ) then
    coefg(1) = fmelm
    coefg(2) = (1.d0-fmelm)
    coefg(3) = zero
    mode     = -1

    !   Converting the mean temperatur boundary conditions into
    !   enthalpy values
    call cothht                                                   &
    !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                       &
        npo    , npot   , th     , ehgazg ,                       &
        hinit  , tentm )
  endif

  do iel = 1, ncel

    ! b. Initialisation

    ! Mass fraction of Unburned (fresh) Gas
    rtp(iel,isca(iyfm))  = 0.0d0*fmelm
    ! Mean Mixture Fraction
    rtp(iel,isca(ifm))   = 0.d0*fmelm
    ! Variance of fuel Mass fraction
    rtp(iel,isca(iyfp2m)) = zero
    ! Variance of Mixture Fraction
    rtp(iel,isca(ifp2m))  = zero

    ! Covariance for NDIRAC >= 3

    if ( ippmod(icolwc).ge. 2 ) then
      rtp(iel,isca(icoyfp))   = zero
    endif

    ! Enthalpy

    if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3             &
                             .or. ippmod(icolwc).eq.5 ) then
      rtp(iel,isca(ihm)) = hinit
    endif

  enddo

  ! ---> Control Output of the user defined initialisation values

  write(nfecra,2000)

  do ii  = 1, nscapp
    iscal = iscapp(ii)
    ivar  = isca(iscal)
    valmax = -grand
    valmin =  grand
    do iel = 1, ncel
      valmax = max(valmax,rtp(iel,ivar))
      valmin = min(valmin,rtp(iel,ivar))
    enddo
    if ( irangp.ge.0 ) then
      call parmax(valmax)
      call parmin(valmin)
    endif

    chaine = nomvar(ipprtp(ivar))
    write(nfecra,2010)chaine(1:8),valmin,valmax
  enddo
  write(nfecra,2020)

endif

!--------
! Formats
!--------

 9001 format(                                                   /,&
'  cs_user_initialization: variables initialisation by user'   ,/,&
                                                                /)

 2000 format(                                                   /,&
                                                                /,&
' -----------------------------------------------------------' ,/,&
                                                                /,&
' ** INITIALISATION OF VARIABLES FOR Libby-Williams model'     ,/,&
'    --------------------------------------------------------' ,/,&
'           ONLY ONE PASS'                                     ,/,&
' ---------------------------------'                           ,/,&
'  Variable  Valeur min  Valeur max'                           ,/,&
' ---------------------------------'                             )

 2010 format(                                                     &
 2x,     a8,      e12.4,      e12.4                              )

 2020 format(                                                     &
' ---------------------------------'                           ,/)

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_initialization
