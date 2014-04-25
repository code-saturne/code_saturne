!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine uskpdc &
!================

 ( nvar   , nscal  ,                                              &
   ncepdp , iappel ,                                              &
   icepdc , izcpdc ,                                              &
   dt     , rtpa   , rtp    , propce ,                            &
   ckupdc )

!===============================================================================
! Purpose:
! --------

!    User subroutine.

!    Define Head losses

! The subroutine uskpdc is called at three different stages in the code
!  (iappel = 1, 2 or 3)

! iappel = 1:
!    Calculation of the number of cells where a head loss term is
!    imposed: ncepdp
!    Called once at the beginning of the calculation

! iappel = 2
!    Identification of the cells where a head loss term is imposed:
!    array icepdc(ncepdc)
!    Called once at the beginning of the calculation

! iappel = 3
!    Calculation of the values of the head loss term
!    Called at each time step

! Note that calling this subroutine completely overwrites head losses
! defined using the GUI.

! ckupdc is the local head loss term

! It appears on the momentum as follows:
!    rho du/dt = - grad p + head_loss        (+ other terms)
!                      with head_loss = - rho ckupdc u ( en kg/(m2 s))

! For a distributed head loss,

!    let ksil = dhl/(0.5 rho u**2) given by the litterature
!    (dhl est is the head loss per unit length)

!    the source term tspdc is equal to dhl = - ksil *(0.5 rho u**2)

!    we have ckupdc = 0.5 ksil abs(U)


! For a singular head loss,

!    let ksil = dhs/(0.5 rho u**2) given by the litterature
!    (dhs est is the singular head loss)

!    the source term tspdc is equal to dhs/l = - ksil/l *(0.5 rho u**2)

!    we have ckupdc = 0.5 ksis/l abs(u)

!    where l is the length over whic we have chosen to represent the
!    singular head loss


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
! icepdc(ncepdp    ! te ! <-- ! numbers of ncepdp cells with head loss         !
! izcpdc(ncelet)   ! ia ! <-- ! cells zone for head loss definition            !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc           ! tr ! <-- ! work array for head loss                       !
!  (ncepdp,6)      !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp
integer          iappel

integer          icepdc(ncepdp)
integer          izcpdc(ncel)

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6)

! Local variables

integer          iel, ielpdc, ikpdc
integer          ilelt, nlelt
integer          izone

double precision alpha, cosalp, sinalp, vit, ck1, ck2

double precision, dimension(:,:), pointer :: vela

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================

! Map field arrays
call field_get_val_prev_v(ivarfl(iu), vela)

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))


if (iappel.eq.1.or.iappel.eq.2) then

  !=============================================================================

  ! 2 calls:

  ! iappel = 1:
  !    Calculation of the number of cells where a head loss term is
  !    imposed: ncepdp
  !    Called once at the beginning of the calculation

  ! iappel = 2
  !    Identification of the cells where a head loss term is imposed:
  !    array icepdc(ncepdc)
  !    Called once at the beginning of the calculation

  ! Notes:

  !    - Do not use ckupdc in this section (it is defined with iappel = 3)
  !    - Use icepdc in this section only with (iappel = 2)

  !=============================================================================

  ! To be completed by the user: cell selection
  ! -------------------------------------------

  ! Example 1: No head loss (default)

  ielpdc = 0


  ! Example 2: head losses define by coodinates for zone
  !            (4 <= x <=6; 2 <= y <= 8).
  !            No head losses else

  if (1.eq.0) then ! This test allows deactivating the example

    izone = 0
    ielpdc = 0

    call getcel('X <= 6.0 and X >= 4.0 and Y >= 2.0 and'//      &
                'Y <= 8.0',nlelt,lstelt)

    izone = izone + 1

    do ilelt = 1, nlelt
      iel = lstelt(ilelt)
      izcpdc(iel) = izone
      ielpdc = ielpdc + 1
      if (iappel.eq.2) icepdc(ielpdc) = iel
    enddo

  endif


! Generic subsection that should not be modified
! ----------------------------------------------

! For iappel = 1,
!    Define ncepdp, the number of cells with head losses
!    This is valid for both examples above

  if (iappel.eq.1) then
    ncepdp = ielpdc
  endif

!-------------------------------------------------------------------------------

else if (iappel.eq.3) then

  !=============================================================================

  ! Third call, at each time step

  ! iappel = 3:

  !    ckupdc: compute head loss coefficients in the calculation coordinates,
  !            organized in order k11, k22, k33, k12, k13, k23

  ! Note:
  !
  !    - make sure diagonal coefficients are positive. The calculation
  !      may crash if this is not the case, and no further check will
  !      be done

  !      ===========================================================

  ! To be completed by the user: coefficient values
  ! -----------------------------------------------

  ! --- Caution
  !   It is important that all ckupdc values are defined (by zero values if
  !   necessary), as they will be used to compute a source term in cells
  !   identified previously.

  !   They are initialized by zero values,
  !   and the user must keep this initialization.
  !                     ====

  do ikpdc = 1, 6
    do ielpdc = 1, ncepdp
      ckupdc(ielpdc,ikpdc) = 0.d0
    enddo
  enddo

  ! --- Diagonal tensor
  !     Example of head losses in direction x

  if (.true.) return  ! (replace .true. with .false. or remove test to activate)

  do ielpdc = 1, ncepdp
    iel=icepdc(ielpdc)
    vit = sqrt(vela(1,iel)**2 + vela(2,iel)**2 + vela(3,iel)**2)
    ckupdc(ielpdc,1) = 10.d0*vit
    ckupdc(ielpdc,2) =  0.d0*vit
    ckupdc(ielpdc,3) =  0.d0*vit
  enddo

  ! ---  3x3 tensor
  !      Example of head losses at alpha = 45 degres x,y
  !      direction x resists by ck1 and y by ck2
  !      ck2 = 0 represents vanes as follows:  ///////
  !      in coordinate system x y

  !                 Y|    /y
  !                  |  /
  !                  |/
  !                  \--------------- X
  !                   \ / ALPHA
  !                    \
  !                     \ x

  if (.true.) return  ! (replace .true. with .false. or remove test to activate)

  alpha  = pi/4.d0
  cosalp = cos(alpha)
  sinalp = sin(alpha)
  ck1 = 10.d0
  ck2 =  0.d0

  do ielpdc = 1, ncepdp
    iel = icepdc(ielpdc)
    vit = sqrt(vela(1,iel)**2 + vela(2,iel)**2 + vela(3,iel)**2)
    ckupdc(ielpdc,1) = (cosalp**2*ck1 + sinalp**2*ck2)*vit
    ckupdc(ielpdc,2) = (sinalp**2*ck1 + cosalp**2*ck2)*vit
    ckupdc(ielpdc,3) =  0.d0
    ckupdc(ielpdc,4) = cosalp*sinalp*(-ck1+ck2)*vit
    ckupdc(ielpdc,5) =  0.d0
    ckupdc(ielpdc,6) =  0.d0
  enddo

  !-----------------------------------------------------------------------------

endif

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine uskpdc
