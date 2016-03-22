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
! --------

!> \file cs_user_head_losses.f90
!>
!> \brief Define Head losses
!>
!> The subroutine cs_user_head_losses is called at three different stages
!> in the code (iappel = 1, 2 or 3)
!>
!> iappel = 1:
!>    Calculation of the number of cells where a head loss term is
!>    imposed: ncepdp
!>    Called once at the beginning of the calculation
!>
!> iappel = 2
!>    Identification of the cells where a head loss term is imposed:
!>    array icepdc(ncepdc)
!>    Called once at the beginning of the calculation
!>
!> iappel = 3
!>    Calculation of the values of the head loss term
!>    Called at each time step
!>
!> Note that calling this subroutine completely overwrites head losses
!> defined using the GUI.
!>
!> ckupdc is the local head loss term
!>
!> It appears on the momentum as follows:
!>    rho du/dt = - grad p + head_loss        (+ other terms)
!>                      with head_loss = - rho ckupdc u (in kg/(m2 s2))
!>
!> For a distributed head loss,
!>
!>    let ksil = dhl/(0.5 rho u**2) given by the litterature
!>    (dhl is the head loss per unit length)
!>
!>    the source term tspdc is equal to dhl = - ksil *(0.5 rho u**2)
!>
!>    we have ckupdc = 0.5 ksil abs(U)
!>
!>
!> For a singular head loss,
!>
!>    let ksil = dhs/(0.5 rho u**2) given by the litterature
!>    (dhs is the singular head loss)
!>
!>    the source term tspdc is equal to dhs/L = - ksil/L *(0.5 rho u**2)
!>
!>    we have ckupdc = 0.5 ksil/L abs(u)
!>
!>    where L is the length over which we have chosen to represent the
!>    singular head loss
!>
!>
!> Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     iappel        stage in the code
!> \param[in]     icepdc        numbers of ncepdp cells with head loss
!> \param[in]     izcpdc        cells zone for head loss definition
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head loss
!_______________________________________________________________________________


subroutine cs_user_head_losses &
 ( ncepdp , iappel ,                                              &
   icepdc , izcpdc ,                                              &
   dt     ,                                                       &
   ckupdc )

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

!< [arg]
integer          ncepdp
integer          iappel

integer          icepdc(ncepdp)
integer          izcpdc(ncel)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6)
!< [arg]

! Local variables

!< [loc_var_dec]
integer          iel, ielpdc, ikpdc
integer          ilelt, nlelt
integer          izone

double precision alpha, cosalp, sinalp, vit, ck1, ck2

double precision, dimension(:,:), pointer :: cvara_vel

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================
!< [map_field_array]
! Map field arrays
call field_get_val_prev_v(ivarfl(iu), cvara_vel)
!< [map_field_array]

!< [allocate]
! Allocate a temporary array for cells selection
allocate(lstelt(ncel))
!< [allocate]

!< [start_1]
if (iappel.eq.1.or.iappel.eq.2) then
!< [start_1]
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
!< [example_1]
  ielpdc = 0
!< [example_1]


  ! Example 2: head losses defined by coodinates for zone
  !            (4 <= x <=6; 2 <= y <= 8).
  !            No head losses else

!< [example_2]
  if (1.eq.0) then ! This test allows deactivating the example

    izone = 0
    ielpdc = 0

    call getcel('X <= 6.0 and X >= 4.0 and Y >= 2.0 and'//      &
                'Y <= 8.0',nlelt,lstelt)
!< [example_2]

!< [start_2]
    izone = izone + 1

    do ilelt = 1, nlelt
      iel = lstelt(ilelt)
      izcpdc(iel) = izone
      ielpdc = ielpdc + 1
      if (iappel.eq.2) icepdc(ielpdc) = iel
    enddo

  endif
!< [start_2]



! Generic subsection that should not be modified
! ----------------------------------------------

! For iappel = 1,
!    Define ncepdp, the number of cells with head losses
!    This is valid for both examples above

!< [generic_subsection_1]
  if (iappel.eq.1) then
    ncepdp = ielpdc
  endif
!< [generic_subsection_1]
!-------------------------------------------------------------------------------

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

  ! --- Diagonal tensor
  !     Example of head losses in direction x

!< [example_3]

  do ielpdc = 1, ncepdp
    iel=icepdc(ielpdc)
    vit = sqrt(cvara_vel(1,iel)**2 + cvara_vel(2,iel)**2 + cvara_vel(3,iel)**2)
    ckupdc(ielpdc,1) = 10.d0*vit
    ckupdc(ielpdc,2) =  0.d0*vit
    ckupdc(ielpdc,3) =  0.d0*vit
  enddo
!< [example_3]

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
!< [example_4]

  alpha  = pi/4.d0
  cosalp = cos(alpha)
  sinalp = sin(alpha)
  ck1 = 10.d0
  ck2 =  0.d0
!< [example_4]

!< [filling]
  do ielpdc = 1, ncepdp
    iel = icepdc(ielpdc)
    vit = sqrt(cvara_vel(1,iel)**2 + cvara_vel(2,iel)**2 + cvara_vel(3,iel)**2)
    ckupdc(ielpdc,1) = (cosalp**2*ck1 + sinalp**2*ck2)*vit
    ckupdc(ielpdc,2) = (sinalp**2*ck1 + cosalp**2*ck2)*vit
    ckupdc(ielpdc,3) =  0.d0
    ckupdc(ielpdc,4) = cosalp*sinalp*(-ck1+ck2)*vit
    ckupdc(ielpdc,5) =  0.d0
    ckupdc(ielpdc,6) =  0.d0
  enddo

  !-----------------------------------------------------------------------------

endif
!< [filling]


!< [deallocate]
! Deallocate the temporary array
deallocate(lstelt)
!< [deallocate]

return
end subroutine cs_user_head_losses