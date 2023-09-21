!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file cscpce.f90
!> \brief Preparation of sending velocity variables for coupling between
!> two instances of code_saturne via boundary faces.
!> Received indformation will be transformed into boundary condition
!> in subroutine \ref csc2cl.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     n_elts        number of distant points
!> \param[in]     f_id          field index
!> \param[in]     f_dim         field dimension
!> \param[in]     reverse       reverse mode
!> \param[in]     elt_ids       connectivity
!> \param[in]     coopts        coordinates of the distants points
!> \param[out]    cw1_to_send   work array for the variable to be exchanged
!> \param[out]    cw2_to_send   work array for the variable to be exchanged
!______________________________________________________________________________

subroutine cscpce &
 ( n_elts  , f_id , f_dim ,                                      &
   reverse , elt_ids ,                                                       &
   coopts  , cw1_to_send , cw2_to_send )

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
use cplsat
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          n_elts
integer          f_id
integer          f_dim
integer          reverse
integer          elt_ids(n_elts)

double precision coopts(3,n_elts), cw1_to_send(f_dim,n_elts)
double precision cw2_to_send(f_dim,f_dim,n_elts)

! Local variables

integer          ipt    , iel    , isou   , jsou, iprev
integer          inc

double precision dx     , dy     , dz
double precision xtau   , rovtau

double precision, dimension(:), pointer ::  crom
double precision, dimension(:,:), allocatable :: grads
double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:), pointer :: cvara_s
double precision, dimension(:,:), pointer :: cvara_v


!===============================================================================

inc    = 1
iprev  = 1 ! previsous time step value

call field_get_val_s(icrom,crom)
xtau = 100.d0*dtref
! For scalars
if (reverse.eq.0) then

  if (f_dim.eq.1) then
    call field_get_val_prev_s(f_id, cvara_s)

    ! Allocate a temporary array
    allocate(grads(3, ncelet))

    call field_gradient_scalar(f_id, iprev, inc, grads)

    ! --- Interpolation
    do ipt = 1, n_elts

      iel = elt_ids(ipt)

      dx = coopts(1,ipt) - xyzcen(1,iel)
      dy = coopts(2,ipt) - xyzcen(2,iel)
      dz = coopts(3,ipt) - xyzcen(3,iel)

      cw1_to_send(1,ipt) = cvara_s(iel) + grads(1,iel)*dx       &
                                         + grads(2,iel)*dy       &
                                         + grads(3,iel)*dz
    enddo

    ! Free memory
    deallocate(grads)

  else if (f_dim.eq.3) then
    call field_get_val_prev_v(f_id, cvara_v)

    ! Allocate a temporary array
    allocate(gradv(3, 3, ncelet))

    call field_gradient_vector(f_id, iprev, inc, gradv)
    ! --- Interpolation
    do ipt = 1, n_elts

      iel = elt_ids(ipt)

      dx = coopts(1,ipt) - xyzcen(1,iel)
      dy = coopts(2,ipt) - xyzcen(2,iel)
      dz = coopts(3,ipt) - xyzcen(3,iel)

      do isou = 1, f_dim
        cw1_to_send(isou,ipt) = cvara_v(isou,iel) + gradv(1,isou,iel)*dx &
                                                  + gradv(2,isou,iel)*dy &
                                                  + gradv(3,isou,iel)*dz
      enddo

    enddo

    ! Free memory
    deallocate(gradv)

  else
    !TODO for tensors
    call csexit(1)
  endif

! mode reverse = 1
else
  if (f_dim.eq.1) then
    call field_get_val_prev_s(f_id, cvara_s)

    do ipt = 1, n_elts

      iel = elt_ids(ipt)

      ! Mass weighted averaged value
      cw1_to_send(1,ipt) = cell_f_vol(iel)*crom(iel)*cvara_s(iel)
      ! Mass
      cw2_to_send(1,1,ipt) = cell_f_vol(iel)*crom(iel)

    enddo

  else if (f_dim.eq.3) then
    call field_get_val_prev_v(f_id, cvara_v)

    do ipt = 1, n_elts

      iel = elt_ids(ipt)

      do isou = 1, f_dim
        ! Mass weighted averaged value
        cw1_to_send(isou,ipt) = cell_f_vol(iel)*crom(iel)*cvara_v(isou,iel)
        ! Implicit part
        do jsou = 1, f_dim
          if (isou.eq.jsou) then
            cw2_to_send(jsou,isou,ipt) = cell_f_vol(iel)*crom(iel)
          else
            cw2_to_send(jsou,isou,ipt) = 0.d0
          endif
        enddo
      enddo

    enddo

  else
    !TODO for tensors
    call csexit(1)
  endif

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
