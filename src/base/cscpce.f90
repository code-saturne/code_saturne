!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!> \param[in]     nptdis        number of distant points
!> \param[in]     f_id          field index
!> \param[in]     f_dim         field dimension
!> \param[in]     locpts        connectivity
!> \param[in]     coopts        coordinates of the distants points
!> \param[out]    rvdis         work array for the variable to be exchanged
!______________________________________________________________________________

subroutine cscpce &
 ( nptdis , f_id   , f_dim ,                                      &
   locpts ,                                                       &
   coopts , rvdis  )

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

integer          nptdis
integer          f_id
integer          f_dim

integer          locpts(nptdis)

double precision coopts(3,nptdis), rvdis(f_dim,nptdis)

! Local variables

integer          ipt    , iel    , isou   , iprev
integer          inc    , iccocg

double precision dx     , dy     , dz

double precision, dimension(:,:), allocatable :: grads
double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:), pointer :: cvara_s
double precision, dimension(:,:), pointer :: cvara_v


!===============================================================================

inc    = 1
iprev  = 1 ! previsous time step value
iccocg = 1

! For scalars
if (f_dim.eq.1) then
  call field_get_val_prev_s(f_id, cvara_s)

  ! Allocate a temporary array
  allocate(grads(3, ncelet))

  call field_gradient_scalar(f_id, iprev, 0, inc, iccocg, grads)

  ! --- Interpolation
  do ipt = 1, nptdis

    iel = locpts(ipt)

    dx = coopts(1,ipt) - xyzcen(1,iel)
    dy = coopts(2,ipt) - xyzcen(2,iel)
    dz = coopts(3,ipt) - xyzcen(3,iel)

    ! FIXME remove reconstruction ?
    rvdis(1,ipt) = cvara_s(iel) + grads(1,iel)*dx       &
                                + grads(2,iel)*dy       &
                                + grads(3,iel)*dz
  enddo

  ! Free memory
  deallocate(grads)

else if (f_dim.eq.3) then
  call field_get_val_prev_v(f_id, cvara_v)

  ! Allocate a temporary array
  allocate(gradv(3, 3, ncelet))

  call field_gradient_vector(f_id, iprev, 0, inc, gradv)

  ! --- Interpolation
  do ipt = 1, nptdis

    iel = locpts(ipt)

    dx = coopts(1,ipt) - xyzcen(1,iel)
    dy = coopts(2,ipt) - xyzcen(2,iel)
    dz = coopts(3,ipt) - xyzcen(3,iel)

    ! FIXME remove reconstruction ?
    do isou = 1, f_dim
      rvdis(isou,ipt) = cvara_v(isou,iel) + gradv(1,isou,iel)*dx       &
                                          + gradv(2,isou,iel)*dy       &
                                          + gradv(3,isou,iel)*dz
    enddo

  enddo

  ! Free memory
  deallocate(gradv)

else
  !TODO for tensors
  call csexit(1)
endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
