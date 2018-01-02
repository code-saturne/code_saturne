!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
! Function:
! ---------
!> \file tsepls.f90
!> \brief Calculation of the E term of the \f$ \varepsilon\f$ equation
!>        (BL-V2/K model)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in,out] w1            work array to store the E-term
!______________________________________________________________________________!

subroutine tsepls &
 ( w1     )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

double precision w1(ncelet)

! Local variables

integer          iel, ifac, inc, iprev
integer          isou, ii, jj, j ,k

double precision w_temp, pnd

double precision, allocatable, dimension(:,:,:) :: w7
double precision, allocatable, dimension(:,:,:) :: gradv

double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu

double precision :: duidxk(3),njsj(3)
!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary for the gradient calculation
allocate(gradv(3, 3, ncelet))

! Allocate work arrays
allocate(w7(ncelet, 3, 3))

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

!===============================================================================
! 2. Calculation of the term d2Ui/dxkdxj*d2Ui/dxkdxj
!===============================================================================

do iel = 1, ncel
  w1(iel) = 0.0d0
enddo

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

! Loop over u, v, w components
do isou = 1, 3

  do iel = 1, ncel
    do k = 1,3
      do j = 1,3
        w7(iel,j,k) = .0d0
      enddo
    enddo
  enddo

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    pnd = pond(ifac)
    duidxk(1) = pnd * gradv(1, isou, ii) + (1.0d0 - pnd) * gradv(1, isou, jj)
    duidxk(2) = pnd * gradv(2, isou, ii) + (1.0d0 - pnd) * gradv(2, isou, jj)
    duidxk(3) = pnd * gradv(3, isou, ii) + (1.0d0 - pnd) * gradv(3, isou, jj)
    njsj(1)   = surfac(1,ifac)
    njsj(2)   = surfac(2,ifac)
    njsj(3)   = surfac(3,ifac)
    do k = 1,3
      do j = 1,3
        w7(ii,j,k) =  w7(ii,j,k) + duidxk(k)*njsj(j)
        w7(jj,j,k) =  w7(jj,j,k) - duidxk(k)*njsj(j)
      end do
    end do

  enddo

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    duidxk(1) = gradv(1, isou, ii)
    duidxk(2) = gradv(2, isou, ii)
    duidxk(3) = gradv(3, isou, ii)
    njsj(1)   = surfbo(1,ifac)
    njsj(2)   = surfbo(2,ifac)
    njsj(3)   = surfbo(3,ifac)
    do k = 1,3
      do j = 1,3
        w7(ii,j,k) =  w7(ii,j,k) + duidxk(k)*njsj(j)
      end do
    end do

  enddo

  do iel = 1, ncel
    w_temp = .0d0
    do k = 1,3
      do j = 1,3
        w_temp = w_temp + (w7(iel,j,k)/volume(iel))**2.d0
      end do
    end do
    w1(iel) = w1(iel) + w_temp
  enddo

enddo

! Free memory
deallocate(gradv)
deallocate(w7)

!----
! End
!----

return

end subroutine
