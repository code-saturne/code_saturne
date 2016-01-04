!-------------------------------------------------------------------------------

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
integer          isou, ii, jj

double precision w1f, w2f, w3f
double precision pnd, flux, somsur

double precision, allocatable, dimension(:) :: w7
double precision, allocatable, dimension(:,:,:) :: gradv

double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary for the gradient calculation
allocate(gradv(3, 3, ncelet))

! Allocate work arrays
allocate(w7(ncelet))

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
    w7(iel) = 0.0d0
  enddo

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    pnd = pond(ifac)
    w1f = pnd * gradv(1, isou, ii) + (1.0d0 - pnd) * gradv(1, isou, jj)
    w2f = pnd * gradv(2, isou, ii) + (1.0d0 - pnd) * gradv(2, isou, jj)
    w3f = pnd * gradv(3, isou, ii) + (1.0d0 - pnd) * gradv(3, isou, jj)

    somsur = surfac(1,ifac) + surfac(2,ifac) + surfac(3,ifac)

    flux = (w1f + w2f + w3f)*somsur

    w7(ii) = w7(ii) + flux
    w7(jj) = w7(jj) - flux

  enddo

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    w1f = gradv(1, isou, ii)
    w2f = gradv(2, isou, ii)
    w3f = gradv(3, isou, ii)
    somsur = surfbo(1,ifac) + surfbo(2,ifac) + surfbo(3,ifac)
    flux = (w1f + w2f + w3f)*somsur
    w7(ii) = w7(ii) + flux

  enddo

  do iel = 1, ncel
    w1(iel) = w1(iel) + (w7(iel)/volume(iel))**2
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
