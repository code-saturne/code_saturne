!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
! --------
!>
!> \file cnlevm.f90
!>
!> \brief Calculation of non linear terms of the quadratic k-epsilon model
!>        (Baglietto et al.)
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     rij          array containing the non linear terms of
!>                             quadratic Boussinesq approximation
!______________________________________________________________________________!

subroutine cnlevm &
 ( rij)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

double precision rij(6,ncelet)

! Local variables

integer iel
integer ii, jj, kk
integer inc, iprev

double precision xstrai(3,3), xrotac(3,3)
double precision sikskj(3,3), wikskj(3,3)
double precision skiwjk(3,3), wikwjk(3,3)
double precision xrij(3,3)
double precision sijsij
double precision xqc1, xqc2, xqc3
double precision xk , xeps, xvisct
double precision xttke, xcmu, xss
double precision d1s2, d2s3

double precision, allocatable, dimension(:,:,:) :: gradv

double precision, dimension(:), pointer :: cvar_k, cvar_ep
double precision, dimension(:), pointer :: visct

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_val_s(ivisct, visct)
call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iep), cvar_ep)

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, 0, inc, gradv)

d1s2 = .5d0
d2s3 = 2.d0/3.d0

!===============================================================================
! 2. Computation
!===============================================================================

do iel = 1, ncelet
  xvisct = visct(iel)
  xeps   = cvar_ep(iel)
  xk     = cvar_k(iel)
  xttke  = xk/xeps

  ! Sij
  xstrai(1,1) = gradv(1, 1, iel)
  xstrai(1,2) = d1s2*(gradv(2, 1, iel)+gradv(1, 2, iel))
  xstrai(1,3) = d1s2*(gradv(3, 1, iel)+gradv(1, 3, iel))
  xstrai(2,1) = xstrai(1,2)
  xstrai(2,2) = gradv(2, 2, iel)
  xstrai(2,3) = d1s2*(gradv(3, 2, iel)+gradv(2, 3, iel))
  xstrai(3,1) = xstrai(1,3)
  xstrai(3,2) = xstrai(2,3)
  xstrai(3,3) = gradv(3, 3, iel)
  ! omegaij
  xrotac(1,1) = 0.d0
  xrotac(1,2) = d1s2*(gradv(2, 1, iel)-gradv(1, 2, iel))
  xrotac(1,3) = d1s2*(gradv(3, 1, iel)-gradv(1, 3, iel))
  xrotac(2,1) = -xrotac(1,2)
  xrotac(2,2) = 0.d0
  xrotac(2,3) = d1s2*(gradv(3, 2, iel)-gradv(2, 3, iel))
  xrotac(3,1) = -xrotac(1,3)
  xrotac(3,2) = -xrotac(2,3)
  xrotac(3,3) = 0.d0

  sijsij    = 0.d0
  do ii = 1,3
    do jj = 1,3
      sijsij        = sijsij + xstrai(ii,jj)*xstrai(ii,jj)
      sikskj(ii,jj) = 0.d0
      wikskj(ii,jj) = 0.d0
      skiwjk(ii,jj) = 0.d0
      wikwjk(ii,jj) = 0.d0
      do kk = 1,3
        sikskj(ii,jj) = sikskj(ii,jj) + xstrai(ii,kk)*xstrai(kk,jj)
        wikskj(ii,jj) = wikskj(ii,jj) + xrotac(ii,kk)*xstrai(kk,jj)
        skiwjk(ii,jj) = skiwjk(ii,jj) + xstrai(kk,ii)*xrotac(jj,kk)
        wikwjk(ii,jj) = wikwjk(ii,jj) + xrotac(ii,kk)*xrotac(jj,kk)
      end do
    end do
  end do

  xss  = xttke*sqrt(.5d0*sijsij)
  xcmu = d2s3/(3.9d0 + xss)

  ! Evaluating "constants"
  xqc1 = cnl1/((cnl4 + cnl5*xss**3.d0)*xcmu)
  xqc2 = cnl2/((cnl4 + cnl5*xss**3.d0)*xcmu)
  xqc3 = cnl3/((cnl4 + cnl5*xss**3.d0)*xcmu)

  do ii = 1,3
    do jj= 1,3
      xrij(ii,jj) = xqc1*xvisct*xttke*sikskj(ii,jj)                 &
                  + xqc2*xvisct*xttke*(wikskj(ii,jj)+skiwjk(ii,jj)) &
                  + xqc3*xvisct*xttke*wikwjk(ii,jj)
    enddo
  enddo

  rij(1,iel) = xrij(1,1)
  rij(2,iel) = xrij(2,2)
  rij(3,iel) = xrij(3,3)
  rij(4,iel) = xrij(1,2)
  rij(5,iel) = xrij(2,3)
  rij(6,iel) = xrij(1,3)

enddo

! Free memory
deallocate(gradv)

!----
! Format
!----


!----
! End
!----

return
end subroutine

