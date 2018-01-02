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

!> \file predvv.f90
!>
!> \brief This subroutine solves a linear system with the gauss method.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nn            matrix columns number
!> \param[in]     mm            matrix rows number
!> \param[in]     aa            linear system matrix
!> \param[out]    xx            linear system solution
!> \param[in]     bb            linear system right hand side
!_______________________________________________________________________________

subroutine gauss(nn ,mm, aa, xx, bb)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use entsor, only:nfecra

!===============================================================================

implicit none

integer          nn, mm
double precision aa(nn,mm), xx(nn), bb(nn)
double precision, allocatable, dimension (:,:) :: ap
double precision, allocatable, dimension (:)   :: bp

integer          ii, jj, kk, ll
integer          ipass, ipivot, ipermu
double precision akk, aik, aij, bi, sum, rr, residu

!===============================================================================

if (nn.ne.mm) then
  write(nfecra,*) &
       "Gaussian elimination (gauss.f90): the matrix is not triangular, ", &
       "stop the calculation."
  call csexit(1)
endif

! --- save the matrix which will be modified by the triangularisation

allocate (ap(nn,mm), bp(nn))

do jj = 1, nn
  do ii = 1, mm
    ap(ii,jj) = aa(ii,jj)
  enddo
  bp(jj) = bb(jj)
enddo

! --- triangulise the matrix

ipass  = 0
ipivot = 0
ipermu = 0

kk = 1
do while (kk.lt.nn)

  ipass = ipass + 1
  if( ipass .gt. (nn) )then
    write(nfecra,*) "Gaussian elimnitation (gauss.f90):"
    write(nfecra,*) "ipass > n',ipass,'>',nn,'-> stop"
    write(nfecra,*) "verify the Gauss pivot"
    write(nfecra,*) "gauss diagnostic:"
    write(nfecra,*) "pivoting on ",ipivot,"cells"
    write(nfecra,*) "permuting on ",ipermu,"cells"
    call csexit(1)
  endif

  ! looking for the greatest coefficient
  ii = kk
  ll = 0
  akk = 0.d0
  do while (ii.le.nn)
    if (abs(ap(ii,kk)).gt.akk) then
      akk = abs(ap(ii,kk))
      ll = ii
    endif
    ii = ii + 1
  enddo

  if (ll.eq.0) then
    write(nfecra,*) &
         "Gaussian elimination (gauss.f90): no non zero pivot => stop"
    call csexit(1)
  endif

  ! pivoting line kk with line ll except if kk = ll
  if (ll.ne.kk) then
    do jj = 1, mm
      aij = ap(kk,jj)
      ap(kk,jj) = ap(ll,jj)
      ap(ll,jj) = aij
    enddo
    bi = bp(kk)
    bp(kk) = bp(ll)
    bp(ll) = bi
  endif

  ! scaling the pivot line
  akk = ap(kk,kk)
  do jj = kk, mm
    ap(kk,jj) = ap(kk,jj)/akk
  enddo
  bp(kk) = bp(kk)/akk

  ! elimination on the last lines
  do ii = kk+1, nn
    aik = ap(ii,kk)
    do jj = kk, mm
      ap(ii,jj) = ap(ii,jj) - aik*ap(kk,jj)
    enddo
    bp(ii) = bp(ii) - aik*bp(kk)
  enddo

  kk = kk + 1

enddo

! --- solves the triangular linear system

xx(nn) = bp(nn)/ap(nn,nn)
do ii = nn-1, 1, -1
  sum = 0.d0
  do jj = ii+1, mm
    sum = sum + ap(ii,jj)*xx(jj)
  enddo
  xx(ii) = (bp(ii)-sum)/ap(ii,ii)
enddo

! --- Computes residu
residu = 0.d0
do ii = 1, nn
  rr = bb(ii)
  do jj = 1, nn
    rr = rr - aa(ii,jj)*xx(jj)
  enddo
  residu = residu + rr**2
enddo

deallocate (ap, bp)

return

end subroutine gauss
