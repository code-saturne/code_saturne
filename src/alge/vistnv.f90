!-------------------------------------------------------------------------------

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

!===============================================================================
! Function:
! ---------

!> \file vistnv.f90
!>
!> \brief This function computes the equivalent tensor viscosity at faces for
!> a 3x3 symetric tensor.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     imvisf        method to compute the viscosity at faces:
!>                               - 0: arithmetic
!>                               - 1: harmonic
!> \param[in]     w1            cell viscosity symmetric tensor
!> \param[out]    viscf         inner face tensor viscosity
!>                               (times surface divided distance)
!> \param[out]    viscb         inner face viscosity
!>                               (surface, must be consistent with flux BCs)
!_______________________________________________________________________________

subroutine vistnv &
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use optcal, only: iporos
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          imvisf

double precision, target :: w1(6,ncelet)
double precision viscf(3,3,nfac), viscb(nfabor)

! Local variables

integer          ifac, iel, ii, jj, isou, jsou
double precision visci(3,3), viscj(3,3), s1(6), s2(6)
double precision pnd, srfddi

double precision, pointer, dimension(:,:) :: viscce => null()
double precision, dimension(:,:), allocatable, target :: w2

!===============================================================================

! ---> Periodicity and parallelism treatment

! Without porosity
if (iporos.eq.0) then
  viscce => w1(:,:)

! With porosity
else if (iporos.eq.1) then
  allocate(w2(6, ncelet))
  do iel = 1, ncel
    do isou = 1, 6
      w2(isou, iel) = porosi(iel)*w1(isou, iel)
    enddo
  enddo
  viscce => w2(:,:)

! With tensorial porosity
else if (iporos.eq.2) then
  allocate(w2(6,ncelet))
  do iel = 1, ncel
    call symmetric_matrix_product(w2(1, iel), porosf(1, iel), w1(1, iel))
  enddo
  viscce => w2(:,:)
endif

! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call syntis(viscce)
endif

! Arithmetic mean
if (imvisf.eq.0) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    visci(1,1) = viscce(1,ii)
    visci(2,2) = viscce(2,ii)
    visci(3,3) = viscce(3,ii)
    visci(1,2) = viscce(4,ii)
    visci(2,1) = viscce(4,ii)
    visci(2,3) = viscce(5,ii)
    visci(3,2) = viscce(5,ii)
    visci(1,3) = viscce(6,ii)
    visci(3,1) = viscce(6,ii)

    viscj(1,1) = viscce(1,jj)
    viscj(2,2) = viscce(2,jj)
    viscj(3,3) = viscce(3,jj)
    viscj(1,2) = viscce(4,jj)
    viscj(2,1) = viscce(4,jj)
    viscj(2,3) = viscce(5,jj)
    viscj(3,2) = viscce(5,jj)
    viscj(1,3) = viscce(6,jj)
    viscj(3,1) = viscce(6,jj)

    do isou = 1, 3
      do jsou = 1, 3
        viscf(isou,jsou,ifac) = 0.5d0*(visci(isou,jsou)+viscj(isou,jsou)) &
                              * surfan(ifac)/dist(ifac)
      enddo
    enddo

  enddo

! Harmonic mean: Kf = Ki . (pnd Ki +(1-pnd) Kj)^-1 . Kj
else

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd  = pond(ifac)

    do isou = 1, 6
      s1(isou) = pnd*viscce(isou, ii) + (1.d0-pnd)*viscce(isou, jj)
    enddo

    call symmetric_matrix_inverse(s2, s1)

    call symmetric_matrix_product(s1, s2, viscce(:, jj))

    call symmetric_matrix_product(s2, viscce(:, ii), s1)

    srfddi = surfan(ifac)/dist(ifac)

    viscf(1,1,ifac) = s2(1)*srfddi
    viscf(2,2,ifac) = s2(2)*srfddi
    viscf(3,3,ifac) = s2(3)*srfddi
    viscf(1,2,ifac) = s2(4)*srfddi
    viscf(2,1,ifac) = s2(4)*srfddi
    viscf(2,3,ifac) = s2(5)*srfddi
    viscf(3,2,ifac) = s2(5)*srfddi
    viscf(1,3,ifac) = s2(6)*srfddi
    viscf(3,1,ifac) = s2(6)*srfddi

  enddo

endif

! Without porosity
if (iporos.eq.0) then

  do ifac = 1, nfabor
    ii = ifabor(ifac)
    viscb(ifac) = surfbn(ifac)
  enddo

! With porosity
else if (iporos.eq.1) then

  do ifac = 1, nfabor
    ii = ifabor(ifac)
    viscb(ifac) = surfbn(ifac)*porosi(ii)
  enddo

! With anisotropic porosity
else if (iporos.eq.2) then

  do ifac = 1, nfabor
    ii = ifabor(ifac)
    viscb(ifac) = surfbn(ifac)*porosi(ii)
  enddo

endif

if (allocated(w2)) deallocate(w2)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    sout          sout = s1 * s2
!> \param[in]     s             symmetric matrix
!_______________________________________________________________________________

subroutine symmetric_matrix_inverse &
 ( sout , s )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision sout(6), s(6)

! Local variables

double precision detinv

!===============================================================================

sout(1) = s(2)*s(3) - s(5)**2
sout(2) = s(1)*s(3) - s(6)**2
sout(3) = s(1)*s(2) - s(4)**2
sout(4) = s(5)*s(6) - s(4)*s(3)
sout(5) = s(4)*s(6) - s(1)*s(5)
sout(6) = s(4)*s(5) - s(2)*s(6)

detinv = 1.d0 / (s(1)*sout(1) + s(4)*sout(4) + s(6)*sout(6))

sout(1) = sout(1) * detinv
sout(2) = sout(2) * detinv
sout(3) = sout(3) * detinv
sout(4) = sout(4) * detinv
sout(5) = sout(5) * detinv
sout(6) = sout(6) * detinv

return
end subroutine
