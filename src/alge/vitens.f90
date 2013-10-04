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

!> \file vitens.f90
!>
!> \brief This function computes the equivalent viscosity at faces for
!> a 3x3 symetric tensor.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     w1            cell viscosity symmetric tensor
!> \param[in]     iwarnp        verbosity
!> \param[out]    weighf        inner face weight between cells i and j
!>                              \f$ \frac{\vect{IF} \cdot \tens{K}_\celli}
!>                                       {\norm{\tens{K}_\celli \cdot \vect{S}}^2} \f$
!>                              and
!>                              \f$ \frac{\vect{FJ} \cdot \tens{K}_\cellj}
!>                                       {\norm{\tens{K}_\cellj \cdot \vect{S}}^2} \f$
!> \param[out]    weighb        boundary face weight
!>                              \f$ \frac{\vect{IF} \cdot \tens{K}_\celli}
!>                                       {\norm{\tens{K}_\celli \cdot \vect{S}}^2} \f$
!> \param[out]    viscf         inner face viscosity
!>                               (times surface divided distance)
!> \param[out]    viscb         inner face viscosity
!>                               (surface, must be consistent with flux BCs)
!_______________________________________________________________________________

subroutine vitens &
 ( w1     , iwarnp,                                               &
   weighf , weighb,                                               &
   viscf  , viscb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use optcal, only: iporos
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          iwarnp

double precision, target :: w1(6,ncelet)
double precision weighf(2,nfac), weighb(nfabor)
double precision viscf(nfac), viscb(nfabor)

! Local variables

integer          ifac, iel, ii, jj, isou
integer          nclipf, nclipb
double precision visci(3,3), viscj(3,3)
double precision viscis, viscjs, fikis, fjkjs, distfi, distfj
double precision temp, eps

double precision, pointer, dimension(:,:) :: viscce => null()
double precision, dimension(:,:), allocatable, target :: w2

!===============================================================================

nclipf = 0
nclipb = 0

eps = 1.d-1

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

! Always Harmonic mean

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

  ! ||Ki.S||^2
  viscis = ( visci(1,1)*surfac(1,ifac)       &
           + visci(1,2)*surfac(2,ifac)       &
           + visci(1,3)*surfac(3,ifac))**2   &
         + ( visci(2,1)*surfac(1,ifac)       &
           + visci(2,2)*surfac(2,ifac)       &
           + visci(2,3)*surfac(3,ifac))**2   &
         + ( visci(3,1)*surfac(1,ifac)       &
           + visci(3,2)*surfac(2,ifac)       &
           + visci(3,3)*surfac(3,ifac))**2

  ! IF.Ki.S
  fikis = ( (cdgfac(1,ifac)-xyzcen(1,ii))*visci(1,1)   &
          + (cdgfac(2,ifac)-xyzcen(2,ii))*visci(2,1)   &
          + (cdgfac(3,ifac)-xyzcen(3,ii))*visci(3,1)   &
          )*surfac(1,ifac)                             &
        + ( (cdgfac(1,ifac)-xyzcen(1,ii))*visci(1,2)   &
          + (cdgfac(2,ifac)-xyzcen(2,ii))*visci(2,2)   &
          + (cdgfac(3,ifac)-xyzcen(3,ii))*visci(3,2)   &
          )*surfac(2,ifac)                             &
        + ( (cdgfac(1,ifac)-xyzcen(1,ii))*visci(1,3)   &
          + (cdgfac(2,ifac)-xyzcen(2,ii))*visci(2,3)   &
          + (cdgfac(3,ifac)-xyzcen(3,ii))*visci(3,3)   &
          )*surfac(3,ifac)

  distfi = (1.d0 - pond(ifac))*dist(ifac)

  ! Take I" so that I"F= eps*||FI||*Ki.n when I" is in cell rji
  temp = eps*sqrt(viscis)*distfi
  if (fikis.lt.temp) then
    fikis = temp
    nclipf = nclipf + 1
  endif

  viscj(1,1) = viscce(1,jj)
  viscj(2,2) = viscce(2,jj)
  viscj(3,3) = viscce(3,jj)
  viscj(1,2) = viscce(4,jj)
  viscj(2,1) = viscce(4,jj)
  viscj(2,3) = viscce(5,jj)
  viscj(3,2) = viscce(5,jj)
  viscj(1,3) = viscce(6,jj)
  viscj(3,1) = viscce(6,jj)

  ! ||Kj.S||^2
  viscjs = ( viscj(1,1)*surfac(1,ifac)       &
           + viscj(1,2)*surfac(2,ifac)       &
           + viscj(1,3)*surfac(3,ifac))**2   &
         + ( viscj(2,1)*surfac(1,ifac)       &
           + viscj(2,2)*surfac(2,ifac)       &
           + viscj(2,3)*surfac(3,ifac))**2   &
         + ( viscj(3,1)*surfac(1,ifac)       &
           + viscj(3,2)*surfac(2,ifac)       &
           + viscj(3,3)*surfac(3,ifac))**2

  ! FJ.Kj.S
  fjkjs = ( (xyzcen(1,jj)-cdgfac(1,ifac))*viscj(1,1)   &
          + (xyzcen(2,jj)-cdgfac(2,ifac))*viscj(2,1)   &
          + (xyzcen(3,jj)-cdgfac(3,ifac))*viscj(3,1)   &
          )*surfac(1,ifac)                             &
        + ( (xyzcen(1,jj)-cdgfac(1,ifac))*viscj(1,2)   &
          + (xyzcen(2,jj)-cdgfac(2,ifac))*viscj(2,2)   &
          + (xyzcen(3,jj)-cdgfac(3,ifac))*viscj(3,2)   &
          )*surfac(2,ifac)                             &
        + ( (xyzcen(1,jj)-cdgfac(1,ifac))*viscj(1,3)   &
          + (xyzcen(2,jj)-cdgfac(2,ifac))*viscj(2,3)   &
          + (xyzcen(3,jj)-cdgfac(3,ifac))*viscj(3,3)   &
          )*surfac(3,ifac)

  distfj = pond(ifac)*dist(ifac)

  ! Take J" so that FJ"= eps*||FJ||*Kj.n when J" is in cell i
  temp = eps*sqrt(viscjs)*distfj
  if (fjkjs.lt.temp) then
    fjkjs = temp
    nclipf = nclipf + 1
  endif

  weighf(1,ifac) = fikis/viscis
  weighf(2,ifac) = fjkjs/viscjs

  viscf(ifac) = 1.d0/(weighf(1,ifac) + weighf(2,ifac))

enddo

do ifac = 1, nfabor

  ii = ifabor(ifac)

  visci(1,1) = viscce(1,ii)
  visci(2,2) = viscce(2,ii)
  visci(3,3) = viscce(3,ii)
  visci(1,2) = viscce(4,ii)
  visci(2,1) = viscce(4,ii)
  visci(2,3) = viscce(5,ii)
  visci(3,2) = viscce(5,ii)
  visci(1,3) = viscce(6,ii)
  visci(3,1) = viscce(6,ii)

  ! ||Ki.S||^2
  viscis = ( visci(1,1)*surfbo(1,ifac)       &
           + visci(1,2)*surfbo(2,ifac)       &
           + visci(1,3)*surfbo(3,ifac))**2   &
         + ( visci(2,1)*surfbo(1,ifac)       &
           + visci(2,2)*surfbo(2,ifac)       &
           + visci(2,3)*surfbo(3,ifac))**2   &
         + ( visci(3,1)*surfbo(1,ifac)       &
           + visci(3,2)*surfbo(2,ifac)       &
           + visci(3,3)*surfbo(3,ifac))**2

  ! IF.Ki.S
  fikis = ( (cdgfbo(1,ifac)-xyzcen(1,ii))*visci(1,1)   &
          + (cdgfbo(2,ifac)-xyzcen(2,ii))*visci(2,1)   &
          + (cdgfbo(3,ifac)-xyzcen(3,ii))*visci(3,1)   &
          )*surfbo(1,ifac)                             &
        + ( (cdgfbo(1,ifac)-xyzcen(1,ii))*visci(1,2)   &
          + (cdgfbo(2,ifac)-xyzcen(2,ii))*visci(2,2)   &
          + (cdgfbo(3,ifac)-xyzcen(3,ii))*visci(3,2)   &
          )*surfbo(2,ifac)                             &
        + ( (cdgfbo(1,ifac)-xyzcen(1,ii))*visci(1,3)   &
          + (cdgfbo(2,ifac)-xyzcen(2,ii))*visci(2,3)   &
          + (cdgfbo(3,ifac)-xyzcen(3,ii))*visci(3,3)   &
          )*surfbo(3,ifac)

  distfi = distb(ifac)

  ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
  temp = eps*sqrt(viscis)*distfi
  if (fikis.lt.temp) then
    fikis = temp
    nclipb = nclipb + 1
  endif

  weighb(ifac) = fikis/viscis

enddo

! Without porosity
if (iporos.eq.0) then

  do ifac = 1, nfabor

    ! Warning: hint must be ||Ki.n||/I"F
    viscb(ifac) = surfbn(ifac)
  enddo

! With porosity
else if (iporos.eq.1) then

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    ! Warning: hint must be ||Ki.n||/I"F
    viscb(ifac) = surfbn(ifac)*porosi(ii)

  enddo

! With tensorial porosity
else if (iporos.eq.2) then

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    ! Warning: hint must be ||Ki.n||/I"F
    viscb(ifac) = surfbn(ifac)*porosi(ii)

  enddo

endif

if (irangp.ge.0) then
  call parsom(nclipf)
  call parsom(nclipb)
endif

if (iwarnp.ge.3) then
  write(nfecra,1000) nclipf, nclipb
endif

if (allocated(w2)) deallocate(w2)

!--------
! Formats
!--------

 1000 format ( &
 'Computing the face viscosity from the tensorial viscosity:',/,&
 '   Number of internal clippings: ',I5                      ,/,&
 '   Number of boundary clippings: ',I5)

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
!> \param[in]     s1            symmetric matrix
!> \param[in]     s2            symmetric matrix
!_______________________________________________________________________________

subroutine symmetric_matrix_product &
 ( sout , s1, s2 )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision sout(6), s1(6), s2(6)

! Local variables

!===============================================================================

! S11
sout(1) = s1(1)*s2(1) + s1(4)*s2(4) + s1(6)*s2(6)
! S22
sout(2) = s1(4)*s2(4) + s1(2)*s2(2) + s1(5)*s2(5)
! S33
sout(3) = s1(6)*s2(6) + s1(5)*s2(5) + s1(3)*s2(3)
! S12 = S21
sout(4) = s1(1)*s2(4) + s1(4)*s2(2) + s1(6)*s2(5)
! S23 = S32
sout(5) = s1(4)*s2(6) + s1(2)*s2(5) + s1(5)*s2(3)
! S13 = S31
sout(6) = s1(1)*s2(6) + s1(4)*s2(5) + s1(6)*s2(3)

return
end subroutine
