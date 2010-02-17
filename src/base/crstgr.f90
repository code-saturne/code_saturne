!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine crstgr &
!================

( iappel , isym   , igr    ,                                      &
  ncelf  , ncelg  , ncelfe , ncelge , nfacf  , nfacg ,            &
  iwarnp ,                                                        &
  ifaclf , ifaclg , irscel , irsfac ,                             &
  volumf , xyzfin , surfaf , xaf0   , xaf0ij ,                    &
  daf    , xaf    ,                                               &
  volumg , xyzgro , surfag , xag0   , xag0ij ,                    &
  dag    , xag    ,                                               &
  w1     , w2     , w3     , w4     )


!===============================================================================
! FONCTION :
! ----------

!  MULTIGRILLE ALGEBRIQUE :
!  CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER A PARTIR
!  DU NIVEAU SUPERIEUR

! structure grille fine ==> grille grossiere
! ------------------------------------------
! [ NCELF,NFACF,IFACLF,ROVDTF,XAF0,VOLUMF,XYZFIN,
!        SURFAF,XAF,XAF0,XAF0IJ ]
! [ NCELG,NFACG,IFACLG,ROVDTG,XAG0,VOLUMG,XYZGRO,
!        SURFAG,XAG,XAG0,XAG0IJ ]

!             explicité dans le cartouche


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iappel           ! e  ! <-- ! numero d'appel                                 !
! isym             ! e  ! <-- ! indicateur = 1 matrice sym                     !
!                  !    !     !            = 2 matrice non sym                 !
! igr              ! e  ! <-- ! niveau du maillage grossier                    !
! ncelf            ! e  ! <-- ! nombre d'elements maillage fin                 !
! ncelg            ! e  ! <-- ! nombre d'elements maillage grossier            !
! ncelfe           ! e  ! <-- ! nombre d'elements etendus fin                  !
! ncelge           ! e  ! <-- ! nombre d'elements etendus grossier             !
! nfacf            ! e  ! <-- ! nombre de faces internes maill. fin            !
! nfacg            ! e  ! <-- ! nombre de faces internes maill. gro.           !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! ifacef           ! te ! <-- ! elements voisins d'une face interne            !
!  (2, nfacf)      !    !     ! sur maillage fin                               !
! ifaceg           ! te ! <-- ! elements voisins d'une face interne            !
!  (2, nfacg)      !    !     ! sur maillage grossier                          !
! irscel           ! te ! <-- ! cellule fine -> cellule grossiere              !
!  (ncelfe)        !    !     !                                                !
! irsfac           ! te ! <-- ! face fine -> face grossiere                    !
!  (nfacf)         !    !     !  = 0 : face interne cel. grossiere             !
!                  !    !     !  < 0 : orientation inverse                     !
!                  !    !     !  > 0 : orientation identique                   !
! daf(ncelf)       ! tr ! <-- ! diagonale de la matrice mail fin               !
! xaf              ! tr ! <-- ! extradiagonale matrice maillage fin            !
!  (nfacf,isym)    !    !     !                                                !
! dag(ncelg)       ! tr ! --> ! diagonale matrice maillage grossier            !
! xag              ! tr ! --> ! extradiagonale matrice maillage                !
!  (nfacg,isym)    !    !     !  grossier                                      !
! ivois(ncelf)     ! te ! --- ! indicateur de voisinage cellules               !
! ip(ncelf)        ! te ! --- ! pointeurs sur voisins pour connecti-           !
!                  !    !     ! vite inverse                                   !
! icelfa           ! te ! --- ! connectivite cellules->faces mailla-           !
!  (2*nfacf)       !    !     ! ge fin                                         !
! icelce           ! te ! --- ! connectivite cellules->cellules                !
!  (2*nfacf)       !    !     ! voisines du maillage fin                       !
! rw(ncelf)        ! tr ! --- ! tableau de travail                             !
! ifaclf(2,nfacf   ! te ! <-- ! cell. voisines face intrn maill fin            !
! xaf0(nfacf,isym) ! tr ! <-- ! extradiagonale matrice p0 mailage fin          !
! volumf(ncelf)    ! tr ! <-- ! volume cellule maillage fin                    !
! xyzfin(3,ncelf)  ! tr ! <-- ! coordonnes cellule maillage fin                !
! surfaf(3,nfacf)  ! tr ! <-- ! surface face interne maillage fin              !
! xafxf0(2,nfacf)  ! tr ! <-- ! integ. xaf0*coord.cell adj. mail.fin           !
! ncelg            ! e  ! <-- ! nombre d'elements maillage grossier            !
! nfacg            ! e  ! <-- ! nombre faces internes maill. grossier          !
! ifaclg(2,nfacg)  ! te ! <-- ! cell. voisines face internes maillage grossier !
! xag0(nfacg,isym) ! tr ! <-- ! extradiagonale matrice p0 maillage grossier    !
! volumg(ncelg)    ! tr ! <-- ! volume cellule maillage grossier               !
! xyzgro(3,ncelg)  ! tr ! <-- ! coordonnes cellule maillage grossier           !
! surfag(3,nfacg)  ! tr ! <-- ! surface face interne maillage grossier         !
! xagxg0(2,nfacg)  ! tr ! <-- ! integ. xag0*coord.cell adj. maillage grossier  !
! w1,..,4(ncel)    ! tr ! <-> ! tableaux de travail                            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================


include "paramx.h"
include "entsor.h"
include "optcal.h"
include "cstnum.h"
include "parall.h"

!===============================================================================

! Arguments

integer          iappel, isym, igr
integer          ncelf, ncelfe, nfacf, ncelg, ncelge, nfacg
integer          iwarnp

double precision daf(ncelf), xaf(nfacf,2)
double precision dag(ncelge), xag(nfacg,2)

integer          ifaclf(2, nfacf), ifaclg(2, nfacg)
integer          irscel(ncelf), irsfac(nfacf)

double precision xaf0(nfacf),volumf(ncelfe), xyzfin(3, ncelfe)
double precision xag0(nfacg),volumg(ncelge), xyzgro(3, ncelge)
double precision surfaf(3, nfacf), xaf0ij(3, nfacf)
double precision surfag(3, nfacg), xag0ij(3, nfacg)
double precision w1(ncelfe), w2(ncelfe), w3(ncelfe), w4(ncelfe)

! Local variables

integer          iel, ii, jj, ig, jg
integer          ifacg, imin, imax
integer          ifac
integer          interp

double precision dsigjg, dsxaij
double precision rmin, rmax, cclip
double precision anmin(2), anmax(2)


!===============================================================================

! Compute volume and center of coarse cells: xyzgro(3,ncelg)
!==========================================================================

if (iappel .eq. 1) then

  do ig = 1, ncelge
    volumg(ig) = 0.d0
  enddo
  do ig = 1, ncelge
    xyzgro(1, ig) = 0.d0
    xyzgro(2, ig) = 0.d0
    xyzgro(3, ig) = 0.d0
  enddo

  do ii = 1, ncelf
    ig = irscel(ii)
    volumg(ig) = volumg(ig) + volumf(ii)
    xyzgro(1, ig) = xyzgro(1, ig) + volumf(ii)*xyzfin(1, ii)
    xyzgro(2, ig) = xyzgro(2, ig) + volumf(ii)*xyzfin(2, ii)
    xyzgro(3, ig) = xyzgro(3, ig) + volumf(ii)*xyzfin(3, ii)
  enddo
  do ig = 1, ncelg
    xyzgro(1, ig) = xyzgro(1, ig) / volumg(ig)
    xyzgro(2, ig) = xyzgro(2, ig) / volumg(ig)
    xyzgro(3, ig) = xyzgro(3, ig) / volumg(ig)
  enddo

! Return to calling function for parallel / periodic synchronization
! of xyzgro and volumg

  return

endif

! P0 restriction of matrixes, "internal" surface:
! xag0(nfacg), surfag(3,nfacgl), xagxg0(2,nfacg)
!==========================================================================

imax = 0

do ifacg = 1, nfacg
  xag0(ifacg) = 0.d0
  surfag(1, ifacg) = 0.d0
  surfag(2, ifacg) = 0.d0
  surfag(3, ifacg) = 0.d0
  xag0ij(1, ifacg) = 0.d0
  xag0ij(2, ifacg) = 0.d0
  xag0ij(3, ifacg) = 0.d0
enddo

do ifac = 1, nfacf

  if (irsfac(ifac).gt.0) then

    ifacg = irsfac(ifac)

    xag0(ifacg) = xag0(ifacg) + xaf0(ifac)

    surfag(1, ifacg) = surfag(1, ifacg) + surfaf(1, ifac)
    surfag(2, ifacg) = surfag(2, ifacg) + surfaf(2, ifac)
    surfag(3, ifacg) = surfag(3, ifacg) + surfaf(3, ifac)
    xag0ij(1, ifacg) = xag0ij(1, ifacg) + xaf0ij(1, ifac)
    xag0ij(2, ifacg) = xag0ij(2, ifacg) + xaf0ij(2, ifac)
    xag0ij(3, ifacg) = xag0ij(3, ifacg) + xaf0ij(3, ifac)

  else if (irsfac(ifac).lt.0) then

    ifacg = - irsfac(ifac)

    xag0(ifacg) = xag0(ifacg) + xaf0(ifac)

    surfag(1, ifacg) = surfag(1, ifacg) - surfaf(1, ifac)
    surfag(2, ifacg) = surfag(2, ifacg) - surfaf(2, ifac)
    surfag(3, ifacg) = surfag(3, ifacg) - surfaf(3, ifac)
    xag0ij(1, ifacg) = xag0ij(1, ifacg) - xaf0ij(1, ifac)
    xag0ij(2, ifacg) = xag0ij(2, ifacg) - xaf0ij(2, ifac)
    xag0ij(3, ifacg) = xag0ij(3, ifacg) - xaf0ij(3, ifac)

  endif

enddo


!===============================================================================
! Finalize computation of matrix in dag, xag
!===============================================================================

! interp = 0 : P0 restriction / P0 prolongation => XAG = XAG0
! interp = 1 : P0 restriction / P1 prolongation => XAG = XAG0ij/IgJg

! Initialization

interp = 1

! Initialize non differential fine mesh term saved in w1

do iel = 1, ncelf
  w1(iel) = daf(iel)
enddo
do iel = ncelf + 1, ncelfe
  w1(iel) = 0.d0
enddo
do ifac = 1, nfacf
  ii = ifaclf(1, ifac)
  jj = ifaclf(2, ifac)
  w1(ii) = w1(ii) + xaf(ifac, 1)
  w1(jj) = w1(jj) + xaf(ifac, isym)
enddo

! Initialize coarse matrix storage on (dag, xag)

do iel = 1, ncelge
  dag(iel) = 0.d0
enddo
do ifac = 1, nfacg
  xag(ifac, 1)= 0.d0
  xag(ifac, isym)= 0.d0
enddo

! Extradiagonal terms
! (symmetric matrixes for now, even with non symmetric storage isym=2)

! Matrix intialized to xag0 (interp=0)

do ifacg = 1, nfacg
  xag(ifacg, 1) = xag0(ifacg)
  xag(ifacg, isym) = xag0(ifacg)
enddo

if (interp.eq.1) then

  imin = 0
  imax = 0

  do ifacg = 1, nfacg

    ig = ifaclg(1, ifacg)
    jg = ifaclg(2, ifacg)

    dsigjg =   (xyzgro(1, jg)-xyzgro(1, ig))*surfag(1, ifacg)     &
             + (xyzgro(2, jg)-xyzgro(2, ig))*surfag(2, ifacg)     &
             + (xyzgro(3, jg)-xyzgro(3, ig))*surfag(3, ifacg)

    dsxaij =   xag0ij(1, ifacg)*surfag(1, ifacg)                  &
             + xag0ij(2, ifacg)*surfag(2, ifacg)                  &
             + xag0ij(3, ifacg)*surfag(3, ifacg)

    if (abs(dsigjg) .gt. epzero) then

      ! Standard
      xag(ifacg, 1)    = dsxaij/dsigjg
      xag(ifacg, isym) = dsxaij/dsigjg

      ! Clipped matrix
      cclip = dsxaij/dsigjg
      if (cclip .lt. xag0(ifacg)) imin = imin+1
      if (cclip .gt. 0.d0)  imax = imax +1
      if (cclip .lt. xag0(ifacg) .or. cclip .gt. 0.d0) then
        xag(ifacg, 1) =  xag0(ifacg)
        xag(ifacg, isym) = xag0(ifacg)
      endif

    endif

  enddo

  if(iwarnp.gt.3) then
    if (irangp .ge. 0) then
      call parcpt(imin)
      call parcpt(imax)
    endif
    write(nfecra, 2001) imin, imax
  endif

  ! Possible P1 matrix / P0 matrix relaxation defined by the user in usini1.f90

  do ifacg = 1, nfacg
    xag(ifacg, 1) = rlxp1*xag(ifacg, 1) +(1.d0-rlxp1)*xag0(ifacg)
    xag(ifacg, isym) = rlxp1*xag(ifacg, isym) +(1.d0-rlxp1)*xag0(ifacg)
  enddo

endif

if (interp.ne.0 .and. interp.ne.1) then

  write(nfecra,*) 'interp incorrectly defined in crstgr'
  write(nfecra,*) '--> Stop in crstgr '
  call csexit(1)

endif


! Diagonal term

do ii = 1, ncelf
  ig = irscel(ii)
  dag(ig) = dag(ig) + w1(ii)
enddo

do ifacg = 1, nfacg

  ig = ifaclg(1, ifacg)
  jg = ifaclg(2, ifacg)

  dag(ig) = dag(ig) - xag(ifacg, 1)
  dag(jg) = dag(jg) - xag(ifacg, isym)

enddo

! Check
!======

if (iwarnp .gt. 3) then

  ! Evaluate fine and coarse matrixes anisotropy

  do ii = 1, ncelfe
    w1(ii) =-1.d12
    w2(ii) =+1.d12
    w3(ii) =-1.d12
    w4(ii) =+1.d12
  enddo

  do ifac = 1, nfacf
    ii = ifaclf(1, ifac)
    jj = ifaclf(2, ifac)
    w1(ii) = max(abs(xaf(ifac, 1)), w1(ii))
    w2(ii) = min(abs(xaf(ifac, 1)), w2(ii))
    w1(jj) = max(abs(xaf(ifac, isym)), w1(jj))
    w2(jj) = min(abs(xaf(ifac, isym)), w2(jj))
  enddo

  do ifacg = 1, nfacg
    ig = ifaclg(1, ifacg)
    jg = ifaclg(2, ifacg)
    w3(ig) = max(abs(xag(ifacg,1)), w3(ig))
    w4(ig) = min(abs(xag(ifacg,1)), w4(ig))
    w3(jg) = max(abs(xag(ifacg, isym)), w3(jg))
    w4(jg) = min(abs(xag(ifacg, isym)), w4(jg))
  enddo

  do ii = 1, ncelf
    w1(ii) = w2(ii)/w1(ii)
  enddo

  do ig = 1, ncelg
    w3(ig) = w4(ig)/w3(ig)
  enddo

  anmin(1) = w1(1)
  anmax(1) = w1(1)
  do ii = 2, ncelf
    if (w1(ii) .lt. anmin(1)) then
      anmin(1) = w1(ii)
    else if (w1(ii) .gt. anmax(1)) then
      anmax(1) = w1(ii)
    endif
  enddo

  anmin(2) = w3(1)
  anmax(2) = w3(1)
  do ig = 2, ncelg
    if (w3(ig) .lt. anmin(2)) then
      anmin(2) = w3(ig)
    else if (w3(ig) .gt. anmax(2)) then
      anmax(2) = w3(ig)
    endif
  enddo

  if (irangp .ge. 0) then
    ii = 2
    call parrmn(ii, anmin)
    call parrmx(ii, anmax)
  endif

  write (nfecra, 2002) anmin(1), anmax(1), anmin(2), anmax(2)

  if (interp .eq. 1) then

    rmin = +1.d10
    rmax = -1.d10
    do ifacg=1,nfacg
      rmin = min(rmin, xag(ifacg,1)/xag0(ifacg))
      rmax = max(rmax, xag(ifacg,1)/xag0(ifacg))
    enddo

    if (irangp .ge. 0) then
      call parmin(rmin)
      call parmax(rmax)
    endif

    write(nfecra, 2003) rmin, rmax

  endif

  ! Evaluate fine and coarse matrixes diagonal dominance

  do ii = 1, ncelf
    w1(ii) = abs(daf(ii))
  enddo
  do ii = ncelf+1, ncelfe
    w1(ii) = 0.d0
  enddo
  do ig = 1, ncelg
    w3(ig) = abs(dag(ig))
  enddo
  do ig = ncelg+1, ncelge
    w3(ig) = 0.d0
  enddo

  do ifac = 1, nfacf
    ii = ifaclf(1, ifac)
    jj = ifaclf(2, ifac)
    w1(ii) = w1(ii) - abs(xaf(ifac, 1))
    w1(jj) = w1(jj) - abs(xaf(ifac, isym))
  enddo

  do ifacg = 1, nfacg
    ig = ifaclg(1, ifacg)
    jg = ifaclg(2, ifacg)
    w3(ig) = w3(ig) - abs(xag(ifacg, 1))
    w3(jg) = w3(jg) - abs(xag(ifacg, isym))
  enddo

  do ii = 1, ncelf
    w1(ii) = w1(ii) / abs(daf(ii))
  enddo
  do ig = 1, ncelg
    w3(ig) = w3(ig) / abs(dag(ig))
  enddo

  anmin(1) = w1(1)
  anmax(1) = w1(1)
  do ii = 2, ncelf
    if (w1(ii) .lt. anmin(1)) then
      anmin(1) = w1(ii)
    else if (w1(ii) .gt. anmax(1)) then
      anmax(1) = w1(ii)
    endif
  enddo

  anmin(2) = w3(1)
  anmax(2) = w3(1)
  do ig = 2, ncelg
    if (w3(ig) .lt. anmin(2)) then
      anmin(2) = w3(ig)
    else if (w3(ig) .gt. anmax(2)) then
      anmax(2) = w3(ig)
    endif
  enddo

  if (irangp .ge. 0) then
    ii = 2
    call parrmn(ii, anmin)
    call parrmx(ii, anmax)
  endif

  write (nfecra, 2004) anmin(1), anmax(1), anmin(2), anmax(2)

endif

!--------
! Formats
!--------

 2001 format(&
  '    crstgr: coarse matrix < xag0 at ', i10,' faces', /, &
  '                          > 0    at ', i10,' faces')
 2002 format(&
  '       fine mesh anisotropy:   min = ', e12.5, /, &
  '                               max = ', e12.5, /, &
  '       coarse mesh anisotropy: min = ', e12.5, /, &
  '                               max = ', e12.5, /)
 2003 format(&
  '       minimum xag_p1/xag_p0 = ', e12.5, /, &
  '       maximum xag_p1/xag_p0 = ', e12.5, /)

 2004 format(&
  '       fine mesh diag dominance:   min = ', e12.5, /, &
  '                                   max = ', e12.5, /, &
  '       coarse mesh diag dominance: min = ', e12.5, /, &
  '                                   max = ', e12.5, /)
!----
! End
!----

return
end subroutine
