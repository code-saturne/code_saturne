!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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
!    nom           !type!mode !                   role                         !
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
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
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
!ifaclf(2,nfacf    ! te ! <-- ! cell. voisines face intrn maill fin            !
!xaf0(nfacf,isym tr ! <-- ! extradiagonale matrice p0 mailage fin          !
! volumf(ncelf)    ! tr ! <-- ! volume cellule maillage fin                    !
!xyzfin(3,ncelf) tr ! <-- ! coordonnes cellule maillage fin                !
!surfaf(3,nfacf) tr ! <-- ! surface face interne maillage fin              !
!xafxf0(2,nfacf) tr ! <-- ! integ. xaf0*coord.cell adj. mail.fin           !
! ncelg            ! e  ! <-- ! nombre d'elements maillage grossier            !
! nfacg            ! e  !  <- ! nombre faces internes maill. grossier          !
!ifaclg(2,nfacg    ! te !  <- ! cell. voisines face intrn mail. gros           !
!xag0(nfacg,isym tr !  <- ! extradiagonale matrice p0 maill.gros.          !
! volumg(ncelg)    ! tr !  <- ! volume cellule maillage grossier               !
!xyzgro(3,ncelg) tr !  <- ! coordonnes cellule maillage grossier           !
!surfag(3,nfacg) tr !  <- ! surface face interne maill. grossier           !
!xagxg0(2,nfacg) tr !  <- ! integ. xag0*coord.cell adj. mail.gro           !
!argu !W1,..,4 (NCEL)! TR ! <->           ! TABLEAUX DE TRAVAILS
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
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

! VARIABLES LOCALES

integer          iel, ii, jj, ig, jg
integer          ifacg, imin, imax
integer          ifac
integer          interp

double precision dsigjg, dsxaij
double precision rmin, rmax, cclip
double precision anmin(2), anmax(2)


!===============================================================================

! CREATION VOLUME, CdG CELLULES GROSSIERES : XYZGRO(3,NCELG)
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

!       RETOUR A LA FONCTION APPELANTE POUR SYNCHRONISATION
!       PARALLELE / PERIODIQUE DE XYZGRO ET VOLUMG

  return

endif

! RESTRICTION P0 MATRICES, SURFACE "INTERNE" :
!       XAG0(NFACG), SURFAG(3,NFACG), XAGXG0(2,NFACG)
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
! FINALISATION CALCUL DE LA MATRICE  DANS DAG, XAG
!===============================================================================

! INTERP= 0 : RESTRICTION P0 /PROLONGATION P0 => XAG=XAG0
! INTERP= 1 : RESTRICTION P0 /PROLONGATION P1 => XAG=XAG0IJ/IgJg

!     INITIALISATION

interp = 0
interp = 1


!     INITIALISATION TERME NON DIFFERENTIEL MESH FIN STOCKE DANS W1
!==============================

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


!     INITIALISATION STOCKAGE MATRICE MESH GROSSIER SUR (DAG, XAG)
!=============================

do iel = 1, ncelge
  dag(iel) = 0.d0
enddo
do ifac = 1, nfacg
  xag(ifac, 1)= 0.d0
  xag(ifac, isym)= 0.d0
enddo


!     TERMES EXTRADIAGONAUX
!     (matrices sym. pour l'instant, meme si stockage non syme isym=2)
!========================================

!     MATRICE INTIALISEE A XAG0 (INTERP=0)

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

!           STANDARD
      xag(ifacg, 1)    = dsxaij/dsigjg
      xag(ifacg, isym) = dsxaij/dsigjg

!           MATRICE CLIPPEE
      cclip= dsxaij/dsigjg
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
    write(nfecra, *)                                              &
      '    crstgr.F : matrice grossiere < XAG0 en ',IMIN,' faces'
    write(nfecra, *)                                              &
      '    crstgr.F : matrice grossiere > 0    en ',IMAX,' faces'
  endif

!       RELAXATION EVENTUELLE MATRICE P1 / MATRICE P0 DEFINI PAR UTILISATEUR
!       DANS usini1.F

  do ifacg = 1, nfacg
    xag(ifacg, 1)                                                 &
      = rlxp1*xag(ifacg, 1) +(1.d0-rlxp1)*xag0(ifacg)
    xag(ifacg, isym)                                              &
      = rlxp1*xag(ifacg, isym) +(1.d0-rlxp1)*xag0(ifacg)
  enddo

endif

if (interp.ne.0 .and. interp.ne.1) then

  WRITE(NFECRA,*) 'INTERP MAL DEFINI DANS crstgr.F'
  WRITE(NFECRA,*) '--> ARRET DANS crstgr.F '
  call csexit(1)

endif


!     TERME DIAGONAL
!============================

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

!     CONTROLE
!=============
!     WRITE(NFECRA,*) 'TYPE INTERPOLATION MATRICE = ',INTERP

!     EVALUATION ANISOTROPIE DES MATRICES FINE ET GROSSIERE

if (iwarnp .gt. 3) then

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

  if (iwarnp .gt. 3) then

    write (nfecra, 2000) anmin(1), anmax(1), anmin(2), anmax(2)

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

      WRITE(NFECRA, *) '      minimum XAG_P1/XAG_P0 = ', RMIN
      WRITE(NFECRA, *) '      maximum XAG_P1/XAG_P0 = ', RMAX

    endif

  endif

endif

!--------
! FORMATS
!--------

 2000 format('       anisotropie maillage fin :      min = ', E12.5, /, &
       '                                       max = ', E12.5, /, &
       '       anisotropie maillage grossier : min = ', E12.5, /, &
       '                                       max = ', E12.5, /)

!----
! FIN
!----

return
end
