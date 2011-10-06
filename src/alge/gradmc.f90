!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine gradmc &
!================

 ( ncelet , ncel   , nfac   , nfabor , ncelbr ,                   &
   inc    , iccocg , nswrgp , idimte , itenso , iphydp , imrgra , &
   iwarnp , nfecra , epsrgp , extrap ,                            &
   ifacel , ifabor , icelbr , ipcvse , ielvse , isympa ,          &
   volume , surfac , surfbo , surfbn , pond   ,                   &
   dist   , distbr , dijpf  , diipb  ,                            &
   fextx  , fexty  , fextz  ,                                     &
   xyzcen , cdgfac , cdgfbo , coefap , coefbp , pvar   ,          &
   cocgb  , cocg   ,                                              &
   dpdx   , dpdy   , dpdz   ,                                     &
   bx     , by     , bz     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DU GRADIENT CELLULE PAR RECONSTRUCTION MOINDRES CARRES 99
!  SUR UN SUPPORT ETENDU
! AVEC PRISE EN COMPTE EVENTUELLE D'UN TERME DE FORCE VOLUMIQUE
!  GENERANT UNE COMPOSANTE DE PRESSION HYDROSTATIQUE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au  moins              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! idimte           ! e  ! <-- ! dimension de la varible (maximum 3)            !
!                  !    !     ! 0 : scalaire (var11) ou assimile               !
!                  !    !     !     scalaire                                   !
!                  !    !     ! 1 : vecteur (var11,var22,var33)                !
!                  !    !     ! 2 : tenseur d'ordre 2 (varij)                  !
! itenso           ! e  ! <-- ! pour l'explicitation de la rotation            !
!                  !    !     ! 0 : scalaire (var11)                           !
!                  !    !     ! 1 : composante de vecteur ou de                !
!                  !    !     !     tenseur (var11) implcite pour la           !
!                  !    !     !     translation                                !
!                  !    !     !11 : reprend le traitement itenso=1 et          !
!                  !    !     !     composante de vecteur ou de                !
!                  !    !     !     tenseur (var11) annulee  pour la           !
!                  !    !     !     rotation                                   !
!                  !    !     ! 2 : vecteur (var11 et var22 et var33)          !
!                  !    !     !     implicite pour la translation              !
! iphydp           ! e  ! <-- ! indicateur de prise en compte de la            !
!                  !    !     ! pression hydrostatique                         !
! imrgra           ! e  ! <-- ! methode de reconstruction du gradient          !
!                  !    !     !  0 reconstruction 97                           !
!                  !    !     !  1 moindre carre                               !
!                  !    !     !  2 moindre carre support etendu                !
!                  !    !     !  3 moindre carre support etendu redui          !
!                  !    !     !  4 reconstr avec init moindres carres          !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! ifacel(2,nfac    ! te ! <-- ! no des elts voisins d'une face intern          !
! ifabor(nfabor    ! te ! <-- ! no de l'elt voisin d'une face de bord          !
! icelbr           ! te ! <-- ! numero global des elements ayant au            !
! (ncelbr)         !    !     !  moins une face de bord                        !
!                  !    !     !  ranges par numero croissant                   !
! ipcvse           ! te ! <-- ! position dans ielvse des voisins               !
! (ncel  )         !    !     !  etendus des cellules                          !
! ielvse (*)       ! te ! <-- ! numero des voisins etendus des cellul          !
! isympa(nfabor    ! te ! <-- ! nul sur les symetries                          !
! volume(ncelet    ! tr ! <-- ! volume des elements                            !
! surfac(3,nfac    ! tr ! <-- ! surf vectorielle des surfaces interne          !
! surfbo           ! tr ! <-- ! surf vectorielle des surfaces de bord          !
!   (3,nfabor)     !    !     !                                                !
! surfbn           ! tr ! <-- ! norme de la surface des faces de bord          !
! (nfabor)         !    !     !                                                !
! pond(nfac)       ! tr ! <-- ! ponderation geometrique (entre 0 et 1          !
! dist(nfac)       ! tr ! <-- ! dist entre les projections orthogonal          !
!                  !    !     !  sur la normale a une face des centre          !
!                  !    !     !  volumes voisins                               !
! distbr(nfabor    ! tr ! <-- ! dist du centre a la face de bord               !
! dijpf(3,nfac)    ! tr ! <-- ! vect i'j', i' (resp. j') projection            !
!                  !    !     !  du centre i (resp. j) sur la normale          !
!                  !    !     !  a la face interne                             !
! diipb            ! tr ! <-- ! vect ii', ii projection du centre i            !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
!  (3,ncelet       !    !     !                                                !
! cdgfac           ! tr ! <-- ! point associes aux facettes fluides            !
!  (3,nfac)        !    !     !                                                !
! cdgfbo           ! tr ! <-- ! points associes aux facettes de bord           !
!  (3,nfabor)      !    !     !                                                !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! pvar  (ncelet    ! tr ! <-- ! variable (ncelet + v. etendu eventuel          !
! fextx,y,z        ! tr ! <-- ! force exterieure generant la pression          !
!   (ncelet)       !    !     !  hydrostatique                                 !
! cocgb            ! tr ! <-- ! contribution des faces internes a              !
!  (ncelbr,3,3)    !    !     !  cocg (cellules de bord)                       !
! cocg             ! tr ! <-- ! couplage des composantes du gradient           !
!  ncelet,3,3      !    !     ! modifie eventuellement aux bords               !
! dpdx dpdy        ! tr ! --> ! gradient de pvar                               !
! dpdz (ncelet     ! tr !     !                                                !
! bx,y,z(ncelet    ! tr ! --- ! tableau de travail pour le grad de p           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use albase
use cplsat
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfac   , nfabor , ncelbr
integer          inc    , iccocg , nswrgp
integer          idimte , itenso , iphydp , imrgra
integer          iwarnp , nfecra
double precision epsrgp , extrap

integer          ifacel(2,nfac),ifabor(nfabor),icelbr(ncelbr)
integer          ipcvse(*), ielvse(*), isympa(nfabor)
double precision volume(ncelet), surfac(3,nfac)
double precision surfbo(3,nfabor), surfbn(nfabor)
double precision pond(nfac), dist(nfac), distbr(nfabor)
double precision dijpf(3,nfac), diipb(3,nfabor)
double precision xyzcen(3,*),cdgfac(3,nfac),cdgfbo(3,nfabor)
double precision coefap(nfabor), coefbp(nfabor), pvar(*)
double precision fextx(*),fexty(*),fextz(*)
double precision cocgb(ncelbr,3,3), cocg(ncelet,3,3)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision bx(ncelet),by(ncelet),bz(ncelet)

! Local variables

integer          lbloc
parameter       (lbloc = 1024)

integer          ii    , jj    , iel   , ielb  , ifac  , ipos
integer          ibloc , nbloc , irel  , iiii  , invemx
double precision pfac  , rkij  , dsij(3), aa(lbloc,3,3)
double precision a11   , a22   , a33   , a12   , a13   , a23
double precision cocg11, cocg12, cocg13, cocg21, cocg22, cocg23
double precision cocg31, cocg32, cocg33
double precision usdij2, pfsx  , pfsy  , pfsz
double precision unsdij, unssbn, umcbsd, unsdet
double precision pfac1 , pfac2 , pfac3 , vecfac, unsvol, extrab
double precision ptmidx, ptmidy, ptmidz

! INDICATEUR DE PREMIER PASSAGE
integer          inicoc
data             inicoc /1/
save             inicoc

! VERIFICATION QUE NCEL > NBRE DE VOISINS ETENDUS D'UNE CELLULE
integer          icesve
data             icesve /-1/
save             icesve

!===============================================================================

!===============================================================================
! 0. TRAITEMENT DE EXTRAG POUR LES SYMETRIES
!===============================================================================

! POUR LES CALCULS 2D EN PARTICULIER, SI ON EXTRAPOLE LE GRADIENT DE
!   PRESSION, ON SE TROUVE AVEC UNE MATRICE COCG NON INVERSIBLE, A
!   CAUSE DE LA TROISIEME DIRECTION.

! COMME L'ON N'A EN PRATIQUE QU'UNE SEULE PHASE,
!   ET QUE EXTRAG = 0 ou 1 (1 POUR LA PRESSION UNIQUEMENT)
!   ON SE PLACE IMPLICITEMENT DANS CETTE SITUATION, AVEC UN STOP
!   DANS VERINI AU CAS OU CE NE SERAIT PAS LE CAS.

! LA MODIFICATION CONSISTE A MULTIPLIER EXTRAP PAR ISYMPA QUI EST
!   NUL SUR LES SYMETRIES : ON N'EXTRAPOLE DONC PAS LE GRADIENT SUR CES
!   FACES.

!===============================================================================
! 1. INITIALISATION DU GRADIENT
!===============================================================================

if( nswrgp.le.1 ) then

  do iel = 1, ncelet
    bx  (iel) = 0.d0
    by  (iel) = 0.d0
    bz  (iel) = 0.d0
  enddo

!  CAS STANDARD, SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  ===============================================================
  if (iphydp.eq.0) then

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

    do ifac = 1,nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      pfac = pond(ifac)*pvar(ii) +(1.d0-pond(ifac))*pvar(jj)
      pfac1 = pfac*surfac(1,ifac)
      pfac2 = pfac*surfac(2,ifac)
      pfac3 = pfac*surfac(3,ifac)
      bx(ii) = bx(ii) +pfac1
      by(ii) = by(ii) +pfac2
      bz(ii) = bz(ii) +pfac3
      bx(jj) = bx(jj) -pfac1
      by(jj) = by(jj) -pfac2
      bz(jj) = bz(jj) -pfac3
    enddo

!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

    do ifac = 1,nfabor
      ii = ifabor(ifac)
      pfac = inc*coefap(ifac) +coefbp(ifac)*pvar(ii)
      bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
      by(ii) = by(ii) +pfac*surfbo(2,ifac)
      bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
    enddo

!  CAS AVEC PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  =====================================================
  else

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

    do ifac = 1,nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      pfac   = pond(ifac)*(pvar(ii)                             &
           -(xyzcen(1,ii)-cdgfac(1,ifac))*fextx(ii)             &
           -(xyzcen(2,ii)-cdgfac(2,ifac))*fexty(ii)             &
           -(xyzcen(3,ii)-cdgfac(3,ifac))*fextz(ii))            &
           +(1.d0-pond(ifac))*(pvar(jj)                         &
           -(xyzcen(1,jj)-cdgfac(1,ifac))*fextx(jj)             &
           -(xyzcen(2,jj)-cdgfac(2,ifac))*fexty(jj)             &
           -(xyzcen(3,jj)-cdgfac(3,ifac))*fextz(jj))
      pfac1  = pfac*surfac(1,ifac)
      pfac2  = pfac*surfac(2,ifac)
      pfac3  = pfac*surfac(3,ifac)
      bx(ii) = bx(ii) +pfac1
      by(ii) = by(ii) +pfac2
      bz(ii) = bz(ii) +pfac3
      bx(jj) = bx(jj) -pfac1
      by(jj) = by(jj) -pfac2
      bz(jj) = bz(jj) -pfac3
    enddo

!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

    do ifac = 1,nfabor
      ii = ifabor(ifac)
      pfac = inc*coefap(ifac) +coefbp(ifac)*(pvar(ii)           &
           -(xyzcen(1,ii)-cdgfbo(1,ifac))*fextx(ii)             &
           -(xyzcen(2,ii)-cdgfbo(2,ifac))*fexty(ii)             &
           -(xyzcen(3,ii)-cdgfbo(3,ifac))*fextz(ii) )
      bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
      by(ii) = by(ii) +pfac*surfbo(2,ifac)
      bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
    enddo

  endif

!     DPDX,DPDY,DPDZ = GRADIENT

  do iel = 1, ncel
    unsvol = 1.d0/volume(iel)
    dpdx(iel) = bx(iel)*unsvol
    dpdy(iel) = by(iel)*unsvol
    dpdz(iel) = bz(iel)*unsvol
  enddo

!     TRAITEMENT DU PARALLELISME

  if(irangp.ge.0) then
    call parcom (dpdx)
    !==========
    call parcom (dpdy)
    !==========
    call parcom (dpdz)
    !==========
  endif

!     TRAITEMENT DE LA PERIODICITE

  if(iperio.eq.1) then
    call percom                                                   &
    !==========
  ( idimte , itenso ,                                             &
    dpdx   , dpdx   , dpdx  ,                                     &
    dpdy   , dpdy   , dpdy  ,                                     &
    dpdz   , dpdz   , dpdz  )
  endif

  return

endif


!===============================================================================
! 2. CONSTRUCTION DES GRADIENTS PAR UNE METHODE DE MOINDRE
!         CARRE POUR LES MAILLAGES TORDUS
!===============================================================================


if( (inicoc.eq.1.or.iale.eq.1.or.imobil.eq.1) .and.iccocg.eq.1) then


! --->  2.1 CALCUL COMPLET DE COCG ET
!           SAUVEGARDE DE LA CONTRIBUTION AUX CELLULES DE BORD

!   INITIALISATION

  do ii = 1, 3
    do jj = 1, 3
      do iel = 1, ncelet
        cocg(iel,ii,jj) = 0.d0
      enddo
    enddo
  enddo


!   ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    unsdij = 1.d0/sqrt( (xyzcen(1,ii)-xyzcen(1,jj))**2          &
                       +(xyzcen(2,ii)-xyzcen(2,jj))**2          &
                       +(xyzcen(3,ii)-xyzcen(3,jj))**2 )
    dsij(1) = (xyzcen(1,jj)-xyzcen(1,ii))*unsdij
    dsij(2) = (xyzcen(2,jj)-xyzcen(2,ii))*unsdij
    dsij(3) = (xyzcen(3,jj)-xyzcen(3,ii))*unsdij

    cocg(ii,1,1) = cocg(ii,1,1) +dsij(1)*dsij(1)
    cocg(ii,2,2) = cocg(ii,2,2) +dsij(2)*dsij(2)
    cocg(ii,3,3) = cocg(ii,3,3) +dsij(3)*dsij(3)
    cocg(ii,1,2) = cocg(ii,1,2) +dsij(1)*dsij(2)
    cocg(ii,1,3) = cocg(ii,1,3) +dsij(1)*dsij(3)
    cocg(ii,2,3) = cocg(ii,2,3) +dsij(2)*dsij(3)

    cocg(jj,1,1) = cocg(jj,1,1) +dsij(1)*dsij(1)
    cocg(jj,2,2) = cocg(jj,2,2) +dsij(2)*dsij(2)
    cocg(jj,3,3) = cocg(jj,3,3) +dsij(3)*dsij(3)
    cocg(jj,1,2) = cocg(jj,1,2) +dsij(1)*dsij(2)
    cocg(jj,1,3) = cocg(jj,1,3) +dsij(1)*dsij(3)
    cocg(jj,2,3) = cocg(jj,2,3) +dsij(2)*dsij(3)
  enddo

! ET COMPLEMENT POUR LE VOISINAGE ETENDU
!     PAS DE VECTORISATION PARTICULIERE A IMPOSER A PRIORI

  if (imrgra.eq.2.or.imrgra.eq.3) then

    do ii = 1, ncel
      do ipos = ipcvse(ii), ipcvse(ii+1)-1
        jj = ielvse(ipos)

        unsdij = 1.d0/sqrt( (xyzcen(1,ii)-xyzcen(1,jj))**2        &
                           +(xyzcen(2,ii)-xyzcen(2,jj))**2        &
                           +(xyzcen(3,ii)-xyzcen(3,jj))**2 )
        dsij(1) = (xyzcen(1,jj)-xyzcen(1,ii))*unsdij
        dsij(2) = (xyzcen(2,jj)-xyzcen(2,ii))*unsdij
        dsij(3) = (xyzcen(3,jj)-xyzcen(3,ii))*unsdij

        cocg(ii,1,1) = cocg(ii,1,1) +dsij(1)*dsij(1)
        cocg(ii,2,2) = cocg(ii,2,2) +dsij(2)*dsij(2)
        cocg(ii,3,3) = cocg(ii,3,3) +dsij(3)*dsij(3)
        cocg(ii,1,2) = cocg(ii,1,2) +dsij(1)*dsij(2)
        cocg(ii,1,3) = cocg(ii,1,3) +dsij(1)*dsij(3)
        cocg(ii,2,3) = cocg(ii,2,3) +dsij(2)*dsij(3)

      enddo
    enddo

  endif

!  SAUVEGARDE DE LA CONTRIBUTION DES FACES INTERNES AUX ELEMENTS DE BORD

!CDIR NODEP
  do ii = 1, ncelbr
    iel = icelbr(ii)
    cocgb(ii,1,1) = cocg(iel,1,1)
    cocgb(ii,1,2) = cocg(iel,1,2)
    cocgb(ii,1,3) = cocg(iel,1,3)
    cocgb(ii,2,1) = cocg(iel,1,2)
    cocgb(ii,2,2) = cocg(iel,2,2)
    cocgb(ii,2,3) = cocg(iel,2,3)
    cocgb(ii,3,1) = cocg(iel,1,3)
    cocgb(ii,3,2) = cocg(iel,2,3)
    cocgb(ii,3,3) = cocg(iel,3,3)
  enddo

  inicoc = 0

!   ASSEMBLAGE A PARTIR DES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    extrab = 1.d0-isympa(ifac)*extrap*coefbp(ifac)
    umcbsd = extrab*(1.d0-coefbp(ifac))/distbr(ifac)
    unssbn = extrab/surfbn(ifac)

    dsij(1) = surfbo(1,ifac)*unssbn +umcbsd*diipb(1,ifac)
    dsij(2) = surfbo(2,ifac)*unssbn +umcbsd*diipb(2,ifac)
    dsij(3) = surfbo(3,ifac)*unssbn +umcbsd*diipb(3,ifac)

    cocg(ii,1,1) = cocg(ii,1,1) +dsij(1)*dsij(1)
    cocg(ii,2,2) = cocg(ii,2,2) +dsij(2)*dsij(2)
    cocg(ii,3,3) = cocg(ii,3,3) +dsij(3)*dsij(3)
    cocg(ii,1,2) = cocg(ii,1,2) +dsij(1)*dsij(2)
    cocg(ii,1,3) = cocg(ii,1,3) +dsij(1)*dsij(3)
    cocg(ii,2,3) = cocg(ii,2,3) +dsij(2)*dsij(3)

  enddo

!   SYMETRISATION

  do iel = 1, ncel
    cocg(iel,2,1) = cocg(iel,1,2)
    cocg(iel,3,1) = cocg(iel,1,3)
    cocg(iel,3,2) = cocg(iel,2,3)
  enddo

!   INVERSION

  nbloc = ncel/lbloc
  if (nbloc.gt.0) then
    do ibloc = 1, nbloc
      do ii = 1, lbloc
        iel = (ibloc-1)*lbloc+ii

        cocg11 = cocg(iel,1,1)
        cocg12 = cocg(iel,1,2)
        cocg13 = cocg(iel,1,3)
        cocg21 = cocg(iel,2,1)
        cocg22 = cocg(iel,2,2)
        cocg23 = cocg(iel,2,3)
        cocg31 = cocg(iel,3,1)
        cocg32 = cocg(iel,3,2)
        cocg33 = cocg(iel,3,3)

        a11=cocg22*cocg33-cocg32*cocg23
        a12=cocg32*cocg13-cocg12*cocg33
        a13=cocg12*cocg23-cocg22*cocg13
        a22=cocg11*cocg33-cocg31*cocg13
        a23=cocg21*cocg13-cocg11*cocg23
        a33=cocg11*cocg22-cocg21*cocg12

        unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

        aa(ii,1,1) = a11 *unsdet
        aa(ii,1,2) = a12 *unsdet
        aa(ii,1,3) = a13 *unsdet
        aa(ii,2,2) = a22 *unsdet
        aa(ii,2,3) = a23 *unsdet
        aa(ii,3,3) = a33 *unsdet

      enddo

      do ii = 1, lbloc
        iel = (ibloc-1)*lbloc+ii
        cocg(iel,1,1) = aa(ii,1,1)
        cocg(iel,1,2) = aa(ii,1,2)
        cocg(iel,1,3) = aa(ii,1,3)
        cocg(iel,2,2) = aa(ii,2,2)
        cocg(iel,2,3) = aa(ii,2,3)
        cocg(iel,3,3) = aa(ii,3,3)
      enddo

    enddo

  endif

  irel = mod(ncel,lbloc)
  if (irel.gt.0) then
    ibloc = nbloc + 1
    do ii = 1, irel
      iel = (ibloc-1)*lbloc+ii

      cocg11 = cocg(iel,1,1)
      cocg12 = cocg(iel,1,2)
      cocg13 = cocg(iel,1,3)
      cocg21 = cocg(iel,2,1)
      cocg22 = cocg(iel,2,2)
      cocg23 = cocg(iel,2,3)
      cocg31 = cocg(iel,3,1)
      cocg32 = cocg(iel,3,2)
      cocg33 = cocg(iel,3,3)

      a11=cocg22*cocg33-cocg32*cocg23
      a12=cocg32*cocg13-cocg12*cocg33
      a13=cocg12*cocg23-cocg22*cocg13
      a22=cocg11*cocg33-cocg31*cocg13
      a23=cocg21*cocg13-cocg11*cocg23
      a33=cocg11*cocg22-cocg21*cocg12

      unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

      aa(ii,1,1) = a11 *unsdet
      aa(ii,1,2) = a12 *unsdet
      aa(ii,1,3) = a13 *unsdet
      aa(ii,2,2) = a22 *unsdet
      aa(ii,2,3) = a23 *unsdet
      aa(ii,3,3) = a33 *unsdet

    enddo

    do ii = 1, irel
      iel = (ibloc-1)*lbloc+ii
      cocg(iel,1,1) = aa(ii,1,1)
      cocg(iel,1,2) = aa(ii,1,2)
      cocg(iel,1,3) = aa(ii,1,3)
      cocg(iel,2,2) = aa(ii,2,2)
      cocg(iel,2,3) = aa(ii,2,3)
      cocg(iel,3,3) = aa(ii,3,3)
    enddo
  endif


!         MATRICE SYMETRIQUE

  do iel = 1, ncel
    cocg(iel,2,1) = cocg(iel,1,2)
    cocg(iel,3,1) = cocg(iel,1,3)
    cocg(iel,3,2) = cocg(iel,2,3)
  enddo

elseif(iccocg.eq.1) then

! --->  2.2 CALCUL AUX BORDS DE COCG

  do ii = 1, 3
    do jj = 1, 3
!CDIR NODEP
      do ielb = 1, ncelbr
        iel = icelbr(ielb)
        cocg(iel,ii,jj) = cocgb(ielb,ii,jj)
      enddo
    enddo
  enddo

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    extrab = 1.d0-isympa(ifac)*extrap*coefbp(ifac)
    umcbsd = extrab*(1.d0-coefbp(ifac))/distbr(ifac)
    unssbn = extrab/surfbn(ifac)

    dsij(1) = surfbo(1,ifac)*unssbn+umcbsd*diipb(1,ifac)
    dsij(2) = surfbo(2,ifac)*unssbn+umcbsd*diipb(2,ifac)
    dsij(3) = surfbo(3,ifac)*unssbn+umcbsd*diipb(3,ifac)

    cocg(ii,1,1) = cocg(ii,1,1) +dsij(1)*dsij(1)
    cocg(ii,2,2) = cocg(ii,2,2) +dsij(2)*dsij(2)
    cocg(ii,3,3) = cocg(ii,3,3) +dsij(3)*dsij(3)
    cocg(ii,1,2) = cocg(ii,1,2) +dsij(1)*dsij(2)
    cocg(ii,1,3) = cocg(ii,1,3) +dsij(1)*dsij(3)
    cocg(ii,2,3) = cocg(ii,2,3) +dsij(2)*dsij(3)

  enddo

!     SYMETRISATION

!CDIR NODEP
  do ielb = 1, ncelbr
    iel = icelbr(ielb)
    cocg(iel,2,1) = cocg(iel,1,2)
    cocg(iel,3,1) = cocg(iel,1,3)
    cocg(iel,3,2) = cocg(iel,2,3)
  enddo


  nbloc = ncelbr/lbloc
  if (nbloc.gt.0) then
    do ibloc = 1, nbloc
      do ii = 1, lbloc
        ielb = (ibloc-1)*lbloc+ii
        iel = icelbr(ielb)

        cocg11 = cocg(iel,1,1)
        cocg12 = cocg(iel,1,2)
        cocg13 = cocg(iel,1,3)
        cocg21 = cocg(iel,2,1)
        cocg22 = cocg(iel,2,2)
        cocg23 = cocg(iel,2,3)
        cocg31 = cocg(iel,3,1)
        cocg32 = cocg(iel,3,2)
        cocg33 = cocg(iel,3,3)

        a11=cocg22*cocg33-cocg32*cocg23
        a12=cocg32*cocg13-cocg12*cocg33
        a13=cocg12*cocg23-cocg22*cocg13
        a22=cocg11*cocg33-cocg31*cocg13
        a23=cocg21*cocg13-cocg11*cocg23
        a33=cocg11*cocg22-cocg21*cocg12

        unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

        aa(ii,1,1) = a11*unsdet
        aa(ii,1,2) = a12*unsdet
        aa(ii,1,3) = a13*unsdet
        aa(ii,2,2) = a22*unsdet
        aa(ii,2,3) = a23*unsdet
        aa(ii,3,3) = a33*unsdet

      enddo

!CDIR NODEP
      do ii = 1, lbloc
        ielb = (ibloc-1)*lbloc + ii
        iel = icelbr(ielb)
        cocg(iel,1,1) = aa(ii,1,1)
        cocg(iel,1,2) = aa(ii,1,2)
        cocg(iel,1,3) = aa(ii,1,3)
        cocg(iel,2,2) = aa(ii,2,2)
        cocg(iel,2,3) = aa(ii,2,3)
        cocg(iel,3,3) = aa(ii,3,3)
      enddo

    enddo

  endif

  irel  = mod(ncelbr,lbloc)
  if (irel .gt.0) then
    ibloc = nbloc + 1
    do ii = 1, irel
      ielb = (ibloc-1)*lbloc+ii
      iel = icelbr(ielb)

      cocg11 = cocg(iel,1,1)
      cocg12 = cocg(iel,1,2)
      cocg13 = cocg(iel,1,3)
      cocg21 = cocg(iel,2,1)
      cocg22 = cocg(iel,2,2)
      cocg23 = cocg(iel,2,3)
      cocg31 = cocg(iel,3,1)
      cocg32 = cocg(iel,3,2)
      cocg33 = cocg(iel,3,3)

      a11=cocg22*cocg33-cocg32*cocg23
      a12=cocg32*cocg13-cocg12*cocg33
      a13=cocg12*cocg23-cocg22*cocg13
      a22=cocg11*cocg33-cocg31*cocg13
      a23=cocg21*cocg13-cocg11*cocg23
      a33=cocg11*cocg22-cocg21*cocg12

      unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

      aa(ii,1,1) = a11*unsdet
      aa(ii,1,2) = a12*unsdet
      aa(ii,1,3) = a13*unsdet
      aa(ii,2,2) = a22*unsdet
      aa(ii,2,3) = a23*unsdet
      aa(ii,3,3) = a33*unsdet

    enddo

!CDIR NODEP
    do ii = 1, irel
      ielb = (ibloc-1)*lbloc + ii
      iel = icelbr(ielb)
      cocg(iel,1,1) = aa(ii,1,1)
      cocg(iel,1,2) = aa(ii,1,2)
      cocg(iel,1,3) = aa(ii,1,3)
      cocg(iel,2,2) = aa(ii,2,2)
      cocg(iel,2,3) = aa(ii,2,3)
      cocg(iel,3,3) = aa(ii,3,3)
    enddo

  endif


!     SYMETRISATION DE LA MATRICE INVERSEE

!CDIR NODEP
  do ielb = 1, ncelbr
    iel = icelbr(ielb)
    cocg(iel,2,1) = cocg(iel,1,2)
    cocg(iel,3,1) = cocg(iel,1,3)
    cocg(iel,3,2) = cocg(iel,2,3)
  enddo

endif

!===============================================================================
! 3. CALCUL DU SECOND MEMBRE
!===============================================================================

do iel = 1, ncelet
  bx(iel) = 0.d0
  by(iel) = 0.d0
  bz(iel) = 0.d0
enddo

!  CAS STANDARD, SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  ===============================================================
if (iphydp.eq.0) then

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

  do ifac = 1,nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    usdij2 = 1.d0/( ( xyzcen(1,jj)-xyzcen(1,ii) )**2            &
                   +( xyzcen(2,jj)-xyzcen(2,ii) )**2            &
                   +( xyzcen(3,jj)-xyzcen(3,ii) )**2 )
    vecfac = ( pvar(jj)-pvar(ii) )*usdij2

    pfsx = ( xyzcen(1,jj)-xyzcen(1,ii) )*vecfac
    pfsy = ( xyzcen(2,jj)-xyzcen(2,ii) )*vecfac
    pfsz = ( xyzcen(3,jj)-xyzcen(3,ii) )*vecfac
    bx(ii) = bx(ii) +pfsx
    by(ii) = by(ii) +pfsy
    bz(ii) = bz(ii) +pfsz
    bx(jj) = bx(jj) +pfsx
    by(jj) = by(jj) +pfsy
    bz(jj) = bz(jj) +pfsz
  enddo

! ET COMPLEMENT POUR LE VOISINAGE ETENDU
!     LA REECRITURE COMPLEXE DE LA BOUCLE EST DESTINEE A DIVISER LE CPU
!       TOTAL SUR VPP PAR 5. POUR CELA, ON UTILISE DPDX DPDY DPDZ COMME
!       TABLEAUX DE TRAVAIL LOCAUX
!     SI ON TOMBAIT SUR UN CAS PATHOLOGIQUES DANS LEQUEL NCEL < NOMBRE
!       DE CELLULES DU VOISINAGE ETENDU D'UNE CELLULE (AVEC PERIODICITE
!       MULTIPLE ET PEU DE MAILLES ???) ON FAIT UN TEST ET ON CONSERVE
!       LA BOUCLE INITIALE. ON AJOUTE CEPENDANT UN IF QUI PERMET SUR VPP
!       DE REDUIRE LE TEMPS CALCUL D'UN FACTEUR 4 SI LA BOUCLE N'EST PAS
!       SCINDEE EN DEUX.

  if (imrgra.eq.2.or.imrgra.eq.3) then

!     ON REGARDE SI NCEL >= NBRE DE VOISINS ETENDUS D'UNE CELLULE
    if (icesve.lt.0) then
      invemx = 0
      do ii = 1, ncel
        invemx = max(invemx,ipcvse(ii+1)-ipcvse(ii))
      enddo
      if (ncel.ge.invemx) then
        icesve = 1
      else
        icesve = 0
      endif
    endif

    if (icesve.gt.0) then

      do ii = 1, ncel
        do ipos = ipcvse(ii), ipcvse(ii+1)-1
          jj = ielvse(ipos)
          usdij2 = 1.d0/( ( xyzcen(1,jj)-xyzcen(1,ii) )**2        &
                         +( xyzcen(2,jj)-xyzcen(2,ii) )**2        &
                         +( xyzcen(3,jj)-xyzcen(3,ii) )**2 )
          vecfac = ( pvar(jj)-pvar(ii) )*usdij2

          iiii = ipos-ipcvse(ii)+1
          dpdx(iiii) = ( xyzcen(1,jj)-xyzcen(1,ii) )*vecfac
          dpdy(iiii) = ( xyzcen(2,jj)-xyzcen(2,ii) )*vecfac
          dpdz(iiii) = ( xyzcen(3,jj)-xyzcen(3,ii) )*vecfac
        enddo
        do ipos = ipcvse(ii), ipcvse(ii+1)-1
          iiii = ipos-ipcvse(ii)+1
          bx(ii) = bx(ii) +dpdx(iiii)
          by(ii) = by(ii) +dpdy(iiii)
          bz(ii) = bz(ii) +dpdz(iiii)
        enddo
      enddo

    else

      do ii = 1, ncel
        if(ipcvse(ii+1).gt.ipcvse(ii)) then
          do ipos = ipcvse(ii), ipcvse(ii+1)-1
            jj = ielvse(ipos)
            usdij2 = 1.d0/( ( xyzcen(1,jj)-xyzcen(1,ii) )**2      &
                           +( xyzcen(2,jj)-xyzcen(2,ii) )**2      &
                           +( xyzcen(3,jj)-xyzcen(3,ii) )**2 )
            vecfac = ( pvar(jj)-pvar(ii) )*usdij2

            pfsx = ( xyzcen(1,jj)-xyzcen(1,ii) )*vecfac
            pfsy = ( xyzcen(2,jj)-xyzcen(2,ii) )*vecfac
            pfsz = ( xyzcen(3,jj)-xyzcen(3,ii) )*vecfac
            bx(ii) = bx(ii) +pfsx
            by(ii) = by(ii) +pfsy
            bz(ii) = bz(ii) +pfsz
          enddo
        endif
      enddo

    endif

  endif


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

  do ifac = 1,nfabor

    ii = ifabor(ifac)

    extrab = (1.d0-isympa(ifac)*extrap*coefbp(ifac))**2
    unsdij = 1.d0/distbr(ifac)
    unssbn = 1.d0/surfbn(ifac)
    umcbsd = (1.d0-coefbp(ifac))*unsdij

    dsij(1) = surfbo(1,ifac)*unssbn + umcbsd*diipb(1,ifac)
    dsij(2) = surfbo(2,ifac)*unssbn + umcbsd*diipb(2,ifac)
    dsij(3) = surfbo(3,ifac)*unssbn + umcbsd*diipb(3,ifac)
    rkij = ( inc*coefap(ifac)                                   &
            +(coefbp(ifac)-1.d0)*pvar(ii) )*unsdij*extrab

    bx(ii) = bx(ii) +dsij(1)*rkij
    by(ii) = by(ii) +dsij(2)*rkij
    bz(ii) = bz(ii) +dsij(3)*rkij

  enddo

!  CAS AVEC PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  =====================================================
elseif(iphydp.ne.0) then

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

  do ifac = 1,nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    usdij2 = 1.d0/( ( xyzcen(1,jj)-xyzcen(1,ii) )**2            &
         +( xyzcen(2,jj)-xyzcen(2,ii) )**2                      &
         +( xyzcen(3,jj)-xyzcen(3,ii) )**2 )
    vecfac = ( pvar(jj)-pvar(ii)                                &
         +(xyzcen(1,ii)-cdgfac(1,ifac))*fextx(ii)               &
         +(xyzcen(2,ii)-cdgfac(2,ifac))*fexty(ii)               &
         +(xyzcen(3,ii)-cdgfac(3,ifac))*fextz(ii)               &
         -(xyzcen(1,jj)-cdgfac(1,ifac))*fextx(jj)               &
         -(xyzcen(2,jj)-cdgfac(2,ifac))*fexty(jj)               &
         -(xyzcen(3,jj)-cdgfac(3,ifac))*fextz(jj) )*usdij2

    pfsx = ( xyzcen(1,jj)-xyzcen(1,ii) )*vecfac
    pfsy = ( xyzcen(2,jj)-xyzcen(2,ii) )*vecfac
    pfsz = ( xyzcen(3,jj)-xyzcen(3,ii) )*vecfac
    bx(ii) = bx(ii) +pfsx
    by(ii) = by(ii) +pfsy
    bz(ii) = bz(ii) +pfsz
    bx(jj) = bx(jj) +pfsx
    by(jj) = by(jj) +pfsy
    bz(jj) = bz(jj) +pfsz
  enddo


! ET COMPLEMENT POUR LE VOISINAGE ETENDU
!     ON SUPPOSE QUE LE MILIEU DU SEGMENT JOIGNANT LES CENTRES
!     PEUT REMPLACER LE CENTRE DE GRAVITE DE LA FACE FICTIVE

!     LA REECRITURE COMPLEXE DE LA BOUCLE EST DESTINEE A DIVISER LE CPU
!       TOTAL SUR VPP PAR 5. POUR CELA, ON UTILISE DPDX DPDY DPDZ COMME
!       TABLEAUX DE TRAVAIL LOCAUX
!     SI ON TOMBAIT SUR UN CAS PATHOLOGIQUES DANS LEQUEL NCEL < NOMBRE
!       DE CELLULES DU VOISINAGE ETENDU D'UNE CELLULE (AVEC PERIODICITE
!       MULTIPLE ET PEU DE MAILLES ???) ON FAIT UN TEST ET ON CONSERVE
!       LA BOUCLE INITIALE. ON AJOUTE CEPENDAT UN IF QUI PERMET SUR VPP
!       DE REDUIRE LE TEMPS CALCUL D'UN FACTEUR 4 SI LA BOUCLE N'EST PAS
!       SCINDEE EN DEUX.

  if (imrgra.eq.2.or.imrgra.eq.3) then

!     ON REGARDE SI NCEL >= NBRE DE VOISINS ETENDUS D'UNE CELLULE
    if (icesve.lt.0) then
      invemx = 0
      do ii = 1, ncel
        invemx = max(invemx,ipcvse(ii+1)-ipcvse(ii))
      enddo
      if (ncel.ge.invemx) then
        icesve = 1
      else
        icesve = 0
      endif
    endif

    if (icesve.gt.0) then

      do ii = 1, ncel
        do ipos = ipcvse(ii), ipcvse(ii+1)-1
          jj = ielvse(ipos)
          usdij2 = 1.d0/( ( xyzcen(1,jj)-xyzcen(1,ii) )**2        &
                         +( xyzcen(2,jj)-xyzcen(2,ii) )**2        &
                         +( xyzcen(3,jj)-xyzcen(3,ii) )**2 )
          ptmidx = 0.5d0*(xyzcen(1,jj)+xyzcen(1,ii))
          ptmidy = 0.5d0*(xyzcen(2,jj)+xyzcen(2,ii))
          ptmidz = 0.5d0*(xyzcen(3,jj)+xyzcen(3,ii))
          vecfac = ( pvar(jj)-pvar(ii)                            &
               +(xyzcen(1,ii)-ptmidx)*fextx(ii)                   &
               +(xyzcen(2,ii)-ptmidy)*fexty(ii)                   &
               +(xyzcen(3,ii)-ptmidz)*fextz(ii)                   &
               -(xyzcen(1,jj)-ptmidx)*fextx(jj)                   &
               -(xyzcen(2,jj)-ptmidy)*fexty(jj)                   &
               -(xyzcen(3,jj)-ptmidz)*fextz(jj) )*usdij2

          iiii = ipos-ipcvse(ii)+1
          dpdx(iiii) = ( xyzcen(1,jj)-xyzcen(1,ii) )*vecfac
          dpdy(iiii) = ( xyzcen(2,jj)-xyzcen(2,ii) )*vecfac
          dpdz(iiii) = ( xyzcen(3,jj)-xyzcen(3,ii) )*vecfac
        enddo
        do ipos = ipcvse(ii), ipcvse(ii+1)-1
          iiii = ipos-ipcvse(ii)+1
          bx(ii) = bx(ii) +dpdx(iiii)
          by(ii) = by(ii) +dpdy(iiii)
          bz(ii) = bz(ii) +dpdz(iiii)
        enddo
      enddo

    else

      do ii = 1, ncel
        if(ipcvse(ii+1).gt.ipcvse(ii)) then
          do ipos = ipcvse(ii), ipcvse(ii+1)-1
            jj = ielvse(ipos)
            usdij2 = 1.d0/( ( xyzcen(1,jj)-xyzcen(1,ii) )**2      &
                           +( xyzcen(2,jj)-xyzcen(2,ii) )**2      &
                           +( xyzcen(3,jj)-xyzcen(3,ii) )**2 )
            ptmidx = 0.5d0*(xyzcen(1,jj)+xyzcen(1,ii))
            ptmidy = 0.5d0*(xyzcen(2,jj)+xyzcen(2,ii))
            ptmidz = 0.5d0*(xyzcen(3,jj)+xyzcen(3,ii))
            vecfac = ( pvar(jj)-pvar(ii)                          &
                 +(xyzcen(1,ii)-ptmidx)*fextx(ii)                 &
                 +(xyzcen(2,ii)-ptmidy)*fexty(ii)                 &
                 +(xyzcen(3,ii)-ptmidz)*fextz(ii)                 &
                 -(xyzcen(1,jj)-ptmidx)*fextx(jj)                 &
                 -(xyzcen(2,jj)-ptmidy)*fexty(jj)                 &
                 -(xyzcen(3,jj)-ptmidz)*fextz(jj) )*usdij2

            pfsx = ( xyzcen(1,jj)-xyzcen(1,ii) )*vecfac
            pfsy = ( xyzcen(2,jj)-xyzcen(2,ii) )*vecfac
            pfsz = ( xyzcen(3,jj)-xyzcen(3,ii) )*vecfac
            bx(ii) = bx(ii) +pfsx
            by(ii) = by(ii) +pfsy
            bz(ii) = bz(ii) +pfsz
          enddo
        endif
      enddo

    endif

  endif


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD
!      NOTER QUE EXTRAB NE DOIT PRENDRE QUE LES VALEURS 0 OU 1

  do ifac = 1,nfabor

    ii = ifabor(ifac)

    extrab = (1.d0-isympa(ifac)*extrap*coefbp(ifac))**2
    unsdij = 1.d0/distbr(ifac)
    unssbn = 1.d0/surfbn(ifac)
    umcbsd = (1.d0-coefbp(ifac))*unsdij

    dsij(1) = surfbo(1,ifac)*unssbn + umcbsd*diipb(1,ifac)
    dsij(2) = surfbo(2,ifac)*unssbn + umcbsd*diipb(2,ifac)
    dsij(3) = surfbo(3,ifac)*unssbn + umcbsd*diipb(3,ifac)
    rkij = ( inc*coefap(ifac)                                   &
         +(coefbp(ifac)-1.d0)*pvar(ii) )*unsdij                 &
         +(coefbp(ifac)-1.d0)*unsdij*(                          &
         (cdgfbo(1,ifac)-xyzcen(1,ii))*fextx(ii)                &
         +(cdgfbo(2,ifac)-xyzcen(2,ii))*fexty(ii)               &
         +(cdgfbo(3,ifac)-xyzcen(3,ii))*fextz(ii) )
    rkij = rkij*extrab

    bx(ii) = bx(ii) +dsij(1)*rkij
    by(ii) = by(ii) +dsij(2)*rkij
    bz(ii) = bz(ii) +dsij(3)*rkij

  enddo

endif

!===============================================================================
! 4. RESOLUTION
!===============================================================================


do iel = 1, ncel
  dpdx(iel) = cocg(iel,1,1)*bx(iel)+cocg(iel,1,2)*by(iel)         &
             +cocg(iel,1,3)*bz(iel)
  dpdy(iel) = cocg(iel,2,1)*bx(iel)+cocg(iel,2,2)*by(iel)         &
             +cocg(iel,2,3)*bz(iel)
  dpdz(iel) = cocg(iel,3,1)*bx(iel)+cocg(iel,3,2)*by(iel)         &
             +cocg(iel,3,3)*bz(iel)
enddo

!     AJOUT DE LA COMPOSANTE HYDROSTATIQUE
if (iphydp.eq.1) then
  do iel = 1, ncel
    dpdx(iel) = dpdx(iel) + fextx(iel)
    dpdy(iel) = dpdy(iel) + fexty(iel)
    dpdz(iel) = dpdz(iel) + fextz(iel)
  enddo
endif

!     TRAITEMENT DU PARALLELISME

  if(irangp.ge.0) then
    call parcom (dpdx)
    !==========
    call parcom (dpdy)
    !==========
    call parcom (dpdz)
    !==========
  endif

!     TRAITEMENT DE LA PERIODICITE

  if(iperio.eq.1) then
    call percom                                                   &
    !==========
  ( idimte , itenso ,                                             &
    dpdx   , dpdx   , dpdx  ,                                     &
    dpdy   , dpdy   , dpdy  ,                                     &
    dpdz   , dpdz   , dpdz  )
  endif


!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
