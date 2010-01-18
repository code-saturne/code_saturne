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

subroutine gradrc &
!================

 ( ncelet , ncel   , nfac   , nfabor , ncelbr ,                   &
   imrgra , inc    , iccocg , nswrgp , idimte , itenso , iphydp , &
   iwarnp , nfecra , epsrgp , extrap ,                            &
   ifacel , ifabor , icelbr , ivar   ,                            &
   volume , surfac , surfbo , pond   , xyzcen , cdgfac , cdgfbo , &
   dijpf  , diipb  , dofij  , fextx  , fexty  , fextz  ,          &
   coefap , coefbp , pvar   ,                                     &
   cocgb  , cocg   ,                                              &
   dpdx   , dpdy   , dpdz   ,                                     &
   bx     , by     , bz     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DU GRADIENT D'UNE VARIABLE AVEC UNE TECHNIQUE IERATIVE
! DE RECONSTRUCTION POUR LES MAILLAGES NON ORTHOGONAUX (NSWRGP>1)
! AVEC PRISE EN COMPTE EVENTUELLE D'UN TERME DE FORCE VOLUMIQUE
! GENERANT UNE COMPOSANTE DE PRESSION HYDROSTATIQUE

! CALCUL DE COCG POUR PRENDRE EN COMPTE LES C.L VARIABLES (FLUX)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au  moins              !
!                  !    !     ! face de bord                                   !
! imrgra           ! e  ! <-- ! methode de calcul du gradient                  !
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
! ivar             ! e  ! <-- ! indicateur du numero de la variable            !
!                  !    !     ! (ou 0 si variable non resolue)                 !
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
! iphydp           ! e  ! <-- ! indicateur de prise en compte de la            !
!                  !    !     ! pression hydrostatique                         !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! ifacel(2,nfac    ! te ! <-- ! no des elts voisins d'une face intern          !
! ifabor(nfabor    ! te ! <-- ! no de l'elt voisin d'une face de bord          !
! icelbr           ! te ! <-- ! numero global des elements ayant au            !
! (ncelbr)         !    !     !  moins une face de bord                        !
! volume(ncelet    ! tr ! <-- ! volume des elements                            !
! surfac(3,nfac    ! tr ! <-- ! surf vectorielle des surfaces interne          !
! surfbo           ! tr ! <-- ! surf vectorielle des surfaces de bord          !
!   (3,nfabor)     !    !     !                                                !
! pond(nfac)       ! tr ! <-- ! ponderation geometrique (entre 0 et 1          !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (3,ncelet        !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (3,nfac)         !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (3,nfabor)       !    !     !                                                !
! dijpf(3,nfac)    ! tr ! <-- ! vect i'j', i' (resp. j') projection            !
!                  !    !     !  du centre i (resp. j) sur la normale          !
!                  !    !     !  a la face interne                             !
! diipb            ! tr ! <-- ! vect ii', ii projection du centre i            !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! pvar  (ncelet    ! tr ! <-- ! variable                                       !
! fextx,y,z        ! tr ! <-- ! force exterieure generant la pression          !
!   (ncelet)       !    !     !  hydrostatique                                 !
! cocgb            ! tr ! <-- ! contribution des faces internes a              !
!  (ncelbr,3,3)    !    !     !  cocg pour les cellules de bord                !
! cocg(ncelet,3    ! tr ! <-- ! couplage des composantes du gradient           !
!     ,3)          !    !     ! lors de la reconstruction                      !
! dpdx dpdy        ! tr ! <-- ! gradient de pvar                               !
! dpdz (ncelet     ! tr !     !   (halo rempli pour la periodicite)            !
! bx,y,z(ncelet    ! tr ! --- ! tableau de travail pour le grad de p           !
! dofij            ! tr ! --> ! vecteur of pour les faces internes             !
! (ndim,nfac  )    !    !     ! o : intersection de ij et la face              !
!                  !    !     ! f : centre de la face                          !
! ra(*)            ! tr ! --- ! tableau de travail pour les reels              !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "cstphy.h"
include "cstnum.h"
include "vector.h"
include "albase.h"
include "cplsat.h"
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer          ncelet , ncel   , nfac   , nfabor , ncelbr
integer          imrgra , inc    , iccocg , nswrgp
integer          ivar   , idimte , itenso , iphydp
integer          iwarnp , nfecra
double precision epsrgp , extrap

integer          ifacel(2,nfac),ifabor(nfabor)
integer          icelbr(ncelbr)
double precision volume(ncelet), surfac(3,nfac)
double precision surfbo(3,nfabor)
double precision pond(nfac)
double precision xyzcen(3,ncelet)
double precision cdgfac(3,nfac),cdgfbo(3,nfabor)
double precision dijpf(3,nfac), diipb(3,nfabor)
double precision coefap(nfabor), coefbp(nfabor), pvar(ncelet)
double precision cocgb(ncelbr,3,3), cocg(ncelet,3,3)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision bx   (ncelet),by   (ncelet),bz   (ncelet)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision dofij(3,nfac)

! VARIABLES LOCALES

integer          lbloc
parameter       (lbloc = 1024)
integer          iel, ifac, ii, jj, kk, ll, mm
integer          isqrt, isweep, nswmax
integer          ibloc,nbloc,irel
double precision pfac,pip,deltpx,deltpy,deltpz
double precision rnorx,rnory,rnorz,rnorm,residu
double precision pfaci
double precision dof,fmoyx,fmoyy,fmoyz
double precision aa(lbloc,3,3),unsdet,pfsx,pfsy,pfsz
double precision a11,a12,a13,a21,a22,a23,a31,a32,a33
double precision cocg11,cocg12,cocg13,cocg21,cocg22,cocg23
double precision cocg31,cocg32,cocg33
double precision pfac1, pfac2, pfac3, unsvol, vecfac

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

isqrt = 1
nswmax = nswrgp
isweep = 1


!===============================================================================
! 2. CALCUL SANS RECONSTRUCTION
!===============================================================================

!     SI INITIALISATION PAR MOINDRES CARRES (IMRGRA=4), B EST DEJA REMPLI
!     SINON (IMRGRA=0) ON CALCULE UN GRADIENT SANS RECONSTRUCTION

if (imrgra.eq.0) then

  do iel = 1, ncelet
    bx  (iel) = 0.d0
    by  (iel) = 0.d0
    bz  (iel) = 0.d0
  enddo

!  CAS STANDARD, SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  ===============================================================
  if (iphydp.eq.0) then

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

    if (ivecti.eq.1) then

!CDIR NODEP
      do ifac = 1,nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        pfac   = pond(ifac)*pvar(ii) +(1.d0-pond(ifac))*pvar(jj)
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

    else

! VECTORISATION NON FORCEE
      do ifac = 1,nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        pfac   = pond(ifac)*pvar(ii) +(1.d0-pond(ifac))*pvar(jj)
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

    endif

!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

    if (ivectb.eq.1) then

!CDIR NODEP
      do ifac = 1,nfabor
        ii = ifabor(ifac)
        pfac = inc*coefap(ifac) +coefbp(ifac)*pvar(ii)
        bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
        by(ii) = by(ii) +pfac*surfbo(2,ifac)
        bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
      enddo

    else

! VECTORISATION NON FORCEE
      do ifac = 1,nfabor
        ii = ifabor(ifac)
        pfac = inc*coefap(ifac) +coefbp(ifac)*pvar(ii)
        bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
        by(ii) = by(ii) +pfac*surfbo(2,ifac)
        bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
      enddo

    endif

!  CAS AVEC PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  =====================================================
  else

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

    if (ivecti.eq.1) then

!CDIR NODEP
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

    else

! VECTORISATION NON FORCEE
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

    endif

!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

    if (ivectb.eq.1) then

!CDIR NODEP
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

    else

!  VECTORISATION NON FORCEE
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
    ( idimte , itenso ,                                           &
      dpdx   , dpdx   , dpdx  ,                                   &
      dpdy   , dpdy   , dpdy  ,                                   &
      dpdz   , dpdz   , dpdz  )
  endif

endif

if( nswrgp.le.1 ) return


!     On incremente IPASS quand on calcule COCG pour la premiere fois
ipass = ipass + 1

!===============================================================================
! 3. RECONSTRUCTION DES GRADIENTS POUR LES MAILLAGES TORDUS
!===============================================================================

!     RESOLUTION SEMI-IMPLICITE SUR TOUT LE MAILLAGE
!     DPDX,DY,DZ = GRADIENT

if(ipass.eq.1 .or. iale.eq.1 .or. imobil.eq.1) then

! ---> CALCUL DE COCG

  do ii = 1, 3
    do jj = 1, 3
      do iel =1,ncelet
        cocg(iel,ii,jj) = 0.d0
      enddo
    enddo
  enddo
  do iel=1,ncel
    cocg(iel,1,1) = volume(iel)
    cocg(iel,2,2) = volume(iel)
    cocg(iel,3,3) = volume(iel)
  enddo

! ---> AJOUT DES CONTRIBUTIONS DES FACES INTERNES
  do ll =1,3
    do mm =1,3

      if (ivecti.eq.1) then

!CDIR NODEP
        do ifac=1,nfac
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

!---> DOF = OF
          dof = dofij(mm,ifac)

          pfaci  = -dof*0.5d0
          vecfac = pfaci*surfac(ll,ifac)
          cocg(ii,ll,mm) = cocg(ii,ll,mm) + vecfac
          cocg(jj,ll,mm) = cocg(jj,ll,mm) - vecfac
        enddo

      else

! VECTORISATION NON FORCEE
        do ifac=1,nfac
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

!---> DOF = OF
          dof = dofij(mm,ifac)

          pfaci  = -dof*0.5d0
          vecfac = pfaci*surfac(ll,ifac)
          cocg(ii,ll,mm) = cocg(ii,ll,mm) + vecfac
          cocg(jj,ll,mm) = cocg(jj,ll,mm) - vecfac
        enddo

      endif

    enddo
  enddo

! ---> SAUVEGADE DU COCG PARTIEL AUX FACES INTERNES DES CELLULES DE BORD
  do ii = 1, ncelbr
    iel = icelbr(ii)
    do ll = 1, 3
      do mm = 1, 3
        cocgb(ii,ll,mm) = cocg(iel,ll,mm)
      enddo
    enddo
  enddo

! ---> AJOUT DES CONTRIBUTIONS DES FACES DE BORD
  do ll =1,3
    do mm =1,3

      if (ivectb.eq.1) then

!CDIR NODEP
        do ifac=1,nfabor
          ii = ifabor(ifac)
          cocg(ii,ll,mm) = cocg(ii,ll,mm)                         &
            - coefbp(ifac)*diipb(mm,ifac)*surfbo(ll,ifac)
        enddo

      else

! VECTORISATION NON FORCEE
        do ifac=1,nfabor
          ii = ifabor(ifac)
          cocg(ii,ll,mm) = cocg(ii,ll,mm)                         &
            - coefbp(ifac)*diipb(mm,ifac)*surfbo(ll,ifac)
        enddo

      endif

    enddo
  enddo

! ---> ON INVERSE POUR TOUTE LES CELLULES : LE COCG POUR LES CELLULES INTERNES
!      RESTE ENSUITE LE MEME TANT QUE LE MAILLAGE NE CHANGE PAS
  nbloc = ncel/lbloc
  if (nbloc.gt.0) then
    do ibloc = 1, nbloc
      do ii =1, lbloc
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
        a21=cocg31*cocg23-cocg21*cocg33
        a22=cocg11*cocg33-cocg31*cocg13
        a23=cocg21*cocg13-cocg11*cocg23
        a31=cocg21*cocg32-cocg31*cocg22
        a32=cocg31*cocg12-cocg11*cocg32
        a33=cocg11*cocg22-cocg21*cocg12

        unsdet = 1.d0/(cocg11*a11 +cocg21*a12+cocg31*a13)

        aa(ii,1,1) = a11*unsdet
        aa(ii,1,2) = a12*unsdet
        aa(ii,1,3) = a13*unsdet
        aa(ii,2,1) = a21*unsdet
        aa(ii,2,2) = a22*unsdet
        aa(ii,2,3) = a23*unsdet
        aa(ii,3,1) = a31*unsdet
        aa(ii,3,2) = a32*unsdet
        aa(ii,3,3) = a33*unsdet

      enddo

      do ii = 1, lbloc
        iel = (ibloc-1)*lbloc+ii
        cocg(iel,1,1) = aa(ii,1,1)
        cocg(iel,1,2) = aa(ii,1,2)
        cocg(iel,1,3) = aa(ii,1,3)
        cocg(iel,2,1) = aa(ii,2,1)
        cocg(iel,2,2) = aa(ii,2,2)
        cocg(iel,2,3) = aa(ii,2,3)
        cocg(iel,3,1) = aa(ii,3,1)
        cocg(iel,3,2) = aa(ii,3,2)
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
      a21=cocg31*cocg23-cocg21*cocg33
      a22=cocg11*cocg33-cocg31*cocg13
      a23=cocg21*cocg13-cocg11*cocg23
      a31=cocg21*cocg32-cocg31*cocg22
      a32=cocg31*cocg12-cocg11*cocg32
      a33=cocg11*cocg22-cocg21*cocg12

      unsdet = 1.d0/(cocg11*a11 +cocg21*a12+cocg31*a13)

      aa(ii,1,1) = a11*unsdet
      aa(ii,1,2) = a12*unsdet
      aa(ii,1,3) = a13*unsdet
      aa(ii,2,1) = a21*unsdet
      aa(ii,2,2) = a22*unsdet
      aa(ii,2,3) = a23*unsdet
      aa(ii,3,1) = a31*unsdet
      aa(ii,3,2) = a32*unsdet
      aa(ii,3,3) = a33*unsdet

    enddo

    do ii = 1, irel
      iel = (ibloc-1)*lbloc+ii
      cocg(iel,1,1) = aa(ii,1,1)
      cocg(iel,1,2) = aa(ii,1,2)
      cocg(iel,1,3) = aa(ii,1,3)
      cocg(iel,2,1) = aa(ii,2,1)
      cocg(iel,2,2) = aa(ii,2,2)
      cocg(iel,2,3) = aa(ii,2,3)
      cocg(iel,3,1) = aa(ii,3,1)
      cocg(iel,3,2) = aa(ii,3,2)
      cocg(iel,3,3) = aa(ii,3,3)
    enddo

  endif

endif

! ---> SI ON DOIT RECALCULER COCG ENSUITE, ON NE LE FAIT PLUS
!      QUE POUR LES CELLULES DE BORD, AVEC LE COCGB STOCKE

if(iccocg.eq.1 .and. ipass.gt.1 .and. iale.eq.0 .and. imobil.eq.0) then

  do ll =1,3
    do mm =1,3

      do kk = 1, ncelbr
        iel = icelbr(kk)
        cocg(iel,ll,mm) = cocgb(kk,ll,mm)
      enddo

      if (ivectb.eq.1) then

!CDIR NODEP
        do ifac=1,nfabor
          ii = ifabor(ifac)
          cocg(ii,ll,mm) = cocg(ii,ll,mm)                         &
            - coefbp(ifac)*diipb(mm,ifac)*surfbo(ll,ifac)
        enddo

      else

! VECTORISATION NON FORCEE
        do ifac=1,nfabor
          ii = ifabor(ifac)
          cocg(ii,ll,mm) = cocg(ii,ll,mm)                         &
            - coefbp(ifac)*diipb(mm,ifac)*surfbo(ll,ifac)
        enddo

      endif

    enddo
  enddo

  do ii = 1, ncelbr

    iel = icelbr(ii)

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
    a21=cocg31*cocg23-cocg21*cocg33
    a22=cocg11*cocg33-cocg31*cocg13
    a23=cocg21*cocg13-cocg11*cocg23
    a31=cocg21*cocg32-cocg31*cocg22
    a32=cocg31*cocg12-cocg11*cocg32
    a33=cocg11*cocg22-cocg21*cocg12

    unsdet = 1.d0/(cocg11*a11 +cocg21*a12+cocg31*a13)

    a11 = a11*unsdet
    a12 = a12*unsdet
    a13 = a13*unsdet
    a21 = a21*unsdet
    a22 = a22*unsdet
    a23 = a23*unsdet
    a31 = a31*unsdet
    a32 = a32*unsdet
    a33 = a33*unsdet

    cocg(iel,1,1) = a11
    cocg(iel,1,2) = a12
    cocg(iel,1,3) = a13
    cocg(iel,2,1) = a21
    cocg(iel,2,2) = a22
    cocg(iel,2,3) = a23
    cocg(iel,3,1) = a31
    cocg(iel,3,2) = a32
    cocg(iel,3,3) = a33
  enddo

endif

! ---> CALCUL DU RESIDU DE NORMALISATION

call prods3(ncelet,ncel,isqrt,bx,bx,by,by,bz,bz,                  &
            rnorx,rnory,rnorz)
rnorm = rnorx +rnory +rnorz
if (volmax.gt.1.d0) rnorm = rnorm / volmax
if( rnorm.le.epzero ) return

!  LE VECTEUR OijFij EST CALCULE DANS CLDIJP

! ---> DEBUT DES ITERATIONS

 100  continue

isweep = isweep +1


!     CALCUL DU SECOND MEMBRE

do iel = 1, ncel
  bx(iel) = -dpdx(iel)*volume(iel)
  by(iel) = -dpdy(iel)*volume(iel)
  bz(iel) = -dpdz(iel)*volume(iel)
enddo

!  CAS STANDARD, SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  ===============================================================
if (iphydp.eq.0) then

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

  if (ivecti.eq.1) then

!CDIR NODEP
    do ifac = 1,nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      vecfac = pond(ifac)*pvar(ii)  +(1.d0-pond(ifac))*pvar(jj)   &
           +(dofij(1,ifac)*(dpdx(ii)+dpdx(jj))                    &
           + dofij(2,ifac)*(dpdy(ii)+dpdy(jj))                    &
           + dofij(3,ifac)*(dpdz(ii)+dpdz(jj)))*0.5d0
      pfsx = vecfac*surfac(1,ifac)
      pfsy = vecfac*surfac(2,ifac)
      pfsz = vecfac*surfac(3,ifac)
      bx(ii) = bx(ii) +pfsx
      by(ii) = by(ii) +pfsy
      bz(ii) = bz(ii) +pfsz
      bx(jj) = bx(jj) -pfsx
      by(jj) = by(jj) -pfsy
      bz(jj) = bz(jj) -pfsz
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1,nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      vecfac = pond(ifac)*pvar(ii)  +(1.d0-pond(ifac))*pvar(jj)   &
           +(dofij(1,ifac)*(dpdx(ii)+dpdx(jj))                    &
           + dofij(2,ifac)*(dpdy(ii)+dpdy(jj))                    &
           + dofij(3,ifac)*(dpdz(ii)+dpdz(jj)))*0.5d0
      pfsx = vecfac*surfac(1,ifac)
      pfsy = vecfac*surfac(2,ifac)
      pfsz = vecfac*surfac(3,ifac)
      bx(ii) = bx(ii) +pfsx
      by(ii) = by(ii) +pfsy
      bz(ii) = bz(ii) +pfsz
      bx(jj) = bx(jj) -pfsx
      by(jj) = by(jj) -pfsy
      bz(jj) = bz(jj) -pfsz
    enddo

  endif


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

  if (ivectb.eq.1) then

!CDIR NODEP
    do ifac = 1,nfabor
      ii = ifabor(ifac)
      pip =  pvar(ii)                                             &
           +diipb(1,ifac)*dpdx(ii) +diipb(2,ifac)*dpdy(ii)        &
           +diipb(3,ifac)*dpdz(ii)
      pfac = inc*coefap(ifac) +coefbp(ifac)*pip
      pfac1= pvar(ii) +(cdgfbo(1,ifac)-xyzcen(1,ii))*dpdx(ii)     &
           +(cdgfbo(2,ifac)-xyzcen(2,ii))*dpdy(ii)                &
           +(cdgfbo(3,ifac)-xyzcen(3,ii))*dpdz(ii)
      pfac = coefbp(ifac)*(extrap*pfac1 +(1.d0-extrap)*pfac)      &
           +(1.d0-coefbp(ifac))*pfac
      bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
      by(ii) = by(ii) +pfac*surfbo(2,ifac)
      bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1,nfabor
      ii = ifabor(ifac)
      pip =  pvar(ii)                                             &
           +diipb(1,ifac)*dpdx(ii) +diipb(2,ifac)*dpdy(ii)        &
           +diipb(3,ifac)*dpdz(ii)
      pfac = inc*coefap(ifac) +coefbp(ifac)*pip
      pfac1= pvar(ii) +(cdgfbo(1,ifac)-xyzcen(1,ii))*dpdx(ii)     &
           +(cdgfbo(2,ifac)-xyzcen(2,ii))*dpdy(ii)                &
           +(cdgfbo(3,ifac)-xyzcen(3,ii))*dpdz(ii)
      pfac = coefbp(ifac)*(extrap*pfac1 +(1.d0-extrap)*pfac)      &
           +(1.d0-coefbp(ifac))*pfac
      bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
      by(ii) = by(ii) +pfac*surfbo(2,ifac)
      bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
    enddo

  endif

!  CAS AVEC PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
!  =====================================================
else

!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

  if (ivecti.eq.1) then

!CDIR NODEP
    do ifac = 1,nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      fmoyx=0.5d0*(fextx(ii)+fextx(jj))
      fmoyy=0.5d0*(fexty(ii)+fexty(jj))
      fmoyz=0.5d0*(fextz(ii)+fextz(jj))
      vecfac   = pond(ifac)*(pvar(ii)                             &
           -(xyzcen(1,ii)-cdgfac(1,ifac))*(fextx(ii)-fmoyx)       &
           -(xyzcen(2,ii)-cdgfac(2,ifac))*(fexty(ii)-fmoyy)       &
           -(xyzcen(3,ii)-cdgfac(3,ifac))*(fextz(ii)-fmoyz))      &
           +(1.d0-pond(ifac))*(pvar(jj)                           &
           -(xyzcen(1,jj)-cdgfac(1,ifac))*(fextx(jj)-fmoyx)       &
           -(xyzcen(2,jj)-cdgfac(2,ifac))*(fexty(jj)-fmoyy)       &
           -(xyzcen(3,jj)-cdgfac(3,ifac))*(fextz(jj)-fmoyz))      &
           +(dofij(1,ifac)*(dpdx(ii)+dpdx(jj))                    &
           + dofij(2,ifac)*(dpdy(ii)+dpdy(jj))                    &
           + dofij(3,ifac)*(dpdz(ii)+dpdz(jj)))*0.5d0
      pfsx = vecfac*surfac(1,ifac)
      pfsy = vecfac*surfac(2,ifac)
      pfsz = vecfac*surfac(3,ifac)
      bx(ii) = bx(ii) +pfsx
      by(ii) = by(ii) +pfsy
      bz(ii) = bz(ii) +pfsz
      bx(jj) = bx(jj) -pfsx
      by(jj) = by(jj) -pfsy
      bz(jj) = bz(jj) -pfsz
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1,nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      fmoyx=0.5d0*(fextx(ii)+fextx(jj))
      fmoyy=0.5d0*(fexty(ii)+fexty(jj))
      fmoyz=0.5d0*(fextz(ii)+fextz(jj))
      vecfac   = pond(ifac)*(pvar(ii)                             &
           -(xyzcen(1,ii)-cdgfac(1,ifac))*(fextx(ii)-fmoyx)       &
           -(xyzcen(2,ii)-cdgfac(2,ifac))*(fexty(ii)-fmoyy)       &
           -(xyzcen(3,ii)-cdgfac(3,ifac))*(fextz(ii)-fmoyz))      &
           +(1.d0-pond(ifac))*(pvar(jj)                           &
           -(xyzcen(1,jj)-cdgfac(1,ifac))*(fextx(jj)-fmoyx)       &
           -(xyzcen(2,jj)-cdgfac(2,ifac))*(fexty(jj)-fmoyy)       &
           -(xyzcen(3,jj)-cdgfac(3,ifac))*(fextz(jj)-fmoyz))      &
           +(dofij(1,ifac)*(dpdx(ii)+dpdx(jj))                    &
           + dofij(2,ifac)*(dpdy(ii)+dpdy(jj))                    &
           + dofij(3,ifac)*(dpdz(ii)+dpdz(jj)))*0.5d0
      pfsx = vecfac*surfac(1,ifac)
      pfsy = vecfac*surfac(2,ifac)
      pfsz = vecfac*surfac(3,ifac)
      bx(ii) = bx(ii) +pfsx
      by(ii) = by(ii) +pfsy
      bz(ii) = bz(ii) +pfsz
      bx(jj) = bx(jj) -pfsx
      by(jj) = by(jj) -pfsy
      bz(jj) = bz(jj) -pfsz
    enddo

  endif


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

  if (ivectb.eq.1) then

!CDIR NODEP
    do ifac = 1,nfabor
      ii = ifabor(ifac)
      pip =  pvar(ii)                                             &
           +diipb(1,ifac)*dpdx(ii) +diipb(2,ifac)*dpdy(ii)        &
           +diipb(3,ifac)*dpdz(ii)
      pfac = inc*coefap(ifac) +coefbp(ifac)*(pip                  &
           -(xyzcen(1,ii)-cdgfbo(1,ifac)+diipb(1,ifac))*fextx(ii) &
           -(xyzcen(2,ii)-cdgfbo(2,ifac)+diipb(2,ifac))*fexty(ii) &
           -(xyzcen(3,ii)-cdgfbo(3,ifac)+diipb(3,ifac))*fextz(ii))
      pfac1= pvar(ii) +(cdgfbo(1,ifac)-xyzcen(1,ii))*dpdx(ii)     &
                    +(cdgfbo(2,ifac)-xyzcen(2,ii))*dpdy(ii)       &
                    +(cdgfbo(3,ifac)-xyzcen(3,ii))*dpdz(ii)
      pfac = coefbp(ifac)*(extrap*pfac1 +(1.d0-extrap)*pfac)      &
           +(1.d0-coefbp(ifac))*pfac
      bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
      by(ii) = by(ii) +pfac*surfbo(2,ifac)
      bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1,nfabor
      ii = ifabor(ifac)
      pip =  pvar(ii)                                             &
           +diipb(1,ifac)*dpdx(ii) +diipb(2,ifac)*dpdy(ii)        &
           +diipb(3,ifac)*dpdz(ii)
      pfac = inc*coefap(ifac) +coefbp(ifac)*(pip                  &
           -(xyzcen(1,ii)-cdgfbo(1,ifac)+diipb(1,ifac))*fextx(ii) &
           -(xyzcen(2,ii)-cdgfbo(2,ifac)+diipb(2,ifac))*fexty(ii) &
           -(xyzcen(3,ii)-cdgfbo(3,ifac)+diipb(3,ifac))*fextz(ii))
      pfac1= pvar(ii) +(cdgfbo(1,ifac)-xyzcen(1,ii))*dpdx(ii)     &
                    +(cdgfbo(2,ifac)-xyzcen(2,ii))*dpdy(ii)       &
                    +(cdgfbo(3,ifac)-xyzcen(3,ii))*dpdz(ii)
      pfac = coefbp(ifac)*(extrap*pfac1 +(1.d0-extrap)*pfac)      &
           +(1.d0-coefbp(ifac))*pfac
      bx(ii) = bx(ii) +pfac*surfbo(1,ifac)
      by(ii) = by(ii) +pfac*surfbo(2,ifac)
      bz(ii) = bz(ii) +pfac*surfbo(3,ifac)
    enddo

  endif

endif

!     INCREMENTATION DU GRADIENT

do iel =1,ncel

  deltpx =                                                        &
 cocg(iel,1,1)*bx(iel)+cocg(iel,1,2)*by(iel)+cocg(iel,1,3)*bz(iel)
  deltpy =                                                        &
 cocg(iel,2,1)*bx(iel)+cocg(iel,2,2)*by(iel)+cocg(iel,2,3)*bz(iel)
  deltpz =                                                        &
 cocg(iel,3,1)*bx(iel)+cocg(iel,3,2)*by(iel)+cocg(iel,3,3)*bz(iel)

  dpdx(iel) = dpdx(iel) +deltpx
  dpdy(iel) = dpdy(iel) +deltpy
  dpdz(iel) = dpdz(iel) +deltpz

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
  call percom                                                     &
  !==========
  ( idimte , itenso ,                                             &
    dpdx   , dpdx   , dpdx  ,                                     &
    dpdy   , dpdy   , dpdy  ,                                     &
    dpdz   , dpdz   , dpdz  )
endif


! ---> TEST DE CONVERGENCE

call prods3(ncelet,ncel,isqrt,bx,bx,by,by,bz,bz,                  &
            rnorx,rnory,rnorz)
residu = rnorx +rnory +rnorz
if (volmax.gt.1.d0) residu = residu / volmax

if( residu.le.epsrgp*rnorm) then
  if( iwarnp.ge.2 ) then
    write (nfecra,1000) isweep,residu/rnorm,rnorm,ivar
  endif
  goto 101
elseif( isweep.ge.nswmax ) then
  if( iwarnp.ge.0) then
     write (nfecra,1000)isweep,residu/rnorm,rnorm,ivar
     write (nfecra,1100)
  endif
  goto 101
else
  goto 100
endif

 101  continue



!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1000 format(1X,'GRADRC ISWEEP = ',I4,' RESIDU NORME: ',E11.4,          &
         ' NORME: ',E11.4,/,1X,'PARAMETRE IVAR = ',I4 )
 1100 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION :         NON CONVERGENCE DE GRADRC           ',/,&
'@    =========                                               ',/,&
'@                                                            '  )

#else

 1000 format(1X,'GRADRC ISWEEP = ',I4,' NORMED RESIDUAL: ',E11.4,       &
         ' NORM: ',E11.4,/,1X,'PARAMETER IVAR = ',I4 )
 1100 format(                                                           &
'@'                                                            ,/,&
'@ @@ WARNING:            NON CONVERGENCE OF GRADRC'           ,/,&
'@    ========'                                                ,/,&
'@'                                                              )

#endif

!----
! FIN
!----

return

end subroutine
