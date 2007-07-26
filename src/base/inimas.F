c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2007 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
                        SUBROUTINE INIMAS
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IVAR1  , IVAR2  , IVAR3  , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGU , IMLIGU ,
     &   IWARNU , NFECRA ,
     &   EPSRGU , CLIMGU , EXTRAU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISYMPA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   ROM    , ROMB   ,
     &   UX     , UY     , UZ     ,
     &   COEFAX , COEFAY , COEFAZ , COEFBX , COEFBY , COEFBZ ,
     &   FLUMAS , FLUMAB ,
     &   DPDX   , DPDY   , DPDZ   , DPDXA  , DPDYA  , DPDZA  ,
     &   QDMX   , QDMY   , QDMZ   , COEFQA ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC                                                               -->
CFONC INCREMENTATION DU FLUX DE MASSE A PARTIR DU CHAMP VECTORIEL ROM.U
CFONC   .     .          -->   -->
CFONC   m   = m   +(rom* U ) . n
CFONC    ij    ij         ij    ij
CFONC
CFONC
CFONC  Pour la reconstruction, grad(rho u) est calcule avec des
CFONC  conditions aux limites approchees :
CFONC     COEFA(rho u) = ROMB * COEFA(u)
CFONC     COEFB(rho u) = COEFB (u)
CFONC
CFONC  et pour le flux de masse au bord on ecrit
CFONC   FLUMAB = [ROMB*COEFA(u) + ROMB*COEFB(u)*Ui
CFONC                 + COEFB(u)*II'.grad(rho u) ].Sij
CFONC  ce qui utilise de petites approximations sur les
CFONC  non-orthogonalites (cf. notice)
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! IVAR1        ! E  !  -> ! VARIABLE DE LA DIRECTION 1           !
CARGU ! IVAR2        ! E  !  -> ! VARIABLE DE LA DIRECTION 2           !
CARGU ! IVAR3        ! E  !  -> ! VARIABLE DE LA DIRECTION 3           !
CARGU ! IMASPE       ! E  !  -> ! SUIVANT L'APPEL DE INIMAS            !
CARGU !              !    !     ! = 1 SI APPEL DE NAVSTO RESOLP        !
CARGU !              !    !     ! = 2 SI APPEL DE DIVRIJ               !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE COURANTE             !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFLMB0       ! E  !  -> ! =1 : FLUX DE MASSE ANNULE SYM-PAROI  !
CARGU ! INIT         ! E  !  -> ! > 0 : INITIALISATION DU FLUX DE MASSE!
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! IMRGRA       ! E  !  -> ! INDICATEUR = 0 GRADRC 97             !
CARGU !              ! E  !  -> !            = 1 GRADMC 99             !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! NSWRGU       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGU       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IWARNU       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! EPSRGU       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! CLIMGU       ! R  !  -> ! COEF GRADIENT*DISTANCE/ECART         !
CARGU ! EXTRAU       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ISYMPA       ! TE ! <-> ! ZERO POUR ANNULER LE FLUX DE MASSE   !
CARGU ! (NFABOR     )!    !     !(SYMETRIES ET PAROIS AVEC CL COUPLEES)!
CARGU !              !    !     ! UN SINON                             !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! ROM(NCELET   ! TR !  -> ! MASSE VOLUMIQUE AUX CELLULES         !
CARGU ! ROMB(NFABOR) ! TR !  -> ! MASSE VOLUMIQUE AUX BORDS            !
CARGU ! UX,Y,Z(NCELET! TR !  -> ! VITESSE                              !
CARGU ! COEFAX, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR UX, UY, UZ!
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! FLUMAS(NFAC) ! TR ! <-> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR ! <-> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! DPD. (NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU ! DPD.A(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU !              !    !     !  AVEC DECENTREMENT AMONT             !
CARGU ! QDM.(NCELET) ! TR !  -  ! TABLEAU DE TRAVAIL POUR LA QDM       !
CARGU ! COEFQA(NFAB,3! TR !  -  ! TABLEAU DE TRAVAIL CL DE LA QDM      !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!
c@argue
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IVAR1  , IVAR2  , IVAR3  , IMASPE , IPHAS
      INTEGER          IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG
      INTEGER          NSWRGU , IMLIGU
      INTEGER          IWARNU , NFECRA
      DOUBLE PRECISION EPSRGU , CLIMGU , EXTRAU
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR), ISYMPA(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION ROM(NCELET), ROMB(NFABOR)
      DOUBLE PRECISION UX(NCELET), UY(NCELET), UZ(NCELET)
      DOUBLE PRECISION COEFAX(NFABOR), COEFAY(NFABOR), COEFAZ(NFABOR)
      DOUBLE PRECISION COEFBX(NFABOR), COEFBY(NFABOR), COEFBZ(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION DPDXA(NCELET),DPDYA(NCELET),DPDZA(NCELET)
      DOUBLE PRECISION QDMX(NCELET),QDMY(NCELET),QDMZ(NCELET)
      DOUBLE PRECISION COEFQA(NDIMFB,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, II, JJ, IEL, IOF, III
      INTEGER          IDIMTE, ITENSO, IAPPEL, IPHYDP
      DOUBLE PRECISION PFAC,PIP,UXFAC,UYFAC,UZFAC
      DOUBLE PRECISION DOFX,DOFY,DOFZ,POND
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C ---> CALCUL DE LA QTE DE MOUVEMENT
C
C
      IF( INIT.EQ.1 ) THEN
        DO IFAC = 1, NFAC
          FLUMAS(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          FLUMAB(IFAC) = 0.D0
        ENDDO
C
      ELSEIF(INIT.NE.0) THEN
        WRITE(NFECRA,1000) INIT
        CALL CSEXIT (1)
      ENDIF
C
      DO IEL = 1, NCEL
        QDMX(IEL) = ROM(IEL)*UX(IEL)
        QDMY(IEL) = ROM(IEL)*UY(IEL)
        QDMZ(IEL) = ROM(IEL)*UZ(IEL)
      ENDDO
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (QDMX)
C       ===========
        CALL PARCOM (QDMY)
C       ===========
        CALL PARCOM (QDMZ)
C       ===========
      ENDIF
C
C ---> PERIODICITE SUR QDM
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 1
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    QDMX   , QDMX   , QDMX  ,
     &    QDMY   , QDMY   , QDMY  ,
     &    QDMZ   , QDMZ   , QDMZ  )
      ENDIF
C
      DO IFAC =1, NFABOR
        COEFQA(IFAC,1) = ROMB(IFAC)*COEFAX(IFAC)
        COEFQA(IFAC,2) = ROMB(IFAC)*COEFAY(IFAC)
        COEFQA(IFAC,3) = ROMB(IFAC)*COEFAZ(IFAC)
      ENDDO
C
C=======================================================================
C 2.  CALCUL DU FLUX DE MASSE SANS TECHNIQUE DE RECONSTRUCTION
C=======================================================================
C
      IF( NSWRGU.LE.1 ) THEN
C
C     FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          POND = RA(IPOND-1+IFAC)
C
          FLUMAS(IFAC) =  FLUMAS(IFAC)
     &     +(POND*QDMX(II)+(1.D0-POND)*QDMX(JJ) )*SURFAC(1,IFAC)
     &     +(POND*QDMY(II)+(1.D0-POND)*QDMY(JJ) )*SURFAC(2,IFAC)
     &     +(POND*QDMZ(II)+(1.D0-POND)*QDMZ(JJ) )*SURFAC(3,IFAC)
C
        ENDDO
C
C
C     FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          UXFAC = INC*COEFQA(IFAC,1) +COEFBX(IFAC)*ROMB(IFAC)*UX(II)
          UYFAC = INC*COEFQA(IFAC,2) +COEFBY(IFAC)*ROMB(IFAC)*UY(II)
          UZFAC = INC*COEFQA(IFAC,3) +COEFBZ(IFAC)*ROMB(IFAC)*UZ(II)
C
          FLUMAB(IFAC) = FLUMAB(IFAC)
     &     +( UXFAC*SURFBO(1,IFAC)
     &                    +UYFAC*SURFBO(2,IFAC) +UZFAC*SURFBO(3,IFAC) )
C
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 4.  CALCUL DU FLUX DE MASSE AVEC TECHNIQUE DE RECONSTRUCTION
C        SI LE MAILLAGE EST NON ORTHOGONAL
C=======================================================================
C
      IF( NSWRGU.GT.1 ) THEN
C
C
C
C
C     TRAITEMENT DE LA PERIODICITE SPEFICIQUE A INIMAS AU DEBUT
C     =========================================================
C
        IF(IPEROT.GT.0) THEN
          IAPPEL = 1
C
          CALL PERMAS
C         ===========
     & ( IMASPE , IPHAS  , IAPPEL ,
     &   ROM    ,
     &   RA(IDUDXY) , RA(IDRDXY) , RA(IWDUDX) , RA(IWDRDX) )
C
        ENDIF
C
C     FLUX DE MASSE SUIVANT X
C     =======================
C
C ---> CALCUL DU GRADIENT
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR1  , IMRGRA , INC    , ICCOCG , NSWRGU , IMLIGU , IPHYDP ,
     &   IWARNU , NFECRA , EPSRGU , CLIMGU , EXTRAU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DPDXA  , DPDXA  , DPDXA  ,
     &   QDMX   , COEFQA(1,1) , COEFBX ,
     &   DPDX   , DPDY   , DPDZ   ,
C        ------   ------   ------
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C
C
C ---> FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          POND = RA(IPOND-1+IFAC)
C
          IOF = IDOFIJ-1+3*(IFAC-1)
C---> DOF = OF
          DOFX = RA(IOF+1)
          DOFY = RA(IOF+2)
          DOFZ = RA(IOF+3)
C
          FLUMAS(IFAC) = FLUMAS(IFAC)
     &         +( POND*QDMX(II) +(1.D0-POND)*QDMX(JJ)
     &           +0.5D0*( DPDX(II) +DPDX(JJ) )*DOFX
     &           +0.5D0*( DPDY(II) +DPDY(JJ) )*DOFY
     &           +0.5D0*( DPDZ(II) +DPDZ(JJ) )*DOFZ    )*SURFAC(1,IFAC)
C
        ENDDO
C
C ---> FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
C
          PIP = ROMB(IFAC) * UX(II)
     &          +DPDX(II)*DIIPBX
     &          +DPDY(II)*DIIPBY +DPDZ(II)*DIIPBZ
          PFAC = INC*COEFQA(IFAC,1) +COEFBX(IFAC)*PIP
C
          FLUMAB(IFAC) = FLUMAB(IFAC) +PFAC*SURFBO(1,IFAC)
C
        ENDDO
C
C
C     FLUX DE MASSE SUIVANT Y
C     =======================
C
C ---> CALCUL DU GRADIENT
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR2  , IMRGRA , INC    , ICCOCG , NSWRGU , IMLIGU , IPHYDP ,
     &   IWARNU , NFECRA , EPSRGU , CLIMGU , EXTRAU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DPDXA  , DPDXA  , DPDXA  ,
     &   QDMY   , COEFQA(1,2) , COEFBY ,
     &   DPDX   , DPDY   , DPDZ   ,
C        ------   ------   ------
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C
C
C ---> FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          POND = RA(IPOND-1+IFAC)
C
          IOF = IDOFIJ-1+3*(IFAC-1)
C
C---> DOF = OF
C
          DOFX = RA(IOF+1)
          DOFY = RA(IOF+2)
          DOFZ = RA(IOF+3)
C
          FLUMAS(IFAC) = FLUMAS(IFAC)
     &         +( POND*QDMY(II) +(1.D0-POND)*QDMY(JJ)
     &           +0.5D0*( DPDX(II) +DPDX(JJ) )*DOFX
     &           +0.5D0*( DPDY(II) +DPDY(JJ) )*DOFY
     &           +0.5D0*( DPDZ(II) +DPDZ(JJ) )*DOFZ    )*SURFAC(2,IFAC)
C
        ENDDO
C
C ---> FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
C
          PIP = ROMB(IFAC) * UY(II)
     &        +DPDX(II)*DIIPBX
     &        +DPDY(II)*DIIPBY +DPDZ(II)*DIIPBZ
          PFAC = INC*COEFQA(IFAC,2) +COEFBY(IFAC)*PIP
C
          FLUMAB(IFAC) = FLUMAB(IFAC) +PFAC*SURFBO(2,IFAC)
C
        ENDDO
C
C     FLUX DE MASSE SUIVANT Z
C     =======================
C
C ---> CALCUL DU GRADIENT
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR3  , IMRGRA , INC    , ICCOCG , NSWRGU , IMLIGU , IPHYDP ,
     &   IWARNU , NFECRA , EPSRGU , CLIMGU , EXTRAU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DPDXA  , DPDXA  , DPDXA  ,
     &   QDMZ     , COEFQA(1,3) , COEFBZ ,
     &   DPDX   , DPDY   , DPDZ   ,
C        ------   ------   ------
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C
C     FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          POND = RA(IPOND-1+IFAC)
C
          IOF = IDOFIJ-1+3*(IFAC-1)
C
C---> DOF = OF
C
          DOFX = RA(IOF+1)
          DOFY = RA(IOF+2)
          DOFZ = RA(IOF+3)
C
          FLUMAS(IFAC) = FLUMAS(IFAC)
     &         +( POND*QDMZ(II) +(1.D0-POND)*QDMZ(JJ)
     &           +0.5D0*( DPDX(II) +DPDX(JJ) )*DOFX
     &           +0.5D0*( DPDY(II) +DPDY(JJ) )*DOFY
     &           +0.5D0*( DPDZ(II) +DPDZ(JJ) )*DOFZ    )*SURFAC(3,IFAC)
C
        ENDDO
C
C ---> FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
C
          PIP = ROMB(IFAC) * UZ(II)
     &        +DPDX(II)*DIIPBX
     &        +DPDY(II)*DIIPBY +DPDZ(II)*DIIPBZ
          PFAC = INC*COEFQA(IFAC,3) +COEFBZ(IFAC)*PIP
C
          FLUMAB(IFAC) = FLUMAB(IFAC) +PFAC*SURFBO(3,IFAC)
C
        ENDDO
C
C
C
C     TRAITEMENT DE LA PERIODICITE SPEFICIQUE A INIMAS A LA FIN
C     =========================================================
C
        IF(IPEROT.GT.0) THEN
          IAPPEL = 2
C
          CALL PERMAS
C         ===========
     & ( IMASPE , IPHAS  , IAPPEL ,
     &   ROM    ,
     &   RA(IDUDXY) , RA(IDRDXY) , RA(IWDUDX) , RA(IWDRDX) )
C
        ENDIF
C
C
C
C
      ENDIF
C
C
C=======================================================================
C 6.  POUR S'ASSURER DE LA NULLITE DU FLUX DE MASSE AUX LIMITES
C       SYMETRIES PAROIS COUPLEES
C=======================================================================
C
      IF(IFLMB0.EQ.1) THEN
C FORCAGE DE FLUMAB a 0 pour la vitesse'
        DO IFAC = 1, NFABOR
          IF(ISYMPA(IFAC).EQ.0) THEN
            FLUMAB(IFAC) = 0.D0
          ENDIF
        ENDDO
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT('INIMAS APPELE AVEC INIT =',I10)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
