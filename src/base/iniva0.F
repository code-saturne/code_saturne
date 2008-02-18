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
C     Copyright (C) 1998-2008 EDF S.A., France
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
                        SUBROUTINE INIVA0
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCOFAB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , FRCXT  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC INITIALISATION DES VARIABLES DE CALCUL, DU PAS DE TEMPS
CFONC   ET DU TABLEAU INDICATEUR DU CALCUL DE LA DISTANCE A LA PAROI
CFONC   AUX VALEURS PAR DEFAUT
CFONC                 AVANT LECTURE EVENTUELLE DU FICHIER SUITE ET
CFONC                 AVANT DE PASSER LA MAIN A L'UTILISATEUR
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
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
CARGU ! NCOFAB       ! E  !  -> ! NOMBRE DE COUPLES COEFA/B POUR LES CL!
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
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
CARGU !              !    !     !  A LA FACE INTERNE                   !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR ! <-> ! VALEUR DU PAS DE TEMPS               !
CARGU ! RTP          ! TR ! <-> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA COEFB  ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! FRCXT(NCELET,! TR ! <-> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
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
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NCOFAB
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,NCOFAB), COEFB(NFABOR,NCOFAB)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IIS   , IVAR  , IPHAS , IPHASS, ISCAL , IMOM
      INTEGER          IEL   , IFAC
      INTEGER          ICLIP , II    , JJ    , IDIM
      INTEGER          IIFLUM, IIFLUA
      INTEGER          IIROM , IIROMB, IIROMA
      INTEGER          IIVISL, IIVIST, IIVISA, IIVISM
      INTEGER          IICP  , IICPA
      INTEGER          IIVISS, IIPTOT
      INTEGER          IPTSNA, IPTSTA, IPTSCA
      INTEGER          IKIPH , IEIPH , IPHIPH, IFBIPH, IOMGIP
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          NN
      DOUBLE PRECISION RO0IPH, VISIPH
      DOUBLE PRECISION XXK, XCMU, TRII
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
C En compressible, ISYMPA initialise (= 1) car utile dans le calcul
C     du pas de temps variable avant passage dans les C.L.
C
      IF ( IPPMOD(ICOMPF).GE.0 ) THEN
        DO IFAC = 1, NFABOR*NPHAS
          IA(IISYMP+IFAC-1) = 1
        ENDDO
      ENDIF
C
C=======================================================================
C 2. PAS DE TEMPS
C=======================================================================
C
      DO IEL = 1, NCEL
        DT (IEL) = DTREF
      ENDDO
C
C=======================================================================
C 3.  INITIALISATION DES PROPRIETES PHYSIQUES
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
C     Masse volumique
        IIROM  = IPPROC(IROM  (IPHAS))
        IIROMB = IPPROB(IROM  (IPHAS))
        RO0IPH = RO0   (IPHAS)
C
C     Masse volumique aux cellules (et au pdt precedent si ordre2 ou icalhy)
        DO IEL = 1, NCEL
          PROPCE(IEL,IIROM)  = RO0IPH
        ENDDO
        IF(IROEXT(IPHAS).GT.0.OR.ICALHY.EQ.1) THEN
          IIROMA = IPPROC(IROMA (IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IIROMA) = PROPCE(IEL,IIROM)
          ENDDO
        ENDIF
C     Masse volumique aux faces de bord (et au pdt precedent si ordre2)
        DO IFAC = 1, NFABOR
          PROPFB(IFAC,IIROMB) = RO0IPH
        ENDDO
        IF(IROEXT(IPHAS).GT.0) THEN
          IIROMA = IPPROB(IROMA (IPHAS))
          DO IFAC = 1, NFABOR
            PROPFB(IFAC,IIROMA) = PROPFB(IFAC,IIROMB)
          ENDDO
        ENDIF
C
C     Viscosite moleculaire
        IIVISL = IPPROC(IVISCL(IPHAS))
        IIVIST = IPPROC(IVISCT(IPHAS))
        VISIPH = VISCL0(IPHAS)
C
C     Viscosite moleculaire aux cellules (et au pdt precedent si ordre2)
        DO IEL = 1, NCEL
          PROPCE(IEL,IIVISL) = VISIPH
        ENDDO
        IF(IVIEXT(IPHAS).GT.0) THEN
          IIVISA = IPPROC(IVISLA(IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IIVISA) = PROPCE(IEL,IIVISL)
          ENDDO
        ENDIF
C     Viscosite turbulente aux cellules (et au pdt precedent si ordre2)
        DO IEL = 1, NCEL
          PROPCE(IEL,IIVIST) = 0.D0
        ENDDO
        IF(IVIEXT(IPHAS).GT.0) THEN
          IIVISA = IPPROC(IVISTA(IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IIVISA) = PROPCE(IEL,IIVIST)
          ENDDO
        ENDIF
C
C     Chaleur massique aux cellules (et au pdt precedent si ordre2)
        IF(ICP(IPHAS).GT.0) THEN
          IICP = IPPROC(ICP(IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IICP) = CP0(IPHAS)
          ENDDO
          IF(ICPEXT(IPHAS).GT.0) THEN
            IICPA  = IPPROC(ICPA(IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IICPA ) = PROPCE(IEL,IICP)
            ENDDO
          ENDIF
        ENDIF
C
C La pression totale sera initialisee a P0 + rho.g.r dans INIVAR
C  si l'utilisateur n'a pas fait d'initialisation personnelle
C Non valable en compressible
        IF (IPPMOD(ICOMPF).LT.0) THEN
          IIPTOT = IPPROC(IPRTOT(IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IIPTOT) = - RINFIN
          ENDDO
        ENDIF
C
      ENDDO
C
C
C     Diffusivite des scalaires
      DO ISCAL = 1, NSCAL
        IF(IVISLS(ISCAL).GT.0) THEN
          IIVISS = IPPROC(IVISLS(ISCAL))
C     Diffusivite aux cellules (et au pdt precedent si ordre2)
          DO IEL = 1, NCEL
            PROPCE(IEL,IIVISS) = VISLS0(ISCAL)
          ENDDO
          IF(IVSEXT(IPHSCA(ISCAL)).GT.0) THEN
            IIVISA = IPPROC(IVISSA(ISCAL))
            DO IEL = 1, NCEL
              PROPCE(IEL,IIVISA) = PROPCE(IEL,IIVISS)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C
C
C     Viscosite de maillage en ALE
      IF (IALE.EQ.1) THEN
        NN = 1
        IF (IORTVM.EQ.1) NN = 3
        DO II = 1, NN
          IIVISM = IPPROC(IVISMA(II))
          DO IEL = 1, NCEL
            PROPCE(IEL,IIVISM) = 1.D0
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
C 4. INITIALISATION STANDARD DES VARIABLES DE CALCUL
C     On complete ensuite pour les variables turbulentes et les scalaires
C=======================================================================
C
C     Toutes les variables a 0
      DO IVAR = 1, NVAR
        DO IEL = 1, NCEL
          RTP(IEL,IVAR) = 0.D0
        ENDDO
      ENDDO
C
C     On met la pression P* a PRED0
      DO IPHAS = 1, NPHAS
        DO IEL = 1, NCEL
          RTP(IEL,IPR(IPHAS)) = PRED0(IPHAS)
        ENDDO
      ENDDO
C
C     Couplage U-P
      IF(IPUCOU.EQ.1) THEN
        DO IEL = 1, NCEL
          RA(ITPUCO+IEL-1         ) = 0.D0
          RA(ITPUCO+IEL-1+  NCELET) = 0.D0
          RA(ITPUCO+IEL-1+2*NCELET) = 0.D0
        ENDDO
      ENDIF
C
C=======================================================================
C 5. INITIALISATION DE K, RIJ ET EPS
C=======================================================================
C
C  Si UREF n'a pas ete donnee par l'utilisateur ou a ete mal initialisee
C    (valeur negative), on met les valeurs de k, Rij, eps et omega a
C    -10*GRAND. On testera ensuite si l'utilisateur les a modifiees dans
C    usiniv ou en lisant un fichier suite.
C
      DO IPHAS = 1, NPHAS
C
        IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50) THEN
C
          IKIPH  = IK (IPHAS)
          IEIPH  = IEP(IPHAS)
C
          XCMU = CMU
          IF (ITURB(IPHAS).EQ.50) XCMU = CV2FMU
C
C
          IF (UREF(IPHAS).GE.0.D0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IKIPH) = 1.5D0*(0.02D0*UREF(IPHAS))**2
              RTP(IEL,IEIPH) = RTP(IEL,IKIPH)**1.5D0*XCMU/ALMAX(IPHAS)
            ENDDO
C
            ICLIP = 1
            IPHASS = IPHAS
            CALL CLIPKE(NCELET , NCEL   , NVAR    , NPHAS  ,
     &                  IPHASS , ICLIP  , IWARNI(IKIPH),
     &                  PROPCE , RTP    )
C
          ELSE
            DO IEL = 1, NCEL
              RTP(IEL,IKIPH) = -GRAND
              RTP(IEL,IEIPH) = -GRAND
            ENDDO
          ENDIF
C
          IF (ITURB(IPHAS).EQ.50) THEN
            IPHIPH = IPHI(IPHAS)
            IFBIPH = IFB (IPHAS)
            DO IEL = 1, NCEL
              RTP(IEL,IPHIPH) = 2.D0/3.D0
              RTP(IEL,IFBIPH) = 0.D0
            ENDDO
          ENDIF
C
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
          IR11IP = IR11(IPHAS)
          IR22IP = IR22(IPHAS)
          IR33IP = IR33(IPHAS)
          IR12IP = IR12(IPHAS)
          IR13IP = IR13(IPHAS)
          IR23IP = IR23(IPHAS)
          IEIPH  = IEP (IPHAS)
C
          IF (UREF(IPHAS).GE.0.D0) THEN
C
            TRII   = (0.02D0*UREF(IPHAS))**2
C
            DO IEL = 1, NCEL
              RTP(IEL,IR11IP) = TRII
              RTP(IEL,IR22IP) = TRII
              RTP(IEL,IR33IP) = TRII
              RTP(IEL,IR12IP) = 0.D0
              RTP(IEL,IR13IP) = 0.D0
              RTP(IEL,IR23IP) = 0.D0
              XXK = 0.5D0*(RTP(IEL,IR11IP)+
     &                     RTP(IEL,IR22IP)+RTP(IEL,IR33IP))
              RTP(IEL,IEIPH) = XXK**1.5D0*CMU/ALMAX(IPHAS)
            ENDDO
            ICLIP = 1
            IPHASS = IPHAS
            CALL CLPRIJ(NCELET , NCEL   , NVAR    , NPHAS  ,
     &                  IPHASS , ICLIP  ,
     &                  PROPCE , RTP    , RTP    )
C
          ELSE
C
            DO IEL = 1, NCEL
              RTP(IEL,IR11IP) = -GRAND
              RTP(IEL,IR22IP) = -GRAND
              RTP(IEL,IR33IP) = -GRAND
              RTP(IEL,IR12IP) = -GRAND
              RTP(IEL,IR13IP) = -GRAND
              RTP(IEL,IR23IP) = -GRAND
              RTP(IEL,IEIPH)  = -GRAND
            ENDDO
C
          ENDIF
C
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
          IKIPH   = IK  (IPHAS)
          IOMGIP  = IOMG(IPHAS)
C
          IF (UREF(IPHAS).GE.0.D0) THEN
C
            DO IEL = 1, NCEL
              RTP(IEL,IKIPH ) = 1.5D0*(0.02D0*UREF(IPHAS))**2
C     on utilise la formule classique eps=k**1.5/Cmu/ALMAX et omega=eps/Cmu/k
              RTP(IEL,IOMGIP) = RTP(IEL,IKIPH)**0.5D0/ALMAX(IPHAS)
            ENDDO
C     pas la peine de clipper, les valeurs sont forcement positives
C
          ELSE
C
            DO IEL = 1, NCEL
              RTP(IEL,IKIPH ) = -GRAND
              RTP(IEL,IOMGIP) = -GRAND
            ENDDO
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 6.  CLIPPING DES GRANDEURS SCALAIRES (SF K-EPS VOIR CI DESSUS)
C=======================================================================
C
      IF (NSCAL.GT.0) THEN
C
C    Clipping des scalaires non variance
        DO IIS = 1, NSCAL
          IF(ISCAVR(IIS).EQ.0) THEN
            ISCAL = IIS
            CALL CLPSCA
C           ===========
     &     ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &       PROPCE , RA     , RTP    )
          ENDIF
        ENDDO
C
C     Clipping des variances qui sont clippees sans recours au scalaire
C        associe
        DO IIS = 1, NSCAL
          IF(ISCAVR(IIS).NE.0.AND.ICLVFL(IIS).NE.1) THEN
            ISCAL = IIS
            CALL CLPSCA
C           ===========
     &     ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &       PROPCE , RA     , RTP    )
          ENDIF
        ENDDO
C
C     Clipping des variances qui sont clippees avec recours au scalaire
C        associe s'il est connu
        DO IIS = 1, NSCAL
          IF( ISCAVR(IIS).LE.NSCAL.AND.ISCAVR(IIS).GE.1.AND.
     &        ICLVFL(IIS).EQ.1 ) THEN
            ISCAL = IIS
            CALL CLPSCA
C           ===========
     &     ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &       PROPCE , RTP(1,ISCA(ISCAVR(IIS))) , RTP      )
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 7.  INITIALISATION DE CONDITIONS AUX LIMITES ET FLUX DE MASSE
C      NOTER QUE LES CONDITIONS AUX LIMITES PEUVENT ETRE UTILISEES DANS
C      PHYVAR, PRECLI
C=======================================================================
C
C Conditions aux limites
      DO II = 1, NCOFAB
        DO IFAC = 1, NFABOR
          COEFA(IFAC,II) = 0.D0
          COEFB(IFAC,II) = 1.D0
        ENDDO
      ENDDO
C
      DO IPHAS = 1, NPHAS
        DO IFAC = 1, NFABOR
          IA(IITYPF-1+IFAC+NFABOR*(IPHAS-1)) = 0
          IA(IITRIF-1+IFAC+NFABOR*(IPHAS-1)) = 0
        ENDDO
      ENDDO
C
C Type symétrie : on en a besoin dans le cas du calcul des gradients
C     par moindres carrés étendu avec extrapolation du gradient au bord
C     La valeur 0 permet de ne pas extrapoler le gradient sur les faces.
C     Habituellement, on évite l'extrapolation sur les faces de symétries
C     pour ne pas tomber sur une indétermination et une matrice 3*3 non
C     inversible dans les configurations 2D).
      DO IPHAS = 1, NPHAS
        DO IFAC = 1, NFABOR
          IA(IISYMP-1+IFAC+NFABOR*(IPHAS-1)) = 0
        ENDDO
      ENDDO
C
C Flux de masse (on essaye de ne pas trop faire les choses 2 fois,
C  sans toutefois faire des tests trop compliques)
      IIFLUM = 0
      DO IVAR = 1, NVAR
        IF(IFLUMA(IVAR).GT.0.AND.IFLUMA(IVAR).NE.IIFLUM) THEN
          IIFLUM = IFLUMA(IVAR)
          DO IFAC = 1, NFAC
            PROPFA(IFAC,IPPROF(IIFLUM)) = 0.D0
          ENDDO
          DO IFAC = 1, NFABOR
            PROPFB(IFAC,IPPROB(IIFLUM)) = 0.D0
          ENDDO
        ENDIF
      ENDDO
C
C Flux de masse "ancien" : on utilise le fait  que IFLUMAA = -1 si
C     le flux ancien n'est pas defini (voir le test dans varpos).
C
      IIFLUA = 0
      DO IVAR = 1, NVAR
        IF(IFLUAA(IVAR).GT.0.AND.IFLUAA(IVAR).NE.IIFLUA) THEN
          IIFLUA = IFLUAA(IVAR)
          IIFLUM = IFLUMA(IVAR)
          DO IFAC = 1, NFAC
            PROPFA(IFAC,IPPROF(IIFLUA)) = PROPFA(IFAC,IPPROF(IIFLUM))
          ENDDO
          DO IFAC = 1, NFABOR
            PROPFB(IFAC,IPPROB(IIFLUA)) = PROPFB(IFAC,IPPROB(IIFLUM))
          ENDDO
        ENDIF
      ENDDO
C
C
C=======================================================================
C 8.  INITIALISATION DES TERMES SOURCES SI EXTRAPOLES
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
C     les termes sources de Navier Stokes
        IF(ISNO2T(IPHAS).GT.0) THEN
          IPTSNA = IPPROC(ITSNSA(IPHAS))
          DO II = 1, NDIM
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSNA+II-1) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
C     les termes sources turbulents
        IF(ISTO2T(IPHAS).GT.0) THEN
          IF(ITYTUR(IPHAS).EQ.2) JJ = 2
          IF(ITYTUR(IPHAS).EQ.3) JJ = 7
          IF(ITURB(IPHAS).EQ.50) JJ = 4
          IF(ITURB(IPHAS).EQ.60) JJ = 2
          IPTSTA = IPPROC(ITSTUA(IPHAS))
          DO II = 1, JJ
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSTA+II-1) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
      ENDDO
C
C     les termes sources des scalaires
      DO IIS = 1, NSCAL
        IF(ISSO2T(IIS).GT.0) THEN
          IPTSCA = IPPROC(ITSSCA(IIS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSCA) = 0.D0
          ENDDO
        ENDIF
      ENDDO
C
C=======================================================================
C 9.  INITIALISATION DES MOYENNES
C=======================================================================
C
      DO IMOM = 1, NBMOMT
        DO IEL = 1, NCEL
          PROPCE(IEL,IPPROC(ICMOME(IMOM))) = 0.D0
        ENDDO
      ENDDO
      DO II = 1,  NBDTCM
        DO IEL = 1, NCEL
          PROPCE(IEL,IPPROC(ICDTMO(II))) = 0.D0
        ENDDO
      ENDDO
      DO II = 1,  NBMOMX
        DTCMOM(II) = 0.D0
      ENDDO
C
C=======================================================================
C 10.  INITIALISATION CONSTANTE DE SMAGORINSKY EN MODELE DYNAMIQUE
C=======================================================================
C
      DO IPHAS = 1, NPHAS
        IF(ITURB(IPHAS).EQ.41) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPPROC(ISMAGO(IPHAS))) = 0.D0
          ENDDO
        ENDIF
      ENDDO
C
C=======================================================================
C 11.  INITIALISATION DU NUMERO DE LA FACE DE PAROI 5 LA PLUS PROCHE
C=======================================================================
C
C     Si IFAPA existe,
C     on suppose qu'il faut le (re)calculer : on init le tab a -1.
C
      DO IPHAS = 1, NPHAS
        IF(IIFAPA(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            IA(IIFAPA(IPHAS)-1+ IEL) = -1
          ENDDO
        ENDIF
      ENDDO
C
C=======================================================================
C 12.  INITIALISATION DE LA FORCE EXTERIEURE QUAND IPHYDR=1
C=======================================================================
C
      IF(IPHYDR.EQ.1) THEN
        DO IPHAS = 1, NPHAS
          DO IEL = 1, NCEL
            FRCXT(IEL,1,IPHAS) = 0.D0
            FRCXT(IEL,2,IPHAS) = 0.D0
            FRCXT(IEL,3,IPHAS) = 0.D0
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
C 12.  INITIALISATIONS EN ALE
C=======================================================================
C
      IF (IALE.EQ.1) THEN
        DO II = 1, NNOD
          IA(IIMPAL+II-1) = 0
          DO IDIM = 1, 3
            RA(IDEPAL+(IDIM-1)*NNOD+II-1) = 0.D0
            RA(IXYZN0+(II-1)*NDIM+IDIM-1) = XYZNOD(IDIM,II)
          ENDDO
        ENDDO
      ENDIF
C
C----
C FIN
C----
C
      RETURN
      END
c@z
