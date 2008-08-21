c@a
c@versb
C-----------------------------------------------------------------------
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
                        SUBROUTINE VERICL
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC VERIFICATION DE ICODCL
CFONC
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
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IVAR, ICODE
      INTEGER          NSTONI        , NSTVIT(NPHSMX), NSTOPP(NPHSMX)
      INTEGER          NSTOKE(NPHSMX), NSTOSC, NSTOVF
      INTEGER          NSTUVW(NPHSMX), NSTOUP(NPHSMX), NSTUKE(NPHSMX)
      INTEGER          NSTRIJ(NPHSMX), NSURIJ(NPHSMX), NSTOV2(NPHSMX)
      INTEGER          NSTUV2(NPHSMX), NSTOKW(NPHSMX), NSTUKW(NPHSMX)
      INTEGER          NSTUSC
      INTEGER          IIS, ICODCU, ICODCV, ICODCW, ICODCK, ICODCE
      INTEGER          ICODCP, ICODCF, ICODOM
      INTEGER          ICOR11, ICOR22, ICOR33, ICOR12, ICOR13, ICOR23
      INTEGER          IPP, IOKCOD, IOK, IPHAS
      INTEGER          IPRIPH, IUIPH , IVIPH , IWIPH , IKIPH , IEPIPH
      INTEGER          IPHIPH, IFBIPH, IOMGIP
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IPPPRP, IPPUIP, IPPVIP, IPPWIP, IPPEPP, IPPKIP
      INTEGER          IPP11P, IPP22P, IPP33P, IPP12P, IPP13P, IPP23P
      INTEGER          IPPPHP, IPPFBP, IPPOMG
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2.  VERIFICATION DE LA CONSISTANCE DES CL
C=======================================================================
C
C 2.1  INITIALISATION
C ====================
C
C Dans USCLIM, on se donne une grande liberte pour la specif. des c.l
C  sur les variables. neanmoins pour limiter la plage des tests, on se
C  donne, pour l'instant, les contraintes suivantes :
C
C   - meme type de c.l pour les 3 composantes de vitesse
C   - pas de conditions de frottemt sur la pression
C   - coherence entre les c.l vitesses et pression
C   - coherence entre les c.l vitesses et turbulence
C
      NSTONI = 0
      NSTOSC = 0
      NSTOVF = 0
      NSTUSC = 0
      DO IPHAS = 1, NPHAS
        NSTVIT(IPHAS) = 0
        NSTOPP(IPHAS) = 0
        NSTOKE(IPHAS) = 0
        NSTRIJ(IPHAS) = 0
        NSTOV2(IPHAS) = 0
        NSTOKW(IPHAS) = 0
        NSTUVW(IPHAS) = 0
        NSTOUP(IPHAS) = 0
        NSTUKE(IPHAS) = 0
        NSURIJ(IPHAS) = 0
        NSTUV2(IPHAS) = 0
        NSTUKW(IPHAS) = 0
      ENDDO
C
C
C 2.2 VERIFICATIONS QUE TOUTES LES CL SONT INITIALISEES
C ======================================================
C
C --- Premiere boucle rapide
      IOKCOD = 0
      DO IVAR = 1, NVAR
        DO IFAC = 1, NFABOR
          ICODE = ICODCL(IFAC,IVAR)
          IF(ICODE.EQ. 0) THEN
            IOKCOD = 1
          ENDIF
        ENDDO
      ENDDO
C
C --- Seconde boucle lente si pb plus haut
      IF(IOKCOD.NE.0) THEN
        DO IPP = 2, NVPPMX
          IF (ITRSVR(IPP).GE.1) THEN
            IVAR = ITRSVR(IPP)
            DO IFAC = 1, NFABOR
              ICODE = ICODCL(IFAC,IVAR)
              IF(ICODE.EQ. 0) THEN
                CHAINE=NOMVAR(IPP)
                WRITE(NFECRA,1000)IFAC,IPRFML(IFMFBR(IFAC),1),
     &                            CHAINE(1:8),ICODCL(IFAC,IVAR)
                NSTONI = NSTONI + 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
C
C 2.3 VERIFICATIONS DE L'ADMISSIBILITE DES CONDITIONS
C ====================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C --- Reperage des variables dans RTP
        IPRIPH = IPR(IPHAS)
        IUIPH  = IU (IPHAS)
        IVIPH  = IV (IPHAS)
        IWIPH  = IW (IPHAS)
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IKIPH  = IK (IPHAS)
          IEPIPH = IEP(IPHAS)
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IR11IP = IR11(IPHAS)
          IR22IP = IR22(IPHAS)
          IR33IP = IR33(IPHAS)
          IR12IP = IR12(IPHAS)
          IR13IP = IR13(IPHAS)
          IR23IP = IR23(IPHAS)
          IEPIPH = IEP(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IKIPH  = IK  (IPHAS)
          IEPIPH = IEP (IPHAS)
          IPHIPH = IPHI(IPHAS)
          IFBIPH = IFB(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IKIPH  = IK  (IPHAS)
          IOMGIP = IOMG(IPHAS)
        ENDIF
C
        IPPPRP = IPPRTP(IPRIPH)
        IPPUIP = IPPRTP(IUIPH )
        IPPVIP = IPPRTP(IVIPH )
        IPPWIP = IPPRTP(IWIPH )
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IPPKIP = IPPRTP(IKIPH )
          IPPEPP = IPPRTP(IEPIPH)
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IPP11P = IPPRTP(IR11IP)
          IPP22P = IPPRTP(IR22IP)
          IPP33P = IPPRTP(IR33IP)
          IPP12P = IPPRTP(IR12IP)
          IPP13P = IPPRTP(IR13IP)
          IPP23P = IPPRTP(IR23IP)
          IPPEPP = IPPRTP(IEPIPH)
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IPPKIP = IPPRTP(IKIPH )
          IPPEPP = IPPRTP(IEPIPH)
          IPPPHP = IPPRTP(IPHIPH)
          IPPFBP = IPPRTP(IFBIPH)
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IPPKIP = IPPRTP(IKIPH )
          IPPOMG = IPPRTP(IOMGIP)
        ENDIF
C
C --- Conditions admissibles pour les composantes de vitesse
        DO IFAC = 1, NFABOR
C
          ICODCU = ICODCL(IFAC,IUIPH)
          ICODCV = ICODCL(IFAC,IVIPH)
          ICODCW = ICODCL(IFAC,IWIPH)
C
          IF(ICODCU.NE. 1.AND.                 ICODCU.NE. 3.AND.
     &       ICODCU.NE. 4.AND.ICODCU.NE. 5.AND.ICODCU.NE. 6.AND.
     &                                         ICODCU.NE. 9) THEN
              CHAINE=NOMVAR(IPPUIP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IUIPH)
              NSTVIT(IPHAS) = NSTVIT(IPHAS) + 1
          ENDIF
          IF(ICODCV.NE. 1.AND.                 ICODCV.NE. 3.AND.
     &       ICODCV.NE. 4.AND.ICODCV.NE. 5.AND.ICODCV.NE. 6.AND.
     &                                         ICODCV.NE. 9) THEN
              CHAINE=NOMVAR(IPPVIP )
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IVIPH)
              NSTVIT(IPHAS) = NSTVIT(IPHAS) + 1
          ENDIF
          IF(ICODCW.NE. 1.AND.                 ICODCW.NE. 3.AND.
     &       ICODCW.NE. 4.AND.ICODCW.NE. 5.AND.ICODCV.NE. 6.AND.
     &                                         ICODCW.NE. 9) THEN
              CHAINE=NOMVAR(IPPWIP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IWIPH)
              NSTVIT(IPHAS) = NSTVIT(IPHAS) + 1
          ENDIF
C
C --- verification que la rugosite est initialisee si icodl=6
          IF(ICODCU.EQ.6 .AND. RCODCL(IFAC,IUIPH,3).LT.EPZERO)THEN
              CHAINE='RUGOSITV'
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IUIPH)
              NSTVIT(IPHAS) = NSTVIT(IPHAS) + 1
          ENDIF
C
C --- on interdit les parois rugueuses en compressible
          IF (ICODCU.EQ.6 .AND. IPPMOD(ICOMPF).GT.0) THEN
              CHAINE=NOMVAR(IPPUIP)
            WRITE(NFECRA,1015)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &           ICODCL(IFAC,IUIPH),IPPMOD(ICOMPF)
            NSTVIT(IPHAS) = NSTVIT(IPHAS) + 1
          ENDIF
C
        ENDDO
C
C --- Conditions admissibles pour la pression
        DO IFAC = 1, NFABOR
C
          IF(ICODCL(IFAC,IPRIPH).NE. 1.AND.
     &       ICODCL(IFAC,IPRIPH).NE. 3) THEN
              CHAINE=NOMVAR(IPPPRP)
            WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                        ICODCL(IFAC,IPRIPH)
            NSTOPP(IPHAS) = NSTOPP(IPHAS) + 1
          ENDIF
C
        ENDDO
C
C --- Conditions admissibles pour k et epsilon
        IF (ITYTUR(IPHAS).EQ.2) THEN
C
          DO IFAC = 1, NFABOR
C
            IF((ICODCL(IFAC,IKIPH ).NE. 1.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 3.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 5.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 6     ).OR.
     &         (ICODCL(IFAC,IEPIPH).NE. 1.AND.
     &          ICODCL(IFAC,IEPIPH).NE. 3.AND.
     &          ICODCL(IFAC,IEPIPH).NE. 5.AND.
     &          ICODCL(IFAC,IEPIPH).NE. 6     ) )THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IKIPH )
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IEPIPH)
              NSTOKE(IPHAS) = NSTOKE(IPHAS) + 1
            ENDIF
C
          ENDDO
C
C --- Conditions admissibles pour Rij et epsilon
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
          IVAR = IR11IP
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IVAR)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &         ICODE.NE. 4.AND.ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPP11P)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
          IVAR = IR22IP
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IVAR)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &         ICODE.NE. 4.AND.ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPP22P)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
          IVAR = IR33IP
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IVAR)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &         ICODE.NE. 4.AND.ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPP33P)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
          IVAR = IR12IP
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IVAR)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &         ICODE.NE. 4.AND.ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPP12P)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
          IVAR = IR13IP
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IVAR)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &         ICODE.NE. 4.AND.ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPP13P)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
          IVAR = IR23IP
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IVAR)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &         ICODE.NE. 4.AND.ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPP23P)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
          DO IFAC = 1, NFABOR
            ICODE = ICODCL(IFAC,IEPIPH)
            IF(ICODE.NE. 1.AND.                ICODE.NE. 3.AND.
     &                         ICODE.NE. 5.AND.ICODE.NE. 6     ) THEN
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),ICODE
              NSTRIJ(IPHAS) = NSTRIJ(IPHAS) + 1
            ENDIF
          ENDDO
C
C --- Conditions admissibles pour k, epsilon, phi et f_barre
        ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
          DO IFAC = 1, NFABOR
C
            IF((ICODCL(IFAC,IKIPH ).NE. 1.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 3.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 5.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 6     ).OR.
     &         (ICODCL(IFAC,IEPIPH).NE. 1.AND.
     &          ICODCL(IFAC,IEPIPH).NE. 3.AND.
     &          ICODCL(IFAC,IEPIPH).NE. 5.AND.
     &          ICODCL(IFAC,IEPIPH).NE. 6     ).OR.
     &         (ICODCL(IFAC,IPHIPH).NE. 1.AND.
     &          ICODCL(IFAC,IPHIPH).NE. 3.AND.
     &          ICODCL(IFAC,IPHIPH).NE. 5.AND.
     &          ICODCL(IFAC,IPHIPH).NE. 6     ).OR.
     &         (ICODCL(IFAC,IFBIPH).NE. 1.AND.
     &          ICODCL(IFAC,IFBIPH).NE. 3.AND.
     &          ICODCL(IFAC,IFBIPH).NE. 5.AND.
     &          ICODCL(IFAC,IFBIPH).NE. 6     ) )THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IKIPH )
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IEPIPH)
              CHAINE=NOMVAR(IPPPHP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IPHIPH )
              CHAINE=NOMVAR(IPPFBP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IFBIPH)
              NSTOV2(IPHAS) = NSTOV2(IPHAS) + 1
C
            ENDIF
C
          ENDDO
C
C --- Conditions admissibles pour k et omega
        ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
          DO IFAC = 1, NFABOR
C
            IF((ICODCL(IFAC,IKIPH ).NE. 1.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 3.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 5.AND.
     &          ICODCL(IFAC,IKIPH ).NE. 6     ).OR.
     &         (ICODCL(IFAC,IOMGIP).NE. 1.AND.
     &          ICODCL(IFAC,IOMGIP).NE. 3.AND.
     &          ICODCL(IFAC,IOMGIP).NE. 5.AND.
     &          ICODCL(IFAC,IOMGIP).NE. 6     ) )THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IKIPH )
              CHAINE=NOMVAR(IPPOMG)
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IOMGIP)
              NSTOKW(IPHAS) = NSTOKW(IPHAS) + 1
            ENDIF
C
          ENDDO
C
        ENDIF
C
      ENDDO
C --- Boucle sur les phases : fin
C
C --- Conditions admissibles pour les scalaires
      IF(NSCAL.GE.1) THEN
        DO IIS = 1,NSCAL
          IVAR = ISCA(IIS)
          DO IFAC = 1, NFABOR
            IF(ICODCL(IFAC,IVAR).NE. 1.AND.
     &         ICODCL(IFAC,IVAR).NE. 3.AND.
     &         ICODCL(IFAC,IVAR).NE. 5.AND.
     &         ICODCL(IFAC,IVAR).NE. 6 ) THEN
              CHAINE=NOMVAR(IPPRTP(IVAR))
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &          ICODCL(IFAC,IVAR)
              NSTOSC = NSTOSC + 1
            ENDIF
            IF(ICODCL(IFAC,IVAR).EQ. 5.AND.
     &         ISCAVR(IIS).GT.0        ) THEN
              CHAINE=NOMVAR(IPPRTP(IVAR))
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &          ICODCL(IFAC,IVAR)
              NSTOVF = NSTOVF + 1
            ENDIF
            IF(ICODCL(IFAC,IVAR).EQ. 6.AND.
     &         ISCAVR(IIS).GT.0        ) THEN
              CHAINE=NOMVAR(IPPRTP(IVAR))
              WRITE(NFECRA,1010)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &          ICODCL(IFAC,IVAR)
              NSTOVF = NSTOVF + 1
            ENDIF
C --- verification que la rugosite scalaire est initialisee si icodl=6
            IF(ICODCL(IFAC,IVAR).EQ.6.AND.
     &         RCODCL(IFAC,IVAR,3).LT.EPZERO)THEN
              CHAINE='RUGOSITS'
              WRITE(NFECRA,1010)IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),
     &                          ICODCL(IFAC,IVAR)
              NSTOSC = NSTOSC + 1
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C 2.4 VERIFICATIONS DES COHERENCES INTER VARIABLES INTRA PHASE
C =============================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C --- Reperage des variables dans RTP
        IPRIPH = IPR(IPHAS)
        IUIPH  = IU (IPHAS)
        IVIPH  = IV (IPHAS)
        IWIPH  = IW (IPHAS)
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IKIPH  = IK (IPHAS)
          IEPIPH = IEP(IPHAS)
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IR11IP = IR11(IPHAS)
          IR22IP = IR22(IPHAS)
          IR33IP = IR33(IPHAS)
          IR12IP = IR12(IPHAS)
          IR13IP = IR13(IPHAS)
          IR23IP = IR23(IPHAS)
          IEPIPH = IEP(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IKIPH  = IK  (IPHAS)
          IEPIPH = IEP (IPHAS)
          IPHIPH = IPHI(IPHAS)
          IFBIPH = IFB(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IKIPH  = IK  (IPHAS)
          IOMGIP = IOMG(IPHAS)
        ENDIF
C
        IPPPRP = IPPRTP(IPRIPH)
        IPPUIP = IPPRTP(IUIPH )
        IPPVIP = IPPRTP(IVIPH )
        IPPWIP = IPPRTP(IWIPH )
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IPPKIP = IPPRTP(IKIPH )
          IPPEPP = IPPRTP(IEPIPH)
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IPP11P = IPPRTP(IR11IP)
          IPP22P = IPPRTP(IR22IP)
          IPP33P = IPPRTP(IR33IP)
          IPP12P = IPPRTP(IR12IP)
          IPP13P = IPPRTP(IR13IP)
          IPP23P = IPPRTP(IR23IP)
          IPPEPP = IPPRTP(IEPIPH)
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IPPKIP = IPPRTP(IKIPH )
          IPPEPP = IPPRTP(IEPIPH)
          IPPPHP = IPPRTP(IPHIPH)
          IPPFBP = IPPRTP(IFBIPH)
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IPPKIP = IPPRTP(IKIPH )
          IPPOMG = IPPRTP(IOMGIP)
        ENDIF
C
C --- Coherence pour les composantes de vitesse
        DO IFAC = 1, NFABOR
C
          ICODCU = ICODCL(IFAC,IUIPH)
          ICODCV = ICODCL(IFAC,IVIPH)
          ICODCW = ICODCL(IFAC,IWIPH)
C
          IF(ICODCU.EQ.4.OR.ICODCU.EQ.5.OR.ICODCU.EQ.6.OR.
     &                      ICODCU.EQ.9.OR.
     &       ICODCV.EQ.4.OR.ICODCV.EQ.5.OR.ICODCV.EQ.6.OR.
     &                      ICODCV.EQ.9.OR.
     &       ICODCW.EQ.4.OR.ICODCW.EQ.5.OR.ICODCW.EQ.6.OR.
     &                      ICODCW.EQ.9                )THEN
C
            IF( ICODCU.NE.ICODCV .OR. ICODCU.NE.ICODCW .OR.
     &                                ICODCV.NE.ICODCW ) THEN
              WRITE(NFECRA,1020)IFAC,IPRFML(IFMFBR(IFAC),1),IPHAS,
     &                          ICODCU,ICODCV,ICODCW
              NSTUVW(IPHAS) = NSTUVW(IPHAS) + 1
            ENDIF
          ENDIF
C
C --- Coherence vitesse pression
C
C      Remarques :
C        Pas de regle stricte de coherence vitesse/pression.
C        Avant on imposait un Dirichlet sur la pression pour en
C        entree/sortie, mais cela ne semble pas imperatif. Le test
C        est laisse en commentaire pour etre recupere si necessaire.
C
C          IF( ICODCU.EQ.9 .OR. ICODCV.EQ.9 .OR. ICODCW.EQ.9 ) THEN
C            IF( ICODCL(IFAC,IPRIPH).NE.1                    ) THEN
C              CHAINE=NOMVAR(IPPPRP)
C              WRITE(NFECRA,1030)
C     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
C     &          ICODCL(IFAC,IPRIPH),ICODCU,ICODCV,ICODCW
C              NSTOUP(IPHAS) = NSTOUP(IPHAS) + 1
C            ENDIF
C          ENDIF
C
        ENDDO
C
C --- Coherence vitesse turbulence
C
        IF(ITYTUR(IPHAS).EQ.2) THEN
C
          DO IFAC = 1, NFABOR
C
            ICODCU = ICODCL(IFAC,IUIPH)
            ICODCV = ICODCL(IFAC,IVIPH)
            ICODCW = ICODCL(IFAC,IWIPH)
            ICODCK = ICODCL(IFAC,IKIPH)
            ICODCE = ICODCL(IFAC,IEPIPH)
C
            IF( (ICODCU.EQ.5 .OR. ICODCV.EQ.5 .OR. ICODCW.EQ.5 .OR.
     &           ICODCK.EQ.5 .OR. ICODCE.EQ.5) .AND.
     &          (ICODCU.NE.5 .OR. ICODCV.NE.5 .OR. ICODCW.NE.5 .OR.
     &           ICODCK.NE.5 .OR. ICODCE.NE.5)                    ) THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IKIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IEPIPH),ICODCU,ICODCV,ICODCW
              NSTUKE(IPHAS) = NSTUKE(IPHAS) + 1
            ENDIF
C
            IF( (ICODCU.EQ.6 .OR. ICODCV.EQ.6 .OR. ICODCW.EQ.6 .OR.
     &           ICODCK.EQ.6 .OR. ICODCE.EQ.6) .AND.
     &          (ICODCU.NE.6 .OR. ICODCV.NE.6 .OR. ICODCW.NE.6 .OR.
     &           ICODCK.NE.6 .OR. ICODCE.NE.6)                    ) THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IKIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IEPIPH),ICODCU,ICODCV,ICODCW
              NSTUKE(IPHAS) = NSTUKE(IPHAS) + 1
            ENDIF
C
          ENDDO
C
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
          DO IFAC = 1, NFABOR
C
            ICODCU = ICODCL(IFAC,IUIPH)
            ICODCV = ICODCL(IFAC,IVIPH)
            ICODCW = ICODCL(IFAC,IWIPH)
            ICOR11 = ICODCL(IFAC,IR11IP)
            ICOR22 = ICODCL(IFAC,IR22IP)
            ICOR33 = ICODCL(IFAC,IR33IP)
            ICOR12 = ICODCL(IFAC,IR12IP)
            ICOR13 = ICODCL(IFAC,IR13IP)
            ICOR23 = ICODCL(IFAC,IR23IP)
            ICODCE = ICODCL(IFAC,IEPIPH)
C
            IF(  (ICODCU.EQ.5 .OR. ICODCV.EQ.5 .OR. ICODCW.EQ.5 .OR.
     &            ICOR11.EQ.5 .OR. ICOR22.EQ.5 .OR.
     &            ICOR33.EQ.5 .OR. ICOR12.EQ.5 .OR.
     &            ICOR13.EQ.5 .OR. ICOR23.EQ.5 .OR.
     &            ICODCE.EQ.5                      ) .AND.
     &           (ICODCU.NE.5 .OR. ICODCV.NE.5 .OR. ICODCW.NE.5 .OR.
     &            ICOR11.NE.5 .OR. ICOR22.NE.5 .OR.
     &            ICOR33.NE.5 .OR. ICOR12.NE.5 .OR.
     &            ICOR13.NE.5 .OR. ICOR23.NE.5 .OR.
     &            ICODCE.NE.5                      )      ) THEN
              WRITE(NFECRA,1040)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),IPHAS,
     &          ICOR11,ICOR22,ICOR33,
     &          ICOR12,ICOR13,ICOR23,
     &          ICODCE,ICODCU,ICODCV,ICODCW
              NSURIJ(IPHAS) = NSURIJ(IPHAS) + 1
            ENDIF
C
            IF(  (ICODCU.EQ.6 .OR. ICODCV.EQ.6 .OR. ICODCW.EQ.6 .OR.
     &            ICOR11.EQ.6 .OR. ICOR22.EQ.6 .OR.
     &            ICOR33.EQ.6 .OR. ICOR12.EQ.6 .OR.
     &            ICOR13.EQ.6 .OR. ICOR23.EQ.6 .OR.
     &            ICODCE.EQ.6                      ) .AND.
     &           (ICODCU.NE.6 .OR. ICODCV.NE.6 .OR. ICODCW.NE.6 .OR.
     &            ICOR11.NE.6 .OR. ICOR22.NE.6 .OR.
     &            ICOR33.NE.6 .OR. ICOR12.NE.6 .OR.
     &            ICOR13.NE.6 .OR. ICOR23.NE.6 .OR.
     &            ICODCE.NE.6                      )      ) THEN
              WRITE(NFECRA,1040)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),IPHAS,
     &          ICOR11,ICOR22,ICOR33,
     &          ICOR12,ICOR13,ICOR23,
     &          ICODCE,ICODCU,ICODCV,ICODCW
              NSURIJ(IPHAS) = NSURIJ(IPHAS) + 1
            ENDIF
C
            IF(  (ICODCU.EQ.4 .OR. ICODCV.EQ.4 .OR. ICODCW.EQ.4 .OR.
     &            ICOR11.EQ.4 .OR. ICOR22.EQ.4 .OR.
     &            ICOR33.EQ.4 .OR. ICOR12.EQ.4 .OR.
     &            ICOR13.EQ.4 .OR. ICOR23.EQ.4
     &                                  ) .AND.
     &           (ICODCU.NE.4 .OR. ICODCV.NE.4 .OR. ICODCW.NE.4 .OR.
     &            ICOR11.NE.4 .OR. ICOR22.NE.4 .OR.
     &            ICOR33.NE.4 .OR. ICOR12.NE.4 .OR.
     &            ICOR13.NE.4 .OR. ICOR23.NE.4 .OR.
     &            ICODCE.NE.3) ) THEN
              WRITE(NFECRA,1040)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),IPHAS,
     &          ICOR11,ICOR22,ICOR33,
     &          ICOR12,ICOR13,ICOR23,
     &          ICODCE,ICODCU,ICODCV,ICODCW
              NSURIJ(IPHAS) = NSURIJ(IPHAS) + 1
            ENDIF
C
          ENDDO
C
        ELSEIF(ITURB(IPHAS).EQ.50 ) THEN
C
          DO IFAC = 1, NFABOR
C
            ICODCU = ICODCL(IFAC,IUIPH)
            ICODCV = ICODCL(IFAC,IVIPH)
            ICODCW = ICODCL(IFAC,IWIPH)
            ICODCK = ICODCL(IFAC,IKIPH)
            ICODCE = ICODCL(IFAC,IEPIPH)
            ICODCP = ICODCL(IFAC,IPHIPH)
            ICODCF = ICODCL(IFAC,IFBIPH)
C
            IF( (ICODCU.EQ.5 .OR. ICODCV.EQ.5 .OR. ICODCW.EQ.5 .OR.
     &           ICODCK.EQ.5 .OR. ICODCE.EQ.5 .OR. ICODCP.EQ.5 .OR.
     &           ICODCF.EQ.5 ) .AND.
     &          (ICODCU.NE.5 .OR. ICODCV.NE.5 .OR. ICODCW.NE.5 .OR.
     &           ICODCK.NE.5 .OR. ICODCE.NE.5 .OR. ICODCP.NE.5 .OR.
     &           ICODCF.NE.5 )                    ) THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IKIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IEPIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPPHP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IPHIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPFBP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IFBIPH),ICODCU,ICODCV,ICODCW
              NSTUV2(IPHAS) = NSTUV2(IPHAS) + 1
C
            IF( (ICODCU.EQ.6 .OR. ICODCV.EQ.6 .OR. ICODCW.EQ.6 .OR.
     &           ICODCK.EQ.6 .OR. ICODCE.EQ.6 .OR. ICODCP.EQ.6 .OR.
     &           ICODCF.EQ.6 ) .AND.
     &          (ICODCU.NE.6 .OR. ICODCV.NE.6 .OR. ICODCW.NE.6 .OR.
     &           ICODCK.NE.6 .OR. ICODCE.NE.6 .OR. ICODCP.NE.6 .OR.
     &           ICODCF.NE.6 )                    ) THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IKIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPEPP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IEPIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPPHP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IPHIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPFBP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IFBIPH),ICODCU,ICODCV,ICODCW
              NSTUV2(IPHAS) = NSTUV2(IPHAS) + 1
C
            ENDIF
C
            ENDIF
C
          ENDDO
C
        ELSEIF(ITURB(IPHAS).EQ.60 ) THEN
C
          DO IFAC = 1, NFABOR
C
            ICODCU = ICODCL(IFAC,IUIPH)
            ICODCV = ICODCL(IFAC,IVIPH)
            ICODCW = ICODCL(IFAC,IWIPH)
            ICODCK = ICODCL(IFAC,IKIPH)
            ICODOM = ICODCL(IFAC,IOMGIP)
C
            IF( (ICODCU.EQ.5 .OR. ICODCV.EQ.5 .OR. ICODCW.EQ.5 .OR.
     &           ICODCK.EQ.5 .OR. ICODOM.EQ.5 ) .AND.
     &          (ICODCU.NE.5 .OR. ICODCV.NE.5 .OR. ICODCW.NE.5 .OR.
     &           ICODCK.NE.5 .OR. ICODOM.NE.5 ) ) THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IKIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPOMG)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IOMGIP),ICODCU,ICODCV,ICODCW
              NSTUKW(IPHAS) = NSTUKW(IPHAS) + 1
            ENDIF
C
            IF( (ICODCU.EQ.6 .OR. ICODCV.EQ.6 .OR. ICODCW.EQ.6 .OR.
     &           ICODCK.EQ.6 .OR. ICODOM.EQ.6 ) .AND.
     &          (ICODCU.NE.6 .OR. ICODCV.NE.6 .OR. ICODCW.NE.6 .OR.
     &           ICODCK.NE.6 .OR. ICODOM.NE.6 ) ) THEN
              CHAINE=NOMVAR(IPPKIP)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IKIPH),ICODCU,ICODCV,ICODCW
              CHAINE=NOMVAR(IPPOMG)
              WRITE(NFECRA,1030)
     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8),IPHAS,
     &          ICODCL(IFAC,IOMGIP),ICODCU,ICODCV,ICODCW
              NSTUKW(IPHAS) = NSTUKW(IPHAS) + 1
            ENDIF
C
          ENDDO
C
        ENDIF
C
      ENDDO
C --- Boucle sur les phases : fin
C
C
C 2.5 VERIFICATIONS DES COHERENCES INTER VARIABLES INTRA PHASE
C =============================================================
C
C --- Coherence vitesse scalaires
C
      IF( NSCAL.GE.1 ) THEN
        DO IIS = 1, NSCAL
          IPHAS = IPHSCA(IIS)
          IF(ITYTUR(IPHAS).EQ.2.OR.ITYTUR(IPHAS).EQ.3) THEN
            IVAR  = ISCA(IIS)
            DO IFAC = 1, NFABOR
              ICODCU = ICODCL(IFAC,IU(IPHAS))
              IF(ICODCL(IFAC,IVAR).EQ.5.AND.ICODCU.NE.5) THEN
                CHAINE=NOMVAR(IPPRTP(IVAR))
                WRITE(NFECRA,1050) IFAC,IPRFML(IFMFBR(IFAC),1),
     &                        CHAINE(1:8), IIS, IPHAS,
     &                        ICODCL(IFAC,IVAR), ICODCU
                NSTUSC = NSTUSC + 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 3.  IMPRESSIONS RECAPITULATIVES
C=======================================================================
C
      IOK = 0
C
      IF( NSTONI.GT.0 .OR. NSTOSC.GT.0 .OR. NSTOVF.GT.0 .OR.
     &                                      NSTUSC.GT.0 ) THEN
        WRITE (NFECRA,1901) NSTONI, NSTOSC, NSTOVF, NSTUSC
        IOK = 1
      ENDIF
C
      DO IPHAS = 1, NPHAS
        IF( NSTVIT(IPHAS).GT.0 .OR. NSTOPP(IPHAS).GT.0 .OR.
     &      NSTOKE(IPHAS).GT.0 .OR. NSTRIJ(IPHAS).GT.0 .OR.
     &      NSTOV2(IPHAS).GT.0 .OR.
     &      NSTUVW(IPHAS).GT.0 .OR. NSTOUP(IPHAS).GT.0 .OR.
     &      NSTUKE(IPHAS).GT.0 .OR. NSURIJ(IPHAS).GT.0 .OR.
     &      NSTUV2(IPHAS).GT.0       ) THEN
          WRITE (NFECRA,1902) IPHAS, NSTVIT(IPHAS),NSTOPP(IPHAS),
     &                               NSTOKE(IPHAS),NSTRIJ(IPHAS),
     &                               NSTOV2(IPHAS),
     &                               NSTUVW(IPHAS),NSTOUP(IPHAS),
     &                               NSTUKE(IPHAS),NSURIJ(IPHAS),
     &                               NSTUV2(IPHAS)
          IOK = 1
        ENDIF
      ENDDO
C
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C 3.  FORMATS
C=======================================================================
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@ COND. LIM. NON INITIALISEES                                ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,
     &'@     ICODCL VARIABLE ', I10                                  ,/,
     &'@                                                            '  )
 1010 FORMAT(
     &'@                                                            ',/,
     &'@ COND. LIM. NON PREVUES                                     ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,
     &'@     ICODCL VARIABLE ', I10                                  ,/,
     &'@                                                            '  )
 1015 FORMAT(
     &'@                                                            ',/,
     &'@ CONDITIONS AUX LIMITES DE PAROI RUGUEUSE INCOMPATIBLES     ',/,
     &'@ AVEC LE MODULE COMPRESSIBLE                                ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,
     &'@     ICODCL VARIABLE =',I10                                  ,/,
     &'@     IPPMOD(ICOMPF)  =',I10                                  ,/,
     &'@                                                            '  )
 1020 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCE COND. LIM. COMPOSANTES DE LA VITESSE           ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; PHASE  ',I10     ,/,
     &'@     ICODCL VITESSE  ',3I10                                  ,/,
     &'@                                                            '  )
 1030 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCE COND. LIM. VITESSE-VARIABLE                    ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,
     &'@                                 PHASE ',I10                 ,/,
     &'@     ICODCL VARIABLE ', I10                                  ,/,
     &'@     ICODCL VITESSE  ',3I10                                  ,/,
     &'@                                                            '  )
 1040 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCE COND. LIM. VITESSE-RIJ-EPSILON                 ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; RIJ-EPSILON     ',/,
     &'@                                 PHASE ',I10                 ,/,
     &'@     ICODCL RIJ-EPS ',7I5                                    ,/,
     &'@     ICODCL VITESSE ',3I5                                    ,/,
     &'@                                                            '  )
 1050 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCE COND. LIM. VITESSE-SCALAIRE                    ',/,
     &'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,
     &'@     SCALAIRE NUMERO ',I10   ,'; PHASE ',I10                 ,/,
     &'@     ICODCL SCALAIRE ',I10   ,'; ICODCL VITESSE ',I10        ,/,
     &'@                                                            '  )
 1901 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@         Conditions aux limites non initialisees  : ',I10    ,/,
     &'@         Conditions aux limites non prevues :               ',/,
     &'@             sur les scalaires                    : ',I10    ,/,
     &'@             sur les scalaires representant                 ',/,
     &'@                                    une variance  : ',I10    ,/,
     &'@         Incoherences :                                     ',/,
     &'@             entre vitesse et scalaires           : ',I10    ,/,
     &'@                                                            ',/,
     &'@         Le calcul ne sera pas execute.                     ',/,
     &'@                                                            ',/,
     &'@         Verifier les parametres donnes via l''interface    ',/,
     &'@           ou usclim.                                       ',/,
     &'@                                                            ',/)
 1902 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@                 POUR LA PHASE ',I10                         ,/,
     &'@                                                            ',/,
     &'@         Conditions aux limites non prevues :               ',/,
     &'@             sur la vitesse                       : ',I10    ,/,
     &'@             sur la pression                      : ',I10    ,/,
     &'@             sur k et epsilon                     : ',I10    ,/,
     &'@             sur Rij et epsilon                   : ',I10    ,/,
     &'@             sur k, epsilon, phi et f_barre       : ',I10    ,/,
     &'@         Incoherences :                                     ',/,
     &'@             entre les composantes de la vitesse  : ',I10    ,/,
     &'@             entre vitesse et pression            : ',I10    ,/,
     &'@             entre vitesse et k-epsilon           : ',I10    ,/,
     &'@             entre vitesse et Rij-epsilon         : ',I10    ,/,
     &'@             entre vitesse et v2f                 : ',I10    ,/,
     &'@                                                            ',/,
     &'@         Le calcul ne sera pas execute.                     ',/,
     &'@                                                            ',/,
     &'@         Verifier les parametres donnes via l''interface    ',/,
     &'@           ou usclim.                                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@ UNINITIALIZED BOUNDARY CONDITIONS                          ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,
     &'@     ICODCL VARIABLE ', I10                                  ,/,
     &'@                                                            '  )
 1010 FORMAT(
     &'@                                                            ',/,
     &'@ UNEXPECTED BOUNDARY CONDITIONS                             ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,
     &'@     ICODCL VARIABLE ', I10                                  ,/,
     &'@                                                            '  )
 1015 FORMAT(
     &'@                                                            ',/,
     &'@ ROUGH WALL BOUNDARY CONDITIONS INCOMPATIBLE WITH THE       ',/,
     &'@ COMPRESSIBLE MODULE                                        ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,
     &'@     ICODCL VARIABLE =',I10                                  ,/,
     &'@     IPPMOD(ICOMPF)  =',I10                                  ,/,
     &'@                                                            '  )
 1020 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY COMPONENT         ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; PHASE  ',I10      ,/,
     &'@     ICODCL VELOCITY ',3I10                                  ,/,
     &'@                                                            '  )
 1030 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-VARIABLE          ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,
     &'@                                 PHASE ',I10                 ,/,
     &'@     ICODCL VARIABLE ', I10                                  ,/,
     &'@     ICODCL VELOCITY ',3I10                                  ,/,
     &'@                                                            '  )
 1040 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-RIJ-EPSILON       ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; RIJ-EPSILON      ',/,
     &'@                                 PHASE ',I10                 ,/,
     &'@     ICODCL RIJ-EPS  ',7I5                                   ,/,
     &'@     ICODCL VELOCITY ',3I5                                   ,/,
     &'@                                                            '  )
 1050 FORMAT(
     &'@                                                            ',/,
     &'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-SCALAR            ',/,
     &'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,
     &'@     SCALAR NUMBER ',I10   ,'; PHASE ',I10                   ,/,
     &'@     ICODCL SCALAR ',I10   ,'; ICODCL VELOCITY ',I10         ,/,
     &'@                                                            '  )
 1901 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS VERIF.    ',/,
     &'@    ********                                                ',/,
     &'@                                                            ',/,
     &'@         Uninitialized boundary conditions        : ',I10    ,/,
     &'@         Unexpected  boundary conditions:                   ',/,
     &'@             on the scalars                       : ',I10    ,/,
     &'@             on the scalars representing                    ',/,
     &'@                                      a variance  : ',I10    ,/,
     &'@         Incoherencies:                                     ',/,
     &'@             between velocity and scalars         : ',I10    ,/,
     &'@                                                            ',/,
     &'@         The calculation will not be run.                   ',/,
     &'@                                                            ',/,
     &'@         Verify the parameters given via the interface or   ',/,
     &'@           usclim.                                          ',/,
     &'@                                                            ',/)
 1902 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS VERIF.    ',/,
     &'@    ********                                                ',/,
     &'@                 FOR PHASE ',I10                             ,/,
     &'@                                                            ',/,
     &'@         Unexpeted boundary conditions:                     ',/,
     &'@             on the velocity                      : ',I10    ,/,
     &'@             on the pressure                      : ',I10    ,/,
     &'@             on k and epsilon                     : ',I10    ,/,
     &'@             on Rij and epsilon                   : ',I10    ,/,
     &'@             on k, epsilon, phi and f_barre       : ',I10    ,/,
     &'@         Incoherencies:                                     ',/,
     &'@             between the velocity components      : ',I10    ,/,
     &'@             between velocity and pressure        : ',I10    ,/,
     &'@             between velocity and k-epsilon       : ',I10    ,/,
     &'@             between velocity and Rij-epsilon     : ',I10    ,/,
     &'@             between velocity and v2f             : ',I10    ,/,
     &'@                                                            ',/,
     &'@         The calculation will not be run.                   ',/,
     &'@                                                            ',/,
     &'@         Verify the parameters given via the interface or   ',/,
     &'@           usclim.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
      RETURN
      END
c@z
