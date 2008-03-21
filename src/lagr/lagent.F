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
                        SUBROUTINE LAGENT
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYCEL , ICOCEL ,
     &   ITYPFB , ITRIFB , IFRLAG , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   SURFBN , DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , VAGAUS , AUXL   , W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC    Gestion de l'injection des particules dans le domaine de calcul
CFONC
CFONC      1. initialisation par l'utilisateur via USLAG2
CFONC         des classes de particules et du type d'interaction
CFONC         particule/face de frontiere.
CFONC
CFONC      2. injection des particules dans le domaine : initialisation
CFONC         des tableau ETTP, ITEPA(IP,JISOR) et TEPA(IP,JRPOI).
CFONC
CFONC      3. modification des conditions d'injection des particules :
CFONC         retouche des ETTP, ITEPA(IP,JISOR) et TEPA(IP,JRPOI).
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
CARGU ! LNDNOD       ! E  !  -> ! LONGUEUR DU TABLEAU ICOCEL
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
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
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! ICOCEL       ! TE !  -> ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (LNDNOD)     !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE !  -> ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (NCELET+1)   !    !     !    POINTEUR DU TABLEAU ICOCEL        !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE !  -> ! TAB D'INDIRECTION POUR TRI DES FACES !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! IFRLAG       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE LAGRANGIEN           !
CARGU ! ITEPA        ! TE ! <-  ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! SURFBN(NFABOR! TR !  -> ! SURFACE DES FACES DE BORD            !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT PREC OU         !
CARGU !              !    !     !    INSTANT COURANT SI NTCABS = 1)    !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR ! <-> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! VAGAUS       ! TR ! <-  ! VARIABLES ALEATOIRES GAUSSIENNES     !
CARGU !(NBPMAX,NVGAUS!    !     !                                      !
CARGU ! AUXL(NBPMAX,3! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! W1..W3(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
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
C
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
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "radiat.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
      INTEGER          ICOCEL(LNDNOD) , ITYCEL(NCELET+1)
      INTEGER          ITEPA(NBPMAX,NIVEP) , IFRLAG(NFABOR)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION SURFBN(NFABOR)
      DOUBLE PRECISION DT(NCELET) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION VAGAUS(NBPMAX,*)
      DOUBLE PRECISION AUXL(NBPMAX,3)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFINIA, IFINRA
C
      INTEGER          IEL , IFAC , IPHAS , IP , NB , NC, II, IFVU
      INTEGER          IIWORK , IOK , N1 , ND , ICHA
      INTEGER          NPT , NFIN , NPAR1  , NPAR2 , MODE , IDVAR
      INTEGER          MAXELT, IDBIA1, ILS
C
      DOUBLE PRECISION VN1 , VN2 , VN3 , PIS6 , D3
      DOUBLE PRECISION DMASSE , RD(1) , AA
      DOUBLE PRECISION XXPART , YYPART , ZZPART
      DOUBLE PRECISION TVPART , UUPART , VVPART , WWPART
      DOUBLE PRECISION DDPART , TTPART
      DOUBLE PRECISION SURF   , VOLP , VITP
C
C***********************************************************************
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C Init aberrante pour forcer l'utilisateur a mettre sa valeur
C
      IPHAS = ILPHAS
C
      PIS6 = PI / 6.D0
C
      DO NB = 1, NFLAGM
        IUSNCL(NB) = 0
        IUSCLB(NB) = 0
      ENDDO
C
      DO NB = 1, NFLAGM
        DO NC = 1, NCLAGM
          DO ND = 1,NDLAGM
            RUSLAG(NC,NB,ND) = 0.D0
          ENDDO
          DO ND = 1,NDLAIM
            IUSLAG(NC,NB,ND) = 0
          ENDDO
        ENDDO
      ENDDO
C
      DO NB = 1, NFLAGM
        DO NC = 1, NCLAGM
C
          IUSLAG(NC,NB,IJNBP) =  0
          IUSLAG(NC,NB,IJFRE) =  0
          IUSLAG(NC,NB,ICLST) =  0
          IUSLAG(NC,NB,IJUVW) = -2
          IUSLAG(NC,NB,IJPRTP) = -2
          IUSLAG(NC,NB,IJPRDP) = -2
          IUSLAG(NC,NB,IJPRPD) = -2
C
          RUSLAG(NC,NB,IUNO)  = -GRAND
          RUSLAG(NC,NB,IUPT)  = -GRAND
          RUSLAG(NC,NB,IVPT)  = -GRAND
          RUSLAG(NC,NB,IWPT)  = -GRAND
          RUSLAG(NC,NB,IPOIT) = -GRAND
C
          RUSLAG(NC,NB,IDPT)  = -GRAND
          RUSLAG(NC,NB,IVDPT) = -GRAND
          RUSLAG(NC,NB,IROPT) = -GRAND
C
          IF ( IPHYLA.EQ.1 ) THEN
C
C    Thermique
C
            IF ( ITPVAR.EQ.1 ) THEN
              RUSLAG(NC,NB,ITPT)  = -GRAND
              RUSLAG(NC,NB,ICPT)  = -GRAND
              RUSLAG(NC,NB,IEPSI) = -GRAND
            ENDIF
C
          ELSE IF ( IPHYLA .EQ. 2 ) THEN
C
C    Charbon
C    (REM : le diametre du coeur retrecissant est calcule dans lagich.F)
C
            IUSLAG(NC,NB,INUCHL) = 0
            RUSLAG(NC,NB,IHPT)   = -GRAND
            RUSLAG(NC,NB,IMCHT)  = -GRAND
            RUSLAG(NC,NB,IMCKT)  = -GRAND
            RUSLAG(NC,NB,ICPT)   = -GRAND
C
          ENDIF
C
        ENDDO
      ENDDO
C
      DO IFAC = 1,NFABOR
        IFRLAG(IFAC) = 0
      ENDDO
C
C    Mise a zero des debits pour chaque zone de bord
C
      DO NB = 1,NFLAGM
        DEBLAG(NB) = 0.D0
      ENDDO
C
C=======================================================================
C 2. Initialisation utilisateur par classe et par frontiere
C=======================================================================
C
      MAXELT = MAX(NCELET,NFAC,NFABOR)
      ILS    = IDEBIA
      IDBIA1 = ILS + MAXELT
      CALL IASIZE('LAGENT',IDBIA1)

      CALL USLAG2
C     ===========
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 3. Controles
C=======================================================================
C
      IOK = 0
C
C --> Les faces doivent toutes appartenir a une zone frontiere
C
      DO IFAC = 1, NFABOR
        IF(IFRLAG(IFAC).LE.0 .OR. IFRLAG(IFAC).GT.NFLAGM) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1000) IFAC,NFLAGM,IFRLAG(IFAC)
        ENDIF
      ENDDO
C
      IF (IOK.GT.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C --> On construit une liste des numeros des zones frontieres.
C
      NFRLAG = 0
      DO IFAC = 1, NFABOR
        IFVU = 0
        DO II = 1, NFRLAG
          IF (ILFLAG(II).EQ.IFRLAG(IFAC)) THEN
            IFVU = 1
          ENDIF
        ENDDO
        IF(IFVU.EQ.0) THEN
          NFRLAG = NFRLAG + 1
          IF(NFRLAG.LE.NFLAGM) THEN
            ILFLAG(NFRLAG) = IFRLAG(IFAC)
          ELSE
            WRITE(NFECRA,1001) NFRLAG
            WRITE(NFECRA,'(I10)') (ILFLAG(II),II=1,NFRLAG)
            CALL CSEXIT (1)
C           ===========
          ENDIF
        ENDIF
      ENDDO
C
C --> Nombre de classes.
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        IF (IUSNCL(NB).LT.0 .OR. IUSNCL(NB).GT.NCLAGM) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1010) NB,NCLAGM,IUSNCL(NB)
        ENDIF
      ENDDO
C
C --> Nombre de particules.
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        DO NC = 1,IUSNCL(NB)
          IF ( IUSLAG(NC,NB,IJNBP).LT.0 .OR.
     &         IUSLAG(NC,NB,IJNBP).GT.NBPMAX ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1020) NC,NB,NBPMAX,IUSLAG(NC,NB,IJNBP)
          ENDIF
        ENDDO
      ENDDO
C
C    Verification des classes de particules : Uniquement un warning
C
      IF (NBCLST.GT.0 ) THEN
        DO II = 1,NFRLAG
          NB = ILFLAG(II)
          DO NC = 1,IUSNCL(NB)
            IF (IUSLAG(NC,NB,ICLST).LE.0      .OR.
     &          IUSLAG(NC,NB,ICLST).GT.NBCLST        ) THEN
              WRITE(NFECRA,1021) NC,NB,NBCLST,IUSLAG(NC,NB,ICLST)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C --> Frequence d'injection.
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (IUSLAG(NC,NB,IJFRE).LT.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1030) NB,NC,IUSLAG(NC,NB,IJFRE)
          ENDIF
        ENDDO
      ENDDO
C
C --> Conditions au bord.
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        IF ( IUSCLB(NB).NE.IENTRL .AND. IUSCLB(NB).NE.ISORTL .AND.
     &       IUSCLB(NB).NE.IREBOL .AND. IUSCLB(NB).NE.IDEPO1 .AND.
     &       IUSCLB(NB).NE.IDEPO2 .AND. IUSCLB(NB).NE.IDEPO3 .AND.
     &       IUSCLB(NB).NE.IENCRL .AND. IUSCLB(NB).NE.JBORD1 .AND.
     &       IUSCLB(NB).NE.JBORD2 .AND. IUSCLB(NB).NE.JBORD3 .AND.
     &       IUSCLB(NB).NE.JBORD4 .AND. IUSCLB(NB).NE.JBORD5 .AND.
     &       IUSCLB(NB).NE.IDEPFA ) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1040) NB
        ENDIF
      ENDDO
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        IF (IUSCLB(NB).EQ.IDEPO2 .AND. IENSI2.NE.1) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1041) NB
        ENDIF
      ENDDO
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        IF (IUSCLB(NB).EQ.IENCRL .AND. IPHYLA.NE.2) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1042) NB
        ENDIF
      ENDDO
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        IF (IUSCLB(NB).EQ.IENCRL .AND. IENCRA.NE.1) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1043) NB
        ENDIF
      ENDDO
C
C --> Type de condition pour le taux de presence.
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF ( IUSLAG(NC,NB,IJPRPD).LT.1 .OR.
     &         IUSLAG(NC,NB,IJPRPD).GT.2      ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1053) NB, NC, IUSLAG(NC,NB,IJPRPD)
          ENDIF
        ENDDO
      ENDDO
C
C --> Type de condition pour la vitesse.
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF ( IUSLAG(NC,NB,IJUVW).LT.-1 .OR.
     &         IUSLAG(NC,NB,IJUVW).GT.2      ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1050) NB,NC,IUSLAG(NC,NB,IJUVW)
          ENDIF
        ENDDO
      ENDDO
C
C --> Type de condition pour le diametre.
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF ( IUSLAG(NC,NB,IJPRDP).LT.1 .OR.
     &         IUSLAG(NC,NB,IJPRDP).GT.2      ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1051) NB, NC, IUSLAG(NC,NB,IJPRDP)
          ENDIF
        ENDDO
      ENDDO
C
C --> Type de condition pour le diametre.
C
      IF ( IPHYLA.EQ.1 .AND.
     &       (ITPVAR.EQ.1 .OR. IDPVAR.EQ.1 .OR. IMPVAR.EQ.1) ) THEN
        DO II = 1, NFRLAG
          NB = ILFLAG(II)
          DO NC = 1, IUSNCL(NB)
            IF ( IUSLAG(NC,NB,IJPRTP).LT.1 .OR.
     &           IUSLAG(NC,NB,IJPRTP).GT.2      ) THEN
              IOK = IOK + 1
              WRITE(NFECRA,1052) NB, NC, IUSLAG(NC,NB,IJPRTP)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C --> Poids statistiques
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (RUSLAG(NC,NB,IPOIT).LE.0.D0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1055) NB, NC, RUSLAG(NC,NB,IPOIT)
          ENDIF
        ENDDO
      ENDDO
C
C --> Debit massique de particule
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (RUSLAG(NC,NB,IDEBT).LT.0.D0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1056) NB, NC, RUSLAG(NC,NB,IDEBT)
          ENDIF
        ENDDO
      ENDDO
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (RUSLAG(NC,NB,IDEBT).GT.0.D0 .AND.
     &        IUSLAG(NC,NB,IJNBP).EQ.0         ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1057) NB, NC, RUSLAG(NC,NB,IDEBT),
     &                                 IUSLAG(NC,NB,IJNBP)
          ENDIF
        ENDDO
      ENDDO
C
C --> Proprietes des particules : le diametre, son ecart-type, et rho
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (RUSLAG(NC,NB,IROPT).LT.0.D0 .OR.
     &        RUSLAG(NC,NB,IDPT) .LT.0.D0 .OR.
     &        RUSLAG(NC,NB,IVDPT).LT.0.D0       ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1060) NB, NC,
     &                         RUSLAG(NC,NB,IROPT),
     &                         RUSLAG(NC,NB,IDPT),
     &                         RUSLAG(NC,NB,IVDPT)
          ENDIF
        ENDDO
      ENDDO
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF ( RUSLAG(NC,NB,IDPT).LT.3.D0*RUSLAG(NC,NB,IVDPT) ) THEN
            IOK = IOK + 1
            WRITE(NFECRA,1065) NB, NC,
     &        RUSLAG(NC,NB,IDPT)-3.D0*RUSLAG(NC,NB,IVDPT)
          ENDIF
        ENDDO
      ENDDO
C
C --> Proprietes des particules : Temperature et CP
C
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
C
        DO II = 1, NFRLAG
          NB = ILFLAG(II)
          DO NC = 1, IUSNCL(NB)
            IF (RUSLAG(NC,NB,ICPT)  .LT.0.D0   .OR.
     &          RUSLAG(NC,NB,ITPT)  .LT.TKELVN       ) THEN
              IOK = IOK + 1
              WRITE(NFECRA,1070)
     &        IPHYLA, ITPVAR, NB, NC,
     &        RUSLAG(NC,NB,ITPT), RUSLAG(NC,NB,ICPT)
            ENDIF
          ENDDO
        ENDDO
C
      ENDIF
C
C --> Proprietes des particules : Emissivite
C
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1 .AND. IIRAYO.EQ.1) THEN
C
        DO II = 1, NFRLAG
          NB = ILFLAG(II)
          DO NC = 1, IUSNCL(NB)
            IF (RUSLAG(NC,NB,IEPSI) .LT.0.D0 .OR.
     &          RUSLAG(NC,NB,IEPSI) .GT.1.D0      ) THEN
              IOK = IOK + 1
              WRITE(NFECRA,1075)
     &        IPHYLA, ITPVAR, NB, NC, RUSLAG(NC,NB,IEPSI)
            ENDIF
          ENDDO
        ENDDO
C
      ENDIF
C
C Charbon
C
      IF (IPHYLA.EQ.2) THEN
C
C --> Numero du charbon
C
        DO II = 1, NFRLAG
          NB = ILFLAG(II)
          DO NC = 1, IUSNCL(NB)
            IF (IUSLAG(NC,NB,INUCHL).LT.0.D0
     &          .AND.  IUSLAG(NC,NB,INUCHL).GT.NCHARB) THEN
              IOK = IOK + 1
              WRITE(NFECRA,1080) NB, NC, NCHARB, IUSLAG(NC,NB,INUCHL)
            ENDIF
          ENDDO
        ENDDO
C
C --> Proprietes des particules de Charbon.
C
        DO II = 1, NFRLAG
          NB = ILFLAG(II)
          DO NC = 1, IUSNCL(NB)
            IF (RUSLAG(NC,NB,IHPT)  .LT.0.D0 .OR.
     &          RUSLAG(NC,NB,ICPT)  .LT.0.D0 .OR.
     &          RUSLAG(NC,NB,IMCHT) .LT.0.D0 .OR.
     &          RUSLAG(NC,NB,IMCKT) .LT.0.D0  ) THEN
              IOK = IOK + 1
              WRITE(NFECRA,1090)
     &        IPHYLA, NB, NC,
     &        RUSLAG(NC,NB,IHPT),  RUSLAG(NC,NB,ICPT),
     &        RUSLAG(NC,NB,IMCHT), RUSLAG(NC,NB,IMCKT)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C --> Stop si erreur.
C
      IF(IOK.GT.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C
C=======================================================================
C 4. Transformation des donnees utilisateur
C=======================================================================
C
C --> Injection des part 1ere iter seulement si freq d'injection nulle
C
      DO II = 1, NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (IUSLAG(NC,NB,IJFRE).EQ.0 .AND. IPLAS.EQ.1) THEN
            IUSLAG(NC,NB,IJFRE) = NTCABS
          ENDIF
          IF (IUSLAG(NC,NB,IJFRE).EQ.0 .AND. IPLAS.GT.1) THEN
            IUSLAG(NC,NB,IJFRE) = NTCABS+1
          ENDIF
        ENDDO
      ENDDO
C
C --> Calcul du nombre de particules a injecter pour cette iteration
C
      NBPNEW = 0
      DNBPNW = 0.D0
C
C     Dans le cas ou on a un taux de presence impose dans une zone,
C     on corrige IUSLAG(NC,NB,IJNBP) donne dans USLAG2 qui n'a
C     pas de sens puisque l'on injecte 1 particule par maille
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
C       pour chaque classe :
        DO NC = 1, IUSNCL(NB)
C         si de nouvelles particules doivent entrer :
          IF (MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0) THEN
C
            IF ( IUSLAG(NC,NB,IJPRPD) .EQ. 2 ) THEN
              IUSLAG(NC,NB,IJNBP)=0
              DO IFAC = 1,NFABOR
                IF (IFRLAG(IFAC).EQ.NB) THEN
                  IUSLAG(NC,NB,IJNBP)=IUSLAG(NC,NB,IJNBP)+1
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        DO NC = 1, IUSNCL(NB)
          IF (MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0) THEN
            NBPNEW = NBPNEW + IUSLAG(NC,NB,IJNBP)
          ENDIF
        ENDDO
      ENDDO
C
C --> Limite du nombre de particules a NBPMAX
C
      IF ( (NBPART+NBPNEW).GT.NBPMAX ) THEN
        WRITE(NFECRA,3000) NBPART,NBPNEW,NBPMAX
        NBPNEW = 0
      ENDIF
C
C --> Si pas de new particules alors RETURN
C
      IF (NBPNEW.EQ.0) RETURN
C
C --> Tirage aleatoire des positions des NBPNEW nouvelles particules
C   au niveau des zones de bord et reperage des cellules correspondantes
C
C   initialisation du compteur de nouvelles particules
C
      NPT = NBPART
C
C     On reserve d'abord la memoire si on a de nouvelles particules
C       on garde IIWORK jusqu'a lagnwc
C
      IIWORK = IDEBIA
      IFINIA = IIWORK + NBPNEW
      IFINRA = IDEBRA
      CALL IASIZE('LAGENT',IFINIA)
C     ===========
C
C     Ensuite, on regarde ou on les met
C
C     pour chaque zone de bord :
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
C       pour chaque classe :
        DO NC = 1, IUSNCL(NB)
C         si de nouvelles particules doivent entrer :
          IF (MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0) THEN
C
            IF ( IUSLAG(NC,NB,IJPRPD) .EQ. 1 ) THEN
C
            CALL LAGNEW
C           ===========
     &  ( IFINIA , IFINRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NPT    , NBPNEW , IUSLAG(NC,NB,IJNBP)      ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    NB     ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    IFRLAG , ITEPA(1,JISOR)  , IA(IIWORK)      ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    SURFBN , ETTP   ,
     &    RDEVEL , RTUSER , RA     )
C
            ELSEIF ( IUSLAG(NC,NB,IJPRPD) .EQ. 2 ) THEN
C
              CALL LAGNPR
C             ===========
     &  ( IFINIA , IFINRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NPT    , NBPNEW , IUSLAG(NC,NB,IJNBP)      ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    NB     ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    IFRLAG , ITEPA(1,JISOR)  , IA(IIWORK)      ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    SURFBN , ETTP   ,
     &    RDEVEL , RTUSER , RA     )
            ENDIF
C
          ENDIF
        ENDDO
      ENDDO
C
C-->TEST DE CONTROLE (NE PAS MODIFIER)
C
      IF ( (NBPART+NBPNEW).NE.NPT ) THEN
        WRITE(NFECRA,3010) NBPNEW, NPT-NBPART
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C--> Injection en continu des particules
C
      IF ( INJCON.EQ.1 ) THEN
C
C   reinitialisation du compteur de nouvelles particules
C
        NPT = NBPART
C
C       pour chaque zone de bord:
C
        DO II = 1,NFRLAG
          NB = ILFLAG(II)
C
C         pour chaque classe :
C
          DO NC = 1, IUSNCL(NB)
C
C           si de nouvelles particules doivent entrer :
            IF ( MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0 ) THEN
C
              CALL LAGNWC
C             ===========
     &  ( IFINIA , IFINRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NPT    , NBPNEW , IUSLAG(NC,NB,IJNBP)      ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    ITYCEL , ICOCEL ,
     &    IFRLAG , ITEPA(1,JISOR)  , IA(IIWORK) ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    SURFBN , ETTP   ,
     &    RDEVEL , RTUSER , RA     )
C
            ENDIF
C
          ENDDO
        ENDDO
C
      ENDIF
C
C-->TEST DE CONTROLE (NE PAS MODIFIER)
C
      IF ( (NBPART+NBPNEW).NE.NPT ) THEN
        WRITE(NFECRA,3010) NBPNEW, NPT-NBPART
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C
C   reinitialisation du compteur de nouvelles particules
C
      NPT = NBPART
C
C     pour chaque zone de bord:
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
C
C       pour chaque classe :
        DO NC = 1, IUSNCL(NB)
C
C         si de nouvelles particules doivent entrer :
          IF (MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0) THEN
C
            DO IP = NPT+1 , NPT+IUSLAG(NC,NB,IJNBP)
              IEL = ITEPA(IP,JISOR)
              IFAC = IA(IIWORK+IP-1)
C
C-->COMPOSANTES DE LA VITESSE DES PARTICULES
C
C             si composantes de la vitesse imposee :
              IF (IUSLAG(NC,NB,IJUVW).EQ.1) THEN
                ETTP(IP,JUP) = RUSLAG(NC,NB,IUPT)
                ETTP(IP,JVP) = RUSLAG(NC,NB,IVPT)
                ETTP(IP,JWP) = RUSLAG(NC,NB,IWPT)
C
C             si norme de la vitesse imposee :
              ELSE IF (IUSLAG(NC,NB,IJUVW).EQ.0) THEN
                AA = -1.D0 / SURFBN(IFAC)
                VN1 = SURFBO(1,IFAC) * AA
                VN2 = SURFBO(2,IFAC) * AA
                VN3 = SURFBO(3,IFAC) * AA
                ETTP(IP,JUP) = VN1 * RUSLAG(NC,NB,IUNO)
                ETTP(IP,JVP) = VN2 * RUSLAG(NC,NB,IUNO)
                ETTP(IP,JWP) = VN3 * RUSLAG(NC,NB,IUNO)
C
C             si vitesse du fluide vu :
              ELSE IF (IUSLAG(NC,NB,IJUVW).EQ.-1) THEN
                ETTP(IP,JUP) = RTPA(IEL,IU(IPHAS))
                ETTP(IP,JVP) = RTPA(IEL,IV(IPHAS))
                ETTP(IP,JWP) = RTPA(IEL,IW(IPHAS))
C
C             si profil de vitesse impose :
              ELSE IF (IUSLAG(NC,NB,IJUVW).EQ.2) THEN
C
               IDVAR = 1
               XXPART = ETTP(IP,JXP)
               YYPART = ETTP(IP,JYP)
               ZZPART = ETTP(IP,JZP)
C
               CALL USLAPR
C              ===========
     & ( IDEBIA , IDEBRA ,
     &   IDVAR  , IEL    , NB     , NC     ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG ,
     &   IDEVEL , ITUSER , IA     ,
     &   XXPART , YYPART , ZZPART ,
     &   TVPART , UUPART , VVPART , WWPART , DDPART , TTPART ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C
                ETTP(IP,JUP) = UUPART
                ETTP(IP,JVP) = VVPART
                ETTP(IP,JWP) = WWPART
C
              ENDIF
C
C-->Vitesse du fluide vu
C
              ETTP(IP,JUF) = RTPA(IEL,IU(IPHAS))
              ETTP(IP,JVF) = RTPA(IEL,IV(IPHAS))
              ETTP(IP,JWF) = RTPA(IEL,IW(IPHAS))
C
C--> TEMPS DE SEJOUR
C
              TEPA(IP,JRTSP) = 0.D0
C
C--> Diametre
C
C             si diametre constant imposee :
              IF (IUSLAG(NC,NB,IJPRDP).EQ.1) THEN
                IF (RUSLAG(NC,NB,IVDPT) .GT. 0.D0) THEN
                  N1 = 1
                  CALL NORMALEN(N1,RD)
                  ETTP(IP,JDP) = RUSLAG(NC,NB,IDPT)
     &                         + RD(1) * RUSLAG(NC,NB,IVDPT)
C
C    On verifie qu'on obtient un diamètre dans la gamme des 99,7%
C
                  D3 = 3.D0 * RUSLAG(NC,NB,IVDPT)
                  IF (ETTP(IP,JDP).LT.RUSLAG(NC,NB,IDPT)-D3)
     &              ETTP(IP,JDP)= RUSLAG(NC,NB,IDPT)
                  IF (ETTP(IP,JDP).GT.RUSLAG(NC,NB,IDPT)+D3)
     &              ETTP(IP,JDP)= RUSLAG(NC,NB,IDPT)
                ELSE
                  ETTP(IP,JDP) = RUSLAG(NC,NB,IDPT)
                ENDIF
C
C             si profil pour le diametre  :
              ELSE IF (IUSLAG(NC,NB,IJPRDP).EQ.2) THEN
C
                IDVAR = 2
                XXPART = ETTP(IP,JXP)
                YYPART = ETTP(IP,JYP)
                ZZPART = ETTP(IP,JZP)
C
                CALL USLAPR
C               ===========
     & ( IDEBIA , IDEBRA ,
     &   IDVAR  , IEL    , NB     , NC     ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG ,
     &   IDEVEL , ITUSER , IA     ,
     &   XXPART , YYPART , ZZPART ,
     &   TVPART , UUPART , VVPART , WWPART , DDPART , TTPART ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C
                ETTP(IP,JDP) = DDPART
C
              ENDIF
C
C--> Autres variables : masse, ... en fonction de la physique
C
              D3 = ETTP(IP,JDP) * ETTP(IP,JDP) * ETTP(IP,JDP)
C
              IF (NBCLST.GT.0) THEN
                ITEPA(IP,JCLST) = IUSLAG(NC,NB,ICLST)
              ENDIF

              IF ( IPHYLA.EQ.0 .OR. IPHYLA.EQ.1 ) THEN
C
                ETTP(IP,JMP) = RUSLAG(NC,NB,IROPT) * PIS6 * D3

                IF ( IPHYLA.EQ.1 .AND. ITPVAR.EQ.1 ) THEN
C
C             si Temperature constante imposee :
                  IF (IUSLAG(NC,NB,IJPRTP).EQ.1) THEN
                    ETTP(IP,JTP) = RUSLAG(NC,NB,ITPT)
C             si profil pour la temperature :
                  ELSE IF (IUSLAG(NC,NB,IJPRTP).EQ.2) THEN
C
                    IDVAR = 3
                    XXPART = ETTP(IP,JXP)
                    YYPART = ETTP(IP,JYP)
                    ZZPART = ETTP(IP,JZP)
C
                    CALL USLAPR
C                   ===========
     & ( IDEBIA , IDEBRA ,
     &   IDVAR  , IEL    , NB     , NC     ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG ,
     &   IDEVEL , ITUSER , IA     ,
     &   XXPART , YYPART , ZZPART ,
     &   TVPART , UUPART , VVPART , WWPART , DDPART , TTPART ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C
                    ETTP(IP,JTP) = TTPART
C
                  ENDIF
C
                  IF ( IPPMOD(ICP3PL).GE.0 .OR.
     &                 IPPMOD(ICPL3C).GE.0 .OR.
     &                 IPPMOD(ICFUEL).GE.0      ) THEN
C
                    ETTP(IP,JTF) = PROPCE(IEL,IPPROC(ITEMP1)) -TKELVI
C
                  ELSE IF ( IPPMOD(ICOD3P).GE.0 .OR.
     &                      IPPMOD(ICOEBU).GE.0 .OR.
     &                      IPPMOD(IELARC).GE.0 .OR.
     &                      IPPMOD(IELJOU).GE.0      ) THEN
C
                    ETTP(IP,JTF) = PROPCE(IEL,IPPROC(ITEMP)) -TKELVI
C
C Kelvin
                  ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.1 ) THEN
C
                    ETTP(IP,JTF) = RTPA(IEL,ISCA(ISCALT(IPHAS))) -TKELVI
C
C Celsius
                  ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.-1 ) THEN
C
                    ETTP(IP,JTF) = RTPA(IEL,ISCA(ISCALT(IPHAS)))
C
                  ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.2 ) THEN
C
                    MODE = 1
                    CALL USTHHT(MODE, RTPA(IEL,ISCA(ISCALT(IPHAS))),
     &                          ETTP(IP,JTF))
C
                  ENDIF
C
                  ETTP(IP,JCP) = RUSLAG(NC,NB,ICPT)
                  TEPA(IP,JREPS) = RUSLAG(NC,NB,IEPSI)
C
                ENDIF
C
              ELSE IF ( IPHYLA.EQ.2 ) THEN
C
                ITEPA(IP,JINCH)  = IUSLAG(NC,NB,INUCHL)
                ETTP(IP,JHP) = RUSLAG(NC,NB,IHPT)
                ETTP(IP,JTF) = PROPCE(IEL,IPPROC(ITEMP1)) - TKELVI
                ETTP(IP,JCP) = RUSLAG(NC,NB,ICPT)
C
                ETTP(IP,JMCH) = RUSLAG(NC,NB,IMCHT)
                ETTP(IP,JMCK) = RUSLAG(NC,NB,IMCKT)
C
                TEPA(IP,JRDCK) = ETTP(IP,JDP)
                TEPA(IP,JRD0P) = ETTP(IP,JDP)
C
                ICHA = ITEPA(IP,JINCH)
                ETTP(IP,JMP) = ETTP(IP,JMCH)
     &                       + ETTP(IP,JMCK)
     &             + XASHCH(ICHA) * PIS6 * D3 * RHO0CH(ICHA)
C
              ENDIF
C
C--> POIDS STATISTIQUE
C
              IF (IUSLAG(NC,NB,IJPRPD).EQ.1) THEN
                TEPA(IP,JRPOI) = RUSLAG(NC,NB,IPOIT)
              ELSE IF (IUSLAG(NC,NB,IJPRPD).EQ.2) THEN
C
                IDVAR = 0
                XXPART = ETTP(IP,JXP)
                YYPART = ETTP(IP,JYP)
                ZZPART = ETTP(IP,JZP)
C
                CALL USLAPR
C               ===========
     & ( IDEBIA , IDEBRA ,
     &   IDVAR  , IEL    , NB     , NC     ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG ,
     &   IDEVEL , ITUSER , IA     ,
     &   XXPART , YYPART , ZZPART ,
     &   TVPART , UUPART , VVPART , WWPART , DDPART , TTPART ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C
                VOLP = PIS6*D3
                SURF = SQRT( SURFBO(1,IFAC)*SURFBO(1,IFAC)
     &                      +SURFBO(2,IFAC)*SURFBO(2,IFAC)
     &                      +SURFBO(3,IFAC)*SURFBO(3,IFAC) )
                VITP = SQRT( ETTP(IP,JUP)*ETTP(IP,JUP)
     &                      +ETTP(IP,JVP)*ETTP(IP,JVP)
     &                      +ETTP(IP,JWP)*ETTP(IP,JWP) )
                TEPA(IP,JRPOI) =TVPART*(SURF*VITP*DTP)/VOLP
C
               ENDIF
C
            ENDDO
C
            NPT = NPT + IUSLAG(NC,NB,IJNBP)
C
          ENDIF
C
        ENDDO
      ENDDO
C
C-->TEST DE CONTROLE (NE PAS MODIFIER)
C
      IF ( (NBPART+NBPNEW).NE.NPT ) THEN
        WRITE(NFECRA,3010) NBPNEW, NPT-NBPART
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C
C=======================================================================
C 5. MODIFICATION DES POIDS POUR AVOIR LE DEBIT
C=======================================================================
C
C   reinitialisation du compteur de nouvelles particules
C
        NPT = NBPART
C
C     pour chaque zone de bord :
C
        DO II = 1,NFRLAG
          NB = ILFLAG(II)
C
C         pour chaque classe :
C
        DO NC = 1,IUSNCL(NB)
C
C         si de nouvelles particules sont entrees,
C         et si on a un debit non nul :
C
          IF ( MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0 .AND.
     &         RUSLAG(NC,NB,IDEBT) .GT. 0.D0        .AND.
     &         IUSLAG(NC,NB,IJNBP) .GT. 0                 ) THEN
C
            DMASSE = 0.D0
            DO IP = NPT+1 , NPT+IUSLAG(NC,NB,IJNBP)
              DMASSE = DMASSE + ETTP(IP,JMP)
            ENDDO
C
C        Calcul des Poids
C
            IF ( DMASSE.GT.0.D0 ) THEN
              DO IP = NPT+1 , NPT+IUSLAG(NC,NB,IJNBP)
                TEPA(IP,JRPOI) = ( RUSLAG(NC,NB,IDEBT)*DTP ) / DMASSE
              ENDDO
            ELSE
              WRITE(NFECRA,1057) NB, NC, RUSLAG(NC,NB,IDEBT),
     &                                   IUSLAG(NC,NB,IJNBP)
              CALL CSEXIT (1)
C             ===========
            ENDIF
C
            ENDIF
C
            NPT = NPT + IUSLAG(NC,NB,IJNBP)
C
          ENDDO
C
        ENDDO
C
C-->TEST DE CONTROLE (NE PAS MODIFIER)
C
      IF ( (NBPART+NBPNEW).NE.NPT ) THEN
        WRITE(NFECRA,3010) NBPNEW, NPT-NBPART
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C 6. CALCUL DE LA MASSE TOTALE INJECTES EN CHAQUE ZONE
C    Attention cette valeur est modifie dans USLABO pour tenir compte
C    des particules qui sortent
C    + calcul du nombres physiques de particules qui rentrent (tenant
C       compte des poids)
C=======================================================================
C
C   reinitialisation du compteur de nouvelles particules
C
      NPT     = NBPART
      DNBPNW = 0.D0
C
C     pour chaque zone de bord :
C
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
        DEBLAG(NB) = 0.D0
C
C       pour chaque classe :
C
        DO NC = 1,IUSNCL(NB)
C
C        si de nouvelles particules sont entrees,
C
          IF ( MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0 .AND.
     &         IUSLAG(NC,NB,IJNBP) .GT. 0                 ) THEN
C
            DO IP = NPT+1 , NPT+IUSLAG(NC,NB,IJNBP)
              DEBLAG(NB) = DEBLAG(NB) + TEPA(IP,JRPOI)*ETTP(IP,JMP)
              DNBPNW = DNBPNW + TEPA(IP,JRPOI)
            ENDDO
C
          ENDIF
C
          NPT = NPT + IUSLAG(NC,NB,IJNBP)
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C 7. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANNEES VUES
C    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
C=======================================================================
C
C   si de nouvelles particules doivent entrer :
C
      NPAR1 = NBPART+1
      NPAR2 = NBPART+NBPNEW
C
      CALL LAGIPN
C     ===========
     &  ( IFINIA , IFINRA ,
     &    NCELET , NCEL   ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NPAR1  , NPAR2  ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    ITEPA  ,
     &    IDEVEL , ITUSER , IA     ,
     &    RTPA   ,
     &    ETTP   , TEPA   , VAGAUS ,
     &    W1     , W2     , W3     ,
     &    RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 8. MODIFICATION DES TABLEAUX DE DONNEES PARTICULAIRES
C=======================================================================
C
      CALL USLAIN
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NBPNEW ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG , IA(IIWORK) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , VAGAUS , W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C   reinitialisation du compteur de nouvelles particules
      NPT = NBPART
C
C     pour chaque zone de bord:
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
C
C       pour chaque classe :
        DO NC = 1, IUSNCL(NB)
C
C         si de nouvelles particules doivent entrer :
          IF (MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0) THEN
C
            DO IP = NPT+1 , NPT+IUSLAG(NC,NB,IJNBP)
C
              IF (ETTP(IP,JDP).LT.0.D0 .AND.
     &            RUSLAG(NC,NB,IVDPT).GT.0.D0) THEN
                WRITE(NFECRA,4000) RUSLAG(NC,NB,IDPT),
     &                             RUSLAG(NC,NB,IVDPT),
     &                             ETTP(IP,JDP)
              ENDIF
C
            ENDDO
C
            NPT = NPT + IUSLAG(NC,NB,IJNBP)
C
          ENDIF
C
        ENDDO
      ENDDO
C
C     Ici on peut laisser choir IIWORK (et donc reprendre
C       IDEBIA et IDEBRA comme indicateur de la zone de memoire libre)
C
C=======================================================================
C 8. IMPRESSIONS POUR POST-PROCESSING EN MODE TRAJECTOIRES
C=======================================================================
C
      IF ( IENSI1.EQ.1 ) THEN
C
        NFIN = 0
C
C   reinitialisation du compteur de nouvelles particules
        NPT = NBPART
C
C       pour chaque zone de bord :
        DO II = 1, NFRLAG
          NB = ILFLAG(II)
C
C         pour chaque classe :
          DO NC = 1, IUSNCL(NB)
C
C           si de nouvelles particules doivent entrer :
            IF (MOD(NTCABS,IUSLAG(NC,NB,IJFRE)).EQ.0) THEN
C
              DO IP = NPT+1 , NPT+IUSLAG(NC,NB,IJNBP)
                CALL ENSLAG
C               ===========
     &           ( IFINIA , IFINRA ,
     &             NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &             NFIN   , IP     ,
     &             ITEPA  ,
     &             ETTP   , TEPA   , RA)
              ENDDO
              NPT = NPT + IUSLAG(NC,NB,IJNBP)
C
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
C 9. NOUVEAU NOMBRE DE PARTICULES TOTAL
C=======================================================================
C
C     NBPART : NOMBRE DE PARTICULES PRESENTES DANS LE DOMAINE
C
C     NBPTOT : NOMBRE DE PARTICULES TOTAL INJECTE DANS
C              LE CALCUL DEPUIS LE DEBUT SUITE COMPRISE
C
      NBPART = NBPART + NBPNEW
      DNBPAR = DNBPAR + DNBPNW
C
      NBPTOT = NBPTOT + NBPNEW
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,
     &'@                                                            ',/,
     &'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,
     &'@    un entier strictement positif et inferieur ou egal a    ',/,
     &'@    NFLAGM = ',I10                                           ,/,
     &'@  Ce numero (IFRLAG(IFAC)) vaut ici ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,
     &'@                                                            ',/,
     &'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,
     &'@    definies par l''utilisateur est NFLAGM = ',I10           ,/,
     &'@    Il a ete depasse.                                       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@  Les zones frontieres NFLAGM premieres zones frontieres    ',/,
     &'@    portent ici les numeros suivants :                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,
     &'@                                                            ',/,
     &'@  Le nombre de classes de la zone numero ',I10   ,' doit    ',/,
     &'@    etre un entier positif ou nul et inferieur ou egal      ',/,
     &'@    a NCLAGM = ',I10                                         ,/,
     &'@  Ce nombre (IUSNCL(NB)  ) vaut ici ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le nombre de particules dans la classe         ',I10       ,/,
     &'@                          dans la zone frontiere ',I10       ,/,
     &'@  doit etre un entier strictement positif et                ',/,
     &'@                      inferieur ou egal a NBPMAX = ',I10     ,/,
     &'@                                                            ',/,
     &'@  Ce nombre (IUSLAG(NC,NB,IJNBP)) vaut ici ',I10             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME A L''EXECUTION DU MODULE LAGRANGIEN',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le numero de groupe statistique de particules             ',/,
     &'@     dans la classe         NC= ',I10                        ,/,
     &'@     dans la zone frontiere NB= ',I10                        ,/,
     &'@  doit etre un entier strictement positif et                ',/,
     &'@                      inferieur ou egal a NBCLST = ',I10     ,/,
     &'@                                                            ',/,
     &'@  Ce nombre (IUSLAG(NC,NB,IJNBP)) vaut ici ',I10             ,/,
     &'@                                                            ',/,
     &'@  Le calcul continue mais cette classe statistique sera     ',/,
     &'@  ignoree                                                   ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  La frequence d''injection des particules doit etre un     ',/,
     &'@    entier positif ou nul (nul signifiant que               ',/,
     &'@    les particules ne injectees qu''au debut du calcul).    ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere NB  =',I10                     ,/,
     &'@    et      pour la classe    NC  =',I10                     ,/,
     &'@    vaut ici IUSLAG (NC,NB,IJFRE) =',I10                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Les conditions aux bords sont representees par            ',/,
     &'@    une variable dont les valeurs sont                      ',/,
     &'@    obligatoirement les suivantes :                         ',/,
     &'@                                                            ',/,
     &'@   = IENTRL  zone d''injection de particules                ',/,
     &'@   = ISORTL  sortie du domaine                              ',/,
     &'@   = IREBOL  rebond des particules                          ',/,
     &'@   = IDEPO1  deposition definitive                          ',/,
     &'@   = IDEPO2  deposition definitive mais la particule reste  ',/,
     &'@             en memoire (utile si IENSI2 = 1 uniquement)    ',/,
     &'@   = IDEPO3  deposition et remise en suspension possible    ',/,
     &'@             suivant les condition de l''ecoulement         ',/,
     &'@   = IENCRL  encrassement (Charbon uniquement IPHYLA = 2)   ',/,
     &'@   = JBORD1  interaction particule/frontiere utilisateur    ',/,
     &'@   = JBORD2  interaction particule/frontiere utilisateur    ',/,
     &'@   = JBORD3  interaction particule/frontiere utilisateur    ',/,
     &'@   = JBORD4  interaction particule/frontiere utilisateur    ',/,
     &'@   = JBORD5  interaction particule/frontiere utilisateur    ',/,
     &'@                                                            ',/,
     &'@  Cette valeur pour la frontiere NB = ',I10                  ,/,
     &'@     est erronees.                                          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  La condition a la limite representee par la variable      ',/,
     &'@    IDEPO2 n''est admissible que pour un post-processing    ',/,
     &'@    en mode deplacement IENSI2 = 1 .                        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENSI2 dans USLAG1 et les           ',/,
     &'@  conditions aux limites pour la frontiere NB =',I10         ,/,
     &'@  dans USLAG2.                                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1042 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  La condition a la limite representee par le type          ',/,
     &'@    IENCRL n''est admissible que lorsque les particules     ',/,
     &'@    transportees sont des grains de charbon IPHYLA = 2      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IPHYLA dans USLAG1 et les           ',/,
     &'@  conditions aux limites pour la frontiere NB =',I10         ,/,
     &'@  dans USLAG2.                                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1043 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  La condition a la limite representee par le type          ',/,
     &'@    ENCRAS n''est admissible que lorsque l''option          ',/,
     &'@    encrassement est enclenche IENCRA = 1                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENCRA dans USLAG1 et les           ',/,
     &'@  conditions aux limites pour la frontiere NB =',I10         ,/,
     &'@  dans USLAG2.                                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le type de condition au bord pour la vitesse est          ',/,
     &'@    represente par un entier dont les valeurs sont          ',/,
     &'@    obligatoirement les suivantes                           ',/,
     &'@   =-1 vitesse fluide imposee                               ',/,
     &'@   = 0 vitesse imposee selon la direction normale a la      ',/,
     &'@       face de bord et de norme RUSLAG(NB,NC,IUNO)          ',/,
     &'@   = 1 vitesse imposee : on donne RUSLAG(NB,NC,IUPT)        ',/,
     &'@                                  RUSLAG(NB,NC,IVPT)        ',/,
     &'@                                  RUSLAG(NB,NC,IWPT)        ',/,
     &'@   = 2 profil de vitesse imposee dans USLAPR                ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere NB = ',I10                     ,/,
     &'@    et      pour la classe    NC = ',I10                     ,/,
     &'@    vaut ici IUSLAG(NC,NB,IJUVW) = ',I10                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@  (Si IUSLAG(NC,NB,IJUVW) vaut -2 il n''a par ete renseigne)',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1051 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le type de condition au bord pour le diametre est         ',/,
     &'@    represente par un entier dont les valeurs sont          ',/,
     &'@    obligatoirement les suivantes                           ',/,
     &'@   = 1 diametre imposee dans la zone :                      ',/,
     &'@                  on donne        RUSLAG(NB,NC,IDPT)        ',/,
     &'@                                  RUSLAG(NB,NC,IVDPT)       ',/,
     &'@   = 2 profil de diametre imposee dans USLAPR               ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere  NB = ',I10                    ,/,
     &'@    et      pour la classe     NC = ',I10                    ,/,
     &'@    vaut ici IUSLAG(NC,NB,IJPRDP) = ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@ (Si IUSLAG(NC,NB,IJPRDP) vaut -2 il n''a par ete renseigne)',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1052 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le type de condition au bord pour la temperature est      ',/,
     &'@    represente par un entier dont les valeurs sont          ',/,
     &'@    obligatoirement les suivantes                           ',/,
     &'@   = 1 temperature imposee dans la zone :                   ',/,
     &'@                  on donne        RUSLAG(NB,NC,ITPT))       ',/,
     &'@   = 2 profil de temperature imposee dans USLAPR            ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere  NB = ',I10                    ,/,
     &'@    et      pour la classe     NC = ',I10                    ,/,
     &'@    vaut ici IUSLAG(NC,NB,IJPRTP) = ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@ (Si IUSLAG(NC,NB,IJPRTP) vaut -2 il n''a par ete renseigne)',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1053 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le type de condition au bord pour la taux de presence est ',/,
     &'@    represente par un entier dont les valeurs sont          ',/,
     &'@    obligatoirement les suivantes                           ',/,
     &'@   = 1 distribution uniforme                                ',/,
     &'@                  on donne        RUSLAG(NB,NC,IPOID))      ',/,
     &'@   = 2 profil de taux de presence imposee dans USLAPR       ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere  NB = ',I10                    ,/,
     &'@    et      pour la classe     NC = ',I10                    ,/,
     &'@    vaut ici IUSLAG(NC,NB,IJPRPD) = ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@ (Si IUSLAG(NC,NB,IJPRPD) vaut -2 il n''a par ete renseigne)',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1055 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le poids statistique des particules doit etre un reel     ',/,
     &'@    strictement positif.                                    ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere NB  =',I10                     ,/,
     &'@    et      pour la classe    NC  =',I10                     ,/,
     &'@    vaut ici RUSLAG (NC,NB,IPOIT) =',E14.5                   ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1056 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le debit massique des particules doit etre un reel        ',/,
     &'@    positif ou nul (nul signifiant que le debit n''est pas  ',/,
     &'@    pris en compte).                                        ',/,
     &'@                                                            ',/,
     &'@  Ce nombre pour la frontiere NB  =',I10                     ,/,
     &'@    et      pour la classe    NC  =',I10                     ,/,
     &'@    vaut ici RUSLAG (NC,NB,IDEBT) =',E14.5                   ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1057 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Les proprietes physiques des particules au bord pour      ',/,
     &'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,
     &'@    ne sont pas physiques.                                  ',/,
     &'@  Le debit massique impose                                  ',/,
     &'@    vaut ici RUSLAG (NC,NB,IDEBT) =',E14.5                   ,/,
     &'@  alors que le nombre de particules injectees est nul       ',/,
     &'@             IUSLAG (NC,NB,IJNBP) =',I10                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1060 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Les proprietes physiques des particules au bord pour      ',/,
     &'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,
     &'@    ne sont pas physiques :                                 ',/,
     &'@    Masse volumique : RUSLAG(NC,NB,IROPT) = ',E14.5          ,/,
     &'@    Diametre moyen  : RUSLAG(NC,NB,IDPT)  = ',E14.5          ,/,
     &'@    Ecart type      : RUSLAG(NC,NB,IVDPT) = ',E14.5          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@  Verifier le fichier dp_FCP si l''option Charbon pulverise ',/,
     &'@    est activee.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1065 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  L''ecart-type fourni est trop grand par rapport           ',/,
     &'@    au diametre moyen pour                                  ',/,
     &'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,
     &'@                                                            ',/,
     &'@  Il y a un risque non nul de calcul d''un diametre negatif.',/,
     &'@                                                            ',/,
     &'@  Theoriquement 99,7% des particules se trouvent entre :    ',/,
     &'@    RUSLAG(NC,NB,IDPT)-3*RUSLAG(NC,NB,IVDPT)                ',/,
     &'@    RUSLAG(NC,NB,IDPT)+3*RUSLAG(NC,NB,IVDPT)                ',/,
     &'@  Pour eviter des diametres aberrants, dans le module       ',/,
     &'@    lagrangien, avec un clipping, on impose que 100% des    ',/,
     &'@    particules doivent etre dans cet intervalle.            ',/,
     &'@                                                            ',/,
     &'@  Or on a :                                                 ',/,
     &'@    RUSLAG(NC,NB,IDPT)-3*RUSLAG(NC,NB,IVDPT) = ',E14.5       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@  Verifier le fichier dp_FCP si l''option Charbon pulverise ',/,
     &'@    est activee.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1070 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Une equation sur la temperature est associee              ',/,
     &'@    aux particules (IPHYLA = ',I10,') :                     ',/,
     &'@    ITPVAR = ',I10                                           ,/,
     &'@  Les proprietes physiques des particules au bord pour      ',/,
     &'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,
     &'@    doivent etre renseignees :                              ',/,
     &'@    Temperature  : RUSLAG(NC,NB,ITPT) = ',E14.5              ,/,
     &'@    Cp           : RUSLAG(NC,NB,ICPT) = ',E14.5              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1075 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Une equation sur la temperature est associee              ',/,
     &'@    aux particules (IPHYLA = ',I10,') :                     ',/,
     &'@    ITPVAR = ',I10                                           ,/,
     &'@    avec prise em compte des echanges thermiques radiatifs. ',/,
     &'@  L''emissivite des particules doit etre renseignee et      ',/,
     &'@    comprise entre 0 et 1 (inclus).                         ',/,
     &'@                                                            ',/,
     &'@  L''emissivite pour la frontiere NB = ',I10                 ,/,
     &'@                     et la classe NC = ',I10                 ,/,
     &'@    vaut : RUSLAG(NC,NB,IEPSI) = ',E14.5                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1080 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@    L''INDICATEUR SUR LE NUMERO DU CHARBON                  ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@  Le numero du charbon injecte pour                         ',/,
     &'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,
     &'@    devrait etre compris entre 1 et NCHARB= ',I10            ,/,
     &'@    Le nombre de charbon NCHARB est donne dans dp_FCP.      ',/,
     &'@                                                            ',/,
     &'@    Il vaut ici  : IUSLAG(NC,NB,INUCHL) = ',I10              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2 et le fichier dp_FCP.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1090 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  L''option de transport de particules de charbon pulverise ',/,
     &'@    est active (IPHYLA = ',I10,')                           ',/,
     &'@  Les proprietes physiques des particules au bord pour      ',/,
     &'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,
     &'@    doivent etre renseignees :                              ',/,
     &'@    Temperature   : RUSLAG(NC,NB,IHPT)  = ',E14.5            ,/,
     &'@    Cp            : RUSLAG(NC,NB,ICPT)  = ',E14.5            ,/,
     &'@    Masse de charbon                                        ',/,
     &'@    reactif       : RUSLAG(NC,NB,IMCHT) = ',E14.5            ,/,
     &'@    Masse de coke : RUSLAG(NC,NB,IMCKT) = ',E14.5            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2 et le fichier dp_FCP.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MODULE LAGRANGIEN                           ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le nombre de nouvelles particules injectees conduit a un  ',/,
     &'@    nombre total de particules superieur au maximum prevu : ',/,
     &'@                                                            ',/,
     &'@    Nombre de particules courant   : NBPART = ',I10          ,/,
     &'@    Nombre de nouvelles particules : NBPNEW = ',I10          ,/,
     &'@    Nombre maximal de particules   : NBPMAX = ',I10          ,/,
     &'@                                                            ',/,
     &'@  Le calcul se poursuit, mais on n''injecte aucune particule',/,
     &'@                                                            ',/,
     &'@  Ajuster NBPMAX dans USLAG1 et verifier USLAG2.            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le nombre de particules injectees dans le domaine         ',/,
     &'@    pour cette iteration Lagrangienne ne correspond pas     ',/,
     &'@    a celui specifie dans les conditions aux limites.       ',/,
     &'@                                                            ',/,
     &'@  Nombre de particules specifie pour l''injection :         ',/,
     &'@    NBPNEW = ',I10                                           ,/,
     &'@  Nombre de particules effectivement injectees :            ',/,
     &'@    NPT-NBPART = ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGENT)                                    ',/,
     &'@                                                            ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,
     &'@                                                            ',/,
     &'@  Le reconstruction du duametre d''une particule a partir   ',/,
     &'@    du diametre moyen et de l''ecart type donne une valeur  ',/,
     &'@    de diametre negatif, a cause d''un tirage aleatoire     ',/,
     &'@    dans un des "bord de la Gaussienne".                    ',/,
     &'@                                                            ',/,
     &'@  Diametre moyen : ',E14.5                                   ,/,
     &'@  Ecart type : ',E14.5                                       ,/,
     &'@  Diametre calcule : ',E14.5                                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Mettre en place une verification du diametre en fonction  ',/,
     &'@  des donnees granulometriques dans USLAIN.                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
