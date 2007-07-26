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
                        SUBROUTINE LAGSTA
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , TEPA   , STATIS , STATIV ,
     &   W1     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC        SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC        -----------------------------------
CFONC
CFONC     CALCUL DES STATISTIQUES SUR LES PARTICULES
CFONC
CFONC
CFONC   ISTTIO = 0 : calcul instationnaire pour le lagrangien
CFONC          = 1 : calcul stationnaire   pour le lagrangien
CFONC
CFONC   ISTALA : calcul statistiques       si  >= 1 sinon pas de stat
CFONC
CFONC   ISUIST : suite calcul statistiques si  >= 1 sinon pas de stat
CFONC
CFONC   IDSTNT : Numero du pas de temps pour debut statistque
CFONC
CFONC   NSTIST : iteration Lagrangienne du debut calcul stationnaire
CFONC
CFONC   NPST   : nombre d'iterations de calcul de stats stationnaires
CFONC
CFONC   NPSTT  : nombre d'iterations total des stats depuis le debut
CFONC              du calcul, partie instationnaire comprise
CFONC              (Attention : uniquement pour affichage listing,
CFONC               ne pas faire de test dessus, preferer IDSTNT)
CFONC
CFONC   TSTAT  : temps physique de calcul des statistiques stationnaires
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
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS(NCELET! TR ! <-  ! CUMUL DES STATISTIQUES VOLUMIQUES    !
CARGU !    NVLSTA)   !    !     !                                      !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! W1(NCELET)   ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
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
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION W1(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          NPT , NV , IEL, NV1, IZCL
      INTEGER          ILVX1 , ILVY1  , ILVZ1  , ILPD1  , ILFV1 , ILTS1
      INTEGER          ILTP1 , ILDP1  , ILMP1
      INTEGER          ILHP1 , ILMCH1 , ILMCK1 , ILDCK1
      DOUBLE PRECISION PIS6 , CONCEN
C
C***********************************************************************
C
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
C 1. INCREMENTATION DES COMPTEURS ET INITIALISATION
C=======================================================================
C
C      ISTTIO = 0 : calcul instationnaire pour le lagrangien
C             = 1 : calcul stationnaire   pour le lagrangien
C      NSTIST : iter de depart pour le debut calcul stationnaire
C      NPST   : nombre de pas de temps pour le cumul des stats
C               stationnaires
C      NPSTT  : nombre de pas de temps total des stats depuis le debut
C               du calcul, partie instationnaire comprise (pour listing)
C      TSTAT  : temps physique de calcul des statistiques stationnaires
C
C-> Si on est en instationnaire, ou si le debut des stat stationnaires
C   n'est pas encore enclenchee, on remet les stat a zero a chaque
C   pas de temps
C
      IF (ISTTIO.EQ.0 .OR. (ISTTIO.EQ.1 .AND. IPLAS.LE.NSTIST)) THEN
C
        NPST  = 0
        TSTAT = 0.D0
C
C Statistiques globales
        DO NV = 1,NVLSTA
          DO IEL = 1,NCEL
            STATIS(IEL,NV) = 0.D0
          ENDDO
        ENDDO
C
        DO NV = 1,NVLSTA-1
          DO IEL = 1,NCEL
            STATIV(IEL,NV) = 0.D0
          ENDDO
        ENDDO
C
C Statistiques par groupe
        IF (NBCLST.GT.0) THEN
C
          DO NV = 1, NVLSTA*NBCLST
            NV1 = NVLSTA + NV
            DO IEL = 1,NCEL
              STATIS(IEL,NV1) = 0.D0
            ENDDO
          ENDDO
C
          DO NV = 1, (NVLSTA-1)*NBCLST
            NV1 = (NVLSTA-1) + NV
            DO IEL = 1,NCEL
              STATIV(IEL,NV1) = 0.D0
            ENDDO
          ENDDO
C
        ENDIF
C
      ENDIF
C
      NPST  = NPST  + 1
      TSTAT = TSTAT + DTP
C
      NPSTT = IPLAS - IDSTNT + 1
C
C=======================================================================
C 2 - CALCUL DES STATISTIQUES PARTICULAIRES
C=======================================================================
C
C     * Moyenne et variance des composantes de la vitesse
C     * Moyenne et variance du taux de presence
C           (i.e. concentration volumique)
C     * Moyenne et variance du temps de séjour
C     * Somme du poids statistiques associé aux particules
C          (i.e. nombre de particules par cellules)
C
C     * Moyenne et variance de la temperature
C     * Moyenne et variance du diametre
C     * Moyenne et variance de la masse
C
C     * Moyenne et variance de la temperature
C     * Moyenne et variance de la masse de charbon reactif
C     * Moyenne et variance de la masse de coke
C     * Moyenne et variance du diametre du coeur retrecissant
C
C
      PIS6 = PI / 6.D0
C
C=======================================================================
C 2.1 -  STATISTIQUES GLOBALES
C=======================================================================
C
      DO NPT = 1,NBPART
C
        IF (ITEPA(NPT,JISOR).GT.0) THEN
C
C     * Moyenne et variance des composantes de la vitesse
C
          IEL = ITEPA(NPT,JISOR)
C
          STATIS(IEL,ILVX) = STATIS(IEL,ILVX)
     &                     + TEPA(NPT,JRPOI) * ETTP(NPT,JUP)
          STATIS(IEL,ILVY) = STATIS(IEL,ILVY)
     &                     + TEPA(NPT,JRPOI) * ETTP(NPT,JVP)
          STATIS(IEL,ILVZ) = STATIS(IEL,ILVZ)
     &                     + TEPA(NPT,JRPOI) * ETTP(NPT,JWP)
C
          STATIV(IEL,ILVX) = STATIV(IEL,ILVX)
     &     + TEPA(NPT,JRPOI) * ETTP(NPT,JUP) * ETTP(NPT,JUP)
C
          STATIV(IEL,ILVY) = STATIV(IEL,ILVY)
     &     + TEPA(NPT,JRPOI) * ETTP(NPT,JVP) * ETTP(NPT,JVP)
          STATIV(IEL,ILVZ) = STATIV(IEL,ILVZ)
     &     + TEPA(NPT,JRPOI) * ETTP(NPT,JWP) * ETTP(NPT,JWP)
C
C     * Moyenne et variance du taux de presence
C
          CONCEN = (ETTP(NPT,JDP)**3) *PIS6
C
          STATIS(IEL,ILFV) = STATIS(IEL,ILFV)
     &      + TEPA(NPT,JRPOI) * CONCEN
C
          STATIV(IEL,ILFV) = STATIV(IEL,ILFV)
     &      + TEPA(NPT,JRPOI) * CONCEN * CONCEN
C
C     * Moyenne et variance du temps de séjour
C
          STATIS(IEL,ILTS) = STATIS(IEL,ILTS)
     &      + TEPA(NPT,JRPOI) *  TEPA(NPT,JRTSP)
C
          STATIV(IEL,ILTS) = STATIV(IEL,ILTS)
     &      + TEPA(NPT,JRPOI) *  TEPA(NPT,JRTSP) *  TEPA(NPT,JRTSP)
C
C     * Somme du poids statistiques associé aux particules
C
          STATIS(IEL,ILPD) = STATIS(IEL,ILPD) + TEPA(NPT,JRPOI)
C
          IF (IPHYLA.EQ.1) THEN
C
C     * Moyenne et variance de la temperature
C
            IF ( ITPVAR .EQ. 1 ) THEN
C
              STATIS(IEL,ILTP) = STATIS(IEL,ILTP)
     &                         + TEPA(NPT,JRPOI) * ETTP(NPT,JTP)
C
              STATIV(IEL,ILTP) = STATIV(IEL,ILTP)
     &             + TEPA(NPT,JRPOI)*ETTP(NPT,JTP)*ETTP(NPT,JTP)
C
            ENDIF
C
C     * Moyenne et variance du diametre
C
            IF ( IDPVAR .EQ. 1 ) THEN
C
              STATIS(IEL,ILDP) = STATIS(IEL,ILDP)
     &                         + TEPA(NPT,JRPOI) * ETTP(NPT,JDP)
C
              STATIV(IEL,ILDP) = STATIV(IEL,ILDP)
     &             + TEPA(NPT,JRPOI)*ETTP(NPT,JDP)*ETTP(NPT,JDP)
C
            ENDIF
C
C     * Moyenne et variance de la masse
C
            IF ( IMPVAR .EQ. 1 ) THEN
C
              STATIS(IEL,ILMP) = STATIS(IEL,ILMP)
     &                         + TEPA(NPT,JRPOI) * ETTP(NPT,JMP)
C
              STATIV(IEL,ILMP) = STATIV(IEL,ILMP)
     &             + TEPA(NPT,JRPOI)*ETTP(NPT,JMP)*ETTP(NPT,JMP)
C
            ENDIF
C
          ELSE IF ( IPHYLA .EQ. 2 ) THEN
C
C     * Moyenne et variance de la temperature
C     * Moyenne et variance de la masse de charbon reactif
C     * Moyenne et variance de la masse de coke
C     * Moyenne et variance du diametre du coeur retrecissant
C
            STATIS(IEL,ILHP)  = STATIS(IEL,ILHP)
     &                        + TEPA(NPT,JRPOI) * ETTP(NPT,JHP)
            STATIS(IEL,ILMCH) = STATIS(IEL,ILMCH)
     &                        + TEPA(NPT,JRPOI) * ETTP(NPT,JMCH)
            STATIS(IEL,ILMCK) = STATIS(IEL,ILMCK)
     &                        + TEPA(NPT,JRPOI) * ETTP(NPT,JMCK)
            STATIS(IEL,ILDCK) = STATIS(IEL,ILDCK)
     &                        + TEPA(NPT,JRPOI) * TEPA(NPT,JRDCK)
C
            STATIV(IEL,ILHP)  = STATIV(IEL,ILHP)
     &          + TEPA(NPT,JRPOI)*ETTP(NPT,JHP)*ETTP(NPT,JHP)
            STATIV(IEL,ILMCH) = STATIV(IEL,ILMCH)
     &          + TEPA(NPT,JRPOI)*ETTP(NPT,JMCH)*ETTP(NPT,JMCH)
            STATIV(IEL,ILMCK) = STATIV(IEL,ILMCK)
     &          + TEPA(NPT,JRPOI)*ETTP(NPT,JMCK)*ETTP(NPT,JMCK)
            STATIV(IEL,ILDCK) = STATIV(IEL,ILDCK)
     &          + TEPA(NPT,JRPOI)*TEPA(NPT,JRDCK)*TEPA(NPT,JRDCK)
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 2.2 -  STATISTIQUES PAR GROUPE
C=======================================================================
C
      IF (NBCLST.GT.0) THEN
C
        DO IZCL = 1,NBCLST
C
          DO NPT = 1,NBPART
C
            IF (ITEPA(NPT,JISOR).GT.0 .AND.
     &                         ITEPA(NPT,JCLST).EQ. IZCL ) THEN
C
C     * Moyenne et variance des composantes de la vitesse
C
              IEL = ITEPA(NPT,JISOR)
C
              ILVX1 = ILVX + IZCL*NVLSTA
              ILVY1 = ILVY + IZCL*NVLSTA
              ILVZ1 = ILVZ + IZCL*NVLSTA
C
              STATIS(IEL,ILVX1) = STATIS(IEL,ILVX1)
     &                       + TEPA(NPT,JRPOI) * ETTP(NPT,JUP)
              STATIS(IEL,ILVY1) = STATIS(IEL,ILVY1)
     &                       + TEPA(NPT,JRPOI) * ETTP(NPT,JVP)
              STATIS(IEL,ILVZ1) = STATIS(IEL,ILVZ1)
     &                       + TEPA(NPT,JRPOI)*ETTP(NPT,JWP)
C
              ILVX1 = ILVX + IZCL*(NVLSTA-1)
              ILVY1 = ILVY + IZCL*(NVLSTA-1)
              ILVZ1 = ILVZ + IZCL*(NVLSTA-1)
C
              STATIV(IEL,ILVX1) = STATIV(IEL,ILVX1)
     &       + TEPA(NPT,JRPOI) * ETTP(NPT,JUP) * ETTP(NPT,JUP)
              STATIV(IEL,ILVY1) = STATIV(IEL,ILVY1)
     &       + TEPA(NPT,JRPOI) * ETTP(NPT,JVP) * ETTP(NPT,JVP)
              STATIV(IEL,ILVZ1) = STATIV(IEL,ILVZ1)
     &       + TEPA(NPT,JRPOI) * ETTP(NPT,JWP) * ETTP(NPT,JWP)
C
C     * Moyenne et variance du taux de presence
C
C
              CONCEN = (ETTP(NPT,JDP)**3) *PIS6
C
              ILFV1 = ILFV+IZCL*NVLSTA
              STATIS(IEL,ILFV1) = STATIS(IEL,ILFV1)
     &                          +TEPA(NPT,JRPOI)*CONCEN
C
              ILFV1 = ILFV+IZCL*(NVLSTA-1)
              STATIV(IEL,ILFV1) = STATIV(IEL,ILFV1)
     &                          +TEPA(NPT,JRPOI)*CONCEN*CONCEN
C
C     * Moyenne et variance du temps de séjour
C
              ILTS1 = ILTS+IZCL*NVLSTA
              STATIS(IEL,ILTS1) = STATIS(IEL,ILTS1)
     &          + TEPA(NPT,JRPOI) *  TEPA(NPT,JRTSP)
C
              ILTS1 = ILTS+IZCL*(NVLSTA-1)
              STATIV(IEL,ILTS1) = STATIV(IEL,ILTS1)
     &          + TEPA(NPT,JRPOI)*TEPA(NPT,JRTSP)*TEPA(NPT,JRTSP)
C
C     * Somme du poids statistiques associé aux particules
C
              ILPD1 = ILPD+IZCL*NVLSTA
              STATIS(IEL,ILPD1) = STATIS(IEL,ILPD1) + TEPA(NPT,JRPOI)
C
             IF ( IPHYLA .EQ. 1 ) THEN
C
C     * Moyenne et variance de la temperature
C
                IF ( ITPVAR .EQ. 1 ) THEN
C
                  ILTP1 = ILTP+IZCL*NVLSTA
                  STATIS(IEL,ILTP1) = STATIS(IEL,ILTP1)
     &                         + TEPA(NPT,JRPOI) * ETTP(NPT,JTP)
C
                  ILTP1 = ILTP+IZCL*(NVLSTA-1)
                  STATIV(IEL,ILTP1) = STATIV(IEL,ILTP1)
     &             + TEPA(NPT,JRPOI)*ETTP(NPT,JTP)*ETTP(NPT,JTP)
C
                ENDIF
C
C     * Moyenne et variance du diametre
C
                IF ( IDPVAR .EQ. 1 ) THEN
C
                  ILDP1 = ILDP+IZCL*NVLSTA
                  STATIS(IEL,ILDP1) = STATIS(IEL,ILDP1)
     &                         + TEPA(NPT,JRPOI) * ETTP(NPT,JDP)
C
                  ILDP1 = ILDP+IZCL*(NVLSTA-1)
                  STATIV(IEL,ILDP1) = STATIV(IEL,ILDP1)
     &             + TEPA(NPT,JRPOI)*ETTP(NPT,JDP)*ETTP(NPT,JDP)
C
                ENDIF
C
C     * Moyenne et variance de la masse
C
                IF ( IMPVAR .EQ. 1 ) THEN
C
                  ILMP1 = ILMP+IZCL*NVLSTA
                  STATIS(IEL,ILMP1) = STATIS(IEL,ILMP1)
     &                         + TEPA(NPT,JRPOI) * ETTP(NPT,JMP)
C
                  ILMP1 = ILMP+IZCL*(NVLSTA-1)
                  STATIV(IEL,ILMP1) = STATIV(IEL,ILMP1)
     &             + TEPA(NPT,JRPOI)*ETTP(NPT,JMP)*ETTP(NPT,JMP)
C
                ENDIF
C
              ELSE IF ( IPHYLA .EQ. 2 ) THEN
C
C     * Moyenne et variance de la temperature
C     * Moyenne et variance de la masse de charbon reactif
C     * Moyenne et variance de la masse de coke
C     * Moyenne et variance du diametre du coeur retrecissant
C
                ILHP1  = ILHP  +IZCL*NVLSTA
                ILMCH1 = ILMCH1+IZCL*NVLSTA
                ILMCK1 = ILMCK1+IZCL*NVLSTA
                ILDCK1 = ILDCK1+IZCL*NVLSTA
C
                STATIS(IEL,ILHP1)  = STATIS(IEL,ILHP1)
     &                        + TEPA(NPT,JRPOI) * ETTP(NPT,JHP)
                STATIS(IEL,ILMCH1) = STATIS(IEL,ILMCH1)
     &                        + TEPA(NPT,JRPOI) * ETTP(NPT,JMCH)
                STATIS(IEL,ILMCK1) = STATIS(IEL,ILMCK1)
     &                        + TEPA(NPT,JRPOI) * ETTP(NPT,JMCK)
                STATIS(IEL,ILDCK1) = STATIS(IEL,ILDCK1)
     &                        + TEPA(NPT,JRPOI) * TEPA(NPT,JRDCK)
C
                ILHP1  = ILHP  +IZCL*(NVLSTA-1)
                ILMCH1 = ILMCH1+IZCL*(NVLSTA-1)
                ILMCK1 = ILMCK1+IZCL*(NVLSTA-1)
                ILDCK1 = ILDCK1+IZCL*(NVLSTA-1)
C
                STATIV(IEL,ILHP1)  = STATIV(IEL,ILHP1)
     &          + TEPA(NPT,JRPOI)*ETTP(NPT,JHP)*ETTP(NPT,JHP)
                STATIV(IEL,ILMCH1) = STATIV(IEL,ILMCH1)
     &          + TEPA(NPT,JRPOI)*ETTP(NPT,JMCH)*ETTP(NPT,JMCH)
                STATIV(IEL,ILMCK1) = STATIV(IEL,ILMCK1)
     &          + TEPA(NPT,JRPOI)*ETTP(NPT,JMCK)*ETTP(NPT,JMCK)
                STATIV(IEL,ILDCK1) = STATIV(IEL,ILDCK1)
     &          + TEPA(NPT,JRPOI)*TEPA(NPT,JRDCK)*TEPA(NPT,JRDCK)
C
              ENDIF
C
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C
C====
C FIN
C====
C
      END
c@z
