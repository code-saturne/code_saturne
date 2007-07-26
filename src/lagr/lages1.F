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
                        SUBROUTINE LAGES1
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &   NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , STATIS ,
     &   TAUP   , TLAG   , PIIL   ,
     &   BX     , VAGAUS , GRADPR , GRADVF , ROMP   ,
     &   BRGAUS , TERBRU , FEXTLA ,
     &   RDEVEL , RTUSER , RA    )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 1 (scheO1)
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
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (PAS DE TEMPS PRECEDENT) !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR ! <-  ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR !  -> ! CUMUL DES STATISTIQUES VOLUMIQUES    !
CARGU !(NCELET,NVLSTA!    !     !                                      !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! PIIL(NBPMAX,3! TR !  -> ! TERME DANS L'INTEGRATION DES EDS Up  !
CARGU ! BX(NBPMAX,3,2! TR !  -> ! CARACTERISTIQUES DE LA TURBULENCE    !
CARGU ! VAGAUS       ! TR !  -> ! VARIABLES ALEATOIRES GAUSSIENNES     !
CARGU !(NBPMAX,NVGAUS!    !     !                                      !
CARGU ! GRADPR(NCEL,3! TR !  -> ! GRADIENT DE PRESSION                 !
CARGU ! GRADVF(NCEL,3! TR !  -> ! GRADIENT DE LA VITESSE DU FLUIDE     !
CARGU ! ROMP         ! TR !  -> ! MASSE VOLUMIQUE DES PARTICULES       !
CARGU ! FEXTLA       ! TR !  -> ! CHAMP DE FORCES EXTERIEUR            !
CARGU !(NCELET,3)    !    !     !    UTILISATEUR (M/S2)                !
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   ,LNDFAC , LNDFBR , NCELBR
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
      DOUBLE PRECISION DT(NCELET) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP) , STATIS(NCELET,*)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3)
      DOUBLE PRECISION PIIL(NBPMAX,3) , BX(NBPMAX,3,2)
      DOUBLE PRECISION VAGAUS(NBPMAX,*)
      DOUBLE PRECISION BRGAUS(NBPMAX,*) , TERBRU(NBPMAX)
      DOUBLE PRECISION GRADPR(NCELET,3) , GRADVF(NCELET,9)
      DOUBLE PRECISION ROMP(NBPMAX)
      DOUBLE PRECISION FEXTLA(NBPMAX,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL , IP , ID , I0 , IROMF , IPHAS , MODE
C
      DOUBLE PRECISION AA , BB , CC , DD , EE
      DOUBLE PRECISION AUX1 , AUX2 ,AUX3 , AUX4 , AUX5 , AUX6
      DOUBLE PRECISION AUX7 , AUX8 , AUX9 , AUX10 , AUX11
      DOUBLE PRECISION TER1F , TER2F , TER3F
      DOUBLE PRECISION TER1P , TER2P , TER3P , TER4P , TER5P
      DOUBLE PRECISION TER1X , TER2X , TER3X , TER4X , TER5X
      DOUBLE PRECISION TCI , FORCE
      DOUBLE PRECISION GAMA2 , OMEGAM , OMEGA2
      DOUBLE PRECISION GRGA2 , GAGAM , GAOME
      DOUBLE PRECISION P11 , P21 , P22 , P31 , P32 , P33
      DOUBLE PRECISION GRAV(3) , ROM , VITF
      DOUBLE PRECISION TEMPF
      DOUBLE PRECISION DDBR, TIX2, TIU2, TIXIU
      DOUBLE PRECISION TBRIX1, TBRIX2, TBRIU
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
C 1. INITIALISATIONS
C=======================================================================
C
      IPHAS = ILPHAS
C
      GRAV(1) = GX
      GRAV(2) = GY
      GRAV(3) = GZ
C
C Pointeur sur la masse volumique en fonction de l'ecoulement
C
      IF ( IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICFUEL).GE.0 ) THEN
        IROMF = IPPROC(IROM1)
      ELSE
        IROMF = IPPROC(IROM(IPHAS))
      ENDIF
C
C=======================================================================
C 2. INTEGRATION DES EDS SUR LES PARTICULES
C=======================================================================
C
      DO ID = 1,3
C
        I0 = ID - 1
C
        DO IP = 1,NBPART
C
          IF (ITEPA(IP,JISOR).GT.0) THEN
C
            IEL = ITEPA(IP,JISOR)
C
            ROM = PROPCE(IEL,IROMF)
C
            IF (ID.EQ.1) VITF = RTPA(IEL,IU(IPHAS))
            IF (ID.EQ.2) VITF = RTPA(IEL,IV(IPHAS))
            IF (ID.EQ.3) VITF = RTPA(IEL,IW(IPHAS))
C
C---> (2.1) Calcul preliminaires :
C     ----------------------------
C
C      calcul de II*TL+<u> et [(grad<P>/rhop+g)*tau_p+<Uf>] ?
C
            TCI = PIIL(IP,ID) * TLAG(IP,ID) + VITF
C
            FORCE =
     &        ( ROM * GRADPR(IEL,ID) / ROMP(IP)
     &         +GRAV(ID)+ FEXTLA(IP,ID)        ) * TAUP(IP)
C
C---> (2.2) Calcul des coefficients/termes deterministes :
C     ----------------------------------------------------
C
C  sur HP-UX : l'option de compil +T genere des messages du type :
C  *** MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW
C  au passage de l'exponentiel si l'exposant est < a -7.0D2
C
C
            AUX1 = EXP( -DTP / TAUP(IP))
            AUX2 = EXP( -DTP / TLAG(IP,ID))
            AUX3 = TLAG(IP,ID) / (TLAG(IP,ID)-TAUP(IP))
C
            AUX4 = TLAG(IP,ID) / (TLAG(IP,ID)+TAUP(IP))
            AUX5 = TLAG(IP,ID) * (1.D0-AUX2)
            AUX6 = BX(IP,ID,NOR) * BX(IP,ID,NOR) * TLAG(IP,ID)

            AUX7 = TLAG(IP,ID) - TAUP(IP)
            AUX8 = BX(IP,ID,NOR) * BX(IP,ID,NOR) * AUX3**2
C
C---> termes pour la trajectoire
C
            AA = TAUP(IP) * (1.D0 - AUX1)
            BB = (AUX5 - AA) * AUX3
            CC = DTP - AA - BB
C
            TER1X = AA * ETTPA(IP,JUP+I0)
            TER2X = BB * ETTPA(IP,JUF+I0)
            TER3X = CC * TCI
            TER4X = (DTP - AA) * FORCE
C
C---> termes pour le fluide vu
C
            TER1F = ETTPA(IP,JUF+I0) * AUX2
            TER2F = TCI * (1.D0-AUX2)
C
C---> termes pour la vitesse des particules
C
            DD = AUX3 * (AUX2 - AUX1)
            EE = 1.D0 - AUX1
C
            TER1P = ETTPA(IP,JUP+I0) * AUX1
            TER2P = ETTPA(IP,JUF+I0) * DD
            TER3P = TCI * (EE-DD)
            TER4P = FORCE * EE
C
C---> (2.3) Calcul des coefficients pour les integrales stochastiques :
C     -----------------------------------------------------------------
C
C---> integrale sur la position des particules
C
            GAMA2  = 0.5D0 * (1.D0 - AUX2*AUX2 )
            OMEGAM = 0.5D0 * AUX4 * ( AUX5 - AUX2*AA )
     &              -0.5D0 * AUX2 * BB
            OMEGAM = OMEGAM * SQRT(AUX6)
C
            OMEGA2 = AUX7
     &             * (AUX7*DTP - 2.D0 * (TLAG(IP,ID)*AUX5-TAUP(IP)*AA))
     &           + 0.5D0 * TLAG(IP,ID) * TLAG(IP,ID) * AUX5
     &             * (1.D0 + AUX2)
     &           + 0.5D0 * TAUP(IP) * TAUP(IP) * AA * (1.D0+AUX1)
     &           - 2.D0 * AUX4 * TLAG(IP,ID) * TAUP(IP) * TAUP(IP)
     &             * (1.D0 - AUX1*AUX2)
C
            OMEGA2 = AUX8 * OMEGA2
C

            IF (ABS(GAMA2).GT.EPZERO) THEN
C
              P21 = OMEGAM / SQRT(GAMA2)
              P22 = OMEGA2 - P21**2
C
              P22 = SQRT( MAX(ZERO,P22) )
C
            ELSE
              P21 = 0.D0
              P22 = 0.D0
            ENDIF
C
            TER5X = P21 * VAGAUS(IP,ID) + P22 * VAGAUS(IP,3+ID)
C
C---> integrale sur la vitesse du fluide vu
C
            P11 = SQRT( GAMA2*AUX6 )
            TER3F = P11*VAGAUS(IP,ID)
C
C---> integrale sur la vitesse des particules
C
            AUX9  = 0.5D0 * TLAG(IP,ID) * (1.D0 - AUX2*AUX2)
            AUX10 = 0.5D0 * TAUP(IP) * (1.D0 - AUX1*AUX1)
            AUX11 = TAUP(IP) * TLAG(IP,ID) * (1.D0 - AUX1*AUX2)
     &               / (TAUP(IP) + TLAG(IP,ID))
C
            GRGA2 = (AUX9 - 2.D0*AUX11 + AUX10) * AUX8
            GAGAM = (AUX9 - AUX11) * (AUX8 / AUX3)
            GAOME = ( (TLAG(IP,ID) - TAUP(IP)) * (AUX5 - AA)
     &             - TLAG(IP,ID) * AUX9 - TAUP(IP) * AUX10
     &             + (TLAG(IP,ID) + TAUP(IP)) * AUX11) * AUX8
C
            IF(P11.GT.EPZERO) THEN
              P31 = GAGAM / P11
            ELSE
              P31 = 0.D0
            ENDIF
C
            IF(P22.GT.EPZERO) THEN
              P32 = (GAOME-P31*P21) / P22
            ELSE
              P32 = 0.D0
            ENDIF
C
            P33 = GRGA2 - P31**2 - P32**2
C
            P33 = SQRT( MAX(ZERO,P33) )
C
            TER5P = P31 * VAGAUS(IP,ID)
     &            + P32 * VAGAUS(IP,3+ID)
     &            + P33 * VAGAUS(IP,6+ID)
C
C
C---> (2.3) Calcul des Termes dans le cas du mouvement Brownien :
C     ---------------------------------------------------------
C
            IF ( LAMVBR .EQ. 1 ) THEN
C
C             Calcul de la temperature du fluide en fonction du type
C             d'ecoulement
C
              IF ( IPPMOD(ICP3PL).GE.0 .OR.
     &             IPPMOD(ICPL3C).GE.0      ) THEN
C
                TEMPF = PROPCE(IEL,IPPROC(ITEMP1))
C
              ELSE IF ( IPPMOD(ICOD3P).GE.0 .OR.
     &          IPPMOD(ICOEBU).GE.0 .OR.
     &          IPPMOD(IELARC).GE.0 .OR.
     &          IPPMOD(IELJOU).GE.0      ) THEN
C
                TEMPF = PROPCE(IEL,IPPROC(ITEMP))
C
              ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.-1 ) THEN
                TEMPF = RTPA(IEL,ISCA(ISCALT(IPHAS)))
C
              ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.1 ) THEN
                TEMPF = RTPA(IEL,ISCA(ISCALT(IPHAS)))
C
              ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.2 ) THEN
                MODE = 1
                CALL USTHHT(MODE,RTPA(IEL,ISCA(ISCALT(IPHAS))),TEMPF)
C               ===========
                TEMPF = TEMPF+TKELVI
              ELSE
                TEMPF = T0(IPHAS)
              ENDIF
C
              DDBR  = SQRT( 2.D0*KBOLTZ*TEMPF
     &                     /(ETTP(IP,JMP)*TAUP(IP)) )
              TIX2  = (TAUP(IP)*DDBR)**2
     &                *(DTP - TAUP(IP)*(1.D0-AUX1)
     &                               *(3.D0-AUX1) /2.D0 )
C
              TIU2  = DDBR*DDBR*TAUP(IP)
     &               *(1.D0-EXP(-2.D0*DTP/TAUP(IP)))/2.D0
C
              TIXIU = (DDBR*TAUP(IP)*(1.D0-AUX1))**2/2.D0
C
              TBRIX2 = TIX2-(TIXIU*TIXIU)/TIU2
              IF ( TBRIX2 .GT.0.D0 ) THEN
                TBRIX2 = SQRT(TBRIX2)*BRGAUS(IP,ID)
              ELSE
                TBRIX2 = 0.D0
              ENDIF
C
              IF ( TIU2 .GT. 0.D0 ) THEN
                TBRIX1 = TIXIU/SQRT(TIU2)*BRGAUS(IP,3+ID)
              ELSE
                TBRIX1 = 0.D0
              ENDIF
C
              IF ( TIU2 .GT. 0.D0 ) THEN
                TBRIU = SQRT(TIU2)*BRGAUS(IP,3+ID)
                TERBRU(IP) = SQRT(TIU2)
              ELSE
                TBRIU = 0.D0
                TERBRU(IP) = 0.D0
              ENDIF
C
            ELSE
              TBRIX1 = 0.D0
              TBRIX2 = 0.D0
              TBRIU  = 0.D0
            ENDIF
C
C=======================================================================
C 3. FINALISATION DES ECRITURES
C=======================================================================
C
C---> trajectoire
C
            ETTP(IP,JXP+I0) = ETTPA(IP,JXP+I0)
     &                       + TER1X + TER2X + TER3X + TER4X + TER5X
     &                       + TBRIX1 + TBRIX2
C
C---> vitesse fluide vu
C
            ETTP(IP,JUF+I0) = TER1F + TER2F + TER3F
C
C---> vitesse particules
C
            ETTP(IP,JUP+I0) = TER1P + TER2P + TER3P + TER4P + TER5P
     &                       + TBRIU
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C----
C FIN
C----
C
      END
c@z
