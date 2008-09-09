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
                        SUBROUTINE LAGES2
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &   NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IBORD  , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , STATIS , TAUP   , TLAG   , PIIL   ,
     &   TSUF   , TSUP   , BX     , TSFEXT ,
     &   VAGAUS , AUXL   , GRADPR , GRADVF ,
     &   ROMP   , BRGAUS , TERBRU , FEXTLA ,
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
CFONC     INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 2 (sche2c)
CFONC
CFONC      Lorsqu'il y a eu interaction avec une face de bord,
CFONC      les calculs de la vitesse de la particule et de
CFONC      la vitesse du fluide vu sont forcement a l'ordre 1
CFONC      (meme si on est a l'ordre 2 par ailleurs).
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
CARGU ! IBORD        ! TE ! <-  ! SI NORDRE=2, CONTIENT LE NUMERO DE LA!
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
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
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT ET PREC)!
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
CARGU ! STATIS       ! TR !  -> !  CUMUL DES STATISTIQUES VOLUMIQUES   !
CARGU !(NCELET,NVLSTA!    !     !                                      !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! PIIL(NBPMAX,3! TR !  -> ! TERME DANS L'INTEGRATION DES EDS Up  !
CARGU ! TSUP(NBPMAX,3! TR !  -> ! PREDICTION 1ER SOUS-PAS POUR         !
CARGU !              !    !     !   LA VITESSE DES PARTICULES          !
CARGU ! TSUF(NBPMAX,3! TR !  -> ! PREDICTION 1ER SOUS-PAS POUR         !
CARGU !              !    !     !   LA VITESSE DU FLUIDE VU            !
CARGU ! BX(NBPMAX,3,2! TR !  -> ! CARACTERISTIQUES DE LA TURBULENCE    !
CARGU ! TSFEXT(NBPMAX! TR ! <-  ! INFOS POUR COUPLAGE RETOUR DYNAMIQUE !
CARGU ! VAGAUS       ! TR !  -> ! VARIABLES ALEATOIRES GAUSSIENNES     !
CARGU !(NBPMAX,NVGAUS!    !     !                                      !
CARGU ! AUXL         ! TR !  -  ! TABLEAU AUXILIAIRE                   !
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
      INTEGER          ITEPA(NBPMAX,NIVEP) , IBORD(NBPMAX)
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
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP) , STATIS(NCELET,*)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3)
      DOUBLE PRECISION PIIL(NBPMAX,3) , BX(NBPMAX,3,2)
      DOUBLE PRECISION TSUF(NBPMAX,3) , TSUP(NBPMAX,3)
      DOUBLE PRECISION TSFEXT(NBPMAX)
      DOUBLE PRECISION VAGAUS(NBPMAX,*) , AUXL(NBPMAX,7)
      DOUBLE PRECISION GRADPR(NCELET,3) , GRADVF(NCELET,9)
      DOUBLE PRECISION ROMP(NBPMAX)
      DOUBLE PRECISION BRGAUS(NBPMAX,*) , TERBRU(NBPMAX)
      DOUBLE PRECISION FEXTLA(NBPMAX,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL , IP , ID , I0 , IROMF , IPHAS
      DOUBLE PRECISION AUX0 , AUX1 , AUX2 , AUX3 , AUX4 , AUX5
      DOUBLE PRECISION AUX6 , AUX7 , AUX8 , AUX9 , AUX10 , AUX11
      DOUBLE PRECISION AUX12 , AUX14 , AUX15 , AUX16
      DOUBLE PRECISION AUX17 , AUX18 , AUX19 , AUX20
      DOUBLE PRECISION TER1 , TER2 , TER3 , TER4 , TER5
      DOUBLE PRECISION SIGE , TAPN
      DOUBLE PRECISION GAMMA2 , OMEGAM , OMEGA2
      DOUBLE PRECISION GRGAM2 , GAGAM , GAOME
      DOUBLE PRECISION P11 , P21 , P22 , P31 , P32 , P33
      DOUBLE PRECISION GRAV(3) , ROM , VITF
      DOUBLE PRECISION TBRIU
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
C=======================================================================
C 2.1 CALCUL A CHAQUE SOUS PAS DE TEMPS
C=======================================================================
C
C---> Calcul de tau_p*A_p et de II*TL+<u> :
C     -------------------------------------
C
      DO ID = 1,3
C
        DO IP = 1,NBPART
C
          IF (ITEPA(IP,JISOR).GT.0) THEN
C
            IEL = ITEPA(IP,JISOR)
C
            ROM = PROPCE(IEL,IROMF)
C
            AUXL(IP,ID) =
     &        ( ROM *GRADPR(IEL,ID) /ROMP(IP)
     &         +GRAV(ID)+FEXTLA(IP,ID)       ) * TAUP(IP)
C
            IF (NOR.EQ.1) THEN
              IF (ID.EQ.1) VITF = RTPA(IEL,IU(IPHAS))
              IF (ID.EQ.2) VITF = RTPA(IEL,IV(IPHAS))
              IF (ID.EQ.3) VITF = RTPA(IEL,IW(IPHAS))
            ELSE
              IF (ID.EQ.1) VITF = RTP(IEL,IU(IPHAS))
              IF (ID.EQ.2) VITF = RTP(IEL,IV(IPHAS))
              IF (ID.EQ.3) VITF = RTP(IEL,IW(IPHAS))
            ENDIF
            AUXL(IP,ID+3) = PIIL(IP,ID) * TLAG(IP,ID) + VITF
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C 2.2 ETAPE DE PREDICTION :
C=======================================================================
C
      IF (NOR.EQ.1) THEN
C
C
C---> Sauvegarde de tau_p^n
C     ---------------------
C
        DO IP = 1,NBPART
C
          IF (ITEPA(IP,JISOR).GT.0) THEN
C
            AUXL(IP,7) = TAUP(IP)
C
          ENDIF
C
        ENDDO
C
C---> Sauvegarde couplage
C     -------------------
C
        IF (IILAGR.EQ.2) THEN
C
          DO IP = 1,NBPART
C
            IF (ITEPA(IP,JISOR).GT.0) THEN
C
              AUX0 = -DTP/TAUP(IP)
              AUX1 = EXP(AUX0)
              TSFEXT(IP) = TAUP(IP) * ETTP(IP,JMP)
     &                   * (-AUX1 + (AUX1-1.D0) / AUX0)
C
            ENDIF
C
          ENDDO
C
        ENDIF
C
C
C---> Chargement des termes a t = t_n :
C     ---------------------------------
C
        DO ID = 1,3
C
          I0 = ID - 1
C
          DO IP = 1,NBPART
C
            IF (ITEPA(IP,JISOR).GT.0) THEN
C
              AUX0 = -DTP / TAUP(IP)
              AUX1 = -DTP / TLAG(IP,ID)
              AUX2 = EXP(AUX0)
              AUX3 = EXP(AUX1)
              AUX4 = TLAG(IP,ID) / (TLAG(IP,ID) - TAUP(IP))
              AUX5 = AUX3 - AUX2
C
              TSUF(IP,ID) = 0.5D0 * ETTPA(IP,JUF+I0) * AUX3
     &                    + AUXL(IP,ID+3) * ( -AUX3 + (AUX3-1.D0) /AUX1)
C
              TER1 = 0.5D0 * ETTPA(IP,JUP+I0) * AUX2
              TER2 = 0.5D0 * ETTPA(IP,JUF+I0) * AUX4 * AUX5
              TER3 = AUXL(IP,ID+3)*
     &              ( - AUX2 + ((TLAG(IP,ID) + TAUP(IP)) / DTP)
     &             * (1.D0 - AUX2)
     &             - (1.D0 + TLAG(IP,ID) / DTP) * AUX4 * AUX5)
              TER4 = AUXL(IP,ID) * (-AUX2 + (AUX2 - 1.D0) / AUX0)
C
              TSUP(IP,ID) = TER1 + TER2 + TER3 + TER4
C
            ENDIF
C
          ENDDO
C
        ENDDO
C
C---> Schema d'Euler :
C     ----------------
C
        CALL LAGES1
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &   NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , STATIS , TAUP   , TLAG   , PIIL   ,
     &   BX     , VAGAUS , GRADPR , GRADVF , ROMP   ,
     &   BRGAUS , TERBRU , FEXTLA ,
     &   RDEVEL , RTUSER , RA     )
C
      ELSE
C
C=======================================================================
C 2.2 ETAPE DE CORRECTION :
C=======================================================================
C
C---> Calcul de Us :
C     --------------
C
        DO ID = 1,3
C
          I0 = ID - 1
C
          DO IP = 1,NBPART
C
            IF (ITEPA(IP,JISOR).GT.0 .AND. IBORD(IP).EQ.0) THEN
C
              AUX0 = -DTP / TAUP(IP)
              AUX1 = -DTP / TLAG(IP,ID)
              AUX2 = EXP(AUX0)
              AUX3 = EXP(AUX1)
              AUX4 = TLAG(IP,ID) / (TLAG(IP,ID) - TAUP(IP))
              AUX5 = AUX3 - AUX2
              AUX6 = AUX3 * AUX3
C
              TER1 = 0.5D0 * ETTPA(IP,JUF+I0) * AUX3
              TER2 = AUXL(IP,ID+3) * (1.D0 - (AUX3 - 1.D0) / AUX1)
C
              TER3 = -AUX6 + (AUX6 - 1.D0) / (2.D0 * AUX1)
              TER4 = 1.D0 - (AUX6 - 1.D0) / (2.D0 * AUX1)
              SIGE = (TER3 * BX(IP,ID,NOR-1) + TER4 * BX(IP,ID,NOR))
     &             * (1.D0 / (1.D0 - AUX6) )
C
              TER5 = 0.5D0 * TLAG(IP,ID) * (1.D0 - AUX6)
C
              ETTP(IP,JUF+I0) = TSUF(IP,ID) + TER1 + TER2
     &                            + SIGE * SQRT(TER5) * VAGAUS(IP,ID)
C
C---> Calcul de Up :
C     --------------
C
              TER1 = 0.5D0 * ETTPA(IP,JUP+I0) * AUX2
              TER2 = 0.5D0 * ETTPA(IP,JUF+I0) * AUX4 * AUX5
              TER3 = AUXL(IP,ID+3)*
     &           ( 1.D0 - ((TLAG(IP,ID) + TAUP(IP)) /DTP) *(1.D0-AUX2)
     &             + (TLAG(IP,ID) / DTP) * AUX4 * AUX5)
     &             + AUXL(IP,ID) * (1.D0 - (AUX2 - 1.D0) / AUX0)
C
              TAPN  = AUXL(IP,7)
              AUX7  = EXP(-DTP / TAPN)
              AUX8  = 1.D0 - AUX3 * AUX7
              AUX9  = 1.D0 - AUX6
              AUX10 = 1.D0 - AUX7 * AUX7
C
              AUX11 = TAPN / (TLAG(IP,ID) + TAPN)
              AUX12 = TLAG(IP,ID) / (TLAG(IP,ID) - TAPN)
              AUX14 = TLAG(IP,ID) - TAPN
C
              AUX15 = TLAG(IP,ID) * (1.D0-AUX3)
              AUX16 = TAPN * (1.D0 - AUX7)
              AUX17 = SIGE * SIGE * AUX12 * AUX12
C
              AUX18 = 0.5D0 * TLAG(IP,ID) * AUX9
              AUX19 = 0.5D0 * TAPN * AUX10
              AUX20 = TLAG(IP,ID) * AUX11 * AUX8
C
C ---> calcul de la matrice de correlation
C
              GAMMA2 = SIGE * SIGE * AUX18
C
              GRGAM2 = AUX17 * (AUX18 - 2.D0 * AUX20 + AUX19)
C
              GAGAM = SIGE * SIGE * AUX12 * (AUX18 - AUX20)
C
              OMEGA2 = AUX17 * (AUX14 *
     &               (AUX14 * DTP - 2.D0 * TLAG(IP,ID) * AUX15
     &               + 2.D0 * TAPN * AUX16)
     &               + TLAG(IP,ID) * TLAG(IP,ID) * AUX18
     &               + TAPN * TAPN * AUX19
     &               - 2.D0 * TLAG(IP,ID) * TAPN * AUX20)
C
              OMEGAM = AUX14 * (1.D0-AUX3) - AUX18 + TAPN * AUX11 * AUX8
              OMEGAM = OMEGAM * SIGE * SIGE * AUX12 * TLAG(IP,ID)
C
              GAOME  = AUX17 * (AUX14 * (AUX15 - AUX16) - TLAG(IP,ID)
     &            * AUX18 - TAPN * AUX19 + TAPN * TLAG(IP,ID) * AUX8 )
C
C ---> simulation du vecteur Gaussien
C
              P11 = SQRT( MAX (ZERO,GAMMA2) )

              IF (P11.GT.EPZERO) THEN
                P21 = OMEGAM / P11
                P22 = OMEGA2 - P21 * P21
                P22 = SQRT( MAX (ZERO,P22) )
              ELSE
                P21 = 0.D0
                P22 = 0.D0
              ENDIF
C
              IF (P11.GT.EPZERO) THEN
                P31 = GAGAM / P11
              ELSE
                P31 = 0.D0
              ENDIF
C
              IF (P22.GT.EPZERO) THEN
                P32 = (GAOME - P31 * P21) / P22
              ELSE
                P32 = 0.D0
              ENDIF
C
              P33 = GRGAM2 - P31 * P31 - P32 * P32
              P33 = SQRT( MAX(ZERO,P33) )
C
              TER4 = P31 * VAGAUS(IP,ID)
     &             + P32 * VAGAUS(IP,3+ID) + P33 * VAGAUS(IP,6+ID)
C
C ---> Calcul des Termes dans le cas du mouvement Brownien
C
              IF ( LAMVBR .EQ. 1 ) THEN
                TBRIU = TERBRU(IP)*BRGAUS(IP,3+ID)
              ELSE
                TBRIU = 0.D0
              ENDIF
C
C ---> finalisation de l'ecriture
C
              ETTP(IP,JUP+I0) = TSUP(IP,ID) +TER1 +TER2 +TER3 +TER4
     &                        + TBRIU
C
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C----
C FIN
C----
C
      END
c@z
