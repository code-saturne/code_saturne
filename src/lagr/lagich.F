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
                        SUBROUTINE LAGICH
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   ITEPA  , IBORD  , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT , TSVAR  ,
     &   CPGD1  , CPGD2  , CPGHT  ,
     &   SKP1   , SKP2   , SKGLOB ,
     &   GAMHET , DELTAH , TEMPF  ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC      INTEGRATION DES EDS POUR LE CHARBON
CFONC
CFONC         - Temperature              (JHP)
CFONC         - Masse de charbon reactif (JMCH)
CFONC         - Masse de coke            (JMCK)
CFONC
CFONC      ET CALCUL DU DIAMETRE DU COEUR RETRECISSANT (JRDCK)
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
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IBORD        ! TE !  -> ! CONTIENT LE NUMERO DE LA             !
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
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
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
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
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! TEMPCT       ! TR !  -> ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU !  (NBPMAX,2)  !    !     !                                      !
CARGU ! TSVAR        ! TR ! <-> ! PREDICTION 1ER SOUS-PAS POUR LA      !
CARGU ! (NBPMAX,NVP1)!    !     !   VARIABLE IVAR, UTILISE POUR LA     !
CARGU !              !    !     !   CORRECTION AU 2EME SOUS-PAS        !
CARGU ! CPGD1,CPGD2, ! TR ! <-  ! TERMES DE DEVOLATILISATION 1 ET 2 ET !
CARGU !  CPGHT(NBPMAX!    !     !   DE COMBUSION HETEROGENE (CHARBON   !
CARGU !              !    !     !   AVEC COUPLAGE RETOUR THERMIQUE)    !
CARGU ! SK1,SK2,     ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU ! SKGLOB(NBPMAX!    !     !                                      !
CARGU ! GAMHET(NBPMAX! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! DELTAH(NBPMAX! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TEMPF(NCELET)! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "radiat.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP , NVP1 , NVEP , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          ITEPA(NBPMAX,NIVEP) , IBORD(NBPMAX)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3) , TEMPCT(NBPMAX,2)
      DOUBLE PRECISION TSVAR(NBPMAX,NVP1)
      DOUBLE PRECISION SKP1(NBPMAX) , SKP2(NBPMAX) , SKGLOB(NBPMAX)
      DOUBLE PRECISION GAMHET(NBPMAX) , DELTAH(NBPMAX) , TEMPF(NCELET)
      DOUBLE PRECISION CPGD1(NBPMAX), CPGD2(NBPMAX), CPGHT(NBPMAX)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          NPT , IEL , ICHA , MODE , IGE
      INTEGER          IROMF , IPHAS
      DOUBLE PRECISION AUX1 , AUX2 , AUX3 , AUX4 , AUX5 , AUX6
      DOUBLE PRECISION TER1 , TER2 , TER3 , DIAMP2, D2
      DOUBLE PRECISION TPK , TFK , SKC , SKDD , SE , PO2
      DOUBLE PRECISION HO2TF , HCTP , HCOTP , DEN , SHERW
      DOUBLE PRECISION COEF , MP0 , D6SPI , DPIS6 , D1S3 , D2S3
      DOUBLE PRECISION GAMDV1(NCHARM2) , GAMDV2(NCHARM2)
      DOUBLE PRECISION F1MC(NCHARM2) , F2MC(NCHARM2)
      DOUBLE PRECISION COEFE(NGAZEM)
C
      DOUBLE PRECISION PRECIS
      PARAMETER ( PRECIS = 1.D-15 )
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      IPHAS = ILPHAS
C
      D6SPI = 6.D0 / PI
      DPIS6 = PI / 6.D0
      D1S3  = 1.D0 / 3.D0
      D2S3  = 2.D0 / 3.D0
C
C --- Si couplage retour thermique :
C
      IF ( LTSTHE.EQ.1 .AND. NOR.EQ.1 ) THEN
        DO NPT = 1,NBPART
          CPGD1(NPT) = 0.D0
          CPGD2(NPT) = 0.D0
          CPGHT(NPT) = 0.D0
        ENDDO
      ENDIF
      IF ( LTSTHE.EQ.1 ) THEN
        COEF = 1.D0 / DBLE(NORDRE)
      ENDIF
C
C=======================================================================
C 2. Pointeur sur la masse volumique en fonction de l'ecoulement
C=======================================================================
C
      IF ( IPPMOD(ICP3PL).GE.0 ) THEN
        IROMF = IPPROC(IROM1)
      ELSE
        IROMF = IPPROC(IROM(IPHAS))
      ENDIF
C
C=======================================================================
C 3. Temperature moyenne Fluide en Kelvin dans le cas ou la phase
C    porteuse est une flamme de charbon pulverise
C=======================================================================
C
      IF ( IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICPL3C).GE.0 ) THEN
C
         DO IEL = 1,NCEL
           TEMPF(IEL) = PROPCE(IEL,IPPROC(ITEMP1))
         ENDDO
C
      ELSE
        WRITE(NFECRA,1000) IPHYLA, IPPMOD(ICPL3C), IPPMOD(ICP3PL)
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C 4. Calcul des constantes de vitesses SPK1 et SPK2 du transfert
C    de masse par devolatilisation avec des lois d'Arrhenius
C=======================================================================
C
C     RR --> Constante des gaz parfaits en J/mol/K
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
          ICHA = ITEPA(NPT,JINCH)
          TPK = ETTP(NPT,JHP) + TKELVI
          AUX1 = 1.D0 / (RR*TPK)
          SKP1(NPT) = A1CH(ICHA) * EXP( -E1CH(ICHA) * AUX1)
          SKP2(NPT) = A2CH(ICHA) * EXP( -E2CH(ICHA) * AUX1)
        ENDIF
      ENDDO
C
C=======================================================================
C 5. Calcul de la masse volumique du coke
C    On suppose pour le calcul de la masse volumique du coke que
C    la devolatilisation a lieu a volume constant
C=======================================================================
C
C --- Initialisation
C
      DO ICHA = 1,NCHARM
        GAMDV1(ICHA) = 0.D0
        GAMDV2(ICHA) = 0.D0
        RHOCK(ICHA) = RHO0CH(ICHA)
      ENDDO
C
C --- Calcul de l'integrale de GAMDV1 et GAMDV2 pour chaque charbon
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
          ICHA = ITEPA(NPT,JINCH)
          AUX1 = SKP1(NPT) * Y1CH(ICHA) * ETTP(NPT,JMCH)
          AUX2 = SKP2(NPT) * Y2CH(ICHA) * ETTP(NPT,JMCH)
C
          GAMDV1(ICHA) = GAMDV1(ICHA) + AUX1
          GAMDV2(ICHA) = GAMDV2(ICHA) + AUX2
C
C --- Couplage retour thermique
C
          IF ( LTSTHE.EQ.1 ) THEN
            CPGD1(NPT) = CPGD1(NPT) + COEF*AUX1
            CPGD2(NPT) = CPGD2(NPT) + COEF*AUX2
          ENDIF
C
        ENDIF
      ENDDO
C
C --- Calcul de la masse volumique moyenne du coke
C
      DO ICHA = 1,NCHARB
        DEN = Y2CH(ICHA)*GAMDV1(ICHA) + Y1CH(ICHA)*GAMDV2(ICHA)
        IF ( DEN.GT.PRECIS ) THEN
          RHOCK(ICHA) = RHO0CH(ICHA)
     &      *( Y2CH(ICHA)*GAMDV1(ICHA)+Y1CH(ICHA)*GAMDV2(ICHA)
     &        -Y1CH(ICHA)*Y2CH(ICHA)*(GAMDV1(ICHA)+GAMDV2(ICHA)) )
     &       / DEN
        ENDIF
      ENDDO
C
C=======================================================================
C 6. Calcul du diametre du coeur retrecissant
C=======================================================================
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
          ICHA = ITEPA(NPT,JINCH)
C
          TEPA(NPT,JRDCK) =
     &             ( (D6SPI / ( 1.D0-XASHCH(ICHA)) )
     &              *( ETTP(NPT,JMCH)/RHO0CH(ICHA)
     &              +ETTP(NPT,JMCK)/RHOCK(ICHA) ) )**D1S3
        ENDIF
      ENDDO
C
C=======================================================================
C 7. Calcul de la constante globale de reaction
C=======================================================================
C
C ---  Hypothese Sherwood = 2.
C
      SHERW = 2.D0
C
      DO NPT = 1,NBPART
C
        IF (ITEPA(NPT,JISOR).GT.0) THEN
C
          ICHA = ITEPA(NPT,JINCH)
C
          TPK = ETTP(NPT,JHP) + TKELVI
C
C --- Coefficient de cinetique chimique de formation de CO
C       en (kg.m-2.s-1.atm(-n))
C
          SKC = AHETCH(ICHA)
     &      * EXP(-EHETCH(ICHA)*4185.D0 / (RR*TPK) )
C
C --- Coefficient de diffusion en  (Kg/m2/s/atm) et constante
C     globale de reaction
C
          IF ( TEPA(NPT,JRDCK).GT.EPSICP ) THEN
            SKDD = SHERW * 2.53D-7 * (TPK**0.75D0) / TEPA(NPT,JRDCK)
            SKGLOB(NPT) = (SKC*SKDD) / (SKC+SKDD)
          ELSE
            SKGLOB(NPT) = SKC
          ENDIF
C
        ENDIF
      ENDDO
C
C=======================================================================
C 8. Calcul de la GAMMAhet , GAMMACH et 0.5(MO2/MC)*(HO2(Tp)-HO2(TF))
C=======================================================================
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
C
          ICHA = ITEPA(NPT,JINCH)
          IEL  = ITEPA(NPT,JISOR)
C
          TPK = ETTP(NPT,JHP) + TKELVI
          TFK = ETTP(NPT,JTF) + TKELVI
C
C --- Calcul de la pression partielle en oxygene (atm)
C                                                 ---
C       PO2 = RHO1*RR*T*YO2/MO2
C
          PO2 = PROPCE(IEL,IROMF) * RR * TEMPF(IEL)
     &        * PROPCE(IEL,IPPROC(IYM1(IO2))) / WMOLE(IO2) / PREFTH
C
C --- Calcul de (Surface efficace)/(Mck**2/3) : SE
C
          SE =  ( PI*(1.D0-XASHCH(ICHA)) )**D1S3
     &        * ( 6.D0/RHOCK(ICHA)       )**D2S3
C
C --- Calcul de la GamHET/(Mck**2/3)
C
          GAMHET(NPT) = SE * PO2 * SKGLOB(NPT)
C
C --- Pas de combustion heterogene si Mch/Mp >= 1.D-3
C
          IF ( ETTPA(NPT,JMCH).GE.(1.D-3*ETTPA(NPT,JMP)) ) THEN
            GAMHET(NPT) = 0.D0
          ENDIF
C
C --- Couplage retour thermique
C
          IF ( LTSTHE.EQ.1 ) THEN
            CPGHT(NPT) = CPGHT(NPT)
     &                 + COEF * GAMHET(NPT)
     &                 * ( ETTP(NPT,JMCK)**D2S3 )
          ENDIF
C
C --- Calcul de Hc(Tp)-Mco/Mc Hco2(Tp)+0.5Mo2/Mc Ho2(Tf)
C
C        Calcul de Hcoke(TP)
C
          HCTP = H02CH(ICHA) + ETTP(NPT,JCP)*(TPK-TREFTH)
C
C        Calcul de MCO/MC HCO(TP)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(ICO) = WMOLE(ICO) / WMOLAT(IATC)
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
          MODE      = -1
          CALL CPTHP1 ( MODE , HCOTP , COEFE  , F1MC , F2MC ,  TPK )
C         ===========
C
C        Calcul de MO2/MC/2. HO2(TF)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(IO2) = WMOLE(IO2) / WMOLAT(IATC) / 2.D0
C
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
          MODE      = -1
          CALL CPTHP1 ( MODE  , HO2TF , COEFE , F1MC  , F2MC , TFK )
C         ===========
C
          DELTAH(NPT) = HCOTP - HO2TF - HCTP
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 9. Integration Masse de Charbon reactif
C=======================================================================
C
      IF (NOR.EQ.1) THEN
        DO NPT = 1,NBPART
          IF (ITEPA(NPT,JISOR).GT.0) THEN
C
            AUX1 = EXP(-(SKP1(NPT)+SKP2(NPT))*DTP)
            ETTP(NPT,JMCH) = ETTPA(NPT,JMCH)*AUX1
C
C Clipping
            IF ( ETTP(NPT,JMCH).LT.PRECIS ) THEN
              ETTP(NPT,JMCH) = 0.D0
            ENDIF
C
          ENDIF
        ENDDO
C
      ELSE IF (NOR.EQ.2) THEN
        DO NPT = 1,NBPART
          IF (ITEPA(NPT,JISOR).GT.0 .AND. IBORD(NPT).EQ.0) THEN
C
            AUX1 = EXP(-(SKP1(NPT)+SKP2(NPT))*DTP)
            ETTP(NPT,JMCH) = 0.5D0 * ( ETTP(NPT,JMCH)
     &                         + ETTPA(NPT,JMCH)*AUX1 )
C
C Clipping
            IF ( ETTP(NPT,JMCH).LT.PRECIS ) THEN
              ETTP(NPT,JMCH) = 0.D0
            ENDIF
C
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 10. Integration Masse de Coke
C=======================================================================
C
      IF (NOR.EQ.1) THEN
        DO NPT = 1,NBPART
          IF (ITEPA(NPT,JISOR).GT.0) THEN
C
            ICHA = ITEPA(NPT,JINCH)
C
            AUX1 = SKP1(NPT) * (1.D0-Y1CH(ICHA)) * ETTPA(NPT,JMCH)
     &           + SKP2(NPT) * (1.D0-Y2CH(ICHA)) * ETTPA(NPT,JMCH)
C
            IF ( ETTPA(NPT,JMCK).GT.PRECIS ) THEN
C
              AUX2 = ETTPA(NPT,JMCK)**D1S3
              AUX3 = AUX2 + D2S3 *GAMHET(NPT) *DTP
C
              TER1 = (AUX1*AUX2-GAMHET(NPT)*ETTPA(NPT,JMCK)) *DTP/AUX3
C
              TSVAR(NPT,JMCK) = 0.5D0 * TER1
C
              ETTP(NPT,JMCK) = ETTPA(NPT,JMCK) + TER1
C
            ELSE
              TER1 = AUX1 * DTP
              TSVAR(NPT,JMCK) = 0.5D0 * TER1
              ETTP(NPT,JMCK) = ETTPA(NPT,JMCK) + TER1
            ENDIF
C  Clipping
            IF ( ETTP(NPT,JMCK).LT.0.D0 ) THEN
              ETTP(NPT,JMCK) = 0.D0
            ENDIF
C
          ENDIF
        ENDDO
C
      ELSE IF (NOR.EQ.2) THEN
        DO NPT = 1,NBPART
          IF (ITEPA(NPT,JISOR).GT.0 .AND. IBORD(NPT).EQ.0) THEN
C
            ICHA = ITEPA(NPT,JINCH)
C
            AUX1 = SKP1(NPT) * (1.D0-Y1CH(ICHA)) * ETTP(NPT,JMCH)
     &           + SKP2(NPT) * (1.D0-Y2CH(ICHA)) * ETTP(NPT,JMCH)
C
            IF ( ETTPA(NPT,JMCK).GT.PRECIS ) THEN
C
              AUX2 = ETTPA(NPT,JMCK)**D1S3
              AUX3 = AUX2 + D2S3 *GAMHET(NPT) *DTP
C
              TER1 = ( AUX1 *AUX2 -GAMHET(NPT) *ETTPA(NPT,JMCK))
     &             * DTP / AUX3
C
              ETTP(NPT,JMCK) = ETTPA(NPT,JMCK)
     &                            +TSVAR(NPT,JMCK)+0.5D0*TER1
C
            ELSE
              TER1 = AUX1*DTP
              ETTP(NPT,JMCK) = ETTPA(NPT,JMCK)
     &                           + TSVAR(NPT,JMCK) + 0.5D0*TER1
            ENDIF
C Clipping
            IF ( ETTP(NPT,JMCK).LT.0.D0 ) THEN
              ETTP(NPT,JMCK) = 0.D0
            ENDIF
C
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 11. Integration de la temperature des grains de charbon
C=======================================================================
C
      IF (NOR.EQ.1) THEN
        DO NPT = 1,NBPART
          IF (ITEPA(NPT,JISOR).GT.0) THEN
C
            ICHA = ITEPA(NPT,JINCH)
            IEL   = ITEPA(NPT,JISOR)
C
            D2 = ETTP(NPT,JDP)*ETTP(NPT,JDP)
            DIAMP2 = XASHCH(ICHA)*TEPA(NPT,JRD0P)*TEPA(NPT,JRD0P)
     &              +(1.D0-XASHCH(ICHA))*TEPA(NPT,JRDCK)*TEPA(NPT,JRDCK)
C
            AUX1 = TEMPCT(NPT,1)*DIAMP2/D2
C
C    Combustion heterogene
            AUX2 = ( -GAMHET(NPT) *(ETTP(NPT,JMCK)**D2S3) *DELTAH(NPT) )
C
C    Rayonnement
C
            AUX3 = PI * DIAMP2
     &           * ( PROPCE(IEL,IPPROC(ILUMI))
     &           - 4.D0*STEPHN*((ETTP(NPT,JHP)+TKELVI)**4) )
C
            AUX4 = ETTPA(NPT,JTF) + AUX1*(AUX2+AUX3)
C
            AUX5 = DTP/AUX1
            AUX6 = EXP(-AUX5)
C
            TER1 = ETTPA(NPT,JHP) * AUX6
            TER2 = AUX4 * (1.D0-AUX6)
            TER3 = AUX4 * ( -AUX6+(1.D0-AUX6) / AUX5 )
C
            TSVAR(NPT,JHP) = 0.5D0 * TER1 + TER3
            ETTP(NPT,JHP) = TER1 + TER2
C
          ENDIF
        ENDDO
C
      ELSE IF (NOR.EQ.2) THEN
        DO NPT = 1,NBPART
          IF (ITEPA(NPT,JISOR).GT.0 .AND. IBORD(NPT).EQ.0) THEN
C
C
            ICHA  = ITEPA(NPT,JINCH)
            IEL   = ITEPA(NPT,JISOR)
C
            D2 = ETTP(NPT,JDP)*ETTP(NPT,JDP)
            DIAMP2 = XASHCH(ICHA)*TEPA(NPT,JRD0P)*TEPA(NPT,JRD0P)
     &             + (1.D0-XASHCH(ICHA))*TEPA(NPT,JRDCK)*TEPA(NPT,JRDCK)
C
            AUX1 = TEMPCT(NPT,1)*DIAMP2/D2
C
C    Combustion heterogene
            AUX2 = -GAMHET(NPT) *(ETTP(NPT,JMCK)**D2S3) *DELTAH(NPT)
C
C    Rayonnement
C
            AUX3 = PI * DIAMP2
     &           * ( PROPCE(IEL,IPPROC(ILUMI))
     &           - 4.D0*STEPHN*((ETTP(NPT,JHP)+TKELVI)**4) )
C
            AUX4 = ETTP(NPT,JTF) + AUX1*(AUX2+AUX3)
C
            AUX5 = DTP / AUX1
            AUX6 = EXP(-AUX5)
C
            TER1 = ETTPA(NPT,JHP) * AUX6
            TER2 = AUX4 * ( 1.D0-((1.D0-AUX6)/AUX5) )
C
            ETTP(NPT,JHP) = TSVAR(NPT,JHP) + 0.5D0*TER1 + TER2
C
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 12. Mise a jour du diametre du coeur retrecissant
C=======================================================================
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
          ICHA = ITEPA(NPT,JINCH)
C
          TEPA(NPT,JRDCK) = ( (D6SPI/(1.D0-XASHCH(ICHA)) )
     &                    *   ( ETTP(NPT,JMCH)/RHO0CH(ICHA)
     &                        + ETTP(NPT,JMCK)/RHOCK(ICHA) ) )**D1S3
        ENDIF
      ENDDO
C
C=======================================================================
C 13. Calcul du diametre des grains de charbon
C=======================================================================
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
          ICHA = ITEPA(NPT,JINCH)
          ETTP(NPT,JDP) = (XASHCH(ICHA)*(TEPA(NPT,JRD0P)**2)
     &                         + (1.D0-XASHCH(ICHA))
     &                          *(TEPA(NPT,JRDCK)**2) )**0.5D0
        ENDIF
      ENDDO
C
C=======================================================================
C 14. Calcul de la masse des grains de charbon
C=======================================================================
C
      DO NPT = 1,NBPART
        IF (ITEPA(NPT,JISOR).GT.0) THEN
          ICHA = ITEPA(NPT,JINCH)
          MP0  = DPIS6 * (TEPA(NPT,JRD0P)**3) * RHO0CH(ICHA)
          ETTP(NPT,JMP) = ETTP(NPT,JMCH) + ETTP(NPT,JMCK)
     &                      + XASHCH(ICHA)*MP0
        ENDIF
      ENDDO
C
C***********************************************************************
C
C=======
C FORMAT
C=======
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON        ',/,
     &'@      EST ACTIVE (LAGICH), ALORS QU''AUCUNE PHYSIQUE        ',/,
     &'@      PARTICULIERE SUR LA COMBUSTION DU CHABON PULVERISE    ',/,
     &'@      N''EST PAS ENCLENCHE (USPPMO).                        ',/,
     &'@                                                            ',/,
     &'@       IPHYLA = ', I10                                       ,/,
     &'@       IPPMOD(ICPL3C) = ', I10                               ,/,
     &'@       IPPMOD(ICP3PL) = ', I10                               ,/,
     &'@                                                            ',/,
     &'@  Le transport Lagrangien de particule de charbon doit      ',/,
     &'@   etre couple avec la combustion d''une flamme de charbon  ',/,
     &'@   pulverise en phase continue.                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
