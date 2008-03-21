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
                        SUBROUTINE LAGCAR
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
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &   VOLUME , DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   ,
     &   PIIL   , BX     , TEMPCT , STATIS ,
     &   GRADPR , GRADVF , ENERGI , DISSIP , ROMP   ,
     &   RDEVEL , RTUSER , RA        )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     CALCUL DES CARACTERISTIQUES DES PARTICULES : Tp, TL et PI
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
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
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
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! TAUP(NBPMAX) ! TR ! <-  ! TEMPS CARACTERISTIQUES DYNAMIQUE     !
CARGU ! TLAG(NBPMAX) ! TR ! <-  ! TEMPS CARACTERISTIQUES FLUIDE        !
CARGU ! PIIL(NBPMAX,3! TR ! <-  ! TERME DANS L'INTEGRATION DES EDS Up  !
CARGU ! BX(NBPMAX,3,2! TR ! <-  ! CARACTERISTIQUES DE LA TURBULENCE    !
CARGU ! TEMPCT       ! TR ! <-  ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU !   (NBPMAX,2) !    !     !                                      !
CARGU ! STATIS(NCELET! TR !  -> ! CUMUL DES STATISTIQUES VOLUMIQUES    !
CARGU !    NVLSTA)   !    !     !                                      !
CARGU ! GRADPR       ! TR !  -> ! GRADIENT DE PRESSION                 !
CARGU ! (NCELET,3)   !    !     !                                      !
CARGU ! GRADVF       ! TR !  -> ! GRADIENT DE LA VITESSE DU FLUIDE     !
CARGU ! (NCELET,3)   !    !     !                                      !
CARGU ! ENERGI(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! DISSIP(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! ROMP(NBPMAX) ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstphy.h"
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
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3)
      DOUBLE PRECISION PIIL(NBPMAX,3) , BX(NBPMAX,3,2)
      DOUBLE PRECISION TEMPCT(NBPMAX,2)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION GRADPR(NCELET,3) , GRADVF(NCELET,9)
      DOUBLE PRECISION ENERGI(NCELET) , DISSIP(NCELET), ROMP(NBPMAX)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL , IP , ID , IGVX , IGVY , IGVZ , IVT
      INTEGER          IROMF , IPHAS
C
      DOUBLE PRECISION CD1 , CD2 , REC , CL , C0 , CB , CBCB
      DOUBLE PRECISION UPART , VPART , WPART
      DOUBLE PRECISION UFLUI , VFLUI , WFLUI
      DOUBLE PRECISION UVWDIF , TL , UVWR
      DOUBLE PRECISION REP , D2 , D3 , FDR , D1S3 , D3S444 , D6SPI
      DOUBLE PRECISION BB1 , BB2 , BB3 , KTIL , BX1 , BX2 , BX3
      DOUBLE PRECISION VPMX , VPMY , VPMZ
      DOUBLE PRECISION R11 , R22 , R33
      DOUBLE PRECISION XNUL , ROM , PRT , FNUS , XRKL , XCP
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
      CD1  = 0.15D0
      CD2  = 0.687D0
      REC  = 1000.D0
      C0   = 2.1D0
      CL   = 1.D0 / (0.5D0 + (3.D0/4.D0)*C0 )
      CB   = 0.8D0
      CBCB = 0.64D0
C
      D6SPI = 6.D0 / PI
      D1S3 = 1.D0 / 3.D0
      D3S444 = 0.44D0 * 3.D0 / 4.D0
C
C Pointeur sur la masse volumique en fonction de l'ecoulement
C
      IF ( IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICFUEL).GE.0 ) THEN
        IROMF = IPPROC(IROM1)
      ELSE
        IROMF = IPPROC(IROM(IPHAS))
      ENDIF
C
C Calcul de la masse volumique
C
      DO IP = 1,NBPART
        IF ( ITEPA(IP,JISOR).GT.0 ) THEN
          D3 = ETTP(IP,JDP) * ETTP(IP,JDP) * ETTP(IP,JDP)
          ROMP(IP) = ETTP(IP,JMP) * D6SPI / D3
        ENDIF
      ENDDO
C
C=======================================================================
C 2. CALCUL DE Tp ET DE Tc SI THERMIQUE
C=======================================================================
C
      DO IP = 1,NBPART
C
        IF ( ITEPA(IP,JISOR) .GT.0 ) THEN
C
          IEL = ITEPA(IP,JISOR)
C
          ROM  = PROPCE(IEL,IROMF)
          XNUL = PROPCE(IEL,IPPROC(IVISCL(IPHAS))) / ROM
C
          UVWR = SQRT( ( ETTP(IP,JUF) -ETTP(IP,JUP) )*
     &                 ( ETTP(IP,JUF) -ETTP(IP,JUP) )
     &               + ( ETTP(IP,JVF) -ETTP(IP,JVP) )*
     &                 ( ETTP(IP,JVF) -ETTP(IP,JVP) )
     &               + ( ETTP(IP,JWF) -ETTP(IP,JWP) )*
     &                 ( ETTP(IP,JWF) -ETTP(IP,JWP) )  )
C
C--->  CALCUL DU REYNOLDS LOCAL
C
          REP  = UVWR * ETTP(IP,JDP) / XNUL
C
C--->  CALCUL DU COEFFICIENT DE TRAINEE
C
          D2 = ETTP(IP,JDP) * ETTP(IP,JDP)
C
          IF (REP.LE.REC) THEN
            FDR = 18.D0 * XNUL * (1.D0 + CD1 * REP**CD2) / D2
          ELSE
            FDR = D3S444 * UVWR / ETTP(IP,JDP)
          ENDIF
C
C--->  CALCUL DE Tp
C
          TAUP(IP) = ROMP(IP) / ROM / FDR
C
C--->  CALCUL UTILISATEUR DE Tp
C
          CALL USLATP
C         ===========
     &     ( IDEBIA , IDEBRA ,
     &       NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &       NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &       NVAR   , NSCAL  , NPHAS  ,
     &       NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &       NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &       IP     , ITEPA  , IDEVEL , ITUSER , IA     ,
     &       REP    , UVWR   , ROM    , ROMP(IP) , XNUL , TAUP(IP) ,
     &       XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &       VOLUME , DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &       ETTP   , ETTPA  , TEPA   ,
     &       RDEVEL , RTUSER , RA     )
C
C--->  CALCUL DE Tc
C
          IF ( (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) .OR.
     &         (IPHYLA.EQ.2)                        ) THEN
C
C     CP fluide
C
            IF (ICP(IPHAS).GT.0) THEN
              XCP = PROPCE( 1,IPPROC(ICP(IPHAS)) )
            ELSE
              XCP = CP0(IPHAS)
            ENDIF
C
C     CALCUL DU NUSSELT LOCAL
C
            IF ( IPPMOD(ICP3PL).GE.0 .OR.
     &           IPPMOD(ICPL3C).GE.0 .OR.
     &           IPPMOD(ICOD3P).GE.1 .OR.
     &           IPPMOD(ICOEBU).EQ.1 .OR.
     &           IPPMOD(ICOEBU).EQ.3 .OR.
     &           IPPMOD(ICFUEL).GE.0 .OR.
     &           IPPMOD(IELARC).GE.0 .OR.
     &           IPPMOD(IELJOU).GE.0      ) THEN
              IVT = IHM
            ELSE
              IVT = ISCALT(IPHAS)
            ENDIF
C
C a priori en combustion gaz ou CP, la diffusvite est toujours constante
C
            IF (IPPMOD(ICOEBU).EQ.0 .OR. IPPMOD(ICOEBU).EQ.2) THEN
              XRKL = DIFTL0 / ROM
            ELSE IF (IVISLS(IVT).GE.1) THEN
              XRKL = PROPCE(IEL,IPPROC(IVISLS(IVT))) / ROM
            ELSE
              XRKL = VISLS0(IVT) / ROM
            ENDIF
C
            PRT  = XNUL / XRKL
            FNUS = 2.D0 + 0.55D0 * REP**0.5D0 * PRT**(D1S3)
C
C Calcul du temps caracteristique thermique Tc
C
            TEMPCT(IP,1) = D2 * ROMP(IP) * ETTP(IP,JCP)
     &                   / ( FNUS * 6.D0 * ROM * XCP * XRKL )
C
C--->  CALCUL UTILISATEUR DE Tc
C
          CALL USLATC
C         ===========
     &     ( IDEBIA , IDEBRA ,
     &       NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &       NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &       NVAR   , NSCAL  , NPHAS  ,
     &       NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &       NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &       IP     , ITEPA  , IDEVEL , ITUSER , IA     ,
     &       REP    , UVWR   , ROM    , ROMP(IP) , XNUL ,
     &       XCP    , XRKL   , TEMPCT(IP,1) ,
     &       XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &       VOLUME , DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &       ETTP   , ETTPA  , TEPA   ,
     &       RDEVEL , RTUSER , RA     )
C
C Terme source implicite pour le couplage retour thermique
C
            TEMPCT(IP,2) = FNUS * PI * ETTP(IP,JDP) * XRKL * ROM
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C
C=======================================================================
C 3. CALCUL DE TL
C=======================================================================
C
C--> Calcul de l'energie turbulente et de la dissipation
C      en fonction du modele de turbulence
C
      IF (IDISTU.EQ.1) THEN
C
        IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50) THEN
          DO IEL = 1,NCEL
            ENERGI(IEL) = RTP(IEL,IK(IPHAS))
            DISSIP(IEL) = RTP(IEL,IEP(IPHAS))
          ENDDO
        ELSE IF (ITYTUR(IPHAS).EQ.3) THEN
          DO IEL = 1,NCEL
            ENERGI(IEL) = 0.5D0*( RTP(IEL,IR11(IPHAS))
     &                           +RTP(IEL,IR22(IPHAS))
     &                           +RTP(IEL,IR33(IPHAS)) )
            DISSIP(IEL) = RTP(IEL,IEP(IPHAS))
          ENDDO
        ELSE IF (ITURB(IPHAS).EQ.60) THEN
          DO IEL = 1,NCEL
            ENERGI(IEL) = RTP(IEL,IK(IPHAS))
            DISSIP(IEL) = CMU*ENERGI(IEL)*RTP(IEL,IOMG(IPHAS))
          ENDDO
        ELSE
          WRITE(NFECRA,2000) IILAGR, IDISTU, IPHAS, ITURB(IPHAS)
          CALL CSEXIT (1)
C              ======
        ENDIF
C
C--> Calcul de TL et BX
C
        DO IP = 1,NBPART
C
          IF (ITEPA(IP,JISOR).GT.0) THEN
C
            IEL = ITEPA(IP,JISOR)
C
C
            IF (DISSIP(IEL).GT.0.D0 .AND.
     &          ENERGI(IEL).GT.0.D0              ) THEN
C
            TL = CL * ENERGI(IEL) / DISSIP(IEL)
            TL = MAX(TL,EPZERO)
C
            UPART = ETTP(IP,JUP)
            VPART = ETTP(IP,JVP)
            WPART = ETTP(IP,JWP)
            UFLUI = ETTP(IP,JUF)
            VFLUI = ETTP(IP,JVF)
            WFLUI = ETTP(IP,JWF)
C
            IF (MODCPL.GT.0 .AND. IPLAS.GT.MODCPL) THEN
              IF (STATIS(IEL,ILPD).GT.SEUIL) THEN
                UPART = STATIS(IEL,ILVX) / STATIS(IEL,ILPD)
                VPART = STATIS(IEL,ILVY) / STATIS(IEL,ILPD)
                WPART = STATIS(IEL,ILVZ) / STATIS(IEL,ILPD)
                UFLUI = RTP(IEL,IU(IPHAS))
                VFLUI = RTP(IEL,IV(IPHAS))
                WFLUI = RTP(IEL,IW(IPHAS))
              ENDIF
            ENDIF
C
            UVWDIF = (UFLUI-UPART) * (UFLUI-UPART)
     &             + (VFLUI-VPART) * (VFLUI-VPART)
     &             + (WFLUI-WPART) * (WFLUI-WPART)
C
            UVWDIF = (3.D0 * UVWDIF) / (2.D0 *ENERGI(IEL))
C
            IF (MODCPL.GT.0 .AND. IPLAS.GT.MODCPL) THEN
C
              IF (IDIRLA.EQ.1) THEN
                BB1 = SQRT( 1.D0 + CBCB *UVWDIF )
                TLAG(IP,1)= TL / BB1
                BB2  = SQRT( 1.D0 + 4.D0 *CBCB *UVWDIF )
                TLAG(IP,2)= TL / BB2
                BB3  = SQRT( 1.D0 + 4.D0 *CBCB *UVWDIF )
                TLAG(IP,3)= TL / BB3
C
              ELSE IF (IDIRLA.EQ.2) THEN
                BB1  = SQRT( 1.D0 + 4.D0 *CBCB *UVWDIF )
                TLAG(IP,1)= TL / BB1
                BB2  = SQRT( 1.D0 + CBCB *UVWDIF )
                TLAG(IP,2)= TL / BB2
                BB3  = SQRT( 1.D0 + 4.D0 *CBCB *UVWDIF )
                TLAG(IP,3)= TL / BB3
C
              ELSE IF (IDIRLA.EQ.3) THEN
                BB1  = SQRT( 1.D0 + 4.D0 *CBCB *UVWDIF )
                TLAG(IP,1)= TL / BB1
                BB2  = SQRT( 1.D0 + 4.D0 *CBCB *UVWDIF )
                TLAG(IP,2)= TL / BB2
                BB3  = SQRT( 1.D0 + CBCB *UVWDIF )
                TLAG(IP,3)= TL / BB3
              ELSE
                WRITE(NFECRA,1000) IDIRLA
                CALL CSEXIT(1)
C                    ======
              ENDIF
C
              IF (ITYTUR(IPHAS).EQ.3) THEN
                R11 = RTP(IEL,IR11(IPHAS))
                R22 = RTP(IEL,IR22(IPHAS))
                R33 = RTP(IEL,IR33(IPHAS))
                KTIL = 3.D0 * ( R11*BB1 + R22*BB2 + R33*BB3  )
     &               / (2.D0 * (BB1+BB2+BB3) )
              ELSE IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &                                    .OR. ITURB(IPHAS).EQ.60) THEN
                KTIL = ENERGI(IEL)
              ENDIF
C
              BX1 = DISSIP(IEL) * ( (C0*BB1*KTIL/ENERGI(IEL))
     &                 +(2.D0 *(BB1*KTIL/ENERGI(IEL) -1.D0)/3.D0) )
              BX2 = DISSIP(IEL) * ( (C0*BB2*KTIL/ENERGI(IEL))
     &                 +(2.D0 *(BB2*KTIL/ENERGI(IEL) -1.D0)/3.D0) )
              BX3 = DISSIP(IEL) * ( (C0*BB3*KTIL/ENERGI(IEL))
     &                 +(2.D0 *(BB3*KTIL/ENERGI(IEL) -1.D0)/3.D0) )
C
              IF (BX1.GT.0.D0) THEN
                BX(IP,1,NOR) = SQRT(BX1)
              ELSE
                BX(IP,1,NOR) = 0.D0
              ENDIF
C
              IF (BX2.GT.0.D0) THEN
                BX(IP,2,NOR) = SQRT(BX2)
              ELSE
                BX(IP,2,NOR) = 0.D0
              ENDIF
C
              IF (BX3.GT.0.D0) THEN
                BX(IP,3,NOR) = SQRT(BX3)
              ELSE
                BX(IP,3,NOR) = 0.D0
              ENDIF
C
            ELSE
C
              TLAG(IP,1) = TL
              TLAG(IP,2) = TL
              TLAG(IP,3) = TL
C
              IF (IDIFFL.EQ.0) THEN
                UVWDIF = SQRT(UVWDIF)
                TLAG(IP,1) = TL/(1.D0 + CB*UVWDIF)
                TLAG(IP,2) = TLAG(IP,1)
                TLAG(IP,3) = TLAG(IP,1)
              ENDIF
C
              BX(IP,1,NOR) = SQRT(C0*DISSIP(IEL))
              BX(IP,2,NOR) = BX(IP,1,NOR)
              BX(IP,3,NOR) = BX(IP,1,NOR)
            ENDIF
C
            ELSE
C
              TLAG(IP,1) = EPZERO
              TLAG(IP,2) = EPZERO
              TLAG(IP,3) = EPZERO
              BX(IP,1,NOR) = ZERO
              BX(IP,2,NOR) = ZERO
              BX(IP,3,NOR) = ZERO
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ELSE
C
        DO IP = 1,NBPART
C
          IF ( ITEPA(IP,JISOR) .GT.0 ) THEN
            TLAG(IP,1) = EPZERO
            TLAG(IP,2) = EPZERO
            TLAG(IP,3) = EPZERO
            BX(IP,1,NOR) = ZERO
            BX(IP,2,NOR) = ZERO
            BX(IP,3,NOR) = ZERO
          ENDIF
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4. CALCUL DE PII
C=======================================================================
C
      DO ID = 1,3
C
        IGVX = 3*(ID-1)+1
        IGVY = 3*(ID-1)+2
        IGVZ = 3*(ID-1)+3
C
        DO IP = 1,NBPART
C
          IF ( ITEPA(IP,JISOR).GT.0 ) THEN
C
C--->   Calcul de II = ( -grad(P)/Rom(f)+grad(<Vf>)*(<Up>-<Uf>) )
C       ou
C       Calcul de II = ( -grad(P)/Rom(f) )
C
            IEL = ITEPA(IP,JISOR)
C
            PIIL(IP,ID) = GRADPR(IEL,ID)
C
            IF (MODCPL.GT.0 .AND. IPLAS.GT.MODCPL) THEN
              IF ( STATIS(IEL,ILPD) .GT. SEUIL ) THEN
                VPMX = STATIS(IEL,ILVX) / STATIS(IEL,ILPD)
                VPMY = STATIS(IEL,ILVY) / STATIS(IEL,ILPD)
                VPMZ = STATIS(IEL,ILVZ) / STATIS(IEL,ILPD)
C
                UFLUI = RTP(IEL,IU(IPHAS))
                VFLUI = RTP(IEL,IV(IPHAS))
                WFLUI = RTP(IEL,IW(IPHAS))
C
                PIIL(IP,ID) = GRADPR(IEL,ID)
     &                       +GRADVF(IEL,IGVX) * (VPMX-UFLUI)
     &                       +GRADVF(IEL,IGVY) * (VPMY-VFLUI)
     &                       +GRADVF(IEL,IGVZ) * (VPMZ-WFLUI)
C
              ENDIF
            ENDIF
C
C--->  Terme purement explicite : probleme avec petit diametre
C      NE PAS EFFACER SVP !
C
c            IF (IILAGR.EQ.2 .AND. ISTALA.GE.1) THEN
C
c              IF (STATIS(IEL,ILPD).GT.SEUIL) THEN
C
c                ROM   = PROPCE(IEL,IROMF)
C
c                FF = ROMP(IP) / ROM
c     &             *( STATIS(IEL,ILFV) / (DBLE(NPST)*VOLUME(IEL)) )
c     &             *( ETTP(IP,JUF+(ID-1)) - ETTP(IP,JUP+(ID-1)))
c     &                   /TAUP(IP)
C
c                PIIL(IP,ID) = PIIL(IP,ID) - FF
C
c              ENDIF
C
c            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C**********************************************************************
C
C--------
C FORMATS
C--------
C
 1000   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE CHOIX DE LA DIRECTION DU MODELE COMPLET              ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGCAR).                   ',/,
     &'@                                                            ',/,
     &'@    IDIRLA DEVRAIT ETRE UN ENTIER EGAL A 1 2 OU 3           ',/,
     &'@       (LA VALEUR 1 POUR UN ECOULEMENT SELON L''AXE X,      ',/,
     &'@        LA VALEUR 2 POUR UN ECOULEMENT SELON L''AXE Y,      ',/,
     &'@        LA VALEUR 3 POUR UN ECOULEMENT SELON L''AXE Z)      ',/,
     &'@       IL VAUT ICI IDIRLA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDIRLA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,
     &'@    DE TURBULENCE SELECTIONNE.                              ',/,
     &'@                                                            ',/,
     &'@   Le module Lagrangien a ete active avec IILAGR = ',I10     ,/,
     &'@     et la dispersion turbulente est prise en compte        ',/,
     &'@                                     avec IDISTU = ',I10     ,/,
     &'@   Le modele de turbulence active pour la phase ',I6         ,/,
     &'@     correspond a ITURB  = ',I10                             ,/,
     &'@   Or, les seuls traitements de la turbulence compatibles   ',/,
     &'@     avec le module Lagrangien et la dispersion turbulente  ',/,
     &'@     sont k-epsilon et Rij-epsilon, v2f et k-omega.         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,
     &'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,
     &'@  USINI1.                                                   ',/,
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
