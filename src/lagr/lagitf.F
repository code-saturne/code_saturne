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
                        SUBROUTINE LAGITF
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
     &   AUXL1  , AUXL2  , TEMPF  ,
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
CFONC      INTEGRATION DES EDS POUR LA TEMPERATURE FLUIDE
CFONC       VU PAR LES PARTICULES.
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
CARGU ! AUXL1(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! AUXL2(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      DOUBLE PRECISION AUXL1(NBPMAX) , AUXL2(NBPMAX)
      DOUBLE PRECISION TEMPF(NCELET)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          NPT   , IEL   , IPHAS  , MODE
      DOUBLE PRECISION CT    , AUX1  , AUX2   , TER1   , TER2
      DOUBLE PRECISION ENERG , DISSIP
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      IPHAS = ILPHAS
C
      CT = 1.D0
C
      MODE = 1
C
C=======================================================================
C 2. Temperature moyenne Fluide en degres Celsius
C=======================================================================
C
      IF ( IPPMOD(ICP3PL).GE.0 .OR.
     &     IPPMOD(ICPL3C).GE.0 .OR.
     &     IPPMOD(ICFUEL).GE.0      ) THEN
C
        DO IEL = 1,NCEL
          TEMPF(IEL) = PROPCE(IEL,IPPROC(ITEMP1)) - TKELVI
        ENDDO
C
      ELSE IF ( IPPMOD(ICOD3P).GE.0 .OR.
     &          IPPMOD(ICOEBU).GE.0 .OR.
     &          IPPMOD(IELARC).GE.0 .OR.
     &          IPPMOD(IELJOU).GE.0      ) THEN
C
        DO IEL = 1,NCEL
          TEMPF(IEL) = PROPCE(IEL,IPPROC(ITEMP)) - TKELVI
        ENDDO
C
      ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.-1 ) THEN
        DO IEL = 1,NCEL
          TEMPF(IEL) = RTP(IEL,ISCA(ISCALT(IPHAS)))
        ENDDO
C
      ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.1 ) THEN
        DO IEL = 1,NCEL
          TEMPF(IEL) = RTP(IEL,ISCA(ISCALT(IPHAS))) - TKELVI
        ENDDO
C
      ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.2 ) THEN
        DO IEL = 1,NCEL
          CALL USTHHT (MODE, RTP(IEL,ISCA(ISCALT(IPHAS))), TEMPF(IEL))
C         ===========
        ENDDO
      ENDIF
C
C=======================================================================
C 3. INTEGRATION DE L'EDS SUR LES PARTICULES
C=======================================================================
C
      DO NPT = 1,NBPART
C
        IF ( ITEPA(NPT,JISOR).GT.0 ) THEN
C
          IEL = ITEPA(NPT,JISOR)
C
          IF (ITYTUR(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.3 .OR.
     &         ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60 ) THEN
C
            IF ( ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50 ) THEN
C
              ENERG  = RTP(IEL,IK(IPHAS))
              DISSIP = RTP(IEL,IEP(IPHAS))
C
            ELSE IF ( ITYTUR(IPHAS).EQ.3 ) THEN
C
              ENERG  = 0.5D0 * ( RTP(IEL,IR11(IPHAS))
     &                         + RTP(IEL,IR22(IPHAS))
     &                         + RTP(IEL,IR33(IPHAS)) )
              DISSIP = RTP(IEL,IEP(IPHAS))
C
            ELSE IF (ITURB(IPHAS).EQ.60) THEN
C
              ENERG  = RTP(IEL,IK(IPHAS))
              DISSIP = CMU*RTP(IEL,IK(IPHAS))*RTP(IEL,IOMG(IPHAS))
C
            ENDIF
C
            AUXL1(NPT) = ENERG / ( CT*DISSIP )
            AUXL1(NPT) = MAX( AUXL1(NPT),EPZERO )
C
          ELSE
C
            AUXL1(NPT) = EPZERO
C
          ENDIF
C
        ENDIF
      ENDDO
C
      IF (NOR.EQ.1) THEN

        DO NPT = 1,NBPART
C
          IF (ITEPA(NPT,JISOR).GT.0) THEN
C
            IEL = ITEPA(NPT,JISOR)
C
            AUX1 = -DTP/AUXL1(NPT)
            AUX2 = EXP(AUX1)
C
            TER1 = ETTPA(NPT,JTF) * AUX2
            TER2 = TEMPF(IEL) * ( 1.D0-AUX2 )
C
            ETTP(NPT,JTF) = TER1 + TER2
C
C            Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2
C
            TSVAR(NPT,JTF) = 0.5D0 * TER1
     &                     + TEMPF(IEL) * ( -AUX2 +(AUX2-1.D0) / AUX1 )
          ENDIF
        ENDDO
C
      ELSE IF (NOR.EQ.2) THEN
c
        DO NPT = 1,NBPART
C
          IF ( ITEPA(NPT,JISOR).GT.0 .AND. IBORD(NPT).EQ.0 ) THEN
C
            IEL = ITEPA(NPT,JISOR)
C
            AUX1 = -DTP/AUXL1(NPT)
            AUX2 = EXP(AUX1)
C
            TER1 = 0.5D0 * ETTPA(NPT,JTF) * AUX2
            TER2 = TEMPF(IEL) * (1.D0 - (AUX2-1.D0) / AUX1)
C
            ETTP(NPT,JTF) = TSVAR(NPT,JTF) + TER1 + TER2
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C
C=======
C FORMAT
C=======
C
C----
C FIN
C----
C
      END
c@z
