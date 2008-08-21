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
                        SUBROUTINE CAUSTA
C                       *****************
C     -------------------------------------------------------------
     & ( IFACE , IMP    , XKAPPA , CSTLOG , YPLULI ,
     &   APOW  , BPOW   , DPOW   ,
     &   UU    , DP     , XNU    , UET      )
C     -------------------------------------------------------------
C***********************************************************************
C     FONCTION  :
C     --------
c@foncb
CFONC
CFONC   CALCUL DE LA VITESSE DE FROTTEMENT UET PAR UNE METHODE DE NEWTON
CFONC    EN SUPPOSANT QUE UU VITESSE TANGENTIELLE REPOND A LA LOI
CFONC     UU/UET = 1/XKAPPA LOG (YPLUS) + CSTLOG
CFONC     YPLUS = UET*DP/XNU
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                         ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IFACE        ! E  !  -> ! NUMERO D'ENTITE TRAITEE              !
CARGU ! IMP          ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU !              !    !     ! > 0 EN CAS D'ERREUR FATALE           !
CARGU !              !    !     ! > 1 NON RESPECT D'UN CRITERE         !
CARGU !              !    !     ! > 2 MINI POUR BONNE ANALYSE          !
CARGU !              !    !     ! > 3 INTERVENTION INHABITUELLE DU CODE!
CARGU !              !    !     ! > 4 IMPRESSION POUR INFO             !
CARGU !              !    !     ! > 5 IDEM + COMPLEMENTS               !
CARGU ! XKAPPA       ! R  !  -> ! CST DE KARMAN                        !
CARGU ! CSTLOG       ! R  !  -> ! CST DE LA LOI LOG                    !
CARGU ! YPLULI       ! R  !  -> ! YPLUS LIMITE                         !
CARGU ! *POW         ! R  !  -> ! COEFFICIENTS DE LA LOI DE WERNER     !
CARGU ! UU           ! R  !  -> ! VITESSE TANGENTIELLE A LA PAROI      !
CARGU ! DP           ! R  !  -> ! DIST A LA PAROI DU PT OU EST PRIS UU !
CARGU ! XNU          ! R  !  -> ! VISCOSITE CINEMATIQUE MOLECULAIRE    !
CARGU ! UET          ! R  ! <-  ! VITESSE DE FROTTEMENT                !
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
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IFACE, IMP
      DOUBLE PRECISION XKAPPA, CSTLOG, YPLULI, APOW, BPOW, DPOW
      DOUBLE PRECISION UU, DP, XNU, UET
C
C VARIABLES LOCALES
C
      INTEGER          IPOS, NITM, NIT
      DOUBLE PRECISION EPS, RR, YP, UETWER, UETMIN, UETA, YSURNU
C
C***********************************************************************
C
C     -------------------
C     1 - INITIALISATIONS
C     -------------------
C
      EPS  = 0.001D0
      NITM = 100
      NIT  = 0
C
      IPOS = 0
      IF ( IMP .GT. 4 ) THEN
         WRITE ( NFECRA,1000 ) IFACE
         IPOS = 1
         WRITE ( NFECRA,1100 ) UU, DP
      ENDIF
C
C     --------------------------
C   2 - CALCUL DU REYNOLDS LOCAL
C     --------------------------
C
      YSURNU = DP / XNU
      RR = UU * YSURNU
C
C     ------------------------------------------
C     3 - CALCUL DE LA VITESSE DE FROTTEMENT UET
C     ------------------------------------------
C
C       SOUS-COUCHE VISQUEUSE
C
      IF ( RR .LE. YPLULI**2 ) THEN
        UET = SQRT( UU / YSURNU )
      ELSE
C
C       RESOLUTION ITERATIVE PAR UNE METHODE DE NEWTON
C
C     On initialise par Werner ou la borne inf de uet assurant la cv
C     de la methode
        UETWER = ( ABS(UU)/APOW/(YSURNU**BPOW) )**DPOW
        UETMIN = EXP(-CSTLOG*XKAPPA)/YSURNU
        UET = MAX(UETWER,UETMIN)
C
 300    CONTINUE
        NIT = NIT+1
C
        UETA = UET
        UET = (XKAPPA*UU + UET)/(LOG(YSURNU*UET)+XKAPPA*CSTLOG+1.D0)
C
        IF ( ABS ( UET-UETA ) .LE. EPS * UETA ) THEN
          IF ( IMP .GT. 5 ) WRITE ( NFECRA,3000 ) NIT, EPS
        ELSE IF ( NIT .GE. NITM ) THEN
          IF ( IMP .GT .1 ) THEN
            IF ( IPOS .LE. 0 ) WRITE ( NFECRA,1000 ) IFACE
            IPOS = 1
            WRITE ( NFECRA,3100 ) NITM, EPS
          ENDIF
        ELSE
          GO TO 300
        ENDIF
C
      ENDIF
C
C     ------------------------------------------------------
C     4 -  POSITION DU POINT INTERIEUR DANS LA COUCHE LIMITE
C     ------------------------------------------------------
C
      IF ( IMP .GT. 4 ) THEN
        YP = YSURNU*UET
        IF ( YP .LE. YPLULI ) THEN
          WRITE ( NFECRA,4000 ) YP
        ELSE
          WRITE ( NFECRA,4100 ) YP
        ENDIF
      ENDIF
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT ( 5X,'APPEL DU SOUS-PROGRAMME CAUSTA A LA FACE IFAC=',I9)
 1100 FORMAT ( 10X,'CONDITIONS D''ENTREE : UU=',E10.3,2X,'DP=',E10.3 )
 3000 FORMAT ( 10X,'CONVERGENCE DU CALCUL DE UET APRES NIT=',I4,
     & 1X,'ITERATIONS, PRECISION DEMANDEE EPS=',E10.3 )
 3100 FORMAT ( 10X,'NOMBRE MAXIMAL D''ITERATIONS ATTEINT POUR LE',
     & 1X,'CALCUL DE UET : NITM=',I4,', PRECISION',
     & 1X,'DEMANDEE EPS=',E10.3 )
 4000 FORMAT ( 10X,'LE PREMIER POINT INTERIEUR EST DANS LA SOUS-COUCHE',
     & 1X,'VISQUEUSE, Y+=',E10.3 )
 4100 FORMAT ( 10X,'LE PREMIER POINT INTERIEUR EST A LA DISTANCE DE',
     & 1X,'LA PAROI, Y+=',E10.3 )
C
#else
C
 1000 FORMAT ( 5X,'CAUSTA SUBROUTINE CALLED FOR FACE IFAC=',I9)
 1100 FORMAT ( 10X,'INPUT CONDITIONS: UU=',E10.3,2X,'DP=',E10.3 )
 3000 FORMAT ( 10X,'UET COMPUTATION CONVERGENCE AFTER NIT=',I4,
     & 1X,'ITERATIONS, DESIRED PRECISION EPS=',E10.3 )
 3100 FORMAT ( 10X,'MAXIMUM NUMBER OF ITERATIONS REACHED FOR THE',
     & 1X,'COMPUTATION OF UET: NITM=',I4,', DESIRED PRECISION',
     & 1X,'EPS=',E10.3 )
 4000 FORMAT ( 10X,'THE FIRST POINT IS IN THE VISCOUS SUBLAYER',
     & 1X,' Y+=',E10.3 )
 4100 FORMAT ( 10X,'THE FIRST POINT IS AT A WALL-DISTANCE OF',
     & 1X,' Y+=',E10.3 )
C
#endif
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
