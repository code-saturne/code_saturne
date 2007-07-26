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
                        SUBROUTINE CPPHY2
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   RTP    , PROPCE )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
CFONC  (CLASSES DE PARTICULES)
CFONC  VALEURS CELLULES
CFONC  ----------------
CFONC
CFONC    FRACTION MASSIQUE DE SOLIDE
CFONC      ET CLIPPING EVENTUELS
CFONC    DIAMETRE
CFONC    MASSE VOLUMIQUE
CFONC      ET CLIPPING EVENTUELS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
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
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
C
      DOUBLE PRECISION RTP(NCELET,*) , PROPCE(NCELET,*)
C
C VARIABLES LOCALES
C
      INTEGER          NBRINT
      PARAMETER       (NBRINT=8)
      INTEGER          IEL    , ICLA   , IPCRO2 , IPCDI2
      INTEGER          N1     , N2     , N3     , N4     , N5    , N6
      INTEGER          N7     , N8     , IPCX2C
      INTEGER          INTTMP(NBRINT)
C
      DOUBLE PRECISION XCH    , DCH    , XNP    , XCK    , DCK , D1S3
      DOUBLE PRECISION XASHCL , XUASH
      DOUBLE PRECISION X2MIN  , X2MAX  , DCKMIN , DCKMAX
      DOUBLE PRECISION DCHMIN , DCHMAX , ROMIN  , ROMAX , COEDMI
      DOUBLE PRECISION RO2INI , ROH2O
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      D1S3 = 1.D0/3.D0
C
C=======================================================================
C 2. CALCUL POUR CHAQUE CLASSE
C    DE LA FRACTION MASSIQUE DE SOLIDE
C    DU DIAMETRE DU COKE
C    DE LA MASSE VOLUMIQUE DU CHARBON
C=======================================================================
C
C --> Coefficient relatif au diametre de coke
C
      COEDMI = 1.2D0
C
      DO ICLA = 1, NCLACP
C
        N1 = 0
        N2 = 0
        N3 = 0
        N4 = 0
        N5 = 0
        N6 = 0
        N7 = 0
        N8 = 0
        X2MIN  =  GRAND
        X2MAX  = -GRAND
        DCHMIN =  GRAND
        DCHMAX = -GRAND
        DCKMIN =  GRAND
        DCKMAX = -GRAND
        ROMIN  =  GRAND
        ROMAX  = -GRAND
C
        DO IEL = 1, NCEL
C
          IPCX2C = IPPROC(IX2(ICLA))
          IPCRO2 = IPPROC(IROM2(ICLA))
          IPCDI2 = IPPROC(IDIAM2(ICLA))
          XCK    = RTP(IEL,ISCA(IXCK(ICLA)))
          XCH    = RTP(IEL,ISCA(IXCH(ICLA)))
          XNP    = RTP(IEL,ISCA(INP(ICLA)))
          XASHCL = XASHCH(ICHCOR(ICLA))
          XUASH  = XNP*XMP0(ICLA)*(1.D0-XASHCL)
C
C --- Calcul de la fraction massique de solide
C
          PROPCE(IEL,IPCX2C) = XCH + XCK + XNP*XMASH(ICLA)
C
C         Prise en compte de l'humidite
C
          IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
            PROPCE(IEL,IPCX2C) = PROPCE(IEL,IPCX2C)
     &                          +RTP(IEL,ISCA(IXWT(ICLA)))
          ENDIF

C
C ---- Clipping eventuels pour la fraction massique de solide
C
          IF ( PROPCE(IEL,IPCX2C) .GT. (1.D0+EPSICP) ) THEN
            N1 = N1 + 1
            X2MAX = MAX(PROPCE(IEL,IPCX2C),X2MAX)
            PROPCE(IEL,IPCX2C) = 1.D0
          ELSE IF ( PROPCE(IEL,IPCX2C) .LT. (ZERO-EPSICP) ) THEN
            N2 = N2 + 1
            X2MIN = MIN(PROPCE(IEL,IPCX2C),X2MIN)
            PROPCE(IEL,IPCX2C) = ZERO
          ENDIF
C
C
C --- Initialisation
C
          PROPCE(IEL,IPCRO2) = RHO20(ICLA)
          PROPCE(IEL,IPCDI2) = DIAM20(ICLA)
C
          IF ( XUASH.GT.EPSICP ) THEN
C
C --- Calcul du diametre du charbon reactif : Dch
C
            DCH = DIAM20(ICLA)*(XCH/XUASH)**D1S3
C
C ---- Clipping eventuels pour le diametre du charbon reactif
C
            IF ( DCH .GT. (DIAM20(ICLA)+EPSICP) ) THEN
              N3 = N3 + 1
              DCHMAX = MAX(DCH,DCHMAX)
              DCH = DIAM20(ICLA)
            ELSE IF ( DCH .LT. (ZERO-EPSICP) ) THEN
              N4 = N4 + 1
              DCHMIN = MIN(DCH,DCHMIN)
              DCH = ZERO
            ENDIF
C
C --- Calcul du diametre du coke : Dck stocke ds PROPCE(IEL,IPCDI2)
C
            DCK = ( (XCH/RHO20(ICLA)+XCK/RHOCK(ICHCOR(ICLA)))/
     &              ((1.D0-XASHCL)*PI/6.D0*XNP) )**D1S3
C
C ---- Clipping eventuels pour le diametre du coke
C
            IF ( DCK .GT. COEDMI*DIAM20(ICLA) ) THEN
              N5 = N5 + 1
              DCKMAX = MAX(DCK,DCKMAX)
              DCK = DIAM20(ICLA)*COEDMI
            ELSE IF ( DCK .LT. (ZERO-EPSICP) ) THEN
              N6 = N6 + 1
              DCKMIN = MIN(DCK,DCKMIN)
              DCK = ZERO
            ENDIF
            PROPCE(IEL,IPCDI2) = DCK
C
C --- Masse volumique
C
            RO2INI = RHO20(ICLA)
C         Prise en compte de l'humidite
            IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
C           pour l'instant on suppose que ROH2O est constant
              ROH2O = 998.203
              RO2INI = RHO20(ICLA)+ RTP(IEL,ISCA(IXWT(ICLA)))
     &                             *ROH2O
            ENDIF
C
            PROPCE(IEL,IPCRO2) =
     &        ( XASHCL*DIAM20(ICLA)**3*RHO20(ICLA) +
     &          (1.D0-XASHCL)*(DCK**3-DCH**3)*RHOCK(ICHCOR(ICLA)) +
     &          (1.D0-XASHCL)*DCH**3*RO2INI ) /
     &        ( XASHCL*DIAM20(ICLA)**3 +
     &          (1.D0-XASHCL)*DCK**3 )
C
C ---- Clipping pour la masse volumique
C
            IF ( PROPCE(IEL,IPCRO2) .GT. (RO2INI+EPSICP) ) THEN
              N7 = N7 + 1
              ROMAX = MAX(PROPCE(IEL,IPCRO2),ROMAX)
              PROPCE(IEL,IPCRO2) = RHO20(ICLA)
            ENDIF
            IF ( PROPCE(IEL,IPCRO2) .LT. (RHOCK(ICHCOR(ICLA))-EPSICP) )
     &                              THEN
              N8 = N8 + 1
              ROMIN = MIN(PROPCE(IEL,IPCRO2),ROMIN)
              PROPCE(IEL,IPCRO2) = RHOCK(ICHCOR(ICLA))
            ENDIF
          ENDIF
C
        ENDDO
C
        IF (IRANGP.GE.0) THEN
C
          INTTMP(1) = N1
          INTTMP(2) = N2
          INTTMP(3) = N3
          INTTMP(4) = N4
          INTTMP(5) = N5
          INTTMP(6) = N6
          INTTMP(7) = N7
          INTTMP(8) = N8
          CALL PARISM (NBRINT,INTTMP)
C         ===========
          N1 = INTTMP(1)
          N2 = INTTMP(2)
          N3 = INTTMP(3)
          N4 = INTTMP(4)
          N5 = INTTMP(5)
          N6 = INTTMP(6)
          N7 = INTTMP(7)
          N8 = INTTMP(8)
C
          CALL PARMAX (X2MAX )
C         ===========
          CALL PARMAX (DCHMAX)
C         ===========
          CALL PARMAX (DCKMAX)
C         ===========
          CALL PARMAX (ROMAX )
C         ===========
C
          CALL PARMIN (X2MIN )
C         ===========
          CALL PARMIN (DCHMIN)
C         ===========
          CALL PARMIN (DCKMIN)
C         ===========
          CALL PARMIN (ROMIN )
C         ===========
        ENDIF
C
        IF ( N1 .GT. 0 ) THEN
           WRITE(NFECRA,1001) ICLA, N1, X2MAX
        ENDIF
        IF ( N2 .GT. 0 ) THEN
           WRITE(NFECRA,1002) ICLA, N2, X2MIN
        ENDIF
        IF ( N3 .GT. 0 ) THEN
           WRITE(NFECRA,1003) ICLA, N3, DCHMAX
        ENDIF
        IF ( N4 .GT. 0 ) THEN
           WRITE(NFECRA,1004) ICLA, N4, DCHMIN
        ENDIF
        IF ( N5 .GT. 0 ) THEN
           WRITE(NFECRA,1005) ICLA, N5, DCKMAX
        ENDIF
        IF ( N6 .GT. 0 ) THEN
           WRITE(NFECRA,1006) ICLA, N6, DCKMIN
        ENDIF
        IF ( N7 .GT. 0 ) THEN
           WRITE(NFECRA,1007) ICLA, N7, ROMAX
        ENDIF
        IF ( N8 .GT. 0 ) THEN
           WRITE(NFECRA,1008) ICLA, N8, ROMIN
        ENDIF
C
      ENDDO
C
C----
C FORMATS
C----
C
 1001 FORMAT(/,1X,' CLIPPING EN MAX DE LA FRM SOL. POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Max       : ',G15.7)
 1002 FORMAT(/,1X,' CLIPPING EN MIN DE LA FRM SOL. POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Max       : ',G15.7)
 1003 FORMAT(/,1X,' CLIPPING EN MAX DU DIAMETRE CH POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Max       : ',G15.7)
 1004 FORMAT(/,1X,' CLIPPING EN MIN DU DIAMETRE CH POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Min       : ',G15.7)
 1005 FORMAT(/,1X,' CLIPPING EN MAX DU DIAMETRE CK POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Max       : ',G15.7)
 1006 FORMAT(/,1X,' CLIPPING EN MIN DU DIAMETRE CK POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Min       : ',G15.7)
 1007 FORMAT(/,1X,' CLIPPING EN MAX DE LA MASSE VOL. POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Max       : ',G15.7)
 1008 FORMAT(/,1X,' CLIPPING EN MIN DE LA MASSE VOL. POUR LA CLASSE ',
     &        I3,/,10X,' Nombre de points : ',I8,
     &           /,10X,' Valeur Min       : ',G15.7)
     &
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
