c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USTMGR
C                       *****************
C
C     -------------------------------------------------------------
     &    ( IAPPEL , IGR    , ISYM   ,
     &      NCELF  , NCELFE , NFACF  , IWARNP ,
     &      IUSMGR , NIW    , NRW    ,
     &      IFACEF ,
     &      DAF    , XAF    , SURFAF , VOLUMF , XYZFIN ,
     &      IRSCEL ,
     &      IW     , RW     )
C     -------------------------------------------------------------
C*******************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC   MULTIGRILLE ALGEBRIQUE :
CFONC
CFONC   DETERMINATION DE LA CONNECTIVITE
CFONC   CELLULE FINE -> CELLULE GROSSIERE (IRSCEL)
CFONC   POUR LA CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER
CFONC   A PARTIR DU NIVEAU SUPERIEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IAPPEL       ! E  !  -> ! 1 : SELECTION/DIMENSIONNEMENT        !
CARGU !              !    !     ! 2 : DETERMINATION DE IRSCEL          !
CARGU ! IGR          ! E  !  -> ! NIVEAU DU MAILLAGE GROSSIER          !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! NCELF        ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE FIN       !
CARGU ! NCELFE       ! E  !  -> ! NOMBRE D'ELEMENTS ETENDUS FIN        !
CARGU ! NFACF        ! E  !  -> ! NOMBRE DE FACES INTERNES MAILL. FIN  !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! IUSMGR       ! E  ! <-  ! 0 : AGGLOMERATION AUTOMATIQUE        !
CARGU !              !    !     ! 1 : ON UTILISE CE SOUS-PROGRAMME     !
CARGU ! NIW / NRW    ! E  ! <-  ! TAILLES DES TABLEAUX IW / RW         !
CARGU !              !    !     ! POUR IAPPEL = 2                      !
CARGU !              !    !     ! (DETERMINEES LORSQUE IAPPEL = 1)     !
CARGU ! IFACEF       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFACF)   !    !     !  DU MAILLAGE FIN                     !
CARGU ! DAF(NCELFE)  ! TR !  -> ! DIAGONALE MATRICE MAILLAGE FIN       !
CARGU ! XAF          ! TR !  -> ! EXTRADIAGONALE MATRICE MAILLAGE FIN  !
CARGU ! (NFACF, ISYM)!    !     !                                      !
CARGU ! SURFAF       ! TR !  -> ! SURFACES FACES INTERNES MAILLAGE FIN !
CARGU ! (3, NFACF)   !    !     !                                      !
CARGU ! VOLUMF       ! TR !  -> ! VOLUMES DES CELLULES DU MAILLAGE FIN !
CARGU ! (NCELFE)     !    !     !                                      !
CARGU ! XYZFIN       ! TR !  -> ! CENTRES DES CELLULES DU MAILLAGE FIN !
CARGU ! (3, NCELFE)  !    !     !                                      !
CARGU ! IRSCEL       ! TE ! <-  ! CELLULE FINE -> CELLULE GROSSIERE    !
CARGU !  (NCELFE)    !    !     !                                      !
CARGU ! IW(NIW)      ! TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! RW(NRW)      ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IAPPEL, IGR, ISYM
      INTEGER          NCELF, NCELFE, NFACF
      INTEGER          IWARNP
      INTEGER          IUSMGR, NIW, NRW
C
      INTEGER          IFACEF(2, NFACF)
      INTEGER          IRSCEL(NCELFE)
C
      DOUBLE PRECISION DAF(NCELFE), XAF(NFACF, ISYM)
      DOUBLE PRECISION SURFAF(3, NFACF), VOLUMF(NCELFE)
      DOUBLE PRECISION XYZFIN(3, NCELFE)
      INTEGER          IW(NIW)
      DOUBLE PRECISION RW(NRW)
C
C
C VARIABLES LOCALES
C
      LOGICAL          CRITER
      INTEGER          NCELG
      INTEGER          IEL, IEL1, IEL2, IEL1V, IFAC, IFAC2
      INTEGER          IELG, IEL1G, IEL2G, IFACG, ICELG, NCELT
      INTEGER          INDIC, IRSPR2, IPAIMP, IW1
      INTEGER          ICOMP, IEDIR
C
      INTEGER          IMP, JMP, KMP, I, J, K, ICEL
      INTEGER          NTEST, IP, JP, KP, IJK, IGLOB, I1, IG, JG
      DOUBLE PRECISION DX, DY, DZ, XYZMIN(3), XYZMAX(3)
      DOUBLE PRECISION XMN, XMX, YMN, YMX, ZMN, ZMX
      DOUBLE PRECISION XG, YG, ZG, RAP, REP, EPSLON
      DOUBLE PRECISION NGROS
C
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF (1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
      NCELG = 0
      EPSLON = +1.D-6
C
C=======================================================================
C 0.  INITIALISATION ET DIMENSIONNEMENT
C=======================================================================
C
      IUSMGR = 0
C
C     DIMENSIONNEMENT MEMOIRE
C
      INDIC  = 1
      IRSPR2 = INDIC + NCELFE
      IPAIMP = IRSPR2 + NCELFE
      NIW    = IPAIMP + NCELFE
C
      IW1    = 1
      NRW    = IW1 + NCELFE
C
      IF (IAPPEL.EQ.1) RETURN
C
C=======================================================================
C 1.  CREATION DES CELLULES GROSSIERES : IRSCEL
C=======================================================================
C
C***********************************************************************
C  MAILLAGE DU CARRE CAS 2D
C
C     BORNE MIN/MAX DU MAILLAGE
C
      DO I = 1, 3
        XYZMIN(I) = +1.D12
        XYZMAX(I) = -1.D12
        DO ICEL = 1, NCELF
          XYZMIN(I) = MIN(XYZFIN(I, ICEL), XYZMIN(I))
          XYZMAX(I) = MAX(XYZFIN(I, ICEL), XYZMAX(I))
        ENDDO
        XYZMIN(I) = XYZMIN(I) -EPSLON
        XYZMAX(I) = XYZMAX(I) +EPSLON
      ENDDO
C
C     EVALUATION DE IMP,JMP,KMP AVEC LA REGLE MAILLAGE REGULIER
C      IMP/(Xmax-Xmin) ~ JMP/(Ymax-Ymin) ~ KMP/(Zmax-Zmin)
C                       NCELG ~ IMP.JMP.KMP ~ NCEL/NGROS
C
C
C     CAS 2D
C
      NGROS = 4
C
      RAP    = NCELF/NGROS
      DX = XYZMAX(1) -XYZMIN(1)
      DY = XYZMAX(2) -XYZMIN(2)
      REP = SQRT(RAP*DX/DY)
      KMP   = 1
      IMP = INT(REP)
      JMP = INT(REP*DY/DX)
      IMP = MAX(IMP, 1)
      JMP = MAX(JMP, 1)
      NCELG = IMP*JMP*KMP
C
C     INCLUSION DANS LE MAILLAGE STRUCTURÉ ORTHOGONAL
C
      DX = (XYZMAX(1) -XYZMIN(1))/IMP
      DY = (XYZMAX(2) -XYZMIN(2))/JMP
      DZ = (XYZMAX(3) -XYZMIN(3))/KMP
C
      DO ICEL = 1, NCELF
        IW(INDIC+ICEL-1) = 0
        IW(IRSPR2+ICEL-1) = 0
      ENDDO
C
      DO ICEL = 1, NCELF
C
C       INITIALISATION
        XG = XYZFIN(1, ICEL)
        YG = XYZFIN(2, ICEL)
        ZG = XYZFIN(3, ICEL)
        IP = 0
        JP = 0
        KP = 0
C
        DO I = 1, IMP
          XMN = XYZMIN(1) +(I-1)*DX
          XMX = XYZMIN(1) +   I *DX
          IF (IP.EQ.0) THEN
            IF (XG.GE.XMN .AND. XG.LE.XMX) THEN
              IP = I
            ENDIF
          ENDIF
        ENDDO
C
        DO J = 1, JMP
          YMN = XYZMIN(2) +(J-1)*DY
          YMX = XYZMIN(2) +   J *DY
          IF (JP.EQ.0) THEN
            IF (YG.GE.YMN .AND. YG.LE.YMX) THEN
              JP = J
            ENDIF
          ENDIF
        ENDDO
C
        DO K = 1, KMP
          ZMN = XYZMIN(3) +(K-1)*DZ
          ZMX = XYZMIN(3) +   K *DZ
          IF (KP.EQ.0) THEN
            IF (ZG.GE.ZMN .AND. ZG.LE.ZMX) THEN
              KP = K
            ENDIF
          ENDIF
        ENDDO
C
        IJK = IP*JP*KP
C
        IF (IJK.GT.0) THEN
          IGLOB = IMP*JMP*(KP-1) +IMP*(JP-1)+IP
          IW(IRSPR2+ICEL-1) = IGLOB
          IW(INDIC+IGLOB-1) = +1
        ELSE
          WRITE(NFECRA,*) ' pb de programmation dans ustmgr.F  '
          WRITE(NFECRA,*) ' arret dans ustmgr.F '
          WRITE(NFECRA,*) ' IP, JP, KP = ', IP, JP, KP
          CALL CSEXIT(1)
        ENDIF
C
      ENDDO
C
C
C     DEFINITION DU MAILLAGE PAIR-IMPAIR POUR VISUALISATION
C
      DO IGLOB =1,NCELF
        IW(IPAIMP+IGLOB-1) = 0
      ENDDO
C
      DO K = 1, KMP
        DO J = 1, JMP
          I1 = MOD(K+J,2) +1
          DO I = I1, IMP, 2
            IGLOB = IMP*JMP*(K-1) +IMP*(J-1) +I
            IW(IPAIMP+IGLOB-1) = -1
          ENDDO
          I1 = MOD(K+J+1,2)+1
          DO I = I1, IMP, 2
            IGLOB = IMP*JMP*(K-1) +IMP*(J-1) +I
            IW(IPAIMP+IGLOB-1) = +1
          ENDDO
        ENDDO
      ENDDO
C
C     Compactage du maillage structure si, au moins, une maille
C     grossiere ne contient aucune maille du maillage fin non structure
C
      NTEST = 0
      DO IGLOB = 1, NCELG
        IF (IW(INDIC+IGLOB-1).LE.0) NTEST=NTEST+1
      ENDDO
C
      IF (NTEST.LE.0) THEN
        DO ICEL=1,NCELF
          IRSCEL(ICEL) = IW(IRSPR2+ICEL-1)
        ENDDO
C
      ELSE
C
        I = 0
        DO IGLOB = 1, NCELG
          IF (IW(INDIC+IGLOB-1).GT.0) THEN
            I = I +1
            IW(INDIC+IGLOB-1) = I
          ENDIF
        ENDDO
C
C       NOUVEAU NCELG
        NCELG = I
C
        DO ICEL=1,NCELF
          IRSCEL(ICEL) = IW(INDIC+IW(IRSPR2+ICEL-1)-1)
        ENDDO
C
      ENDIF
C
      I = 0
      J = NCELF
      DO ICEL = 1, NCELF
        I = MIN(IRSCEL(ICEL), I)
        J = MAX(IRSCEL(ICEL), J)
      ENDDO
C
      IEDIR = 0
C     IEDIR = 1
C
      IF (IEDIR.GT.0) THEN
C
C       EXTRACTION DES POINTS DIRICHLETS (DIAGONALE FORTEMENT DOMINANTE)
C
        DO ICEL = 1, NCELF
          IW(INDIC+ICEL-1) = 0
          IW(IRSPR2+ICEL-1) = 0
          IW(IW1+ICEL-1) = DAF(ICEL)
        ENDDO
C
        DO IFAC = 1, NFACF
          I = IFACEF(1,IFAC)
          J = IFACEF(2,IFAC)
          IW(IW1+I-1) = IW(IW1+I-1) +XAF(IFAC,1)
          IW(IW1+J-1) = IW(IW1+J-1) +XAF(IFAC,ISYM)
          IG = IRSCEL(I)
          JG = IRSCEL(J)
          IW(INDIC+IG-1) = IW(INDIC+IG-1) +1
          IW(INDIC+JG-1) = IW(INDIC+JG-1) +1
        ENDDO
C
        REP = 0.1D0
        DO ICEL = 1, NCELF
          XMX = IW(IW1+ICEL-1)/DAF(ICEL)
          IG = IRSCEL(ICEL)
          IF (XMX.GT.REP) THEN
            IF (IW(INDIC+IG-1).GE.2) THEN
              NCELG = NCELG+1
              IRSCEL(ICEL) = NCELG
              IW(INDIC+IG-1) = IW(INDIC+IG-1) -1
              IW(INDIC+NCELG-1) = IW(INDIC+NCELG-1) -1
            ENDIF
          ENDIF
        ENDDO
C
      ENDIF
C
C
C**********************************************************************
C
C--------
C FORMATS
C--------
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
