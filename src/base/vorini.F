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
                        SUBROUTINE VORINI
C                       *****************
C     -------------------------------------------------------------
     & ( NCEVOR , NVOR   , IENT   ,
     &   IVORCE , XYZ    , YZCEL  , XU     , XV     , XW     ,
     &   YZVOR  , SIGNV  , TEMPS  , TPSLIM )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ----------
c@foncb
CFONC
CFONC INITIALISATION DE LA METHODE DES VORTEX POUR LES ENTREES EN L.E.S.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEVOR       ! E  !  -> ! NOMBRE DE FACE A L'ENTREE OU EST     !
CARGU !              !    !     ! UTILISE LA METHODE                   !
CARGU ! NVOR         ! E  !  -> ! NOMBRE DE VORTEX A L'ENTREE          !
CARGU ! IENT         ! E  !  -> ! NUMERO DE L'ENTREE                   !
CARGU ! IVORCE       ! TE ! <-> ! NUMERO DU VORTEX LE PLUS PROCHE D'UNE!
CARGU !     (NVOMAX) !    !     ! FACE DONNEE                          !
CARGU ! XYZ(ICVMAX,3)!    !  -> ! COORDONNEES DES FACES D'ENTREE DANS  !
CARGU !              !    !     ! LE CALCUL                            !
CARGU ! YZCEL        ! TR !  -> ! COORDONNEES DES FACES D'ENTREE DANS  !
CARGU !   (ICVMAX ,2)!    !     ! LE REFERENTIEL LOCAL                 !
CARGU ! XU(ICVMAX)   ! TR !  -  ! COMPOSANTE DE VITESSE PRINCIPALE     !
CARGU ! XV(ICVMAX)   ! TR ! <-> ! COMPOSANTES DE VITESSE TRANSVERSES   !
CARGU ! XW(ICVMAX)   ! TR ! <-> !                                      !
CARGU ! YZVOR        ! TR ! <-> ! COORDONNEES DU CENTRE DES VORTEX     !
CARGU !   (NVOMAX,2) !    !     !                                      !
CARGU ! SIGNV(NVOMAX)! TR ! <-> ! SENS DE ROTATION DES VORTEX          !
CARGU ! TEMPS        ! TR ! <-> ! TEMPS ECOULE DEPUIS LA CREATION      !
CARGU !     (NVOMAX) !    !     ! DU VORTEX                            !
CARGU ! TPSLIM       ! TR ! <-> ! DUREE DE VIE DU VORTEX               !
CARGU !     (NVOMAX) !    !     !                                      !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCEVOR , NVOR   , IENT
      INTEGER          IVORCE(NVOMAX)
C
      DOUBLE PRECISION XYZ(ICVMAX,3)   , YZCEL(ICVMAX,2)
      DOUBLE PRECISION XU(ICVMAX)      , XV(ICVMAX)      , XW(ICVMAX)
      DOUBLE PRECISION YZVOR(NVOMAX,2) , SIGNV(NVOMAX)
      DOUBLE PRECISION TEMPS(NVOMAX)   , TPSLIM(NVOMAX)
C
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, KK, III, IUN, IVORT, IIENT, IOK
C
      DOUBLE PRECISION DD
      DOUBLE PRECISION DRAND(1), PHIDAT, XX, YY, ZZ
      DOUBLE PRECISION UU, VV, WW
      DOUBLE PRECISION U_VOR, EK_VOR, EE_VOR
C
      INTEGER          ILECT
      DATA             ILECT /0/
      SAVE             ILECT
C=======================================================================
C 1. CALCUL DU REPERE LOCAL ET CHANGEMENT DE REPERE
C=======================================================================
C
      ILECT = ILECT + 1
      IF(ILECT.EQ.1) THEN
        WRITE(NFECRA,1000) NNENT, ISUIVO
      ENDIF
C
      IF(ICAS(IENT).EQ.1.OR.ICAS(IENT).EQ.2.OR.ICAS(IENT).EQ.3) THEN
C
        DIR3(1,IENT)=DIR1(2,IENT)*DIR2(3,IENT)-DIR1(3,IENT)*DIR2(2,IENT)
        DIR3(2,IENT)=DIR1(3,IENT)*DIR2(1,IENT)-DIR1(1,IENT)*DIR2(3,IENT)
        DIR3(3,IENT)=DIR1(1,IENT)*DIR2(2,IENT)-DIR1(2,IENT)*DIR2(1,IENT)
C
      ELSEIF(ICAS(IENT).EQ.4) THEN
C
C     On s'aide du vecteur surface d'une face de l'entree (supposee plane)
C     pour definir le repere local
C
        VV = SQRT(SURF(1,IENT)**2 + SURF(2,IENT)**2 + SURF(3,IENT)**2)
C
        DIR3(1,IENT) = -SURF(1,IENT)/VV
        DIR3(2,IENT) = -SURF(2,IENT)/VV
        DIR3(3,IENT) = -SURF(3,IENT)/VV
C
C     On se fixe, par exemple, x1 = 0 et y1 = 1 et on norme le vecteur
C
        DIR1(1,IENT) = 0.D0
        DIR1(2,IENT) = 1.D0
        IF (ABS(DIR3(3,IENT)).GT.EPZERO) THEN
          DIR1(3,IENT) = -DIR3(2,IENT)/DIR3(3,IENT)
        ELSE
          DIR1(3,IENT) = 0.D0
        ENDIF
C
        VV = SQRT(DIR1(1,IENT)**2 + DIR1(2,IENT)**2 + DIR1(3,IENT)**2)
C
        DIR1(1,IENT) = DIR1(1,IENT)/VV
        DIR1(2,IENT) = DIR1(2,IENT)/VV
        DIR1(3,IENT) = DIR1(3,IENT)/VV
C
C     On obtient le dernier vecteur par produit vectoriel des deux autres
C
        DIR2(1,IENT) =
     &    DIR3(2,IENT)*DIR1(3,IENT) - DIR3(3,IENT)*DIR1(2,IENT)
        DIR2(2,IENT) =
     &    DIR3(3,IENT)*DIR1(1,IENT) - DIR3(1,IENT)*DIR1(3,IENT)
        DIR2(3,IENT) =
     &    DIR3(1,IENT)*DIR1(2,IENT) - DIR3(2,IENT)*DIR1(1,IENT)
C
      ENDIF
C
C - Changement de repere (on suppose que les vecteurs sont normes)
C
      DO II = 1, NCEVOR
        XX = XYZ(II,1) - CEN(1,IENT)
        YY = XYZ(II,2) - CEN(2,IENT)
        ZZ = XYZ(II,3) - CEN(3,IENT)
        YZCEL(II,1) = DIR1(1,IENT)*XX+DIR1(2,IENT)*YY+DIR1(3,IENT)*ZZ
        YZCEL(II,2) = DIR2(1,IENT)*XX+DIR2(2,IENT)*YY+DIR2(3,IENT)*ZZ
      ENDDO
C
C - Dimensions min et max de l entree
C
      YMAX(IENT) = -GRAND
      YMIN(IENT) =  GRAND
      ZMAX(IENT) = -GRAND
      ZMIN(IENT) =  GRAND
      DO II = 1, NCEVOR
        YMAX(IENT) = MAX(YMAX(IENT),YZCEL(II,1))
        YMIN(IENT) = MIN(YMIN(IENT),YZCEL(II,1))
        ZMAX(IENT) = MAX(ZMAX(IENT),YZCEL(II,2))
        ZMIN(IENT) = MIN(ZMIN(IENT),YZCEL(II,2))
      ENDDO
C
C - Verification
C
      IF(ICAS(IENT).EQ.1) THEN
        IF(LLY(IENT).LT.YMAX(IENT)-YMIN(IENT).OR.
     &       LLZ(IENT).LT.ZMAX(IENT)-ZMIN(IENT)) THEN
          WRITE(NFECRA,2000) IENT
          CALL CSEXIT(1)
        ENDIF
      ELSEIF(ICAS(IENT).EQ.2) THEN
        IF(LLD(IENT).LT.YMAX(IENT)-YMIN(IENT).OR.
     &       LLD(IENT).LT.ZMAX(IENT)-ZMIN(IENT)) THEN
          WRITE(NFECRA,2000) IENT
          CALL CSEXIT(1)
        ENDIF
      ENDIF
C=======================================================================
C 2. IMPRESSIONS DES PARAMETES A CHAQUE ENTREE
C=======================================================================
C
      CALL VORIMP(IENT)
C     ===========
C=======================================================================
C 3. REMPLISSAGE DES TABLEAUX DE DONNEES EN ENTREE
C=======================================================================
C
      IUN = 1
C
      IF(ICAS(IENT).EQ.1.OR.ICAS(IENT).EQ.2.OR.ICAS(IENT).EQ.3) THEN
        OPEN(FILE=FICVOR(IENT),UNIT=IMPDVO)
        REWIND(IMPDVO)
        DO II = 1, NDAT(IENT)
          READ(IMPDVO,*)
     &         XDAT(II,IENT) , YDAT(II,IENT) , ZDAT(II,IENT)  ,
     &         UDAT(II,IENT) , VDAT(II,IENT) , WDAT(II,IENT)  ,
     &         DUDAT(II,IENT), KDAT(II,IENT) , EPSDAT(II,IENT)
        ENDDO
        CLOSE(IMPDVO)
        WRITE(NFECRA,3000)
      ELSEIF(ICAS(IENT).EQ.4) THEN
        XDAT(1,IENT)   = CEN(1,IENT)
        YDAT(1,IENT)   = CEN(2,IENT)
        ZDAT(1,IENT)   = CEN(3,IENT)
        UDAT(1,IENT)   = UDEBIT(IENT)
        VDAT(1,IENT)   = 0.D0
        WDAT(1,IENT)   = 0.D0
        DUDAT(1,IENT)  = 0.D0
        KDAT(1,IENT)   = KDEBIT(IENT)
        EPSDAT(1,IENT) = EDEBIT(IENT)
      ENDIF
C
C On suppose que les donnees sont fournies
C dans le repere du calcul (et non dans le repere local).
C C'est plus simple pour l'utilisateur, et plus
C naturel pour faire du couplage
C
      DO II = 1, NDAT(IENT)
        XX = XDAT(II,IENT) - CEN(1,IENT)
        YY = YDAT(II,IENT) - CEN(2,IENT)
        ZZ = ZDAT(II,IENT) - CEN(3,IENT)
        UU = UDAT(II,IENT)
        VV = VDAT(II,IENT)
        WW = WDAT(II,IENT)
        YDAT(II,IENT) = DIR1(1,IENT)*XX+DIR1(2,IENT)*YY+DIR1(3,IENT)*ZZ
        ZDAT(II,IENT) = DIR2(1,IENT)*XX+DIR2(2,IENT)*YY+DIR2(3,IENT)*ZZ
        UDAT(II,IENT) = DIR3(1,IENT)*UU+DIR3(2,IENT)*VV+DIR3(3,IENT)*WW
        VDAT(II,IENT) = DIR1(1,IENT)*UU+DIR1(2,IENT)*VV+DIR1(3,IENT)*WW
        WDAT(II,IENT) = DIR2(1,IENT)*UU+DIR2(2,IENT)*VV+DIR2(3,IENT)*WW
      ENDDO
C
C --- Verfication des donnees
C
      IOK = 0
      DO II = 1, NDAT(IENT)
        IF(UDAT(II,IENT).LE.0.D0.OR.KDAT(II,IENT).LE.0.D0.OR.
     &       EPSDAT(II,IENT).LE.0.D0) THEN
          WRITE(NFECRA,3100) IENT
          CALL CSEXIT (1)
        ENDIF
        IF(ICAS(IENT).EQ.1) THEN
          IF(YDAT(II,IENT).LT.-LLY(IENT)/2.D0.OR.
     &         YDAT(II,IENT).GT.LLY(IENT)/2.D0.OR.
     &         ZDAT(II,IENT).LT.-LLZ(IENT)/2.D0.OR.
     &         ZDAT(II,IENT).GT.LLZ(IENT)/2.D0) THEN
            IOK = IOK + 1
          ENDIF
        ELSEIF(ICAS(IENT).EQ.2) THEN
          IF(YDAT(II,IENT).LT.-LLD(IENT)/2.D0.OR.
     &         YDAT(II,IENT).GT.LLD(IENT)/2.D0.OR.
     &         ZDAT(II,IENT).LT.-LLD(IENT)/2.D0.OR.
     &         ZDAT(II,IENT).GT.LLD(IENT)/2.D0) THEN
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
      IF(IOK.GT.0) THEN
        WRITE(NFECRA,3200) IENT
      ENDIF
C
C=======================================================================
C 4. LECTURE DU FICHIER SUITE / INITIALISATION DU CHAMP DE VORTEX
C=======================================================================
C
      IF(ISUIVO.EQ.1) THEN
C
        IF(IENT.EQ.1) THEN
          OPEN(UNIT=IMPMVO,FILE=FICMVO)
          REWIND(IMPMVO)
        ENDIF
C
        READ(IMPMVO,100) IIENT
        READ(IMPMVO,100) IVORT
        IF(IVORT.NE.NVOR.OR.IIENT.NE.IENT) THEN
          WRITE(NFECRA,4500) IENT, IVORT, NVOR
          INITVO(IENT) = 1
        ELSE
          DO II = 1, NVOR
            READ(IMPMVO,200) YZVOR(II,1), YZVOR(II,2),
     &             TEMPS(II), TPSLIM(II), SIGNV(II)
          ENDDO
          INITVO(IENT) = 0
          WRITE(NFECRA,4000)
        ENDIF
C
        IF(IENT.EQ.NNENT) THEN
          CLOSE(IMPMVO)
        ENDIF
C
      ENDIF
C
      IF(ISUIVO.EQ.0.OR.INITVO(IENT).EQ.1) THEN
C
C-----------------------
C  Tirage des positions
C-----------------------
        IUN = 1
        IF(ICAS(IENT).EQ.1)THEN
          DO II = 1, NVOR
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,1) = LLY(IENT) * DRAND(1) - LLY(IENT)/2.D0
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,2) = LLZ(IENT) * DRAND(1) - LLZ(IENT)/2.D0
          ENDDO
        ELSEIF(ICAS(IENT).EQ.2) THEN
          DO II = 1, NVOR
 15         CONTINUE
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,1) = LLD(IENT) * DRAND(1) - LLD(IENT)/2.0D0
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,2) = LLD(IENT) * DRAND(1) - LLD(IENT)/2.0D0
            IF ((YZVOR(II,1)**2+YZVOR(II,2)**2).GT.
     &           (LLD(IENT)/2.D0)**2) THEN
              GOTO 15
            ENDIF
          ENDDO
        ELSEIF(ICAS(IENT).EQ.3.OR.ICAS(IENT).EQ.4) THEN
          DO II = 1, NVOR
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,1) = YMIN(IENT) + LLY(IENT) * DRAND(1)
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,2) = ZMIN(IENT) + LLZ(IENT) * DRAND(1)
          ENDDO
        ENDIF
C--------------
C Duree de vie
C--------------
        IF(ITLIVO(IENT).EQ.1) THEN
          DO II = 1, NVOR
            CALL ZUFALL(IUN,DRAND(1))
            TEMPS(II) = DRAND(1)*TLIMVO(IENT)
C
C on fait cela pour que les vortex ne disparaissent pas tous
C en meme temps.
C
            TPSLIM(II) = TLIMVO(IENT)
          ENDDO
        ELSEIF(ITLIVO(IENT).EQ.2) THEN
          DO II = 1, NVOR
            YY = YZVOR(II,1)
            ZZ = YZVOR(II,2)
            III = 0
            U_VOR  =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &                YDAT(1,IENT),ZDAT(1,IENT),UDAT(1,IENT),III)
            EK_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &                YDAT(1,IENT),ZDAT(1,IENT),KDAT(1,IENT),III)
            EE_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &                YDAT(1,IENT),ZDAT(1,IENT),EPSDAT(1,IENT),III)
            TPSLIM(II) = 5.D0*CMU*EK_VOR**(3.D0/2.D0)/EE_VOR
            TPSLIM(II) = TPSLIM(II)/U_VOR
            TEMPS(II) = 0.D0
          ENDDO
        ENDIF
C------------------
C Sens de rotation
C------------------
        DO II = 1, NVOR
          SIGNV(II) = 1.D0
          CALL ZUFALL(IUN,DRAND(1))
          IF (DRAND(1).LT.0.5D0) SIGNV(II) = -1.D0
        ENDDO
      ENDIF
C
C=======================================================================
C 5. AJOUT DE LA VITESSE MOYENNE POUR U
C=======================================================================
C
      DO II = 1, NCEVOR
        YY = YZCEL(II,1)
        ZZ = YZCEL(II,2)
        III = 0
        U_VOR = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &          YDAT(1,IENT),ZDAT(1,IENT),UDAT(1,IENT),III)
        XU(II)= U_VOR
      ENDDO
C
C=======================================================================
C 6. RECHERCHE DE LA FACE LA PLUS PROCHE DE CHAQUE VORTEX
C=======================================================================
      DO II = 1, NVOR
        KK = 0
        DD = GRAND
        DO JJ = 1, NCEVOR
          IF(((YZCEL(JJ,1)-YZVOR(II,1))**2+
     &        (YZCEL(JJ,2)-YZVOR(II,2))**2).LT.DD)THEN
            DD = (YZCEL(JJ,1)-YZVOR(II,1))**2
     &          +(YZCEL(JJ,2)-YZVOR(II,2))**2
            KK = JJ
          ENDIF
        ENDDO
        IVORCE(II) = KK
      ENDDO
C
C FORMATS
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(
     &'                                                             ',/,
     &' ** METHODE DES VORTEX                                       ',/,
     &'    ------------------                                       ',/,
     &'       NNENT  = ',4X,I10,    ' (Nombre d entrees            )',/,
     &'       ISUIVO = ',4X,I10,    ' (1 : suite de calcul         )'  )
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LES DIMENSIONS MAX DE L''ENTREE SONT INCOMPATIBLES AVEC ',/,
     &'@    LES DONNEES A L''ENTREE ',I10                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3000 FORMAT(
     &'                                                             ',/,
     &' --  Fin de la lecture du fichier de donnees                 ',/)
 3100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    U, K ET EPSILON SONT DES GRANDEURS QUI SONT DEFINIES    ',/,
     &'@    POSITIVES DANS LE REPERE LOCAL DE L''ENTREE             ',/,
     &'@                                                            ',/,
     &'@    VERIFIER LE FICHIER DE DONNEE DE L''ENTREE ',I10         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LES DIMENSIONS MAX DE L''ENTREE SONT INCOMPATIBLES AVEC ',/,
     &'@    CELLES DU FICHIER DE DONNEES                            ',/,
     &'@                                                            ',/,
     &'@    VERIFIER LE FICHIER DE DONNEE DE L''ENTREE ',I10         ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4000 FORMAT(
     &' --  Fin de la lecture du fichier suite                      ',/)
C
 4500   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@  LE NOMBRE DE VORTEX A CHANGE A L''ENTREE',I10              ,/,
     &'@    NVORT VALAIT PRECECEDEMENT ',I10                         ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@    ET VAUT MAINTENANT  ',I10                                ,/,
     &'@                                                            ',/,
     &'@  LA METHODE EST REINITIALISE A CETTE ENTREE                ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 100  FORMAT(I10)
 200  FORMAT(5E13.5)
C
#else
C
 1000 FORMAT(
     &'                                                             ',/,
     &' ** VORTEX METHOD                                            ',/,
     &'    -------------                                            ',/,
     &'       NNENT  = ',4X,I10,    ' (Number of inlets            )',/,
     &'       ISUIVO = ',4X,I10,    ' (1: calculation restart      )'  )
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE MAX DIMENSIONS OF THE INLET ARE INCOMPATIBLE WITH   ',/,
     &'@    THE DATA AT THE INLET ',I10                              ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3000 FORMAT(
     &'                                                             ',/,
     &' --  End reading the data file                               ',/)
 3100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    U, K AND EPSILON ARE QUANTITIES WHICH MUST BE POSITIVE  ',/,
     &'@    IN THE LOCAL FRAME OF THE INLET                         ',/,
     &'@                                                            ',/,
     &'@    VERIFY THE DATA FILE FOR THE INLET ',I10                 ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE MAX DIMENSIONS OF THE INLET ARE INCOMPATIBLE WITH   ',/,
     &'@    THE ONES FROM THE DATA FILE                             ',/,
     &'@                                                            ',/,
     &'@    VERIFY THE DATA FILE FOR THE INLET ',I10                 ,/,
     &'@                                                            ',/,
     &'@  The calculation will be run.                              ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4000 FORMAT(
     &' --  End reading the restart file                            ',/)
C
 4500   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@  THE NUMBER OF VORTICES HAS CHANGED AT THE INLET ',I10      ,/,
     &'@    NVORT WAS PREVIOUSLY ',I10                               ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@    IT IS CURRENTLY ',I10                                    ,/,
     &'@                                                            ',/,
     &'@  THE METHOD IS RE-INITIALIZED AT THIS INLET                ',/,
     &'@                                                            ',/,
     &'@  The calculation will be run.                              ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 100  FORMAT(I10)
 200  FORMAT(5E13.5)
C
#endif
C
      RETURN
      END
c@z
