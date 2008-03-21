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
C                   LIBRAIRIE DE SOUS-PROGRAMMES RMODAK
C                   ***********************************
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
CFONC    -----------------------------------------
CFONC
CFONC    SOUS-PROGRAMMES DE LA PHYSIQUE PARTICULIERE
CFONC     RELATIFS AU CALCUL DU COEFFICIENT D'ABSORPTION AVEC MODAK
CFONC
CFONC References : MODAK A.T.,
CFONC              "Radiation from products of combustion",
CFONC              Fire Research, 1 pp. 339-361, 1978.
CFONC
CFONC              MECHITOUA N.
CFONC              "Modelisation numerique du rayonnment dans les
CFONC               milieux semi-transoparents",
CFONC              Raport EDF, HE/44/87-15, 1987.
CFONC
c@fonce
C
C **********************************************************************
C
c@a
                        SUBROUTINE ABSORB
C                       *****************
C     -------------------------------------------------------------
     & ( TS    , TE     , PATH   , SOOTK  ,
     &   PCO2  , PH2O   , ALPHA  )
C     -------------------------------------------------------------
C
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     ON CALCULE LES ABSORPTIVITES (PAR RAPPORT
CFONC     A UNE SOURCE DE CORPS NOIR) D'UN MELANGE GAZEUX ISOTHERME,
CFONC     HOMOGENE DE SUIE, CO2 ET H2O A LA PRESSION TOTALE D'1 ATM.
CFONC
CFONC     SI LA TEMPERATURE DU CORPS NOIR EST EGALE A LA TEMPERATURE DU
CFONC     MELANGE, L'ABSORPTIVITE EST EGALE A L'EMISSIVITE.
CFONC     LES EMISSIVITES AINSI CALCULEES SONT EN BON ACCORD AVEC LES
CFONC     CALCULS SPECTRAUX ET LES MESURES EXPERIMENTALES
CFONC
CFONC     TS ET TE DOIVENT ETRE COMPRIS ENTRE 300 ET 2000 KELVIN
CFONC
CFONC     LA LONGUEUR D'ONDE 0.94 MICRON.
CFONC     SOOTK EST LIEE A LA FRACTION VOLUMIQUE DE SUIE FV SUIVANT
CFONC     LA FORMULE :
CFONC                    SOOTK=7FV/0.94E-6
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! TS           ! R  !  -> ! TEMPERATURE DU CORPS NOIR (K)        !
CARGU ! TE           ! R  !  -> ! TEMPERATURE DU MELANGE (K)           !
CARGU ! PATH         ! R  !  -> ! PENETRATION DU RAYONNEMENT DANS LE   !
CARGU !              !    !     ! MELANGE (M)                          !
CARGU ! SOOTK        ! R  !  -> ! COEFFICIENT D'ABSORPTION DES SUIES   !
CARGU ! PCO2         ! R  !  -> ! PRESSION PARTIELLE DE CO2 DANS UN    !
CARGU !              !    !     !  MELANGE DE PRESSSION TOTALE 1 ATM.  !
CARGU ! PH2O         ! R  !  -> ! PRESSION PARTIELLE DE H2O DANS UN    !
CARGU !              !    !     !  MELANGE DE PRESSSION TOTALE 1 ATM.  !
CARGU ! ALPHA        ! R  ! <-  ! ABSORPTIVITE
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
      DOUBLE PRECISION TS, TE, PATH, SOOTK, PCO2, PH2O, ALPHA
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION TMAX, TMIN, PTOTAL, RATIO, PATHL, PCL, PWL
      DOUBLE PRECISION AS, TAUS, AG, POWER, ZETA
      DOUBLE PRECISION EMIGAS
C
C=======================================================================
C CALCUL
C=======================================================================
C
      TMAX = 3000.D0
      TMIN = 298.D0
      IF (TS.LT.TMIN  .OR.  TS.GT.TMAX ) GOTO 1
      IF (TE.LT.TMIN  .OR.  TE.GT.TMAX ) GOTO 2
C
C --- Pression totale : PTOTAL
C
      PTOTAL = PCO2 + PH2O
C
      IF ( PTOTAL.GT.1.D0 ) GOTO 3
C
C --- Rapport temeperature melange et temperature source : RATIO
C
      RATIO = TE/TS
C
C --- Loncueur de penetration du rayonment effectif : PATHL
C
      PATHL = PATH/RATIO
      PCL   = PCO2*PATHL
      PWL   = PH2O*PATHL
      IF (PCL.GT.5.98D0 .OR. PWL.GT.5.98D0) GOTO 4
C
C --- Calcul de l'absortivite des suies : AS
C
      AS = 0.D0
      IF (SOOTK.LE.0.D0) GOTO 51
      CALL TASOOT
C     ===========
     & ( SOOTK , PATH , TS , TAUS )
C
      AS = 1.D0-TAUS
C
 51   CONTINUE
C
C --- Calcul de l'absorptivite du gaz : AG
C                 = emissivite du gaz
C
      AG = 0.D0
      IF (PCO2.LT.0.0011D0  .AND.   PH2O.LT.0.0011D0) GOTO 52
      IF (PCL .LT.0.0011D0  .AND.   PWL .LT.0.0011D0) GOTO 52
      AG = EMIGAS( PATHL, PCO2, PH2O, TS )
C
C --- Calcul de la fraction de vapeur d'eau : ZETA
C
      ZETA  = PH2O/PTOTAL
      POWER = 0.65D0-0.2D0*ZETA
      AG    = AG*(RATIO**POWER)
C
 52   CONTINUE
      ALPHA = AS + AG -AS*AG
C     ALPHA = ABS(ALPHA)
      IF (ALPHA.LE.1.D-8) GOTO 8
C
      RETURN
C
 4    CONTINUE
      WRITE(NFECRA,1000)
      GOTO 8
C
 3    CONTINUE
      WRITE(NFECRA,1001)
      GOTO 8
C
 2    CONTINUE
      WRITE(NFECRA,1002)
      GOTO 8
C
 1    CONTINUE
      WRITE(NFECRA,1003)
C
 8    CONTINUE
      ALPHA= 1.D-8
C
C
C========
C FORMATS
C========
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ERREUR MODAK :                                          ',/,
     &'@    ************                                            ',/,
     &'@ LE PRODUIT PATH*TS/T*PCO2 OU PATH*TS/T*PH2O                ',/,
     &'@ DEPASSE LA VALEUR 5.98 ATM.METRE.                          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ERREUR MODAK :                                          ',/,
     &'@    ************                                            ',/,
     &'@ LA SOMME DES PRESSIONS PARTIELLES DES GAZ CO2 ET H2O       ',/,
     &'@ DEPASSE UN ATMOSPHERE.                                     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ERREUR MODAK :                                          ',/,
     &'@    ************                                            ',/,
     &'@ LA TEMPERATURE DU MELANGE TE                               ',/,
     &'@ SORT DES LIMITES DU DOMAINE.                               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ERREUR MODAK :                                          ',/,
     &'@    ************                                            ',/,
     &'@ LA TEMPERATURE DU CORPS NOIR TS                            ',/,
     &'@ SORT DES LIMITES DU DOMAINE.                               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
C
C **********************************************************************
c@a
C
                        SUBROUTINE CHEBYC
C                       *****************
C     -------------------------------------------------------------
     & ( NORPOL , ARGPOL , VALPOL )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC      CALCUL DU POLYNOME DE CHEBYCHEV D'ORDRE NORPOL
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NORPOL       ! E  !  -> ! ORDRE DU POLYNOME DE  CHEBYCHEV      !
CARGU ! ARGPOL       ! R  !  -> ! ARGUMENT DU POLYNOME DE CHEBYCHEV    !
CARGU ! VALPOL       ! R  ! <-  ! VALEUR DU POLYNOME DE CHEBYCHEV      !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NORPOL
      DOUBLE PRECISION ARGPOL , VALPOL
C
C VARIABLES LOCALES
C
      INTEGER          ICT
      DOUBLE PRECISION F, VM2,VM1
C
C=======================================================================
C CALCUL
C=======================================================================
C
      VALPOL = 1.D0
      IF (NORPOL.LE.0) THEN
        GOTO 1
      ELSE
        GOTO 2
      ENDIF
C
 1    RETURN
C
 2    VALPOL = ARGPOL
C
      IF ( (NORPOL-1).LE.0) THEN
        GOTO 1
      ELSE
        GOTO 3
      ENDIF
 3    F   = ARGPOL+ARGPOL
      VM1 = ARGPOL
      VM2 = 1.D0
      DO ICT = 2, NORPOL
        VALPOL = F*VM1-VM2
        VM2    = VM1
        VM1    = VALPOL
      ENDDO
C
      RETURN
      END
c@z
C
C **********************************************************************
C
c@a
                       SUBROUTINE ASYMPT
C                      *****************
C     -------------------------------------------------------------
     & ( ZZ    , ZZV    )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONCC   CALCUL DE L'EXPANSION ASYMPTOTIQUE POUR LA FONCTION PENTAGAMMA
CFONCC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! ZZ           ! R  !  -> !                                      !
CARGU ! ZZV          ! R  ! <-  !                                      !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION ZZ, ZZV
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION ZI1, ZI2, ZI3, D1S3
C
C=======================================================================
C CALCUL
C=======================================================================
C
      D1S3 = 1.D0/3.D0
      ZI1 = 1.D0/ZZ
      ZI2 = ZI1*ZI1
      ZI3 = ZI1*ZI2
      ZZV = ZI3 * (
     &              (2.D0+3.D0*ZI1) +
     &              ZI2*(2.D0 + ZI2*(-1.D0 +
     &                               ZI2*( 1.D0+D1S3
     *                                     + ZI2*(-3.D0+10.D0*ZI2) )
     &                              )
     &                  )
     &            )
C
      RETURN
      END
c@z
C
C **********************************************************************
C
c@a
                       SUBROUTINE TASOOT
C                      *****************
C     -------------------------------------------------------------
     & ( ZKLED  , PATHL  , TBLACK , TAUS   )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC   CALCUL DE LA TRANSMISSIVITE(TAUS) DE PATH
CFONC   A UNE TEMPERATURE DONNEE.
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! ZKLED        ! R  !  -> !                                      !
CARGU ! PATHL        ! R  !  -> ! PENETRATION DU RAYONNEMENT DANS LE   !
CARGU !              !    !     ! MELANGE                              !
CARGU ! TBLACK       ! R  !  -> ! TEMPERATURE SOURCE OU TEMPERATURE    !
CARGU !              !    !     !  DU GAZ                              !
CARGU ! TAUS         ! R  ! <-  ! TRANSMISSIVITE DE PATH               !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION ZKLED, PATHL, TBLACK, TAUS
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION ARG , VAL
C
C=======================================================================
C CALCUL
C=======================================================================
C
      IF ( ZKLED.LE.0.D0 ) GOTO 1
C
      ARG = 1.D0 + ZKLED*PATHL*TBLACK*6.5333D-5
C
      CALL PENTAG
C     ===========
     & ( ARG    , VAL )
C
      TAUS = VAL*.1539897336D0
      RETURN
C
 1    TAUS = 1.D0
C
      RETURN
      END
c@z
C
C **********************************************************************
C
c@a
                       SUBROUTINE PENTAG
C                      *****************
C     -------------------------------------------------------------
     & ( ARGFPE , VALFPE )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    CALCUL DE LA VALEUR VAL DE LA
CFONC    FONCTION PENTAGAMMA D ARGUMENT X.
CFONC    ON UTILISE LES FORMULES ASYMPTOTIQUES ET DE RECURRENCE
CFONC    D'ABRAMOWITZ ET STEGUN.
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! ARGFPE       ! R  !  -> ! ARGUMENT DE LA FONCTION PENTAGAMMA   !
CARGU ! VALFPE       ! R  ! <-  ! VALEUR DE LA FONCTION PENTAGAMMA     !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION ARGFPE, VALFPE
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION ZZ, ZZV, ZS
C
C=======================================================================
C   CALCUL
C=======================================================================
C
      IF (ARGFPE.GE.4.D0) GOTO 1
      IF (ARGFPE.GE.3.D0) GOTO 2
      IF (ARGFPE.GE.2.D0) GOTO 3
C
      ZS = ( 1.D0/(ARGFPE+2.D0)**4 + 1.D0/(ARGFPE+1.D0)**4
     &                             + 1.D0/ARGFPE**4          )*6.D0
      ZZ = ARGFPE+3.D0
      CALL ASYMPT
C     ===========
     & ( ZZ, ZZV )
C
      GOTO 4
C
 3    CONTINUE
      ZS = (1.D0/(ARGFPE+1.D0)**4+1.D0/ARGFPE**4)*6.D0
      ZZ = ARGFPE+2.D0
      CALL ASYMPT
C     ===========
     & ( ZZ, ZZV )
      GOTO 4
C
 2    CONTINUE
      ZS= 6.D0/ARGFPE**4
      ZZ= ARGFPE+1.D0
      CALL ASYMPT
C     ===========
     & ( ZZ, ZZV )
C
      GOTO 4
C
 1    CONTINUE
      ZS = 0.D0
C La ligne suivante est ajoutee simplement pour eviter un warning
C    sous Foresys (var non initialisee) en attendant que la routine
C    soit revue de maniere globale (pour l'instant, l'utilisation
C    de rmodak est bloquee en amont).
      ZZ= ARGFPE
C
      CALL ASYMPT
C     ===========
     & ( ZZ, ZZV )
C
 4    CONTINUE
      VALFPE = ZZV+ZS
C
      RETURN
      END
c@z
C
C **********************************************************************
C
c@a
                       FUNCTION FDLECK
C                      ***************
C     -------------------------------------------------------------
     & ( VAL    , PL     , TE     )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC   CALCUL DE LA CORRECTION APPORTEE POUR LE MELANGE
CFONC   DE CO2 ET H2O LORSQUE LES LONGUEURS D ONDES SONT AU DELA DE 2.7
CFONC   ET 15 MICRONS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! VAL          ! R  ! ->  !                                      !
CARGU ! PL           ! R  ! ->  !                                      !
CARGU ! TE           ! R  ! ->  !                                      !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION VAL, PL, TE, FDLECK
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION TERM, TERM2, TERM3, TT, TT2, AA, BB, CC
C
C=======================================================================
C   CALCUL
C=======================================================================
C
      IF (PL.LT.0.1D0) GOTO 1
C
      TERM  = VAL/(10.7D0+101.D0*VAL) - VAL**10.4D0/111.7D0
      TERM2 = LOG10(101.325D0*PL)
      TERM2 = TERM2**2.76D0
      TT    = TE/1000.D0
      TT2   = TT*TT
      AA    = -1.0204082D0
      BB    =  2.2448979D0
      CC    = -0.23469386D0
      TERM3 = AA*TT2+BB*TT+CC
C
C --- TERM3 represente l'ajustement de temperature
C
      FDLECK = TERM*TERM2*TERM3
C
      RETURN
C
 1    FDLECK= 0.D0
C
      RETURN
      END
c@z
C
C **********************************************************************
C
c@a
                       FUNCTION EMIGAS
C                      ***************
C      -------------------------------------------------------------
     & ( PATHL  , PC     , PW     , TE     )
C      -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     CALCUL DE L EMISSIVITE A UN PATH DONNE
CFONC     D'UN MELANGE DE CO2 ET H2O A LA TEMPERATURE TE
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! PATHL        ! R  ! ->  ! VALEUR DU PATH                       !
CARGU ! PC           ! R  ! ->  ! PRESSION PARTIELLE DE CO2            !
CARGU ! PW           ! R  ! ->  ! PRESSION PARTIELLE DE H2O            !
CARGU ! TE           ! R  ! ->  ! TEMPERATURE                          !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION PATHL, PC, PW, TE, EMIGAS
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION TMIN, TMAX, PCL, EC
      DOUBLE PRECISION PWL, PCWL, DELS, EW, PCPW, XI, FDLECK
C
C=======================================================================
C   CALCUL
C=======================================================================
C
      TMIN   = 298.D0
      TMAX   = 3000.D0
      EMIGAS = 0.D0
C
      IF ( TE.LT.TMIN .OR. TE.GT.TMAX ) RETURN
C
      EC = 0.D0
C
      IF ( PC.LT.0.0011D0 .OR. PC.GT.1.D0 ) GOTO 1
C
      PCL = PC*PATHL
C
      IF ( PCL.LT.0.0011D0 .OR. PCL.GT.5.98D0 ) GOTO 1
C
      CALL SCRTCH
C     ===========
     & ( PC , PCL , TE , 1 , EC)
C
 1    CONTINUE
C
      IF ( PW.LT.0.0011D0 .OR. PW.GT.1.D0 ) GOTO 2
C
      PWL = PW*PATHL
C
      IF ( PWL.LT.0.0011D0 .OR. PWL.GT.5.98D0 ) GOTO 2
C
      CALL SCRTCH
C     ===========
     & ( PW , PWL , TE , 2 , EW )
C
      EMIGAS = EC + EW
C
      IF ( EC.LE.0.0D0 ) RETURN
C
      PCPW = PC + PW
      XI   = PW / PCPW
      IF ( XI.LT.0.01D0 ) RETURN
      PCWL = PCPW*PATHL
C
      IF ( PCWL.LT.0.1D0 ) RETURN
C
      DELS   = FDLECK(XI,PCWL,TE)
      EMIGAS = EMIGAS - DELS
C
      RETURN
C
2     CONTINUE
      EMIGAS = EC
C
      RETURN
      END
c@z
C
C **********************************************************************
C
c@a
                           SUBROUTINE SCRTCH
C                          *****************
C     -------------------------------------------------------------
     & ( PP     , PL     , TE     , INDEX  , VAL    )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! PP           ! R  ! ->  !                                      !
CARGU ! PL           ! R  ! ->  !                                      !
CARGU ! TE           ! R  ! ->  !                                      !
CARGU ! INDEX        ! E  ! ->  !                                      !
CARGU ! VAL          ! R  ! ->  !                                      !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          INDEX
      DOUBLE PRECISION PP, PL, TE, VAL
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, KK, III, JJJ, KKK
      DOUBLE PRECISION CC(3,4,4), CW(3,4,4), SC(3,4,4)
      DOUBLE PRECISION XX, YY, ZZ, V6, V7
      DOUBLE PRECISION TIX, TJY, TKZ
C
C=======================================================================
C   CALCUL
C=======================================================================
C
      IF ( INDEX.EQ.2 ) GOTO 2
C
C --- CC represente untebleau de 48 elements pour CO2
C
      CC(1,1,1) = -.2754568D1
      CC(1,1,2) = -.2997857D0
      CC(1,1,3) = -.1232494D0
      CC(1,1,4) =  .1279287D-1
      CC(1,2,1) =  .1503051D1
      CC(1,2,2) =  .3156449D0
      CC(1,2,3) =  .1058126D-1
      CC(1,2,4) = -.3729625D-1
      CC(1,3,1) = -.247411D0
      CC(1,3,2) = -.3323846D-1
      CC(1,3,3) = -.1819471D-1
      CC(1,3,4) =  .2289789D-1
      CC(1,4,1) =  .4994029D-1
      CC(1,4,2) = -.1986786D-2
      CC(1,4,3) =  .3007898D-2
      CC(1,4,4) = -.1175598D-2
      CC(2,1,1) =  .5737722D-2
      CC(2,1,2) = -.9328458D-2
      CC(2,1,3) =  .2906286D-2
      CC(2,1,4) =  .422752D-3
      CC(2,2,1) = -.3151784D-2
      CC(2,2,2) =  .5632821D-2
      CC(2,2,3) = -.3260295D-2
      CC(2,2,4) =  .7065884D-3
      CC(2,3,1) =  .1668751D-3
      CC(2,3,2) = -.7326533D-3
      CC(2,3,3) =  .3639855D-3
      CC(2,3,4) =  .3228318D-3
      CC(2,4,1) =  .7386638D-3
      CC(2,4,2) = -.7277073D-3
      CC(2,4,3) =  .5925968D-3
      CC(2,4,4) = -.2021413D-3
      CC(3,1,1) =  .3385611D-2
      CC(3,1,2) = -.5439185D-2
      CC(3,1,3) =  .176456D-2
      CC(3,1,4) =  .3036031D-3
      CC(3,2,1) = -.18627D-2
      CC(3,2,2) =  .3236275D-2
      CC(3,2,3) = -.195225D-2
      CC(3,2,4) =  .3474022D-3
      CC(3,3,1) =  .1204807D-3
      CC(3,3,2) = -.4479927D-3
      CC(3,3,3) =  .2497521D-3
      CC(3,3,4) =  .1812996D-3
      CC(3,4,1) =  .4218169D-3
      CC(3,4,2) = -.4046608D-3
      CC(3,4,3) =  .3256861D-3
      CC(3,4,4) = -.9514981D-4
C
      GOTO 4
C
 2    CONTINUE
C
C --- CW represente untebleau de 48 elements pour H2O
C
      CW(1,1,1) = -.2594279D1
      CW(1,1,2) = -.7118472D0
      CW(1,1,3) = -.9956839D-3
      CW(1,1,4) =  .1226560D-1
      CW(1,2,1) =  .2510331D1
      CW(1,2,2) =  .6481808D0
      CW(1,2,3) = -.3330587D-1
      CW(1,2,4) = -.5524345D-2
      CW(1,3,1) = -.4191636D0
      CW(1,3,2) = -.1375180D0
      CW(1,3,3) =  .3877930D-1
      CW(1,3,4) =  .8862328D-3
      CW(1,4,1) = -.322912D-1
      CW(1,4,2) = -.1820241D-1
      CW(1,4,3) = -.2223133D-1
      CW(1,4,4) = -.5940781D-3
      CW(2,1,1) =  .1126869D0
      CW(2,1,2) = -.8133829D-1
      CW(2,1,3) =  .1514940D-1
      CW(2,1,4) =  .1393980D-2
      CW(2,2,1) = -.9298805D-2
      CW(2,2,2) =  .4550660D-1
      CW(2,2,3) = -.2082008D-1
      CW(2,2,4) =  .2013361D-2
      CW(2,3,1) = -.4375032D-1
      CW(2,3,2) =  .1924597D-1
      CW(2,3,3) =  .8859877D-2
      CW(2,3,4) = -.4618414D-2
      CW(2,4,1) =  .7077876D-2
      CW(2,4,2) = -.2096188D-1
      CW(2,4,3) =  .1458262D-2
      CW(2,4,4) =  .3851421D-2
      CW(3,1,1) =  .5341517D-1
      CW(3,1,2) = -.3407693D-1
      CW(3,1,3) =  .4354611D-2
      CW(3,1,4) =  .1492038D-2
      CW(3,2,1) = -.4708178D-2
      CW(3,2,2) =  .2086896D-1
      CW(3,2,3) = -.9477533D-2
      CW(3,2,4) =  .6153272D-3
      CW(3,3,1) = -.2104622D-1
      CW(3,3,2) =  .7515796D-2
      CW(3,3,3) =  .5965509D-2
      CW(3,3,4) = -.2756144D-2
      CW(3,4,1) =  .4318975D-2
      CW(3,4,2) = -.1005744D-1
      CW(3,4,3) =  .4091084D-3
      CW(3,4,4) =  .2550435D-2
C
 4    CONTINUE
C
      XX = LOG(PP)/3.45D0 + 1.D0
      YY  = (LOG(PL)+2.555D0) / 4.345D0
      ZZ  = (TE-1150.D0) / 850.D0
      VAL = 0.D0
C
      DO II =1, 3
        III = II-1
        CALL CHEBYC
C       ===========
     &  ( III , XX , TIX )
C
        V6 = 0.D0
        DO JJ = 1, 4
          JJJ = JJ-1
          CALL CHEBYC
C         ===========
     &    ( JJJ, YY , TJY )
C
          V7 = 0.D0
          DO KK = 1, 4
            KKK = KK-1
            CALL CHEBYC
C           ===========
     &      ( KKK , ZZ , TKZ )
C
            IF (INDEX.EQ.1) SC(II,JJ,KK) = CC(II,JJ,KK)
            IF (INDEX.EQ.2) SC(II,JJ,KK) = CW(II,JJ,KK)
            V7 = V7 + TKZ*SC(II,JJ,KK)
          ENDDO
C
            V6 = V6 + V7*TJY
        ENDDO
C
        VAL = VAL + V6*TIX
C
      ENDDO
C
      VAL = EXP(VAL)
C
      RETURN
      END
c@z
