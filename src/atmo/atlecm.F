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
                        SUBROUTINE ATLECM
C                       *****************
C     -------------------------------------------------------------
     & ( IMODE  , TMPROM , ZTPROM , ZDPROM ,
     &   XMET   , YMET   , PMER   ,
     &   TTPROM , QVPROM ,
     &   UPROM  , VPROM  , EKPROM , EPPROM,
     &   RPROM  , TPPROM , PHPROM   )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC
CFONC              READING THE METEO PROFILE DATA
CFONC
CFONC              SPECIFIC PHYSICS ATMOSPHERIC MODULE
CFONC
CFONC              IMODE = 0 : READING OF DIMENSIONS ONLY
CFONC              IMODE = 1 : READING OF DATA
CFONC
CFONC              WARNING : ARGUMENTS ARE DEFINED ONLY ON SECOND CALL
CFONC                OF THE ROUTINE (IMODE=1)
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "atincl.h"
C
C***********************************************************************
C
C ARGUMENTS
      INTEGER           IMODE
      DOUBLE PRECISION  TMPROM(NBMETM)
      DOUBLE PRECISION  ZTPROM(NBMETT) , ZDPROM(NBMETD)
      DOUBLE PRECISION  XMET(NBMETM)   , YMET(NBMETM)  , PMER(NBMETM)
      DOUBLE PRECISION  TTPROM(NBMETT,NBMETM) , QVPROM(NBMETT,NBMETM)
      DOUBLE PRECISION  UPROM(NBMETD,NBMETM)  , VPROM(NBMETD,NBMETM)
      DOUBLE PRECISION  EKPROM(NBMETD,NBMETM) , EPPROM(NBMETD,NBMETM)
      DOUBLE PRECISION  RPROM(NBMETT,NBMETM)  , TPPROM(NBMETT,NBMETM)
      DOUBLE PRECISION  PHPROM(NBMETT,NBMETM)
C
C VARIABLES LOCALES
C
      INTEGER ITP, II, IOS, K, IPHAS
      DOUBLE PRECISION ANNEE,QUANT,HEURE,MINUTE,SECOND
      DOUBLE PRECISION RAIR,PSOL,RAP,RSCP,TMOY
      CHARACTER*80     CCOMNT
      CHARACTER*1      CSAUTE
C
C***********************************************************************
C
      IF (IMODE.EQ.0) THEN
        WRITE(NFECRA,*) 'reading of dimensions for meteo profiles'
      ELSE
        WRITE(NFECRA,*) 'reading of meteo profiles data'
       ENDIF
C
C=======================================================================
C 0. INITIALISATION
C=======================================================================

      CSAUTE = '/'
C
C --> Ouverture du fichier
      OPEN ( UNIT=IMPMET, FILE=FICMET,
     &       STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &                                          IOSTAT=IOS, ERR=99 )
      REWIND(UNIT=IMPMET,ERR=99)
C
      ITP=0
C
C=======================================================================
C 1. LOOP ON TIME
C=======================================================================
C
 100  CONTINUE
      ITP = ITP+1
C
C ---> READING THE COMMENTS
C ---------------------------

 101  READ(IMPMET,'(A80)',ERR=999,END=906) CCOMNT
C
      IF(CCOMNT(1:1).EQ.CSAUTE) GO TO 101
      BACKSPACE(IMPMET)
C
C=======================================================================
C 2. READING THE TIME OF THE PROFILE
C=======================================================================
C
C --> annee, quantieme, heure, minute, seconde  du profil (TU)
C
      IF (IMODE.EQ.0) THEN
        READ(IMPMET,*,ERR=999,END=906)
      ELSE
        READ(IMPMET,*,ERR=999,END=906) ANNEE, QUANT,HEURE,MINUTE,SECOND
C
C --> calcul du temps relatif au debut de la simulation
C     TTPROM(NBMETM) = instants correspondant aux profils
C TD!!!! calcul par rapport ANNDEB, QUADEB, HEUDEB, MINDEB, SCDDEB
C + tests d'incohérence
C
        TMPROM(ITP)= SECOND
C
C
C --> verification de l'ordre chronologique des profils
C
        IF(ITP.GT.1) THEN
          IF(TMPROM(ITP).LT.TMPROM(ITP-1)) THEN
            WRITE(NFECRA,8000)
            CALL CSEXIT (1)
C           ===========
          ENDIF
        ENDIF
C
      ENDIF
C
C=======================================================================
C 3. READING THE POSITION OF THE PROFILE
C=======================================================================

 102  READ(IMPMET,'(A80)',ERR=999,END=999) CCOMNT
C
      IF(CCOMNT(1:1).EQ.CSAUTE) GO TO 102
      BACKSPACE(IMPMET)
C
C
      IF (IMODE.EQ.0) THEN
        READ(IMPMET,*,ERR=999,END=999)
      ELSE
        READ(IMPMET,*,ERR=999,END=999) XMET(ITP), YMET(ITP)
      ENDIF
C
C=======================================================================
C 4. READING THE SEA-LEVEL PRESSURE
C=======================================================================

 103  READ(IMPMET,'(A80)',ERR=999,END=999) CCOMNT
C
      IF(CCOMNT(1:1).EQ.CSAUTE) GO TO 103
      BACKSPACE(IMPMET)
C
C
      IF (IMODE.EQ.0) THEN
        READ(IMPMET,*,ERR=999,END=999)
      ELSE
        READ(IMPMET,*,ERR=999,END=999) PMER(ITP)
      ENDIF
C
C=========================================================================
C 5. READING THE TEMPERATURE AND HUMIDITY PROFILES
C=========================================================================
C
 104  READ(IMPMET,'(A80)',ERR=999,END=999) CCOMNT
C
      IF(CCOMNT(1:1).EQ.CSAUTE) GO TO 104
      BACKSPACE(IMPMET)
C
C
      IF (IMODE.EQ.0) THEN
        READ(IMPMET,*,ERR=999,END=999) NBMETT
C
        IF(NBMETT.LE.1) THEN
          WRITE(NFECRA,8001)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        DO II=1,NBMETT
          READ (IMPMET,*,ERR=999,END=999)
        ENDDO
C
      ELSE
C
        READ(IMPMET,*,ERR=999,END=999)
C
        DO II=1,NBMETT
C
C     Altitude, temperature, humidite
          READ (IMPMET,*,ERR=999,END=999) ZTPROM(II),
     &                                    TTPROM(II,ITP), QVPROM(II,ITP)
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 6. COMPUTING HYDRO PRESSURE PROFILE  (LAPLACE INTEGRATION)
C=======================================================================
C  NOTE : BASEE SUR LA PRESSION PMER
C         ET UNE INTEGRATION DE LAPLACE DE BAS EN HAUT
C
      IF (IMODE.EQ.1) THEN
        IPHAS = 1
        PHPROM(1,ITP)=PMER(ITP)
        RAIR=287.D0
        PSOL=P0(IPHAS)
        RSCP=RAIR/CP0(IPHAS)

        DO K=2,NBMETT
          TMOY=0.5*(TTPROM(K-1,ITP)+TTPROM(K,ITP))+TKELVI
C         RHMOY=RAIR*(1.+(RVSRA-1.)*
C    &         (QVPROM(K-1,ITP)+QVPROM(K,ITP))/2.*IH2O)
C         RAP=-ABS(GZ)*(ZK-ZKM1)/RHMOY/TMOY
          RAP=-ABS(GZ)*(ZTPROM(K)-ZTPROM(K-1))/RAIR/TMOY
          PHPROM(K,ITP)=PHPROM(K-1,ITP)*EXP(RAP)
        ENDDO
C
      ENDIF
C
C=======================================================================
C 7. COMPUTING THE POT. TEMPERATURE PROFILE AND THE DENSITY PROFILE
C=======================================================================
C
      IF (IMODE.EQ.1) THEN
        DO K=1,NBMETT
C         RHUM=RAIR*(1.D0+(RVSRA-1.D0)*QVPROM(K,ITP)*IH2O)

C         IF ((ISCALT(IPHAS).EQ.-1).OR.(IPHYSI.EQ.0)) THEN
C     MASSE VOLUMIQUE CONSTANTE
C           RPROM(K,ITP)=PMER/(TPROM(K,ITP)+TKELVI)/RHUM
C         ELSE
C     MASSE VOLUMIQUE VARIABLE
C           RPROM(K,ITP)=PHPROM(K,ITP)/(TTPROM(K,ITP)+TKELVI)/RHUM
C         ENDIF
          RPROM(K,ITP)=PHPROM(K,ITP)/(TTPROM(K,ITP)+TKELVI)/RAIR
C    &            *(1.D0+(RVSRA-CPVCPA)*QVPROM(K,ITP)*IH2O)
          TPPROM(K,ITP)=(TTPROM(K,ITP)+TKELVI)*
     &         ((PSOL/PHPROM(K,ITP))**RSCP)
        ENDDO
C
      ENDIF
C
C========================================================================
C 8. READING THE VELOCITY PROFILE
C========================================================================
C
C
 105  READ(IMPMET,'(A80)',ERR=999,END=999) CCOMNT
C
      IF(CCOMNT(1:1).EQ.CSAUTE) GO TO 105
      BACKSPACE(IMPMET)
C
C
      IF (IMODE.EQ.0) THEN
C
        READ(IMPMET,*,ERR=999,END=999) NBMETD
C
        IF(NBMETD.LE.1) THEN
          WRITE(NFECRA,8002)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        DO II=1,NBMETD
          READ (IMPMET,*,ERR=999,END=999)
        ENDDO
C
      ELSE
C
        READ(IMPMET,*,ERR=999,END=999)
C
        DO II=1,NBMETD
C     Altitude,  u,  v, k, epsilon
          READ (IMPMET,*,ERR=999,END=999) ZDPROM(II),
     &                                    UPROM(II,ITP),  VPROM(II,ITP),
     &                                    EKPROM(II,ITP), EPPROM(II,ITP)
C
        ENDDO
C
      ENDIF
C
C========================================================================
C 9. PRINTINGS
C========================================================================
C
      IF (IMODE.EQ.1) THEN
        IF(ITP.EQ.1) THEN
          WRITE(NFECRA,*)
          WRITE(NFECRA,*)
     &         '==================================================='
          WRITE(NFECRA,*)
     &         'printing meteo profiles'
        ENDIF
        WRITE(NFECRA,*)
     &       'annee,quantieme, heure, minute, seconde'
        WRITE(NFECRA,7995)
     &       ANNEE, QUANT,HEURE,MINUTE,SECOND
 7995   FORMAT(1X, 5F8.2)
        WRITE(NFECRA,*)
     &       'tmprom(itp)'
        WRITE(NFECRA,7996)
     &       TMPROM(ITP)
 7996   FORMAT(1X, F8.2)
        WRITE(NFECRA,*)
     &       'zdprom,uprom,vprom,ekprom,epprom'
        DO II=1,NBMETD
          WRITE(NFECRA,7997)
     &         ZDPROM(II), UPROM(II,ITP), VPROM(II,ITP),
     &         EKPROM(II,ITP), EPPROM(II,ITP)
 7997     FORMAT(1X, 3F8.2,2E10.3)
        ENDDO
        WRITE(NFECRA,*)
     &       'ztprom,ttprom,tpprom,rprom,phprom,qvprom'
        DO II=1,NBMETT
          WRITE(NFECRA,7998)
     &         ZTPROM(II), TTPROM(II,ITP), TPPROM(II,ITP),
     &         RPROM(II,ITP), PHPROM(II,ITP), QVPROM(II,ITP)
 7998     FORMAT(1X, 3F8.2,F8.4,F12.3,E10.3)
        ENDDO
C
      ENDIF
C
C========================================================================
C 10. END OF THE LOOP ON TIME
C========================================================================
C
      GOTO 100
C
 906  CONTINUE
C
      IF (IMODE.EQ.0) NBMETM= ITP-1
C
      CLOSE(UNIT=IMPMET)
C
C ---
C END
C ---
      RETURN
C
C====================
C XX. ERROR OUTPUTS
C====================
C
 99   CONTINUE
      WRITE ( NFECRA,9998 )
      CALL CSEXIT (1)
C     ===========
C
 999  CONTINUE
      WRITE ( NFECRA,9999 )
      CALL CSEXIT (1)
C     ===========
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
 8000 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ATDIMM)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE               ',/,
     &'@                                                            ',/,
     &'@              Erreur  dans le fichier meteo:                ',/,
     &'@      ordre chronologique des profils non respecte          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ATDIMM)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,
     &'@                                                            ',/,
     &'@              Erreur  dans le fichier meteo:                ',/,
     &'@  le nombre de mesures de temperatures doit être            ',/,
     &'@  superieur a 2                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8002 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ATDIMM)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,
     &'@                                                            ',/,
     &'@              Erreur  dans le fichier meteo:                ',/,
     &'@  le nombre de mesures de vitesses doit être                ',/,
     &'@  superieur a 2                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9998 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ATDIMM)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,
     &'@                                                            ',/,
     &'@  Erreur a l''ouverture du fichier meteo                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ATDIMM)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,
     &'@                                                            ',/,
     &'@  Erreur a la lecture du fichier meteo.                     ',/,
     &'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,
     &'@    ou son format inadapte.                                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 8000 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA (ATDIMM)      ',/,
     &'@     *******                                                ',/,
     &'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,
     &'@                                                            ',/,
     &'@              Error in the meteo profile file:              ',/,
     &'@      non chronological meteo profiles                      ',/,
     &'@                                                            ',/,
     &'@  The computation will not be run                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA (ATDIMM)      ',/,
     &'@     *******                                                ',/,
     &'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,
     &'@                                                            ',/,
     &'@              Error in the meteo profile file:              ',/,
     &'@  the number of temperature measurements must be larger     ',/,
     &'@  than 2                                                    ',/,
     &'@                                                            ',/,
     &'@  The computation will not be run                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8002 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA (ATDIMM)      ',/,
     &'@     *******                                                ',/,
     &'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,
     &'@                                                            ',/,
     &'@              Error in the meteo profile file:              ',/,
     &'@  the number of velocity measurements must be larger        ',/,
     &'@  than 2                                                    ',/,
     &'@                                                            ',/,
     &'@  The computation will not be run                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9998 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA (ATDIMM)      ',/,
     &'@     *******                                                ',/,
     &'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,
     &'@                                                            ',/,
     &'@  Error opening the meteo profile file                      ',/,
     &'@                                                            ',/,
     &'@  The computation will not be run                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA (ATDIMM)      ',/,
     &'@     *******                                                ',/,
     &'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,
     &'@                                                            ',/,
     &'@  Error opening the meteo profile file                      ',/,
     &'@    The meteo profile file has been opened but its content  ',/,
     &'@    is incomplete or under a wrong format                   ',/,
     &'@                                                            ',/,
     &'@  The computation will not be run                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
      END
C
C
c@z
