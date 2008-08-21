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
                        SUBROUTINE INITI1
C                       *****************
C     -------------------------------------------------------------
     & ( IVERIF )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DES COMMONS
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! LONGIB       ! E  ! <-> ! LONGUEUR DU TABLEAU  DES ENTIERS  IA !
CARGU ! LONGRB       ! E  ! <-> ! LONGUEUR DU TABLEAU  DES REELS    RA !
CARGU ! IDBIA0       ! E  ! <-  ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  ! <-  ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "matiss.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IVERIF
C
C VARIABLES LOCALES
C
      INTEGER          IOK
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION DES VARIABLES EN COMMON
C      AVANT INTERVENTION UTILISATEUR.
C      . mots-cles dynamiques
C      . entsor.h
C      . dimens.h
C      . numvar.h
C      . pointe.h
C      . optcal.h
C      . mltgrd.h
C      . cstphy.h
C=======================================================================
C
      CALL INIINI
C     ===========
C
C=======================================================================
C 2. ENTREE DES DONNEES PAR L'UTILISATEUR
C      ET POSITIONNEMENT DES VARIABLES (VARPOS)
C=======================================================================
C
      CALL INIUSI(IVERIF)
C     ===========
C
      CALL PPINI1
C     ===========
C
      CALL RAYOPT
C     ===========
C
      CALL LAGOPT
C     ===========
C
      CALL MTINI1
C     ===========
C
C En mode verification, on change IMRGRA en consequence, avant les
C   impressions dans le listing (on en a besoin dans MODINI pour ANOMAX)
C Si d'autres tests elementaires viennent se greffer, on pourra faire
C   une routine specifique d'initialisation des parametres des tests.
C On met aussi la sortie en ASCII, pour le cas ou on utilise autovalid
        IF (IVERIF.EQ.1) IMRGRA = 0
        IF (IVERIF.EQ.2) IMRGRA = 1
        IF (IVERIF.EQ.3) IMRGRA = 2
        IF (IVERIF.EQ.4) IMRGRA = 3
        IF (IVERIF.EQ.5) IMRGRA = 4
        IF (IVERIF.GT.0) OPTCHR = "text"
C
C=======================================================================
C 3. MODIFS APRES USINI1
C=======================================================================
C
      CALL MODINI
C     ===========
C
C=======================================================================
C 4. VERIFS APRES USINI1
C=======================================================================
C
      IOK = 0
C
      CALL VERINI (IOK)
C     ===========
C
      IF(IOK.GT.0) THEN
        WRITE(NFECRA,9999)IOK
        CALL CSEXIT (1)
      ELSE
        WRITE(NFECRA,9998)
      ENDIF
C
C
#if defined(_CS_LANG_FR)
C
 9998 FORMAT(
     &'                                                             ',/,
     &' Pas d erreur detectee lors de la verification des donnees   ',/,
     &'                               (interface, usini1 et autres).',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,
     &'@                                                            ',/,
     &'@  Se reporter aux impressions precedentes pour plus de      ',/,
     &'@    renseignements.                                         ',/,
     &'@  Verifier les donnees entrees dans l''interface, usini1 ou ',/,
     &'@    les autres sous-programmes d''initialisation.           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 9998 FORMAT(
     &''                                                             ,/,
     &' No error detected during the data verification'              ,/,
     &'                              (interface, usini1 and others).',/)
 9999 FORMAT(
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION'                ,/,
     &'@    ********'                                                ,/,
     &'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLET'  ,/,
     &'@'                                                            ,/,
     &'@  The calculation will not be run (',I10,' errors).'         ,/,
     &'@'                                                            ,/,
     &'@  See previous impressions for more informations.'           ,/,
     &'@  Verify the provided data in the interface, usini1 or'      ,/,
     &'@    the other initialization subroutines.'                   ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/)
C
#endif
C
C=======================================================================
C 5. IMPRESSIONS
C=======================================================================
C
      CALL IMPINI
C     ===========
C
      IF (IMATIS.EQ.1) THEN
        CALL MTIMPI
C       ===========
      ENDIF
C
C
      RETURN
      END
c@z
