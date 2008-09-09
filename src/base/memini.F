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
                        SUBROUTINE MEMINI
C                       *****************
C     -------------------------------------------------------------
     & ( LONGIA , LONGRA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC          GESTION MEMOIRE INITIALE
CFONC          MACROS TABLEAUX, DEVELOPPEUR, UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! LONGIA LONGRA! E  !  -> ! LONGUEUR DE IA     RA                !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
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
      INCLUDE "dimens.h"
      INCLUDE "entsor.h"
      INCLUDE "ihmpre.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER LONGIA , LONGRA
      INTEGER NIDEVE , NRDEVE , NITUSE , NRTUSE
C
C VARIABLES LOCALES
C
      INTEGER II, IOK
      INTEGER ICOFTU(16)
C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
C --- Dimension des macros tableaux entier IA et reel RA
C
      LONGIA = 0
      LONGRA = 0
C
C --- Dimension des tableaux utilisateurs et developpeurs
C      supplementaires (entiers et reels) generalement un multiple
C      de NCEL, NFABOR ...
C
      NIDEVE = 0
      NRDEVE = 0
C
      NITUSE = 0
      NRTUSE = 0
C
      DO II = 1, 16
        ICOFTU(II) = 0
      ENDDO
C
C=======================================================================
C 1. MEMOIRE INITIALE POUR UTILISATEUR, DEVELOPPEUR
C=======================================================================
C
C   - Interface Code_Saturne
C     ======================
C
      IF(IIHMPR.EQ.1) THEN
C
        CALL UIUSAR(ICOFTU)
C       ===========
C
        NITUSE = ICOFTU(1)*NCELET + ICOFTU(2)*NFAC + ICOFTU(3)*NFABOR
     &         + ICOFTU(4)
        NRTUSE = ICOFTU(5)*NCELET + ICOFTU(6)*NFAC + ICOFTU(7)*NFABOR
     &         + ICOFTU(8)
        LONGIA = ICOFTU( 9)*NCELET + ICOFTU(10)*NFAC + ICOFTU(11)*NFABOR
     &         + ICOFTU(12)
        LONGRA = ICOFTU(13)*NCELET + ICOFTU(14)*NFAC + ICOFTU(15)*NFABOR
     &         + ICOFTU(16)
C
      ENDIF
C
C   - Sous-programme utilisateur
C     ==========================
C
      CALL USTBTR
C     ===========
     & ( NCEL   , NCELET , NFAC   , NFABOR , NNOD  ,
     &   LONGIA , LONGRA ,
     &   NIDEVE , NITUSE , NRDEVE , NRTUSE )
C
C
C=======================================================================
C 2. VERIFICATIONS
C=======================================================================
C
      IOK = 0
C
      IF (LONGIA.LT.0) THEN
        WRITE(NFECRA,1000)'LONGIA',LONGIA
        IOK = IOK + 1
      ENDIF
      IF (LONGRA.LT.0) THEN
        WRITE(NFECRA,1000)'LONGRA',LONGIA
        IOK = IOK + 1
      ENDIF
C
      IF (NIDEVE.LT.0) THEN
        WRITE(NFECRA,1000)'NIDEVE',NIDEVE
        IOK = IOK + 1
      ENDIF
      IF (NRDEVE.LT.0) THEN
        WRITE(NFECRA,1000)'NRDEVE',NRDEVE
        IOK = IOK + 1
      ENDIF
      IF (NITUSE.LT.0) THEN
        WRITE(NFECRA,1000)'NITUSE',NITUSE
        IOK = IOK + 1
      ENDIF
      IF (NRTUSE.LT.0) THEN
        WRITE(NFECRA,1000)'NRTUSE',NRTUSE
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.GT.0) CALL CSEXIT(1)
C                   ===========
C
C=======================================================================
C 3. DIMENSIONNEMENT MEMOIRE PAR DEFAUT
C=======================================================================
C
      IF (LONGIA.EQ.0) THEN
        LONGIA =  30 * MAX(NCELET, NFABOR)
      ENDIF
C
      IF (LONGRA.EQ.0) THEN
        LONGRA = 120 * MAX(NCELET, NFABOR)
      ENDIF
C
C----
C FORMATS
C----
C
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER POSITIF OU NUL               ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou ustbtr.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
