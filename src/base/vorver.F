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
                        SUBROUTINE VORVER
C                       *****************
C     -------------------------------------------------------------
     & ( NFABOR , IREPVO , IAPPEL )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  VERIFICATION DES PARAMETRES DE CALCUL DE LA METHODE DES VORTEX
CFONC    APRES INTERVENTION UTILISATEUR
CFONC    (COMMONS)
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! IREPVO       ! TE !  -> ! NUMERO DE L'ENTREE ASSOCIE A CHAQUE  !
CARGU !     (NFABOR) !    !     ! FACE DE BORD (=0 SI PAS DE VORTEX)   !
CARGU ! IAPPEL       ! E  !  -> ! INDIQUE LES DONNES A VERIFIER        !
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
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "optcal.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NFABOR , IREPVO(NFABOR)  , IAPPEL
C
C VARIABLES LOCALES
C
      INTEGER          IOK, NMAX , II, JJ
      DOUBLE PRECISION NORM1, NORM2, CROSP
C
C***********************************************************************
C
      IOK = 0
C
C=======================================================================
C 1. IAPPEL = 1
C=======================================================================
C
      IF(IAPPEL.EQ.1) THEN
C
C --- Nombre d'entree
C
        IF (NNENT.GT.NENTMX.OR.NNENT.LE.0) THEN
          WRITE(NFECRA,1000) NENTMX, NNENT
          IOK = IOK + 1
        ENDIF
C
C --- Nombre de vortex
C
        DO II = 1, NNENT
          IF(NVORT(II).LE.0) THEN
            WRITE(NFECRA,2000) II, NVORT(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C --- Reperage des entrees
C
        NMAX = 0
        DO II = 1, NFABOR
          IF(IREPVO(II).LT.0.OR.IREPVO(II).GT.NNENT) THEN
            WRITE(NFECRA,2200) NNENT,IREPVO(II)
            IOK = IOK + 1
          ELSE
            NMAX = MAX(NMAX,IREPVO(II))
          ENDIF
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARCMX(NMAX)
        ENDIF
        IF(NMAX.EQ.0) THEN
          WRITE(NFECRA,2200) NNENT,NMAX
          IOK = IOK + 1
        ENDIF
C
C=======================================================================
C 2. IAPPEL = 2
C=======================================================================
C
      ELSEIF(IAPPEL.EQ.2) THEN
C
C --- Dimensions des Tableaux
C
        IF (NNENT.GT.NENTMX.OR.NNENT.LE.0) THEN
          WRITE(NFECRA,1000) NENTMX, NNENT
          IOK = IOK + 1
        ENDIF
C
        NMAX = 0
C
        DO II = 1, NNENT
          IF(ICAS(II).EQ.4) THEN
            NDAT(II) = 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(NDAT(II).LE.0) THEN
            WRITE(NFECRA,1100) NDATMX, NDAT(II), II
            IOK = IOK + 1
          ELSE
            NMAX = MAX(NMAX,NDAT(II))
          ENDIF
          IF(NMAX.GT.NDATMX) THEN
            WRITE(NFECRA,1100) NDATMX, NMAX, II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C --- Nombre de vortex
C
        DO II = 1, NNENT
          IF(NVORT(II).LE.0) THEN
            WRITE(NFECRA,2000) II, NVORT(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C --- Suite de calcul
C
        IF(ISUIVO.NE.0.AND.ISUIVO.NE.1) THEN
          WRITE(NFECRA,2100) ISUIVO
          IOK = IOK + 1
        ENDIF
C
C --- Reperage des entrees
C
        NMAX = 0
        DO II = 1, NFABOR
          IF(IREPVO(II).LT.0.OR.IREPVO(II).GT.NNENT) THEN
            WRITE(NFECRA,2200) NNENT,IREPVO(II)
            IOK = IOK + 1
          ELSE
            NMAX = MAX(NMAX,IREPVO(II))
          ENDIF
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARCMX(NMAX)
        ENDIF
        IF(NMAX.EQ.0) THEN
          WRITE(NFECRA,2200) NNENT,NMAX
          IOK = IOK + 1
        ENDIF
C
C --- Definition du cas
C
        DO II = 1, NNENT
          IF(ICAS(II).NE.1.AND.ICAS(II).NE.2.AND.
     &         ICAS(II).NE.3.AND.ICAS(II).NE.4) THEN
            WRITE(NFECRA,3000) II, ICAS(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C --- Repere locale
C
        DO II = 1, NNENT
          IF (ICAS(II).NE.4) THEN
            NORM1 = DIR1(1,II)**2.D0+DIR1(2,II)**2.D0+DIR1(3,II)**2.D0
            NORM2 = DIR1(1,II)**2.D0+DIR1(2,II)**2.D0+DIR1(3,II)**2.D0
            CROSP = (DIR1(2,II)*DIR2(3,II)-DIR1(3,II)*DIR2(2,II))**2.D0
     &            + (DIR1(1,II)*DIR2(3,II)-DIR1(3,II)*DIR2(1,II))**2.D0
     &            + (DIR1(1,II)*DIR2(2,II)-DIR1(2,II)*DIR2(1,II))**2.D0
            IF(ABS(NORM1*NORM2*CROSP).LT.EPZERO) THEN
              IOK = IOK + 1
              WRITE(NFECRA,3100)II,SQRT(NORM1),SQRT(NORM2),SQRT(CROSP)
            ENDIF
            IF(ABS(NORM1-1.D0).GT.EPZERO.OR.
     &         ABS(NORM2-1.D0).GT.EPZERO.AND.
     &         ABS(NORM1*NORM2*CROSP).GT.EPZERO) THEN
              WRITE(NFECRA,3150) II, SQRT(NORM1),SQRT(NORM2)
            ENDIF
          ENDIF
        ENDDO
C
C --- Conditions aux limites
C
        DO II = 11, NNENT
          IF(ICAS(II).EQ.1) THEN
            DO JJ = 1, 4
              IF(ICLVOR(JJ,II).NE.1
     &             .AND.ICLVOR(JJ,II).NE.2
     &             .AND.ICLVOR(JJ,II).NE.3) THEN
                WRITE(NFECRA,3200) ICLVOR(JJ,II), JJ, II
                IOK = IOK + 1
              ENDIF
            ENDDO
            IF((ICLVOR(1,II).EQ.3.AND.ICLVOR(3,II).NE.3).OR.
     &           (ICLVOR(1,II).NE.3.AND.ICLVOR(3,II).EQ.3).OR.
     &           (ICLVOR(2,II).EQ.3.AND.ICLVOR(4,II).NE.3).OR.
     &           (ICLVOR(2,II).NE.3.AND.ICLVOR(4,II).EQ.3)) THEN
              WRITE(NFECRA,3250) ICLVOR(JJ,II), JJ, II
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
C --- Dimension carateristiques
C
        DO II = 1, NNENT
          IF(ICAS(II).EQ.1.AND.(LLZ(II).LE.0.D0.OR.
     &         LLY(II).LE.0.D0)) THEN
            WRITE(NFECRA,3400) LLZ(II), LLY(II), II
            IOK = IOK + 1
          ENDIF
          IF(ICAS(II).EQ.2.AND.LLD(II).LE.0.D0) THEN
            WRITE(NFECRA,3500) LLD(II), II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C --- Parametres pour la duree de vie, la taille, et le deplacement des vortex
C
        DO II = 1, NNENT
          IF(ITLIVO(II).NE.1.AND.ITLIVO(II).NE.2) THEN
            WRITE(NFECRA,4000) ITLIVO(II), II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(ITLIVO(II).EQ.1.AND.TLIMVO(II).LE.0.D0) THEN
            WRITE(NFECRA,4100) TLIMVO(II), II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(ISGMVO(II).NE.1.AND.ISGMVO(II).NE.2.AND.
     &         ISGMVO(II).NE.3) THEN
            WRITE(NFECRA,4200) ISGMVO(II), II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(ISGMVO(II).EQ.1.AND.XSGMVO(II).LE.0.D0) THEN
            WRITE(NFECRA,4300) XSGMVO(II), II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(IDEPVO(II).NE.1.AND.IDEPVO(II).NE.2.AND.
     &         IDEPVO(II).NE.3) THEN
            WRITE(NFECRA,4400) IDEPVO(II), II
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(IDEPVO(II).EQ.1.AND.UD(II).LE.0.D0) THEN
            WRITE(NFECRA,4500) UD(II), II
C     IOK = IOK + 1
          ENDIF
        ENDDO
C
C --- Donnees utilisateur
C
        DO II = 1, NNENT
          IF(ICAS(II).EQ.4.AND.UDEBIT(II).LE.0.D0) THEN
            WRITE(NFECRA,5000) II, UDEBIT(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(ICAS(II).EQ.4.AND.KDEBIT(II).LE.0.D0) THEN
            WRITE(NFECRA,5100) II, KDEBIT(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO II = 1, NNENT
          IF(ICAS(II).EQ.4.AND.EDEBIT(II).LE.0.D0) THEN
            WRITE(NFECRA,5200) II, EDEBIT(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
      ENDIF
C=======================================================================
C 3. FORMATS
C=======================================================================
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    NNENT DOIT ETRE UN ENTIER POSITIF INFERIEUR A ',I10      ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort ou augmenter la taille de NENTMX dans     ',/,
     &'@  vortex.h                                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    NDAT DOIT ETRE UN ENTIER POSITIF INFERIEUR A ',I10       ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort ou agmenter la taille de NDAT dans        ',/,
     &'@  vortex.h                                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE NOMBRE DE VORTEX NVORT A L''ENTREE',I10               ,/,
     &'@    DOIT ETRE UN ENTIER POSITIF                             ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ISUIVO DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    IREPVO INDIQUE LE NUMERO DE L''ENTREE                   ',/,
     &'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A ',I10  ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ICAS INDIQUE LE CAS TRAITE                              ',/,
     &'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 4     ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LES VECTEURS DIR1 ET DIR2 NE PEUVENT ETRE NI NULS       ',/,
     &'@    NI COLINEAIRES                                          ',/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@    ||DIR1|| VAUT ',E14.5                                    ,/,
     &'@    ||DIR2|| VAUT ',E14.5                                    ,/,
     &'@    ||DIR1 ^ DIR2 || VAUT ',E14.5                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3150 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LES VECTEURS DIR1 ET DIR2 NE SONT PAS NORMES            ',/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@    ||DIR1|| VAUT ',E14.5                                    ,/,
     &'@    ||DIR2|| VAUT ',E14.5                                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ICLVOR INDIQUE LE TYPE DE CONDITION AUX LIMITES         ',/,
     &'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 4     ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    SUR LE COTE ',I10                                        ,/,
     &'@    DE L''ENTREE ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3250 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VERFIFIER LA COHERENCE DES CONDITION AUX LIMITES        ',/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@    NOMBRE IMPAIR DE CONDITIONS PERIODIQUES                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LLY ET LLZ SONT LES DIMENSIONS DE L''ENTREE             ',/,
     &'@    ELLES VALENT ICI ',E14.5, ' ET ',E14.5                   ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LLD EST LE DIAMETRE DE LA CONDUITE                      ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ITLIVO INDIQUE LE TYPE DE MODELE UTILISE POUR CALCULER  ',/,
     &'@    LA DUREE DE VIE DES VORTEX                              ',/,
     &'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 2     ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    TLIMVO EST LA DUREE DE VIE MAXIMALE DES VORTEX          ',/,
     &'@    ELLE VAUT ICI ',E14.5                                    ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ISGMVO INDIQUE LE TYPE DE MODELE UTILISE POUR CALCULER  ',/,
     &'@    LA TAILLE DES VORTEX                                    ',/,
     &'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 3     ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    XSGMVO EST LA TAILLE MAXIMALE DES VORTEX                ',/,
     &'@    ELLE VAUT ICI ',E14.5                                    ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    IDEPVO INDIQUE LE TYPE DE MODELE UTILISE POUR CALCULER  ',/,
     &'@    LE DEPLACEMENT DES VORTEX                               ',/,
     &'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 3     ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@    UD EST LA VITESSE MAXIMALE DE DEPLCAMENT DES VORTEX     ',/,
     &'@    ELLE VAUT ICI ',E14.5                                    ,/,
     &'@    A L''ENTREE ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA VITESSE DEBITANTE UDEBIT A L''ENTREE ',I10            ,/,
     &'@    VAUT ',E14.5                                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L ENERGIE CINETIQUE KDEBIT A L''ENTREE ',I10             ,/,
     &'@    VAUT ',E14.5                                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA DISSIPATION EDEBIT A L''ENTREE ',I10                  ,/,
     &'@    VAUT ',E14.5                                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    NNENT MUST BE A POSITIVE INTEGER LOWER THAN',I10         ,/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort or increase the size of NENTMX in vortex.h  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    NDAT MUST BE A POSITIVE INTEGER LOWER THAN ',I10         ,/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort or increase the size of NDAT in vortex.h    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE NUMBER OF VORTICES NVORT AT THE INLET ',I10          ,/,
     &'@    MUST BE A POSITIVE INTEGER                              ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    ISUIVO MUST BE AN INTEGER EQUAL TO 0 OR 1               ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    IREPVO GIVES THE INLET NUMBER                           ',/,
     &'@    IT IS A POSITIVE INTEGER LOWER THAN ',I10                ,/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    ICAS GIVES THE CHOSEN CASE                              ',/,
     &'@    IT IS A POSITIVE INTEGER LOWER THAN 4                   ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE VECTORS DIR1 AND DIR2 CANNOT BE NULL NOR COLINEAR   ',/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@    ||DIR1|| IS ',E14.5                                      ,/,
     &'@    ||DIR2|| IS ',E14.5                                      ,/,
     &'@    ||DIR1 ^ DIR2 || IS ',E14.5                              ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3150 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE VECTORS DIR1 AND DIR2 ARE NOT NORMED                ',/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@    ||DIR1|| IS ',E14.5                                      ,/,
     &'@    ||DIR2|| IS ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  The calculation will be run.                              ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    ICLVOR GIVES THE TYPE OF THE BOUNDARY CONDITIONS        ',/,
     &'@    IT IS A POSITIVE INTEGER LOWER THAN 4                   ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@    ON THE SIDE ',I10                                        ,/,
     &'@    OF THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3250 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    VERIFY THE COHERENCY OF THE BOUNDARY CONDITIONS         ',/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@    ODD NUMBER OF PERIODIC CONDITIONS                       ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    LLY AND LLZ ARE THE INLET DIMENSIONS                    ',/,
     &'@    THEIR VALUES ARE ',E14.5, ' AND ',E14.5                  ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    LLD IS THE PIPE DIAMETER                                ',/,
     &'@    ITS VALUE IS ',E14.5                                     ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    ITLIVO GIVES THE MODEL TYPE USED TO COMPUTE THE         ',/,
     &'@    VORTICES LIFE TIME                                      ',/,
     &'@    IT IS A POSITIVE INTEGER LOWER THAN OR EQUAL TO 2       ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    TLIMVO IS THE MAXIMUM LIFE TIME OF THE VORTICES         ',/,
     &'@    ITS VALUE IS ',E14.5                                     ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    ISGMVO GIVES THE MODEL TYPE USED TO COMPUTE THE         ',/,
     &'@    VORTICES SIZE                                           ',/,
     &'@    IT IS A POSITIVE INTEGER LOWER THAN OR EQUAL TO 3       ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    XSGMVO IS THE MAXIMUM SIZE OF THE VORTICES              ',/,
     &'@    ITS VALUE IS ',E14.5                                     ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    IDEPVO GIVES THE MODEL TYPE USED TO COMPUTE THE         ',/,
     &'@    VORTICES DISPLACEMENT                                   ',/,
     &'@    IT IS A POSITIVE INTEGER LOWER THAN OR EQUAL TO 3       ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    UD IS THE MAXIMUM DISPLACEMENT VELOCITY OF THE VORTICES ',/,
     &'@    ITS VALUE IS ',E14.5                                     ,/,
     &'@    AT THE INLET ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will be run.                              ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE INLET VELOCITY UDEBIT AT THE INLET ',I10             ,/,
     &'@    IS ',E14.5                                               ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE INLET KINETIC ENERGY KDEBIT AT THE INLET ',I10       ,/,
     &'@    IS ',E14.5                                               ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE INLET DISSIPATION EDEBIT AT THE INLET ',I10          ,/,
     &'@    IS ',E14.5                                               ,/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
#endif
C
C=======================================================================
C 3. SORTIE ET IMPRESSIONS FINALES
C=======================================================================
C
      IF(IOK.GT.0) THEN
        WRITE(NFECRA,9999)IOK
        CALL CSEXIT (1)
      ELSE
        WRITE(NFECRA,9998)
      ENDIF
C
#if defined(_CS_LANG_FR)
C
 9998 FORMAT(
     &'                                                             ',/,
     &' Pas d''erreur detectee lors de la verification des donnees  ',/,
     &'    pour la methode des vortex (usvort).                     ',/)
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
     &'@  Verifier usvort.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 9998 FORMAT(
     &'                                                             ',/,
     &' No error detected during the verification of the parameters ',/,
     &'    for the vortex method (usvort).                          ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLETE ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run (',I10,' errors).         ',/,
     &'@                                                            ',/,
     &'@  Refer to the previous warnings for further information.   ',/,
     &'@  Verify usvort.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
      RETURN
      END
c@z
