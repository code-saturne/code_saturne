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
                        SUBROUTINE RAYOPT
C                       *****************
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC   1) Initialisation par defaut du parametrage du module de
CFONC      transferts thermiques radiatifs
CFONC   2) Lecture du parametrage utilisateur
CFONC   3) Controle de coherence avec les physiques particulieres
CFONC   4) Verifications du parametrage utilisateur
CFONC
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "radiat.h"
      INCLUDE "ihmpre.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, IPHAS, IPHAOK, IOK , IISCAL, ISCAOK
      CHARACTER        CAR4*4
C
C***********************************************************************
C=======================================================================
C 1. INITIALISATIONS PAR DEFAUT DU MODULE DE TRANSFERTS RADIATIFS
C                        ^^^^^^
C=======================================================================
C
C
C-->  IIRAYO = 0 : PAS DE TRANSFERTS RADIATIFS
C            = 1 : TRANSFERTS RADIATIFS
C
      IIRAYO = 0
C
C-->  POUR CHAQUE PHASE:
C     IRAYON = 0 : PAS DE TRANSFERTS RADIATIFS
C            = 1 : TRANSFERTS RADIATIFS, METHODE DES ORDONNEES DISCRETES
C            = 2 : TRANSFERTS RADIATIFS, APPROXIMATION P-1
C     On initialise a -1 pour montrer que ce n'est pas initialise ...
C        (on fera un test apres usray1)
      DO IPHAS = 1,NPHSMX
        IRAYON(IPHAS) = -1
      ENDDO
C
C
C-->  CALCUL DU COEFFICIENT D'ABSORPTION
C     (REMPLI AUTOMATIQUEMENT ENSUITE POUR LA PHYSIQUE PARTICULIERE)
C      IMODAK = 0 : sans utiliser modak
C               1 : a l'aide modak
C
      IMODAK = 0
C
C-->  INDICATEUR SUITE DE CALCUL (LECTURE DU FICHIER SUITE)
C
      ISUIRD = -1
C
C-->  FREQUENCE DE PASSAGE DANS LE MODULE RAYONNEMENT
C
      NFREQR = -1
C
C-->  NOMBRE DE DIRECTIONS : 32 OU 128
C
      NDIREC = -1
C
C-->  NOMBRE DE BANDES SPECTRALES (PAS UTILISE)
C
      NBANDE = 1
C
C-->  POURCENTAGE DE CELLULES OU L'ON ADMET QUE LA LONGUEUR OPTIQUE DEPASSE
C       L'UNITE POUR LE MODELE P-1
C
      XNP1MX = 10.D0
C
C-->  INITIALISATION DU MODE DE CALCUL DU TERME SOURCE RADIATIF EXPLICITE
C     IDIVER = 0 => CALCUL SEMI-ANALYTIQUE (OBLIGATOIRE SI TRANSPARENT)
C     IDIVER = 1 => CALCUL CONSERVATIF
C     IDIVER = 2 => CALCUL SEMI-ANALYTIQUE CORRIGE POUR ETRE CONSERVATIF
C     REMARQUE : SI TRANSPARENT IDIVER = -1 AUTOMATIQUEMENT DANS RAYDOM
C
      IDIVER = -1
C
C--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS TEMPERATURE DE PAROI
C
      IIMPAR = -1
C
C--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS RESOLUTION LUMINANCE
C
      IIMLUM = -1
C
C--> INITIALISATION DES DONNEES POST-PROCESSING
C    NBRAYP : NOMBRE MAX DES SORTIES DE VARIABLES CELLULES
C    NBRAYF : NOMBRE MAX DES SORTIES DE VARIABLES FACETTES DE BORD
C
      DO IPHAS = 1,NPHSMX
        DO II = 1,NBRAYP
          WRITE(CAR4,'(I4.4)') II + NBRAYP * (IPHAS -1)
          NBRVAP(II,IPHAS) = 'RAYTCE'//CAR4
          IRAYVP(II,IPHAS) = -1
        ENDDO
        DO II = 1,NBRAYF
          WRITE(CAR4,'(I4.4)') II+ NBRAYF * (IPHAS -1)
          NBRVAF(II,IPHAS) = 'RAYTFB'//CAR4
          IRAYVF(II,IPHAS) = -1
        ENDDO
      ENDDO
C
C=======================================================================
C 2. INITIALISATIONS UTILISATEURS
C                    ^^^^^^^^^^^^
C=======================================================================
C
C   - Interface Code_Saturne
C     ======================
C
      IF (IIHMPR.EQ.1) THEN
C
        DO IPHAS = 1, NPHAS
          DO JJ = 1,NBRAYP
            II = JJ + NBRAYP * (IPHAS -1)
            CALL FCNMRA(NBRVAP(JJ,IPHAS), LEN(NBRVAP(JJ,IPHAS)), II)
          ENDDO
        ENDDO
        DO IPHAS = 1, NPHAS
          DO JJ = 1, NBRAYF
            II = JJ + NPHAS * NBRAYP + NBRAYF* (IPHAS -1)
            CALL FCNMRA(NBRVAF(JJ,IPHAS), LEN(NBRVAF(JJ,IPHAS)), II)
          ENDDO
        ENDDO
C
        CALL UIRAY1
C       ===========
     &     ( NBRAYP, NBRAYF, NPHAS,
     &       IRAYON, ISUIRD, NDIREC, NFREQR, IDIVER, IIMPAR, IIMLUM,
     &       IRAYVP, IRAYVF )
C
        DO IPHAS = 1, NPHAS
          DO JJ = 1,NBRAYP
            II = JJ + NBRAYP * (IPHAS -1)
            CALL CFNMRA(NBRVAP(JJ,IPHAS), LEN(NBRVAP(JJ,IPHAS)), II)
          ENDDO
        ENDDO
        DO IPHAS = 1, NPHAS
          DO JJ = 1, NBRAYF
            II = JJ + NPHAS * NBRAYP + NBRAYF* (IPHAS -1)
            CALL CFNMRA(NBRVAF(JJ,IPHAS), LEN(NBRVAF(JJ,IPHAS)), II)
          ENDDO
        ENDDO
C
      ENDIF
C
      CALL USRAY1
C     ===========
C
C=======================================================================
C 3. VERIFICATION LA COHERENCE D'UTILISATION DU MODULE DE RAYONNEMENT
C    AVEC LA THERMIQUE OU LES PHYSIQUES PARTICULIERES (COMBUSTION)
C=======================================================================
C
      IOK = 0
C
C--> IRAYON
C
C     Si physique classique ou particuliere sans fichier paramétrique
C       et IRAYON pas modifie dans USRAY1
C       IRAYON(IPHAS) = 0 (pas de rayonnement).
C
      IF( IPPMOD(IPHPAR).LE.1 ) THEN
C
        DO IPHAS = 1, NPHAS
          IF(IRAYON(IPHAS).EQ.-1) THEN
            IRAYON(IPHAS) = 0
          ENDIF
        ENDDO
C
C
C     Sinon, IRAYON doit etre issu du fichier
C       parametrique ; IMODAK aussi.
C
      ELSE
C
        IPHAS = 1
C
        IF (IRAYON(IPHAS).EQ.-1) THEN
          IF (IRAYPP.EQ.0) THEN
            IRAYON(IPHAS) = 0
            IMODAK        = 0
          ELSEIF (IRAYPP.EQ.1) THEN
            IRAYON(IPHAS) = 1
            IMODAK        = 0
          ELSEIF (IRAYPP.EQ.2) THEN
            IRAYON(IPHAS) = 1
            IMODAK        = 1
          ELSEIF (IRAYPP.EQ.3) THEN
            IRAYON(IPHAS) = 2
            IMODAK        = 0
          ELSEIF (IRAYPP.EQ.4) THEN
            IRAYON(IPHAS) = 2
            IMODAK        = 1
          ELSE
            WRITE(NFECRA,1010)
            IOK = IOK + 1
          ENDIF
C
        ELSE
          WRITE(NFECRA,2000)
          IOK = IOK + 1
        ENDIF
C
      ENDIF
C
C--> ISCSTH
C
C     Si physique particuliere avec fichier parametrique,
C       ISCSTH a ete renseigne
C       dans ppini1 ou coini1 (a verifier en elec).
C     Si physique classique et ISCSTH pas modifie dans USRAY1
C                           ou pas de variable thermique
C       STOP.
C
      IF(IPPMOD(IPHPAR).LE.1) THEN
C
C     Il y a une seule phase ; si on rayonne
        DO IPHAS = 1, NPHAS
          IF (IRAYON(IPHAS).EQ.1 .OR. IRAYON(IPHAS).EQ.2) THEN
C        On cherche s'il y a un scalaire thermique
            ISCAOK = 0
            DO IISCAL = 1, NSCAL
              IF (IISCAL.EQ.ISCALT(IPHAS)) THEN
                ISCAOK = 1
C           Et on regarde si on a dit enthalpie
                IF (ISCSTH(IISCAL).NE.2) THEN
                  WRITE(NFECRA,3000) IPHAS,IISCAL,IISCAL
                  IOK = IOK + 1
                ENDIF
              ENDIF
            ENDDO
            IF (ISCAOK.EQ.0) THEN
              WRITE(NFECRA,3001) IPHAS,IPHAS
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
      ELSE
C
C     Pour toutes les phases qui rayonnent
        DO IPHAS = 1, NPHAS
          IF (IRAYON(IPHAS).EQ.1 .OR. IRAYON(IPHAS).EQ.2) THEN
C
C        On cherche s'il y a un scalaire thermique
            ISCAOK = 0
            DO IISCAL = 1, NSCAL
              IF (IISCAL.EQ.ISCALT(IPHAS)) THEN
                ISCAOK = 1
C
C           Et on regarde si on a dit temp C, K ou enthalpie
                IF (ABS(ISCSTH(IISCAL)).NE.1.AND.
     &                  ISCSTH(IISCAL) .NE.2      ) THEN
                  WRITE(NFECRA,3010) IPHAS,IISCAL,IISCAL
                  IOK = IOK + 1
                ENDIF
C
              ENDIF
            ENDDO
            IF(ISCAOK.EQ.0)THEN
              WRITE(NFECRA,3011)IPHAS,IPHAS
              IOK = IOK + 1
            ENDIF
C
          ENDIF
        ENDDO
C
      ENDIF
C
C--> Stop si erreur.
C
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    IRAYPP NE PEUT PRENDRE POUR VALEURS QUE 0 1 2 3 OU 4    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Arret dans rayopt.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    EN PHYSIQUE PARTICULIERE IRAYON EST DONNE PAR FICHIER   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1 (IRAYON ne doit pas etre modifie).        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE ACTIVEE : ENTHALPIE NECESSAIRE    ',/,
     &'@                                                            ',/,
     &'@  Avec rayonnement, pour la phase ',I10   ,', il faut       ',/,
     &'@    preciser la variable energetique representee par le     ',/,
     &'@    scalaire ',I10   ,' en renseignant ISCSTH(',I10   ,')   ',/,
     &'@    dans usini1 : soit                                      ',/,
     &'@               -1 temperature en C                          ',/,
     &'@                1 temperature en K                          ',/,
     &'@                2 enthalpie                                 ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Avec physique particuliere, cette initialisation aurait   ',/,
     &'@    du etre automatique.                           ~~~~~~   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE ACTIVEE : ENTHALPIE NECESSAIRE    ',/,
     &'@                                                            ',/,
     &'@  Lorsque le rayonnement est utilise (phase ',I10   ,', il  ',/,
     &'@    faut indiquer qu''un scalaire represente la variable    ',/,
     &'@    energetique (enthalpie) en renseignant                  ',/,
     &'@    ISCALT(',I10   ,') dans usini1 (numero du scalaire).    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Avec physique particuliere, cette initialisation aurait   ',/,
     &'@    du etre automatique.                           ~~~~~~   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    ISCSTH DOIT ETRE RENSEIGNE OBLIGATOIREMENT DANS USINI1  ',/,
     &'@                                                            ',/,
     &'@  Avec rayonnement, pour la phase ',I10   ,', il faut       ',/,
     &'@    preciser la variable energetique representee par le     ',/,
     &'@    scalaire ',I10   ,' en renseignant ISCSTH(',I10   ,')   ',/,
     &'@    dans usini1 : soit                                      ',/,
     &'@               -1 temperature en C                          ',/,
     &'@                1 temperature en K                          ',/,
     &'@                2 enthalpie                                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    IL FAUT UTILISER UNE VARIABLE ENERGETIQUE.              ',/,
     &'@                                                            ',/,
     &'@  Lorsque le rayonnement est utilise (phase ',I10   ,', il  ',/,
     &'@    faut indiquer qu''un scalaire represente la variable    ',/,
     &'@    energetique (temperature ou enthalpie) en renseignant   ',/,
     &'@    ISCALT(',I10   ,') dans usini1 (numero du scalaire).    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C 4. FAIT-ON DU RAYONNEMENT SEMI-TRANSPARENT ?
C=======================================================================
C
      IIRAYO = 0
      DO IPHAS = 1,NPHAS
        IF (IRAYON(IPHAS).EQ.1 .OR. IRAYON(IPHAS).EQ.2) IIRAYO = 1
      ENDDO
C
C=======================================================================
C 5. VERIFICATIONS (uniquement s'il y a du rayonnement)
C=======================================================================
C
      IOK = 0
C
      IF (IIRAYO.EQ.1) THEN
C
C --> ISUIRD
C
        IF (ISUIRD.NE.0 .AND. ISUIRD.NE.1 ) THEN
          WRITE(NFECRA,4000)ISUIRD
          IOK = IOK + 1
        ENDIF
C
C --> NFREQR
C
        IF (NFREQR.LE.0) THEN
          WRITE(NFECRA,4010)NFREQR
          IOK = IOK + 1
        ENDIF
C
C --> NDIREC
C     Choix entre 32 et 128 directions (cf raysol)
        IF (NDIREC.NE.32 .AND. NDIREC.NE.128 ) THEN
          WRITE(NFECRA,4020) NDIREC
          IOK = IOK + 1
        ENDIF
C
C --> IDIVER
C     Choix entre 0  1 et 2
        IF (IDIVER.NE.0 .AND. IDIVER.NE.1 .AND. IDIVER.NE.2) THEN
          WRITE(NFECRA,4030) IDIVER
          IOK = IOK + 1
        ENDIF
C
C --> IIMPAR
C     Choix entre 0  1 et 2
        IF (IIMPAR.NE.0 .AND. IIMPAR.NE.1 .AND. IIMPAR.NE.2) THEN
          WRITE(NFECRA,4040) IIMPAR
          IOK = IOK + 1
        ENDIF
C
C --> IIMLUM
C     Choix entre 0  1 et 2
        IF (IIMLUM.NE.0 .AND. IIMLUM.NE.1 .AND. IIMLUM.NE.2) THEN
          WRITE(NFECRA,4050) IIMLUM
          IOK = IOK + 1
        ENDIF
C
C --> IRAYVP & IRAYVF
C     Choix entre -1 et 1
        DO IPHAS = 1,NPHAS
          IF (IRAYON(IPHAS).EQ.1 .OR. IRAYON(IPHAS).EQ.2) THEN
            DO II = 1,NBRAYP
              IF (IRAYVP(II,IPHAS).NE.1 .AND.
     &            IRAYVP(II,IPHAS).NE.-1) THEN
                WRITE(NFECRA,4060) NBRVAP(II,IPHAS), IRAYVP(II,IPHAS)
                IOK = IOK + 1
              ENDIF
            ENDDO
            DO II = 1,NBRAYF
              IF (IRAYVF(II,IPHAS).NE.1 .AND.
     &            IRAYVF(II,IPHAS).NE.-1) THEN
                WRITE(NFECRA,4070) NBRVAF(II,IPHAS), IRAYVF(II,IPHAS)
                IOK = IOK + 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C
C --> CONNECTIVITE
        IPHAOK = 0
        DO IPHAS = 1,NPHAS
          DO II = 1,NBRAYF
            IF (IRAYVF(II,IPHAS).EQ.1) IPHAOK = IPHAOK + 1
          ENDDO
        ENDDO
        IF (IPHAOK.NE.0) THEN
          IF (NNOD.EQ.0 .OR. (LNDFAC.EQ.0 .AND. LNDFBR.EQ.0)) THEN
            WRITE(NFECRA,4080) NNOD, LNDFAC, LNDFBR
            IOK = IOK + 1
          ENDIF
        ENDIF
C
      ENDIF
C
C--> Stop si erreur.
C
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    INDICATEUR DE SUITE DE CALCUL NON ADMISSIBLE            ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de suite de calcul doit etre 0 ou 1 (ISUIRD)',/,
     &'@    Il vaut ici ISUIRD = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@     FREQUENCE DE PASSAGE DANS LE MODULE DE RAYONNEMENT     ',/,
     &'@     NON ADMISSIBLE                                         ',/,
     &'@                                                            ',/,
     &'@  La frequence de passage doit etre superieure ou egale a 1 ',/,
     &'@    Elle vaut ici NFREQR = ',I10                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    NOMBRE DE DIRECTIONS NON ADMISSIBLE                     ',/,
     &'@                                                            ',/,
     &'@  Le nombre de directions doit etre 32 ou 128 (NDIREC)      ',/,
     &'@    Il vaut ici NDIREC = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,
     &'@    *********                                               ',/,
     &'@    INDICATEUR DU MODE DE CALCUL DU TERME SOURCE RADIATIF   ',/,
     &'@    EXPLICITE NON ADMISSIBLE                                ',/,
     &'@                                                            ',/,
     &'@  L''indicateur du mode de calcul doit etre 0, 1 ou 2       ',/,
     &'@    Il vaut ici IDIVER = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,
     &'@    *********                                               ',/,
     &'@    NIVEAU D''AFFICHAGE DES RENSIGNEMENTS DES               ',/,
     &'@    TEMPERATURE DE PAROI NON ADMISSIBLE                     ',/,
     &'@                                                            ',/,
     &'@  Le niveau d''affichage doit etre 0, 1 ou 2  (IIMPAR)      ',/,
     &'@    Il vaut ici IIMPAR = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,
     &'@    *********                                               ',/,
     &'@    NIVEAU D''AFFICHAGE DES RENSIGNEMENTS SUR LA            ',/,
     &'@    RESOLUTION DE LA LUMINANCE NON ADMISSIBLE               ',/,
     &'@                                                            ',/,
     &'@  Le niveau d''affichage doit etre 0, 1 ou 2  (IIMLUM)      ',/,
     &'@    Il vaut ici IIMLUM = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4060 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,
     &'@    *********                                               ',/,
     &'@    INDICATEUR DE SORTIE EN POSTPROCESSING                  ',/,
     &'@    POUR ',A40                                               ,/,
     &'@    NON ADMISSIBLE                                          ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de postprocessing doit etre -1 ou 1         ',/,
     &'@    Il vaut ici IRAYVP = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4070 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,
     &'@    *********                                               ',/,
     &'@    INDICATEUR DE SORTIE EN POSTPROCESSING                  ',/,
     &'@    POUR ',A40                                               ,/,
     &'@    NON ADMISSIBLE                                          ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de postprocessing doit etre -1 ou 1         ',/,
     &'@    Il vaut ici IRAYVF = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage.                                    ',/,
     &'@                                                            ',/,
     &'@  On conseille donc de verifier usray1.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4080 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,
     &'@    *********                                               ',/,
     &'@    LES CONNECTIVITES FACES -> SOMMETS OU LES COORDONNEES   ',/,
     &'@    DES NOEUDS SOMMETS NE SONT PAS FOURNIS                  ',/,
     &'@                                                            ',/,
     &'@  Le module rayonnement est active mais les noeuds sommets  ',/,
     &'@   ou les connectivites faces -> sommets n''ont pas ete     ',/,
     &'@   fournis au noyau.                                        ',/,
     &'@                                                            ',/,
     &'@    NNDOD  = ', I10                                          ,/,
     &'@    LNDFAC = ', I10                                          ,/,
     &'@    LNDFBR = ', I10                                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier si le fichier geomet contient bien les           ',/,
     &'@    connectivites faces -> sommets si IFOENV = 0.           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C=======================================================================
C 6. NOMBRE DE PHASES POUR LESQUELLES ON FAIT DU RAYONNEMENT
C=======================================================================
C
C--> NPHAST : NOMBRE DE PHASES POUR LESQUELLES ON FAIT DU RAYONNEMENT
C--> IRAPHA : NUMEROS DES PHASES POUR LESQUELLES ON FAIT DU RAYONNEMENT
C
C--> NPHASC : EN CP, NOMBRE DE PHASE (= 1 FORCEMENT) + NOMBRE DE CLASSES
C             SERT A STOCKER LE COEFF D'ABSORPTION PAR EXEMPLE POUR LE
C             MELANGE ET LES CLASSES DE PARTICULES (VOIR MEMRA1)
C
      NPHAST = 0
      DO IPHAS = 1,NPHAS
        IF (IRAYON(IPHAS).GE.1) THEN
          NPHAST = NPHAST +1
          IRAPHA(NPHAST) = IPHAS
        ENDIF
      ENDDO
C
      NPHASC = NPHAST
C
C--> Pour le Charbon Pulverise :
      IF ( IPPMOD(ICP3PL) .GE. 0 ) THEN
        NPHASC = NPHASC + NCLACP
      ENDIF
C
C--> Pour le Fuel :
      IF ( IPPMOD(ICFUEL) .GE. 0 ) THEN
        NPHASC = NPHASC + NCLAFU
      ENDIF
C
      RETURN
C
      END
c@z
