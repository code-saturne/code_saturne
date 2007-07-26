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
                        SUBROUTINE ELLECD
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
CFONC             RELATIF AU MODULE ELECTRIQUE
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
C VARIABLES LOCALES
C
      INTEGER          IT, IOS , IESP, II , I
C
C***********************************************************************
C=======================================================================
C 0. INITIALISATION
C=======================================================================
C
      NGAZG = 1
C
C=======================================================================
C 1. LECTURE DU FICHIER DONNEES SPECIFIQUES
C=======================================================================
C
C --> Ouverture du fichier
C
      IF ( IPPMOD(IELARC).GE.1 .OR. IPPMOD(IELJOU).EQ.3
     &                         .OR. IPPMOD(IELJOU).EQ.4 ) THEN
        OPEN ( UNIT=IMPFPP, FILE=FICFPP,
     &        STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &                                        IOSTAT=IOS, ERR=99 )
        REWIND ( UNIT=IMPFPP,ERR=99 )
      ENDIF
C
C==========================================
C 2. LECTURE D'UN FICHIER ARC ELECTRIQUE
C==========================================
C     Il n'y a pas de lecture de fichier prevue pour l'effet Joule
C
C
C
C  on transporte l'enthalpie, on suppose qu'il y a du rayonnement
C  la conductivite, la masse volumique, l'emissivite et
C  la viscosite dependent de la temperature
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
        NGAZG = 0
C
C ----- NB de constituants et Nb de points de tabulation
C     on saute les lignes de commentaires
        DO II = 1, 7
          READ (IMPFPP,*)
        ENDDO
        READ ( IMPFPP,*,ERR=999,END=999 ) NGAZG,NPO
C
        IF ( NPO.GT.NPOT ) THEN
          WRITE(NFECRA,8000) NPOT
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        IF ( NGAZG.GT. NGAZGM) THEN
          WRITE(NFECRA,8001) NGAZGM, NGAZG
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        IF ( NGAZG.LT. 1) THEN
          WRITE(NFECRA,8002)NGAZG
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C -----  Lecture de l'indicateur pour savoir ce que represente XKABEL
C      on saute les lignes de commentaires
        DO II = 1, 5
          READ (IMPFPP,*)
        ENDDO
        READ ( IMPFPP,*,ERR=999,END=999 ) IXKABE
        IF ( IXKABE .LT. 0 .OR. IXKABE .GE. 3 ) THEN
          WRITE(NFECRA,8003) IXKABE
          CALL CSEXIT (1)
C         ===========
        ENDIF
C     Si IXKABE=1, il faut activer explicitement le rayonnement
        IF (IXKABE.EQ.1) IRAYPP = 1
C
C ----- En fonction de la temperature pour chaque espece courante
C          Enthalpie massique
C          Masse volumique
C          Chaleur massique
C          Conductivite electrique
C          Viscosite laminaire
C          Conductivite Thermique
C          Coefficent d'absorption (rayonnement)
C     on saute les lignes de commentaires au debut
C
        IF(NGAZG.GT.0.AND.NPO.GT.0) THEN
          DO II = 1, 7
            READ (IMPFPP,*)
          ENDDO
          DO IESP = 1, NGAZG
            DO IT = 1, NPO
              READ (IMPFPP,*,ERR=999,END=999 )
     &             TH(IT)         ,EHGAZG(IESP,IT),RHOEL(IESP,IT),
     &             CPEL(IESP,IT)  ,SIGEL(IESP,IT) ,
     &             VISEL(IESP,IT) ,XLABEL(IESP,IT),
     &             XKABEL(IESP,IT)
            ENDDO
          ENDDO
        ENDIF
C
      ENDIF
C
C==========================================
C 3. LECTURE D'UN FICHIER EFFET JOULE
C==========================================
C
      IF ( IPPMOD(IELJOU).EQ. 3 .OR. IPPMOD(IELJOU).EQ. 4 ) THEN
C
C ----- Lecture du transfo de reference
C
        READ (IMPFPP,*,ERR=999,END=999 ) NTFREF
C
C ----- Nombre de transfo
C
C       on saute 2 lignes de commentaires
        READ (IMPFPP,*)
        READ (IMPFPP,*)
        READ (IMPFPP,*,ERR=999,END=999 ) NBTRF
C
C     Boucle sur le nombre de transfo
C
        DO I=1,NBTRF
C
C         on saute la ligne de commentaire
          READ (IMPFPP,*)
C
C         Tension primaire
          READ (IMPFPP,*,ERR=999,END=999 ) TENSPR(I)
C
C         Rapport du nombre de spires
          READ (IMPFPP,*,ERR=999,END=999 ) RNBS(I)
C
C         Impedances complexes
          READ (IMPFPP,*,ERR=999,END=999 ) ZR(I),ZI(I)
C
C         Type de branchement primaire
          READ (IMPFPP,*,ERR=999,END=999 ) IBRPR(I)
C
C         Type de branchement secondaire
          READ (IMPFPP,*,ERR=999,END=999 ) IBRSEC(I)
C
        ENDDO
C
C ----- Nombre d'electrodes
C
        READ (IMPFPP,*)
        READ (IMPFPP,*)
        READ (IMPFPP,*,ERR=999,END=999 ) NBELEC
C
C       Boucle sur le nombre d'electrodes
C
        DO I=1,NBELEC
C
C         Tension primaire
          READ (IMPFPP,*,ERR=999,END=999 ) IELECC(I),IELECT(I),
     &                                     IELECB(I)
C
        ENDDO
C
      ENDIF
C
C==========================================
C 4. LECTURE D'UN FICHIER MIGRATION IONIQUE
C==========================================
C
C    c'est plus complique (chaque espece peut avoir sa propre mobilite)
C    mais il n'y a pas forcement de rayonnement, ni meme de
C    chauffage significatif par effet Joule
C
C======================================
C 5. CALCULS DE DONNEES COMPLEMENTAIRES
C======================================
C
C
      RETURN
C
C
C====================
C 3. SORTIE EN ERREUR
C====================
C
  99  CONTINUE
      WRITE ( NFECRA,9998 )
      CALL CSEXIT (1)
C     ===========
C
  999 CONTINUE
      WRITE ( NFECRA,9999 )
      CALL CSEXIT (1)
C     ===========
C
C
C--------
C FORMATS
C--------
C
C
 8000 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,
     &'@                                                            ',/,
     &'@  Le nombre de points de tabulation lu dans le fichier de   ',/,
     &'@    donnees doit etre un entier inferieur ou egal           ',/,
     &'@              a NPOT   = ',I10                               ,/,
     &'@    Il vaut ici NPO    = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes courantes lu dans le fichier de      ',/,
     &'@    doit etre un entier inferieur ou egal                   ',/,
     &'@              a NGAZGM = ',I10                               ,/,
     &'@    Il vaut ici NGAZG  = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8002 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes courantes lu dans le fichier de      ',/,
     &'@    doit etre un entier superieur ou egal a 1.              ',/,
     &'@    Il vaut ici NGAZG = ',I10                                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8003 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,
     &'@                                                            ',/,
     &'@  La valeur de l''indicateur pour le rayonnement            ',/,
     &'@    doit etre comprise entre 0 et 2                         ',/,
     &'@    elle vaut ici IXKABE = ',I10                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9998 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,
     &'@                                                            ',/,
     &'@  Erreur a l''ouverture du fichier parametrique.            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,
     &'@                                                            ',/,
     &'@  Erreur a la lecture du fichier parametrique.              ',/,
     &'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,
     &'@    ou son format inadapte.                                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      END
C
C
c@z
