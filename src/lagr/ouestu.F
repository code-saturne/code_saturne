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
                        SUBROUTINE OUESTU
C                       *****************
C     -------------------------------------------------------------
     & ( NFECRA , NDIM   , NPOINT ,
     &   IERROR ,
     &   PX     , PY     , PZ     ,
     &   QX     , QY     , QZ     ,
     &   CDGX   , CDGY   , CDGZ   ,
     &   CELX   , CELY   , CELZ   ,
     &   ITYPF7 , ICONF7 , XYZNOD ,
     &   INDIAN )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC   reperage d'une particule par rapport a une face :
CFONC   soient P le point de depart de la particule
CFONC   et Q le point d'arrivee ; on obtient :
CFONC
CFONC   INDIAN =  0 le rayon PQ ne sort pas de la cellule par cette face
CFONC   INDIAN = -1 meme cellule
CFONC   INDIAN =  1 sortie de la cellule par cette face
CFONC
CFONC   ATTENTION : POUR UTILISER CES SUBROUTINE DE REPERAGE, IL EST
CFONC   NECESSAIRE QUE LE TYPE DOUBLE PRECISION SOIT CODE SUR 8 OCTETS,
CFONC   i.e. DOUBLE PRECISION = REAL*8
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER DE SORTIE LISTING   !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE (=3)           !
CARGU ! NPOINT       ! E  !  -> ! NOMBRE DE NOEUDS                     !
CARGU ! IERROR       ! E  ! <-  ! INDICATEUR D'ERREUR                  !
CARGU ! PX,PY,PZ     ! R  !  -> ! POINT DE DEPART DE LA PARTICULE      !
CARGU ! QX,QY,QZ     ! R  !  -> ! POINT DE D'ARRIVE DE LA PARTICULE    !
CARGU ! CDGX,..,CDGZ ! R  !  -> ! CENTRE DE GRAVITE DE LA FACE         !
CARGU ! CELX,..,CELZ ! R  !  -> ! CENTRE DE GRAVITE DE LA CELLULE      !
CARGU !              !    !     !   CONTENANT LE POINT DE DEPART       !
CARGU ! ITYPF7       ! E  !  -> ! NOMBRE DE POINTS SUPPORT DE LA FACE  !
CARGU ! ICONF7(ITYPF7! E  !  -> ! NUMERO DES POINTS SUPPORT            !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! INDIAN       ! E  ! <-  ! INDICATEUR D'ORIENTATION             !
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
C
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "cstnum.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NFECRA , NDIM , NPOINT
      INTEGER          INDIAN , IERROR
      INTEGER          ITYPF7 , ICONF7(ITYPF7)
C
      DOUBLE PRECISION PX , PY , PZ
      DOUBLE PRECISION QX , QY , QZ
      DOUBLE PRECISION CDGX , CDGY , CDGZ
      DOUBLE PRECISION CELX , CELY , CELZ
      DOUBLE PRECISION XYZNOD(NDIM,NPOINT)
C
C VARIABLES LOCALES
C
      INTEGER          IN , IS , ISIGN , IL , IT , II
      INTEGER          IJKLUG , ITTOUR , IPTURB , ISENSF , IPOS
C
      DOUBLE PRECISION VALMAX
      DOUBLE PRECISION CR1X , CR1Y , CR1Z
      DOUBLE PRECISION CR2X , CR2Y , CR2Z
C
C***********************************************************************
C
C
C=======================================================================
C -1.  MACRO DE DEBUGGAGE DEVELOPPEUR
C=======================================================================
C
C            ATTENTION INTERVENTION DEVELOPPEUR UNIQUEMENT.
C
C     Si cette macro est vrai elle permet les impressions
C     dans le listing du detail des erreurs de reperage.
C
C     PAR DEFAUT : DEBUG_OUESTU = 0
C
C
#define DEBUG_OUESTU 0
C
C=======================================================================
C 0. Initialisations
C=======================================================================
C
      IERROR = 0
C
      ITTOUR = 0
C
      INDIAN = 0
C
      IJKLUG = 0
C
      IPTURB = 0
C
C Calcul de ValMax en local :  VALEUR MAX POUR L'ARRONDI
C
      VALMAX = EPZERO
C
C ---> Face
      DO IS = 1,ITYPF7
C
        II = ICONF7(IS)
        DO IL = 1,3
          VALMAX = MAX( VALMAX,ABS( XYZNOD(IL,II) ) )
        ENDDO
C
      ENDDO
C
C ---> Centre de gravite de la face
      VALMAX = MAX(VALMAX,ABS(CDGX))
      VALMAX = MAX(VALMAX,ABS(CDGY))
      VALMAX = MAX(VALMAX,ABS(CDGZ))
C
C ---> Centre de gravite de la cellule de depart
      VALMAX = MAX(VALMAX,ABS(CELX))
      VALMAX = MAX(VALMAX,ABS(CELY))
      VALMAX = MAX(VALMAX,ABS(CELZ))
C
C ---> Point de depart de la particule
      VALMAX = MAX(VALMAX,ABS(PX))
      VALMAX = MAX(VALMAX,ABS(PY))
      VALMAX = MAX(VALMAX,ABS(PZ))
C
C ---> Point d'arrivee de la particule
      VALMAX = MAX(VALMAX,ABS(QX))
      VALMAX = MAX(VALMAX,ABS(QY))
      VALMAX = MAX(VALMAX,ABS(QZ))
C
C=======================================================================
C 1. Position relative de P et Q : les 2 points sont-ils confondus ?
C=======================================================================
C
      CALL COLOCA
C     ===========
     &  ( VALMAX ,
     &    PX     , PY     , PZ     ,
     &    QX     , QY     , QZ     ,
     &    IPOS   )
C
C--> Si P et Q sont confondus, la particule est dans le meme cellule
C
      IF (IPOS.EQ.1) THEN
        INDIAN = -1
        RETURN
      ENDIF
C
C=======================================================================
C 2. Verification de l'orientation de la face
C=======================================================================
C
C   Ici on verifie dans quel sens sont lus les points
C   de la face par rapport a la cellule dans laquelle
C   se trouve la particule (au point P)
C   "bien" orientee : ISENSF = 1
C   "mal" orientee : ISENSF = -1
C
C--> Coordonnees du premier point support S(1) de la face
C
      II = ICONF7(1)
      CR1X = XYZNOD(1,II)
      CR1Y = XYZNOD(2,II)
      CR1Z = XYZNOD(3,II)
C
C--> Coordonnees du deuxieme point support S(2) de la face
C
      II = ICONF7(2)
      CR2X = XYZNOD(1,II)
      CR2Y = XYZNOD(2,II)
      CR2Z = XYZNOD(3,II)
C
C--> Orientation de PgS(1)S(2)
C
      CALL COTURN
C     ===========
     &  ( VALMAX ,
     &    PX     , PY    , PZ    ,
     &    CDGX   , CDGY  , CDGZ  ,
     &    CR1X   , CR1Y  , CR1Z  ,
     &    CR2X   , CR2Y  , CR2Z  ,
     &    ISENSF , IPTURB        )
C
C--> Si l'orientation precedente a echouee, on en essaie une autre,
C    mais elle est dangereuse dans le cas des cellules convaves.
C
      IF (ISENSF.EQ.0) THEN
C
C--> Orientation de GgS(1)S(2)
C
        CALL COTURN
C       ===========
     &  ( VALMAX ,
     &    CELX   , CELY  , CELZ  ,
     &    CDGX   , CDGY  , CDGZ  ,
     &    CR1X   , CR1Y  , CR1Z  ,
     &    CR2X   , CR2Y  , CR2Z  ,
     &    ISENSF , IPTURB        )
C
      ENDIF
C
      IF (ISENSF.EQ.0) THEN
#if DEBUG_OUESTU
        WRITE(NFECRA,1010)
#endif
        IERROR = 1
        RETURN
      ENDIF
C
C--------
C FORMATS
C--------
C
#if DEBUG_OUESTU
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     LA RECHERCHE DE L''ORIENTATION DE LA FACE A ECHOUEE    ',/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
#endif
C
C=======================================================================
C 2. Test sur le 1er point support de la face
C=======================================================================
C
C-->orientation de PQgS(1)
C
      CALL COTURN
C     ===========
     &  ( VALMAX ,
     &    PX     , PY    , PZ    ,
     &    QX     , QY    , QZ    ,
     &    CDGX   , CDGY  , CDGZ  ,
     &    CR1X   , CR1Y  , CR1Z  ,
     &    ISIGN  , IPTURB        )
C
      ISIGN = ISIGN * ISENSF
C
      IF (ISIGN.EQ.0) THEN
#if DEBUG_OUESTU
        WRITE(NFECRA,1020)
#endif
        IERROR = 1
        RETURN
      ENDIF
C
C--------
C FORMATS
C--------
C
#if DEBUG_OUESTU
 1020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     L''ORIENTATION DE PQgS(1) A ECHOUEE                    ',/,
     &'@                                                            ',/,
     &'@     Cause probable : les points PQgS(1) sont coplanaires.  ',/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
#endif
C
C=======================================================================
C 3. Boucle sur tous les points supports S(i) de la face
C    On cherche dans quel est le triangle construit sur les sommets et
C    le CDG de la face, qui est perce par la droite PQ
C=======================================================================
C
      DO IS = 2, ITYPF7
C
C--> Coordonnees des points supports
C
        II = ICONF7(IS)
        CR1X = XYZNOD(1,II)
        CR1Y = XYZNOD(2,II)
        CR1Z = XYZNOD(3,II)
C
C--> Orientation de PQgS(i) avec 2 =< i =< ITYPF7
C
        CALL COTURN
C       ===========
     &    ( VALMAX ,
     &      PX     , PY     , PZ     ,
     &      QX     , QY     , QZ     ,
     &      CDGX   , CDGY   , CDGZ   ,
     &      CR1X   , CR1Y   , CR1Z   ,
     &      IN     , IPTURB        )
C
        IN = IN * ISENSF
C
        IF (IN.EQ.0) THEN
#if DEBUG_OUESTU
          WRITE(NFECRA,1030) IS
#endif
          IERROR = 1
          RETURN
        ENDIF
C
C--> Si inversion de l'orientation de PQgS(i) alors orientation
C    de PQS(i-1)S(i)
C
        IF (ISIGN.EQ.-IN) THEN
C
C--> Le triangle gS(i-1)S(i) est t-il traverse par le rayon PQ ?
C    si oui, on repere le triangle en question
C
          IF (ISIGN.EQ.1) IJKLUG = IS
C
          ISIGN = IN
C
C--> Orientation de PQS(i-1)S(i)
C
          II = ICONF7(IS-1)
          CR2X = XYZNOD(1,II)
          CR2Y = XYZNOD(2,II)
          CR2Z = XYZNOD(3,II)
C
          CALL COTURN
C         ===========
     &      ( VALMAX ,
     &        PX     , PY     , PZ     ,
     &        QX     , QY     , QZ     ,
     &        CR2X   , CR2Y   , CR2Z   ,
     &        CR1X   , CR1Y   , CR1Z   ,
     &        IT     , IPTURB          )
C
          IT = IT * ISENSF
C
          IF (IT.EQ.0) THEN
#if DEBUG_OUESTU
            WRITE(NFECRA,1040)
#endif
            IERROR = 1
            RETURN
          ENDIF
C
          ITTOUR = ITTOUR + IT
C
        ENDIF
C
      ENDDO
C
      IF (ITTOUR.NE.-2 .AND. ITTOUR.NE.2 .AND. ITTOUR.NE.0) THEN
#if DEBUG_OUESTU
        WRITE(NFECRA,2010) ITTOUR
#endif
        IERROR = 1
        RETURN
      ELSE IF ((ITTOUR.EQ.-2 .OR. ITTOUR.EQ.2) .AND. IJKLUG.EQ.0) THEN
#if DEBUG_OUESTU
        WRITE(NFECRA,2020)
#endif
        IERROR = 1
        RETURN
      ENDIF
C
C--> Si la droite PQ ne traverse pas la face, ou s'elle rentre dans
C    la cellule par cette face, on retourne dans LAGCEL
C
      IF (ITTOUR.EQ.0 .OR. ITTOUR.EQ.-2) RETURN
C
C--------
C FORMATS
C--------
C
#if DEBUG_OUESTU
 1030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     L''ORIENTATION DE PQgS(i) A ECHOUEE, i = ',I2           ,/,
     &'@                                                            ',/,
     &'@     Cause probable : les points PQgS(i) sont coplanaires.  ',/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     L''ORIENTATION DE PQS(i-1)S(i) A ECHOUEE               ',/,
     &'@                                                            ',/,
     &'@     Cause Probable : points PQS(i-1)S(i) coplanaires.      ',/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     L''INDICE DE PASSAGE DOIT ESTRE EGALE A -2, 0 ou 2     ',/,
     &'@     IL EST CALCULE A : ',I2                                 ,/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     LA DETERMINATION DU TRIANGLE DE PASSAGE                ',/,
     &'@     DE LA PARTICULE A ECHOUEE                              ',/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
#endif
C
C=======================================================================
C 4. Position relative des points d'arrive et de depart et de la face
C=======================================================================
C
C--> Dans la cas ou la droite PQ construite sur les points d'arrive Q
C    et de depart P de la particule, sort de la cellule par la face
C    courante (ie ITTOUR = 2), on veut savoir si P et Q se trouve
C    de part et d'autre de la face, ou du meme cote :
C
C    ATTENTION : test perturbe (IPTURB=1) : si Q est sur la face,
C                la detection fonctionnera (a moins que QgS(i-1)s(i)
C                ne soient coplanaires).
C
      IPTURB = 1
C
      II = ICONF7(IJKLUG-1)
      CR1X = XYZNOD(1,II)
      CR1Y = XYZNOD(2,II)
      CR1Z = XYZNOD(3,II)
C
      II = ICONF7(IJKLUG)
      CR2X = XYZNOD(1,II)
      CR2Y = XYZNOD(2,II)
      CR2Z = XYZNOD(3,II)
C
C--> Orientation de QgS(i-1)S(i)
C
          CALL COTURN
C         ===========
     &      ( VALMAX ,
     &        QX     , QY     , QZ     ,
     &        CDGX   , CDGY   , CDGZ   ,
     &        CR1X   , CR1Y   , CR1Z   ,
     &        CR2X   , CR2Y   , CR2Z   ,
     &        ISIGN  , IPTURB          )
C
      IF (ISIGN.EQ.0) THEN
#if DEBUG_OUESTU
        WRITE(NFECRA,3010)
#endif
        IERROR = 1
        RETURN
      ENDIF
C
C--> Position relative entre P Q et la face
C
      INDIAN = -ISIGN * ISENSF
C
C--------
C FORMATS
C--------
C
#if DEBUG_OUESTU
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,
     &'@     LA POSITION RELATIVE DE Q PAR RAPPORT A P              ',/,
     &'@     A LA FACE A ECHOUEE                                    ',/,
     &'@                                                            ',/,
     &'@     Cause Probable : points QgS(i-1)S(i) coplanaires.      ',/,
     &'@                                                            ',/,
     &'@  La particule est perdue.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
#endif
C
C----
C FIN
C----
C
      END
c@z
