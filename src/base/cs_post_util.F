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
                        SUBROUTINE PSTCWR
C                       ******************
C     -------------------------------------------------------------
     & ( NUMGEP , NOMCAS , NOMREP , NOMFMT , OPTFMT ,
     &   INDMOD , NTCHR )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  CREATION D'UN "WRITER" A PARTIR DES DONNEES FOURNIES PAR LA
CFONC  COUCHE FORTRAN : ENCAPSULATION COUCHE C POUR LA TRANSMISSION
CFONC  DES LONGUEURS DES CHAINES DE CARACTERES
CFONC
CFONC  UN WRITER CORRESPOND AU CHOIX D'UN NOM DE CAS, DE REPERTOIRE,
CFONC  ET DE FORMAT, AINSI QU'UN INDICATEUR PRECISANT SI LES MAILLAGES
CFONC  ASSOCIES DOIVENT DEPENDRE OU NON DU TEMPS, ET LA FREQUENCE
CFONC  DE SORTIE PAR DEFAUT POUR LES VARIABLES ASSOCIEES

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NUMGEP       ! E  !  -> ! IDENTIFICATEUR DU GESTIONNAIRE       !
CARGU !              !    !     ! (< 0 POUR GESTIONNAIRE RESERVE,      !
CARGU !              !    !     !  > 0 POUR GESTIONNAIRE UTILISATEUR)  !
CARGU ! NOMCAS       ! A  !  -> ! NOM DU CAS ASSOCIE                   !
CARGU ! NOMREP       ! A  !  -> ! NOM DU REPERTOIRE ASSOCIE            !
CARGU ! NOMFMT       ! A  !  -> ! NOM DE FORMAT ASSOCIE                !
CARGU ! OPTFMT       ! E  !  -> ! OPTIONS ASSOCIEES AU FORMAT          !
CARGU ! INDMOD       ! E  !  -> ! 0 : MAILLAGES FIGES                  !
CARGU !              !    !     ! 1 : MAILLAGES DEFORMABLES            !
CARGU !              !    !     ! 2 : MAILLAGES MODIFIABLES            !
CARGU ! NTCHR        ! E  !  -> ! FREQUENCE DE SORTIE PAR DEFAUT       !
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
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*32     NOMCAS , NOMFMT
      CHARACTER*96     NOMREP , OPTFMT
      INTEGER          NUMGEP , INDMOD , NTCHR
C
C VARIABLES LOCALES
C
      INTEGER          LNMCAS , LNMREP , LNMFMT , LOPFMT
C
C***********************************************************************
C
      LNMCAS = LEN(NOMCAS)
      LNMREP = LEN(NOMREP)
      LNMFMT = LEN(NOMFMT)
      LOPFMT = LEN(OPTFMT)
C
      CALL PSTCW1 (NUMGEP, NOMCAS, NOMREP, NOMFMT, OPTFMT,
C     ===========
     &             LNMCAS, LNMREP, LNMFMT, LOPFMT,
     &             INDMOD, NTCHR)
C
      RETURN
C
      END
c@z
c@a
                        SUBROUTINE PSTCMA
C                       ******************
C     -------------------------------------------------------------
     & ( NUMMAI , NOMMAI ,
     &   NBRCEL , NBRFAC , NBRFBR , LSTCEL , LSTFAC , LSTFBR )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  CREATION D'UN MAILLAGE DE POST TRAITEMENT A PARTIR DES DONNEES
CFONC  FOURNIES PAR LA COUCHE FORTRAN : ENCAPSULATION COUCHE C
CFONC  POUR LA TRANSMISSION DES LONGUEURS DES CHAINES DE CARACTERES
CFONC
CFONC  LES LISTES DE CELLULES OU FACES A EXTRAIRE SONT TRIEES EN SORTIE,
CFONC  QU'ELLES LE SOIENT DEJA EN ENTREE OU NON.
CFONC
CFONC  LA LISTE DES CELLULES ASSOCIEES N'EST NECESSAIRE QUE SI LE NOMBRE
CFONC  DE CELLULES A EXTRAIRE EST STRICTEMENT SUPERIEUR A 0 ET INFERIEUR
CFONC  AU NOMBRE DE CELLULES DU MAILLAGE.
CFONC
CFONC  LES LISTES DE FACES NE SONT PRISES EN COMPTE QUE SI LE NOMBRE DE
CFONC  CELLULES A EXTRAIRE EST NUL ; SI LE NOMBRE DE FACES DE BORD A
CFONC  EXTRAIRE EST EGAL AU NOMBRE DE FACES DE BORD DU MAILLAGE GLOBAL,
CFONC  ET LE NOMBRE DE FACES INTERNES A EXTRAIRE EST NUL, ALORS ON
CFONC  EXTRAIT PAR DEFAUT LE MAILLAGE DE BORD, ET LA LISTE DES FACES DE
CFONC  BORD ASSOCIEES N'EST DONC PAS NECESSAIRE.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NUMMAI       ! E  !  -> ! IDENTIFICATEUR DU MAILLAGE           !
CARGU !              !    !     ! (< 0 POUR MAILLAGE RESERVE,   ,      !
CARGU !              !    !     !  > 0 POUR MAILLAGE UTILISATEUR)      !
CARGU ! NOMMAI       ! A  !  -> ! NOM DU MAILLAGE ASSOCIE              !
CARGU ! NBRCEL       ! E  !  -> ! NOMBRE DE CELLULES ASSOCIEES         !
CARGU ! NBRFAC       ! E  !  -> ! NOMBRE DE FACES INTERNES ASSOCIEES   !
CARGU ! NBRFBR       ! E  !  -> ! NOMBRE DE FACES DE BORD ASSOCIEES    !
CARGU ! LSTCEL       ! E  !  -> ! LISTE DES CELLULES ASSOCIEES         !
CARGU !              ! E  !     ! (INUTILE SI NBRCEL >= NCEL)          !
CARGU ! LSTFAC       ! E  !  -> ! LISTE DES FACES INTERNES ASSOCIEES   !
CARGU ! LSTFBR       ! E  !  -> ! LISTE DES FACES DE BORD ASSOCIEES    !
CARGU !              ! E  !     ! (INUTILE SI    NBRFBR = NFABOR       !
CARGU !              ! E  !     !             ET NBRFAC = 0     )      !
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
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*32     NOMMAI
      INTEGER          NUMMAI, NBRCEL, NBRFAC, NBRFBR
C
      INTEGER          LSTCEL(NBRCEL), LSTFAC(NBRFAC), LSTFBR(NBRFBR)
C
C VARIABLES LOCALES
C
      INTEGER          LNMMAI
C
C***********************************************************************
C
      LNMMAI = LEN(NOMMAI)
C
      CALL PSTCM1 (NUMMAI, NOMMAI, LNMMAI,
C     ===========
     &             NBRCEL, NBRFAC, NBRFBR, LSTCEL, LSTFAC, LSTFBR)
C
      RETURN
C
      END
c@z
c@a
                        SUBROUTINE PSTEVA
C                       *****************
C     -------------------------------------------------------------
     & ( NUMMAI , NOMVAR , DIMVAR , IENTLA , IVARPR , NTCABS , TTCABS ,
     &   VARCEL , VARFAC , VARFBO )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  ECRITURE D'UN CHAMP DE POST TRAITEMENT ASSOCIE AUX CELLULES
CFONC  OU FACES D'UN MAILLAGE A PARTIR DES DONNEES FOURNIES PAR LA
CFONC  COUCHE FORTRAN :
CFONC  ENCAPSULATION COUCHE C POUR LA TRANSMISSION DES LONGUEURS DES
CFONC  CHAINES DE CARACTERES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NUMMAI       ! A  !  -> ! NUMERO DU MAILLAGE ASSOCIE           !
CARGU ! NOMVAR       ! E  !  -> ! NOM DE LA VARIABLE ASSOCIEE          !
CARGU ! DIMVAR       ! E  !  -> ! 1 POUR SCALAIRE, 3 POUR VECTEUR      !
CARGU ! IENTLA       ! E  !  -> ! SI VECTEUR, 1 SI VALEURS ENTRELACEES !
CARGU !              !    !     ! (X1, Y1, Z1, X2, Y2, ..., YN, ZN),   !
CARGU !              !    !     ! 0 SINON (X1, X2, ...XN, Y1, Y2, ...) !
CARGU ! IVARPR       ! E  !  -> ! 1 SI VARIABLE DEFINIE SUR MAILLAGE   !
CARGU !              !    !     ! "PARENT", 0 SI VARIABLE RESTREINTE   !
CARGU !              !    !     ! AU MAILLAGE NUMMAI                   !
CARGU ! NTCABS       ! E  !  -> ! NUMERO DE PAS DE TEMPS (-1 POUR UNE  !
CARGU !              !    !     ! VARIABLE INDEPENDANTE DU TEMPS)      !
CARGU ! TTCABS       ! R  !  -> ! TEMPS PHYSIQUE ASSOCIE               !
CARGU ! VARCEL(*)    ! R  !  -> ! VALEURS AUX CELLULES ASSOCIEES       !
CARGU ! VARFAC(*)    ! R  !  -> ! VALEURS AUX FACES INTERNES ASSOCIEES !
CARGU ! VARFBO(*)    ! R  !  -> ! VALEURS AUX FACES DE BORD ASSOCIEES  !
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
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*32     NOMVAR
      INTEGER          NUMMAI, DIMVAR, IENTLA, IVARPR, NTCABS
      DOUBLE PRECISION TTCABS, VARCEL(*), VARFAC(*), VARFBO(*)
C
C VARIABLES LOCALES
C
      INTEGER          LNMVAR
C
C***********************************************************************
C
      LNMVAR = LEN(NOMVAR)
C
      CALL PSTEV1 (NUMMAI, NOMVAR, LNMVAR, DIMVAR, IENTLA, IVARPR,
C     ===========
     &             NTCABS, TTCABS, VARCEL, VARFAC, VARFBO)
C
      RETURN
C
      END
c@z
c@a
                        SUBROUTINE PSTSNV
C                       *****************
C     -------------------------------------------------------------
     & ( NOMVAR , NOMVA2 , NOMVA3 )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  SUPPRESSION DU CARACTERE X, x, OU 1 D'UNE CHAINE DE CARACTERES
CFONC  FORTRAN SI LES CHAINES COMPAREES SONT IDENTIQUES AU DERNIER
CFONC  CARACTERE PRES, RESPECTIVEMENT Y, y, OU 2 ET Z, z, OU 3
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NOMVAR       ! E  !  -> ! NOM DE LA VARIABLE ASSOCIEE          !
CARGU ! NOMVA2       ! E  !  -> ! NOM DE LA VARIABLE 2 ASSOCIEE        !
CARGU ! NOMVA3       ! E  !  -> ! NOM DE LA VARIABLE 3 ASSOCIEE        !
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
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*32     NOMVAR, NOMVA2, NOMVA3
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ
      INTEGER          LNMVAR, LNMVA2, LNMVA3
C
C***********************************************************************
C
      LNMVAR = LEN(NOMVAR)
      LNMVA2 = LEN(NOMVA2)
      LNMVA3 = LEN(NOMVA3)
C
      IF ((LNMVAR .EQ. LNMVA2) .AND. (LNMVAR .EQ. LNMVA3)) THEN
C
        DO 10 II = LNMVAR, 1, -1
          IF (     NOMVAR(II:II) .NE. ' '
     &        .OR. NOMVA2(II:II) .NE. ' '
     &        .OR. NOMVA3(II:II) .NE. ' ') THEN
            GOTO 20
          ENDIF
 10     CONTINUE
C
 20     CONTINUE
C
        IF (II .GT. 1) THEN
C
          JJ = II
C
C         On prevoit le cas ou c'est l'avant-dernier caractere
C         qui change, comme avec VitesX1, VitesX2, ... en
C         cas de calcul avec plusieurs phases
C
          IF (      (II .GT. 2)
     &        .AND. (NOMVAR(II:II) .EQ. NOMVA2(II:II))
     &        .AND. (NOMVAR(II:II) .EQ. NOMVA3(II:II))) THEN
            II = JJ-1
          ENDIF
C
C         On supprime le caractere lie a la dimension
C
          IF (      NOMVAR(II:II) .EQ. 'X'
     &        .AND. NOMVA2(II:II) .EQ. 'Y'
     &        .AND. NOMVA3(II:II) .EQ. 'Z') THEN
            NOMVAR(II:II) = ' '
          ELSE IF (      NOMVAR(II:II) .EQ. 'x'
     &             .AND. NOMVA2(II:II) .EQ. 'y'
     &             .AND. NOMVA3(II:II) .EQ. 'z') THEN
            NOMVAR(II:II) = ' '
          ELSE IF (      NOMVAR(II:II) .EQ. '1'
     &             .AND. NOMVA2(II:II) .EQ. '2'
     &             .AND. NOMVA3(II:II) .EQ. '3') THEN
            NOMVAR(II:II) = ' '
          ENDIF
C
C         Si l'on a supprime l'avant-dernier caractere, on
C         decale le dernier caractere
C
          IF (II .EQ. JJ+1) THEN
            NOMVAR(II:II) = NOMVAR(JJ:JJ)
            NOMVAR(JJ:JJ) = ' '
          ENDIF
C
        ENDIF
C
      ENDIF
C
      RETURN
C
      END
c@z

