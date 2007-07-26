c@a
                        SUBROUTINE GETFBR
C                       ******************
C     -------------------------------------------------------------
     & ( FSTR , FACNB, FACES)
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
c@versb
CVERS    CODE_SATURNE VERSION 1.3              COPYRIGHT EDF 1999-2006 *
c@verse
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA LISTE DES FACES DE BORD CORRESPONDANT A
CFONC A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! FSTR         ! A  !  -> ! CHAINE A INTERPRETER                 !
CARGU ! FACES        ! E  !  <- ! FACES SELECTIONNEES                  !
CARGU ! FACNB        ! E  !  <- ! NOMBRE DE FACES SELECTIONNEES        !
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
      CHARACTER*(*)    FSTR
      INTEGER      FACES(*), FACNB
C
C VARIABLES LOCALES
C
      INTEGER          LENSTR
C
C***********************************************************************
C
      LENSTR=LEN(FSTR)
      CALL CSGFBR(FSTR, LENSTR, FACNB, FACES)
      END



                        SUBROUTINE GETCEL
C                       ******************
C     -------------------------------------------------------------
     & ( FSTR , CELLNB, CELLS)
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
c@versb
CVERS    CODE_SATURNE VERSION 1.3              COPYRIGHT EDF 1999-2006 *
c@verse
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA LISTE DES CELLULES CORRESPONDANT A
CFONC A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! FSTR         ! A  !  -> ! CHAINE A INTERPRETER                 !
CARGU ! CELLS        ! E  !  <- ! FACES SELECTIONNEES                  !
CARGU ! CELLNB       ! E  !  <- ! NOMBRE DE FACES SELECTIONNEES        !
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
      CHARACTER*(*) FSTR
      INTEGER       CELLS(*), CELLNB
C
C VARIABLES LOCALES
C
      INTEGER          LENSTR
C
C***********************************************************************
C
      LENSTR=LEN(FSTR)
      CALL CSGCEL(FSTR, LENSTR, CELLNB, CELLS)
      END


                        SUBROUTINE GETFAC
C                       ******************
C     -------------------------------------------------------------
     & ( FSTR , FACNB, FACES)
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
c@versb
CVERS    CODE_SATURNE VERSION 1.3              COPYRIGHT EDF 1999-2006 *
c@verse
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA LISTE DES FACES INTERNES CORRESPONDANT A
CFONC A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! FSTR         ! A  !  -> ! CHAINE A INTERPRETER                 !
CARGU ! CELLS        ! E  !  <- ! FACES SELECTIONNEES                  !
CARGU ! CELLNB       ! E  !  <- ! NOMBRE DE FACES SELECTIONNEES        !
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
      CHARACTER*(*) FSTR
      INTEGER       FACES(*), FACNB
C
C VARIABLES LOCALES
C
      INTEGER          LENSTR
C
C***********************************************************************
C
      LENSTR=LEN(FSTR)
      CALL CSGFAC(FSTR, LENSTR, FACNB, FACES)
      END
