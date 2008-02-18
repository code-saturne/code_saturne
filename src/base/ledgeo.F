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
                        SUBROUTINE LEDGEO
C                       *****************
C     -------------------------------------------------------------
     &     ( NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &       NPRFML , NFML   , NSOM   , LNDFAC , LNDFBR )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC LECTURE DES DIMENSIONS DES TABLEAUX ENTITES GEOMETRIQUES
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NDIM         ! E  ! <-  ! DIMENSION DE L'ESPACE (=3)           !
CARGU ! NCELET       ! E  ! <-  ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  ! <-  ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  ! <-  ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  ! <-  ! NOMBRE DE FACES DE BORD              !
CARGU ! NPRFML       ! E  ! <-  ! NOMBRE DE PROPIETES DES FAMILLES     !
CARGU !              !    !     ! DE FACES DE BORD                     !
CARGU ! NFML         ! E  ! <-  ! NOMBRE DE FAMILLES DE FACES DE BORD  !
CARGU ! NSOM         ! E  ! <-  ! NOMBRE DE SOMMETS DU MAILLAGE        !
CARGU ! LNDFAC       ! E  ! <-  ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  ! <-  ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
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
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NPRFML , NFML   , NSOM   , LNDFAC , LNDFBR
C
C VARIABLES LOCALES
C
      INTEGER          IEL , IFAC , IFML , NN , N1 , IP
      INTEGER          NTETRA , NPYRAM , NPRISM , NHEXAE
C
C***********************************************************************
C
      NDIM = 0
      NCEL = 0
      NFAC = 0
      NFABOR = 0
      NPRFML = 0
      NFML   = 0
      NSOM   = 0
C
      LNDFAC = 0
      LNDFBR = 0
C
C=======================================================================
C 1.  Ouverture de FICGEO et lecture des entetes principales
C=======================================================================
C
C---> Ouverture
C
      OPEN (FILE=FICGEO,UNIT=IMPGEO,FORM='formatted',ERR=9000)
C
C---> Ecriture avec rubrique ... avant
C
      NDIM = 3
C
      READ (IMPGEO,   *,ERR=9002,END=9003)
      READ (IMPGEO,   *,ERR=9002,END=9003)
      READ (IMPGEO,1100,ERR=9002,END=9003) NCEL, NFAC, NFABOR, NSOM
C
      READ (IMPGEO,   *,ERR=9002,END=9003)
      READ (IMPGEO,   *,ERR=9002,END=9003)
      READ (IMPGEO,1100,ERR=9002,END=9003) NTETRA,NPYRAM,NPRISM,NHEXAE
C
      READ (IMPGEO,   *,ERR=9002,END=9003)
      READ (IMPGEO,   *,ERR=9002,END=9003)
      READ (IMPGEO,1100,ERR=9002,END=9003) NPRFML, NFML
C
C=======================================================================
C 2.  Positionnement dans FICGEO pour la lecture des tableaux optionnels
C     de connectivites faces sommets (IPNFAC, IPNFBR, NODFAC, NODFBR)
C=======================================================================
C
C-->  Positionnement dans le fichier avant lecture
C
C     Connectivite faces internes - cellules
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1,NFAC
        READ(IMPGEO,*)
      ENDDO
C
C     Connectivite faces de bord - cellule
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1, NFABOR
        READ(IMPGEO,*)
      ENDDO
C
C     Coordonnees du centre des cellules
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IEL = 1,NCEL
        READ(IMPGEO,*)
      ENDDO
C
C     Surfaces des faces internes
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1,NFAC
        READ(IMPGEO,*)
      ENDDO
C
C     Surfaces des faces de bord
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1,NFABOR
        READ(IMPGEO,*)
      ENDDO
C
C     Coordonnees du centre des faces
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1,NFAC
        READ(IMPGEO,*)
      ENDDO
C
C     Coordonnees du centre des faces de bord
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1,NFABOR
        READ(IMPGEO,*)
      ENDDO
C
C     Familles des faces de bord
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFAC = 1,NFABOR
        READ(IMPGEO,*)
      ENDDO
C
C     Proprietes des familles
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO IFML = 1,NFML
        READ(IMPGEO,*)
      ENDDO
C
C     Coordonnees des noeuds
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO N1 = 1,NSOM
        READ(IMPGEO,*)
      ENDDO
C
C     Connectivite cellules points
      READ(IMPGEO,*)
      READ(IMPGEO,*)
      DO N1 = 1,NCEL
        READ(IMPGEO,*)
      ENDDO
C
C=======================================================================
C 3.  Calcul des dimensions LNDFAC et LNDFBR si possible
C=======================================================================
C
C     Remarque : si les tableaux ne sont pas disponibles,
C                on aura une erreur de fin de fichier.
C                Dans ce cas, on sort en fermant normalement
C                le fichier, en mettant LNDFAC et LNDFBR à 0.
C
C-->  Dimension du tableau de connectivite faces de bord -> points
C
      READ(IMPGEO,*,ERR=9000,END=100)
      READ(IMPGEO,*,ERR=9000,END=100)
      LNDFBR = 0
      DO N1 = 1,NFABOR
        READ(IMPGEO,1200,ERR=9000,END=100) NN , IP
        LNDFBR = LNDFBR + IP
      ENDDO
C
C-->  Dimension du tableau de connectivite faces internes -> points
C
      READ(IMPGEO,*,ERR=9000,END=100)
      READ(IMPGEO,*,ERR=9000,END=100)
      LNDFAC = 0
      DO N1 = 1,NFAC
        READ(IMPGEO,1200,ERR=9000,END=100) NN , IP
        LNDFAC = LNDFAC + IP
      ENDDO
C
C-->  Si l'on a pu lire les tableaux de connectivite faces -> points,
C     on conserve les dimensions LNDFBR et LNDFAC calculees ;
C     sinon, on les remet a zero
C
      GOTO 200
C
 100  LNDFBR = 0
      LNDFAC = 0
C
 200  CONTINUE
C
C=======================================================================
C 3.  Fermeture du fichier et mise a jour des structures C
C=======================================================================
C
      CLOSE(IMPGEO)
C
C---> MISE A JOUR STRUCTURES C
C
      NCELET = NCEL
C
      CALL DIMGEO
C     ===========
     &     (NDIM  , NCELET, NCEL  , NFAC  , NFABOR , NSOM  ,
     &      LNDFAC, LNDFBR, NFML  , NPRFML,
     &      NTETRA, NPYRAM, NPRISM, NHEXAE)
C
C---> FORMATS
C
 1100 FORMAT(20I10)
 1200 FORMAT(I10,I10)
C
      RETURN
C
C ERREURS
C
 9000 CONTINUE
      WRITE(NFECRA,8000)FICGEO
      CALL CSEXIT (1)
 9002 CONTINUE
      WRITE(NFECRA,8002)
     &            FICGEO,NDIM,NCEL,NFAC,NFABOR,NPRFML,NFML  ,NSOM
      CALL CSEXIT (1)
 9003 CONTINUE
      WRITE(NFECRA,8003)
     &            FICGEO,NDIM,NCEL,NFAC,NFABOR,NPRFML,NFML  ,NSOM
      CALL CSEXIT (1)
C
C FORMATS
C
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''OUVERTURE DU FICHIER GEOMETRIE   ',/,
     &'@    *********                                               ',/,
     &'@      ERREUR DANS ledgeo POUR LE FICHIER ',A6                ,/,
     &'@                                                            ',/,
     &'@      Le calcul ne peut etre execute.                       ',/,
     &'@                                                            ',/,
     &'@      Verifier le fichier geometrie (existence, droits      ',/,
     &'@        d''acces, recopie correcte dans le repertoire de    ',/,
     &'@        travail).                                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,
     &'@    *********                                               ',/,
     &'@      ERREUR DANS ledgeo POUR LE FICHIER ',A6                ,/,
     &'@                                                            ',/,
     &'@      L''etat des dimensions lues est le suivant :          ',/,
     &'@        NDIM      NCEL      NFAC    NFABOR                  ',/,
     &'@  ',4I10                                                     ,/,
     &'@        NPRFML    NFML      NSOM                            ',/,
     &'@  ',3I10                                                     ,/,
     &'@                                                            ',/,
     &'@      Le calcul ne peut etre execute.                       ',/,
     &'@                                                            ',/,
     &'@      Verifier le fichier geometrie (existence, droits      ',/,
     &'@        d''acces, format...).                               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,
     &'@    *********                                               ',/,
     &'@      FIN PREMATUREE DANS ledgeo POUR LE FICHIER ',A6        ,/,
     &'@                                                            ',/,
     &'@      L''etat des dimensions lues est le suivant :          ',/,
     &'@        NDIM      NCEL      NFAC    NFABOR                  ',/,
     &'@  ',4I10                                                     ,/,
     &'@        NPRFML    NFML      NSOM                            ',/,
     &'@  ',3I10                                                     ,/,
     &'@                                                            ',/,
     &'@      Le calcul ne peut etre execute.                       ',/,
     &'@                                                            ',/,
     &'@      Verifier le fichier geometrie (existence, droits      ',/,
     &'@        d''acces, format...).                               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
      END
c@z
