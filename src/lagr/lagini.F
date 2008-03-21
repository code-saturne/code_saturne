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
                        SUBROUTINE LAGINI
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   , NFAC   , NFABOR ,
     &   LNDNOD ,
     &   IFACEL , IFABOR ,
     &   NBRFAC ,
     &   IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC    Remplissage de LNDNOD pour l'allocation de memoire dans MEMLA1
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS                    !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! LNDNOD       ! E  ! <-  ! DIM. CONNECT. CELLULES->FACES        !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! NBRFAC(NCEL) ! TE !  -  ! TABLEAU DE TRAVAIL ENTIER            !
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
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
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "period.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
      INTEGER          NFAC   , NFABOR
      INTEGER          LNDNOD
      INTEGER          IFACEL(2,NFAC)  , IFABOR(NFABOR)
      INTEGER          NBRFAC(NCELET)
      INTEGER          IA(*)
C
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL , IFAC , IP
C***********************************************************************
C=======================================================================
C 0. GESTION DE LA MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. Calcul de la dimension du tableau de connectivites cellules->faces
C=======================================================================
C
C-->Initialisation
C
      DO IEL = 1,NCELET
        NBRFAC(IEL) = 0
      ENDDO
C
C-->Pour chaque cellule on compte le nombre de faces internes qui l'entourent
C
      DO IFAC = 1,NFAC
        DO IP = 1,2
          IEL = IFACEL(IP,IFAC)
          NBRFAC(IEL) = NBRFAC(IEL) + 1
        ENDDO
      ENDDO
C
C-->Pour chaque cellule on compte le nombre de faces de bord qui l'entourent
C
      DO IFAC = 1,NFABOR
        IEL = IFABOR(IFAC)
        NBRFAC(IEL) = NBRFAC(IEL) + 1
      ENDDO
C
C-->verif : chaque cellule ne peut avoir moins de 4 faces en regard
C   Rem : d'apres la routine LEDGEO on a toujours NDIM = 3
C
      IP = 0
      DO IEL = 1,NCEL
        IF (NBRFAC(IEL).LT.4) THEN
          IP = IP + 1
        ENDIF
      ENDDO
      IF (IP.GT.0) THEN
        WRITE(NFECRA,9000) IP
        CALL CSEXIT (1)
      ENDIF
C
C--->Traitement de la periodicite
C
      IF (IPERIO.EQ.1) THEN
C
        IP = 0
        DO IEL = NCEL+1,NCELET
          IF (NBRFAC(IEL).NE.1) THEN
            IP = IP + 1
          ENDIF
        ENDDO
        IF (IP.GT.0) THEN
          WRITE(NFECRA,9001) IP
          CALL CSEXIT (1)
        ENDIF
C
      ENDIF
C
C-->Calcul de la dimension du tableau de connectivite
C
      LNDNOD = 0
      DO IEL = 1,NCELET
        LNDNOD = LNDNOD + NBRFAC(IEL)
      ENDDO
C
C=======================================================================
C 2. Ouverture du fichier listing specifique lagrangien
C=======================================================================
C
      OPEN ( UNIT=IMPLAL, FILE=FICLAL,
     &     STATUS='UNKNOWN', FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      REWIND(IMPLAL)
C
C--------
C FORMATS
C--------
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@  Il y a ',I10,' cellules qui ont moins de 4 faces.         ',/,
     &'@   Erreur rencontree dans LAGINI (module Lagrangien).       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@  Il y a ',I10,' cellules du halo periodique qui            ',/,
     &'@   ne comportent pas qu''une unique face.                   ',/,
     &'@   Erreur rencontree dans LAGINI (module Lagrangien).       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
