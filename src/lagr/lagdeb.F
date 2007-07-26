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
                        SUBROUTINE LAGDEB
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICOCEL , ITYCEL ,
     &   IDEVEL , ITUSER , IA     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC  Ce sous-programme est lu lors du premier passage dans le module
CFONC  Lagrangien (IPLAR=1).
CFONC
CFONC  Construction du stockage compact (tableaux "positions"/"valeurs")
CFONC  des connectivites cellules -> faces
CFONC  (pour tous formats de maillage).
CFONC  On part des connectivites faces->cellule IFACEL et IFABOR.
CFONC  Les faces de bord sont stockees avec un numero negatif.
CFONC
CFONC
CFONC  On rappelle que la lecture des connectivites face->sommets
CFONC  pour un maillage .slc est faite dans LAGGET, appelee par LETGEO.
CFONC
CFONC
CFONC
CFONC ----------------------------------------------------------------
CFONC  Schema d'association entre les tableaux ITYCEL "positions"
CFONC  et ICOCEL "valeurs"
CFONC ----------------------------------------------------------------
CFONC
CFONC
CFONC          Tableau des valeurs ICOCEL (de dimension LNDNOD)
CFONC
CFONC
CFONC .---.------..------.------.------..------.------..------.------.
CFONC |   | ...  ||IFACn | ...  |      || IFACm| ...  || ...  |      |
CFONC `---'------'`------'------'------'`------'------'`------'------'
CFONC   1           iVal         jVal-1   jVal                  LNDNOD
CFONC                 |                     |                      |
CFONC                 |                     |                      |
CFONC                 `----------.      .---'         .------------'
CFONC                            |      |             |
CFONC                            |      |             |
CFONC         .-----.-------.------.------.-------.------.--------.
CFONC         |  1  |  ...  | iVal | jVal |  ...  |LNDNOD|LNDNOD+1|
CFONC         `-----'-------'------'------'-------'------'--------'
CFONC            1            iPos  iPos+1         nPos-1 nPos=NCELET+1
CFONC
CFONC
CFONC          Tableau des positions ITYCEL (de dimension NCELET+1)
CFONC
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
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDNOD       ! E  !  -> ! DIM. CONNECTIVITE CELLULES->FACES    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! ICOCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (LNDNOD)     !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (NCELET+1)   !    !     !    POINTEUR DU TABLEAU ICOCEL        !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ICOCEL(LNDNOD) , ITYCEL(NCELET+1)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE) , RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
C
      INTEGER          IEL , IFAC , IP , N1
      INTEGER          NBFAC1 , NBFAC2
C
C***********************************************************************
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. Connectivite cellules -> faces
C    Les faces de bord sont stockees avec des valeurs negatives
C=======================================================================
C
C-->Initialisation
C
      DO N1 = 1,NCELET+1
        ITYCEL(N1) = 0
      ENDDO
      DO N1 = 1,LNDNOD
        ICOCEL(N1) = 0
      ENDDO
C
C-->Premier passage : on compte le nombre de faces par cellules
C
      DO IFAC = 1,NFAC
        DO IP = 1,2
          IEL = IFACEL(IP,IFAC)
          ITYCEL(IEL) = ITYCEL(IEL) + 1
        ENDDO
      ENDDO
C
      DO IFAC = 1,NFABOR
        IEL = IFABOR(IFAC)
        ITYCEL(IEL) = ITYCEL(IEL) + 1
      ENDDO
C
C-->Verification
C
      IP = 0
      DO IEL = 1,NCELET
        IP = IP + ITYCEL(IEL)
      ENDDO
      IF (IP.NE.LNDNOD) THEN
        WRITE(NFECRA,5000) LNDNOD , IP
        CALL CSEXIT (1)
      ENDIF
C
C-->Deuxieme passage : on calcule les positions du tableau ICOCEL
C
      NBFAC2 = ITYCEL(1)
      ITYCEL(1) = 1
      DO IEL = 2, NCELET+1
        NBFAC1 = NBFAC2
        NBFAC2 = ITYCEL(IEL)
        ITYCEL(IEL) = ITYCEL(IEL-1) + NBFAC1
      ENDDO
C
C-->Verification
C
      IF (ITYCEL(NCELET+1).NE.LNDNOD+1) THEN
        WRITE(NFECRA,1000) ITYCEL(NCELET+1) , LNDNOD+1
        CALL CSEXIT (1)
      ENDIF
C
C-->Troisieme passage : on remplit le tableau de connectivites ICOCEL
C   Rem : les faces de bord sont stockees avec des valeurs negatives
C
      DO IFAC = 1,NFAC
        DO IP = 1,2
          IEL = IFACEL(IP,IFAC)
          N1 = ITYCEL(IEL)
          DO WHILE (ICOCEL(N1).NE.0)
            N1 = N1 + 1
            IF (N1.EQ.ITYCEL(IEL+1)) THEN
              WRITE(NFECRA,2000)
              CALL CSEXIT (1)
            ENDIF
          ENDDO
          ICOCEL(N1) = IFAC
        ENDDO
      ENDDO
C
      DO IFAC = 1,NFABOR
        IEL = IFABOR(IFAC)
        N1 = ITYCEL(IEL)
        DO WHILE (ICOCEL(N1).NE.0)
          N1 = N1 + 1
          IF (N1.EQ.ITYCEL(IEL+1)) THEN
            WRITE(NFECRA,3000)
            CALL CSEXIT (1)
          ENDIF
        ENDDO
        ICOCEL(N1) = -IFAC
      ENDDO
C
C-->Archive : ancienne version du 3e passage NE PAS SUPPRIMER SVP
C   Attention ne marche pas dans le cas de la TURBINE,
C   peut-etre a cause du fait que deux cellules ont en commun
C   plusieurs faces
C
c     N1 = 0
c     DO IEL = 1,NCEL
c       DO IFAC = 1,NFAC
c         DO IP = 1,2
c           IF (IFACEL(IP,IFAC).EQ.IEL) THEN
c             N1 = N1 + 1
c             ICOCEL(N1) = ICOCEL(N1) - IFAC
c           ENDIF
c         ENDDO
c       ENDDO
c       DO IFAC = 1,NFABOR
c         IF (IFABOR(IFAC).EQ.IEL) THEN
c           N1 = N1 + 1
c           ICOCEL(N1) = ICOCEL(N1)  - IFAC
c         ENDIF
c       ENDDO
c     ENDDO
C
C
C-->Verification
C
      DO N1 = 1,LNDNOD
        IF (ICOCEL(N1).EQ.0) THEN
          WRITE(NFECRA,4000)
          CALL CSEXIT (1)
        ENDIF
      ENDDO
C
C***********************************************************************
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    Erreur dans la construction de la connectivite          ',/,
     &'@    cellules -> faces (LAGDEB)                              ',/,
     &'@                                                            ',/,
     &'@    Mauvais remplissage de ITYCEL : le dernier              ',/,
     &'@    element de ITYCEL ne vaut pas LNDNOD+1                  ',/,
     &'@                                                            ',/,
     &'@    ITYCEL(NCELET+1) = ',I10                                 ,/,
     &'@    LNDNOD+1         = ',I10                                 ,/,
     &'@                                                            ',/,
     &'@    Verifier le remplissage de ITYCEL dans LAGDEB           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    Erreur dans la construction de la connectivite          ',/,
     &'@    cellules -> faces (LAGDEB)                              ',/,
     &'@                                                            ',/,
     &'@    Mauvais remplissage de ITYCEL avec les faces internes   ',/,
     &'@                                                            ',/,
     &'@    Verifier le remplissage de ITYCEL dans LAGDEB           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    Erreur dans la construction de la connectivite          ',/,
     &'@    cellules -> faces (LAGDEB)                              ',/,
     &'@                                                            ',/,
     &'@    Mauvais remplissage de ITYCEL avec les faces de bord    ',/,
     &'@                                                            ',/,
     &'@    Verifier le remplissage de ITYCEL dans LAGDEB           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    Erreur dans la construction de la connectivite          ',/,
     &'@    cellules -> faces (LAGDEB)                              ',/,
     &'@                                                            ',/,
     &'@    Mauvais remplissage de ICOCEL                           ',/,
     &'@    (ICOCEL contient au moins une valeur nulle)             ',/,
     &'@                                                            ',/,
     &'@    Verifier le remplissage de ICOCEL dans LAGDEB           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 5000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    Erreur dans la construction de la connectivite          ',/,
     &'@    cellules -> faces (LAGDEB)                              ',/,
     &'@                                                            ',/,
     &'@    La taille du tableau ICOCEL n''est plus identique a     ',/,
     &'@    celle calculee dans LAGINI                              ',/,
     &'@                                                            ',/,
     &'@     - LAGINI : LNDNOD = ',I10                               ,/,
     &'@     - LAGDEB : LNDNOD = ',I10                               ,/,
     &'@                                                            ',/,
     &'@    Verifier le calcul de LNDNOD dans LAGINI et LAGDEB      ',/,
     &'@    ou l''integrite des tableaux IFACEL et IFABOR           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
