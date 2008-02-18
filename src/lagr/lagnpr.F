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
                        SUBROUTINE LAGNPR
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NPT    , NPTNEW , NEW    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IZONE  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFRLAG , ISORTI , IWORKP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   SURFBN , ETTP   ,
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
CFONC    TIRAGE ALEATOIRE DE LA POSITION DES PARTICULES A INJECTER
CFONC    DANS LE DOMAINE DE CALCUL AU NIVEAU DES FACES
CFONC    D'UNE ZONE DE COULEUR : POSITIONS + REPERAGE DE LA CELLULE
CFONC    ATTENTION : CE TIRAGE ALEATOIRE NE FONCTIONNE QU'AVEC DES FACES
CFONC    TRIANGULAIRE OU QUADRANGULAIRES.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU !              !    !     !                                      !
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
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! NPT          ! E  ! <-  ! NOMBRE COURANT DE PARTICULES         !
CARGU ! NPTNEW       ! E  !  -> ! NOMBRE TOTAL DE NOUVELLES PARTICULES !
CARGU !              !    !     ! POUR TOUTES LES ZONES D'INJECTION    !
CARGU ! NEW          ! E  !  -> ! NOMBRE DE NOUVELLES PART A INJECTER  !
CARGU !              !    !     ! POUR LA ZONE D'INJECTION COURANTE    !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IZONE        ! E  !  -> ! NUMERO  DE LA ZONE D'INJECTION       !
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
CARGU ! IFRLAG       ! TE !  -> ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU !   (NFABOR)   !    !     !  POUR LE MODULE LAGRANGIEN           !
CARGU ! ISORTI       ! TE ! <-  ! POUR CHAQUE PARTICULE :              !
CARGU !   (NBPMAX)   !    !     !    * NUMERO DE SA CELLULE            !
CARGU !              !    !     !    * 0 SI SORTIE DU DOMAINE          !
CARGU ! IWORKP(NPTNEW! TE ! <-  ! NUMERO DE LA FACE D'INJECTION        !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! SURFBN(NFABOR! TR !  -> ! SURFACE DES FACES DE BORD            !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
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
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C**********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NPT    , NPTNEW , NEW
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IZONE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFRLAG(NFABOR) , ISORTI(NBPMAX) , IWORKP(NPTNEW)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION SURFBN(NFABOR)
      DOUBLE PRECISION ETTP(NBPMAX,NVP)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          NNN , NNNN , IFAC , II
      INTEGER          IFRBR , MINFAC , MAXFAC
      INTEGER          ICONFO(100)
C
      DOUBLE PRECISION SURFM , RDA(1) , RD , EPS , PM1 , PM6
      DOUBLE PRECISION CTR(6,3) , VEC(3) , ARE(2,3) , SURFTR(2)
C
C***********************************************************************
C
      EPS = 1.D-3
C
C     CALCUL DE LA SURFACE MAX DE L'ENTREE :
C
      SURFM = -10.D0
      MINFAC = NFABOR+1
      MAXFAC = 0
C
      DO IFAC = 1,NFABOR
        IFRBR = IFRLAG(IFAC)
        IF (IFRBR.EQ.IZONE) THEN
          SURFM = MAX( SURFM,SURFBN(IFAC) )
          MINFAC = MIN (IFAC,MINFAC)
          MAXFAC = MAX (IFAC,MAXFAC)
        ENDIF
      ENDDO
C
C     STOP PAR SECURITE
C
      IF (MAXFAC.EQ.0.OR.MINFAC.EQ.NFABOR+1) THEN
        WRITE(NFECRA,9000) IZONE
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C     BOUCLE SUR LES FACES DE LA ZONE
C
      DO IFAC = 1,NFABOR
        IFRBR = IFRLAG(IFAC)
        IF (IFRBR.EQ.IZONE) THEN
C
C       incrementation du pointeur sur les particules
C
          NPT = NPT + 1
C
C    ATTENTION :
C
C         type de face : 3 ou 4 points supports
C         pour l'instant je ne sais pas traiter les autres
C         avec plus de points supports...
C
          II = IPNFBR(IFAC+1)-IPNFBR(IFAC)
          IF (II.GT.4) GOTO 100
C
C        si face a 4 points, on choisit l'un des deux triangles
C
          IF (II.EQ.4) THEN
C
            II = 0
            DO NNNN = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
              II = II+1
              ICONFO(II) = NODFBR(NNNN)
            ENDDO
C
C        longueur des arretes 1 et 2 du premier triangle :
C
            ARE(1,1) = XYZNOD(1,ICONFO(2)) - XYZNOD(1,ICONFO(1))
            ARE(1,2) = XYZNOD(2,ICONFO(2)) - XYZNOD(2,ICONFO(1))
            ARE(1,3) = XYZNOD(3,ICONFO(2)) - XYZNOD(3,ICONFO(1))
            ARE(2,1) = XYZNOD(1,ICONFO(3)) - XYZNOD(1,ICONFO(1))
            ARE(2,2) = XYZNOD(2,ICONFO(3)) - XYZNOD(2,ICONFO(1))
            ARE(2,3) = XYZNOD(3,ICONFO(3)) - XYZNOD(3,ICONFO(1))
C
C        surface du premier triangle
C
            VEC(1) = ARE(1,2)*ARE(2,3)-ARE(1,3)*ARE(2,2)
            VEC(2) = ARE(1,3)*ARE(2,1)-ARE(1,1)*ARE(2,3)
            VEC(3) = ARE(1,1)*ARE(2,2)-ARE(1,2)*ARE(2,1)
            SURFTR(1) = SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2)
C
C        longueur des arretes 1 et 2 du deuxieme triangle :
C
            ARE(1,1) = XYZNOD(1,ICONFO(3)) - XYZNOD(1,ICONFO(1))
            ARE(1,2) = XYZNOD(2,ICONFO(3)) - XYZNOD(2,ICONFO(1))
            ARE(1,3) = XYZNOD(3,ICONFO(3)) - XYZNOD(3,ICONFO(1))
            ARE(2,1) = XYZNOD(1,ICONFO(4)) - XYZNOD(1,ICONFO(1))
            ARE(2,2) = XYZNOD(2,ICONFO(4)) - XYZNOD(2,ICONFO(1))
            ARE(2,3) = XYZNOD(3,ICONFO(4)) - XYZNOD(3,ICONFO(1))
C
C        surface du deuxieme triangle
C
            VEC(1) = ARE(1,2)*ARE(2,3) - ARE(1,3)*ARE(2,2)
            VEC(2) = ARE(1,3)*ARE(2,1) - ARE(1,1)*ARE(2,3)
            VEC(3) = ARE(1,1)*ARE(2,2) - ARE(1,2)*ARE(2,1)
            SURFTR(2) = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
C
C        tirage d'un nombre aleatoire entre 0 et 1
C        pour determiner quel triangle choisir
C
            NNN = 1
            CALL ZUFALL(NNN,RDA)
            RD=RDA(1)
C
C        si le deuxieme triangle est choisit, on reorganise
C        les points : 4 <--> 2
C
            IF (RD.LE.( SURFTR(2) / (SURFTR(1)+SURFTR(2)) )) THEN
              NNNN = ICONFO(4)
              ICONFO(4) = ICONFO(2)
              ICONFO(2) = NNNN
            ENDIF
C
C dans le cas ou la face est un triangle...
C
          ELSE IF (II.EQ.3) THEN
C
            II = 0
            DO NNNN = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
              II = II+1
              ICONFO(II) = NODFBR(NNNN)
            ENDDO
C
          ENDIF
C
C        constitution des coordonnees du triangle
C
          DO NNNN = 1,3
            CTR(NNNN,1) = XYZNOD(1,ICONFO(NNNN))
            CTR(NNNN,2) = XYZNOD(2,ICONFO(NNNN))
            CTR(NNNN,3) = XYZNOD(3,ICONFO(NNNN))
          ENDDO
C
C        tirage aleatoire d'un point dans le triangle
C        constitue des points 1,2,3
C
 200      CONTINUE
C
C        1) tirage du point 4 sur l'arete 12
C
 300      CONTINUE
          NNN = 1
          CALL ZUFALL(NNN,RDA)
          RD=RDA(1)
          IF (RD.EQ.0.D0 .OR. RD.EQ.1.D0) GOTO 300
C
          DO NNNN = 1,3
            CTR(4,NNNN) = RD*CTR(1,NNNN) + (1.D0-RD)*CTR(2,NNNN)
          ENDDO
C
C        2) tirage du point 5 sur l'arete 13
C
 400      CONTINUE
          NNN = 1
          CALL ZUFALL(NNN,RDA)
          RD=RDA(1)
          IF (RD.EQ.0.D0 .OR. RD.EQ.1.D0) GOTO 400
C
          DO NNNN = 1,3
            CTR(5,NNNN) = RD*CTR(1,NNNN) +( 1.D0-RD)*CTR(3,NNNN)
          ENDDO
C
C        3) le point 6 est le sommet du parallelogramme 1465
C
          DO NNNN = 1,3
            CTR(6,NNNN) = CTR(4,NNNN) + CTR(5,NNNN) - CTR(1,NNNN)
          ENDDO
C
C        4) reste a verifier que le point 6 appartient au triangle 123
C
C        4.1) vecteur normal au triangle : 12^13
C
          VEC(1) = (CTR(2,2)-CTR(1,2))*(CTR(3,3)-CTR(1,3))
     &           - (CTR(2,3)-CTR(1,3))*(CTR(3,2)-CTR(1,2))
          VEC(2) = (CTR(2,3)-CTR(1,3))*(CTR(3,1)-CTR(1,1))
     &           - (CTR(2,1)-CTR(1,1))*(CTR(3,3)-CTR(1,3))
          VEC(3) = (CTR(2,1)-CTR(1,1))*(CTR(3,2)-CTR(1,2))
     &           - (CTR(2,2)-CTR(1,2))*(CTR(3,1)-CTR(1,1))
C
C        4.2) produit mixte pour le point 1 :
C
          PM1 = 0.D0
          PM1 = PM1 + VEC(1) *
     &        ( (CTR(2,2)-CTR(1,2))*(CTR(3,3)-CTR(2,3))
     &         -(CTR(2,3)-CTR(1,3))*(CTR(3,2)-CTR(2,2)) )
          PM1 = PM1 + VEC(2) *
     &        ( (CTR(2,3)-CTR(1,3))*(CTR(3,1)-CTR(2,1))
     &         -(CTR(2,1)-CTR(1,1))*(CTR(3,3)-CTR(2,3)) )
          PM1 = PM1 + VEC(3) *
     &        ( (CTR(2,1)-CTR(1,1))*(CTR(3,2)-CTR(2,2))
     &         -(CTR(2,2)-CTR(1,2))*(CTR(3,1)-CTR(2,1)) )
C
C        4.3) produit mixte pour le point 6 :
C
          PM6 = 0.D0
          PM6 = PM6 + VEC(1) *
     &        ( (CTR(2,2)-CTR(6,2))*(CTR(3,3)-CTR(2,3))
     &         -(CTR(2,3)-CTR(6,3))*(CTR(3,2)-CTR(2,2)) )
          PM6 = PM6 + VEC(2) *
     &        ( (CTR(2,3)-CTR(6,3))*(CTR(3,1)-CTR(2,1))
     &         -(CTR(2,1)-CTR(6,1))*(CTR(3,3)-CTR(2,3)) )
          PM6 = PM6 + VEC(3) *
     &        ( (CTR(2,1)-CTR(6,1))*(CTR(3,2)-CTR(2,2))
     &         -(CTR(2,2)-CTR(6,2))*(CTR(3,1)-CTR(2,1)) )
C
C        4.4) 6 est dans le triangle si PM1*PM6>=0
C
            IF (PM1*PM6.LT.0.D0) GOTO 200
C
C        5) POUR PLUS DE SECURITE, ON DEPLACE LE POINT
C           D'UN EPSILON EN DIRECTION DU POINT PRESSION
C
C ATTENTION : CE DECALAGE PEUT ETRE DANGEREUX DANS LE CAS
C             DE CELULES CONCAVES
C
          CTR(6,1) = CTR(6,1) + (XYZCEN(1,IFABOR(IFAC))-CTR(6,1))*EPS
          CTR(6,2) = CTR(6,2) + (XYZCEN(2,IFABOR(IFAC))-CTR(6,2))*EPS
          CTR(6,3) = CTR(6,3) + (XYZCEN(3,IFABOR(IFAC))-CTR(6,3))*EPS
C
C        LE TRAITEMENT EST TERMINE POUR LE POINT NPT,
C        ON REMPLIT LES TABLEAUX POUR LE LAGRANGIEN :
C
          ETTP(NPT,JXP) = CTR(6,1)
          ETTP(NPT,JYP) = CTR(6,2)
          ETTP(NPT,JZP) = CTR(6,3)
C
          ISORTI(NPT) = IFABOR(IFAC)
          IWORKP(NPT) = IFAC
C
        ENDIF
C
 100    CONTINUE
C
      ENDDO
C
C***********************************************************************
C
C-------
C FORMAT
C-------
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGNPR).                                   ',/,
     &'@                                                            ',/,
     &'@    PROBLEME DANS LA GESTION DE NOUVELLES PARTICULES        ',/,
     &'@                                                            ',/,
     &'@  Le nombre de faces de la zone ',I10                        ,/,
     &'@    est egal a zero.                                        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
