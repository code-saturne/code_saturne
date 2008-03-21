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
                        SUBROUTINE LAGNWC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NPT    , NPTNEW , NEW    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYCEL , ICOCEL ,
     &   IFRLAG , ISORTI , IWORKP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   SURFBN , ETTP   ,
     &   RDEVEL , RTUSER , RA     )
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
CFONC        INJECTION EN CONTINUE DES PARTICULES
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
CARGU ! NPT          ! E  ! <-  ! NOMBRE COURANT DE PARTICULES         !
CARGU ! NPTNEW       ! E  !  -> ! NOMBRE TOTAL DE NOUVELLES PARTICULES !
CARGU !              !    !     ! POUR TOUTES LES ZONES D'INJECTION    !
CARGU ! NEW          ! E  !  -> ! NOMBRE DE NOUVELLES PART A INJECTER  !
CARGU !              !    !     ! POUR LA ZONE D'INJECTION COURANTE    !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU !   (NFABOR)   !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU !   (NFABOR)   !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU !   (NCELET)   !    !     !                                      !
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
CARGU ! ICOCEL       ! TE !  -> ! CONNECTIVITE CELLULES -> FACES       !
CARGU !   (LNDNOD)   !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE !  -> ! CONNECTIVITE CELLULES -> FACES       !
CARGU !   (NCELET+1) !    !     !    POINTEUR DU TABLEAU ICOCEL        !
CARGU ! ISORTI       ! TE ! <-> ! POUR CHAQUE PARTICULE :              !
CARGU !   (NBPMAX)   !    !     !    * NUMERO DE SA CELLULE            !
CARGU !              !    !     !    * 0 SI SORTIE DU DOMAINE          !
CARGU ! IWORKP(NPTNEW! TE !  -> ! NUMERO DE LA FACE D'INJECTION        !
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
      INCLUDE "period.h"
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
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ICOCEL(LNDNOD) , ITYCEL(NCELET+1)
      INTEGER          ISORTI(NBPMAX)
      INTEGER          IFRLAG(NFABOR) , IWORKP(NPTNEW)
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
      INTEGER          IDEBIA, IDEBRA, IFINIA
C
      INTEGER          NP , IEL , N1  , IFAC , KFAC , NBP
      INTEGER          II , JJ  , IN  , ISORT
      INTEGER          INDIAN , IFAOLD , IFANEW
      INTEGER          ITYPFO , ICONFO(100)
      INTEGER          IDEHOR , IERRIE , ICECPT
      INTEGER          ICELCR, IPERCR ,ITEPAS, IPER
C
      DOUBLE PRECISION RD(1)
      DOUBLE PRECISION XF, YF, ZF
      DOUBLE PRECISION UP, VP, WP, UF, VF, WF
      DOUBLE PRECISION PTA(3), PTB(3), VECT(3), VECTN(3)
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C Traitement de la periodicite
C
      IF (IPERIO.EQ.1) THEN
C
        ICELCR = IDEBIA
        IPERCR = ICELCR + NCELET-NCEL
        IFINIA = IPERCR + NCELET-NCEL
        CALL IASIZE('LAGNWC', IFINIA)
C       ===========
C
        DO IEL = 1,NCELET-NCEL
          IA(ICELCR+IEL-1) = 0
          IA(IPERCR+IEL-1) = -1
        ENDDO
        CALL PERLOC(IA(ICELCR), IA(IPERCR))
C       ===========
C
      ENDIF
C
C=======================================================================
C 2. Injection des particules
C=======================================================================
C
C     BOUCLE SUR LES NOUVELLES PARTICULES :
C
      DO NP = 1,NEW
C
C   RAZ du compteur de cellules traversees
C
        ICECPT = 0
C
C-->IDEHOR = 1 si exterieur du domaine de calcul
C
        IDEHOR = 0
C
C      Incrementation du pointeur sur les particules
C
        NPT = NPT + 1
C
C      On sauvegarde la maille de depart
C
        ISORT = ISORTI(NPT)
C
C      Tirage aleatoire
C
        N1 = 1
        CALL ZUFALL(N1,RD)
C
C      Nouvelle position : on le fait dans la
C      direction de la vitesse mais possibilite
C      de le faire par rapport a la normale
C
        XF = ETTP(NPT,JXP) + RD(1) *DTP *ETTP(NPT,JUP)
        YF = ETTP(NPT,JYP) + RD(1) *DTP *ETTP(NPT,JVP)
        ZF = ETTP(NPT,JZP) + RD(1) *DTP *ETTP(NPT,JWP)
C
        UP = ETTP(NPT,JUP)
        VP = ETTP(NPT,JVP)
        WP = ETTP(NPT,JWP)
C
        UF = ETTP(NPT,JUF)
        VF = ETTP(NPT,JVF)
        WF = ETTP(NPT,JWF)
C
        IF (  XF.EQ. ETTP(NPT,JXP)
     &  .AND. YF.EQ. ETTP(NPT,JYP)
     &  .AND. ZF.EQ. ETTP(NPT,JZP) ) GOTO 300
C
C Est-ce que le point est dans le domaine (Clone de LAGCEL)
C
C      Numero de face d 'injection dans l'element
C
        IFANEW = IWORKP(NPT)
C
 100    CONTINUE
C
          IEL    = ISORT
          IFAOLD = IFANEW
          INDIAN = 0
          ICECPT = ICECPT + 1
C
C         ---> Elimination des particules qui posent problemes
C              la particule reste au niveau de la face d'entree
C              (boucles infinies)
C
          IF (ICECPT.GT.30) THEN
            IDEHOR = 1
            GOTO 200
          ENDIF
C
C      --> balayage des KFAC faces entourant la cellule IEL
C           Elles sont stockees entre ITYCEL(IEL) et ITYCEL(IEL+1)-1
C          (donc KFAC ne peut pas valoir ITYCEL(IEL+1)...)
C
          KFAC = ITYCEL(IEL)-1
C
          DO WHILE (INDIAN.EQ.0)
C
            KFAC = KFAC + 1
C
C           --> Erreur : la particule reste au niveau de la face
C                        d'entree
C
            IF (KFAC.EQ.ITYCEL(IEL+1)) THEN
              IDEHOR = 1
              GOTO 200
            ENDIF
C
            IFAC = ICOCEL(KFAC)
C
C-->boucle sur les faces internes
C
C-->si la face interne a deja ete traitee dans la cellule precedente
C   son numero est dans IFAOLD et on ne la retraite pas une seconde fois
C
            IF (IFAC.GT.0 .AND. IFAC.NE.IFAOLD) THEN
C
              IN = 0
              DO NBP = IPNFAC(IFAC),IPNFAC(IFAC+1)-1
                IN = IN + 1
                ICONFO(IN) = NODFAC(NBP)
              ENDDO
              ITYPFO = IPNFAC(IFAC+1) - IPNFAC(IFAC) + 1
              ICONFO(ITYPFO) = ICONFO(1)
C
              CALL OUESTU
C             ===========
     &   (    NFECRA , NDIM   , NNOD ,
     &        IERRIE ,
     &        ETTP(NPT,JXP)  , ETTP(NPT,JYP)  , ETTP(NPT,JZP)  ,
     &        XF                , YF                , ZF               ,
     &        CDGFAC(1,IFAC)    , CDGFAC(2,IFAC)    , CDGFAC(3,IFAC)   ,
     &        XYZCEN(1,IEL)     , XYZCEN(2,IEL)     , XYZCEN(3,IEL)    ,
     &        ITYPFO , ICONFO , XYZNOD ,
     &        INDIAN )
C
              IF (IERRIE.EQ.1) THEN
                IDEHOR = 1
                GOTO 200
C
C-->si la particule passe dans la cellule voisine
C
              ELSE IF (INDIAN.EQ.1) THEN
C
C-->si la particule NPT est dans la cellule II alors le voisin ne
C   peut etre que la cellule JJ, et vice versa, et inversement, et
C   ainsi de suite.
C
                IFANEW = IFAC
C
                II = IFACEL(1,IFAC)
                JJ = IFACEL(2,IFAC)
C
                IF (IEL.EQ.II) THEN
                  ISORT = JJ
                ELSE IF (IEL.EQ.JJ) THEN
                  ISORT = II
                ENDIF
C
C Traitement de la periodicite
C Meme commentaire que dans lagcel.F
C
                IF (ISORT.GT.NCEL) THEN
C
                  ITEPAS = ISORT
                  ISORT = IA(ICELCR+ITEPAS-NCEL-1)
C
C                On recupere les informations sur la periodicite
C
                  IPER  = IA(IPERCR+ITEPAS-NCEL-1)
C
C                 MODIFICATION DE LA POSITION
C
                  PTA(1)= XF
                  PTA(2)= YF
                  PTA(3)= ZF
C
                  CALL LAGPER(IPER,PTA, PTB)
C                 ===========
C
                  XF = PTB(1)
                  YF = PTB(2)
                  ZF = PTB(3)
C
C                 MODIFICATION DE LA VITESSE
C
                  VECT(1) = UP
                  VECT(2) = VP
                  VECT(3) = WP
C
                  CALL LAGVEC(IPER, VECT, VECTN)
C                 ===========
C
                  UP = VECTN(1)
                  VP = VECTN(2)
                  WP = VECTN(3)
C
C                 MODIFICATION DE LA VITESSE FLUIDE VUE
C
                  VECT(1) = UF
                  VECT(2) = VF
                  VECT(3) = WF
C
                  CALL LAGVEC(IPER, VECT, VECTN)
C                 ===========
C
                  UF = VECTN(1)
                  VF = VECTN(2)
                  WF = VECTN(3)
C
                  IFANEW = 0
C
                ENDIF
C
C--> Retour pour balayage des face de la cellule suivante
C
                GOTO 100
C
              ENDIF
C
C--> Balayage des faces de bord (reperees par leur valeur negative
C    dans ICOCEL)
C
C    resultat : INDIAN =  0 le rayon PQ ne sort pas de la cellule par
C    ~~~~~~~~               cette face
C               INDIAN = -1 meme cellule
C               INDIAN =  1 interaction avec la frontiere
C
            ELSE IF (IFAC.LT.0 .AND. IFAC.NE.IFAOLD) THEN
C
              IFAC = -IFAC
C
              IN = 0
              DO NBP = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                IN = IN + 1
                ICONFO(IN) = NODFBR(NBP)
              ENDDO
              ITYPFO = IPNFBR(IFAC+1) - IPNFBR(IFAC) + 1
              ICONFO(ITYPFO) = ICONFO(1)
C
              CALL OUESTU
C             ===========
     &   (    NFECRA , NDIM   , NNOD ,
     &        IERRIE ,
     &        ETTP(NPT,JXP)  , ETTP(NPT,JYP)  , ETTP(NPT,JZP)  ,
     &        XF                , YF                , ZF               ,
     &        CDGFBO(1,IFAC)    , CDGFBO(2,IFAC)    , CDGFBO(3,IFAC)   ,
     &        XYZCEN(1,IEL)     , XYZCEN(2,IEL)     , XYZCEN(3,IEL)    ,
     &        ITYPFO , ICONFO , XYZNOD ,
     &        INDIAN )
C
C-->si la trajectoire de la particule traverse la face de bord
C    alors on laisse la particule au niveau de la face d'entree
C
              IF (IERRIE.EQ.1 .OR. INDIAN.EQ.1) THEN
                IDEHOR = 1
                GOTO 200
              ENDIF
C
            ENDIF
C
C fin de DO WHILE
          ENDDO
C
C-->Fin de la boucle principale sur les particules
C
 200      CONTINUE
C
C Si le point est dans le domaine, alors on injecte la
C particule du point sinon on laisse la particule au niveau
C de l'entree
C
        IF (IDEHOR.EQ.0) THEN
          ETTP(NPT,JXP) = XF
          ETTP(NPT,JYP) = YF
          ETTP(NPT,JZP) = ZF
          ISORTI(NPT)   = ISORT
        ENDIF
C
 300      CONTINUE
C
      ENDDO
C
C**********************************************************************
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
