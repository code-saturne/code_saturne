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
                        SUBROUTINE LAGCEL
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ICOCEL , ITYCEL , IFRLAG , ITEPA  , IBORD  ,
     &   INDEP  , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   SURFBN ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , PARBOR , AUXL   ,
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
CFONC    Trajectographie des particules : le propos de ce sous-programme
CFONC    est de donner le numero de la cellule d'arrive d'une particule
CFONC    connaissant les coordonnees du point de depart, celles
CFONC    du point d'arrive, ainsi que le numero de la cellule de depart.
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
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE ! <-  ! TAB D'INDIRECTION POUR TRI DES FACES !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ICOCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU !   (LNDNOD)   !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU !   (NCELET+1) !    !     !                                      !
CARGU ! IFRLAG       ! TE !  -> ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU !   (NFABOR)   !    !     !  POUR LE MODULE LAGRANGIEN           !
CARGU ! ITEPA        ! TE ! <-  ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IBORD        ! TE ! <-  ! SI NORDRE=2, CONTIENT LE NUMERO DE LA!
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
CARGU ! INDEP        ! TE ! <-  ! POUR CHAQUE PARTICULE :              !
CARGU !   (NBPMAX)   !    !     !   NUMERO DE LA CELLULE DE DEPART     !
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
CARGU ! SURFBN(NFABOR! TR ! ->  ! SURFACE DES FACES DE BORD            !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT ET PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! PARBOR(NFABOR! TR ! <-> ! CUMUL DES STATISTIQUES AUX FRONTIERES!
CARGU !    NVISBR)   !    !     !                                      !
CARGU ! AUXL(NBPMAX,3! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstphy.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
      INTEGER          ICOCEL(LNDNOD) , ITYCEL(NCELET+1)
      INTEGER          IFRLAG(NFABOR) , ITEPA(NBPMAX,NIVEP)
      INTEGER          IBORD(NBPMAX)
      INTEGER          INDEP(NBPMAX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION SURFBN(NFABOR)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR) , AUXL(NBPMAX,3)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
C
      INTEGER          IEL, IFAC, KFAC, NBP, ICECPT
      INTEGER          II, JJ, IN, IP
      INTEGER          INDIAN, IFAOLD, IFANEW
      INTEGER          ISUIVI, IERROR, IERRIE
      INTEGER          ITYPFO, ICONFO(100)
C
      INTEGER          ICELCR , IPERCR, ITEPAS, IPER
      DOUBLE PRECISION PTA(3), PTB(3), VECT(3), VECTN(3)
C
C=======================================================================
C -1.  MACRO DE DEBUGGAGE DEVELOPPEUR
C=======================================================================
C
C            ATTENTION INTERVENTION DEVELOPPEUR UNIQUEMENT.
C
C     Si cette macro est vrai elle permet les impressions listing
C     de l'avance des particules (si IMPLTG = 1), et permet la sortie
C     des fichers Ensight de debuggage en cas de perte de particule
C     en mode sans erreur (si IERRIE = 1).
C     Lorsque cette macro est vrai, il est conseille d'avoir un NBPART
C     petit, sans quoi les temps de calcul seront enormes a cause
C     des ecritures disque du au fichier SCRATCH4.lag.
C
C     PAR DEFAUT : DEBUG_LAGCEL = 0
C
C
#define DEBUG_LAGCEL 0
C
C
C
#if DEBUG_LAGCEL
      INTEGER          IMPLTG
      INTEGER          IPT    , KPT    , NPT    , IPART  , IELOLD
      INTEGER          NQUAD4 , NTRIA3 , NSIDED
#endif
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
C 1.  INITIALISATION
C=======================================================================
C
      DO IP = 1,NBPMAX
        IBORD(IP) = 0
      ENDDO
C
#if DEBUG_LAGCEL
C--> IMPLTG =1 ecriture de la progression des particules dans le listing
C    (utile si  DEBUG_LAGCEL = 1)
C
C mettre IMPLTG a 1 ici
      IMPLTG = 1
C
#endif
C
C
C--> Si IERRIE est non nul alors on se place en mode sans erreur de
C    trajectoire : il y a arret du calcul a la 1ere erreur sur le
C    reperage.
C
      IERRIE = 0
C
C--> La valeur de IERR bascule a 1 s'il y a une erreur et si IERRIE = 1
C
C ne pas modifier cette initialisation SVP
      IERR = 0
C
      NBPERR = 0
      DNBPER = 0.D0
C
C
C--> Si on est en instationnaire, RAZ des statistiques aux frontieres
C
      IF (IENSI3.EQ.1) THEN
C
        IF (ISTTIO.EQ.0 .OR. (ISTTIO.EQ.1 .AND. IPLAS.LE.NSTBOR)) THEN
          TSTATP = 0.D0
          NPSTF = 0
          DO II = 1,NVISBR
            DO IFAC = 1,NFABOR
              PARBOR(IFAC,II) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
        TSTATP = TSTATP + DTP
        NPSTF  = NPSTF  + 1
        NPSTFT = NPSTFT + 1
C
      ENDIF
C
      IF (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) THEN
        NPENCR = 0
        DNPENC = 0.D0
      ENDIF
C
C
C Traitement de la periodicite
C
      IF (IPERIO.EQ.1) THEN
C
        ICELCR = IDEBIA
        IPERCR = ICELCR + NCELET-NCEL
        IFINIA = IPERCR + NCELET-NCEL
        CALL IASIZE('LAGCEL', IFINIA)
C       ===========
C
        DO IEL = 1,NCELET-NCEL
          IA(ICELCR+IEL-1) = 0
          IA(IPERCR+IEL-1) = -1
        ENDDO
C
        CALL PERLOC(IA(ICELCR), IA(IPERCR))
C       ===========
C
      ENDIF
C
C=======================================================================
C 2. Reperage des particules dans le maillage
C=======================================================================
C
C--> Boucle principale : les particules sont traitees une par une
C
      DO IP = 1,NBPART
C
C utile car IDEPO2
        IF (ITEPA(IP,JISOR).GT.0) THEN
C
          IFANEW = 0
          ICECPT = 0
C
C--> Ouverture du fichier SCRATCH qui permettra l'ecriture des infos
C    en mode DEBUG
C
#if DEBUG_LAGCEL
          OPEN (IMPLA4,FILE='SCRATCH4.lag',
     &          STATUS='UNKNOWN',FORM='UNFORMATTED',
     &          ACCESS='SEQUENTIAL')
          NQUAD4 = 0
          NTRIA3 = 0
          NSIDED = 0
#endif
C
C--> Etiquette de retour pour le suivi de la particule dans la cellule voisine
C
 100      CONTINUE
C
#if DEBUG_LAGCEL
          IF (IMPLTG.EQ.1) WRITE(NFECRA,9001) IP, ITEPA(IP,JISOR)
#endif
C
C--> Attention :
C
C 1) IEL est cst jusqu'au GOTO 100, ITEPA(IP,JISOR) est modifie
C    jusqu'au GOTO 100
C
C 2) IFAOLD contient le numero de la derniere face traversee par la
C    trajectoire de la particule IP, si > 0 face interne,
C    si < 0 face de bord, IFAOLD est cst jusqu'au GOTO 100,
C    IFANEW est modifie jusqu'au GOTO 100
C
          IEL    = ITEPA(IP,JISOR)
          ISUIVI = -999
          IFAOLD = IFANEW
          INDIAN = 0
          ICECPT = ICECPT + 1
C
C         ---> Elimination des particules qui posent problemes
C              (boucles infinies)
C
          IF (ICECPT.GT.30) THEN
            ITEPA(IP,JISOR) = 0
            NBPERR = NBPERR + 1
            DNBPER = DNBPER + TEPA(IP,JRPOI)
            IF (IERRIE.EQ.0) THEN
#if DEBUG_LAGCEL
              IF (IMPLTG.EQ.1) WRITE(NFECRA,9103) IP
#endif
              GOTO 200
            ELSE
              WRITE(NFECRA,9103) IP
              IERR = 1
              GOTO 300
            ENDIF
          ENDIF
C
C--> Balayage des KFAC faces entourant la cellule IEL,
C    elles sont stockees entre ITYCEL(IEL) et ITYCEL(IEL+1)-1
C    (donc KFAC ne peut pas valoir ITYCEL(IEL+1)...)
C
C    Remarque : si on veut supprimer le DO WHILE il faut utiliser
C               la ligne suivante
C         DO KFAC = ITYCEL(IEL),ITYCEL(IEL+1)-1
C
          KFAC = ITYCEL(IEL)-1
C
          DO WHILE (INDIAN.EQ.0)
C
            KFAC = KFAC + 1
C
C--> Cas ou aucune face a INDIAN= -1 ou 1 sur la cellule n'a ete
C    detectee, gestion de l'erreur :
C
            IF (KFAC.EQ.ITYCEL(IEL+1)) THEN
              IF (IERRIE.EQ.0) THEN
                ITEPA(IP,JISOR) = 0
                NBPERR = NBPERR + 1
                DNBPER = DNBPER + TEPA(IP,JRPOI)
#if DEBUG_LAGCEL
                IF (IMPLTG.EQ.1)  WRITE(NFECRA,9102) IEL,IP
#endif
                GOTO 200
              ELSE
                WRITE(NFECRA,9102) IEL,IP
                IERR = 1
                GOTO 300
               ENDIF
            ENDIF
C
            IFAC = ICOCEL(KFAC)
C
C--> Boucle sur les faces internes (numero positif dans ICOCEL)
C    resultat : INDIAN =  0 le rayon PQ ne sort pas de la cellule
C    ~~~~~~~~               par cette face
C               INDIAN = -1 meme cellule
C               INDIAN =  1 sortie de la cellule par cette face
C
C--> Si la face interne a deja ete traitee dans la cellule precedente
C    son numero est dans IFAOLD et on ne la retraite pas une 2eme fois
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
#if DEBUG_LAGCEL
              NBP = IPNFAC(IFAC+1) - IPNFAC(IFAC)
              WRITE(IMPLA4)
     &        NBP,IFAC,IEL,(NODFAC(IN),IN=IPNFAC(IFAC),IPNFAC(IFAC+1)-1)
              IF (NBP.EQ.4) THEN
                NQUAD4 = NQUAD4 + 1
              ELSE IF (NBP.EQ.3) THEN
                NTRIA3 = NTRIA3 + 1
              ELSE IF (NBP.GE.5) THEN
                NSIDED = NSIDED + 1
              ENDIF
#endif
C
              CALL OUESTU
C             ===========
     &   (    NFECRA , NDIM   , NNOD ,
     &        IERROR ,
     &        ETTPA(IP,JXP)  , ETTPA(IP,JYP)  , ETTPA(IP,JZP)  ,
     &        ETTP(IP,JXP)   , ETTP(IP,JYP)   , ETTP(IP,JZP)   ,
     &        CDGFAC(1,IFAC)    , CDGFAC(2,IFAC)    , CDGFAC(3,IFAC)   ,
     &        XYZCEN(1,IEL)     , XYZCEN(2,IEL)     , XYZCEN(3,IEL)    ,
     &        ITYPFO , ICONFO , XYZNOD ,
     &        INDIAN )
C
#if DEBUG_LAGCEL
              IF (IMPLTG.EQ.1) WRITE(NFECRA,9004) IP , IFAC , INDIAN
#endif
C
C         ---> Elimination des particules qui posent problemes
C
              IF (IERROR.EQ.1) THEN
                IERROR = 0
                ITEPA(IP,JISOR) = 0
                NBPERR = NBPERR + 1
                DNBPER = DNBPER + TEPA(IP,JRPOI)
                IF (IERRIE.EQ.0) THEN
#if DEBUG_LAGCEL
                  IF (IMPLTG.EQ.1) WRITE(NFECRA,9101) IFAC,IP
#endif
                  GOTO 200
                ELSE
                  WRITE(NFECRA,9101) IFAC,IP
                  IERR = 1
                  GOTO 300
                ENDIF
C
C--> Si la particule passe dans la cellule voisine
C
              ELSE IF (INDIAN.EQ.1) THEN
C
C--> Si la particule IP est dans la cellule II alors le voisin ne
C    peut etre que la cellule JJ, et vice versa, et inversement, et
C    ainsi de suite.
C
                IFANEW = IFAC
C
                II = IFACEL(1,IFAC)
                JJ = IFACEL(2,IFAC)
C
                IF (IEL.EQ.II) THEN
                  ITEPA(IP,JISOR) = JJ
                ELSE IF (IEL.EQ.JJ) THEN
                  ITEPA(IP,JISOR) = II
                ENDIF
C
C Traitement de la periodicite (debut).
C  Ne marchera pas si on passe en parallele. Le fait d'être dans le halo
C  n'est pas suffisant pour conclure si la cellule est periodique.
C
                IF (ITEPA(IP,JISOR).GT.NCEL) THEN
C
C                 Si on est sur un plan de periodicite, on passe
C                 a l'ordre 1 sur les schemas
C
                  IBORD(IP) = -1
C
                  ITEPAS = ITEPA(IP,JISOR)
                  ITEPA(IP,JISOR) = IA(ICELCR+ITEPAS-NCEL-1)
C
C                 On recupere les informations sur la peridodicite
C
                  IPER  = IA(IPERCR+ITEPAS-NCEL-1)
C
C Faire un test si IPER != -1 pour ne traiter que les cellules periodiques
C Finir l'implémentation dans PERLOC
C
C                 POINT DE DEPART
C
                  PTA(1) = ETTPA(IP,JXP)
                  PTA(2) = ETTPA(IP,JYP)
                  PTA(3) = ETTPA(IP,JZP)
C
                  CALL LAGPER(IPER, PTA, PTB)
C                 ===========
C
                  ETTPA(IP,JXP) = PTB(1)
                  ETTPA(IP,JYP) = PTB(2)
                  ETTPA(IP,JZP) = PTB(3)
C
C                 POINT D'ARRIVEE
C
                  PTA(1) = ETTP(IP,JXP)
                  PTA(2) = ETTP(IP,JYP)
                  PTA(3) = ETTP(IP,JZP)
C
                  IPER  = IA(IPERCR+ITEPAS-NCEL-1)
C
                  CALL LAGPER(IPER, PTA, PTB)
C                 ===========
C
                  ETTP(IP,JXP) = PTB(1)
                  ETTP(IP,JYP) = PTB(2)
                  ETTP(IP,JZP) = PTB(3)
C
C                 MODIFICATION DES VITESSES PARTICULES
C
                  VECT(1) = ETTPA(IP,JUP)
                  VECT(2) = ETTPA(IP,JVP)
                  VECT(3) = ETTPA(IP,JWP)
C
                  CALL LAGVEC(IPER, VECT, VECTN)
C                 ===========
C
                  ETTPA(IP,JUP) = VECTN(1)
                  ETTPA(IP,JVP) = VECTN(2)
                  ETTPA(IP,JWP) = VECTN(3)
C
                  VECT(1) = ETTP(IP,JUP)
                  VECT(2) = ETTP(IP,JVP)
                  VECT(3) = ETTP(IP,JWP)
C
                  CALL LAGVEC(IPER, VECT, VECTN)
C                 ===========
C
                  ETTP(IP,JUP) = VECTN(1)
                  ETTP(IP,JVP) = VECTN(2)
                  ETTP(IP,JWP) = VECTN(3)
C
C                 MODIFICATION DES VITESSES FLUIDES VUES
C
                  VECT(1) = ETTPA(IP,JUF)
                  VECT(2) = ETTPA(IP,JVF)
                  VECT(3) = ETTPA(IP,JWF)
C
                  CALL LAGVEC(IPER, VECT, VECTN)
C                 ===========
C
                  ETTPA(IP,JUF) = VECTN(1)
                  ETTPA(IP,JVF) = VECTN(2)
                  ETTPA(IP,JWF) = VECTN(3)
C
                  VECT(1) = ETTP(IP,JUF)
                  VECT(2) = ETTP(IP,JVF)
                  VECT(3) = ETTP(IP,JWF)
C
                  CALL LAGVEC(IPER, VECT, VECTN)
C                 ===========
C
                  ETTP(IP,JUF) = VECTN(1)
                  ETTP(IP,JVF) = VECTN(2)
                  ETTP(IP,JWF) = VECTN(3)
C
                  IFANEW = 0
C
                ENDIF
C
C Traitement de la periodicite (fin)
C
C--> Retour pour balayage des faces de la cellule suivante
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
#if DEBUG_LAGCEL
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              WRITE(IMPLA4)
     &        NBP, -IFAC, IEL,
     &        (NODFBR(IN),IN=IPNFBR(IFAC),IPNFBR(IFAC+1)-1)
              IF (NBP.EQ.4) THEN
                NQUAD4 = NQUAD4 + 1
              ELSE IF (NBP.EQ.3) THEN
                NTRIA3 = NTRIA3 + 1
              ELSE IF (NBP.GE.5) THEN
                NSIDED = NSIDED + 1
              ENDIF
#endif
C
              CALL OUESTU
C             ===========
     &   (    NFECRA , NDIM   , NNOD ,
     &        IERROR ,
     &        ETTPA(IP,JXP)  , ETTPA(IP,JYP)  , ETTPA(IP,JZP)  ,
     &        ETTP(IP,JXP)   , ETTP(IP,JYP)   , ETTP(IP,JZP)   ,
     &        CDGFBO(1,IFAC)    , CDGFBO(2,IFAC)    , CDGFBO(3,IFAC)   ,
     &        XYZCEN(1,IEL)     , XYZCEN(2,IEL)     , XYZCEN(3,IEL)    ,
     &        ITYPFO , ICONFO , XYZNOD ,
     &        INDIAN )
C
#if DEBUG_LAGCEL
              IF (IMPLTG.EQ.1) WRITE(NFECRA,9002) IP , IFAC , INDIAN
#endif
C
C         ---> Elimination des particules qui posent problemes
C
              IF (IERROR.EQ.1) THEN
                IERROR = 0
                ITEPA(IP,JISOR) = 0
                NBPERR = NBPERR + 1
                DNBPER = DNBPER + TEPA(IP,JRPOI)
                IF (IERRIE.EQ.0) THEN
#if DEBUG_LAGCEL
                  IF (IMPLTG.EQ.1) WRITE(NFECRA,9101) IFAC,IP
#endif
                  GOTO 200
                ELSE
                  WRITE(NFECRA,9101) IFAC,IP
                  IERR = 1
                  GOTO 300
                ENDIF
C
C--> Si la trajectoire de la particule traverse la face de bord
C
              ELSE IF (INDIAN.EQ.1) THEN
                IF (NORDRE.EQ.2) IBORD(IP) = IFAC
C
#if DEBUG_LAGCEL
                IF (IMPLTG.EQ.1) WRITE(NFECRA,9003) IFAC
#endif
C
C--> Traitement de l'interaction particule/frontiere
C
C   1) modification du numero de la cellule de depart
C      (ITEPA(IP,JISOR) = IEL ou 0)
C
C   2) P devient K, intersection entre la rayon PQ et le plan
C      de la face de bord
C
C   3) Q est a determiner selon la nature de l'interaction
C      particule/frontiere
C
C   resultat : ISUIVI = 0 -> on ne suit plus la particule
C   ~~~~~~~~                 dans le maillage apres USLABO
C              ISUIVI = 1 -> la particule continue a etre suivie
C              ISUIVI = -999 valeur initiale aberrante
C
C   Blindage dans USLABO : ISUIVI ne peut que valoir 0 ou 1 apres.
C
                CALL USLABO
C               ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   IFAC   , IP     , ISUIVI ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , IFRLAG , ITEPA  , INDEP  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   SURFBN , DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , PARBOR , ETTP(1,JUP) ,
     &                                       ETTP(1,JUF) , AUXL   ,
     &   RDEVEL , RTUSER , RA     )
C
C--> Si la particule continue sa route (ex : rebond) apres l'interaction
C    avec la frontiere, il faut continuer a la suivre donc GOTO 100
C
                IF (ISUIVI.EQ.1) THEN
                  IFANEW = -IFAC
                  GOTO 100
                ELSE IF (ISUIVI.NE.0 .AND. ISUIVI.NE.1) THEN
                  WRITE (NFECRA,8010) IFRLAG(IFAC) , ISUIVI
                  CALL CSEXIT (1)
C                 ===========
                ENDIF
C
              ENDIF
C
C  fin de IF (IFAC.GT.0 .AND. IFAC.NE.IFAOLD) THEN
            ENDIF
C
C fin de DO WHILE
          ENDDO
C
C--> Fin de la boucle sur les faces entourant la cellule courante
C
 200      CONTINUE
C
#if DEBUG_LAGCEL
          IF (IMPLTG.EQ.1) WRITE(NFECRA,9005) IP , ITEPA(IP,JISOR)
#endif
C
C-->Fin de la boucle principale sur les particules
C
        ENDIF
C
      ENDDO
C
      RETURN
C
C
C=======================================================================
C 3. Gestion des erreurs et ecriture des fichiers debugX
C=======================================================================
C
C--> debugX.CASE sont une aide graphique au debuggage, ils contiennent :
C   1) le dernier morceau de la trajectoire de la particule
C      qui a un pepin,
C   2) les faces (avec leur numero global) entourant les volumes
C      qui contiennent le dernier morceau de trajectoire,
C   3) le fichier debug2.geom contient en plus les CDG faces et cellules
C      qui ont pour numero d'element Ensight celui de l'element
C      de maillage (face ou cellule) associe.
C
C--> On ecrit 4 fichiers :
C   1) debug1.geom et debug1.CASE ecris au format ensight,
C   2) debug2.geom et debug2.CASE ecris au format ensight gold.
C   Les deux fichiers .CASE peuvent etre lus avec Ensight7 indifferemment.
C
C--> L'ecriture des fichiers debugX est declenchee par defaut
C    par une erreur sur la detection de la trajectoire (IERR=1).
C
C   ATTENTION : seules les faces de forme triangles et quadrangles sont
C   =========   reconnues par le format ensight et donc enregistrees
C               dans debug1.geom !
C               Le format ensight gold reconnait toutes formes de faces
C               et doit etre utilise par defaut...
C
C
 300  CONTINUE
C
#if DEBUG_LAGCEL
C
      IF (IERR.EQ.1) THEN
C
C
C...FORMAT ENSIGHT
C
C
            WRITE(NFECRA,9050)
            OPEN (IMPLA1,FILE='debug1.geom',
     &            STATUS='UNKNOWN',FORM='FORMATTED',
     &            ACCESS='SEQUENTIAL')
C
            WRITE(IMPLA1,'(A)') 'geometrie debuggage'
            WRITE(IMPLA1,'(A)') 'au format ensight6'
            WRITE(IMPLA1,'(A)') 'node id given'
            WRITE(IMPLA1,'(A)') 'element id given'
            WRITE(IMPLA1,'(A)') 'coordinates'
            WRITE(IMPLA1,'(I8)') 4*NQUAD4 + 3*NTRIA3 + 2
C
            KPT = 0
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
              IF (NBP.LT.5) THEN
                DO IN = 1,NBP
                  KPT = MAX(KPT,ICONFO(IN))
                  WRITE(IMPLA1,'(I8,3E12.5)') ICONFO(IN),
     &                                        XYZNOD(1,ICONFO(IN)),
     &                                        XYZNOD(2,ICONFO(IN)),
     &                                        XYZNOD(3,ICONFO(IN))
                ENDDO
              ENDIF
            ENDDO
C
            WRITE(IMPLA1,'(I8,3E12.5)') KPT + 1,
     &                                  ETTPA(IP,JXP),
     &                                  ETTPA(IP,JYP),
     &                                  ETTPA(IP,JZP)
C
            WRITE(IMPLA1,'(I8,3E12.5)') KPT + 2,
     &                                  ETTP(IP,JXP),
     &                                  ETTP(IP,JYP),
     &                                  ETTP(IP,JZP)
C
            WRITE(IMPLA1,'(A)')    'part 1'
            WRITE(IMPLA1,'(A,I9)') 'detail particule ',IP
            WRITE(IMPLA1,'(A)')    'bar2'
            NPT = 1
            WRITE(IMPLA1,'(I8)')   NPT
            WRITE(IMPLA1,'(3I8)')  NPT , KPT + 1 , KPT + 2
C
            WRITE(IMPLA1,'(A)') 'part 2'
            WRITE(IMPLA1,'(A)') 'faces'
C
            IF (NTRIA3.GT.0) THEN
              WRITE(IMPLA1,'(A)')  'tria3'
              WRITE(IMPLA1,'(I8)') NTRIA3
            ENDIF
C
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
              IF (IFAC.LT.0) IFAC = -IFAC
              IF (NBP.EQ.3) THEN
                WRITE(IMPLA1,'(4I8)')
     &          IFAC , ( ICONFO(IN) , IN=1,NBP )
              ENDIF
            ENDDO
C
            IF (NQUAD4.GT.0) THEN
              WRITE(IMPLA1,'(A)')  'quad4'
              WRITE(IMPLA1,'(I8)') NQUAD4
            ENDIF
C
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
              IF (IFAC.LT.0) IFAC = -IFAC
              IF (NBP.EQ.4) THEN
                WRITE(IMPLA1,'(5I8)')
     &          IFAC , ( ICONFO(IN) , IN=1,NBP )
              ENDIF
            ENDDO
C
            CLOSE(IMPLA1)
C

            WRITE(NFECRA,9055)
            OPEN (IMPLA1,FILE='debug1.CASE',
     &            STATUS='UNKNOWN',FORM='FORMATTED',
     &            ACCESS='SEQUENTIAL')
              WRITE(IMPLA1,'(A)') 'FORMAT'
              WRITE(IMPLA1,'(A)') 'type:     ensight'
              WRITE(IMPLA1,'(A)') 'GEOMETRY'
              WRITE(IMPLA1,'(A)') 'model:    debug1.geom'
            CLOSE(IMPLA1)
C
C
C...FORMAT ENSIGHT GOLD
C
C-> Ouverture
C
            WRITE(NFECRA,9060)
            OPEN (IMPLA1,FILE='debug2.geom',
     &            STATUS='UNKNOWN',FORM='FORMATTED',
     &            ACCESS='SEQUENTIAL')
C
C-> Entete
C
            WRITE(IMPLA1,'(A)') 'geometrie debuggage'
            WRITE(IMPLA1,'(A)') 'au format ensight gold'
            WRITE(IMPLA1,'(A)') 'node id given'
            WRITE(IMPLA1,'(A)') 'element id given'
C
C-> Id des points du segment de trajectoire
C
            KPT = 0
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
              DO IN = 1,NBP
                KPT = MAX(KPT,ICONFO(IN))
              ENDDO
            ENDDO
C
C-> Part de la trajectoire
C
            IPART = 1
            WRITE(IMPLA1,'(A)')     'part'
            WRITE(IMPLA1,'(I10)')   IPART
            WRITE(IMPLA1,'(A,I9)')  'detail particule ',IP
            WRITE(IMPLA1,'(A)')     'coordinates'
            NPT = 2
            WRITE(IMPLA1,'(I10)') NPT
            WRITE(IMPLA1,'(I10)') KPT + 1
            WRITE(IMPLA1,'(I10)') KPT + 2
            DO IN   = 0,2
              WRITE(IMPLA1,'(E12.5)') ETTPA(IP,JXP+IN)
              WRITE(IMPLA1,'(E12.5)') ETTP(IP,JXP+IN)
            ENDDO
C
            WRITE(IMPLA1,'(A)')    'bar2'
            NPT = 1
            KPT = 2
            WRITE(IMPLA1,'(I10)')  NPT
            WRITE(IMPLA1,'(I10)')  IP
            WRITE(IMPLA1,'(2I10)') NPT,KPT
C
C-> Part des triangles
C
            IF (NTRIA3.GT.0) THEN
              IPART = IPART+1
              WRITE(IMPLA1,'(A)')    'part'
              WRITE(IMPLA1,'(I10)')  IPART
              WRITE(IMPLA1,'(A)')    'faces triangles'
              WRITE(IMPLA1,'(A)')    'coordinates'
              WRITE(IMPLA1,'(I10)')  NTRIA3*3
C
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.EQ.3) THEN
                  DO IN = 1,NBP
                    WRITE(IMPLA1,'(I10)') ICONFO(IN)
                  ENDDO
                ENDIF
              ENDDO
              DO IN   = 1,3
                REWIND(IMPLA4)
                DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                  READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
                  IF (NBP.EQ.3) THEN
                    DO II = 1,NBP
                      WRITE(IMPLA1,'(E12.5)') XYZNOD(IN,ICONFO(II))
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
C
              WRITE(IMPLA1,'(A)')  'tria3'
              WRITE(IMPLA1,'(I10)') NTRIA3
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (IFAC.LT.0) IFAC = -IFAC
                IF (NBP.EQ.3) WRITE(IMPLA1,'(I10)') IFAC
              ENDDO
              KPT = 0
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.EQ.3) THEN
                  WRITE(IMPLA1,'(3I10)') KPT+1,KPT+2,KPT+3
                  KPT = KPT+3
                ENDIF
              ENDDO
            ENDIF
C
C-> Part des quadrangles
C
            IF (NQUAD4.GT.0) THEN
              IPART = IPART+1
              WRITE(IMPLA1,'(A)')    'part'
              WRITE(IMPLA1,'(I10)')  IPART
              WRITE(IMPLA1,'(A)')    'faces quadrangles'
              WRITE(IMPLA1,'(A)')    'coordinates'
              WRITE(IMPLA1,'(I10)')  NQUAD4*4
C
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.EQ.4) THEN
                  DO IN = 1,NBP
                    WRITE(IMPLA1,'(I10)') ICONFO(IN)
                  ENDDO
                ENDIF
              ENDDO
              DO IN   = 1,3
                REWIND(IMPLA4)
                DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                  READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
                  IF (NBP.EQ.4) THEN
                    DO JJ = 1,NBP
                      WRITE(IMPLA1,'(E12.5)') XYZNOD(IN,ICONFO(JJ))
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
C
              WRITE(IMPLA1,'(A)')  'quad4'
              WRITE(IMPLA1,'(I10)') NQUAD4
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (IFAC.LT.0) IFAC = -IFAC
                IF (NBP.EQ.4) WRITE(IMPLA1,'(I10)') IFAC
              ENDDO
              KPT = 0
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.EQ.4) THEN
                  WRITE(IMPLA1,'(4I10)') KPT+1,KPT+2,KPT+3,KPT+4
                  KPT = KPT+4
                ENDIF
              ENDDO
            ENDIF
C
C-> Part des Polygones
C
            IF (NSIDED.GT.0) THEN
              IPART = IPART+1
              WRITE(IMPLA1,'(A)')    'part'
              WRITE(IMPLA1,'(I10)')  IPART
              WRITE(IMPLA1,'(A)')    'faces polygones'
              WRITE(IMPLA1,'(A)')    'coordinates'
C
              KPT = 0
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.GE.5) KPT = KPT + NBP
              ENDDO
              WRITE(IMPLA1,'(I10)')  KPT
C
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.GE.5) THEN
                  DO IN = 1,NBP
                    WRITE(IMPLA1,'(I10)') ICONFO(IN)
                  ENDDO
                ENDIF
              ENDDO
              DO IN   = 1,3
                REWIND(IMPLA4)
                DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                  READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
                  IF (NBP.GE.5) THEN
                    DO JJ = 1,NBP
                      WRITE(IMPLA1,'(E12.5)') XYZNOD(IN,ICONFO(JJ))
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
C
              WRITE(IMPLA1,'(A)')  'nsided'
              WRITE(IMPLA1,'(I10)') NSIDED
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (IFAC.LT.0) IFAC = -IFAC
                IF (NBP.GE.5) WRITE(IMPLA1,'(I10)') IFAC
              ENDDO
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.GE.5) WRITE(IMPLA1,'(I10)') NBP
              ENDDO
              KPT = 0
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(IN),IN = 1,NBP)
                IF (NBP.GE.5) THEN
                  WRITE(IMPLA1,'(100I10)') (KPT+II,II = 1,NBP)
                  KPT = KPT+NBP
                ENDIF
              ENDDO
            ENDIF
C
C-> Part des points supplementaires (CDG faces et CDG cellules)
C
C
            IPART = IPART + 1
            WRITE(IMPLA1,'(A)')   'part'
            WRITE(IMPLA1,'(I10)') IPART
            WRITE(IMPLA1,'(A)')   'CDG faces et CDG cellules'
            WRITE(IMPLA1,'(A)')   'coordinates'
C
            KPT    = 0
            IELOLD = 0
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
              IF (IEL.NE.IELOLD) THEN
                IELOLD  = IEL
                KPT     = KPT + 1
              ENDIF
            ENDDO
            WRITE(IMPLA1,'(I10)') NQUAD4 + NTRIA3 + NSIDED + KPT
C
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
              IF (IFAC.GT.0) THEN
                WRITE(IMPLA1,'(I10)') IFAC
              ELSE IF (IFAC.LT.0) THEN
                WRITE(IMPLA1,'(I10)') -IFAC
              ENDIF
            ENDDO
C
            IELOLD  = 0
            REWIND(IMPLA4)
            DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
              READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
              IF (IEL.NE.IELOLD) THEN
                IELOLD = IEL
                WRITE(IMPLA1,'(I10)') IEL
              ENDIF
            ENDDO
C
            DO IN = 1,3
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
                IF (IFAC.GT.0) THEN
                  WRITE(IMPLA1,'(E12.5)') CDGFAC(IN,IFAC)
                ELSE IF (IFAC.LT.0) THEN
                  WRITE(IMPLA1,'(E12.5)') CDGFBO(IN,-IFAC)
                ENDIF
              ENDDO
              IELOLD  = 0
              REWIND(IMPLA4)
              DO IPT = 1, NQUAD4 + NTRIA3 + NSIDED
                READ(IMPLA4) NBP, IFAC, IEL, (ICONFO(II),II = 1,NBP)
                IF (IEL.NE.IELOLD) THEN
                  IELOLD = IEL
                  WRITE(IMPLA1,'(E12.5)') XYZCEN(IN,IEL)
                ENDIF
              ENDDO
            ENDDO
C
C-> Fermeture du fichier geometrique
C
            CLOSE(IMPLA1)
C
C-> Ecriture du fichier CASE
C
            WRITE(NFECRA,9065)
            OPEN (IMPLA1,FILE='debug2.CASE',
     &            STATUS='UNKNOWN',FORM='FORMATTED',
     &            ACCESS='SEQUENTIAL')
              WRITE(IMPLA1,'(A)') 'FORMAT'
              WRITE(IMPLA1,'(A)') 'type:     ensight gold'
              WRITE(IMPLA1,'(A)') 'GEOMETRY'
              WRITE(IMPLA1,'(A)') 'model:    debug2.geom'
            CLOSE(IMPLA1)
C
      ENDIF
C
      CLOSE(IMPLA4,STATUS='DELETE')
C
      IF (IERR.EQ.1) THEN
        WRITE(NFECRA,9100)
        RETURN
      ENDIF
C
      WRITE(NFECRA,9999) IERRIE , IERR
      CALL CSEXIT (1)
C     ===========
C
#endif
C
C--------
C FORMATS
C--------
C
 8010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    L''INDICATEUR DE SUIVI DE PARTICULE DANS LE MAILLAGE    ',/,
     &'@       APRES INTERACTION AVEC LA FRONTIERE NB = ',I10        ,/,
     &'@       A UNE VALEUR NON PERMISE (LAGCEL).                   ',/,
     &'@                                                            ',/,
     &'@   ISUIVI DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1              ',/,
     &'@       IL VAUT ICI ISUIVI = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ISUIVI dans la subroutine USLABO.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#if DEBUG_LAGCEL
 9001 FORMAT(3X,'@@ Particule N°',I6,' >> cellule de depart : ',I7)
C
 9002 FORMAT(3X,'@@ Particule N°',I6,'  > face de bord : ',I7
     &           ,' REPERE = ',I2)
C
 9003 FORMAT(3X,'@@ Interaction frontiere -> face de bord : ',I7)
C
 9004 FORMAT(3X,'@@ Particule N°',I6,'  > face interne : ',I7
     &           ,' REPERE = ',I2)
C
 9005 FORMAT(3X,'@@ Particule N°',I6,' >> cellule d''arrive : ',I7,/)
C
C
 9050 FORMAT(/3X,'** ECRITURE DU FICHIER debug1.geom ',
     &        /3X,'   AU FORMAT ENSIGHT6')
C
 9055 FORMAT(/3X,'** ECRITURE DU FICHIER debug1.CASE ',
     &        /3X,'   AU FORMAT ENSIGHT6')
C
 9060 FORMAT(/3X,'** ECRITURE DU FICHIER debug2.geom ',
     &        /3X,'   AU FORMAT ENSIGHT GOLD')
C
 9065 FORMAT(/3X,'** ECRITURE DU FICHIER debug2.CASE ',
     &       /3X,'   AU FORMAT ENSIGHT GOLD',/)
C
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE REPERAGE D''UNE PARTICULE     ',/,
     &'@    *********    (LAGCEL)                                   ',/,
     &'@                                                            ',/,
     &'@  Les fichiers debug1.CASE et debug2.CASE peuvent fournir   ',/,
     &'@    une explication de la cause de cet echec de             ',/,
     &'@    l''algorithme.                                          ',/,
     &'@  La cause la plus probable est la presence d''une face     ',/,
     &'@    non plane dans le maillage (due a un recollement        ',/,
     &'@    par exemple).                                           ',/,
     &'@                                                            ',/,
     &'@  Choisir de preference le fichier debug2.CASE              ',/,
     &'@    au format ensight gold, celui-ci contient plus          ',/,
     &'@    d''inforamtion que debug1.CASE au format ensight6.      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute en mode sans erreur.       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
#endif
C
 9101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE REPERAGE D''UNE PARTICULE     ',/,
     &'@    *********    (LAGCEL)                                   ',/,
     &'@                                                            ',/,
     &'@  ECHEC DE REPERAGE SUR LA FACE ',I10                        ,/,
     &'@    POUR LA PARTICULE ',I10                                  ,/,
     &'@                                                            ',/,
     &'@  Explication possible : le reperage echoue lorsque la      ',/,
     &'@    position d''arrive de la particule fait partie du       ',/,
     &'@    plan de la face.                                        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute en mode sans erreur.       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9102 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE REPERAGE D''UNE PARTICULE     ',/,
     &'@    *********    (LAGCEL)                                   ',/,
     &'@                                                            ',/,
     &'@  AUCUNE FACE ENTOURANT LA CELLULE ',I10                     ,/,
     &'@    N''EST TRAVERSEE PAR LA DROITE PASSANT PAR LES POINTS   ',/,
     &'@    DE DEPART ET D''ARRIVEE DE LA PARTICULE ',I10            ,/,
     &'@                                                            ',/,
     &'@  Explication possible : mauvais traitement d''une          ',/,
     &'@    cellule concave (amelioration a venir).                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute en mode sans erreur.       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9103 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE REPERAGE D''UNE PARTICULE     ',/,
     &'@    *********    (LAGCEL)                                   ',/,
     &'@                                                            ',/,
     &'@  ECHEC DE REPERAGE POUR LA PARTICULE ',I10                  ,/,
     &'@    TROP DE CELLULES PARCOURUES                             ',/,
     &'@                                                            ',/,
     &'@  Explication possible : mauvais traitement d''une          ',/,
     &'@    cellule concave (amelioration a venir).                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute en mode sans erreur.       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#if DEBUG_LAGCEL
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE REPERAGE D''UNE PARTICULE     ',/,
     &'@    *********    (LAGCEL)                                   ',/,
     &'@                                                            ',/,
     &'@  LA DETECTION DE LA TRAJECTOIRE DE LA PARTICULE EST SORTIE ',/,
     &'@    EN ERREUR (MODE SANS ERREUR ENCLENCHE AVEC              ',/,
     &'@    IERRIE = ',I10,'),                                      ',/,
     &'@    SANS QUE L''INDICATEUR D''ECRITURE DES FICHIERS DEBUG   ',/,
     &'@    AIT ETE ACTIVE (IERR = ',I10,').                        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
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
