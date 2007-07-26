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
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  GESTION DES STRUCTURES MOBILES EN ALE AVEC COUPLAGE INTERNE
CFONC
CFONC  ON TROUVERA ICI DEUX ROUTINES DIFFERENTES :
CFONC
CFONC  - USSTR1 : APPELEE A L'INITIALISATION, POUR DEFINIR LES STRUCTURES
CFONC                ET LEURS PARAMETRES INITIAUX (VITESSE, DEPLACEMENT)
CFONC
CFONC  - USSTR2 : APPELE A CHAQUE PAS DE TEMPS POUR DEFINIR LES
CFONC                CARACTERISTIQUES (POTENTIELLEMENT VARIABLES) DES
CFONC                STRUCTURES
CFONC
c@fonce
C-----------------------------------------------------------------------
c@a
                        SUBROUTINE USSTR1
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IDFSTR ,
     &   IDEVEL , ITUSER , IA     ,
     &   AEXXST , BEXXST , CFOPRE ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   XSTR0  , VSTR0  , XSTREQ ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------

C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  GESTION DES STRUCTURES MOBILES EN ALE AVEC COUPLAGE INTERNE
CFONC
CFONC    DEFINITION DES STRUCTURES
CFONC    DEPLACEMENT ET VITESSE INITIAUX
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
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
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
CARGU ! IDFSTR(NFABOR! TE !  -> ! DEFINITION DES STRUCTURES            !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! AEXXST,BEXXST! R  !  -> ! COEFFICIENTS DE PREDICTION DU DEPLACT!
CARGU ! CFOPRE       ! R  !  -> ! COEFF DE PREDICTION DES EFFORTS      !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! XSTR0(NDIM,  ! TR !  -> ! DEPLACEMENT INITIAL DES STRUCTURES   !
CARGU !       NBSTRU)!    !     !                                      !
CARGU ! VSTR0(NDIM,  ! TR !  -> ! VITESSE INITIALE DES STRUCTURES      !
CARGU !       NBSTRU)!    !     !                                      !
CARGU ! XSTREQ(NDIM, ! TR !  -> ! DEPLACEMENT DU MAILLAGE INITIAL PAR  !
CARGU !       NBSTRU)!    !     ! RAPPORT A L'EQUILIBRE                !
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
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NBSTRU
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDFSTR(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION AEXXST, BEXXST, CFOPRE
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION XSTR0(3,NSTRMX), XSTREQ(3,NSTRMX)
      DOUBLE PRECISION VSTR0(3,NSTRMX)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IFML, ICOUL
      INTEGER          ILELT, NLELT
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C
C
      IF(1.EQ.1) THEN
        NBSTRU = 0
        RETURN
      ENDIF
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  INITIALISATIONS
C
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2.  DEFINITION DES STRUCTURES (appel a l'initialisation)
C=======================================================================
C
C     On remplit le tableau IDFSTR(NFABOR)
C     IDFSTR(IFAC) est le numero de la structure a laquelle appartient
C       la face de bord IFAC (0 si elle n'appartient a aucune structure)
C     Le nombre de structure est automatiquement determine a partir du
C       plus grand element de IDFSTR (les numeros des structures doivent
C       donc etre affectes de maniere sequentielle sans trou en commencant
C       par 1).

C Dans l'exemple ci-dessous la structure 1 est bordee par les faces de
C   couleur 4, la structure 2 par les faces de couleur 6
C
C
      CALL GETFBR('4',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IDFSTR(IFAC) = 1
C
      ENDDO
C
      CALL GETFBR('6',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IDFSTR(IFAC) = 2
C
      ENDDO
C
C
C Pour chaque structure, on definit si necessaire
C - une vitesse initiale VSTR0
C - un deplacement initial XSTR0 (i.e. la valeur du deplacement XSTR
C     au temps t=0 par rapport au maillage initial)
C - un deplacement par rapport a l'equilibre XSTREQ (i.e. le deplacement
C     du maillage initial par rapport a la position d'equilibre de la
C     structure ; la force de retour exercee par le ressort a un temps
C     donne pour un deplacement XSTR sera donc -k*(XSTR+XSTREQ) ).
C Toutes les composantes de XSTR0, XSTREQ et VSTR0 sont initialisees a 0
C
C     Exemple : deplacement initial y=2 pour la structure 1
C               deplacement par rapport a l'equilibre yeq=1 pour la
C                  structure 1
C               vitesse initiale uz=0.5 pour la structure 2
C
C Dans le cas d'un calcul initial, ou d'une suite d'un calcul sans ALE,
C   une iteration 0 est automatiquement realisee pour gerer un eventuel
C   deplacement initial des structures. Si necessaire, positionner
C   ITALIN a 1 dans usalin pour activer une iteration 0 dans les autres
C   cas.
C
      XSTR0(2,1)  = 2.D0
      XSTREQ(2,1) = 1.D0
      VSTR0(3,2)  =-0.5D0
C
C --- Si necessaire on definit les coefficients d'extrapolation utiles
C       en couplage explicite :
C       deplacement predit = X(n) + AEXXST.DT.X'(n)
C                                 + BEXXST.DT.( X'(n)-X'(n-1) )
C       force envoyee a la structure = CFOPRE.F(n) + (1.D0-CFOPRE).F(n-1)
C
      AEXXST =  0.5D0
      BEXXST =  0.0D0
      CFOPRE =  2.D0
C
C --- Ecriture des fichiers historiques des structures mobiles
C     (deplacement, vitesse, acceleration, force)
C     La periodicite de sortie est la meme que pour les historiques
C     standards (NTHIST)
      IHISTR = 1
C
      RETURN
C
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USSTR2
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NBSTRU ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IDFSTR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DTCEL  ,
     &   XMSTRU , XCSTRU , XKSTRU , XSTREQ , XSTR   , VSTR   , FORSTR ,
     &   DTSTR  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     -------------------------------------------------------------

C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  GESTION DES STRUCTURES MOBILES EN ALE AVEC COUPLAGE INTERNE
CFONC
CFONC  SPECIFICATION DES CARACTERISTIQUES DES STRUCTURES
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NBSTRU       ! E  !  -> ! NOMBRE DE STRUCTURES DEFINIES        !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
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
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
CARGU ! IDFSTR(NFABOR! TE !  -> ! DEFINITION DES STRUCTURES            !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DTCEL(NCELET)! TR !  -> ! PAS DE TEMPS DANS LES CELLULES       !
CARGU ! XMSTRU(NDIM, ! TR ! <-  ! MATRICE DE MASSE DES STRUCTURES      !
CARGU !  NDIM,NBSTRU)!    !     !                                      !
CARGU ! XCSTRU(NDIM, ! TR ! <-  ! MATRICE DE FRICTION DES STRUCTURES   !
CARGU !  NDIM,NBSTRU)!    !     !                                      !
CARGU ! XKSTRU(NDIM, ! TR ! <-  ! MATRICE DE RAIDEUR DES STRUCTURES    !
CARGU !  NDIM,NBSTRU)!    !     !                                      !
CARGU ! XSTREQ(NDIM, ! TR !  -> ! DEPLACEMENT DU MAILLAGE INITIAL PAR  !
CARGU !       NBSTRU)!    !     ! RAPPORT A L'EQUILIBRE                !
CARGU ! XSTR(NDIM,   ! TR !  -> ! DEPLACEMENT DES STRUCTURES           !
CARGU !       NBSTRU)!    !     !                                      !
CARGU ! VSTR(NDIM,   ! TR !  -> ! VITESSE  DES STRUCTURES              !
CARGU !       NBSTRU)!    !     !                                      !
CARGU ! FORSTR(NDIM  ! TR ! <-> ! EFFORT SUR LES STRUCTURES (CONTIENT  !
CARGU !       NBSTRU)!    !     !            LES EFFORTS DUS AU FLUIDE)!
CARGU ! DTSTR(NBSTRU)! TR ! <-  ! PAS DE TEMPS DES STRUCTURES          !
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
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NBSTRU
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDFSTR(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DTCEL(NCELET)
      DOUBLE PRECISION XMSTRU(3,3,NSTRMX)
      DOUBLE PRECISION XCSTRU(3,3,NSTRMX)
      DOUBLE PRECISION XKSTRU(3,3,NSTRMX)
      DOUBLE PRECISION XSTREQ(3,NSTRMX)
      DOUBLE PRECISION XSTR(3,NSTRMX)
      DOUBLE PRECISION VSTR(3,NSTRMX)
      DOUBLE PRECISION FORSTR(3,NSTRMX)
      DOUBLE PRECISION DTSTR(NSTRMX)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          II, JJ, ISTR
      DOUBLE PRECISION THETA, SINT, COST, XM, XC, XK, FX, FY, FZ
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C
C
      IF(1.EQ.1) RETURN
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  INITIALISATIONS
C
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2.  COEFFICIENTS DES STRUCTURES (appel a chaque pas de temps)
C=======================================================================
C
C
C     On remplit ici les coefficients de definissant la structure.
C      - sa masse M                     (XMSTRU)
C      - son coefficient de friction C  (XCSTRU)
C      - son coefficient de raideur K   (XKSTRU)
C
C     Le tableau FORSTR contient les efforts exerces par le fluide
C       sur chacune des structures. Il est possible d'y ajouter une
C       composante de force exterieure (gravite par exemple)
C
C     Le tableau XSTR contient le deplacement des structures par rapport
C       au maillage initial
C     Le tableau XSTR0 contient le deplacement des structures dans le
C       maillage initial par rapport a leur position d'equilibre
C     Le tableau XPSTR contient la vitesse des structures
C
C     XSTR, XSTR0 et VSTR sont des DONNEES pouvant servir eventuellement
C     a calculer M, C et K. ILS NE DOIVENT PAS ETRE MODIFIES.
C
C     L'evolution du systeme est resolue de maniere tridimensionnelle,
C       ces trois coefficients sont donc en fait des matrices 3x3.
C
C     L'equation resolue est
C
C       M.X'' + C.X' + K.(X+X0) = F
C       = -     = -    =  - --    -
C
C       (X est le deplacement par rapport a la position initiale du maillage,
C        X0 est le deplacement de la position dans le maillage initial par
C        rapport a l aposition d'equilibre)
C
C     La resolution est effectuee par la methode de Newmark HHT.
C     Le pas de temps utilise peut etre different du pas de temps
C       fluide (a definir dans le tableau DTSTR, initialise par defaut
C       au pas de temps fluide).
C
C
C     On met a zero tous les coefficients
      DO ISTR = 1, NBSTRU
C
        DO II = 1, 3
          DO JJ = 1, 3
            XMSTRU(II,JJ,ISTR) = 0.D0
            XCSTRU(II,JJ,ISTR) = 0.D0
            XKSTRU(II,JJ,ISTR) = 0.D0
          ENDDO
        ENDDO
C
      ENDDO
C
C Dans l'exemple ci-dessous, la structure 1, de masse 5 kg est retenue par
C   ressort isotrope de raideur 2 N/m et de coefficient de friction 3 kg.s
C
      DO II = 1, 3
        XMSTRU(II,II,1) = 5.D0
        XCSTRU(II,II,1) = 2.D0
        XKSTRU(II,II,1) = 3.D0
      ENDDO
C
C
C Dans l'exemple ci-dessous, la structure 2 est contrainte a un mouvement
C   decompose en deux :
C  - dans le plan xOy, le mouvement est force dans une direction X, avec
C    une masse xm, une friction xc et une raideur xk (et la composante
C    normale Y est donc forcee a 0). L'axe X est incline d'un angle THETA
C    par rapport a l'axe x du repere global.
C  - dans la direction z, le mouvement est un mouvement d'oscillation
C    harmonique de masse 1 et de raideur 1 (et de friction nulle) avec un
C    forcage externe en 3.cos(4.t) (en plus des efforts fluides).
C
      THETA = PI/6.D0
      COST = COS(THETA)
      SINT = SIN(THETA)
C
C
C Dans le repere local on a donc
C      xm.X'' + xc.X' + xk.X = FX
C                          Y = 0
C         Z''         +    Z = FZ + 3.cos(4.t)
C
C   FX, FY et FZ sont les composantes des efforts fluides dans le repere
C     local. Soit, a partir des composantes dans le repere global :
C     FX =  COST*Fx + SINT*Fy
C     FY = -SINT*Fx + COST*Fy
C     FZ = Fz
C
C Apres changement de repere, on obtient donc :
C
      XM = 1.D0
      XC = 3.D-1
      XK = 2.D0
      FX = FORSTR(1,2)
      FY = FORSTR(2,2)
C
      XMSTRU(1,1,2) = XM*COST**2
      XMSTRU(1,2,2) = XM*COST*SINT
      XMSTRU(2,1,2) = XM*COST*SINT
      XMSTRU(2,2,2) = XM*SINT**2
      XMSTRU(3,3,2) = 1.D0
C
      XCSTRU(1,1,2) = XC*COST**2
      XCSTRU(1,2,2) = XC*COST*SINT
      XCSTRU(2,1,2) = XC*COST*SINT
      XCSTRU(2,2,2) = XC*SINT**2
C
      XKSTRU(1,1,2) = (XK-1.D0)*COST**2 + 1.D0
      XKSTRU(1,2,2) = (XK-1.D0)*COST*SINT
      XKSTRU(2,1,2) = (XK-1.D0)*COST*SINT
      XKSTRU(2,2,2) = (XK-1.D0)*SINT**2 + 1.D0
      XKSTRU(3,3,2) = 1.D0
C
      FORSTR(1,2) = FX*COST**2   + FY*SINT*COST
      FORSTR(2,2) = FX*SINT*COST + FY*SINT**2
      FORSTR(3,2) = FORSTR(3,2) + 3.D0*COS(4.D0*TTCABS)
C
      DO ISTR = 1, NBSTRU
        DTSTR(ISTR) = DTCEL(1)
      ENDDO
C
C
      RETURN
C
      END
c@z
