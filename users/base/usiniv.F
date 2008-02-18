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
                        SUBROUTINE USINIV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB , COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE UTILISATEUR : INITIALISATION DES VARIABLES DE CALCUL
CFONC
CFONC  Cette routine est appelee en debut de calcul (suite ou non)
CFONC      avant le debut de la boucle en temps
CFONC
CFONC  Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
CFONC      les variables de calcul,
CFONC      les valeurs du pas de temps
CFONC
CFONC
CFONC  On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
CFONC      ou relues d'un fichier suite
CFONC  On ne dispose des variables VISCLS, CP (quand elles sont
CFONC      definies) que si elles ont pu etre relues dans un fichier
CFONC      suite de calcul
CFONC
CFONC  Les proprietes physiaues sont accessibles dans le tableau
CFONC      PROPCE (prop au centre), PROPFA (aux faces internes),
CFONC      PROPFB (prop aux faces de bord)
CFONC      Ainsi,
CFONC       PROPCE(IEL,IPPROC(IROM  (IPHAS))) designe ROM   (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(IVISCL(IPHAS))) designe VISCL (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(ICP   (IPHAS))) designe CP    (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)
CFONC
CFONC       PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)
CFONC
CFONC       PROPFB(IFAC,IPPROB(IROM  (IPHAS))) designe ROMB  (IFAC,IPHAS)
CFONC       PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)
CFONC
CFONC
CFONC
CFONC
CFONC
CFONC  LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
CFONC      SE FERA EN STANDARD DANS LE SOUS PROGRAMME USPHYV
CFONC      ET PAS ICI
CFONC
CFONC
CFONC
CFONC  POUR LA COMBUSTION, LE CHARBON, L'INITIALISATION EST FAITE
CFONC    PAR DEFAUT : IL EST CONSEILLE DE LA CONSERVER.
CFONC
CFONC L'identification des cellules concernees peut s'appuyer sur la
CFONC commande GETCEL.
CFONC
CFONC   GETCEL(CHAINE,NLELT,LSTELT) :
CFONC   - CHAINE est une chaine de caractere fournie par l'utilisateur
CFONC     qui donne les criteres de selection
CFONC   - NLTELT est renvoye par la commande. C'est un entier qui
CFONC     correspond au nombre de cellules trouveees repondant au
CFONC     critere
CFONC   - LSTELT est renvoye par la commande. C'est un tableau d'entiers
CFONC     de taille NLTELT donnant la liste des cellules trouvees
CFONC     repondant au critere.
CFONC
CFONC   CHAINE peut etre constitue de :
CFONC   - references de couleurs (ex. : 1, 8, 26, ...
CFONC   - references de groupes (ex. : entrees, groupe1, ...)
CFONC   - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
CFONC   Ces criteres peuvent etre combines par des operateurs logiques
CFONC   (AND et OR) et des parentheses
CFONC   ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
CFONC   les cellules de couleur 1, appartenant aux groupes 'groupe2'
CFONC   ou 'groupe3' et de coordonnee Y inferieure a 1.
CFONC
c@fonce
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
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
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
CARGU ! DT(NCELET)   ! TR ! <-> ! VALEUR DU PAS DE TEMPS               !
CARGU ! RTP          ! TR ! <-> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA COEFB  ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
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
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML), MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL, IUTILE
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) THEN
C       Indicateur de non passage dans le sous-programme
        IUSINI = 0
        RETURN
      ENDIF
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 1.  INITIALISATION VARIABLES LOCALES
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2. INITIALISATION DES INCONNUES :
C      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
C=======================================================================
C
C --- INCONNUES
C     Exemple : 1 est le numero du scalaire utilisateur 1
C                 parmi tous les scalaires (utilisateurs + physiques
C                 particulieres).
C               ISCA(1) est le numero de la variable de calcul
C                 associee.
C               RTP(IEL,ISCA(1)) est la valeur de la variable
C                 de calcul dans la cellule IEL.
C
      IF (ISUITE.EQ.0) THEN
C
        DO IEL = 1, NCEL
          RTP(IEL,ISCA(1)) = 25.D0
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 3. PAS DE TEMPS : EXEMPLE
C       SI ON FAIT UNE SUITE EN PAS DE TEMPS NON UNIFORME ET CONSTANT
C       ET QU'ON DESIRE MODIFIER LE PAS DE TEMPS STOCKE DANS LE SUITE
C         (POUR S'AFFRANCHIR D'UNE EVOLUTION TROP PROGRESSIVE PAR EX)
C=======================================================================
C
C ----------------------------------------------
C
C Il est assez courant que l'on oublie d'eliminer cet exemple
C   de ce sous-programme.
C On a donc prevu le test suivant pour eviter les mauvaises surprises
C
      IUTILE = 0
C
      IF(IUTILE.EQ.0) RETURN
C
C ----------------------------------------------
C
C
      IF(ISUITE.EQ.1.AND.(IDTVAR.EQ.1.OR.IDTVAR.EQ.2)) THEN
        DO IEL = 1, NCEL
          DT (IEL) = 10.D0*DT(IEL)
        ENDDO
      ENDIF
C
C----
C FORMATS
C----
C
C----
C FIN
C----
C
      RETURN
      END
c@z
