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
                        SUBROUTINE USPROJ
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVEP   , NIVEP  , NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , STATIS , STATIV , TSLAGR , PARBOR ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC MODIFICATION UTILISATEUR EN FIN DE PAS DE TEMPS
CFONC    TOUT EST POSSIBLE,
CFONC
CFONC
CFONC ON DONNE ICI PLUSIEURS EXEMPLES :
CFONC
CFONC   - CALCUL DE BILAN THERMIQUE
CFONC     (au besoin, voir "ADAPTATION A UN SCALAIRE QUELCONQUE")
CFONC
CFONC   - CALCUL DES EFFORTS GLOBAUX SUR UN SOUS-ENSEMBLE DE FACES
CFONC
CFONC   - MODIFICATION ARBITRAIRE D'UNE VARIABLE DE CALCUL
CFONC
CFONC   - EXTRACTION D'UN PROFIL 1D
CFONC
CFONC   - IMPRESSION D'UN MOMENT
CFONC
CFONC   - EXEMPLES D'UTILISATION DES ROUTINES DE PARALLELISME
CFONC
CFONC CES EXEMPLES SONT DONNES EN SUPPOSANT UN CAS AVEC PERIODICITE
CFONC   (IPERIO.GT.0) ET PARALLELISME (IRANGP.GE.0).
CFONC
CFONC
CFONC LE CALCUL DE BILAN THERMIQUE FOURNIT EN OUTRE UNE TRAME POUR
CFONC   PLUSIEURS CHOSES
CFONC   - CALCUL DE GRADIENT (AVEC LES PRECAUTIONS UTILES EN PARALLELE ET
CFONC     PERIODIQUE)
CFONC   - CALCUL DE GRANDEUR DEPENDANT DES VALEURS AUX CELLULES VOISINES
CFONC     D'UNE FACE (AVEC LES PRECAUTIONS A PRENDRE EN PARALLELE ET
CFONC     PERIODIQUE : VOIR L'ECHANGE DE DT ET DE CP)
CFONC   - CALCUL D'UNE SOMME SUR LES PROCESSEURS LORS D'UN CALCUL
CFONC     PARALLELE (PARSOM)
CFONC
CFONC
CFONC IDENTIFICATION DES CELLULES/FACES DE BORD/FACES INTERNES
CFONC ========================================================
CFONC
CFONC   Les commandes GETCEL, GETFBR et GETFAC permettent d'identifier
CFONC   respectivement les cellules, faces ou faces de bord en fonction
CFONC   de differents criteres.
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
CFONC   La syntaxe des commandes GETFBR et GETFAC est identique.
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULES AUTORISE    !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
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
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !  (NFABOR+1)  !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR  ) !    !     !  (OPTIONNEL)                         !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
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
CARGU ! STATIS       ! TR !  -> ! MOYENNES STATISTIQUES                !
CARGU !(NCELET,NVLSTA!    !     !                                      !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,NTERSL!    !     !   LAGRANGIEN SUR LA PHASE PORTEUSE   !
CARGU ! PARBOR       ! TR !  -> ! INFOS SUR INTERACTION DES PARTICULES !
CARGU !(NFABOR,NVISBR!    !     !   AUX FACES DE BORD                  !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA), STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL    , IELG   , IFAC   , IFACG  , IVAR
      INTEGER          IEL1   , IEL2   , IELTSM
      INTEGER          IORTHO , IMPOUT
      INTEGER          IFINIA , IFINRA
      INTEGER          IGRADX , IGRADY , IGRADZ
      INTEGER          ITRAVX , ITRAVY , ITRAVZ , ITRECO
      INTEGER          INC    , ICCOCG
      INTEGER          NSWRGP , IMLIGP , IPHYDP , IWARNP
      INTEGER          IUTILE , IPHAS  , ICLVAR , III
      INTEGER          IPCROM , IPCVST , IFLMAS , IFLMAB , IPCCP, IPCVSL
      INTEGER          IDIMTE , ITENSO , ISCAL
      INTEGER          II     , NBR    , IRANGV , IRANG1 , NPOINT
      INTEGER          IMOM   , IPCMOM , IDTCM
      INTEGER          ITAB(3), IUN
      INTEGER          NCESMP , ICESMP , ISMACP , ITPSMP
      INTEGER          ILELT  , NLELT
C
      DOUBLE PRECISION XRTPA  , XRTP
      DOUBLE PRECISION XBILAN , XBILVL , XBILPA , XBILPT
      DOUBLE PRECISION XBILSY , XBILEN , XBILSO , XBILDV
      DOUBLE PRECISION XBILMI , XBILMA
      DOUBLE PRECISION EPSRGP , CLIMGP , EXTRAP
      DOUBLE PRECISION XFLUXF , XGAMMA
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ, SURFBN, DISTBR
      DOUBLE PRECISION VISCT, FLUMAB , XCP , XVSL, CP0IPH, RRR
      DOUBLE PRECISION XFOR(3), XYZ(3), XABS, XU, XV, XW, XK, XEPS
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C ---> Gestion memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C=======================================================================
C 2. EXEMPLE : CALCUL DE BILAN D'ENERGIE RELATIF A LA TEMPERATURE
C    ------------------------------------------------------------
C
C   On suppose que l'on souhaite faire des bilans (convectifs et
C     diffusifs) aux frontieres du domaine de calcul represente
C     ci-dessous (frontieres reperees par leurs couleurs).
C
C   Le scalaire considere est la temperature. On fera egalement
C     intervenir la chaleur massique (pour obtenir des bilans en Joules).
C
C
C
C   Domaine et couleurs
C   -------------------
C                  6
C      --------------------------
C      |                        |
C      |                        |
C   7  |           1            | 5
C      |     ^                  |
C      |     |                  |
C      --------------------------
C
C         2  3             4
C
C
C    2, 4, 7 : parois adiabatiques
C    6       : paroi  a temperature imposee
C    3       : entree
C    5       : sortie
C    1       : symetrie
C
C
C  ---------------------------------------------------------------------
C
C     POUR LA SIGNIFICATION PHYSIQUE DE CALCULS, IL EST UTILE D'ADOPTER
C       UN PAS DE TEMPS UNIFORME EN ESPACE (IDTVAR = 0 OU 1)
C     EN OUTRE, MEME EN SUITE DE CALCUL, LE BILAN EST FAUX SI INPDT0=1
C       (VISCT NON INITIALISE ET T(n-1) NON CONNU)
C
C  ---------------------------------------------------------------------
C
C     VARIABLE TEMPERATURE : IVAR = ISCA(ISCALT) (utiliser RTP(IEL,IVAR))
C
C  ---------------------------------------------------------------------
C
C     LE BILAN A L'INSTANT N VAUT :
C
C
C       n   iel=ncelet                n-1
C  BILAN  = SOMME { VOLUME(iel)*CP*ROM(iel)
C           iel=1
C                                      *(RTPA(iel,ivar)-RTP(iel,ivar)) }
C
C           ifac=nfabor
C         + SOMME { SURFBN(ifac)*DT(IFABOR(ifac))*CP
C           ifac=1
C
C              * [ VISLS0(ISCALT) + VISCT(IFABOR(ifac))/SIGMAS(ISCALT) ]
C              / DISTBR(ifac)
C              * [ COEFA(ifac,ICLVAR)
C                   + (COEFB(ifac,ICLVAR)-1.D0)*RTP(IFABOR(ifac,ivar)) ] }
C
C           ifac=nfabor
C         + SOMME { DT(IFABOR(ifac))*CP
C           ifac=1
C                              *RTP(IFABOR(ifac,ivar))*(-FLUMAB(ifac)) }
C
C
C  Le premier terme (nul en stationnaire) est negatif si la quantite
C    d'energie a decru dans le volume.
C  Les autres termes (convection, diffusion) sont positifs si la
C    quantite d'energie a augmente dans le volume par les
C    apports des conditions aux limites.
C
C  En regime stationnaire, un bilan positif indique donc
C   un gain d'energie.
C
C  ---------------------------------------------------------------------
C
C
C   AVEC ROM CALCULE PAR LA LOI REGISSANT LA MASSE VOLUMIQUE DANS USPHYV,
C     SOIT, PAR EXEMPLE :
C           n-1
C        ROM (iel) = P0 / [ RR * ( RTPA(iel,IVAR) + TKELV ) ]
C
C
C  ---------------------------------------------------------------------
C
C    CP ET LAMBDA/CP PEUVENT ETRE VARIABLES
C
C  ---------------------------------------------------------------------
C  ---------------------------------------------------------------------
C
C
C
C
C    ADAPTATION A UN SCALAIRE QUELCONQUE
C    -----------------------------------
C
C    L'approche peut s'utiliser pour realiser le bilan d'un scalaire
C      quelconque (mais les bilans ne sont plus attendus en Joules et
C      la chaleur massique n'intervient donc plus)
C
C    Pour cela :
C
C      - remplacer ISCALT(IPHAS) par le numero ISCAL du scalaire souhaite
C        ISCAL pouvant varier de 1 a NSCAL
C
C      - positionner IPCCP a 0 independamment de la valeur de ICP(IPHAS)
C        et affecter la valeur 1 a  CP0IPH (au lieu de CP0(IPHAS)).
C
C=======================================================================
C
C
C  Le bilan n'est pas valable si INPDT0=1
      IF (INPDT0.EQ.0) THEN
C
C
C 2.1 INITIALISATION
C ==================
C
C --> Variables locales
C     -----------------
C
C       XBILVL : bilan volumique des termes instationnaires
C       XBILDV : bilan volumique du au terme en div(rho u)
C       XBILPA : bilan en paroi adiabatique
C       XBILPT : bilan en paroi a temperature imposee
C       XBILSY : bilan en symetrie
C       XBILEN : bilan en entree
C       XBILSO : bilan en sortie
C       XBILMI : bilan lie aux injections de masse
C       XBILMA : bilan lie aux aspirations de masse
C       XBILAN : bilan total
C
      XBILVL = 0.D0
      XBILDV = 0.D0
      XBILPA = 0.D0
      XBILPT = 0.D0
      XBILSY = 0.D0
      XBILEN = 0.D0
      XBILSO = 0.D0
      XBILMI = 0.D0
      XBILMA = 0.D0
      XBILAN = 0.D0
C
C --- On travaillera sur la phase 1 uniquement
      IPHAS = 1
C
C --- Le numero du scalaire temperature est ISCAL
      ISCAL = ISCALT(IPHAS)
C
C --- Le numero de variable temperature est IVAR
      IVAR = ISCA(ISCAL)
C
C --- Le numero pour les conditions aux limites est
      ICLVAR = ICLRTP(IVAR,ICOEF)
C
C --- Numero des grandeurs physiques
      IPCROM = IPPROC(IROM(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IVAR))
      IFLMAB = IPPROB(IFLUMA(IVAR))
C
C --- On stocke dans IPCCP un indicateur permettant de determiner si
C       la Chaleur massique est constante (=CP0) ou variable. Elle sera
C       utilisee pour realiser les bilans (XBILVL est en Joules).
      IF(ICP(IPHAS).GT.0) THEN
        IPCCP  = IPPROC(ICP   (IPHAS))
      ELSE
        IPCCP  = 0
        CP0IPH = CP0(IPHAS)
      ENDIF
C
C --- On stocke dans IPCVSL un indicateur permettant de determiner si
C       la diffusivite est constante (=VISLS0) ou variable. Elle sera
C       utilisee pour les termes diffusifs.
      IF(IVISLS(ISCAL).GT.0) THEN
        IPCVSL = IPPROC(IVISLS(ISCAL))
      ELSE
        IPCVSL = 0
      ENDIF
C
C
C --> Echange de Cp et de Dt
C     ----------------------
C
C      Pour calculer les valeurs de flux aux faces internes, il est
C        necessaire d'avoir acces aux variables dans les cellules voisines.
C        En particulier, il faut connaitre la chaleur massique et la valeur
C        du pas de temps. Pour cela,
C        - dans les calculs paralleles, il est necessaire que
C          les cellules situees sur un bord de sous-domaine connaissent
C          la valeur de ces variables dans les cellules situees en
C          vis-a-vis sur le sous-domaine voisin.
C        - dans les calculs periodiques, il est necessaire que
C          les cellules periodiques aient acces a la valeur de ces0
C          variables dans les cellules periodiques correspondantes
C
C      Pour cela, il est necessaire d'appeler les routines de
C        communication PARCOM (parallelisme) et PERCOM (periodicite)
C        pour echanger les valeurs de Cp et de Dt avant de calculer le
C        gradient. L'appel de ces routines doit etre dait dans cet ordre
C        PARCOM puis PERCOM (si parallelisme et periodicite coexistent).
C
C      Si le calcul n'est ni periodique, ni parallele, on peut conserver
C        appels (les tests sur IPERIO et IRANGP assurent la generalite)
C
C
C
C    - Echange pour le parallelisme
C
        IF(IRANGP.GE.0) THEN
C
C       Echange de Dt
          CALL PARCOM (DT)
C         ===========
C
C       Echange de Cp si variable (sinon CP0(IPHAS) est utilise)
          IF(IPCCP.GT.0) THEN
            CALL PARCOM (PROPCE(1,IPCCP))
C           ===========
          ENDIF
C
        ENDIF
C
C
C    - Echange pour la periodicite
C
        IF(IPERIO.EQ.1) THEN
C
          IDIMTE = 0
          ITENSO = 0
C
C       Echange de Dt
          CALL PERCOM
C         ===========
     &      ( IDIMTE , ITENSO ,
     &        DT     , DT     , DT     ,
     &        DT     , DT     , DT     ,
     &        DT     , DT     , DT     )
C
C       Echange de Cp si variable (sinon CP0(IPHAS) est utilise)
          IF(IPCCP.GT.0) THEN
            CALL PERCOM
C           ===========
     &      ( IDIMTE , ITENSO ,
     &        PROPCE(1,IPCCP) , PROPCE(1,IPCCP) , PROPCE(1,IPCCP) ,
     &        PROPCE(1,IPCCP) , PROPCE(1,IPCCP) , PROPCE(1,IPCCP) ,
     &        PROPCE(1,IPCCP) , PROPCE(1,IPCCP) , PROPCE(1,IPCCP) )
          ENDIF
C
        ENDIF
C
C
C
C
C
C --> Calcul de la valeur reconstruite en I' pour les mailles de bord
C     ---------------------------------------------------------------
C
C     Pour les maillages orthogonaux, elle doit etre egale
C        a la valeur au centre de la cellule
C     Cette valeur est calculee dans RA(ITRECO+IFAC-1)
C                                                   (avec IFAC=1,NFABOR)
C
C     Dans le cas de maillages orthogonaux, on peut simplifier :
C       il suffit d'affecter RTP(IEL,IVAR) a RA(ITRECO+IFAC-1),
C                                                  avec IEL=IFABOR(IFAC)
C       (cette option correspond a la deuxieme branche du if ci dessous,
C        avec IORTHO different de 0)
C
C
      IORTHO = 0
C
C
C --> Cas des maillages non orthogonaux
C
      IF(IORTHO.EQ.0) THEN
C
C --- Reservation de la memoire
C
        IFINIA = IDEBIA
C
        IGRADX = IDEBRA
        IGRADY = IGRADX+NCELET
        IGRADZ = IGRADY+NCELET
        ITRAVX = IGRADZ+NCELET
        ITRAVY = ITRAVX+NCELET
        ITRAVZ = ITRAVY+NCELET
        ITRECO = ITRAVZ+NCELET
        IFINRA = ITRECO+NFABOR
C
C --- Verification de la disponibilite de la memoire
C
        CALL IASIZE('USPROJ',IFINIA)
        CALL RASIZE('USPROJ',IFINRA)
C
C
C --- Calcul du gradient de la temperature
C
C
C      Pour calculer le gradient de Temperature
C        - dans les calculs paralleles, il est necessaire que
C          les cellules situees sur un bord de sous-domaine connaissent
C          la valeur de temperature dans les cellules situees en
C          vis-a-vis sur le sous-domaine voisin.
C        - dans les calculs periodiques, il est necessaire que
C          les cellules periodiques aient acces a la valeur de la
C          temperature des cellules periodiques correspondantes
C
C      Pour cela, il est necessaire d'appeler les routines de
C        communication PARCOM (parallelisme) et PERCOM (periodicite)
C        pour echanger les valeurs de temperature avant de calculer le
C        gradient. L'appel a ces routines doit etre fait dans cet ordre
C        PARCOM puis PERCOM (pour les cas ou parallelisme et periodicite
C        coexistent).
C      En effet, on se situe ici a la fin du pas de temps n. Or,
C        les variables RTP ne seront echangees qu'en debut du pas de
C        temps n+1. Ici, seules les variables RTPA (obtenues a la fin
C        du pas de temps n-1) ont deja ete echangees.
C
C      Si le calcul n'est ni periodique, ni parallele, on peut conserver
C        appels (les tests sur IPERIO et IRANGP assurent la generalite)
C
C
C
C    - Echange pour le parallelisme
C
        IF(IRANGP.GE.0) THEN
C
          CALL PARCOM (RTP(1,IVAR))
C         ===========
C
        ENDIF
C
C    - Echange pour la periodicite
C
        IF(IPERIO.EQ.1) THEN
C
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &      ( IDIMTE , ITENSO ,
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
C
        ENDIF
C
C
C    - Calcul du gradient
C
        INC = 1
        ICCOCG = 1
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IWARNP = IWARNI(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(ITRAVX) , RA(ITRAVX) , RA(ITRAVX) ,
     &   RTP(1,IVAR) , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   RA(IGRADX) , RA(IGRADY) , RA(IGRADZ) ,
C        ----------   ----------   ----------
     &   RA(ITRAVX) , RA(ITRAVY) , RA(ITRAVZ) ,
     &   RDEVEL , RTUSER , RA     )
C
C
C    - Calcul de la valeur reconstruite dans les cellules de bord
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
          RA(ITRECO+IFAC-1) = RTP(IEL,IVAR)
     &         + DIIPBX*RA(IGRADX+IEL-1)
     &         + DIIPBY*RA(IGRADY+IEL-1)
     &         + DIIPBZ*RA(IGRADZ+IEL-1)
        ENDDO
C
C
C
C
C --> Cas des maillages orthogonaux
C
      ELSE
C
C --- Reservation de la memoire
C
        IFINIA = IDEBIA
C
        ITRECO = IDEBRA
        IFINRA = ITRECO+NFABOR
C
C --- Verification de la disponibilite de la memoire
C
        CALL IASIZE('USPROJ',IFINIA)
        CALL RASIZE('USPROJ',IFINRA)
C
C ---  Calcul de la valeur reconstruite (en fait, ici, affectation
C     de la valeur non reconstruite)
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          RA(ITRECO+IFAC-1) = RTP(IEL,IVAR)
        ENDDO
C
      ENDIF
C
C
C
C
C
C 2.2 CALCUL DU BILAN A L'INSTANT N
C =================================
C
C --> Bilan sur les volumes internes
C     ------------------------------
C
C       La masse volumique ROM est celle qui a ete calculee en debut de
C       pas de temps a partir de la temperature au pas de temps
C       precedent (dans usphyv ou, si elle est constante RO0)
C
C
      IF(IPCCP.GT.0) THEN
        DO IEL = 1, NCEL
          XRTPA = RTPA(IEL,IVAR)
          XRTP  = RTP (IEL,IVAR)
          XBILVL = XBILVL
     &        + VOLUME(IEL) * PROPCE(IEL,IPCCP) * PROPCE(IEL,IPCROM)
     &                                               * ( XRTPA - XRTP )
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          XRTPA = RTPA(IEL,IVAR)
          XRTP  = RTP (IEL,IVAR)
          XBILVL = XBILVL
     &        + VOLUME(IEL) * CP0IPH * PROPCE(IEL,IPCROM)
     &                                               * ( XRTPA - XRTP )
        ENDDO
      ENDIF
C
C
C --> Bilan sur toutes les faces, internes et de bord, pour prendre
C     en compte le terme en div(rho u)-----------------------------
C     --------------------------------
C
C     Attention, on fait intervenir les valeurs de Cp et de Dt dans
C     les cellules voisines des faces internes, ce qui necessite d'avoir
C     pris ses precautions en periodicite et parallelisme.
C
C --- Si Cp est variable
C     Noter que si Cp est variable, ecrire ici un bilan sur l'equation
C     de la temerature n'est pas absolument correct
C
      IF(IPCCP.GT.0) THEN
        DO IFAC = 1, NFAC
          IEL1 = IFACEL(1,IFAC)
          IEL2 = IFACEL(2,IFAC)
          XBILDV = XBILDV + PROPFA(IFAC,IFLMAS)
     &          *(DT(IEL1)*PROPCE(IEL1,IPCCP)*RTP(IEL1,IVAR)
     &           -DT(IEL2)*PROPCE(IEL2,IPCCP)*RTP(IEL2,IVAR))
        ENDDO
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          XBILDV = XBILDV + DT(IEL)
     &         * PROPCE(IEL,IPCCP)
     &         * PROPFB(IFAC,IFLMAB)
     &         * RTP(IEL,IVAR)
        ENDDO
C
C --- Si Cp est constant
C
      ELSE
        DO IFAC = 1, NFAC
          IEL1 = IFACEL(1,IFAC)
          IEL2 = IFACEL(2,IFAC)
          XBILDV = XBILDV + (DT(IEL1)+ DT(IEL2))*0.5D0
     &         * CP0IPH
     &         * PROPFA(IFAC,IFLMAS)
     &         * (RTP(IEL1,IVAR) - RTP(IEL2,IVAR))
        ENDDO
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          XBILDV = XBILDV + DT(IEL)
     &         * CP0IPH
     &         * PROPFB(IFAC,IFLMAB)
     &         * RTP(IEL,IVAR)
        ENDDO
      ENDIF
C
C     En cas de terme source de masse, on ajoute la contribution en Gamma*Tn+1
C
      NCESMP=NCETSM(IPHAS)
      IF (NCESMP.GT.0) THEN
        ICESMP = IICESM(IPHAS)
        ISMACP = ISMACE(IPHAS)
        ITPSMP = IITPSM(IPHAS)
        DO IELTSM = 1, NCESMP
          IEL = IA(ICESMP+IELTSM-1)
          XRTP  = RTP (IEL,IVAR)
          XGAMMA = RA( ISMACP+IELTSM+NCESMP*(IPR(IPHAS)-1)-1)
          IF(IPCCP.GT.0) THEN
            XBILDV = XBILDV
     &           - VOLUME(IEL) * PROPCE(IEL,IPCCP) * DT(IEL)
     &           * XGAMMA * XRTP
          ELSE
            XBILDV = XBILDV
     &           - VOLUME(IEL) * CP0IPH * DT(IEL)
     &           * XGAMMA * XRTP
          ENDIF
        ENDDO
      ENDIF
C
C
C
C --> Bilan sur les faces de bord
C     ---------------------------
C
C      On distingue ici les differents types de faces de bord
C        pour mieux analyser l'information, mais ce n'est pas oblige.
C
C
C     - Calcul de la contribution des parois de couleur 2, 4, 7
C         (ici adiabatiques, donc flux nul a priori)
C
      CALL GETFBR('2 or 4 or 7',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C ---   Element de bord
C
        IEL    = IFABOR(IFAC)
C
C ---   Variables geometriques
C
        SURFBN = RA(ISRFBN-1+IFAC)
        DISTBR = RA(IDISTB-1+IFAC)
C
C ---   Variables physiques
C
        VISCT  = PROPCE(IEL,IPCVST)
        FLUMAB = PROPFB(IFAC,IFLMAB)
C
        IF(IPCCP.GT.0) THEN
          XCP = PROPCE(IEL,IPCCP)
        ELSE
          XCP    = CP0IPH
        ENDIF
C
        IF(IPCVSL.GT.0) THEN
          XVSL = PROPCE(IEL,IPCVSL)
        ELSE
          XVSL = VISLS0(ISCAL)
        ENDIF
C
C ---   Calcul de la contribution au flux sur la facette courante
C         (flux de diffusion et de convection, negatif si entrant)
C
          XFLUXF =          SURFBN * DT(IEL) * XCP *
     &     (XVSL+VISCT/SIGMAS(ISCAL))/DISTBR *
     &     (COEFA(IFAC,ICLVAR)+(COEFB(IFAC,ICLVAR)-1.D0)
     &                                               *RA(ITRECO+IFAC-1))
     &                    - FLUMAB * DT(IEL) * XCP *
     &     (COEFA(IFAC,ICLVAR)+ COEFB(IFAC,ICLVAR)
     &                                               *RA(ITRECO+IFAC-1))
C
          XBILPA = XBILPA + XFLUXF
C
      ENDDO
C
C
C     - Calcul de la contribution des parois de couleur 6
C         (ici a temperature imposee ;
C                                   le flux convectif est nul a priori)
C
      CALL GETFBR('6',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C ---   Element de bord
C
        IEL    = IFABOR(IFAC)
C
C ---   Variables geometriques
C
        SURFBN = RA(ISRFBN-1+IFAC)
        DISTBR = RA(IDISTB-1+IFAC)
C
C ---   Variables physiques
C
        VISCT  = PROPCE(IEL,IPCVST)
        FLUMAB = PROPFB(IFAC,IFLMAB)
C
        IF(IPCCP.GT.0) THEN
          XCP = PROPCE(IEL,IPCCP)
        ELSE
          XCP    = CP0IPH
        ENDIF
C
        IF(IPCVSL.GT.0) THEN
          XVSL = PROPCE(IEL,IPCVSL)
        ELSE
          XVSL = VISLS0(ISCAL)
        ENDIF
C
C ---   Calcul de la contribution au flux sur la facette courante
C         (flux de diffusion et de convection, negatif si entrant)
C
          XFLUXF =          SURFBN * DT(IEL) * XCP *
     &     (XVSL+VISCT/SIGMAS(ISCAL))/DISTBR *
     &     (COEFA(IFAC,ICLVAR)+(COEFB(IFAC,ICLVAR)-1.D0)
     &                                               *RA(ITRECO+IFAC-1))
     &                    - FLUMAB * DT(IEL) * XCP *
     &     (COEFA(IFAC,ICLVAR)+ COEFB(IFAC,ICLVAR)
     &                                               *RA(ITRECO+IFAC-1))
C
          XBILPT = XBILPT + XFLUXF
C
      ENDDO
C
C
C     - Calcul de la contribution des symetries (couleur 1)
C         (a priori nul).
C
      CALL GETFBR('1',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C ---   Element de bord
C
        IEL    = IFABOR(IFAC)
C
C ---   Variables geometriques
C
        SURFBN = RA(ISRFBN-1+IFAC)
        DISTBR = RA(IDISTB-1+IFAC)
C
C ---   Variables physiques
C
        VISCT  = PROPCE(IEL,IPCVST)
        FLUMAB = PROPFB(IFAC,IFLMAB)
C
        IF(IPCCP.GT.0) THEN
          XCP = PROPCE(IEL,IPCCP)
        ELSE
          XCP    = CP0IPH
        ENDIF
C
        IF(IPCVSL.GT.0) THEN
          XVSL = PROPCE(IEL,IPCVSL)
        ELSE
          XVSL = VISLS0(ISCAL)
        ENDIF
C
C ---   Calcul de la contribution au flux sur la facette courante
C         (flux de diffusion et de convection, negatif si entrant)
C
          XFLUXF =          SURFBN * DT(IEL) * XCP *
     &     (XVSL+VISCT/SIGMAS(ISCAL))/DISTBR *
     &     (COEFA(IFAC,ICLVAR)+(COEFB(IFAC,ICLVAR)-1.D0)
     &                                               *RA(ITRECO+IFAC-1))
     &                    - FLUMAB * DT(IEL) * XCP *
     &     (COEFA(IFAC,ICLVAR)+ COEFB(IFAC,ICLVAR)
     &                                               *RA(ITRECO+IFAC-1))
C
          XBILSY = XBILSY + XFLUXF
C
      ENDDO
C
C
C     - Calcul de la contribution en entree (couleur 3)
C         (flux de diffusion et de convection)
C
      CALL GETFBR('3',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C ---   Element de bord
C
        IEL    = IFABOR(IFAC)
C
C ---   Variables geometriques
C
        SURFBN = RA(ISRFBN-1+IFAC)
        DISTBR = RA(IDISTB-1+IFAC)
C
C ---   Variables physiques
C
        VISCT  = PROPCE(IEL,IPCVST)
        FLUMAB = PROPFB(IFAC,IFLMAB)
C
        IF(IPCCP.GT.0) THEN
          XCP = PROPCE(IEL,IPCCP)
        ELSE
          XCP    = CP0IPH
        ENDIF
C
        IF(IPCVSL.GT.0) THEN
          XVSL = PROPCE(IEL,IPCVSL)
        ELSE
          XVSL = VISLS0(ISCAL)
        ENDIF
C
C ---   Calcul de la contribution au flux sur la facette courante
C         (flux de diffusion et de convection, negatif si entrant)
C
          XFLUXF =          SURFBN * DT(IEL) * XCP *
     &     (XVSL+VISCT/SIGMAS(ISCAL))/DISTBR *
     &     (COEFA(IFAC,ICLVAR)+(COEFB(IFAC,ICLVAR)-1.D0)
     &                                               *RA(ITRECO+IFAC-1))
     &                    - FLUMAB * DT(IEL) * XCP *
     &     (COEFA(IFAC,ICLVAR)+ COEFB(IFAC,ICLVAR)
     &                                               *RA(ITRECO+IFAC-1))
C
          XBILEN = XBILEN + XFLUXF
C
      ENDDO
C
C     - Calcul de la contribution en sortie (couleur 5)
C         (flux de diffusion et de convection)
C
      CALL GETFBR('5',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C ---   Element de bord
C
        IEL    = IFABOR(IFAC)
C
C ---   Variables geometriques
C
        SURFBN = RA(ISRFBN-1+IFAC)
        DISTBR = RA(IDISTB-1+IFAC)
C
C ---   Variables physiques
C
        VISCT  = PROPCE(IEL,IPCVST)
        FLUMAB = PROPFB(IFAC,IFLMAB)
C
        IF(IPCCP.GT.0) THEN
          XCP = PROPCE(IEL,IPCCP)
        ELSE
          XCP    = CP0IPH
        ENDIF
C
        IF(IPCVSL.GT.0) THEN
          XVSL = PROPCE(IEL,IPCVSL)
        ELSE
          XVSL = VISLS0(ISCAL)
        ENDIF
C
C ---   Calcul de la contribution au flux sur la facette courante
C         (flux de diffusion et de convection, negatif si entrant)
C
          XFLUXF =          SURFBN * DT(IEL) * XCP *
     &     (XVSL+VISCT/SIGMAS(ISCAL))/DISTBR *
     &     (COEFA(IFAC,ICLVAR)+(COEFB(IFAC,ICLVAR)-1.D0)
     &                                               *RA(ITRECO+IFAC-1))
     &                    - FLUMAB * DT(IEL) * XCP *
     &     (COEFA(IFAC,ICLVAR)+ COEFB(IFAC,ICLVAR)
     &                                               *RA(ITRECO+IFAC-1))
C
          XBILSO = XBILSO + XFLUXF
C
      ENDDO
C
C
C --> Bilan sur les termes sources de masse
C     -------------------------------------
C     On va séparer les injections de masse des aspirations
C     pour plus de generalite
C
      NCESMP=NCETSM(IPHAS)
      IF (NCESMP.GT.0) THEN
        ICESMP = IICESM(IPHAS)
        ISMACP = ISMACE(IPHAS)
        ITPSMP = IITPSM(IPHAS)
        DO IELTSM = 1, NCESMP
c     suivant le type d'injection, on utilise la valeur SMACEL ou la valeur ambiante
c     de la température
          IEL = IA(ICESMP+IELTSM-1)
          XGAMMA = RA( ISMACP+IELTSM+NCESMP*(IPR(IPHAS)-1)-1)
          IF ( IA( ITPSMP+IELTSM+NCESMP*(IVAR-1)-1).EQ.0
     &         .OR. XGAMMA.LT.0.D0 ) THEN
            XRTP  = RTP (IEL,IVAR)
          ELSE
            XRTP  = RA( ISMACP+IELTSM+NCESMP*(IVAR-1)-1)
          ENDIF
          IF(IPCCP.GT.0) THEN
            IF (XGAMMA.LT.0.D0) THEN
              XBILMA = XBILMA
     &             + VOLUME(IEL) * PROPCE(IEL,IPCCP) * DT(IEL)
     &             * XGAMMA * XRTP
            ELSE
              XBILMI = XBILMI
     &             + VOLUME(IEL) * PROPCE(IEL,IPCCP) * DT(IEL)
     &             * XGAMMA * XRTP
            ENDIF
          ELSE
            IF (XGAMMA.LT.0.D0) THEN
              XBILMA = XBILMA
     &             + VOLUME(IEL) * CP0IPH * DT(IEL)
     &             * XGAMMA * XRTP
            ELSE
              XBILMI = XBILMI
     &             + VOLUME(IEL) * CP0IPH * DT(IEL)
     &             * XGAMMA * XRTP
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C
C     - Somme des grandeurs sur tous les processeurs (calculs paralleles)
C
      IF (IRANGP.GE.0) THEN
        CALL PARSOM (XBILVL)
        CALL PARSOM (XBILDV)
        CALL PARSOM (XBILPA)
        CALL PARSOM (XBILPT)
        CALL PARSOM (XBILSY)
        CALL PARSOM (XBILEN)
        CALL PARSOM (XBILSO)
        CALL PARSOM (XBILMI)
        CALL PARSOM (XBILMA)
      ENDIF
C
C
C
C
C --> Bilan total
C     -----------
C
C      On ajoute les differentes contributions calculees plus haut.
C
      XBILAN = XBILVL + XBILDV + XBILPA + XBILPT + XBILSY + XBILEN
     &     + XBILSO+ XBILMI + XBILMA
C
C
C
C
C 2.3 ECRITURE DU BILAN A L'INSTANT N
C ===================================
C
      WRITE(NFECRA,2000)
     & NTCABS, XBILVL, XBILDV, XBILPA, XBILPT, XBILSY, XBILEN, XBILSO,
     &     XBILMI, XBILMA, XBILAN
C
C
 2000 FORMAT(/,
     & 3X,'** BILAN THERMIQUE **',/,
     & 3X,'   ---------------',/,
     & '---','------',
     & '------------------------------------------------------------',/,
     & 'bt ','  ITER',
     & '   Volumique  Divergence  Paroi Adia  Paroi Timp    Symetrie',
     & '      Entree      Sortie  Masse inj.  Masse asp.  Bil. total',/,
     & 'bt ',I6,10E12.4,/,
     & '---','------',
     & '------------------------------------------------------------')
C
C- Fin du test sur INPDT0
      ENDIF
C
C=======================================================================
C 3. EXEMPLE : CALCUL DES EFFORTS GLOBAUX SUR UN SOUS-ENSEMBLE DE FACES
C
C           A FAIRE AVEC PRECAUTIONS ...
C           L'UTILISATEUR PREND SES RESPONSABILITES.
C=======================================================================
C
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
C     Si les efforts ont bien ete calcules :
      IF (INEEDF.EQ.1) THEN
C
        DO II = 1, NDIM
          XFOR(II) = 0.D0
        ENDDO
C
        CALL GETFBR('2 or 3',NLELT,LSTELT)
C       ===========
C
        DO IELT = 1, NLELT
C
          IFAC = LSTELT(IELT)
C
          DO II = 1, NDIM
            XFOR(II) = XFOR(II) + RA(IFORBR + (IFAC-1)*NDIM + II-1)
          ENDDO
C
        ENDDO
C
        IF (IRANGP.GE.0) THEN
          CALL PARRSM(NDIM,XFOR)
        ENDIF
C
      ENDIF
C
C=======================================================================
C 4. EXEMPLE : MISE A 20 DE LA TEMPERATURE DANS UNE ZONE DONNEE
C              A PARTIR DU TEMPS 12s
C
C           A FAIRE AVEC PRECAUTIONS ...
C           L'UTILISATEUR PREND SES RESPONSABILITES.
C=======================================================================
C
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
      IPHAS = 1
      ISCAL = ISCALT(IPHAS)
C
      IF (TTCABS.GE.12.D0) THEN
C
        IF (ISCAL.GT.0.AND.ISCAL.LE.NSCAL) THEN
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(ISCAL)) = 20.D0
          ENDDO
        ENDIF
C
        WRITE(NFECRA,3000)
C
      ENDIF
C
 3000 FORMAT(/,
     & ' MODIFICATION UTILISATEUR DES VARIABLES EN FIN DE PAS DE TEMPS',
     & /)
C
C=======================================================================
C 5. EXEMPLE : EXTRACTION D'UN PROFIL 1D
C    -----------------------------------
C
C    On cherche ici a extraire le profil de U, V, W, k et eps sur une
C     courbe 1D quelconque en fonction de l'abscisse curviligne.
C    Le profil est ecrit dans le fichier "profil.dat" (ne pas oublier de
C     prevoir son rapatriement par le script de lancement).
C
C    - la courbe utilisee ici est le segment [(0;0;0),(0;0.1;0)], mais la
C      generalisation a une courbe quelconque est simple.
C    - la routine est adaptee au parallelisme et a la periodicite, ainsi
C      qu'au differents modeles de turbulence.
C    - la courbe 1D est discretisee en NPOINT points. Pour chacun de ces
C      points, on calcule le centre de cellule le plus proche et on
C      imprime la valeur des variables au centre de cette cellule. Pour plus
C      de coherence, l'ordonnee imprimee est celle du centre de la cellule.
C    - on evite d'imprimer deux fois la meme cellule (si plusieurs points de
C      la courbe sont associes a la meme cellule).
C=======================================================================
C
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
      IF (NTCABS.EQ.NTMABS) THEN
C
C     Seul le processeur de rang 0 (parallele) ou -1 (scalaire) ecrit dans le
C     fichier. On utilise les unites "utilisateur".
        IMPOUT = IMPUSR(1)
        IF (IRANGP.LE.0) THEN
          OPEN(IMPOUT,FILE="profil.dat")
          WRITE(IMPOUT,*)
     &         '# z(m) U(m/s) V(m/s) W(m/s) k(m2/s2) eps(m2/s3)'
        ENDIF
C
        IPHAS  = 1
        NPOINT = 200
        IEL1   = -999
        IRANG1 = -999
        DO II = 1, NPOINT
C
          XYZ(1) = 0.D0
          XYZ(2) = FLOAT(II-1)/FLOAT(NPOINT-1)*0.1D0
          XYZ(3) = 0.D0
C
          CALL FINDPT
C         ===========
     &      (NCELET, NCEL, XYZCEN,
     &       XYZ(1), XYZ(2), XYZ(3), IEL, IRANGV)
C
          IF ((IEL.NE.IEL1).OR.(IRANGV.NE.IRANG1)) THEN
            IEL1   = IEL
            IRANG1 = IRANGV
C
C     On remplit les variables temporaires XU, XV, ... pour le processeur
C      qui contient le point et on envoie ensuite l'info aux autres
C      processeurs.
            IF (IRANGP.EQ.IRANGV) THEN
              XABS = XYZCEN(2,IEL)
              XU   = RTP(IEL,IU(IPHAS))
              XV   = RTP(IEL,IV(IPHAS))
              XW   = RTP(IEL,IW(IPHAS))
              XK   = 0.D0
              XEPS = 0.D0
              IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &             .OR. ITURB(IPHAS).EQ.60) THEN
                XK = RTP(IEL,IK(IPHAS))
              ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
                XK = ( RTP(IEL,IR11(IPHAS))+RTP(IEL,IR22(IPHAS))+
     &               RTP(IEL,IR33(IPHAS)) )/2.D0
              ENDIF
              IF (ITYTUR(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.3
     &             .OR. ITURB(IPHAS).EQ.50) THEN
                XEPS = RTP(IEL,IEP(IPHAS))
              ELSEIF (ITURB(IPHAS).EQ.60) THEN
                XEPS = CMU*RTP(IEL,IK(IPHAS))*RTP(IEL,IOMG(IPHAS))
              ENDIF
            ELSE
              XABS = 0.D0
              XU   = 0.D0
              XV   = 0.D0
              XW   = 0.D0
              XK   = 0.D0
              XEPS = 0.D0
            ENDIF
C
C           Envoi aux autres processeurs si parallele
            IF (IRANGP.GE.0) THEN
              IUN = 1
              CALL PARBCR(IRANGV,IUN,XABS)
              CALL PARBCR(IRANGV,IUN,XU  )
              CALL PARBCR(IRANGV,IUN,XV  )
              CALL PARBCR(IRANGV,IUN,XW  )
              CALL PARBCR(IRANGV,IUN,XK  )
              CALL PARBCR(IRANGV,IUN,XEPS)
            ENDIF
C
            IF (IRANGP.LE.0)
     &           WRITE(IMPOUT,99) XABS,XU,XV,XW,XK,XEPS
C
 99         FORMAT(6G17.9)
C
          ENDIF
C
        ENDDO
C
        IF (IRANGP.LE.0) CLOSE(IMPOUT)
C
      ENDIF
C
C
C=======================================================================
C 6. EXEMPLE : IMPRESSION DU PREMIER MOMENT CALCULE
C=======================================================================
C
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
      IF(NBMOMT.GT.0) THEN
C
C     Numero du moment : IMOM
        IMOM = 1
C
C     Position dans PROPCE du tableau de cumul temporel des moments
C       PROPCE(IEL,IPCMOM)
        IPCMOM = IPPROC(ICMOME(IMOM))
C
C     Le cumul temporel des moments doit etre divise par la variable
C       de cumul du temps qui est un tableau NCEL ou un reel :
C             un tableau NCEL   si IDTMOM(IMOM) > 0 : PROPCE(IEL,IDTCM)
C             ou un simple reel si IDTMOM(IMOM) < 0 : DTCMOM(IDTCM)
C
        IF(IDTMOM(IMOM).GT.0) THEN
          IDTCM = IPPROC(ICDTMO(IDTMOM(IMOM)))
          DO IEL = 1, NCEL
            WRITE(NFECRA,4000) IEL,PROPCE(IEL,IPCMOM)/
     &           MAX(PROPCE(IEL,IDTCM),EPZERO)
          ENDDO
        ELSEIF(IDTMOM(IMOM).LT.0) THEN
          IDTCM = -IDTMOM(IMOM)
          DO IEL = 1, NCEL
            WRITE(NFECRA,4000) IEL,PROPCE(IEL,IPCMOM)/
     &           MAX(DTCMOM(IDTCM),EPZERO)
          ENDDO
        ENDIF
C
      ENDIF
C
 4000 FORMAT(' Cellule ',I10,'   Premier moment ',E14.5)
C
C
C=======================================================================
C 7. EXEMPLE : UTILISATION DES ROUTINES DE CALCUL PARALLELE
C              POUR LES OPERATIONS SUIVANTES~:
C=======================================================================
C
C   Cet exemple n'a pas d'autre utilite que de fournir la liste des
C     routines utilisables pour simplifier certaines operations
C     globales en parallele.
C
C   ATTENTION, ces routines modifient leur argument
C
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
C Maximum d'un compteur entier II, ici le nombre de cellules par processeur
C
C     Valeur locale
      II = NCEL
C     Calcul du maximum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARCMX(II)
      ENDIF
C     Ecriture de la valeur maximale renvoyee
      WRITE(NFECRA,5010)II
 5010 FORMAT(' USPROJ: Nombre de cellules max par processeur = ', I10)
C
C
C Somme d'un compteur entier II, ici le nombre de cellules
C
C     Valeur locale
      II = NCEL
C     Calcul de la somme sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARCPT(II)
      ENDIF
C     Ecriture de la valeur somme renvoyee
      WRITE(NFECRA,5020)II
 5020 FORMAT(' USPROJ: Nombre de cellules total = ', I10)
C
C Somme d'un reel RRR, ici le volume
C
C     Valeur locale
      RRR = 0.D0
      DO IEL = 1, NCEL
        RRR = RRR + VOLUME(IEL)
      ENDDO
C     Calcul de la somme sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARSOM(RRR)
      ENDIF
C     Ecriture de la valeur somme renvoyee
      WRITE(NFECRA,5030)RRR
 5030 FORMAT(' USPROJ: Volume du domaine total = ', E14.5)
C
C Maximum d'un reel RRR, ici le volume par processeur
C
C     Valeur locale
      RRR = 0.D0
      DO IEL = 1, NCEL
        RRR = RRR + VOLUME(IEL)
      ENDDO
C     Calcul du maximum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARMAX(RRR)
      ENDIF
C     Ecriture de la valeur maximale renvoyee
      WRITE(NFECRA,5040)RRR
 5040 FORMAT(' USPROJ: Volume max par processeur = ', E14.5)
C
C Minimum d'un reel RRR, ici le volume par processeur
C
C     Valeur locale
      RRR = 0.D0
      DO IEL = 1, NCEL
        RRR = RRR + VOLUME(IEL)
      ENDDO
C     Calcul du minimum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARMIN(RRR)
      ENDIF
C     Ecriture de la valeur minimale renvoyee
      WRITE(NFECRA,5050)RRR
 5050 FORMAT(' USPROJ: Volume min par processeur = ', E14.5)
C
C  Maximum d'un réel et valeurs reelles associées
C       ici le volume et sa localisation
C
C     Maximum local et sa localisation (NBR=3 coordonnees)
      NBR = 3
      RRR  = -1.D0
      XYZ(1) = 0.D0
      XYZ(2) = 0.D0
      XYZ(3) = 0.D0
      DO IEL = 1, NCEL
        IF (RRR.LT.VOLUME(IEL)) THEN
          RRR = VOLUME(IEL)
          XYZ(1) = XYZCEN(1,IEL)
          XYZ(2) = XYZCEN(2,IEL)
          XYZ(3) = XYZCEN(3,IEL)
        ENDIF
      ENDDO
C     Calcul du maximum et localisation sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARMXL(NBR,RRR,XYZ)
      ENDIF
C     Ecriture de la valeur maximale et localisation renvoyees
      WRITE(NFECRA,5060)RRR,XYZ(1),XYZ(2),XYZ(3)
 5060 FORMAT(' USPROJ: Volume max = ', E14.5,/,
     &       '         Localisation x,y,z = ',3E14.5)
C
C  Minimum d'un réel et valeurs reelles associées
C       ici le volume et sa localisation
C
C     Minimum local et sa localisation (NBR=3 coordonnees)
      NBR = 3
      RRR  = 1.D+30
      XYZ(1) = 0.D0
      XYZ(2) = 0.D0
      XYZ(3) = 0.D0
      DO IEL = 1, NCEL
        IF (RRR.GT.VOLUME(IEL)) THEN
          RRR = VOLUME(IEL)
          XYZ(1) = XYZCEN(1,IEL)
          XYZ(2) = XYZCEN(2,IEL)
          XYZ(3) = XYZCEN(3,IEL)
        ENDIF
      ENDDO
C     Calcul du minimum et localisation sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARMNL(NBR,RRR,XYZ)
      ENDIF
C     Ecriture de la valeur minimale et localisation renvoyees
      WRITE(NFECRA,5070)RRR,XYZ(1),XYZ(2),XYZ(3)
 5070 FORMAT(' USPROJ: Volume min = ', E14.5,/,
     &       '         Localisation x,y,z = ',3E14.5)
C
C  Somme d'un tableau d'entiers ici le nombre
C         de cellules, de faces et de faces de bord
C
C     Valeurs locales
      NBR = 3
      ITAB(1) = NCEL
      ITAB(2) = NFAC
      ITAB(3) = NFABOR
C     Calcul de la somme sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARISM(NBR,ITAB)
      ENDIF
C     Ecriture de la valeur somme renvoyee
      WRITE(NFECRA,5080)ITAB(1),ITAB(2),ITAB(3)
 5080 FORMAT(' USPROJ: Nombre de cellules       = ',I10,/,
     &       '         Nombre de faces internes = ',I10,/,
     &       '         Nombre de faces de bord  = ',I10)
C
C  Maximum d'un tableau d'entiers ici le nombre
C         de cellules, de faces et de faces de bord
C
C     Valeurs locales
      NBR = 3
      ITAB(1) = NCEL
      ITAB(2) = NFAC
      ITAB(3) = NFABOR
C     Calcul du maximum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARIMX(NBR,ITAB)
      ENDIF
C     Ecriture de la valeur maximale renvoyee
      WRITE(NFECRA,5090)ITAB(1),ITAB(2),ITAB(3)
 5090 FORMAT(' USPROJ: Nombre de cellules max       par proc = ',I10,/,
     &       '         Nombre de faces internes max par proc = ',I10,/,
     &       '         Nombre de faces de bord  max par proc = ',I10)
C
C  Minimum d'un tableau d'entiers ici le nombre
C         de cellules, et de faces de bord
C     attention, une somme similaire pour compter les faces internes
C         compte 2 fois les faces internes situes sur les bords des
C         processeurs
C
C     Valeurs locales
      NBR = 2
      ITAB(1) = NCEL
      ITAB(2) = NFABOR
C     Calcul du minimum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARIMN(NBR,ITAB)
      ENDIF
C     Ecriture de la valeur minimale renvoyee
      WRITE(NFECRA,5100)ITAB(1),ITAB(2)
 5100 FORMAT(' USPROJ: Nombre de cellules       min par proc = ',I10,/,
     &       '         Nombre de faces de bord  min par proc = ',I10)
C
C  Somme d'un tableau de reels ici les trois composantes de la vitesse
C         (dans le but de faire une moyenne ensuite par exemple)
C
C     Valeurs locales
      NBR = 3
      XYZ(1) = 0.D0
      XYZ(2) = 0.D0
      XYZ(3) = 0.D0
      DO IEL = 1, NCEL
        XYZ(1) = XYZ(1)+RTP(IEL,IU(1))
        XYZ(2) = XYZ(2)+RTP(IEL,IV(1))
        XYZ(3) = XYZ(3)+RTP(IEL,IW(1))
      ENDDO
C     Calcul de la somme sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARRSM(NBR,XYZ)
      ENDIF
C     Ecriture de la valeur somme renvoyee
      WRITE(NFECRA,5110)XYZ(1),XYZ(2),XYZ(3)
 5110 FORMAT(' USPROJ: Somme de U sur le domaine = ',E14.5,/,
     &       '         Somme de V sur le domaine = ',E14.5,/,
     &       '         Somme de W sur le domaine = ',E14.5)
C
C  Maximum d'un tableau de reels ici les trois composantes de la vitesse
C
C     Valeurs locales
      NBR = 3
      XYZ(1) = RTP(1,IU(1))
      XYZ(2) = RTP(1,IV(1))
      XYZ(3) = RTP(1,IW(1))
      DO IEL = 1, NCEL
        XYZ(1) = MAX(XYZ(1),RTP(IEL,IU(1)))
        XYZ(2) = MAX(XYZ(2),RTP(IEL,IV(1)))
        XYZ(3) = MAX(XYZ(3),RTP(IEL,IW(1)))
      ENDDO
C     Calcul du maximum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARRMX(NBR,XYZ)
      ENDIF
C     Ecriture de la valeur somme renvoyee
      WRITE(NFECRA,5120)XYZ(1),XYZ(2),XYZ(3)
 5120 FORMAT(' USPROJ: Maximum de U sur le domaine = ',E14.5,/,
     &       '         Maximum de V sur le domaine = ',E14.5,/,
     &       '         Maximum de W sur le domaine = ',E14.5)
C
C  Minimum d'un tableau de reels ici les trois composantes de la vitesse
C
C     Valeurs locales
      NBR = 3
      XYZ(1) = RTP(1,IU(1))
      XYZ(2) = RTP(1,IV(1))
      XYZ(3) = RTP(1,IW(1))
      DO IEL = 1, NCEL
        XYZ(1) = MIN(XYZ(1),RTP(IEL,IU(1)))
        XYZ(2) = MIN(XYZ(2),RTP(IEL,IV(1)))
        XYZ(3) = MIN(XYZ(3),RTP(IEL,IW(1)))
      ENDDO
C     Calcul du minimum sur les processeurs
      IF (IRANGP.GE.0) THEN
        CALL PARRMN(NBR,XYZ)
      ENDIF
C     Ecriture de la valeur somme renvoyee
      WRITE(NFECRA,5130)XYZ(1),XYZ(2),XYZ(3)
 5130 FORMAT(' USPROJ: Minimum de U sur le domaine = ',E14.5,/,
     &       '         Minimum de V sur le domaine = ',E14.5,/,
     &       '         Minimum de W sur le domaine = ',E14.5)
C
C  Envoi d'un tableau de valeurs entieres locales aux autres processeurs
C     par exemple le nombre de cellules, de faces internes et de faces
C     de bord du processeur numero IRANGV=0.
C
C     Valeurs locales
      IRANGV = 0
      NBR = 3
      ITAB(1) = NCEL
      ITAB(2) = NFAC
      ITAB(3) = NFABOR
C     Envoi aux autres
      IF (IRANGP.GE.0) THEN
        CALL PARBCI(IRANGV,NBR,ITAB)
      ENDIF
C     Ecriture par chaque processeur de la valeur envoyee par le proc 0
      WRITE(NFECRA,5140)IRANGV,ITAB(1),ITAB(2),ITAB(3)
 5140 FORMAT(' USPROJ: Sur le processeur ', I10 ,/,
     &       '         Nombre de cellules       = ',I10,/,
     &       '         Nombre de faces internes = ',I10,/,
     &       '         Nombre de faces de bord  = ',I10)
C
C  Envoi d'un tableau de valeurs reelles locales aux autres processeurs
C     par exemple trois valeurs de la vitesse
C     du processeur numero IRANGV=0.
C
C     Valeurs locales
      IRANGV = 0
      NBR = 3
      XYZ(1) = RTP(1,IU(1))
      XYZ(2) = RTP(1,IV(1))
      XYZ(3) = RTP(1,IW(1))
C     Envoi aux autres
      IF (IRANGP.GE.0) THEN
        CALL PARBCR(IRANGV,NBR,XYZ)
      ENDIF
C     Ecriture par chaque processeur de la valeur envoyee par le proc 0
      WRITE(NFECRA,5150)IRANGV,XYZ(1),XYZ(2),XYZ(3)
 5150 FORMAT(' USPROJ: Sur le processeur ', I10 ,/,
     &       '         Vitesse U dans la premiere cellule = ',E14.5,/,
     &       '         Vitesse V dans la premiere cellule = ',E14.5,/,
     &       '         Vitesse W dans la premiere cellule = ',E14.5)
C
C  Recuperation par tous les processeurs du numero global d'une cellule
C     appartenant a un processeur donne (exemple de la cellule IEL
C     du processeur IRANGV)
      IEL = 1
      IRANGV = 0
      IF (IRANGP.GE.0) THEN
        CALL PARCEL(IEL,IRANGV,IELG)
      ELSE
        IELG = -1
      ENDIF
C     Ecriture par chaque processeur du numero global de la cellule IEL
C       du processeur IRANGV
      WRITE(NFECRA,5160)IEL,IRANGV,IELG
 5160 FORMAT(' USPROJ: La cellule de numero local IEL    = ',I10,/,
     &       '            sur le processeur       IRANGV = ',I10,/,
     &       '            a pour numero global    IELG   = ',I10)
C
C  Recuperation par le processeur courant du numero global d'une cellule
C     IEL lui appartenant (exemple avec IEL=1, tous les processeurs ayant
C     au moins une cellule)
      IEL = 1
      CALL PARCLG(IEL,IRANGP,IELG)
C     Ecriture par chaque processeur du numero global de sa cellule IEL
C       (si le processeur a moins de IEL cellules, affiche 0)
C       (en sequentiel, affiche IEL)
      WRITE(NFECRA,5170)IEL,IRANGP,IELG
 5170 FORMAT(' USPROJ: La cellule de numero local IEL    = ',I10,/,
     &       '            sur le processeur       IRANGP = ',I10,/,
     &       '            a pour numero global    IELG   = ',I10)
C
C  Recuperation par le processeur courant du numero global d'une face
C     interne IFAC lui appartenant (exemple avec IFAC=1)
      IFAC = 1
      CALL PARFIG(IFAC,IRANGP,IFACG)
C     Ecriture par chaque processeur du numero global de sa face IFAC
C       (si le processeur a moins de IFAC faces internes, affiche 0)
C       (en sequentiel, affiche IFAC)
      WRITE(NFECRA,5180)IFAC,IRANGP,IFACG
 5180 FORMAT(' USPROJ: La face interne de numero local IFAC = ',I10,/,
     &       '            sur le processeur          IRANGP = ',I10,/,
     &       '            a pour numero global        IFACG = ',I10)
C
C  Recuperation par le processeur courant du numero global d'une face
C     de bord IFAC lui appartenant (exemple avec IFAC=1)
      IFAC = 1
      CALL PARFBG(IFAC,IRANGP,IFACG)
C     Ecriture par chaque processeur du numero global de sa face IFAC
C       (si le processeur a moins de IFAC faces de bord , affiche 0)
C       (en sequentiel, affiche IFAC)
      WRITE(NFECRA,5190)IFAC,IRANGP,IFACG
 5190 FORMAT(' USPROJ: La face de bord de numero local IFAC = ',I10,/,
     &       '            sur le processeur          IRANGP = ',I10,/,
     &       '            a pour numero global        IFACG = ',I10)
C
      RETURN
      END
c@z
