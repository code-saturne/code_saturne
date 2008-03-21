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
                        SUBROUTINE MTPHYV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE MATISSE : REMPLISSAGE DES VARIABLES PHYSIQUES
CFONC              COPIE DE LA ROUTINE UTILISATEUR USPHYV
CFONC
CFONC
CFONC
CFONC  ATTENTION : (on conserve pour memoire les mises en garde de usphyv)
CFONC  =========
CFONC
CFONC
CFONC  Il est INTERDIT de modifier la viscosite turbulente VISCT ici
CFONC         ========
CFONC   (une routine specifique est dediee a cela : usvist)
CFONC
CFONC
CFONC   Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
CFONC      ==================
CFONC     dans usini1 si on souhaite imposer une chaleur specifique
CFONC     CP variable pour la phase IPHAS (sinon: ecrasement memoire).
CFONC
CFONC
CFONC   Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
CFONC      ==================
CFONC      dans usini1 si on souhaite une diffusivite VISCLS variable
CFONC      pour le scalaire considere (sinon: ecrasement memoire).
CFONC
CFONC
CFONC
CFONC
CFONC  Remarques :
CFONC  ---------
CFONC
CFONC  Cette routine est appelee au debut de chaque pas de temps
CFONC
CFONC     Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
CFONC     grandeurs initialisees avant appel sont celles donnees
CFONC       - dans usini1 :
CFONC          . la masse volumique (initialisee a RO0(IPHAS))
CFONC          . la viscosite       (initialisee a VISCL0(IPHAS))
CFONC       - dans usiniv :
CFONC          . les variables de calcul  (initialisees a 0 par defaut
CFONC              ou a la valeur donnee dans usiniv)
CFONC
CFONC  On peut donner ici les lois de variation aux cellules
CFONC      - de la masse volumique                      ROM    kg/m3
CFONC          (et eventuellememt aux faces de bord     ROMB   kg/m3)
CFONC      - de la viscosite moleculaire                VISCL  kg/(m s)
CFONC      - de la chaleur specifique associee          CP     J/(kg degres)
CFONC      - des "diffusivites" associees aux scalaires VISCLS kg/(m s)
CFONC
CFONC
CFONC  On dispose des types de faces de bord au pas de temps
CFONC    precedent (sauf au premier pas de temps, ou les tableaux
CFONC    ITYPFB et ITRIFB n'ont pas ete renseignes)
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
CARGU ! NPHMX        ! E  !  -> ! NPHSMX                               !
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
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! IBROM        ! TE !  -> ! INDICATEUR DE REMPLISSAGE DE ROMB    !
CARGU !   (NPHMX   ) !    !     !                                      !
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
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! W1...8(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstnum.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "matiss.h"
C
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR), IBROM(NPHMX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET),W4(NCELET)
      DOUBLE PRECISION W5(NCELET),W6(NCELET),W7(NCELET),W8(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL   , ICOUL , IFML
      INTEGER          IPHAS , IVARTA, IVARTC
      INTEGER          IPCROM, IPCVIS, IPCVSL, IPCVSC, IPCVSP
      DOUBLE PRECISION XRTP  , UN0   , VISCT1
      DOUBLE PRECISION TRFMTK, USTMAT, DRFMAT
      DOUBLE PRECISION HRAY  , ZETA  , F12   , CXY
C
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
C --- Initialisation memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. REPERAGE DES VARIABLES
C=======================================================================
C
C --- Une seule phase
      IPHAS  = 1
C
C --- Reperage des scalaires (température air et de peau des colis)
C     ISCALT est complete en retour de l'IHM
C
      IVARTA  = ISCA(ISCALT(IPHAS))
      IVARTC  = ISCA(ITPCMT)
C
C --- Rang de la masse volumique de la phase courante IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCROM
C
      IPCROM = IPPROC(IROM(IPHAS))
C
C --- Rang de la viscosite dynamique de la phase courante IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCVIS
C
      IPCVIS = IPPROC(IVISCL(IPHAS))
C
C --- Rang des diffusivites des scalaires
C     dans PROPCE, prop. physiques au centre des elements       : IPCVS*
C     - IPCVSL : scalaire ITAAMT = T air ambiant (scalaire temperature)
C     - IPCVSC : scalaire ITPCMT = T PeauColis
C     - IPCVSP : scalaire ITPPMT = T PeauParoi
C
      IPCVSL = IPPROC(IVISLS(ISCALT(IPHAS)))
      IPCVSC = IPPROC(IVISLS(ITPCMT))
      IPCVSP = IPPROC(IVISLS(ITPPMT))
C
C
C=======================================================================
C 2. GRANDEURS GEOMETRIQUES
C=======================================================================
C
C --- CXY et F12 (facteur de forme rayonnement)
      CXY  = (PTRRES + PLGRES) * 0.5D0
      ZETA = CXY/DMCONT
      F12  = 2.D0/PI*( (ZETA**2 - 1.D0)**0.5D0
     &     - ZETA + ASIN(1.D0/ZETA) )
C
C
C=======================================================================
C 3. MASSE VOLUMIQUE
C=======================================================================
C
C --- Masse volumique au centre des cellules
C
C     Utilisation de la loi des gaz parfait
C     (attention, il faut des temperatures en Kelvin)
      TRFMTK = TRFMAT+TKELVI
      DO IEL = 1, NCEL
        XRTP = RTP(IEL,IVARTA)
        PROPCE(IEL,IPCROM) = TRFMTK*RRFMAT /(XRTP + TKELVI)
      ENDDO
C
C=======================================================================
C 4. VISCOSITE DYNAMIQUE ET DIFFUSIVITE DES SCALAIRES
C=======================================================================
C
C --- Remarque : par securite, on s'assure que la viscosite et la
C       diffusivite ne sont jamais nulles en ajoutant si besoins la
C       viscosite moleculaire de l'air.
C
C
      DO IEL = 1, NCEL
C
C
C --- Vitesse de reference locale : max de la vitesse instantanee et
C       locale et de la vitesse de reference calculee dans mtini1
C
        UN0 = SQRT(RTPA(IEL,IU(IPHAS))**2 + RTPA(IEL,IV(IPHAS))**2
     &       + RTPA(IEL,IW(IPHAS))**2)
        UN0 = MAX(VITREF,UN0)
C
C
C --- Viscosite totale de reference : VISCL0 calculee dans mtini1
C       MAIS en convection naturelle, VITREF est modifie dans mttsns
C         il faut donc alors recalculer VISCL0 ici.
C       Par securite on le recalcule aussi en convection forcee.
C
C     La viscosite dynamique turbulente est modelisee sous la forme
C       mu_t = rho * u*_ref * D_ref, formule dans laquelle :
C     . u*_ref est la vitesse de frottement deduite de la vitesse de
C              reference VITREF calculee precedemment et d'une intensite
C              turbulente imposee RTURB0
C     . D_ref  est une distance de reference basee sur l'encombrement
C              du milieu
C
C     La valeur de u*_ref est calculee a partir de l'energie cinetique
C       turbulente k par u*_ref = Cmu**(1/4) k**(1/2). La valeur de k
C       se deduit de l'intensite turbulente I, supposee connue, par
C       k = 3/2 (I V_ref)**2 avec V_ref la vitesse de reference
C       VITREF calculee precedemment.
C
C     La valeur de D_ref est calculee par D_ref = 0.2(Pt - d) ou Pt est
C       le pas transverse du reseau et d le diametre des conteneurs.
C
C     Ainsi, mu_t = rho * u*_ref * L_ref
C            . u*_ref = Cmu**(1/4) * (3/2)**(1/2) * I * V_ref
C            . L_ref  = 0.2 * (Pt - d)
C
C     On ajoute la viscosite moleculaire XMUMAT a mu_t pour obtenir
C       la viscosite dynamique totale VISCL0
C
        USTMAT = CMU**0.25D0 * SQRT(1.5D0) * (RTURB0/100.D0) * VITREF
        DRFMAT = 0.2D0 * (PTRRES-DMCONT)
C
        VISCL0(IPHAS) = RO0(IPHAS) * USTMAT * DRFMAT  + XMUMAT
C
C
C --- Viscosite totale VISCT1 calculee comme ci-dessus mais
C     . avec la vitesse de reference UN0 calculee ci-dessus
C     . avec la masse volumique instantanee et locale
C
C     La viscosite dynamique turbulente est modelisee sous la forme
C       mu_t = rho * u*_ref * D_ref, formule dans laquelle :
C     . u*_ref est la vitesse de frottement deduite de la vitesse de
C              reference UN0    calculee ci-dessus et d'une intensite
C              turbulente imposee RTURB0
C     . D_ref  est un distance de reference basee sur l'encombrement
C              du milieu
C
C     La valeur de u*_ref est calculee à partir de l'energie cinetique
C       turbulente k par u*_ref = Cmu**(1/4) k**(1/2). La valeur de k
C       se deduit de l'intensite turbulente I, supposee connue, par
C       k = 3/2 (I V_ref)**2 avec V_ref la vitesse de reference
C       UN0    calculee precedemment.
C
C     La valeur de D_ref est calculee par D_ref = 0.2(Pt - d) ou Pt est
C       le pas transverse du reseau et d le diametre des conteneurs.
C
C     Ainsi, mu_t = rho * u*_ref * L_ref
C            . u*_ref = Cmu**(1/4) * (3/2)**(1/2) * I * V_ref
C            . L_ref  = 0.2 * (Pt - d)
C
C     On ajoute la viscosite moleculaire XMUMAT a mu_t pour obtenir
C       la viscosite dynamique totale VISCL0
C
        USTMAT = CMU**0.25D0 * SQRT(1.5D0) * (RTURB0/100.D0) * UN0
        DRFMAT = 0.2D0 * (PTRRES-DMCONT)
C
        VISCT1 = PROPCE(IEL,IPCROM) * USTMAT * DRFMAT  + XMUMAT
C
C
C --- Reperage des elements selon les zones du maillage
C
C     Couleur de l'element courant
C     (hypotheses : une couleur et une seule, pas de groupe,
C                   et donc pas de boucle sur les familles)
        IFML  = IFMCEL(IEL   )
        ICOUL = IPRFML(IFML,1)
C
C
C --- Zone de stockage
        IF(ICOUL.EQ.ICMTST) THEN
C
C     Viscosite totale
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),VISCT1)
C
C     Conductivite totale pour l'air (Nb de Prandlt total = 1.)
          PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVIS)
C
C     Diffusivite equivalente pour la temperature de peau des colis
C     . Au dessus de HRESO ou entreposage en alveole : zero
C     . En dessous de HRESO sans alveole : rayonnement
C       phi lineique = HRAY*RF*(CXY-DMCONT)*gradT
C       avec RF=DMCONT/CXY rapport surfacique
C     On ajoute la viscosite moleculaire
C       (pour eviter zero dans les conditions aux limites)
          IF(XYZCEN(3,IEL).GE.HRESO.OR.IALVEO.EQ.1) THEN
            PROPCE(IEL,IPCVSC) = 0.D0
          ELSE
            HRAY = STEPHN*F12*EMICON*(RTP(IEL,IVARTC)+TKELVI)**3
            PROPCE(IEL,IPCVSC) = HRAY*(CXY-DMCONT)*DMCONT/CXY
          ENDIF
          PROPCE(IEL,IPCVSC) = PROPCE(IEL,IPCVSC) + XMUMAT
C
C     Diffusivite equivalente pour la temperature de peau des murs
C       Pas de diffusion pour ce scalaire : traitement par terme source
C       (viscosite moleculaire pour eviter zero)
          PROPCE(IEL,IPCVSP) = XMUMAT
C
C
C --- Obstacles (registres amont)
C
        ELSEIF(ICOUL.EQ.ICMTRI)THEN
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),VISCT1)
          PROPCE(IEL,IPCVSL) = XMUMAT
          PROPCE(IEL,IPCVSC) = XMUMAT
          PROPCE(IEL,IPCVSP) = XMUMAT
C
C --- Obstacles (registres aval)
C
        ELSEIF(ICOUL.EQ.ICMTRO) THEN
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),VISCT1)
          PROPCE(IEL,IPCVSL) = XMUMAT
          PROPCE(IEL,IPCVSC) = XMUMAT
          PROPCE(IEL,IPCVSP) = XMUMAT
C
C --- Defaut (les autres zones)
C
        ELSEIF(ICOUL.EQ.ICMTDF) THEN
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),VISCT1)
          PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVIS)
          PROPCE(IEL,IPCVSC) = XMUMAT
          PROPCE(IEL,IPCVSP) = XMUMAT
C
C --- Cheminee d'alimentation (au dessus du convergent eventuel)
C     . on dope la diffusion pour homogeneiser T et U
C       et diminuer les risques de refoulement
C
        ELSEIF(ICOUL.EQ.ICMTCI) THEN
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),20.D0*VISCT1)
          PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVIS)
          PROPCE(IEL,IPCVSC) = XMUMAT
          PROPCE(IEL,IPCVSP) = XMUMAT
C
C --- Cheminee d'evacuation (au dessus du convergent eventuel)
C     . on dope la diffusion pour homogeneiser T et U
C       et diminuer les risques de refoulement
C
        ELSEIF(ICOUL.EQ.ICMTCO) THEN
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),20.D0*VISCT1)
          PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVIS)
          PROPCE(IEL,IPCVSC) = XMUMAT
          PROPCE(IEL,IPCVSP) = XMUMAT
C
C --- Sinon ...
C
        ELSE
          PROPCE(IEL,IPCVIS) = MAX(VISCL0(IPHAS),VISCT1)
          PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVIS)
          PROPCE(IEL,IPCVSC) = XMUMAT
          PROPCE(IEL,IPCVSP) = XMUMAT
        ENDIF
C
      ENDDO
C
C
C=======================================================================
C FORMATS
C=======================================================================
C
C----
C FIN
C----
C
      RETURN
      END
c@z
