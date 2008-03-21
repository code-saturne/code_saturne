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
                        SUBROUTINE MTTSNS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   CRVEXP , CRVIMP ,
     &   DAM    , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
CC FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES TERMES SOURCES POUR LA COMPOSANTE DE VITESSE IVAR
CFONC
CFONC     SOUS-PROGRAMME SPECIFIQUE A MATISSE (COPIE DE USTSNS)
CFONC
CFONC
CFONC LES TERMES SOURCES SONT LES PERTES QUI N'ONT PAS ETE TRAITEES
CFONC     DANS MTKPDC (I.E. TOUT SAUF LES REGISTRES)
CFONC
CFONC IL SERAIT POSSIBLE DE REGROUPER MTTSNS ET MTKPDC EN UN SEUL
CFONC     SOUS-PROGRAMME.
CFONC
CFONC
CFONC ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP
CFONC
CFONC ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
CFONC     CRVEXP en kg variable/s :
CFONC      ex : pour la vitesse            kg m/s2
CFONC     CRVIMP en kg /s :
CFONC
CFONC VEILLER A UTILISER UN CRVIMP NEGATIF
CFONC  (ON IMPLICITERA CRVIMP
CFONC   IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
CFONC    MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
CFONC        -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2
CFONC
CFONC CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
CFONC    PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
CFONC    POUR LA RESOLUTION
CFONC
CFONC  EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
CFONC    FOURNIR CRVEXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
CFONC            CRVIMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
CFONC                                      ON LE SUPPOSE NEGATIF)
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
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
CARGU ! NCEPDP       ! E  !  -> ! NOMBRE DE CELLULES AVEC PDC          !
CARGU ! NCKPDP       ! E  !  -> ! NBR DE COEF DU TENSEUR DE PDC (3 OU 6!
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IVAR         ! E  !  -> ! NUMERO DE VARIABLE                   !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
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
CARGU ! ICEPDC(NCELET! TE !  -> ! NUMERO DES NCEPDP CELLULES AVEC PDC  !
CARGU ! ICETSM(NCESMP! TE !  -> ! NUMERO DES CELLULES A SOURCE DE MASSE!
CARGU ! ITYPSM       ! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LES     !
CARGU ! (NCESMP,NVAR)!    !     !  VARIABLES (cf. USTSMA)              !
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
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT            PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! CRVEXP(NCELET! TR ! <-  ! TABLEAU DE TRAVAIL POUR PART EXPLICIT!
CARGU ! CRVIMP(NCELET! TR ! <-  ! TABLEAU DE TRAVAIL POUR PART IMPLICIT!
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! W1...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
C-----------------------------------------------------------------------
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "matiss.h"
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
      INTEGER          NCEPDP , NCKPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IVAR   , IPHAS
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION CRVEXP(NCELET), CRVIMP(NCELET)
      DOUBLE PRECISION DAM(NCELET ),XAM(NFAC ,2)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IPCROM, IPP
      INTEGER          IEL   , IFML  , ICOUL
      INTEGER          INDIXY, IZONE , NZONES
      INTEGER          MODNTL
      DOUBLE PRECISION XLIMIN, XLIMAX, YRGMIN, YRGMAX, ZALMIN, ZALMAX
      DOUBLE PRECISION XZOMIN, XZOMAX, YZOMIN, YZOMAX
      DOUBLE PRECISION XTZC  , XLZC  , XHZC  , XLCE  , XLCS
      DOUBLE PRECISION SALV  , XLALV , DHPLEN
      DOUBLE PRECISION CXSD  , CYSD  , CTSD  , RS0   , GAPSCX
      DOUBLE PRECISION RAPSRF, SMESR , SMSSR
      DOUBLE PRECISION FPFL  , QLQP
      DOUBLE PRECISION ROE   , WE0   , VE0
      DOUBLE PRECISION ROS   , WS0   , VS0   , TS0
      DOUBLE PRECISION ROM   , W0    , WALV  , UN0   , UNH   , UREEL
      DOUBLE PRECISION REYE0 , REYS0 , REYALV, REYRES, REYHOR
      DOUBLE PRECISION RICHAR, BETMAT, GMAT
      DOUBLE PRECISION TETA0 , SDE   , SDS   , XXAA  , XXEM
      DOUBLE PRECISION ACOEFF, APRIM
      DOUBLE PRECISION CPE   , CPS   , CPC   , PDCE  , PDCS
      DOUBLE PRECISION VAMCAA, VAMHAA, PDCFAA, PDCBAA, PDCRAA, CMPTAA
C
C***********************************************************************
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
C --- Gestion memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Reperage des variables
C
      IPCROM = IPPROC(IROM  (IPHAS))
C
C --- Reperage du pas de temps pour les impressions listing
C
      IF(NTLIST.GT.0) THEN
        MODNTL = MOD(NTCABS,NTLIST)
      ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
        MODNTL = 0
      ELSE
        MODNTL = 1
      ENDIF
C
C --- Impressions
C
      IPP    = IPPRTP(IVAR)
      CHAINE = NOMVAR(IPP)
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1000) CHAINE(1:8)
      ENDIF
      IF ((IRANGP.LE.0).AND.(MODNTL.EQ.0)) THEN
        WRITE(NFECRA,2001)CHAINE(1:8)
      ENDIF
C
C --- Compteurs pour moyennes et affichages
C
C     Vitesse alveole moyenne corrigee (reelle)
      VAMCAA = 0.D0
C     Perte de charge moyenne alveole par frottement (Pa)
      PDCFAA = 0.D0
C     Perte de charge moyenne alveole par bifurcation laterale (Pa)
      PDCBAA = 0.D0
C     Perte de charge moyenne alveole par alimentation (reunification) (Pa)
      PDCRAA = 0.D0
C     Vitesse alveole moyenne homogeneisee
      VAMHAA = 0.D0
C     Compteur pour moyenne
      CMPTAA = 0.D0
C
C
C=======================================================================
C 1. DIMENSIONS GEOMETRIQUES
C=======================================================================
C
C
C --- Dimensions de la zone de stockage
C     (Transverse (x), Longitudinal (y), Hauteur (z))
C
      XTZC = PTRRES*NPTRAN
      XLZC = PLGRES*NPLGRS
      XHZC = EPCHEL*NCHEST
C
C
C --- Section debitante des cheminees (SDE : Entree Amont et SDS : Sortie)
C       FRDTRA, le facteur de reduction transverse, est le facteur
C         par lequel il faut diviser les sections d'entree ou les
C         debits si l'on represente un entrepot (reel, 3D) par un
C         calcul 2D.
C       SDCHEA et SDCHES sont les sections passantes reelles en entree
C         et sortie, au dessus des convergents eventuels.
C       SDE    et SDS    sont les sections passantes correspondantes
C         a prendre en compte dans le calcul (sections reelles divisees
C         par FRDTRA, rendant compte d'une reduction eventuelle du domaine
C         dans la direction transverse par rapport a la realite).
C         Neanmoins, SDE et SDS ne sont pas necessairement les sections
C         passantes geometriquement representees dans le maillage (voir
C         les rapports SMESR et SMSSR).
C
      SDE = SDCHEA/FRDTRA
      SDS = SDCHES/FRDTRA
C
C
C --- Longueur cheminee
C     (XLCE : cheminee d'Entree sans le convergent,
C      XLCS : cheminee de Sortie sans le convergent ou sans le toit)
C
      XLCE = HCHALI - HCONVE
      IF (ITYPEN.EQ.1) THEN
        XLCS = HCHEVA - HFTTOI
      ELSE
        XLCS = HCHEVA - HCONVE
      ENDIF
C
C
C --- Surface debitante d'une alveole (entre chemise et colis)
C
      SALV = PI * ABS(DHALVE**2 - DMCONT**2)*0.25D0
C
C
C --- Rapport des sections modele/voulu : SMESR et SMSSR
C
C       On multiplie les vitesses issues directement du calcul
C         dans les cheminees (au dessus des convergents eventuels)
C         par SMESR ou SMSRS pour obtenir des vitesses reelles.
C       SMESR est associe a la cheminee d'entree et SMSSR a la
C         cheminee de sortie.
C       SMESR et SMSSR different de FRDTRA (voir FRDTRA plus haut).
C
C       SMESR (resp. SMSSR) est le rapport entre la section passante
C         de la cheminee d'entree (resp. de sortie) telle qu'elle est
C         representee geometriquement dans le modele et la section
C         passante que l'on souhaite a modeliser (les sections sont
C         prises au dessus du convergent eventuel).
C       La section representee geometriquement peut etre differente de
C         la section que l'on souhaite modeliser par exemple si les
C         cheminees reelles sont de section circulaire alors qu'on les
C         modelise de section carree.
C       La section que l'on souhaite modeliser peut etre differente de
C         la section reelle si l'on realise un calcul 2D pour
C         representer une configuration reelle 3D (rapport FRDTRA).
C
C      Pour un cas modelise en 3D :
C         FRDTRA = 1, i.e. SDE=SDCHEA : la section que l'on souhaite
C           modeliser est egale a la section reelle,
C         on a, pour la cheminee d'entree (et de meme pour la cheminee
C           de sortie) :
C             SMESR = EPCHEM*XTZC/(RCONVE*SDCHEA)
C         ceci est bien le rapport des sections passantes
C         "(representee geometriquement)/(que l'on souhaite modeliser)"
C         (les sections sont prises au dessus du convergent eventuel).
C
C       Pour un cas modelise en 2D :
C         FRDTRA non unite, i.e. SDE=SDCHEA/FRDTRA : la section que
C           l'on souhaite modeliser est differente de la section reelle
C           (en general la section que l'on souhaite modeliser est une
C            fraction de la section reelle, dans la mesure ou la
C            puissance des colis modelises en 2D est une fraction de
C            la puissance totale reelle)
C         RCONV = 1 puisqu'il n'est pas possible de representer
C           geometriquement un convergent sur un maillage 2D,
C         on a, pour la cheminee d'entree (et de meme pour la cheminee
C           de sortie) :
C             SMESR = EPCHEM*XTZC/SDE
C         ceci est bien le rapport des sections passantes
C         "(representee geometriquement)/(que l'on souhaite modeliser)"
C         (les sections etant prises au dessus du convergent eventuel).
C
      SMESR = EPCHEM*XTZC/(RCONVE*SDE)
      SMSSR = EPCHEM*XTZC/(RCONVE*SDS)
C
C
C --- Rapport pas du reseau / diametre d'un conteneur
C
      CXSD = PTRRES/DMCONT
      CYSD = PLGRES/DMCONT
      CTSD = SQRT((CXSD*0.5D0)**2 + CYSD**2)
C
C
C --- Porosite de gap entre tubes
C
      GAPSCX = (CXSD - 1.D0) / CXSD
C
C=======================================================================
C 2. GRANDEURS PHYSIQUES
C=======================================================================
C
C --- Temperature et vitesse verticale en entree et sortie
C       (la vitesse servira pour le calcul des pertes de charges)
C
C     W moyenne dans la cheminee d'entree et volume correspondant
      WE0 = 0.D0
      VE0 = 0.D0
C     W moyenne dans la cheminee de sortie et volume correspondant
      WS0 = 0.D0
      VS0 = 0.D0
C     T moyenne dans la cheminee de sortie
      TS0 = 0.D0
C
      DO IEL = 1, NCEL
        W0 = RTPA(IEL,IW(IPHAS))
        IFML  = IFMCEL(IEL   )
        ICOUL = IPRFML(IFML,1)
C     Cheminee d'entree
        IF(ICOUL.EQ.ICMTCI) THEN
          VE0 = VE0 + VOLUME(IEL)
          WE0 = WE0 + VOLUME(IEL)*W0
C     Cheminee de sortie
        ELSEIF(ICOUL.EQ.ICMTCO) THEN
          VS0 = VS0 + VOLUME(IEL)
          WS0 = WS0 + VOLUME(IEL)*W0
          TS0 = TS0 + VOLUME(IEL)*RTPA(IEL,ISCA(ITAAMT))
        ENDIF
      ENDDO
C
      WE0 = ABS(WE0)/MAX(VE0,EPZERO)
      WS0 = ABS(WS0)/MAX(VS0,EPZERO)
      TS0 =     TS0 /MAX(VS0,EPZERO)
C
C
C --- Masse volumique entree a TINIT et sortie a TS0
C
      ROE = (TRFMAT+TKELVI)*RRFMAT/(TINIT+TKELVI)
      ROS = (TRFMAT+TKELVI)*RRFMAT/(TS0  +TKELVI)
C
C
C --- Calcul des coefficients de perte de charge des cheminees
C
C   - Entree
C
C     Nombre de Reynolds (on impose 1 au minimum)
      REYE0 = ROE*WE0*SMESR*DHCHEA/XMUMAT + 1.D0
C
C     Coefficient de perte de charge : frottement + diffuseur + filtre
      CPE = 0.3164D0 * REYE0**(-0.25D0) * (XLCE/DHCHEA)
     &     + PDCCHA + PDCFCH
C
C
C   - Sortie
C
C     Nombre de Reynolds (on impose 1 au minimum)
      REYS0 = ROS*WS0*SMSSR*DHCHES/XMUMAT + 1.D0
C
C     Coefficient de perte de charge : frottement + diffuseur + clapet
      CPS = 0.3164D0 * REYS0**(-0.25D0) * (XLCS/DHCHES)
     &     + PDCCHE + PDCCCH
C
C     Pertes de charge (en Pa) et affichage
      IF (IVAR.EQ.IW(IPHAS)) THEN
C
        PDCE = CPE * 0.5D0*ROE*(WE0*SMESR)**2
        PDCS = CPS * 0.5D0*ROS*(WS0*SMSSR)**2
C
        IF ((IRANGP.LE.0).AND.(NTCABS.EQ.NTMABS)) THEN
          WRITE(IMPMAT,5008) PDCE
          WRITE(IMPMAT,5009) PDCS
        ENDIF
C
        IF ((IRANGP.LE.0).AND.(MODNTL.EQ.0)) THEN
          WRITE(NFECRA,2002)
     &         REYE0    , REYS0    ,
     &         WE0*SMESR, WS0*SMSSR,
     &         WE0      , WS0      ,
     &         SMESR    , SMSSR    ,
     &         PDCE     , PDCS     ,
     &         TINIT    , TS0      ,
     &         ROE      , ROS      ,
     &         CPE      , CPS
        ENDIF
C
      ENDIF
C
C
C  --- Perte de charge reseau Diagramme 8.11 et 8.12 IdelCik
C     Definition de XXEM et XXAA
C
C    - Reseau de conteneurs en ligne (reseau a pas carre)
C
      IF(ICONLG.EQ.1) THEN
C
        XXEM = -0.2D0
        RS0 = (CXSD - 1.D0)/(CYSD - 1.D0)
        XXAA = NPLGRS * (CXSD - 1.D0)**(-0.5D0)
        IF(PTRRES.LE.PLGRES) THEN
          XXAA = XXAA * 1.52D0 * RS0**(-0.2D0)
        ELSE
          XXEM = XXEM * RS0**(-2.D0)
          XXAA = XXAA * 0.32D0 * (RS0-0.9D0)**(-0.68D0)
        ENDIF
C
C    - Reseau de conteneurs en quinconce (reseau a pas triangulaire)
C
      ELSE
C
        XXEM =-0.27D0
        RS0 = (CXSD - 1.D0)/(CTSD - 1.D0)
        XXAA = (4.6D0 - 2.7D0*RS0)*(2.D0 - CXSD) + 3.2D0
        XXAA = MIN(XXAA,3.2D0)
C       Pour N RANGEE de tubes
        XXAA = XXAA * (NPLGRS+1)
C
      ENDIF
C
C     Le Reynolds est calcule pour un ecoulement de circulation naturelle
C       pour un diametre de DH m
C       Il sert dans les pertes de charges plus bas
      REYRES = (ROE+ROS)*0.5D0 * VITREF * DMCONT / (GAPSCX * XMUMAT)
C
C     La vitesse sert a l'affichage
      UN0  = VITREF/(GAPSCX**2)
C
C     Affichage
      IF( IVAR.EQ.IW(IPHAS).AND.(MODNTL.EQ.0)
     &                     .AND.(IRANGP.LE.0) ) THEN
        WRITE(NFECRA,2003)
     &       0.5D0*RRFMAT*UN0*VITREF*AMPPDC*XXAA*(REYRES**XXEM),
     &       XXAA*(REYRES**XXEM)
      ENDIF
C
C
C=======================================================================
C 3. DEFINITION DES PERTES DE CHARGE
C=======================================================================
C
C === Boucle sur les cellules et selection selon leur couleur
C ===========================================================
C
      DO IEL = 1, NCEL
C
        ROM = PROPCE(IEL,IPCROM)
C
        IFML  = IFMCEL(IEL   )
        ICOUL = IPRFML(IFML,1)
C
C === Cheminee d'alimentation (au dessus du convergent eventuel)
C ===========================================================
C
        IF(ICOUL.EQ.ICMTCI) THEN
          CRVIMP(IEL) =
     &         -0.5D0*VOLUME(IEL)*ROM*WE0*SMESR*SMESR/XLCE*CPE
C
C === Cheminee d'evacuation (au dessus du convergent eventuel)
C ===========================================================
C
        ELSEIF(ICOUL.EQ.ICMTCO) THEN
          CRVIMP(IEL) =
     &         -0.5D0*VOLUME(IEL)*ROM*WS0*SMSSR*SMSSR/XLCS*CPS
C
C === Sortie de la zone de stockage
C ===========================================================
C
        ELSEIF((ICOUL.EQ.ICMTRO).OR.
     &         (ICOUL.EQ.ICMTJO.AND.IALVEO.EQ.1)) THEN
C
C     Pour representer les murs en sortie, on impose des pertes de
C       charge "infinies". Les portes de sortie doivent rester libres.
C
C     La position des murs dans la direction X est definie par des
C       zones dont les bornes XLIMIN et XLIMAX sont stockees dans
C       VIZCAR (avec l'indicateur ICPDCS faisant reference aux pertes de
C       charges a la sortie de la zone de stockage)
          NZONES = NZOCAR(ILIGNE,ICPDCS)
          DO IZONE = 1, NZONES
            XLIMIN = PTRRES*VIZCAR(1,IZONE,ILIGNE,ICPDCS)
            XLIMAX = PTRRES*VIZCAR(2,IZONE,ILIGNE,ICPDCS)
            IF ( (XYZCEN(1,IEL).GE.XLIMIN).AND.
     &           (XYZCEN(1,IEL).LT.XLIMAX) ) THEN
              CRVIMP(IEL) = CRVIMP(IEL) - 1.D5*VOLUME(IEL)
            ENDIF
          ENDDO
C
C     La position des murs dans la direction Z est definie par des
C       zones dont les bornes ZALMIN et ZALMAX sont stockees dans
C       VIZCAR (avec l'indicateur ICPDCS faisant reference aux pertes de
C       charges a la sortie de la zone de stockage)
          NZONES = NZOCAR(IALTIT,ICPDCS)
          DO IZONE = 1, NZONES
            ZALMIN = EPCHEL*VIZCAR(1,IZONE,IALTIT,ICPDCS)
            ZALMAX = EPCHEL*VIZCAR(2,IZONE,IALTIT,ICPDCS)
            IF ( (XYZCEN(3,IEL).GE.ZALMIN).AND.
     &           (XYZCEN(3,IEL).LT.ZALMAX) ) THEN
              CRVIMP(IEL) = CRVIMP(IEL) - 1.D5*VOLUME(IEL)
            ENDIF
          ENDDO
C
C
C === Reseau dans la zone de stockage
C ===========================================================
C
C     On traite successivement :
C       . pertes de charge de base en reseau
C       . pertes de charge verticales dans les alveoles
C       . pertes de charge horizontales pour le stockage en alveoles
C       . pertes de charge en zone libre de colis, stockage sans alveole
C
C     En resume :
C
C     - par defaut,   partout,              pour x, y, z : pdc=    reseau
C
C     - en alveole,   zone de colis, pour       z : pdc=    frot., bif., reu.
C     - en alveole,   plenums,              pour x, y    : pdc=pdc+frot., bif., reu.
C     - en alveole,   zone de colis, pour x, y    : pdc=    plenum inf*1.D4
C
C     - sans alveole, plenum sup,    pour x, y, z : pdc=    zero
C     - sans alveole, zone de colis, pour x, y, z : pdc=pdc*cartes (*0 ou *1)
C
C
        ELSEIF(ICOUL.EQ.ICMTST) THEN
C
C --- Vitesse de reference (corrigee : vitesse en reseau)
C
          UN0 = SQRT(RTPA(IEL,IU(IPHAS))**2
     &         + RTPA(IEL,IV(IPHAS))**2
     &         + RTPA(IEL,IW(IPHAS))**2)
          UN0 = MAX(VITREF,UN0)
          UN0 = UN0/(GAPSCX**2)
C
C
C --- Pertes de charge de base en reseau (cas standard)
C
          IF(IVAR.EQ.IW(IPHAS)) THEN
            TETA0 = 0.53D0 * XXAA * REYRES**XXEM
          ELSE
            TETA0 = 1.00D0 * XXAA * REYRES**XXEM
          ENDIF
C
C       Perte de charge en reseau de base
C       Attention, cette valeur est ecrasee ensuite
C         - pour le traitement des alveoles eventuelles
C         - dans les zones libres de colis
C
          CRVIMP(IEL) =
     &         -0.5D0*VOLUME(IEL)*ROM*TETA0*UN0*AMPPDC/XLZC
C
C
C --- Pertes de charges verticales dans les alveoles
C
C     On traite successivement :
C       . la perte de charge de frottement
C       . la perte de charge par bifurcation branche laterale
C       . la perte de charge d'alimentation apres alveole (reunification)
C
          IF(IALVEO.EQ.1) THEN
C
            IF ( XYZCEN(3,IEL).LT.HRESO.AND.
     &           XYZCEN(3,IEL).GT.HPLEN ) THEN
C
              IF(IVAR.EQ.IW(IPHAS)) THEN
C
C
C   - Perte de charge verticale par frottement
C
C     Vitesse verticale homogeneisee dans l'alveole
                WALV    = ABS (RTPA(IEL,IW(IPHAS)))
C     Rapport surface d'un pas / surface passante reelle d'une alveole
                RAPSRF = PTRRES * PLGRES / SALV
C     Reynolds en alveole base sur une vitesse corrigee (minimum a 1)
                REYALV  = 1.D0 +
     &               ROM*WALV*RAPSRF*ABS(DHALVE-DMCONT)/XMUMAT
C     Hauteur d'alveole
                XLALV   = HRESO - HPLEN
C     Coefficient de perte de charge reguliere par frottement par metre
                CPC     =
     &               0.3164D0*REYALV**(-0.25D0)/ABS(DHALVE-DMCONT)
C
C     Perte de charge verticale par frottement
C       On ecrase la perte de charge standard en reseau
C       calculee precedemment.
                CRVIMP(IEL) =
     &               - 0.5D0*VOLUME(IEL)*ROM*WALV*(RAPSRF**2)*CPC
C
C
C     Compteurs pour impressions
C
                CMPTAA = CMPTAA + 1.D0
                PDCFAA = PDCFAA
     &               + 0.5D0*ROM*(WALV*RAPSRF)**2 * CPC*XLALV
                VAMCAA = VAMCAA + WALV * RAPSRF
                VAMHAA = VAMHAA + WALV
C
C
C   - Perte de charge verticale par bifurcation branche laterale
C     (cf Idel'cik, page 263, diagramme 7.21)
C
C     Rapport de surface (FPFL) et de debit (QLQP)
C       entre la branche laterale et la branche rectiligne de
C       la bifurcation (cf Idel'cik)
                FPFL = HPLEN * PTRRES / SALV
                QLQP = 1.D0 / NPLGRS
C
                IF(QLQP*FPFL.GT.0.8D0) THEN
                  APRIM = 0.9D0
                ELSE
                  APRIM = 1.0D0
                ENDIF
C
C     Coefficient de perte de charge singuliere par bifurcation
                CPC = APRIM*(1.D0+(QLQP*FPFL)**2)/((QLQP*FPFL)**2)
C
C     Perte de charge verticale par bifurcation laterale
C       On l'ajoute a la perte de charge par frottement
C       calculee precedemment.
C       Comme la perte de charge est repartie sur la hauteur de l'alveole
C       il faut diviser par XLALV
                CRVIMP(IEL) = CRVIMP(IEL)
     &               - 0.5D0*VOLUME(IEL)*ROM*WALV*(RAPSRF**2)
     &               * CPC/XLALV
C
C
C     Compteurs pour impressions
                PDCBAA = PDCBAA + 0.5D0*ROM*(WALV*RAPSRF)**2 * CPC
C
C
C   - Perte de charge verticale par alimentation apres alveole (reunification)
C     (cf Idel'cik, page 249, diagramme 7.7)
C
C     Rapport de surface (FPFL) et de debit (QLQP)
C       entre la branche laterale et la branche rectiligne de
C       l'alimention (reunification) (cf Idel'cik)
                FPFL = (XHZC - HRESO) * PTRRES / SALV
                QLQP = 1.D0 / NPLGRS
C
                ACOEFF = 1.D0
C
C     Coefficient de perte de charge singuliere par alimentation
                CPC = ACOEFF
     &               * ( 1.D0 + (QLQP*FPFL)**2 - 2.D0*(1.D0-QLQP) )
     &               / ( (QLQP*FPFL)**2 )
C
C
C     Perte de charge verticale par alimentation (reunification)
C       On l'ajoute a la perte de charge par frottement et a la
C       perte de charge singuliere calculee precedemment.
C       Comme la perte de charge est repartie sur la hauteur de l'alveole
C       il faut diviser par XLALV
                CRVIMP(IEL) = CRVIMP(IEL)
     &               - 0.5D0*VOLUME(IEL)*ROM*WALV*(RAPSRF**2)
     &               * CPC/XLALV
C
C     Compteurs pour impressions
                PDCRAA = PDCRAA + 0.5D0*ROM*(WALV*RAPSRF)**2 * CPC
C
C
C     Fin du traitement des pertes de charges verticales en alveole
              ENDIF
            ENDIF
          ENDIF
C
C
C --- Pertes de charges horizontales pour le stockage en alveoles
C
C     On traite successivement :
C       . la perte de charge de reunion (au dessus des colis) et
C                            de bifurcation (au dessous des colis)
C       . la perte de charge de frottement
C                             au dessus et au dessous des colis
C       . la perte de charge dans la zone de colis
C
          IF(IALVEO.EQ.1) THEN
C
            IF(IVAR.EQ.IU(IPHAS).OR.IVAR.EQ.IV(IPHAS)) THEN
C
C
C   - Perte de charge horizontale de reunion et de bifurcation
C
C     La perte de charge calculee dans la zone des colis sera modifiee
C       plus bas (pour simplifier les tests, les formules valables
C       au dessus et au dessous des colis sont egalement appliquees
C       dans la zone des colis).
C
C     Rapport de debit (QLQP) entre la branche laterale et
C       la branche rectiligne (cf Idel'cik)
              QLQP = 1.D0 / NPLGRS
C
              IF(XYZCEN(3,IEL).GT.HRESO) THEN
C
C     Coefficient de perte de charge horizontale de reunion
C       pour le stockage en alveoles (apres alveole, branche rectiligne)
C       (cf Idel'cik, page 249, diagramme 7.7)
C
                CPC = (1.55D0*QLQP-QLQP**2) / (1.D0-QLQP)**2
C
              ELSE
C
C     Coefficient de perte de charge horizontale de bifurcation
C       pour le stockage en alveoles (avant alveole, branche rectiligne)
C       (cf Idel'cik, page 265, diagramme 7.23)
C
C     Ici, on traite la zone XYZCEN(3,IEL).LE.HRESO et donc non
C       seulement le dessous des colis mais aussi la zone des colis.
C       Neanmoins, la perte de charge dans la zone des colis sera
C       modifiee plus bas.
C
                CPC = 0.4D0 * QLQP**2 / (1.D0-QLQP)**2
C
              ENDIF
C
C
C     Perte de charge horizontale de reunion ou de bifurcation
C       On l'ajoute a CRVIMP
C       Comme la perte de charge est repartie sur la hauteur de l'alveole
C         il faut diviser par XLZC
C       On utilise la vitesse de gap UN0 (les puits sont supposes traverser
C         les plenums superieur et inferieur)
              CRVIMP(IEL) = CRVIMP(IEL)
     &             - 0.5D0*VOLUME(IEL)*ROM*UN0 * CPC/XLZC
C
C
C   - Perte de charge horizontale de frottement
C       au dessus et au dessous des colis
C
C     La perte de charge calculee dans la zone des colis sera modifiee
C       plus bas (pour simplifier les tests, les formules valables
C       au dessus et au dessous des colis sont egalement appliques
C       dans la zone des colis).
C
C     Vitesse homogeneisee locale
              UNH = SQRT(RTPA(IEL,IU(IPHAS))**2
     &             + RTPA(IEL,IV(IPHAS))**2
     &             + RTPA(IEL,IW(IPHAS))**2)
C
C     Reynolds au dessus et au dessous des colis (minimum 1)
C       seul le diametre hydraulique considere differe
              IF(XYZCEN(3,IEL).GT.HRESO) THEN
                DHPLEN = 2.D0*(XHZC - HRESO)
              ELSE
                DHPLEN = 2.D0*HPLEN
              ENDIF
              REYHOR = ROM*UNH*DHPLEN/XMUMAT + 1.D0
C
C     Coefficient de perte de charge reguliere par frottement par metre
              CPC =
     &             0.3164D0 * REYHOR**(-0.25D0) / DHPLEN
C
C     Perte de charge horizontale par frottement
C       On l'ajoute a la perte de charge par reunion et bifurcation
C       On ne doit pas diviser par XLZC car CPC est deja calcule par
C         metre : on n'a pas multiplie par XLZC dans CPC
              CRVIMP(IEL) = CRVIMP(IEL)
     &             - 0.5D0*VOLUME(IEL)*ROM*UN0 * CPC
C
C
C   - Perte de charge horizontale dans la zone des colis
C     On augmente la valeur de perte de charge calculee
C       precedemment, pour guider l'ecoulement sur la verticale
C
              IF ( XYZCEN(3,IEL).LT.HRESO.AND.
     &             XYZCEN(3,IEL).GT.HPLEN) THEN
                CRVIMP(IEL) = CRVIMP(IEL) * 1.D4
              ENDIF
C
C     Fin du traitement des pertes de charges horizontales en alveole
            ENDIF
          ENDIF
C
C
C --- Perte de charge en zone libre de colis en stockage sans alveole
C     Dans les zones libres, on annule la perte de charge
C
          IF(IALVEO.NE.1) THEN
C
C     Au dessus des colis
            IF (XYZCEN(3,IEL).GE.HRESO) THEN
              CRVIMP(IEL) = 0.D0
            ENDIF
C
C     Les entrepots ne sont pas uniformement remplis.
C     Pour représenter la carte d'encombrement des rangees de colis, on
C       repere les rangees occupees (coordonnee Y, indicateur IRANGE)
C       par leurs bornes YZOMIN et YZOMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPDCR faisant reference aux pertes de
C       charge en reseau)
C
            YRGMIN = 0.D0
            YRGMAX = NPLGRS*PLGRES
C
            IF ( (XYZCEN(2,IEL).GE.YRGMIN).AND.
     &           (XYZCEN(2,IEL).LT.YRGMAX) ) THEN
C
              INDIXY = 0
              NZONES = NZOCAR(IRANGE,ICPDCR)
              DO IZONE = 1, NZONES
                YZOMIN = PLGRES*VIZCAR(1,IZONE,IRANGE,ICPDCR)
                YZOMAX = PLGRES*VIZCAR(2,IZONE,IRANGE,ICPDCR)
                IF ( (XYZCEN(2,IEL).GE.YZOMIN).AND.
     &               (XYZCEN(2,IEL).LT.YZOMAX) ) THEN
                  INDIXY = 1
                ENDIF
              ENDDO
C
              IF (INDIXY.EQ.0) THEN
                CRVIMP(IEL) = 0.D0
              ENDIF
C
            ENDIF
C
C     Les entrepots ne sont pas uniformement remplis.
C     Pour représenter la carte d'encombrement des lignes de colis, on
C       repere les lignes occupees (coordonnee X, indicateur ILIGNE)
C       par leurs bornes XZOMIN et XZOMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPDCR faisant reference aux pertes de
C       charge en reseau)
C
            XLIMIN = 0.D0
            XLIMAX = NPTRAN*PTRRES
C
            IF ( (XYZCEN(1,IEL).GE.XLIMIN).AND.
     &           (XYZCEN(1,IEL).LT.XLIMAX) ) THEN
              INDIXY = 0
              NZONES = NZOCAR(ILIGNE,ICPDCR)
              DO IZONE = 1, NZONES
                XZOMIN = PTRRES*VIZCAR(1,IZONE,ILIGNE,ICPDCR)
                XZOMAX = PTRRES*VIZCAR(2,IZONE,ILIGNE,ICPDCR)
                IF ( (XYZCEN(1,IEL).GE.XZOMIN).AND.
     &               (XYZCEN(1,IEL).LT.XZOMAX) ) THEN
                  INDIXY = 1
                ENDIF
              ENDDO
C
              IF (INDIXY.EQ.0) THEN
                CRVIMP(IEL) = 0.D0
              ENDIF
C
            ENDIF
C
C     Fin du traitement des pertes de charges en zone libre sans alveole
          ENDIF
C
C
C === Entree de la zone de stockage
C ===========================================================
C
        ELSEIF((ICOUL.EQ.ICMTRI).OR.
     &         (ICOUL.EQ.ICMTJI.AND.IALVEO.EQ.1)) THEN
C
C     Pour représenter les murs en entree, on impose des pertes de
C       charge "infinies". Les portes d'entree doivent rester libres.
C
C     La position des murs dans la direction X est definie par des
C       zones dont les bornes XLIMIN et XLIMAX sont stockees dans
C       VIZCAR (avec l'indicateur ICPDCE faisant reference a la sortie)
          NZONES = NZOCAR(ILIGNE,ICPDCE)
          DO IZONE = 1, NZONES
            XLIMIN = PTRRES*VIZCAR(1,IZONE,ILIGNE,ICPDCE)
            XLIMAX = PTRRES*VIZCAR(2,IZONE,ILIGNE,ICPDCE)
            IF ( (XYZCEN(1,IEL).GE.XLIMIN).AND.
     &           (XYZCEN(1,IEL).LT.XLIMAX) ) THEN
              CRVIMP(IEL) = CRVIMP(IEL) - 1.D5*VOLUME(IEL)
            ENDIF
          ENDDO
C
C     La position des murs dans la direction Z est definie par des
C       zones dont les bornes ZALMIN et ZALMAX sont stockees dans
C       VIZCAR (avec l'indicateur ICPDCE faisant reference a la sortie)
          NZONES = NZOCAR(IALTIT,ICPDCE)
          DO IZONE = 1, NZONES
            ZALMIN = EPCHEL*VIZCAR(1,IZONE,IALTIT,ICPDCE)
            ZALMAX = EPCHEL*VIZCAR(2,IZONE,IALTIT,ICPDCE)
            IF ( (XYZCEN(3,IEL).GE.ZALMIN).AND.
     &           (XYZCEN(3,IEL).LT.ZALMAX) ) THEN
              CRVIMP(IEL) = CRVIMP(IEL) - 1.D5*VOLUME(IEL)
            ENDIF
          ENDDO
C
C     Fin du test sur les couleurs et de la boucle sur les elements
        ENDIF
      ENDDO
C
C=======================================================================
C 4. CALCULS ET AFFICHAGES FINALS
C=======================================================================
C
C --- Moyennes relatives au cas alveoles
C
      IF(CMPTAA.GT.0.D0) THEN
C
        VAMCAA = VAMCAA/CMPTAA
        PDCFAA = PDCFAA/CMPTAA
        PDCBAA = PDCBAA/CMPTAA
        PDCRAA = PDCRAA/CMPTAA
        VAMHAA = VAMHAA/CMPTAA
C
        IF ((MODNTL.EQ.0).AND.(IRANGP.LE.0)) THEN
          WRITE(NFECRA,2004) VAMCAA, VAMHAA,
     &         PDCFAA, PDCBAA, PDCRAA
        ENDIF
      ENDIF
C
C
C --- En convection naturelle,
C       . actualisation et impression de la vitesse de reference
C         (utilisee dans mtphyv en particulier)
C       . calcul et impression du Richardson
C
      IF(ICOFOR.EQ.0.AND.IVAR.EQ.IW(IPHAS))THEN
C
C   - Vitesse de reference
C
C     Calcul de la moyenne des vitesses homogenes moyennes
C       dans les cheminees d'entree et de sortie :
C       0.5D0*(WE0+WS0)
C     Cette vitesse est relative a une section de cheminee
C       sur le maillage de EPCHEM*XTZC/RCONVE (section dans
C       la partie haute, au dessus du convergent eventuel) :
C       le debit associe, represente sur le maillage, est donc
C       Q=0.5D0*(WE0+WS0)*EPCHEM*XTZC/RCONVE
C     Ce debit traverse ensuite la zone de stockage de section
C       passante XTZC*XHZC. La vitesse qui s'en deduit dans
C       cette region est donc Q/(XTZC*XHZC), soit
C       0.5D0*(WE0+WS0)*EPCHEM*XTZC/(RCONVE*XTZC*XHZC)
C     On la relaxe a partir de l'ancienne valeur de VITREF
C
        VITREF = ( 0.5D0*(WE0+WS0)*EPCHEM/(RCONVE*XHZC)
     &       + VITREF )*0.5D0
C
C     Impressions
        IF ((MODNTL.EQ.0).AND.(IRANGP.LE.0)) THEN
          WRITE(NFECRA,2005) VITREF, VISCL0(IPHAS)
        ENDIF
C
C
C   - Richardson
C
C     Calcul selon la formule Ri = g beta DeltaT H_ref / U**2 avec :
C       . beta = 1/((TS0+TINIT)*0.5D0 + TKELVI)
C       . DeltaT = TS0-TINIT (temperature sortie-entree)
C       . H_ref = XHZC (hauteur zone de stockage)
C       . U = VITREF/GAPSCX (U gap)
C
        GMAT   = SQRT(GX**2+GY**2+GZ**2)
        BETMAT = 1.D0/((TS0+TINIT)*0.5D0 + TKELVI)
        UREEL  = VITREF/GAPSCX
C
        RICHAR = GMAT*BETMAT*(TS0-TINIT)*XHZC/(UREEL**2)
C
C     Impressions
        IF ((MODNTL.EQ.0).AND.(IRANGP.LE.0)) THEN
          WRITE(NFECRA,2006) UREEL, RICHAR
        ENDIF
        IF ((IRANGP.LE.0).AND.(NTCABS.EQ.NTMABS)) THEN
          WRITE(IMPMAT,5011) RICHAR
        ENDIF
C
      ENDIF
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(' TERMES SOURCES MATISSE POUR LA VARIABLE ',A8,/)
 2001 FORMAT (/,3X,'** INFORMATIONS SUR MATISSE, VARIABLE ',A8,/,
     &          3X,'   -------------------------------------------')
 2002 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Cheminees (au dessus des convergents eventuels)         ',/,
     &'mati ---------------------------   entree      sortie   unite',/,
     &'mati --------------------------------------------------------',/,
     &'mati Nombre de Reynolds        ', E12.5   , E12.5    ,'     -',/,
     &'mati Vitesse reelle            ', E12.5   , E12.5    ,'   m/s',/,
     &'mati Vitesse homogeneisee      ', E12.5   , E12.5    ,'   m/s',/,
     &'mati Sct° maillage/a representer', E11.5  , E12.5    ,'     -',/,
     &'mati Perte de charge           ', E12.5   , E12.5    ,'    Pa',/,
     &'mati T° (TINIT, moyenne sortie)', E12.5   , E12.5    ,'    °C',/,
     &'mati Masse volumique associee  ', E12.5   , E12.5    ,' kg/m3',/,
     &'mati Coeff. de perte de charge ', E12.5   , E12.5    ,'     -',/,
     &'mati --------------------------------------------------------',/)
 2003 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Reseau                                                  ',/,
     &'mati ------                                                  ',/,
     &'mati Perte de charge reseau    (/0.53)     ', E12.5  ,'    Pa',/,
     &'mati Coeff. de perte de charge (/0.53)     ', E12.5  ,'     -',/,
     &'mati --------------------------------------------------------',/)
 2004 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Alveoles                                                ',/,
     &'mati --------                                                ',/,
     &'mati Vitesse moyenne reelle                ', E12.5  ,'   m/s',/,
     &'mati Vitesse moyenne homogeneisee          ', E12.5  ,'   m/s',/,
     &'mati Pdc reguliere  moyenne de frottement  ', E12.5  ,'    Pa',/,
     &'mati Pdc singuliere moyenne de bifurcation ', E12.5  ,'    Pa',/,
     &'mati Pdc singuliere moyenne de reunification', E11.5 ,'    Pa',/,
     &'mati --------------------------------------------------------',/)
 2005 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Vitesse de reference (zone de stockage)                 ',/,
     &'mati --------------------                                    ',/,
     &'mati Vitesse debitante homogeneisee       ', E12.5  ,'    m/s',/,
     &'mati Viscosite dynamique totale           ', E12.5  ,' kg/m/s',/,
     &'mati --------------------------------------------------------',/)
 2006 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Richardson (zone de stockage)                           ',/,
     &'mati ----------                                              ',/,
     &'mati Vitesse reelle de reference           ', E12.5  ,'   m/s',/,
     &'mati Richardson base sur la vitesse reelle ', E12.5  ,'     -',/,
     &'mati --------------------------------------------------------',/)
C
 5008 FORMAT(' Perte de charge de cheminee d''entree                  ',
     & ' :', E12.5, ' Pa')
 5009 FORMAT(' Perte de charge de cheminee de sortie                  ',
     & ':', E12.5, ' Pa')
 5011 FORMAT(' Nombre de Richardson                                   ',
     & ':', E12.5)
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
