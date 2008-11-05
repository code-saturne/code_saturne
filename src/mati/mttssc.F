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
                        SUBROUTINE MTTSSC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   ICONRA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   CRVEXP , CRVIMP ,
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    , W11    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC
CFONC  CALCUL DES TERMES SOURCES POUR LE SCALAIRE ISCAL
CFONC
CFONC     SOUS-PROGRAMME SPECIFIQUE A MATISSE (COPIE DE USTSSC)
CFONC
CFONC     LES SEULS SCALAIRES TRAITES SONT
CFONC       ISCA(ITAAMT) TEMPERATURE DE L'AIR
CFONC       ISCA(ITPCMT) TEMPERATURE DE PEAU DES COLIS
CFONC       ISCA(ITPPMT) TEMPERATURE DE PEAU DES MURS
CFONC
CFONC
CFONC ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP
CFONC
CFONC ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
CFONC     CRVEXP en kg variable/s :
CFONC      ex : pour les temperatures      kg degres/s
CFONC           pour les enthalpies        Joules/s
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
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ISCAL        ! E  !  -> ! NUMERO DU SCALAIRE                   !
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
CARGU ! ICONRA       ! TE ! <-> ! TAB DE CONNECTIVITE POUR             !
CARGU ! (NCELET+1)   !    !     !  LE RAYONNEMENT ET LES PANACHES      !
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
CARGU ! CKUPDC       ! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !  (NCEPDP,6)  !    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! CRVEXP(NCELET! TR ! <-  ! TABLEAU DE TRAVAIL POUR PART EXPLICIT!
CARGU ! CRVIMP(NCELET! TR ! <-  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! XAM(NFAC,2)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..11(NCELET! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
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
      INCLUDE "cstnum.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
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
      INTEGER          NCEPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ISCAL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
      INTEGER          ICONRA(NCELET)
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
      DOUBLE PRECISION CKUPDC(NCEPDP,6), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION CRVEXP(NCELET), CRVIMP(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION XAM(NFAC,2)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR  , IISCVR, IPCROM, IPHAS
      INTEGER          ICOUL , IFML  , IEL   , JEL
      INTEGER          MODNTL
      INTEGER          NZONES, IZONE
      DOUBLE PRECISION XNU0  , XLAMB0, PRAND0, PHI0
      DOUBLE PRECISION GMAT  , UNTIER
      DOUBLE PRECISION FFORM , EMIS
      DOUBLE PRECISION HPLUS , HMOIN , HMAX  , DIFX  , DIFY  , DHMUR
      DOUBLE PRECISION FRAGH , XXAV  , TBULK , UN0   , PTOT
      DOUBLE PRECISION RALEIZ, REYNOL
      DOUBLE PRECISION XNSRAZ, XNSREY, XNSALV
      DOUBLE PRECISION HRAZ  , HREY  , H0    , HRAY
      DOUBLE PRECISION HMREY , HMRAZ , HM0   , HMRAY , HNB
      DOUBLE PRECISION CPCONT, CPPARO
      DOUBLE PRECISION XLIMIN, XLIMAX, YRGMIN, YRGMAX, ZALMIN, ZALMAX
      DOUBLE PRECISION SQZ
      DOUBLE PRECISION DHE   , DHS   , PTCN
      DOUBLE PRECISION FPTSTO, PVCCSC, PCCSC
C
C***********************************************************************
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C --- Gestion memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C --- Numero du scalaire a traiter : ISCAL (argument)
C
C --- Numero de la variable associee au scalaire a traiter ISCAL
      IVAR = ISCA(ISCAL)
C
C --- Nom de la variable associee au scalaire a traiter ISCAL
      CHAINE = NOMVAR(IPPRTP(IVAR))
C
C --- Indicateur de variance
C         Si ISCAVR = 0 :
C           le scalaire ISCAL n'est pas une variance
C         Si ISCAVR > 0 et ISCAVR < NSCAL + 1 :
C           le scalaire ISCAL est une variance associee
C           au scalaire ISCAVR
      IISCVR = ISCAVR(ISCAL)
C
C --- Numero de phase associee au scalaire ISCAL
      IPHAS = IPHSCA(ISCAL)
C
C --- Masse volumique
      IPCROM = IPPROC(IROM(IPHAS))
C
C
C --- Reperage du pas de temps pour les impressions listing
      IF(NTLIST.GT.0) THEN
        MODNTL = MOD(NTCABS,NTLIST)
      ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
        MODNTL = 0
      ELSE
        MODNTL = 1
      ENDIF
C
C --- Impression
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1000) CHAINE(1:8)
      ENDIF
      IF((IRANGP.LE.0).AND.(MODNTL.EQ.0)) THEN
        WRITE(NFECRA,2001) CHAINE(1:8)
      ENDIF
C
C --- Numerique
      UNTIER = 1.D0/3.D0
C
C=======================================================================
C 1. DIMENSIONS GEOMETRIQUES
C=======================================================================
C
C
C=======================================================================
C 2. PROPRIETES PHYSIQUES
C=======================================================================
C
C --- Proprietes de l'air
C
C     Nombre de Prandtl moleculaire
      PRAND0 = 0.7D0
C     Viscosite cinematique moleculaire (nu en m2.s-1)
      XNU0   = XMUMAT/RO0(IPHAS)
C     Condutivite thermique moleculaire (lambda en W.m-1.K-1)
      XLAMB0 = (XNU0/PRAND0)*CP0(IPHAS)*RO0(IPHAS)
C
C
C --- Capacite thermique fictive pour le solide
C
C     On cherche un etat stationnaire comme limite d'un transitoire
C       Les equations portant sur la temperature T des solides
C         comportent un terme en CP dT/dt qui tend vers zero a
C         convergence. La valeur de CP est utilisee pour accelerer
C         la convergence : elle correspond a l'utilisation d'un
C         pas de temps adapte pour la temperature des solides (comme
C         on pourrait le faire avec CDTVAR)
C
      CPCONT = 1.D0
      CPPARO = 1.D0
C
C
C=======================================================================
C 2. EROSION DES PANACHES DE CONVECTION NATURELLE
C=======================================================================
C
C     Le but est de reporter dans la zone de ciel une partie de la
C       puissance transmise a l'air dans la zone de stockage.
C       Ceci permet de simuler l'effet des panaches, le debit
C       enthalpique de ceux-ci etant prescrit par l'utilisateur.
C
C     Calcul de DHS, DHE
C
C --- Debit enthalpique de convection naturelle calcule (initialisation)
      DHS = 0.D0
C
C --- Debit enthalpique de convection naturelle impose pour les panaches
C
C     Nul par defaut (pas de modelisation des panaches)
      DHE = 0.D0
C     Calcul, si modelisation des panaches
      IF(IMDCNT.EQ.1) THEN
        DHE = DHPCNT/FRDTRA
      ENDIF
C
C=======================================================================
C 2. RAYONNEMENT
C=======================================================================
C
C --- Flux d'un conteneur (puissance ramenee a la surface en W.m-2)
C
C       En cas de stockage en alveole, il y a un plenum inferieur libre
      IF(IALVEO.EQ.1) THEN
        PHI0 = PUICON/((HRESO-HPLEN)*DMCONT*PI)
      ELSE
        PHI0 = PUICON/(HRESO*DMCONT*PI)
      ENDIF
C
C --- Facteur de forme (alveole ou mur) / conteneur
C
C     En alveole, le conteneur voit l'alveole entiere et c'est tout
      IF(IALVEO.EQ.1) THEN
C
        FFORM = 1.D0
C
C     Sans alveole, le facteur de forme depend de la structure du reseau
      ELSE
C
C       Pour un reseau en ligne, considerer une rangee suffit
C         (valable jusqu'a PTRRES/DMCONT=4)
        FFORM =  SQRT( SQRT(DMCONT/PTRRES) + (DMCONT/PTRRES) )
C
C       Pour un reseau en quinconce (reseau a pas triangulaire), on
C         complete par la contribution d'une seconde rangee
        IF(ICONLG.EQ.0) THEN
          FFORM = MAX( (FFORM + 0.22D0), 1.D0 )
        ENDIF
C
      ENDIF
C
C
C  --- Emissivite equivalente (conteneur / mur)
      EMIS = 1.D0/(1.D0/EMICON + 1.D0/EMIMUR)
C
C  --- Densite de surface d'echange (par unite de hauteur)
      XXAV = PI*DMCONT/(PTRRES*PLGRES)
C
C
C  --- Connectivite pour le rayonnement colis / plafond
C                et pour la modelisation des panaches
C        Une connectivite est necessaire pour relier la couche
C          superieure de la zone occupee par des colis et les
C          mailles situees directement sous le plafond
C
C     Cote du haut des mailles du haut des conteneurs
      HPLUS = HRESO
C     Cote du bas des mailles du haut des conteneurs
      HMOIN = HRESO - EPCHEL
C     Cote du bas des mailles touchant le plafond
      HMAX  = EPCHEL*NCHEST - EPCHEL
C
C     Calcul de la connectivite
C       Attention, c'est quadratique ...
C       mais on ne le fait qu'au premier passage
C
C     Calcul uniquement si pas deja fait et si utile
C       (utile si les colis ne touchent pas le plafond)
      IF( (ICNROK.EQ.0).AND.(HRESO.LT.EPCHEL*NCHEST) ) THEN
C     On repere les mailles de plafond dans la zone de stockage
        DO IEL = 1, NCEL
          IFML  = IFMCEL(IEL   )
          ICOUL = IPRFML(IFML,1)
          IF(ICOUL.EQ.ICMTST.AND.XYZCEN(3,IEL).GT.HMAX)THEN
C     On repere la maille du haut des colis situee dessous
            DO JEL = 1, NCEL
              DIFX = ABS(XYZCEN(1,JEL)-XYZCEN(1,IEL))
              DIFY = ABS(XYZCEN(2,JEL)-XYZCEN(2,IEL))
              IF(DIFX.LT.PTRRES*1.D-2.AND.DIFY.LT.PLGRES*1.D-2) THEN
                IF ( XYZCEN(3,JEL).LT.HPLUS.AND.
     &               XYZCEN(3,JEL).GT.HMOIN ) THEN
                  ICONRA(IEL) = JEL
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C     On indique que la connectivite a ete calculee
        ICNROK = 1
      ENDIF
C
C=======================================================================
C
C  ON VA COMPLETER CI-DESSOUS LES TERMES SOURCES
C     - PUREMENT EXPLICITES (PUISSANCE)
C     - PARTIELLEMENT EXPLICITES ET IMPLICITES (ECHANGES)
C
C
C  ON UTILISE LE MODELE DE TERME SOURCE SUIVANT, POUR UN SCALAIRE F :
C
C                             S = A * F + B
C
C     APPARAISSANT DANS LES EQUATIONS DE TRANSPORT DE F SOUS LA FORME :
C
C                       RHO VOLUME D(F)/Dt = VOLUME*S
C
C
C   CE TERME A UNE PARTIE QU'ON VEUT IMPLICITER         : A
C           ET UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B
C
C
C   PAR EXEMPLE SI ON A :
C     A = - RHO / TAUF
C     B =   RHO * PRODF
C        AVEC
C     TAUF   = 10.D0  [secondes  ] (TEMPS DE DISSIPATION DE F)
C     PRODF  = 100.D0 [variable/s] (PRODUCTION DE F PAR UNITE DE TEMPS)
C
C   ON A ALORS
C     CRVIMP(IEL) = VOLUME(IEL)* A = - VOLUME(IEL) (RHO / TAUF )
C     CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL) (RHO * PRODF)
C
C=======================================================================
C
C === Initialisation
C ===========================================================
C
C --- Compteurs
C       pour moyenne et affichage, mais pas seulement
C
C     Coefficient d'echange moyen de convection forcee
      HMREY = 0.D0
C     Coefficient d'echange moyen de convection naturelle
      HMRAZ = 0.D0
C     Coefficient d'echange convectif moyen efficace
      HM0   = 0.D0
C     Coefficient d'echange radiatif moyen pour les parois
      HMRAY = 0.D0
C     Compteur pour moyenne
      HNB   = 0.D0
C     Puissance totale sans panache
      PTOT  = 0.D0
C     Puissance totale avec panache
      PTCN  = 0.D0
C
C     Gravite
      GMAT = SQRT(GX**2+GY**2+GZ**2)
C
C
C === Boucle sur les cellules et selection selon couleur
C ===========================================================
C
      DO IEL = 1, NCEL
C
C     Vitesse horizontale
        UN0 = SQRT(RTP(IEL,IU(IPHAS))**2+RTP(IEL,IV(IPHAS))**2)
C
C     Temperature de reference Kelvin
        TBULK = (TINIT+RTP(IEL,ISCA(ITAAMT)))*0.5D0 + TKELVI
C
C     Raleigh/DeltaT = FRAGH = g*H**3*beta/(nu a)
C       . g    = GMAT
C       . H    = XYZCEN(3,IEL) : on suppose que le sol est a z=0
C       . beta = 1/TBULK
C       . nu   = XNU0
C       . a    = nu/Pr = XNU0/PRAND0
C
        FRAGH = GMAT*XYZCEN(3,IEL)**3*PRAND0/(TBULK*XNU0**2)
C
C     Couleur pour selection de la zone de stockage
        IFML  = IFMCEL(IEL   )
        ICOUL = IPRFML(IFML,1)
C
        IF(ICOUL.EQ.ICMTST) THEN
C
C
C === Calcul du coefficient d'echange par rayonnement
C ===========================================================
C
C     Coeff d'echange = flux / (S (Tcolis-Tparoi))
C       = Sigma Epsilon FactForme (Tcolis+Tparoi)(Tcolis**2+Tparoi**2)
C
C     La paroi consideree est le plafond ou les murs ou les alveoles
C
C --- Zone de ciel
C       On utilise la connectivite calculee precedemment
          IF(ICNROK.EQ.1.AND.XYZCEN(3,IEL).GT.HMAX) THEN
C
C     En presence d'alveoles,
C       le plafond echange avec les alveoles face a face
C       (seule l'alveole directement en face est vue,
C        l'alveole voit un disque de plafond)
C       La temperature du plafond et des alveoles est representee
C       par le meme scalaire
C       Ce coefficient est utilise pour le calcul des termes sources
C       permettant de determiner la temperature de peau des parois
C       (plafond, mur, alveoles)
C       On le verra aussi intervenir pour un calcul de temperature
C       de colis, mais dans la zone de ciel ou la temperature de colis
C       n'a pas de signification (et comme il n'y a pas de diffusion
C       en alveoles, c'est ok)
            IF(IALVEO.EQ.1) THEN
              HRAY = STEPHN*EMIMUR
     &             * (0.25D0*PI*DMCONT**2)/(PTRRES*PLGRES)
     &             * ( RTP(ICONRA(IEL),ISCA(ITPPMT))+TKELVI
     &                +RTP(       IEL ,ISCA(ITPPMT))+TKELVI)
     &             * ((RTP(ICONRA(IEL),ISCA(ITPPMT))+TKELVI)**2
     &               +(RTP(       IEL ,ISCA(ITPPMT))+TKELVI)**2)
C
C     Sans alveoles,
C       le plafond echange avec les colis face a face
C       Ce coefficient sera utilise pour le calcul des termes sources
C       permettant de determiner la temperature de peau du plafond
            ELSE
              HRAY = STEPHN*EMIS
     &             * (0.25D0*PI*DMCONT**2)/(PTRRES*PLGRES)
     &             * ( RTP(ICONRA(IEL),ISCA(ITPCMT))+TKELVI
     &                +RTP(       IEL ,ISCA(ITPPMT))+TKELVI)
     &             * ((RTP(ICONRA(IEL),ISCA(ITPCMT))+TKELVI)**2
     &               +(RTP(       IEL ,ISCA(ITPPMT))+TKELVI)**2)
            ENDIF
C
C --- Zones de murs ou d'alveoles
C       les murs et les alveoles echangent avec les colis
C       le facteur de forme calcule precedemment s'applique
C       il prend en compte le cas avec ou sans alveoles
C       Ce coefficient sera utilise pour le calcul des termes sources
C       permettant de determiner
C         - la temperature de peau des parois
C         - la temperature de peau des colis, en alveole (FFORM=1)
C
          ELSE
            HRAY = STEPHN*EMIS
     &           * FFORM
     &           * ( RTP(IEL,ISCA(ITPCMT))+TKELVI
     &              +RTP(IEL,ISCA(ITPPMT))+TKELVI)
     &           * ((RTP(IEL,ISCA(ITPCMT))+TKELVI)**2
     &             +(RTP(IEL,ISCA(ITPPMT))+TKELVI)**2)
          ENDIF
C
C
C === Puissance transmise a l'air et a la peau des colis par conduction
C ===========================================================
C
C     Pour l'air :
C       en stationnaire, toute la puissance des colis se retrouve tot ou
C       tard dans l'air (elle peut transiter par la peau des colis,
C       par les alveoles, par les murs, mais elle est extraite par
C       l'air).
C
C     Pour la peau des colis :
C       toute la puissance des colis est transmise par conduction a
C       la peau des colis.
C
C
          IF (ISCAL.EQ.ITAAMT.OR.ISCAL.EQ.ITPCMT) THEN
C
C --- Puissance des colis par unite de surface au sol
C       On prend en compte la carte de remplissage ensuite
C
C     On doit fournit un CRVEXP en Watt/CP, soit de la dimension de
C       PUICON/CP.
C
C     Comme la puissance integree sur la hauteur du colis doit etre
C       la puissance totale du colis, PUICON, on doit imposer dans
C       chaque maille une fraction de PUICON.
C     Pour une repartition uniforme en hauteur, cette fraction est
C       la hauteur de la maille VOLUME/(PTRRES*PLGRES) divisee par
C       la hauteur totale du colis.
C     Les repartitions non uniformes en hauteur sont traitees par
C       les cartes de repartition en altitude, plus bas.
C     On impose donc CRVEXP = PUICON/CP * VOLUME/(PTRRES*PLGRES)
C       puis on divise par la hauteur de colis (plus bas, lors du
C       traitement des cartes de repartition en altitude).
C
C     Noter que comme on multiplie par VOLUME/(PTRRES*PLGRES) ici
C       et qu'on divise plus bas par EPCHEL, on pourrait simplifier,
C       puisque, a priori, EPCHEL = VOLUME/(PTRRES*PLGRES)
C
C
C     Air ambiant
            IF (ISCAL.EQ.ITAAMT) THEN
              CRVEXP(IEL) = PUICON/CP0(IPHAS)
     &             * VOLUME(IEL)/(PTRRES*PLGRES)
C     Peau des colis
            ELSEIF(ISCAL.EQ.ITPCMT) THEN
              CRVEXP(IEL) = PUICON/CPCONT
     &             * VOLUME(IEL)/(PTRRES*PLGRES)
            ENDIF
C
C --- Prise en compte des cartes de remplissage (adimensionnelles)
C     Les entrepots ne sont pas uniformement remplis.
C
C     Pour représenter la carte d'encombrement des lignes de colis, on
C       repere les lignes occupees (coordonnee X, indicateur ILIGNE)
C       par leurs bornes XLIMIN et XLIMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPUIS faisant reference a la puissance)
            NZONES = NZOCAR(ILIGNE,ICPUIS)
            DO IZONE = 1, NZONES
              XLIMIN = PTRRES*VIZCAR(1,IZONE,ILIGNE,ICPUIS)
              XLIMAX = PTRRES*VIZCAR(2,IZONE,ILIGNE,ICPUIS)
              IF ((XYZCEN(1,IEL).GE.XLIMIN).AND.
     &             (XYZCEN(1,IEL).LT.XLIMAX)) THEN
                CRVEXP(IEL) = CRVEXP(IEL) * VCARTH(IZONE,ILIGNE)
              ENDIF
            ENDDO
C
C     Pour représenter la carte d'encombrement des rangees de colis, on
C       repere les rangees occupees (coordonnee Y, indicateur IRANGE)
C       par leurs bornes YRGMIN et YRGMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPUIS faisant reference a la puissance)
            NZONES = NZOCAR(IRANGE,ICPUIS)
            DO IZONE = 1, NZONES
              YRGMIN = PLGRES*VIZCAR(1,IZONE,IRANGE,ICPUIS)
              YRGMAX = PLGRES*VIZCAR(2,IZONE,IRANGE,ICPUIS)
              IF ((XYZCEN(2,IEL).GE.YRGMIN).AND.
     &             (XYZCEN(2,IEL).LT.YRGMAX)) THEN
                CRVEXP(IEL) = CRVEXP(IEL) * VCARTH(IZONE,IRANGE)
              ENDIF
            ENDDO
C
C     Pour représenter la carte d'encombrement en hauteur de colis, on
C       repere les couches occupees (coordonnee Z, indicateur IALTIT)
C       par leurs bornes ZALMIN et ZALMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPUIS faisant reference a la puissance)
C     SQZ permet d'integrer VCARTH sur les zones afin que
C       de normer VCARTH de telle sorte que la multiplication par
C       VCARTH/SQZ permette de retrouver, en integrant sur la hauteur,
C       toute la puissance du colis.
            SQZ = 0.D0
            NZONES = NZOCAR(IALTIT,ICPUIS)
            DO IZONE = 1, NZONES
              ZALMIN = EPCHEL*VIZCAR(1,IZONE,IALTIT,ICPUIS)
              ZALMAX = EPCHEL*VIZCAR(2,IZONE,IALTIT,ICPUIS)
              SQZ = SQZ + VCARTH(IZONE,IALTIT) *
     &             ( VIZCAR(2,IZONE,IALTIT,ICPUIS) -
     &             VIZCAR(1,IZONE,IALTIT,ICPUIS) )
            ENDDO
            DO IZONE = 1, NZONES
              ZALMIN = EPCHEL*VIZCAR(1,IZONE,IALTIT,ICPUIS)
              ZALMAX = EPCHEL*VIZCAR(2,IZONE,IALTIT,ICPUIS)
              IF ((XYZCEN(3,IEL).GE.ZALMIN).AND.
     &             (XYZCEN(3,IEL).LT.ZALMAX)) THEN
                CRVEXP(IEL) = CRVEXP(IEL)*
     &               VCARTH(IZONE,IALTIT)/(SQZ*EPCHEL)
              ENDIF
            ENDDO
C
C       On calcule la puissance totale injectee (Watt)
C
C     Air ambiant
            IF (ISCAL.EQ.ITAAMT) THEN
              PTOT = PTOT + CRVEXP(IEL)*CP0(IPHAS)
            ENDIF
C
          ENDIF
C
C
C === Echange convectif conteneur / air et radiatif conteneur / mur
C ===========================================================
C
C     Pour la determination de la temperature de peau des colis
C
C     On calcule le coefficient d'echange convectif efficace selon
C       la configuration (convection naturelle, forcee, mixte).
C     On ajoute le terme de rayonnement a partir du coefficient
C       equivalent HRAY calcule precedemment.
C
          IF (ISCAL.EQ.ITPCMT) THEN
C
C --- Convection naturelle verticale le long des colis
C       Van Vliet&Ross et Vliet&Liu
C
C     Rayleigh a partir du flux d'un conteneur
            RALEIZ = FRAGH*PHI0*XYZCEN(3,IEL)/XLAMB0
C     Nusselt associe
            XNSRAZ = 0.17D0 * RALEIZ**0.25D0
            XNSRAZ = MAX(XNSRAZ,(0.6D0 * RALEIZ**0.2D0))
C     Coefficient d'echange de convection naturelle       conteneur / air
C       Hypothese que le zero est a la base des colis
            HRAZ   = XNSRAZ*XLAMB0/XYZCEN(3,IEL)
C     Coefficient d'echange moyen de convection naturelle
            HMRAZ  = HMRAZ + HRAZ
C
C
C --- Convection forcee transversale dans le reseau de colis
C       Zukauskas
C
C     Reynolds (vitesse horizontale gap, diametre
            REYNOL = (UN0*PTRRES/(PTRRES-DMCONT)) * DMCONT / XNU0
C     Nusselt associe
            XNSREY = 0.35D0*(PTRRES/PLGRES)**0.2D0
     &           * 0.89D0*REYNOL**0.6D0
C     Coefficient d'echange de convection forcee          conteneur / air
            HREY   = XNSREY*XLAMB0/DMCONT
C     Coefficient d'echange moyen de convection forcee
            HMREY = HMREY + HREY
C
C --- Compteur pour les moyennes
C
            HNB = HNB + 1.D0
C
C
C --- Prise en compte du type de reseau (en ligne ou en quinconce)
C       et de la presence d'alveoles
C
C   - Alveoles
C       convection naturelle uniquement
C
            IF(IALVEO.EQ.1) THEN
C
C     Nusselt associe
              XNSALV = 0.17D0 * RALEIZ**0.25D0
C     Coefficient d'echange de convection naturelle       conteneur / air
C       Hypothese que le zero est a la base des colis
              H0 = XNSALV*XLAMB0/XYZCEN(3,IEL)
C
C
C   - Reseau en ligne (pas carre) sans alveoles
C       2/3 conv.force (face+cotes) + 1/3 conv.nat (sillage arriere)
C
            ELSEIF(ICONLG.EQ.1)THEN
              H0 = (HRAZ+2.D0*HREY)/3.D0
C
C   - Reseau en quinconce (pas triangulaire) sans alveoles
C       conv.force (pas de protection du colis amont)
C
            ELSE
              H0 = HREY
            ENDIF
C
C   - Sans alveoles
C       le coefficient d'echange est toujours au moins egal a celui
C       de la convection naturelle
            IF(IALVEO.NE.1) THEN
              H0 = MAX(HRAZ,H0)
            ENDIF
C
C    - Coefficient d'echange convectif moyen efficace
            HM0 = HM0 + H0
C
C
C --- Ajout du terme conductif et du terme radiatif
C
C       CRVEXP contient deja la puissance des colis transmise
C         par conduction
C       On veut ajouter dans CRVEXP des Watt/CP soit H0*Surface*T/CP.
C         On a ici XXAV*VOLUME = PI*DMCONT*EPCHEL, qui est la surface
C         d'échange dans la cellule.
C
C
C   - Alveoles
C       echange par convection naturelle et
C       les colis rayonnent avec le plafond et les alveoles
C
            IF(IALVEO.EQ.1) THEN
C
              CRVEXP(IEL) = CRVEXP(IEL)
     &             +H0  *(XXAV*VOLUME(IEL))
     &                  *RTP(IEL,ISCA(ITAAMT))/CPCONT
     &             +HRAY*(XXAV*VOLUME(IEL))
     &                  *RTP(IEL,ISCA(ITPPMT))/CPCONT
              CRVIMP(IEL) = -(H0+HRAY)*(XXAV*VOLUME(IEL))/CPCONT
C
C
C   - Sans alveole
C       echange par convection mixte seul
C       les colis ne rayonnent qu'entre eux, et ceci est pris en compte
C         par diffusion
C
            ELSE
C
              CRVEXP(IEL) = CRVEXP(IEL)
     &             +H0  *(XXAV*VOLUME(IEL))
     &                  *RTP(IEL,ISCA(ITAAMT))/CPCONT
              CRVIMP(IEL) = - H0      *(XXAV*VOLUME(IEL))/CPCONT
C
            ENDIF
C
          ENDIF
C
C
C
C === Echange convectif mur / air et radiatif mur / conteneur
C ===========================================================
C
C     Pour la determination de la temperature de peau des parois
C       (murs et alveoles)
C
          IF (ISCAL.EQ.ITPPMT) THEN
C
C --- Echange par convection paroi / air
C
C   - Alveoles
C       convection naturelle uniquement
C
            IF(IALVEO.EQ.1) THEN
C
C     Rayleigh a partir du flux d'un conteneur
              RALEIZ = FRAGH*PHI0*XYZCEN(3,IEL)/XLAMB0
C     Nusselt associe
              XNSALV = 0.17D0 * RALEIZ**0.25D0
C     Coefficient d'echange de convection naturelle       alveole / air
C       Hypothese que le zero est a la base des colis
              HRAZ = XNSALV*XLAMB0/XYZCEN(3,IEL)
C     Coefficient d'echange efficace alveole / air
C       Convection naturelle seule
              H0 = HRAZ
C
C
C   - Sans alveole
C       convection naturelle (Mac Adams) ou convection forcee
C
            ELSE
C
C    .Convection naturelle
C
C     Rayleigh a partir de l'ecart de temperature paroi - air
              RALEIZ =
     &             FRAGH*(RTP(IEL,ISCA(ITPPMT))-RTP(IEL,ISCA(ITAAMT)))
C     Nusselt associe (Mac Adams)
              XNSRAZ = 0.12D0*(ABS(RALEIZ))**UNTIER
C     Coefficient d'echange de convection naturelle       paroi / air
              HRAZ   = XNSRAZ*XLAMB0/XYZCEN(3,IEL)
C     Coefficient d'echange moyen de convection naturelle
              HMRAZ  = HMRAZ + HRAZ
C
C
C    .Convection forcee
C
C     Espace mur reseau
              DHMUR  = PTRRES-DMCONT
C     Reynolds associé (vitesse horizontale hors reseau, pres des murs)
              REYNOL = DHMUR*UN0/XNU0
C     Nusselt associe (Colburn)
              XNSREY = 0.023D0*REYNOL**0.8D0*PRAND0**UNTIER
C     Coefficient d'echange de convection forcee          paroi / air
              HREY   = XNSREY*XLAMB0/DHMUR
C     Coeficient d'echange moyen de convection forcee
              HMREY  = HMREY + HREY
C
C
C     .Coefficient d'echange efficace (le plus efficace est dominant)
              H0     = MAX(HRAZ,HREY)
C
            ENDIF
C
C
C   - Coefficient d'echange convectif moyen efficace
            HM0   = HM0 + H0
C   - Coefficient d'echange moyen de convection naturelle
            HMRAZ = HMRAZ + HRAZ
C   - Compteur pour les moyennes
            HNB   = HNB + 1.D0
C
C
C --- Ajout du terme conductif et du terme radiatif
C
C       On veut ajouter dans CRVEXP des Watt/CP soit H0*Surface*T/CP.
C         On calcule H0*VOLUME*T/CP et on divise par une distance
C         plus bas. On pourrait faire les choses plus directement.
C
C       Les echanges radiatifs sont les echanges avec les colis.
C
C
C   - Coefficient d'echange radiatif moyen pour les parois
            HMRAY = HMRAY + HRAY
C
C
C   - Echanges convectifs
            CRVEXP(IEL) = H0*VOLUME(IEL)*RTP(IEL,ISCA(ITAAMT))/CPPARO
            CRVIMP(IEL) =-H0*VOLUME(IEL)/CPPARO
C
C
C   - Echange radiatif dans le ciel
            IF(ICNROK.EQ.1.AND.XYZCEN(3,IEL).GT.HMAX) THEN
C
C    .Alveoles : le plafond echange avec les alveoles
              IF(IALVEO.EQ.1) THEN
                CRVEXP(IEL) = CRVEXP(IEL)
     &               +HRAY*VOLUME(IEL)
     &               *RTP(ICONRA(IEL),ISCA(ITPPMT))/CPPARO
C
C    .Sans alveole : le plafond echange avec le dessus des colis
              ELSE
                CRVEXP(IEL) = CRVEXP(IEL)
     &               +HRAY*VOLUME(IEL)
     &               *RTP(ICONRA(IEL),ISCA(ITPCMT))/CPPARO
              ENDIF
C
C
C   - Echange radiatif avec les colis pour les autres parois
C       (murs et alveoles)
            ELSE
              CRVEXP(IEL) = CRVEXP(IEL)
     &             +HRAY*VOLUME(IEL)*RTP(IEL,ISCA(ITPCMT))/CPPARO
C
            ENDIF
C
C
C   - La partie implicite est la meme pour tous
            CRVIMP(IEL) = CRVIMP(IEL) - HRAY*VOLUME(IEL)/CPPARO
C
C
C   - On divise par une distance
C       Plus exactement, on multiplie par la surface d'echange et
C       on divise par le volume (on aurait donc pu eviter de multiplier
C       par le volume).
C
C
C    .Au dessus ou au dessous des colis, le plafond ou le sol echangent
C       sur une surface de PTRRES*PLGRES
C       dans des cellules de volume EPCHEL*PTRRES*PLGRES.
C       Le rapport Surface/Volume est donc 1/EPCHEL
C
            IF(XYZCEN(3,IEL).LT.EPCHEL.OR.XYZCEN(3,IEL).GT.HMAX) THEN
              CRVEXP(IEL) = CRVEXP(IEL)/EPCHEL
              CRVIMP(IEL) = CRVIMP(IEL)/EPCHEL
C
C
C    .En stockage en alveoles, les alveoles echangent
C       sur une surface de PI*DHALVE*EPCHEL
C       dans des cellules de volume EPCHEL*PTRRES*PLGRES.
C       Le rapport Surface/Volume est donc PI*DHALVE/(PTRRES*PLGRES)
C     A priori DHALVE est positif, mais on laisse le ABS (pour le
C       cas, ou il serait defini comme la difference entre le
C       diametre des colis et des alveoles)
C
            ELSEIF(IALVEO.EQ.1) THEN
              CRVEXP(IEL) =
     &             CRVEXP(IEL)*PI*ABS(DHALVE)/(PTRRES*PLGRES)
              CRVIMP(IEL) =
     &             CRVIMP(IEL)*PI*ABS(DHALVE)/(PTRRES*PLGRES)
C
C
C    .Sans alveoles, les murs echangent
C       sur une surface de PLGRES*EPCHEL
C       dans des cellules de volume EPCHEL*PTRRES*PLGRES.
C       Le rapport Surface/Volume est donc 1/PTRRES
C       La ou il n'y a pas de mur, la temperature de mur ne signifie
C       rien, mais on peut quand meme l'y calculer (il n'y a pas de
C       phenomene de diffusion, tout est local)
C
            ELSE
              CRVEXP(IEL) = CRVEXP(IEL)/PTRRES
              CRVIMP(IEL) = CRVIMP(IEL)/PTRRES
C
            ENDIF
C
          ENDIF
C
C
C     Fin du test sur la zone de stockage
        ENDIF
      ENDDO
C
C
C === Pour l'air, sans modelisation des panaches de convection naturelle
C ===========================================================
C
      IF(ISCAL.EQ.ITAAMT) THEN
C
        IF(IMDCNT.EQ.0)THEN
C
C --- Calcul du debit enthalpique de convection naturelle
C       sortant en haut de la zone des colis : rho*W*S*Cp*DeltaT
C       (il aurait ete plus correct d'utiliser FLUMAS)
C
          DO IEL = 1, NCEL
C
            IFML  = IFMCEL(IEL   )
            ICOUL = IPRFML(IFML,1)
            IF(ICOUL.EQ.ICMTST) THEN
C
              IF(ICNROK.EQ.1.AND.XYZCEN(3,IEL).GT.HMAX) THEN
                JEL = ICONRA(IEL)
                DHS = DHS +
     &               PROPCE(JEL,IPCROM)*RTP(JEL,IW(IPHAS))
     &               *(VOLUME(JEL)/EPCHEL)
     &               *CP0(IPHAS)*(RTP(JEL,ISCA(ITAAMT))-TINIT)
              ENDIF
C
            ENDIF
          ENDDO
C
C
C --- Calcul de la puissance en kW,
C       avec correction si representation 3D en 2D
C
          PUITOT = PTOT*FRDTRA*1.D-3
C
C
C --- Calcul du debit enthalpique de convection naturelle en kW,
C       avec correction si representation 3D en 2D
C
          DEBCON = DHS*FRDTRA*1.D-3
C
C --- Impression de PUITOT et DEBCON
          IF (IRANGP.LE.0.AND.MODNTL.EQ.0) THEN
            WRITE(NFECRA,2002) PUITOT, DEBCON
          ENDIF
C
C
C === Pour l'air, modelisation des panaches de convection naturelle
C ===========================================================
C
        ELSEIF(IMDCNT.EQ.1) THEN
C
C
C --- Le debit enthalpique des panaches ne peut pas etre plus grand
C       que la puissance totale disponible ; sinon, c'est qu'il
C       a ete mal predit ou que la puissance a ete mal imposee :
C       les cartes ou la puissance d'un conteneur peuvent etre en
C       cause
C     On fait ce test ici, dans la mesure ou PTOT est calcule juste
C       au dessus.
C
          IF(DHPCNT.GE.PTOT*FRDTRA)THEN
            WRITE(NFECRA,9001) IMDCNT, DHPCNT, PTOT*FRDTRA,
     &           DHPCNT, PTOT*FRDTRA
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
C --- Modification du terme source pour la modelisation des panaches
C
C       On redistribue la puissance transmise a l'air uniquement :
C         une fraction de la puissance totale est retiree de la zone
C         de stockage et injectee directement dans le ciel au dessus
C         des conteneurs. La repartition se fait a priori de maniere
C         homogene (le panache de tous les conteneurs presents recoit
C         la meme fraction de la puissance totale), puis on applique
C         la carte de repartition thermique XY.
C
C   - Calcul des grandeurs independantes de la maille
C
C     Fraction de la puissance totale a conserver
C       dans la zone de stockage
          FPTSTO = 1.D0-DHE/PTOT
C
C     Puissance Volumique a ajouter dans la zone de Ciel au dessus
C       d'un Conteneur de puissance PUICON, divise par CP
          PVCCSC = ( (DHE/PTOT)*PUICON / (HERCNT*PTRRES*PLGRES) )
     &         / CP0(IPHAS)
C
          DO IEL = 1, NCEL
C
            IFML  = IFMCEL(IEL   )
            ICOUL = IPRFML(IFML,1)
C
            IF(ICOUL.EQ.ICMTST) THEN
C
C   - Reduction dans la zone de stockage
C
C     La puissance est reduite pour tous les conteneurs de maniere
C       homogene, en fonction du rapport entre le debit enthalpique
C       des panaches impose par l'utilisateur et la puissance totale
C       presente dans l'entrepot. Les cartes de repartition thermique
C       ayant deja ete appliquees a CRVEXP, c'est bien le rapport de
C       reduction (et non pas la valeur absolue de la reduction) qui
C       est homogene.
C
              CRVEXP(IEL) = CRVEXP(IEL)*FPTSTO
C
C
              IF ( XYZCEN(3,IEL).LT.(HERCNT+HRESO).AND.
     &             XYZCEN(3,IEL).GT.HRESO) THEN
C
C   - Report dans la zone des panaches (1/2)
C
C     Calcul de la puissance a ajouter dans la zone de ciel, au dessus
C       d'un conteneur de puissance PUICON (deja divise par CP)
C
                PCCSC = PVCCSC * VOLUME(IEL)
C
C     Les entrepots ne sont pas uniformement remplis.
C     Application de la carte thermique horizontale XY pour
C       selectionner les positions auxquelles un conteneur (ou
C       une fraction de conteneur est presente, puisque l'on
C       peut avoir des demi-conteneurs, par exemple sur un plan
C       de symetrie)
C
C     Pour représenter la carte d'encombrement des lignes de colis, on
C       repere les lignes occupees (coordonnee X, indicateur ILIGNE)
C       par leurs bornes XLIMIN et XLIMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPUIS faisant reference a la puissance)
                NZONES = NZOCAR(ILIGNE,ICPUIS)
                DO IZONE = 1, NZONES
                  XLIMIN = PTRRES*VIZCAR(1,IZONE,ILIGNE,ICPUIS)
                  XLIMAX = PTRRES*VIZCAR(2,IZONE,ILIGNE,ICPUIS)
                  IF ((XYZCEN(1,IEL).GE.XLIMIN).AND.
     &                 (XYZCEN(1,IEL).LT.XLIMAX)) THEN
                    PCCSC = PCCSC *
     &                   VCARTH(IZONE,ILIGNE)
                  ENDIF
                ENDDO
C
C     Pour représenter la carte d'encombrement des rangees de colis, on
C       repere les rangees occupees (coordonnee Y, indicateur IRANGE)
C       par leurs bornes YRGMIN et YRGMAX qui sont stockees dans VIZCAR
C       (avec l'indicateur ICPUIS faisant reference a la puissance)
                NZONES = NZOCAR(IRANGE,ICPUIS)
                DO IZONE = 1, NZONES
                  YRGMIN = PLGRES*VIZCAR(1,IZONE,IRANGE,ICPUIS)
                  YRGMAX = PLGRES*VIZCAR(2,IZONE,IRANGE,ICPUIS)
                  IF ((XYZCEN(2,IEL).GE.YRGMIN).AND.
     &                 (XYZCEN(2,IEL).LT.YRGMAX)) THEN
                    PCCSC = PCCSC *
     &                   VCARTH(IZONE,IRANGE)
                  ENDIF
                ENDDO
C
C     Ajout de la puissance des panaches a la puissance deja existante
C
                CRVEXP(IEL) = CRVEXP(IEL) + PCCSC
C
C     Fin si zone de panache
              ENDIF
C
C     Fin si zone de stockage
            ENDIF
C
C
C --- Calcul de la puissance totale incluant la modelisation
C       des panaches de convection natuelle
C       (normalement la puissance totale n'a pas ete modifiee)
C
            PTCN = PTCN + CRVEXP(IEL)*CP0(IPHAS)
C
C     Fin boucle NCEL
          ENDDO
C
C --- Calcul et impression de la puissance totale incluant la modelisation
C       des panaches de convection naturelle
C       en kW, avec correction si representation 3D en 2D
C
          PUITOT = PTCN*FRDTRA*1.D-3
          IF ((IRANGP.LE.0).AND.(MODNTL.EQ.0)) THEN
            WRITE(NFECRA,2003) PUITOT
          ENDIF
C
C     Fin si on modelise les panaches
        ENDIF
C     Fin si scalaire = air ambiant
      ENDIF
C
C
C
C === Impressions pour la temperature de peau des colis et des parois
C ===========================================================
C
C --- Valeur fictive : pas de convection forcee en alveole
C     (valeur bizarre, mais uniquement pour affichage)
      IF (IALVEO.EQ.1) THEN
        HMREY = -1.00001D10
      ENDIF
C
C --- Temperature de peau des colis
C
      IF(ISCAL.EQ.ITPCMT) THEN
C
C     Si rien a moyenner, on prend une valeur negative
C       (mais les numerateurs seront nuls, normalement)
        IF (HNB.LT.EPZERO) THEN
          HNB = -EPZERO
        ENDIF
C
C     Stockage en common pour impression ulterieure
        CFECCA = HM0/HNB
C
C     Impressions des coefficients d'echange moyens
        IF ((IRANGP.LE.0).AND.(MODNTL.EQ.0)) THEN
          WRITE(NFECRA,2004) HMREY/HNB, HMRAZ/HNB, HM0/HNB
        ENDIF
C
      ENDIF
C
C --- Temperature de peau des parois
C
      IF(ISCAL.EQ.ITPPMT) THEN
C
C     Si rien a moyenner, on prend une valeur negative
C       (mais les numerateurs seront nuls, normalement)
        IF (HNB.LT.EPZERO) THEN
          HNB = -EPZERO
        ENDIF
C
C     Stockage en common pour impression ulterieure
        CFECMA = HM0/HNB
C
C     Impressions des coefficients d'echange moyens
        IF ((IRANGP.LE.0).AND.(MODNTL.EQ.0)) THEN
          WRITE(NFECRA,2005)
     &         HMREY/HNB, HMRAZ/HNB, HMRAY/HNB, HM0/HNB
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
C
 2001 FORMAT (/,3X,'** INFORMATIONS SUR MATISSE, VARIABLE ',A8,/,
     &          3X,'   -------------------------------------------')
 2002 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Puissance (sans modelisation panaches)                  ',/,
     &'mati --------------------------------------                  ',/,
     &'mati Puissance totale                      ', E12.5  ,'    kW',/,
     &'mati Debit enthalpique de conv. naturelle  ', E12.5  ,'    kW',/,
     &'mati --------------------------------------------------------',/)
 2003 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Puissance (avec modelisation panaches)                  ',/,
     &'mati --------------------------------------                  ',/,
     &'mati Puissance totale avec erosion panaches', E12.5  ,'    kW',/,
     &'mati --------------------------------------------------------',/)
 2004 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Coeff. d echange moyens (calcul T° de peau des colis)   ',/,
     &'mati -----------------------------------------------------   ',/,
     &'mati Convection forcee                   ', E12.5  ,' W/m2/°C',/,
     &'mati Convection naturelle                ', E12.5  ,' W/m2/°C',/,
     &'mati Efficace                            ', E12.5  ,' W/m2/°C',/,
     &'mati --------------------------------------------------------',/)
 2005 FORMAT (
     &'mati --------------------------------------------------------',/,
     &'mati Coeff. d echange moyens (calcul T° de peau des parois)  ',/,
     &'mati ------------------------------------------------------  ',/,
     &'mati Convection forcee                   ', E12.5  ,' W/m2/°C',/,
     &'mati Convection naturelle                ', E12.5  ,' W/m2/°C',/,
     &'mati Rayonnement                         ', E12.5  ,' W/m2/°C',/,
     &'mati Efficace                            ', E12.5  ,' W/m2/°C',/,
     &'mati --------------------------------------------------------',/)
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTTSSC)                      ',/,
     &'@    *********                                               ',/,
     &'@      DEBIT ENTHALPIQUE DES PANACHES INADAPTE.              ',/,
     &'@                                                            ',/,
     &'@    La modelisation des panaches est activee                ',/,
     &'@      avec IMDCNT = ',I10   ,'                              ',/,
     &'@    Le debit enthalpique des panaches de convection         ',/,
     &'@      naturelle doit etre strictement inferieur a la        ',/,
     &'@      puissance totale des colis presents dans l''entrepot  ',/,
     &'@      puisqu''il en represente une partie.                  ',/,
     &'@                                                            ',/,
     &'@    Le debit enthalpique prescrit est  DHPCNT  =',E12.5      ,/,
     &'@    La puissance totale reelle est PTOT*FRDTRA =',E12.5      ,/,
     &'@                                                            ',/,
     &'@    L''inegalite suivante n''est pas verifiee :             ',/,
     &'@          DHPCNT    <    PTOT*FRDTRA                        ',/,
     &'@       ',E12.5   ,'     ',E12.5   ,'                        ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Desactiver la modelisation des panaches ou              ',/,
     &'@      modifier la valeur du debit enthalpique prescrit ou   ',/,
     &'@      modifier la puissance totale (carte de remplissage ou ',/,
     &'@      puissance des colis).                                 ',/,
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
