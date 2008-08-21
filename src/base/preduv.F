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
                        SUBROUTINE PREDUV
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IAPPEL ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ITERNS , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM , IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   FLUMAS , FLUMAB ,
     &   TSLAGR , COEFA  , COEFB  ,
     &   CKUPDC , SMACEL , FRCXT  ,
     &   TRAVA  , XIMPA  , UVWK   , DFRCXT , TPUCOU , TRAV   ,
     &   VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , XNORMP , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS N-S 1 PHASE INCOMPRESSIBLE OU RO VARIABLE
CFONC SUR UN PAS DE TEMPS (CONVECTION/DIFFUSION - PRESSION /CONTINUITE)
CFONC
CFONC AU PREMIER APPEL,  ON EFFECTUE LA PREDICITION DES VITESSES
CFONC                ET  ON CALCULE UN ESTIMATEUR SUR LA VITESSE PREDITE
CFONC
CFONC AU DEUXIEME APPEL, ON CALCULE UN ESTIMATEUR GLOBAL
CFONC                POUR NAVIER-STOKES :
CFONC    ON UTILISE TRAV, SMBR ET LES TABLEAUX DE TRAVAIL
CFONC    ON APPELLE BILSC2 AU LIEU DE CODITS
CFONC    ON REMPLIT LE PROPCE ESTIMATEUR IESTOT
CFONC    CE DEUXIEME APPEL INTERVIENT DANS NAVSTO APRES RESOLP
CFONC    LORS DE CE DEUXIEME APPEL
CFONC      RTPA ET RTP SONT UN UNIQUE TABLEAU (= RTP)
CFONC      LE FLUX DE MASSE EST LE FLUX DE MASSE DEDUIT DE LA VITESSE
CFONC       AU CENTRE CONTENUE DANS RTP
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
CARGU ! IAPPEL       ! E  !  -> ! NUMERO D'APPEL DU SOUS PROGRAMME     !
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
CARGU ! ITERNS       ! E  !  -> ! NUMERO D'ITERATION SUR NAVSTO        !
CARGU ! NCEPDP       ! E  !  -> ! NOMBRE DE CELLULES AVEC PDC          !
CARGU ! NCKPDP       ! E  !  -> ! NBR DE COEF DU TENSEUR DE PDC (3 OU 6!
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
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
CARGU ! IFACLG(2,NFAC! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! IRESPR(NCELET! TE !  -  ! TAB ENTIER MULTIGRILLE               !
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
CARGU ! FLUMAS       ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU !  (NFAC  )    !    !     !   (DISTINCTION IAPPEL=1 OU 2)        !
CARGU ! FLUMAB       ! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU !  (NFABOR  )  !    !     !    (DISTINCTION IAPPEL=1 OU 2)       !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! TRAVA,XIMPA  ! TR ! <-> ! TABLEAU DE TRAVAIL POUR COUPLAGE     !
CARGU ! UVWK         ! TR ! <-> ! TABLEAU DE TRAVAIL POUR COUPLAGE U/P !
CARGU !              !    !     ! SERT A STOCKER LA VITESSE DE         !
CARGU !              !    !     ! L'ITERATION PRECEDENTE               !
CARGU !DFRCXT(NCELET,! TR ! <-  ! VARIATION DE FORCE EXTERIEURE        !
CARGU !   3,NPHAS)   !    !     !  GENERANT LA PRESSION HYDROSTATIQUE  !
CARGU ! TPUCOU       ! TR ! <-  ! COUPLAGE VITESSE PRESSION            !
CARGU ! (NCELEL,NDIM)!    !     !                                      !
CARGU ! TRAV(NCELET,3! TR ! <-  ! SMB QUI SERVIRA POUR NORMALISATION   !
CARGU !              !    !     !  DANS RESOLP                         !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! VISCFI(NFAC) ! TR !  -  ! IDEM VISCF POUR INCREMENTS           !
CARGU ! VISCBI(NFABOR! TR !  -  ! IDEM VISCB POUR INCREMENTS           !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR  (NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! XNORMP(NCELET! TR !  -> ! TABLEAU REEL POUR NORME RESOLP       !
CARGU ! COEFU(NFAB,3)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TSLAGR(NCELET! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !   NTERSL)    !    !     !   LAGRANGIEN                         !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "matiss.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , IAPPEL
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , ITERNS
      INTEGER          NCEPDP , NCKPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
      INTEGER          IFACLG(2,NFAC), IRESPR(NCELET)
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
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS), DFRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION TRAVA(NCELET,NDIM,NPHAS)
      DOUBLE PRECISION XIMPA(NCELET,NDIM,NPHAS),UVWK(NCELET,NDIM,NPHAS)
      DOUBLE PRECISION TPUCOU(NCELET,NDIM), TRAV(NCELET,3)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION VISCFI(NFAC), VISCBI(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET)
      DOUBLE PRECISION SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION XNORMP(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IEL   , IELPDC, IFAC  , IVAR  , ISOU
      INTEGER          ICCOCG, INC   , INIT  , IPHYDP, II    , ISQRT
      INTEGER          IRESLP, NSWRGP, IMLIGP, IWARNP, IPPT  , IPP
      INTEGER          IPRIPH, IKIPH , IUIPH , IVIPH , IWIPH
      INTEGER                  ICLIPR, ICLIUP, ICLIVP, ICLIWP
      INTEGER                          ICLIK , ICLVAR, ICLVAF
      INTEGER          IPCROM, IPCROA, IPCROO , IPCVIS, IPCVST
      INTEGER          ICONVP, IDIFFP, NDIRCP, NITMAP, NSWRSP
      INTEGER          IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IDIMTE, ITENSO
      INTEGER          IESPRP, IESTOP
      INTEGER          IPTSNA
      INTEGER          IFLMB0, NSWRP , IMASPE, IISMPH, IPBROM
      INTEGER          IDIAEX, IDTVA0
      INTEGER          MAXELT, ILS
      DOUBLE PRECISION RNORM , VITNOR
      DOUBLE PRECISION ROMVOM, RO0IPH, DROM
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, RELAXP, BLENCP, EPSILP
      DOUBLE PRECISION VIT1  , VIT2  , VIT3, XKB, PIP, PFAC, PFAC1
      DOUBLE PRECISION CPDC11, CPDC22, CPDC33, CPDC12, CPDC13, CPDC23
      DOUBLE PRECISION D2S3  , THETAP, THETP1, THETS , DTSROM
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPRIPH = IPR(IPHAS)
      IUIPH  = IU(IPHAS)
      IVIPH  = IV(IPHAS)
      IWIPH  = IW(IPHAS)
      IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &     .OR. ITURB(IPHAS).EQ.60) THEN
        IKIPH  = IK(IPHAS)
      ENDIF
C
      ICLIPR = ICLRTP(IPRIPH,ICOEF)
      ICLIUP = ICLRTP(IUIPH ,ICOEF)
      ICLIVP = ICLRTP(IVIPH ,ICOEF)
      ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
      IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &     .OR. ITURB(IPHAS).EQ.60) THEN
        ICLIK  = ICLRTP(IKIPH ,ICOEF)
      ENDIF
C
      RO0IPH = RO0(IPHAS)
C
C     Reperage de rho au bord
      IPBROM = IPPROB(IROM  (IPHAS))
C     Reperage de rho courant (ie en cas d'extrapolation rho^n+1/2)
      IPCROM = IPPROC(IROM  (IPHAS))
C     Reperage de rho^n en cas d'extrapolation
      IF(IROEXT(IPHAS).GT.0) THEN
        IPCROA = IPPROC(IROMA(IPHAS))
      ELSE
        IPCROA = 0
      ENDIF
C
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C     Theta relatif aux termes sources explicites
      THETS  = THETSN(IPHAS)
      IF(ISNO2T(IPHAS).GT.0) THEN
        IPTSNA = IPPROC(ITSNSA(IPHAS))
      ELSE
        IPTSNA = 0
      ENDIF
C
C
C=======================================================================
C 2.  GRADIENT DE PRESSION ET GRAVITE
C=======================================================================
C
C ---> PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE :
C       UNIQUEMENT AU PREMIER APPEL (I.E. POUR LE CALCUL NORMAL,
C       LE DEUXIEME APPEL EST DESTINE A UN CALCUL D'ESTIMATEUR)
C
      IF (IAPPEL.EQ.1.AND.IPHYDR.EQ.1) THEN
C
c force ext au pas de temps precedent :
C     FRCXT a ete initialise a zero
C     (est deja utilise dans typecl, et est mis a jour a la fin
C     de navsto)
C
        DO IEL = 1, NCEL
C
C variation de force (utilise dans resolp)
          DROM = (PROPCE(IEL,IPCROM)-RO0IPH)
          DFRCXT(IEL,1,IPHAS) = DROM*GX - FRCXT(IEL,1,IPHAS)
          DFRCXT(IEL,2,IPHAS) = DROM*GY - FRCXT(IEL,2,IPHAS)
          DFRCXT(IEL,3,IPHAS) = DROM*GZ - FRCXT(IEL,3,IPHAS)
        ENDDO
C     Ajout eventuel des pertes de charges
        IF (NCEPDP.GT.0) THEN
          DO IELPDC = 1, NCEPDP
            IEL=ICEPDC(IELPDC)
            VIT1   = RTPA(IEL,IUIPH)
            VIT2   = RTPA(IEL,IVIPH)
            VIT3   = RTPA(IEL,IWIPH)
            CPDC11 = CKUPDC(IELPDC,1)
            CPDC22 = CKUPDC(IELPDC,2)
            CPDC33 = CKUPDC(IELPDC,3)
            IF (NCKPDP.EQ.6) THEN
              CPDC12 = CKUPDC(IELPDC,4)
              CPDC13 = CKUPDC(IELPDC,5)
              CPDC23 = CKUPDC(IELPDC,6)
            ELSE
              CPDC12 = 0.D0
              CPDC13 = 0.D0
              CPDC23 = 0.D0
            ENDIF
            DFRCXT(IEL,1,IPHAS) = DFRCXT(IEL,1,IPHAS)
     &           -PROPCE(IEL,IPCROM)*(
     &           CPDC11*VIT1+CPDC12*VIT2+CPDC13*VIT3)
            DFRCXT(IEL,2,IPHAS) = DFRCXT(IEL,2,IPHAS)
     &           -PROPCE(IEL,IPCROM)*(
     &           CPDC12*VIT1+CPDC22*VIT2+CPDC23*VIT3)
            DFRCXT(IEL,3,IPHAS) = DFRCXT(IEL,3,IPHAS)
     &           -PROPCE(IEL,IPCROM)*(
     &           CPDC13*VIT1+CPDC23*VIT2+CPDC33*VIT3)
          ENDDO
        ENDIF
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (DFRCXT(1,1,IPHAS))
C         ===========
          CALL PARCOM (DFRCXT(1,2,IPHAS))
C         ===========
          CALL PARCOM (DFRCXT(1,3,IPHAS))
C         ===========
        ENDIF
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 1
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    DFRCXT(1,1,IPHAS),DFRCXT(1,1,IPHAS),DFRCXT(1,1,IPHAS),
     &    DFRCXT(1,2,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,2,IPHAS),
     &    DFRCXT(1,3,IPHAS),DFRCXT(1,3,IPHAS),DFRCXT(1,3,IPHAS))
        ENDIF
C
      ENDIF
C
C
C
C ---> PRISE EN COMPTE DU GRADIENT DE PRESSION
C
      ICCOCG = 1
      INC    = 1
      NSWRGP = NSWRGR(IPRIPH)
      IMLIGP = IMLIGR(IPRIPH)
      IWARNP = IWARNI(IPRIPH)
      EPSRGP = EPSRGR(IPRIPH)
      CLIMGP = CLIMGR(IPRIPH)
      EXTRAP = EXTRAG(IPRIPH)
C
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPRIPH , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCXT(1,1,IPHAS), FRCXT(1,2,IPHAS), FRCXT(1,3,IPHAS),
     &   RTPA(1,IPRIPH)  , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C    Calcul des efforts aux parois (partie 2/5), si demande
C    La pression a la face est calculee comme dans gradrc/gradmc
C    et on la transforme en pression totale
C    On se limite a la premiere iteration (pour faire simple par
C      rapport a la partie issue de condli, hors boucle)
      IF (INEEDF.EQ.1 .AND. ITERNS.EQ.1) THEN
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          II = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(II+1)
          DIIPBY = RA(II+2)
          DIIPBZ = RA(II+3)
          PIP =  RTPA(IEL,IPRIPH)
     &         +DIIPBX*W1(IEL) +DIIPBY*W2(IEL)
     &         +DIIPBZ*W3(IEL)
          PFAC = COEFA(IFAC,ICLIPR) +COEFB(IFAC,ICLIPR)*PIP
          PFAC1= RTPA(IEL,IPRIPH)
     &         +(CDGFBO(1,IFAC)-XYZCEN(1,IEL))*W1(IEL)
     &         +(CDGFBO(2,IFAC)-XYZCEN(2,IEL))*W2(IEL)
     &         +(CDGFBO(3,IFAC)-XYZCEN(3,IEL))*W3(IEL)
          PFAC = COEFB(IFAC,ICLIPR)*(EXTRAG(IPRIPH)*PFAC1
     &         +(1.D0-EXTRAG(IPRIPH))*PFAC)
     &         +(1.D0-COEFB(IFAC,ICLIPR))*PFAC
     &         + RO0(IPHAS)*(GX*(CDGFBO(1,IFAC)-XYZP0(1,IPHAS))
     &         + GY*(CDGFBO(2,IFAC)-XYZP0(2,IPHAS))
     &         + GZ*(CDGFBO(3,IFAC)-XYZP0(3,IPHAS)) )
     &         - PRED0(IPHAS)
c on ne rajoute pas P0, pour garder un maximum de precision
c     &         + P0(IPHAS)
          DO ISOU = 1, 3
            RA(IFORBR+(IFAC-1)*NDIM + ISOU-1) =
     &           RA(IFORBR+(IFAC-1)*NDIM + ISOU-1)
     &           + PFAC*SURFBO(ISOU,IFAC)
          ENDDO
        ENDDO
      ENDIF
C
C ---> RESIDU DE NORMALISATION POUR RESOLP
C     Test d'un residu de normalisation de l'etape de pression
C       plus comprehensible = div(rho u* + dt gradP^(n))-Gamma
C       i.e. second membre du systeme en pression hormis la partie
C       pression (sinon a convergence, on tend vers 0)
C       Represente les termes que la pression doit equilibrer
C     On calcule ici div(rho dt/rho gradP^(n)) et on complete a la fin
C       avec  div(rho u*)
C     Pour grad P^(n) on suppose que des CL de Neumann homogenes
C       s'appliquent partout : on peut donc utiliser les CL de la
C       vitesse pour u*+dt/rho gradP^(n). Comme on calcule en deux fois,
C       on utilise les CL de vitesse homogenes pour dt/rho gradP^(n)
C       ci-dessous et les CL de vitesse completes pour u* a la fin.
C
      IF(IAPPEL.EQ.1.AND.IRNPNW.EQ.1) THEN
C
C     Calcul de dt/rho*grad P
        DO IEL = 1, NCEL
          DTSROM = DT(IEL)/PROPCE(IEL,IPCROM)
          TRAV(IEL,1) = W1(IEL)*DTSROM
          TRAV(IEL,2) = W2(IEL)*DTSROM
          TRAV(IEL,3) = W3(IEL)*DTSROM
        ENDDO
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (TRAV(1,1))
C         ===========
          CALL PARCOM (TRAV(1,2))
C         ===========
          CALL PARCOM (TRAV(1,3))
C         ===========
        ENDIF
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 1
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    TRAV(1,1),TRAV(1,1),TRAV(1,1),
     &    TRAV(1,2),TRAV(1,2),TRAV(1,2),
     &    TRAV(1,3),TRAV(1,3),TRAV(1,3))
        ENDIF
C
C     Calcul de rho dt/rho*grad P.n aux faces
C       Pour gagner du temps, on ne reconstruit pas.
        INIT   = 1
        INC    = 0
        ICCOCG = 1
        IFLMB0 = 1
        NSWRP  = 1
        IMLIGP = IMLIGR(IUIPH )
        IWARNP = IWARNI(IPRIPH)
        EPSRGP = EPSRGR(IUIPH )
        CLIMGP = CLIMGR(IUIPH )
        EXTRAP = EXTRAG(IUIPH )
C
        IMASPE = 1
C
        IISMPH = IISYMP+NFABOR*(IPHAS-1)
C
        CALL INIMAS
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRP  , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   TRAV(1,1) , TRAV(1,2) , TRAV(1,3) ,
     &   COEFA(1,ICLIUP), COEFA(1,ICLIVP), COEFA(1,ICLIWP),
     &   COEFB(1,ICLIUP), COEFB(1,ICLIVP), COEFB(1,ICLIWP),
     &   VISCF  , VISCB  ,
     &   W4     , W5     , W6     , W7     , W8     , W9     ,
     &   SMBR   , DRTP   , ROVSDT , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
C     Calcul de div(rho dt/rho*grad P)
        INIT = 1
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                                IFACEL,IFABOR,VISCF,VISCB,XNORMP)
C
C     Ajout de -Gamma
        IF (NCESMP.GT.0) THEN
          DO II = 1, NCESMP
            IEL = ICETSM(II)
            XNORMP(IEL) = XNORMP(IEL)-VOLUME(IEL)*SMACEL(II,IPRIPH)
          ENDDO
        ENDIF
C
C     On conserve XNORMP, on complete avec u* a la fin et
C       on le transfere a resolp
C
      ENDIF
C
C
C     Au premier appel, TRAV est construit directement ici.
C     Au second  appel (estimateurs), TRAV contient deja
C       l'increment temporel.
C     On pourrait fusionner en initialisant TRAV a zero
C       avant le premier appel, mais ca fait des operations en plus.
C
C     Remarques :
C       rho g sera a l'ordre 2 s'il est extrapole.
C       si on itere sur navsto, ca ne sert a rien de recalculer rho g a
C         chaque fois (ie on pourrait le passer dans trava) mais ce n'est
C         pas cher.
      IF(IAPPEL.EQ.1) THEN
C
        IF (IPHYDR.EQ.1) THEN
          DO IEL = 1, NCEL
            TRAV(IEL,1) = (FRCXT(IEL,1,IPHAS) - W1(IEL)) * VOLUME(IEL)
            TRAV(IEL,2) = (FRCXT(IEL,2,IPHAS) - W2(IEL)) * VOLUME(IEL)
            TRAV(IEL,3) = (FRCXT(IEL,3,IPHAS) - W3(IEL)) * VOLUME(IEL)
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            DROM = (PROPCE(IEL,IPCROM)-RO0IPH)
            TRAV(IEL,1) = ( DROM*GX - W1(IEL) ) * VOLUME(IEL)
            TRAV(IEL,2) = ( DROM*GY - W2(IEL) ) * VOLUME(IEL)
            TRAV(IEL,3) = ( DROM*GZ - W3(IEL) ) * VOLUME(IEL)
          ENDDO
        ENDIF
C
      ELSEIF(IAPPEL.EQ.2) THEN
C
        IF (IPHYDR.EQ.1) THEN
          DO IEL = 1, NCEL
            TRAV(IEL,1) =
     &           TRAV(IEL,1) + ( FRCXT(IEL,1,IPHAS) -
     &                           W1(IEL) )*VOLUME(IEL)
            TRAV(IEL,2) =
     &           TRAV(IEL,2) + ( FRCXT(IEL,2,IPHAS) -
     &                           W2(IEL) )*VOLUME(IEL)
            TRAV(IEL,3) =
     &           TRAV(IEL,3) + ( FRCXT(IEL,3,IPHAS) -
     &                           W3(IEL) )*VOLUME(IEL)
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            DROM = (PROPCE(IEL,IPCROM)-RO0IPH)
            TRAV(IEL,1) =
     &           TRAV(IEL,1) + ( DROM*GX - W1(IEL) )*VOLUME(IEL)
            TRAV(IEL,2) =
     &           TRAV(IEL,2) + ( DROM*GY - W2(IEL) )*VOLUME(IEL)
            TRAV(IEL,3) =
     &           TRAV(IEL,3) + ( DROM*GZ - W3(IEL) )*VOLUME(IEL)
          ENDDO
        ENDIF
C
      ENDIF
C
C   Pour IAPPEL = 1 (ie appel standard sans les estimateurs)
C     TRAV rassemble les termes sources  qui seront recalcules
C       a toutes les iterations sur navsto
C     Si on n'itere pas sur navsto et qu'on n'extrapole pas les
C       termes sources, TRAV contient tous les termes sources
C       jusqu'au basculement dans SMBR
C     A ce niveau, TRAV contient -grad P et rho g
C       P est suppose pris a n+1/2
C       rho est eventuellement interpole a n+1/2
C
C
C ---> INITIALISATION DU TABLEAU TRAVA et PROPCE AU PREMIER PASSAGE
C     (A LA PREMIERE ITER SUR NAVSTO)
C
C     TRAVA rassemble les termes sources qu'il suffit de calculer
C       a la premiere iteration sur navsto quand il y a plusieurs iter.
C     Quand il n'y a qu'une iter, on cumule directement dans TRAV
C       ce qui serait autrement alle dans TRAVA
C     PROPCE rassemble les termes sources explicites qui serviront
C       pour le pas de temps suivant en cas d'extrapolation (plusieurs
C       iter sur navsto ou pas)
C
C     A la premiere iter sur navsto
      IF(ITERNS.EQ.1) THEN
C
C       Si on   extrapole     les T.S. : -theta*valeur precedente
          IF(ISNO2T(IPHAS).GT.0) THEN
C         S'il n'y a qu'une    iter : TRAV  incremente
            IF(NTERUP.EQ.1) THEN
              DO II = 1, NDIM
                DO IEL = 1, NCEL
                  TRAV (IEL,II)        = TRAV (IEL,II)
     &                 - THETS*PROPCE(IEL,IPTSNA+II-1)
                ENDDO
              ENDDO
C         S'il   y a plusieurs iter : TRAVA initialise
            ELSE
              DO II = 1, NDIM
                DO IEL = 1, NCEL
                  TRAVA(IEL,II,IPHAS)  =
     &                 - THETS*PROPCE(IEL,IPTSNA+II-1)
                ENDDO
              ENDDO
            ENDIF
C         Et on initialise PROPCE pour le remplir ensuite
            DO II = 1, NDIM
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSNA+II-1) = 0.D0
              ENDDO
            ENDDO
C
C       Si on n'extrapole pas les T.S. : pas de PROPCE
          ELSE
C         S'il   y a plusieurs iter : TRAVA initialise
C           sinon TRAVA n'existe pas
            IF(NTERUP.GT.1) THEN
              DO II = 1, NDIM
                DO IEL = 1, NCEL
                  TRAVA(IEL,II,IPHAS)  = 0.D0
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C
        ENDIF
C
C ---> 2/3 RHO * GRADIENT DE K SI k-epsilon ou k-omega
C      NB : ON NE PREND PAS LE GRADIENT DE (RHO K), MAIS
C           CA COMPLIQUERAIT LA GESTION DES CL ...
C     On peut se demander si l'extrapolation en temps sert a
C       quelquechose
C
C     Ce terme explicite est calcule une seule fois,
C       a la premiere iter sur navsto : il va dans PROPCE si on
C       doit l'extrapoler en temps ; il va dans TRAVA si on n'extrapole
C       pas mais qu'on itere sur navsto. Il va dans TRAV si on
C       n'extrapole pas et qu'on n'itere pas sur navsto.
      IF( (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &     .OR. ITURB(IPHAS).EQ.60)
     &       .AND.IGRHOK(IPHAS).EQ.1.AND.ITERNS.EQ.1)THEN
        ICCOCG = 1
        INC    = 1
        NSWRGP = NSWRGR(IKIPH)
        IMLIGP = IMLIGR(IKIPH)
        EPSRGP = EPSRGR(IKIPH)
        CLIMGP = CLIMGR(IKIPH)
        EXTRAP = EXTRAG(IKIPH)
C
        IWARNP = IWARNI(IUIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IKIPH  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W6     , W6     , W6     ,
     &   RTPA(1,IKIPH)   , COEFA(1,ICLIK)  , COEFB(1,ICLIK)  ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        D2S3 = 2.D0/3.D0
C
C     Si on extrapole les termes source en temps : PROPCE
        IF(ISNO2T(IPHAS).GT.0) THEN
C       Calcul de rho^n grad k^n      si rho non extrapole
C                 rho^n grad k^n si rho     extrapole
          IPCROO = IPCROM
          IF(IROEXT(IPHAS).GT.0) IPCROO = IPCROA
          DO IEL = 1, NCEL
            ROMVOM = -PROPCE(IEL,IPCROO)*VOLUME(IEL)*D2S3
            PROPCE(IEL,IPTSNA  )=PROPCE(IEL,IPTSNA  )+W1(IEL)*ROMVOM
            PROPCE(IEL,IPTSNA+1)=PROPCE(IEL,IPTSNA+1)+W2(IEL)*ROMVOM
            PROPCE(IEL,IPTSNA+2)=PROPCE(IEL,IPTSNA+2)+W3(IEL)*ROMVOM
          ENDDO
C     Si on n'extrapole pas les termes sources en temps : TRAV ou TRAVA
        ELSE
          IF(NTERUP.EQ.1) THEN
            DO IEL = 1, NCEL
              ROMVOM = -PROPCE(IEL,IPCROM)*VOLUME(IEL)*D2S3
              TRAV (IEL,1)       =
     &        TRAV (IEL,1)       + W1(IEL) * ROMVOM
              TRAV (IEL,2)       =
     &        TRAV (IEL,2)       + W2(IEL) * ROMVOM
              TRAV (IEL,3)       =
     &        TRAV (IEL,3)       + W3(IEL) * ROMVOM
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              ROMVOM = -PROPCE(IEL,IPCROM)*VOLUME(IEL)*D2S3
              TRAVA(IEL,1,IPHAS) =
     &        TRAVA(IEL,1,IPHAS) + W1(IEL) * ROMVOM
              TRAVA(IEL,2,IPHAS) =
     &        TRAVA(IEL,2,IPHAS) + W2(IEL) * ROMVOM
              TRAVA(IEL,3,IPHAS) =
     &        TRAVA(IEL,3,IPHAS) + W3(IEL) * ROMVOM
            ENDDO
          ENDIF
        ENDIF
C
C    Calcul des efforts aux parois (partie 3/5), si demande
        IF (INEEDF.EQ.1) THEN
          DO IFAC = 1, NFABOR
            IEL = IFABOR(IFAC)
            II = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(II+1)
            DIIPBY = RA(II+2)
            DIIPBZ = RA(II+3)
            XKB = RTPA(IEL,IKIPH) + DIIPBX*W1(IEL)
     &           + DIIPBY*W2(IEL) + DIIPBZ*W3(IEL)
            XKB = COEFA(IFAC,ICLIK)+COEFB(IFAC,ICLIK)*XKB
            XKB = D2S3*PROPCE(IEL,IPCROM)*XKB
            DO ISOU = 1, 3
              RA(IFORBR+(IFAC-1)*NDIM + ISOU-1) =
     &             RA(IFORBR+(IFAC-1)*NDIM + ISOU-1)
     &             + XKB*SURFBO(ISOU,IFAC)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C
C ---> TERMES DE GRADIENT TRANSPOSE
C
C     Ce terme explicite est calcule une seule fois,
C       a la premiere iter sur navsto :
C       si on extrapole il va dans PROPCE,
C       sinon si on itere sur navsto dans TRAVA
C             sinon                  dans TRAV
      IF (IVISSE(IPHAS).EQ.1.AND.ITERNS.EQ.1) THEN
C
C     On utilise temporairement TRAV comme tableau de travail.
C     Son contenu est stocke dans W7, W8 et W9 jusqu'apres vissec
        DO IEL = 1,NCEL
          W7(IEL) = TRAV(IEL,1)
          W8(IEL) = TRAV(IEL,2)
          W9(IEL) = TRAV(IEL,3)
          TRAV(IEL,1) = 0.D0
          TRAV(IEL,2) = 0.D0
          TRAV(IEL,3) = 0.D0
        ENDDO
C
        CALL VISSEC
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   TRAV   ,
C        ------
     &   VISCF  , VISCB  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si on extrapole les termes source en temps :
C       PROPCE recoit les termes de gradient transpose et
C       TRAV retrouve sa valeur
        IF(ISNO2T(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSNA  ) = PROPCE(IEL,IPTSNA  ) + TRAV(IEL,1)
            PROPCE(IEL,IPTSNA+1) = PROPCE(IEL,IPTSNA+1) + TRAV(IEL,2)
            PROPCE(IEL,IPTSNA+2) = PROPCE(IEL,IPTSNA+2) + TRAV(IEL,3)
            TRAV(IEL,1)  = W7(IEL)
            TRAV(IEL,2)  = W8(IEL)
            TRAV(IEL,3)  = W9(IEL)
          ENDDO
C
C     Si on n'extrapole pas les termes source en temps :
C       Si on itere sur navsto
C         TRAVA recoit les termes de gradient transpose et
C         TRAV retrouve sa valeur
        ELSE
          IF(NTERUP.GT.1) THEN
            DO IEL = 1, NCEL
              TRAVA(IEL,1,IPHAS) = TRAVA(IEL,1,IPHAS) + TRAV(IEL,1)
              TRAVA(IEL,2,IPHAS) = TRAVA(IEL,2,IPHAS) + TRAV(IEL,2)
              TRAVA(IEL,3,IPHAS) = TRAVA(IEL,3,IPHAS) + TRAV(IEL,3)
              TRAV(IEL,1)  = W7(IEL)
              TRAV(IEL,2)  = W8(IEL)
              TRAV(IEL,3)  = W9(IEL)
            ENDDO
C       Si on n'itere pas sur navsto
C         TRAV retrouve sa valeur augmentee des termes de gradient transpose
          ELSE
            DO IEL = 1, NCEL
              TRAV(IEL,1)  = W7(IEL) + TRAV(IEL,1)
              TRAV(IEL,2)  = W8(IEL) + TRAV(IEL,2)
              TRAV(IEL,3)  = W9(IEL) + TRAV(IEL,3)
            ENDDO
          ENDIF
        ENDIF
C
      ENDIF
C
C
C ---> TERMES DE PERTES DE CHARGE
C     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT
C
      IF((NCEPDP.GT.0).AND.(IPHYDR.EQ.0)) THEN
C
C     Les termes diagonaux sont places dans TRAV ou TRAVA,
C       La prise en compte de UVWK a partir de la seconde iteration
C       est faite directement dans codits.
        IF(ITERNS.EQ.1) THEN
C
C       On utilise temporairement TRAV comme tableau de travail.
C       Son contenu est stocke dans W7, W8 et W9 jusqu'apres tsepdc
          DO IEL = 1,NCEL
            W7(IEL) = TRAV(IEL,1)
            W8(IEL) = TRAV(IEL,2)
            W9(IEL) = TRAV(IEL,3)
            TRAV(IEL,1) = 0.D0
            TRAV(IEL,2) = 0.D0
            TRAV(IEL,3) = 0.D0
          ENDDO
C
          IDIAEX = 1
          CALL TSEPDC
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDP , NCKPDP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  , IDIAEX ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , TRAV   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si on itere sur navsto, on utilise TRAVA ; sinon TRAV
          IF(NTERUP.GT.1) THEN
            DO IEL = 1, NCEL
              TRAVA(IEL,1,IPHAS) = TRAVA(IEL,1,IPHAS) + TRAV(IEL,1)
              TRAVA(IEL,2,IPHAS) = TRAVA(IEL,2,IPHAS) + TRAV(IEL,2)
              TRAVA(IEL,3,IPHAS) = TRAVA(IEL,3,IPHAS) + TRAV(IEL,3)
              TRAV(IEL,1)  = W7(IEL)
              TRAV(IEL,2)  = W8(IEL)
              TRAV(IEL,3)  = W9(IEL)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              TRAV(IEL,1)  = W7(IEL) + TRAV(IEL,1)
              TRAV(IEL,2)  = W8(IEL) + TRAV(IEL,2)
              TRAV(IEL,3)  = W9(IEL) + TRAV(IEL,3)
            ENDDO
          ENDIF
        ENDIF
C
C     Les termes extradiagonaux ne sont calcules qu'au premier passage
C       si on extrapole il va dans PROPCE,
C       sinon si on itere sur navsto dans TRAVA
C             sinon                  dans TRAV
C
        IF(ITERNS.EQ.1) THEN
C
C       On utilise temporairement TRAV comme tableau de travail.
C       Son contenu est stocke dans W7, W8 et W9 jusqu'apres tsepdc
          DO IEL = 1,NCEL
            W7(IEL) = TRAV(IEL,1)
            W8(IEL) = TRAV(IEL,2)
            W9(IEL) = TRAV(IEL,3)
            TRAV(IEL,1) = 0.D0
            TRAV(IEL,2) = 0.D0
            TRAV(IEL,3) = 0.D0
          ENDDO
C
          IDIAEX = 2
          CALL TSEPDC
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDP , NCKPDP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  , IDIAEX ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , TRAV   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Si on extrapole les termes source en temps :
C       PROPCE recoit les termes extradiagonaux et
C       TRAV retrouve sa valeur
          IF(ISNO2T(IPHAS).GT.0) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSNA  ) = PROPCE(IEL,IPTSNA  ) + TRAV(IEL,1)
              PROPCE(IEL,IPTSNA+1) = PROPCE(IEL,IPTSNA+1) + TRAV(IEL,2)
              PROPCE(IEL,IPTSNA+2) = PROPCE(IEL,IPTSNA+2) + TRAV(IEL,3)
              TRAV(IEL,1)  = W7(IEL)
              TRAV(IEL,2)  = W8(IEL)
              TRAV(IEL,3)  = W9(IEL)
            ENDDO
C
C     Si on n'extrapole pas les termes source en temps :
C       Si on itere sur navsto
C         TRAVA recoit les termes extradiagonaux et
C         TRAV retrouve sa valeur
          ELSE
            IF(NTERUP.GT.1) THEN
              DO IEL = 1, NCEL
                TRAVA(IEL,1,IPHAS) = TRAVA(IEL,1,IPHAS) + TRAV(IEL,1)
                TRAVA(IEL,2,IPHAS) = TRAVA(IEL,2,IPHAS) + TRAV(IEL,2)
                TRAVA(IEL,3,IPHAS) = TRAVA(IEL,3,IPHAS) + TRAV(IEL,3)
                TRAV(IEL,1)  = W7(IEL)
                TRAV(IEL,2)  = W8(IEL)
                TRAV(IEL,3)  = W9(IEL)
              ENDDO
C       Si on n'itere pas sur navsto
C         TRAV retrouve sa valeur augmentee des termes extradiagonaux
            ELSE
              DO IEL = 1, NCEL
                TRAV(IEL,1)  = W7(IEL) + TRAV(IEL,1)
                TRAV(IEL,2)  = W8(IEL) + TRAV(IEL,2)
                TRAV(IEL,3)  = W9(IEL) + TRAV(IEL,3)
              ENDDO
            ENDIF
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C
C ---> - DIVERGENCE DE RIJ
C
      IF(ITYTUR(IPHAS).EQ.3.AND.ITERNS.EQ.1) THEN
C
        DO ISOU = 1, 3
C
          IF(ISOU.EQ.1) IVAR = IUIPH
          IF(ISOU.EQ.2) IVAR = IVIPH
          IF(ISOU.EQ.3) IVAR = IWIPH
C
          CALL DIVRIJ
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISOU   , IVAR   , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   VISCF  , VISCB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
          INIT = 1
          CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                                   IFACEL,IFABOR,VISCF,VISCB,W1)
C
C     Si on extrapole les termes source en temps :
C       PROPCE recoit les termes de divergence
          IF(ISNO2T(IPHAS).GT.0) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSNA+ISOU-1 ) =
     &        PROPCE(IEL,IPTSNA+ISOU-1 ) - W1(IEL)
            ENDDO
C     Si on n'extrapole pas les termes source en temps :
          ELSE
C       si on n'itere pas sur navsto : TRAV
            IF(NTERUP.EQ.1) THEN
              DO IEL = 1, NCEL
                TRAV(IEL,ISOU) = TRAV(IEL,ISOU) - W1(IEL)
              ENDDO
C       si on itere sur navsto       : TRAVA
            ELSE
              DO IEL = 1, NCEL
                TRAVA(IEL,ISOU,IPHAS) = TRAVA(IEL,ISOU,IPHAS) - W1(IEL)
              ENDDO
            ENDIF
          ENDIF
C
        ENDDO
C
      ENDIF
C
C
C ---> "VITESSE" DE DIFFUSION FACETTE
C      SI ON FAIT AUTRE CHOSE QUE DU K EPS, IL FAUDRA LA METTRE
C        DANS LA BOUCLE
C
      IF( IDIFF(IUIPH).GE. 1 ) THEN
C
C --- Si la vitesse doit etre diffusee, on calcule la viscosite
C       pour le second membre (selon Rij ou non)
C
        IF (ITYTUR(IPHAS).EQ.3) THEN
          DO IEL = 1, NCEL
            W1(IEL) = PROPCE(IEL,IPCVIS)
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            W1(IEL) = PROPCE(IEL,IPCVIS)
     &                            + IDIFFT(IUIPH)*PROPCE(IEL,IPCVST)
          ENDDO
        ENDIF
C
        CALL VISCFA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
C
C     Quand on n'est pas en Rij ou que irijnu = 0, les tableaux
C       VISCFI, VISCBI se trouvent remplis par la meme occasion
C       (ils sont confondus avec VISCF, VISCB)
C     En Rij avec irijnu = 1, on calcule la viscosite increment
C       de la matrice dans VISCFI, VISCBI
C
        IF(ITYTUR(IPHAS).EQ.3.AND.IRIJNU(IPHAS).EQ.1) THEN
          DO IEL = 1, NCEL
            W1(IEL) = PROPCE(IEL,IPCVIS)
     &                            + IDIFFT(IUIPH)*PROPCE(IEL,IPCVST)
          ENDDO
          CALL VISCFA
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     ,
     &   VISCFI , VISCBI ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
C
      ELSE
C
C --- Si la vitesse n'a pas de diffusion, on annule la viscosite
C      (matrice et second membre sont dans le meme tableau,
C       sauf en Rij avec IRIJNU = 1)
C
        DO IFAC = 1, NFAC
          VISCF(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
C
        IF(ITYTUR(IPHAS).EQ.3.AND.IRIJNU(IPHAS).EQ.1) THEN
          DO IFAC = 1, NFAC
            VISCFI(IFAC) = 0.D0
          ENDDO
          DO IFAC = 1, NFABOR
            VISCBI(IFAC) = 0.D0
          ENDDO
        ENDIF
C
      ENDIF
C
C
C 2.2  RESOLUTION IMPLICITE NON COUPLEE DES 3 COMPO. DE VITESSES
C ==============================================================
C
C
C ---> AU PREMIER APPEL,
C      MISE A ZERO DE L'ESTIMATEUR POUR LA VITESSE PREDITE
C      S'IL DOIT ETRE CALCULE
C
      IF (IAPPEL.EQ.1) THEN
        IF(IESCAL(IESPRE,IPHAS).GT.0) THEN
          IESPRP = IPPROC(IESTIM(IESPRE,IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IESPRP) =  0.D0
          ENDDO
        ENDIF
      ENDIF
C
C ---> AU DEUXIEME APPEL,
C      MISE A ZERO DE L'ESTIMATEUR TOTAL POUR NAVIER-STOKES
C      (SI ON FAIT UN DEUXIEME APPEL, ALORS IL DOIT ETRE CALCULE)
C
      IF(IAPPEL.EQ.2) THEN
        IESTOP = IPPROC(IESTIM(IESTOT,IPHAS))
        DO IEL = 1, NCEL
          PROPCE(IEL,IESTOP) = 0.D0
        ENDDO
      ENDIF
C
C
C ---> BOUCLE SUR LES DIRECTIONS DE L'ESPACE (U, V, W)
C
C
C Remarque : On suppose que le couplage vitesse pression
C  n'est valable que pour une seule phase.
C
      DO ISOU = 1, 3
C
        IF(ISOU.EQ.1) THEN
          IVAR = IUIPH
          IPPT = IPPTX
        ENDIF
        IF(ISOU.EQ.2) THEN
          IVAR = IVIPH
          IPPT = IPPTY
        ENDIF
        IF(ISOU.EQ.3) THEN
          IVAR = IWIPH
          IPPT = IPPTZ
        ENDIF
        IPP  = IPPRTP(IVAR)
C
        ICLVAR = ICLRTP(IVAR,ICOEF)
        ICLVAF = ICLRTP(IVAR,ICOEFF)
C
C
C ---> TERMES SOURCES UTILISATEURS
C
        DO IEL = 1, NCEL
          W7    (IEL) = 0.D0
          DRTP  (IEL) = 0.D0
        ENDDO
C
C     Le calcul des parties implicite et explicite des termes sources
C       utilisateurs est faite uniquement a la premiere iter sur navsto.
        IF(ITERNS.EQ.1) THEN
C
          IF (IMATIS.EQ.1) THEN
C
            CALL MTTSNS
C           ===========
     & ( IDEBIA , IDEBRA ,
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
     &   W7     , DRTP   ,
C        ------   ------
     &   DAM    , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSE
C
            MAXELT = MAX(NCELET, NFAC, NFABOR)
            ILS    = IDEBIA
            IFINIA = ILS + MAXELT
            CALL IASIZE('PREDUV',IFINIA)
C
            CALL USTSNS
C           ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W7     , DRTP   ,
C        ------   ------
     &   DAM    , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
          ENDIF
C
        ENDIF
C
C     On conserve la partie implicite pour les autres iter sur navsto
        IF(ITERNS.EQ.1.AND.NTERUP.GT.1) THEN
          DO IEL = 1, NCEL
            XIMPA(IEL,ISOU,IPHAS) = DRTP(IEL)
          ENDDO
        ENDIF
C
C     On ajoute a TRAV ou TRAVA la partie issue des termes implicites
C       en utilisant DRTP
C       La prise en compte de UVWK a partir de la seconde iteration
C       est faite directement dans codits.
C     En schema std en temps, on continue a mettre MAX(-DRTP,0) dans la matrice
C     Avec termes sources a l'ordre 2, on implicite DRTP quel que soit son signe
C       (si on le met dans la matrice ou non selon son signe, on risque de ne pas
C        avoir le meme traitement d'un pas de temps au suivant)
        IF(ITERNS.EQ.1) THEN
          IF(NTERUP.GT.1) THEN
            DO IEL = 1, NCEL
              TRAVA(IEL,ISOU,IPHAS) = TRAVA(IEL,ISOU,IPHAS)
     &               + DRTP(IEL)*RTPA(IEL,IVAR)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              TRAV(IEL,ISOU) = TRAV(IEL,ISOU)
     &               + DRTP(IEL)*RTPA(IEL,IVAR)
            ENDDO
          ENDIF
        ENDIF
C
C     A la premiere iter sur navsto, on ajoute la partie issue des
C       termes explicites
        IF(ITERNS.EQ.1) THEN
C     Si on extrapole les termes source en temps :
C       PROPCE recoit les termes explicites
          IF(ISNO2T(IPHAS).GT.0) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSNA+ISOU-1 ) =
     &        PROPCE(IEL,IPTSNA+ISOU-1 ) + W7(IEL)
            ENDDO
C     Si on n'extrapole pas les termes source en temps :
          ELSE
C       si on n'itere pas sur navsto : TRAV
            IF(NTERUP.EQ.1) THEN
              DO IEL = 1, NCEL
                TRAV(IEL,ISOU) = TRAV(IEL,ISOU) + W7(IEL)
              ENDDO
C       si on itere sur navsto : TRAVA
            ELSE
              DO IEL = 1, NCEL
                TRAVA(IEL,ISOU,IPHAS) =
     &          TRAVA(IEL,ISOU,IPHAS) + W7(IEL)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
C
C
C ---> TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume
C
        INIT = 1
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                           IFACEL,IFABOR,FLUMAS,FLUMAB,W1)
C
C
C     On ajoute a TRAV ou TRAVA la partie issue des termes implicites
        IF(ITERNS.EQ.1) THEN
          IF(NTERUP.GT.1) THEN
            DO IEL = 1, NCEL
              TRAVA(IEL,ISOU,IPHAS) = TRAVA(IEL,ISOU,IPHAS)
     &               +ICONV(IVAR)*W1(IEL)*RTPA(IEL,IVAR)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              TRAV(IEL,ISOU) = TRAV(IEL,ISOU)
     &               +ICONV(IVAR)*W1(IEL)*RTPA(IEL,IVAR)
            ENDDO
          ENDIF
        ENDIF
C
        IF(IAPPEL.EQ.1) THEN
C     Extrapolation ou non, meme forme par coherence avec bilsc2
          DO IEL = 1, NCEL
            ROVSDT(IEL) =
     &         ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &         -ICONV(IVAR)*W1(IEL)*THETAV(IVAR)
          ENDDO
C
C     Le remplissage de ROVSDT est toujours indispensable,
C       meme si on peut se contenter de n'importe quoi pour IAPPEL=2.
        ELSE
          DO IEL = 1, NCEL
            ROVSDT(IEL) = 0.D0
          ENDDO
        ENDIF
C
C ---> TERMES SOURCES UTILISATEUR
C
        IF(IAPPEL.EQ.1) THEN
          IF(ISNO2T(IPHAS).GT.0) THEN
            THETAP = THETAV(IVAR)
            IF(ITERNS.GT.1) THEN
              DO IEL = 1, NCEL
                ROVSDT(IEL) = ROVSDT(IEL) -XIMPA(IEL,ISOU,IPHAS)*THETAP
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                ROVSDT(IEL) = ROVSDT(IEL) -DRTP(IEL)*THETAP
              ENDDO
            ENDIF
          ELSE
            IF(ITERNS.GT.1) THEN
              DO IEL = 1, NCEL
                ROVSDT(IEL) = ROVSDT(IEL)
     &               + MAX(-XIMPA(IEL,ISOU,IPHAS),ZERO)
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                ROVSDT(IEL) = ROVSDT(IEL)
     &               + MAX(-DRTP(IEL),ZERO)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
C
C
C ---> PERTES DE CHARGE
C
C  Au second appel, on n'a pas besoin de rovsdt
        IF(IAPPEL.EQ.1) THEN
          IF (NCEPDP.GT.0) THEN
            IF(ISNO2T(IPHAS).GT.0) THEN
              THETAP = THETAV(IVAR)
            ELSE
              THETAP = 1.D0
            ENDIF
            DO IELPDC = 1, NCEPDP
              IEL = ICEPDC(IELPDC)
              ROVSDT(IEL) = ROVSDT(IEL) +
     &        PROPCE(IEL,IPCROM)*VOLUME(IEL)*CKUPDC(IELPDC,ISOU)*THETAP
            ENDDO
          ENDIF
        ENDIF
C
C
C --->  TERMES DE SOURCE DE MASSE
C
        IF (NCESMP.GT.0) THEN
C
C     On calcule les termes Gamma (uinj - u)
C       -Gamma u a la premiere iteration est mis dans
C          TRAV ou TRAVA selon qu'on itere ou non sur navsto
C       Gamma uinj a la premiere iteration est placee dans W1
C       ROVSDT a chaque iteration recoit Gamma
          IF(NTERUP.EQ.1) THEN
            CALL CATSMA
C           ===========
     &  ( NCELET , NCEL , NCESMP , ITERNS , ISNO2T(IPHAS), THETAV(IVAR),
     &    ICETSM , ITYPSM(1,IVAR) ,
     &    VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) ,SMACEL(1,IPR(IPHAS)) ,
     &    TRAV(1,ISOU)        , ROVSDT , W1 )
          ELSE
            CALL CATSMA
C           ===========
     &  ( NCELET , NCEL , NCESMP , ITERNS , ISNO2T(IPHAS), THETAV(IVAR),
     &    ICETSM , ITYPSM(1,IVAR) ,
     &    VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) ,SMACEL(1,IPR(IPHAS)) ,
     &    TRAVA(1,ISOU,IPHAS) , ROVSDT , W1 )
          ENDIF
C
C     A la premiere iter sur navsto, on ajoute la partie Gamma uinj
          IF(ITERNS.EQ.1) THEN
C     Si on extrapole les termes source en temps :
C       PROPCE recoit les termes explicites
            IF(ISNO2T(IPHAS).GT.0) THEN
              DO IEL = 1,NCEL
                PROPCE(IEL,IPTSNA+ISOU-1 ) =
     &          PROPCE(IEL,IPTSNA+ISOU-1 ) + W1(IEL)
              ENDDO
C     Si on n'extrapole pas les termes source en temps :
            ELSE
C       si on n'itere pas sur navsto : TRAV
              IF(NTERUP.EQ.1) THEN
                DO IEL = 1,NCEL
                  TRAV(IEL,ISOU)  = TRAV(IEL,ISOU) + W1(IEL)
                ENDDO
C       si on itere sur navsto : TRAVA
              ELSE
                DO IEL = 1,NCEL
                  TRAVA(IEL,ISOU,IPHAS) =
     &            TRAVA(IEL,ISOU,IPHAS) + W1(IEL)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
C
        ENDIF
C
C
C ---> INITIALISATION DU SECOND MEMBRE
C
C     Si on extrapole les TS
        IF(ISNO2T(IPHAS).GT.0) THEN
          THETP1 = 1.D0 + THETS
C       Si on n'itere pas sur navsto : TRAVA n'existe pas
          IF(NTERUP.EQ.1) THEN
            DO IEL = 1, NCEL
              SMBR(IEL) =  TRAV(IEL,ISOU)
     &             + THETP1*PROPCE(IEL,IPTSNA+ISOU-1)
            ENDDO
C       Si on   itere     sur navsto : tout existe
          ELSE
            DO IEL = 1, NCEL
              SMBR(IEL) =  TRAV(IEL,ISOU) + TRAVA(IEL,ISOU,IPHAS)
     &             + THETP1*PROPCE(IEL,IPTSNA+ISOU-1)
            ENDDO
          ENDIF
C     Si on n'extrapole pas les TS : PROPCE n'existe pas
        ELSE
C       Si on n'itere pas sur navsto : TRAVA n'existe pas
          IF(NTERUP.EQ.1) THEN
            DO IEL = 1, NCEL
              SMBR(IEL) =  TRAV(IEL,ISOU)
            ENDDO
C       Si on   itere     sur navsto : TRAVA existe
          ELSE
            DO IEL = 1, NCEL
              SMBR(IEL) =  TRAV(IEL,ISOU) + TRAVA(IEL,ISOU,IPHAS)
            ENDDO
          ENDIF
        ENDIF
C
C
C ---> LAGRANGIEN : COUPLAGE RETOUR
C
C     L'ordre 2 sur les termes issus du lagrangien necessiterait de
C       decomposer TSLAGR(IEL,ISOU) en partie implicite et
C       explicite, comme c'est fait dans ustsns.
C     Pour le moment, on n'y touche pas.
        IF (IILAGR.EQ.2 .AND. LTSDYN.EQ.1 .AND. IPHAS.EQ.ILPHAS)  THEN
C
          DO IEL = 1, NCEL
            SMBR(IEL)   = SMBR(IEL) + TSLAGR(IEL,ITSVX+ISOU-1)
          ENDDO
C  Au second appel, on n'a pas besoin de rovsdt
          IF(IAPPEL.EQ.1) THEN
            DO IEL = 1, NCEL
              ROVSDT(IEL) = ROVSDT(IEL) + MAX(-TSLAGR(IEL,ITSLI),ZERO)
            ENDDO
          ENDIF
C
        ENDIF
C
C ---> VERSIONS ELECTRIQUES : Arc Electrique (Force de Laplace)
C     Pour le moment, pas d'ordre 2 en temps.
C
        IF ( IPPMOD(IELARC) .GE. 1 ) THEN
          DO IEL = 1,NCEL
            SMBR(IEL) = SMBR(IEL)
     &                 + VOLUME(IEL)*PROPCE(IEL,IPPROC(ILAPLA(ISOU)))
          ENDDO
        ENDIF
C
C
C ---> PARAMETRES POUR LA RESOLUTION DU SYSTEME OU LE CALCUL DE l'ESTIMATEUR
C
        ICONVP = ICONV (IVAR)
        IDIFFP = IDIFF (IVAR)
        IRESLP = IRESOL(IVAR)
        NDIRCP = NDIRCL(IVAR)
        NITMAP = NITMAX(IVAR)
CMO        IMRGRA
        NSWRSP = NSWRSM(IVAR)
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IRCFLP = IRCFLU(IVAR)
        ISCHCP = ISCHCV(IVAR)
        ISSTPP = ISSTPC(IVAR)
        IMGRP  = IMGR  (IVAR)
        NCYMXP = NCYMAX(IVAR)
        NITMFP = NITMGF(IVAR)
CMO        IPP
        IWARNP = IWARNI(IVAR)
        BLENCP = BLENCV(IVAR)
        EPSILP = EPSILO(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
        RELAXP = RELAXV(IVAR)
        THETAP = THETAV(IVAR)
C
        IF(IAPPEL.EQ.1) THEN
C
          IESCAP = IESCAL(IESPRE,IPHAS)
C
C ---> FIN DE LA CONSTRUCTION ET DE LA RESOLUTION DU SYSTEME
C
          IF(ITERNS.EQ.1) THEN
C
C  Attention, dans le cas des estimateurs, DAM fournit l'estimateur
C     des vitesses predites
            CALL CODITS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     FLUMAS , FLUMAB ,
     &   VISCFI , VISCBI , VISCF  , VISCB  ,
     &   ROVSDT , SMBR   , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSEIF(ITERNS.GT.1) THEN
C
            CALL CODITS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , RELAXP, THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , UVWK(1,ISOU,IPHAS) ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     FLUMAS , FLUMAB ,
     &   VISCFI , VISCBI , VISCF  , VISCB  ,
     &   ROVSDT , SMBR   , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
C     DANS LE CAS DE PERTES DE CHARGE, ON UTILISE LES TABLEAUX
C     TPUCOU POUR LA PHASE D'IMPLICITATION
C     Attention, il faut regarder s'il y a des pdc sur un proc quelconque,
C       pas uniquement en local.
          IF((NCPDCT(IPHAS).GT.0).AND.(IPUCOU.EQ.0)) THEN
            DO IEL = 1,NCEL
              TPUCOU(IEL,ISOU) = DT(IEL)
            ENDDO
            DO IELPDC = 1, NCEPDP
              IEL=ICEPDC(IELPDC)
              TPUCOU(IEL,ISOU) = 1.D0/(
     &             1.D0/DT(IEL)+CKUPDC(IELPDC,ISOU))
            ENDDO
          ENDIF
C
C     COUPLAGE P/U RENFORCE : CALCUL DU VECTEUR T, STOCKE DANS TPUCOU
C       ON PASSE DANS CODITS, EN NE FAISANT QU'UN SEUL SWEEP, ET EN
C       INITIALISANT TPUCOU A 0 POUR QUE LA PARTIE CV/DIFF AJOUTEE
C       PAR BILSC2 SOIT NULLE
C     NSWRSP = -1 INDIQUERA A CODITS QU'IL NE FAUT FAIRE QU'UN SWEEP
C       ET QU'IL FAUT METTRE INC A 0 (POUR OTER LES DIRICHLETS DANS LES
C       CL DES MATRICES POIDS)
C
          IF (IPUCOU.EQ.1) THEN
            NSWRSP = -1
            DO IEL = 1, NCEL
              SMBR(IEL) = VOLUME(IEL)
            ENDDO
            DO IEL = 1, NCELET
              TPUCOU(IEL,ISOU) = 0.D0
            ENDDO
            IESCAP = 0
C
            CALL CODITS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGRP  , NCYMXP , NITMFP , IPPT   , IWARNP ,
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   TPUCOU(1,ISOU)  , TPUCOU(1,ISOU)  ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     FLUMAS , FLUMAB ,
     &   VISCFI , VISCBI , VISCF  , VISCB  ,
     &   ROVSDT , SMBR   , TPUCOU(1,ISOU)  ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
            DO IEL = 1, NCELET
              TPUCOU(IEL,ISOU) = PROPCE(IEL,IPCROM)*TPUCOU(IEL,ISOU)
            ENDDO
C
          ENDIF
C
C
C --->  ESTIMATEUR SUR LA VITESSE PREDITE : ON SOMME SUR LES COMPOSANTES
C
          IF(IESCAL(IESPRE,IPHAS).GT.0) THEN
            IESPRP = IPPROC(IESTIM(IESPRE,IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IESPRP) =  PROPCE(IEL,IESPRP) + DAM(IEL)
            ENDDO
          ENDIF
C
        ELSEIF(IAPPEL.EQ.2) THEN
C
C ---> FIN DE LA CONSTRUCTION DE L'ESTIMATEUR
C        RESIDU SECOND MEMBRE(Un+1,Pn+1) + RHO*VOLUME*( Un+1 - Un )/DT
C
          INC = 1
          ICCOCG = 1
C     Pas de relaxation en stationnaire
          IDTVA0 = 0
C
          CALL BILSC2
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVA0 , IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP(1,IVAR)     , RTP(1,IVAR)     ,
     &   COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &   FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          IESTOP = IPPROC(IESTIM(IESTOT,IPHAS))
          DO IEL = 1, NCEL
            PROPCE(IEL,IESTOP) =
     &           PROPCE(IEL,IESTOP)+ (SMBR(IEL)/VOLUME(IEL))**2
          ENDDO
C
        ENDIF
C
      ENDDO
C
C
C --->  APRES LA BOUCLE SUR U, V, W,
C        FIN DU CALCUL DE LA NORME POUR RESOLP
C
      IF(IAPPEL.EQ.1.AND.IRNPNW.EQ.1) THEN
C
C     Calcul de div(rho u*)
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTP(1,IUIPH))
C         ===========
          CALL PARCOM (RTP(1,IVIPH))
C         ===========
          CALL PARCOM (RTP(1,IWIPH))
C         ===========
        ENDIF
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 1
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    RTP(1,IUIPH),RTP(1,IUIPH),RTP(1,IUIPH),
     &    RTP(1,IVIPH),RTP(1,IVIPH),RTP(1,IVIPH),
     &    RTP(1,IWIPH),RTP(1,IWIPH),RTP(1,IWIPH))
        ENDIF
C
C       Pour gagner du temps, on ne reconstruit pas.
        INIT   = 1
        INC    = 1
        ICCOCG = 1
        IFLMB0 = 1
        NSWRP  = 1
        IMLIGP = IMLIGR(IUIPH )
        IWARNP = IWARNI(IPRIPH)
        EPSRGP = EPSRGR(IUIPH )
        CLIMGP = CLIMGR(IUIPH )
        EXTRAP = EXTRAG(IUIPH )
C
        IMASPE = 1
C
        IISMPH = IISYMP+NFABOR*(IPHAS-1)
C
        CALL INIMAS
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRP  , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   RTP(1,IUIPH) , RTP(1,IVIPH) , RTP(1,IWIPH) ,
     &   COEFA(1,ICLIUP), COEFA(1,ICLIVP), COEFA(1,ICLIWP),
     &   COEFB(1,ICLIUP), COEFB(1,ICLIVP), COEFB(1,ICLIWP),
     &   VISCF  , VISCB  ,
     &   W4     , W5     , W6     , W7     , W8     , W9     ,
     &   SMBR   , DRTP   , ROVSDT , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
        INIT = 0
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                                IFACEL,IFABOR,VISCF,VISCB,XNORMP)
C
C     Calcul de la norme
C       RNORMP qui servira dans resolp
        ISQRT = 1
        CALL PRODSC(NCELET,NCEL,ISQRT,XNORMP,XNORMP,RNORMP(IPHAS))
C
      ENDIF
C
C --->  APRES LA BOUCLE SUR U, V, W,
C        FIN DU CALCUL DES ESTIMATEURS ET IMPRESSION
C
      IF(IAPPEL.EQ.1) THEN
C
C --->  ESTIMATEUR SUR LA VITESSE PREDITE : ON PREND LA RACINE (NORME)
C         SANS OU AVEC VOLUME (ET DANS CE CAS C'EST LA NORME L2)
C
        IF(IESCAL(IESPRE,IPHAS).GT.0) THEN
          IESPRP = IPPROC(IESTIM(IESPRE,IPHAS))
          IF(IESCAL(IESPRE,IPHAS).EQ.1) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IESPRP) =  SQRT(PROPCE(IEL,IESPRP)            )
            ENDDO
          ELSEIF(IESCAL(IESPRE,IPHAS).EQ.2) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IESPRP) =  SQRT(PROPCE(IEL,IESPRP)*VOLUME(IEL))
            ENDDO
          ENDIF
        ENDIF
C
C ---> IMPRESSION DE NORME
C
        IF (IWARNI(IUIPH).GE.2) THEN
          RNORM = -1.D0
          DO IEL = 1, NCEL
            VITNOR =
     &       SQRT(RTP(IEL,IUIPH)**2+RTP(IEL,IVIPH)**2+RTP(IEL,IWIPH)**2)
            RNORM = MAX(RNORM,VITNOR)
          ENDDO
          IF (IRANGP.GE.0) CALL PARMAX (RNORM)
C                          ===========
          WRITE(NFECRA,1100) IPHAS, RNORM
        ENDIF
C
      ELSEIF (IAPPEL.EQ.2) THEN
C
C --->  ESTIMATEUR SUR NAVIER-STOKES TOTAL : ON PREND LA RACINE (NORME)
C         SANS OU AVEC VOLUME (ET DANS CE CAS C'EST LA NORME L2)
C
        IESTOP = IPPROC(IESTIM(IESTOT,IPHAS))
        IF(IESCAL(IESTOT,IPHAS).EQ.1) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IESTOP) = SQRT(PROPCE(IEL,IESTOP)            )
          ENDDO
        ELSEIF(IESCAL(IESTOT,IPHAS).EQ.2) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IESTOP) = SQRT(PROPCE(IEL,IESTOP)*VOLUME(IEL))
          ENDDO
        ENDIF
C
      ENDIF
C
C--------
C FORMATS
C--------
#if defined(_CS_LANG_FR)
C
 1100 FORMAT(/,
     & 1X,'Phase ',I4,' : Vitesse maximale apres prediction ',E12.4)
C
#else
C
 1100 FORMAT(/,
     & 1X,'Phase ',I4,' : Maximum velocity after prediction ',E12.4)
C
#endif
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
