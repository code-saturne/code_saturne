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
                        SUBROUTINE TURBKW
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR , COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   S2KW   , DIVUKW , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRK  , SMBRW  , TINSTK , TINSTW ,
     &   XF1    , W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS K-OMEGA SST 1 PHASE INCOMPRESSIBLE OU
CFONC  RHO VARIABLE SUR UN PAS DE TEMPS
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
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,*)    !    !     !     LAGRANGIEN                       !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! S2KW(NCELET) ! TR !  -> ! TABLEAU CONTENANT S2=2.SijSij        !
CARGU ! DIVUKW(NCELET! TR !  -> ! DIVERGENCE DE U (CALCULEE PAR GRDCEL)!
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR.(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! TINST.(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! XF1(NCELET)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR COEF F1      !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "optcal.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "pointe.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
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
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION S2KW(NCELET), DIVUKW(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET),   SMBRK(NCELET),  SMBRW(NCELET)
      DOUBLE PRECISION TINSTK(NCELET), TINSTW(NCELET), XF1(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IEL   , IFAC  , INIT  , INC   , ICCOCG, IVAR
      INTEGER          II, IIVAR , IIUN  , IFACPT
      INTEGER          ICLIPK, ICLIPW, ISQRT
      INTEGER          NSWRGP, IMLIGP, IPHYDP
      INTEGER          IPRIPH, IKIPH , IOMGIP
      INTEGER          ICLIKP, ICLOMG
      INTEGER          ICLVAR, ICLVAF
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP
      INTEGER          NITMAP, NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IPCROM, IPBROM, IPCVST, IPCVIS, IFLMAS, IFLMAB
      INTEGER          IWARNP, IPP
      INTEGER          IPTSTA
      INTEGER          IPCROO, IPBROO, IPCVTO, IPCVLO
      INTEGER          MAXELT, ILS
      DOUBLE PRECISION RNORM , D2S3, DIVP23, EPZ2
      DOUBLE PRECISION DELTK , DELTW, A11, A12, A22, A21
      DOUBLE PRECISION UNSDET, ROMVSD
      DOUBLE PRECISION PRDTUR, XK, XW, XEPS, XNU
      DOUBLE PRECISION VISCT , ROM
      DOUBLE PRECISION BLENCP, EPSILP, EPSRGP, CLIMGP, EXTRAP, RELAXP
      DOUBLE PRECISION THETP1, THETAK, THETAW, THETS, THETAP, EPSRSP
      DOUBLE PRECISION TUEXPK, TUEXPW
      DOUBLE PRECISION CDKW, XARG1, XXF1, XGAMMA, XBETA, SIGMA, PRODUC
      DOUBLE PRECISION VAR, VRMIN, VRMAX
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      EPZ2 = EPZERO**2
C
      IPRIPH = IPR (IPHAS)
      IKIPH  = IK  (IPHAS)
      IOMGIP = IOMG(IPHAS)
C
      ICLIKP = ICLRTP(IKIPH ,ICOEF)
      ICLOMG = ICLRTP(IOMGIP,ICOEF)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IKIPH))
      IFLMAB = IPPROB(IFLUMA(IKIPH))
      IPBROM = IPPROB(IROM  (IPHAS))
C
      THETS  = THETST(IPHAS)
C
      IPCROO = IPCROM
      IPBROO = IPBROM
      IPCVTO = IPCVST
      IPCVLO = IPCVIS
      IF(ISTO2T(IPHAS).GT.0) THEN
        IF (IROEXT(IPHAS).GT.0) THEN
          IPCROO = IPPROC(IROMA(IPHAS))
          IPBROO = IPPROB(IROMA(IPHAS))
        ENDIF
        IF(IVIEXT(IPHAS).GT.0) THEN
          IPCVTO = IPPROC(IVISTA(IPHAS))
          IPCVLO = IPPROC(IVISLA(IPHAS))
        ENDIF
      ENDIF
C
      IF(ISTO2T(IPHAS).GT.0) THEN
        IPTSTA = IPPROC(ITSTUA(IPHAS))
      ELSE
        IPTSTA = 0
      ENDIF
C
      IF(IWARNI(IKIPH).GE.1) THEN
        WRITE(NFECRA,1000)IPHAS
      ENDIF
C
C
C=======================================================================
C 2. CALCUL DE dk/dxj.dw/dxj
C      Le terme est stocke dans         W1
C      En sortie de l'etape on conserve W1
C=======================================================================
C
      ICCOCG = 1
      INC = 1
C
C
      NSWRGP = NSWRGR(IKIPH)
      IMLIGP = IMLIGR(IKIPH)
      IWARNP = IWARNI(IKIPH)
      EPSRGP = EPSRGR(IKIPH)
      CLIMGP = CLIMGR(IKIPH)
      EXTRAP = EXTRAG(IKIPH)
      IPHYDP = 0
C
      CALL GRDCEL
C     ===========
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
     &   W7     , W7     , W7     ,
     &   RTPA(1,IKIPH)   , COEFA(1,ICLIKP) , COEFB(1,ICLIKP) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W7     , W8     , XF1    ,
     &   RDEVEL , RTUSER , RA     )
C
C
      NSWRGP = NSWRGR(IOMGIP)
      IMLIGP = IMLIGR(IOMGIP)
      IWARNP = IWARNI(IOMGIP)
      EPSRGP = EPSRGR(IOMGIP)
      CLIMGP = CLIMGR(IOMGIP)
      EXTRAP = EXTRAG(IOMGIP)
      IPHYDP = 0
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IOMGIP , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W7     , W7     , W7     ,
     &   RTPA(1,IOMGIP)  , COEFA(1,ICLOMG) , COEFB(1,ICLOMG) ,
     &   W4     , W5     , W6     ,
C        ------   ------   ------
     &   W7     , W8     , XF1    ,
     &   RDEVEL , RTUSER , RA     )
C
      DO IEL = 1, NCEL
        W1(IEL) = W1(IEL)*W4(IEL)+W2(IEL)*W5(IEL)+W3(IEL)*W6(IEL)
      ENDDO
C
C============================================
C 3. CALCUL DU COEFFICIENT DE PONDERATION F1
C      Le terme est stocke dans         XF1
C      En sortie de l'etape on conserve W1,XF1
C============================================
C
C
      IF(ABS(ICDPAR).EQ.2) THEN
        DO IEL = 1, NCEL
          IFACPT = IA(IIFAPA(IPHAS)-1+IEL)
          W2(IEL) =
     &         (CDGFBO(1,IFACPT)-XYZCEN(1,IEL))**2
     &         +(CDGFBO(2,IFACPT)-XYZCEN(2,IEL))**2
     &         +(CDGFBO(3,IFACPT)-XYZCEN(3,IEL))**2
          W2(IEL) = SQRT(W2(IEL))
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          W2(IEL) =  MAX(RA(IDIPAR+IEL-1),EPZERO)
        ENDDO
      ENDIF
C
C     En cas d'ordre 2 on utilise les valeurs en n car le terme en (1-F1)*W1
C     sera dans PROPCE. Du coup, on aura quand meme certaines "constantes"
C     intervenant dans des termes en n+1/2 (ex sigma_k pour la diffusion) calcules
C     a partir de F1 en n -> mais l'effet sur les "constantes" est faible
C     -> a garder en tete si on fait vraiment de l'ordre 2 en temps en k-omega
      DO IEL = 1, NCEL
        ROM = PROPCE(IEL,IPCROO)
        XNU = PROPCE(IEL,IPCVLO)/ROM
        XK = RTPA(IEL,IKIPH)
        XW  = RTPA(IEL,IOMGIP)
        CDKW = 2*ROM/CKWSW2/XW*W1(IEL)
        CDKW = MAX(CDKW,1.D-20)
        XARG1 = MAX( SQRT(XK)/CMU/XW/W2(IEL),
     &       500.D0*XNU/XW/W2(IEL)**2 )
        XARG1 = MIN(XARG1,
     &       4.D0*ROM*XK/CKWSW2/CDKW/W2(IEL)**2)
        XF1(IEL) = TANH(XARG1**4)
      ENDDO
C
C=======================================================================
C 4. CALCUL DU TERME DE PRODUCTION
C      Les termes sont stockes dans     TINSTK,TINSTW
C      En sortie de l'etape on conserve W1,XF1,TINSTK,TINSTW
C=======================================================================
C
      D2S3 = 2.D0/3.D0
      DO IEL = 1, NCEL
        XK   = RTPA(IEL,IKIPH)
        XEPS = CMU*RTPA(IEL,IOMGIP)*XK
        TINSTK(IEL) = S2KW(IEL) - D2S3*DIVUKW(IEL)*DIVUKW(IEL)
        TINSTW(IEL) = PROPCE(IEL,IPCVTO)*TINSTK(IEL)
     &       -D2S3*PROPCE(IEL,IPCROO)*XK*DIVUKW(IEL)
        TINSTK(IEL) = MIN(TINSTW(IEL),CKWC1*PROPCE(IEL,IPCROO)*XEPS)
      ENDDO
C
C=======================================================================
C 4. CALCUL DU TERME DE GRAVITE
C      Les termes sont stockes dans     TINSTK,TINSTW,W2
C      En sortie de l'etape on conserve W1,W2,XF1,TINSTK,TINSTW
C=======================================================================
C
      IF(IGRAKE(IPHAS).EQ.1) THEN
C
C --- Terme de gravite G = BETA*G*GRAD(SCA)/PRDTUR/RHO
C     Ici on calcule   G =-G*GRAD(RHO)/PRDTUR/RHO
C
        ICCOCG = 1
        INC = 1
C
C     Le choix ci dessous a l'avantage d'etre simple
C
        NSWRGP = NSWRGR(IKIPH)
        EPSRGP = EPSRGR(IKIPH)
        IMLIGP = IMLIGR(IKIPH)
        IWARNP = IWARNI(IKIPH)
        CLIMGP = CLIMGR(IKIPH)
        EXTRAP = EXTRAG(IKIPH)
C
C     Conditions aux limites sur ROM : Dirichlet ROMB
C       On utilise VISCB pour stocker le COEFB relatif a ROM
C       On impose en Dirichlet (COEFA) la valeur ROMB
C
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
C
        IIVAR = 0
C
        IPHYDP = 0
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IIVAR  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W5     , W5     , W5     ,
     &   PROPCE(1,IPCROO), PROPFB(1,IPBROO), VISCB  ,
     &   W2     , W3     , W4     ,
C        ------   ------   ------
     &   W5     , W6     , W7     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C      Production et terme de gravite
C        TINSTK=MIN(P,C1*EPS)+G et TINSTW=P+(1-CE3)*G
C        On conserve G dans W2 pour la phase de couplage des termes sources
C
        IF(ISCALT(IPHAS).GT.0.AND.NSCAL.GE.ISCALT(IPHAS)) THEN
          PRDTUR = SIGMAS(ISCALT(IPHAS))
        ELSE
          PRDTUR = 1.D0
        ENDIF
C
        DO IEL = 1, NCEL
          W2(IEL) = -(W2(IEL)*GX+W3(IEL)*GY+W4(IEL)*GZ)/
     &         (PROPCE(IEL,IPCROO)*PRDTUR)
          TINSTW(IEL)=TINSTW(IEL)+PROPCE(IEL,IPCVTO)*MAX(W2(IEL),ZERO)
          TINSTK(IEL)=TINSTK(IEL)+PROPCE(IEL,IPCVTO)*W2(IEL)
        ENDDO
C
      ENDIF
C
C
C======================================================================
C 5. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
C      Les termes sont stockes dans     SMBRK,SMBRW,DAM,W3
C      En sortie de l'etape on conserve W1-3,XF1,TINSTK,TINSTW,SMBRK,SMBRW,DAM
C=======================================================================
C
      DO IEL = 1, NCEL
        SMBRK(IEL) = 0.D0
        SMBRW(IEL) = 0.D0
        DAM  (IEL) = 0.D0
        W3   (IEL) = 0.D0
      ENDDO
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('TURBKW',IFINIA)
C
      CALL USTSKW
C     ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL , S2KW   , DIVUKW ,
     &   W1     , W2     , XF1    ,
     &   SMBRK  , SMBRW  , DAM    , W3     ,
C        ------   ------   ------   ------
     &   VISCF  , VISCB  , XAM    , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , DRTP   ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 6. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume
C
C      Le terme est stocke dans         W4
C      En sortie de l'etape on conserve W1-4,XF1,TINSTK,TINSTW,SMBRK,SMBRW,DAM
C=======================================================================
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W4)
C
C=======================================================================
C 7. ON FINALISE LE CALCUL DES TERMES SOURCES
C
C      Les termes sont stockes dans     SMBRK, SMBRW
C      En sortie de l'etape on conserve SMBRK,SMBRW,DAM,W1-4
C=======================================================================
C
      DO IEL = 1, NCEL
C
        VISCT  = PROPCE(IEL,IPCVTO)
        ROM    = PROPCE(IEL,IPCROO)
        XK     = RTPA(IEL,IKIPH)
        XW     = RTPA(IEL,IOMGIP)
        XXF1   = XF1(IEL)
        XGAMMA = XXF1*CKWGM1 + (1.D0-XXF1)*CKWGM2
        XBETA  = XXF1*CKWBT1 + (1.D0-XXF1)*CKWBT2
C
        SMBRK(IEL) = SMBRK(IEL) + VOLUME(IEL)*(
     &       TINSTK(IEL)
     &       -CMU*ROM*XW*XK )
C
        SMBRW(IEL) = SMBRW(IEL) + VOLUME(IEL)*(
     &       ROM*XGAMMA/VISCT*TINSTW(IEL)
     &       -XBETA*ROM*XW**2
     &       +2.D0*ROM/XW*(1.D0-XXF1)/CKWSW2*W1(IEL) )
C
      ENDDO
C
C======================================================================
C 8. PRISE EN COMPTE DES TERMES D'ACCUMULATION DE MASSE ET
C         DE LA DEUXIEME PARTIE DES TS UTILISATEURS (PARTIE EXPLICITE)
C         STOCKAGE POUR EXTRAPOLATION EN TEMPS
C      On utilise                       SMBRK,SMBRW
C      En sortie de l'etape on conserve SMBRK,SMBRW,DAM,W2-4
C
C    Remarque : l'extrapolation telle qu'elle est ecrite n'a pas grand
C               sens si IKECOU=1
C=======================================================================
C
C     Si on extrapole les T.S.
      IF(ISTO2T(IPHAS).GT.0) THEN
C
        DO IEL = 1, NCEL
C
C       Sauvegarde pour echange
          TUEXPK = PROPCE(IEL,IPTSTA)
C       Pour la suite et le pas de temps suivant
          PROPCE(IEL,IPTSTA) = SMBRK(IEL)
C       Termes dependant de la variable resolue et theta PROPCE
          SMBRK(IEL) = ICONV(IKIPH)*W4(IEL)*RTPA(IEL,IKIPH)
     &         - THETS*TUEXPK
C       On suppose -DAM > 0 : on implicite
C         le terme utilisateur dependant de la variable resolue
          SMBRK(IEL) = DAM(IEL)*RTPA(IEL,IKIPH) + SMBRK(IEL)
C
C       Sauvegarde pour echange
          TUEXPW = PROPCE(IEL,IPTSTA+1)
C       Pour la suite et le pas de temps suivant
          PROPCE(IEL,IPTSTA+1) = SMBRW(IEL)
C       Termes dependant de la variable resolue et theta PROPCE
          SMBRW(IEL) = ICONV(IOMGIP)*W4(IEL)*RTPA(IEL,IOMGIP)
     &         - THETS*TUEXPW
C       On suppose -W3 > 0 : on implicite
C         le terme utilisateur dependant de la variable resolue
          SMBRW(IEL) =  W3(IEL)*RTPA(IEL,IOMGIP) + SMBRW(IEL)
C
        ENDDO
C
C     Si on n'extrapole pas les T.S.
      ELSE
        DO IEL = 1, NCEL
          SMBRK(IEL) = SMBRK(IEL) + DAM(IEL)*RTPA(IEL,IKIPH)
     &         +ICONV(IKIPH)*W4(IEL)*RTPA(IEL,IKIPH)
          SMBRW(IEL) = SMBRW(IEL) + W3 (IEL)*RTPA(IEL,IOMGIP)
     &         +ICONV(IOMGIP)*W4(IEL)*RTPA(IEL,IOMGIP)
        ENDDO
      ENDIF
C
C======================================================================
C 8.1 PRISE EN COMPTE DES TERMES SOURCES LAGRANGIEN : PARTIE EXPLICITE
C     COUPLAGE RETOUR
C=======================================================================
C
C     Ordre 2 non pris en compte
      IF (IILAGR.EQ.2 .AND. LTSDYN.EQ.1 .AND. IPHAS.EQ.ILPHAS) THEN
C
        DO IEL = 1,NCEL
C
C Termes sources explicte et implicte sur k
C
          SMBRK(IEL)  = SMBRK(IEL) + TSLAGR(IEL,ITSKE)
C
C Termes sources explicte sur omega : on reprend la constante CE4 directement
C    du k-eps sans justification ... a creuser si necessaire !
C
          SMBRW(IEL)  = SMBRW(IEL)
     &                + CE4 *TSLAGR(IEL,ITSKE) * PROPCE(IEL,IPCROO)
     &                /PROPCE(IEL,IPCVTO)
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 9. PRISE EN COMPTE DES TERMES DE CONV/DIFF DANS LE SECOND MEMBRE
C
C      Tableaux de travail              W7, W8, W1, TINSTK, TINSTW, DRTP
C      Les termes sont stockes dans     W5 et W6, puis ajoutes a SMBRK, SMBRW
C      En sortie de l'etape on conserve W2-6,SMBRK,SMBRW,DAM
C=======================================================================
C
C     Ceci ne sert a rien si IKECOU n'est pas egal a 1
C
      IF (IKECOU(IPHAS).EQ.1) THEN
C
        DO IEL = 1, NCEL
          W5 (IEL) = 0.D0
          W6 (IEL) = 0.D0
        ENDDO
C
C ---> Traitement de k
C
        IVAR   = IKIPH
C
        IPP    = IPPRTP(IVAR)
C
        ICLVAR = ICLRTP(IVAR,ICOEF )
        ICLVAF = ICLRTP(IVAR,ICOEFF)
        CHAINE = NOMVAR(IPP)
C
        IF( IDIFF(IVAR).GE. 1 ) THEN
C
          DO IEL = 1, NCEL
            XXF1 = XF1(IEL)
            SIGMA = XXF1*CKWSK1 + (1.D0-XXF1)*CKWSK2
            W7(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMA
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
     &   W7     ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )

        ELSE
C
          DO IFAC = 1, NFAC
            VISCF(IFAC) = 0.D0
          ENDDO
          DO IFAC = 1, NFABOR
            VISCB(IFAC) = 0.D0
          ENDDO
C
        ENDIF
C
        ICCOCG = 1
        INC    = 1
        ICONVP = ICONV (IVAR)
        IDIFFP = IDIFF (IVAR)
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IRCFLP = IRCFLU(IVAR)
        ISCHCP = ISCHCV(IVAR)
        ISSTPP = ISSTPC(IVAR)
        IWARNP = IWARNI(IVAR)
        BLENCP = BLENCV(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
        RELAXP = RELAXV(IVAR)
        THETAP = THETAV(IVAR)
C
        CALL BILSC2
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &   COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   W5  ,
C        --
     &   W7     , W8     , W1     , TINSTK , TINSTW , DRTP ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (IWARNI(IVAR).GE.2) THEN
          ISQRT = 1
          CALL PRODSC(NCELET,NCEL,ISQRT,SMBRK,SMBRK,RNORM)
          WRITE(NFECRA,1100) CHAINE(1:8) ,RNORM
        ENDIF
C
C
C ---> Traitement de omega
C
        IVAR   = IOMGIP
C
        IPP    = IPPRTP(IVAR)
C
        ICLVAR = ICLRTP(IVAR,ICOEF )
        ICLVAF = ICLRTP(IVAR,ICOEFF)
        CHAINE = NOMVAR(IPP)
C
        IF( IDIFF(IVAR).GE. 1 ) THEN
          DO IEL = 1, NCEL
            XXF1 = XF1(IEL)
            SIGMA = XXF1*CKWSW1 + (1.D0-XXF1)*CKWSW2
            W7(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMA
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
     &   W7     ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSE
C
          DO IFAC = 1, NFAC
            VISCF(IFAC) = 0.D0
          ENDDO
          DO IFAC = 1, NFABOR
            VISCB(IFAC) = 0.D0
          ENDDO
C
        ENDIF
C
        ICCOCG = 1
        INC    = 1
        ICONVP = ICONV (IVAR)
        IDIFFP = IDIFF (IVAR)
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IRCFLP = IRCFLU(IVAR)
        ISCHCP = ISCHCV(IVAR)
        ISSTPP = ISSTPC(IVAR)
        IWARNP = IWARNI(IVAR)
        BLENCP = BLENCV(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
        RELAXP = RELAXV(IVAR)
        THETAP = THETAV(IVAR)
C
        CALL BILSC2
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &   COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   W6  ,
C        --
     &   W7     , W8     , W1     , TINSTK , TINSTW , DRTP ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (IWARNI(IVAR).GE.2) THEN
          ISQRT = 1
          CALL PRODSC(NCELET,NCEL,ISQRT,SMBRW,SMBRW,RNORM)
          WRITE(NFECRA,1100) CHAINE(1:8) ,RNORM
        ENDIF
C
        DO IEL = 1,NCEL
          SMBRK(IEL) = SMBRK(IEL) + W5(IEL)
          SMBRW(IEL) = SMBRW(IEL) + W6(IEL)
        ENDDO
C
      ENDIF
C
C=======================================================================
C 10. AJOUT DES TERMES SOURCES DE MASSE EXPLICITES
C
C       Les parties implicites eventuelles sont conservees dans W7 et W8
C         et utilisees dans la phase d'implicitation cv/diff
C
C       Les termes sont stockes dans     SMBRK, SMBRW, W7, W8
C       En sortie de l'etape on conserve W2-8,SMBRK,SMBRW,DAM
C=======================================================================
C
      IF (NCESMP.GT.0) THEN
C
        DO IEL = 1, NCEL
          W7(IEL) = 0.D0
          W8(IEL) = 0.D0
        ENDDO
C
C       Entier egal a 1 (pour navsto : nb de sur-iter)
        IIUN = 1
C
C       On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
        IVAR = IKIPH
        CALL CATSMA
C       ===========
     & ( NCELET , NCEL   , NCESMP , IIUN   ,
     &                                 ISTO2T(IPHAS) , THETAV(IVAR) ,
     &   ICETSM , ITYPSM(1,IVAR) ,
     &   VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) , SMACEL(1,IPRIPH) ,
     &   SMBRK  , W7     , TINSTK )
        IVAR = IOMGIP
        CALL CATSMA
C       ===========
     & ( NCELET , NCEL   , NCESMP , IIUN   ,
     &                                 ISTO2T(IPHAS) , THETAV(IVAR) ,
     &   ICETSM , ITYPSM(1,IVAR) ,
     &   VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) , SMACEL(1,IPRIPH) ,
     &   SMBRW  , W8     , TINSTW )
C
C       Si on extrapole les TS on met Gamma Pinj dans PROPCE
        IF(ISTO2T(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSTA  ) = PROPCE(IEL,IPTSTA  ) + TINSTK(IEL)
            PROPCE(IEL,IPTSTA+1) = PROPCE(IEL,IPTSTA+1) + TINSTW(IEL)
          ENDDO
C       Sinon on le met directement dans SMBR
        ELSE
          DO IEL = 1, NCEL
            SMBRK(IEL) = SMBRK(IEL) + TINSTK(IEL)
            SMBRW(IEL) = SMBRW(IEL) + TINSTW(IEL)
          ENDDO
        ENDIF
C
      ENDIF
C
C     Finalisation des termes sources
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          SMBRK(IEL) = SMBRK(IEL) + THETP1 * PROPCE(IEL,IPTSTA)
          SMBRW(IEL) = SMBRW(IEL) + THETP1 * PROPCE(IEL,IPTSTA+1)
        ENDDO
      ENDIF
C
C=======================================================================
C 11. INCREMENTS DES TERMES SOURCES DANS LE SECOND MEMBRE
C
C       Les termes sont stockes dans     SMBRK, SMBRW
C       En sortie de l'etape on conserve W3-8,SMBRK,SMBRW
C=======================================================================
C
C     Ordre 2 non pris en compte
      IF(IKECOU(IPHAS).EQ.1) THEN
C
        DO IEL = 1, NCEL
C
          ROM = PROPCE(IEL,IPCROM)
C
C   RESOLUTION COUPLEE
C
          ROMVSD     = 1.D0/(ROM*VOLUME(IEL))
          SMBRK(IEL) = SMBRK(IEL)*ROMVSD
          SMBRW(IEL) = SMBRW(IEL)*ROMVSD
          DIVP23     = D2S3*MAX(DIVUKW(IEL),ZERO)
          PRODUC     = S2KW(IEL)-D2S3*DIVUKW(IEL)**2+W2(IEL)
          XK         = RTPA(IEL,IKIPH)
          XW         = RTPA(IEL,IOMGIP)
          XXF1       = XF1(IEL)
          XGAMMA     = XXF1*CKWGM1 + (1.D0-XXF1)*CKWGM2
          XBETA      = XXF1*CKWBT1 + (1.D0-XXF1)*CKWBT2
C
          A11 = 1.D0/DT(IEL)
     &         - 1.D0/XW*MIN(PRODUC,ZERO)+DIVP23+CMU*XW
          A12 = CMU*XK
          A21 = 0.D0
          A22 = 1.D0/DT(IEL)+XGAMMA*DIVP23+2.D0*XBETA*XW
C
          UNSDET = 1.D0/(A11*A22 -A12*A21)
C
          DELTK = ( A22*SMBRK(IEL) -A12*SMBRW(IEL) )*UNSDET
          DELTW = (-A21*SMBRK(IEL) +A11*SMBRW(IEL) )*UNSDET
C
C     NOUVEAU TERME SOURCE POUR CODITS
C
          ROMVSD = ROM*VOLUME(IEL)/DT(IEL)
C
          SMBRK(IEL) = ROMVSD*DELTK
          SMBRW(IEL) = ROMVSD*DELTW
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 12. TERMES INSTATIONNAIRES
C
C     Les termes sont stockes dans     TINSTK, TINSTW
C     En sortie de l'etape on conserve SMBRK, SMBRW,  TINSTK, TINSTW
C=======================================================================
C
C --- PARTIE EXPLICITE
C
C     on enleve la convection/diffusion au temps n a SMBRK et SMBRW
C     s'ils ont ete calcules
      IF (IKECOU(IPHAS).EQ.1) THEN
        DO IEL = 1, NCEL
          SMBRK(IEL) = SMBRK(IEL) - W5(IEL)
          SMBRW(IEL) = SMBRW(IEL) - W6(IEL)
        ENDDO
      ENDIF
C
C --- RHO/DT et DIV
C     Extrapolation ou non, meme forme par coherence avec bilsc2
C
      DO IEL = 1, NCEL
        ROM = PROPCE(IEL,IPCROM)
        ROMVSD = ROM*VOLUME(IEL)/DT(IEL)
        TINSTK(IEL) = ISTAT(IKIPH)*ROMVSD
     &               -ICONV(IKIPH)*W4(IEL)*THETAV(IKIPH)
        TINSTW(IEL) = ISTAT(IOMGIP)*ROMVSD
     &               -ICONV(IOMGIP)*W4(IEL)*THETAV(IOMGIP)
      ENDDO
C
C --- Source de masse (le theta est deja inclus par catsma)
      IF (NCESMP.GT.0) THEN
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) + W7(IEL)
          TINSTW(IEL) = TINSTW(IEL) + W8(IEL)
        ENDDO
      ENDIF
C
C --- Termes sources utilisateurs
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETAK = THETAV(IKIPH)
        THETAW = THETAV(IOMGIP)
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) -DAM(IEL)*THETAK
          TINSTW(IEL) = TINSTW(IEL) -W3 (IEL)*THETAW
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) + MAX(-DAM(IEL),ZERO)
          TINSTW(IEL) = TINSTW(IEL) + MAX(-W3 (IEL),ZERO)
        ENDDO
      ENDIF
C
C --- PRISE EN COMPTE DES TERMES LAGRANGIEN : COUPLAGE RETOUR
C
C     Ordre 2 non pris en compte
      IF (IILAGR.EQ.2 .AND. LTSDYN.EQ.1 .AND. IPHAS.EQ.1) THEN
C
        DO IEL = 1,NCEL
C
C Termes sources implicite sur k
C
          TINSTK(IEL) = TINSTK(IEL) + MAX(-TSLAGR(IEL,ITSLI),ZERO)
C
C Termes sources implicte sur omega
C
          TINSTW(IEL) = TINSTW(IEL)
     &          + MAX( (-CE4*TSLAGR(IEL,ITSKE)/RTPA(IEL,IKIPH)) , ZERO)
C
        ENDDO
C
      ENDIF
C
C Si IKECOU=0, on implicite plus fortement k et omega
C
      IF(IKECOU(IPHAS).EQ.0)THEN
        DO IEL=1,NCEL
          XW    = RTPA(IEL,IOMGIP)
          XXF1  = XF1(IEL)
          XBETA = XXF1*CKWBT1 + (1.D0-XXF1)*CKWBT2
          ROM = PROPCE(IEL,IPCROM)
            TINSTK(IEL) = TINSTK(IEL) +
     &         VOLUME(IEL)*CMU*ROM*XW
            TINSTW(IEL) = TINSTW(IEL) +
     &         VOLUME(IEL)*XBETA*ROM*XW
        ENDDO
      ENDIF
C
C
C=======================================================================
C 13. RESOLUTION
C
C       On utilise                      SMBRK, SMBRW,  TINSTK, TINSTW
C=======================================================================
C
C ---> Traitement de k
C
      IVAR = IKIPH
      ICLVAR = ICLRTP(IVAR,ICOEF )
      ICLVAF = ICLRTP(IVAR,ICOEFF)
C
      IPP    = IPPRTP(IVAR)
C
C     "VITESSE" DE DIFFUSION FACETTE
C
      IF( IDIFF(IVAR).GE. 1 ) THEN
C
        DO IEL = 1, NCEL
          XXF1 = XF1(IEL)
          SIGMA = XXF1*CKWSK1 + (1.D0-XXF1)*CKWSK2
          W1(IEL) = PROPCE(IEL,IPCVIS)
     &                        + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMA
        ENDDO
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
      ELSE
C
        DO IFAC = 1, NFAC
          VISCF(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
C
      ENDIF
C
C     RESOLUTION POUR K
C
      ICONVP = ICONV (IVAR)
      IDIFFP = IDIFF (IVAR)
      IRESLP = IRESOL(IVAR)
      NDIRCP = NDIRCL(IVAR)
      NITMAP = NITMAX(IVAR)
      NSWRSP = NSWRSM(IVAR)
      NSWRGP = NSWRGR(IVAR)
      IMLIGP = IMLIGR(IVAR)
      IRCFLP = IRCFLU(IVAR)
      ISCHCP = ISCHCV(IVAR)
      ISSTPP = ISSTPC(IVAR)
      IESCAP = 0
      IMGRP  = IMGR  (IVAR)
      NCYMXP = NCYMAX(IVAR)
      NITMFP = NITMGF(IVAR)
CMO      IPP    =
      IWARNP = IWARNI(IVAR)
      BLENCP = BLENCV(IVAR)
      EPSILP = EPSILO(IVAR)
      EPSRSP = EPSRSM(IVAR)
      EPSRGP = EPSRGR(IVAR)
      CLIMGP = CLIMGR(IVAR)
      EXTRAP = EXTRAG(IVAR)
      RELAXP = RELAXV(IVAR)
      THETAP = THETAV(IVAR)
C
      CALL CODITS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRSP , EPSRGP , CLIMGP , EXTRAP ,
     &   RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   TINSTK , SMBRK  , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C ---> Traitement de omega
C
      IVAR = IOMGIP
      ICLVAR = ICLRTP(IVAR,ICOEF )
      ICLVAF = ICLRTP(IVAR,ICOEFF)
C
      IPP    = IPPRTP(IVAR)
C
C
C     "VITESSE" DE DIFFUSION FACETTE
C
      IF( IDIFF(IVAR).GE. 1 ) THEN
        DO IEL = 1, NCEL
          XXF1 = XF1(IEL)
          SIGMA = XXF1*CKWSW1 + (1.D0-XXF1)*CKWSW2
          W1(IEL) = PROPCE(IEL,IPCVIS)
     &                        + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMA
        ENDDO
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
      ELSE
C
        DO IFAC = 1, NFAC
          VISCF(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
C
      ENDIF
C
C RESOLUTION POUR OMEGA
C
      ICONVP = ICONV (IVAR)
      IDIFFP = IDIFF (IVAR)
      IRESLP = IRESOL(IVAR)
      NDIRCP = NDIRCL(IVAR)
      NITMAP = NITMAX(IVAR)
      NSWRSP = NSWRSM(IVAR)
      NSWRGP = NSWRGR(IVAR)
      IMLIGP = IMLIGR(IVAR)
      IRCFLP = IRCFLU(IVAR)
      ISCHCP = ISCHCV(IVAR)
      ISSTPP = ISSTPC(IVAR)
      IESCAP = 0
      IMGRP  = IMGR  (IVAR)
      NCYMXP = NCYMAX(IVAR)
      NITMFP = NITMGF(IVAR)
CMO      IPP    =
      IWARNP = IWARNI(IVAR)
      BLENCP = BLENCV(IVAR)
      EPSILP = EPSILO(IVAR)
      EPSRSP = EPSRSM(IVAR)
      EPSRGP = EPSRGR(IVAR)
      CLIMGP = CLIMGR(IVAR)
      EXTRAP = EXTRAG(IVAR)
      RELAXP = RELAXV(IVAR)
      THETAP = THETAV(IVAR)
C
      CALL CODITS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVAR , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRSP , EPSRGP , CLIMGP , EXTRAP ,
     &   RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   TINSTW , SMBRW  , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C=======================================================================
C 14. CLIPPING
C=======================================================================
C
C     Calcul des Min/Max avant clipping, pour affichage
      DO II = 1, 2
        IF(II.EQ.1) THEN
          IVAR = IKIPH
        ELSEIF(II.EQ.2) THEN
          IVAR = IOMGIP
        ENDIF
        IPP  = IPPRTP(IVAR)
C
        VRMIN =  GRAND
        VRMAX = -GRAND
        DO IEL = 1, NCEL
          VAR = RTP(IEL,IVAR)
          VRMIN = MIN(VRMIN,VAR)
          VRMAX = MAX(VRMAX,VAR)
        ENDDO
        IF (IRANGP.GE.0) THEN
          CALL PARMAX (VRMAX)
C         ===========
          CALL PARMIN (VRMIN)
C         ===========
        ENDIF
        VARMNA(IPP) = VRMIN
        VARMXA(IPP) = VRMAX
C
      ENDDO
C
C     On clippe simplement k et omega par valeur absolue
      ICLIPK = 0
      ICLIPW = 0
      DO IEL = 1, NCEL
        XK = RTP(IEL,IKIPH )
        XW = RTP(IEL,IOMGIP)
        IF (ABS(XK).LE.EPZ2) THEN
          ICLIPK = ICLIPK + 1
          RTP(IEL,IKIPH) = MAX(RTP(IEL,IKIPH),EPZ2)
        ELSEIF(XK.LE.0.D0) THEN
          ICLIPK = ICLIPK + 1
          RTP(IEL,IKIPH) = -XK
        ENDIF
        IF (ABS(XW).LE.EPZ2) THEN
          ICLIPW = ICLIPW + 1
          RTP(IEL,IOMGIP) = MAX(RTP(IEL,IOMGIP),EPZ2)
        ELSEIF(XW.LE.0.D0) THEN
          ICLIPW = ICLIPW + 1
          RTP(IEL,IOMGIP) = -XW
        ENDIF
      ENDDO
C
      IF (IRANGP.GE.0) THEN
        CALL PARCPT (ICLIPK)
C       ===========
        CALL PARCPT (ICLIPW)
C       ===========
      ENDIF
C
C ---  Stockage nb de clippings pour listing
C
      ICLPMN(IPPRTP(IKIPH )) = ICLIPK
      ICLPMN(IPPRTP(IOMGIP)) = ICLIPW
C
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU K-OMEGA                     ',/,
     &'      ----------------------------------                     ',/)
 1100 FORMAT(1X,A8,' : BILAN EXPLICITE = ',E14.5)
C
#else
C
 1000 FORMAT(/,
     &'   ** PHASE ',I4,' SOLVING K-OMEGA'                           ,/,
     &'      ----------------------------'                           ,/)
 1100 FORMAT(1X,A8,' : EXPLICIT BALANCE = ',E14.5)
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
