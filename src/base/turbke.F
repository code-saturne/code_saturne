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
                        SUBROUTINE TURBKE
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
     &   VISCF  , VISCB  , PRDV2F ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRK  , SMBRE  , TINSTK , TINSTE ,
     &   DIVU   , W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS K-EPS 1 PHASE INCOMPRESSIBLE OU
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
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! PRDV2F(NCELET! TR !  -> ! TABLEAU DE STOCKAGE DU TERME DE      !
CARGU !              !    !     ! PROD DE TURBULENCE POUR LE V2F       !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR.(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! TINST.(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! DIVU(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION PRDV2F(NCELET)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET),   SMBRK(NCELET),  SMBRE(NCELET)
      DOUBLE PRECISION TINSTK(NCELET), TINSTE(NCELET), DIVU(NCELET)
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
      INTEGER          IIVAR , IIUN
      INTEGER          ICLIP , ISQRT
      INTEGER          NSWRGP, IMLIGP, IPHYDP
      INTEGER          IPRIPH, IUIPH , IVIPH , IWIPH
      INTEGER          IKIPH , IEIPH , IPHIPH
      INTEGER          ICLIUP, ICLIVP, ICLIWP
      INTEGER          ICLVAR, ICLVAF
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP
      INTEGER          NITMAP, NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IPCROM, IPBROM, IPCVST, IPCVIS, IFLMAS, IFLMAB
      INTEGER          IWARNP, IPP
      INTEGER          IPTSTA
      INTEGER          IPCROO, IPBROO, IPCVTO, IPCVLO
      INTEGER          MAXELT, ILS
      DOUBLE PRECISION RNORM , D2S3, DIVP23
      DOUBLE PRECISION DELTK , DELTE, A11, A12, A22, A21
      DOUBLE PRECISION GRAVKE, EPSSUK, UNSDET, ROMVSD
      DOUBLE PRECISION PRDTUR, XK, XEPS, XPHI, XNU, TTKE, TTMIN, TT
      DOUBLE PRECISION VISCT , ROM   , CEPS1 , CTSQNU
      DOUBLE PRECISION BLENCP, EPSILP, EPSRGP, CLIMGP, EXTRAP, RELAXP
      DOUBLE PRECISION EPSRSP
      DOUBLE PRECISION THETP1, THETAK, THETAE, THETS, THETAP
      DOUBLE PRECISION TUEXPK, TUEXPE
      DOUBLE PRECISION CMUETA, SQRCMU, XS
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
      IPRIPH = IPR (IPHAS)
      IUIPH  = IU  (IPHAS)
      IVIPH  = IV  (IPHAS)
      IWIPH  = IW  (IPHAS)
      IKIPH  = IK  (IPHAS)
      IEIPH  = IEP (IPHAS)
      IPHIPH = IPHI(IPHAS)
C
      ICLIUP = ICLRTP(IUIPH,ICOEF)
      ICLIVP = ICLRTP(IVIPH,ICOEF)
      ICLIWP = ICLRTP(IWIPH,ICOEF)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IUIPH))
      IFLMAB = IPPROB(IFLUMA(IUIPH))
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
        IF (ITURB(IPHAS).EQ.20) THEN
          WRITE(NFECRA,1000)IPHAS
        ELSE IF (ITURB(IPHAS).EQ.21) THEN
          WRITE(NFECRA,1001)IPHAS
        ELSE
          WRITE(NFECRA,1002)IPHAS
        ENDIF
      ENDIF
C
C     Pour la prod lineaire on a besoin de SQRT(Cmu)
      SQRCMU = SQRT(CMU)
C
C=======================================================================
C 2. CALCUL DE SIJ SIJ ET DE DIVU
C
C      Tableaux de travail              SMBRK,TINSTE,W1,W2,W3,W4,W5,W6
C      SijSij est stocke dans           TINSTK
C      DivU est stocke dans             DIVU
C      En sortie de l'etape on conserve TINSTK, DIVU
C=======================================================================
C
      ICCOCG = 1
      INC = 1
C
C ON S'APPUIE SUR LES TABLEAUX DE TRAVAIL W4,W5,W6
C     ON UTILISE EGALEMENT SMBRK  ET TINSTE COMME
C                           TABLEAUX DE TRAVAIL TEMPORAIRES
C     W1,W2,W3 SONT UTILISES DANS GRDCEL
C     TINSTK RECOIT LA PRODUCTION
C        TINSTK = {(2 (S11)**2 + 2 (S22)**2 +2 (S33)**2 )
C           +         ((2 S12)**2 + (2 S13)**2 +(2 S23)**2 )}
C               = 2 Sij.Sij
C
C
C SMBRK  = DUDX ,W4 = DUDY ,W5 = DUDZ
C
      NSWRGP = NSWRGR(IUIPH)
      IMLIGP = IMLIGR(IUIPH)
      IWARNP = IWARNI(IKIPH)
      EPSRGP = EPSRGR(IUIPH)
      CLIMGP = CLIMGR(IUIPH)
      EXTRAP = EXTRAG(IUIPH)
      IPHYDP = 0
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IUIPH  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IUIPH)   , COEFA(1,ICLIUP) , COEFB(1,ICLIUP) ,
     &   SMBRK  , W4     , W5     ,
C        ------   ------   ------
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C TINSTK = (S11)**2
C DIVU = DUDX
C
      DO IEL = 1, NCEL
        TINSTK(IEL)   = SMBRK(IEL)**2
        DIVU  (IEL)   = SMBRK(IEL)
      ENDDO
C
C               ,W4     = DUDY ,W5     = DUDZ
C TINSTE = DVDX ,SMBRK  = DVDY ,W6     = DVDZ
C
      NSWRGP = NSWRGR(IVIPH)
      IMLIGP = IMLIGR(IVIPH)
      IWARNP = IWARNI(IKIPH)
      EPSRGP = EPSRGR(IVIPH)
      CLIMGP = CLIMGR(IVIPH)
      EXTRAP = EXTRAG(IVIPH)
      IPHYDP = 0
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVIPH  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IVIPH)   , COEFA(1,ICLIVP) , COEFB(1,ICLIVP) ,
     &   TINSTE , SMBRK  , W6     ,
C        ------   ------   ------
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C TINSTK = 2 (S11)**2 + 2 (S22)**2
C      + (2 S12)**2
C DIVU = DUDX + DVDY
C
      DO IEL = 1, NCEL
        TINSTK (IEL)   = 2.D0*(TINSTK(IEL) + SMBRK(IEL)**2 )
     &           +  (TINSTE(IEL)+W4(IEL))**2
        DIVU  (IEL)   = DIVU(IEL) + SMBRK(IEL)
      ENDDO
C
C               ,              ,W5     = DUDZ
C               ,              ,W6     = DVDZ
C W4     = DWDX ,TINSTE = DWDY ,SMBRK  = DWDZ
C
      NSWRGP = NSWRGR(IWIPH)
      IMLIGP = IMLIGR(IWIPH)
      IWARNP = IWARNI(IKIPH)
      EPSRGP = EPSRGR(IWIPH)
      CLIMGP = CLIMGR(IWIPH)
      EXTRAP = EXTRAG(IWIPH)
      IPHYDP = 0
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IWIPH  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IWIPH)   , COEFA(1,ICLIWP) , COEFB(1,ICLIWP) ,
     &   W4     , TINSTE , SMBRK  ,
C        ------   ------   ------
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C TINSTK = PRODUCTION = (
C      2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
C      + (2 S12)**2 + (2 S13)**2 + (2 S23)**2 )
C DIVU = DUDX + DVDY + DWDZ
C
      DO IEL = 1, NCEL
        TINSTK (IEL)   = TINSTK(IEL) + 2.D0*SMBRK(IEL)**2
     &           +  (W4(IEL)+W5(IEL))**2
     &           +  (TINSTE(IEL)+W6(IEL))**2
        DIVU  (IEL)   = DIVU(IEL) + SMBRK(IEL)
      ENDDO
C
C On libere   SMBRK,TINSTE,W1,W2,W3,W4,W5,W6
C
C======================================================================
C 3. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
C
C      On passe 2 Sij.Sij = TINSTK et la divergence DIVU
C      Tableaux de travail                        W1, W2, W3, W4, W5, W6
c                                VISCF VISCB XAM DRTP SMBRK SMBRE TINSTE
C      La partie a expliciter est stockee dans    W7, W8
C      La partie a impliciter est stockee dans    DAM, W9
C      En sortie de l'etape on conserve           TINSTK, DIVU,
C                                                 W7 , W8, DAM, W9
C=======================================================================
C viscf viscb xam drtp smbrk smbre tinste
      DO IEL = 1, NCEL
        DAM(IEL) = 0.D0
        W9 (IEL) = 0.D0
        W7 (IEL) = 0.D0
        W8 (IEL) = 0.D0
      ENDDO
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('TURBKE',IFINIA)
C
      CALL USTSKE
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
     &   COEFA  , COEFB  , CKUPDC , SMACEL , TINSTK , DIVU   ,
     &   W7     , W8     , DAM    , W9     ,
C        ------   ------   ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   DRTP   , SMBRK  , SMBRE  , TINSTE ,
     &   RDEVEL , RTUSER , RA     )
C
C On libere W1, W2, W3, W4, W5, W6
C           VISCF, VISCB, XAM, DRTP, SMBRK, SMBRE, TINSTE
C=======================================================================
C 4. AJOUT DE - 2/3 DIV(U) * DIV(U)
C
C      En sortie de l'etape on conserve TINSTK, DIVU,
C                                       W7 , W8, DAM, W9
C=======================================================================
C
C Dans le cas de la production lineaire, seul le terme en divu est
C multiplie par VISCT. Pour les autres modeles, la multiplication par
C VISCT sera faite ulterieurement.
C A ce stade, TINSTK contient S**2
      D2S3 = 2.D0/3.D0
      IF (ITURB(IPHAS).EQ.21) THEN
        DO IEL = 1, NCEL
          ROM   = PROPCE(IEL,IPCROO)
          VISCT = PROPCE(IEL,IPCVTO)
          XS = SQRT(TINSTK(IEL))
          CMUETA = CMU*RTPA(IEL,IKIPH)/RTPA(IEL,IEIPH)*XS
          CMUETA = MIN(CMUETA,SQRCMU)
          TINSTK(IEL) = ROM*CMUETA*XS*RTPA(IEL,IKIPH)
     &         - D2S3*VISCT*DIVU(IEL)*DIVU(IEL)
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) - D2S3*DIVU(IEL)*DIVU(IEL)
        ENDDO
      ENDIF
C
C=======================================================================
C 5. CALCUL DU TERME DE GRAVITE
C
C      Les s.m. recoivent production et termes de gravite
C      Tableaux de travail              W1,W2,W3,W4,W5,W6,VISCB
C      Les s.m. sont stockes dans       TINSTK TINSTE
C      En sortie de l'etape on conserve TINSTK, TINSTE,
C                                       DIVU,
C                                       W7 , W8, DAM, W9
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
     &   W1     , W1     , W1     ,
     &   PROPCE(1,IPCROO), PROPFB(1,IPBROO), VISCB  ,
     &   W4     , W5     , W6     ,
C        ------   ------   ------
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C      Production et terme de gravite
C        TINSTK=P+G et TINSTE=P+(1-CE3)*G
C
        IF(ISCALT(IPHAS).GT.0.AND.NSCAL.GE.ISCALT(IPHAS)) THEN
          PRDTUR = SIGMAS(ISCALT(IPHAS))
        ELSE
          PRDTUR = 1.D0
        ENDIF
C
C     En production lineaire, on multiplie tout de suite le terme
C     de gravite par VISCT, car le reste est deja multiplie.
C     Dans les autres cas, la multiplication est faite plus tard.
        IF (ITURB(IPHAS).EQ.21) THEN
          DO IEL = 1, NCEL
            GRAVKE = -(W4(IEL)*GX+W5(IEL)*GY+W6(IEL)*GZ)/
     &           (PROPCE(IEL,IPCROO)*PRDTUR)
            TINSTE(IEL)=TINSTK(IEL)+PROPCE(IEL,IPCVTO)*MAX(GRAVKE,ZERO)
            TINSTK(IEL)=TINSTK(IEL)+PROPCE(IEL,IPCVTO)*GRAVKE
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            GRAVKE = -(W4(IEL)*GX+W5(IEL)*GY+W6(IEL)*GZ)/
     &           (PROPCE(IEL,IPCROO)*PRDTUR)
            TINSTE(IEL) = TINSTK(IEL) + MAX( GRAVKE,ZERO )
            TINSTK(IEL) = TINSTK(IEL) + GRAVKE
          ENDDO
        ENDIF
C
      ELSE
C
C
C --- Production sans termes de gravite
C       TINSTK=TINSTE=P
C
        DO IEL = 1, NCEL
          TINSTE(IEL) = TINSTK(IEL)
        ENDDO
C
      ENDIF
C
C     En V2F, on stocke TINSTK dans PRDV2F qui sera complete plus loin pour
C     contenir le terme de production complet
      IF (ITURB(IPHAS).EQ.50) THEN
        DO IEL = 1, NCEL
          PRDV2F(IEL) = TINSTK(IEL)
        ENDDO
      ENDIF
C
C On libere W1,W2,W3,W4,W5,W6,VISCB
C
C=======================================================================
C 6. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume
C
C      Le terme est stocke dans         W1
C      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
C                                       W7 , W8, DAM, W9
C=======================================================================
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W1)
C
C=======================================================================
C 7. ON FINALISE LE CALCUL DES TERMES SOURCES
C
C      Les termes sont stockes dans     SMBRK, SMBRE
C      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
C                                       SMBRK, SMBRE
C                                       W7 , W8, DAM, W9
C=======================================================================
C
C SMBRE = CEPS1 EPSILON/K (PROD + G ) - RO0 VOLUME EPSILON EPSILON/K
C SMBRK = PROD + G                  - RO0 VOLUME EPSILON
C
C     Si on extrapole les termes sources et rho  , il faut ici rho^n
C                                        et visct, il faut ici visct^n
C
      IF (ITURB(IPHAS).EQ.20) THEN
C
        DO IEL = 1, NCEL
C
          VISCT = PROPCE(IEL,IPCVTO)
          ROM   = PROPCE(IEL,IPCROO)
C
          SMBRK(IEL) = VOLUME(IEL)*(
     &         VISCT*TINSTK(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL)
     &         -ROM*RTPA(IEL,IEIPH) )
C
          SMBRE(IEL) = VOLUME(IEL)*RTPA(IEL,IEIPH)/RTPA(IEL,IKIPH)*(
     &         CE1*( VISCT*TINSTE(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL) )
     &         -CE2*ROM*RTPA(IEL,IEIPH) )
C
        ENDDO
C
      ELSE IF (ITURB(IPHAS).EQ.21) THEN
C
        DO IEL = 1, NCEL
C
          ROM   = PROPCE(IEL,IPCROO)
C
          SMBRK(IEL) = VOLUME(IEL)*(
     &         TINSTK(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL)
     &         -ROM*RTPA(IEL,IEIPH) )
C
          SMBRE(IEL) = VOLUME(IEL)*RTPA(IEL,IEIPH)/RTPA(IEL,IKIPH)*(
     &         CE1*(TINSTE(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL) )
     &         -CE2*ROM*RTPA(IEL,IEIPH) )
C
        ENDDO
C
      ELSE IF (ITURB(IPHAS).EQ.50) THEN
C
        DO IEL = 1, NCEL
C
          VISCT = PROPCE(IEL,IPCVTO)
          ROM   = PROPCE(IEL,IPCROO)
          XEPS = RTPA(IEL,IEIPH )
          XK   = RTPA(IEL,IKIPH )
          XPHI = RTPA(IEL,IPHIPH)
          XPHI = MAX(XPHI,EPZERO)
          XNU  = PROPCE(IEL,IPCVLO)/ROM
          CEPS1= 1.4D0*(1.D0+CV2FA1*SQRT(1.D0/XPHI))
          TTKE = XK / XEPS
          TTMIN = CV2FCT*SQRT(XNU/XEPS)
          TT = MAX(TTKE,TTMIN)
C
          SMBRK(IEL) = VOLUME(IEL)*(
     &         VISCT*TINSTK(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL)
     &         -ROM*RTPA(IEL,IEIPH) )
C
          SMBRE(IEL) = VOLUME(IEL)/TT*(
     &         CEPS1*( VISCT*TINSTE(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL) )
     &         -CV2FE2*ROM*RTPA(IEL,IEIPH) )
C
C     On stocke la partie en Pk dans PRDV2F pour etre reutilise dans RESV2F
          PRDV2F(IEL) = VISCT*PRDV2F(IEL)
     &         -D2S3*ROM*RTPA(IEL,IKIPH)*DIVU(IEL)
C
        ENDDO
C
      ENDIF
C
C
C======================================================================
C 8. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
C                        ET ACCUMULATION DE MASSE    : PARTIE EXPLICITE
C      On utilise                       W1,  W7, W8, DAM, W9
C      Les termes sont stockes dans     SMBRK, SMBRE
C      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
C                                       SMBRK, SMBRE
C                                       DAM, W9
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
          PROPCE(IEL,IPTSTA) = SMBRK(IEL) + W7(IEL)
C       Termes dependant de la variable resolue et theta PROPCE
          SMBRK(IEL) = ICONV(IKIPH)*W1(IEL)*RTPA(IEL,IKIPH)
     &         - THETS*TUEXPK
C       On suppose -DAM > 0 : on implicite
C         le terme utilisateur dependant de la variable resolue
          SMBRK(IEL) = DAM(IEL)*RTPA(IEL,IKIPH) + SMBRK(IEL)
C
C       Sauvegarde pour echange
          TUEXPE = PROPCE(IEL,IPTSTA+1)
C       Pour la suite et le pas de temps suivant
          PROPCE(IEL,IPTSTA+1) = SMBRE(IEL) + W8(IEL)
C       Termes dependant de la variable resolue et theta PROPCE
          SMBRE(IEL) = ICONV(IEIPH)*W1(IEL)*RTPA(IEL,IEIPH)
     &         - THETS*TUEXPE
C       On suppose -W9 > 0 : on implicite
C         le terme utilisateur dependant de la variable resolue
          SMBRE(IEL) =  W9(IEL)*RTPA(IEL,IEIPH) + SMBRE(IEL)
C
        ENDDO
C
C     Si on n'extrapole pas les T.S.
      ELSE
        DO IEL = 1, NCEL
          SMBRK(IEL) = SMBRK(IEL) + DAM(IEL)*RTPA(IEL,IKIPH) + W7(IEL)
     &         +ICONV(IKIPH)*W1(IEL)*RTPA(IEL,IKIPH)
          SMBRE(IEL) = SMBRE(IEL) + W9 (IEL)*RTPA(IEL,IEIPH) + W8(IEL)
     &         +ICONV(IEIPH)*W1(IEL)*RTPA(IEL,IEIPH)
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
C Termes sources explicte sur Eps
C
          SMBRE(IEL)  = SMBRE(IEL)
     &                + CE4 *TSLAGR(IEL,ITSKE) *RTPA(IEL,IEIPH)
     &                                         /RTPA(IEL,IKIPH)
C
        ENDDO
C
      ENDIF
C
C ON LIBERE W7, W8
C
C=======================================================================
C 9. PRISE EN COMPTE DES TERMES DE CONV/DIFF DANS LE SECOND MEMBRE
C
C      Tableaux de travail              W2, W3, W4, W5, W6, DRTP
C      Les termes sont stockes dans     W7 et W8, puis ajoutes a SMBRK, SMBRE
C      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
C                                       SMBRK, SMBRE
C                                       DAM, W7, W8, W9
C=======================================================================
C
C     Ceci ne sert a rien si IKECOU n'est pas egal a 1
C
      IF (IKECOU(IPHAS).EQ.1) THEN
C
        DO IEL = 1, NCEL
          W7 (IEL) = 0.D0
          W8 (IEL) = 0.D0
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
            W4(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMAK
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
     &   W4     ,
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
     &   W7  ,
C        --
     &   W2     , W3     , W4     , W5     , W6     , DRTP ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (IWARNI(IVAR).GE.2) THEN
          ISQRT = 1
          CALL PRODSC(NCELET,NCEL,ISQRT,SMBRK,SMBRK,RNORM)
          WRITE(NFECRA,1100) CHAINE(1:8) ,RNORM
        ENDIF
C
C
C ---> Traitement de epsilon
C
        IVAR   = IEIPH
C
        IPP    = IPPRTP(IVAR)
C
        ICLVAR = ICLRTP(IVAR,ICOEF )
        ICLVAF = ICLRTP(IVAR,ICOEFF)
        CHAINE = NOMVAR(IPP)
C
        IF( IDIFF(IVAR).GE. 1 ) THEN
          DO IEL = 1, NCEL
            W4(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMAE
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
     &   W4     ,
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
     &   W8  ,
C        --
     &   W2     , W3     , W4     , W5     , W6     , DRTP ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (IWARNI(IVAR).GE.2) THEN
          ISQRT = 1
          CALL PRODSC(NCELET,NCEL,ISQRT,SMBRE,SMBRE,RNORM)
          WRITE(NFECRA,1100) CHAINE(1:8) ,RNORM
        ENDIF
C
        DO IEL = 1,NCEL
          SMBRK(IEL) = SMBRK(IEL) + W7(IEL)
          SMBRE(IEL) = SMBRE(IEL) + W8(IEL)
        ENDDO
C
      ENDIF
C
C=======================================================================
C 10. AJOUT DES TERMES SOURCES DE MASSE EXPLICITES
C
C       Les parties implicites eventuelles sont conservees dans W2 et W3
C         et utilisees dans la phase d'implicitation cv/diff
C
C       Les termes sont stockes dans     SMBRK, SMBRE, W2, W3
C       En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
C                                        SMBRK, SMBRE
C                                        DAM, W9, W2, W3
C=======================================================================
C
      IF (NCESMP.GT.0) THEN
C
        DO IEL = 1, NCEL
          W2(IEL) = 0.D0
          W3(IEL) = 0.D0
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
     &   SMBRK  , W2     , W4 )
        IVAR = IEIPH
        CALL CATSMA
C       ===========
     & ( NCELET , NCEL   , NCESMP , IIUN   ,
     &                                 ISTO2T(IPHAS) , THETAV(IVAR) ,
     &   ICETSM , ITYPSM(1,IVAR) ,
     &   VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) , SMACEL(1,IPRIPH) ,
     &   SMBRE  , W3     , W5 )
C
C       Si on extrapole les TS on met Gamma Pinj dans PROPCE
        IF(ISTO2T(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSTA  ) = PROPCE(IEL,IPTSTA  ) + W4(IEL)
            PROPCE(IEL,IPTSTA+1) = PROPCE(IEL,IPTSTA+1) + W5(IEL)
          ENDDO
C       Sinon on le met directement dans SMBR
        ELSE
          DO IEL = 1, NCEL
            SMBRK(IEL) = SMBRK(IEL) + W4(IEL)
            SMBRE(IEL) = SMBRE(IEL) + W5(IEL)
          ENDDO
        ENDIF
C
      ENDIF
C
C     ON LIBERE                       W4, W5, W6, TINSTE
C
C     Finalisation des termes sources
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          SMBRK(IEL) = SMBRK(IEL) + THETP1 * PROPCE(IEL,IPTSTA)
          SMBRE(IEL) = SMBRE(IEL) + THETP1 * PROPCE(IEL,IPTSTA+1)
        ENDDO
      ENDIF
C
C=======================================================================
C 11. INCREMENTS DES TERMES SOURCES DANS LE SECOND MEMBRE
C
C       On utilise                       TINSTK, TINSTE, DIVU
C       Les termes sont stockes dans     SMBRK, SMBRE
C       En sortie de l'etape on conserve W1, SMBRK, SMBRE,
C                                        DAM, W9, W2, W3, W7, W8
C=======================================================================
C
C     Ordre 2 non pris en compte
      IF(IKECOU(IPHAS).EQ.1) THEN
C
        IF (ITURB(IPHAS).EQ.20) THEN
C
          DO IEL = 1, NCEL
C
            ROM = PROPCE(IEL,IPCROM)
C
C   RESOLUTION COUPLEE
C
            ROMVSD=1.D0/(ROM*VOLUME(IEL))
            SMBRK(IEL)=SMBRK(IEL)*ROMVSD
            SMBRE(IEL)=SMBRE(IEL)*ROMVSD
            DIVP23= D2S3*MAX(DIVU(IEL),ZERO)
C
            EPSSUK = RTPA(IEL,IEIPH)/RTPA(IEL,IKIPH)
C
            A11 = 1.D0/DT(IEL)
     &           -2.D0*RTPA(IEL,IKIPH)/RTPA(IEL,IEIPH)
     &           *CMU*MIN(TINSTK(IEL),ZERO)+DIVP23
            A12 = 1.D0
            A21 = -CE1*CMU*TINSTE(IEL)-CE2*EPSSUK*EPSSUK
            A22 = 1.D0/DT(IEL)+CE1*DIVP23
     &           +2.D0*CE2*EPSSUK
C
            UNSDET = 1.D0/(A11*A22 -A12*A21)
C
            DELTK = ( A22*SMBRK(IEL) -A12*SMBRE(IEL) )*UNSDET
            DELTE = (-A21*SMBRK(IEL) +A11*SMBRE(IEL) )*UNSDET
C
C     NOUVEAU TERME SOURCE POUR CODITS
C
            ROMVSD = ROM*VOLUME(IEL)/DT(IEL)
C
            SMBRK(IEL) = ROMVSD*DELTK
            SMBRE(IEL) = ROMVSD*DELTE
C
          ENDDO
C
C     Dans verini on bloque la combinaison ITURB=21/IKECOU=1
        ELSE IF (ITURB(IPHAS).EQ.21) THEN
C
          WRITE(NFECRA,*)'IKECOU=1 NON VALIDE EN K-EPS PROD LIN'
          CALL CSEXIT (1)
C  Section non totalement validee (a priori ca marche, mais pas trop stable) :
C  en fait le v2f est meilleur avec IKECOU=0, on bloque donc la combinaison
C  ITURB=50/IKECOU=1 au niveau de verini. Ces lignes sont donc inaccessibles.
C  On les laisse au cas ou .....
        ELSE IF (ITURB(IPHAS).EQ.50) THEN
C
          DO IEL = 1, NCEL
C
            ROM = PROPCE(IEL,IPCROM)
C
C   RESOLUTION COUPLEE
C
            ROMVSD=1.D0/(ROM*VOLUME(IEL))
            SMBRK(IEL)=SMBRK(IEL)*ROMVSD
            SMBRE(IEL)=SMBRE(IEL)*ROMVSD
            DIVP23= D2S3*MAX(DIVU(IEL),ZERO)
C
            XEPS = RTPA(IEL,IEIPH )
            XK   = RTPA(IEL,IKIPH )
            XPHI = RTPA(IEL,IPHIPH)
            XPHI = MAX(XPHI,EPZERO)
            XNU  = PROPCE(IEL,IPCVIS)/PROPCE(IEL,IPCROM)
            CTSQNU= CV2FCT*SQRT(XNU)
            CEPS1= 1.4D0*(1.D0+CV2FA1*SQRT(1.D0/XPHI))
            EPSSUK = XEPS/XK
            TTKE = XK / XEPS
            TTMIN = CV2FCT*SQRT(XNU/XEPS)
C
            IF(TTKE.GT.TTMIN) THEN
              A11 = 1.D0/DT(IEL)
     &             -2.D0*XK/XEPS*XPHI
     &             *CV2FMU*MIN(TINSTK(IEL),ZERO)+DIVP23
C     Pour A12 on fait comme en k-eps standard pour l'instant,
C     on ne prend pas le terme en P+G ... est-ce judicieux ?
              A12 = 1.D0
              A21 = -CEPS1*CV2FMU*XPHI*TINSTE(IEL)-CV2FE2*EPSSUK*EPSSUK
              A22 = 1.D0/DT(IEL)+CEPS1*DIVP23
     &             +2.D0*CV2FE2*EPSSUK
            ELSE
              A11 = 1.D0/DT(IEL)
     &             -CV2FMU*XPHI*CTSQNU*MIN(TINSTK(IEL),ZERO)/SQRT(XEPS)
     &             +DIVP23
C     Pour A12 on fait comme en k-eps standard pour l'instant,
C     on ne prend pas le terme en P+G ... est-ce judicieux ?
              A12 = 1.D0
C     Le terme en DIVP23 dans A21 n'est pas forcement judicieux
C     (a-t-on besoin du MAX ?)
              A21 = -CEPS1*CV2FMU*XPHI*TINSTE(IEL)
     &             +CEPS1*SQRT(XEPS)/CTSQNU*DIVP23
              A22 = 1.D0/DT(IEL)+1.D0/2.D0*CEPS1*DIVP23*XK
     &             /CTSQNU/SQRT(XEPS)
     &             +3.D0/2.D0*CV2FE2/CTSQNU*SQRT(XEPS)
            ENDIF
C
            UNSDET = 1.D0/(A11*A22 -A12*A21)
C
            DELTK = ( A22*SMBRK(IEL) -A12*SMBRE(IEL) )*UNSDET
            DELTE = (-A21*SMBRK(IEL) +A11*SMBRE(IEL) )*UNSDET
C
C     NOUVEAU TERME SOURCE POUR CODITS
C
            ROMVSD = ROM*VOLUME(IEL)/DT(IEL)
C
            SMBRK(IEL) = ROMVSD*DELTK
            SMBRE(IEL) = ROMVSD*DELTE
C
          ENDDO

        ENDIF
C
      ENDIF
C
C     ON LIBERE                       TINSTK, TINSTE, DIVU
C
C=======================================================================
C 12. TERMES INSTATIONNAIRES
C
C     On utilise                       W1, W2, W3, W7, W8
C                                      DAM, W9
C     Les termes sont stockes dans     TINSTK, TINSTE
C     En sortie de l'etape on conserve SMBRK, SMBRE,  TINSTK, TINSTE
C=======================================================================
C
C --- PARTIE EXPLICITE
C
C     on enleve la convection/diffusion au temps n a SMBRK et SMBRE
C     si on les avait calcules
      IF (IKECOU(IPHAS).EQ.1) THEN
        DO IEL = 1, NCEL
          SMBRK(IEL) = SMBRK(IEL) - W7(IEL)
          SMBRE(IEL) = SMBRE(IEL) - W8(IEL)
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
     &               -ICONV(IKIPH)*W1(IEL)*THETAV(IKIPH)
        TINSTE(IEL) = ISTAT(IEIPH)*ROMVSD
     &               -ICONV(IEIPH)*W1(IEL)*THETAV(IEIPH)
      ENDDO
C
C --- Source de masse (le theta est deja inclus par catsma)
      IF (NCESMP.GT.0) THEN
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) + W2(IEL)
          TINSTE(IEL) = TINSTE(IEL) + W3(IEL)
        ENDDO
      ENDIF
C
C --- Termes sources utilisateurs
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETAK = THETAV(IKIPH)
        THETAE = THETAV(IEIPH)
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) -DAM(IEL)*THETAK
          TINSTE(IEL) = TINSTE(IEL) -W9 (IEL)*THETAE
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          TINSTK(IEL) = TINSTK(IEL) + MAX(-DAM(IEL),ZERO)
          TINSTE(IEL) = TINSTE(IEL) + MAX(-W9 (IEL),ZERO)
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
C Termes sources implicte sur Eps
C
          TINSTE(IEL) = TINSTE(IEL)
     &          + MAX( (-CE4*TSLAGR(IEL,ITSKE)/RTPA(IEL,IKIPH)) , ZERO)
C
        ENDDO
C
      ENDIF
C
C Si IKECOU=0, on implicite plus fortement k et eps
C
      IF(IKECOU(IPHAS).EQ.0)THEN
        IF(ITYTUR(IPHAS).EQ.2)THEN
          DO IEL=1,NCEL
            XEPS = RTPA(IEL,IEIPH )
            XK   = RTPA(IEL,IKIPH )
            ROM = PROPCE(IEL,IPCROM)
            TTKE = XK / XEPS
            IF(XK.GT.1.D-12) THEN
              TINSTK(IEL) = TINSTK(IEL) +
     &             ROM*VOLUME(IEL)/TTKE
            ENDIF
            TINSTE(IEL) = TINSTE(IEL) +
     &           CE2*ROM*VOLUME(IEL)/TTKE
          ENDDO
        ELSE IF(ITURB(IPHAS).EQ.50)THEN
          DO IEL=1,NCEL
            XEPS = RTPA(IEL,IEIPH )
            XK   = RTPA(IEL,IKIPH )
            ROM = PROPCE(IEL,IPCROM)
            XNU  = PROPCE(IEL,IPCVIS)/ROM
            TTKE = XK / XEPS
            TTMIN = CV2FCT*SQRT(XNU/XEPS)
            TT = MAX(TTKE,TTMIN)
            IF(XK.GT.1.D-12) THEN
              TINSTK(IEL) = TINSTK(IEL) +
     &             ROM*VOLUME(IEL)/TTKE
            ENDIF
            TINSTE(IEL) = TINSTE(IEL) +
     &           CV2FE2*ROM*VOLUME(IEL)/TT
          ENDDO
C
        ENDIF
      ENDIF
C
C ON LIBERE W1, W2, W3, DAM, W9
C
C=======================================================================
C 13. RESOLUTION
C
C       On utilise                      SMBRK, SMBRE,  TINSTK, TINSTE
C       Tableaux de travail             W1, W2, W3, W4, W5, W6
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
          W1(IEL) = PROPCE(IEL,IPCVIS)
     &                        + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMAK
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
C ---> Traitement de epsilon
C
      IVAR = IEIPH
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
          W1(IEL) = PROPCE(IEL,IPCVIS)
     &                        + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMAE
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
C RESOLUTION POUR EPSILON
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
     &   TINSTE , SMBRE  , RTP(1,IVAR)     ,
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
      ICLIP = 1
      IWARNP = IWARNI(IKIPH)
      CALL CLIPKE
C     ===========
     & ( NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IPHAS  , ICLIP  , IWARNP ,
     &   PROPCE , RTP    )
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU K-EPSILON                   ',/,
     &'      ------------------------------------                   ',/)
 1001 FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU K-EPSILON A PROD LINEAIRE   ',/,
     &'      ----------------------------------------------------   ',/)
 1002 FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU V2F (K ET EPSILON)          ',/,
     &'      ---------------------------------------------          ',/)
 1100 FORMAT(1X,A8,' : BILAN EXPLICITE = ',E14.5)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
