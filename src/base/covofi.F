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
                        SUBROUTINE COVOFI
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  , ITSPDV ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , TSLAGR ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
CFONC    POUR UN SCALAIRE SUR UN PAS DE TEMPS
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
CARGU ! ISCAL        ! E  !  -> ! NUMERO DU SCALAIRE                   !
CARGU ! ITSPDV       ! E  !  -> ! CALCUL TERMES SOURCES PROD ET DISSIP !
CARGU !              !    !     !  (0 : NON , 1 : OUI)                 !
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
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBRS(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
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
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "radiat.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
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
      INTEGER          ISCAL  , ITSPDV
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
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR  , IPHAS
      INTEGER          IFAC  , IEL
      INTEGER          INIT  , INC   , ICCOCG, ISQRT, III, IIUN, IBCL
      INTEGER          IVARSC, ISCALA
      INTEGER          IISCAV, IICP
      INTEGER          IKIPH , IEIPH , IOMGIP
      INTEGER          IR11IP, IR22IP, IR33IP
      INTEGER          ICLVAR, ICLVAF
      INTEGER          IPCROM, IPCVST, IPCVSL, IFLMAS, IFLMAB
      INTEGER          IPPVAR, IPP   , IPHYDP, IPTSCA, IPCVSO
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP, NITMAP
      INTEGER          NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          MAXELT, ILS   , IDBIA1
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, RELAXP, BLENCP, EPSILP
      DOUBLE PRECISION EPSRSP
      DOUBLE PRECISION RHOVST, XK    , XE    , SCLNOR
      DOUBLE PRECISION THETV , THETS , THETAP, THETP1
      DOUBLE PRECISION SMBEXP
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
C ---  Numero de phase associee au scalaire traite
      IPHAS  = IPHSCA(ISCAL)
C
C --- Numero de variable de calcul et de post associe au scalaire traite
      IVAR   = ISCA(ISCAL)
      IPPVAR = IPPRTP(IVAR)
C
C --- Numero du scalaire eventuel associe dans le cas fluctuation
C         et numero de variable de calcul
      IISCAV = ISCAVR(ISCAL)
      IF(IISCAV.GT.0.AND.IISCAV.LE.NSCAL) THEN
        IVARSC = ISCA(IISCAV)
      ELSE
        IVARSC = 0
      ENDIF
C
C --- Numero des variables de calcul
      IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50) THEN
        IKIPH  = IK  (IPHAS)
        IEIPH  = IEP (IPHAS)
      ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
        IR11IP = IR11(IPHAS)
        IR22IP = IR22(IPHAS)
        IR33IP = IR33(IPHAS)
        IEIPH  = IEP (IPHAS)
      ELSEIF(ITURB(IPHAS).EQ.60) THEN
        IKIPH  = IK  (IPHAS)
        IOMGIP = IOMG(IPHAS)
      ENDIF
C
C --- Numero des conditions aux limites
      ICLVAR = ICLRTP(IVAR,ICOEF)
      ICLVAF = ICLRTP(IVAR,ICOEFF)
C
C --- Numero des grandeurs physiques
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IVAR ))
      IFLMAB = IPPROB(IFLUMA(IVAR ))
      IF(IVISLS(ISCAL).GT.0) THEN
        IPCVSL = IPPROC(IVISLS(ISCAL))
      ELSE
        IPCVSL = 0
      ENDIF
C
C --- Numero du terme source dans PROPCE si extrapolation
      IF(ISSO2T(ISCAL).GT.0) THEN
        IPTSCA = IPPROC(ITSSCA(ISCAL))
      ELSE
        IPTSCA = 0
      ENDIF
C
C     S pour Source, V pour Variable
      THETS  = THETSS(ISCAL)
      THETV  = THETAV(IVAR )
C
      CHAINE = NOMVAR(IPPVAR)
C
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1000) CHAINE(1:8)
      ENDIF
C
C=======================================================================
C 2. TERMES SOURCES
C=======================================================================
C
C --> Initialisation
C
C
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
C
      DO IEL = 1, NCEL
        SMBRS(IEL) = 0.D0
      ENDDO
C
      ISCALA = ISCAL
C
C
      IF (IMATIS.EQ.1) THEN
C
C     Matisse
        CALL MTTSSC
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCALA ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IA(IICONR) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBRS  , ROVSDT ,
C        ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , DRTP   , DAM    ,
     &   RDEVEL , RTUSER , RA     )

      ELSE
C
C     Utilisateur standard
C
        MAXELT = MAX(NCELET,NFAC,NFABOR)
        ILS    = IDEBIA
        IDBIA1 = ILS + MAXELT
        CALL IASIZE('COVOFI',IDBIA1)
C
        CALL USTSSC
C       ===========
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCALA ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBRS  , ROVSDT ,
C        ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , DRTP   , DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C     Si on extrapole les TS :
C       SMBRS recoit -theta PROPCE du pas de temps precedent
C         (on aurait pu le faire avant ustssc, mais avec le risque que
C          l'utilisateur l'ecrase)
C       SMBRS recoit la partie du terme source qui depend de la variable
C       A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
C         on implicite le terme (donc ROVSDT*RTPA va dans SMBRS)
C       En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
C         quand meme dans SMBRS (pas d'autre choix)
      IF(ISSO2T(ISCAL).GT.0) THEN
        DO IEL = 1, NCEL
C        Stockage temporaire pour economiser un tableau
          SMBEXP = PROPCE(IEL,IPTSCA)
C        Terme source utilisateur explicite
          PROPCE(IEL,IPTSCA) = SMBRS(IEL)
C        Terme source du pas de temps precedent et
C        On suppose -ROVSDT > 0 : on implicite
C           le terme source utilisateur (le reste)
          SMBRS(IEL) = ROVSDT(IEL)*RTPA(IEL,IVAR) - THETS*SMBEXP
C        Diagonale
          ROVSDT(IEL) = - THETV*ROVSDT(IEL)
        ENDDO
C     Si on n'extrapole pas les TS :
      ELSE
        DO IEL = 1, NCEL
C        Terme source utilisateur
          SMBRS(IEL) = SMBRS(IEL) + ROVSDT(IEL)*RTPA(IEL,IVAR)
C        Diagonale
          ROVSDT(IEL) = MAX(-ROVSDT(IEL),ZERO)
        ENDDO
      ENDIF
C
C --> Physique particulieres
C     Ordre 2 non pris en compte
C
      IF ( NSCAPP.GT.0 ) THEN
        CALL PPTSSC
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCALA ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBRS  , ROVSDT , TSLAGR ,
C        ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , DRTP   , DAM    ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
C
C --> Rayonnement
C     Ordre 2 non pris en compte
C
      IF (IRAYON(IPHAS).GE.1 .AND. ISCAL.EQ.ISCALT(IPHAS)) THEN
        CALL RAYSCA
C       ===========
     &  ( ISCALT(IPHAS),IPHAS,NCELET,NCEL,
     &    SMBRS,ROVSDT,VOLUME,RA(ITSRE),RA(ITSRI))
      ENDIF
C
C     Charbon pulverise
C     Ordre 2 non pris en compte
C
      IF ( IRAYON(IPHAS).GE.1 ) THEN
        IF ( IPPMOD(ICP3PL) .GT. 0 .AND.
     &     ( ISCA(ISCAL).GE.ISCA(IH2(1)) .AND.
     &       ISCA(ISCAL).LE.ISCA(IH2(NCLACP)) ) ) THEN
C
          CALL CPRAYS
C         ===========
     &    ( IVAR  ,NCELET, NCEL  ,
     &      VOLUME,PROPCE,RA(ITSRE),RA(ITSRI),SMBRS,ROVSDT)
        ENDIF
      ENDIF
C
C --> Fuel
C     Ordre 2 non pris en compte
C     Pour l'instant rayonnement non compatible avec Fuel
C
      IF ( IRAYON(IPHAS).GE.1 .AND.
     &     IPPMOD(ICFUEL) .GE. 0  .AND. ISCAL.EQ. IHLF ) THEN
C
        WRITE(NFECRA,9001)
        CALL CSEXIT (1)
C       ===========
C
      ENDIF

C --> Lagrangien (couplage retour thermique)
C     Ordre 2 non pris en compte
C
      IF (IILAGR.EQ.2 .AND. LTSTHE.EQ.1 .AND. IPHAS.EQ.ILPHAS)  THEN
C
        IF (ISCSTH(ISCAL).EQ.2) THEN
C
C    --> Enthalpie
C
          DO IEL = 1,NCEL
            SMBRS (IEL) = SMBRS(IEL)  + TSLAGR(IEL,ITSTE)
            ROVSDT(IEL) = ROVSDT(IEL) + MAX(TSLAGR(IEL,ITSTI),ZERO)
          ENDDO
C
        ELSE IF (ISCSTH(ISCAL).EQ.1 .OR. ISCSTH(ISCAL).EQ.-1 ) THEN
C
C    --> Temperature  :
C
          IF (ICP(IPHAS).EQ.0) THEN
C
C       --> Cp constant
C
            DO IEL = 1,NCEL
              SMBRS (IEL) = SMBRS(IEL)  + TSLAGR(IEL,ITSTE)/CP0(IPHAS)
              ROVSDT(IEL) = ROVSDT(IEL) + MAX(TSLAGR(IEL,ITSTI),ZERO)
            ENDDO
C
          ELSE IF (ICP(IPHAS).GT.0) THEN
C
C       --> Cp variable
C
            IICP = IPPROC(ICP(IPHAS))
            DO IEL = 1,NCEL
              SMBRS (IEL) = SMBRS(IEL) + TSLAGR(IEL,ITSTE)
     &                                 / PROPCE(IEL,IICP)
              ROVSDT(IEL) = ROVSDT(IEL) + MAX(TSLAGR(IEL,ITSTI),ZERO)
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C     TERMES DE SOURCE DE MASSE
C
      IF (NCESMP.GT.0) THEN
C
C       Entier egal a 1 (pour navsto : nb de sur-iter)
        IIUN = 1
C
C       On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
        CALL CATSMA
C       ===========
     & ( NCELET , NCEL   , NCESMP , IIUN   , ISSO2T(ISCAL) , THETV  ,
     &   ICETSM , ITYPSM(1,IVAR) ,
     &   VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) , SMACEL(1,IPR(IPHAS)),
     &   SMBRS , ROVSDT , W1)
C
C       Si on extrapole les TS on met Gamma Pinj dans PROPCE
        IF(ISSO2T(ISCAL).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSCA) = PROPCE(IEL,IPTSCA) + W1(IEL)
          ENDDO
C       Sinon on le met directement dans SMBRS
        ELSE
          DO IEL = 1, NCEL
            SMBRS(IEL) = SMBRS(IEL) + W1(IEL)
          ENDDO
        ENDIF
C
      ENDIF
C
C
C     TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W1)
C
C     Extrapolation ou non, le terme d'accumulation de masse va dans SMBRS
      DO IEL = 1, NCEL
        SMBRS(IEL) = SMBRS(IEL)
     &              + ICONV(IVAR)*W1(IEL)*RTPA(IEL,IVAR)
      ENDDO
C
C     Extrapolation ou non le terme d'accumulation de masse a la meme forme
C       par coherence avec bilsc2
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL) - ICONV(IVAR)*W1(IEL)*THETV
      ENDDO
C
C
C     TERME INSTATIONNAIRE
C
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL)
     &           + ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
      ENDDO
C
C
C
C    SI ON CALCULE LA VARIANCE DES FLUCTUATIONS D'UN SCALAIRE,
C      ON RAJOUTE LES TERMES DE PRODUCTION ET DE DISSIPATION
C
      IF (ITSPDV.EQ.1) THEN
C
        IF(ITYTUR(IPHAS).EQ.2.OR.ITYTUR(IPHAS).EQ.3
     &       .OR.ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60) THEN
C
C Remarque : on a prevu la possibilite de scalaire associe non
C  variable de calcul, mais des adaptations sont requises
C
          IF(IVARSC.GT.0) THEN
            III = IVARSC
          ELSE
            WRITE(NFECRA,9000)IVARSC
            CALL CSEXIT(1)
          ENDIF
C
          INC = 1
          ICCOCG = 1
          NSWRGP = NSWRGR(III)
          IMLIGP = IMLIGR(III)
          IWARNP = IWARNI(III)
          EPSRGP = EPSRGR(III)
          CLIMGP = CLIMGR(III)
          EXTRAP = EXTRAG(III)
          IPHYDP = 0
C
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   III    , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,III) , COEFA(1,ICLRTP(III,ICOEF)) ,
     &                 COEFB(1,ICLRTP(III,ICOEF)) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Traitement de la production
C     On utilise MAX(PROPCE,ZERO) car en LES dynamique on fait un clipping
C     tel que (mu + mu_t)>0, donc mu_t peut etre negatif et donc
C     potentiellement (lambda/Cp + mu_t/sigma) aussi
C     Ceci ne pose probleme que quand on resout une equation de variance
C     de scalaire avec un modele LES ... ce qui serait curieux mais n'est
C     pas interdit par le code.
C       Si extrapolation : dans PROPCE
          IF (ISSO2T(ISCAL).GT.0) THEN
C         On prend la viscosite a l'instant n, meme si elle est extrapolee
            IPCVSO = IPCVST
            IF(IVIEXT(IPHAS).GT.0) IPCVSO = IPPROC(IVISTA(IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSCA) = PROPCE(IEL,IPTSCA)
     &             + 2.D0*MAX(PROPCE(IEL,IPCVSO),ZERO)
     &             *VOLUME(IEL)/SIGMAS(ISCAL)
     &               *(W1(IEL)**2 + W2(IEL)**2 + W3(IEL)**2)
            ENDDO
C       Sinon : dans SMBRS
          ELSE
            IPCVSO = IPCVST
            DO IEL = 1, NCEL
              SMBRS(IEL) = SMBRS(IEL)
     &             + 2.D0*MAX(PROPCE(IEL,IPCVSO),ZERO)
     &             *VOLUME(IEL)/SIGMAS(ISCAL)
     &             *(W1(IEL)**2 + W2(IEL)**2 + W3(IEL)**2)
            ENDDO
          ENDIF
C
C     Traitement de la dissipation
          IF (ISSO2T(ISCAL).GT.0) THEN
            THETAP = THETV
          ELSE
            THETAP = 1.D0
          ENDIF
          DO IEL = 1, NCEL
            IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50) THEN
              XK     = RTPA(IEL,IKIPH)
              XE     = RTPA(IEL,IEIPH)
            ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
              XK     =
     &        0.5D0*(RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
              XE     = RTPA(IEL,IEIPH)
            ELSEIF(ITURB(IPHAS).EQ.60) THEN
              XK     = RTPA(IEL,IKIPH)
              XE     = CMU*XK*RTPA(IEL,IOMGIP)
            ENDIF
            RHOVST = PROPCE(IEL,IPCROM)*XE/
     &                 (XK * RVARFL(ISCAL))*VOLUME(IEL)
C     La diagonale recoit eps/Rk, (*theta eventuellement)
            ROVSDT(IEL) = ROVSDT(IEL) + RHOVST*THETAP
C     SMBRS recoit la dissipation
            SMBRS(IEL) = SMBRS(IEL) - RHOVST*RTPA(IEL,IVAR)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C
      IF(ISSO2T(ISCAL).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          SMBRS(IEL) = SMBRS(IEL) + THETP1 * PROPCE(IEL,IPTSCA)
        ENDDO
      ENDIF
C
C     "VITESSE" DE DIFFUSION FACETTE
C
C     On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
C     (clipping sur (mu + mu_t)). On aurait pu prendre
C     MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
C     considere ici comme non physique.
      IF( IDIFF(IVAR).GE. 1 ) THEN
        IF(IPCVSL.EQ.0)THEN
          DO IEL = 1, NCEL
            W1(IEL) = VISLS0(ISCAL)
     &         + IDIFFT(IVAR)*MAX(PROPCE(IEL,IPCVST),ZERO)/SIGMAS(ISCAL)
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            W1(IEL) = PROPCE(IEL,IPCVSL)
     &         + IDIFFT(IVAR)*MAX(PROPCE(IEL,IPCVST),ZERO)/SIGMAS(ISCAL)
          ENDDO
        ENDIF
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
C
C=======================================================================
C 3. RESOLUTION
C=======================================================================
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
      IPP    = IPPVAR
      IWARNP = IWARNI(IVAR)
      BLENCP = BLENCV(IVAR)
      EPSILP = EPSILO(IVAR)
      EPSRSP = EPSRSM(IVAR)
      EPSRGP = EPSRGR(IVAR)
      CLIMGP = CLIMGR(IVAR)
      EXTRAP = EXTRAG(IVAR)
      RELAXP = RELAXV(IVAR)
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
     &   RELAXP , THETV  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   ROVSDT , SMBRS  , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 3. IMPRESSIONS ET CLIPPINGS
C=======================================================================
C
      IF(IVARSC.GT.0) THEN
        III = IVARSC
      ELSE
C Valeur bidon
        III = 1
      ENDIF
C
      CALL CLPSCA
C     ===========
     & ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &   PROPCE , RTP(1,III)      , RTP    )
C
C
C BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
C Ceci devrait etre valable avec le theta schema sur les Termes source
C
      IF (IWARNI(IVAR).GE.2) THEN
        IF(NSWRSM(IVAR).GT.1) THEN
          IBCL = 1
        ELSE
          IBCL = 0
        ENDIF
        DO IEL = 1, NCEL
          SMBRS(IEL) = SMBRS(IEL)
     &            - ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &                *(RTP(IEL,IVAR)-RTPA(IEL,IVAR))*IBCL
        ENDDO
        ISQRT = 1
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBRS,SMBRS,SCLNOR)
        WRITE(NFECRA,1200)CHAINE(1:8) ,SCLNOR
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,
     &'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,
     &'      ---------------------------                            ',/)
 1200 FORMAT(1X,A8,' : BILAN EXPLICITE = ',E14.5)
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS COVOFI                          ',/,
     &'@    *********                                               ',/,
     &'@    IVARSC DOIT ETRE UN ENTIER POSITIF STRICTEMENT          ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS COVOFI                          ',/,
     &'@    *********                                               ',/,
     &'@    RAYONNEMENT ET COMBUSTION FUEL NON COMPATIBLE DANS      ',/,
     &'@    LA VERSION ACTUELLE.                                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
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
