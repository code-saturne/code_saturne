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
                        SUBROUTINE TURRIJ
C                       *****************
C     -------------------------------------------------------------
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
     &   TSLAGR ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   VISCF  , VISCB  , COEFAX ,
     &   DAM    , XAM    , DRTP   ,
     &   SMBR   , ROVSDT , GRDVIT , PRODUC , GRAROX , GRAROY , GRAROZ ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS Rij-EPS 1 PHASE INCOMPRESSIBLE OU
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
CARGU ! COEFAX(NFABOR! TR !  -  ! TAB DE TRAV POUR COND.LIM. PAROI     !
CARGU !              ! TR !  -  !   ATTENTION : UNIQUEMENT AVEC ECHO   !
CARGU !              ! TR !  -  !   DE PAROI ET ABS(ICDPAR) = 1        !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR (NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! DRTP(NCELET) ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR?(NCELET)! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! GRDVIT       ! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME GRAD   !
CARGU !  (NCELET,3,3)!    !     !    DE VITESSE     UNIQT POUR ITURB=31!
CARGU ! PRODUC       ! TR !  -> ! TABLEAU DE TRAVAIL POUR PRODUCTION   !
CARGU !  (6,NCELET)  !    !     ! (SANS RHO VOLUME) UNIQT POUR ITURB=30!
CARGU ! GRAROX,Y,Z   ! TR !  -  ! TABLEAU DE TRAVAIL POUR GRAD ROM     !
CARGU !  (NCELET)    !    !     !                                      !
CARGU ! W?(NCELET)   ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR), COEFAX(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION GRDVIT(NCELET,3,3), PRODUC(6,NCELET)
      DOUBLE PRECISION GRAROX(NCELET), GRAROY(NCELET), GRAROZ(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC  , IEL   , IVAR  , ISOU  , II
      INTEGER          INC   , ICCOCG
      INTEGER          IPP   , IWARNP, ICLIP
      INTEGER          IPRIPH, IUIPH , IVIPH , IWIPH
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEIPH
      INTEGER          ICLIUP, ICLIVP, ICLIWP
      INTEGER          NSWRGP, IMLIGP, IPHYDP
      INTEGER          IPCROM, IPBROM, IPCROO, IPBROO, IIVAR
      INTEGER          IITSLA
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
C
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
      IR11IP = IR11(IPHAS)
      IR22IP = IR22(IPHAS)
      IR33IP = IR33(IPHAS)
      IR12IP = IR12(IPHAS)
      IR13IP = IR13(IPHAS)
      IR23IP = IR23(IPHAS)
      IEIPH  = IEP (IPHAS)
C
      ICLIUP = ICLRTP(IUIPH,ICOEF)
      ICLIVP = ICLRTP(IVIPH,ICOEF)
      ICLIWP = ICLRTP(IWIPH,ICOEF)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPBROM = IPPROB(IROM  (IPHAS))
C
      IF(IWARNI(IEIPH).GE.1) THEN
        IF (ITURB(IPHAS).EQ.30) THEN
          WRITE(NFECRA,1000) IPHAS
        ELSE
          WRITE(NFECRA,1001) IPHAS
        ENDIF
      ENDIF
C
C
C     SI ITURB=30 (RIJ STD) ON STOCKE DIRECTEMENT LA PRODUCTION DANS
C     LE TABLEAU PRODUC
C     SI ITURB=31 (SSG) ON STOCKE LE GRADIENT DE VITESSE DANS GRDVIT
C
C=======================================================================
C 2.a CALCUL DU TENSEUR DE PRODUCTION POUR LE RIJ STANDARD
C     W7 = P11 , W8 = P22 , W9 = P33
C     W10 = P12 , W11 = P13 , W9 = P23
C=======================================================================
C
      IF (ITURB(IPHAS).EQ.30) THEN
C INITIALISATIONS DE W7 ... W12
C
        DO II = 1 , 6
          DO IEL = 1, NCEL
            PRODUC(II,IEL) = 0.0D0
          ENDDO
        ENDDO
C
C CALCUL DU GRADIENT DES 3 COMPOSANTES DE LA VITESSE
C
        ICCOCG = 1
        INC    = 1
C
C GRADIENT SUIVANT X
C
        NSWRGP = NSWRGR(IUIPH)
        IMLIGP = IMLIGR(IUIPH)
        IWARNP = IWARNI(IUIPH)
        EPSRGP = EPSRGR(IUIPH)
        CLIMGP = CLIMGR(IUIPH)
        EXTRAP = EXTRAG(IUIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
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
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
        DO IEL = 1 , NCEL
C
          PRODUC(1,IEL) = PRODUC(1,IEL)
     &         - 2.0D0*(RTPA(IEL,IR11IP)*W1(IEL) +
     &         RTPA(IEL,IR12IP)*W2(IEL) +
     &         RTPA(IEL,IR13IP)*W3(IEL) )
C
          PRODUC(4,IEL) = PRODUC(4,IEL)
     &         - (RTPA(IEL,IR12IP)*W1(IEL) +
     &         RTPA(IEL,IR22IP)*W2(IEL) +
     &         RTPA(IEL,IR23IP)*W3(IEL) )
C
          PRODUC(5,IEL) = PRODUC(5,IEL)
     &         - (RTPA(IEL,IR13IP)*W1(IEL) +
     &         RTPA(IEL,IR23IP)*W2(IEL) +
     &         RTPA(IEL,IR33IP)*W3(IEL) )
C
        ENDDO
C
C Gradient suivant Y
C
        NSWRGP = NSWRGR(IVIPH)
        IMLIGP = IMLIGR(IVIPH)
        IWARNP = IWARNI(IVIPH)
        EPSRGP = EPSRGR(IVIPH)
        CLIMGP = CLIMGR(IVIPH)
        EXTRAP = EXTRAG(IVIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
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
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        DO IEL = 1 , NCEL
C
          PRODUC(2,IEL) = PRODUC(2,IEL)
     &         - 2.0D0*(RTPA(IEL,IR12IP)*W1(IEL) +
     &         RTPA(IEL,IR22IP)*W2(IEL) +
     &         RTPA(IEL,IR23IP)*W3(IEL) )
C
          PRODUC(4,IEL) = PRODUC(4,IEL)
     &         - (RTPA(IEL,IR11IP)*W1(IEL) +
     &         RTPA(IEL,IR12IP)*W2(IEL) +
     &         RTPA(IEL,IR13IP)*W3(IEL) )
C
          PRODUC(6,IEL) = PRODUC(6,IEL)
     &         - (RTPA(IEL,IR13IP)*W1(IEL) +
     &         RTPA(IEL,IR23IP)*W2(IEL) +
     &         RTPA(IEL,IR33IP)*W3(IEL) )
C
        ENDDO
C
C Gradient suivant Z
C
        NSWRGP = NSWRGR(IWIPH)
        IMLIGP = IMLIGR(IWIPH)
        IWARNP = IWARNI(IWIPH)
        EPSRGP = EPSRGR(IWIPH)
        CLIMGP = CLIMGR(IWIPH)
        EXTRAP = EXTRAG(IWIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
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
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        DO IEL = 1 , NCEL
C
          PRODUC(3,IEL) = PRODUC(3,IEL)
     &         - 2.0D0*(RTPA(IEL,IR13IP)*W1(IEL) +
     &         RTPA(IEL,IR23IP)*W2(IEL) +
     &         RTPA(IEL,IR33IP)*W3(IEL) )
C
          PRODUC(5,IEL) = PRODUC(5,IEL)
     &         - (RTPA(IEL,IR11IP)*W1(IEL) +
     &         RTPA(IEL,IR12IP)*W2(IEL) +
     &         RTPA(IEL,IR13IP)*W3(IEL) )
C
          PRODUC(6,IEL) = PRODUC(6,IEL)
     &         - (RTPA(IEL,IR12IP)*W1(IEL) +
     &         RTPA(IEL,IR22IP)*W2(IEL) +
     &         RTPA(IEL,IR23IP)*W3(IEL) )
C
        ENDDO
C
      ELSE
C
C=======================================================================
C 2.b CALCUL DU GRADIENT DE VITESSE POUR LE RIJ SSG
C     GRDVIT(IEL,I,J) = dUi/dxj(IEL)
C=======================================================================
C
C CALCUL DU GRADIENT DES 3 COMPOSANTES DE LA VITESSE
C
        ICCOCG = 1
        INC    = 1
C
C GRADIENT SUIVANT X
C
        NSWRGP = NSWRGR(IUIPH)
        IMLIGP = IMLIGR(IUIPH)
        IWARNP = IWARNI(IUIPH)
        EPSRGP = EPSRGR(IUIPH)
        CLIMGP = CLIMGR(IUIPH)
        EXTRAP = EXTRAG(IUIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
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
     &   GRDVIT(1,1,1)   , GRDVIT(1,1,2)   , GRDVIT(1,1,3)   ,
C        -------------     -------------     -------------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C Gradient suivant Y
C
        NSWRGP = NSWRGR(IVIPH)
        IMLIGP = IMLIGR(IVIPH)
        IWARNP = IWARNI(IVIPH)
        EPSRGP = EPSRGR(IVIPH)
        CLIMGP = CLIMGR(IVIPH)
        EXTRAP = EXTRAG(IVIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
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
     &   GRDVIT(1,2,1)   , GRDVIT(1,2,2)   , GRDVIT(1,2,3)   ,
C        -------------     -------------     -------------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C Gradient suivant Z
C
        NSWRGP = NSWRGR(IWIPH)
        IMLIGP = IMLIGR(IWIPH)
        IWARNP = IWARNI(IWIPH)
        EPSRGP = EPSRGR(IWIPH)
        CLIMGP = CLIMGR(IWIPH)
        EXTRAP = EXTRAG(IWIPH)
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
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
     &   GRDVIT(1,3,1)   , GRDVIT(1,3,2)   , GRDVIT(1,3,3)   ,
C        -------------     -------------     -------------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C=======================================================================
C 3.  CALCUL DU GRADIENT DE ROM POUR LES TERMES DE GRAVITE
C=======================================================================
C
      IF(IGRARI(IPHAS).EQ.1) THEN
C
C Conditions aux limites : Dirichlet ROMB
C   On utilise VISCB pour stocker le coefb relatif a ROM
C   On impose en Dirichlet (COEFA) la valeur ROMB
C
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
C
C Le choix ci dessous a l'avantage d'etre simple
C
        NSWRGP = NSWRGR(IR11IP)
        IMLIGP = IMLIGR(IR11IP)
        IWARNP = IWARNI(IR11IP)
        EPSRGP = EPSRGR(IR11IP)
        CLIMGP = CLIMGR(IR11IP)
        EXTRAP = EXTRAG(IR11IP)
        IPHYDP = 0
C
        IIVAR = 0
C
C     Si on extrapole les termes sources et rho, on utilise cpdt rho^n
        IPCROO = IPCROM
        IPBROO = IPBROM
        IF(ISTO2T(IPHAS).GT.0.AND.IROEXT(IPHAS).GT.0) THEN
          IPCROO = IPPROC(IROMA(IPHAS))
          IPBROO = IPPROB(IROMA(IPHAS))
        ENDIF
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IIVAR  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP ,  IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   PROPCE(1,IPCROO), PROPFB(1,IPBROO), VISCB           ,
     &   GRAROX , GRAROY , GRAROZ ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C=======================================================================
C 4.  Boucle sur les variables Rij (6 variables)
C     L'ordre est R11 R22 R33 R12 R13 R23 (La place de ces variables
C     est IR11...
C     On resout les equation dans une routine semblable a covofi.F
C=======================================================================
C
C
C
      DO ISOU = 1, 6
        IF    (ISOU.EQ.1) THEN
          IVAR   = IR11IP
        ELSEIF(ISOU.EQ.2) THEN
          IVAR   = IR22IP
        ELSEIF(ISOU.EQ.3) THEN
          IVAR   = IR33IP
        ELSEIF(ISOU.EQ.4) THEN
          IVAR   = IR12IP
        ELSEIF(ISOU.EQ.5) THEN
          IVAR   = IR13IP
        ELSEIF(ISOU.EQ.6) THEN
          IVAR   = IR23IP
        ENDIF
        IPP    = IPPRTP(IVAR)
C
        IF (IILAGR.EQ.2 .AND. IPHAS.EQ.1) THEN
          IITSLA = ITSR11 + (ISOU-1)
        ENDIF

C
C     Rij-epsilon standard (LRR)
        IF (ITURB(IPHAS).EQ.30) THEN
          CALL RESRIJ
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM(1,IVAR)  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , PRODUC , GRAROX , GRAROY , GRAROZ ,
     &   CKUPDC , SMACEL(1,IVAR)  , SMACEL(1,IPRIPH),
     &   VISCF  , VISCB  , COEFAX ,
     &   TSLAGR(1,IITSLA) , TSLAGR(1,ITSLI) ,
     &   DAM    , XAM    , DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSE
C     Rij-epsilon SSG
          CALL RESSSG
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM(1,IVAR)  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , GRDVIT , GRAROX , GRAROY , GRAROZ ,
     &   CKUPDC , SMACEL(1,IVAR)  , SMACEL(1,IPRIPH),
     &   VISCF  , VISCB  , COEFAX ,
     &   TSLAGR(1,IITSLA) , TSLAGR(1,ITSLI) ,
     &   DAM    , XAM    , DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
C
      ENDDO
C
C=======================================================================
C 5.  RESOLUTION DE EPSILON
C=======================================================================
C
         IVAR   = IEIPH
         IPP    = IPPRTP(IVAR)
         ISOU   = 7
C
         CALL RESEPS
C        ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM(1,IVAR)  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , GRDVIT , PRODUC ,GRAROX , GRAROY , GRAROZ ,
     &   CKUPDC , SMACEL(1,IVAR)  , SMACEL(1,IPRIPH),
     &   VISCF  , VISCB  ,
     &   TSLAGR ,
     &   DAM    , XAM    , DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 6. CLIPPING
C=======================================================================
C
         ICLIP  = 2
         CALL CLPRIJ
C        ===========
     & ( NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IPHAS  , ICLIP  ,
     &   PROPCE , RTPA   , RTP    )
C
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU Rij-EPSILON LRR             ',/,
     &'      ------------------------------------------             ',/)
 1001 FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU Rij-EPSILON SSG             ',/,
     &'      ------------------------------------------             ',/)
C
#else
C
 1000 FORMAT(/,
     &'   ** PHASE ',I4,' SOLVING Rij-EPSILON LRR'                   ,/,
     &'      ------------------------------------'                   ,/)
 1001 FORMAT(/,
     &'   ** PHASE ',I4,' SOLVING Rij-EPSILON SSG'                   ,/,
     &'      ------------------------------------'                   ,/)
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
