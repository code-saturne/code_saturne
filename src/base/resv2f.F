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
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE RESV2F
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM , IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   VISCF  , VISCB  , PRDV2F ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER , RA     )
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
CFONC    POUR PHI ET DE DIFFUSION POUR F_BARRE DANS LE CADRE DU
CFONC    MODELE V2F PHI-MODEL
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
CARGU ! IVAR         ! E  !  -> ! NUMERO DE VARIABLE                   !
CARGU ! ISOU         ! E  !  -> ! NUMERO DE PASSAGE                    !
CARGU ! IPP          ! E  !  -> ! NUMERO DE VARIABLE POUR SORTIES POST !
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
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! W1..10(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IPHAS  , IVAR   , ISOU   , IPP
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
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION PRDV2F(NCELET)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          INIT  , IFAC  , IEL   , INC   , ICCOCG
      INTEGER          IIUN
      INTEGER          IPRIPH,IUIPH
      INTEGER          IKIPH, IEIPH, IPHIPH, IFBIPH
      INTEGER          ICLVAR, ICLVAF
      INTEGER          ICLIKP, ICLPHI, ICLFBP
      INTEGER          IPCROM, IPCROO, IPCVIS, IPCVLO, IPCVST, IPCVSO
      INTEGER          IFLMAS, IFLMAB
      INTEGER          NSWRGP, IMLIGP, IWARNP, IPHYDP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP
      INTEGER          NITMAP, NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IPTSTA
      DOUBLE PRECISION BLENCP, EPSILP, EPSRGP, CLIMGP, EXTRAP, RELAXP
      DOUBLE PRECISION TUEXPE, THETS , THETV , THETAP, THETP1
      DOUBLE PRECISION D2S3, D1S4, D3S2
      DOUBLE PRECISION XK, XE, XNU, XROM, TTKE, TTMIN, LLKE, LLMIN
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
      IKIPH  = IK(IPHAS)
      IEIPH  = IEP(IPHAS)
      IPHIPH = IPHI(IPHAS)
      IFBIPH = IFB(IPHAS)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IUIPH))
      IFLMAB = IPPROB(IFLUMA(IUIPH))
C
      ICLIKP = ICLRTP(IKIPH,ICOEF)
      ICLPHI = ICLRTP(IPHIPH,ICOEF)
      ICLFBP = ICLRTP(IFBIPH,ICOEF)
C
      IF(ISTO2T(IPHAS).GT.0) THEN
        IPTSTA = IPPROC(ITSTUA(IPHAS))
      ELSE
        IPTSTA = 0
      ENDIF
C
      D2S3 = 2.0D0/3.0D0
      D1S4 = 1.0D0/4.0D0
      D3S2 = 3.0D0/2.0D0
C
      IF(IWARNI(IPHIPH).GE.1) THEN
        WRITE(NFECRA,1000)IPHAS
      ENDIF
C
C=======================================================================
C 2. CALCUL DU TERME EN GRAD PHI.GRAD K
C=======================================================================
C
        ICCOCG = 1
        INC = 1
        IVAR = IPHIPH
C
        NSWRGP = NSWRGR(IVAR )
        IMLIGP = IMLIGR(IVAR )
        IWARNP = IWARNI(IVAR )
        EPSRGP = EPSRGR(IVAR )
        CLIMGP = CLIMGR(IVAR )
        EXTRAP = EXTRAG(IVAR )
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHIPH , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IPHIPH ) , COEFA(1,ICLPHI) , COEFB(1,ICLPHI) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        ICCOCG = 1
        INC = 1
        IVAR = IKIPH
C
        NSWRGP = NSWRGR(IVAR )
        IMLIGP = IMLIGR(IVAR )
        IWARNP = IWARNI(IVAR )
        EPSRGP = EPSRGR(IVAR )
        CLIMGP = CLIMGR(IVAR )
        EXTRAP = EXTRAG(IVAR )
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
     &   W4     , W4     , W4     ,
     &   RTPA(1,IKIPH )  , COEFA(1,ICLIKP) , COEFB(1,ICLIKP) ,
     &   W4     , W5     , W6     ,
C        ------   ------   ------
     &   W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
        DO IEL = 1, NCEL
          W1(IEL) = W1(IEL)*W4(IEL) + W2(IEL)*W5(IEL) + W3(IEL)*W6(IEL)
        ENDDO
C
C=======================================================================
C 3. RESOLUTION DE L'EQUATION DE F_BARRE
C=======================================================================
C
      IVAR = IFBIPH
      ISOU = 3
      ICLVAR = ICLFBP
      ICLVAF = ICLFBP
      IPP    = IPPRTP(IVAR)
C
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1100) NOMVAR(IPP)
      ENDIF
C
C     S pour Source, V pour Variable
      THETS  = THETST(IPHAS)
      THETV  = THETAV(IVAR )
C
      IPCROO = IPCROM
      IPCVLO = IPCVIS
      IF(ISTO2T(IPHAS).GT.0) THEN
        IF (IROEXT(IPHAS).GT.0) THEN
          IPCROO = IPPROC(IROMA(IPHAS))
        ENDIF
        IF(IVIEXT(IPHAS).GT.0) THEN
          IPCVLO = IPPROC(IVISLA(IPHAS))
        ENDIF
      ENDIF
C
      DO IEL = 1, NCEL
        SMBR(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
C
C=======================================================================
C 3.1 TERMES SOURCES  UTILISATEURS
C=======================================================================
C
      CALL USTSV2
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL , PRDV2F , W1     ,
     &   SMBR   , ROVSDT ,
C        ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W2     , W3     , W4     , W5     , W6     , W7     ,
     &   W8     , W9     , W10    , DAM    , DRTP   ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si on extrapole les T.S.
      IF(ISTO2T(IPHAS).GT.0) THEN
        DO IEL = 1, NCEL
C       Sauvegarde pour echange
          TUEXPE = PROPCE(IEL,IPTSTA+ISOU-1)
C       Pour la suite et le pas de temps suivant
C       On met un signe "-" car on résout en fait "-div(grad fb) = ..."
          PROPCE(IEL,IPTSTA+ISOU-1) = - SMBR(IEL)
C       Second membre du pas de temps precedent
C       on implicite le terme source utilisateur (le reste)
          SMBR(IEL) = - ROVSDT(IEL)*RTPA(IEL,IVAR) - THETS*TUEXPE
C       Diagonale
          ROVSDT(IEL) = THETV*ROVSDT(IEL)
        ENDDO
      ELSE
        DO IEL = 1, NCEL
C       On met un signe "-" car on résout en fait "-div(grad fb) = ..."
C       On resout par gradient conjugue, donc on n'impose pas le signe
C          de ROVSDT
          SMBR(IEL)   = -ROVSDT(IEL)*RTPA(IEL,IVAR) - SMBR(IEL)
c          ROVSDT(IEL) =  ROVSDT(IEL)
        ENDDO
      ENDIF
C
C
C=======================================================================
C 3.2 TERME SOURCE DE F_BARRE
C     SMBR=1/L^2*(f_b + 1/T(C1-1)(phi-2/3) - C2*Pk/k/rho
C     -2*nu/k*grad_phi*grad_k -nu*div(grad(phi)) )
C  En fait on met un signe "-" car l'eq resolue est
C    -div(grad f_b) = SMBR
C=======================================================================
C
C     On calcule le terme en -VOLUME*div(grad(phi)) par itrgrp,
C     et on le stocke dans W2
C     Attention, les VISCF et VISCB calcules ici servent a ITRGRP mais
C     aussi a CODITS qui suit
C
      DO IEL = 1, NCEL
        W3(IEL) = 1.D0
      ENDDO
      CALL VISCFA
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W3     ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
C
C
      ICCOCG = 1
      INC = 1
      INIT = 1
C
      NSWRGP = NSWRGR(IPHIPH)
      IMLIGP = IMLIGR(IPHIPH)
      IWARNP = IWARNI(IPHIPH)
      EPSRGP = EPSRGR(IPHIPH)
      CLIMGP = CLIMGR(IPHIPH)
      EXTRAP = EXTRAG(IPHIPH)
      IPHYDP = 0

      CALL ITRGRP
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W2     , W2     , W2     ,
     &   RTPA(1,IPHIPH)   , COEFA(1,ICLPHI) , COEFB(1,ICLPHI) ,
     &   VISCF  , VISCB  ,
     &   W3     , W3     , W3     ,
     &   W2     ,
C        --
     &   W4     , W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C      On stocke T dans W3 et L^2 dans W4
C      Dans le cas de l'ordre 2 en temps, T est calcule en n
C      (il sera extrapole) et L^2 en n+theta (meme si k et eps restent en n)
      DO IEL=1,NCEL
        XK = RTPA(IEL,IKIPH)
        XE = RTPA(IEL,IEIPH)
        XNU  = PROPCE(IEL,IPCVLO)/PROPCE(IEL,IPCROO)
        TTKE = XK / XE
        TTMIN = CV2FCT*SQRT(XNU/XE)
        W3(IEL) = MAX(TTKE,TTMIN)
C
        XNU  = PROPCE(IEL,IPCVIS)/PROPCE(IEL,IPCROM)
        LLKE = XK**D3S2/XE
        LLMIN = CV2FET*(XNU**3/XE)**D1S4
        W4(IEL) = ( CV2FCL*MAX(LLKE,LLMIN) )**2
      ENDDO
C
C     Terme explicite, stocke temporairement dans W5
C     W2 est deja multiplie par le volume et contient deja
C     un signe "-" (issu de ITRGRP)
      DO IEL = 1, NCEL
          XROM = PROPCE(IEL,IPCROO)
          XNU  = PROPCE(IEL,IPCVLO)/XROM
          XK = RTPA(IEL,IKIPH)
          XE = RTPA(IEL,IEIPH)
          W5(IEL) = - VOLUME(IEL)*
     &         ( (CV2FC1-1.D0)*(RTPA(IEL,IPHIPH)-D2S3)/W3(IEL)
     &           -CV2FC2*PRDV2F(IEL)/XROM/XK
     &           -2.0D0*XNU/XE/W3(IEL)*W1(IEL) ) - XNU*W2(IEL)
      ENDDO
C     Si on extrapole les T.S : PROPCE
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          PROPCE(IEL,IPTSTA+ISOU-1) =
     &    PROPCE(IEL,IPTSTA+ISOU-1) + W5(IEL)
          SMBR(IEL) = SMBR(IEL) + THETP1*PROPCE(IEL,IPTSTA+ISOU-1)
        ENDDO
C     Sinon : SMBR
      ELSE
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR(IEL) + W5(IEL)
        ENDDO
      ENDIF
C
C     Terme implicite
      DO IEL = 1, NCEL
        SMBR(IEL) = ( - VOLUME(IEL)*RTPA(IEL,IFBIPH) + SMBR(IEL) )
     &       /W4(IEL)
      ENDDO
C
C ---> Matrice
C
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETAP = THETV
      ELSE
        THETAP = 1.D0
      ENDIF
      DO IEL = 1, NCEL
        ROVSDT(IEL) = (ROVSDT(IEL) + VOLUME(IEL)*THETAP)/W4(IEL)
      ENDDO
C
C
C
C=======================================================================
C 3.3 RESOLUTION EFFECTIVE DE L'EQUATION DE F_BARRE
C=======================================================================
C
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
      IWARNP = IWARNI(IVAR)
      BLENCP = BLENCV(IVAR)
      EPSILP = EPSILO(IVAR)
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
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETV  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   ROVSDT , SMBR   , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   ,
     &   W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER , RA     )
C
C
C=======================================================================
C 4. RESOLUTION DE L'EQUATION DE PHI
C=======================================================================
C
      IVAR = IPHIPH
      ISOU = 4
      ICLVAR = ICLPHI
      ICLVAF = ICLPHI
      IPP    = IPPRTP(IVAR)
C
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1100) NOMVAR(IPP)
      ENDIF
C
C     S pour Source, V pour Variable
      THETS  = THETST(IPHAS)
      THETV  = THETAV(IVAR )
C
      IPCROO = IPCROM
      IPCVSO = IPCVST
      IF(ISTO2T(IPHAS).GT.0) THEN
        IF (IROEXT(IPHAS).GT.0) THEN
          IPCROO = IPPROC(IROMA(IPHAS))
        ENDIF
        IF(IVIEXT(IPHAS).GT.0) THEN
          IPCVSO = IPPROC(IVISTA(IPHAS))
        ENDIF
      ENDIF
C
      DO IEL = 1, NCEL
        SMBR(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
C
C=======================================================================
C 4.1 TERMES SOURCES  UTILISATEURS
C=======================================================================
C
      CALL USTSV2
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL , PRDV2F , W1     ,
     &   SMBR   , ROVSDT ,
C        ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W2     , W3     , W4     , W5     , W6     , W7     ,
     &   W8     , W9     , W10    , DAM    , DRTP   ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si on extrapole les T.S.
      IF(ISTO2T(IPHAS).GT.0) THEN
        DO IEL = 1, NCEL
C       Sauvegarde pour echange
          TUEXPE = PROPCE(IEL,IPTSTA+ISOU-1)
C       Pour la suite et le pas de temps suivant
          PROPCE(IEL,IPTSTA+ISOU-1) = SMBR(IEL)
C       Second membre du pas de temps precedent
C       On suppose -ROVSDT > 0 : on implicite
C          le terme source utilisateur (le reste)
          SMBR(IEL) = ROVSDT(IEL)*RTPA(IEL,IVAR) - THETS*TUEXPE
C       Diagonale
          ROVSDT(IEL) = - THETV*ROVSDT(IEL)
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          SMBR(IEL)   = ROVSDT(IEL)*RTPA(IEL,IVAR) + SMBR(IEL)
          ROVSDT(IEL) = MAX(-ROVSDT(IEL),ZERO)
        ENDDO
      ENDIF
C
C=======================================================================
C 4.2 TERME SOURCE DE MASSE
C=======================================================================
C
C
      IF (NCESMP.GT.0) THEN
C
C       Entier egal a 1 (pour navsto : nb de sur-iter)
        IIUN = 1
C
C       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
        CALL CATSMA
C       ===========
     & ( NCELET , NCEL   , NCESMP , IIUN   , ISTO2T(IPHAS) , THETV ,
     &   ICETSM , ITYPSM(1,IVAR) ,
     &   VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) , SMACEL(1,IPRIPH) ,
     &   SMBR   ,  ROVSDT , W2 )
C
C       Si on extrapole les TS on met Gamma Pinj dans PROPCE
        IF(ISTO2T(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSTA+ISOU-1) =
     &      PROPCE(IEL,IPTSTA+ISOU-1) + W2(IEL)
          ENDDO
C       Sinon on le met directement dans SMBR
        ELSE
          DO IEL = 1, NCEL
            SMBR(IEL) = SMBR(IEL) + W2(IEL)
          ENDDO
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C 4.3 TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
C    ET TERME INSTATIONNAIRE
C=======================================================================
C
C ---> Calcul de mij
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W2)
C
C ---> Ajout au second membre
C
      DO IEL = 1, NCEL
        SMBR(IEL) = SMBR(IEL)
     &              + ICONV(IVAR)*W2(IEL)*RTPA(IEL,IVAR)
      ENDDO
C
C ---> Ajout dans la diagonale de la matrice
C     Extrapolation ou non, meme forme par coherence avec bilsc2
C
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL)
     &           + ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &           - ICONV(IVAR)*W2(IEL)*THETV
      ENDDO
C
C=======================================================================
C 4.4 TERME SOURCE DE PHI
C     SMBR=rho*f_barre - phi/k*Pk +2/k*mu_t/sigmak*grad_phi*grad_k
C=======================================================================
C
C     Terme explicite, stocke temporairement dans W2
C
      DO IEL = 1, NCEL
c   Le terme en f_barre est pris en RTP et pas en RTPA ... a priori meilleur
c  Rq : si on reste en RTP, il faut modifier le cas de l'ordre 2 (qui
c       necessite RTPA pour l'extrapolation).
        W2(IEL)   =  VOLUME(IEL)*
     &       ( PROPCE(IEL,IPCROO)*RTP(IEL,IFBIPH)
     &         +2.D0/RTPA(IEL,IKIPH)*PROPCE(IEL,IPCVSO)/SIGMAK*W1(IEL) )
      ENDDO
C
C     Si on extrapole les T.S : PROPCE
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          PROPCE(IEL,IPTSTA+ISOU-1) =
     &    PROPCE(IEL,IPTSTA+ISOU-1) + W2(IEL)
          SMBR(IEL) = SMBR(IEL) + THETP1*PROPCE(IEL,IPTSTA+ISOU-1)
        ENDDO
C     Sinon : SMBR
      ELSE
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR(IEL) + W2(IEL)
        ENDDO
      ENDIF
C
C     Terme implicite
      DO IEL = 1, NCEL
        SMBR(IEL) = SMBR(IEL)
     &       - VOLUME(IEL)*PRDV2F(IEL)*RTPA(IEL,IPHIPH)/RTPA(IEL,IKIPH)
      ENDDO
C
C ---> Matrice
C
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETAP = THETV
      ELSE
        THETAP = 1.D0
      ENDIF
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL)
     &       + VOLUME(IEL)*PRDV2F(IEL)/RTPA(IEL,IKIPH)*THETAP
      ENDDO
C

C
C=======================================================================
C 4.5 TERMES DE DIFFUSION
C=======================================================================
C ---> Viscosite
C Normalement, dans les equations du phi-model, seul la viscosite
C  turbulente intervient dans la diffusion de phi (le terme en mu
C  a disparu passant de f a f_barre). Mais tel
C  quel, cela rend le calcul instable (car mu_t tend vers 0 a la paroi
C  ce qui decouple phi de sa condition a la limite et le terme de diffusion
C  moleculaire etant integre dans f_barre, c'est comme s'il etait traite
C  en explicite).
C  -> on rajoute artificiellement de la diffusion (sachant que comme k=0 a
C  la paroi, on se moque de la valeur de phi).
C
        IF( IDIFF(IVAR).GE. 1 ) THEN
          DO IEL = 1, NCEL
            W2(IEL) = PROPCE(IEL,IPCVIS) + PROPCE(IEL,IPCVST)/SIGMAK
          ENDDO
C
          CALL VISCFA
C        ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W2     ,
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
C
C=======================================================================
C 4.6 RESOLUTION EFFECTIVE DE L'EQUATION DE PHI
C=======================================================================
C
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR(IEL) + THETP1*PROPCE(IEL,IPTSTA+ISOU-1)
        ENDDO
      ENDIF
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
      IWARNP = IWARNI(IVAR)
      BLENCP = BLENCV(IVAR)
      EPSILP = EPSILO(IVAR)
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
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETV  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , RTPA(1,IVAR)    ,
     &                     COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   ROVSDT , SMBR   , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 10. CLIPPING
C=======================================================================
C
         CALL CLPV2F
C        ===========
     & ( NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IPHAS  , IWARNI(IPHI(IPHAS)) ,
     &   PROPCE , RTP    )
C
C--------
C FORMATS
C--------
C
 1000    FORMAT(/,
     &'   ** PHASE ',I4,' RESOLUTION DU V2F (PHI ET F_BARRE)        ',/,
     &'      -----------------------------------------------        ',/)
 1100    FORMAT(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)
C
C12345678 : MAX: 12345678901234 MIN: 12345678901234 NORM: 12345678901234
C----
C FIN
C----
C
      RETURN
C
      END
c@z
