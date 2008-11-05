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
                        SUBROUTINE CFENER
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
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
CFONC    POUR L'ENERGIE TOTALE SUR UN PAS DE TEMPS
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
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cfpoin.h"
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
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
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
      INTEGER          IFINIA, IFINRA
      INTEGER          IVAR  , IPHAS
      INTEGER          IFAC  , IEL
      INTEGER          INIT  , ISQRT , III
      INTEGER          ICLVAR, ICLVAF
      INTEGER          IPCROM, IPCVST, IPCVSL, IFLMAS, IFLMAB
      INTEGER          IPPVAR, IPP
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP, NITMAP
      INTEGER          NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, BLENCP, EPSILP
      DOUBLE PRECISION SCLNOR, THETAP, EPSRSP
C
      INTEGER          IWB    , INC    , ICCOCG , ICOEFA , ICOEFB
      INTEGER          IVAR0  , IPHYDP , IIJ , II , JJ
      INTEGER          ICCFTH , IMODIF
      INTEGER          IEL1  , IEL2, IIFRU, IIFBE
      INTEGER          IDIMTE, ITENSO , ITERNS
      INTEGER          MAXELT, ILS, IDBIA1
      DOUBLE PRECISION FLUX
      DOUBLE PRECISION DIJPFX, DIJPFY, DIJPFZ, POND  , PIP   , PJP
      DOUBLE PRECISION DIIPFX, DIIPFY, DIIPFZ, DJJPFX, DJJPFY, DJJPFZ
c      DOUBLE PRECISION FLUI  , FLUJ
C
C***********************************************************************
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
C --- Indicateur flux de bord Rusanov
      IF(IIFBRU.GT.0) THEN
        IIFRU = IIFBRU+(IPHAS-1)*NFABOR
      ELSE
        IIFRU = 1
      ENDIF
C
C --- Indicateur flux conductif de bord imposé
      IF(IIFBET.GT.0) THEN
        IIFBE = IIFBET+(IPHAS-1)*NFABOR
      ELSE
        IIFBE = 1
      ENDIF
C
C --- Impressions
      CHAINE = NOMVAR(IPPVAR)
C
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1000) CHAINE(1:8)
      ENDIF
C
C --- Reservation de la memoire
C
      CALL MEMCFE
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IWB    ,
     &   IFINIA , IFINRA )
C
      IDEBIA = IFINIA
      IDEBRA = IFINRA
C
C=======================================================================
C 2. TERMES SOURCES
C=======================================================================
C
C --> Theta-schema de resolution
C
C Pour l'instant on prend THETA=1 et on ne code pas le theta-schema
C
C --> Initialisation
C
      DO IEL = 1, NCEL
        SMBRS(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
C
C
C     TERME SOURCE VOLUMIQUE DE CHALEUR : RHO*PHI *VOLUME
C     =================================          v
C
      MAXELT = MAX(NCELET,NFAC,NFABOR)
      ILS = IDEBIA
      IDBIA1 = ILS + MAXELT
      CALL IASIZE('CFENER',IDBIA1)
C
      CALL USTSSC
C     ===========
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
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
      DO IEL = 1, NCEL
        SMBRS(IEL) = SMBRS(IEL) + ROVSDT(IEL)*RTP(IEL,IVAR)
        ROVSDT(IEL) = MAX(-ROVSDT(IEL),ZERO)
      ENDDO
C
C
C     TERMES DE SOURCE DE MASSE
C     =========================
C
C     GAMMA(IEL) = SMACEL(IEL,IPR(IPHAS))
C
C     Terme implicite : GAMMA*VOLUME
C                                                        n
C     Terme explicite : GAMMA*VOLUME*e   - GAMMA*VOLUME*e
C                                     inj
      IF (NCESMP.GT.0) THEN
        ITERNS = 1
        CALL CATSMA ( NCELET , NCEL , NCESMP , ITERNS ,
     &                ISNO2T(IPHAS), THETAV(IVAR),
     &                ICETSM , ITYPSM(1,IVAR) ,
     &                VOLUME , RTPA(1,IVAR) , SMACEL(1,IVAR) ,
     &                SMACEL(1,IPR(IPHAS)) , SMBRS , ROVSDT , W1)
      ENDIF
C
C                                     __    n+1
C     TERME D'ACCUMULATION DE MASSE : >  (Q    .n)  *S
C     =============================   --   pr     ij  ij
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W1)
C
C                                      __    n+1             n
C     TERME INSTATIONNAIRE EXPLICITE : >  (Q    .n)  *S   * e
C     ==============================   --   pr     ij  ij
C
      DO IEL = 1, NCEL
        SMBRS(IEL) = SMBRS(IEL) + ICONV(IVAR)*W1(IEL)*RTPA(IEL,IVAR)
      ENDDO
C
C                                      RHO*VOLUME   __    n+1
C     TERME INSTATIONNAIRE IMPLICITE : ---------- - >  (Q    .n)  *S
C     ==============================       DT       --   pr     ij  ij
C
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL) - ICONV(IVAR)*W1(IEL)
     &           + ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
      ENDDO
C
C                                       __        v
C     TERME DE DISSIPATION VISQUEUSE  : >  ((SIGMA *U).n)  *S
C     ==============================    --               ij  ij
C
      IF( IDIFF(IU(IPHAS)).GE. 1 ) THEN
C                             ^^^
        CALL CFDIVS
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBRS  , RTP(1,IU(IPHAS)), RTP(1,IV(IPHAS)), RTP(1,IW(IPHAS)),
C        ------
     &   W9     ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C                                         __   P        n+1
C     TERME DE TRANSPORT DE PRESSION  : - >  (---)  *(Q    .n)  *S
C     ==============================      --  RHO ij   pr     ij  ij
C
C
      IF(IGRDPP(IPHAS).GT.0) THEN
        DO IEL = 1, NCEL
          W9(IEL) = RTP(IEL,ISCA(IRHO(IPHAS)))
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          W9(IEL) = RTPA(IEL,ISCA(IRHO(IPHAS)))
        ENDDO
      ENDIF
C
C     Avec Reconstruction : ca pose probleme pour l'instant
C
C
C   Calcul du gradient de P/RHO
C
c      DO IEL = 1, NCEL
c        W7(IEL) = RTP(IEL,IPR(IPHAS))/W9(IEL)
c      ENDDO
C
C Rq : A defaut de connaitre les parametres pour P/RHO on prend ceux de P
C
c      III = IPR(IPHAS)
c      INC = 1
c      ICCOCG = 1
c      NSWRGP = NSWRGR(III)
c      IMLIGP = IMLIGR(III)
c      IWARNP = IWARNI(III)
c      EPSRGP = EPSRGR(III)
c      CLIMGP = CLIMGR(III)
c      EXTRAP = EXTRAG(III)
C
C       On alloue localement 2 tableaux de NFABOR pour le calcul
C       de COEFA et COEFB de P/RHO
C
c      ICOEFA = IDEBRA
c      ICOEFB = ICOEFA + NFABOR
c      IFINRA = ICOEFB + NFABOR
c      CALL RASIZE ('CFENER',IFINRA)
C     ===========
C
c      DO IFAC = 1, NFABOR
c        RA(ICOEFA+IFAC-1) = ZERO
c        RA(ICOEFB+IFAC-1) = 1.D0
c      ENDDO
C
C En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
c      IF(IRANGP.GE.0) THEN
c        CALL PARCOM(W7)
C       ===========
c      ENDIF
C
C    Periodique
c      IF(IPERIO.EQ.1) THEN
c        IDIMTE = 0
c        ITENSO = 0
c        CALL PERCOM
C       ===========
c     &( IDIMTE , ITENSO ,
c     &  W7     , W7     , W7     ,
c     &  W7     , W7     , W7     ,
c     &  W7     , W7     , W7     )
c      ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
c      IVAR0 = 0
c      IPHYDP = 0
c      CALL GRDCEL
C     ===========
c     & ( IDEBIA , IFINRA ,
c     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
c     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
c     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
c     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
c     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
c     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
c     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
c     &   IDEVEL , ITUSER , IA     ,
c     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
c     &   W7     , W7     , W7     ,
c     &   W7     , RA(ICOEFA) , RA(ICOEFB)  ,
c     &   W1     , W2     , W3     ,
c     &   W4     , W5     , W6     ,
c     &   RDEVEL , RTUSER , RA     )
C
C       On libere la place dans RA
C
c      IFINRA = IDEBRA
C
C     Faces internes
c      DO IFAC = 1, NFAC
C
c        II = IFACEL(1,IFAC)
c        JJ = IFACEL(2,IFAC)
C
c        IIJ = IDIJPF-1+3*(IFAC-1)
c        DIJPFX = RA(IIJ+1)
c        DIJPFY = RA(IIJ+2)
c        DIJPFZ = RA(IIJ+3)
C
c        POND   = RA(IPOND-1+IFAC)
C
C        Calcul II' et JJ'
C
c        DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
c     &           (1.D0-POND) * DIJPFX)
c        DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
c     &           (1.D0-POND) * DIJPFY)
c        DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
c     &           (1.D0-POND) * DIJPFZ)
c        DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
c     &               POND  * DIJPFX
c        DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
c     &               POND  * DIJPFY
c        DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
c     &               POND  * DIJPFZ
C
c        PIP = W7(II)
c     &       +W1(II)*DIIPFX+W2(II)*DIIPFY+W3(II)*DIIPFZ
C
c        PJP = W7(JJ)
c     &       +W1(JJ)*DJJPFX+W2(JJ)*DJJPFY+W3(JJ)*DJJPFZ
C
c        FLUI = (PROPFA(IFAC,IFLMAS)+ABS(PROPFA(IFAC,IFLMAS)))
c        FLUJ = (PROPFA(IFAC,IFLMAS)-ABS(PROPFA(IFAC,IFLMAS)))
C
c        VISCF(IFAC) = -(POND*PIP*FLUI+POND*PJP*FLUJ)
C
c      ENDDO
C
C     Sans Reconstruction
C
C     En periodique et parallele, echange avant utilisation
C       des valeurs aux faces
C
C    Parallele
      IF(IRANGP.GE.0) THEN
        CALL PARCOM(W9)
C       ===========
      ENDIF
C
C    Periodique
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &( IDIMTE , ITENSO ,
     &  W9     , W9     , W9     ,
     &  W9     , W9     , W9     ,
     &  W9     , W9     , W9     )
      ENDIF
C
C     Faces internes
      DO IFAC = 1, NFAC
        IEL1 = IFACEL(1,IFAC)
        IEL2 = IFACEL(2,IFAC)
        VISCF(IFAC) =
     &     - RTP(IEL1,IPR(IPHAS))/W9(IEL1)
     &     *0.5D0*( PROPFA(IFAC,IFLMAS) +ABS(PROPFA(IFAC,IFLMAS)) )
     &     - RTP(IEL2,IPR(IPHAS))/W9(IEL2)
     &     *0.5D0*( PROPFA(IFAC,IFLMAS) -ABS(PROPFA(IFAC,IFLMAS)) )
      ENDDO
C
C     Faces de bord : pour les faces ou on a calcule un flux de Rusanov,
C       on remplace la contribution standard par le flux de Rusanov qui
C       contient tous les flux convectifs (et il faudra donc eliminer le
C       flux convectif dans cfbsc2)
C
      IF(IIFBRU.GT.0) THEN
C
        DO IFAC = 1, NFABOR
          IF(IA(IIFRU+IFAC-1).EQ.0) THEN
C
            IEL = IFABOR(IFAC)
            VISCB(IFAC) = - PROPFB(IFAC,IFLMAB)
     &    * ( COEFA(IFAC,ICLRTP(IPR(IPHAS),ICOEF))
     &      + COEFB(IFAC,ICLRTP(IPR(IPHAS),ICOEF))*RTP(IEL,IPR(IPHAS)) )
     &    / ( COEFA(IFAC,ICLRTP(ISCA(IRHO(IPHAS)),ICOEF))
     &      + COEFB(IFAC,ICLRTP(ISCA(IRHO(IPHAS)),ICOEF))*W9(IEL) )
C
          ELSE
            VISCB(IFAC) = - PROPFB(IFAC,IPPROB(IFBENE(IPHAS)))
          ENDIF
        ENDDO
C
      ELSE
        DO IFAC = 1, NFABOR
C
          IEL = IFABOR(IFAC)
          VISCB(IFAC) = - PROPFB(IFAC,IFLMAB)
     &    * ( COEFA(IFAC,ICLRTP(IPR(IPHAS),ICOEF))
     &      + COEFB(IFAC,ICLRTP(IPR(IPHAS),ICOEF))*RTP(IEL,IPR(IPHAS)) )
     &    / ( COEFA(IFAC,ICLRTP(ISCA(IRHO(IPHAS)),ICOEF))
     &      + COEFB(IFAC,ICLRTP(ISCA(IRHO(IPHAS)),ICOEF))*W9(IEL) )
C
        ENDDO
      ENDIF
C
C     Divergence
      INIT = 0
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &              IFACEL,IFABOR,VISCF,VISCB,SMBRS)
C
C
C     TERME DE FORCES DE PESANTEUR : RHO*g.U *VOLUME
C     ============================
C
      DO IEL = 1, NCEL
        SMBRS(IEL) = SMBRS(IEL) + W9(IEL)*VOLUME(IEL)
     &                           *( GX*RTP(IEL,IU(IPHAS))
     &                            + GY*RTP(IEL,IV(IPHAS))
     &                            + GZ*RTP(IEL,IW(IPHAS)) )
      ENDDO
C
C                                       Kij*Sij           LAMBDA   Cp   MUT
C     "VITESSE" DE DIFFUSION FACETTE : --------- avec K = ------ + -- .------
C     ==============================    IJ.nij              Cv     Cv  SIGMAS
C
      IF( IDIFF(IVAR).GE. 1 ) THEN
C
C     MUT/SIGMAS
        DO IEL = 1, NCEL
          W1(IEL) = PROPCE(IEL,IPCVST)/SIGMAS(ISCAL)
        ENDDO
C     CP*MUT/SIGMAS
        IF(ICP(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL)*PROPCE(IEL,IPPROC(ICP(IPHAS)))
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL)*CP0(IPHAS)
          ENDDO
        ENDIF
C     (CP/CV)*MUT/SIGMAS
        IF(ICV(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL)/PROPCE(IEL,IPPROC(ICV(IPHAS)))
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL)/CV0(IPHAS)
          ENDDO
        ENDIF
C     (CP/CV)*MUT/SIGMAS+LAMBDA/CV
        IF(IPCVSL.EQ.0)THEN
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL) + VISLS0(ISCAL)
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL) + PROPCE(IEL,IPCVSL)
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
C
C     TERME DIFFUSIF COMPLEMENTAIRE : - div( K grad ( epsilon - Cv.T ) )
C     =============================                   1  2
C                                     - div( K grad ( -.u  ) )
C                                                     2
C
C     Terme complementaire au centre des cellules
        ICCFTH = 7
        IMODIF = 0
        CALL USCFTH
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W9     , RA(IWB), W8     , W1     ,
C        ------   -------
     &   RDEVEL , RTUSER , RA     )
C
C     Calcul de la divergence avec reconstruction
C
C
C   Calcul du gradient de (0.5*u*u+EPSILONsup)
C
C
      DO IEL = 1, NCEL
        W7(IEL) =0.5D0*( RTP(IEL,IU(IPHAS))**2
     &                  +RTP(IEL,IV(IPHAS))**2
     &                  +RTP(IEL,IW(IPHAS))**2 ) + W9(IEL)
      ENDDO
C
C Rq : A defaut de connaitre les parametres, on prend ceux de la Vitesse
C
      III = IU(IPHAS)
      INC = 1
      ICCOCG = 1
      NSWRGP = NSWRGR(III)
      IMLIGP = IMLIGR(III)
      IWARNP = IWARNI(III)
      EPSRGP = EPSRGR(III)
      CLIMGP = CLIMGR(III)
      EXTRAP = EXTRAG(III)
C
C       On alloue localement 2 tableaux de NFABOR pour le calcul
C       de COEFA et COEFB
C
      ICOEFA = IDEBRA
      ICOEFB = ICOEFA + NFABOR
      IFINRA = ICOEFB + NFABOR
      CALL RASIZE ('CFENER',IFINRA)
C     ===========
C
      DO IFAC = 1, NFABOR
        RA(ICOEFA+IFAC-1) = ZERO
        RA(ICOEFB+IFAC-1) = 1.D0
      ENDDO
C
C En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
      IF(IRANGP.GE.0) THEN
        CALL PARCOM(W7)
C       ===========
      ENDIF
C
C    Periodique
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &( IDIMTE , ITENSO ,
     &  W7     , W7     , W7     ,
     &  W7     , W7     , W7     ,
     &  W7     , W7     , W7     )
      ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
      IVAR0 = 0
      IPHYDP = 0
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W7     , W7     , W7     ,
     &   W7     , RA(ICOEFA) , RA(ICOEFB)  ,
     &   W1     , W2     , W3     ,
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C       On libere la place dans RA
C
      IFINRA = IDEBRA
C
C     Faces internes
C
      DO IFAC = 1, NFAC
C
        II = IFACEL(1,IFAC)
        JJ = IFACEL(2,IFAC)
C
        IIJ = IDIJPF-1+3*(IFAC-1)
        DIJPFX = RA(IIJ+1)
        DIJPFY = RA(IIJ+2)
        DIJPFZ = RA(IIJ+3)
C
        POND   = RA(IPOND-1+IFAC)
C
C        Calcul II' et JJ'
C
        DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &           (1.D0-POND) * DIJPFX)
        DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &           (1.D0-POND) * DIJPFY)
        DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &           (1.D0-POND) * DIJPFZ)
        DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &               POND  * DIJPFX
        DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &               POND  * DIJPFY
        DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &               POND  * DIJPFZ
C
        PIP = W7(II)
     &       +W1(II)*DIIPFX+W2(II)*DIIPFY+W3(II)*DIIPFZ
C
        PJP = W7(JJ)
     &       +W1(JJ)*DJJPFX+W2(JJ)*DJJPFY+W3(JJ)*DJJPFZ
C
        FLUX = VISCF(IFAC)*(PIP-PJP)
C
        SMBRS(II) = SMBRS(II) - FLUX
        SMBRS(JJ) = SMBRS(JJ) + FLUX
C
      ENDDO
C
C
C       Assemblage a partir des facettes de bord
C     Pour les faces à flux imposé ou temperature imposée, tout est
C       pris par le terme de diffusion de l'energie. On ne doit donc
C       pas prendre en compte la contribution des termes en u2 et e-CvT
C       quand IA(IIFBE+IFAC-1).NE.0
C
        IF(IIFBET.GT.0) THEN
          DO IFAC = 1, NFABOR
C
            IF(IA(IIFBE+IFAC-1).EQ.0) THEN
              IEL = IFABOR(IFAC)
C
              FLUX = VISCB(IFAC)*( W9(IEL) - RA(IWB +IFAC-1)
     &           + 0.5D0*( RTP(IEL,IU(IPHAS))**2 -
     &   ( COEFA(IFAC,ICLRTP(IU(IPHAS),ICOEF))
     &   + COEFB(IFAC,ICLRTP(IU(IPHAS),ICOEF))*RTP(IEL,IU(IPHAS)) )**2
     &                   + RTP(IEL,IV(IPHAS))**2 -
     &   ( COEFA(IFAC,ICLRTP(IV(IPHAS),ICOEF))
     &   + COEFB(IFAC,ICLRTP(IV(IPHAS),ICOEF))*RTP(IEL,IV(IPHAS)) )**2
     &                   + RTP(IEL,IW(IPHAS))**2 -
     &   ( COEFA(IFAC,ICLRTP(IW(IPHAS),ICOEF))
     &   + COEFB(IFAC,ICLRTP(IW(IPHAS),ICOEF))*RTP(IEL,IW(IPHAS)) )**2))
C
              SMBRS(IEL) = SMBRS(IEL) - FLUX
C
            ENDIF
C
          ENDDO
C
C     Sinon : meme code, mais sans le test
        ELSE
C
          DO IFAC = 1, NFABOR
C
            IEL = IFABOR(IFAC)
C
            FLUX = VISCB(IFAC)*( W9(IEL) - RA(IWB +IFAC-1)
     &           + 0.5D0*( RTP(IEL,IU(IPHAS))**2 -
     &   ( COEFA(IFAC,ICLRTP(IU(IPHAS),ICOEF))
     &   + COEFB(IFAC,ICLRTP(IU(IPHAS),ICOEF))*RTP(IEL,IU(IPHAS)) )**2
     &                   + RTP(IEL,IV(IPHAS))**2 -
     &   ( COEFA(IFAC,ICLRTP(IV(IPHAS),ICOEF))
     &   + COEFB(IFAC,ICLRTP(IV(IPHAS),ICOEF))*RTP(IEL,IV(IPHAS)) )**2
     &                   + RTP(IEL,IW(IPHAS))**2 -
     &   ( COEFA(IFAC,ICLRTP(IW(IPHAS),ICOEF))
     &   + COEFB(IFAC,ICLRTP(IW(IPHAS),ICOEF))*RTP(IEL,IW(IPHAS)) )**2))
C
            SMBRS(IEL) = SMBRS(IEL) - FLUX
C
          ENDDO
        ENDIF
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
C=======================================================================
C 4. RESOLUTION
C=======================================================================
C
      ICONVP = ICONV (IVAR)
      IDIFFP = IDIFF (IVAR)
      IRESLP = IRESOL(IVAR)
      NITMAP = NITMAX(IVAR)
      NDIRCP = NDIRCL(IVAR)
      NSWRSP = NSWRSM(IVAR)
      NSWRGP = NSWRGR(IVAR)
      IMLIGP = IMLIGR(IVAR)
      IRCFLP = IRCFLU(IVAR)
      ISCHCP = ISCHCV(IVAR)
      ISSTPP = ISSTPC(IVAR)
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
      THETAP = THETAV(IVAR)
      IESCAP = 0
C
      CALL CFCDTS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISCAL  , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP , IIFBRU ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRSP , EPSRGP , CLIMGP , EXTRAP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , IA(IIFRU) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
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
C 5. IMPRESSIONS ET CLIPPINGS
C=======================================================================
C
C Valeur bidon
        III = 1
C
      CALL CLPSCA
C     ===========
     & ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &   PROPCE , RTP(1,III)      , RTP    )
C
C --- Traitement utilisateur pour gestion plus fine des bornes
C       et actions correctives éventuelles.
        ICCFTH = -4
        IMODIF = 0
        CALL USCFTH
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C --- Bilan explicite (voir codits : on enleve l'increment)
C
      IF (IWARNI(IVAR).GE.2) THEN
        DO IEL = 1, NCEL
          SMBRS(IEL) = SMBRS(IEL)
     &            - ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &                *(RTP(IEL,IVAR)-RTPA(IEL,IVAR))
     &                * MAX(0,MIN(NSWRSM(IVAR)-2,1))
        ENDDO
        ISQRT = 1
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBRS,SMBRS,SCLNOR)
        WRITE(NFECRA,1200)CHAINE(1:8) ,SCLNOR
      ENDIF
C
C=======================================================================
C 6. ACTUALISATION FINALE DE LA PRESSION (et calcul de la température)
C=======================================================================
C                               n+1      n+1  n+1
C On utilise l'equation d'etat P   =P(RHO   ,H   )
C
C --- Calcul de P et T au centre des cellules
        ICCFTH = 24
        IMODIF = 0
        CALL USCFTH
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RTP(1,IPR(IPHAS)) , RTP(1,ISCA(ITEMPK(IPHAS))) , W8     , W9 ,
C        -----------------   --------------------------
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 7. COMMUNICATION DE LA PRESSION, DE L'ENERGIE ET DE LA TEMPERATURE
C=======================================================================
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (RTP(1,IPR(IPHAS)))
C       ===========
        CALL PARCOM (RTP(1,IVAR))
C       ===========
        CALL PARCOM (RTP(1,ISCA(ITEMPK(IPHAS))))
C       ===========
      ENDIF
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     & ( IDIMTE , ITENSO ,
     &   RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)),
     &   RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)),
     &   RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)))
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     & ( IDIMTE , ITENSO ,
     &   RTP(1,IVAR)     , RTP(1,IVAR)     , RTP(1,IVAR),
     &   RTP(1,IVAR)     , RTP(1,IVAR)     , RTP(1,IVAR),
     &   RTP(1,IVAR)     , RTP(1,IVAR)     , RTP(1,IVAR) )
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     & ( IDIMTE , ITENSO ,
     &   RTP(1,ISCA(ITEMPK(IPHAS))) , RTP(1,ISCA(ITEMPK(IPHAS))) ,
     &   RTP(1,ISCA(ITEMPK(IPHAS))) , RTP(1,ISCA(ITEMPK(IPHAS))) ,
     &   RTP(1,ISCA(ITEMPK(IPHAS))) , RTP(1,ISCA(ITEMPK(IPHAS))) ,
     &   RTP(1,ISCA(ITEMPK(IPHAS))) , RTP(1,ISCA(ITEMPK(IPHAS))) ,
     &   RTP(1,ISCA(ITEMPK(IPHAS))) )
      ENDIF
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,
     &'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,
     &'      ---------------------------                            ',/)
 1200 FORMAT(1X,A8,' : BILAN EXPLICITE = ',E14.5)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
