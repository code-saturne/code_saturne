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
                        SUBROUTINE DTTVAR
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IWARNP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   VISCF  , VISCB  , DAM    , COFBDT , W1     , W2     , W3     ,
     &   COEFBR , GRAROX , GRAROY , GRAROZ , WCF    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL DU PAS DE TEMPS LOCAL
CFONC AFFICHAGE DES NOMBRES DE COURANT + FOURIER MINIMUM, MAXIMUM
CFONC  On dispose des types de faces de bord au pas de temps
CFONC    precedent (sauf au premier pas de temps, ou les tableaux
CFONC    ITYPFB et ITRIFB n'ont pas ete renseignes)
CFONC
CFONC Sous programme utilise dans le cas une seule phase (ou
CFONC  si seule la phase 1 pilote le pas de temps)
c@fonce
C-----------------------------------------------------------------------
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
CARGU ! NCEPDP       ! E  !  -> ! NOMBRE DE CELLULES AVEC PDC          !
CARGU ! NCKPDP       ! E  !  -> ! NBR DE COEF DU TENSEUR DE PDC (3 OU 6!
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES AVEC TSM          !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
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
CARGU ! DT(NCELET)   ! TR ! <-> ! PAS DE TEMPS                         !
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
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,NVAR)!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     ! POUR IVAR=IPR, SMACEL=FLUX DE MASSE  !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! COFBDT(NFABOR! TR !  -  ! CONDITION LIMITE PAS DE TEMPS        !
CARGU ! W1,2,3(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU ! GRARO.(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL (IPTLRO=1)       !
CARGU ! COEFBR(NFABOR! TR !  -  ! TABLEAU DE TRAVAIL (IPTLRO=1)        !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , IWARNP
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
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET ), COFBDT(NFABOR)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
C   Attention, COEFBR n'est defini  que pour IPTLRO = 1
      DOUBLE PRECISION COEFBR(NFABOR)
C   Attention, GRAROX, GRAROY, GRAROZ ne sont
C   definis que pour IPTLRO = 1 ou en compressible
      DOUBLE PRECISION GRAROX(NCELET),GRAROY(NCELET),GRAROZ(NCELET)
C   Attention, WCF n'est defini qu'en compressible
      DOUBLE PRECISION WCF(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*8      CNOM
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IEL, ICFMAX, ICFMIN, IDIFF0, ICONV0, ISYM
      INTEGER          MODNTL
      INTEGER          IPHAS, IUIPH, IPCVIS, IPCVST
      INTEGER          IFLMAS, IFLMAB
      INTEGER          ICOU, IFOU , ICOUCF
      INTEGER          INC, ICCOCG
      INTEGER          NSWRGP, IMLIGP , IPHYDP
      INTEGER          IPCROM, IPBROM, IIVAR
      INTEGER          NBRVAL
      INTEGER          IPCCOU, IPCFOU
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
      DOUBLE PRECISION CFMAX,CFMIN, COUFOU, W1MIN, W2MIN, W3MIN
      DOUBLE PRECISION UNPVDT, ROM
      DOUBLE PRECISION XYZMAX(3), XYZMIN(3)
      DOUBLE PRECISION DTSDTM,DTSDT0
C
C***********************************************************************
C
C=======================================================================
C 0.  INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPHAS   = 1
      IUIPH   = IU(IPHAS)
      IFLMAS  = IPPROF(IFLUMA(IUIPH))
      IFLMAB  = IPPROB(IFLUMA(IUIPH))
      IPCVIS  = IPPROC(IVISCL(IPHAS))
      IPCVST  = IPPROC(IVISCT(IPHAS))
      IPCROM  = IPPROC(IROM  (IPHAS))
      IPBROM  = IPPROB(IROM  (IPHAS))
      IPCCOU  = IPPROC(ICOUR (IPHAS))
      IPCFOU  = IPPROC(IFOUR (IPHAS))
C
      IF(NTLIST.GT.0) THEN
        MODNTL = MOD(NTCABS,NTLIST)
      ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
        MODNTL = 0
      ELSE
        MODNTL = 1
      ENDIF
C
      IF (
     &   .NOT. ( ICONV(IUIPH).GE.1.AND.
     &           (IWARNP.GE.2.OR.MODNTL.EQ.0) ) .AND.
     &   .NOT. ( IDIFF(IUIPH).GE.1.AND.
     &           (IWARNP.GE.2.OR.MODNTL.EQ.0) ) .AND.
     &   .NOT. ( IPPMOD(ICOMPF).GE.0.AND.
     &           (IWARNP.GE.2.OR.MODNTL.EQ.0) ) .AND.
     &   .NOT. ( IDTVAR.EQ.1.OR.IDTVAR.EQ.2.OR.
     &           ( (IWARNP.GE.2.OR.MODNTL.EQ.0).AND.
     &             (IDIFF(IUIPH).GE.1.OR.ICONV(IUIPH).GE.1
     &                               .OR.IPPMOD(ICOMPF).GE.0)  ) )
     &   ) THEN
C
        RETURN
C
      ENDIF
C
C=======================================================================
C 1.  CONDITION LIMITE POUR MATRDT
C=======================================================================
C

      DO IFAC = 1, NFABOR
        IF(PROPFB(IFAC,IFLMAB).LT.0.D0) THEN
          COFBDT(IFAC) = 0.D0
        ELSE
          COFBDT(IFAC) = 1.D0
        ENDIF
      ENDDO
C
C=======================================================================
C 2.  CALCUL DE LA LIMITATION EN COMPRESSIBLE
C=======================================================================
C
C     On commence par cela afin de disposer de VISCF VISCB comme
C       tableaux de travail.
C
        IF(IPPMOD(ICOMPF).GE.0) THEN
C
          CALL CFDTTV
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IWARNP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   WCF    ,
C        ---
     &   VISCF  , VISCB  , COFBDT , W1     , W2     , DAM    ,
     &   GRAROX , GRAROY , GRAROZ ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF

C
C=======================================================================
C 3.  CALCUL DE LA VISCOSITE FACETTES
C=======================================================================
C
C
C     On s'en sert dans les divers matrdt suivants
C
C     "VITESSE" DE DIFFUSION FACETTE
C
      IF( IDIFF(IUIPH).GE. 1 ) THEN
        DO IEL = 1, NCEL
          W1    (IEL) = PROPCE(IEL,IPCVIS)
     &                                +IDIFFT(IUIPH)*PROPCE(IEL,IPCVST)
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
        DO IFAC = 1, NFAC
          VISCF(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
      ENDIF
C
C=======================================================================
C 4.  ALGORITHME INSTATIONNAIRE
C=======================================================================
C
      IF (IDTVAR.GE.0) THEN
C
C=======================================================================
C 4.1  PAS DE TEMPS VARIABLE A PARTIR DE COURANT ET FOURIER IMPOSES
C=======================================================================
C
C     On calcule le pas de temps thermique max (meme en IDTVAR=0, pour affichage)
C     DTTMAX = 1/SQRT(MAX(0+,gradRO.g/RO) -> W3
C
        IF (IPTLRO.EQ.1) THEN
C
          DO IFAC = 1, NFABOR
            COEFBR(IFAC) = 0.D0
          ENDDO
C
          NSWRGP = NSWRGR(IPR(IPHAS))
          IMLIGP = IMLIGR(IPR(IPHAS))
          IWARNP = IWARNI(IPR(IPHAS))
          EPSRGP = EPSRGR(IPR(IPHAS))
          CLIMGP = CLIMGR(IPR(IPHAS))
          EXTRAP = 0.D0
          IPHYDP = 0
C
          IIVAR = 0
          INC   = 1
          ICCOCG = 1
C
          CALL GRDCEL
C         ===========
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
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM), COEFBR ,
     &   GRAROX , GRAROY , GRAROZ ,
C        ------   ------   ------
     &   W1     , W2     , DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1, NCEL
            W3(IEL) = (GRAROX(IEL)*GX + GRAROY(IEL)*GY + GRAROZ(IEL)*GZ)
     &           /PROPCE(IEL,IPCROM)
            W3(IEL) = 1.D0/SQRT(MAX(EPZERO,W3(IEL)))
C
          ENDDO
C
C     On met le nombre de clippings a 0 (il le restera pour IDTVAR=0)
          NCLPTR = 0
C
        ENDIF
C
C
        IF (IDTVAR.EQ.1.OR.IDTVAR.EQ.2) THEN
C
          ICOU = 0
          IFOU = 0
C
C 4.1.1 LIMITATION PAR LE COURANT
C =============================
C
          IF ( COUMAX.GT.0.D0.AND.ICONV(IUIPH).GE.1 ) THEN
C
C     ICOU = 1 marque l'existence d'une limitation par le COURANT
            ICOU = 1
C
C ---> CONSTRUCTION DE U/DX           (COURANT         ) =W1
C
            IDIFF0 = 0
C
C     Matrice a priori non symetrique
            ISYM = 2
C
            CALL MATRDT
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONV(IUIPH)    , IDIFF0          , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COFBDT , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
            DO IEL = 1, NCEL
              ROM = PROPCE(IEL,IPCROM)
              W1    (IEL) = DAM(IEL)/(ROM*VOLUME(IEL))
            ENDDO
C
C ---> CALCUL DE W1     = PAS DE TEMPS VARIABLE VERIFIANT
C       LE NOMBRE DE COURANT         MAXIMUM PRESCRIT PAR L'UTILISATEUR
C
            DO IEL = 1, NCEL
              W1    (IEL) = COUMAX/MAX( W1    (IEL), EPZERO)
            ENDDO
C
C ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE
C
            IF (IDTVAR.EQ.1) THEN
              W1MIN = GRAND
              DO IEL = 1, NCEL
                W1MIN = MIN(W1MIN,W1(IEL))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMIN (W1MIN)
C               ===========
              ENDIF
              DO IEL = 1, NCEL
                W1(IEL) = W1MIN
              ENDDO
            ENDIF
C
          ENDIF
C
C 4.1.2 LIMITATION PAR LE FOURIER
C =============================
C
          IF ( FOUMAX.GT.0.D0.AND.IDIFF(IUIPH).GE.1 ) THEN
C
C     IFOU = 1 marque l'existence d'une limitation par le FOURIER
            IFOU = 1
C
            ICONV0 = 0
C                                   2
C ---> CONSTRUCTION DE      +2.NU/DX  (         FOURIER) =W2
C
C     Matrice a priori symetrique
            ISYM = 1
C
            CALL MATRDT
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONV0          , IDIFF(IUIPH)    , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COFBDT , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
            DO IEL = 1, NCEL
              ROM = PROPCE(IEL,IPCROM)
              W2    (IEL) = DAM(IEL)/(ROM*VOLUME(IEL))
            ENDDO
C
C ---> CALCUL DE W2     = PAS DE TEMPS VARIABLE VERIFIANT
C       LE NOMBRE DE         FOURIER MAXIMUM PRESCRIT PAR L'UTILISATEUR
C
            DO IEL = 1, NCEL
              W2    (IEL) = FOUMAX/MAX( W2    (IEL), EPZERO)
            ENDDO
C
C ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE
C
            IF (IDTVAR.EQ.1) THEN
              W2MIN = GRAND
              DO IEL = 1, NCEL
                W2MIN = MIN(W2MIN,W2(IEL))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMIN (W2MIN)
C               ===========
              ENDIF
              DO IEL = 1, NCEL
                W2(IEL) = W2MIN
              ENDDO
            ENDIF
C
          ENDIF
C
C 4.1.3 LIMITATION POUR L'ALGORITHME COMPRESSIBLE
C =============================================
C     Il est important de conserver WCF intact : on le reutilise
C     plus bas pour l'affichage
C
          ICOUCF = 0
          IF ( COUMAX.GT.0.D0.AND.IPPMOD(ICOMPF).GE.0 ) THEN
C
            ICOUCF = 1
C
C ---> CALCUL DE DAM     = PAS DE TEMPS VARIABLE VERIFIANT
C       LA CONTRAINTE CFL MAXIMUM PRESCRITE PAR L'UTILISATEUR
C
            DO IEL = 1, NCEL
              DAM(IEL) = COUMAX/MAX( WCF(IEL), EPZERO)
            ENDDO
C
C ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE
C
            IF (IDTVAR.EQ.1) THEN
              W3MIN = GRAND
              DO IEL = 1, NCEL
                W3MIN = MIN(W3MIN,DAM(IEL))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMIN (W3MIN)
C               ===========
              ENDIF
              DO IEL = 1, NCEL
                DAM(IEL) = W3MIN
              ENDDO
            ENDIF
C
          ENDIF
C
C 4.1.4 ON PREND LA PLUS CONTRAIGNANTE DES LIMITATIONS
C ==================================================
C    (le minimum des deux si elles existent et
C     celle qui existe s'il n'en existe qu'une)
C
          IF(ICOU.EQ.1.AND.IFOU.EQ.1) THEN
            DO IEL = 1, NCEL
              W1(IEL) = MIN(W1(IEL),W2(IEL))
            ENDDO
          ELSEIF(ICOU.EQ.0.AND.IFOU.EQ.1) THEN
            DO IEL = 1, NCEL
              W1(IEL) = W2(IEL)
            ENDDO
          ENDIF
C
C
C     En compressible, on prend obligatoirement
C     en compte la limitation associee à la masse volumique.
C
          IF(ICOUCF.EQ.1) THEN
            DO IEL = 1, NCEL
              W1(IEL) = MIN(W1(IEL),DAM(IEL))
            ENDDO
          ENDIF
C
C 4.1.5 ON CALCULE EFFECTIVEMENT LE PAS DE TEMPS
C ============================================
C
C --->  MONTEE           PROGRESSIVE DU PAS DE TEMPS
C              DESCENTE  IMMEDIATE   DU PAS DE TEMPS
C
          DO IEL = 1, NCEL
            IF( W1    (IEL).GE.DT(IEL) ) THEN
              UNPVDT = 1.D0+VARRDT
              DT(IEL) = MIN( UNPVDT*DT(IEL), W1    (IEL) )
            ELSE
              DT(IEL) =                      W1    (IEL)
            ENDIF
          ENDDO
C
C
C 4.1.6 ON LIMITE PAR LE PAS DE TEMPS "THERMIQUE" MAX
C =================================================
C     DTTMAX = W3 = 1/SQRT(MAX(0+,gradRO.g/RO)
C     on limite le pas de temps a DTTMAX
C
          IF (IPTLRO.EQ.1) THEN
C
C
C  On clippe le pas de temps a DTTMAX
C     (affiche dans ecrlis)
C
            NCLPTR = 0
C
            DO IEL = 1, NCEL
              IF ( DT(IEL).GT.W3(IEL) ) THEN
                NCLPTR = NCLPTR +1
                DT(IEL) = W3(IEL)
              ENDIF
            ENDDO
C
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (NCLPTR)
C             ===========
            ENDIF
C
C ---> PAS DE TEMPS UNIFORME : on reuniformise le pas de temps
C
            IF (IDTVAR.EQ.1) THEN
              W3MIN = GRAND
              DO IEL = 1, NCEL
                W3MIN = MIN(W3MIN,DT(IEL))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMIN (W3MIN)
C               ===========
              ENDIF
              DO IEL = 1, NCEL
                DT(IEL) = W3MIN
              ENDDO
            ENDIF
C
          ENDIF
C
C 4.1.7 ON CLIPPE LE PAS DE TEMPS PAR RAPPORT A DTMIN ET DTMAX
C ==========================================================
C
          ICFMIN = 0
          ICFMAX = 0
C
          DO IEL = 1, NCEL
C
            IF( DT(IEL).GT.DTMAX      ) THEN
              ICFMAX = ICFMAX +1
              DT(IEL) = DTMAX
            ENDIF
            IF( DT(IEL).LT.DTMIN      ) THEN
              ICFMIN = ICFMIN +1
              DT(IEL) = DTMIN
            ENDIF
C
          ENDDO
C
          IF (IRANGP.GE.0) THEN
            CALL PARCPT (ICFMIN)
C           ===========
            CALL PARCPT (ICFMAX)
C           ===========
          ENDIF
C
          ICLPMX(IPPDT) = ICFMAX
          ICLPMN(IPPDT) = ICFMIN
C
          IF( IWARNP.GE.2) THEN
            WRITE (NFECRA,1003) ICFMIN,DTMIN,ICFMAX,DTMAX
          ENDIF
C
        ENDIF
C
C     Rapport DT sur DTmax lie aux effets de densite
C       (affichage dans ecrlis)
        IF (IPTLRO.EQ.1) THEN
C
          DTSDTM = 0.D0
          DO IEL = 1, NCEL
            DTSDT0 = DT(IEL)/W3(IEL)
            IF ( DTSDT0 .GT. DTSDTM ) THEN
              DTSDTM = DTSDT0
              ICFMAX = IEL
            ENDIF
          ENDDO
          XYZMAX(1) = XYZCEN(1,ICFMAX)
          XYZMAX(2) = XYZCEN(2,ICFMAX)
          XYZMAX(3) = XYZCEN(3,ICFMAX)
C
          IF (IRANGP.GE.0) THEN
            NBRVAL = 3
            CALL PARMXL (NBRVAL, DTSDTM, XYZMAX)
C           ===========
          ENDIF
          RPDTRO(1) = DTSDTM
          RPDTRO(2) = XYZMAX(1)
          RPDTRO(3) = XYZMAX(2)
          RPDTRO(4) = XYZMAX(3)
C
        ENDIF
C
C=======================================================================
C 4.2  CALCUL DU NOMBRE DE COURANT POUR AFFICHAGE
C=======================================================================
C
        IF ( ICONV(IUIPH).GE.1.AND.
     &       (IWARNP.GE.2.OR.MODNTL.EQ.0) ) THEN
C
          IDIFF0 = 0
          CNOM   =' COURANT'
C
C     CONSTRUCTION DE U/DX           (COURANT         ) =W1
C
C MATRICE A PRIORI NON SYMETRIQUE
C
          ISYM = 2
C
          CALL MATRDT
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONV(IUIPH)    , IDIFF0          , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COFBDT , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1, NCEL
            ROM = PROPCE(IEL,IPCROM)
            W1    (IEL) = DAM(IEL)/(ROM*VOLUME(IEL))
          ENDDO
C
C     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM
C
          CFMAX = -GRAND
          CFMIN =  GRAND
          ICFMAX= 0
          ICFMIN= 0
C
          DO IEL = 1, NCEL
C
            COUFOU = W1(IEL)*DT(IEL)
            PROPCE(IEL,IPCCOU) = COUFOU
C
            IF( COUFOU.LE.CFMIN ) THEN
              CFMIN  = COUFOU
              ICFMIN = IEL
            ENDIF
C
            IF( COUFOU.GE.CFMAX ) THEN
              CFMAX  = COUFOU
              ICFMAX = IEL
            ENDIF
C
          ENDDO
C
          XYZMIN(1) = XYZCEN(1,ICFMIN)
          XYZMIN(2) = XYZCEN(2,ICFMIN)
          XYZMIN(3) = XYZCEN(3,ICFMIN)
          XYZMAX(1) = XYZCEN(1,ICFMAX)
          XYZMAX(2) = XYZCEN(2,ICFMAX)
          XYZMAX(3) = XYZCEN(3,ICFMAX)
C
          IF (IRANGP.GE.0) THEN
            NBRVAL = 3
            CALL PARMNL (NBRVAL, CFMIN, XYZMIN)
C           ===========
            CALL PARMXL (NBRVAL, CFMAX, XYZMAX)
C           ===========
          ENDIF
C
          IF(IWARNP.GE.2) THEN
            WRITE(NFECRA,1001) CNOM,CFMAX,XYZMAX(1),XYZMAX(2),XYZMAX(3)
            WRITE(NFECRA,1002) CNOM,CFMIN,XYZMIN(1),XYZMIN(2),XYZMIN(3)
          ENDIF
C
C       -> pour listing
          PTPLOC(1,1) = CFMIN
          PTPLOC(1,2) = XYZMIN(1)
          PTPLOC(1,3) = XYZMIN(2)
          PTPLOC(1,4) = XYZMIN(3)
          PTPLOC(2,1) = CFMAX
          PTPLOC(2,2) = XYZMAX(1)
          PTPLOC(2,3) = XYZMAX(2)
          PTPLOC(2,4) = XYZMAX(3)
C
        ENDIF
C
C=======================================================================
C 4.3  CALCUL DU NOMBRE DE FOURIER POUR AFFICHAGE
C=======================================================================
C
        IF ( IDIFF(IUIPH).GE.1.AND.
     &       (IWARNP.GE.2.OR.MODNTL.EQ.0) ) THEN
C
          ICONV0 = 0
          CNOM   =' FOURIER'
C                                   2
C     CONSTRUCTION DE      +2.NU/DX  (         FOURIER) =W1
C
C MATRICE A PRIORI SYMETRIQUE
C
          ISYM = 1
C
          CALL MATRDT
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONV0          , IDIFF(IUIPH)    , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COFBDT , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1, NCEL
            ROM = PROPCE(IEL,IPCROM)
            W1    (IEL) = DAM(IEL)/(ROM*VOLUME(IEL))
          ENDDO
C
C     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM
C
          CFMAX  = -GRAND
          CFMIN  =  GRAND
          ICFMAX = 0
          ICFMIN = 0
C
          DO IEL = 1, NCEL
C
            COUFOU = W1(IEL)*DT(IEL)
            PROPCE(IEL,IPCFOU) = COUFOU
C
            IF( COUFOU.LE.CFMIN ) THEN
              CFMIN  = COUFOU
              ICFMIN = IEL
            ENDIF
C
            IF( COUFOU.GE.CFMAX ) THEN
              CFMAX  = COUFOU
              ICFMAX = IEL
            ENDIF
C
          ENDDO
C
          XYZMIN(1) = XYZCEN(1,ICFMIN)
          XYZMIN(2) = XYZCEN(2,ICFMIN)
          XYZMIN(3) = XYZCEN(3,ICFMIN)
          XYZMAX(1) = XYZCEN(1,ICFMAX)
          XYZMAX(2) = XYZCEN(2,ICFMAX)
          XYZMAX(3) = XYZCEN(3,ICFMAX)
C
          IF (IRANGP.GE.0) THEN
            NBRVAL = 3
            CALL PARMNL (NBRVAL, CFMIN, XYZMIN)
C           ===========
            CALL PARMXL (NBRVAL, CFMAX, XYZMAX)
C           ===========
          ENDIF
C
          IF(IWARNP.GE.2) THEN
            WRITE(NFECRA,1001) CNOM,CFMAX,XYZMAX(1),XYZMAX(2),XYZMAX(3)
            WRITE(NFECRA,1002) CNOM,CFMIN,XYZMIN(1),XYZMIN(2),XYZMIN(3)
          ENDIF
C
C       -> pour listing
          PTPLOC(3,1) = CFMIN
          PTPLOC(3,2) = XYZMIN(1)
          PTPLOC(3,3) = XYZMIN(2)
          PTPLOC(3,4) = XYZMIN(3)
          PTPLOC(4,1) = CFMAX
          PTPLOC(4,2) = XYZMAX(1)
          PTPLOC(4,3) = XYZMAX(2)
          PTPLOC(4,4) = XYZMAX(3)
C
        ENDIF
C
C=======================================================================
C 4.4  CALCUL DU NOMBRE DE COURANT/FOURIER POUR AFFICHAGE
C=======================================================================
C
C     En incompressible uniquement (en compressible, on preferera
C       afficher la contrainte liee a la masse volumique)
C
        IF ( (IWARNP.GE.2.OR.MODNTL.EQ.0).AND.
     &      (IDIFF(IUIPH).GE.1.OR.ICONV(IUIPH).GE.1)
     &      .AND.(IPPMOD(ICOMPF).LT.0)               ) THEN
C
          CNOM   =' COU/FOU'
C                                   2
C     CONSTRUCTION DE U/DX +2.NU/DX  (COURANT +FOURIER) =W1
C
C MATRICE A PRIORI NON SYMETRIQUE
C
          ISYM = 1
          IF (ICONV(IUIPH).GT.0) ISYM = 2
C
          CALL MATRDT
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONV(IUIPH)    , IDIFF(IUIPH)    , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COFBDT , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   DAM    ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1, NCEL
            ROM = PROPCE(IEL,IPCROM)
            W1    (IEL) = DAM(IEL)/(ROM*VOLUME(IEL))
          ENDDO
C
C     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM
C
          CFMAX  = -GRAND
          CFMIN  =  GRAND
          ICFMAX = 0
          ICFMIN = 0
C
          DO IEL = 1, NCEL
C
            COUFOU = W1(IEL)*DT(IEL)
C
            IF( COUFOU.LE.CFMIN ) THEN
              CFMIN  = COUFOU
              ICFMIN = IEL
            ENDIF
C
            IF( COUFOU.GE.CFMAX ) THEN
              CFMAX  = COUFOU
              ICFMAX = IEL
            ENDIF
C
          ENDDO
C
          XYZMIN(1) = XYZCEN(1,ICFMIN)
          XYZMIN(2) = XYZCEN(2,ICFMIN)
          XYZMIN(3) = XYZCEN(3,ICFMIN)
          XYZMAX(1) = XYZCEN(1,ICFMAX)
          XYZMAX(2) = XYZCEN(2,ICFMAX)
          XYZMAX(3) = XYZCEN(3,ICFMAX)
C
          IF (IRANGP.GE.0) THEN
            NBRVAL = 3
            CALL PARMNL (NBRVAL, CFMIN, XYZMIN)
C           ===========
            CALL PARMXL (NBRVAL, CFMAX, XYZMAX)
C           ===========
          ENDIF
C
          IF(IWARNP.GE.2) THEN
            WRITE(NFECRA,1001) CNOM,CFMAX,XYZMAX(1),XYZMAX(2),XYZMAX(3)
            WRITE(NFECRA,1002) CNOM,CFMIN,XYZMIN(1),XYZMIN(2),XYZMIN(3)
          ENDIF
C
C       -> pour listing
          PTPLOC(5,1) = CFMIN
          PTPLOC(5,2) = XYZMIN(1)
          PTPLOC(5,3) = XYZMIN(2)
          PTPLOC(5,4) = XYZMIN(3)
          PTPLOC(6,1) = CFMAX
          PTPLOC(6,2) = XYZMAX(1)
          PTPLOC(6,3) = XYZMAX(2)
          PTPLOC(6,4) = XYZMAX(3)
C
        ENDIF
C
C=======================================================================
C 4.5  CALCUL DE LA CONTRAINTE CFL DE LA MASSE VOL. POUR AFFICHAGE
C=======================================================================
C
C En Compressible uniquement
C
        IF ( (IWARNP.GE.2.OR.MODNTL.EQ.0).AND.
     &       (IPPMOD(ICOMPF).GE.0)                        ) THEN
C
          CNOM   =' CFL/MAS'
C
C
C     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM
C
          CFMAX  = -GRAND
          CFMIN  =  GRAND
          ICFMAX = 0
          ICFMIN = 0
C
          DO IEL = 1, NCEL
C
            COUFOU = WCF(IEL)*DT(IEL)
C
            IF( COUFOU.LE.CFMIN ) THEN
              CFMIN  = COUFOU
              ICFMIN = IEL
            ENDIF
C
            IF( COUFOU.GE.CFMAX ) THEN
              CFMAX  = COUFOU
              ICFMAX = IEL
            ENDIF
C
          ENDDO
C
          XYZMIN(1) = XYZCEN(1,ICFMIN)
          XYZMIN(2) = XYZCEN(2,ICFMIN)
          XYZMIN(3) = XYZCEN(3,ICFMIN)
          XYZMAX(1) = XYZCEN(1,ICFMAX)
          XYZMAX(2) = XYZCEN(2,ICFMAX)
          XYZMAX(3) = XYZCEN(3,ICFMAX)
C
          IF (IRANGP.GE.0) THEN
            NBRVAL = 3
            CALL PARMNL (NBRVAL, CFMIN, XYZMIN)
C           ===========
            CALL PARMXL (NBRVAL, CFMAX, XYZMAX)
C           ===========
          ENDIF
C
          IF(IWARNP.GE.2) THEN
            WRITE(NFECRA,1001) CNOM,CFMAX,XYZMAX(1),XYZMAX(2),XYZMAX(3)
            WRITE(NFECRA,1002) CNOM,CFMIN,XYZMIN(1),XYZMIN(2),XYZMIN(3)
          ENDIF
C
C       -> pour listing
          PTPLOC(5,1) = CFMIN
          PTPLOC(5,2) = XYZMIN(1)
          PTPLOC(5,3) = XYZMIN(2)
          PTPLOC(5,4) = XYZMIN(3)
          PTPLOC(6,1) = CFMAX
          PTPLOC(6,2) = XYZMAX(1)
          PTPLOC(6,3) = XYZMAX(2)
          PTPLOC(6,4) = XYZMAX(3)
C
        ENDIF
C
C=======================================================================
C 5.   ALGORITHME STATIONNAIRE
C=======================================================================
      ELSE
C
        ISYM = 1
        IF (ICONV(IUIPH).GT.0) ISYM = 2
C
        CALL MATRDT
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONV(IUIPH)    , IDIFF(IUIPH)    , ISYM,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COEFB(1,IUIPH)  , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &                                                VISCF  , VISCB  ,
     &   DT     , RDEVEL , RTUSER , RA )
C
        DO IEL = 1, NCEL
          DT(IEL) = RELAXV(IUIPH)*PROPCE(IEL,IPCROM)
     &         *VOLUME(IEL)/MAX(DT(IEL),EPZERO)
        ENDDO

      ENDIF
C
C--------
C FORMATS
C--------
C
 1001 FORMAT ( /,A8,' MAX= ',E11.4,
     & ' EN ',E11.4,' ',E11.4,' ',E11.4)
 1002 FORMAT (   A8,' MIN= ',E11.4,
     & ' EN ',E11.4,' ',E11.4,' ',E11.4)
 1003 FORMAT ( /,'CLIPPINGS DE DT : ',
     &                             I10,' A ',E11.4,', ',I10,' A ',E11.4)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
