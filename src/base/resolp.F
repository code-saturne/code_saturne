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
                        SUBROUTINE RESOLP
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM , ISOSTD , IFACLG , IRESPR , IDTSCA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   FRCXT  , DFRCXT , TPUCOU , TRAV   ,
     &   VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   , ROVSDT , TSLAGR ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , FRCHY  , DFRCHY , COEFU  , TRAVA  ,
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
CARGU ! ISOSTD       ! TE !  -> ! INDICATEUR DE SORTIE STANDARD        !
CARGU !    (NFABOR+1)!    !     !  +NUMERO DE LA FACE DE REFERENCE     !
CARGU ! IFACLG(2,NFAC! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! IRESPR(NCELET! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! IDTSCA       ! E  !  -> ! INDICATEUR DE PAS DE TEMPS NON SCALAI!
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
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU !DFRCXT(NCELET,! TR !  -> ! VARIATION DE FORCE EXTERIEURE        !
CARGU !   3,NPHAS)   !    !     !  GENERANT LAPRESSION HYDROSTATIQUE   !
CARGU ! TPUCOU       ! TR !  -> ! COUPLAGE VITESSE PRESSION            !
CARGU ! (NCELEL,NDIM)!    !     !                                      !
CARGU ! TRAV(NCELET,3! TR !  -> ! SMB POUR NORMALISATION DE RESIDU     !
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
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !  (NCELET,*)  !    !     !   LAGRANGIEN                         !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! FRCHY(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !  NDIM  )     !    !     !                                      !
CARGU ! DFRCHY(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !  NDIM  )     !    !     !                                      !
CARGU ! COEFU(NFAB,3)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TRAVA        ! TR !  -> ! TABLEAU DE TRAVAIL POUR COUPLAGE     !
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
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
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
      INTEGER          ISOSTD(NFABOR+1,NPHAS)
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
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS), DFRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION TPUCOU(NCELET,NDIM), TRAV(NCELET,3)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION VISCFI(NFAC), VISCBI(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET)
      DOUBLE PRECISION SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION FRCHY(NCELET,NDIM), DFRCHY(NCELET,NDIM)
      DOUBLE PRECISION COEFU(NFABOR,3), TRAVA(NCELET,NDIM,NPHAS)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ICCOCG, INC   , INIT  , ISYM  , IPOL  , ISQRT
      INTEGER          IEL   , IFAC  , IFAC0 , IEL0
      INTEGER          IRESLP, NSWRP , NSWMPR
      INTEGER          ISWEEP, NITERF, ICYCLE, NGR   , IFINIA , IFINRA
      INTEGER          IFLMB0, IFCSOR
      INTEGER          NSWRGP, IMLIGP, IWARNP, IWARMG
      INTEGER          IPRIPH, IUIPH , IVIPH , IWIPH ,ICLIPF
      INTEGER                  ICLIPR, ICLIUP, ICLIVP, ICLIWP
      INTEGER          IPCROM, IPCROA, IPBROM, IFLMAS, IFLMAB
      INTEGER          IPP
      INTEGER          IISMPH
      INTEGER          IDIFFP, ICONVP, NDIRCP
      INTEGER          NITMAP, IMGRP , NCYMAP, NITMGP
      INTEGER          IINVPE, IMASPE, INDHYD
      INTEGER          IDIMTE, ITENSO, IESDEP
      INTEGER          IDTSCA
      DOUBLE PRECISION RESIDU, PHYDR0
      DOUBLE PRECISION ARDTSR, ARSR  , ARAKPH, UNSARA, THETAP
      DOUBLE PRECISION DTSROM, UNSVOM, ROMRO0, RO0IPH
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, EPSILP
      DOUBLE PRECISION DROM  , DRONM1
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C --- Memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Impressions
      IPP    = IPPRTP(IPR(IPHAS))
C
C --- Variables
      IPRIPH = IPR(IPHAS)
      IUIPH  = IU (IPHAS)
      IVIPH  = IV (IPHAS)
      IWIPH  = IW (IPHAS)
C
C --- Conditions aux limites
C     (ICLRTP(IPRIPH,ICOEFF) pointe vers ICLRTP(IPRIPH,ICOEF) si IPHYDR=0)
      ICLIPR = ICLRTP(IPRIPH,ICOEF)
      ICLIPF = ICLRTP(IPRIPH,ICOEFF)
      ICLIUP = ICLRTP(IUIPH ,ICOEF)
      ICLIVP = ICLRTP(IVIPH ,ICOEF)
      ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
      IISMPH = IISYMP+NFABOR*(IPHAS-1)
C
C --- Grandeurs physiques
      IPCROM = IPPROC(IROM  (IPHAS ))
      IF(ICALHY.EQ.1) THEN
        IPCROA = IPPROC(IROMA(IPHAS))
      ELSE
        IPCROA = 0
      ENDIF
      IPBROM = IPPROB(IROM  (IPHAS ))
      IFLMAS = IPPROF(IFLUMA(IPRIPH))
      IFLMAB = IPPROB(IFLUMA(IPRIPH))
C
      RO0IPH = RO0(IPHAS)
C
C --- Options de resolution
      ISYM  = 1
      IF( ICONV (IPRIPH).GT.0 ) THEN
        ISYM  = 2
      ENDIF
C
      IF (IRESOL(IPRIPH).EQ.-1) THEN
        IRESLP = 0
        IPOL   = 0
        IF( ICONV(IPRIPH).GT.0 ) THEN
          IRESLP = 1
          IPOL   = 0
        ENDIF
      ELSE
        IRESLP = MOD(IRESOL(IPRIPH),1000)
        IPOL   = (IRESOL(IPRIPH)-IRESLP)/1000
      ENDIF
C
      ARAKPH = ARAK(IPHAS)
C
      ISQRT = 1
C
C=======================================================================
C 2.  RESIDU DE NORMALISATION
C=======================================================================
C
      IF(IRNPNW.NE.1) THEN
C
        IF (IPHYDR.EQ.1) THEN
          DO IEL = 1, NCEL
            UNSVOM = -1.D0/VOLUME(IEL)
            TRAV(IEL,1) = TRAV(IEL,1)*UNSVOM
     &           + FRCXT(IEL,1,IPHAS)
     &           + DFRCXT(IEL,1,IPHAS)
            TRAV(IEL,2) = TRAV(IEL,2)*UNSVOM
     &           + FRCXT(IEL,2,IPHAS)
     &           + DFRCXT(IEL,2,IPHAS)
            TRAV(IEL,3) = TRAV(IEL,3)*UNSVOM
     &           + FRCXT(IEL,3,IPHAS)
     &           + DFRCXT(IEL,3,IPHAS)
          ENDDO
        ELSE
          IF(ISNO2T(IPHAS).GT.0) THEN
            DO IEL = 1, NCEL
              UNSVOM = -1.D0/VOLUME(IEL)
              ROMRO0 = PROPCE(IEL,IPCROM)-RO0IPH
              TRAV(IEL,1) = (TRAV(IEL,1)+TRAVA(IEL,1,IPHAS))*UNSVOM
     &             + ROMRO0*GX
              TRAV(IEL,2) = (TRAV(IEL,2)+TRAVA(IEL,2,IPHAS))*UNSVOM
     &             + ROMRO0*GY
              TRAV(IEL,3) = (TRAV(IEL,3)+TRAVA(IEL,3,IPHAS))*UNSVOM
     &             + ROMRO0*GZ
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              UNSVOM = -1.D0/VOLUME(IEL)
              ROMRO0 = PROPCE(IEL,IPCROM)-RO0IPH
              TRAV(IEL,1) = TRAV(IEL,1)*UNSVOM
     &             + ROMRO0*GX
              TRAV(IEL,2) = TRAV(IEL,2)*UNSVOM
     &             + ROMRO0*GY
              TRAV(IEL,3) = TRAV(IEL,3)*UNSVOM
     &             + ROMRO0*GZ
            ENDDO
          ENDIF
        ENDIF
        DO IEL = 1, NCEL
          DTSROM = DT(IEL)/PROPCE(IEL,IPCROM)
          TRAV(IEL,1) = RTP(IEL,IUIPH) +DTSROM*TRAV(IEL,1)
          TRAV(IEL,2) = RTP(IEL,IVIPH) +DTSROM*TRAV(IEL,2)
          TRAV(IEL,3) = RTP(IEL,IWIPH) +DTSROM*TRAV(IEL,3)
        ENDDO
C
C ---> TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (TRAV(1,1))
C          ===========
          CALL PARCOM (TRAV(1,2))
C          ===========
          CALL PARCOM (TRAV(1,3))
C          ===========
        ENDIF
C
C ON IMPOSE LA PERIODICITE SUR TRAV
C
        IF(IPERIO.EQ.1) THEN
C
          IDIMTE = 1
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    TRAV(1,1) , TRAV(1,1) , TRAV(1,1) ,
     &    TRAV(1,2) , TRAV(1,2) , TRAV(1,2) ,
     &    TRAV(1,3) , TRAV(1,3) , TRAV(1,3) )
C
        ENDIF
C
C ON NE RECONSTRUIT PAS POUR GAGNER DU TEMPS
C   EPSRGR N'EST DONC PAS UTILISE
C
        INIT   = 1
        INC    = 1
        ICCOCG = 1
        IFLMB0 = 1
        IF (IALE.EQ.1) IFLMB0 = 0
        NSWRP  = 1
        IMLIGP = IMLIGR(IUIPH )
        IWARNP = IWARNI(IPRIPH)
        EPSRGP = EPSRGR(IUIPH )
        CLIMGP = CLIMGR(IUIPH )
        EXTRAP = EXTRAG(IUIPH )
C
        IMASPE = 1
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
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
        INIT = 1
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &       IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W1)
C
        IF (NCESMP.GT.0) THEN
          DO IEL = 1, NCESMP
            W1(ICETSM(IEL)) = W1(ICETSM(IEL))
     &       -VOLUME(ICETSM(IEL))*SMACEL(IEL,IPRIPH)/PROPCE(IEL,IPCROM)
          ENDDO
        ENDIF
C
C ---> LAGRANGIEN : COUPLAGE RETOUR
C
        IF (IILAGR.EQ.2 .AND. LTSMAS.EQ.1 .AND. IPHAS.EQ.ILPHAS) THEN
C
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL) -TSLAGR(IEL,ITSMAS) / PROPCE(IEL,IPCROM)
          ENDDO
C
        ENDIF
C
        CALL PRODSC(NCELET,NCEL,ISQRT,W1,W1,RNORMP(IPHAS))
C
        IF(IWARNI(IPRIPH).GE.2) THEN
          CHAINE = NOMVAR(IPP)
          WRITE(NFECRA,1300)CHAINE(1:8) ,RNORMP(IPHAS)
        ENDIF
        DERVAR(IPP) = RNORMP(IPHAS)
        NBIVAR(IPP) = 0
C
      ELSE
C
        IF(IWARNI(IPRIPH).GE.2) THEN
          CHAINE = NOMVAR(IPP)
          WRITE(NFECRA,1300)CHAINE(1:8) ,RNORMP(IPHAS)
        ENDIF
        DERVAR(IPP) = RNORMP(IPHAS)
        NBIVAR(IPP) = 0
C
      ENDIF
C
C=======================================================================
C 3.  CALCUL DE L'INCREMENT DE PRESSION HYDROSTATIQUE (SI NECESSAIRE)
C=======================================================================
C
      IF (IPHYDR.EQ.1) THEN
C     L'INCREMENT EST STOCKE PROVISOIREMENT DANS RTP(.,IPRIPH)
C     on resout une equation de Poisson avec des conditions de
C     flux nul partout
C     Ce n'est utile que si on a des faces de sortie
        IFCSOR = ISOSTD(NFABOR+1,IPHAS)
        IF (IRANGP.GE.0) THEN
          CALL PARCMX (IFCSOR)
        ENDIF
C
        IF (IFCSOR.LE.0) THEN
          INDHYD = 0
        ELSE
          DO IFAC=1,NFABOR
            COEFA(IFAC,ICLIPF) = 0.D0
            COEFB(IFAC,ICLIPF) = 1.D0
          ENDDO
C
          IF (ICALHY.EQ.1) THEN
C
C
C     Il serait necessaire de communiquer pour periodicite et parallelisme
C      avec PARCOM et PERCOM sur le vecteur
C      DFRCHY(IEL,1) DFRCHY(IEL,2) DFRCHY(IEL,3)
C     On peut economiser la communication tant que DFRCHY ne depend que de
C      RHO et RHO n-1 qui ont ete communiques auparavant.
C     Exceptionnellement, on fait donc le calcul sur NCELET.
            DO IEL = 1, NCELET
              DRONM1 = (PROPCE(IEL,IPCROA)-RO0IPH)
              DROM   = (PROPCE(IEL,IPCROM)-RO0IPH)
              FRCHY(IEL,1)  = DRONM1*GX
              FRCHY(IEL,2)  = DRONM1*GY
              FRCHY(IEL,3)  = DRONM1*GZ
              DFRCHY(IEL,1) = DROM  *GX - FRCHY(IEL,1)
              DFRCHY(IEL,2) = DROM  *GY - FRCHY(IEL,2)
              DFRCHY(IEL,3) = DROM  *GZ - FRCHY(IEL,3)
            ENDDO
C
            CALL CALHYD
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR , INDHYD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCHY (1,1) , FRCHY (1,2) , FRCHY (1,3) ,
     &   DFRCHY(1,1) , DFRCHY(1,2) , DFRCHY(1,3) ,
     &   RTP(1,IPRIPH)   , PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   COEFA(1,ICLIPF) , COEFB(1,ICLIPF) ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , ROVSDT ,
     &   RDEVEL , RTUSER , RA     )
          ELSE
            INDHYD = 0
          ENDIF
C
        ENDIF
      ENDIF
C
C
C=======================================================================
C 4.  PREPARATION DE LA MATRICE DU SYSTEME A RESOUDRE
C=======================================================================
C
C ---> TERME INSTATIONNAIRE
C
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
C
C ---> "VITESSE" DE DIFFUSION FACETTE
C
      IF( IDIFF(IPRIPH).GE. 1 ) THEN
        IF (IDTSCA.EQ.0) THEN
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
     &   DT     ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
        ELSE
          CALL VISORT
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   TPUCOU(1,1) , TPUCOU(1,2) , TPUCOU(1,3) ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
      ELSE
        DO IFAC = 1, NFAC
          VISCF(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = 0.D0
        ENDDO
      ENDIF
C
      ICONVP = ICONV (IPRIPH)
      IDIFFP = IDIFF (IPRIPH)
      NDIRCP = NDIRCL(IPRIPH)
C
      THETAP = 1.D0
      CALL MATRIX
C     ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   ICONVP , IDIFFP , NDIRCP ,
     &   ISYM   , NFECRA ,
     &   THETAP ,
     &   IFACEL , IFABOR ,
     &   COEFB(1,ICLIPR) , ROVSDT ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB), VISCF  , VISCB  ,
     &   DAM    , XAM    )
C
C=======================================================================
C 5.  INITIALISATION DU FLUX DE MASSE
C=======================================================================
C
C --- Flux de masse predit et premiere composante Rhie et Chow
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
     &   RTPA(1,IPRIPH)  , COEFA(1,ICLIPR) , COEFB(1,ICLIPR)  ,
     &   TRAV(1,1) , TRAV(1,2) , TRAV(1,3) ,
C        ---------   ---------   ---------
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C
      IF (IPHYDR.EQ.1) THEN
        DO IEL = 1, NCEL
          TRAV(IEL,1) = TRAV(IEL,1) - FRCXT(IEL,1,IPHAS)
          TRAV(IEL,2) = TRAV(IEL,2) - FRCXT(IEL,2,IPHAS)
          TRAV(IEL,3) = TRAV(IEL,3) - FRCXT(IEL,3,IPHAS)
        ENDDO
      ENDIF
C
      IF (IDTSCA.EQ.0) THEN
        DO IEL = 1, NCEL
          ARDTSR  = ARAKPH*(DT(IEL)/PROPCE(IEL,IPCROM))
          TRAV(IEL,1) = RTP(IEL,IUIPH) + ARDTSR*TRAV(IEL,1)
          TRAV(IEL,2) = RTP(IEL,IVIPH) + ARDTSR*TRAV(IEL,2)
          TRAV(IEL,3) = RTP(IEL,IWIPH) + ARDTSR*TRAV(IEL,3)
        ENDDO
      ELSE
        DO IEL=1,NCEL
          ARSR  = ARAKPH/PROPCE(IEL,IPCROM)
          TRAV(IEL,1) = RTP(IEL,IUIPH)+ARSR*TPUCOU(IEL,1)*TRAV(IEL,1)
          TRAV(IEL,2) = RTP(IEL,IVIPH)+ARSR*TPUCOU(IEL,2)*TRAV(IEL,2)
          TRAV(IEL,3) = RTP(IEL,IWIPH)+ARSR*TPUCOU(IEL,3)*TRAV(IEL,3)
        ENDDO
      ENDIF
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (TRAV(1,1))
C       ===========
        CALL PARCOM (TRAV(1,2))
C       ===========
        CALL PARCOM (TRAV(1,3))
C       ===========
      ENDIF
C
C ON IMPOSE LA PERIODICITE SUR TRAV
C
      IF(IPERIO.EQ.1) THEN
C
        IDIMTE = 1
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    TRAV(1,1) , TRAV(1,1) , TRAV(1,1) ,
     &    TRAV(1,2) , TRAV(1,2) , TRAV(1,2) ,
     &    TRAV(1,3) , TRAV(1,3) , TRAV(1,3) )
C
      ENDIF
C
      INIT   = 1
      INC    = 1
      ICCOCG = 1
      IFLMB0 = 1
      IF (IALE.EQ.1) IFLMB0 = 0
      NSWRGP = NSWRGR(IUIPH )
      IMLIGP = IMLIGR(IUIPH )
      IWARNP = IWARNI(IPRIPH)
      EPSRGP = EPSRGR(IUIPH )
      CLIMGP = CLIMGR(IUIPH )
      EXTRAP = EXTRAG(IUIPH )
C
      IMASPE = 1
C
      CALL INIMAS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
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
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
C
C --- Projection aux faces des forces exterieures
C
      IF (IPHYDR.EQ.1) THEN
        INIT   = 0
        INC    = 0
        ICCOCG = 1
        NSWRGP = NSWRGR(IPRIPH)
        IMLIGP = IMLIGR(IPRIPH)
        IWARNP = IWARNI(IPRIPH)
        EPSRGP = EPSRGR(IPRIPH)
        CLIMGP = CLIMGR(IPRIPH)
C
        IF (IDTSCA.EQ.0) THEN
          CALL PROJTS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   COEFB(1,ICLIPR) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   VISCF  , VISCB  ,
     &   DT     , DT     , DT     ,
     &   RDEVEL , RTUSER , RA     )
        ELSE
          CALL PROJTS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   COEFB(1,ICLIPR) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1)     , TPUCOU(1,2)     , TPUCOU(1,3)     ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
      ENDIF
C
      INIT   = 0
      INC    = 1
      ICCOCG = 1
C
      IF(ARAKPH.GT.0.D0) THEN
C
C --- Prise en compte de Arak : la viscosite face est multipliee
C       Le pas de temps aussi. On retablit plus bas.
        DO IFAC = 1, NFAC
          VISCF(IFAC) = ARAKPH*VISCF(IFAC)
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = ARAKPH*VISCB(IFAC)
        ENDDO

        IF (IDTSCA.EQ.0) THEN
          DO IEL = 1, NCEL
            DT(IEL) = ARAKPH*DT(IEL)
          ENDDO
C
          NSWRGP = NSWRGR(IPRIPH )
          IMLIGP = IMLIGR(IPRIPH )
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH )
          CLIMGP = CLIMGR(IPRIPH )
          EXTRAP = EXTRAG(IPRIPH )
          CALL ITRMAS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCXT(1,1,IPHAS), FRCXT(1,2,IPHAS), FRCXT(1,3,IPHAS),
     &   RTPA(1,IPRIPH)  , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   DT     , DT     , DT     ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Projection du terme source pour oter la partie hydrostat de la pression
          IF (IPHYDR.EQ.1) THEN
            INIT   = 0
            INC    = 0
            ICCOCG = 1
            NSWRGP = NSWRGR(IPRIPH)
            IMLIGP = IMLIGR(IPRIPH)
            IWARNP = IWARNI(IPRIPH)
            EPSRGP = EPSRGR(IPRIPH)
            CLIMGP = CLIMGR(IPRIPH)
C     on passe avec un pseudo coefB=1, pour avoir 0 aux faces de bord
            DO IFAC = 1,NFABOR
              COEFB(IFAC,ICLIPF) = 1.D0
            ENDDO
C
            CALL PROJTS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCXT(1,1,IPHAS), FRCXT(1,2,IPHAS), FRCXT(1,3,IPHAS),
     &   COEFB(1,ICLIPF) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   VISCF  , VISCB  ,
     &   DT     , DT     , DT     ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C --- Correction du pas de temps
          UNSARA = 1.D0/ARAKPH
          DO IEL = 1, NCEL
            DT(IEL) = DT(IEL)*UNSARA
          ENDDO
C
        ELSE
C
          DO IEL = 1, NCEL
            TPUCOU(IEL,1) = ARAKPH*TPUCOU(IEL,1)
            TPUCOU(IEL,2) = ARAKPH*TPUCOU(IEL,2)
            TPUCOU(IEL,3) = ARAKPH*TPUCOU(IEL,3)
          ENDDO
C
          NSWRGP = NSWRGR(IPRIPH )
          IMLIGP = IMLIGR(IPRIPH )
          IWARNP = IWARNI(IPRIPH )
          EPSRGP = EPSRGR(IPRIPH )
          CLIMGP = CLIMGR(IPRIPH )
          EXTRAP = EXTRAG(IPRIPH )
          CALL ITRMAS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCXT(1,1,IPHAS), FRCXT(1,2,IPHAS), FRCXT(1,3,IPHAS),
     &   RTPA(1,IPRIPH)  , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1)     , TPUCOU(1,2)     , TPUCOU(1,3)     ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Projection du terme source pour oter la partie hydrostat de la pression
          IF (IPHYDR.EQ.1) THEN
            INIT   = 0
            INC    = 0
            ICCOCG = 1
            NSWRGP = NSWRGR(IPRIPH)
            IMLIGP = IMLIGR(IPRIPH)
            IWARNP = IWARNI(IPRIPH)
            EPSRGP = EPSRGR(IPRIPH)
            CLIMGP = CLIMGR(IPRIPH)
C     on passe avec un pseudo coefB=1, pour avoir 0 aux faces de bord
            DO IFAC = 1,NFABOR
              COEFB(IFAC,ICLIPF) = 1.D0
            ENDDO
C
            CALL PROJTS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCXT(1,1,IPHAS), FRCXT(1,2,IPHAS), FRCXT(1,3,IPHAS),
     &   COEFB(1,ICLIPF) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1)     , TPUCOU(1,2)     , TPUCOU(1,3)     ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
C --- Correction du pas de temps
          UNSARA = 1.D0/ARAKPH
          DO IEL = 1, NCEL
            TPUCOU(IEL,1) = UNSARA*TPUCOU(IEL,1)
            TPUCOU(IEL,2) = UNSARA*TPUCOU(IEL,2)
            TPUCOU(IEL,3) = UNSARA*TPUCOU(IEL,3)
          ENDDO
C
        ENDIF

C --- Correction de la viscosite aux faces
        DO IFAC = 1, NFAC
          VISCF(IFAC) = VISCF(IFAC)*UNSARA
        ENDDO
        DO IFAC = 1, NFABOR
          VISCB(IFAC) = VISCB(IFAC)*UNSARA
        ENDDO
C
      ENDIF
C
C     Calcul des CL pour l'increment de pression
C     On commence par affecter les CL classiques
C     (COEFA=0 et COEFB=COEFB(P), puis on change
C     les CL de sortie en mettant COEFA a l'increment
C     de pression hydrostatique, decale pour valoir 0
C     sur la face de reference
      IF (IPHYDR.EQ.1) THEN
        DO IFAC=1,NFABOR
          COEFA(IFAC,ICLIPF) = 0.D0
        ENDDO
        IF (INDHYD.EQ.1) THEN
          IFAC0 = ISOSTD(NFABOR+1,IPHAS)
          IF (IFAC0.LE.0) THEN
            PHYDR0 = 0.D0
          ELSE
            IEL0 = IFABOR(IFAC0)
            PHYDR0 = RTP(IEL0,IPRIPH)
     &           +(CDGFBO(1,IFAC0)-XYZCEN(1,IEL0))*DFRCXT(IEL0,1,IPHAS)
     &           +(CDGFBO(2,IFAC0)-XYZCEN(2,IEL0))*DFRCXT(IEL0,2,IPHAS)
     &           +(CDGFBO(3,IFAC0)-XYZCEN(3,IEL0))*DFRCXT(IEL0,3,IPHAS)
          ENDIF
          IF (IRANGP.GE.0) THEN
            CALL PARSOM (PHYDR0)
          ENDIF
          DO IFAC=1,NFABOR
            IF (ISOSTD(IFAC,IPHAS).EQ.1) THEN
              IEL=IFABOR(IFAC)
              COEFA(IFAC,ICLIPF) = RTP(IEL,IPRIPH)
     &             +(CDGFBO(1,IFAC)-XYZCEN(1,IEL))*DFRCXT(IEL,1,IPHAS)
     &             +(CDGFBO(2,IFAC)-XYZCEN(2,IEL))*DFRCXT(IEL,2,IPHAS)
     &             +(CDGFBO(3,IFAC)-XYZCEN(3,IEL))*DFRCXT(IEL,3,IPHAS)
     &             - PHYDR0
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C
C=======================================================================
C 6.  PREPARATION DU MULTIGRILLE ALGEBRIQUE
C=======================================================================
C
      IF (IMGR(IPRIPH).GT.0) THEN
C
C --- Initialisations
C
        IFINIA = IDEBIA
        IFINRA = IDEBRA
        NGR = 0
C
C --- Creation de la hierarchie de maillages
C
        IWARMG = IWARNI(IPRIPH)
        CALL CLMLGA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IWARMG , NFECRA , IFINIA , IFINRA , NGR    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   RDEVEL , RTUSER , RA     )
C
        IDEBIA = IFINIA
        IDEBRA = IFINRA
C
      ENDIF
C
C=======================================================================
C 7.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
C=======================================================================
C
C --- Nombre de sweeps
      NSWMPR = NSWRSM(IPRIPH)
C
C --- Mise a zero des variables
C       RTP(.,IPR) sera l'increment de pression cumule
C       DRTP       sera l'increment d'increment a chaque sweep
C       W7         sera la divergence du flux de masse predit
      DO IEL = 1,NCEL
        RTP(IEL,IPRIPH) = 0.D0
        DRTP(IEL) = 0.D0
        SMBR(IEL) = 0.D0
      ENDDO
C
C --- Divergence initiale
      INIT = 1
      CALL DIVMAS
     &  (NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W7)
C
C --- Termes sources de masse
      IF (NCESMP.GT.0) THEN
        DO IEL = 1, NCESMP
          W7(ICETSM(IEL)) = W7(ICETSM(IEL))
     &           -VOLUME(ICETSM(IEL))*SMACEL(IEL,IPRIPH)
        ENDDO
      ENDIF
C
C ---> Terme sources de masse
      IF (IILAGR.EQ.2 .AND. LTSMAS.EQ.1 .AND. IPHAS.EQ.ILPHAS) THEN
        DO IEL = 1, NCEL
          W7(IEL) = W7(IEL) -TSLAGR(IEL,ITSMAS)
        ENDDO
      ENDIF
C
C --- Boucle de reconstruction : debut
      DO 100 ISWEEP = 1, NSWMPR
C
C --- Mise a jour du second membre
C     (signe "-" a cause de celui qui est implicitement dans la matrice)
        DO IEL = 1, NCEL
          SMBR(IEL) = - W7(IEL) - SMBR(IEL)
        ENDDO
C
C --- Test de convergence du calcul
C
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBR,SMBR,RESIDU)
        IF (IWARNI(IPRIPH).GE.2) THEN
           CHAINE = NOMVAR(IPP)
           WRITE(NFECRA,1400)CHAINE(1:8),ISWEEP,RESIDU
        ENDIF
        IF (ISWEEP.EQ.1) RNSMBR(IPP) = RESIDU
C
C  Test a modifier eventuellement
C (il faut qu'il soit plus strict que celui de gradco)
        IF( RESIDU .LE. 10.D0*EPSILO(IPRIPH)*RNORMP(IPHAS) ) THEN
C --- Si convergence, calcul de l'indicateur
C                     mise a jour du flux de masse et sortie
C
C
C --- Calcul d'indicateur, avec prise en compte
C       du volume (norme L2) ou non
C
          IF(IESCAL(IESDER,IPHAS).GT.0) THEN
            IESDEP = IPPROC(IESTIM(IESDER,IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IESDEP) = ABS(SMBR(IEL))/VOLUME(IEL)
            ENDDO
            IF(IESCAL(IESDER,IPHAS).EQ.2) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IESDEP) =
     &            PROPCE(IEL,IESDEP)*SQRT(VOLUME(IEL))
              ENDDO
            ENDIF
          ENDIF
C
C
C
          ICCOCG = 1
          INIT = 0
          INC  = 0
          IF (IPHYDR.EQ.1) INC = 1
C --- en cas de prise en compte de Phydro, on met INC=1 pour prendre en
C     compte les CL de COEFA(.,ICLIPF)
          NSWRGP = NSWRGR(IPRIPH)
          IMLIGP = IMLIGR(IPRIPH)
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH)
          CLIMGP = CLIMGR(IPRIPH)
          EXTRAP = EXTRAG(IPRIPH)
          IF (IDTSCA.EQ.0) THEN
            CALL ITRMAS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   RTP(1,IPRIPH)   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   DT     , DT     , DT     ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSE
C
            CALL ITRMAS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   RTP(1,IPRIPH)   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1) , TPUCOU(1,2) , TPUCOU(1,3) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
          GOTO 101
C
        ENDIF
C
C --- Resolution implicite sur l'increment d'increment DRTP
        DO IEL = 1, NCEL
          DRTP(IEL) = 0.D0
        ENDDO
C
        CHAINE = NOMVAR(IPP)
        NITMAP = NITMAX(IPRIPH)
        IMGRP  = IMGR  (IPRIPH)
        NCYMAP = NCYMAX(IPRIPH)
        NITMGP = NITMGF(IPRIPH)
        IWARNP = IWARNI(IPRIPH)
        EPSILP = EPSILO(IPRIPH)
C
C ---> TRAITEMENT PERIODICITE
C     (La pression est un scalaire,
C                 pas de pb pour la rotation: IINVPE=1)
       IINVPE = 1
C
        CALL INVERS
C       ===========
     & ( CHAINE(1:8)     , IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IPOL   , IRESLP , NITMAP , IMGRP  , NGR    ,
     &   NCYMAP , NITMGP ,
     &   IWARNP , NFECRA , NITERF , ICYCLE , IINVPE ,
     &   EPSILP , RNORMP(IPHAS)   , RESIDU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , SMBR   , DRTP   ,
     &   DAG    , XAG    , W1     , W2     ,
     &   W3     , W4     , W5     , W6     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
        NBIVAR(IPP) = NBIVAR(IPP) + NITERF
        IF(ABS(RNORMP(IPHAS)).GT.0.D0) THEN
          RESVAR(IPP) = RESIDU/RNORMP(IPHAS)
        ELSE
          RESVAR(IPP) = 0.D0
        ENDIF
C
        IF( ISWEEP.EQ.NSWMPR ) THEN
C
C --- Si dernier sweep :
C       Calcul d'estimateur
C       Incrementation du flux de masse
C         avec reconstruction a partir de (dP)^(NSWMPR-1)
C       Puis on rajoute la correction en (d(dP))^(NSWMPR)
C         sans reconstruction pour assurer la divergence nulle
C
C
C --- Calcul d'indicateur, avec prise en compte
C       du volume (norme L2) ou non
C
          IF(IESCAL(IESDER,IPHAS).GT.0) THEN
            IESDEP = IPPROC(IESTIM(IESDER,IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IESDEP) = ABS(SMBR(IEL))/VOLUME(IEL)
            ENDDO
            IF(IESCAL(IESDER,IPHAS).EQ.2) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IESDEP) =
     &            PROPCE(IEL,IESDEP)*SQRT(VOLUME(IEL))
              ENDDO
            ENDIF
          ENDIF
C
C --- Incrementation du flux de masse et correction
C
          ICCOCG = 1
          INIT = 0
          INC  = 0
          IF (IPHYDR.EQ.1) INC = 1
          NSWRGP = NSWRGR(IPRIPH)
          IMLIGP = IMLIGR(IPRIPH)
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH)
          CLIMGP = CLIMGR(IPRIPH)
          EXTRAP = EXTRAG(IPRIPH)
C
          IF (IDTSCA.EQ.0) THEN
            CALL ITRMAS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   RTP(1,IPRIPH)   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   DT          , DT          , DT          ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     pour le dernier increment, on ne reconstruit pas, on n'appelle donc
C     pas GRDCEL. La valeur des DFRCXT (qui doit normalement etre nul)
C     est donc sans importance
          ICCOCG = 0
          NSWRP = 0
          INC = 0
C
            CALL ITRMAS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRP  , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   DRTP            , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   DT          , DT          , DT          ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSE
C
            CALL ITRMAS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   RTP(1,IPRIPH)   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1) , TPUCOU(1,2) , TPUCOU(1,3) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     pour le dernier increment, on ne reconstruit pas, on n'appelle donc
C     pas GRDCEL. La valeur des DFRCXT (qui doit normalement etre nul)
C     est donc sans importance
          ICCOCG = 0
          NSWRP = 0
          INC = 0
C
            CALL ITRMAS
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRP  , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   DRTP            , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1) , TPUCOU(1,2) , TPUCOU(1,3) ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
C     Mise a jour de l'increment de pression
          DO IEL = 1, NCEL
            RTP(IEL,IPRIPH) = RTP(IEL,IPRIPH) + DRTP(IEL)
          ENDDO
C
        ELSE
C
C --- Si ce n'est pas le dernier sweep
C       Mise a jour de l'increment de pression et calcul direct de la
C       partie en gradient d'increment de pression du second membre
C       (avec reconstruction)
C
          IF (IDTVAR.GE.0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IPRIPH) = RTP(IEL,IPRIPH)
     &                        + RELAXV(IPRIPH)*DRTP(IEL)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              RTP(IEL,IPRIPH) = RTP(IEL,IPRIPH)
     &                        + DRTP(IEL)
            ENDDO
          ENDIF
C
          ICCOCG = 1
          INIT = 1
          INC  = 0
          IF (IPHYDR.EQ.1) INC = 1
          NSWRGP = NSWRGR(IPRIPH)
          IMLIGP = IMLIGR(IPRIPH)
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH)
          CLIMGP = CLIMGR(IPRIPH)
          EXTRAP = EXTRAG(IPRIPH)
C
          IF (IDTSCA.EQ.0) THEN
C
            CALL ITRGRP
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   RTP(1,IPRIPH)   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   DT          , DT          , DT          ,
     &   SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSE
C
            CALL ITRGRP
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   RTP(1,IPRIPH)   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR) ,
     &   VISCF  , VISCB  ,
     &   TPUCOU(1,1) , TPUCOU(1,2) , TPUCOU(1,3) ,
     &   SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
        ENDIF
C
 100  CONTINUE
C --- Boucle de reconstruction : fin
C
      IF(IWARNI(IPRIPH).GE.2) THEN
         CHAINE = NOMVAR(IPP)
         WRITE( NFECRA,1600)CHAINE(1:8),NSWMPR
      ENDIF
C
 101  CONTINUE
C
C     REACTUALISATION DE LA PRESSION
C
      IF (IDTVAR.LT.0) THEN
        DO IEL = 1, NCEL
          RTP(IEL,IPRIPH) = RTPA(IEL,IPRIPH)
     &         + RELAXV(IPRIPH)*RTP(IEL,IPRIPH)
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          RTP(IEL,IPRIPH) = RTPA(IEL,IPRIPH) + RTP(IEL,IPRIPH)
        ENDDO
      ENDIF
C
C--------
C FORMATS
C--------
C
 1300 FORMAT(1X,A8,' : RESIDU DE NORMALISATION =', E14.6)
 1400 FORMAT(1X,A8,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)
 1600 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' ETAPE DE PRESSION                   ',/,
     &'@    *********                                               ',/,
     &'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,
     &'@                                                            '  )
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
