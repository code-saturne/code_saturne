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
                        SUBROUTINE CFQDMV
C                       *****************
C     --------------------------------------------------------------
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
     &   FLUMAS , FLUMAB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL , FRCXT  , DFRCXT ,
     &   TPUCOU , TRAV   , VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
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
CARGU !DFRCXT(NCELET,! TR !  -> ! VARIATION DE FORCE EXTERIEURE        !
CARGU !   3,NPHAS)   !    !     !  GENERANT LAPRESSION HYDROSTATIQUE   !
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
CARGU ! COEFU(NFAB,3)! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
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
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS), DFRCXT(NCELET,3,NPHAS)
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
      DOUBLE PRECISION COEFU(NFABOR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL   , IELPDC, IFAC  , IVAR  , ISOU  , III
      INTEGER          ICCOCG, INC   , INIT  , IPHYDP, II
      INTEGER          IRESLP, NSWRGP, IMLIGP, IWARNP, IPP
      INTEGER          IPRIPH, IKIPH , IUIPH , IVIPH , IWIPH
      INTEGER          ICLIK , ICLVAR, ICLVAF, ICLIPR
      INTEGER          IPCROM, IPCVIS, IPCVST
      INTEGER          ICONVP, IDIFFP, NDIRCP, NITMAP, NSWRSP
      INTEGER          IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IDIMTE, ITENSO
      INTEGER          IDIAEX, ITERNS
      INTEGER          IIFRU
      DOUBLE PRECISION RNORM , VITNOR
      DOUBLE PRECISION ROMVOM, RTPROM
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, BLENCP, EPSILP
      DOUBLE PRECISION VIT1  , VIT2  , VIT3  , THETAP, PFAC, PFAC1
      DOUBLE PRECISION CPDC11, CPDC22, CPDC33, CPDC12, CPDC13, CPDC23
      DOUBLE PRECISION D2S3  , PBORD , DIIPBX, DIIPBY, DIIPBZ, PIP, XKB
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
      IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &     .OR. ITURB(IPHAS).EQ.60) THEN
        ICLIK  = ICLRTP(IKIPH ,ICOEF)
      ENDIF
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C
C     Indicateur flux de bord Rusanov
      IF(IIFBRU.GT.0) THEN
        IIFRU = IIFBRU+(IPHAS-1)*NFABOR
      ELSE
        IIFRU = 1
      ENDIF
C
C=======================================================================
C 2.  GRADIENT DE PRESSION ET GRAVITE
C=======================================================================
C
C ---> PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE :
C
      IF (IPHYDR.EQ.1) THEN
C
C     on doit pouvoir adapter l'option iphydr au compressible,
C       mais noter plusieurs points
C       - on dispose de la masse volumique au pas de temps courant et au
C         pas de temps précédent (au pas de temps précédent dans propce
C         en particulier)
C       - la correction de pression est ici généree par la resolution de
C         l'energie (la pression change sans que rho ne change)
C       - si l'objectif se limite a adapter le calcul de grad p pour
C         qu'il compense le rho0 g, noter quand meme que l'on ne resout
C         pas en rho-rho0 et que P est tjrs cohérent avec rho (par la
C         thermo)
C
        DO IEL = 1, NCEL
C
C variation de force (utilise dans resolp)
C
          IF(IGRDPP(IPHAS).GT.0) THEN
            RTPROM = RTP(IEL,ISCA(IRHO(IPHAS)))
          ELSE
            RTPROM = RTPA(IEL,ISCA(IRHO(IPHAS)))
          ENDIF
C
          DFRCXT(IEL,1,IPHAS) = RTPROM*GX - FRCXT(IEL,1,IPHAS)
          DFRCXT(IEL,2,IPHAS) = RTPROM*GY - FRCXT(IEL,2,IPHAS)
          DFRCXT(IEL,3,IPHAS) = RTPROM*GZ - FRCXT(IEL,3,IPHAS)
        ENDDO
C     Ajout eventuel des pertes de charges
        IF (NCEPDP.GT.0) THEN
          DO IELPDC = 1, NCEPDP
            IEL=ICEPDC(IELPDC)
            VIT1   = RTP(IEL,IUIPH)
            VIT2   = RTP(IEL,IVIPH)
            VIT3   = RTP(IEL,IWIPH)
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
     & -RTP(IEL,ISCA(IRHO(IPHAS)))*(CPDC11*VIT1+CPDC12*VIT2+CPDC13*VIT3)
            DFRCXT(IEL,2,IPHAS) = DFRCXT(IEL,2,IPHAS)
     & -RTP(IEL,ISCA(IRHO(IPHAS)))*(CPDC12*VIT1+CPDC22*VIT2+CPDC23*VIT3)
            DFRCXT(IEL,3,IPHAS) = DFRCXT(IEL,3,IPHAS)
     & -RTP(IEL,ISCA(IRHO(IPHAS)))*(CPDC13*VIT1+CPDC23*VIT2+CPDC33*VIT3)
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
C       Fin du test sur IPHYDR
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
     &   RTP(1,IPRIPH)   , COEFA(1,ICLRTP(IPRIPH,ICOEF))  ,
     &                     COEFB(1,ICLRTP(IPRIPH,ICOEF))  ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
      IF (IPHYDR.EQ.1) THEN
        DO IEL = 1, NCEL
          TRAV(IEL,1) = ( FRCXT(IEL,1,IPHAS) - W1(IEL) )*VOLUME(IEL)
          TRAV(IEL,2) = ( FRCXT(IEL,2,IPHAS) - W2(IEL) )*VOLUME(IEL)
          TRAV(IEL,3) = ( FRCXT(IEL,3,IPHAS) - W3(IEL) )*VOLUME(IEL)
        ENDDO
      ELSE
        DO IEL = 1, NCEL
C
          IF(IGRDPP(IPHAS).GT.0) THEN
            RTPROM = RTP(IEL,ISCA(IRHO(IPHAS)))
          ELSE
            RTPROM = RTPA(IEL,ISCA(IRHO(IPHAS)))
          ENDIF
C
          TRAV(IEL,1) = ( RTPROM*GX - W1(IEL) )*VOLUME(IEL)
          TRAV(IEL,2) = ( RTPROM*GY - W2(IEL) )*VOLUME(IEL)
          TRAV(IEL,3) = ( RTPROM*GZ - W3(IEL) )*VOLUME(IEL)
        ENDDO
      ENDIF
C
C    Calcul des efforts aux parois (partie 2/5), si demande
C    La pression a la face est calculee comme dans gradrc/gradmc
      IF (INEEDF.EQ.1) THEN
        ICLIPR = ICLRTP(IPRIPH,ICOEF)
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
          DO ISOU = 1, 3
            RA(IFORBR+(IFAC-1)*NDIM + ISOU-1) =
     &           RA(IFORBR+(IFAC-1)*NDIM + ISOU-1)
     &           + PFAC*SURFBO(ISOU,IFAC)
          ENDDO
        ENDDO
      ENDIF
C
C
C     Elimination du flux au bord associé au gradient de pression :
C       il est pris en compte par les conditions aux limites dans
C       le flux de Rusanov
C
      IF(IIFBRU.GT.0) THEN
C
        DO IFAC = 1, NFABOR
C
          IF(IA(IIFRU+IFAC-1).EQ.1) THEN
C
            IEL = IFABOR(IFAC)
C
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
C
            PIP = RTP(IEL,IPRIPH)
     &           +(W1(IEL)*DIIPBX+W2(IEL)*DIIPBY+W3(IEL)*DIIPBZ)
C
            PBORD = COEFA(IFAC,ICLRTP(IPRIPH,ICOEF))
     &           + COEFB(IFAC,ICLRTP(IPRIPH,ICOEF))*PIP
C
            TRAV(IEL,1) = TRAV(IEL,1) + PBORD*SURFBO(1,IFAC)
            TRAV(IEL,2) = TRAV(IEL,2) + PBORD*SURFBO(2,IFAC)
            TRAV(IEL,3) = TRAV(IEL,3) + PBORD*SURFBO(3,IFAC)
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
C     Flux de C.L. associé à Rusanov (PROPFB contient la contribution
C       de - div(rho u u) - grad P si on est passé dans cfrusb
C       ou 0 sinon).
C     Pour ne pas ajouter le flux div(rho u u) deux fois, on a remplace
C       codits et bilsc2 par cfcdts et cfbsc2 qui ne different des
C       precedents que par les indicateurs qui permettent de
C       ne pas prendre en compte le flux convectif aux faces de bord
C       pour lesquelles on est passe dans cfrusb
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
        TRAV(IEL,1) =  TRAV(IEL,1) - PROPFB(IFAC,IPPROB(IFBRHU(IPHAS)))
        TRAV(IEL,2) =  TRAV(IEL,2) - PROPFB(IFAC,IPPROB(IFBRHV(IPHAS)))
        TRAV(IEL,3) =  TRAV(IEL,3) - PROPFB(IFAC,IPPROB(IFBRHW(IPHAS)))
      ENDDO
C
C
C ---> 2/3 RHO * GRADIENT DE K SI k epsilon
C      NB : ON NE PREND PAS LE GRADIENT DE (RHO K), MAIS
C           CA COMPLIQUERAIT LA GESTION DES CL ...
C
      IF( (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &     .OR. ITURB(IPHAS).EQ.60) .AND.IGRHOK(IPHAS).EQ.1) THEN
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
     &   RTP(1,IKIPH)    , COEFA(1,ICLIK)  , COEFB(1,ICLIK)  ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        D2S3 = 2.D0/3.D0
        DO IEL = 1, NCEL
          ROMVOM = -RTP(IEL,ISCA(IRHO(IPHAS)))*VOLUME(IEL)*D2S3
          TRAV(IEL,1) = TRAV(IEL,1) + W1(IEL) * ROMVOM
          TRAV(IEL,2) = TRAV(IEL,2) + W2(IEL) * ROMVOM
          TRAV(IEL,3) = TRAV(IEL,3) + W3(IEL) * ROMVOM
        ENDDO
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
C
      ENDIF
C
C
C
C ---> TERMES DE GRADIENT TRANSPOSE
C
      IF (IVISSE(IPHAS).EQ.1) THEN
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
      ENDIF
C
C
C ---> TERMES DE PERTES DE CHARGE
C     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT
C
      IF((NCEPDP.GT.0).AND.(IPHYDR.EQ.0)) THEN
C
        IDIAEX = 1
        CALL TSEPDC
C       ===========
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
     &   RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , TRAV   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C ---> - DIVERGENCE DE RIJ
C
      IF(ITYTUR(IPHAS).EQ.3 ) THEN
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
     &   RTP    , PROPCE , PROPFA , PROPFB ,
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
          DO IEL = 1, NCEL
             TRAV(IEL,ISOU) = TRAV(IEL,ISOU) -  W1(IEL)
          ENDDO
C
        ENDDO
C
      ENDIF
C
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
            W1(IEL) = PROPCE(IEL,IPCVIS) + PROPCE(IEL,IPCVST)
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
        IF(ITYTUR(IPHAS).EQ.3 .AND. IRIJNU(IPHAS).EQ.1) THEN
          DO IEL = 1, NCEL
            W1(IEL) = PROPCE(IEL,IPCVIS) + PROPCE(IEL,IPCVST)
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
        ENDIF
        IF(ISOU.EQ.2) THEN
          IVAR = IVIPH
        ENDIF
        IF(ISOU.EQ.3) THEN
          IVAR = IWIPH
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
          SMBR  (IEL) = 0.D0
          DRTP  (IEL) = 0.D0
        ENDDO
C
        CALL USTSNS
C       ===========
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
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBR   , DRTP   ,
C        ------   ------
     &   DAM    , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
        DO IEL = 1, NCEL
          ROVSDT(IEL) = MAX(-DRTP(IEL),ZERO)
          SMBR  (IEL) = SMBR(IEL) + DRTP(IEL) * RTP(IEL,IVAR)
        ENDDO
C
C
C ---> TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume
C
        INIT = 1
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                           IFACEL,IFABOR,FLUMAS,FLUMAB,W1)

C
C
C ---> AJOUT DANS LE TERME SOURCE ET DANS LE TERME INSTATIONNAIRE
C
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR  (IEL) +
     &         TRAV(IEL,ISOU)+ICONV(IVAR)*W1(IEL)*RTPA(IEL,IVAR)
        ENDDO
C
        DO IEL = 1, NCEL
          ROVSDT(IEL) = ROVSDT(IEL) +
     &         ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &         -ICONV(IVAR)*W1(IEL)
        ENDDO
C
C
C ---> PERTES DE CHARGE
C
        IF (NCEPDP.GT.0) THEN
          DO IELPDC = 1, NCEPDP
            IEL = ICEPDC(IELPDC)
            ROVSDT(IEL) = ROVSDT(IEL) +
     &           PROPCE(IEL,IPCROM)*VOLUME(IEL)*CKUPDC(IELPDC,ISOU)
          ENDDO
        ENDIF
C
C
C --->  TERMES DE SOURCE DE MASSE
C
        IF (NCESMP.GT.0) THEN
          ITERNS = 1
          CALL CATSMA ( NCELET , NCEL , NCESMP , ITERNS ,
     &                  ISNO2T(IPHAS), THETAV(IVAR),
     &                  ICETSM , ITYPSM(1,IVAR)  ,
     &                  VOLUME , RTP(1,IVAR) , SMACEL(1,IVAR) ,
     &                  SMACEL(1,IPR(IPHAS)) , SMBR , ROVSDT , W1)
        ENDIF
C
C
C
C ---> PARAMETRES POUR LA RESOLUTION DU SYSTEME
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
        THETAP = THETAV(IVAR)
        IESCAP = 0
C
C
C ---> FIN DE LA CONSTRUCTION ET DE LA RESOLUTION DU SYSTEME
C
        CALL CFCDTS
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP ,  NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP , IIFBRU ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , IA(IIFRU) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA(1,IVAR)    , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     FLUMAS , FLUMAB ,
     &   VISCFI , VISCBI , VISCF  , VISCB  ,
     &   ROVSDT , SMBR   , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     PAS DE COUPLAGE INSTATIONNAIRE EN COMPRESSIBLE
C
      ENDDO
C
C --->  FIN DE LA BOUCLE SUR U, V, W,
C
C
C ---> IMPRESSION DE NORME
C
      IF (IWARNI(IUIPH).GE.2) THEN
        RNORM = -1.D0
        DO IEL = 1, NCEL
          VITNOR =
     &     SQRT(RTP(IEL,IUIPH)**2+RTP(IEL,IVIPH)**2+RTP(IEL,IWIPH)**2)
          RNORM = MAX(RNORM,VITNOR)
        ENDDO
        IF (IRANGP.GE.0) CALL PARMAX (RNORM)
C                        ===========
        WRITE(NFECRA,1100) IPHAS, RNORM
      ENDIF
C
C
C--------
C FORMATS
C--------
C
 1100 FORMAT(/,
     & 1X,'Phase ',I4,' : Vitesse maximale apres qdm ',E12.4)
C----
C FIN
C----
C
      RETURN
      END
c@z
