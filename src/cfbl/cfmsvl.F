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
                        SUBROUTINE CFMSVL
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
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    , W11    , W12    ,
     &   WFLMAS , WFLMAB ,
     &   COEFU  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
CFONC    POUR LA MASSE VOLUMIQUE SUR UN PAS DE TEMPS
CFONC    (ALGORITHME COMPRESSIBLE EN RHO, U, E)
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
CARGU ! W1..12(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! WFLMAS(NFAC) ! TR !  -  ! TABLEAU DE W FLUX DE MASSE AUX FACES !
CARGU ! WFLMAB(NFABOR! TR !  -  ! TABLEAU DE W FLUX DE MASSE AUX BORDS !
CARGU ! COEFU(NFABO,3! TR !  -  ! TABLEAU DE TRAVAIL CL DE LA QDM      !
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
      DOUBLE PRECISION W1(NCELET) , W2(NCELET) , W3(NCELET)
      DOUBLE PRECISION W4(NCELET) , W5(NCELET) , W6(NCELET)
      DOUBLE PRECISION W7(NCELET) , W8(NCELET) , W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET), W12(NCELET)
      DOUBLE PRECISION WFLMAS(NFAC), WFLMAB(NFABOR)
      DOUBLE PRECISION COEFU(NFABOR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR  , IPHAS
      INTEGER          IFAC  , IEL
      INTEGER          INIT  , INC   , ICCOCG, ISQRT , II, JJ, III
      INTEGER          ICLVAR, ICLVAF
      INTEGER          IFLMAS, IFLMAB
      INTEGER          IPPVAR, IPP   , IPHYDP
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          ISTATP, ICONVP, IDIFFP, IRESLP, NDIRCP, NITMAP
      INTEGER          NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, BLENCP, EPSILP
      DOUBLE PRECISION EPSRSP
      DOUBLE PRECISION SCLNOR
C
      INTEGER          ICCFTH, IMODIF
      INTEGER          IDIMTE, ITENSO
      INTEGER          IIJ
      INTEGER          IWFABG, IWFBBG
      DOUBLE PRECISION DIJPFX, DIJPFY, DIJPFZ, POND
      DOUBLE PRECISION DIIPFX, DIIPFY, DIIPFZ, DJJPFX, DJJPFY, DJJPFZ
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
      DOUBLE PRECISION PIP   , PJP   , THETV, RELAXP
C
C***********************************************************************
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IWFABG = IDEBRA
      IWFBBG = IWFABG+NFAC
      IDEBRA = IWFBBG+NFABOR
      CALL RASIZE('CFMSVL',IDEBRA)
C
C --- Numero de phase associee au scalaire traite
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
C --- Flux de masse associe a l'energie
      IFLMAS = IPPROF(IFLUMA(ISCA(IENERG(IPHAS))))
      IFLMAB = IPPROB(IFLUMA(ISCA(IENERG(IPHAS))))
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
      DO IEL = 1, NCEL
        SMBRS(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
C
C
C     TERME SOURCE DE MASSE
C     =====================
C
      IF (NCESMP.GT.0) THEN
        DO II = 1, NCESMP
          IEL = ICETSM(II)
          SMBRS(IEL) = SMBRS(IEL) + SMACEL(IEL,IPR(IPHAS))*VOLUME(IEL)
        ENDDO
      ENDIF
C
C
C     TERME INSTATIONNAIRE
C     ====================
C
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL) + ISTAT(IVAR)*(VOLUME(IEL)/DT(IEL))
      ENDDO
C
C=======================================================================
C 3. CALCUL DU "FLUX DE MASSE" ET DE LA "VISCOSITE" AUX FACES
C=======================================================================
C
C     Ici VISCF et VISCB sont deux tableaux de travail.
C     On calcule WFLMAS et WFLMAB, WFABGS , WFBBGS
C
      CALL CFMSGS
C     ===========
     & ( IDEBIA , IDEBRA ,
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
     &   WFLMAS , WFLMAB , RA(IWFABG) , RA(IWFBBG) ,
C        ------   ------   ------   ------
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    , W11    , W12    ,
     &   VISCF  , VISCB  , COEFU  , XAM    ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C
C     Calcul du gradient de rho pour la reconstruction de rho
C       (contribution du terme convectif)
C
      IRCFLP = IRCFLU(IVAR)
C
      IF(IRCFLP.GT.0) THEN
C
        INC    = 1
        ICCOCG = 1
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IPHYDP = 0
        IWARNP = IWARNI(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IVAR)    ,
     &   COEFA(1,ICLRTP(IVAR,ICOEF)) , COEFB(1,ICLRTP(IVAR,ICOEF)) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSE
          DO II = 1, NCELET
            W1(II) = 0.D0
            W2(II) = 0.D0
            W3(II) = 0.D0
          ENDDO
        ENDIF
C
C
C     Au bord, on écrase WFLMAB pour imposer le débit souhaité.
C     Si on explicite le terme de convection (choix retenu par défaut,
C       seul testé), pas de problème.
C     Si on implicite le terme de convection, on n'impose pas
C       nécessairement le bon débit (cela dépend du traitement de
C       la convection au bord dans bilsc2 et de la condition à la
C       limite sur rho : ainsi, si l'on se place sur une sortie pour
C       laquelle la condition n'est pas valeur bord = valeur interne,
C       la convection étant traitée en upwind, c'est la valeur interne
C       de rho qui intervient dans le flux de bord et non pas
C       la valeur de bord comme supposé ci-dessous.
C
      IF(ICONV(IVAR).LE.0) THEN
C
        IF(IRCFLP.GT.0) THEN
C
          DO IFAC = 1, NFAC
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
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            PIP = RTPA(II,IVAR)
     &           + IRCFLP*(W1(II)*DIIPFX+W2(II)*DIIPFY+W3(II)*DIIPFZ)
            PJP = RTPA(JJ,IVAR)
     &           + IRCFLP*(W1(JJ)*DJJPFX+W2(JJ)*DJJPFY+W3(JJ)*DJJPFZ)
C
            WFLMAS(IFAC) = -0.5D0*
     &           ( PIP          *(WFLMAS(IFAC)+ABS(WFLMAS(IFAC)))
     &           + PJP          *(WFLMAS(IFAC)-ABS(WFLMAS(IFAC))))
          ENDDO
C
        ELSE
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            WFLMAS(IFAC) = -0.5D0*
     &           ( RTPA(II,IVAR)*(WFLMAS(IFAC)+ABS(WFLMAS(IFAC)))
     &           + RTPA(JJ,IVAR)*(WFLMAS(IFAC)-ABS(WFLMAS(IFAC))))
          ENDDO
        ENDIF
C
        DO IFAC = 1, NFABOR
          WFLMAB(IFAC) = -PROPFB(IFAC,IFLMAB)
        ENDDO
C
        INIT = 0
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &              IFACEL,IFABOR,WFLMAS,WFLMAB,SMBRS)
C
        DO IFAC = 1, NFAC
          RA(IWFABG+IFAC-1) = - RA(IWFABG+IFAC-1)
        ENDDO
        DO IFAC = 1, NFABOR
          RA(IWFBBG+IFAC-1) = - RA(IWFBBG+IFAC-1)
        ENDDO
        INIT = 0
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &              IFACEL,IFABOR,RA(IWFABG),RA(IWFBBG),SMBRS)
C
      ELSE
C
        IF(IRCFLP.GT.0) THEN
C
          DO IFAC = 1, NFABOR
            II = IFABOR(IFAC)
C
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
C
            PIP = RTPA(II,IVAR)
     &           +IRCFLP*(W1(II)*DIIPBX+W2(II)*DIIPBY+W3(II)*DIIPBZ)
C
            WFLMAB(IFAC) = -PROPFB(IFAC,IFLMAB)/
     &           ( COEFA(IFAC,ICLRTP(IVAR,ICOEF))
     &           + COEFB(IFAC,ICLRTP(IVAR,ICOEF))*PIP           )
          ENDDO
C
        ELSE
          DO IFAC = 1, NFABOR
            II = IFABOR(IFAC)
            WFLMAB(IFAC) = -PROPFB(IFAC,IFLMAB)/
     &           ( COEFA(IFAC,ICLRTP(IVAR,ICOEF))
     &           + COEFB(IFAC,ICLRTP(IVAR,ICOEF))*RTPA(II,IVAR) )
          ENDDO
        ENDIF
C
      ENDIF
C
C
C     On calcule VISCF et VISCB
C
      CALL CFMSVS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   VISCF  , VISCB  ,
C        ------   ------
     &   W1     , W9     , W10    ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C
C     On annule la viscosité au bord afin que la contribution
C       au flux de bord soit exactement WFLMAB (dans lequel on a mis
C       le flux de masse souhaité). Si on n'annule pas, on risque
C       d'obtenir des contributions non nulles de la partie diffusive,
C       sauf si on a imposé une condition de Neumann homogene sur rho.
C       Pour le moment, on prefere prendre des precautions (la
C       modification de VISCB se traduit simplement par la modification
C       de la matrice portant sur les incréments, ou encore par une
C       une condition de Neumann homogène sur les incréments).
C
      DO IFAC = 1, NFABOR
        VISCB (IFAC) = 0.D0
      ENDDO
C
C=======================================================================
C 4. RESOLUTION
C=======================================================================
C
      ISTATP = ISTAT (IVAR)
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
      THETV  = THETAV(IVAR)
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
     &   COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     WFLMAS          , WFLMAB          ,
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
C --- Clipping aux bornes définies par l'utilisateur ou par defaut
C       (par défaut, pas de borne contraignate)
C
C     Valeur bidon
      III = 1
C
      CALL CLPSCA
C     ===========
     & ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &   PROPCE , RTP(1,III)      , RTP    )
C
C --- Traitement utilisateur pour gestion plus fine des bornes
C       et actions correctives éventuelles.
        ICCFTH = -2
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
     &   W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER , RA     )
C
C
C --- Bilan explicite (voir codits : on enleve l'increment)
C
      IF (IWARNI(IVAR).GE.2) THEN
        DO IEL = 1, NCEL
          SMBRS(IEL) = SMBRS(IEL)
     &            - ISTAT(IVAR)*(VOLUME(IEL)/DT(IEL))
     &                *(RTP(IEL,IVAR)-RTPA(IEL,IVAR))
     &                * MAX(0,MIN(NSWRSM(IVAR)-2,1))
        ENDDO
        ISQRT = 1
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBRS,SMBRS,SCLNOR)
        WRITE(NFECRA,1200)CHAINE(1:8) ,SCLNOR
      ENDIF
C
C=======================================================================
C 6. COMMUNICATION DE RHO
C=======================================================================
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (RTP(1,IVAR))
C       ===========
      ENDIF
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     & ( IDIMTE , ITENSO ,
     &   RTP(1,IVAR)     , RTP(1,IVAR)     , RTP(1,IVAR),
     &   RTP(1,IVAR)     , RTP(1,IVAR)     , RTP(1,IVAR),
     &   RTP(1,IVAR)     , RTP(1,IVAR)     , RTP(1,IVAR) )
      ENDIF
C
C     On ne remplit pas PROPCE et PROPFB ici, car on veut disposer de
C       rho à l'instant précédent pour résoudre en qdm et energie
C     On modifiera PROPCE et PROPFB après resolution de l'energie.
C
C
C=======================================================================
C 7. CALCUL DU FLUX DE MASSE ACOUSTIQUE AUX FACES
C=======================================================================
C
C Ce flux est stocke en tant que flux de masse associe a l'energie
      INC    = 1
      ICCOCG = 1
C
      CALL CFBSC3
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP(1,IVAR)     , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &                     COEFA(1,ICLVAF) , COEFB(1,ICLVAF) ,
     &                     WFLMAS          , WFLMAB          ,
     &   VISCF  , VISCB  ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
C        ----------------  ----------------
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Si ICONV = 0, le terme convectif n'est pas traite par CFBSC3
C       il faut donc le rajouter a la main
C     Noter egalement que si ICONV = 0, PROPFB contient zero au sortir
C       de cfbsc3 et qu'on lui ajoute donc le flux de masse WFLMAB
C       (cohérent avec le flux imposé au bord et avec un signe négatif,
C        car WFLMAB etait, ci-dessus, utilise au second membre)
      IF(ICONV(IVAR).LE.0) THEN
        DO IFAC = 1, NFAC
          PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMAS) - WFLMAS(IFAC)
     &                                              - RA(IWFABG+IFAC-1)
        ENDDO
        DO IFAC = 1, NFABOR
          PROPFB(IFAC,IFLMAB) = PROPFB(IFAC,IFLMAB) - WFLMAB(IFAC)
     &                                              - RA(IWFBBG+IFAC-1)
        ENDDO
      ENDIF
C
C
C=======================================================================
C 8. ACTUALISATION DE LA PRESSION
C=======================================================================
C                               Pred      n+1  n
C On utilise l'equation d'etat P    =P(rho   ,e )
C
C --- Calcul de P au centre des cellules et actualisation de RTP
      IF(IGRDPP(IPHAS).GT.0) THEN
C
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
     &   RTP(1,IPR(IPHAS))        , W8     , W9     , W10    ,
C        -----------------
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 9. COMMUNICATION DE LA PRESSION
C=======================================================================
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTP(1,IPR(IPHAS)))
C         ===========
        ENDIF
C
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     & ( IDIMTE , ITENSO ,
     &   RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)),
     &   RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)),
     &   RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)), RTP(1,IPR(IPHAS)))
        ENDIF
C
      ENDIF
C
C     Pas de vérification ni de clipping de la pression, car on
C       considere que la masse volumique et l'énergie ont été vérifiees,
C       sont correctes et donc que la pression l'est aussi (elle
C       vient d'être calculée avec un sous-pgm utilisateur).
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
      END
c@z
