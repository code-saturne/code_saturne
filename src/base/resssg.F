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
                        SUBROUTINE RESSSG
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITPSMP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , GRDVIT , GRAROX , GRAROY , GRAROZ ,
     &   CKUPDC , SMCELP , GAMMA  ,
     &   VISCF  , VISCB  , COEFAX ,
     &   TSLAGE , TSLAGI ,
     &   DAM    , XAM    , DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
CFONC    POUR Rij modele SSG
CFONC VAR  = R11 R22 R33 R12 R13 R23
CFONC ISOU =  1   2   3   4   5   6
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
CARGU ! ITPSMP       ! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LA      !
CARGU ! (NCESMP)     !    !     !  VARIABLES (cf. USTSMA)              !
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
CARGU ! GRDVIT       ! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME GRAD   !
CARGU !  (NCELET,3,3)!    !     !                         DE VITESSE   !
CARGU ! GRAROX,Y,Z   ! TR !  -> ! TABLEAU DE TRAVAIL POUR GRAD ROM     !
CARGU !  (NCELET)    !    !     !                                      !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMCELP(NCESMP! TR !  -> ! VALEUR DE LA VARIABLE ASSOCIEE A LA  !
CARGU !              !    !     !  SOURCE DE MASSE                     !
CARGU ! GAMMA(NCESMP)! TR !  -> ! VALEUR DU FLUX DE MASSE              !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! COEFAX(NFABOR! TR !  -  ! TAB DE TRAV POUR COND.LIM. PAROI     !
CARGU !              ! TR !  -  !   ATTENTION : UNIQUEMENT AVEC ECHO   !
CARGU !              ! TR !  -  !   DE PAROI ET ABS(ICDPAR) = 1        !
CARGU ! TSLAGE(NCELET! TR !  -> ! TS EXPLICITE COUPLAGE RETOUR LAGR.   !
CARGU ! TSLAGI(NCELET! TR !  -> ! TS IMPLICITE COUPLAGE RETOUR LAGR.   !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
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
      INCLUDE "pointe.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IPHAS  , IVAR   , ISOU   , IPP
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITPSMP(NCESMP)
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
      DOUBLE PRECISION GRDVIT(NCELET,3,3)
      DOUBLE PRECISION GRAROX(NCELET), GRAROY(NCELET), GRAROZ(NCELET)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP)
      DOUBLE PRECISION SMCELP(NCESMP), GAMMA(NCESMP)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR), COEFAX(NFABOR)
      DOUBLE PRECISION TSLAGE(NCELET),TSLAGI(NCELET)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          INIT  , IFAC  , IEL
      INTEGER          II    , JJ    , KK    , IIUN  , III   , JJJ
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEIPH , IUIPH
      INTEGER          IPCROM, IPCVIS, IFLMAS, IFLMAB, IPCROO
      INTEGER          ICLVAR, ICLVAF
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP
      INTEGER          NITMAP, NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IPTSTA
      INTEGER          MAXELT, ILS
      DOUBLE PRECISION BLENCP, EPSILP, EPSRGP, CLIMGP, EXTRAP, RELAXP
      DOUBLE PRECISION EPSRSP
      DOUBLE PRECISION TRPROD, TRRIJ , RCTSE , DELTIJ
      DOUBLE PRECISION TUEXPR, THETS , THETV , THETP1
      DOUBLE PRECISION AIKSJK, AIKRJK, AII ,AKLSKL, AIKAKJ
      DOUBLE PRECISION XANISO(3,3), XSTRAI(3,3), XROTAC(3,3), XPROD(3,3)
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
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1000) NOMVAR(IPP)
      ENDIF
C
      IUIPH  = IU  (IPHAS)
      IR11IP = IR11(IPHAS)
      IR22IP = IR22(IPHAS)
      IR33IP = IR33(IPHAS)
      IR12IP = IR12(IPHAS)
      IR13IP = IR13(IPHAS)
      IR23IP = IR23(IPHAS)
      IEIPH  = IEP (IPHAS)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IUIPH))
      IFLMAB = IPPROB(IFLUMA(IUIPH))
C
      ICLVAR = ICLRTP(IVAR,ICOEF)
      ICLVAF = ICLRTP(IVAR,ICOEFF)
C
      DELTIJ = 1.0D0
      IF(ISOU.GT.3) THEN
        DELTIJ = 0.0D0
      ENDIF
C
C     S pour Source, V pour Variable
      THETS  = THETST(IPHAS)
      THETV  = THETAV(IVAR )
C
      IPCROO = IPCROM
      IF(ISTO2T(IPHAS).GT.0.AND.IROEXT(IPHAS).GT.0) THEN
        IPCROO = IPPROC(IROMA(IPHAS))
      ENDIF
      IPTSTA = 0
      IF(ISTO2T(IPHAS).GT.0) THEN
        IPTSTA = IPPROC(ITSTUA(IPHAS))
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
C 2. TERMES SOURCES  UTILISATEURS
C=======================================================================
C(le deuxieme argument GRDVIT est lu en PRODUC dans ustsri, mais ce
C tableau n'est dimensionne et utilise qu'en modele Rij standard)
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('RESSSG',IFINIA)
C
      CALL USTSRI
C     ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITPSMP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMCELP , GAMMA  , GRDVIT , GRDVIT ,
     &   SMBR   , ROVSDT ,
C        ------   ------
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , DAM    , DRTP   ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si on extrapole les T.S.
      IF(ISTO2T(IPHAS).GT.0) THEN
        DO IEL = 1, NCEL
C       Sauvegarde pour echange
          TUEXPR = PROPCE(IEL,IPTSTA+ISOU-1)
C       Pour la suite et le pas de temps suivant
          PROPCE(IEL,IPTSTA+ISOU-1) = SMBR(IEL)
C       Second membre du pas de temps precedent
C       On suppose -ROVSDT > 0 : on implicite
C          le terme source utilisateur (le reste)
          SMBR(IEL) = ROVSDT(IEL)*RTPA(IEL,IVAR)  - THETS*TUEXPR
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
C 2. TERMES SOURCES  LAGRANGIEN : COUPLAGE RETOUR
C=======================================================================
C
C     Ordre 2 non pris en compte
       IF (IILAGR.EQ.2 .AND. LTSDYN.EQ.1 .AND. IPHAS.EQ.ILPHAS) THEN
         DO IEL = 1,NCEL
           SMBR(IEL)   = SMBR(IEL)   + TSLAGE(IEL)
           ROVSDT(IEL) = ROVSDT(IEL) + MAX(-TSLAGI(IEL),ZERO)
         ENDDO
       ENDIF
C
C=======================================================================
C 3. TERME SOURCE DE MASSE
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
     & ( NCELET , NCEL   , NCESMP , IIUN   , ISTO2T(IPHAS) , THETV  ,
     &   ICETSM , ITPSMP ,
     &   VOLUME , RTPA(1,IVAR) , SMCELP , GAMMA  ,
     &   SMBR   ,  ROVSDT , W1 )
C
C       Si on extrapole les TS on met Gamma Pinj dans PROPCE
        IF(ISTO2T(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSTA+ISOU-1) =
     &      PROPCE(IEL,IPTSTA+ISOU-1) + W1(IEL)
          ENDDO
C       Sinon on le met directement dans SMBR
        ELSE
          DO IEL = 1, NCEL
            SMBR(IEL) = SMBR(IEL) + W1(IEL)
          ENDDO
        ENDIF
C
      ENDIF
C
C=======================================================================
C 4. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
C    ET TERME INSTATIONNAIRE
C=======================================================================
C
C ---> Calcul de mij
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &               IFACEL,IFABOR,PROPFA(1,IFLMAS),PROPFB(1,IFLMAB),W1)
C
C ---> Ajout au second membre
C
      DO IEL = 1, NCEL
        SMBR(IEL) = SMBR(IEL)
     &              + ICONV(IVAR)*W1(IEL)*RTPA(IEL,IVAR)
      ENDDO
C
C ---> Ajout dans la diagonale de la matrice
C     Extrapolation ou non, meme forme par coherence avec bilsc2
C
      DO IEL=1,NCEL
        ROVSDT(IEL) = ROVSDT(IEL)
     &            + ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &            - ICONV(IVAR)*W1(IEL)*THETV
      ENDDO
C
C
C=======================================================================
C 5. PRODUCTION, PHI1, PHI2, ET DISSIPATION
C=======================================================================
C
C ---> Terme source
C     -rho*epsilon*( Cs1*aij + Cs2*(aikajk -1/3*aijaij*deltaij))
C     -Cr1*P*aij + Cr2*rho*k*sij - Cr3*rho*k*sij*sqrt(aijaij)
C     +Cr4*rho*k(aik*sjk+ajk*sik-2/3*akl*skl*deltaij)
C     +Cr5*rho*k*(aik*rjk + ajk*rik)
C     -2/3*epsilon*deltaij
C
      IF(ISOU.EQ.1)THEN
        III = 1
        JJJ = 1
      ELSEIF(ISOU.EQ.2)THEN
        III = 2
        JJJ = 2
      ELSEIF(ISOU.EQ.3)THEN
        III = 3
        JJJ = 3
      ELSEIF(ISOU.EQ.4)THEN
        III = 1
        JJJ = 2
      ELSEIF(ISOU.EQ.5)THEN
        III = 1
        JJJ = 3
      ELSEIF(ISOU.EQ.6)THEN
        III = 2
        JJJ = 3
      ENDIF
C
      DO IEL=1,NCEL
C
        XPROD(1,1) = -2.0D0*(RTPA(IEL,IR11IP)*GRDVIT(IEL,1,1) +
     &                       RTPA(IEL,IR12IP)*GRDVIT(IEL,1,2) +
     &                       RTPA(IEL,IR13IP)*GRDVIT(IEL,1,3) )
        XPROD(1,2) = -(      RTPA(IEL,IR11IP)*GRDVIT(IEL,2,1) +
     &                       RTPA(IEL,IR12IP)*GRDVIT(IEL,2,2) +
     &                       RTPA(IEL,IR13IP)*GRDVIT(IEL,2,3) )
     &               -(      RTPA(IEL,IR12IP)*GRDVIT(IEL,1,1) +
     &                       RTPA(IEL,IR22IP)*GRDVIT(IEL,1,2) +
     &                       RTPA(IEL,IR23IP)*GRDVIT(IEL,1,3) )
        XPROD(1,3) = -(      RTPA(IEL,IR11IP)*GRDVIT(IEL,3,1) +
     &                       RTPA(IEL,IR12IP)*GRDVIT(IEL,3,2) +
     &                       RTPA(IEL,IR13IP)*GRDVIT(IEL,3,3) )
     &               -(      RTPA(IEL,IR13IP)*GRDVIT(IEL,1,1) +
     &                       RTPA(IEL,IR23IP)*GRDVIT(IEL,1,2) +
     &                       RTPA(IEL,IR33IP)*GRDVIT(IEL,1,3) )
        XPROD(2,2) = -2.0D0*(RTPA(IEL,IR12IP)*GRDVIT(IEL,2,1) +
     &                       RTPA(IEL,IR22IP)*GRDVIT(IEL,2,2) +
     &                       RTPA(IEL,IR23IP)*GRDVIT(IEL,2,3) )
        XPROD(2,3) = -(      RTPA(IEL,IR12IP)*GRDVIT(IEL,3,1) +
     &                       RTPA(IEL,IR22IP)*GRDVIT(IEL,3,2) +
     &                       RTPA(IEL,IR23IP)*GRDVIT(IEL,3,3) )
     &               -(      RTPA(IEL,IR13IP)*GRDVIT(IEL,2,1) +
     &                       RTPA(IEL,IR23IP)*GRDVIT(IEL,2,2) +
     &                       RTPA(IEL,IR33IP)*GRDVIT(IEL,2,3) )
        XPROD(3,3) = -2.0D0*(RTPA(IEL,IR13IP)*GRDVIT(IEL,3,1) +
     &                       RTPA(IEL,IR23IP)*GRDVIT(IEL,3,2) +
     &                       RTPA(IEL,IR33IP)*GRDVIT(IEL,3,3) )
        XPROD(2,1) = XPROD(1,2)
        XPROD(3,1) = XPROD(1,3)
        XPROD(3,2) = XPROD(2,3)
C
        TRPROD = 0.5D0 * (XPROD(1,1) + XPROD(2,2) + XPROD(3,3) )
        TRRIJ = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                   RTPA(IEL,IR33IP))
C-----> aII = aijaij
        AII    = 0.D0
        AKLSKL = 0.D0
        AIKSJK = 0.D0
        AIKRJK = 0.D0
        AIKAKJ = 0.D0
        XANISO(1,1) = RTPA(IEL,IR11IP)/TRRIJ-2.D0/3.D0
        XANISO(2,2) = RTPA(IEL,IR22IP)/TRRIJ-2.D0/3.D0
        XANISO(3,3) = RTPA(IEL,IR33IP)/TRRIJ-2.D0/3.D0
        XANISO(1,2) = RTPA(IEL,IR12IP)/TRRIJ
        XANISO(1,3) = RTPA(IEL,IR13IP)/TRRIJ
        XANISO(2,3) = RTPA(IEL,IR23IP)/TRRIJ
        XANISO(2,1) = XANISO(1,2)
        XANISO(3,1) = XANISO(1,3)
        XANISO(3,2) = XANISO(2,3)
C
        XSTRAI(1,1) = GRDVIT(IEL,1,1)
        XSTRAI(1,2) = 0.5D0*(GRDVIT(IEL,1,2)+GRDVIT(IEL,2,1))
        XSTRAI(1,3) = 0.5D0*(GRDVIT(IEL,1,3)+GRDVIT(IEL,3,1))
        XSTRAI(2,1) = XSTRAI(1,2)
        XSTRAI(2,2) = GRDVIT(IEL,2,2)
        XSTRAI(2,3) = 0.5D0*(GRDVIT(IEL,2,3)+GRDVIT(IEL,3,2))
        XSTRAI(3,1) = XSTRAI(1,3)
        XSTRAI(3,2) = XSTRAI(2,3)
        XSTRAI(3,3) = GRDVIT(IEL,3,3)
C
        XROTAC(1,1) = 0.D0
        XROTAC(1,2) = 0.5D0*(GRDVIT(IEL,1,2)-GRDVIT(IEL,2,1))
        XROTAC(1,3) = 0.5D0*(GRDVIT(IEL,1,3)-GRDVIT(IEL,3,1))
        XROTAC(2,1) = -XROTAC(1,2)
        XROTAC(2,2) = 0.D0
        XROTAC(2,3) = 0.5D0*(GRDVIT(IEL,2,3)-GRDVIT(IEL,3,2))
        XROTAC(3,1) = -XROTAC(1,3)
        XROTAC(3,2) = -XROTAC(2,3)
        XROTAC(3,3) = 0.D0
C
        DO II=1,3
         DO JJ = 1,3
           AII    = AII+XANISO(II,JJ)*XANISO(II,JJ)
           AKLSKL = AKLSKL + XANISO(II,JJ)*XSTRAI(II,JJ)
         ENDDO
        ENDDO
C
        DO KK = 1,3
           AIKSJK = AIKSJK + XANISO(III,KK)*XSTRAI(JJJ,KK)
     &              +XANISO(JJJ,KK)*XSTRAI(III,KK)
           AIKRJK = AIKRJK + XANISO(III,KK)*XROTAC(JJJ,KK)
     &              +XANISO(JJJ,KK)*XROTAC(III,KK)
           AIKAKJ = AIKAKJ + XANISO(III,KK)*XANISO(KK,JJJ)
        ENDDO
C
C     Si on extrapole les TS (rarissime), on met tout dans PROPCE.
C     On n'implicite pas le terme en Cs1*aij ni le terme en Cr1*P*aij.
C     Sinon, on met tout dans SMBR et on peut impliciter Cs1*aij
C     et Cr1*P*aij. Ici on stocke le second membre et le terme implicite
C     dans W1 et W2, pour eviter d'avoir un test IF(ISTO2T(IPHAS).GT.0)
C     dans la boucle NCEL
C     Dans le terme en W1, qui a vocation a etre extrapole, on utilise
C     naturellement IPCROO.
C     L'implicitation des deux termes pourrait se faire aussi en cas
C     d'extrapolation, en isolant ces deux termes et les mettant dans
C     SMBR et pas PROPCE et en utilisant IPCROM ... a modifier si le
C     besoin s'en fait vraiment sentir !
        W1(IEL) = PROPCE(IEL,IPCROO) * VOLUME(IEL) *(
     &       XPROD(III,JJJ) -RTPA(IEL,IEIPH)
     &       *(CSSGS1*XANISO(III,JJJ)+CSSGS2*(
     &       AIKAKJ -1.D0/3.D0*DELTIJ*AII))
     &       -CSSGR1*TRPROD*XANISO(III,JJJ)
     &       +TRRIJ*XSTRAI(III,JJJ)*(CSSGR2-CSSGR3*SQRT(AII))
     &       +CSSGR4*TRRIJ*(AIKSJK-2.D0/3.D0*DELTIJ*AKLSKL)
     &       +CSSGR5*TRRIJ*(AIKRJK)-2.D0/3.D0*RTPA(IEL,IEIPH)*DELTIJ)
C
        W2(IEL) = VOLUME(IEL)/TRRIJ*PROPCE(IEL,IPCROM)*(
     &         CSSGS1*RTPA(IEL,IEIPH) + CSSGR1*MAX(TRPROD,0.D0) )
C
      ENDDO
C
      IF(ISTO2T(IPHAS).GT.0) THEN
C
        DO IEL = 1, NCEL
          PROPCE(IEL,IPTSTA+ISOU-1) = PROPCE(IEL,IPTSTA+ISOU-1)
     &         + W1(IEL)
        ENDDO
C
      ELSE
C
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR(IEL) + W1(IEL)
          ROVSDT(IEL) = ROVSDT(IEL) + W2(IEL)
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 7. TERMES DE GRAVITE
C=======================================================================
C
      IF(IGRARI(IPHAS).EQ.1) THEN
C
        DO IEL = 1, NCEL
          W7(IEL) = 0.D0
        ENDDO
C
        CALL RIJTHE
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , GRAROX , GRAROY , GRAROZ , W7     ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si on extrapole les T.S. : PROPCE
      IF(ISTO2T(IPHAS).GT.0) THEN
        DO IEL = 1, NCEL
           PROPCE(IEL,IPTSTA+ISOU-1) =
     &     PROPCE(IEL,IPTSTA+ISOU-1) + W7(IEL)
         ENDDO
C     Sinon SMBR
       ELSE
         DO IEL = 1, NCEL
           SMBR(IEL) = SMBR(IEL) + W7(IEL)
         ENDDO
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C 9. TERMES DE DIFFUSION
C=======================================================================
C ---> Viscosite
C
      IF( IDIFF(IVAR).GE. 1 ) THEN
        DO IEL = 1, NCEL
          TRRIJ = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                   RTPA(IEL,IR33IP))
          RCTSE = 2.D0/3.D0 * PROPCE(IEL,IPCROM) * CSRIJ * TRRIJ**2
     &                  / RTPA(IEL,IEIPH)
          W1(IEL) = PROPCE(IEL,IPCVIS)
     &                            + IDIFFT(IVAR)*RCTSE
        ENDDO
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
C 10. RESOLUTION
C=======================================================================
C
      IF(ISTO2T(IPHAS).GT.0) THEN
        THETP1 = 1.D0 + THETS
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR(IEL) + THETP1*PROPCE(IEL,IPTSTA+ISOU-1)
        ENDDO
      ENDIF
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
CMO      IPP    =
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
     &   ROVSDT , SMBR   , RTP(1,IVAR)     ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C=======================================================================
C 11. IMPRESSIONS
C=======================================================================
C
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)
C
#else
C
 1000 FORMAT(/,'           SOLVING VARIABLE ',A8           ,/)
C
#endif
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
