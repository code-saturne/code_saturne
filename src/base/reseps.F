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
                        SUBROUTINE RESEPS
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
     &   ICEPDC , ICETSM , ITPSMP , IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , GRDVIT , PRODUC , GRAROX , GRAROY , GRAROZ ,
     &   CKUPDC , SMCELP , GAMMA  ,
     &   VISCF  , VISCB  ,
     &   TSLAGR ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
CFONC    POUR EPSILON DANS LE CAS DU MODELE Rij-epsilon
CFONC
CFONC    On a ici ISOU   = 7 (7 ieme variable du Rij-epsilon)
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
CARGU ! GRDVIT       ! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME GRAD   !
CARGU !  (NCELET,3,3)!    !     !    DE VITESSE     UNIQT POUR ITURB=31!
CARGU ! PRODUC       ! TR !  -> ! TABLEAU DE TRAVAIL POUR PRODUCTION   !
CARGU !  (6,NCELET)  !    !     ! (SANS RHO VOLUME) UNIQT POUR ITURB=30!
CARGU ! GRAROX,Y,Z   ! TR !  -> ! TABLEAU DE TRAVAIL POUR GRAD ROM     !
CARGU !  (NCELET)    !    !     !                                      !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMCELP(NCESMP! TR !  -> ! VALEUR DE LA VARIABLE ASSOCIEE A LA  !
CARGU !              !    !     !  SOURCE DE MASSE                     !
CARGU ! GAMMA(NCESMP)! TR !  -> ! VALEUR DU FLUX DE MASSE              !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !  (NCELET,*)  !    !     !   LAGRANGIEN                         !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
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
      DOUBLE PRECISION PRODUC(6,NCELET), GRDVIT(NCELET,3,3)
      DOUBLE PRECISION GRAROX(NCELET), GRAROY(NCELET), GRAROZ(NCELET)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP)
      DOUBLE PRECISION SMCELP(NCESMP), GAMMA(NCESMP)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          INIT  , IFAC  , IEL   , INC   , ICCOCG
      INTEGER          II    ,JJ     , IIUN
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEIPH , IUIPH
      INTEGER          ICLVAR, ICLVAF
      INTEGER          IPCROM, IPCROO, IPCVIS, IPCVST, IFLMAS, IFLMAB
      INTEGER          NSWRGP, IMLIGP, IWARNP, IPHYDP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP
      INTEGER          NITMAP, NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IDIMTE, ITENSO
      INTEGER          IPTSTA
      INTEGER          MAXELT, ILS
      DOUBLE PRECISION BLENCP, EPSILP, EPSRGP, CLIMGP, EXTRAP, RELAXP
      DOUBLE PRECISION TRPROD , TRRIJ ,CSTEPS, RCTSE
      DOUBLE PRECISION GRDPX,GRDPY,GRDPZ,GRDSN
      DOUBLE PRECISION SURFN2
      DOUBLE PRECISION TSEPS , KSEPS , CEPS2
      DOUBLE PRECISION TUEXPE, THETS , THETV , THETAP, THETP1
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
      IUIPH  = IU(IPHAS)
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
      IPCVST = IPPROC(IVISCT(IPHAS))
      IFLMAS = IPPROF(IFLUMA(IUIPH))
      IFLMAB = IPPROB(IFLUMA(IUIPH))
C
      ICLVAR = ICLRTP(IVAR,ICOEF)
      ICLVAF = ICLRTP(IVAR,ICOEFF)
C
C     Constante Ce2, qui vaut CE2 pour ITURB=30 et CSSGE2 pour ITRUB=31
      IF (ITURB(IPHAS).EQ.30) THEN
        CEPS2 = CE2
      ELSE
        CEPS2 = CSSGE2
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
      IF(ISTO2T(IPHAS).GT.0) THEN
        IPTSTA = IPPROC(ITSTUA(IPHAS))
      ELSE
        IPTSTA = 0
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
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('RESEPS',IFINIA)
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
     &   COEFA  , COEFB  , CKUPDC , SMCELP , GAMMA  , GRDVIT , PRODUC ,
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
C 3.  TERMES SOURCES LAGRANGIEN : COUPLAGE RETOUR
C=======================================================================
C
C     Ordre 2 non pris en compte
       IF (IILAGR.EQ.2 .AND. LTSDYN.EQ.1 .AND. IPHAS.EQ.ILPHAS) THEN
C
         DO IEL = 1,NCEL
C Ts sur eps
           TSEPS = -0.5D0 * ( TSLAGR(IEL,ITSR11)
     &                      + TSLAGR(IEL,ITSR22)
     &                      + TSLAGR(IEL,ITSR33) )
C rapport k/eps
           KSEPS = 0.5D0 * ( RTPA(IEL,IR11IP)
     &                     + RTPA(IEL,IR22IP)
     &                     + RTPA(IEL,IR33IP) )
     &                     / RTPA(IEL,IEIPH)
C
           SMBR(IEL)   = SMBR(IEL) + CE4 *TSEPS *RTPA(IEL,IEIPH) /KSEPS
           ROVSDT(IEL) = ROVSDT(IEL) + MAX( (-CE4*TSEPS/KSEPS) , ZERO)
         ENDDO
C
       ENDIF
C
C=======================================================================
C 4. TERME SOURCE DE MASSE
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
C 5. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
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
      DO IEL = 1, NCEL
        ROVSDT(IEL) = ROVSDT(IEL)
     &           + ISTAT(IVAR)*(PROPCE(IEL,IPCROM)/DT(IEL))*VOLUME(IEL)
     &           - ICONV(IVAR)*W1(IEL)*THETV
      ENDDO
C
C
C=======================================================================
C 6. PRODUCTION RHO * Ce1 * epsilon / k * P
C    DISSIPATION RHO*Ce2.epsilon/k*epsilon
C=======================================================================
C
C
C ---> Calcul de k pour la suite du sous-programme
C       on utilise un tableau de travail puisqu'il y en a...
      DO IEL = 1, NCEL
        W8(IEL) = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                     RTPA(IEL,IR33IP))
      ENDDO
C ---> Calcul de la trace de la production, suivant qu'on est en
C     Rij standard ou en SSG (utilisation de PRODUC ou GRDVIT)
      IF (ITURB(IPHAS).EQ.30) THEN
        DO IEL = 1, NCEL
          W9(IEL) = 0.5D0*(PRODUC(1,IEL)+PRODUC(2,IEL)+PRODUC(3,IEL))
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          W9(IEL) = -( RTPA(IEL,IR11IP)*GRDVIT(IEL,1,1) +
     &                 RTPA(IEL,IR12IP)*GRDVIT(IEL,1,2) +
     &                 RTPA(IEL,IR13IP)*GRDVIT(IEL,1,3) +
     &                 RTPA(IEL,IR12IP)*GRDVIT(IEL,2,1) +
     &                 RTPA(IEL,IR22IP)*GRDVIT(IEL,2,2) +
     &                 RTPA(IEL,IR23IP)*GRDVIT(IEL,2,3) +
     &                 RTPA(IEL,IR13IP)*GRDVIT(IEL,3,1) +
     &                 RTPA(IEL,IR23IP)*GRDVIT(IEL,3,2) +
     &                 RTPA(IEL,IR33IP)*GRDVIT(IEL,3,3) )
        ENDDO
      ENDIF
C
C
C     Terme explicite (Production)
C
      DO IEL = 1, NCEL
C     Demi-traces
        TRPROD = W9(IEL)
        TRRIJ  = W8(IEL)
        W1(IEL)   =             PROPCE(IEL,IPCROO)*VOLUME(IEL)*
     &       CE1*RTPA(IEL,IEIPH)/TRRIJ*TRPROD
      ENDDO
C
C     Si on extrapole les T.S : PROPCE
      IF(ISTO2T(IPHAS).GT.0) THEN
        DO IEL = 1, NCEL
          PROPCE(IEL,IPTSTA+ISOU-1) =
     &    PROPCE(IEL,IPTSTA+ISOU-1) + W1(IEL)
        ENDDO
C     Sinon : SMBR
      ELSE
        DO IEL = 1, NCEL
          SMBR(IEL) = SMBR(IEL) + W1(IEL)
        ENDDO
      ENDIF
C
C     Terme implicite (Dissipation)
      DO IEL = 1, NCEL
        TRRIJ  = W8(IEL)
        SMBR(IEL) = SMBR(IEL) - PROPCE(IEL,IPCROM)*VOLUME(IEL)*
     &                         CEPS2*RTPA(IEL,IEIPH)**2/TRRIJ
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
        TRRIJ  = W8(IEL)
        ROVSDT(IEL) = ROVSDT(IEL) + CEPS2*PROPCE(IEL,IPCROM)*VOLUME(IEL)
     &                     *RTPA(IEL,IEIPH)/TRRIJ*THETAP
      ENDDO
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
C 8.a TERMES DE DIFFUSION  A.grad(Eps) : PARTIE EXTRADIAGONALE EXPLICITE
C     RIJ STANDARD
C=======================================================================
C
      IF (ITURB(IPHAS).EQ.30) THEN
C
C ---> Calcul du grad(Eps)
C
C
        ICCOCG = 1
        INC = 1
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
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IVAR )   , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C ---> Calcul des termes extradiagonaux de A.grad(Eps)
C
        DO IEL = 1 , NCEL
          TRRIJ  = W8(IEL)
          CSTEPS  = PROPCE(IEL,IPCROO) * CRIJEP *TRRIJ / RTPA(IEL,IEIPH)
          W4(IEL) = CSTEPS * ( RTPA(IEL,IR12IP) * W2(IEL) +
     &         RTPA(IEL,IR13IP) * W3(IEL) )
          W5(IEL) = CSTEPS * ( RTPA(IEL,IR12IP) * W1(IEL) +
     &         RTPA(IEL,IR23IP) * W3(IEL) )
          W6(IEL) = CSTEPS * ( RTPA(IEL,IR13IP) * W1(IEL) +
     &         RTPA(IEL,IR23IP) * W2(IEL) )
        ENDDO
C
C ---> Assemblage de { A.grad(Eps) } .S aux faces
C
        CALL VECTDS
C       ===========
     & (NDIM   , NCELET , NCEL   , NFAC   , NFABOR   ,
     &  IFACEL , IFABOR , IA     ,
     &  SURFAC , SURFBO , RA(IPOND) ,
     &  W4     , W5     , W6     ,
     &  VISCF  , VISCB  , RA     )
C
        INIT = 1
        CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                                   IFACEL,IFABOR,VISCF,VISCB,W4)
C
C     Si on extrapole les termes sources
        IF(ISTO2T(IPHAS).GT.0) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPTSTA+ISOU-1) =
     &           PROPCE(IEL,IPTSTA+ISOU-1) + W4(IEL)
          ENDDO
C     Sinon
        ELSE
          DO IEL = 1, NCEL
            SMBR(IEL) = SMBR(IEL) + W4(IEL)
          ENDDO
        ENDIF
C
C
C=======================================================================
C 8.b TERMES DE DIFFUSION  A.grad(Eps) : PARTIE DIAGONALE
C     RIJ STANDARD
C=======================================================================
C     Implicitation de (grad(eps).n)n en gradient facette
C     Si IDIFRE=1, terme correctif explicite
C        grad(eps)-(grad(eps).n)n calcule en gradient cellule
C     Les termes de bord sont uniquement pris en compte dans la partie
C        en (grad(eps).n)n
c     (W1,W2,W3) contient toujours le gradient de la variable traitee
C
C     La parcom/percom-isation du gradient de epsilon a ete faite dans
C       grdcel. Pas utile de recommencer.
C
        IF (IDIFRE(IPHAS).EQ.1) THEN
C
          DO IEL = 1, NCEL
            TRRIJ  = W8(IEL)
            CSTEPS = PROPCE(IEL,IPCROO) * CRIJEP *TRRIJ/RTPA(IEL,IEIPH)
            W4(IEL)=CSTEPS*RTPA(IEL,IR11IP)
            W5(IEL)=CSTEPS*RTPA(IEL,IR22IP)
            W6(IEL)=CSTEPS*RTPA(IEL,IR33IP)
          ENDDO
C
C --->  TRAITEMENT DU PARALLELISME (MEMES DOUTES QUE PERIODICITE ?)
C
          IF(IRANGP.GE.0) THEN
            CALL PARCOM (W4)
C           ===========
            CALL PARCOM (W5)
C           ===========
            CALL PARCOM (W6)
C           ===========
          ENDIF
C
C -->   TRAITEMENT DE LA PERIODICITE
          IF(IPERIO.EQ.1) THEN
C
            IDIMTE = 21
            ITENSO = 0
C
            CALL PERCOM
C           ===========
     &    ( IDIMTE , ITENSO ,
     &      W4     , W4     , W4    ,
     &      W5     , W5     , W5    ,
     &      W6     , W6     , W6    )
C
          ENDIF
C
          DO IFAC = 1, NFAC
C
            II=IFACEL(1,IFAC)
            JJ=IFACEL(2,IFAC)
C
            SURFN2 = RA(ISRFAN-1+IFAC)**2
C
            GRDPX=0.5D0*(W1(II)+W1(JJ))
            GRDPY=0.5D0*(W2(II)+W2(JJ))
            GRDPZ=0.5D0*(W3(II)+W3(JJ))
            GRDSN=GRDPX*SURFAC(1,IFAC)+GRDPY*SURFAC(2,IFAC)
     &           +GRDPZ*SURFAC(3,IFAC)
            GRDPX=GRDPX-GRDSN*SURFAC(1,IFAC)/SURFN2
            GRDPY=GRDPY-GRDSN*SURFAC(2,IFAC)/SURFN2
            GRDPZ=GRDPZ-GRDSN*SURFAC(3,IFAC)/SURFN2
C
            VISCF(IFAC)= 0.5D0*(
     &           (W4(II)+W4(JJ))*GRDPX*SURFAC(1,IFAC)
     &           +(W5(II)+W5(JJ))*GRDPY*SURFAC(2,IFAC)
     &           +(W6(II)+W6(JJ))*GRDPZ*SURFAC(3,IFAC))
C
          ENDDO
C
          DO IFAC = 1, NFABOR
            VISCB(IFAC) = 0.D0
          ENDDO
C
          INIT = 1
          CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &       IFACEL,IFABOR,VISCF,VISCB,W1)
C
C     Si on extrapole les termes sources
          IF(ISTO2T(IPHAS).GT.0) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IPTSTA+ISOU-1) =
     &             PROPCE(IEL,IPTSTA+ISOU-1) + W1(IEL)
            ENDDO
C     Sinon
          ELSE
            DO IEL = 1, NCEL
              SMBR(IEL) = SMBR(IEL) + W1(IEL)
            ENDDO
          ENDIF
C
        ENDIF
C
C ---> Viscosite orthotrope pour partie implicite
C
        IF( IDIFF(IVAR).GE. 1 ) THEN
          DO IEL = 1, NCEL
            TRRIJ  = W8(IEL)
            RCTSE = PROPCE(IEL,IPCROM) * CRIJEP * TRRIJ/RTPA(IEL,IEIPH)
            W1(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*RCTSE*RTPA(IEL,IR11IP)
            W2(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*RCTSE*RTPA(IEL,IR22IP)
            W3(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*RCTSE*RTPA(IEL,IR33IP)
          ENDDO
C
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
     &   W1     , W2     , W3     ,
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
      ELSE
C=======================================================================
C 8.c TERMES DE DIFFUSION
C     RIJ SSG
C=======================================================================
C ---> Viscosite
C
        IF( IDIFF(IVAR).GE. 1 ) THEN
          DO IEL = 1, NCEL
            W1(IEL) = PROPCE(IEL,IPCVIS)
     &           + IDIFFT(IVAR)*PROPCE(IEL,IPCVST)/SIGMAE
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
      ENDIF
C
C
C=======================================================================
C 9. RESOLUTION
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
      NDIRCP = NDIRCL(IVAR)
      IRESLP = IRESOL(IVAR)
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
C 10. IMPRESSIONS
C=======================================================================
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)
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
