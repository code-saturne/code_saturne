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
                        SUBROUTINE RESRIJ
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
     &   COEFA  , COEFB  , PRODUC , GRAROX , GRAROY , GRAROZ ,
     &   CKUPDC , SMCELP , GAMMA  ,
     &   VISCF  , VISCB  , COEFAX ,
     &   TSLAGE , TSLAGI ,
     &   DAM    , XAM    , DAG    , XAG    , DRTP   , SMBR   , ROVSDT ,
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
CFONC    POUR Rij (modele standard LRR)
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
CARGU ! PRODUC       ! TR !  -> ! TABLEAU DE TRAVAIL POUR PRODUCTION   !
CARGU !  (6,NCELET)  !    !     ! (SANS RHO VOLUME)                    !
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
      DOUBLE PRECISION PRODUC(6,NCELET)
      DOUBLE PRECISION GRAROX(NCELET), GRAROY(NCELET), GRAROZ(NCELET)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP)
      DOUBLE PRECISION SMCELP(NCESMP), GAMMA(NCESMP)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR), COEFAX(NFABOR)
      DOUBLE PRECISION TSLAGE(NCELET),TSLAGI(NCELET)
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
      INTEGER          IDEBIA, IDEBRA
      INTEGER          INIT  , IFAC  , IEL   , INC   , ICCOCG
      INTEGER          II    , JJ    , IIUN
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEIPH , IUIPH
      INTEGER          IPCROM, IPCVIS, IFLMAS, IFLMAB, IPCROO
      INTEGER          ICLVAR, ICLVAF
      INTEGER          NSWRGP, IMLIGP, IWARNP, IPHYDP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP
      INTEGER          NITMAP, NSWRSP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          IMGRP , NCYMXP, NITMFP
      INTEGER          IDIMTE, ITENSO
      INTEGER          IPTSTA
      INTEGER          ISOLUC
      DOUBLE PRECISION BLENCP, EPSILP, EPSRGP, CLIMGP, EXTRAP, RELAXP
      DOUBLE PRECISION TRPROD, TRRIJ , CSTRIJ, RCTSE , DELTIJ
      DOUBLE PRECISION GRDPX , GRDPY , GRDPZ , GRDSN
      DOUBLE PRECISION SURFN2
      DOUBLE PRECISION TUEXPR, THETS , THETV , THETP1
      DOUBLE PRECISION D1S3  , D2S3
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
      D1S3 = 1.D0/3.D0
      D2S3 = 2.D0/3.D0
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
C(le premier argument PRODUC est lu en GRDVIT dans ustsri, mais ce
C tableau n'est dimensionne et utilise qu'en modele Rij SSG)
C
      CALL USTSRI
C     ===========
     & ( IDEBIA , IDEBRA ,
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
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMCELP , GAMMA  , PRODUC , PRODUC ,
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
C
C ---> Calcul de k pour la suite du sous-programme
C       on utilise un tableau de travail puisqu'il y en a...
      DO IEL = 1, NCEL
        W8(IEL) = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                     RTPA(IEL,IR33IP))
      ENDDO
C
C ---> Terme source
C
C      (1-CRIJ2) Pij (pour toutes les composantes de Rij)
C
C      DELTAIJ*(2/3.CRIJ2.P+2/3.CRIJ1.EPSILON)
C                    (termes diagonaux pour R11, R22 et R33)
C
C      -DELTAIJ*2/3*EPSILON
C
C     Si on extrapole les TS
C       On modifie la partie implicite :
C         Dans PHI1, on ne prendra que RHO CRIJ1 EPSILON/K et non pas
C                                  RHO CRIJ1 EPSILON/K (1-2/3 DELTAIJ)
C         Cela permet de conserver k^n au lieu de (R11^(n+1)+R22^n+R33^n)
C         Ce choix est discutable. C'est la solution ISOLUC = 1
C       Si on veut tout prendre en implicite (comme c'est fait
C         en ordre 1 std), c'est la solution ISOLUC = 2
C       -> a tester plus precisement si necessaire
C
C
C     Si on extrapole les TS
      IF(ISTO2T(IPHAS).GT.0) THEN
C
        ISOLUC = 1
C
        DO IEL = 1, NCEL
C
C     Demi-traces de Prod et R
          TRPROD = 0.5D0*(PRODUC(1,IEL)+PRODUC(2,IEL)+PRODUC(3,IEL))
          TRRIJ  = W8(IEL)
C
C     Calcul de Prod+Phi1+Phi2-Eps
C       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
C       Dans PROPCE :
C       = rhoPij-C1rho eps/k(   -2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
C       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij           }
          PROPCE(IEL,IPTSTA+ISOU-1) = PROPCE(IEL,IPTSTA+ISOU-1)
     &                          + PROPCE(IEL,IPCROO) * VOLUME(IEL)
     &      *(   DELTIJ*D2S3*
     &           (  CRIJ2*TRPROD
     &            +(CRIJ1-1.D0)* RTPA(IEL,IEIPH)  )
     &         +(1.0D0-CRIJ2)*PRODUC(ISOU,IEL)               )
C       Dans SMBR
C       =       -C1rho eps/k(Rij         )
C       = rho{                                     -C1eps/kRij}
          SMBR(IEL) = SMBR(IEL) + PROPCE(IEL,IPCROM) * VOLUME(IEL)
     &      *( -CRIJ1*RTPA(IEL,IEIPH)/TRRIJ * RTPA(IEL,IVAR)  )
C
C     Calcul de la partie implicite issue de Phi1
C       = C1rho eps/k(1        )
          ROVSDT(IEL) = ROVSDT(IEL) + PROPCE(IEL,IPCROM) * VOLUME(IEL)
     &                            *CRIJ1*RTPA(IEL,IEIPH)/TRRIJ*THETV
C
        ENDDO
C
C     Si on veut impliciter un bout de -C1rho eps/k(   -2/3k dij)
        IF(ISOLUC.EQ.2) THEN
C
          DO IEL = 1, NCEL
C
            TRRIJ  = W8(IEL)
C
C    On enleve a PROPCE (avec IPCROO)
C       =       -C1rho eps/k(   -1/3Rij dij)
            PROPCE(IEL,IPTSTA+ISOU-1) = PROPCE(IEL,IPTSTA+ISOU-1)
     &                          - PROPCE(IEL,IPCROO) * VOLUME(IEL)
     &      *(DELTIJ*D1S3*CRIJ1*RTPA(IEL,IEIPH)/TRRIJ * RTPA(IEL,IVAR))
C    On ajoute a SMBR (avec IPCROM)
C       =       -C1rho eps/k(   -1/3Rij dij)
            SMBR(IEL)                 = SMBR(IEL)
     &                          + PROPCE(IEL,IPCROM) * VOLUME(IEL)
     &      *(DELTIJ*D1S3*CRIJ1*RTPA(IEL,IEIPH)/TRRIJ * RTPA(IEL,IVAR))
C    On ajoute a ROVSDT (avec IPCROM)
C       =        C1rho eps/k(   -1/3    dij)
            ROVSDT(IEL) = ROVSDT(IEL) + PROPCE(IEL,IPCROM) * VOLUME(IEL)
     &      *(DELTIJ*D1S3*CRIJ1*RTPA(IEL,IEIPH)/TRRIJ                 )
          ENDDO
C
        ENDIF
C
C     Si on n'extrapole pas les termes sources
      ELSE
C
        DO IEL = 1, NCEL
C
C     Demi-traces de Prod et R
          TRPROD = 0.5D0*(PRODUC(1,IEL)+PRODUC(2,IEL)+PRODUC(3,IEL))
          TRRIJ  = W8(IEL)
C
C     Calcul de Prod+Phi1+Phi2-Eps
C       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
C       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij-C1eps/kRij}
          SMBR(IEL) = SMBR(IEL) + PROPCE(IEL,IPCROM) * VOLUME(IEL)
     &      *(   DELTIJ*D2S3*
     &           (  CRIJ2*TRPROD
     &            +(CRIJ1-1.D0)* RTPA(IEL,IEIPH)  )
     &         +(1.0D0-CRIJ2)*PRODUC(ISOU,IEL)
     &         -CRIJ1*RTPA(IEL,IEIPH)/TRRIJ * RTPA(IEL,IVAR)  )
C
C     Calcul de la partie implicite issue de Phi1
C       = C1rho eps/k(1-1/3 dij)
          ROVSDT(IEL) = ROVSDT(IEL) + PROPCE(IEL,IPCROM) * VOLUME(IEL)
     &         *(1.D0-D1S3*DELTIJ)*CRIJ1*RTPA(IEL,IEIPH)/TRRIJ
        ENDDO
C
      ENDIF
C
C=======================================================================
C 6. TERMES D'ECHO DE PAROI
C=======================================================================
C
      IF(IRIJEC(IPHAS).EQ.1) THEN
C
        DO IEL = 1, NCEL
          W7(IEL) = 0.D0
        ENDDO
C
        CALL RIJECH
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
     &   COEFA  , COEFB  , PRODUC , W7   ,
     &   COEFAX , VISCB  ,
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
C 8. TERMES DE DIFFUSION  A.grad(Rij) : PARTIE EXTRADIAGONALE EXPLICITE
C=======================================================================
C
C ---> Calcul du grad(Rij)
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
C     ===========
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
C ---> Calcul des termes extradiagonaux de A.grad(Rij)
C
       DO IEL = 1, NCEL
        TRRIJ = W8(IEL)
        CSTRIJ = PROPCE(IEL,IPCROO) * CSRIJ *TRRIJ / RTPA(IEL,IEIPH)
        W4(IEL) = CSTRIJ * ( RTPA(IEL,IR12IP) * W2(IEL)
     &                      +RTPA(IEL,IR13IP) * W3(IEL) )
        W5(IEL) = CSTRIJ * ( RTPA(IEL,IR12IP) * W1(IEL)
     &                      +RTPA(IEL,IR23IP) * W3(IEL) )
        W6(IEL) = CSTRIJ * ( RTPA(IEL,IR13IP) * W1(IEL)
     &                      +RTPA(IEL,IR23IP) * W2(IEL) )
       ENDDO
C
C
C ---> Assemblage de { A.grad(Rij) } .S aux faces
C
       CALL VECTDS
C      ===========
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
     &    PROPCE(IEL,IPTSTA+ISOU-1) + W4(IEL)
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
C 9. TERMES DE DIFFUSION  A.grad(Rij) : PARTIE DIAGONALE
C=======================================================================
C     Implicitation de (grad(Rij).n)n en gradient facette
C     Si IDIFRE=1, terme correctif explicite
C        grad(Rij)-(grad(Rij).n)n calcule en gradient cellule
C     Les termes de bord sont uniquement pris en compte dans la partie
C        en (grad(Rij).n)n
c     (W1,W2,W3) contient toujours le gradient de la variable traitee
C
C     Attention en periodicite on percom-ise le gradient comme si c'etait
C       un vecteur (alors que dans grdcel on l'a fait comme si c'etait
C       un tenseur ...).
C     A modifier eventuellement. Pour le moment on conserve donc
C       PERCOM et PARCOM.
C
      IF (IDIFRE(IPHAS).EQ.1) THEN
C
        DO IEL = 1, NCEL
          TRRIJ = W8(IEL)
          CSTRIJ = PROPCE(IEL,IPCROO) * CSRIJ *TRRIJ / RTPA(IEL,IEIPH)
          W4(IEL) = CSTRIJ*RTPA(IEL,IR11IP)
          W5(IEL) = CSTRIJ*RTPA(IEL,IR22IP)
          W6(IEL) = CSTRIJ*RTPA(IEL,IR33IP)
        ENDDO
C
C --->  TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (W1)
C         ===========
          CALL PARCOM (W2)
C         ===========
          CALL PARCOM (W3)
C         ===========
          CALL PARCOM (W4)
C         ===========
          CALL PARCOM (W5)
C         ===========
          CALL PARCOM (W6)
C         ===========
        ENDIF
C
C -->   TRAITEMENT DE LA PERIODICITE
C ->    IL RESTE DES DOUTES LA-DESSUS
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 1
          ITENSO = 0
          CALL PERCOM
C         ===========
     &    ( IDIMTE , ITENSO ,
     &      W1     , W1     , W1    ,
     &      W2     , W2     , W2    ,
     &      W3     , W3     , W3    )
C
          IDIMTE = 21
          ITENSO = 0
          CALL PERCOM
C         ===========
     &    ( IDIMTE , ITENSO ,
     &      W4     , W4     , W4    ,
     &      W5     , W5     , W5    ,
     &      W6     , W6     , W6    )
C
        ENDIF
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          SURFN2 = RA(ISRFAN-1+IFAC)**2
C
          GRDPX = 0.5D0*(W1(II)+W1(JJ))
          GRDPY = 0.5D0*(W2(II)+W2(JJ))
          GRDPZ = 0.5D0*(W3(II)+W3(JJ))
          GRDSN = GRDPX*SURFAC(1,IFAC)+GRDPY*SURFAC(2,IFAC)
     &           +GRDPZ*SURFAC(3,IFAC)
          GRDPX = GRDPX-GRDSN*SURFAC(1,IFAC)/SURFN2
          GRDPY = GRDPY-GRDSN*SURFAC(2,IFAC)/SURFN2
          GRDPZ = GRDPZ-GRDSN*SURFAC(3,IFAC)/SURFN2
C
          VISCF(IFAC)= 0.5D0*(
     &          (W4(II)+W4(JJ))*GRDPX*SURFAC(1,IFAC)
     &         +(W5(II)+W5(JJ))*GRDPY*SURFAC(2,IFAC)
     &         +(W6(II)+W6(JJ))*GRDPZ*SURFAC(3,IFAC))
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
     &      PROPCE(IEL,IPTSTA+ISOU-1) + W1(IEL)
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
C
C ---> Viscosite orthotrope pour partie implicite
C
      IF( IDIFF(IVAR).GE. 1 ) THEN
        DO IEL = 1, NCEL
          TRRIJ = W8(IEL)
          RCTSE = PROPCE(IEL,IPCROM) * CSRIJ * TRRIJ / RTPA(IEL,IEIPH)
          W1(IEL) = PROPCE(IEL,IPCVIS)
     &                            + IDIFFT(IVAR)*RCTSE*RTPA(IEL,IR11IP)
          W2(IEL) = PROPCE(IEL,IPCVIS)
     &                            + IDIFFT(IVAR)*RCTSE*RTPA(IEL,IR22IP)
          W3(IEL) = PROPCE(IEL,IPCVIS)
     &                            + IDIFFT(IVAR)*RCTSE*RTPA(IEL,IR33IP)
        ENDDO
C
        CALL VISORT
C       ===========
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
