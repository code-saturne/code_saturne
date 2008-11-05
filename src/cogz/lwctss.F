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
                        SUBROUTINE LWCTSS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IZFPPP , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBRS  , ROVSDT ,
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    , W11    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME PREMELANGE MODELE LWC
CFONC    ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
CFONC    SUR UN PAS DE TEMPS
CFONC
CFONC ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
CFONC ---------   DE CELUI DE USTSSC.F
CFONC
CFONC ON RESOUT ROVSDT*D(VAR) = SMBRS
CFONC
CFONC ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
CFONC   UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
CFONC   ECRASER
CFONC
CFONC POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
CFONC   QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
CFONC   SMBRS
CFONC
CFONC DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
CFONC ECRIRE :
CFONC           SMBRS  = SMBRS  + CEXP + CIMP*VAR
CFONC           ROVSDT = ROVSDT + MAX(-CIMP,ZERO)
CFONC
CFONC ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
CFONC     SMBRS en kg variable/s :
CFONC      ex : pour la vitesse            kg m/s2
CFONC           pour les temperatures      kg degres/s
CFONC           pour les enthalpies        Joules/s
CFONC     ROVSDT en kg /s
CFONC
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
CARGU ! IZFPPP       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE PHYS. PART.          !
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
CARGU ! SMBRS(NCELET)! TR ! <-  ! SECOND MEMBRE EXPLICITE              !
CARGU ! ROVSDT(NCELET! TR ! <-  ! PARTIE DIAGONALE IMPLICITE           !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! XAM(NFAC,2)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..11(NCELET! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
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
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
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
      INTEGER          IZFPPP(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
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
      DOUBLE PRECISION SMBRS(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION XAM(NFAC,2)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR, IEL, IPHAS, IDIRAC, IVAR0
      INTEGER          IPHYDP , ITENSO , IDIMTE
      INTEGER          INC , ICCOCG
      INTEGER          IPCVST
      INTEGER          IPCROM, II
C
      INTEGER          IPTSCL(NDRACM), IPFMAL(NDRACM)
      INTEGER          IPFMEL(NDRACM), IPRHOL(NDRACM)
      DOUBLE PRECISION SUM, EPSI
      DOUBLE PRECISION TSGRAD, TSCHIM, TSDISS
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
      EPSI   = 1.0D-10
C
C --- Numero du scalaire a traiter : ISCAL
C
C --- Numero de la variable associee au scalaire a traiter ISCAL
      IVAR = ISCA(ISCAL)
C
C --- Numero de phase associee au scalaire ISCAL
      IPHAS = IPHSCA(ISCAL)
C
C ---
      IPCROM = IPPROC(IROM(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C
C --- Numero des grandeurs physiques (voir usclim)
      DO IDIRAC = 1, NDIRAC
        IPTSCL(IDIRAC) = IPPROC(ITSCL(IDIRAC))
        IPFMAL(IDIRAC) = IPPROC(IFMAL(IDIRAC))
        IPFMEL(IDIRAC) = IPPROC(IFMEL(IDIRAC))
        IPRHOL(IDIRAC) = IPPROC(IRHOL(IDIRAC))
      ENDDO
C
C=======================================================================
C 2. PRISE EN COMPTE DES TERMES SOURCES
C=======================================================================
C
      IF ( IVAR.EQ.ISCA(IYFM) ) THEN
C
C ---> Terme source pour la fraction massique moyenne de fuel
C
        DO IEL = 1, NCEL
            SUM = ZERO
            DO IDIRAC = 1, NDIRAC
              SUM  = SUM + PROPCE(IEL,IPRHOL(IDIRAC))
     &           *PROPCE(IEL,IPTSCL(IDIRAC))*VOLUME(IEL)
            ENDDO
C
C terme implicite
C
            IF (RTPA(IEL,IVAR).GT.EPSI) THEN
              ROVSDT(IEL) = ROVSDT(IEL) + MAX(-SUM/RTPA(IEL,IVAR),ZERO)
            ENDIF
C
C terme explicite
C
             SMBRS(IEL) =  SMBRS(IEL) + SUM
C
        ENDDO
C
      ENDIF
C
C ---> Terme source pour la variance de la fraction massique moyenne de fuel
C
      IF (IVAR.EQ.ISCA(IYFP2M)) THEN
C
        DO IEL = 1, NCEL
          SUM = ZERO
          DO IDIRAC = 1, NDIRAC
            SUM  = SUM + (PROPCE(IEL,IPTSCL(IDIRAC))*VOLUME(IEL)
     &        *(PROPCE(IEL,IPFMAL(IDIRAC)) - RTPA(IEL,ISCA(IYFM)))
     &             *PROPCE(IEL,IPRHOL(IDIRAC)))
          ENDDO
          SMBRS(IEL) = SMBRS(IEL) + SUM
        ENDDO
C
      ENDIF
C
C ---> Terme source pour la covariance
C
      IF ( IVAR.EQ.ISCA(ICOYFP)) THEN
C
C --- Calcul du gradient de F
C     =======================
C
        II = ISCA(IFM)
        DO IEL = 1, NCEL
          W10(IEL) = RTPA(IEL,II)
        ENDDO
C
C En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
        IF(IRANGP.GE.0) THEN
          CALL PARCOM(W10)
C         ===========
        ENDIF
C
C    Periodique
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    W10     , W10     , W10     ,
     &    W10     , W10     , W10     ,
     &    W10     , W10     , W10     )
        ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
        IVAR0 = 0
        IPHYDP = 0
        INC = 1
        ICCOCG = 1
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGR(II) , IMLIGR(II) ,
     &   IPHYDP , IWARNI(II) , NFECRA ,
     &   EPSRGR(II) , CLIMGR(II) , EXTRAG(II) ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W10    , W10    , W10    ,
     &   W10    , COEFA(1,ICLRTP(II,ICOEF))  ,
     &            COEFB(1,ICLRTP(II,ICOEF))  ,
     &   W1              , W2              , W3     ,
C        d./dx1          , d./dx2          , d./dx3 ,
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C --- Calcul du gradient de Yfuel
C     ===========================
C
        II = ISCA(IYFM)
        DO IEL = 1, NCEL
          W11(IEL) = RTPA(IEL,II)
        ENDDO

C En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
        IF(IRANGP.GE.0) THEN
          CALL PARCOM(W11)
C         ===========
        ENDIF
C
C    Periodique
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    W11     , W11     , W11     ,
     &    W11     , W11     , W11     ,
     &    W11     , W11     , W11     )
        ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
        IVAR0 = 0
        IPHYDP = 0
        INC = 1
        ICCOCG = 1
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGR(II) , IMLIGR(II) ,
     &   IPHYDP , IWARNI(II) , NFECRA ,
     &   EPSRGR(II) , CLIMGR(II) , EXTRAG(II) ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W11    , W11    , W11    ,
     &   W11    , COEFA(1,ICLRTP(II,ICOEF))  ,
     &            COEFB(1,ICLRTP(II,ICOEF))  ,
     &   W7              , W8              , W9     ,
C        d./dx1          , d./dx2          , d./dx3 ,
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )

C
C --- Calcul du terme source
C     ======================
C
C
C ---> Calcul de K et Epsilon en fonction du modele de turbulence
C
C
C ---- TURBULENCE
C
        IF (ITYTUR(IPHAS).EQ.2) THEN
C
          DO IEL = 1, NCEL
            W10(IEL) = RTPA(IEL,IK(IPHAS))
            W11(IEL) = RTPA(IEL,IEP(IPHAS))
          ENDDO
C
        ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
C
          DO IEL = 1, NCEL
            W10(IEL) = ( RTPA(IEL,IR11(IPHAS))
     &                  +RTPA(IEL,IR22(IPHAS))
     &                  +RTPA(IEL,IR33(IPHAS)) ) / 2.D0
            W11(IEL) = RTPA(IEL,IEP(IPHAS))
          ENDDO
C
        ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
          DO IEL = 1, NCEL
            W10(IEL) = RTPA(IEL,IK(IPHAS))
            W11(IEL) = RTPA(IEL,IEP(IPHAS))
          ENDDO
C
        ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
          DO IEL = 1, NCEL
            W10(IEL) = RTPA(IEL,IK(IPHAS))
            W11(IEL) = CMU*RTPA(IEL,IK(IPHAS))*RTPA(IEL,IOMG(IPHAS))
          ENDDO
C
        ENDIF
C
        DO IEL=1,NCEL
C
C  A confirmer :
C   Le terme de dissipation devrait etre implicite
C   Dans le terme de dissipation, il manque une constante Cf
C   Peut-elle etre consideree egale a 1 ?
C   Verifier le signe du terme de production
C-
C terme implicite
C
C
          W11(IEL) = W11(IEL)/(W10(IEL)*RVARFL(ISCAL))
     &         *VOLUME(IEL)*PROPCE(IEL,IPCROM)
          ROVSDT(IEL) = ROVSDT(IEL) + MAX(W11(IEL),ZERO)
C
C terme de gradient
C
          TSGRAD =  (2.0D0
     &         * PROPCE(IEL,IPCVST)/(SIGMAS(ISCAL))
     &         *(W1(IEL)*W7(IEL)+W2(IEL)*W8(IEL)+W3(IEL)*W9(IEL)))
     &         *VOLUME(IEL)

C
C terme de dissipation
C
          TSDISS = -W11(IEL) * RTPA(IEL,IVAR)
C
C terme de chimique
C
          TSCHIM = ZERO
          DO IDIRAC = 1, NDIRAC
            TSCHIM =   TSCHIM
     &           + (PROPCE(IEL,IPTSCL(IDIRAC))
     &           *(PROPCE(IEL,IPFMEL(IDIRAC))-RTPA(IEL,ISCA(IFM)))
     &           *VOLUME(IEL))*PROPCE(IEL,IPRHOL(IDIRAC))
          ENDDO
C
C --> Somme des termes
C
          SMBRS(IEL) = SMBRS(IEL) + TSCHIM + TSGRAD + TSDISS
C
         ENDDO
C
       ENDIF
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
