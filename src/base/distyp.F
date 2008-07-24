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
                        SUBROUTINE DISTYP
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ISYMPA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DISTPA , PROPCE , UETBOR , DISTY  ,
     &   DAM    , XAM    , SMBDP  , ROVSDP ,
     &   RTPDP  , DRTP   ,
     &   QX     , QY     , QZ     , COEFQ  , ROM    , ROMB   ,
     &   FLUMAS , FLUMAB ,
     &   COEFAX , COEFBX , COEFAY , COEFBY , COEFAZ , COEFBZ ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , W7     ,
     &   W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL DE LA DISTANCE A LA PAROI PAR UTILISATION D'UNE EQUATION
CFONC DE TRANSPORT
CFONC
CFONC  On resout dT/dt + grad (T U) = Gamma Ti
CFONC
CFONC     avec U = grad(distance a la paroi)/||grad(distance a la paroi)||
CFONC          Gamma = div U
CFONC          Ti = T (puits ou injection a la valeur locale)
CFONC     et aux parois
CFONC          T = u*/nu
CFONC     et ailleurs
CFONC          grad(T).n = 0
CFONC
CFONC     on cherche un etat stationnaire
CFONC
c@fonce
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE COURANTE          !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ISYMPA       ! TE !  -> ! ZERO POUR ANNULER LE FLUX DE MASSE   !
CARGU ! (NFABOR     )!    !     ! (TRANSMIS MAIS NON UTILISE)          !
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
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! DISTPA(NCELET! TR !  -> ! TAB DES DISTANCES A LA PAROI         !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! UETBOR       ! TR !  -> ! VITESSE DE FROTTEMENT AU BORD        !
CARGU ! (NFABOR,NPHAS!    !     !  POUR VAN DRIEST EN LES              !
CARGU ! DISTY(NCELET)! TR ! <-  ! DISTANCE Y+                          !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! SMBDP        ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! ROVSDP       ! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! RTPDP(NCELET)! TR !  -  ! VAR DE TRAVAIL DU SCLAIRE DIFFUSE    !
CARGU ! DRTP(NCELET) ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! QX,Y,Z(NCELET! TR !  -  ! NORMALE AUX ISO-DISTANCE A LA PAROI  !
CARGU ! COEFQ        ! TR !  -  ! TABLEAU DE TRAVAIL CL INIMAS         !
CARGU ! (NFABOR,NDIM)!    !     !                                      !
CARGU ! ROM(NCELET   ! TR !  -  ! MASSE VOLUMIQUE AUX CELLULES         !
CARGU ! ROMB(NFABOR) ! TR !  -  ! MASSE VOLUMIQUE AUX BORDS            !
CARGU ! FLUMAS(NFAC) ! TR !  -  ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -  ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! COEFAX,Y,Z, B! TR !  -  ! TABLEAUX DES COND LIM POUR QX, QY, QZ!
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
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
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "coincl.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , IPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS),ISYMPA(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DISTPA(NCELET),PROPCE(NCELET,*)
      DOUBLE PRECISION UETBOR(NFABOR,NPHAS)
      DOUBLE PRECISION DISTY(NCELET)
      DOUBLE PRECISION DAM(NCELET)  , XAM(NFAC,2)
      DOUBLE PRECISION SMBDP(NCELET), ROVSDP(NCELET)
      DOUBLE PRECISION RTPDP(NCELET), DRTP(NCELET)
      DOUBLE PRECISION QX(NCELET)   , QY(NCELET)  , QZ(NCELET)
      DOUBLE PRECISION COEFQ(NFABOR,3)
      DOUBLE PRECISION ROM(NCELET)   ,ROMB(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC)  ,FLUMAB(NFABOR)
      DOUBLE PRECISION COEFAX(NFABOR),COEFBX(NFABOR)
      DOUBLE PRECISION COEFAY(NFABOR),COEFBY(NFABOR)
      DOUBLE PRECISION COEFAZ(NFABOR),COEFBZ(NFABOR)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IDIMTE, ITENSO
      INTEGER          IDTVA0, IVAR  , ICONVP, IDIFFP
      INTEGER          NDIRCP, IRESLP
      INTEGER          IESCAP, IFLMB0, IMASPE, IPHYDP
      INTEGER          NCYMXP, NITMFP, IPP
      INTEGER          IFAC  , IEL   , IPCVIS, INIT  , IPCROM
      INTEGER          INC   , ICCOCG, ISYM  , NTCONT, INFPAR
C
      DOUBLE PRECISION XNORME, DTMINY, DTMAXY, RELAXP, THETAP, TIMEY
      DOUBLE PRECISION XUSNMX, XUSNMN, XNORM0
      DOUBLE PRECISION DISMAX, DISMIN, USNA
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPASS  = IPASS + 1
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
C
C=======================================================================
C 2. AU PREMIER PAS DE TEMPS
C=======================================================================
C
C     Au premier pas de temps, on a en general u* = 0 (ou faux)
C       on ne calcule pas y+
C
C     En effet ca prend du temps, d'autant plus que u* est petit, car il
C       alors calculer y+ jusqu'a une grande distance des parois
C
      IF(NTCABS.EQ.1) THEN
C
        DO IEL = 1, NCEL
          DISTY(IEL) = GRAND
        ENDDO
C
        IF(IWARNY.GE.1) THEN
          WRITE(NFECRA,7000)
        ENDIF
C
        RETURN
C
      ENDIF
C
C=======================================================================
C 3. CALCUL DE Q=Grad(DISTPA)/|Grad(DISTPA)|
C=======================================================================
C
C     La distance a la paroi vaut 0 en paroi
C       par definition et obeit a un flux nul ailleurs
C
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &     ITYPFB(IFAC,IPHAS).EQ.IPARUG ) THEN
          COEFAX(IFAC) = 0.0D0
          COEFBX(IFAC) = 0.0D0
        ELSE
          COEFAX(IFAC) = 0.0D0
          COEFBX(IFAC) = 1.0D0
        ENDIF
      ENDDO
C
C     Calcul du gradient
C
      IF(IRANGP.GE.0) CALL PARCOM (DISTPA)
C                     ===========
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &( IDIMTE , ITENSO ,
     &  DISTPA , DISTPA , DISTPA ,
     &  DISTPA , DISTPA , DISTPA ,
     &  DISTPA , DISTPA , DISTPA )
      ENDIF
C
      INC    = 1
      ICCOCG = 1
      IPHYDP = 0
      IVAR   = 0
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGY , IMLIGY , IPHYDP ,
     &   IWARNY , NFECRA , EPSRGY , CLIMGY , EXTRAY ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   DISTPA , COEFAX , COEFBX ,
     &   QX     , QY     , QZ     ,
C        ------   ------   ------
     &   W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Normalisation (attention, le gradient peut etre nul, parfois)
C
      DO IEL = 1, NCEL
        XNORME = MAX(SQRT(QX(IEL)**2+QY(IEL)**2+QZ(IEL)**2),EPZERO)
        QX(IEL) = QX(IEL)/XNORME
        QY(IEL) = QY(IEL)/XNORME
        QZ(IEL) = QZ(IEL)/XNORME
      ENDDO
C
C=======================================================================
C 4. CALCUL DU FLUX DE Q ET DE GAMMA = div(Q)
C=======================================================================
C
C
      DO IFAC = 1, NFABOR
        ROMB(IFAC) = 1.D0
      ENDDO
      DO IEL = 1, NCELET
        ROM(IEL)  = 1.D0
      ENDDO

C
C     Le gradient normal de la distance a la paroi vaut -1 en paroi
C       par definition et obeit a un flux nul ailleurs
C
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &     ITYPFB(IFAC,IPHAS).EQ.IPARUG) THEN
          XNORME = MAX(RA(ISRFBN+IFAC-1),EPZERO**2)
          COEFAX(IFAC) = -SURFBO(1,IFAC)/XNORME
          COEFBX(IFAC) = 0.D0
          COEFAY(IFAC) = -SURFBO(2,IFAC)/XNORME
          COEFBY(IFAC) = 0.D0
          COEFAZ(IFAC) = -SURFBO(3,IFAC)/XNORME
          COEFBZ(IFAC) = 0.D0
        ELSE
          COEFAX(IFAC) = 0.D0
          COEFBX(IFAC) = 1.D0
          COEFAY(IFAC) = 0.D0
          COEFBY(IFAC) = 1.D0
          COEFAZ(IFAC) = 0.D0
          COEFBZ(IFAC) = 1.D0
        ENDIF
      ENDDO
C
C
C     Parallelisme en preparation du calcul du flux
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (QX)
C       ===========
        CALL PARCOM (QY)
C       ===========
        CALL PARCOM (QZ)
C       ===========
      ENDIF
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 1
        ITENSO = 0
        CALL PERCOM
C       ===========
     &( IDIMTE , ITENSO ,
     &  QX , QX , QX ,
     &  QY , QY , QY ,
     &  QZ , QZ , QZ )
      ENDIF
C
C
C     Calcul du flux de masse
C
C     On ne le met pas a zero en paroi (justement)
      IFLMB0 = 0
C     On l'initialise a 0
      INIT   = 1
C     On prend en compte les Dirichlet
      INC    = 1
C     On recalcule les gradients complets
      ICCOCG = 1
C     Calcul de flux std (pas de divrij)
      IMASPE = 1
C     Il ne s'agit ni de U ni de R
      IVAR = 0
C
      CALL INIMAS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IVAR   , IVAR   , IVAR   , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGY , IMLIGY ,
     &   IWARNY , NFECRA ,
     &   EPSRGY , CLIMGY , EXTRAY ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISYMPA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   ROM    , ROMB   ,
     &   QX     , QY     , QZ     ,
     &   COEFAX , COEFAY , COEFAZ , COEFBX , COEFBY , COEFBZ ,
     &   FLUMAS , FLUMAB ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFQ  ,
     &   RDEVEL , RTUSER , RA     )
C
C
C
C     A partir d'ici, QX, QY et QZ sont des tableaux de travail.
C
C
C
C     Calcul de la divergence du flux de masse
C
      INIT = 1
C
      CALL DIVMAS
C     ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR , INIT   , NFECRA ,
     &   IFACEL , IFABOR , FLUMAS , FLUMAB , QX)
C
C=======================================================================
C 5.CONDITIONS LIMITES
C=======================================================================
C
C     Dirichlet en u*/nu aux parois, et flux nul ailleurs
C
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &     ITYPFB(IFAC,IPHAS).EQ.IPARUG) THEN
          IEL = IFABOR(IFAC)
          COEFAX(IFAC) = UETBOR(IFAC,IPHAS)
     &                  *PROPCE(IEL,IPCROM)/PROPCE(IEL,IPCVIS)
          COEFBX(IFAC) = 0.0D0
        ELSE
          COEFAX(IFAC) = 0.0D0
          COEFBX(IFAC) = 1.0D0
        ENDIF
      ENDDO
C
C=======================================================================
C 6.CALCUL DU PAS DE TEMPS
C=======================================================================
C
C     On vise un Courant infini (de l'ordre de 1000).
C
C     On calcule avec MATRDT DAM = Sigma a S/d
      ICONVP = 1
      IDIFFP = 0
C     La matrice est non symetrique
      ISYM   = 2
C
      CALL MATRDT
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONVP , IDIFFP , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COEFBX , FLUMAS , FLUMAB , FLUMAS , FLUMAB , W2     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Le Courant est COUMXY = DT W2 / VOLUME
C         d'ou DTMINY = MIN(COUMXY * VOLUME/W2)
C     Au cas ou une cellule serait a W2(IEL)=0,
C       on ne la prend pas en compte
C
C     On prend dans QZ un pas de temps variable,
C       si on ne trouve pas de pas de temps, on prend le minimum
C       (cellules sources)
      DTMINY = GRAND
      DTMAXY = -GRAND
      DO IEL = 1, NCEL
        QZ(IEL) = -GRAND
        IF(W2(IEL).GT.EPZERO) THEN
          QZ(IEL) = COUMXY*VOLUME(IEL)/W2(IEL)
          DTMINY  = MIN(QZ(IEL),DTMINY)
          DTMAXY  = MAX(QZ(IEL),DTMAXY)
        ENDIF
      ENDDO
      IF(IRANGP.GE.0) THEN
        CALL PARMIN (DTMINY)
        CALL PARMAX (DTMAXY)
      ENDIF
      DTMINY = MAX(DTMINY,EPZERO)
C
      DO IEL = 1, NCEL
        IF(QZ(IEL).LE.0.D0) THEN
          QZ(IEL) = DTMINY
        ENDIF
      ENDDO
C
      IF(IWARNY.GE.2) THEN
        WRITE(NFECRA,2000)DTMINY,DTMAXY
      ENDIF
C
C=======================================================================
C 7. DIAGONALE DE LA MATRICE
C=======================================================================
C
      DO IEL = 1, NCEL
        ROVSDP(IEL) = VOLUME(IEL)*ROM(IEL)/QZ(IEL)-QX(IEL)
      ENDDO
C
C=======================================================================
C 8. BOUCLE EN TEMPS
C=======================================================================
C
C     Initialisations
C    =================
C
C     Iterations
      NTCONT = 0
C
C     Temps
      TIMEY = 0.D0
C
C
C
C     Inconnue
C       Au cas ou on n'atteint pas tout a fait l'etat stationnaire,
C       il faut que le yplus ne soit pas nul dans la zone ou les
C       conditions aux limites n'ont pas ete convectees. On voudrait
C       plutot que yplus y soit maximum.
C       Si on utilise zero ou une valeur negative comme initialisation,
C       on risque de se retrouver avec des valeurs proches de
C       zero issues de la diffusion due au schema upwind au voisinage
C       du front convecte et donc avec des yplus proches de zero
C       n'importe ou.
C       On va donc utiliser la valeur max de u*/nu.
C
C       A partir du second pas de temps, on a egalement le yplus du pas
C         de temps precedent
C
C     On calcule le min et le max
      XUSNMX = -GRAND
      XUSNMN =  GRAND
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &     ITYPFB(IFAC,IPHAS).EQ.IPARUG) THEN
          XUSNMX = MAX(XUSNMX,COEFAX(IFAC))
          XUSNMN = MIN(XUSNMN,COEFAX(IFAC))
        ENDIF
      ENDDO
      IF(IRANGP.GE.0) THEN
        CALL PARMAX (XUSNMX)
        CALL PARMIN (XUSNMN)
      ENDIF
C
      IF(IPASS.EQ.1) THEN
        DO IEL = 1, NCELET
          RTPDP(IEL) = XUSNMX
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          USNA = DISTY(IEL)/MAX(DISTPA(IEL),EPZERO)
          USNA = MAX(USNA,XUSNMN)
          USNA = MIN(USNA,XUSNMX)
          RTPDP(IEL) = USNA
        ENDDO
      ENDIF
C
C     Norme de reference (moyenne des u*/nu)
C       (on divise par le nombre de faces)
      XNORM0 = 0.D0
      INFPAR = 0
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &     ITYPFB(IFAC,IPHAS).EQ.IPARUG) THEN
          INFPAR = INFPAR+1
          XNORM0 = XNORM0 + COEFAX(IFAC)**2
        ENDIF
      ENDDO
      IF(IRANGP.GE.0) THEN
        CALL PARCPT (INFPAR)
        CALL PARSOM (XNORM0)
      ENDIF
      XNORM0 = XNORM0/DBLE(INFPAR)
C
C     Pour ne pas diviser par zero
      IF(XNORM0.LE.EPZERO**2) GOTO 100
C
C
C     Debut des iterations
C    ======================
C
      DO NTCONT = 1, NTCMXY
C
C     Instant (arbitrairement +DTMINY)
        TIMEY = TIMEY + DTMINY
C
C -- Echange pour les cas paralleles
C       a la premiere iteration, c'est inutile (on a fait l'init sur NCELET)
C
        IF(NTCONT.GT.1.OR.IPASS.GT.1) THEN
C
          IF(IRANGP.GE.0) CALL PARCOM (RTPDP)
C                         ===========
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            CALL PERCOM
C           ===========
     &      ( IDIMTE , ITENSO ,
     &        RTPDP, RTPDP, RTPDP,
     &        RTPDP, RTPDP, RTPDP,
     &        RTPDP, RTPDP, RTPDP)
          ENDIF
C
        ENDIF
C
C
C     Sauvegarde de la solution pour test de convergence
C    ====================================================
C
        DO IEL = 1, NCEL
          QY(IEL) = RTPDP(IEL)
        ENDDO

C     Second membre
C    ===============
C       Obligatoirement a tous les pas de temps
C
        DO IEL = 1, NCEL
          SMBDP(IEL) = QX(IEL)*RTPDP(IEL)
        ENDDO
C
C     Resolution
C    ============
C
C     La variable n'est pas la vitesse ou une composante de Rij
        IVAR = 0
C     Le cas est convectif, non diffusif
        ICONVP = 1
        IDIFFP = 0
C     Il y a des Dirichlet (car il y a des parois)
        NDIRCP = 1
C     On resout par la methode automatique
        IRESLP = -1
C     Pas d'estimateurs, ni de multigrille (100 et 10 sont arbitraires)
        IESCAP = 0
        NCYMXP = 100
        NITMFP = 10
C     La case 1 est une poubelle
        IPP    = 1
        NOMVAR(IPP) = 'YplusPar'
C     Ordre 1 en temps (etat stationnaire cherche)
        THETAP = 1.D0
C     Pas de stationnaire ni de relaxation -> a modifier eventuellement
        IDTVA0 = 0
        RELAXP = 1.D0
C
        CALL CODITS
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVA0 , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAY ,
     &   IMRGRA , NSWRSY , NSWRGY , IMLIGY , IRCFLY ,
     &   ISCHCY , ISSTPY , IESCAP ,
     &   IMGRPY , NCYMXP , NITMFP , IPP    , IWARNY ,
     &   BLENCY , EPSILY , EPSRSY , EPSRGY , CLIMGY , EXTRAY ,
     &   RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPDP  , RTPDP  ,
     &   COEFAX , COEFBX , COEFAX , COEFBX , FLUMAS , FLUMAB ,
     &   FLUMAS , FLUMAB , FLUMAS , FLUMAB ,
     &   ROVSDP , SMBDP  , RTPDP  ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Clipping (indispensable si on initialise par u*/nu du pas de
C    ==========                                   temps precedent)
C
        DO IEL = 1, NCEL
          RTPDP(IEL) = MAX(RTPDP(IEL),XUSNMN)
          RTPDP(IEL) = MIN(RTPDP(IEL),XUSNMX)
        ENDDO
C
C     Test d'arret
C    ==============
C
C     on utilise QY dans lequel la solution precedente a ete sauvee
C     on souhaite que la variation de l'inconnue sur chaque cellule
C       soit inferieure a un pourcentage de la valeur moyenne de
C       u*/nu calculee sur les faces de bord
C     on limite le test aux cellules qui pourront avoir un interet pour
C       VanDriest, c'est a dire qu'on ignore celles a y+ > YPLMXY
C       comme on ne connait pas encore y+, on se base sur min(u*/nu) :
C       on ignore les cellules a y min(u*/nu) > YPLMXY
C
C        XNORME = 0.D0
C        DO IEL = 1, NCEL
C          XNORME = XNORME + (RTPDP(IEL)-QY(IEL))**2
C        ENDDO
C        IF(IRANGP.GE.0) THEN
C          CALL PARSOM (XNORME)
C        ENDIF
C        XNORME = XNORME/DBLE(NCELGB)
C
        XNORME = -GRAND
        DO IEL = 1, NCEL
          IF(DISTPA(IEL)*XUSNMN.LE.YPLMXY) THEN
            XNORME = MAX(XNORME,(RTPDP(IEL)-QY(IEL))**2)
          ENDIF
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARMAX (XNORME)
        ENDIF
C
        IF(IWARNY.GE.2) THEN
          WRITE(NFECRA,3000)NTCONT,XNORME,XNORM0,XNORME/XNORM0
        ENDIF
C
        IF(XNORME.LE.EPSCVY*XNORM0) GOTO 100
C
      ENDDO
C
      WRITE(NFECRA,8000)XNORME,XNORM0,XNORME/XNORM0,NTCMXY
C
 100  CONTINUE
C
C=======================================================================
C 9. CALCUL DE YPLUS ET IMPRESSIONS
C=======================================================================
C
C
      DO IEL = 1, NCEL
        DISTY(IEL) = RTPDP(IEL)*DISTPA(IEL)
      ENDDO
C
      DISMAX = -GRAND
      DISMIN =  GRAND
C
      DO IEL = 1, NCEL
        DISMIN = MIN(DISTY(IEL),DISMIN)
        DISMAX = MAX(DISTY(IEL),DISMAX)
      ENDDO
C
      IF (IRANGP.GE.0) THEN
        CALL PARMIN(DISMIN)
        CALL PARMAX(DISMAX)
      ENDIF
C
      IF(IWARNY.GE.1) THEN
        WRITE(NFECRA,1000)DISMIN, DISMAX, MIN(NTCONT,NTCMXY)
      ENDIF
C
C

 1000 FORMAT(
     &'                                                             ',/,
     &' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,
     &'    ------------------------------------                     ',/,
     &'                                                             ',/,
     &'  Distance+ min = ',E14.5    ,' Distance+ max = ',E14.5       ,/,
     &'                                                             ',/,
     &'     (Calcul de la distance realise en ',I10   ,' iterations)',/)
 2000 FORMAT(
     &'                                                             ',/,
     &' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,
     &'    ------------------------------------                     ',/,
     &'                                                             ',/,
     &' Yplus:  Dt min = ',E14.5    ,'        Dt max = ',E14.5       ,/)
 3000 FORMAT(
     &'                                                             ',/,
     &' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,
     &'    ------------------------------------                     ',/,
     &'                                                             ',/,
     &' Yplus:  iteration   residu abs.     reference   residu rel. ',/,
     &' Yplus: ',I10    ,E14.5        ,E14.5        ,E14.5           ,/)
 7000 FORMAT(
     &'                                                             ',/,
     &' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,
     &'    ------------------------------------                     ',/,
     &'                                                             ',/,
     &'  Elle n''est pas calculee au premier pas de temps           ',/)
 8000   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,
     &'@    *********                                               ',/,
     &'@     Le systeme iteratif pour le calcul de la distance a la ',/,
     &'@       paroi adimensionnelle est imparfaitement converge.   ',/,
     &'@                                                            ',/,
     &'@          Residu     Reference            Residu relatif    ',/,
     &'@  ',2E14.5,12X,E14.5                                         ,/,
     &'@                                                            ',/,
     &'@     Augmenter la valeur de NTCMXY dans usini1 peut resoudre',/,
     &'@       le probleme.                                         ',/,
     &'@     La valeur actuelle de cet entier est ',I10              ,/,
     &'@                                                            ',/,
     &'@     Le calcul se poursuit.                                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C
      RETURN
      END
c@z
