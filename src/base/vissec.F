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
                        SUBROUTINE VISSEC
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
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
     &   VISCF  , VISCB  , VISTOT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC AJOUT AU SECOND MEMBRE DES TERMES
CFONC
CFONC          GRAD( (K -2/3 MU) DIV(U) ) + DIV( MU (GRAD_TRANSPOSE(U)) )
CFONC
CFONC AVEC MU = MU_LAMINAIRE + MU_TURBULENT
CFONC   ET  K = VISCOSITE EN VOLUME (NULLE EN GENERAL)
CFONC
CFONC DIV(U)            EST CALCULE A PARTIR DE FLUMAS/RHO
CFONC GRAD_TRANSPOSE(U) EST UN GRADIENT CELLULE
CFONC
CFONC REMARQUES :
CFONC   - Theoriquement le terme en div(u) devrait plutot etre calcule
CFONC       par un gradient cellule, pour correspondre exactement au
CFONC       terme en dUj/dxi. Mais comme la partie en dUi/dxj est
CFONC       calculee completement autrement (gradient facette et implicitation)
CFONC       de toute facon on n'aura jamais Trace(tau_ij)=0 exactement.
CFONC   - Pour la meme raison, comme le terme en dUi/dxj est calcule sur les
CFONC       elements de bord et pas celui en dUj/dxi, il est difficile de
CFONC       traiter le terme en div(u) de maniere rigoureuse. Il est donc
CFONC       conserve sur les elements de bord.
CFONC   - En LES, le tenseur <(u-<u>)(u-<u>)> est modelise par mut <S>
CFONC       et non pas par mut <S> - 2/3 mut Tr(<S>) Id + 2/3 k Id
CFONC       de sorte que il n'apparait pas ici de mut div<u>
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
CARGU ! IPHAS        ! E  !  -> ! PHASE COURANTE                       !
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
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT PREC)           !
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
CARGU ! TRAV(NCELET,3! TR ! <-> ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! VISTOT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR MU           !
CARGU ! W1...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "vector.h"
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
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION TRAV(NCELET,3)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION VISTOT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ICCOCG, INC, IEL, IFAC, IVAR, ISOU, II, JJ, INIT
      INTEGER          IDIM
      INTEGER          IUIPH, IVIPH, IWIPH
      INTEGER          ICLVAR
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          IPCROM, IPBROM, IPCVIS, IPCVST, IFLMAS, IFLMAB
      INTEGER          IDIMTE, ITENSO, IPHYDP
      INTEGER          IPCVSV
C
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
      DOUBLE PRECISION ROMF, D2S3M, VECFAC
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IUIPH  = IU(IPHAS)
      IVIPH  = IV(IPHAS)
      IWIPH  = IW(IPHAS)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C
C
      IF(IPPMOD(ICOMPF).GE.0) THEN
        IF(IVISCV(IPHAS).GT.0) THEN
          IPCVSV = IPPROC(IVISCV(IPHAS))
        ELSE
          IPCVSV = 0
        ENDIF
      ELSE
        IPCVSV = -1
      ENDIF
C
C
      IFLMAS = IPPROF(IFLUMA(IUIPH))
C
      IPBROM = IPPROB(IROM  (IPHAS))
      IFLMAB = IPPROB(IFLUMA(IUIPH))
C
C
C     Si on extrapole les termes sources, on prend les prop a l'instant n
      IF(ISNO2T(IPHAS).GT.0) THEN
        IF(IROEXT(IPHAS).GT.0) THEN
          IPCROM = IPPROC(IROMA (IPHAS))
          IPBROM = IPPROB(IROMA (IPHAS))
        ENDIF
        IF(IVIEXT(IPHAS).GT.0) THEN
          IPCVIS = IPPROC(IVISLA(IPHAS))
          IPCVST = IPPROC(IVISTA(IPHAS))
        ENDIF
C     Il faudrait aussi faire quelque chose pour le flux de masse, non ?
      ENDIF
C
C --- Calcul de la viscosite totale
C
      IF (ITYTUR(IPHAS).EQ.3) THEN
        DO IEL = 1, NCEL
          VISTOT(IEL) = PROPCE(IEL,IPCVIS)
        ENDDO
      ELSE
        DO IEL = 1, NCEL
          VISTOT(IEL) = PROPCE(IEL,IPCVIS) + PROPCE(IEL,IPCVST)
        ENDDO
      ENDIF
C
C    Pour la periodicite de rotation, il faut avoir calcule
C      le gradient avec grdcel. La seule solution consiste donc a
C      echanger VISTOT puis a faire le produit, y compris sur les
C      cellules halo (calcul sur le halo, exceptionnellement).
C    Pour le parallelisme, on s'aligne sur la sequence ainsi definie.
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (VISTOT)
C       ===========
      ENDIF
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &( IDIMTE , ITENSO ,
     &  VISTOT , VISTOT , VISTOT ,
     &  VISTOT , VISTOT , VISTOT ,
     &  VISTOT , VISTOT , VISTOT )
      ENDIF
C
C
C=======================================================================
C 2.  CALCUL DES TERMES EN GRAD_TRANSPOSE
C=======================================================================
C
      DO ISOU = 1, 3
C
        IF (ISOU.EQ.1) IVAR = IUIPH
        IF (ISOU.EQ.2) IVAR = IVIPH
        IF (ISOU.EQ.3) IVAR = IWIPH
C
C Ceci pointe eventuellement sur ICLRTP(IVAR,ICOEF)
        ICLVAR = ICLRTP(IVAR,ICOEFF)
C
C --- Calcul du gradient de la vitesse
C
        ICCOCG = 1
        INC    = 1
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IWARNP = IWARNI(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
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
     &   W4     , W4     , W4     ,
     &   RTPA(1,IVAR)    , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
        DO IEL = 1, NCELET
          W6(IEL) = 1.D0
        ENDDO
        DO IFAC = 1, NFABOR
          W6(IFABOR(IFAC)) = 0.D0
        ENDDO
C
C --- Assemblage sur les faces internes
C
        DO IDIM = 1, 3
C
C On a echange le gradient dans grdcel et vistot plus haut
C
          IF(IDIM.EQ.1) THEN
            DO IEL = 1, NCELET
              W4(IEL) = VISTOT(IEL)*W1(IEL)
            ENDDO
          ELSEIF(IDIM.EQ.2) THEN
            DO IEL = 1, NCELET
              W4(IEL) = VISTOT(IEL)*W2(IEL)
            ENDDO
          ELSEIF(IDIM.EQ.3) THEN
            DO IEL = 1, NCELET
              W4(IEL) = VISTOT(IEL)*W3(IEL)
            ENDDO
          ENDIF
C
C
C
C
C On initialise TRAV(NCEL+1, NCELET)
C     (valeur bidon, mais pas NaN : les calculs sur le halo sont
C      par principe denue de sens, sauf exception)
          IF(NCELET.GT.NCEL) THEN
            DO IEL = NCEL+1, NCELET
              TRAV(IEL,IDIM) = 0.D0
            ENDDO
          ENDIF
C
C
C
          IF(IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
            DO IFAC = 1, NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
CMO             VECFAC = SURFAC(ISOU,IFAC)
CMO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
              VECFAC = SURFAC(ISOU,IFAC)*(W4(II)+W4(JJ))*0.5D0
              TRAV(II,IDIM) = TRAV(II,IDIM) + VECFAC*W6(II)
              TRAV(JJ,IDIM) = TRAV(JJ,IDIM) - VECFAC*W6(JJ)
            ENDDO
C
          ELSE
C
C VECTORISATION NON FORCEE
            DO IFAC = 1, NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
CMO             VECFAC = SURFAC(ISOU,IFAC)
CMO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
              VECFAC = SURFAC(ISOU,IFAC)*(W4(II)+W4(JJ))*0.5D0
              TRAV(II,IDIM) = TRAV(II,IDIM) + VECFAC*W6(II)
              TRAV(JJ,IDIM) = TRAV(JJ,IDIM) - VECFAC*W6(JJ)
            ENDDO
C
          ENDIF
C
C
C --- Assemblage sur les faces de bord
C
CMO          IF(IVECTB.EQ.1) THEN
CMOC
CMO!OCL NOVREC,VRL(16)
CMO            DO IFAC = 1, NFABOR
CMO             II = IFABOR(IFAC)
CMO             TRAV(II,IDIM) = TRAV(II,IDIM) + SURFBO(ISOU,IFAC)*W4(II)
CMO            ENDDO
CMOC
CMO          ELSE
CMOC
CMOC VECTORISATION NON FORCEE
CMO            DO IFAC = 1, NFABOR
CMO             II = IFABOR(IFAC)
CMO             TRAV(II,IDIM) = TRAV(II,IDIM) + SURFBO(ISOU,IFAC)*W4(II)
CMO            ENDDO
CMOC
CMO          ENDIF
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C 3.  CALCUL DES TERMES EN DIV
C=======================================================================
C  Pour periodicite et parallelisme, ROM est echange dans phyvar.
C     ou apres avoir ete calcule dans cfmsvl en compressible
C
C  Ici pour l'ordre 2 en temps, il faudrait tout prendre en n...
C
C
      DO IFAC = 1, NFAC
        II = IFACEL(1,IFAC)
        JJ = IFACEL(2,IFAC)
CMO        ROMF = POND(IFAC)*PROPCE(II,IPCROM)
CMO     &      + (1.D0-POND(IFAC))*PROPCE(JJ,IPCROM)
        ROMF = (PROPCE(II,IPCROM)+PROPCE(JJ,IPCROM))*0.5D0
        VISCF(IFAC) = PROPFA(IFAC,IFLMAS)/ROMF
      ENDDO
      DO IFAC = 1, NFABOR
        VISCB(IFAC) = PROPFB(IFAC,IFLMAB)/PROPFB(IFAC,IPBROM)
      ENDDO
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                                   IFACEL,IFABOR,VISCF ,VISCB ,W1)
C
      D2S3M = -2.D0/3.D0
C
C
      IF(IPCVSV.GT.0) THEN
        DO IEL = 1, NCEL
          W4(IEL) = ( PROPCE(IEL,IPCVSV) + D2S3M*VISTOT(IEL) )
     &            * W1(IEL)/VOLUME(IEL)
        ENDDO
      ELSEIF(IPCVSV.EQ.0) THEN
        DO IEL = 1, NCEL
          W4(IEL) = ( VISCV0(IPHAS)      + D2S3M*VISTOT(IEL) )
     &            * W1(IEL)/VOLUME(IEL)
        ENDDO
      ELSE
C
        IF( ITYTUR(IPHAS).EQ.4) THEN
          DO IEL = 1, NCEL
            W4(IEL) = D2S3M*PROPCE(IEL,IPCVIS)*W1(IEL)/VOLUME(IEL)
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            W4(IEL) = D2S3M*VISTOT(IEL)*W1(IEL)/VOLUME(IEL)
          ENDDO
        ENDIF
      ENDIF
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) CALL PARCOM (W4)
C                     ===========
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    W4     , W4     , W4    ,
     &    W4     , W4     , W4    ,
     &    W4     , W4     , W4    )
      ENDIF
C
C
      DO IFAC = 1, NFAC
        II = IFACEL(1,IFAC)
        JJ = IFACEL(2,IFAC)
CMO        VISCF (IFAC) = (POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
        VISCF (IFAC) = (W4(II)+W4(JJ))*0.5D0
      ENDDO
C
      DO ISOU = 1, 3
C
        IDIM = ISOU
C
C
C --- Assemblage sur les faces internes
C
        IF(IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VECFAC = SURFAC(ISOU,IFAC)*VISCF(IFAC)
            TRAV(II,IDIM) = TRAV(II,IDIM) + VECFAC
            TRAV(JJ,IDIM) = TRAV(JJ,IDIM) - VECFAC
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VECFAC = SURFAC(ISOU,IFAC)*VISCF(IFAC)
            TRAV(II,IDIM) = TRAV(II,IDIM) + VECFAC
            TRAV(JJ,IDIM) = TRAV(JJ,IDIM) - VECFAC
          ENDDO
C
        ENDIF
C
C
C --- Assemblage sur les faces de bord
C
        IF(IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFABOR
           II = IFABOR(IFAC)
           TRAV(II,IDIM) = TRAV(II,IDIM) + SURFBO(ISOU,IFAC)*W4(II)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFABOR
           II = IFABOR(IFAC)
           TRAV(II,IDIM) = TRAV(II,IDIM) + SURFBO(ISOU,IFAC)*W4(II)
          ENDDO
C
        ENDIF
C
C --- Calcul des efforts aux bords (partie 4/5)
C
        IF (INEEDF.EQ.1) THEN
          DO IFAC = 1, NFABOR
            II = IFABOR(IFAC)
            RA(IFORBR+(IFAC-1)*NDIM+ISOU-1) =
     &             RA(IFORBR+(IFAC-1)*NDIM+ISOU-1)
     &           + SURFBO(ISOU,IFAC)*W4(II)
          ENDDO
        ENDIF
C
      ENDDO
C
C
      RETURN
C
      END
c@z
