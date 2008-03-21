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
                        SUBROUTINE VISDYN
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMAGOR ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     , W10    , XMIJ   ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DE LA VISCOSITE "TURBULENTE" POUR
CFONC UN MODELE LES SMAGORINSKI DYNAMIQUE
CFONC
CFONC  SMAGO = LijMij/MijMij
CFONC
CFONC  PROPCE(1,IVISCT(IPHAS)) = ROM * SMAGO  * L**2 * SQRT ( 2 * Sij.Sij )
CFONC        Sij = (DUi/Dxj + DUj/Dxi)/2
CFONC
CFONC  On dispose des types de faces de bord au pas de temps
CFONC    precedent (sauf au premier pas de temps, ou les tableaux
CFONC    ITYPFB et ITRIFB n'ont pas ete renseignes)
CFONC
c@fonce
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
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
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
CARGU ! SMAGOR(NCELET! TR !  -> ! CONSTANTE DE SMAGORINSKY DANS LE CAS !
CARGU ! , NPHAS)     !    !     ! D'UN MODLELE DYNAMIQUE               !
CARGU ! W1..10(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! XMIJ(NCELET,6! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
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
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION SMAGOR(NCELET)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET),W4(NCELET)
      DOUBLE PRECISION W5(NCELET),W6(NCELET),W7(NCELET),W8(NCELET)
      DOUBLE PRECISION W9(NCELET),W10(NCELET),XMIJ(NCELET,6)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          II, IEL, ICCOCG, INC
      INTEGER          IUIPH, IVIPH, IWIPH
      INTEGER          IPCLIU, IPCLIV, IPCLIW
      INTEGER          IPCROM, IPCVST, IPHYDP
      INTEGER          ICLIPC
      INTEGER          IDIMTE, ITENSO
      DOUBLE PRECISION COEF, RADEUX, DEUX, DELTA, DELTAF
      DOUBLE PRECISION S11, S22, S33, S11F, S22F, S33F
      DOUBLE PRECISION DUDY, DUDZ, DVDX, DVDZ, DWDX, DWDY
      DOUBLE PRECISION DUDYF, DUDZF, DVDXF, DVDZF, DWDXF, DWDYF
      DOUBLE PRECISION XFIL, XA, XB, XFIL2, XSMGMX
      DOUBLE PRECISION AIJ, BIJ
      DOUBLE PRECISION XL11, XL22, XL33, XL12, XL13, XL23
      DOUBLE PRECISION XM11, XM22, XM33, XM12, XM13, XM23
      DOUBLE PRECISION SMAGMA, SMAGMN, SMAGMY
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
C --- Memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Numero des variables (dans RTP)
      IUIPH = IU(IPHAS)
      IVIPH = IV(IPHAS)
      IWIPH = IW(IPHAS)
C
C --- Rang des variables dans PROPCE (prop. physiques au centre)
      IPCVST = IPPROC(IVISCT(IPHAS))
      IPCROM = IPPROC(IROM  (IPHAS))
C
C --- Rang des c.l. des variables dans COEFA COEFB
C        (c.l. std, i.e. non flux)
      IPCLIU = ICLRTP(IUIPH,ICOEF)
      IPCLIV = ICLRTP(IVIPH,ICOEF)
      IPCLIW = ICLRTP(IWIPH,ICOEF)
C
C --- Pour le calcul de la viscosite de sous-maille
      XFIL   = XLESFL(IPHAS)
      XFIL2  = XLESFD(IPHAS)
      XA     = ALES(IPHAS)
      XB     = BLES(IPHAS)
      DEUX   = 2.D0
      RADEUX = SQRT(DEUX)
      XSMGMX = SMAGMX(IPHAS)
C
C=======================================================================
C 2.  CALCUL DES GRADIENTS DE VITESSE ET DE
C       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
C=======================================================================
C
C     Les RTPA ont ete echange pour les calculs en parallele,
C       au debut du pas de temps (donc pas utile de le refaire ici)
C
      ICCOCG = 1
      INC = 1
      IPHYDP = 0
C
C W1 = DUDX, W2 = DUDY, W3=DUDZ
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IUIPH  , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IUIPH) , IMLIGR(IUIPH) , IPHYDP , IWARNI(IUIPH) ,
     &   NFECRA , EPSRGR(IUIPH) , CLIMGR(IUIPH) , EXTRAG(IUIPH) ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W6     , W6     , W6     ,
     &   RTPA(1,IUIPH) , COEFA(1,IPCLIU) , COEFB(1,IPCLIU) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
C Filtrage de W1, LE RESULTAT EST DANS W6
C
      CALL CFILTR
C     ===========
     & ( W1     , W6     , W7     , W8     )
C
      DO IEL = 1, NCEL
        S11   = W1(IEL)
        S11F  = W6(IEL)
        XMIJ(IEL,1)        = S11
        PROPCE(IEL,IPCVST) = S11**2
        W9(IEL) = S11F**2
      ENDDO
C
C
C            W2 = DUDY, W3=DUDZ
C W4 = DVDX, W1 = DVDY, W5=DVDZ
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVIPH  , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IVIPH) , IMLIGR(IVIPH) , IPHYDP , IWARNI(IVIPH) ,
     &   NFECRA , EPSRGR(IVIPH) , CLIMGR(IVIPH) , EXTRAG(IVIPH) ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W6     , W6     , W6     ,
     &   RTPA(1,IVIPH) , COEFA(1,IPCLIV) , COEFB(1,IPCLIV) ,
     &   W4     , W1     , W5     ,
C        ------   ------   ------
     &   W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      CALL CFILTR
C     ===========
     & ( W1     , W6     , W7     , W8    )
C
      DO IEL = 1, NCEL
        S22  = W1(IEL)
        S22F = W6(IEL)
        XMIJ(IEL,2) = S22
        PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST) + S22**2
        W9(IEL) = W9(IEL) + S22F**2
      ENDDO
C
      CALL CFILTR
C     ===========
     & ( W2     , W6     , W8     , W1     )
C
      CALL CFILTR
C     ===========
     & ( W4     , W7     , W8     , W1     )
C
      DO IEL = 1, NCEL
        DUDY  = W2(IEL)
        DVDX  = W4(IEL)
        DUDYF = W6(IEL)
        DVDXF = W7(IEL)
        XMIJ(IEL,4) = 0.5D0*(DUDY+DVDX)
        PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST) + 0.5D0*(DUDY+DVDX)**2
        W9(IEL) = W9(IEL) + 0.5D0*(DUDYF+DVDXF)**2
      ENDDO
C
C                       W3=DUDZ
C                       W5=DVDZ
C W2 = DWDX, W4 = DWDY, W1=DWDZ
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IWIPH  , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IWIPH) , IMLIGR(IWIPH) , IPHYDP , IWARNI(IWIPH) ,
     &   NFECRA , EPSRGR(IWIPH) , CLIMGR(IWIPH) , EXTRAG(IWIPH) ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W6     , W6     , W6     ,
     &   RTPA(1,IWIPH) , COEFA(1,IPCLIW) , COEFB(1,IPCLIW) ,
     &   W2     , W4     , W1     ,
C        ------   ------   ------
     &   W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      CALL CFILTR
C     ===========
     & ( W1     , W6     , W7     , W8     )
C
      DO IEL = 1, NCEL
        S33  = W1(IEL)
        S33F = W6(IEL)
        XMIJ(IEL,3) = S33
        PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST) + S33**2
        W9(IEL) = W9(IEL) + S33F**2
      ENDDO
C
      CALL CFILTR
C     ===========
     & ( W2     , W1     , W7     , W8     )
C
      CALL CFILTR
C     ===========
     & ( W3     , W6     , W7     , W8     )
C
      DO IEL = 1, NCEL
        DUDZ  = W3(IEL)
        DWDX  = W2(IEL)
        DUDZF = W6(IEL)
        DWDXF = W1(IEL)
        XMIJ(IEL,5) = 0.5D0*(DUDZ+DWDX)
        PROPCE(IEL,IPCVST) =
     &    PROPCE(IEL,IPCVST) + 0.5D0*(DUDZ+DWDX)**2
        W9(IEL) = W9(IEL) + 0.5D0*(DUDZF+DWDXF)**2
      ENDDO
C
      CALL CFILTR
C     ===========
     & ( W4     , W1     , W7     , W8     )
C
      CALL CFILTR
C     ===========
     & ( W5     , W6     , W7     , W8     )
C
      DO IEL = 1, NCEL
        DVDZ  = W5(IEL)
        DWDY  = W4(IEL)
        DVDZF = W6(IEL)
        DWDYF = W1(IEL)
        XMIJ(IEL,6) = 0.5D0*(DVDZ+DWDY)
        PROPCE(IEL,IPCVST) =
     &    PROPCE(IEL,IPCVST) + 0.5D0*(DVDZ+DWDY)**2
        W9(IEL) = W9(IEL) + 0.5D0*(DVDZF+DWDYF)**2
      ENDDO
C
      DO IEL = 1, NCEL
        PROPCE(IEL,IPCVST) = RADEUX*SQRT(PROPCE(IEL,IPCVST))
        W9(IEL)            = RADEUX*SQRT(W9(IEL)           )
      ENDDO
C
C     Ici XMIJ contient Sij
C         PROPCE(IEL,IPCVST) contient ||S||
C            SQRT(2)*SQRT(S11^2+S22^2+S33^2+2(S12^2+S13^2+S23^2))
C         W9                 contient ||SF||
C            SQRT(2)*SQRT(S11F^2+S22F^2+S33F^2+2(S12F^2+S13F^2+S23F^2))
C
C=======================================================================
C 3.  CALCUL DE Mij
C=======================================================================
C
      DO IEL = 1, NCEL
        W7(IEL) = XFIL *(XA*VOLUME(IEL))**XB
      ENDDO
C
      DO II = 1, 6
C
        CALL CFILTR
C       ===========
     & ( XMIJ(1,II) , W1     , W2     , W3     )
C
        DO IEL = 1, NCEL
          DELTA = W7(IEL)
          W2(IEL) = -DEUX*DELTA**2*PROPCE(IEL,IPCVST)*XMIJ(IEL,II)
        ENDDO
C
        CALL CFILTR
C       ===========
     & ( W2     , W3     , W4     , W5     )
C
        DO IEL = 1, NCEL
          DELTA = W7(IEL)
          DELTAF = XFIL2*DELTA
          AIJ    = -DEUX*DELTAF**2*W9(IEL)*W1(IEL)
          BIJ    =  W3(IEL)
          XMIJ(IEL,II) = AIJ - BIJ
        ENDDO
C
      ENDDO
C
C     Ici Aij contient alpha_ij, Bij contient beta_ij tilde
C        et XMIJ contient M_ij
C
C=======================================================================
C 4.  CALCUL DE LA CONSTANTE DE SMAGORINSKY DYNAMIQUE
C=======================================================================
C
C FILTRAGE DE LA VITESSE ET DE SON CARRE
C
C
C U**2
      DO IEL = 1,NCEL
        W9(IEL) = RTP(IEL,IUIPH)*RTP(IEL,IUIPH)
      ENDDO
      CALL CFILTR
C     ===========
     & ( W9     , W1     , W7     , W8     )
C
C V**2
      DO IEL = 1,NCEL
        W9(IEL) = RTP(IEL,IVIPH)*RTP(IEL,IVIPH)
      ENDDO
      CALL CFILTR
C     ===========
     & ( W9     , W2     , W7     , W8     )
C
C W**2
      DO IEL = 1,NCEL
        W9(IEL) = RTP(IEL,IWIPH)*RTP(IEL,IWIPH)
      ENDDO
      CALL CFILTR
C     ===========
     & ( W9     , W3     , W7     , W8     )
C
C UV
      DO IEL = 1,NCEL
        W9(IEL) = RTP(IEL,IUIPH)*RTP(IEL,IVIPH)
      ENDDO
      CALL CFILTR
C     ===========
     & ( W9     , W4     , W7     , W8     )
C
C UW
      DO IEL = 1,NCEL
        W9(IEL) = RTP(IEL,IUIPH)*RTP(IEL,IWIPH)
      ENDDO
      CALL CFILTR
C     ===========
     & ( W9     , W5     , W7     , W8     )
C
C VW
      DO IEL = 1,NCEL
        W9(IEL) = RTP(IEL,IVIPH)*RTP(IEL,IWIPH)
      ENDDO
      CALL CFILTR
C     ===========
     & ( W9     , W6     , W7     , W8     )
C
C U
      CALL CFILTR
C     ===========
     & ( RTP(1,IUIPH)    , W7     , W8     , W9     )
C
C V
      CALL CFILTR
C     ===========
     & ( RTP(1,IVIPH)    , W8     , W9     , SMAGOR )
C
C W
      CALL CFILTR
C     ===========
     & ( RTP(1,IWIPH)    , W9     , SMAGOR , W10    )
C
      DO IEL = 1, NCEL
C
C --- Calcul de Lij
        XL11 = W1(IEL) - W7(IEL) * W7(IEL)
        XL22 = W2(IEL) - W8(IEL) * W8(IEL)
        XL33 = W3(IEL) - W9(IEL) * W9(IEL)
        XL12 = W4(IEL) - W7(IEL) * W8(IEL)
        XL13 = W5(IEL) - W7(IEL) * W9(IEL)
        XL23 = W6(IEL) - W8(IEL) * W9(IEL)
C
        XM11 = XMIJ(IEL,1)
        XM22 = XMIJ(IEL,2)
        XM33 = XMIJ(IEL,3)
        XM12 = XMIJ(IEL,4)
        XM13 = XMIJ(IEL,5)
        XM23 = XMIJ(IEL,6)
C ---Calcul de Mij :: Lij
        W1(IEL) = XM11 * XL11 + 2.D0* XM12 * XL12 + 2.D0* XM13 * XL13 +
     &                                XM22 * XL22 + 2.D0* XM23 * XL23 +
     &                                                    XM33 * XL33
C ---Calcul de Mij :: Mij
        W2(IEL) = XM11 * XM11 + 2.D0* XM12 * XM12 + 2.D0* XM13 * XM13 +
     &                                XM22 * XM22 + 2.D0* XM23 * XM23 +
     &                                                    XM33 * XM33
C
      ENDDO
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM(W1)
        CALL PARCOM(W2)
      ENDIF
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    W1     , W1     , W1    ,
     &    W1     , W1     , W1    ,
     &    W1     , W1     , W1             )
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    W2     , W2     , W2    ,
     &    W2     , W2     , W2    ,
     &    W2     , W2     , W2             )
      ENDIF
C
C     Par defaut on fait une moyenne locale du numerateur et du
C     denominateur, puis seulement on fait le rapport.
C     L'utilisateur peut faire autrement dans USSMAG
C
      CALL CFILTR
C     ===========
     & ( W1     , W3     , W5     , W6     )
C
      CALL CFILTR
C     ===========
     & ( W2     , W4     , W5     , W6     )
C
      DO IEL = 1, NCEL
        IF(ABS(W4(IEL)).LE.EPZERO) THEN
          SMAGOR(IEL) = XSMGMX**2
        ELSE
          SMAGOR(IEL) = W3(IEL)/W4(IEL)
        ENDIF
      ENDDO
C
      CALL USSMAG
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMAGOR , W1     , W2     ,
     &   W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
      ICLIPC = 0
      DO IEL = 1, NCEL
        IF(SMAGOR(IEL).GE.XSMGMX**2) THEN
          SMAGOR(IEL) = XSMGMX**2
          ICLIPC = ICLIPC + 1
        ELSEIF(SMAGOR(IEL).LE.-XSMGMX**2) THEN
          SMAGOR(IEL) = -XSMGMX**2
          ICLIPC = ICLIPC + 1
        ENDIF
      ENDDO
C
C=======================================================================
C 3.  CALCUL DE LA VISCOSITE (DYNAMIQUE)
C=======================================================================
C
C On clippe en (mu + mu_t)>0 dans phyvar
C
      DO IEL = 1, NCEL
        COEF = SMAGOR(IEL)
        DELTA  = XFIL * (XA*VOLUME(IEL))**XB
        PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCROM)
     &       * COEF * DELTA**2 * PROPCE(IEL,IPCVST)
      ENDDO
C
C     Quelques impressions
      IF(IWARNI(IUIPH).GE.1) THEN
C
        SMAGMA = -1.0D12
        SMAGMN =  1.0D12
        SMAGMY =  0.D0
        DO IEL = 1, NCEL
          SMAGMA = MAX(SMAGMA,SMAGOR(IEL))
          SMAGMN = MIN(SMAGMN,SMAGOR(IEL))
          SMAGMY = SMAGMY + SMAGOR(IEL)*VOLUME(IEL)
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARMAX(SMAGMA)
C         ===========
          CALL PARMIN(SMAGMN)
C         ===========
          CALL PARSOM(SMAGMY)
C         ===========
          CALL PARCPT(ICLIPC)
C         ===========
        ENDIF
        SMAGMY = SMAGMY / VOLTOT
        WRITE(NFECRA,1000) ICLIPC
        WRITE(NFECRA,2001) IPHAS
        WRITE(NFECRA,2002) SMAGMA, SMAGMN, SMAGMY
        WRITE(NFECRA,2003)
C
      ENDIF
C
C----
C FORMAT
C----
C
 1000 FORMAT(
     &' Nb Clipping Constante Smagorinsky par valeurs maximales ',I10,/)
 2001 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &' --- Informations sur la constante de Smagorinsky^2          ',/,
     &' ----------------------------------                          ',/,
     &' Valeur moy  Valeur min  Valeur max                          ',/,
     &' ----------------------------------                          '  )
 2002 FORMAT(
     & E12.4    ,      E12.4,      E12.4                               )
 2003 FORMAT(
     &' ----------------------------------                          ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
