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
                        SUBROUTINE LAGGRA
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP    , PROPCE , COEFA  , COEFB  ,
     &   GRADPR , GRADVF ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     1)   Calcul de (- (GRADIENT DE PRESSION) / ROM )
CFONC
CFONC     2)   Calcul du gradient de Vitesse
CFONC
c@fonce
C-----------------------------------------------------------------------
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
CARGU ! LNDNOD       ! E  !  -> ! DIM. CONNECTIVITE CELLULES->FACES    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
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
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! COEFA,COEFB  ! TR !  -> ! TABLEAUX DES COND LIM POUR PVAR      !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! GRADPR(NCEL,3! TR ! <-  ! GRADIENT DE PRESSION                 !
CARGU ! GRADVF(NCEL,9! TR ! <-  ! GRADIENT DE VITESSE FLUIDE           !
CARGU ! W1...W3(NCEL)! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
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
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION COEFA(NDIMFB,*) , COEFB(NDIMFB,*)
      DOUBLE PRECISION RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION GRADPR(NCELET,3) , GRADVF(NCELET,9)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
C
      INTEGER          INC , ICCOCG , IPRIPH , ICLIPR
      INTEGER          IUIPH , IVIPH , IWIPH
      INTEGER          IPCLIU , IPCLIV , IPCLIW
      INTEGER          IROMF  , IDIMTE , ITENSO
      INTEGER          IEL    , IPHYDP , IPHAS
      DOUBLE PRECISION UNSRHO
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 0. PARAMETRES
C=======================================================================
C
C-->Rappel de differents paramtrage possibles pour le calcul du gradient
C
C     INC     = 1     ! 0 RESOL SUR INCREMENT 1 SINON
C     ICCOCG  = 1     ! 1 POUR RECALCUL DE COCG 0 SINON
C     NSWRGR  = 100   ! mettre 1 pour un maillage regulie (par defaut : 100)
C     IMLIGR  = -1    ! LIMITATION DU GRADIENT : < 0 PAS DE LIMITATION
C     IWARNI  = 1     ! NIVEAU D'IMPRESSION
C     EPSRGR  = 1.D-8 ! PRECISION RELATIVE POUR LA REC GRA 97
C     CLIMGR  = 1.5D0 ! COEF GRADIENT*DISTANCE/ECART
C     EXTRAG  = 0     ! COEF GRADIENT*DISTANCE/ECART
C     IPHYDP  = 0     ! =1 NE CONCERNE QUE LE GRADIENT DE PRESSION
C
C
C-->Parametrage des calculs des gradients
C
      INC     = 1
      ICCOCG  = 1
      IPHYDP  = 0
C
C-->Numero de la phase continue qui est la phase porteuse
C
      IPHAS = ILPHAS
C
C=======================================================================
C 1. CALCUL DE :  - (GRADIENT DE PRESSION)/ROM
C=======================================================================
C
      IPRIPH = IPR(IPHAS)
      ICLIPR = ICLRTP(IPRIPH,ICOEF)
C
C Calcul du gradient de pression
C
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF (IRANGP.GE.0) CALL PARCOM (RTP(1,IPRIPH))
C                      ===========
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF (IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  (IDIMTE , ITENSO ,
     &   RTP(1,IPRIPH)  , RTP(1,IPRIPH) , RTP(1,IPRIPH) ,
     &   RTP(1,IPRIPH)  , RTP(1,IPRIPH) , RTP(1,IPRIPH) ,
     &   RTP(1,IPRIPH)  , RTP(1,IPRIPH) , RTP(1,IPRIPH)  )
      ENDIF
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPRIPH , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IPRIPH)  , IMLIGR(IPRIPH)  , IPHYDP ,
     &   IWARNI(IPRIPH)  , NFECRA ,
     &   EPSRGR(IPRIPH)  , CLIMGR(IPRIPH)  , EXTRAG(IPRIPH)  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTP(1,IPRIPH)  , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   GRADPR(1,1)     , GRADPR(1,2)     , GRADPR(1,3)     ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C Pointeur sur la masse volumique en fonction de l'ecoulement
C
      IF ( IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICFUEL).GE.0 ) THEN
        IROMF = IPPROC(IROM1)
      ELSE
        IROMF = IPPROC(IROM(IPHAS))
      ENDIF
C
C Calcul de -Grad P / Rom
C
      DO IEL = 1,NCEL
        UNSRHO = 1.D0 / PROPCE(IEL,IROMF)
        GRADPR(IEL,1) = -GRADPR(IEL,1) * UNSRHO
        GRADPR(IEL,2) = -GRADPR(IEL,2) * UNSRHO
        GRADPR(IEL,3) = -GRADPR(IEL,3) * UNSRHO
      ENDDO
C
C=======================================================================
C 2. CALCUL DU GRADIENT DE VITESSE :
C=======================================================================
C
      IF (MODCPL.GT.0 .AND. IPLAS.GE.MODCPL) THEN
C
        IUIPH = IU(IPHAS)
        IVIPH = IV(IPHAS)
        IWIPH = IW(IPHAS)
        IPCLIU = ICLRTP(IUIPH,ICOEF)
        IPCLIV = ICLRTP(IVIPH,ICOEF)
        IPCLIW = ICLRTP(IWIPH,ICOEF)
C
C     COMPOSANTE X
C     ============
C
C ---> TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTP(1,IUIPH))
C         ===========
        ENDIF
C
C ---> TRAITEMENT DE LA PERIODICITE
C
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &    (IDIMTE , ITENSO ,
     &     RTP(1,IUIPH),RTP(1,IUIPH),RTP(1,IUIPH),
     &     RTP(1,IUIPH),RTP(1,IUIPH),RTP(1,IUIPH),
     &     RTP(1,IUIPH),RTP(1,IUIPH),RTP(1,IUIPH) )
        ENDIF
C
C    Sans prise en compte de la pression hydrostatique
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IUIPH  , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IUIPH)   , IMLIGR(IUIPH)  , IPHYDP ,
     &   IWARNI(IUIPH)   , NFECRA ,
     &   EPSRGR(IUIPH)   , CLIMGR(IUIPH)  , EXTRAG(IUIPH)  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTP(1,IUIPH)   , COEFA(1,IPCLIU) , COEFB(1,IPCLIU) ,
     &   GRADVF(1,1)     , GRADVF(1,2)     , GRADVF(1,3)     ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C     COMPOSANTE Y
C     ============
C
C ---> TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTP(1,IVIPH))
C         ===========
        ENDIF
C
C ---> TRAITEMENT DE LA PERIODICITE
C
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &    (IDIMTE , ITENSO ,
     &     RTP(1,IVIPH),RTP(1,IVIPH),RTP(1,IVIPH),
     &     RTP(1,IVIPH),RTP(1,IVIPH),RTP(1,IVIPH),
     &     RTP(1,IVIPH),RTP(1,IVIPH),RTP(1,IVIPH) )
        ENDIF
C
C    Sans prise en compte de la pression hydrostatique
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVIPH  , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IVIPH)   , IMLIGR(IVIPH)  , IPHYDP ,
     &   IWARNI(IVIPH)   , NFECRA ,
     &   EPSRGR(IVIPH)   , CLIMGR(IVIPH)  , EXTRAG(IVIPH)  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTP(1,IVIPH)   , COEFA(1,IPCLIV) , COEFB(1,IPCLIV) ,
     &   GRADVF(1,4)     , GRADVF(1,5)     , GRADVF(1,6)     ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C     COMPOSANTE Z
C     ============
C
C ---> TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTP(1,IWIPH))
C         ===========
        ENDIF
C
C ---> TRAITEMENT DE LA PERIODICITE
C
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &    (IDIMTE , ITENSO ,
     &     RTP(1,IWIPH),RTP(1,IWIPH),RTP(1,IWIPH),
     &     RTP(1,IWIPH),RTP(1,IWIPH),RTP(1,IWIPH),
     &     RTP(1,IWIPH),RTP(1,IWIPH),RTP(1,IWIPH) )
        ENDIF
C
C
C    Sans prise en compte de la pression hydrostatique
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IWIPH  , IMRGRA , INC    , ICCOCG ,
     &   NSWRGR(IWIPH)   , IMLIGR(IWIPH)  , IPHYDP ,
     &   IWARNI(IWIPH)   , NFECRA ,
     &   EPSRGR(IWIPH)   , CLIMGR(IWIPH)  , EXTRAG(IWIPH)  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTP(1,IWIPH)   , COEFA(1,IPCLIW) , COEFB(1,IPCLIW) ,
     &   GRADVF(1,7)     , GRADVF(1,8)     , GRADVF(1,9)     ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
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
