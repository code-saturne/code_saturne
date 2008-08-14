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
                        SUBROUTINE DIVERV
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     ,
     &   DIV    , UX     , VY     , WZ     ,
     &   COEFAX , COEFAY , COEFAZ ,
     &   COEFBX , COEFBY , COEFBZ ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC         CALCUL DE LA DIVERGENCE D'UN VECTEUR
CFONC
CFONC    (On ne s'embete pas, on appelle 3 fois le gradient)
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! DIV(NCELET)  ! TR ! <-  ! DIVERGENCE DU VECTEUR                !
CARGU ! UX,UY,UZ     ! TR ! <-  ! COMPOSANTE DU VECTEUR                !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! COEFAX,...   ! TR ! ->  ! CONDITIONS AUX LIMITES POUR LES      !
CARGU ! COEFBZ       !    !     ! FACES DE BORD                        !
CARGU ! (NFABOR)     !    !     !                                      !
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
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET)
      DOUBLE PRECISION DIV(NCELET)
      DOUBLE PRECISION UX(NCELET) , VY(NCELET) , WZ(NCELET)
      DOUBLE PRECISION COEFAX(NFABOR) , COEFAY(NFABOR) , COEFAZ(NFABOR)
      DOUBLE PRECISION COEFBX(NFABOR) , COEFBY(NFABOR) , COEFBZ(NFABOR)
      DOUBLE PRECISION W1(NCELET) , W2(NCELET) , W3(NCELET)
      DOUBLE PRECISION W4(NCELET) , W5(NCELET) , W6(NCELET)
      DOUBLE PRECISION W7(NCELET) , W8(NCELET) , W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR0
      INTEGER          IPHYDP,IDIMTE,ITENSO
      INTEGER          IEL
      INTEGER          INC, ICCOCG
      INTEGER          NSWRGP, IMLIGP, IWARNP
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
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
C En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
      IF(IRANGP.GE.0) THEN
        CALL PARCOM(UX)
C       ===========
        CALL PARCOM(VY)
C       ===========
        CALL PARCOM(WZ)
C       ===========
      ENDIF
C
C    Periodique
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 1
        ITENSO = 0
        CALL PERCOM
C       ===========
     &( IDIMTE , ITENSO ,
     &  UX     , UX     , UX     ,
     &  VY     , VY     , VY     ,
     &  WZ     , WZ     , WZ     )
      ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
        IVAR0 = 0
C
C    Sans prise en compte de la pression hydrostatique
C
      IPHYDP = 0
C
      INC = 1
      ICCOCG = 1
      NSWRGP = 100
      IMLIGP = -1
      IWARNP = 2
      EPSRGP = 1.D-8
      CLIMGP = 1.5D0
      EXTRAP = 0.D0
C
C=======================================================================
C 1. Calcul du gradient de UX DANS W1
C=======================================================================
C
      CALL GRDCEL
C     ===========
     &( IDEBIA , IDEBRA ,
     &  NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &  NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &  NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &  IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &  IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &  IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &  IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &  IDEVEL , ITUSER , IA     ,
     &  XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &  UX     , UX     , UX     ,
     &  UX     , COEFAX , COEFBX ,
     &  W1     , W4     , W5     ,
     &  W6     , W7     , W8     ,
     &  RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 2. Calcul du gradient de VY DANS W2
C=======================================================================
C
      CALL GRDCEL
C     ===========
     &( IDEBIA , IDEBRA ,
     &  NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &  NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &  NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &  IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &  IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &  IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &  IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &  IDEVEL , ITUSER , IA     ,
     &  XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &  VY     , VY     , VY     ,
     &  VY     , COEFAY , COEFBY ,
     &  W4     , W2     , W5     ,
     &  W6     , W7     , W8     ,
     &  RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 3. Calcul du gradient de VZ DANS W3
C=======================================================================
C
      CALL GRDCEL
C     ===========
     &( IDEBIA , IDEBRA ,
     &  NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &  NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &  NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &  IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &  IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &  IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &  IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &  IDEVEL , ITUSER , IA     ,
     &  XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &  WZ     , WZ     , WZ     ,
     &  WZ     , COEFAZ , COEFBZ ,
     &  W5     , W6     , W3     ,
     &  W7     , W8     , W9     ,
     &  RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 4. Calcul de la divergence du vecteur (UX,VY,WZ)
C=======================================================================
C
      DO IEL = 1,NCEL
        DIV(IEL) = W1(IEL) + W2(IEL) + W3(IEL)
      ENDDO
C
C----
C FIN
C----
C
      END
c@z
