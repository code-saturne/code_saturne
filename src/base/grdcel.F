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
                        SUBROUTINE GRDCEL
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   PVAR   , COEFAP , COEFBP ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC APPEL DES DIFFERENTES ROUTINES DE CALCUL DE GRADIENT CELLULE
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IVAR         ! E  !  -> ! NUMERO DE LA VARIABLE                !
CARGU !              !    !     !   DESTINE A ETRE UTILISE POUR LA     !
CARGU !              !    !     !   PERIODICITE UNIQUEMENT (PERING)    !
CARGU !              !    !     !   ON POURRA DONNER IVAR=0 SI LA      !
CARGU !              !    !     !   VARIABLE N'EST NI UNE COMPOSANTE DE!
CARGU !              !    !     !   LA VITESSE, NI UNE COMPOSANTE DU   !
CARGU !              !    !     !   TENSEUR DES CONTRAINTES RIJ        !
CARGU ! IMRGRA       ! E  !  -> ! METHODE DE RECONSTRUCTION DU GRADIENT!
CARGU !              !    !     !  0 RECONSTRUCTION 97                 !
CARGU !              !    !     !  1 MOINDRES CARRES                   !
CARGU !              !    !     !  2 MOINDRES CARRES SUPPORT ETENDU    !
CARGU !              !    !     !    COMPLET                           !
CARGU !              !    !     !  3 MOINDRES CARRES AVEC SELECTION DU !
CARGU !              !    !     !    SUPPORT ETENDU                    !
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! NSWRGP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGP       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! IPHYDP       ! E  !  -> ! INDICATEUR DE PRISE EN COMPTE DE LA  !
CARGU !              !    !     ! PRESSION HYDROSTATIQUE               !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! EPSRGP       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! CLIMGP       ! R  !  -> ! COEF GRADIENT*DISTANCE/ECART         !
CARGU ! EXTRAP       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
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
CARGU ! PVAR  (NCELET! TR !  -> ! VARIABLE (PRESSION)                  !
CARGU ! COEFAP,COEFBP! TR !  -> ! TABLEAUX DES COND LIM POUR PVAR      !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! FEXTX,Y,Z    ! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   (NCELET)   !    !     !  HYDROSTATIQUE                       !
CARGU ! DPDX,DPDY    ! TR ! <-  ! GRADIENT DE PVAR                     !
CARGU ! DPDZ (NCELET !    !     !                                      !
CARGU ! DPDXA (NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
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
      INCLUDE "pointe.h"
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
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP
      INTEGER          IMLIGP ,IWARNP  , IPHYDP , NFECRA
      DOUBLE PRECISION EPSRGP , CLIMGP , EXTRAP
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION FEXTX(NCELET),FEXTY(NCELET),FEXTZ(NCELET)
      DOUBLE PRECISION PVAR(NCELET), COEFAP(NFABOR), COEFBP(NFABOR)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION DPDXA(NCELET),DPDYA(NCELET),DPDZA(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IDIMTE , ITENSO
      INTEGER          IIU(NPHSMX),IIV(NPHSMX),IIW(NPHSMX)
      INTEGER          IITYTU(NPHSMX)
      INTEGER          IIR11(NPHSMX),IIR22(NPHSMX),IIR33(NPHSMX)
      INTEGER          IIR12(NPHSMX),IIR13(NPHSMX),IIR23(NPHSMX)
      INTEGER          IMLINI
C
      DOUBLE PRECISION CLIMIN
C
C***********************************************************************
C
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 0. PREPARATION POUR PERIODICITE DE ROTATION
C=======================================================================
C
C Par defaut, on traitera le gradient comme un vecteur ...
C   (i.e. on suppose que c'est le gradient d'une grandeurs scalaire)
C
C S'il n'y a pas de rotation, les echanges d'informations seront
C   faits par percom (implicite)
C
C S'il y a une ou des periodicites de rotation,
C   on determine si la variables est un vecteur (vitesse)
C   ou un tenseur (de Reynolds)
C   pour lui appliquer dans percom le traitement adequat.
C   On positionne IDIMTE et ITENSO
C   et on recupere le gradient qui convient.
C Notons que si on n'a pas, auparavant, calcule et stocke les gradients
C   du halo on ne peut pas les recuperer ici (...).
C   Aussi ce sous programme est-il appele dans phyvar (dans perinu perinr)
C   pour calculer les gradients au debut du pas de temps et les stocker
C   dans DUDXYZ et DRDXYZ
C
C Il est necessaire que ITENSO soit toujours initialise, meme hors
C   periodicite, donc on l'initialise au prealable a sa valeur par defaut.
C
      IDIMTE = 1
      ITENSO = 0
C
      IF(IPERIO.EQ.1) THEN
C
C       On recupere d'abord certains pointeurs necessaires a PERING
C
        CALL PERGRA
C         ===========
     &  ( NPHSMX , NPHAS  ,
     &    IIU    , IIV    , IIW    ,
     &    IITYTU ,
     &    IIR11  , IIR22  , IIR33  , IIR12  , IIR13  , IIR23  )
C
        CALL PERING
C       ===========
     &  ( NPHAS  , IVAR   ,
     &    IDIMTE , ITENSO , IPEROT , IGUPER , IGRPER ,
     &    IIU    , IIV    , IIW    , IITYTU ,
     &    IIR11  , IIR22  , IIR33  , IIR12  , IIR13  , IIR23  ,
     &    DPDX   , DPDY   , DPDZ   ,
     &    RA(IDUDXY) , RA(IDRDXY)  )
      ENDIF
C
C=======================================================================
C 1. CALCUL DU GRADIENT
C=======================================================================
C
C     CALCUL VOLUME FINIS PUIS ITERATIONS DE RECONSTRUCTION
      IF (IMRGRA.EQ.0) THEN
C
        CALL GRADRC
C       ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
     &   IMRGRA , INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , EXTRAP ,
     &   IFACEL , IFABOR , IA(IICELB) , IVAR ,
     &   VOLUME , SURFAC , SURFBO , RA(IPOND), XYZCEN , CDGFAC , CDGFBO,
     &   RA(IDIJPF) , RA(IDIIPB) , RA(IDOFIJ) , FEXTX , FEXTY , FEXTZ  ,
     &   COEFAP , COEFBP , PVAR   ,
     &   RA(ICOCGB) , RA(ICOCG)   ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   DPDXA  , DPDYA  , DPDZA  )
C
C     MOINDRES CARRES
      ELSEIF(IMRGRA.EQ.1.OR.IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
C
       CALL CGRDMC
C       ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
     &   INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP , IMRGRA ,
     &   IWARNP , NFECRA , EPSRGP , EXTRAP ,
     &   IFACEL , IFABOR , IA(IICELB) , IA(IISYMP) ,
     &   VOLUME , SURFAC , SURFBO , RA(ISRFBN) , RA(IPOND)   ,
     &   RA(IDIST)   , RA(IDISTB) ,
     &                 RA(IDIJPF) , RA(IDIIPB)  ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   XYZCEN , CDGFAC , CDGFBO , COEFAP , COEFBP , PVAR   ,
     &   RA(ICOCGB)  , RA(ICOCG)  ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   DPDXA  , DPDYA  , DPDZA  )
C
C     MOINDRES CARRES PUIS ITERATIONS DE RECONSTRUCTION
      ELSEIF(IMRGRA.EQ.4) THEN
C
       CALL CGRDMC
C       ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
     &   INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP , IMRGRA ,
     &   IWARNP , NFECRA , EPSRGP , EXTRAP ,
     &   IFACEL , IFABOR , IA(IICELB) , IA(IISYMP) ,
     &   VOLUME , SURFAC , SURFBO , RA(ISRFBN) , RA(IPOND)   ,
     &   RA(IDIST)   , RA(IDISTB) ,
     &                 RA(IDIJPF) , RA(IDIIPB)  ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   XYZCEN , CDGFAC , CDGFBO , COEFAP , COEFBP , PVAR   ,
     &   RA(ICOCGB)  , RA(ICOCG)  ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   DPDXA  , DPDYA  , DPDZA  )
C
C on force la limitation de la solution initiale avec les options par defaut
C pour toutes les variables, quel que soit le choix de l'utilisateur sur la
C limitation du gradient final
C
        IMLINI = 1
        CLIMIN = 1.5D0
C
        CALL CLMGRD
C       ===========
     & ( IMRGRA , IMLINI , IWARNP , ITENSO , CLIMIN ,
     &   PVAR   , DPDX   , DPDY   , DPDZ   )
C
        CALL GRADRC
C       ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
     &   IMRGRA , INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , EXTRAP ,
     &   IFACEL , IFABOR , IA(IICELB) , IVAR ,
     &   VOLUME , SURFAC , SURFBO , RA(IPOND), XYZCEN , CDGFAC , CDGFBO,
     &   RA(IDIJPF) , RA(IDIIPB) , RA(IDOFIJ) , FEXTX , FEXTY , FEXTZ  ,
     &   COEFAP , COEFBP , PVAR   ,
     &   RA(ICOCIB) , RA(ICOCI)   ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   DPDXA  , DPDYA  , DPDZA  )
C
      ENDIF
C
C
C=======================================================================
C 2. LIMITATION DU GRADIENT (EVENTUELLE)
C=======================================================================
C
      CALL CLMGRD
C     ===========
     & ( IMRGRA , IMLIGP , IWARNP , ITENSO , CLIMGP ,
     &   PVAR   , DPDX   , DPDY   , DPDZ   )
C
C
      RETURN
      END
c@z
