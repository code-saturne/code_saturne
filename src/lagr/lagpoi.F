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
                        SUBROUTINE LAGPOI
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC ,LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICOCEL , ITYCEL , IFRLAG , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , STATIS ,
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
CFONC      RESOLUTION DE L'EQUATION DE POISSON POUR LES VITESSE MOYENNES
CFONC                  DES PARTICULES
CFONC        ET CORRECTION DES VITESSES INSTANTANNEES
CFONC                  DES PARTICULES
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
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
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
CARGU ! ICOCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (LNDNOD)     !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (NCELET+1)   !    !     !    POINTEUR DU TABLEAU ICOCEL        !
CARGU ! IFRLAG       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE LAGRANGIEN           !
CARGU ! ITEPA        ! TE ! <-  ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
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
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR !  -> ! MOYENNES STATISTIQUES                !
CARGU !(NCELET,NVLSTA!    !     !                                      !
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
C
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "radiat.h"
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
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ICOCEL(LNDNOD) , ITYCEL(NCELET+1)
      INTEGER          IFRLAG(NFABOR) ,  ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFINIA, IFINRA
      INTEGER          NPT , IEL , IFAC
      INTEGER          IMLGRD
      INTEGER          IPHILA , IPHIL
      INTEGER          IW1   , IW2   , IW3   , IW4 , IW5
      INTEGER          IW6   , IW7   , IW8   , IW9
      INTEGER          IIFACL , IIRESP
      INTEGER          IDTR   , IFMALA , IFMALB
      INTEGER          IVISCF , IVISCB , IDAM   , IXAM   , IDAG   , IXAG
      INTEGER          IDRTP  , ISMBR  , IROVSD
      INTEGER          ICOEFAP , ICOEFBP
      INTEGER          IVAR0
      INTEGER          INC, ICCOCG
      INTEGER          NSWRGP , IMLIGP , IWARNP
      INTEGER          IDIMTE , ITENSO , IPHYDP
      DOUBLE PRECISION EPSRGP , CLIMGP , EXTRAP
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C   Multigrille algebrique  IMLGRD = 0 --> NON
C                           IMLGRD = 1 --> OUI
C   Modifier aussi sa valeur dans LAGEQP
C
      IMLGRD = 0
C
      IIFACL =       IDEBIA
      IIRESP =       IIFACL + NFAC*2*IMLGRD
      IFINIA =       IIRESP + NCELET*IMLGRD
      CALL IASIZE('LAGPOI',IFINIA)
C     ===========
C
      IDTR   = IDEBRA
      IVISCF = IDTR   + NCELET
      IVISCB = IVISCF + NFAC
      IDAM   = IVISCB + NFABOR
      IXAM   = IDAM   + NCELET
      IDAG   = IXAM   + NFAC*2
      IXAG   = IDAG   + NCELET*IMLGRD
      IDRTP  = IXAG   + NFAC*2*IMLGRD
      ISMBR  = IDRTP  + NCELET
      IROVSD = ISMBR  + NCELET
      IFMALA = IROVSD + NCELET
      IFMALB = IFMALA + NFAC
C
      IPHILA  = IFMALB + NFABOR
      IPHIL   = IPHILA  + NCELET
      IW1    = IPHIL   + NCELET
      IW2    = IW1    + NCELET
      IW3    = IW2    + NCELET
      IW4    = IW3    + NCELET
      IW5    = IW4    + NCELET
      IW6    = IW5    + NCELET
      IW7    = IW6    + NCELET
      IW8    = IW7    + NCELET
      IW9    = IW8    + NCELET
      IFINRA = IW9    + NCELET
      CALL RASIZE('LAGPOI',IFINRA)
C     ==========
C
      DO IEL=1,NCEL
        IF ( STATIS(IEL,ILPD) .GT. SEUIL ) THEN
          STATIS(IEL,ILVX) = STATIS(IEL,ILVX)
     &                      /STATIS(IEL,ILPD)
          STATIS(IEL,ILVY) = STATIS(IEL,ILVY)
     &                      /STATIS(IEL,ILPD)
          STATIS(IEL,ILVZ) = STATIS(IEL,ILVZ)
     &                      /STATIS(IEL,ILPD)
          STATIS(IEL,ILFV) = STATIS(IEL,ILFV)
     &                      /( DBLE(NPST) * VOLUME(IEL) )
        ELSE
          STATIS(IEL,ILVX) = 0.D0
          STATIS(IEL,ILVY) = 0.D0
          STATIS(IEL,ILVZ) = 0.D0
          STATIS(IEL,ILFV) = 0.D0
        ENDIF
      ENDDO
C
      CALL LAGEQP
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIFACL) , IA(IIRESP) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , PROPCE , PROPFA , PROPFB ,
     &   RA(IVISCF) , RA(IVISCB) ,
     &   RA(IDAM) , RA(IXAM) , RA(IDAG) , RA(IXAG) ,
     &   RA(IDRTP) , RA(ISMBR) , RA(IROVSD) ,
     &   RA(IFMALA) , RA(IFMALB) ,
     &   STATIS(1,ILVX) , STATIS(1,ILVY) , STATIS(1,ILVZ) ,
     &   STATIS(1,ILFV) ,
     &   RA(IPHILA) , RA(IPHIL) ,
     &   W1     , W2     , W3     , RA(IW1) , RA(IW2) ,
     &   RA(IW3) , RA(IW4) , RA(IW5) , RA(IW6) ,
     &   RA(IW7) , RA(IW8) , RA(IW9) ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C Calcul du gradient du Correcteur PHI
C ====================================
C
C
C       On alloue localement 2 tableaux de NFABOR pour le calcul
C         de COEFA et COEFB de W1,W2,W3
C
      ICOEFAP = IFINRA
      ICOEFBP = ICOEFAP + NFABOR
      IFINRA  = ICOEFBP + NFABOR
      CALL RASIZE ('LAGEQP',IFINRA)
C     ===========
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
        RA(ICOEFAP+IFAC-1) = RA(IPHIL+IEL-1)
        RA(ICOEFBP+IFAC-1) = ZERO
      ENDDO
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
C
C En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
      IF(IRANGP.GE.0) THEN
        CALL PARCOM(RA(IPHIL))
C       ===========
      ENDIF
C
C    Periodique
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    RA(IPHIL) , RA(IPHIL) , RA(IPHIL) ,
     &    RA(IPHIL) , RA(IPHIL) , RA(IPHIL) ,
     &    RA(IPHIL) , RA(IPHIL) , RA(IPHIL)  )
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
      CALL GRDCEL
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IPHIL) , RA(IPHIL) , RA(IPHIL)    ,
     &   RA(IPHIL) , RA(ICOEFAP) , RA(ICOEFBP) ,
     &   W1       , W2    , W3 ,
     &   RA(IW1)  , RA(IW2) , RA(IW3) ,
     &   RDEVEL , RTUSER , RA     )
C
C CORRECTION DES VITESSES MOYENNES ET RETOUR AU CUMUL
C
      DO IEL = 1,NCEL
        IF ( STATIS(IEL,ILPD) .GT. SEUIL ) THEN
          STATIS(IEL,ILVX) = STATIS(IEL,ILVX) - W1(IEL)
          STATIS(IEL,ILVY) = STATIS(IEL,ILVY) - W2(IEL)
          STATIS(IEL,ILVZ) = STATIS(IEL,ILVZ) - W3(IEL)
        ENDIF
      ENDDO
C
      DO IEL = 1,NCEL
        IF ( STATIS(IEL,ILPD) .GT. SEUIL ) THEN
          STATIS(IEL,ILVX) = STATIS(IEL,ILVX)*STATIS(IEL,ILPD)
          STATIS(IEL,ILVY) = STATIS(IEL,ILVY)*STATIS(IEL,ILPD)
          STATIS(IEL,ILVZ) = STATIS(IEL,ILVZ)*STATIS(IEL,ILPD)
          STATIS(IEL,ILFV) = STATIS(IEL,ILFV)
     &                      *( DBLE(NPST) * VOLUME(IEL) )
        ENDIF
      ENDDO
C
C CORRECTION DES VITESSES INSTANTANNES
C
      DO NPT = 1,NBPART
        IF ( ITEPA(NPT,JISOR).GT.0 ) THEN
          IEL = ITEPA(NPT,JISOR)
          ETTP(NPT,JUP) = ETTP(NPT,JUP) - W1(IEL)
          ETTP(NPT,JVP) = ETTP(NPT,JVP) - W2(IEL)
          ETTP(NPT,JWP) = ETTP(NPT,JWP) - W3(IEL)
        ENDIF
      ENDDO
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
C----
C FIN
C----
C
      END
c@z
