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
                        SUBROUTINE LAGUNE
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICOCEL , ITYCEL , IFRLAG , ITEPA  , INDEP  , IBORD  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , STATIS , STATIV , TSLAGR , PARBOR ,
     &   TAUP   , TLAG   , PIIL   , BX     , VAGAUS , TSUF   , TSUP   ,
     &   TSVAR  , TEMPCT , TSFEXT , CPGD1  , CPGD2  , CPGHT  ,
     &   GRADPR , GRADVF , CROULE , BRGAUS , TERBRU ,
     &   W1     , W2     , W3     , AUXL   , AUXL2  ,
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
CFONC    Sous-programme principal du module de modelisation Lagrangienne
CFONC    des ecoulements diphasiques a inclusions dispersees.
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
CARGU !   (LNDNOD)   !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE ! <-  ! CONNECTIVITE CELLULES -> FACES       !
CARGU !   (NCELET+1) !    !     !    POINTEUR DU TABLEAU ICOCEL        !
CARGU ! IFRLAG       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU !   (NFABOR)   !    !     !  POUR LE MODULE LAGRANGIEN           !
CARGU ! ITEPA        ! TE ! <-  ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! INDEP        ! TE ! <-  ! POUR CHAQUE PARTICULE :              !
CARGU !   (NBPMAX)   !    !     !   NUMERO DE LA CELLULE DE DEPART     !
CARGU ! IBORD        ! TE ! <-  ! CONTIENT LE NUMERO DE LA             !
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
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
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT ET PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR ! <-  ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR ! <-  ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR ! <-  ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR ! <-  ! MOYENNES STATISTIQUES                !
CARGU !(NCELET,NVLSTA!    !     !                                      !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! TSLAGR       ! TR ! <-  ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,NTERSL!    !     !   LAGRANGIEN SUR LA PHASE PORTEUSE   !
CARGU ! PARBOR       ! TR ! <-  ! INFOS SUR INTERACTION DES PARTICULES !
CARGU !(NFABOR,NVISBR!    !     !   AUX FACES DE BORD                  !
CARGU ! TAUP(NBPMAX) ! TR ! <-  ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR ! <-  ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! PIIL(NBPMAX,3! TR ! <-  ! TERME DANS L'INTEGRATION DES EDS Up  !
CARGU ! BX(NBPMAX,3,2! TR ! <-  ! CARACTERISTIQUES DE LA TURBULENCE    !
CARGU ! VAGAUS       ! TR ! <-  ! VARIABLES ALEATOIRES GAUSSIENNES     !
CARGU !(NBPMAX,NVGAUS!    !     !                                      !
CARGU ! TSUP(NBPMAX,3! TR ! <-  ! PREDICTION 1ER SOUS-PAS POUR         !
CARGU !              !    !     !   LA VITESSE DES PARTICULES          !
CARGU ! TSUF(NBPMAX,3! TR ! <-  ! PREDICTION 1ER SOUS-PAS POUR         !
CARGU !              !    !     !   LA VITESSE DU FLUIDE VU            !
CARGU ! TSVAR        ! TR ! <-  ! PREDICTION 1ER SOUS-PAS POUR LA      !
CARGU ! (NBPMAX,NVP1)!    !     !   VARIABLE COURANTE, UTILISE POUR LA !
CARGU ! TEMPCT       ! TR ! <-  ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU ! (NBPMAX,2)   !    !     !                                      !
CARGU ! TSFEXT(NBPMAX! TR ! <-  ! FORCES EXTERNES                      !
CARGU ! CPGD1,CPGD2, ! TR ! <-  ! TERMES DE DEVOLATILISATION 1 ET 2 ET !
CARGU !  CPGHT(NBPMAX!    !     !   DE COMBUSION HETEROGENE (CHARBON   !
CARGU !              !    !     !   AVEC COUPLAGE RETOUR THERMIQUE)    !
CARGU ! GRADPR(NCEL,3! TR ! <-  ! GRADIENT DE PRESSION                 !
CARGU ! GRADVF(NCEL,9! TR ! <-  ! GRADIENT DE VITESSE FLUIDE           !
CARGU ! CROULE       ! TR ! <-  ! FONCTION D'IMPORTANCE POUR ROULETTE  !
CARGU !   (NCELET)   !    !     !   RUSSE                              !
CARGU ! W1..W3(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU ! AUXL(NBPMAX,3! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! AUXL2        ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !    (NBPMAX,7)!    !     !                                      !
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
      INTEGER          IFRLAG(NFABOR) , ITEPA(NBPMAX,NIVEP)
      INTEGER          INDEP(NBPMAX) , IBORD(NBPMAX)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3) , PIIL(NBPMAX,3)
      DOUBLE PRECISION VAGAUS(NBPMAX,*) , BX(NBPMAX,3,2)
      DOUBLE PRECISION TSUF(NBPMAX,3) , TSUP(NBPMAX,3)
      DOUBLE PRECISION TSVAR(NBPMAX,NVP1)
      DOUBLE PRECISION TEMPCT(NBPMAX,2) , TSFEXT(NBPMAX)
      DOUBLE PRECISION CPGD1(NBPMAX) , CPGD2(NBPMAX) , CPGHT(NBPMAX)
      DOUBLE PRECISION BRGAUS(NBPMAX,*) , TERBRU(NBPMAX)
      DOUBLE PRECISION GRADPR(NCELET,3) , GRADVF(NCELET,9)
      DOUBLE PRECISION CROULE(NCELET)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION AUXL(NBPMAX,3) , AUXL2(NBPMAX,7)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFINIA, IFINRA
C
      INTEGER          IP     , NPT    , IOK
      INTEGER          NFIN   , NPARS  , IEL    , IVF
      INTEGER          NPAR1  , NPAR2
      INTEGER          IFORCE , IITSLG
      INTEGER          MODNTL
C
      DOUBLE PRECISION DNPARS
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C=======================================================================
C 0.  GESTION MEMOIRE ET COMPTEUR DE PASSAGE
C=======================================================================
C
      IPASS = IPASS + 1
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IPLAR = IPLAR + 1
      IPLAS = IPLAS + 1
C
      NBPNEW = 0
      NPCSUP = 0
      NPCLON = 0
      NPKILL = 0
      NPENCR = 0
      NBPOUT = 0
      NBPERR = 0
C
      DNBPNW = 0.D0
      DNPCSU = 0.D0
      DNPCLO = 0.D0
      DNPKIL = 0.D0
      DNPENC = 0.D0
      DNBPOU = 0.D0
      DNBPER = 0.D0
C
C-->Sur Champ fige Lagrangien : RTPA = RTP
C   Rem : cette boucle pourrait etre faite au 1er passage
C         mais la presence de usproj incite a la prudence...
C
      IF (IILAGR.EQ.3) THEN
        DO IVF = 1,NVAR
          DO IEL = 1,NCEL
            RTPA(IEL,IVF) = RTP(IEL,IVF)
          ENDDO
        ENDDO
      ENDIF
C
C-->au premier passage relatif :
C
      IF (IPLAR.EQ.1) THEN
C
C       Connectivite cellules -> faces
C
        CALL LAGDEB
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICOCEL , ITYCEL ,
     &   IDEVEL , ITUSER , IA     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C=======================================================================
C 2.  MISE A JOUR DES NOUVELLES PARTICULES ENTREES DANS LE DOMAINE
C=======================================================================
C
C Au premier pas de temps on initalise les particules avec RTP et
C non RTPA car RTPA = initialisation
C
      IF ( NTCABS.EQ.1 ) THEN
C
        CALL LAGENT
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYCEL     , ICOCEL      ,
     &   IA(IITYPF) , IA(IITRIF)  , IFRLAG , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(ISRFBN) , DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , VAGAUS , AUXL   , W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
      ELSE
C
        CALL LAGENT
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYCEL     , ICOCEL      ,
     &   IA(IITYPF) , IA(IITRIF)  , IFRLAG , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(ISRFBN) , DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , VAGAUS , AUXL   , W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
C
C=======================================================================
C 2.1 CALCUL DE LA FONCTION D'IMPORTANCE POUR LA ROULETTE RUSSE
C=======================================================================
C
      IF (IROULE.GE.1) THEN
C
        CALL USLARU
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF) , IA(IITRIF)  , ITEPA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(ISRFBN)  , DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , VAGAUS , CROULE , AUXL ,
     &   RA(IDIPAR) , RA(IYPPAR) ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
        IOK = 0
        DO IEL = 1,NCEL
          IF (CROULE(IEL).LE.0.D0) IOK = IOK + 1
        ENDDO
        IF (IOK.NE.0) THEN
          WRITE(NFECRA,9001)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
      ENDIF
C
C=======================================================================
C 3.  GESTION DU TEMPS QUI PASSE...
C=======================================================================
C
C-->Gestion du pas de temps Lagrangien
C
      DTP = DTREF
C
C-->Incrementation du TEMPS COURANT LAGRANGIEN
C
      TTCLAG = TTCLAG + DTP
C
C-->Test pour savoir si le domaine contient des particules
C
      IF (NBPART.EQ.0) GOTO 20
C
C-->On enregistre l'element de depart de la particule
C
      DO IP = 1,NBPART
        INDEP(IP) = ITEPA(IP,JISOR)
      ENDDO
C
C=======================================================================
C 4.  GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
C=======================================================================
C
C Au premier pas de temps on calcul les gradient avec RTP et
C non RTPA car RTPA = initialisation (gradients nuls)
C
      IF ( NTCABS.EQ.1 ) THEN
C
        CALL LAGGRA
C       ===========
     & ( IDEBIA , IDEBRA ,
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
C
      ELSE
C
        CALL LAGGRA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , COEFA  , COEFB  ,
     &   GRADPR , GRADVF ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C=======================================================================
C 4.  Initialisation des variables aleatoires gaussiennes
C=======================================================================
C
C---> CALCUL DES TIRAGES ALEATOIRES
C     remarque : NORMALEN est dans le fichier ZUFALL.F
C     ^^^^^^^^
C
      IF (IDISTU.EQ.1) THEN
        DO IVF = 1,NVGAUS
          CALL NORMALEN(NBPART, VAGAUS(1,IVF))
        ENDDO
      ELSE
        DO IVF = 1,NVGAUS
          DO IP = 1,NBPMAX
            VAGAUS(IP,IVF) = 0.D0
          ENDDO
        ENDDO
      ENDIF
C
C---> CALCUL DES TIRAGES ALEATOIRES POUR LE MVT BROWNIEN
C
      IF ( LAMVBR .EQ. 1 ) THEN
C
        DO IVF = 1,NBRGAU
          CALL NORMALEN(NBPART, BRGAUS(1,IVF))
        ENDDO
C
      ENDIF
C
C=======================================================================
C 5. PROGRESSION DES PARTICULES
C=======================================================================
C
 10   CONTINUE
C
      NOR = MOD(NOR,NORDRE)
      NOR = NOR + 1
C
C---> Recopie des resultats de l'etape precedente :
C
      IF (NOR.EQ.1) THEN
C
        DO IVF = 1,NVP
          DO IP = 1,NBPART
            ETTPA(IP,IVF) = ETTP(IP,IVF)
          ENDDO
        ENDDO
C
      ENDIF
C
C-----> CALCUL GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
C       EN N+1 (avec RTP)
C
      IF (NOR.EQ.2 .AND. IILAGR.NE.3) THEN
C
        CALL LAGGRA
C       ===========
     & ( IDEBIA , IDEBRA ,
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
C
      ENDIF
C
C-----> CALCUL DES CARACTERISTIQUES DES PARTICULES
C
      IF (NOR.EQ.1) THEN
C
C      sous pas de temps n (avec RTPA)
C
        CALL LAGCAR
C       ===========
     &   ( IDEBIA , IDEBRA ,
     &     NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &     NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &     NVAR   , NSCAL  , NPHAS  ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NTERSL , NVLSTA , NVISBR ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &     ITEPA  , IDEVEL , ITUSER , IA     ,
     &     XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &     VOLUME , DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &     ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   ,
     &     PIIL   , BX     , TEMPCT , STATIS ,
     &     GRADPR , GRADVF , W1     , W2     , AUXL(1,1)  ,
     &     RDEVEL , RTUSER , RA     )
C
      ELSE
C
C     sous pas de temps n+1 (avec RTP)
C
        CALL LAGCAR
C       ===========
     &   ( IDEBIA , IDEBRA ,
     &     NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &     NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &     NVAR   , NSCAL  , NPHAS  ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NTERSL , NVLSTA , NVISBR ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &     ITEPA  , IDEVEL , ITUSER , IA     ,
     &     XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &     VOLUME , DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &     ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   ,
     &     PIIL   , BX     , TEMPCT , STATIS ,
     &     GRADPR , GRADVF , W1     , W2     , AUXL(1,1) ,
     &     RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
C     POSITION, VITESSE FLUIDE, VITESSE PARTICULE
C
      CALL LAGESP
C     ===========
     &   ( IDEBIA , IDEBRA ,
     &     NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &     NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &     NVAR   , NSCAL  , NPHAS  ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NTERSL , NVLSTA , NVISBR ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &     ITEPA  , IBORD  , IDEVEL , ITUSER , IA     ,
     &     XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &     VOLUME ,
     &     DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &     ETTP   , ETTPA  , TEPA   ,
     &     STATIS , STATIV , TAUP   , TLAG   , PIIL   ,
     &     TSUF   , TSUP   , BX     , TSFEXT ,
     &     VAGAUS , GRADPR , GRADVF , BRGAUS , TERBRU ,
     &     AUXL(1,1) , AUXL2 , RDEVEL , RTUSER , RA     )
C
C---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
C     LIEES AUX PHYSIQUES PARTICULIERES PARTICULAIRES
C
      IF ( IPHYLA.EQ.1 .OR. IPHYLA.EQ.2 ) THEN
C
        IF ( NOR.EQ.1 ) THEN
          CALL LAGPHY
C         ===========
     &    ( IDEBIA , IDEBRA ,
     &      NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &      NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &      NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &      NTERSL , NVLSTA , NVISBR ,
     &      NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &      ITEPA  , IBORD  ,
     &      IDEVEL , ITUSER , IA     ,
     &      XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &      DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &      ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &      TSVAR  , AUXL   , CPGD1  , CPGD2  , CPGHT  ,
     &      W1     , W2     , W3     ,
     &      RDEVEL , RTUSER , RA     )
        ELSE
          CALL LAGPHY
C         ===========
     &    ( IDEBIA , IDEBRA ,
     &      NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &      NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &      NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &      NTERSL , NVLSTA , NVISBR ,
     &      NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &      ITEPA  , IBORD  ,
     &      IDEVEL , ITUSER , IA     ,
     &      XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &      DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &      ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &      TSVAR  , AUXL   , CPGD1  , CPGD2  , CPGHT  ,
     &      W1     , W2     , W3     ,
     &      RDEVEL , RTUSER , RA     )
        ENDIF
C
      ENDIF
C
C=======================================================================
C 6.  Couplage Retour - Calcul des termes sources
C=======================================================================
C
      IF (IILAGR.EQ.2 .AND. NOR.EQ.NORDRE) THEN
C
        IFINIA = IDEBIA
        IITSLG = IDEBRA
        IFINRA = IITSLG + NTERSL*NBPMAX
        CALL RASIZE('LAGUNE',IFINRA)
C       ===========
C
        CALL LAGCOU
C       ===========
     &   ( IFINIA , IFINRA ,
     &     NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &     NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NTERSL , NVLSTA , NVISBR ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &     ITEPA  , INDEP  , IBORD  ,
     &     IDEVEL , ITUSER , IA     ,
     &     VOLUME , RTP    , PROPCE ,
     &     ETTP   , ETTPA  , TEPA   , TAUP   ,
     &     TEMPCT , TSFEXT , TSLAGR ,
     &     CPGD1  , CPGD2  , CPGHT  ,
     &     RA(IITSLG)      , W1     , W2   ,
     &     AUXL(1,1) , AUXL(1,2)   , AUXL(1,3) ,
     &     RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C=======================================================================
C 7.  Reperage des particules - Traitement des conditions aux limites
C     pour la position des particules
C=======================================================================
C
      IF (NOR.EQ.1) THEN
C
        CALL LAGCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF)      , IA(IITRIF)      ,
     &   ICOCEL , ITYCEL , IFRLAG , ITEPA  , IBORD  , INDEP  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(ISRFBN)      ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , PARBOR , AUXL   ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (IERR.EQ.1) THEN
          CALL LAGERR
C         ===========
          GOTO 20
        ENDIF
C
      ENDIF
C
C=======================================================================
C 9.  ELIMINATION DES PARTICULES QUI SONT SORTIES DU DOMAINE
C=======================================================================
C
C     ATTENTION : NBPOUT contient les particules sorties de facon
C                 normal + les particules sorties en erreur de reperage.
C
      IF (NOR.EQ.NORDRE) THEN
C
        CALL LAGELI
C       ===========
     & ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NPARS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   DNPARS ,
     &   ETTP   , ETTPA  , TEPA   ,
     &   RDEVEL , RTUSER , RA )
C
        NBPOUT = NPARS
        DNBPOU = DNPARS
C
      ENDIF
C
C=======================================================================
C 10.  TEMPS DE SEJOUR
C=======================================================================
C
      IF (NOR.EQ.NORDRE) THEN
C
        DO NPT = 1,NBPART
          IF ( ITEPA(NPT,JISOR).GT.0 ) THEN
            TEPA(NPT,JRTSP) = TEPA(NPT,JRTSP) + DTP
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 11.  CALCUL STATISTIQUES
C=======================================================================
C
      IF (NOR.EQ.NORDRE .AND. ISTALA.EQ.1 .AND. IPLAS.GE.IDSTNT) THEN
C
        CALL LAGSTA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , TEPA   , STATIS , STATIV ,
     &   W1     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C=======================================================================
C 12.  Equation de Poisson
C=======================================================================
C
      IF (NOR.EQ.NORDRE .AND. ILAPOI.EQ.1) THEN
C
        CALL LAGPOI
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
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
C
      ENDIF
C
C=======================================================================
C 13.  Methode de reduction de variances : Clonage/Fusion des particules
C=======================================================================
C
      IF ( NOR.EQ.NORDRE .AND. IROULE.GE.1 ) THEN
C
        CALL LAGRUS
C       ===========
     &   ( IDEBIA , IDEBRA ,
     &     NCELET , NCEL   ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &     ITEPA  , INDEP  ,
     &     IDEVEL , ITUSER , IA     ,
     &     ETTP   , ETTPA  , TEPA   , CROULE ,
     &     RDEVEL , RTUSER , RA     )
C
        IF (NPCLON.GT.0) THEN
C
          NPAR1 = NBPART - NPCLON + 1
          NPAR2 = NBPART
C
          CALL LAGIPN
C         ===========
     &    ( IDEBIA , IDEBRA ,
     &      NCELET , NCEL   ,
     &      NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &      NPAR1  , NPAR2  ,
     &      NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &      ITEPA  ,
     &      IDEVEL , ITUSER , IA     ,
     &      RTP    ,
     &      ETTP   , TEPA   , VAGAUS ,
     &      W1     , W2     , W3     ,
     &      RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 14. UN AUTRE TOUR ?
C=======================================================================
C
      IF (NORDRE.EQ.2 .AND. NOR.EQ.1) GOTO 10
C
C=======================================================================
C 15. BRANCHEMENT UTILISATEUR POUR MODIF DES VARIABLES EVENTUELLES
C     EN FIN D'ITERATION LAGRANGIENNE
C=======================================================================
C
      CALL USLAST
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &   STATIS , STATIV ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 16. Visualisations
C=======================================================================
C
 20   CONTINUE
C
      NFIN = 0
C
C-->Stockage des trajectoires au format Ensight Gold
C
      IF (IENSI1.EQ.1) THEN
C
        IFORCE = 0
C
        CALL ENSLAG
C       ===========
     &   ( IDEBIA , IDEBRA ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NFIN   , IFORCE ,
     &     ITEPA  ,
     &     ETTP   , TEPA   , RA )
      ENDIF
C
      IF (IENSI2.EQ.1) THEN
        CALL ENSWAF
C       ===========
     &   ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NFIN   ,
     &     ITEPA  ,
     &     ETTP   , TEPA , AUXL   )
      ENDIF
C
C=======================================================================
C 17. NOMBRE DE PARITICULES PERDUES (SUITES COMPRISES)
C=======================================================================
C
      NBPERT = NBPERT + NBPERR
C
C=======================================================================
C 18. ECRITURE SUR FICHIERS DES INFORMATIONS SUR LE NOMBRE DE PARTICULES
C        - nombre de particules dans le domaine
C        - nombre de particules entrantes
C        - nombre de particules sorties
C        - ...
C
C=======================================================================
C
      IF (IPASS.EQ.1) THEN
         MODNTL = 0
      ELSEIF(NTLAL.GT.0) THEN
         MODNTL = MOD(NTCABS,NTLAL)
      ELSEIF(NTLAL.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
         MODNTL = 0
      ELSE
         MODNTL = 1
      ENDIF
C
      IF (MODNTL.EQ.0) THEN
         CALL LAGAFF
C        ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT , STATIS ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )

      ENDIF
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA TECHNIQUE DE CLONAGE/FUSION DES PARTICULES           ',/,
     &'@      EST ENCLENCHEE AVEC UNE FONCTION D''IMPORTANCE        ',/,
     &'@      COMPORTANT DES VALEURS NEGATIVES OU NULLES            ',/,
     &'@      (LAGUNE).                                             ',/,
     &'@                                                            ',/,
     &'@    LES ELEMENTS DU TABLEAU CROULE DOIVENT STRICTEMENT      ',/,
     &'@      POSITIFS.                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les valeurs de CROULE dans la subroutine USLARU. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
