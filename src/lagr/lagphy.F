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
                        SUBROUTINE LAGPHY
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IBORD  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   ,
     &   TEMPCT , TSVAR  , AUXL   ,
     &   CPGD1  , CPGD2  , CPGHT  ,
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
CFONC      INTEGRATION DES EDS CONCERNANT LES PHYSIQUES PARTICULIERES
CFONC        LIEES AUX PARTICULES :
CFONC
CFONC          - Temperature du fluide vu par les particules,
CFONC          - Temperature des particules,
CFONC          - Diametre des particules
CFONC          - Masse des particules
CFONC          - Variables liees aux grains de charbon (Temp,MCH,MCK),
CFONC          - Variables Utilisateur supplementaires.
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
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IBORD        ! TE !  -> ! CONTIENT LE NUMERO DE LA             !
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
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR ! <-> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! TEMPCT       ! TR !  -> ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU !  (NBPMAX,2)  !    !     !                                      !
CARGU ! TSVAR        ! TR ! <-> ! PREDICTION 1ER SOUS-PAS POUR LA      !
CARGU ! (NBPMAX,NVP1)!    !     !   VARIABLE COURANTE, UTILISE POUR LA !
CARGU !              !    !     !   CORRECTION AU 2EME SOUS-PAS        !
CARGU ! AUXL(NBPMAX,3! TR !  -  ! TABLEAU DE TRAVAIL LAGRANGIEN        !
CARGU ! CPGD1,CPGD2, ! TR ! <-  ! TERMES DE DEVOLATILISATION 1 ET 2 ET !
CARGU !  CPGHT(NBPMAX!    !     !   DE COMBUSION HETEROGENE (CHARBON   !
CARGU !              !    !     !   AVEC COUPLAGE RETOUR THERMIQUE)    !
CARGU ! W1..W3(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP) , IBORD(NBPMAX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3) , TEMPCT(NBPMAX,2)
      DOUBLE PRECISION TSVAR(NBPMAX,NVP1) , AUXL(NBPMAX,3)
      DOUBLE PRECISION CPGD1(NBPMAX) , CPGD2(NBPMAX) , CPGHT(NBPMAX)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFINIA , IFINRA
      INTEGER          IWL1   , IWL2
C
C***********************************************************************
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
C
C=======================================================================
C 2. INTEGRATION DE LA TEMPERATURE FLUIDE VU PAR LES PARTICULES
C=======================================================================
C
      IF ( IPHYLA.EQ.2 .OR. (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) ) THEN
C
        CALL LAGITF
C       ===========
     &  ( IDEBIA , IDEBRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &    NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NTERSL , NVLSTA , NVISBR ,
     &    ITEPA  , IBORD  , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &    DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &    ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &    TSVAR  , AUXL(1,1) , AUXL(1,2)  , W1     ,
     &    RA     )
C
      ENDIF
C
C=======================================================================
C 3. INTEGRATION DE LA TEMPERATURE DES PARTICULES
C=======================================================================
C
      IF ( IPHYLA.EQ.1 .AND. ITPVAR.EQ.1 ) THEN
C
        CALL LAGITP
C       ===========
     &  ( IDEBIA , IDEBRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &    NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NTERSL , NVLSTA , NVISBR ,
     &    ITEPA  , IBORD  , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &    DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &    ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &    TSVAR  , AUXL(1,1) , AUXL(1,2)  ,
     &    RA     )
C
      ENDIF
C
C=======================================================================
C 4. INTEGRATION DU DIAMETRE DES PARTICULES
C=======================================================================
C
      IF ( IPHYLA.EQ.1 .AND. IDPVAR.EQ.1 ) THEN
C
        CALL LAGIDP
C       ===========
     &  ( IDEBIA , IDEBRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &    NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NTERSL , NVLSTA , NVISBR ,
     &    ITEPA  , IBORD  , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &    DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &    ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &    TSVAR  , AUXL(1,1) , AUXL(1,2)  ,
     &    RA     )
C
      ENDIF
C
C=======================================================================
C 5. INTEGRATION DE LA MASSE DES PARTICULES
C=======================================================================
C
      IF (IPHYLA.EQ.1 .AND. IMPVAR.EQ.1) THEN
C
        CALL LAGIMP
C       ===========
     &  ( IDEBIA , IDEBRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &    NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NTERSL , NVLSTA , NVISBR ,
     &    ITEPA  , IBORD  , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &    DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &    ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &    TSVAR  , AUXL(1,1) , AUXL(1,2)  ,
     &    RA     )
C
      ENDIF
C
C=======================================================================
C 6. INTEGRATION DES EQUATIONS DU CHARBON : HP, MCH, MCK
C=======================================================================
C
      IF (IPHYLA.EQ.2) THEN
C
        IFINIA = IDEBIA
        IWL1   = IDEBRA
        IWL2   = IWL1 + NBPMAX
        IFINRA = IWL2 + NBPMAX
        CALL RASIZE('LAGUNE',IFINRA)
C       ==========
C
        CALL LAGICH
C       ===========
     &  ( IFINIA , IFINRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &    NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NTERSL , NVLSTA , NVISBR ,
     &    ITEPA  , IBORD  , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &    DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &    ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT , TSVAR  ,
     &    CPGD1  , CPGD2  , CPGHT  ,
     &    AUXL(1,1) , AUXL(1,2) , AUXL(1,3) ,
     &    RA(IWL1)  , RA(IWL2)  , W1        , RA     )
C
      ENDIF
C
C=======================================================================
C 7. INTEGRATION DES VARIABLES UTILISATEURS SUPPLEMENTAIRES
C=======================================================================
C
      IF (NVLS.GE.1) THEN
C
        CALL USLAED
C       ===========
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
     &      ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   ,
     &      TEMPCT , TSVAR  ,
     &      AUXL(1,1) ,  AUXL(1,2) ,  AUXL(1,3) ,
     &      W1     , W2     , W3     ,
     &      RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C***********************************************************************
C
C----
C FIN
C----
C
      END
c@z
