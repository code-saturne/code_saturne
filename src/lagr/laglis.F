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
                        SUBROUTINE LAGLIS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
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
     &   ETTP   , TEPA   , STATIS , STATIV , TSLAGR , PARBOR ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC      ECITURE DES INFOS DANS LE LISTING
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
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
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
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR !  -> ! CUMUL POUR LES MOYENNES DES          !
CARGU !(NCELET,NVLSTA!    !     !   STATISTIQUES VOLUMIQUES            !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,NTERSL!    !     !   LAGRANGIEN SUR LA PHASE PORTEUSE   !
CARGU ! PARBOR       ! TR !  -> ! INFOS SUR INTERACTION DES PARTICULES !
CARGU !(NFABOR,NVISBR!    !     !   AUX FACES DE BORD                  !
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
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
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
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE) , IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFINIA, IFINRA
      INTEGER          IFAC , IEL , IVF , ITABVR , NBRCEL
      INTEGER          IVFF , IFLU , ICLA , II , NB
      DOUBLE PRECISION AA , BB , GMAX , GMIN , GMOY
      CHARACTER        CHCOND*16
C
C***********************************************************************
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IF (NBPART.NE.0) THEN
        AA = 100.D0 / DBLE(NBPART)
      ELSE
        AA = 0.D0
      ENDIF
C
C=======================================================================
C 2. AFFICHAGE LISTING
C=======================================================================
C
      WRITE (NFECRA,1000)
C
C NOMBRE DE PARTICULES
C
      WRITE(NFECRA,1001)
      WRITE(NFECRA,1010) IPLAS , IPLAR
      WRITE(NFECRA,1001)
      WRITE(NFECRA,1020)
      WRITE(NFECRA,1002)
      WRITE(NFECRA,1031) NBPNEW, DNBPNW
      IF (IROULE.GE.1) THEN
        WRITE(NFECRA,1037) NPCSUP, DNPCSU
        WRITE(NFECRA,1032) NPCLON, DNPCLO
        WRITE(NFECRA,1034) NPKILL, DNPKIL
      ENDIF
      IF (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) THEN
        WRITE(NFECRA,1038) NPENCR, DNPENC
      ENDIF
      WRITE(NFECRA,1033) NBPOUT-NBPERR, (DNBPOU-DNBPER)
      WRITE(NFECRA,1035) NBPERR, DNBPER
      WRITE(NFECRA,1036) NBPART, DNBPAR
      WRITE(NFECRA,1001)
      IF (NBPTOT.GT.0) THEN
        WRITE(NFECRA,1050) (NBPERT*100.D0)/DBLE(NBPTOT)
        WRITE(NFECRA,1001)
      ENDIF
C
C DEBIT SUR CHAQUE ZONE
C
      WRITE(NFECRA,7000)
      WRITE(NFECRA,1002)
      DO II = 1,NFRLAG
        NB = ILFLAG(II)
C
        IF ( IUSCLB(NB) .EQ. IENTRL) THEN
           CHCOND = 'ENTREE'
        ELSE IF ( IUSCLB(NB) .EQ. IREBOL) THEN
           CHCOND = 'REBOND'
        ELSE IF ( IUSCLB(NB) .EQ. ISORTL) THEN
           CHCOND = 'SORTIE'
        ELSE IF ( IUSCLB(NB) .EQ. IDEPO1 .OR.
     &       IUSCLB(NB) .EQ. IDEPO2 .OR.
     &       IUSCLB(NB) .EQ. IDEPO3     ) THEN
          CHCOND = 'DEPOSITION'
        ELSE IF ( IUSCLB(NB) .EQ. IENCRL) THEN
          CHCOND = 'ENCRASSEMENT'
        ELSE IF ( IUSCLB(NB) .EQ. IDEPFA) THEN
          CHCOND = 'FORCES_CHIMIQUES'
        ELSE
          CHCOND = 'UTILISATEUR'
        ENDIF
C
        WRITE(NFECRA,7001) NB,DEBLAG(NB)/DTP,CHCOND
      ENDDO
      WRITE(NFECRA,1001)
C
C STATISTIQUES VOLUMIQUES
C
      IF (ISTALA.EQ.1) THEN
        WRITE(NFECRA,2000)
        WRITE(NFECRA,1002)
        WRITE(NFECRA,2005) IDSTNT
C
        IF (IPLAS.GE.IDSTNT) THEN
C
          IF (ISTTIO.EQ.0) THEN
            WRITE(NFECRA,2010) NPSTT
          ENDIF
          IF (ISTTIO.EQ.1 .AND. IPLAS.LT.NSTIST) THEN
            WRITE(NFECRA,2020) NPSTT
            WRITE(NFECRA,2030) NSTIST
          ELSE IF (ISTTIO.EQ.1 .AND. IPLAS.GE.NSTIST) THEN
            WRITE(NFECRA,2020) NPSTT
            WRITE(NFECRA,2040) NPST
          ENDIF
          WRITE(NFECRA,1001)
C
          IF (NVLSTA.GT.0) THEN
            WRITE(NFECRA,3010)
            WRITE(NFECRA,1002)
C
            IFINIA = IDEBIA
            ITABVR = IDEBRA
            IFINRA = ITABVR + NCELET
            CALL RASIZE('LAGLIS',IFINRA)
C           ===========
C
C     MOYENNE
C
            DO IVF = 1, NVLSTA
C
              IVFF = IVF
              ICLA = 0
              IFLU = 0
C
              CALL USLAEN
C             ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVFF   , IVFF   , IVFF   , IFLU   , ILPD   , ICLA   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATIS , STATIV , RA(ITABVR),
     &   RDEVEL , RTUSER , RA     )
C
              NBRCEL = 0
C
              GMAX = -GRAND
              GMIN =  GRAND
              GMOY =  0.D0
C
              DO IEL = 1,NCEL
                BB = RA(ITABVR+IEL-1)
                IF (BB.GT.SEUIL) THEN
                  NBRCEL = NBRCEL + 1
                  GMAX = MAX (GMAX, BB)
                  GMIN = MIN (GMIN, BB)
                  GMOY = GMOY + BB
                ENDIF
              ENDDO
C
              IF (NBRCEL.GT.0) THEN
                GMOY = GMOY /DBLE(NBRCEL)
              ELSE
                GMAX =  0.D0
                GMIN =  0.D0
                GMOY =  0.D0
              ENDIF
              WRITE(NFECRA,3020) NOMLAG(IVF),  GMIN, GMAX, GMOY
C
            ENDDO
          ENDIF
        ENDIF
C
        WRITE(NFECRA,1001)
C
      ENDIF
C
C STATISTIQUES PARIETALES
C
      IF (IENSI3.EQ.1) THEN
C
        WRITE(NFECRA,5000)
        WRITE(NFECRA,1002)
        IF (ISTTIO.EQ.1) THEN
          IF (IPLAS.GE.NSTBOR) THEN
            WRITE(NFECRA,5020) NPSTF
          ELSE
            WRITE(NFECRA,5010) NSTBOR
          ENDIF
        ENDIF
        WRITE(NFECRA,5030) NPSTFT
        WRITE(NFECRA,1001)
C
        IF (NVISBR.GT.0) THEN
          WRITE(NFECRA,6000)
          WRITE(NFECRA,1002)
C
          IF (NVISBR.GT.1) THEN
C
            ITABVR = IDEBRA
            IFINRA = ITABVR + NFABOR
            CALL RASIZE('LAGLIS',IFINRA)
C           ===========

            DO IFAC = 1,NFABOR
              IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
                RA(ITABVR+IFAC-1) = 1.D0 / PARBOR(IFAC,INBR)
              ELSE
                RA(ITABVR+IFAC-1) = 0.D0
              ENDIF
            ENDDO
C
          ENDIF
C
          DO IVF = 1, NVISBR
C
            IVFF = IVF
            CALL LAGSTF
C           ===========
     &       ( NCELET , NFABOR , NVISBR ,
     &         IVFF   ,
     &         GMIN   , GMAX   , GMOY   ,
     &         PARBOR , RA(ITABVR)      )
C
            WRITE(NFECRA,6010) NOMBRD(IVF),  GMIN, GMAX, GMOY
C
          ENDDO
C
          WRITE(NFECRA,1001)
C
        ENDIF
C
      ENDIF
C
C INFO SUR LE COUPLAGE RETOUR
C
      IF (IILAGR.EQ.2) THEN
C
        IF (ISTTIO.EQ.0) THEN
          WRITE(NFECRA,4000)
          WRITE(NFECRA,1002)
C
        ELSE IF (ISTTIO.EQ.1) THEN
          WRITE(NFECRA,4010)
          WRITE(NFECRA,1002)
C
          IF (IPLAS.LT.NSTITS) THEN
            WRITE(NFECRA,4020) NSTITS
          ELSE IF (IPLAS.GE.NSTIST) THEN
            WRITE(NFECRA,4030) NPTS
          ENDIF
C
        ENDIF
C
        WRITE(NFECRA,4050) VMAX
        WRITE(NFECRA,4060) TMAMAX
        WRITE(NFECRA,4070) NTXERR
C
        WRITE(NFECRA,1001)
C
      ENDIF
C
C***********************************************************************
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN',/3X,
     &          '   -------------------------------------')
C
 1001 FORMAT('-----------------------------------------------------',
     &       '----------')
C
 1002 FORMAT('   ---------------------------------------------------',
     &       '-----')
C
 1010 FORMAT('Iters Lagrangiennes absolues/relatives : ',
     &         I10,' /',I10)
C
 1020 FORMAT('   Pour cette iteration, nombre de particules',/,
     &       '   (sans et avec leur poids statistique) :')
 1031 FORMAT('ln  nouvelles injectees                ',I8,3X,E14.5)
 1032 FORMAT('ln  nouvelles par clonage              ',I8,3X,E14.5)
 1033 FORMAT('ln  sorties ou deposees definitivement ',I8,3X,E14.5)
 1034 FORMAT('ln  eliminees par roulette russe       ',I8,3X,E14.5)
 1035 FORMAT('ln  perdues par erreur de reperage     ',I8,3X,E14.5)
 1036 FORMAT('ln  total restantes en fin de passage  ',I8,3X,E14.5)
 1037 FORMAT('ln  qui ont subit le clonage           ',I8,3X,E14.5)
 1038 FORMAT('ln  de charbon encrassees              ',I8,3X,E14.5)
C
 1050 FORMAT('% de particules perdues (suites comprises) : ',E10.4)
C
 2000 FORMAT('   Statistiques volumiques :')
C
 2005 FORMAT('Debut des statistiques a l''iteration absolue :   ',I10)
C
 2010 FORMAT('Nombre d''iterations dans les stats instationnaires : ',
     &        I10)
C
 2020 FORMAT('Nombre d''iterations total dans les statistiques :',I10)
C
 2030 FORMAT('RAZ des stats (debut du calcul stationnaire : ',I10,' )')
C
 2040 FORMAT('Nombre d''iterations dans les stats stationnaires :',I9)
C
 3010 FORMAT('   Stats volumiques  Valeur min    Valeur max    ',
     &       'Valeur moy')
C
 3020 FORMAT('lc  ',A13,2X,E12.5,2X,E12.5,2X,E12.5)
C
 4000 FORMAT('   Termes sources de couplage-retour instationnaires :')
C
 4010 FORMAT('   Termes sources de couplage-retour :')
C
 4020 FORMAT('RAZ des TS (debut du calcul stationnaire : ',I10,')')
C
 4030 FORMAT('Nombre d''iterations dans les TS stationnaires :',I10)
C
 4050 FORMAT('Taux volumiques max de particules : ',E14.5)
C
 4060 FORMAT('Taux massiques max de particules :  ',E14.5)
C
 4070 FORMAT('Nbr de cellules qui ont un taux volumique > 0.8 :',I10)
C
 5000 FORMAT('   Statistiques aux frontieres :')
C
 5010 FORMAT('RAZ des stats aux frontieres (debut stationnaire : ',
     &        I8,')')
C
 5020 FORMAT('Nbr d''iters des stats frontieres stationnaires : ',I10)
C
 5030 FORMAT('Nbr d''iters total dans les stats aux frontieres :',I10)
C
 6000 FORMAT('   Stats frontieres  Valeur min    Valeur max    ',
     &       'Valeur moy')
C
 6010 FORMAT('lp  ', A13, 2X, E12.5, 2X, E12.5, 2X, E12.5)
C
 7000 FORMAT(3X,'ZONE          DEBIT(kg/s)       TYPE CL    ')
C
 7001 FORMAT(2X, I3, 10X, E12.5, 9X, A16)
C
C====
C FIN
C====
C
      END
c@z
C
C---------------------------------------------------------------
CIterations Lagrangiennes absolues/relatives : 1234567 /1234567
C---------------------------------------------------------------
C   Pour cette iteration, nombre particules :
C   ---------------------------------------------------------
Cln   nouvelles injectees                12345678  0.3067e+04
Cln   qui ont ete clonees                12345678
Cln   nouvelles par clonage              12345678
Cln   sortantes du domaine               12345678
Cln   eliminees par roulette russe       12345678
Cln   perdues par erreur de reperage     12345678
Cln   total restantes en fin de passage  12345678
C---------------------------------------------------------------
CPoucentage de particules perdues pour ce calcul : 1234567890
C---------------------------------------------------------------
C   Statistiques volumiques :
C   ---------------------------------------------------------
CDebut des stat a l'iteration Lagrangienne absolue : 1234567
CNombre de passages dans les stat instationnaires : 1234567
CNombre de passages total dans les stat : 1234567
CRAZ des statistiques stationnaires (debut du cumul : 1234567)
CNombre de passages dans les stat stationnaires : 1234567
C---------------------------------------------------------------
C   Stat volumiques   Valeur min    Valeur max    Valeur moy
C   ---------------------------------------------------------
Clc   1234567890123  -0.47302e+03  -0.30677e+04  -0.30677e+04
C---------------------------------------------------------------
C   Termes sources de couplage-retour
C   ---------------------------------------------------------
CTaux volumiques max de particules :
CTaux massiques max de particules :
CNbr de cellules qui ont un taux volumique > 0.8 :
C---------------------------------------------------------------
C   Statistiques aux frontieres :
C   ---------------------------------------------------------
CRAZ des stat aux frontieres (debut stationnaire :   1234567890)
CNbr d''iter des stat aux frontieres stationnaires : 1234567890
CNbr d''iter total dans les stat aux frontieres :    1234567890
C---------------------------------------------------------------
C   Stat volumiques   Valeur min    Valeur max    Valeur moy
C   ---------------------------------------------------------
Clp   nombreImpacts       -0.47302e+03  -0.30677e+04
C---------------------------------------------------------------
C
