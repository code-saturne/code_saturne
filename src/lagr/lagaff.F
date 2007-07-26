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
                        SUBROUTINE LAGAFF
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
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT , STATIS ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC        SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC        -----------------------------------
CFONC
CFONC   ECRITURE SUR FICHIERS DES INFORMATIONS SUR LE NOMBRE DE PARTICULES
CFONC        - nombre de particules dans le domaine
CFONC        - nombre de particules entrantes
CFONC        - nombre de particules sorties
CFONC
CFONC
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
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !  (NFABOR+1)  !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR  ) !    !     !  (OPTIONNEL)                         !
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
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA(NBPMAX, ! TR !  -> ! CARACTERISTIQUES DES PARTICULES      !
CARGU !       NVEP)  !    !     !  AUX PARTICULES (poids, ...)         !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! TEMPCT       ! TR !  -> ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU !  (NBPMAX,2)  !    !     !                                      !
CARGU ! STATIS(NCELET! TR ! <-> ! CUMUL DES STATISTIQUES VOLUMIQUES    !
CARGU !    NVLSTA)   !    !     !                                      !
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
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "radiat.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
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
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
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
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3) , TEMPCT(NBPMAX,2)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFINIA, IFINRA
C
      INTEGER          IPHAS , II1 , II2
      CHARACTER        NOMFIC*300
      DOUBLE PRECISION DNBPR
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
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
C 1. INITIALISATION
C=======================================================================
C
      IPHAS = ILPHAS
      IPASS = IPASS + 1
C
C=======================================================================
C 2. OUVERTURE DU FICHIER DE STOCKAGE
C=======================================================================
C
C     Seul le premier processeur ecrit les informations
      IF (IRANGP.LE.0) THEN
C
        IF (IPASS.EQ.1 ) THEN
C
          IF ( IROULE .GE. 1 .AND.
     &         (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) ) THEN
            WRITE(IMPLAL,1000)
          ELSEIF ( IROULE .GE. 1 .AND.
     &           (IPHYLA.NE.2 .OR. IENCRA.NE.1) ) THEN
            WRITE(IMPLAL,1001)
          ELSEIF ( IROULE .NE. 1 .AND.
     &           (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) ) THEN
            WRITE(IMPLAL,1002)
          ELSE
            WRITE(IMPLAL,1003)
          ENDIF
C
        ENDIF
C
C=======================================================================
C 2 - Ecriture des INFORMATIONS
C=======================================================================
C
        IF (NBPTOT.GT.0) THEN
          DNBPR = (NBPERT*100.D0)/DBLE(NBPTOT)
        ELSE
          DNBPR = 0
        ENDIF
C
        IF ( IROULE.GE.1 .AND.
     &       (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) ) THEN
C
          WRITE(IMPLAL,2000) IPLAS,(DTP*IPLAS),
     &         NBPART        , DNBPAR        ,
     &         NBPNEW        , DNBPNW        ,
     &         NBPOUT-NBPERR , DNBPOU-DNBPER ,
     &         NBPERR        , DNBPER        ,
     &         DNBPR         ,
     &         NPCSUP        , DNPCSU        ,
     &         NPCLON        , DNPCLO        ,
     &         NPKILL        , DNPKIL        ,
     &         NPENCR        , DNPENC
C
        ELSEIF ( IROULE.GE.1 .AND.
     &         (IPHYLA.NE.2 .OR. IENCRA.NE.1) ) THEN
C
          WRITE(IMPLAL,2001) IPLAS,(DTP*IPLAS),
     &         NBPART        , DNBPAR        ,
     &         NBPNEW        , DNBPNW        ,
     &         NBPOUT-NBPERR , DNBPOU-DNBPER ,
     &         NBPERR        , DNBPER        ,
     &         DNBPR         ,
     &         NPCSUP        , DNPCSU        ,
     &         NPCLON        , DNPCLO        ,
     &         NPKILL        , DNPKIL
C
        ELSEIF ( IROULE.LT.1 .AND.
     &         (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) ) THEN
C
          WRITE(IMPLAL,2002) IPLAS,(DTP*IPLAS),
     &         NBPART        , DNBPAR        ,
     &         NBPNEW        , DNBPNW        ,
     &         NBPOUT-NBPERR , DNBPOU-DNBPER ,
     &         NBPERR        , DNBPER        ,
     &         DNBPR         ,
     &         NPENCR        , DNPENC
C
        ELSE
C
          WRITE(IMPLAL,2003) IPLAS,(DTP*IPLAS),
     &         NBPART        , DNBPAR        ,
     &         NBPNEW        , DNBPNW        ,
     &         NBPOUT-NBPERR , DNBPOU-DNBPER ,
     &         NBPERR        , DNBPER        ,
     &         DNBPR
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/,
     &       '#    -------------------------------------     ',/,
     &       '#                                         ',/,
     &       '# colonne  1 : numero de pas de temps      ',/,
     &       '# colonne  2 : temps physique              ',/,
     &       '# colonne  3 : nbre inst. de part.    ',/,
     &       '# colonne  4 : nbre inst. de part. (avec poids)',/,
     &       '# colonne  5 : nbre inst. de part. injectees',/,
     &       '# colonne  6 : nbre inst. de part. injectees',
     &       ' (avec poids)',/,
     &       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,
     &       '# colonne  8 : nbre inst. de part. sorties ou deposees',
     &       ' (avec poids)',/,
     &       '# colonne  9 : nbre inst. de part. perdues (reperage) ',/,
     &       '# colonne 10 : nbre inst. de part. perdues',
     &        ' (reperage, avec poids)',/,
     &       '# colonne 11 : % de part. perdues',/,
     &       '# colonne 12 : nbre inst. de part. qui ont subi le',
     &       ' clonage',/,
     &       '# colonne 13 : nbre inst. de part. qui ont subi le',
     &       ' clonage (avec poids)',/,
     &       '# colonne 14 : nbre inst. de nouvel. part. par clonage',/,
     &       '# colonne 15 : nbre inst. de nouvel. part. par clonage',
     &       ' (avec poids)',/,
     &       '# colonne 16 : nbre inst. de nouvel. part. eliminees par',
     &       ' roulette russe ',/,
     &       '# colonne 17 : nbre inst. de nouvel. part. eliminees par',
     &       ' roulette russe (avec poids)',/,
     &       '# colonne 18 : nbre inst. de part encrassees',
     &       ' (Charbon)) '/,
     &       '# colonne 19 : nbre inst. de part encrassees',
     &       ' (Charbon, avec poids))',/,
     &       '# ')
C
 1001 FORMAT('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/,
     &       '#    -------------------------------------     ',/,
     &       '#                                         ',/,
     &       '# colonne  1 : numero de pas de temps      ',/,
     &       '# colonne  2 : temps physique              ',/,
     &       '# colonne  3 : nbre inst. de part.    ',/,
     &       '# colonne  4 : nbre inst. de part. (avec poids)',/,
     &       '# colonne  5 : nbre inst. de part. injectees',/,
     &       '# colonne  6 : nbre inst. de part. injectees',
     &       ' (avec poids)',/,
     &       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,
     &       '# colonne  8 : nbre inst. de part. sorties ou deposees',
     &       ' (avec poids)',/,
     &       '# colonne  9 : nbre inst. de part. perdues (reperage) ',/,
     &       '# colonne 10 : nbre inst. de part. perdues',
     &        ' (reperage, avec poids)',/,
     &       '# colonne 11 : % de part. perdues',/,
     &       '# colonne 12 : nbre inst. de part. qui ont subi le',
     &       ' clonage',/,
     &       '# colonne 13 : nbre inst. de part. qui ont subi le',
     &       ' clonage (avec poids)',/,
     &       '# colonne 14 : nbre inst. de nouvel. part. par clonage',/,
     &       '# colonne 15 : nbre inst. de nouvel. part. par clonage',
     &       ' (avec poids)',/,
     &       '# colonne 16 : nbre inst. de nouvel. part. eliminees par',
     &       ' roulette russe ',/,
     &       '# colonne 17 : nbre inst. de nouvel. part. eliminees par',
     &       ' roulette russe (avec poids)',/,
     &       '# ')
C
 1002 FORMAT('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/,
     &       '#    -------------------------------------     ',/,
     &       '#                                         ',/,
     &       '# colonne  1 : numero de pas de temps      ',/,
     &       '# colonne  2 : temps physique              ',/,
     &       '# colonne  3 : nbre inst. de part.    ',/,
     &       '# colonne  4 : nbre inst. de part. (avec poids)',/,
     &       '# colonne  5 : nbre inst. de part. injectees',/,
     &       '# colonne  6 : nbre inst. de part. injectees',
     &       ' (avec poids)',/,
     &       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,
     &       '# colonne  8 : nbre inst. de part. sorties ou deposees',
     &       ' (avec poids)',/,
     &       '# colonne  9 : nbre inst. de part. perdues (reperage)',/,
     &       '# colonne 10 : nbre inst. de part. perdues',
     &        ' (reperage, avec poids)',/,
     &       '# colonne 11 : % de part. perdues ',/,
     &       '# colonne 12 : nbre inst. de part. encrassees',
     &       ' (Charbon)) '/,
     &       '# colonne 13 : nbre inst. de part. encrassees',
     &       ' (Charbon, avec poids))',/,
     &       '# ')
C
 1003 FORMAT('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/,
     &       '#    -------------------------------------     ',/,
     &       '#                                         ',/,
     &       '# colonne  1 : numero de pas de temps      ',/,
     &       '# colonne  2 : temps physique              ',/,
     &       '# colonne  3 : nbre inst. de part.    ',/,
     &       '# colonne  4 : nbre inst. de part. (avec poids)',/,
     &       '# colonne  5 : nbre inst. de part. injectees',/,
     &       '# colonne  6 : nbre inst. de part. injectees',
     &       ' (avec poids)',/,
     &       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,
     &       '# colonne  8 : nbre inst. de part. sorties ou deposees',
     &       ' (avec poids)',/,
     &       '# colonne  9 : nbre inst. de part. perdues (reperage)',/,
     &       '# colonne 10 : nbre inst. de part. perdues',
     &        ' (reperage, avec poids)',/,
     &       '# colonne 11 : % de part. perdues ',/,
     &       '# ')
C
 2000 FORMAT(1X,I8,2X,E10.4,2X,4(I8,2X,E10.4),2X,E10.4,4(I8,2X,E10.4))
 2001 FORMAT(1X,I8,2X,E10.4,2X,4(I8,2X,E10.4),2X,E10.4,3(I8,2X,E10.4))
 2002 FORMAT(1X,I8,2X,E10.4,2X,4(I8,2X,E10.4),2X,E10.4,1(I8,2X,E10.4))
 2003 FORMAT(1X,I8,2X,E10.4,2X,4(I8,2X,E10.4),2X,E10.4)
C
C====
C FIN
C====
C
      RETURN
C
      END
c@z
