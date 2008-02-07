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
                        SUBROUTINE LAGOUT
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
     &   ICOCEL , ITYCEL , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   , PARBOR , STATIS , STATIV , TSLAGR ,
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
CFONC  1. Ecriture du fichier suite 'lagava' :
CFONC      * variables sur les particules (ETTP)
CFONC      * informations sur les particules (ITEPA, TEPA)
CFONC
CFONC  2. Ecriture du fichier suite statistiques et termes sources
CFONC      'lasava' :
CFONC      * statistiques volumiques (STATIS)
CFONC      * statistiques aux frontieres (PARBOR)
CFONC      * termes sources de couplage retour (TSLAGR)
CFONC
CFONC  3. Finalisation des sorties graphiques
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
CARGU ! ICOCEL       ! TE !  -> ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (LNDNOD)     !    !     !    FACE DE BORD SI NUMERO NEGATIF    !
CARGU ! ITYCEL       ! TE !  -> ! CONNECTIVITE CELLULES -> FACES       !
CARGU ! (NCELET+1)   !    !     !    POINTEUR DU TABLEAU ICOCEL        !
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
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "radiat.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
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
      INTEGER          ITEPA(NBPMAX,NIVEP)
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
      DOUBLE PRECISION ETTP(NBPMAX,NVP) ,  TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
C
      CHARACTER        RUBRIQ*64 , CAR4*4
      CHARACTER        NOMNVL(NVPLMX)*60 , NOMTSL(NVPLMX)*60
      CHARACTER        NOMITE(NVPLMX)*64 , NOMRTE(NVPLMX)*64
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFINIA , IFINRA
      INTEGER          IERROR , IRTYP  , ITYSUP , NBVAL
      INTEGER          IVERS  , ILECEC , IPHAS
      INTEGER          NFIN   , IFORCE , ICHA   , II
      INTEGER          ITRAV1
      INTEGER          IPAS   , JJ
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
C 1. ECRITURE DU FICHIER SUITE : VARIABLES LIEES AUX PARTICULES
C=======================================================================
C
      IPHAS = ILPHAS
C
C ---> Ouverture (et on saute si erreur)
C     IFOAVL   =0 --> non formate       =1 --> formate
C     ILECEC = 2 : ecriture
C
      IF (IFOAVL.EQ.0) THEN
        WRITE(NFECRA,6010)
      ELSE IF (IFOAVL.EQ.1) THEN
        WRITE(NFECRA,6020)
      ELSE
        WRITE(NFECRA,9011) FICAVL, 'IFOAVL', IFOAVL
      ENDIF
C
      ILECEC = 2
      CALL OPNSUI(FICAVL, LEN(FICAVL), ILECEC, IFOAVL, IMPAVL, IERROR)
C     ===========
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9010) FICAVL
        GOTO 9998
      ENDIF
C
      WRITE(NFECRA,6011)
C
C
C Entete et Infos sur le calcul ou on saute si erreur
C     On inclut une rubrique destinee a distinguer ce fichier
C       d'un autre fichier suite
C     Pour le moment, IVERS n'est pas utilise
C
      ITYSUP = 0
      NBVAL  = 1
C
      IVERS  = 111
      RUBRIQ = 'version_fichier_suite_Lagrangien_variables'
      IRTYP  = 1
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IVERS,
     &            IERROR)
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9010)
        GOTO 9998
      ENDIF
C
C Temps (par securite)
C
      RUBRIQ = 'nombre_iterations_Lagrangiennes'
      IRTYP  = 1
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IPLAS,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
     &  'nombre_iterations_Lagrangiennes                             ',
     &  'IPLAS', IPLAS
      ENDIF
C
      RUBRIQ = 'temps_physique_Lagrangien'
      IRTYP  = 2
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TTCLAG,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9021)
     &  'temps_physique_Lagrangien                                   ',
     &  'TTCLAG', TTCLAG
      ENDIF
C
C Infos sur le suivi du calcul
C
      IRTYP  = 1
C
      RUBRIQ = 'nombre_courant_particules'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NBPART,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9030)
     &  'nombre_courant_particules                                   ',
     &  'NBPART', NBPART
        GOTO 9998
      ENDIF
C
      RUBRIQ = 'nombre_total_particules'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NBPTOT,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
     &  'nombre_total_particules                                     ',
     &  'NBPTOT', NBPTOT
      ENDIF
C
      RUBRIQ = 'nombre_particules_perdues'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NBPERT,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
     &  'nombre_particules_perdues                                   ',
     &  'NBPERT', NBPERT
      ENDIF
C
      RUBRIQ = 'indicateur_physique_particules'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IPHYLA,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9030)
     &  'indicateur_physique_particules                              ',
     &  'IPHYLA', IPHYLA
        GOTO 9998
      ENDIF
C
      RUBRIQ = 'indicateur_temperature_particules'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,ITPVAR,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9030)
     &  'indicateur_temperature_particules                           ',
     &  'ITPVAR', ITPVAR
        GOTO 9998
      ENDIF
C
      RUBRIQ = 'indicateur_diametre_particules'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IDPVAR,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
     &  'indicateur_diametre_particules                              ',
     &  'IDPVAR', IDPVAR
      ENDIF
C
      RUBRIQ = 'indicateur_masse_particules'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IMPVAR,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
     &  'indicateur_masse_particules                                 ',
     &  'IMPVAR', IMPVAR
      ENDIF
C
      RUBRIQ = 'nombre_variables_utilisateur'
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NVLS,
     &            IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
     &  'nombre_variables_utilisateur                                ',
     &  'NVLS', NVLS
      ENDIF
C
      WRITE(NFECRA,6012)
C
C Variables particulaires
C
      NOMNVL(JXP) = 'variable_positionX_particule'
      NOMNVL(JYP) = 'variable_positionY_particule'
      NOMNVL(JZP) = 'variable_positionZ_particule'
      NOMNVL(JUP) = 'variable_vitesseU_particule'
      NOMNVL(JVP) = 'variable_vitesseV_particule'
      NOMNVL(JWP) = 'variable_vitesseW_particule'
      NOMNVL(JUF) = 'variable_vitesseU_fluide_vu'
      NOMNVL(JVF) = 'variable_vitesseV_fluide_vu'
      NOMNVL(JWF) = 'variable_vitesseW_fluide_vu'
      NOMNVL(JMP) = 'variable_masse_particule'
      NOMNVL(JDP) = 'variable_diametre_particule'
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
        NOMNVL(JTP) = 'variable_temperature_particule'
        NOMNVL(JTF) = 'variable_temperature_fluide_vu'
        NOMNVL(JCP) = 'variable_chaleur_specifique_particule'
      ELSEIF (IPHYLA.EQ.2) THEN
        NOMNVL(JHP) = 'variable_temperature_particule'
        NOMNVL(JTF) = 'variable_temperature_fluide_vu'
        NOMNVL(JMCH) = 'variable_masse_charbon_reactif'
        NOMNVL(JMCK) = 'variable_masse_coke'
        NOMNVL(JCP) = 'variable_chaleur_specifique_particule'
      ENDIF
      IF (NVLS.GT.0) THEN
        DO II = 1,NVLS
          WRITE(CAR4,'(I4.4)') II
          NOMNVL(JVLS(II)) = 'variable_supplementaire_'//CAR4
        ENDDO
      ENDIF
C
      ITYSUP = 0
      NBVAL  = NBPART
      IRTYP  = 2
C
      DO II = JMP,JWF
        RUBRIQ = NOMNVL(II)
        CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ETTP(1,II),IERROR)
        IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
          WRITE(NFECRA,9100) RUBRIQ
          GOTO 9998
        ENDIF
      ENDDO
C
      DO II = 1,JMP-1
        RUBRIQ = NOMNVL(II)
        CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ETTP(1,II),IERROR)
        IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
          WRITE(NFECRA,9101) RUBRIQ
        ENDIF
      ENDDO
C
C Caracteristiques et infos particulaires (ENTIERS)
C
      NOMITE(JISOR) = 'numero_cellule_particules'
      IF (IPHYLA.EQ.2) THEN
        NOMITE(JINCH) = 'numero_charbon'
      ENDIF
C
      ITYSUP = 0
      NBVAL  = NBPART
      IRTYP  = 1
C
      RUBRIQ = NOMITE(JISOR)
      CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            ITEPA(1,JISOR),IERROR)
      IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
        WRITE(NFECRA,9100) RUBRIQ
        GOTO 9998
      ENDIF
C
C groupe statistique particules
C
      IF (NBCLST .GT. 0 ) THEN
        NOMITE(JCLST) = 'numero_groupe_statistiques'
C
        ITYSUP = 0
        NBVAL  = NBPART
        IRTYP  = 1
C
        RUBRIQ = NOMITE(JCLST)
        CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ITEPA(1,JCLST),IERROR)
        IF(IERROR.NE.0) THEN
C           advienne que pourra sur le format
          WRITE(NFECRA,9100) RUBRIQ
          GOTO 9998
        ENDIF
      ENDIF
C
C Numero du charbon des particules
C
      IF (IPHYLA.EQ.2) THEN
        RUBRIQ = NOMITE(JINCH)
        CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ITEPA(1,JINCH),IERROR)
        IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
          WRITE(NFECRA,9101) RUBRIQ
        ENDIF
      ENDIF
C
C Caracteristiques et infos particulaires (REELS)
C
      NOMRTE(JRTSP) = 'temps_sejour_particules'
      NOMRTE(JRPOI) = 'poids_statistiques_particules'
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1 .AND.IIRAYO.EQ.1) THEN
        NOMRTE(JREPS) = 'emissivite_particules'
      ENDIF
      IF (IPHYLA.EQ.2) THEN
        NOMRTE(JRDCK) = 'diametre_coeur_retrecissant_charbon'
        NOMRTE(JRD0P) = 'diametre_initial_charbon'
        NOMRTE(JRR0P) = 'masse_volumique_initial_charbon'
      ENDIF
C
      ITYSUP = 0
      NBVAL  = NBPART
      IRTYP  = 2
C
      DO II = 1, NVEP
        RUBRIQ = NOMRTE(II)
        CALL ECRSUI(IMPAVL,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              TEPA(1,II),IERROR)
        IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
          WRITE(NFECRA,9101) RUBRIQ
        ENDIF
      ENDDO
C
      WRITE(NFECRA,6013)
C
C ---> Fermeture du fichier suite
      CALL CLSSUI(IMPAVL,IERROR)
C
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9140) FICAVL
      ENDIF
C
C ---> En cas d'erreur, on continue quand meme
 9998 CONTINUE
C
      WRITE(NFECRA,6014)
C
C=======================================================================
C 2. ECRITURE DU FICHIER SUITE STATISTIQUES ET TERMES SOURCES
C    DE COUPLAGE RETOUR
C=======================================================================
C
C
      IF ( (ISTALA.EQ.1 .AND. IPLAS.GE.IDSTNT) .OR.
     &      IILAGR.EQ.2                        .OR.
     &     (IENSI3.EQ.1 .AND. NVISBR.GT.0)          ) THEN
C
C
C ---> Ouverture (et on saute si erreur)
C     IFOVLS    =0 --> non formate       =1 --> formate
C     ILECEC = 2 : ecriture
C
        IF (IFOVLS.EQ.0) THEN
          WRITE(NFECRA,7010)
        ELSE IF (IFOVLS.EQ.1) THEN
          WRITE(NFECRA,7020)
        ELSE
          WRITE(NFECRA,9011) FICVLS, 'IFOVLS', IFOVLS
        ENDIF
C
        ILECEC = 2
        CALL OPNSUI(FICVLS, LEN(FICVLS), ILECEC, IFOVLS, IMPVLS, IERROR)
C       ===========
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9510) FICVLS
          GOTO 9999
        ENDIF
C
        WRITE(NFECRA,7011)
C
C
C Entete et Infos sur le calcul ou on saute si erreur
C     On inclut une rubrique destinee a distinguer ce fichier
C       d'un autre fichier suite
C     Pour le moment, IVERS n'est pas utilise
C
        ITYSUP = 0
        NBVAL  = 1
C
        IVERS  = 111
        RUBRIQ = 'version_fichier_suite_Lagrangien_statistiques'
        IRTYP  = 1
        CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IVERS,
     &              IERROR)
C
        IF(IERROR.NE.0) THEN
          WRITE(NFECRA,9510)
          GOTO 9999
        ENDIF
C
C ---> On ecrit ISTTIO c'est utile dans tous les cas
C
        RUBRIQ = 'indicateur_ecoulement_stationnaire'
        IRTYP  = 1
        CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ISTTIO, IERROR)
C
        IF(IERROR.NE.0) THEN
          WRITE(NFECRA,9510)
          GOTO 9999
        ENDIF
C
C --> En premier, on ecrit les statistiques volumiques
C
        IF (ISTALA.EQ.1 .AND. IPLAS.GE.IDSTNT) THEN
C
          RUBRIQ = 'iteration_debut_statistiques'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IDSTNT,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'iteration_debut_statistiques                                ',
     &      'IDSTNT', IDSTNT
          ENDIF
C
          RUBRIQ = 'iteration_debut_statistiques_stationnaires'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NSTIST,IERROR)
          IF(IERROR.NE.0)  THEN
            WRITE(NFECRA,9520)
     &  'iteration_debut_statistiques_stationnaires                  ',
     &      'NSTIST', NSTIST
          ENDIF
C
          RUBRIQ = 'nombre_iterations_statistiques_stationnaires'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPST,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'nombre_iterations_statistiques_stationnaires                ',
     &      'NPST', NPST
          ENDIF
C
          RUBRIQ = 'temps_statistiques_stationnaires'
          IRTYP  = 2
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TSTAT,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9521)
     &  'temps_statistiques_stationnaires                            ',
     &      'TSTAT', TSTAT
          ENDIF
C
          RUBRIQ = 'classe_statistique_particules'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NBCLST,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'classes_statistiques                                        ',
     &      'NBCLST', NBCLST
          ENDIF
C
          RUBRIQ = 'nombre_statistiques_utilisateur'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NVLSTS,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'nombre_statistiques_utilisateur                             ',
     &     'NVLSTS', NVLSTS
          ENDIF
C
C  Statistiques volumiques
C
          ITYSUP = 1
          IRTYP  = 2
          NBVAL  = 1
C
          DO IPAS  = 0,NBCLST
            DO JJ = 1,NVLSTA
C
              II = IPAS*NVLSTA +JJ
              IF (IPAS.GT.0) THEN
                WRITE(CAR4,'(I4.4)') IPAS
                RUBRIQ = 'moy_stat_vol_groupe_'//CAR4//'_'//NOMLAG(II)
              ELSE
                RUBRIQ = 'moy_stat_vol_'//NOMLAG(II)
              ENDIF
              CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  STATIS(1,II),IERROR)
C
              IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
                WRITE(NFECRA,9550) RUBRIQ
              ENDIF
            ENDDO
C
            DO JJ = 1,NVLSTA-1
              IF (IPAS.GT.0) THEN
                WRITE(CAR4,'(I4.4)') IPAS
                RUBRIQ = 'var_stat_vol_groupe_'//CAR4//'_'//NOMLAV(II)
              ELSE
                RUBRIQ = 'var_stat_vol_'//NOMLAV(II)
              ENDIF
              CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  STATIV(1,II),IERROR)
C
              IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
                WRITE(NFECRA,9550) RUBRIQ
              ENDIF
C
            ENDDO
C
          ENDDO
C
        ENDIF
C
C --> En second, c'est le tour des statistiques aux frontieres
C
        IF (IENSI3.EQ.1 .AND. NVISBR.GT.0) THEN
C
          ITYSUP = 0
          NBVAL  = 1
C
          RUBRIQ = 'iteration_debut_stats_frontieres_stationnaires'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NSTBOR,IERROR)
          IF(IERROR.NE.0)  THEN
            WRITE(NFECRA,9520)
     &  'iteration_debut_stats_frontieres_stationnaires              ',
     &      'NSTBOR', NSTBOR
          ENDIF
C
          RUBRIQ = 'nombre_iterations_stats_frontieres'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPSTFT,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'nombre_iterations_stats_frontieres                          ',
     &      'NPSTFT', NPSTFT
          ENDIF
C
          RUBRIQ = 'nombre_iterations_stats_frontieres_stationnaires'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPSTF,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'nombre_iterations_stats_frontieres_stationnaires            ',
     &      'NPSTF', NPSTF
          ENDIF
C
          RUBRIQ = 'temps_stats_frontieres_stationnaires'
          IRTYP  = 2
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TSTATP,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9521)
     &  'temps_stats_frontieres_stationnaires                        ',
     &      'TSTATP', TSTATP
          ENDIF
C
          RUBRIQ = 'nombre_stats_frontieres_utilisateur'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NUSBOR,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9521)
     &  'nombre_stats_frontieres_utilisateur                         ',
     &      'NUSBOR', NUSBOR
          ENDIF
C
C  Statistiques aux frontieres
C
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 2
C
          DO II = 1,NVISBR
            RUBRIQ = 'stat_bord_'//NOMBRD(II)
            CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PARBOR(1,II),IERROR)
            IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
              WRITE(NFECRA,9550) RUBRIQ
            ENDIF
          ENDDO
C
        ENDIF
C
C --> Enfin, en cas de couplage retour, on ecrit les termes sources
C
        IF (IILAGR.EQ.2) THEN
C
          ITYSUP = 0
          NBVAL  = 1
C
          RUBRIQ = 'iteration_debut_termes_sources_stationnaires'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NSTITS,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'iteration_debut_termes_sources_stationnaires                ',
     &      'NSTITS', NSTITS
          ENDIF
C
          RUBRIQ = 'nombre_iterations_termes_sources_stationnaires'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPTS,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'nombre_iterations_termes_sources_stationnaires              ',
     &      'NPTS', NPTS
          ENDIF
C
          RUBRIQ = 'modele_turbulence_termes_sources'
          IRTYP  = 1
          CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                ITURB(IPHAS),IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9520)
     &  'modele_turbulence_termes_sources                            ',
     &      'ITURB',ITURB(IPHAS)
          ENDIF
C
C       On donne des labels au different TS pour les noms de rubriques
C       On donne le meme label au keps, au v2f et au k-omega (meme variable k)
C
          IF (LTSDYN.EQ.1) THEN
            NOMTSL(ITSVX) = 'terme_source_vitesseX'
            NOMTSL(ITSVY) = 'terme_source_vitesseY'
            NOMTSL(ITSVZ) = 'terme_source_vitesseZ'
            NOMTSL(ITSLI) = 'terme_source_vitesse_implicite'
            IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &           .OR. ITURB(IPHAS).EQ.60) THEN
              NOMTSL(ITSKE) = 'terme_source_turbulence_keps'
            ELSE IF (ITYTUR(IPHAS).EQ.3) THEN
              NOMTSL(ITSR11) = 'terme_source_turbulence_R11'
              NOMTSL(ITSR12) = 'terme_source_turbulence_R12'
              NOMTSL(ITSR13) = 'terme_source_turbulence_R13'
              NOMTSL(ITSR22) = 'terme_source_turbulence_R22'
              NOMTSL(ITSR23) = 'terme_source_turbulence_R23'
              NOMTSL(ITSR33) = 'terme_source_turbulence_R33'
            ENDIF
          ENDIF
          IF (LTSMAS.EQ.1) THEN
            NOMTSL(ITSMAS) = 'terme_source_masse'
          ENDIF
          IF (LTSTHE.EQ.1) THEN
            IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
              NOMTSL(ITSTE) = 'terme_source_thermique_explicite'
              NOMTSL(ITSTI) = 'terme_source_thermique_implicite'
            ELSE IF (IPHYLA.EQ.2) THEN
              NOMTSL(ITSTE) = 'terme_source_thermique_explicite'
              NOMTSL(ITSTI) = 'terme_source_thermique_implicite'
              DO ICHA = 1,NCHARB
                WRITE(CAR4,'(I4.4)') ICHA
                NOMTSL(ITSMV1(ICHA)) = 'terme_source_legeres_F1_'//CAR4
                NOMTSL(ITSMV2(ICHA)) = 'terme_source_lourdes_F2_'//CAR4
              ENDDO
              NOMTSL(ITSCO) = 'terme_source_F3'
              NOMTSL(ITSFP4) = 'terme_source_variance_traceur_air'
            ENDIF
          ENDIF
C
C  Termes source de couplage retour
C
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
C
          DO II = 1,NTERSL
            RUBRIQ = NOMTSL(II)
            CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  TSLAGR(1,II),IERROR)
            IF(IERROR.NE.0) THEN
C         advienne que pourra sur le format
              WRITE(NFECRA,9550) RUBRIQ
            ENDIF
          ENDDO
C
C  Dans le cas specifique de la combustion de grains de charbon
C  avec un couplage retour sur une combustion gaz en phase porteuse
C
C      --> A verifier l'utilite de cette sauvegarde pour une suite...
C
          IF (IPPMOD(ICPL3C).GE.0) THEN
            DO II = 1, NSALPP
              ICHA = NSALTO-NSALPP+II
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              WRITE(CAR4,'(I4.4)') II
              RUBRIQ = 'scalaires_physiques_pariculieres_charbon'//CAR4
              CALL ECRSUI(IMPVLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPPROC(ICHA)),IERROR)
              IF(IERROR.NE.0) WRITE(NFECRA,9550) RUBRIQ
            ENDDO

          ENDIF
C
        ENDIF
C
        WRITE(NFECRA,7013)
C
C ---> Fermeture du fichier suite
        CALL CLSSUI(IMPVLS,IERROR)
C
        IF(IERROR.NE.0) THEN
          WRITE(NFECRA,9700) FICVLS
        ENDIF
C
C ---> En cas d'erreur, on continue quand meme
 9999   CONTINUE
C
        WRITE(NFECRA,7014)
C
      ENDIF
C
C=======================================================================
C 3. Visualisations
C=======================================================================
C
      IF (NTCABS.LT.NTMABS) RETURN
C
      NFIN = 1
C
C-->Stockage des trajectoires au format Ensight Gold
C
      IF (IENSI1.EQ.1) THEN
C
        IFORCE = 0
        CALL ENSLAG
C       ===========
     &   ( IDEBIA , IDEBRA ,
     &     NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NFIN   , IFORCE ,
     &     ITEPA  ,
     &     ETTP   , TEPA   , RA )
C
      ENDIF
C
C-->Stockage des deplacements au format Ensight Gold
C
      IF (IENSI2.EQ.1) THEN
C
        IFINIA = IDEBIA
        ITRAV1 = IDEBRA
        IFINRA = ITRAV1 + 3*NBPMAX
        CALL RASIZE ('LAGOUT',IFINRA)
C       ===========
C
        CALL ENSWAF
C       ===========
     &   ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     NFIN   ,
     &     ITEPA  ,
     &     ETTP   , TEPA , RA(ITRAV1)   )
C
      ENDIF
C
C=======================================================================
C 4. FIN
C=======================================================================
C
      RETURN
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
 6010 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN       ',/,
     &           3X,'   -------------------------------------       ',/,
     &           3X,' Ecriture d''un fichier suite binaire          ',/,
     &           3X,'   sur les variables liees aux particules      ',/)
C
 6020 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN       ',/,
     &           3X,'   -------------------------------------       ',/,
     &           3X,' Ecriture d''un fichier suite formate (ASCII)  ',/,
     &           3X,'   sur les variables liees aux particules      ',/)
C
 6011 FORMAT (   3X,'   Debut de l''ecriture                        '  )
 6012 FORMAT (   3X,'   Fin de l''ecriture des infos sur le calcul  '  )
 6013 FORMAT (   3X,'   Fin de l''ecriture des infos particulaires  '  )
 6014 FORMAT (   3X,' Fin de l''ecriture du fichier suite           ',/,
     &           3X,'   sur les variables liees aux particules      ',/)
C
 7010 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN       ',/,
     &           3X,'   -------------------------------------       ',/,
     &           3X,' Ecriture d''un fichier suite binaire          ',/,
     &           3X,'   sur les statistiques volumiques et aux      ',/,
     &           3X,'   fontieres, ainsi que les termes sources     ',/,
     &           3X,'   de couplage retour                          ',/)
C
 7020 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN       ',/,
     &           3X,'   -------------------------------------       ',/,
     &           3X,' Ecriture d''un fichier suite formate (ASCII)  ',/,
     &           3X,'   sur les statistiques volumiques et aux      ',/,
     &           3X,'   fontieres, ainsi que les termes sources     ',/,
     &           3X,'   de couplage retour                          ',/)
C
 7011 FORMAT (   3X,'   Debut de l''ecriture des stats et TS        '  )
 7013 FORMAT (   3X,'   Fin de l''ecriture des statistiques et TS   '  )
 7014 FORMAT (   3X,' Fin de l''ecriture du fichier suite           ',/,
     &           3X,'   sur les statistiques et TS couplage retour  ',/)
C
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@    ERREUR A L''OUVERTURE EN ECRITURE DU FICHIER SUITE      ',/,
     &'@      (',A13,')                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul se termine mais ne fournira pas de fichier      ',/,
     &'@    suite sur les caracteristiques des particules.          ',/,
     &'@                                                            ',/,
     &'@  Verifier que le repertoire de travail est accessible en   ',/,
     &'@    ecriture et que le fichier suite peut y etre cree.      ',/,
     &'@  Voir le sous-programme LAGOUT.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@    ERREUR A L''OUVERTURE EN ECRITURE DU FICHIER SUITE      ',/,
     &'@      (',A13,')                                             ',/,
     &'@                                                            ',/,
     &'@  Le type du fichier ',A6,' doit valoir 0 ou 1              ',/,
     &'@    il vaut ici  ',I3                                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul se termine mais ne fournira pas de fichier      ',/,
     &'@    suite sur les caracteristiques des particules.          ',/,
     &'@                                                            ',/,
     &'@  Verifier que le repertoire de travail est accessible en   ',/,
     &'@    ecriture et que le fichier suite peut y etre cree.      ',/,
     &'@  Voir le sous-programme LAGOUT.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,
     &'@        ',A10   ,'  = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,
     &'@        ',A10   ,'  = ',E14.5                                ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,
     &'@        ',A10   ,'  = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue mais ne fournira pas de fichier      ',/,
     &'@      suite sur les caracteristiques des particules.        ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue mais ne fournira pas de fichier      ',/,
     &'@      suite sur les caracteristiques des particules.        ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9140 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A LA FERMETURE DU FICHIER SUITE LAGRANGIEN   ',/,
     &'@    *********    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom : ',A13                   ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9510 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,
     &'@                 DE COUPLAGE RETOUR                         ',/,
     &'@                                                            ',/,
     &'@    ERREUR A L''OUVERTURE EN ECRITURE DU FICHIER SUITE      ',/,
     &'@      (',A13,')                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul continue mais ne fournira pas de fichier        ',/,
     &'@    suite sur les statistiques volumiques et aux frontieres ',/,
     &'@    ainsi que sur les termes sources de couplage retour     ',/,
     &'@                                                            ',/,
     &'@  Verifier que le repertoire de travail est accessible en   ',/,
     &'@    ecriture et que le fichier suite peut y etre cree.      ',/,
     &'@  Voir le sous-programme LAGOUT.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9520 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,
     &'@                 DE COUPLAGE RETOUR                         ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,
     &'@        ',A10   ,'  = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9521 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,
     &'@                 DE COUPLAGE RETOUR                         ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,
     &'@        ',A10   ,'  = ',E14.5                                ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9550 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,
     &'@    *********    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,
     &'@                 DE COUPLAGE RETOUR                         ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@    Voir le sous-programme LAGOUT.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9700 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A LA FERMETURE DU FICHIER SUITE LAGRANGIEN   ',/,
     &'@    *********    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,
     &'@                 DE COUPLAGE RETOUR                         ',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom : ',A13                   ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
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
