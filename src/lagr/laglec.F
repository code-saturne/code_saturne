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
                        SUBROUTINE LAGLEC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   ITEPA  , IA     ,
     &   RTPA   , PROPCE ,
     &   ETTP   , TEPA   , STATIS , STATIV , PARBOR , TSLAGR , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     Lecture des fichiers suite Lagrangien "lagamo" et "lasamo"
CFONC     contenant les informations sur les particule, les statistiques
CFONC     volumiques et aux frontieres, ainsi que les termes sources
CFONC     de couplage retour.
CFONC
CFONC     Tous les tableaux sont initialise a zero avant d'être remplis
CFONC     dans le cas d'une suite (sinon ils restent a zero).
CFONC     On realise donc ici l'initialisation des tableaux ouverts
CFONC     dans MEMLA1, ce qui termine l'etape d'initialisation debutee
CFONC     dans LAGOPT.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CTODO cartouche
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
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES INSTANT PRECEDENT        !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR !  -> ! CUMUL POUR LES MOYENNES DES          !
CARGU !(NCELET,NVLSTA!    !     !   STATISTIQUES VOLUMIQUES            !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! PARBOR       ! TR !  -> ! INFOS SUR INTERACTION DES PARTICULES !
CARGU !(NFABOR,NVISBR!    !     !   AUX FACES DE BORD                  !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,NTERSL!    !     !   LAGRANGIEN SUR LA PHASE PORTEUSE   !
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
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
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          ITEPA(NBPMAX,NIVEP)  , IA(*)
C
      DOUBLE PRECISION RTPA(NCELET,*) , PROPCE(NCELET,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        RUBRIQ*64 , CAR4*4, CAR8*8, KAR8*8
      CHARACTER        NOMNVL(NVPLMX)*60 , NOMTSL(NVPLMX)*60
      CHARACTER        NOMITE(NVPLMX)*64 , NOMRTE(NVPLMX)*64
      INTEGER          IDEBIA , IDEBRA
      INTEGER          NCELOK , NFAIOK , NFABOK , NSOMOK
      INTEGER          IERROR , IRTYP  , ITYSUP , NBVAL
      INTEGER          ILECEC , NBERRO , IVERS  , IBIDON
      INTEGER          MVLS   , IVAR   , IP     , ICHA
      INTEGER          IFAC   , IEL    , IOK    , IPHAS
      INTEGER          JPHYLA , JTPVAR , JDPVAR , JMPVAR
      INTEGER          JSTTIO , JDSTNT , MSTIST , MVLSTS
      INTEGER          MSTBOR , MUSBOR , MSTITS , JTURB, JTYTUR
      INTEGER          MODE   , IPAS   , IVL    , NCLSTO
C
C***********************************************************************
C=======================================================================
C 0. Gestion memoire
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. Initialisations par defaut
C=======================================================================
C
C---> Il faut faire dans cette routine les initialisations des
C     tableaux lagrangiens ouverts dans la routine MEMLA1
C     (sauf ITYCEL et ICOCEL qui sont initialises dans LAGDEB),
C
      IPHAS = ILPHAS
C
      DO IVAR = 1,NVP
        DO IP = 1,NBPMAX
          ETTP(IP,IVAR) = 0.D0
        ENDDO
      ENDDO
C
      DO IVAR = 1,NIVEP
        DO IP = 1,NBPMAX
          ITEPA(IP,IVAR) = 0
        ENDDO
      ENDDO
C
      DO IVAR = 1,NVEP
        DO IP = 1,NBPMAX
          TEPA(IP,IVAR) = 0.D0
        ENDDO
      ENDDO
C
      IF (ISTALA.EQ.1) THEN
        DO IPAS  = 0,NBCLST
          DO IVL = 1,NVLSTA
            IVAR = IPAS*NVLSTA +IVL
            DO IEL = 1,NCEL
              STATIS(IEL,IVAR) = 0.D0
            ENDDO
          ENDDO
        ENDDO
        DO IPAS  = 0,NBCLST
          DO IVL = 1,NVLSTA-1
            IVAR = IPAS*(NVLSTA-1) +IVL
            DO IEL = 1,NCEL
              STATIV(IEL,IVAR) = 0.D0
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      IF (IILAGR.EQ.2) THEN
        DO IVAR = 1,NTERSL
          DO IEL = 1,NCEL
            TSLAGR(IEL,IVAR) = 0.D0
          ENDDO
        ENDDO
      ENDIF
C
      IF (IENSI3.EQ.1 .AND. NVISBR.GT.0) THEN
        DO IVAR = 1,NVISBR
          DO IFAC = 1,NFABOR
            PARBOR(IFAC,IVAR) = 0.D0
          ENDDO
        ENDDO
      ENDIF
C
      IF (ISUILA.EQ.0) RETURN
C
C=======================================================================
C 2. LECTURE DU FICHIER SUITE : VARIABLES LIEES AUX PARTICULES
C=======================================================================
C
C  ---> Ouverture
C
      WRITE(NFECRA,6000)
C
C     (ILECEC=1:lecture)
      ILECEC = 1
      CALL OPNSUI(FICAML,LEN(FICAML),ILECEC,IBIDON,IMPAML,IERROR)
C     ===========
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9010) FICAML, FICAML
        CALL CSEXIT (1)
      ENDIF
C
      WRITE(NFECRA,6010)
C
C  ---> Type de fichier suite
C        Pourrait porter le numero de version si besoin.
C        On ne se sert pas de IVERS pour le moment
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      RUBRIQ = 'version_fichier_suite_Lagrangien_variables'
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            IVERS,IERROR)
C
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9020) FICAML, FICAML
        CALL CSEXIT (1)
      ENDIF
C
C  ---> Tests
C
      IOK = 0
C
C     Dimensions des supports
C
      CALL TSTSUI(IMPAML,NCELOK,NFAIOK,NFABOK,NSOMOK)
C     ===========
      IF(NCELOK.EQ.0) THEN
        WRITE(NFECRA,9030) FICAML
        IOK = IOK + 1
      ENDIF
C
      IF(NFAIOK.EQ.0) WRITE(NFECRA,9031) FICAML,'internes','internes'
C
      IF(NFABOK.EQ.0) WRITE(NFECRA,9031) FICAML,'de bord ','de bord '
C
C     Nombre de particules dans le domaine du calcul
C
      ITYSUP = 0
      NBVAL  = 1
C
      RUBRIQ = 'nombre_courant_particules'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            NBPART,IERROR)
C
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9040) FICAML,
     &  'nombre_courant_particules                                   ',
     &  FICAML
        IOK = IOK + 1
      ENDIF
      IF(NBPART.GT.NBPMAX) THEN
        WRITE(NFECRA,9050) FICAML, NBPART, NBPMAX
        IOK = IOK + 1
      ENDIF
C
C     Physique associee aux particules
C
      RUBRIQ = 'indicateur_physique_particules'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            JPHYLA,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9040) FICAML, RUBRIQ, FICAML
        IOK = IOK + 1
      ENDIF
C
      RUBRIQ = 'indicateur_temperature_particules'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            JTPVAR,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9040) FICAML,
     &  'indicateur_temperature_particules                           ',
     &  FICAML
        IOK = IOK + 1
      ENDIF
C
C     Arret
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
      RUBRIQ = 'indicateur_diametre_particules'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            JDPVAR,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9062) FICAML,
     &  'indicateur_diametre_particules                              '
        JDPVAR = IDPVAR
      ENDIF
C
      RUBRIQ = 'indicateur_masse_particules'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            JMPVAR,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9062) FICAML,
     &  'indicateur_masse_particules                                 '
        JMPVAR = IMPVAR
      ENDIF
C
C ---> On previent si des parametres sont differents
C
      IF ( JPHYLA.NE.IPHYLA .OR.
     &     JTPVAR.NE.ITPVAR .OR.
     &     JDPVAR.NE.IDPVAR .OR.
     &     JMPVAR.NE.IMPVAR      ) THEN
        WRITE(NFECRA,9070) FICAML,
     &                     JPHYLA, JTPVAR, JDPVAR, JMPVAR,
     &                     IPHYLA, ITPVAR, IDPVAR, IMPVAR
      ENDIF
C
C ---> Verification de la compatibilite si changement de thermique
C
      IF (JPHYLA.NE.0 .AND. IPHYLA.EQ.0) THEN
        WRITE(NFECRA,9071) FICAML
      ENDIF
C
      IF (ITPVAR.EQ.1 .AND. JTPVAR.EQ.0) THEN
        WRITE(NFECRA,9072) FICAML, TPART, CPPART
      ENDIF
C
      IF (IPHYLA.EQ.2 .AND. JPHYLA.NE.2) THEN
        WRITE(NFECRA,9073) FICAML, FICAML
        CALL CSEXIT (1)
      ENDIF
C
      IF ( (JPHYLA.EQ.2 .AND. IPHYLA.EQ.1) .OR.
     &     (JPHYLA.EQ.1 .AND. IPHYLA.EQ.2)      ) THEN
        WRITE(NFECRA,9074) FICAML, FICAML
        CALL CSEXIT (1)
      ENDIF
C
C ---> Infos suivi du calcul
C
      RUBRIQ = 'nombre_iterations_Lagrangiennes'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            IPLAS,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9060) FICAML,
     &  'nombre_iterations_Lagrangiennes                             ',
     &  'IPLAS',IPLAS
      ENDIF
C
      IF(ISTALA.EQ.1 .AND. ISUIST.EQ.0 .AND. IPLAS.GE.IDSTNT) THEN
        WRITE(NFECRA,9065) FICAML, ISUIST, IPLAS +1, IDSTNT
        CALL CSEXIT (1)
      ENDIF
C
      IF(IENSI3.EQ.1 .AND. ISUIST.EQ.0 .AND. IPLAS.GE.NSTBOR) THEN
        WRITE(NFECRA,9066) FICAML, ISUIST, IPLAS +1, NSTBOR
        CALL CSEXIT (1)
      ENDIF
C
      RUBRIQ = 'temps_physique_Lagrangien'
      IRTYP  = 2
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            TTCLAG,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9061) FICAML,
     &  'temps_physique_Lagrangien                                   ',
     &  'TTCLAG',TTCLAG
      ENDIF
C
      RUBRIQ = 'nombre_total_particules'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            NBPTOT,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9060) FICAML,
     &  'nombre_total_particules                                     ',
     &  'NBPTOT',NBPTOT
      ENDIF
C
      RUBRIQ = 'nombre_particules_perdues'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            NBPERT,IERROR)
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9060) FICAML,
     &  'nombre_particules_perdues                                   ',
     &  'NBPERT',NBPERT
      ENDIF
C
      RUBRIQ = 'nombre_variables_utilisateur'
      IRTYP  = 1
      CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            MVLS,IERROR)
      IF(IERROR.NE.0) THEN
        MVLS = 0
        IF (NVLS.GT.0) THEN
          WRITE(NFECRA,9062) FICAML,
     &  'nombre_variables_utilisateur                                '
        ENDIF
      ENDIF
C
      IF (NVLS.LT.MVLS) THEN
        WRITE(NFECRA,9080) FICAML, MVLS, NVLS, NVLS, NVLS
        MVLS = NVLS
      ELSEIF (NVLS.GT.MVLS ) THEN
        WRITE(NFECRA,9080) FICAML, MVLS, NVLS, NVLS, NVLS
      ENDIF
C
C --> Caracteristiques et infos particulaires (ENTIERS)
C
      NBERRO = 0
C
      NOMITE(JISOR) = 'numero_cellule_particules'
      IF (NBCLST.GT.0) THEN
        NOMITE(JCLST) = 'numero_groupe_statistiques'
      ENDIF
      IF (IPHYLA.EQ.2) THEN
        NOMITE(JINCH) = 'numero_charbon'
      ENDIF
C
      NBVAL  = NBPART
      IRTYP  = 1
C
      DO IVAR = 1, NIVEP
        RUBRIQ = NOMITE(IVAR)
        CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ITEPA(1,IVAR),IERROR)
        NBERRO = NBERRO+IERROR
      ENDDO
C
C --> Caracteristiques et infos particulaires (REELS)
C
      NOMRTE(JRTSP) = 'temps_sejour_particules'
      NOMRTE(JRPOI) = 'poids_statistiques_particules'
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1 .AND. IIRAYO.EQ.1) THEN
        NOMRTE(JREPS) = 'emissivite_particules'
      ENDIF
      IF (IPHYLA.EQ.2) THEN
        NOMRTE(JRDCK) = 'diametre_coeur_retrecissant_charbon'
        NOMRTE(JRD0P) = 'diametre_initial_charbon'
        NOMRTE(JRR0P) = 'masse_volumique_initial_charbon'
      ENDIF
C
      NBVAL  = NBPART
      IRTYP  = 2
C
      DO IVAR = 1, NVEP
        RUBRIQ = NOMRTE(IVAR)
        CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              TEPA(1,IVAR),IERROR)
        NBERRO = NBERRO+IERROR
      ENDDO
C
      DNBPAR = 0.D0
      DO IP  = 1, NBPART
        DNBPAR = DNBPAR + TEPA(IP,JRPOI)
      ENDDO
C
C --> Variables particulaires
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
      IF (MVLS.GT.0) THEN
        DO IP = 1,MVLS
          WRITE(CAR4,'(I4.4)') IP
          NOMNVL(JVLS(IP)) = 'variable_supplementaire_'//CAR4
        ENDDO
      ENDIF
C
      NBVAL  = NBPART
      IRTYP  = 2
C
      DO IVAR = JMP, JWF
        RUBRIQ = NOMNVL(IVAR)
        CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ETTP(1,IVAR),IERROR)
        NBERRO = NBERRO+IERROR
      ENDDO
C
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN

        IF (ITPVAR.EQ.1 .AND. JTPVAR.EQ.0) THEN
C
          DO IP = 1,NBPART
            IEL = ITEPA(IP,JISOR)
            ETTP(IP,JTP) = TPART
            ETTP(IP,JCP) = CPPART
            IF ( IPPMOD(ICP3PL).GE.0 .OR.
     &           IPPMOD(ICPL3C).GE.0 .OR.
     &           IPPMOD(ICFUEL).GE.0      ) THEN
               ETTP(IP,JTF) = PROPCE(IEL,IPPROC(ITEMP1)) -TKELVI
            ELSE IF ( IPPMOD(ICOD3P).GE.0 .OR.
     &                IPPMOD(ICOEBU).GE.0 .OR.
     &                IPPMOD(IELARC).GE.0 .OR.
     &                IPPMOD(IELJOU).GE.0      ) THEN
               ETTP(IP,JTF) = PROPCE(IEL,IPPROC(ITEMP)) -TKELVI
            ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.1 ) THEN
               ETTP(IP,JTF) = RTPA(IEL,ISCA(ISCALT(IPHAS))) -TKELVI
            ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.-1 ) THEN
               ETTP(IP,JTF) = RTPA(IEL,ISCA(ISCALT(IPHAS)))
            ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.2 ) THEN
               MODE = 1
               CALL USTHHT(MODE, RTPA(IEL,ISCA(ISCALT(IPHAS))),
     &                     ETTP(IP,JTF))
            ENDIF
          ENDDO
C
        ELSE IF (ITPVAR.EQ.1 .AND. JTPVAR.EQ.1) THEN
          DO IVAR = JTP,JCP
            RUBRIQ = NOMNVL(IVAR)
            CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  ETTP(1,IVAR),IERROR)
            NBERRO = NBERRO+IERROR
          ENDDO
        ENDIF
C
      ELSE IF (IPHYLA.EQ.2) THEN
        DO IVAR = JHP,JCP
          RUBRIQ = NOMNVL(IVAR)
          CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                ETTP(1,IVAR),IERROR)
          NBERRO = NBERRO+IERROR
        ENDDO
      ENDIF
C
      IF (MVLS.GT.0) THEN
        DO IVAR = 1,MVLS
          RUBRIQ = NOMNVL(JVLS(IVAR))
          CALL LECSUI(IMPAML,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                ETTP(1,JVLS(IVAR)),IERROR)
          NBERRO = NBERRO+IERROR
        ENDDO
      ENDIF
C
C  ---> Si pb : arret
C
C
      IF(NBERRO.NE.0) THEN
        WRITE(NFECRA,9041) FICAML, FICAML
        CALL CSEXIT (1)
      ENDIF
C
      WRITE(NFECRA,6011)

C  ---> Fermeture du fichier suite
C
      CALL CLSSUI(IMPAML,IERROR)
C
      IF(IERROR.NE.0) THEN
        WRITE(NFECRA,9090) FICAML, FICAML
      ENDIF
C
      WRITE(NFECRA,6099)
C
C
C=======================================================================
C 3. LECTURE DU FICHIER SUITE STATISTIQUES ET TERMES SOURCES
C    DE COUPLAGE RETOUR
C=======================================================================
C
      IF (ISUIST.EQ.1) THEN
C
C  ---> Ouverture
C
        WRITE(NFECRA,7000)
C
C     (ILECEC=1:lecture)
        ILECEC = 1
        CALL OPNSUI(FICMLS,LEN(FICMLS),ILECEC,IBIDON,IMPMLS,IERROR)
C       ===========
        IF(IERROR.NE.0) THEN
          WRITE(NFECRA,9010) FICMLS, FICMLS
          CALL CSEXIT (1)
        ENDIF
C
        WRITE(NFECRA,7010)
C
C  ---> Type de fichier suite
C        Pourrait porter le numero de version si besoin.
C        On ne se sert pas de IVERS pour le moment
C
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
C
        RUBRIQ = 'version_fichier_suite_Lagrangien_statistiques'
        CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              IVERS,IERROR)
        IF(IERROR.NE.0) THEN
          WRITE(NFECRA,9020) FICAML, FICAML
          CALL CSEXIT (1)
        ENDIF
C
        RUBRIQ = 'indicateur_ecoulement_stationnaire'
        CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              JSTTIO,IERROR)
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9040) FICAML,
     &  'indicateur_ecoulement_stationnaire                          ',
     &    FICAML
          CALL CSEXIT (1)
        ENDIF
C
C     Dimensions des supports
C
        CALL TSTSUI(IMPMLS,NCELOK,NFAIOK,NFABOK,NSOMOK)
C       ===========
        IF(NCELOK.EQ.0) THEN
          WRITE(NFECRA,9030) FICMLS
          CALL CSEXIT (1)
        ENDIF
C
        IF(NFAIOK.EQ.0) WRITE(NFECRA,9031) FICMLS,'internes','internes'
C
        IF(NFABOK.EQ.0) WRITE(NFECRA,9031) FICMLS,'de bord ','de bord '
C
C
C --> Est-on cense lire une suite de stats volumiques ?
C
        IF (ISTALA.EQ.1 .AND. IPLAS.GE.IDSTNT) THEN
C
          NBERRO = 0
          ITYSUP = 0
          NBVAL  = 1
C
          RUBRIQ = 'iteration_debut_statistiques'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                JDSTNT,IERROR)
          NBERRO = NBERRO+IERROR
C
          RUBRIQ = 'iteration_debut_statistiques_stationnaires'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                MSTIST,IERROR)
          NBERRO = NBERRO+IERROR
C
C  ---> S'il y a des erreurs, on suppose que c'est parce que le fichier
C         suite ne contient pas d'infos sur les stats volumiques.
C         Dans ce cas, si on est en instationnaire on se dit que c'est
C         pas grave, on saute l'etape et on continue. Par contre si on
C         est dans une configuration de calcul de stats volumiques
C         en stationnaire on stoppe.
C
          IF(NBERRO.NE.0) THEN
            IF ( ISTTIO.EQ.0 .OR.
     &          (ISTTIO.EQ.1 .AND. IPLAS.LT.NSTIST) ) THEN
              WRITE(NFECRA,9110) FICMLS, ISTTIO, IDSTNT, NSTIST, IPLAS+1
              GOTO 9991
            ELSE
              WRITE(NFECRA,9120) FICMLS, ISTTIO, IDSTNT, NSTIST, IPLAS+1
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
C --> A partir d'ici on considere que le fichier suite contient
C       des stats volumiques
C
          IF ( JSTTIO.NE.ISTTIO .OR.
     &         JDSTNT.NE.IDSTNT .OR.
     &         MSTIST.NE.NSTIST     ) THEN
            WRITE (NFECRA,9130) FICMLS,
     &                          JSTTIO, JDSTNT, MSTIST,
     &                          ISTTIO, IDSTNT, NSTIST
          ENDIF
C
C  --> Lecture de l'avancement du calcul stats volumiques
C
          RUBRIQ = 'nombre_iterations_statistiques_stationnaires'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPST,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9060) FICMLS,
     &  'nombre_iterations_statistiques_stationnaires                ',
     &      'NPST',NPST
          ENDIF
C
          RUBRIQ = 'temps_statistiques_stationnaires'
          IRTYP  = 2
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TSTAT,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9061) FICMLS,
     &  'temps_statistiques_stationnaires                            ',
     &      'TSTAT',TSTAT
          ENDIF
C
          RUBRIQ = 'classe_statistique_particules'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NCLSTO,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9061) FICMLS,
     &  'classes_statistiques                                        ',
     &      'NBCLST',NCLSTO
          ENDIF
C
C  --> Verif de coherence de l'avancement du calcul avec les
C       indicateurs de calcul de la suite actuelle :
C
C    1) Amont Instationnaire -> Actuel Instationnaire : OK
C         (NPST = 0)
C                            -> Actuel Stationnaire : Exit sauf debut
C
C    2) Amont Stationnaire   -> Actuel Instationnaire : OK
C         (NPST > 0)                            (pertes Stats amont)
C                            -> Actuel Stationnaire : OK si IDSTNT et
C                                 NSTIST n'ont pas change, sinon Exit)
C
          IF (NPST.EQ.0 .AND. (ISTTIO.EQ.1 .AND. NSTIST.LE.IPLAS)) THEN
            WRITE(NFECRA,9140) FICMLS, IPLAS+1, NSTIST
            CALL CSEXIT (1)
          ENDIF
C
          IF ( NPST.GT.0 .AND.
     &        ( (ISTTIO.EQ.1 .AND. IPLAS.LE.NSTIST) .OR.
     &           ISTTIO.EQ.0)                              ) THEN
            WRITE(NFECRA,9141) FICMLS
          ENDIF
C
          IF (NPST.GT.0 .AND. (ISTTIO.EQ.1 .AND. IPLAS.GE.NSTIST)) THEN
           IF (  JDSTNT.NE.IDSTNT .OR.
     &           MSTIST.NE.NSTIST      ) THEN
              WRITE(NFECRA,9142) FICMLS
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
          IF ( NBCLST .NE. NCLSTO ) THEN
            WRITE(NFECRA,9143) FICMLS
            CALL CSEXIT (1)
          ENDIF
C
C --> Stats supplementaires utilisateurs
C
          RUBRIQ = 'nombre_statistiques_utilisateur'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                MVLSTS,IERROR)
C
          IF (NVLSTS.LT.MVLSTS) THEN
            WRITE(NFECRA,9150) FICMLS, MVLSTS, NVLSTS, NVLSTS, NVLSTS
          ENDIF
C
C  --> Lecture des Statistiques volumiques. Pas de traitement d'erreurs,
C        on suppose qu'elles sont dues a un changement de physique.
C
          ITYSUP = 1
          IRTYP  = 2
          NBVAL  = 1
C
          DO IPAS  = 0,NBCLST
            DO IVL = 1,NVLSTA
              IVAR  = IPAS*NVLSTA +IVL
              IF (IPAS.GT.0) THEN
                WRITE(CAR4,'(I4.4)') IPAS
                RUBRIQ = 'moy_stat_vol_groupe_'//CAR4//'_'//NOMLAG(IVAR)
              ELSE
                RUBRIQ = 'moy_stat_vol_'//NOMLAG(IVAR)
              ENDIF
              CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    STATIS(1,IVAR),IERROR)
            ENDDO
C
            DO IVL = 1,NVLSTA-1
              IF (IPAS.GT.0) THEN
                WRITE(CAR4,'(I4.4)') IPAS
                RUBRIQ = 'var_stat_vol_groupe_'//CAR4//'_'//NOMLAV(IVAR)
              ELSE
                RUBRIQ = 'var_stat_vol_'//NOMLAV(IVAR)
              ENDIF
              CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    STATIV(1,IVAR),IERROR)
            ENDDO
          ENDDO
C
        ENDIF

C
 9991   CONTINUE
C
        IF (IENSI3.EQ.1 .AND. NVISBR.GT.0 .AND. NFABOK.NE.0) THEN
C
          ITYSUP = 0
          NBVAL  = 1
C
          RUBRIQ = 'iteration_debut_stats_frontieres_stationnaires'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                MSTBOR,IERROR)
C
C  ---> S'il y a une erreur, on suppose que c'est parce que le fichier
C         suite ne contient pas d'infos sur les stats aux frontieres.
C         Dans ce cas, si on est en instationnaire on se dit que c'est
C         pas grave, on saute l'etape et on continue. Par contre si on
C         est dans une configuration de calcul de stats aux frontieres
C         en stationnaire on stoppe.
C
          IF(IERROR.NE.0) THEN
            IF ( ISTTIO.EQ.0 .OR.
     &          (ISTTIO.EQ.1 .AND. IPLAS.LT.NSTBOR) ) THEN
              WRITE(NFECRA,9210) FICMLS, ISTTIO, NSTBOR, IPLAS+1
              GOTO 9992
            ELSE
              WRITE(NFECRA,9220) FICMLS, ISTTIO, NSTBOR, IPLAS+1
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
C --> A partir d'ici on considere que le fichier suite contient
C       des stats volumiques
C
          IF ( JSTTIO.NE.ISTTIO .OR.
     &         MSTBOR.NE.NSTBOR     ) THEN
            WRITE (NFECRA,9230) FICMLS,
     &                          JSTTIO, MSTBOR,
     &                          ISTTIO, NSTBOR
          ENDIF
C
C  --> Lecture de l'avancement du calcul stats aux frontieres
C
          RUBRIQ = 'nombre_iterations_stats_frontieres'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPSTFT,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9060) FICMLS,
     &  'nombre_iterations_stats_frontieres                          ',
     &      'NPSTFT',NPSTFT
          ENDIF
C
          RUBRIQ = 'nombre_iterations_stats_frontieres_stationnaires'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPSTF,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9060) FICMLS,
     &  'nombre_iterations_stats_frontieres_stationnaires            ',
     &      'NPSTF',NPSTF
          ENDIF
C
          RUBRIQ = 'temps_stats_frontieres_stationnaires'
          IRTYP  = 2
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TSTATP,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9060) FICMLS,
     &  'temps_stats_frontieres_stationnaires                        ',
     &      'TSTATP',TSTATP
          ENDIF
C
C  --> Verif de coherence de l'avancement du calcul avec les
C       indicateurs de calcul de la suite actuelle :
C
C
          IF (NPSTF.EQ.0 .AND. (ISTTIO.EQ.1 .AND. NSTBOR.LE.IPLAS)) THEN
            WRITE(NFECRA,9240) FICMLS, IPLAS+1, NSTBOR
            CALL CSEXIT (1)
          ENDIF
C
          IF ( NPSTF.GT.0 .AND.
     &        ( (ISTTIO.EQ.1 .AND. IPLAS.LE.NSTBOR) .OR.
     &           ISTTIO.EQ.0)                             ) THEN
            WRITE(NFECRA,9241) FICMLS
          ENDIF
C
          IF (NPSTF.GT.0 .AND. (ISTTIO.EQ.1 .AND. IPLAS.GE.NSTBOR)) THEN
           IF (MSTBOR.NE.NSTBOR) THEN
              WRITE(NFECRA,9242) FICMLS
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
C --> Stats supplementaires utilisateurs
C
          RUBRIQ = 'nombre_stats_frontieres_utilisateur'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                MUSBOR,IERROR)
C
          IF (NUSBOR.LT.MUSBOR) THEN
            WRITE(NFECRA,9250) FICMLS, MUSBOR, NUSBOR, NUSBOR, NUSBOR
          ENDIF
C
C  --> Lecture des stats aux frontieres. Pas de traitement d'erreurs,
C        on suppose qu'elles sont dues a un changement de physique.
C
          ITYSUP = 3
          IRTYP  = 2
          NBVAL  = 1
C
          DO IVAR = 1,NVISBR
            RUBRIQ = 'stat_bord_'//NOMBRD(IVAR)
            CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PARBOR(1,IVAR),IERROR)
          ENDDO
C
        ENDIF
C
 9992   CONTINUE
C
        IF (IILAGR.EQ.2) THEN
C
          ITYSUP = 0
          NBVAL  = 1
C
          RUBRIQ = 'iteration_debut_termes_sources_stationnaires'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                MSTITS,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9020) FICMLS,
     &  'iteration_debut_termes_sources_stationnaires                ',
     &      'NSTITS',MSTITS
          ENDIF
C
C  ---> S'il y a une erreur, on suppose que c'est parce que le fichier
C         suite ne contient pas d'infos sur les TS de couplage retour.
C         Dans ce cas, si on est en instationnaire on se dit que c'est
C         pas grave, on saute l'etape et on continue. Par contre si on
C         est dans une configuration de calcul de stats aux frontieres
C         en stationnaire on stoppe.
C
          IF(IERROR.NE.0) THEN
            IF ( ISTTIO.EQ.0 .OR.
     &          (ISTTIO.EQ.1 .AND. IPLAS.LT.NSTITS) ) THEN
              WRITE(NFECRA,9310) FICMLS, ISTTIO, NSTITS, IPLAS+1
              GOTO 9993
            ELSE
              WRITE(NFECRA,9320) FICMLS, ISTTIO, NSTITS, IPLAS+1
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
C --> A partir d'ici on considere que le fichier suite contient
C       des stats volumiques
C
          RUBRIQ = 'modele_turbulence_termes_sources'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                JTURB,IERROR)
C
          JTYTUR = JTURB/10
C
          IF ( JSTTIO.NE.ISTTIO .OR.
     &         MSTITS.NE.NSTITS     ) THEN
            IF (JTYTUR.EQ.2) CAR8 = 'k-eps'
            IF (JTYTUR.EQ.3) CAR8 = 'Rij-eps'
            IF (JTURB.EQ.50) CAR8 = 'v2f'
            IF (JTURB.EQ.60) CAR8 = 'k-omega'
            IF (ITYTUR(IPHAS).EQ.2) KAR8 = 'k-eps'
            IF (ITYTUR(IPHAS).EQ.3) KAR8 = 'Rij-eps'
            IF (ITURB(IPHAS).EQ.50) KAR8 = 'v2f'
            IF (ITURB(IPHAS).EQ.60) KAR8 = 'k-omega'
            WRITE (NFECRA,9330) FICMLS,
     &                          JSTTIO, MSTITS, CAR8,
     &                          ISTTIO, NSTITS, KAR8
          ENDIF
C
C
C  --> Lecture de l'avancement du couplage retour
C
          RUBRIQ = 'nombre_iterations_termes_sources_stationnaires'
          IRTYP  = 1
          CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NPTS,IERROR)
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9060) FICMLS,
     &  'nombre_iterations_termes_sources_stationnaires              ',
     &      'NPTS',NPTS
          ENDIF
C
C  --> Verif de coherence de l'avancement du calcul avec les
C       indicateurs de calcul de la suite actuelle :
C
          IF (NPTS.EQ.0 .AND. (ISTTIO.EQ.1 .AND. NSTITS.LE.IPLAS)) THEN
            WRITE(NFECRA,9340) FICMLS, IPLAS+1, NSTITS
            CALL CSEXIT (1)
          ENDIF
C
          IF ( NPTS.GT.0 .AND.
     &        ( (ISTTIO.EQ.1 .AND. IPLAS.LE.NSTITS) .OR.
     &           ISTTIO.EQ.0)                             ) THEN
            WRITE(NFECRA,9341) FICMLS
          ENDIF
C
          IF (NPTS.GT.0 .AND. (ISTTIO.EQ.1 .AND. IPLAS.GE.NSTITS)) THEN
           IF (MSTITS.NE.NSTITS) THEN
              WRITE(NFECRA,9342) FICMLS
              CALL CSEXIT (1)
            ENDIF
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
          DO IVAR = 1,NTERSL
            RUBRIQ = NOMTSL(IVAR)
            CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  TSLAGR(1,IVAR),IERROR)
          ENDDO
C
C
C  Dans le cas specifique de la combustion de grains de charbon
C  avec un couplage retour sur une combustion gaz en phase porteuse
C
C      --> A verifier l'utilite de cette lecture pour une suite...
C
          IF (IPPMOD(ICPL3C).GE.0) THEN
            DO IVAR = 1, NSALPP
              ICHA = NSALTO-NSALPP+IVAR
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              WRITE(CAR4,'(I4.4)') IVAR
              RUBRIQ = 'scalaires_physiques_pariculieres_charbon'//CAR4
              CALL LECSUI(IMPMLS,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPPROC(ICHA)),IERROR)
            ENDDO
          ENDIF
C
        ENDIF
C
 9993   CONTINUE
C
        WRITE(NFECRA,6011)

C  ---> Fermeture du fichier suite
C
        CALL CLSSUI(IMPMLS,IERROR)
C
        IF(IERROR.NE.0) WRITE(NFECRA,9090) FICMLS, FICMLS
C
C ---> En cas d'erreur, on continue quand meme
C
        WRITE(NFECRA,7099)
C
      ENDIF
C
      WRITE(NFECRA,2000)
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
C
 2000 FORMAT(
     &'                                                             ',/,
     &'-------------------------------------------------------------',/)
C
 6000 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/,
     &           3X,'   -------------------------------------     ',/,
     &           3X,' Lecture d''un fichier suite                 ',/,
     &           3X,'   sur les variables liees aux particules    '  )
 6010 FORMAT (   3X,'   Debut de la lecture                       '  )
 6011 FORMAT (   3X,'   Fin   de la lecture                       '  )
 6099 FORMAT (   3X,' Fin de la lecture du fichier suite          ',/,
     &           3X,'   sur les variables liees aux particules    ',/)
C
 7000 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/,
     &           3X,'   -------------------------------------     ',/,
     &           3X,' Lecture d''un fichier suite                 ',/,
     &           3X,'   sur les statistiques et TS couplage retour'  )
 7010 FORMAT (   3X,'   Debut de la lecture                       '  )
 7099 FORMAT (   3X,' Fin de la lecture du fichier suite          ',/,
     &           3X,'   sur les statistiques et TS couplage retour'  )
C
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE D''UN FICHIER SUITE      ',/,
     &'@    *********                                    LAGRANGIEN ',/,
     &'@                                                            ',/,
     &'@      ERREUR A L''OUVERTURE DU FICHIER SUITE ',A13           ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier l''existence et le nom (',A13      ,') du      ',/,
     &'@        fichier suite dans le repertoire de travail.        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      TYPE DE FICHIER INCORRECT                             ',/,
     &'@                                                            ',/,
     &'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,
     &'@      suite Lagrangien.                                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        a un fichier suite Lagrangien.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de cellules a ete modifie                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Le nombre de faces ',A8  ,' a ete modifie.              ',/,
     &'@                                                            ',/,
     &'@    Le calcul peut etre execute mais les donnees            ',/,
     &'@      sur les faces ',A8  ,' ne seront pas relues           ',/,
     &'@      dans le fichier suite.                                ',/,
     &'@    Elles seront initialisees par des valeurs par defaut.   ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite ',A13      ,' utilise     ',/,
     &'@        correspond bien a un fichier suite Lagrangien,      ',/,
     &'@        et qu''il n''a pas ete endommage.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DES VARIABLES PARTICULAIRES       ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite ',A13      ,' utilise     ',/,
     &'@        correspond bien a un fichier suite Lagrangien,      ',/,
     &'@        et qu''il n''a pas ete endommage.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de particules contenu dans le fichier suite   ',/,
     &'@      est superieur au nombre de particule maximal autorise ',/,
     &'@      NBPART = ',I10                                         ,/,
     &'@      NBPMAX = ',I10                                         ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Augmenter la valeur de NBPMAX dans USLAG1.              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9060 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le mot cle est initialise avec sa valeur par defaut     ',/,
     &'@      ou celle donnee dans le sous-programme USLAG1 :       ',/,
     &'@        ',A10   ,'  = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9061 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le mot cle est initialise avec sa valeur par defaut     ',/,
     &'@      ou celle donnee dans le sous-programme USLAG1 :       ',/,
     &'@        ',A10   ,'  = ',E14.5                                ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9062 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9065 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    L''INDICATEUR DE CALCUL DES STATISTIQUES VOLUMIQUES     ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGLEC).                   ',/,
     &'@                                                            ',/,
     &'@    LORSQU''IL N''Y A PAS DE SUITE DE CALCUL SUR LES        ',/,
     &'@    STATISTIQUES VOLUMIQUES (ISUIST = ',  I3, '),'           ,/,
     &'@    IDSTNT DOIT ETRE UN ENTIER SUPERIEUR AU NUMERO          ',/,
     &'@       DE L''ITERATION LAGRANGIENNE DE REDEMARRAGE ', I10    ,/,
     &'@       IL VAUT ICI IDSTNT = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDSTNT dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9066 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    L''INDICATEUR DE CALCUL STATIONNAIRES DES STATISTIQUES  ',/,
     &'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGLEC).    ',/,
     &'@                                                            ',/,
     &'@    LORSQU''IL N''Y A PAS DE SUITE DE CALCUL SUR LES        ',/,
     &'@    STATISTIQUES AUX FRONTIERES (ISUIST = ',  I3, '),'       ,/,
     &'@    NSTBOR DOIT ETRE UN ENTIER SUPERIEUR AU NUMERO          ',/,
     &'@       DE L''ITERATION LAGRANGIENNE DE REDEMARRAGE ', I10    ,/,
     &'@       IL VAUT ICI NSTBOR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NSTBOR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9070 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant la physique associee         ',/,
     &'@      aux particules sont modifies :                        ',/,
     &'@                                                            ',/,
     &'@              IPHYLA    ITPVAR    IDPVAR    IMPVAR          ',/,
     &'@  AMONT : ',4I10                                             ,/,
     &'@  ACTUEL: ',4I10                                             ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9071 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Aucune selection de physique associee aux particules    ',/,
     &'@      n''est active. Les donnees amont sont perdues.        ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9072 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Une equation sur la temperature des particules est      ',/,
     &'@      enclenchee en cours de calcul.                        ',/,
     &'@    Initialisation par defaut :                             ',/,
     &'@       Temperature TPART = ', E14.5                          ,/,
     &'@       Chaleur massique CPPART = ', E14.5                    ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9073 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    L''indicateur d''un calcul Lagrangien de grains         ',/,
     &'@      de charbon est enclenche (IPHYLA = 2).                ',/,
     &'@    Le fichier suite ',A13      ,' ne correspond pas        ',/,
     &'@      a un calcul Lagrangien de grains de charbon.          ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier IPHYLA dans le sous-programme USLAG1.          ',/,
     &'@    Verifier le fichier suite utilise.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9074 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le fichier suite ',A13      ,' correspond               ',/,
     &'@      a un calcul Lagrangien de grains de charbon.          ',/,
     &'@    L''indicateur de physique actuel associee aux particules',/,
     &'@      a une valeur non permise dans le cadre d''une suite   ',/,
     &'@      d''un calcul Lagrangien de grains de charbon.         ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier IPHYLA dans le sous-programme USLAG1.          ',/,
     &'@    Verifier le fichier suite utilise.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9080 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    L''indicateur du  nombre de variables supplementaires   ',/,
     &'@      utilisateur est modifie, ou n''a pas pu etre relu.    ',/,
     &'@                                                            ',/,
     &'@              NVLS                                          ',/,
     &'@    AMONT : ',I10   ,'      ACTUEL : ',I10                   ,/,
     &'@                                                            ',/,
     &'@    Si ACTUEL > AMONT, on initialise les ',I10   ,' 1eres   ',/,
     &'@      variables supplementaires actuelles avec celles       ',/,
     &'@      du fichier suite, les autres sont initialisees a zero.',/,
     &'@                                                            ',/,
     &'@    Si ACTUEL < AMONT, on initialise les ',I10   ,' 1eres   ',/,
     &'@      variables supplementaires actuelles avec les 1eres    ',/,
     &'@      du fichier suite, le reste des variables du fichier   ',/,
     &'@      suite sont perdues.                                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9090 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA FERMETURE DU FICHIER SUITE             ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom : ',A13                   ,/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9110 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      LES STATISTIQUES VOLUMIQUES DU CALCUL AMONT           ',/,
     &'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,
     &'@        DU FICHIER SUITE                                    ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des                       ',/,
     &'@      statistiques volumiques sont positionnes en mode      ',/,
     &'@      instationnaire ou en debut de calcul stationnaire :   ',/,
     &'@                                                            ',/,
     &'@          ISTTIO    IDSTNT    NSTIST    Iter de redemarrage ',/,
     &'@      ',4I10                                                 ,/,
     &'@                                                            ',/,
     &'@    Elles seront initialisees par des valeurs par defaut.   ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9120 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      LES STATISTIQUES VOLUMIQUES DU CALCUL AMONT           ',/,
     &'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,
     &'@        DU FICHIER SUITE                                    ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des                       ',/,
     &'@      statistiques volumiques sont positionnes              ',/,
     &'@      en mode stationnaire :                                ',/,
     &'@                                                            ',/,
     &'@          ISTTIO    IDSTNT    NSTIST    Iter de redemarrage ',/,
     &'@      ',4I10                                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9130 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des                       ',/,
     &'@      statistiques volumiques sont modifies :               ',/,
     &'@                                                            ',/,
     &'@              ISTTIO    IDSTNT    NSTIST                    ',/,
     &'@  AMONT : ',3I10                                             ,/,
     &'@  ACTUEL: ',3I10                                             ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9140 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL ACTUEL DES STATISTIQUES VOLUMIQUES            ',/,
     &'@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER        ',/,
     &'@      SUITE CONTIENT DES STATISTIQUES INSTATIONNAIRES.      ',/,
     &'@                                                            ',/,
     &'@    NSTIST devrait etre un entier superieur ou egal         ',/,
     &'@      a l''iteration Lagrangienne absolue de redemarrage    ',/,
     &'@      du calcul (iteration : ',I10   ,')                    ',/,
     &'@                                                            ',/,
     &'@      Il vaut ici NSTIST = ',I10                             ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9141 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL ACTUEL DES STATISTIQUES VOLUMIQUES            ',/,
     &'@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER      ',/,
     &'@      SUITE CONTIENT DES STATISTIQUES STATIONNAIRES.        ',/,
     &'@                                                            ',/,
     &'@    Les statistiques volumiques stationnaires amont         ',/,
     &'@      seront remises a zero.                                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9142 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE                 ',/,
     &'@      STATISTIQUES VOLUMIQUES EN MODE STATIONNAIRE          ',/,
     &'@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES    ',/,
     &'@      ON ETE MODIFIEES.                                     ',/,
     &'@                                                            ',/,
     &'@    Pour eviter les incoherences dans le calcul             ',/,
     &'@      IDSTNT et NSTIST ne devraient pas etre modifies entre ',/,
     &'@      deux calculs de statistiques volumiques stationnaires.',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9143 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE                 ',/,
     &'@      STATISTIQUES VOLUMIQUES EN MODE STATIONNAIRE          ',/,
     &'@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES    ',/,
     &'@      ON ETE MODIFIEES.                                     ',/,
     &'@                                                            ',/,
     &'@    Pour eviter les incoherences dans le calcul             ',/,
     &'@      NBCLST ne devrait pas etre modifies entre             ',/,
     &'@      deux calculs de statistiques volumiques.              ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9150 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    L''indicateur du  nombre de statistiques volumiques     ',/,
     &'@      supplementaires utilisateur est modifie,              ',/,
     &'@      ou n''a pas pu etre relu.                             ',/,
     &'@                                                            ',/,
     &'@              NVLSTS                                        ',/,
     &'@    AMONT : ',I10   ,'      ACTUEL : ',I10                   ,/,
     &'@                                                            ',/,
     &'@    Si ACTUEL > AMONT, on initialise les ',I10   ,' 1eres   ',/,
     &'@      statistiques supplementaires actuelles avec celles    ',/,
     &'@      du fichier suite, les autres sont initialisees a zero.',/,
     &'@                                                            ',/,
     &'@    Si ACTUEL < AMONT, on initialise les ',I10   ,' 1eres   ',/,
     &'@      statistiques supplementaires actuelles avec les 1eres ',/,
     &'@      du fichier suite, le reste des statistiques du fichier',/,
     &'@      suite sont perdues.                                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      LES STATISTIQUES AUX FRONTIERES DU CALCUL AMONT       ',/,
     &'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,
     &'@        DU FICHIER SUITE                                    ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des                       ',/,
     &'@      statistiques aux frontieres sont positionnes en mode  ',/,
     &'@      instationnaire ou en debut de calcul stationnaire :   ',/,
     &'@                                                            ',/,
     &'@          ISTTIO    NSTBOR    Iter de redemarrage           ',/,
     &'@      ',3I10                                                 ,/,
     &'@                                                            ',/,
     &'@    Elles seront initialisees par des valeurs par defaut.   ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9220 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      LES STATISTIQUES AUX FRONTIERES DU CALCUL AMONT       ',/,
     &'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,
     &'@        DU FICHIER SUITE                                    ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des                       ',/,
     &'@      statistiques aux frontieres sont positionnes          ',/,
     &'@      en mode stationnaire :                                ',/,
     &'@                                                            ',/,
     &'@          ISTTIO    NSTBOR    Iter de redemarrage           ',/,
     &'@      ',3I10                                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9230 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des                       ',/,
     &'@      statistiques aux frontieres sont modifies :           ',/,
     &'@                                                            ',/,
     &'@                        ISTTIO    NSTBOR                    ',/,
     &'@            AMONT : ',2I10                                   ,/,
     &'@            ACTUEL: ',2I10                                   ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9240 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL ACTUEL DES STATISTIQUES AUX FRONTIERES        ',/,
     &'@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER        ',/,
     &'@      SUITE CONTIENT DES STATISTIQUES INSTATIONNAIRES.      ',/,
     &'@                                                            ',/,
     &'@    NSTBOR devrait etre un entier superieur ou egal         ',/,
     &'@      a l''iteration Lagrangienne absolue de redemarrage    ',/,
     &'@      du calcul (iteration : ',I10   ,')                    ',/,
     &'@                                                            ',/,
     &'@      Il vaut ici NSTBOR = ',I10                             ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9241 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL ACTUEL DES STATISTIQUES AUX FRONTIERES        ',/,
     &'@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER      ',/,
     &'@      SUITE CONTIENT DES STATISTIQUES STATIONNAIRES.        ',/,
     &'@                                                            ',/,
     &'@    Les statistiques aux frontieres stationnaires amont     ',/,
     &'@      seront remises a zero.                                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9242 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE                 ',/,
     &'@      STATISTIQUES AUX FRONTIERES EN MODE STATIONNAIRE      ',/,
     &'@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES    ',/,
     &'@      ON ETE MODIFIEES.                                     ',/,
     &'@                                                            ',/,
     &'@    Pour eviter les incoherences dans le calcul             ',/,
     &'@      NSTBOR ne devrait pas etre modifie entre deux calculs ',/,
     &'@      de statistiques aux frontieres stationnaires.         ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9250 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    L''indicateur du  nombre de statistiques aux frontieres ',/,
     &'@      supplementaires utilisateur est modifie,              ',/,
     &'@      ou n''a pas pu etre relu.                             ',/,
     &'@                                                            ',/,
     &'@              NUSBOR                                        ',/,
     &'@    AMONT : ',I10   ,'      ACTUEL : ',I10                   ,/,
     &'@                                                            ',/,
     &'@    Si ACTUEL > AMONT, on initialise les ',I10   ,' 1eres   ',/,
     &'@      statistiques supplementaires actuelles avec celles    ',/,
     &'@      du fichier suite, les autres sont initialisees a zero.',/,
     &'@                                                            ',/,
     &'@    Si ACTUEL < AMONT, on initialise les ',I10   ,' 1eres   ',/,
     &'@      statistiques supplementaires actuelles avec les 1eres ',/,
     &'@      du fichier suite, le reste des statistiques du fichier',/,
     &'@      suite sont perdues.                                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9310 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      LES TERMES SOURCES DE COUPLAGE RETOUR DU CALCUL AMONT ',/,
     &'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTS          ',/,
     &'@        DU FICHIER SUITE                                    ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des termes sources        ',/,
     &'@      de couplage retour sont positionnes en mode           ',/,
     &'@      instationnaire ou en debut de calcul stationnaire :   ',/,
     &'@                                                            ',/,
     &'@          ISTTIO    NSTITS    Iter de redemarrage           ',/,
     &'@      ',3I10                                                 ,/,
     &'@                                                            ',/,
     &'@    Ils seront initialisees par des valeurs par defaut.     ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@                                                            ',/,
     &'@      LES TERMES SOURCES DE COUPLAGE RETOUR DU CALCUL AMONT ',/,
     &'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTS          ',/,
     &'@        DU FICHIER SUITE                                    ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des termes sources        ',/,
     &'@      de couplage retour sont positionnes                   ',/,
     &'@      en mode stationnaire :                                ',/,
     &'@                                                            ',/,
     &'@          ISTTIO    NSTITS    Iter de redemarrage           ',/,
     &'@      ',3I10                                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9330 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Les indicateurs concernant le calcul                    ',/,
     &'@      instationnaire/stationnaire des termes sources        ',/,
     &'@      de couplage retour sont modifies :                    ',/,
     &'@                                                            ',/,
     &'@                   ISTTIO    NSTITS    Turbulence           ',/,
     &'@       AMONT : ',2I10,A13                                    ,/,
     &'@       ACTUEL: ',2I10,A13                                    ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1 et USINI1                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9340 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL ACTUEL DES TERMES SOURCES DE COUPLAGE RETOUR  ',/,
     &'@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER        ',/,
     &'@      SUITE CONTIENT DES TERMES SOURCES INSTATIONNAIRES.    ',/,
     &'@                                                            ',/,
     &'@    NSTITS devrait etre un entier superieur ou egal         ',/,
     &'@      a l''iteration Lagrangienne absolue de redemarrage    ',/,
     &'@      du calcul (iteration : ',I10   ,')                    ',/,
     &'@                                                            ',/,
     &'@      Il vaut ici NSTITS = ',I10                             ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9341 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL ACTUEL DES TERMES SOURCES DE COUPLAGE RETOUR  ',/,
     &'@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER      ',/,
     &'@      SUITE CONTIENT DES TERMES SOURCES STATIONNAIRES.      ',/,
     &'@                                                            ',/,
     &'@    Les termes sources de couplage retour stationnaires     ',/,
     &'@      amont seront remises a zero.                          ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@    Il est conseille de verifier ces indicateurs dans       ',/,
     &'@      le sous-programme USLAG1.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9342 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********     LAGRANGIEN ',A13                           ,/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DES TERMES         ',/,
     &'@      SOURCES DE COUPLAGE RETOUR EN MODE STATIONNAIRE       ',/,
     &'@      MAIS LES INDICATEURS DE CONTROLES DES TERMES SOURCES  ',/,
     &'@      ON ETE MODIFIEES.                                     ',/,
     &'@                                                            ',/,
     &'@    Pour eviter les incoherences dans le calcul             ',/,
     &'@      NSTITS ne devrait pas etre modifie entre deux calculs ',/,
     &'@      de termes sources de couplage retour stationnaires.   ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,
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
