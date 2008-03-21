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
                        SUBROUTINE MEMLA1
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   LNDNOD ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   IIITEP , IICOCE , IITYCE ,
     &   IETTP  , IETTPA , IITEPA , ISTATC , ISTATV ,
     &   ITSLAG , ISTATF ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC          Reservation de la memoire pour les tableaux qui doivent
CFONC          etre conserves en dehors de la boucle en temps, notamment
CFONC          pour le post-processing.
CFONC
CFONC          Remarque : tous les tableaux ouverts ici doivent
CFONC                     recevoir une initialisation de valeurs
CFONC                     par defaut dans le sous-programme LAGLEC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA EN ENTREE                   !
CARGU ! IDBRA0       ! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS RA EN ENTREE                   !
CARGU ! NDIM         ! E  !  -> ! DIMENSION (3)                        !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! LNDNOD       ! E  !  -> ! DIM. CONNECTIVITE CELLULES->FACES    !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! IIITEP       ! E  ! <-  ! POINTEUR SUR ITEPA                   !
CARGU ! IICOCE       ! E  ! <-  ! POINTEUR SUR ICOCEL                  !
CARGU ! IITYCE       ! E  ! <-  ! POINTEUR SUR ITYCEL                  !
CARGU ! IETTP        ! E  ! <-  ! POINTEUR SUR ETTP                    !
CARGU ! IETTPA       ! E  ! <-  ! POINTEUR SUR ETTPA (simple init. ici)!
CARGU ! IITEPA       ! E  ! <-  ! POINTEUR SUR TEPA                    !
CARGU ! ISTATC       ! E  ! <-  ! POINTEUR SUR STATIS                  !
CARGU ! ITSLAG       ! E  ! <-  ! POINTEUR SUR TSLAGR                  !
CARGU ! ISTATF       ! E  ! <-  ! POINTEUR SUR PARBOR                  !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS RA EN SORTIE                   !
CARGU .______________.____._____.______________________________________.
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
      INTEGER          LNDNOD
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          IIITEP , IICOCE , IITYCE
      INTEGER          IETTP  , IETTPA
      INTEGER          IITEPA , ISTATC , ISTATV , ITSLAG , ISTATF
      INTEGER          IFINIA , IFINRA
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
C
C***********************************************************************
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
C
C     IETTPA est simplement initialise ici, pour eviter un arret avec les
C       options de compilation check bounds. La reservation memoire
C       correspondante est faite dans memla2.
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
      IF (IILAGR.EQ.0) THEN
C
        IIITEP =       IDEBIA
        IICOCE =       IIITEP
        IITYCE =       IICOCE
        IIFRLA =       IITYCE
        IFINIA =       IIFRLA
C
        IETTPA =       IDEBRA
        IETTP  =       IETTPA
        IITEPA =       IETTP
        ISTATC =       IITEPA
        ISTATV =       ISTATC
        ITSLAG =       ISTATV
        ISTATF =       ITSLAG
        IFINRA =       ISTATF
C
      ELSE
C
        IIITEP =       IDEBIA
        IICOCE =       IIITEP + NBPMAX * NIVEP
        IITYCE =       IICOCE + LNDNOD
        IIFRLA =       IITYCE + NCELET + 1
        IFINIA =       IIFRLA + NFABOR
C
        IETTPA =       IDEBRA
        IETTP  =       IETTPA
        IITEPA =       IETTP  + NBPMAX * NVP
        ISTATC =       IITEPA + NBPMAX * NVEP
        ISTATV =       ISTATC + NCELET * NVLSTA
     &                        + NCELET * NVLSTA * NBCLST
        ITSLAG =       ISTATV + NCELET * MAX((NVLSTA-1),0)
     &                        + NCELET * MAX((NVLSTA-1),0) * NBCLST
        ISTATF =       ITSLAG + NTERSL * NCELET
        IFINRA =       ISTATF + NFABOR * NVISBR
C
      ENDIF
C
C---> VERIFICATION
C
      CALL IASIZE('MEMLA1',IFINIA)
C     ===========
C
      CALL RASIZE('MEMLA1',IFINRA)
C     ===========
C
      RETURN
      END
c@z
