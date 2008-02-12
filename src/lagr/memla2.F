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
                        SUBROUTINE MEMLA2
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFABOR , NCELET , NFAC   ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   IINDEP , IIBORD , IETTPA , IAUXL  , IAUXL2 ,
     &   ITAUP  , IITLAG , IPIIL  ,
     &   IVAGAU , ITSUF  , ITSUP  , IBX    ,
     &   IGRADP , IGRADV , ICROUL ,
     &   ITEPCT , ITSFEX , ITSVAR ,
     &   ICPGD1 , ICPGD2 , ICPGHT ,
     &   IBRGAU , ITEBRU ,
     &   IW1    , IW2    , IW3    ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC  Reservation de la memoire pour les tableaux qui ne doivent pas
CFONC  etre conserves en dehors de la boucle en temps.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/IDBRA0! TR !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DES!
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS                    !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! IINDEP       ! E  ! <-  ! POINTEUR SUR INDEP (NUM CEL DE DEPART!
CARGU ! IIBORD       ! E  ! <-  ! POINTEUR SUR IBORD                   !
CARGU ! IETTPA       ! E  ! <-  ! POINTEUR SUR ETTPA                   !
CARGU ! IAUXL        ! E  ! <-  ! POINTEUR SUR AUXL                    !
CARGU ! IAUXL2       ! E  ! <-  ! POINTEUR SUR AUXL2                   !
CARGU ! ITAUP        ! E  ! <-  ! POINTEUR SUR TAUP                    !
CARGU ! IITLAG       ! E  ! <-  ! POINTEUR SUR TLAG                    !
CARGU ! IPIIL        ! E  ! <-  ! POINTEUR SUR PIIL                    !
CARGU ! IVAGAU       ! E  ! <-  ! POINTEUR SUR VAGAUS                  !
CARGU ! ITSUF        ! E  ! <-  ! POINTEUR SUR TSUF                    !
CARGU ! ITSUP        ! E  ! <-  ! POINTEUR SUR TSUP                    !
CARGU ! IBX          ! E  ! <-  ! POINTEUR SUR BX                      !
CARGU ! IGRADP       ! E  ! <-  ! POINTEUR SUR GRADPR                  !
CARGU ! IGRADV       ! E  ! <-  ! POINTEUR SUR GRADVF                  !
CARGU ! ICROUL       ! E  ! <-  ! POINTEUR SUR CROULE                  !
CARGU ! ITEPCT       ! E  ! <-  ! POINTEUR SUR TEMPCT                  !
CARGU ! ITSFEX       ! E  ! <-  ! POINTEUR SUR TSFEXT                  !
CARGU ! ITSVAR       ! E  ! <-  ! POINTEUR SUR TSVAR                   !
CARGU ! ICPGD1       ! E  ! <-  ! POINTEUR SUR CPGD1                   !
CARGU ! ICPGD2       ! E  ! <-  ! POINTEUR SUR CPGD2                   !
CARGU ! ICPGHT       ! E  ! <-  ! POINTEUR SUR CPGHT                   !
CARGU ! IW1          ! E  ! <-  ! POINTEUR SUR W1                      !
CARGU ! IW2          ! E  ! <-  ! POINTEUR SUR W2                      !
CARGU ! IW3          ! E  ! <-  ! POINTEUR SUR W3                      !
CARGU !              !    !     !                                      !
CARGU !              ! TR !     !  TABLEAUX IA/RA                      !
CARGU ! IFINIA       ! TR ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              ! TR !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! TR ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              ! TR !     !  DANS IA EN SORTIE                   !
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
      INTEGER          NFABOR , NCELET , NFAC
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          IINDEP , IIBORD , IETTPA , IAUXL , IAUXL2
      INTEGER          ITAUP  , IITLAG , IPIIL
      INTEGER          IVAGAU , ITSUF  , ITSUP  , IBX
      INTEGER          IGRADP , IGRADV , ICROUL
      INTEGER          ITEPCT , ITSFEX , ITSVAR
      INTEGER          ICPGD1 , ICPGD2 , ICPGHT
      INTEGER          IBRGAU , ITEBRU
      INTEGER          IW1     , IW2    , IW3
      INTEGER          IFINIA , IFINRA
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
C
C***********************************************************************
C
C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
C
      IIBORD =       IDEBIA
      IFINIA =       IIBORD + NBPMAX
C
      IINDEP = IFINIA
      IFINIA = IINDEP + NBPMAX
C
      IETTPA =       IDEBRA
      IAUXL  =       IETTPA + NBPMAX * NVP
      ITAUP  =       IAUXL  + NBPMAX * 3
      IITLAG =       ITAUP  + NBPMAX
      IPIIL  =       IITLAG + NBPMAX * 3
      IVAGAU =       IPIIL  + NBPMAX * 3
      ITSUF  =       IVAGAU + NBPMAX * NVGAUS
      ITSUP  =       ITSUF  + NBPMAX * 3
      IBX    =       ITSUP  + NBPMAX * 3
      ITSVAR =       IBX    + NBPMAX * 3 * 2
      IGRADP =       ITSVAR + NBPMAX * NVP1
      IW1    =       IGRADP + NCELET * 3
      IW2    =       IW1    + NCELET
      IW3    =       IW2    + NCELET
      IFINRA =       IW3    + NCELET
C
      IF ( (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) .OR.
     &      IPHYLA.EQ.2                         ) THEN
        ITEPCT = IFINRA
        IFINRA = ITEPCT + 2*NBPMAX
      ELSE
        ITEPCT = 1
      ENDIF
C
      IF (IILAGR.EQ.2) THEN
        ITSFEX = IFINRA
        IFINRA = ITSFEX + NBPMAX
      ELSE
        ITSFEX = 1
      ENDIF
C
      IF ( IILAGR.EQ.2 .AND. IPHYLA.EQ.2
     &                 .AND. LTSTHE .EQ.1 ) THEN
        ICPGD1 = IFINRA
        ICPGD2 = ICPGD1 + NBPMAX
        ICPGHT = ICPGD2 + NBPMAX
        IFINRA = ICPGHT + NBPMAX
      ELSE
        ICPGD1 = 1
        ICPGD2 = 1
        ICPGHT = 1
      ENDIF
C
      IF (MODCPL.GT.0) THEN
        IGRADV =     IFINRA
        IFINRA =     IGRADV + NCELET * 9
      ELSE
        IGRADV =  1
      ENDIF
C
      IF (IROULE.EQ.1) THEN
        ICROUL = IFINRA
        IFINRA = ICROUL + NCELET
      ELSE
        ICROUL = 1
      ENDIF
C
      IF ( LAMVBR .EQ. 1 ) THEN
        IBRGAU = IFINRA
        ITEBRU = IBRGAU + NBPMAX * NBRGAU
        IFINRA = ITEBRU + NBPMAX
      ELSE
        IBRGAU = 1
      ENDIF
C
      IF (NORDRE.EQ.2) THEN
        IAUXL2 = IFINRA
        IFINRA = IAUXL2 + NBPMAX*7
      ELSE
        IAUXL2 = 1
      ENDIF
C
C---> VERIFICATION
C
      CALL IASIZE('MEMLA2',IFINIA)
C     ==========
C
      CALL RASIZE('MEMLA2',IFINRA)
C     ==========
C
      RETURN
      END
c@z
