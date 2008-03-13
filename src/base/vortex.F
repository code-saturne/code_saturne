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
                        SUBROUTINE VORTEX
C                       *****************
C     -------------------------------------------------------------
     & ( IVORCE , VISCO  , XYZ    ,
     &   YZCEL  , XU     , XV     , XW     ,
     &   YZVOR  , YZVORA , SIGNV  ,
     &   SIGMA  , GAMMA  , TEMPS  , TPSLIM )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ----------
c@foncb
CFONC
CFONC  GESTION DES ENTREES L.E.S. PAR LA METHODE DES VORTEX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IVORCE       ! TE ! <-> ! NUMERO DU VORTEX LE PLUS PROCHE D'UN !
CARGU !     (NVOMAX) !    !     ! VORTEX DONNE                         !
CARGU ! VISCO        ! TR !  -> ! VISCOSITE CINEMATIQUE SUR LES FACES  !
CARGU !(ICVMAX,NNENT)!    !     ! D'ENTREE                             !
CARGU ! XYZ(ICVMAX,3)!    !  -> ! COORDONNEES DES CELLULES D'ENTREE    !
CARGU !              !    !     ! DANS LE REFERENTIEL GLOBAL           !
CARGU ! YZCEL        ! TR !  -> ! COORDONNEES DES FACES D'ENTREE DANS  !
CARGU !   (NELVMX ,2)!    !     ! LE REFERENTIEL LOCAL                 !
CARGU ! XU(NELVMX)   ! TR ! <-> ! COMPOSANTE DE VITESSE PRINCIPALE     !
CARGU ! XV(NELVMX)   ! TR ! <-> ! COMPOSANTES DE VITESSE TRANSVERSES   !
CARGU ! XW(NELVMX)   ! TR ! <-> !                                      !
CARGU ! YZVOR        ! TR ! <-> ! COORDONNEES DU CENTRE DES VORTEX     !
CARGU !   (NVOMAX,2) !    !     !                                      !
CARGU ! YZVORA       ! TR ! <-> ! ANCIENNES COORDONNEES DU CENTRE      !
CARGU !   (NVOMAX,2) !    !     ! DES VORTEX                           !
CARGU ! SIGNV(NVOMAX)! TR ! <-> ! SENS DE ROTATION DES VORTEX          !
CARGU ! SIGMA        ! TR ! <-> ! TAILLE DES VORTEX                    !
CARGU !(NVOMAX,NNENT)!    !     !                                      !
CARGU ! GAMMA        ! TR ! <-> ! INTENSITE DES VORTEX                 !
CARGU !(NVOMAX,2,NNEN!    !     ! (DANS LES DEUX DIRECTIONS DU PLAN)   !
CARGU ! TEMPS        ! TR ! <-> ! TEMPS ECOULE DEPUIS LA CREATION      !
CARGU !     (NVOMAX) !    !     ! DU VORTEX                            !
CARGU ! TPSLIM       ! TR ! <-> ! DUREE DE VIE DU VORTEX               !
CARGU !(NVOMAX,NNENT)!    !     !                                      !
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
      INCLUDE "optcal.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IVORCE(NVOMAX,NNENT)
C
      DOUBLE PRECISION YZCEL(ICVMAX,2,NNENT)    , VISCO(ICVMAX,NNENT)
      DOUBLE PRECISION XYZ(ICVMAX,3,NNENT)      , XU(ICVMAX,NNENT)
      DOUBLE PRECISION XV(ICVMAX,NNENT)         , XW(ICVMAX,NNENT)
      DOUBLE PRECISION YZVOR(NVOMAX,2,NNENT)    , YZVORA(NVOMAX,2,NNENT)
      DOUBLE PRECISION SIGNV(NVOMAX,NNENT)
      DOUBLE PRECISION SIGMA(NVOMAX,NNENT)      , GAMMA(NVOMAX,2,NNENT)
      DOUBLE PRECISION TEMPS(NVOMAX,NNENT)      , TPSLIM(NVOMAX,NNENT)
C
C     VARIABLES LOCALES
C
      INTEGER          II, IENT
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C
C L'EQUATION DE LANGEVIN RESTE A TRAVAILLER POUR UN CAS 3D QUELCONQUE
C OU A MODIFIER. VU LE PEU D'IMPORTANCE QU'ELLE A SUR CE QUI SE PASSE
C EN AVAL DE L'ENTREE (L'IMPORTANT ETANT D'IMPOSER V' ET W'), ON NE VA
C PAS PLUS LOIN (ON ANNULE CES CONTRIBUTION POUR LE MOMENT POUR LES
C CAS 3 ET 4)
C
      IPASS = IPASS + 1
C
      DO IENT = 1, NNENT
C
        IF (IPASS.EQ.1)THEN
C
          CALL VORINI
C         ===========
     & ( ICVOR(IENT)     , NVORT(IENT)     ,
     &   IENT   , IVORCE(1,IENT)  ,
     &   XYZ(1,1,IENT)   , YZCEL(1,1,IENT) ,
     &   XU(1,IENT)      , XV(1,IENT)      , XW(1,IENT)      ,
     &   YZVOR(1,1,IENT) , SIGNV(1,IENT)   , TEMPS(1,IENT)   ,
     &   TPSLIM(1,IENT)  )
C
        ENDIF
C
C=======================================================================
C 2. DEPLACEMENT DU VORTEX
C=======================================================================
C
        CALL VORDEP
C       ===========
     & ( ICVOR(IENT)     , NVORT(IENT)     , IENT   , DTREF  ,
     &   IVORCE(1,IENT)  , YZCEL(1,1,IENT) ,
     &   XU(1,IENT)      , XV(1,IENT)      , XW(1,IENT)      ,
     &   YZVOR(1,1,IENT) , YZVORA(1,1,IENT), SIGNV(1,IENT)   ,
     &   TEMPS(1,IENT)   , TPSLIM(1,IENT)  )
C
C=======================================================================
C 3. CALCUL DE LA VITESSE
C=======================================================================
C
        CALL VORVIT
C       ===========
     & ( ICVOR(IENT)     , NVORT(IENT)     , IENT   ,
     &   IVORCE(1,IENT)  , VISCO(1,IENT)   ,
     &   YZCEL(1,1,IENT) , XU(1,IENT)      , XV(1,IENT)       ,
     &   XW(1,IENT)      , YZVOR(1,1,IENT) , SIGNV(1,IENT)    ,
     &   SIGMA(1,IENT)   , GAMMA(1,1,IENT) , TEMPS(1,IENT)    )
C
C=======================================================================
C 4. CALCUL DES FLUCTUATIONS DANS LE SENS DE L'ECOULEMENT
C=======================================================================
C
        CALL VORLGV
C       ===========
     & ( ICVOR(IENT)     , IENT   , DTREF  ,
     &   YZCEL(1,1,IENT) , XU(1,IENT)      ,
     &   XV(1,IENT)      , XW(1,IENT)      )
C
      ENDDO
C
C=======================================================================
C 5. ECRITURE DU FICHIER SUITE
C=======================================================================
C
C on ecrit a tous les pas de temps pour eviter
C les mauvaises surprises en cas de fin prematuree.
C Il ne faut pas mettre cette partie dans la boucle sur IENT
C car on accede deja a l'unite IMPMVO(=IMPDVO) dans VORINI.
C Seul le premier processeur ecrit (test avant l'appel à VORTEX)
C
      OPEN(UNIT=IMPVVO,FILE=FICVVO)
      REWIND(IMPVVO)
      DO IENT = 1, NNENT
        WRITE(IMPVVO,100) IENT
        WRITE(IMPVVO,100) NVORT(IENT)
        DO II = 1, NVORT(IENT)
          WRITE(IMPVVO,200) YZVOR(II,1,IENT),YZVOR(II,2,IENT),
     &         TEMPS(II,IENT), TPSLIM(II,IENT), SIGNV(II,IENT)
        ENDDO
      ENDDO
      CLOSE(IMPVVO)
C
C=======================================================================
C 6. FIN
C=======================================================================
C
C FORMATS
C
 100  FORMAT(I10)
 200  FORMAT(5E13.5)
C
      RETURN
C
      END
c@z
