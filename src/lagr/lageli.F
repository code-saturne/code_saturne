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
                        SUBROUTINE LAGELI
C                       *****************
C
C      -------------------------------------------------
     & ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NPARS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IDEVEL , ITUSER , IA     ,
     &   DNPARS ,
     &   ETTP   , ETTPA  , TEPA   , RDEVEL , RTUSER , RA )
C      -------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC    ELIMINATION DES PARTICULES QUI SONT SORTIES DU DOMAIME
CFONC      --> on gere la memoire pour eviter les places libres
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NPARS        ! E  ! <-  ! NOMBRE MAX DE PARTICULES SORTIES     !
CARGU !              !    !     !   ELIMINEES                          !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ITEPA        ! TE ! <-> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! DNPARS       ! E  ! <-  ! NOMBRE MAX DE PARTICULES SORTIES     !
CARGU !              !    !     !   ELIMINEES  (POIDS STAT INCLUS)     !
CARGU ! ETTP         ! TR ! <-> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR ! <-> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
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
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C**********************************************************************
C
C ARGUMENTS
C
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NPARS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION DNPARS
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C  VARIABLES LOCALES
C
      INTEGER NBP , NPT , I , IVAR
      DOUBLE PRECISION DNBP
C
C***********************************************************************
C
      NBP  = NBPART
      DNBP = DNBPAR
C
      NPARS  = 0
      DNPARS = 0.D0
C
      DO NPT = NBPART,1,-1
C
       IF (NBPART.LT.1) THEN
         WRITE(NFECRA,*) ' erreur lageli '
       ENDIF
C
        IF (ITEPA(NPT,JISOR).EQ.0) THEN
C
          NPARS  = NPARS  + 1
          DNPARS = DNPARS + TEPA(NBP,JRPOI)
C
C      ---> la particule est sortie du domaine
C
          IF (NPT.EQ.NBP) THEN
C
C        ---> c'est la derniere particule, on la supprime seulement
C
            NBP  = NBP - 1
            DNBP = DNBP - TEPA(NBP,JRPOI)
C
            DO I = 1,NLISTE
              IF ( LISTE(I).EQ.NPT ) THEN
                LISTE(I) = -1
              ENDIF
            ENDDO
C
          ELSE
C
C        ---> la particule NPT est supprime et on met a la place la
C             particule NBP
C
            DO IVAR = 1,NVP
              ETTP(NPT,IVAR) = ETTP(NBP,IVAR)
              ENDDO
C
            DO IVAR = 1,NVP
              ETTPA(NPT,IVAR) = ETTPA(NBP,IVAR)
            ENDDO
C
            DO IVAR = 1,NVEP
              TEPA(NPT,IVAR) = TEPA(NBP,IVAR)
            ENDDO
C
            DO IVAR = 1,NIVEP
              ITEPA(NPT,IVAR) = ITEPA(NBP,IVAR)
            ENDDO
C
            DO I = 1,NLISTE
              IF (LISTE(I).EQ.NPT) THEN
                LISTE(I) = -1
              ENDIF
            ENDDO
C
            DO I = 1,NLISTE
              IF (LISTE(I).EQ.NBP) THEN
                LISTE(I) = NPT
              ENDIF
            ENDDO
C
            NBP  = NBP - 1
            DNBP = DNBP - TEPA(NBP,JRPOI)
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C     ---> On met NBPART a la bonne valeur
C
      NBPART = NBP
      DNBPAR = DNBP
C
C***********************************************************************
C
C====
C FIN
C====
C
      END
c@z
