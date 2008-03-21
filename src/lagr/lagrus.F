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
                        SUBROUTINE LAGRUS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , INDEP  ,
     &   IDEVEL , ITUSER , IA     ,
     &   ETTP   , ETTPA  , TEPA  , CROULE ,
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
CFONC      Roulette russe et clonage applique aux particules
CFONC      suivant un critere d'importance (CROULE)
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
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! INDEP(NBPMAX)! TE ! <-> ! NUMERO DE SA CELLULE DE DEPART       !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! ETTP         ! TR ! <-> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR ! <-> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! CROULE(NCELET! TR !  -> ! CRITERE D'IMPORTANCE                 !
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
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP) , INDEP(NBPMAX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION CROULE(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL    , IELD    , NCLO    , NPARS
      INTEGER          NPT    , N       , N1      , IVA    , NC
      DOUBLE PRECISION AUX(1) , COEFF  , PNEW     , DNPARS
C
C***********************************************************************
C
C=======================================================================
C 0. Initialisation
C=======================================================================
C
C     NPCLON : NOMBRE DE NOUVELLES PARTICULES PAR CLONNAGE
C
C     NPKILL : NOMBRE DE PARTICULES VICTIMES DE LA ROULETTE RUSSE
C
C     NPCSUP : NOMBRE DE PARTICULES QUI ON SUBIT LE CLONNAGE
C
      NPCLON = 0
      NPCSUP = 0
      NPKILL = 0
C
      DNPCLO = 0.D0
      DNPCSU = 0.D0
      DNPKIL = 0.D0
C
C=======================================================================
C 1. Clonage / Fusion (ou "Roulette Russe")
C=======================================================================
C
C
C Boucle sur les particules
C
      DO NPT = 1,NBPART
C
        IF (ITEPA(NPT,JISOR).NE.INDEP(NPT)) THEN
C
          IEL  = ITEPA(NPT,JISOR)
          IELD = INDEP(NPT)
C
C Rapport des fonction d'importance entre la cellule de depart
C et celle d'arrivee
C
          COEFF = CROULE(IEL) / CROULE(IELD)
C
          IF (COEFF.LT.1.D0) THEN
C
C---------------
C ROULETTE RUSSE
C---------------
C
            N1 = 1
            CALL ZUFALL(N1,AUX(1))
C
            IF (AUX(1).LT.COEFF) THEN
C
C La particule survit avec une probabilite COEFF
C
              TEPA(NPT,JRPOI) = TEPA(NPT,JRPOI)/COEFF
C
            ELSE
C
C La particule est supprimee avec une probabilite (1-COEFF)
C
              ITEPA(NPT,JISOR) = 0
              NPKILL = NPKILL + 1
              DNPKIL = DNPKIL + TEPA(NPT,JRPOI)
            ENDIF
C
          ELSE IF (COEFF.GT.1.D0) THEN
C
C--------
C CLONAGE
C--------
C
            N = INT(COEFF)
            N1 = 1
            CALL ZUFALL(N1,AUX(1))
C
            IF (AUX(1).LT.(COEFF-DBLE(N))) THEN
C
C Clonage en N+1 particules
C
              NCLO = N + 1
C
            ELSE
C
C Clonage en N particules
C
              NCLO = N
C
            ENDIF
C
            IF ((NBPART+NPCLON+NCLO+1).GT.NBPMAX) THEN
              WRITE(NFECRA,5000) NBPART, NPCLON+NCLO+1, NBPMAX
              GOTO 1000
            ENDIF
C
            NPCSUP = NPCSUP + 1
            DNPCSU = DNPCSU + TEPA(NPT,JRPOI)
            PNEW = TEPA(NPT,JRPOI) / DBLE(NCLO)
C
            DO NC = 1,NCLO
C
              NPCLON = NPCLON + 1
              DNPCLO = DNPCLO + PNEW
C
              DO IVA = 1,NVP
                ETTP(NBPART+NPCLON,IVA) = ETTP(NPT,IVA)
              ENDDO
C
              DO IVA = 1,NVP
                ETTPA(NBPART+NPCLON,IVA) = ETTPA(NPT,IVA)
              ENDDO
C
              DO IVA = 1,NVEP
                TEPA(NBPART+NPCLON,IVA) = TEPA(NPT,IVA)
              ENDDO
C
              TEPA(NBPART+NPCLON,JRPOI) = PNEW
C
              DO IVA = 1,NIVEP
                ITEPA(NBPART+NPCLON,IVA) = ITEPA(NPT,IVA)
              ENDDO
C
            ENDDO
C
C Modif de la particule elle meme
C
            ITEPA(NPT,JISOR) = 0
C
          ENDIF
        ENDIF
      ENDDO
C
 1000 CONTINUE
C
C Actualisation du nouveau nombre de particules
C
      NBPART = NBPART + NPCLON
      DNBPAR = DNBPAR + DNPCLO
C
C=======================================================================
C 2. On elimine les particules qui ont perdu à la Roulette Russe
C    et celles qui ont subit le clonage.
C=======================================================================
C
      CALL LAGELI
C     ===========
     & ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NPARS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  ,
     &   IDEVEL , ITUSER , IA    ,
     &   DNPARS ,
     &   ETTP   , ETTPA  , TEPA   ,
     &   RDEVEL , RTUSER , RA    )
C
      IF ( NPARS.NE.(NPKILL+NPCSUP) ) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT(1)
C       ===========
      ENDIF
C
C-------
C FORMAT
C-------
C
 5000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MODULE LAGRANGIEN                           ',/,
     &'@    *********                                               ',/,
     &'@                CLONAGE / FUSION DES PARTICULES             ',/,
     &'@                                                            ',/,
     &'@  Le nombre de nouvelles particules clonees conduit a un    ',/,
     &'@    nombre total de particules superieur au maximum prevu : ',/,
     &'@    Nombre de particules courant   : NBPART = ',I10          ,/,
     &'@    Nombre de particules clonnees  : NPCLON = ',I10          ,/,
     &'@    Nombre maximal de particules   : NBPMAX = ',I10          ,/,
     &'@                                                            ',/,
     &'@  On ne clone plus de particules por cette iteration.       ',/,
     &'@                                                            ',/,
     &'@  Verifier NBPMAX dans USLAG1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MODULE LAGRANGIEN                           ',/,
     &'@    *********                                               ',/,
     &'@                CLONAGE / FUSION DES PARTICULES             ',/,
     &'@                                                            ',/,
     &'@  La somme des particules detruites a la Roulette Russe     ',/,
     &'@    avec celles qui ont subit le clonage                    ',/,
     &'@    est different de celui des particules eliminees.        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier LAGRUS et LAGELI.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
