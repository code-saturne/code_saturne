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
                        SUBROUTINE LAGITG
C                       *****************
C     -----------------------------------------------------------------
     & ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   IVAR   ,
     &   ISORTI , IBORD  ,
     &   ETTP   , ETTPA  , TCARAC , PIP    , TSVAR  )
C     -----------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC              INTEGRATION DE L'EDS POUR LA VARIABLE IVAR
CFONC
CFONC
CFONC                  d V     V - PIP
CFONC          EDS :   --- = - -------
CFONC                  d t      TCARAC
CFONC
CFONC                Lorsqu'il y a eu interaction avec une face de bord,
CFONC                  l'integration est degeneree a l'ordre 1
CFONC                  (meme si on est a l'ordre 2 par ailleurs).
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
CARGU ! IVAR         ! E  !  -> ! NUMERO DE LA VARIABLE A INTEGRER     !
CARGU !              !    !     ! DANS LE TABLEAU ETTP                 !
CARGU ! ISORTI(NBPMAX! TE !  -> ! POUR CHAQUE PARTICULE :              !
CARGU !              !    !     !    * NUMERO DE SA CELLULE            !
CARGU !              !    !     !    * 0 SI SORTIE DU DOMAINE          !
CARGU ! IBORD        ! TE !  -> ! CONTIENT LE NUMERO DE LA             !
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
CARGU ! ETTP         ! TR ! <-  ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TCARAC(NBPMAX! TR !  -> ! TEMPS CARACTERISTIQUE ASSOCIE A L'EDS!
CARGU ! PIP(NBPMAX)  ! TR !  -> ! SECOND MEMBRE ASSOCIE A L'EDS        !
CARGU ! TSVAR        ! TR ! <-  ! PREDICTION 1ER SOUS-PAS POUR LA      !
CARGU ! (NBPMAX,NVP1)!    !     !   VARIABLE IVAR, UTILISE POUR LA     !
CARGU !              !    !     !   CORRECTION AU 2EME SOUS-PAS        !
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
      INTEGER          NBPMAX , NVP , NVP1 , NVEP , NIVEP
      INTEGER          IVAR
      INTEGER          ISORTI(NBPMAX) , IBORD(NBPMAX)
C
      DOUBLE PRECISION ETTP(NBPMAX,NVP) ,  ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TCARAC(NBPMAX) , PIP(NBPMAX)
      DOUBLE PRECISION TSVAR(NBPMAX,NVP1)
C
C VARIABLES LOCALES
C
      INTEGER          NPT
      DOUBLE PRECISION AUX1 , AUX2 , TER1 , TER2 , TER3
C
C***********************************************************************
C
C
        IF (NOR.EQ.1) THEN
C
          DO NPT = 1, NBPART
            IF (ISORTI(NPT).GT.0) THEN
C
              IF (TCARAC(NPT).LE.0.D0) THEN
                WRITE(NFECRA,2000) IVAR, TCARAC(NPT), NPT
                CALL CSEXIT (1)
              ENDIF
C
              AUX1 = DTP/TCARAC(NPT)
              AUX2 = EXP(-AUX1)
C
              TER1 = ETTPA(NPT,IVAR) * AUX2
              TER2 = PIP(NPT) * (1.D0-AUX2)
C
C        Pour le cas NORDRE= 1 ou s'il y a rebond,
C          le ETTP suivant est le resultat final
C
              ETTP(NPT,IVAR) = TER1 + TER2
C
C        Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2
C
              TER3 = ( -AUX2 + (1.D0-AUX2) / AUX1 ) * PIP(NPT)
              TSVAR(NPT,IVAR) = 0.5D0 * TER1 + TER3
            ENDIF
          ENDDO
C
        ELSE IF (NOR.EQ.2) THEN
C
          DO NPT = 1, NBPART
            IF (ISORTI(NPT).GT.0 .AND. IBORD(NPT).EQ.0)  THEN
C
              IF (TCARAC(NPT).LE.0.D0) THEN
                WRITE(NFECRA,2000) IVAR, TCARAC(NPT), NPT
                CALL CSEXIT (1)
              ENDIF
C
              AUX1 = DTP/TCARAC(NPT)
              AUX2 = EXP(-AUX1)
C
              TER1 = 0.5D0 * ETTPA(NPT,IVAR) * AUX2
              TER2 = PIP(NPT) * ( 1.D0 - (1.D0-AUX2) / AUX1 )
C
C        Pour le cas NORDRE= 2, le ETTP suivant est le resultat final
C
              ETTP(NPT,IVAR) = TSVAR(NPT,IVAR) + TER1 + TER2
            ENDIF
          ENDDO
C
        ELSE
          WRITE(NFECRA,1000) NOR
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C***********************************************************************
C
C-------
C FORMAT
C-------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    L''INDICATEUR SUR L''ORDRE D''INTEGRATION               ',/,
     &'@       DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES          ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGITG).                   ',/,
     &'@                                                            ',/,
     &'@    NORDRE DEVRAIT ETRE UN ENTIER EGAL A 1 OU 2             ',/,
     &'@       IL VAUT ICI NORDRE = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NORDRE dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    LE TEMPS CARACTERISTIQUE LIE A L''EQUATION              ',/,
     &'@      DIFFERENTIELLE STOCHASTIQUE DE LA VARIABLE            ',/,
     &'@      NUMERO ',I10 ,'UNE VALEUR NON PERMISE (LAGITG).       ',/,
     &'@                                                            ',/,
     &'@    TCARAC DEVRAIT ETRE UN ENTIER STRICTEMENT POSITIF       ',/,
     &'@       IL VAUT ICI TCARAC = ', E10.4                         ,/,
     &'@       POUR LA PARTICULE NUMERO ',I10                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
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
