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
                        SUBROUTINE INITI2
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   JCELBR ,
     &   IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC FIN INITIALISATION DES COMMONS
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! JCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! IA(*)        ! TR !  -  ! TABLEAU DE TRAVAIL POUR LES ENTIERS  !
CARGU ! RA(*)        ! TR !  -  ! TABLEAU DE TRAVAIL POUR LES REELS    !
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
C-----------------------------------------------------------------------
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "dimens.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          JCELBR
      INTEGER          IA(*)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IPHAS
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2. DIMENSIONS DE dimens.h :    MLG ET NON ORTH
C=======================================================================
C
C---> COMMON GEOMET
C
C --- DEDUCTION DES AUTRES DIMENSIONS MLG ET NON ORTH
C
      NCELBR   = JCELBR
C
C=======================================================================
C 3. TABLEAUX DE cstphy.h
C=======================================================================
C
C---> COMMON TURBUL
C
C --- SI ALMAX < 0 , IL EST RECALCULE
C
      WRITE(NFECRA,1000)
C
      DO IPHAS = 1, NPHAS
        IF(ALMAX(IPHAS).LE.0.D0) THEN
          ALMAX(IPHAS) = VOLTOT**.333D0
          WRITE(NFECRA,1100) IPHAS,ALMAX(IPHAS)
          WRITE(NFECRA,1102)
          IF(ITYTUR(IPHAS).EQ.2.OR.ITYTUR(IPHAS).EQ.3
     &         .OR. ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60) THEN
            WRITE(NFECRA,1101)
          ENDIF
        ENDIF
      ENDDO
C
C
C
 1000 FORMAT(
     &'                                                             '  )
 1100 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )'  )
 1101 FORMAT(
     &'       ALMAX est la longueur utilisee pour initialiser       ',/,
     &'                                               la turbulence.'  )
 1102 FORMAT(
     &'       ALMAX est la racine cubique du volume du domaine.     ',/)
C
C=======================================================================
C 4. FIN
C=======================================================================
C
      RETURN
      END
c@z
