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
                        SUBROUTINE PERGRA
C                       *****************
C     -------------------------------------------------------------
     & ( NPHJMX , NPHAS  ,
     &   JU     , JV     , JW     ,
     &   JTYTUR ,
     &   JR11   , JR22   , JR33   , JR12   , JR13   , JR23   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  Recuperation de certains COMMON necessaires a PERING
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NPHJMX       ! E  !  -> ! NOMBRE DE PHASES MAX                 !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! JU, JV, JW   ! TE ! <-  ! NUMERO DE VARIABLE POUR U, V, W      !
CARGU ! JTYTUR       ! TE ! <-  ! INDICATEUR MODELE DE TURBULENCE      !
CARGU ! JR11...JR23  ! TE ! <-  ! NUMERO DE VARIABLE POUR RIJ          !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NPHJMX, NPHAS
      INTEGER          JU(NPHJMX),JV(NPHJMX),JW(NPHJMX)
      INTEGER          JTYTUR(NPHJMX)
      INTEGER          JR11(NPHJMX),JR22(NPHJMX),JR33(NPHJMX)
      INTEGER          JR12(NPHJMX),JR13(NPHJMX),JR23(NPHJMX)
C
C VARIABLES LOCALES
C
      INTEGER          IPHAS
C
C***********************************************************************
C
C
C Recuperation des COMMON de "optcal"
C
      DO IPHAS = 1, NPHAS
        JTYTUR(IPHAS) = ITYTUR(IPHAS)
      ENDDO
C
C Recuperation des COMMON de "numvar.h"
C
      DO IPHAS = 1, NPHAS
        JU  (IPHAS) = IU  (IPHAS)
        JV  (IPHAS) = IV  (IPHAS)
        JW  (IPHAS) = IW  (IPHAS)
        IF(ITYTUR(IPHAS).EQ.3) THEN
          JR11(IPHAS) = IR11(IPHAS)
          JR22(IPHAS) = IR22(IPHAS)
          JR33(IPHAS) = IR33(IPHAS)
          JR12(IPHAS) = IR12(IPHAS)
          JR13(IPHAS) = IR13(IPHAS)
          JR23(IPHAS) = IR23(IPHAS)
        ELSE
          JR11(IPHAS) = 0
          JR22(IPHAS) = 0
          JR33(IPHAS) = 0
          JR12(IPHAS) = 0
          JR13(IPHAS) = 0
          JR23(IPHAS) = 0
        ENDIF
      ENDDO
C
C----
C FIN
C----
C
      RETURN
      END
c@z
