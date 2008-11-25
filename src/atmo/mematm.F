c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE MEMATM
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IFINIA , IFINRA , RA )
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   MEMORY MANAGEMENT FOR THE ATMOSPHERIC SPECIFIC PHYSICS MODULE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU !______________.____._____.______________________________________!
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
C     COMMON DATA
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "atincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0, IDBRA0, IFINIA, IFINRA
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IMODE
      DOUBLE PRECISION W1
C
C=======================================================================
C  1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C  2. MEMORY MANAGEMENT FOR THE METEO PROFILES
C=======================================================================
C
C  2.1 READING THE METEO FILE IN ORDER TO ESTIMATE THE REQUIRED SIZE OF
C  THE TABLES
C----------------------------------------------------------------------
      IF (IMETEO.GT.0) THEN
C
        IMODE = 0
C
C     Nb les arguments ne sont pas utilises quand IMODE=0
        CALL ATLECM
C       ===========
     & ( IMODE      ,
     &   W1 , W1 , W1 ,
     &   W1 , W1 , W1 ,
     &   W1 , W1 ,
     &   W1 , W1 ,
     &   W1 , W1 ,
     &   W1 , W1 , W1 )
C
        ENDIF

C  2.2 MEMORY MANAGEMENT :
C ------------------------
C
C --> Integer
C
      IFINIA = IDEBIA
C
C --> Real
C
      ITMMET = IDEBRA
      IZDMET = ITMMET + NBMETM
      IZTMET = IZDMET + NBMETD
      IUMET  = IZTMET + NBMETT
      IVMET  = IUMET  + NBMETD*NBMETM
      IWMET  = IVMET  + NBMETD*NBMETM
      IEKMET = IWMET  + NBMETD*NBMETM
      IEPMET = IEKMET + NBMETD*NBMETM
      ITTMET = IEPMET + NBMETD*NBMETM
      IQVMET = ITTMET + NBMETT*NBMETM
      IPMER  = IQVMET + NBMETT*NBMETM
      IXMET  = IPMER  + NBMETM
      IYMET  = IXMET  + NBMETM
      IRMET  = IYMET  + NBMETM
      ITPMET = IRMET  + NBMETT*NBMETM
      IPHMET = ITPMET + NBMETT*NBMETM
      IFINRA = IPHMET + NBMETT*NBMETM
C
C --> Verification
C
      CALL IASIZE('MEMATM',IFINIA)
C     ===========
C
      CALL RASIZE('MEMATM',IFINRA)
C     ===========
C
C----
C FIN
C----
C
      RETURN
      END
C

