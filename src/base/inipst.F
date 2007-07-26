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
                        SUBROUTINE  INIPST
C                       ******************
C     -------------------------------------------------------------
     & ( IPSTVL , IPSTBO , IPSTSY ,
     &   IPSTMD , NTPST  , FMTPST , OPTPST )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  RENVOIE LES PARAMETRES CONTENUS DANS LES COMMONS ET UTILES
CFONC  A L'INITIALISATION DU POST-TRAITEMENT
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IPSTVL       ! E  ! <-  ! INDICATEUR POUR LE MAILLAGE VOLUMIQUE!
CARGU ! IPSTBO       ! E  ! <-  ! INDICATEUR POUR LE MAILLAGE DE PEAU  !
CARGU ! IPSTSY       ! E  ! <-  ! INDICATEUR POUR LES MAILLAGES DE     !
CARGU !              !    !     !  PEAU COUPLES AVEC SYRTHES           !
CARGU ! IPSTMD       ! E  ! <-  ! INDICATEUR DE MAILLAGE DEFORMABLE :  !
CARGU !              !    !     !  0 : PAS DE DEFORMATION ;            !
CARGU !              !    !     !  1 : DEFORMATION A TOPOLOGIE CSTE    !
CARGU !              !    !     !  2 : TOPOLOGIE VARIABLE              !
CARGU !              !    !     ! 10 : IDEM QUE 0 AVEC ECRITURE D'UN   !
CARGU !              !    !     !                      CHAMP DE DEPLACT!
CARGU !              !    !     ! 11 : IDEM QUE 1 AVEC ECRITURE D'UN   !
CARGU !              !    !     !                      CHAMP DE DEPLACT!
CARGU !              !    !     ! 12 : IDEM QUE 2 AVEC ECRITURE D'UN   !
CARGU !              !    !     !                      CHAMP DE DEPLACT!
CARGU ! NTPST        ! E  ! <-  ! FREQUENCE DES SORTIES POST-TRAITEMENT!
CARGU ! FMTPST       ! A  ! <-  ! NOM DU FORMAT DE POST-TRAITEMENT     !
CARGU ! FMTPST       ! A  ! <-  ! OPTIONS DU FORMAT DE POST-TRAITEMENT !
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
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IPSTVL , IPSTBO , IPSTSY , IPSTMD , NTPST
C
      CHARACTER        FMTPST(32)
      CHARACTER        OPTPST(64)
C
C VARIABLES LOCALES
C
      INTEGER          II
C
C***********************************************************************
C
C Transmission des paramètres du common Fortran à la partie C.
C
      IPSTVL = ICHRVL
      IPSTBO = ICHRBO
      IPSTSY = ICHRSY
C
      IPSTMD = ICHRMD
      NTPST = NTCHR
C
      DO II = 1, LEN(FMTCHR)
        FMTPST(II) = FMTCHR(II:II)
      ENDDO
      DO II = 1, LEN(OPTCHR)
        OPTPST(II) = OPTCHR(II:II)
      ENDDO
C
C***********************************************************************
C
      RETURN
C
      END
c@z
