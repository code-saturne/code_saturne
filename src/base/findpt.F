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
                        SUBROUTINE FINDPT
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , XYZCEN ,
     &   XX     , YY     , ZZ     , NODE   , NDRANG)
C     -------------------------------------------------------------
C***********************************************************************
c@foncb
CFONC
CFONC  FONCTION
CFONC  --------
CFONC   RECHERCHE DU NOEUD LE PLUS PROCHE DU POINT XX,YY,ZZ
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! XYZCEN(NDIM  ! TR !  -->! TABLE DES COORDONNEES DU             !
CARGU !        NCELET!    !     !           CENTRE DES VOLUMES         !
CARGU ! XX,YY,ZZ     ! TR !  -->! COORDONNEES DU NOEUD CHERCHE         !
CARGU ! NODE         ! E  ! <-- ! NOEUD CHERCHE (NUMEROTATION GLOBALE) !
CARGU !              !    !     !  ZERO SI PLANTAGE                    !
CARGU ! NDRANG       ! E  ! <-- ! RANG DU PROCESSUS ASSOCIE            !
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
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
      INTEGER          NCELET, NCEL, NODE, NDRANG
      DOUBLE PRECISION XYZCEN(3,NCELET)
      DOUBLE PRECISION XX, YY, ZZ
C
      INTEGER          II
      DOUBLE PRECISION XX1, YY1, ZZ1, DIS2, DIS2MN
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      NODE = INT((NCEL+1)/2)
C
      XX1 = XYZCEN(1,NODE)
      YY1 = XYZCEN(2,NODE)
      ZZ1 = XYZCEN(3,NODE)
      DIS2MN = (XX-XX1)**2+(YY-YY1)**2+(ZZ-ZZ1)**2
C
      DO II = 1, NCEL
         XX1 = XYZCEN(1,II)
         YY1 = XYZCEN(2,II)
         ZZ1 = XYZCEN(3,II)
         DIS2 = (XX-XX1)**2+(YY-YY1)**2+(ZZ-ZZ1)**2
         IF (DIS2.LT.DIS2MN) THEN
            NODE = II
            DIS2MN = DIS2
         ENDIF
      ENDDO
C
      IF (IRANGP.GE.0) THEN
         CALL PARFPT (NODE, NDRANG, DIS2MN)
C        ===========
      ELSE
        NDRANG = -1
      ENDIF
C
      RETURN
      END
c@z
