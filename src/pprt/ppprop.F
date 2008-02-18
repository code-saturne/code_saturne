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
                        SUBROUTINE PPPROP
C                       *****************
C     -------------------------------------------------------------
     & ( IPROPP , IPPPST )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DES POSITIONS DES VARIABLES D'ETAT SELON
CFONC    LE TYPE DE PHYSIQUE PARTICULIERE
CFONC    (DANS VECTEURS PROPCE, PROPFA, PROPFB)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IPROPP       ! E  ! <-> ! NUMERO DE LA DERNIERE PROPRIETE      !
CARGU !              !    !     !  (LES PROPRIETES SONT DANS PROPCE,   !
CARGU !              !    !     !   PROPFA ou PRPFB)                   !
CARGU ! IPPPST       ! E  ! <-> ! POINTEUR INDIQUANT LE RANG DE LA     !
CARGU !              !    !     !  DERNIERE GRANDEUR DEFINIE AUX       !
CARGU !              !    !     !  CELLULES (RTP,PROPCE...) POUR LE    !
CARGU !              !    !     !  POST TRAITEMENT                     !
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER       IPROPP, IPPPST
C
C***********************************************************************
C
C ---> Physique particuliere : Combustion Gaz
C
      IF ( IPPMOD(ICOD3P).GE.0 .OR. IPPMOD(ICOEBU).GE.0
     &                         .OR. IPPMOD(ICOLWC).GE.0 ) THEN
        CALL COPROP(IPROPP,IPPPST)
C       ===========
      ENDIF
C
C ---> Physique particuliere :  Combustion Charbon Pulverise
C
      IF ( IPPMOD(ICP3PL).GE.0 ) THEN
        CALL CPPROP(IPROPP,IPPPST)
C       ===========
      ENDIF
C
C ---> Physique particuliere :  Combustion Charbon Pulverise
C      Couplee Transport Lagrangien des particules de charbon
C
      IF ( IPPMOD(ICPL3C).GE.0 ) THEN
        CALL CPLPRO (IPROPP,IPPPST)
C       ===========
      ENDIF
C
C ---> Physique particuliere : Combustion Fuel
C
      IF ( IPPMOD(ICFUEL).GE.0 ) THEN
        CALL FUPROP(IPROPP,IPPPST)
C       ===========
      ENDIF
C
C ---> Physique particuliere : Compressible
C
      IF ( IPPMOD(ICOMPF).GE.0 ) THEN
        CALL CFPROP(IPROPP,IPPPST)
C       ===========
      ENDIF
C
C ---> Physique particuliere : Versions electriques
C
      IF ( IPPMOD(IELJOU).GE.1 .OR.
     &     IPPMOD(IELARC).GE.1 .OR.
     &     IPPMOD(IELION).GE.1       ) THEN
        CALL ELPROP(IPROPP,IPPPST)
C       ===========
      ENDIF
C
      END
c@z
