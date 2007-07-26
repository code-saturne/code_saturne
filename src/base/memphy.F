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
                        SUBROUTINE MEMPHY
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NVAR   , NCELET , NCEL   , NFAC   , NFABOR , NPHAS  ,
     &   IW1    , IW2    , IW3    , IW4    ,
     &   IW5    , IW6    , IW7    , IW8    ,
     &   IW9    , IW10   , IXMIJ  ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE CALCUL DES PROPRIETES PHYSIQUES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/IDBRA0! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DES!
CARGU !              !    !     !  TABLEAUX IA/RA                      !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NVAR         ! E  !  -> ! NOMBRE DE VARIABLES                  !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! IW1,2,3,4,5,6! E  ! <-  ! "POINTEUR" SUR W1,2,3,4,5,6,7,8,9,10 !
CARGU ! IXMIJ        ! E  ! <-  ! "POINTEUR" SUR XMIJ                  !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
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
      INCLUDE "optcal.h"
C
C***********************************************************************
C
C
      INTEGER          IDBIA0 ,IDBRA0
      INTEGER          NVAR
      INTEGER          NCELET , NCEL   , NFAC   , NFABOR, NPHAS
      INTEGER          IW1    , IW2    , IW3    , IW4
      INTEGER          IW5    , IW6    , IW7    , IW8
      INTEGER          IW9    , IW10   , IXMIJ
      INTEGER          IFINIA , IFINRA
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IMEMPH , IPHAS
C
C***********************************************************************
C
C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
      IFINIA =       IDEBIA
C
      IW1    =       IDEBRA
      IW2    =       IW1    + NCELET
      IW3    =       IW2    + NCELET
      IW4    =       IW3    + NCELET
      IW5    =       IW4    + NCELET
      IW6    =       IW5    + NCELET
      IW7    =       IW6    + NCELET
      IW8    =       IW7    + NCELET
      IW9    =       IW8    + NCELET
C
      IMEMPH = 0
      DO IPHAS =1, NPHAS
        IF(ITURB(IPHAS).EQ.41)  IMEMPH = 1
      ENDDO
      IF(IMEMPH.EQ.1) THEN
        IW10   = IW9  + NCELET
        IXMIJ  = IW10 + NCELET
        IFINRA = IXMIJ+ 6*NCELET
      ELSE
        IW10   = IW9
        IXMIJ  = IW10
        IFINRA = IXMIJ
      ENDIF
C
C
C---> VERIFICATION
C
      CALL IASIZE('MEMPHY',IFINIA)
C     ===========
C
      CALL RASIZE('MEMPHY',IFINRA)
C     ===========
C
      RETURN
      END
c@z
