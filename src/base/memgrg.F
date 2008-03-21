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
                        SUBROUTINE MEMGRG
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCEL   , NFAC   , ISYM   ,
     &   JIVOIS , JIP    , JIW    , JINDIF , JIRESF ,
     &   JICELF , JICELC , JRW    ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE MULTIGRILLE POUR CREATION GRILLE GROSSIERE
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
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! JIVOIS       ! E  ! <-  ! POINTEUR SUR IVOIS                   !
CARGU ! JIP          ! E  ! <-  ! POINTEUR SUR IP                      !
CARGU ! JIW          ! E  ! <-  ! POINTEUR SUR IW                      !
CARGU ! JINDIF       ! E  ! <-  ! POINTEUR SUR INDICF                  !
CARGU ! JIRESF       ! E  ! <-  ! POINTEUR SUR IRESFA                  !
CARGU ! JICELF       ! E  ! <-  ! POINTEUR SUR ICELFA                  !
CARGU ! JICELC       ! E  ! <-  ! POINTEUR SUR ICELCE                  !
CARGU ! JIRW         ! E  ! <-  ! POINTEUR SUR RW                      !
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
C
C***********************************************************************
C
C
      INTEGER          IDBIA0 ,IDBRA0
      INTEGER          NCEL   , NFAC   , ISYM
      INTEGER          JIVOIS , JIP    , JIW    , JINDIF , JIRESF
      INTEGER          JICELF , JICELC , JRW
      INTEGER          IFINIA , IFINRA
C
      INTEGER          IDEBIA, IDEBRA
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
      JIVOIS = IDEBIA
      JIP    = JIVOIS + NCEL
      JIW    = JIP    + NCEL
      JINDIF = JIW    + NCEL
      JIRESF = JINDIF + NFAC
      JICELF = JIRESF + NFAC
      JICELC = JICELF + NFAC*2
      IFINIA = JICELC + NFAC*2
C
      JRW    = IDEBRA
      IFINRA = JRW    + NCEL
C
C---> VERIFICATION
C
      CALL IASIZE('MEMGRG',IFINIA)
C     ===========
C
      CALL RASIZE('MEMGRG',IFINRA)
C     ===========
C
      RETURN
      END
c@z
