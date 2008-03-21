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
                        SUBROUTINE MEMINI
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IIDEVE , IITUSE ,
     &   IRDEVE , IRTUSE ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC          GESTION MEMOIRE INITIALE
CFONC          COUPLAGE, DEVELOPPEUR, UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/RA0   ! E  !  -> ! DEBUT DE LA PLACE LIBRE DANS IA ER RA!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IIDEVE       ! E  ! <-  ! POINTEUR DE IDEVEL DANS IA           !
CARGU ! IITUSE       ! E  ! <-  ! POINTEUR DE ITUSER DANS IA           !
CARGU ! IRDEVE       ! E  ! <-  ! POINTEUR DE RDEVEL DANS RA           !
CARGU ! IRTUSE       ! E  ! <-  ! POINTEUR DE RTUSER DANS RA           !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS RA EN SORTIE                   !
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
C***********************************************************************
C
C     - REMARQUE : INTERDICTION DE TOUCHER A IDBIA0 IDBRA0
C
      INTEGER IDBIA0 , IDBRA0
      INTEGER NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER NFML   , NPRFML
      INTEGER NNOD   , LNDFAC , LNDFBR
      INTEGER NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER IIDEVE , IITUSE
      INTEGER IRDEVE , IRTUSE
      INTEGER IFINIA , IFINRA
C
      INTEGER IDEBIA, IDEBRA
C
C=======================================================================
C 1. GESTION MEMOIRE DEVELOPPEUR ET UTILISATEUR
C=======================================================================
C
C---> MEMOIRE CONSERVEE JUSQU'A LA FIN
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IIDEVE = IDEBIA
      IITUSE = IIDEVE+NIDEVE
      IFINIA = IITUSE+NITUSE
C
      IRDEVE = IDEBRA
      IRTUSE = IRDEVE+NRDEVE
      IFINRA = IRTUSE+NRTUSE
C
C---> VERIFICATION
C
      CALL IASIZE('MEMINI',IFINIA)
C     ===========
C
      CALL RASIZE('MEMINI',IFINRA)
C     ===========
C
      RETURN
      END
c@z
