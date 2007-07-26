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
                        SUBROUTINE VECTDS
C                       *****************
C     -------------------------------------------------------------
     & ( NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   IFACEL , IFABOR , IA     ,
     &   SURFAC , SURFBO , POND   ,
     &   VECTX  , VECTY  , VECTZ  ,
     &   VALF   , VALB   , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL A LA FACE DE (VECT)ij . S
CFONC   A PARTIR DU VECTEUR VECTX, VECTY, VECTZ
CFONC  UTILISE POUR LE CALCUL DU TERME DE DIFFUSION DE Rij ET Epsilon
CFONC   EN Rij-Epsilon LRR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IFABOR(NFABOR! TE !  -> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! SURFAC(3,NFAC! TR !  -> ! SURF VECTORIELLE DES SURFACES INTERNE!
CARGU ! SURFBO(      ! TR !  -> ! SURF VECTORIELLE DES SURFACES        !
CARGU !  (3,NFABOR   !    !     !                               DE BORD!
CARGU ! POND(NFAC    ! TR !  -> ! PONDERATION POUR INTERPOLATION FACES !
CARGU ! VECTX (NCELET! TR !  -> ! COMPOSANTE X DU VECTEUR   ENTRE      !
CARGU ! VECTY (NCELET! TR !  -> ! COMPOSANTE Y DU VECTEUR   ENTRE      !
CARGU ! VECTZ (NCELET! TR !  -> ! COMPOSANTE Z DU VECTEUR   ENTRE      !
CARGU ! VALF (NFAC)  ! TR ! <-  ! VECT*SURFACE      AUX FACES INTERNES !
CARGU ! VALB (NFABOR ! TR ! <-  ! VECT*SURFACE      AUX FACES DE BORD  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
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
C-----------------------------------------------------------------------
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NDIM, NCELET, NCEL, NFAC, NFABOR
      INTEGER          IFACEL(2,NFAC), IFABOR(NFABOR)
      INTEGER          IA(*)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION POND(NFAC)
      DOUBLE PRECISION VECTX(NCELET), VECTY(NCELET), VECTZ(NCELET)
      DOUBLE PRECISION VALF(NFAC), VALB(NFABOR)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IFAC, IEL1, IEL2
      INTEGER          IDIMTE, ITENSO
      DOUBLE PRECISION VALFX, VALFY, VALFZ
C
C***********************************************************************
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (VECTX)
C       ===========
        CALL PARCOM (VECTY)
C       ===========
        CALL PARCOM (VECTZ)
C       ===========
      ENDIF
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 1
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    VECTX  , VECTX  , VECTX ,
     &    VECTY  , VECTY  , VECTY ,
     &    VECTZ  , VECTZ  , VECTZ )
      ENDIF
C
C
      DO IFAC = 1 , NFAC
C
        IEL1 = IFACEL(1,IFAC)
        IEL2 = IFACEL(2,IFAC)
C
        VALFX =       POND(IFAC)  * VECTX(IEL1) +
     &          (1.D0-POND(IFAC)) * VECTX(IEL2)
        VALFY =       POND(IFAC)  * VECTY(IEL1) +
     &          (1.D0-POND(IFAC)) * VECTY(IEL2)
        VALFZ =       POND(IFAC)  * VECTZ(IEL1) +
     &          (1.D0-POND(IFAC)) * VECTZ(IEL2)
C
        VALF(IFAC) = VALFX*SURFAC(1,IFAC) +
     &               VALFY*SURFAC(2,IFAC) +
     &               VALFZ*SURFAC(3,IFAC)
       ENDDO
C
       DO IFAC = 1 , NFABOR
C
C     On met VALB a zero, ce qui revient a negliger la partie
C       extradiagonale du tenseur de diffusion au bord.
CMO          IEL1 = IFABOR(IFAC)
CMOC
CMO          VALB(IFAC) = VECTX(IEL1)*SURFBO(1,IFAC) +
CMO     &                 VECTY(IEL1)*SURFBO(2,IFAC) +
CMO     &                 VECTZ(IEL1)*SURFBO(3,IFAC)
          VALB(IFAC) = 0.D0
C
       ENDDO
C
       RETURN
       END
c@z
