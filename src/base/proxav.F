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
                        SUBROUTINE PROXAV
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , ISYM   , IINVPE ,
     &   IFACEL , XA     , VX     , VY     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC PRODUIT MATRICE VECTEUR Y = (XA).X
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYMETRIQUE    !
CARGU !              !    !     !            = 2 MATRICE NON SYMETRIQUE!
CARGU ! IINVPE       ! E  !  -> ! INDICATEUR POUR ANNULER LES INCREMENT!
CARGU !              !    !     ! EN PERIODICITE DE ROTATION (=2) OU   !
CARGU !              !    !     ! POUR LES ECHANGER NORMALEMENT DE     !
CARGU !              !    !     ! MANIERE SCALAIRE (=1)                !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! XA(NFAC,ISYM)! TR !  -> ! EXTRA DIAGONALE DE LA MATRICE        !
CARGU ! VX(NCELET    ! TR !  -> ! VECTEUR A MULTIPLIER                 !
CARGU ! VY(NCELET    ! TR ! <-  ! VECTEUR RESULTAT                     !
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
      INCLUDE "vector.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL   , NFAC   , ISYM , IINVPE
      INTEGER          IFACEL(2,NFAC)
      DOUBLE PRECISION XA(NFAC,ISYM),VX(NCELET),VY(NCELET)
C
C VARAIBLES LOCALES
C
      INTEGER IFAC,II,JJ,IEL,IDIMTE,ITENSO
C
C***********************************************************************
C
C     1 - PRODUIT MATRICE/VECTEUR PAR LA DIAGONALE
C     --------------------------------------------
C
      DO IEL = 1, NCELET
        VY(IEL) = 0.D0
      ENDDO
C
C
C     2 - PRODUIT MATRICE/VECTEUR TERMES X-TRADIAGONAUX
C     -------------------------------------------------
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) CALL PARCOM (VX)
C                     ===========
C
C --> TRAITEMENT DE LA PERIODICITE
      IF(IPERIO.EQ.1) THEN
        IF(IINVPE.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    VX     , VX     , VX    ,
     &    VX     , VX     , VX    ,
     &    VX     , VX     , VX    )
         ELSEIF(IINVPE.EQ.2) THEN
            IDIMTE = 0
            ITENSO = 11
            CALL PERCOM
C           ===========
     &  ( IDIMTE , ITENSO ,
     &    VX     , VX     , VX    ,
     &    VX     , VX     , VX    ,
     &    VX     , VX     , VX    )
C
C        Utile a Codits (produit avec une variable non en increment)
          ELSEIF(IINVPE.EQ.3) THEN
            IDIMTE = 0
            ITENSO = 1
            CALL PERCOM
C           ===========
     &  ( IDIMTE , ITENSO ,
     &    VX     , VX     , VX    ,
     &    VX     , VX     , VX    ,
     &    VX     , VX     , VX    )
          ENDIF
      ENDIF
C
C
      IF( ISYM.EQ.1 ) THEN
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VY(II) = VY(II) +XA(IFAC,1)*VX(JJ)
            VY(JJ) = VY(JJ) +XA(IFAC,1)*VX(II)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VY(II) = VY(II) +XA(IFAC,1)*VX(JJ)
            VY(JJ) = VY(JJ) +XA(IFAC,1)*VX(II)
          ENDDO
C
        ENDIF
C
      ELSE
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VY(II) = VY(II) +XA(IFAC,1)*VX(JJ)
            VY(JJ) = VY(JJ) +XA(IFAC,2)*VX(II)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VY(II) = VY(II) +XA(IFAC,1)*VX(JJ)
            VY(JJ) = VY(JJ) +XA(IFAC,2)*VX(II)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C--------
C FORMATS
C--------
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
