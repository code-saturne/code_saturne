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
                        SUBROUTINE PRCPOL
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , IINVPE , IPOL   , ISYM   , NFECRA ,
     &   IFACEL ,
     &   DAM    , XAM    , RK     , GK     , WK     )
C     -------------------------------------------------------------
C***********************************************************************
C
C FONCTION :
C --------
c@foncb
CFONC                                                    -1
CFONC                PRECONDITIONNEMENT DU RESIDU  GK = C  RK
CFONC AFIN D'AUGMENTER L'ORDRE DU POLYNOME SANS SURCOUT IMPORTANT, ON
CFONC "IMPLICITERA" LA RESOLUTION VIA UN MAILLAGE DAMIER (2 COULEURS)
CFONC                                          -1
CFONC IPOL = 0 : POLYNOMIAL ORDRE 0 ==> GK = DA  RK
CFONC                                            -1   -1      -1
CFONC IPOL = 1 : POLYNOMIAL ORDRE 1 ==> GK = ( DA  -DA. XA. DA  ) RK
CFONC
CFONC IPOL = 2 : POLYNOMIAL ORDRE 2 ==>
CFONC                 -1   -1      -1    -1      -1      -1
CFONC        GK = ( DA  -DA. XA. DA   +DA. XA. DA. XA. DA  ) RK
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
CARGU ! IINVPE       ! E  !  -> ! INDICATEUR POUR ANNULER LES INCREMENT!
CARGU !              !    !     ! EN PERIODICITE DE ROTATION (=2) OU   !
CARGU !              !    !     ! POUR LES ECHANGER NORMALEMENT DE     !
CARGU !              !    !     ! MANIERE SCALAIRE (=1)                !
CARGU ! IPOL         ! E  !  -> ! DEGRE DU POLYNOME                    !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYMETRIQUE    !
CARGU !              !    !     !            = 2 MATRICE NON SYMETRIQUE!
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! DAM(NCELET   ! TR !  -> ! DIAGONALE DE MATRICE A PRECONDITIONN !
CARGU ! XAM(NFAC,*)  ! TR !  -> ! EXTRADIAGONALE DE MATRICE A PRECOND  !
CARGU ! RK(NCELET    ! TR !  -> ! VECTEUR RESIDU DONNE                 !
CARGU ! GK(NCELET    ! TR ! <-  ! VECTEUR RESULTANT                    !
CARGU ! WK(NCELET    ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL   , NFAC
      INTEGER          IINVPE , IPOL   , ISYM
      INTEGER          NFECRA
C
      INTEGER          IFACEL(2,NFAC)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,ISYM)
      DOUBLE PRECISION RK(NCELET), GK(NCELET), WK(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL    , IIPOL
C
C***********************************************************************
C
C
      IF(ISYM.NE.1.AND.ISYM.NE.2) THEN
        WRITE(NFECRA,1000)ISYM
        CALL CSEXIT (1)
      ENDIF
C
C     0 - POLYNOME DE DEGRE 0 (PRECONDITIONNEMENT DIAGONAL)
C       ---------------------------------------------------
C
      IF(IPOL.EQ.0) THEN
C
        DO IEL = 1, NCEL
          GK(IEL) = RK(IEL)/DAM(IEL)
        ENDDO
C
C     N - POLYNOME DE DEGRE N
C       ---------------------
C                        N=1                    N=2
C                 -1    -1       -1     -1       -1       -1
C        GK = ( DA  - DA  .XA. DA   + DA  .XA. DA  .XA. DA  + ... ) RK
C
C                         N=1            N=2         N=3
C                -1            -1           -1           -1
C        GK =  DA  (Id - XA. DA  (1 - XA. DA  (1 - XA. DA  (1...) ))) RK
C
      ELSEIF(IPOL.GT.0) THEN
C
      DO IEL = 1, NCEL
        GK(IEL) = RK(IEL)/DAM(IEL)
      ENDDO
C
        DO IIPOL = 1, IPOL
C
C        CALCUL DE WK = XA. GK
          CALL PROXAV(NCELET, NCEL, NFAC, ISYM, IINVPE,
     &                IFACEL, XAM, GK, WK)
C                         -1
C        CALCUL DE GK = DA  (RK - XA.GK)
          DO IEL = 1, NCEL
            GK(IEL) = (RK(IEL) - WK(IEL))/DAM(IEL)
          ENDDO
        ENDDO
C
      ENDIF

C
C--------
C FORMATS
C--------
C
 1000 FORMAT(' PRCPOL APPELE AVEC ISYM =',I10)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
