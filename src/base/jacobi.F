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
                        SUBROUTINE JACOBI
C                       *****************
C     -------------------------------------------------------------
     & ( CNOM   , NCELET , NCEL   , NFAC   ,
     &   ISYM   , NITMAP , IINVPE , IWARNP , NFECRA , NITERF ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL ,
     &   DAM    , XAM    , SMBRP  , VX     , RK     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DU SYSTEME (DAM+ XAM).VX = SMBRP PAR UNE METHODE DE
CFONC                          JACOBI OU GAUSS SEIDEL RELAXEE
CFONC
CFONC VX SUPPOSE INITIALISE EN ENTREE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                       ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! CNOM         ! A  !  -> ! NOM DE LA VARIABLE                   !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! NITMAP       ! E  !  -> ! NOMBRE MAX D'ITER POUR RESOL ITERATIV!
CARGU ! IINVPE       ! E  !  -> ! INDICATEUR POUR ANNULER LES INCREMENT!
CARGU !              !    !     ! EN PERIODICITE DE ROTATION (=2) OU   !
CARGU !              !    !     ! POUR LES ECHANGER NORMALEMENT DE     !
CARGU !              !    !     ! MANIERE SCALAIRE (=1)                !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! NITERF       ! E  ! <-  ! NOMBRE D'ITERATIONS EFFECTUEES       !
CARGU ! EPSILP       ! R  !  -> ! PRECISION POUR RESOL ITER            !
CARGU ! RNORM        ! R  !  -> ! NORMALISATION DU RESIDU              !
CARGU ! RESIDU       ! R  ! <-  ! RESIDU FINAL NON NORME               !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! DAM(NCELET   ! TR !  -> ! DIAGONALE DE LA MATRICE              !
CARGU ! XAM(NFAC,ISYM! TR !  -> ! EXTRADIAGONALE DE LA MATRICE         !
CARGU ! SMBRP(NCELET ! TR !  -> ! SECOND MEMBRE DU SYSTEME             !
CARGU ! VX   (NCELET ! TR ! <-> ! SOLUTION      DU SYSTEME             !
CARGU ! RK   (NCELET ! TR !  -  ! AUXILIAIRE DE TRAVAIL                !
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
      INCLUDE "cstnum.h"
      INCLUDE "vector.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*8      CNOM
      INTEGER          NCELET , NCEL   , NFAC
      INTEGER          ISYM   , NITMAP , IINVPE , IWARNP , NFECRA
      INTEGER          NITERF
      INTEGER          IDIMTE , ITENSO
      DOUBLE PRECISION EPSILP , RNORM  , RESIDU
C
      INTEGER          IFACEL(2,NFAC)
      DOUBLE PRECISION DAM(NCELET ),XAM(NFAC ,ISYM),SMBRP(NCELET)
      DOUBLE PRECISION VX(NCELET),RK(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IFAC, IEL, II, JJ, ISQRT
C
C***********************************************************************
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
      IF(ISYM.NE.2.AND.ISYM.NE.1) THEN
        WRITE(NFECRA,9000)CNOM,ISYM
        CALL CSEXIT (1)
      ENDIF
C
      NITERF = 0
C
      ISQRT = 1
      CALL PRODSC(NCELET,NCEL,ISQRT,SMBRP,SMBRP,RESIDU)
C
      IF( RNORM.LE.EPZERO .OR. RESIDU.LE.EPZERO ) THEN
        IF(IWARNP.GE.2) THEN
          WRITE (NFECRA,1000) CNOM,RNORM,RESIDU
        ENDIF
        RETURN
      ENDIF
C
C
C=======================================================================
C 2.  ITERATION COURANTE (STD OU MAILLAGE PAIR/IMPAIR)
C=======================================================================
C
 100  CONTINUE
      NITERF = NITERF +1
C
C
      DO IEL=1,NCEL
        RK(IEL) = VX(IEL)
      ENDDO
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) CALL PARCOM (RK)
C                     ===========
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
C
        IF(IINVPE.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    RK     , RK     , RK    ,
     &    RK     , RK     , RK    ,
     &    RK     , RK     , RK    )
        ELSEIF(IINVPE.EQ.2) THEN
          IDIMTE = 0
          ITENSO = 11
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    RK     , RK     , RK    ,
     &    RK     , RK     , RK    ,
     &    RK     , RK     , RK    )
        ENDIF
C
      ENDIF
C
C ---> CAS STANDARD
C
      DO IEL = 1, NCEL
        VX(IEL) = SMBRP(IEL)
      ENDDO
      IF(NCELET.GT.NCEL) THEN
        DO IEL = NCEL+1, NCELET
          VX(IEL) = 0.D0
        ENDDO
      ENDIF
C
C
C On distingue les deux cas pour permettre au compilateur
C   d'optimiser (ISYM = 1 et ISYM = 2).
C   La seule difference entre les deux blocs est l'indice de
C     XAM(IFAC,.) dans la ligne
C     VX(JJ) = VX(JJ) -XAM(IFAC,ISYM)*RK(II)
C
      IF(ISYM.EQ.1) THEN
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VX(II) = VX(II) -XAM(IFAC,1)*RK(JJ)
            VX(JJ) = VX(JJ) -XAM(IFAC,1)*RK(II)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VX(II) = VX(II) -XAM(IFAC,1)*RK(JJ)
            VX(JJ) = VX(JJ) -XAM(IFAC,1)*RK(II)
          ENDDO
C
        ENDIF
C
      ELSE
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VX(II) = VX(II) -XAM(IFAC,1)*RK(JJ)
            VX(JJ) = VX(JJ) -XAM(IFAC,2)*RK(II)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VX(II) = VX(II) -XAM(IFAC,1)*RK(JJ)
            VX(JJ) = VX(JJ) -XAM(IFAC,2)*RK(II)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C
      DO IEL = 1, NCEL
        VX(IEL) = VX(IEL)/DAM(IEL)
      ENDDO
C
C=======================================================================
C 3.  TEST DE CONVERGENCE ET IMPRESSIONS
C=======================================================================
C
      RESIDU = 0.D0
      DO IEL=1,NCEL
        RESIDU = RESIDU +( DAM(IEL)*(VX(IEL)-RK(IEL)) )**2
      ENDDO
      IF (IRANGP.GE.0) CALL PARSOM (RESIDU)
C                      ===========
      RESIDU = SQRT(RESIDU)
C
      IF( RESIDU.GT.EPSILP*RNORM ) THEN
        IF( NITERF.LT.NITMAP ) THEN
          IF(IWARNP.GE.3) THEN
            WRITE (NFECRA,1100)CNOM,NITERF,RESIDU,RESIDU/RNORM
          ENDIF
          GOTO 100
        ELSE
          IF(IWARNP.GE.1) THEN
            WRITE (NFECRA,1100)CNOM,NITERF,RESIDU,RESIDU/RNORM
            WRITE (NFECRA,1200)CNOM
          ENDIF
        ENDIF
      ELSE
        IF(IWARNP.GE.2) THEN
          WRITE (NFECRA,1100)CNOM,NITERF,RESIDU,RESIDU/RNORM
        ENDIF
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(1X,A8,
     & ' JACOBI: SORTIE IMMEDIATE: RNORM:', E11.4,' RESIDU:',E11.4)
 1100 FORMAT (1X,A8,
     &  ' : JACOBI NITER: ',I5,' RES ABS:',E11.4,' RES NOR:',E11.4 )
 1200 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' NON CONVERGENCE DE JACOBI           ',/,
     &'@    *********                                               ',/,
     &'@                                                            '  )
 9000 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS jacobi                           ',/,
     &'@    *********                                               ',/,
     &'@     APPEL DE jacobi POUR ',A8 ,' AVEC ISYM   = ',I10        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C12345678 : JACOBI NITER: 10000 RES ABS: 1234567890 RES NOR: 124567890
C----
C FIN
C----
C
      RETURN
C
      END
c@z
