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
                        SUBROUTINE GRADCO
C                       *****************
C     -------------------------------------------------------------
     & ( CNOM   , NCELET , NCEL   , NFAC   ,
     &   ISYM   , IPOL   , NITMAP , IINVPE ,
     &                              IWARNP , NFECRA , NITERF ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL ,
     &   DAM    , XAM    , SMBRP  , VX     ,
     &   RK     , DK     , GK     , ZK     , WK     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DU SYSTEME (DAM+ XAM).VX = SMBRP PAR UNE METHODE DE
CFONC                          GRADIENT CONJUGUE PRECONDITIONNE
CFONC   VX SUPPOSE INITIALISE EN ENTREE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
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
CARGU ! IPOL         ! E  !  -> ! DEGRE DU POLYNOME POUR PRECOND       !
CARGU !              !    !     !         (0 -> DIAGONAL)              !
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
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! DAM(NCELET   ! TR !  -> ! DIAGONALE DE LA MATRICE              !
CARGU ! XAM(NFAC,ISYM! TR !  -> ! EXTRADIAGONALE DE LA MATRICE         !
CARGU ! SMBRP(NCELET ! TR !  -> ! SECOND MEMBRE DU SYSTEME             !
CARGU ! VX   (NCELET ! TR ! <-> ! SOLUTION      DU SYSTEME             !
CARGU ! RK,DK,GK,ZK  ! TR !  -  ! AUXILIAIRES DE TRAVAIL               !
CARGU !      (NCELET !    !     !                                      !
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
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*8      CNOM
      INTEGER          NCELET , NCEL   , NFAC
      INTEGER          ISYM   , IPOL   , NITMAP , IINVPE
      INTEGER          IWARNP , NFECRA
      INTEGER          NITERF
      DOUBLE PRECISION EPSILP , RNORM  , RESIDU
C
      INTEGER          IFACEL(2,NFAC)
      DOUBLE PRECISION DAM(NCELET ),XAM(NFAC ,ISYM),SMBRP(NCELET)
      DOUBLE PRECISION VX(NCELET)
      DOUBLE PRECISION RK(NCELET),DK(NCELET),GK(NCELET),ZK(NCELET)
      DOUBLE PRECISION WK(NCELET)
C
C
C VARIABLES LOCALES
C
#if defined(_CS_HAVE_BLAS)
      INTEGER          UN
#endif
      INTEGER          IEL, ISQRT, IGRPPS
      DOUBLE PRECISION RO1,RO2,ROK,RKGKM1,RKGK,RGSRG
C
C***********************************************************************
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
      IF (IRANGP.GE.0) THEN
        IGRPPS = 1
      ELSE
        IGRPPS = 0
      ENDIF
C
#if defined(_CS_HAVE_BLAS)
      UN = 1
#endif
      NITERF = 0
C
      ISQRT = 1
      CALL PRODSC(NCELET,NCEL,ISQRT,SMBRP,SMBRP,RESIDU)
C
      IF( RNORM.LE.EPZERO .OR. RESIDU.LE.EPZERO ) THEN
        IF(IWARNP.GE.2) THEN
          WRITE (NFECRA,1000)CNOM, RNORM , RESIDU
        ENDIF
        RETURN
      ENDIF
C
      DO IEL=1,NCEL
        RK(IEL) = 0.D0
        DK(IEL) = 0.D0
        GK(IEL) = 0.D0
        ZK(IEL) = 0.D0
      ENDDO
C
C
C=======================================================================
C 2.  INITIALISATION DU CALCUL ITERATIF
C=======================================================================
C
C     RESIDU ET DIRECTION DE DESCENTE
C
      CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPE,
     &            IFACEL,DAM,XAM,VX,RK)
C
      DO IEL=1,NCEL
        RK(IEL) = RK(IEL) -SMBRP(IEL)
        DK(IEL) = RK(IEL)
      ENDDO
C
C     PRECONDITIONNEMENT POLYNOMIAL D'ORDRE IPOL
C
      CALL PRCPOL(NCELET,NCEL,NFAC,IINVPE,
     &            IPOL,ISYM,NFECRA,IFACEL,DAM,XAM,RK,GK,WK)
C
C     DIRECTION DE DESCENTE
C
#if defined(_CS_HAVE_BLAS)
      CALL DCOPY(NCEL,GK,UN,DK,UN)
#else
      DO IEL=1,NCEL
        DK(IEL) = GK(IEL)
      ENDDO
#endif
C
      ISQRT = 0
      CALL PRODSC(NCELET,NCEL,ISQRT,RK,GK,RKGKM1)
      CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPE,
     &            IFACEL,DAM,XAM,DK,ZK)
C
C     PARAMETRE DE DESCENTE
      ISQRT = 0
      CALL PRODS2(NCELET,NCEL,ISQRT,RK,DK,DK,ZK,RO1,RO2)
      ROK = -RO1/RO2
#if defined(_CS_HAVE_BLAS)
      CALL DAXPY(NCEL,ROK,DK,UN,VX,UN)
      CALL DAXPY(NCEL,ROK,ZK,UN,RK,UN)
#else
      DO IEL=1,NCEL
        VX(IEL) = VX(IEL) +ROK*DK(IEL)
        RK(IEL) = RK(IEL) +ROK*ZK(IEL)
      ENDDO
#endif
C
C     TEST DE CONVERGENCE
C
      ISQRT = 1
      CALL PRODSC(NCELET,NCEL,ISQRT,RK,RK,RESIDU)
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
      GOTO 900
C
C=======================================================================
C 3.  ITERATION COURANTE
C=======================================================================
C
 100   CONTINUE
       NITERF = NITERF +1
C
      IF (IGRPPS.EQ.0) THEN
C
C       CALCUL DU RESIDU
        ISQRT = 1
        CALL PRODSC(NCELET,NCEL,ISQRT,RK,RK,RESIDU)
C
      ELSE
C
C       PRECONDITIONNEMENT POLYNOMIAL
C
        CALL PRCPOL(NCELET,NCEL,NFAC,IINVPE,
     &              IPOL,ISYM,NFECRA,IFACEL,DAM,XAM,RK,GK,WK)
C
C       CALCUL DU RESIDU AVEC PREPARATION DU PARAMETRE DE DESCENTE
        ISQRT = 0
        CALL PRODS2(NCELET,NCEL,ISQRT,RK,RK,RK,GK,RESIDU,RKGK)
        RESIDU = SQRT(RESIDU)
C
      ENDIF
C
C     TEST DE CONVERGENCE
C
 200  CONTINUE
C
      IF( RESIDU.GT.EPSILP*RNORM ) THEN
        IF( NITERF.LT.NITMAP ) THEN
          IF(IWARNP.GE.3) THEN
            WRITE (NFECRA,1100)CNOM,NITERF,RESIDU,RESIDU/RNORM
          ENDIF
          GOTO 300
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
      GOTO 900
C
 300  CONTINUE
C
      IF (IGRPPS.EQ.0) THEN
C
C       PRECONDITIONNEMENT POLYNOMIAL
C
        CALL PRCPOL(NCELET,NCEL,NFAC,IINVPE,
     &              IPOL,ISYM,NFECRA,IFACEL,DAM,XAM,RK,GK,WK)
C
C       PREPARATION DU PARAMETRE DE DESCENTE
        ISQRT = 0
        CALL PRODSC(NCELET,NCEL,ISQRT,RK,GK,RKGK)
C
      ENDIF
C
C     FIN DU CALCUL DU PARAMETRE DE DESCENTE ET PRODUIT MATRICE/VECTEUR
C
      RGSRG = RKGK/RKGKM1
      RKGKM1 = RKGK
C
#if defined(_CS_HAVE_ESSL)
      CALL DZAXPY(NCEL,RGSRG,DK,UN,GK,UN,DK,UN)
#else
      DO IEL=1,NCEL
        DK(IEL) = GK(IEL) +RGSRG*DK(IEL)
      ENDDO
#endif
      CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPE,
     &            IFACEL,DAM,XAM,DK,ZK)
C
      ISQRT = 0
      CALL PRODS2(NCELET,NCEL,ISQRT,RK,DK,DK,ZK,RO1,RO2)
      ROK = -RO1/RO2
C
#if defined(_CS_HAVE_BLAS)
      CALL DAXPY(NCEL,ROK,DK,UN,VX,UN)
#else
      DO IEL=1,NCEL
        VX(IEL) = VX(IEL) +ROK*DK(IEL)
      ENDDO
#endif
C
#if defined(_CS_HAVE_BLAS)
      CALL DAXPY(NCEL,ROK,ZK,UN,RK,UN)
#else
C Sans ce test, le Fujitsu VPP 5000 ne vectorise pas la seconde boucle
C    qu'il trouve trop compliquee ...
      IF(0.LE.1) THEN
      DO IEL=1,NCEL
        RK(IEL) = RK(IEL) +ROK*ZK(IEL)
      ENDDO
      ENDIF
#endif
C
      GOTO 100
C
C=======================================================================
C 4.  FIN
C=======================================================================
C
 900  CONTINUE
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(1X,A8,
     & ' GRADCO: SORTIE IMMEDIATE: RNORM:', E11.4,' RESIDU:',E11.4)
 1100 FORMAT (1X,A8,
     &  ': GRADCO NITER: ',I5,' RES ABS:',E11.4,' RES NOR:',E11.4 )
 1200 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' NON CONVERGENCE DE GRADCO           ',/,
     &'@    *********                                               ',/,
     &'@                                                            '  )
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
