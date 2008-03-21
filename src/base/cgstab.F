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
                        SUBROUTINE CGSTAB
C                       *****************
C     -------------------------------------------------------------
     & ( CNOM   , NCELET , NCEL   , NFAC   ,
     &   ISYM   , IPOL   , NITMAP , IINVPE ,
     &                              IWARNP , NFECRA , NITERF ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL ,
     &   DAM    , XAM    , SMBRP  , VX     ,
     &   RES0   , RK     , PK     , UK     , VK     , ZK     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DU SYSTEME (DAM+ XAM).VX = SMBRP PAR UNE METHODE DE
CFONC   GRADIENT CONJUGUE CARRE STABILISE PRECONDITIONNE
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
CARGU ! DAM(NCELET   ! TR !  -> ! DIAGONALE DE LA MATRICE              !
CARGU ! XAM(NFAC,ISYM! TR !  -> ! EXTRADIAGONALE DE LA MATRICE         !
CARGU ! SMBRP(NCELET ! TR !  -> ! SECOND MEMBRE DU SYSTEME             !
CARGU ! VX   (NCELET ! TR ! <-> ! SOLUTION      DU SYSTEME             !
CARGU ! RES0,RK,PK,UK! TR !  -  ! AUXILIAIRES DE TRAVAIL               !
CARGU ! VK,ZK        !    !     !                                      !
CARGU !      (NCELET !    !     !                                      !
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
      DOUBLE PRECISION RK(NCELET),PK(NCELET),ZK(NCELET)
      DOUBLE PRECISION UK(NCELET),VK(NCELET),RES0(NCELET)
C
C
C VARIABLES LOCALES
C
      INTEGER          IEL, ISQRT
      DOUBLE PRECISION RO1,RO2,ALFA
      DOUBLE PRECISION BETA,BETAM1,GAMA,OMEG,UKRES0
      DOUBLE PRECISION EPZLOC
C
C***********************************************************************
C
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
      IF (NITMAP.LE.0) RETURN
C
C    Precision pour l'arret du calcul
C     (dans l'esprit, doit etre inferieure a EPZER0).
      EPZLOC = 1.D-30
C
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
        RES0(IEL) = 0.D0
        RK(IEL) = 0.D0
        PK(IEL) = 0.D0
        ZK(IEL) = 0.D0
        UK(IEL) = 0.D0
        VK(IEL) = 0.D0
      ENDDO
C
C
C=======================================================================
C 2.  INITIALISATION DU CALCUL ITERATIF
C=======================================================================
C
      CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPE,
     &            IFACEL,DAM,XAM,VX,RES0)
C
      DO IEL = 1, NCEL
        RES0(IEL) = -RES0(IEL) +SMBRP(IEL)
      ENDDO
      DO IEL = 1, NCEL
        RK(IEL) = RES0(IEL)
      ENDDO
C
      ALFA   = 1.D0
      BETAM1 = 1.D0
      GAMA   = 1.D0
C
C=======================================================================
C 3.  ITERATION COURANTE
C=======================================================================
C
 100   CONTINUE
       NITERF = NITERF +1
C
C CALCUL DE BETA ET OMEGA
C
      ISQRT = 0
      CALL PRODSC(NCELET,NCEL,ISQRT,RES0,RK,BETA)
      IF (ABS(BETA).LE.EPZLOC) THEN
        IF (IWARNP.GE.2) THEN
          WRITE (NFECRA,1100)CNOM,NITERF,RESIDU,RESIDU/RNORM
        ENDIF
        RETURN
      ENDIF
C
      IF(ABS(ALFA).LT.EPZLOC)THEN
        WRITE(NFECRA,1400) CNOM, EPZLOC
        CALL CSEXIT (1)
      ENDIF
C
      OMEG = BETA*GAMA/(ALFA*BETAM1)
      BETAM1 = BETA
C
C CALCUL DE PK
C
      DO IEL = 1, NCEL
        PK(IEL) = RK(IEL) + OMEG*(PK(IEL)-ALFA*UK(IEL))
      ENDDO
C
C CALCUL DE ZK = C * PK
C
      CALL PRCPOL(NCELET,NCEL,NFAC,IINVPE,
     &            IPOL,ISYM,NFECRA,IFACEL,DAM,XAM,PK,ZK,UK)
C
C CALCUL DE UK = A * ZK
C
      CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPE,
     &            IFACEL,DAM,XAM,ZK,UK)
C
C CALCUL DE UKRES0 = UK.RES0 ET DE GAMA
C
      ISQRT = 0
      CALL PRODSC(NCELET,NCEL,ISQRT,UK,RES0,UKRES0)
      GAMA = BETA/UKRES0
C
C PREMIERE ACTUALISATION DE VX ET DE RK
C
      DO IEL = 1, NCEL
        VX(IEL) = VX(IEL) + GAMA*ZK(IEL)
        RK(IEL) = RK(IEL) - GAMA*UK(IEL)
      ENDDO
C
C CALCUL DE ZK = C * RK        (ZK PRECEDENT EST ECRASE)
C                              (VK EST UN TABLEAU DE TRAVAIL)
      CALL PRCPOL(NCELET,NCEL,NFAC,IINVPE,
     &            IPOL,ISYM,NFECRA,IFACEL,DAM,XAM,RK,ZK,VK)
C
C CALCUL DE VK = A * ZK
C
      CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPE,
     &            IFACEL,DAM,XAM,ZK,VK)
C
C CALCUL DE ALFA
C
      ISQRT = 0
      CALL PRODS2(NCELET,NCEL,ISQRT,VK,RK,VK,VK,RO1,RO2)
C
C TEST de NULLITE DE RO2
C
      IF(RO2.LT.EPZLOC) THEN
        WRITE(NFECRA,1300) CNOM, EPZLOC
        RETURN
      ENDIF
C
      ALFA = RO1/RO2
C
C ACTUALISATION FINALE DE VX ET DE RK
C
      DO IEL=1,NCEL
        VX(IEL) = VX(IEL) + ALFA*ZK(IEL)
        RK(IEL) = RK(IEL) - ALFA*VK(IEL)
      ENDDO
C
C=======================================================================
C 4.  TEST DE CONVERGENCE ET IMPRESSIONS
C=======================================================================
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
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(1X,A8,
     & ' CGSTAB: SORTIE IMMEDIATE: RNORM:', E11.4,' RESIDU:',E11.4)
 1100 FORMAT (1X,A8,
     &  ': CGSTAB NITER: ',I5,' RES ABS:',E11.4,' RES NOR:',E11.4 )
 1200 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' NON CONVERGENCE DE CGSTAB           ',/,
     &'@    *********                                               ',/,
     &'@                                                            '  )
 1300 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' NON CONVERGENCE DE CGSTAB           ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    Le carre de la norme du vecteur de descente du Bi-CGSTAB',/,
     &'@      est inferieur a ',E12.4                                ,/,
     &'@    La resolution ne progresse plus.                        ',/,
     &'@                                                            '  )
 1400 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' NON CONVERGENCE DE CGSTAB           ',/,
     &'@    *********                                               ',/,
     &'@     ARRET DANS CGSTAB                                      ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@    Le coefficient ALPHA du Bi-CGSTAB est                   ',/,
     &'@      inferieur a ',E12.4                                    ,/,
     &'@                                                            ',/,
     &'@    La matrice ne peut plus etre consideree comme inversible',/,
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
