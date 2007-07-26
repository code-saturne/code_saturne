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
                        SUBROUTINE RESMGR
C                       *****************
C     -------------------------------------------------------------
     & ( CNOM   , IDBIA0 , IDBRA0 , NGR    , NCELET , NCEL   , NFAC   ,
     &   ISYM   , IPOL   , IRESMO , NITSMO , NCYMXP , IINVPE ,
     &   IWARNP , NFECRA , ICYCLE ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL , IRESPR , IFACEW ,
     &   DAF    , XAF    , DAW    , XAW    , BF     , BW     ,
     &   XF     , XW     , W1     , W2     , W3     , W4     ,
     &   W5     , W6     ,
     &   IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC   MULTIGRILLE :
CFONC RESOLUTION DU SYSTEME (DAF+ XAF).XF = BF PAR UN ALGORITHME
CFONC MULTIGRILLE ALGEBRIQUE A PLUSIEURS NIVEAUX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! CNOM         ! A  !  -> ! NOM DE LA VARIABLE                   !
CARGU ! NGR          ! E  !  -> ! NOMBRE DE NIVEAUX DE GRILLES         !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! IPOL         ! E  !  -> ! DEGRE DU POLYNOME POUR PRECOND       !
CARGU !              !    !     !         (0 -> DIAGONAL)              !
CARGU ! IRESMO(0:NGR)! E  !  -> ! INDICATEUR = 0 GRADCO                !
CARGU !              !    !     !            = 1 JACOBI                !
CARGU ! NITSMO       ! E  !  -> ! NOMBRE MAX D'ITER POUR RESOL ITERATIV!
CARGU !  (2,0:NGR)   !    !     !                                      !
CARGU ! IINVPE       ! E  !  -> ! INDICATEUR POUR ANNULER LES INCREMENT!
CARGU !              !    !     ! EN PERIODICITE DE ROTATION (=2) OU   !
CARGU !              !    !     ! POUR LES ECHANGER NORMALEMENT DE     !
CARGU !              !    !     ! MANIERE SCALAIRE (=1)                !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! ICYCLE       ! E  ! <-  ! NOMBRE DE CYCLES EFFECTUES MOINS 1   !
CARGU ! EPSILP       ! R  !  -> ! PRECISION POUR RESOL ITER            !
CARGU ! RNORM        ! R  !  -> ! NORMALISATION DU RESIDU              !
CARGU ! RESIDU       ! R  ! <-  ! RESIDU FINAL NON NORME               !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IRESPR(NCEL) ! TE !  -> ! RESTR. DES RESIDUS ET PROL. DES SOL  !
CARGU ! IFACEW       ! TE !  -> ! NUMERO DES MACROS VOISINS D'UNE      !
CARGU !      (2,NFAC)!    !     !  MACRO FACE                          !
CARGU ! DAF(NCEL)    ! TR !  -> ! DIAGONALE MAILLAGE FIN               !
CARGU ! XAF(NFAC,ISYM! TR !  -> ! EXTRADIAGONALE MAILLAGE FIN          !
CARGU ! DAW(NCEG)    ! TR !  -  ! DIAGONALE MAILLAGE GROSSIER          !
CARGU ! XAW(NFAC,ISYM! TR !  -  ! EXTRADIAGONALE MAILLAGE GROSSIER     !
CARGU ! BF   (NCEL)  ! TR !  -> ! SECOND MEMBRE MAILLAGE FIN           !
CARGU ! BW   (NCEL)  ! TR !  -  ! SECOND MEMBRE MAILLAGE GROSSIER      !
CARGU ! XF   (NCEL)  ! TR ! <-> ! SOLUTION      DU SYSTEME             !
CARGU ! XW   (NCEL)  ! TR !  -  ! TABLEAU DE TRAVAIL MAILLAGE GROSSIER !
CARGU ! W?(NCEL)     ! TR !  -  ! AUXILIAIRES DE TRAVAIL               !
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
C
C***********************************************************************
C
      CHARACTER*8      CNOM
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NGR , NCELET , NCEL , NFAC
      INTEGER          ISYM , IPOL , NCYMXP , IINVPE, IWARNP , NFECRA
      INTEGER          ICYCLE
      DOUBLE PRECISION EPSILP , RNORM  , RESIDU
C
      INTEGER          IRESMO(0:NGR) , NITSMO(2,0:NGR)
      INTEGER          IFACEL(2,NFAC) , IRESPR(NCEL) , IFACEW(2,NFAC)
      INTEGER          IA(*)
      DOUBLE PRECISION DAF(NCEL), XAF(NFAC,ISYM), BF(NCEL)
      DOUBLE PRECISION DAW(NCEL), XAW(NFAC,ISYM), BW(NCEL)
      DOUBLE PRECISION XF(NCEL) , XW(NCEL)
      DOUBLE PRECISION W1(NCEL) , W2(NCEL) , W3(NCEL) , W4(NCEL)
      DOUBLE PRECISION W5(NCEL) , W6(NCEL)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFIN , IEL , ISQRT , IWARMG
      INTEGER          NITFIN
C
C***********************************************************************
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      ICYCLE = 0
      NITFIN = 0
C
      ISQRT = 1
      CALL PRODSC(NCELET,NCEL,ISQRT,BF,BF,RESIDU)
      IF( RNORM.LE.EPZERO .OR. RESIDU.LE.EPZERO ) THEN
        IF(IWARNP.GE.2) THEN
          WRITE (NFECRA,1000) RNORM, RESIDU
        ENDIF
        RETURN
      ENDIF
C
      IWARMG = IWARNP
C
      DO IEL = 1, NCEL
        W1(IEL) = 0.D0
        W2(IEL) = 0.D0
        W3(IEL) = 0.D0
        W4(IEL) = 0.D0
        W5(IEL) = 0.D0
        W6(IEL) = 0.D0
        XW(IEL) = 0.D0
        BW(IEL) = 0.D0
      ENDDO
C
      IF (IWARNP.GE.2) THEN
        WRITE( NFECRA,1100) CNOM
      ENDIF
C
C=======================================================================
C 2.  CYCLES EN V
C=======================================================================
C
 100  CONTINUE
      ICYCLE = ICYCLE +1
      IFIN = 0
C
      IF (IWARNP.GE.2) THEN
        WRITE( NFECRA,1110) ICYCLE
      ENDIF
C
      CALL CYCMGR
C     ===========
     & ( CNOM   , IDEBIA , IDEBRA ,
     &   NGR    , NCELET , NCEL   , NFAC   , ISYM   , IPOL   ,
     &   IINVPE , IWARMG , NFECRA , ICYCLE , IFIN   ,
     &   IRESMO , NITSMO , NCYMXP ,
     &   NITFIN ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL , IRESPR , IFACEW ,
     &   DAF    , XAF    , DAW    , XAW    , BF     , BW     ,
     &   XF     , XW     , W1     , W2     , W3     , W4     ,
     &   W5     , W6     ,
     &   IA     , RA      )
C
      IF (IFIN.NE.0) RETURN
C
      GOTO 100
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(1X,A8,
     & ' RESMGR : SORTIE IMMEDIATE: RNORM: ', E11.4,' RESIDU: ',E11.4)
C
 1100  FORMAT(/,1X,A8,' : RESMGR')
 1110  FORMAT(/,10X,'CYCLE NUMERO : ',I5)
C
C----
C FIN
C----
C
      END
c@z
