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
                        SUBROUTINE CYCMGR
C                       *****************
C     -------------------------------------------------------------
     & ( CNOM   , IDBIA0 , IDBRA0 ,
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
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC   MULTIGRILLE :
CFONC UN V-CYCLE
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
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NGR          ! E  !  -> ! NOMBRE DE NIVEAUX DE GRILLES         !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! IPOL         ! E  !  -> ! DEGRE DU POLYNOME POUR PRECOND       !
CARGU !              !    !     !         (0 -> DIAGONAL)              !
CARGU ! IINVPE       ! E  !  -> ! INDICATEUR POUR ANNULER LES INCREMENT!
CARGU !              !    !     ! EN PERIODICITE DE ROTATION (=2) OU   !
CARGU !              !    !     ! POUR LES ECHANGER NORMALEMENT DE     !
CARGU !              !    !     ! MANIERE SCALAIRE (=1)                !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! ICYCLE       ! E  !  -> ! NOMBRE DE CYCLES EFFECTUES MOINS 1   !
CARGU ! IFIN         ! E  ! <-  ! INDICATEUR DE CONVERGENCE            !
CARGU ! IRESMO(0:NGR)! E  !  -> ! INDICATEUR = 0 GRADCO                !
CARGU !              !    !     !            = 1 JACOBI                !
CARGU ! NITSMO       ! E  !  -> ! NOMBRE MAX D'ITER POUR RESOL ITERATIV!
CARGU !  (2,0:NGR)   !    !     ! (1,IGR) : DESCENTE                   !
CARGU !              !    !     ! (2,IGR) : REMONTEE                   !
CARGU ! NITFIN       ! E  ! <-> ! NOMBRE ITERATION MAILLAGE FIN        !
CARGU ! EPSILP       ! R  !  -> ! PRECISION POUR RESOL ITER            !
CARGU ! RNORM        ! R  !  -> ! NORMALISATION DU RESIDU              !
CARGU ! RESIDU       ! R  ! <-  ! RESIDU FINAL NON NORME               !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IRESPR(NCEL) ! TE !  -> ! RESTR. DES RESIDUS ET PROL. DES SOL  !
CARGU ! IFACEW       ! TE !  -  ! NUMERO DES MACROS VOISINS D'UNE      !
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
      INCLUDE "mltgrd.h"
C
C***********************************************************************
C
      CHARACTER*8      CNOM
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NGR , NCELET, NCEL , NFAC , ISYM , IPOL
      INTEGER          IWARMG , NFECRA , ICYCLE , IFIN , IINVPE
      INTEGER          NCYMXP , NITFIN
      DOUBLE PRECISION EPSILP , RNORM  , RESIDU
C
      INTEGER          IRESMO(0:NGR) , NITSMO(2,0:NGR)
      INTEGER          IFACEL(2,NFAC) , IRESPR(NCEL), IFACEW(2,NFAC)
      INTEGER          IA(*)
      DOUBLE PRECISION DAF(NCEL ) , XAF(NFAC,ISYM), BF(NCEL )
      DOUBLE PRECISION DAW(NCEL) , XAW(NFAC,ISYM), BW(NCEL)
      DOUBLE PRECISION XF(NCEL) , XW(NCEL)
      DOUBLE PRECISION W1(NCEL),W2(NCEL),W3(NCEL),W4(NCEL)
      DOUBLE PRECISION W5(NCEL),W6(NCEL)
      DOUBLE PRECISION RA(*)
C
C
C VARIABLES LOCALES
C
      CHARACTER*8      CNOMMG
      CHARACTER*4      NUM
      INTEGER          IDEBIA , IDEBRA , IFINIA , IFINRA
      INTEGER          IEL , IFAC , NCELW, NCELW2 , NFACW
      INTEGER          IGR , NITERF , IPROL , ISQRT
      INTEGER          IWARSO
      INTEGER          IPX(NGRMMX) , IPB(NGRMMX)
      DOUBLE PRECISION EPSPGR
C
C***********************************************************************
C
C=======================================================================
C 0.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IWARSO = IWARMG - 1
C
C PRECISION ACCRUE SUR MAILLAGE LE PLUS GROSSIER
C
      EPSPGR = EPSILP*1D-5
C
C RESERVATION MEMOIRE POUR LE STOCKAGE DES CORRECTIONS ET DES RESIDUS
C POINTEURS IPX ET IPB
C
      IFINIA = IDEBIA
      NCELW  = IA(IPGREN(JNCEL,1))
      IPX(1) = IDEBRA
      IPB(1) = IPX(1) + NCELW
      DO IGR = 2, NGR
        NCELW2 = IA(IPGREN(JNCEL,IGR))
        NCELW  = IA(IPGREN(JNCEL,IGR-1))
        IPX(IGR) = IPB(IGR-1) + NCELW
        IPB(IGR) = IPX(IGR)   + NCELW2
      ENDDO
      IFINRA = IPB(NGR) + NCELW
C
      CALL RASIZE('CYCMGR',IFINRA)
C     ===========
C
      DO IEL = 1, NCEL
        XW(IEL) = XF(IEL)
        BW(IEL) = BF(IEL)
      ENDDO
C
C=======================================================================
C 1.  BOUCLE DE DESCENTE
C=======================================================================
C
      IF (IWARMG.GE.2) THEN
        WRITE(NFECRA,1200)
      ENDIF
C
      DO IGR = 0, NGR-1
C
C 1.1 ON PASSE LES DESCRIPTEURS DE MAILLAGE DANS DES TABLEAUX DE TRAVAIL
C
        IF (IGR.EQ.0) THEN
          NCELW  = NCEL
          NFACW  = NFAC
          NCELW2 = IA(IPGREN(JNCEL,IGR+1))
          DO IEL = 1, NCELW
            IRESPR(IEL) = IA(IPGREN(JIRESP,IGR+1)+IEL-1)
            DAW   (IEL) = DAF(IEL)
          ENDDO
          DO IFAC = 1, NFACW
            IFACEW(1,IFAC) = IFACEL(1,IFAC)
            IFACEW(2,IFAC) = IFACEL(2,IFAC)
            XAW(IFAC,1)    = XAF(IFAC,1)
            IF (ISYM.EQ.2) THEN
              XAW(IFAC,2)  = XAF(IFAC,2)
            ENDIF
          ENDDO
        ELSE
          NCELW  = IA(IPGREN(JNCEL,IGR))
          NFACW  = IA(IPGREN(JNFAC,IGR))
          NCELW2 = IA(IPGREN(JNCEL,IGR+1))
          DO IEL = 1, NCELW
            IRESPR(IEL) = IA(IPGREN(JIRESP,IGR+1)+IEL-1)
            DAW   (IEL) = RA(IPGRRE(JDA,IGR)     +IEL-1)
          ENDDO
          DO IFAC = 1, NFACW
            IFACEW(1,IFAC) = IA(IPGREN(JIFACE,IGR)+2*(IFAC-1)  )
            IFACEW(2,IFAC) = IA(IPGREN(JIFACE,IGR)+2*(IFAC-1)+1)
            XAW(IFAC,1)    = RA(IPGRRE(JXA,IGR)+IFAC-1)
            IF (ISYM.EQ.2) THEN
              XAW(NFACW+IFAC,1)  = RA(IPGRRE(JXA,IGR)+NFACW+IFAC-1)
            ENDIF
          ENDDO
        ENDIF
C
C 1.2 PASSAGE DANS LE LISSEUR
C
        IF (IWARMG.GE.2) THEN
          WRITE(NFECRA,1220)
        ENDIF
C
        NITERF = 0
        CNOMMG(1:4) = 'MAIL'
        WRITE(NUM,'(I4)') IGR
        CNOMMG(5:8) = NUM
C
        IF( IRESMO(IGR).EQ.0 ) THEN
          CALL GRADCO
C         ===========
     &   ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &     ISYM   , IPOL   , NITSMO(1,IGR) , IINVPE ,
     &                                IWARSO , NFECRA , NITERF ,
C                                                       ------
     &     EPSILP , RNORM  , RESIDU ,
C                            ------
     &     IFACEW ,
     &     DAW    , XAW    , BW     , XW     ,
C                                     ------
     &     W1     , W2     , W3     , W4     , W5    )
        ELSEIF ( IRESMO(IGR).EQ.1 ) THEN
          CALL JACOBI
C         ===========
     & (   CNOMMG , NCELW  , NCELW  , NFACW  ,
     &     ISYM   , NITSMO(1,IGR) , IINVPE ,
     &                                IWARSO , NFECRA , NITERF ,
C                                                       ------
     &     EPSILP , RNORM  , RESIDU ,
C                            ------
     &     IFACEW ,
     &     DAW    , XAW    , BW     , XW     , W1     )
C                                              ------
        ELSEIF ( IRESMO(IGR).EQ.2 ) THEN
          CALL CGSTAB
C         ===========
     &   ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &     ISYM   , IPOL   , NITSMO(1,IGR) , IINVPE ,
     &                                IWARSO , NFECRA , NITERF ,
C                                                       ------
     &     EPSILP , RNORM  , RESIDU ,
C                            ------
     &     IFACEW ,
     &     DAW    , XAW    , BW     , XW     ,
C                                     ------
     &     W1     , W2     , W3     , W4     ,
     &     W5     , W6     )
        ENDIF
C
        IF (IGR.EQ.0) NITFIN = NITFIN + NITERF
C
C 1.3 RESTRICTION DU RESIDU
C
        CALL PROMAV(NCELW,NCELW,NFACW,ISYM,IINVPE,IFACEW,
     &              DAW,XAW,XW,W2)
        DO IEL = 1, NCELW
          W2(IEL) = BW(IEL) - W2(IEL)
        ENDDO
        IPROL = 2
        CALL PROLON (NCELW,NCELW,NCELW2,IPROL,NFECRA,IRESPR,W2,BW)
C
C 1.4 TEST DE CONVERGENCE SI DEBUT DU CYCLE (MAILLAGE FIN)
C
      IF (IGR.EQ.0) THEN
C
        ISQRT = 1
        CALL PRODSC(NCELW,NCELW,ISQRT,W2,W2,RESIDU)
        IF (IWARMG.GE.3) THEN
          WRITE( NFECRA,1120) ICYCLE,NITERF,RESIDU/RNORM
        ENDIF
        IF( RESIDU.LE.EPSILP*RNORM) THEN
          IF (IWARMG.GE.1) THEN
            WRITE( NFECRA,1110)
            WRITE( NFECRA,1120) ICYCLE,NITFIN,RESIDU/RNORM
            WRITE( NFECRA,1130)
          ENDIF
          IFIN = 1
          RETURN
        ELSEIF(ICYCLE.GT.NCYMXP) THEN
          IF (IWARMG.GE.1) THEN
            WRITE( NFECRA,1110)
            WRITE( NFECRA,1120) ICYCLE,NITFIN,RESIDU/RNORM
            WRITE( NFECRA,1130)
            WRITE( NFECRA,1140) NCYMXP
          ENDIF
          IFIN = 1
          RETURN
        ENDIF
C
      ENDIF
C
C 1.5 STOCKAGE DU RESIDU ET DE LA CORRECTION
C
        DO IEL = 1, NCELW2
          RA(IPB(IGR+1)+IEL-1) = BW(IEL)
        ENDDO
        IF (IGR.EQ.0) THEN
          DO IEL = 1, NCELW
            XF(IEL) = XW(IEL)
          ENDDO
        ELSE
          DO IEL = 1, NCELW
            RA(IPX(IGR)+IEL-1) = XW(IEL)
          ENDDO
        ENDIF
C
C 1.6 INITIALISATION DE LA CORRECTION
C
        DO IEL = 1, NCELW2
          XW(IEL) = 0.D0
        ENDDO
C
      ENDDO
C
C=======================================================================
C 2.  RESOLUTION JUSQU'A CONVERGENCE SUR NIVEAU LE + GROSSIER : NGR
C=======================================================================
C
      IF (IWARMG.GE.2) THEN
        WRITE(NFECRA,1230)
      ENDIF
C
      NCELW  = IA(IPGREN(JNCEL,NGR))
      NFACW  = IA(IPGREN(JNFAC,NGR))
      DO IEL = 1, NCELW
        DAW(IEL) = RA(IPGRRE(JDA,NGR)+IEL-1)
      ENDDO
      DO IFAC = 1, NFACW
        IFACEW(1,IFAC) = IA(IPGREN(JIFACE,NGR)+2*(IFAC-1)  )
        IFACEW(2,IFAC) = IA(IPGREN(JIFACE,NGR)+2*(IFAC-1)+1)
        XAW(IFAC,1)    = RA(IPGRRE(JXA,NGR)+IFAC-1)
        IF (ISYM.EQ.2) THEN
          XAW(NFACW+IFAC,1)  = RA(IPGRRE(JXA,NGR)+NFACW+IFAC-1)
        ENDIF
      ENDDO
C
      NITERF = 0
      CNOMMG(1:4) = 'MAIL'
      WRITE(NUM,'(I4)') NGR
      CNOMMG(5:8) = NUM
C
      IF( IRESMO(NGR).EQ.0 ) THEN
        CALL GRADCO
C       ===========
     & ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &   ISYM   , IPOL   , NITSMO(1,NGR) , IINVPE ,
     &                              IWARSO , NFECRA , NITERF ,
C                                                     ------
     &   EPSPGR , RNORM  , RESIDU ,
C                          ------
     &   IFACEW ,
     &   DAW    , XAW    , BW     , XW     ,
C                                   ------
     &   W1     , W2     , W3     , W4     , W5     )
      ELSEIF ( IRESMO(IGR).EQ.1 ) THEN
        CALL JACOBI
C       ===========
     & ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &   ISYM   , NITSMO(1,NGR) , IINVPE ,
     &                              IWARSO , NFECRA , NITERF ,
C                                                     ------
     &   EPSPGR , RNORM  , RESIDU ,
C                          ------
     &   IFACEW ,
     &   DAW    , XAW    , BW     , XW     , W1     )
C                                            ------
      ELSEIF ( IRESMO(IGR).EQ.2 ) THEN
        CALL CGSTAB
C       ===========
     & ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &   ISYM   , IPOL   , NITSMO(1,NGR) , IINVPE ,
     &                              IWARSO , NFECRA , NITERF ,
C                                                     ------
     &   EPSPGR , RNORM  , RESIDU ,
C                          ------
     &   IFACEW ,
     &   DAW    , XAW    , BW     , XW     ,
C                                   ------
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     )
      ENDIF
C
C STOCKAGE DE LA CORRECTION
C
      DO IEL = 1, NCELW
        RA(IPX(NGR)+IEL-1) = XW(IEL)
      ENDDO
C
C=======================================================================
C 3.  BOUCLE DE REMONTEE
C=======================================================================
C
      IF (IWARMG.GE.2) THEN
        WRITE(NFECRA,1210)
      ENDIF
C
      DO 100 IGR = NGR-1, 0, -1
C
C 3.1  ON PASSE LES DESCRIPTEURS DE FICHIER DANS DES TABLEAUX DE TRAVAIL
C
        IF (IGR.EQ.0) THEN
          NCELW  = NCEL
          NFACW  = NFAC
          NCELW2 = IA(IPGREN(JNCEL,IGR+1))
          DO IEL = 1, NCELW
            IRESPR(IEL) = IA(IPGREN(JIRESP,IGR+1)+IEL-1)
            DAW   (IEL) = DAF(IEL)
          ENDDO
          DO IFAC = 1, NFACW
            IFACEW(1,IFAC) = IFACEL(1,IFAC)
            IFACEW(2,IFAC) = IFACEL(2,IFAC)
            XAW(IFAC,1)    = XAF(IFAC,1)
            IF (ISYM.EQ.2) THEN
              XAW(IFAC,2)  = XAF(IFAC,2)
            ENDIF
          ENDDO
        ELSE
          NCELW  = IA(IPGREN(JNCEL,IGR))
          NFACW  = IA(IPGREN(JNFAC,IGR))
          NCELW2 = IA(IPGREN(JNCEL,IGR+1))
          DO IEL = 1, NCELW
            IRESPR(IEL) = IA(IPGREN(JIRESP,IGR+1)+IEL-1)
            DAW   (IEL) = RA(IPGRRE(JDA,IGR)     +IEL-1)
          ENDDO
          DO IFAC = 1, NFACW
            IFACEW(1,IFAC) = IA(IPGREN(JIFACE,IGR)+2*(IFAC-1)  )
            IFACEW(2,IFAC) = IA(IPGREN(JIFACE,IGR)+2*(IFAC-1)+1)
            XAW(IFAC,1)    = RA(IPGRRE(JXA,IGR)+IFAC-1)
            IF (ISYM.EQ.2) THEN
              XAW(NFACW+IFAC,1)  = RA(IPGRRE(JXA,IGR)+NFACW+IFAC-1)
            ENDIF
          ENDDO
        ENDIF
C
C 3.2 PROLONGATION DE LA CORRECTION
C
        IPROL = 1
        CALL PROLON (NCELW,NCELW,NCELW2,IPROL,NFECRA,IRESPR,W2,XW)
        IF (IGR.EQ.0) THEN
          DO IEL = 1, NCELW
            XW(IEL) = XF(IEL) + W2(IEL)
          ENDDO
        ELSE
          DO IEL = 1, NCELW
            XW(IEL) = RA(IPX(IGR)+IEL-1) + W2(IEL)
          ENDDO
        ENDIF
C
C 3.3 PASSAGE DANS LE LISSEUR
C
C     RECUPERATION DU RESIDU DANS BW
C
        IF (IGR.EQ.0) THEN
          DO IEL = 1, NCELW
            BW(IEL) = BF(IEL)
          ENDDO
        ELSE
          DO IEL = 1, NCELW
            BW(IEL) = RA(IPB(IGR)+IEL-1)
          ENDDO
        ENDIF
C
        IF (IWARMG.GE.2) THEN
          WRITE(NFECRA,1220)
        ENDIF
C
        NITERF = 0
        CNOMMG(1:4) = 'MAIL'
        WRITE(NUM,'(I4)') IGR
        CNOMMG(5:8) = NUM
C
        IF( IRESMO(IGR).EQ.0 ) THEN
          CALL GRADCO
C         ===========
     &   ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &     ISYM   , IPOL   , NITSMO(2,IGR) ,  IINVPE ,
     &                                IWARSO , NFECRA , NITERF ,
C                                                       ------
     &     EPSILP , RNORM  , RESIDU ,
C                            ------
     &     IFACEW ,
     &     DAW    , XAW    , BW     , XW     ,
C                                     ------
     &     W1     , W2     , W3     , W4     , W5     )
         ELSEIF ( IRESMO(IGR).EQ.1 ) THEN
          CALL JACOBI
C         ===========
     & (   CNOMMG , NCELW  , NCELW  , NFACW  ,
     &     ISYM   , NITSMO(2,IGR) , IINVPE ,
     &                                IWARSO , NFECRA , NITERF ,
C                                                       ------
     &     EPSILP , RNORM  , RESIDU ,
C                            ------
     &     IFACEW ,
     &     DAW    , XAW    , BW     , XW     , W1     )
C                                              ------
         ELSEIF ( IRESMO(IGR).EQ.2 ) THEN
          CALL CGSTAB
C         ===========
     &   ( CNOMMG , NCELW  , NCELW  , NFACW  ,
     &     ISYM   , IPOL   , NITSMO(2,IGR) , IINVPE ,
     &                                IWARSO , NFECRA , NITERF ,
C                                                       ------
     &     EPSILP , RNORM  , RESIDU ,
C                            ------
     &     IFACEW ,
     &     DAW    , XAW    , BW     , XW     ,
C                                     ------
     &     W1     , W2     , W3     , W4     ,
     &     W5     , W6     )
        ENDIF
C
        IF (IGR.EQ.0) NITFIN = NITFIN + NITERF
C
 100  CONTINUE
C
C=======================================================================
C 4.  FINALISATIONS
C=======================================================================
C
      DO IEL = 1, NCEL
        XF(IEL) = XW(IEL)
      ENDDO
C
C--------
C FORMATS
C--------
C
 1110  FORMAT
     &(1X,
     &'-------------------------------------------------'/,1X,
     &' NB DE | CUMUL DES ITERATIONS SUR | RESIDU NORME',/,
     & 1X,
     &'CYCLES | LE MAILLAGE FIN          | MAILLAGE FIN',/,
     & 1X,
     &'-------|--------------------------|--------------')
C
 1120 FORMAT(2X,I4,'  |',9X,I5,12X,'|  ',E11.4)
 1130 FORMAT(1X,
     &'-------------------------------------------------')
 1140 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' MULTIGRILLE ALGEBRIQUE              ',/,
     &'@    *********                                               ',/,
     &'@  Nombre de cycles     maximal ',I10   ,' atteint           ',/,
     &'@                                                            '  )
 1200 FORMAT(' CYCMGR : DESCENTE')
 1210 FORMAT(' CYCMGR : REMONTEE')
 1220 FORMAT('   PASSAGE DANS LE LISSEUR')
 1230 FORMAT('   RESOLUTION SUR NIVEAU LE + GROSSIER')
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
