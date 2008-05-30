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
                        SUBROUTINE CRSTGR
C                       *****************
C     -------------------------------------------------------------
     &( IAPPEL , ISYM   , IGR    ,
     &  NCELF  , NCELG  , NCELFE , NCELGE , NFACF  , NFACG ,
     &  IWARNP ,
     &  IFACLF , IFACLG , IRSCEL , IRSFAC ,
     &  VOLUMF , XYZFIN , SURFAF , XAF0   , XAF0IJ ,
     &  DAF    , XAF    ,
     &  VOLUMG , XYZGRO , SURFAG , XAG0   , XAG0IJ ,
     &  DAG    , XAG    ,
     &  W1     , W2     , W3     , W4     )
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC   MULTIGRILLE ALGEBRIQUE :
CFONC   CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER A PARTIR
CFONC   DU NIVEAU SUPERIEUR
CFONC
CFONC  structure grille fine ==> grille grossiere
CFONC  ------------------------------------------
CFONC  [ NCELF,NFACF,IFACLF,ROVDTF,XAF0,VOLUMF,XYZFIN,
CFONC         SURFAF,XAF,XAF0,XAF0IJ ]
CFONC  [ NCELG,NFACG,IFACLG,ROVDTG,XAG0,VOLUMG,XYZGRO,
CFONC         SURFAG,XAG,XAG0,XAG0IJ ]
CFONC
CFONC              explicité dans le cartouche
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IAPPEL       ! E  !  -> ! NUMERO D'APPEL                       !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! IGR          ! E  !  -> ! NIVEAU DU MAILLAGE GROSSIER          !
CARGU ! NCELF        ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE FIN       !
CARGU ! NCELG        ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE GROSSIER  !
CARGU ! NCELFE       ! E  !  -> ! NOMBRE D'ELEMENTS ETENDUS FIN        !
CARGU ! NCELGE       ! E  !  -> ! NOMBRE D'ELEMENTS ETENDUS GROSSIER   !
CARGU ! NFACF        ! E  !  -> ! NOMBRE DE FACES INTERNES MAILL. FIN  !
CARGU ! NFACG        ! E  !  -> ! NOMBRE DE FACES INTERNES MAILL. GRO. !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! IFACEF       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU !  (2, NFACF)  !    !     ! SUR MAILLAGE FIN                     !
CARGU ! IFACEG       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU !  (2, NFACG)  !    !     ! SUR MAILLAGE GROSSIER                !
CARGU ! IRSCEL       ! TE ! <-  ! CELLULE FINE -> CELLULE GROSSIERE    !
CARGU !  (NCELFE)    !    !     !                                      !
CARGU ! IRSFAC       ! TE !  -> ! FACE FINE -> FACE GROSSIERE          !
CARGU !  (NFACF)     !    !     !  = 0 : FACE INTERNE CEL. GROSSIERE   !
CARGU !              !    !     !  < 0 : ORIENTATION INVERSE           !
CARGU !              !    !     !  > 0 : ORIENTATION IDENTIQUE         !
CARGU ! DAF(NCELF)   ! TR !  -> ! DIAGONALE DE LA MATRICE MAIL FIN     !
CARGU ! XAF          ! TR ! <-  ! EXTRADIAGONALE MATRICE MAILLAGE FIN  !
CARGU !  (NFACF,ISYM)!    !     !                                      !
CARGU ! DAG(NCELG)   ! TR ! <-  ! DIAGONALE MATRICE MAILLAGE GROSSIER  !
CARGU ! XAG          ! TR ! <-  ! EXTRADIAGONALE MATRICE MAILLAGE      !
CARGU !  (NFACG,ISYM)!    !     !  GROSSIER                            !
CARGU ! IVOIS(NCELF) ! TE !  -  ! INDICATEUR DE VOISINAGE CELLULES     !
CARGU ! IP(NCELF)    ! TE !  -  ! POINTEURS SUR VOISINS POUR CONNECTI- !
CARGU !              !    !     ! VITE INVERSE                         !
CARGU ! ICELFA       ! TE !  -  ! CONNECTIVITE CELLULES->FACES MAILLA- !
CARGU !  (2*NFACF)   !    !     ! GE FIN                               !
CARGU ! ICELCE       ! TE !  -  ! CONNECTIVITE CELLULES->CELLULES      !
CARGU !  (2*NFACF)   !    !     ! VOISINES DU MAILLAGE FIN             !
CARGU ! RW(NCELF)    ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !IFACLF(2,NFACF! TE !  -> ! CELL. VOISINES FACE INTRN MAILL FIN  !
CARGU !XAF0(NFACF,ISYM TR !  -> ! EXTRADIAGONALE MATRICE P0 MAILAGE FIN!
CARGU ! VOLUMF(NCELF)! TR !  -> ! VOLUME CELLULE MAILLAGE FIN          !
CARGU !XYZFIN(3,NCELF) TR !  -> ! COORDONNES CELLULE MAILLAGE FIN      !
CARGU !SURFAF(3,NFACF) TR !  -> ! SURFACE FACE INTERNE MAILLAGE FIN    !
CARGU !XAFXF0(2,NFACF) TR !  -> ! INTEG. XAF0*COORD.CELL ADJ. MAIL.FIN !
CARGU ! NCELG        ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE GROSSIER  !
CARGU ! NFACG        ! E  !  <- ! NOMBRE FACES INTERNES MAILL. GROSSIER!
CARGU !IFACLG(2,NFACG! TE !  <- ! CELL. VOISINES FACE INTRN MAIL. GROS !
CARGU !XAG0(NFACG,ISYM TR !  <- ! EXTRADIAGONALE MATRICE P0 MAILL.GROS.!
CARGU ! VOLUMG(NCELG)! TR !  <- ! VOLUME CELLULE MAILLAGE GROSSIER     !
CARGU !XYZGRO(3,NCELG) TR !  <- ! COORDONNES CELLULE MAILLAGE GROSSIER !
CARGU !SURFAG(3,NFACG) TR !  <- ! SURFACE FACE INTERNE MAILL. GROSSIER !
CARGU !XAGXG0(2,NFACG) TR !  <- ! INTEG. XAG0*COORD.CELL ADJ. MAIL.GRO !
cargu !W1,..,4 (NCEL)! TR ! <-> ! TABLEAUX DE TRAVAILS
CARGU !______________!____!_____!______________________________________!
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
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IAPPEL, ISYM, IGR
      INTEGER          NCELF, NCELFE, NFACF, NCELG, NCELGE, NFACG
      INTEGER          IWARNP
C
      DOUBLE PRECISION DAF(NCELF), XAF(NFACF,2)
      DOUBLE PRECISION DAG(NCELGE), XAG(NFACG,2)
C
      INTEGER          IFACLF(2, NFACF), IFACLG(2, NFACG)
      INTEGER          IRSCEL(NCELF), IRSFAC(NFACF)
C
      DOUBLE PRECISION XAF0(NFACF),VOLUMF(NCELFE), XYZFIN(3, NCELFE)
      DOUBLE PRECISION XAG0(NFACG),VOLUMG(NCELGE), XYZGRO(3, NCELGE)
      DOUBLE PRECISION SURFAF(3, NFACF), XAF0IJ(3, NFACF)
      DOUBLE PRECISION SURFAG(3, NFACG), XAG0IJ(3, NFACG)
      DOUBLE PRECISION W1(NCELFE), W2(NCELFE), W3(NCELFE), W4(NCELFE)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, II, JJ, IG, JG, IIG, JJG
      INTEGER          IFACG, IMIN, IMAX
      INTEGER          IFAC, IFAC1, IFAC2, N2FAC, IFACI, NFACI
      INTEGER          INTERP, IIMIN
C
      DOUBLE PRECISION DSIGJG, DSXAIJ
      DOUBLE PRECISION RMIN, RMAX, SX, SY, SZ, CCLIP
      DOUBLE PRECISION ANMIN(2), ANMAX(2)
C
C
C***********************************************************************
C
C CREATION VOLUME, CdG CELLULES GROSSIERES : XYZGRO(3,NCELG)
C==================================================================
C
      IF (IAPPEL .EQ. 1) THEN
C
        DO IG = 1, NCELGE
          VOLUMG(IG) = 0.D0
        ENDDO
        DO IG = 1, NCELGE
          XYZGRO(1, IG) = 0.D0
          XYZGRO(2, IG) = 0.D0
          XYZGRO(3, IG) = 0.D0
        ENDDO
C
        DO II = 1, NCELF
          IG = IRSCEL(II)
          VOLUMG(IG) = VOLUMG(IG) + VOLUMF(II)
          XYZGRO(1, IG) = XYZGRO(1, IG) + VOLUMF(II)*XYZFIN(1, II)
          XYZGRO(2, IG) = XYZGRO(2, IG) + VOLUMF(II)*XYZFIN(2, II)
          XYZGRO(3, IG) = XYZGRO(3, IG) + VOLUMF(II)*XYZFIN(3, II)
        ENDDO
        DO IG = 1, NCELG
          XYZGRO(1, IG) = XYZGRO(1, IG) / VOLUMG(IG)
          XYZGRO(2, IG) = XYZGRO(2, IG) / VOLUMG(IG)
          XYZGRO(3, IG) = XYZGRO(3, IG) / VOLUMG(IG)
        ENDDO
C
C       RETOUR A LA FONCTION APPELANTE POUR SYNCHRONISATION
C       PARALLELE / PERIODIQUE DE XYZGRO ET VOLUMG
C
        RETURN
C
      ENDIF
C
C RESTRICTION P0 MATRICES, SURFACE "INTERNE" :
C       XAG0(NFACG), SURFAG(3,NFACG), XAGXG0(2,NFACG)
C==================================================================
C
      IMAX = 0

      DO IFACG = 1, NFACG
        XAG0(IFACG) = 0.D0
        SURFAG(1, IFACG) = 0.D0
        SURFAG(2, IFACG) = 0.D0
        SURFAG(3, IFACG) = 0.D0
        XAG0IJ(1, IFACG) = 0.D0
        XAG0IJ(2, IFACG) = 0.D0
        XAG0IJ(3, IFACG) = 0.D0
      ENDDO
C
      DO IFAC = 1, NFACF
C
        IF (IRSFAC(IFAC).GT.0) THEN
C
          IFACG = IRSFAC(IFAC)

          XAG0(IFACG) = XAG0(IFACG) + XAF0(IFAC)

          SURFAG(1, IFACG) = SURFAG(1, IFACG) + SURFAF(1, IFAC)
          SURFAG(2, IFACG) = SURFAG(2, IFACG) + SURFAF(2, IFAC)
          SURFAG(3, IFACG) = SURFAG(3, IFACG) + SURFAF(3, IFAC)
          XAG0IJ(1, IFACG) = XAG0IJ(1, IFACG) + XAF0IJ(1, IFAC)
          XAG0IJ(2, IFACG) = XAG0IJ(2, IFACG) + XAF0IJ(2, IFAC)
          XAG0IJ(3, IFACG) = XAG0IJ(3, IFACG) + XAF0IJ(3, IFAC)

        ELSE IF (IRSFAC(IFAC).LT.0) THEN

          IFACG = - IRSFAC(IFAC)

          XAG0(IFACG) = XAG0(IFACG) + XAF0(IFAC)

          SURFAG(1, IFACG) = SURFAG(1, IFACG) - SURFAF(1, IFAC)
          SURFAG(2, IFACG) = SURFAG(2, IFACG) - SURFAF(2, IFAC)
          SURFAG(3, IFACG) = SURFAG(3, IFACG) - SURFAF(3, IFAC)
          XAG0IJ(1, IFACG) = XAG0IJ(1, IFACG) - XAF0IJ(1, IFAC)
          XAG0IJ(2, IFACG) = XAG0IJ(2, IFACG) - XAF0IJ(2, IFAC)
          XAG0IJ(3, IFACG) = XAG0IJ(3, IFACG) - XAF0IJ(3, IFAC)

        ENDIF
C
      ENDDO
C
C
C=======================================================================
C FINALISATION CALCUL DE LA MATRICE  DANS DAG, XAG
C=======================================================================
C
C INTERP= 0 : RESTRICTION P0 /PROLONGATION P0 => XAG=XAG0
C INTERP= 1 : RESTRICTION P0 /PROLONGATION P1 => XAG=XAG0IJ/IgJg
C
C     INITIALISATION
C
      INTERP = 0
      INTERP = 1
C
C
C     INITIALISATION TERME NON DIFFERENTIEL MESH FIN STOCKE DANS W1
C======================
C
      DO IEL = 1, NCELF
        W1(IEL) = DAF(IEL)
      ENDDO
      DO IEL = NCELF + 1, NCELFE
        W1(IEL) = 0.D0
      ENDDO
      DO IFAC = 1, NFACF
        II = IFACLF(1, IFAC)
        JJ = IFACLF(2, IFAC)
        W1(II) = W1(II) + XAF(IFAC, 1)
        W1(JJ) = W1(JJ) + XAF(IFAC, ISYM)
      ENDDO
C
C
C     INITIALISATION STOCKAGE MATRICE MESH GROSSIER SUR (DAG, XAG)
C=====================
C
      DO IEL = 1, NCELGE
        DAG(IEL) = 0.D0
      ENDDO
      DO IFAC = 1, NFACG
        XAG(IFAC, 1)= 0.D0
        XAG(IFAC, ISYM)= 0.D0
      ENDDO
C
C
C     TERMES EXTRADIAGONAUX
C     (matrices sym. pour l'instant, meme si stockage non syme isym=2)
C================================
C
C     MATRICE INTIALISEE A XAG0 (INTERP=0)
C
      DO IFACG = 1, NFACG
        XAG(IFACG, 1) = XAG0(IFACG)
        XAG(IFACG, ISYM) = XAG0(IFACG)
      ENDDO
C
      IF (INTERP.EQ.1) THEN
C
        IMIN = 0
        IMAX = 0
C
        DO IFACG = 1, NFACG
C
          IG = IFACLG(1, IFACG)
          JG = IFACLG(2, IFACG)

          DSIGJG =   (XYZGRO(1, JG)-XYZGRO(1, IG))*SURFAG(1, IFACG)
     &             + (XYZGRO(2, JG)-XYZGRO(2, IG))*SURFAG(2, IFACG)
     &             + (XYZGRO(3, JG)-XYZGRO(3, IG))*SURFAG(3, IFACG)
C
          DSXAIJ =   XAG0IJ(1, IFACG)*SURFAG(1, IFACG)
     &             + XAG0IJ(2, IFACG)*SURFAG(2, IFACG)
     &             + XAG0IJ(3, IFACG)*SURFAG(3, IFACG)
C
          IF (ABS(DSIGJG) .GT. EPZERO) THEN
C
C           STANDARD
            XAG(IFACG, 1)    = DSXAIJ/DSIGJG
            XAG(IFACG, ISYM) = DSXAIJ/DSIGJG

C           MATRICE CLIPPEE
            CCLIP= DSXAIJ/DSIGJG
            IF (CCLIP .LT. XAG0(IFACG)) IMIN = IMIN+1
            IF (CCLIP .GT. 0.D0)  IMAX = IMAX +1
            IF (CCLIP .LT. XAG0(IFACG) .OR. CCLIP .GT. 0.D0) THEN
              XAG(IFACG, 1) =  XAG0(IFACG)
              XAG(IFACG, ISYM) = XAG0(IFACG)
            ENDIF
C
          ENDIF
C
        ENDDO
C
        IF(IWARNP.GT.3) THEN
          IF (IRANGP. GE. 0) THEN
            CALL PARCPT(IMIN)
            CALL PARCPT(IMAX)
          ENDIF
          WRITE(NFECRA, *)
     &      '    crstgr.F : matrice grossiere < XAG0 en ',IMIN,' faces'
          WRITE(NFECRA, *)
     &      '    crstgr.F : matrice grossiere > 0    en ',IMAX,' faces'
        ENDIF
C
C       RELAXATION EVENTUELLE MATRICE P1 / MATRICE P0 DEFINI PAR UTILISATEUR
C       DANS usini1.F
C
        DO IFACG = 1, NFACG
          XAG(IFACG, 1)
     &      = RLXP1*XAG(IFACG, 1) +(1.D0-RLXP1)*XAG0(IFACG)
          XAG(IFACG, ISYM)
     &      = RLXP1*XAG(IFACG, ISYM) +(1.D0-RLXP1)*XAG0(IFACG)
        ENDDO
C
      ENDIF
C
      IF (INTERP.NE.0 .AND. INTERP.NE.1) THEN
C
        WRITE(NFECRA,*) 'INTERP MAL DEFINI DANS crstgr.F'
        WRITE(NFECRA,*) '--> ARRET DANS crstgr.F '
        CALL CSEXIT(1)
C
      ENDIF
C
C
C     TERME DIAGONAL
C====================
C
      DO II = 1, NCELF
        IG = IRSCEL(II)
        DAG(IG) = DAG(IG) + W1(II)
      ENDDO
C
      DO IFACG = 1, NFACG
C
        IG = IFACLG(1, IFACG)
        JG = IFACLG(2, IFACG)
C
        DAG(IG) = DAG(IG) - XAG(IFACG, 1)
        DAG(JG) = DAG(JG) - XAG(IFACG, ISYM)
C
      ENDDO
C
C     CONTROLE
C=============
c     WRITE(NFECRA,*) 'TYPE INTERPOLATION MATRICE = ',INTERP
C
C     EVALUATION ANISOTROPIE DES MATRICES FINE ET GROSSIERE
C
      IF (IWARNP .GT. 3) THEN
C
        DO II = 1, NCELFE
          W1(II) =-1.D12
          W2(II) =+1.D12
          W3(II) =-1.D12
          W4(II) =+1.D12
        ENDDO
C
        DO IFAC = 1, NFACF
          II = IFACLF(1, IFAC)
          JJ = IFACLF(2, IFAC)
          W1(II) = MAX(ABS(XAF(IFAC, 1)), W1(II))
          W2(II) = MIN(ABS(XAF(IFAC, 1)), W2(II))
          W1(JJ) = MAX(ABS(XAF(IFAC, ISYM)), W1(JJ))
          W2(JJ) = MIN(ABS(XAF(IFAC, ISYM)), W2(JJ))
        ENDDO
C
        DO IFACG = 1, NFACG
          IG = IFACLG(1, IFACG)
          JG = IFACLG(2, IFACG)
          W3(IG) = MAX(ABS(XAG(IFACG,1)), W3(IG))
          W4(IG) = MIN(ABS(XAG(IFACG,1)), W4(IG))
          W3(JG) = MAX(ABS(XAG(IFACG, ISYM)), W3(JG))
          W4(JG) = MIN(ABS(XAG(IFACG, ISYM)), W4(JG))
        ENDDO
C
        DO II = 1, NCELF
          W1(II) = W2(II)/W1(II)
        ENDDO
C
        DO IG = 1, NCELG
          W3(IG) = W4(IG)/W3(IG)
        ENDDO
C
        ANMIN(1) = W1(1)
        ANMAX(1) = W1(1)
        DO II = 2, NCELF
          IF (W1(II) .LT. ANMIN(1)) THEN
            ANMIN(1) = W1(II)
          ELSE IF (W1(II) .GT. ANMAX(1)) THEN
            ANMAX(1) = W1(II)
          ENDIF
        ENDDO
C
        ANMIN(2) = W3(1)
        ANMAX(2) = W3(1)
        DO IG = 2, NCELG
          IF (W3(IG) .LT. ANMIN(2)) THEN
            ANMIN(2) = W3(IG)
          ELSE IF (W3(IG) .GT. ANMAX(2)) THEN
            ANMAX(2) = W3(IG)
          ENDIF
        ENDDO
C
        IF (IRANGP. GE. 0) THEN
          II = 2
          CALL PARRMN(II, ANMIN)
          CALL PARRMX(II, ANMAX)
        ENDIF
C
        IF (IWARNP .GT. 3) THEN
C
          WRITE (NFECRA, 2000) ANMIN(1), ANMAX(1), ANMIN(2), ANMAX(2)
C
          IF (INTERP .EQ. 1) THEN
C
            RMIN = +1.D10
            RMAX = -1.D10
            DO IFACG=1,NFACG
              RMIN = MIN(RMIN, XAG(IFACG,1)/XAG0(IFACG))
              RMAX = MAX(RMAX, XAG(IFACG,1)/XAG0(IFACG))
            ENDDO
C
            IF (IRANGP. GE. 0) THEN
              CALL PARMIN(RMIN)
              CALL PARMAX(RMAX)
            ENDIF

            WRITE(NFECRA, *) '      minimum XAG_P1/XAG_P0 = ', RMIN
            WRITE(NFECRA, *) '      maximum XAG_P1/XAG_P0 = ', RMAX
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 2000 FORMAT('       anisotropie maillage fin :      min = ', E12.5, /,
     &       '                                       max = ', E12.5, /,
     &       '       anisotropie maillage grossier : min = ', E12.5, /,
     &       '                                       max = ', E12.5, /)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
