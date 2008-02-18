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
                        SUBROUTINE CALGRG
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   IDBIGR , IDBRGR ,
     &   NFECRA ,
     &   ISYM   , IGR    , NTGREN , NTGRRE ,
     &   JNCEL  , JNFAC  , JIFACE , JIRESP , JDA    , JXA    ,
     &   NCELF  , NFACF  ,
     &   IFACEF ,
     &   DAF    , XAF    ,
     &   NCELG  , NFACG  ,
     &   JPGREN , JPGRRE ,
     &   IRESPR , IFACEG ,
     &   DAG    , XAG    ,
     &   IVOIS  , IP     , IW     ,
     &   INDICF , IRESFA , ICELFA , ICELCE , RW )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC        MULTIGRILLE ALGEBRIQUE :
CFONC   CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER A PARTIR
CFONC   DU NIVEAU SUPERIEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IDBIGR       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU !              !    !     !                         APRES CALGRG !
CARGU ! IDBRGR       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU !              !    !     !                         APRES CALGRG !
CARGU ! NFECRA       ! E  !  -> ! NUMERO DU FICHIER D'IMPRESSIONS      !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! IGR          ! E  !  -> !NIVEAU DU MAILLAGE GROSSIER           !
CARGU ! NTGREN       ! E  !  -> ! DIMENSION DU TABLEAUX DE POINTEURS   !
CARGU !              !    !     !  ENTIERS SUR LES MAILLAGES GROSSIERS !
CARGU ! NTGRRE       ! E  !  -> ! DIMENSION DU TABLEAUX DE POINTEURS   !
CARGU !              !    !     !  REELS SUR LES MAILLAGES GROSSIERS   !
CARGU ! JNCEL        ! E  !  -> ! ORDRE DE STOCKAGE DE NCELG           !
CARGU !              !    !     !  DANS JPGREN                         !
CARGU ! JNFAC        ! E  !  -> ! ORDRE DE STOCKAGE DE NFACG           !
CARGU !              !    !     !  DANS JPGREN                         !
CARGU ! JIFACE       ! E  !  -> ! ORDRE DE STOCKAGE DE IFACEG          !
CARGU !              !    !     !  DANS JPGREN                         !
CARGU ! JIRESP       ! E  !  -> ! ORDRE DE STOCKAGE DE IRESPR          !
CARGU !              !    !     !  DANS JPGREN                         !
CARGU ! JDA          ! E  !  -> ! ORDRE DE STOCKAGE DE JDA             !
CARGU !              !    !     !  DANS JPGRRE                         !
CARGU ! JXA          ! E  !  -> ! ORDRE DE STOCKAGE DE JXA             !
CARGU !              !    !     !  DANS JPGRRE                         !
CARGU ! NCELF        ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE FIN       !
CARGU ! NFACF        ! E  !  -> ! NOMBRE DE FACES INTERNES MAILL. FIN  !
CARGU ! IFACEF       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFACF)   !    !     ! SUR MAIL FIN                         !
CARGU ! DAF(NCELF)   ! TR !  -> ! DIAGONALE DE LA MATRICE MAIL FIN     !
CARGU !XAF(NFACF,ISYM! TR !  -> ! EXTRADIAGONALE DE LA MATRICE MAIL FIN!
CARGU ! NCELG        ! E  ! <-  ! NOMBRE D'ELEMENTS MAILLAGE GROSSIER  !
CARGU ! NFACG        ! E  ! <-  ! NOMBRE DE FACES INTERNES             !
CARGU !              !    !     !  MAILLAGE GROSSIER                   !
CARGU ! JPGREN       ! TE ! <-  ! POINTEURS SUR DESCRIPTEURS ENTIERS   !
CARGU !  (NTGREN)    !    !     !  DU MAILLAGE GROSSIER                !
CARGU ! JPGRRE       ! TE ! <-  ! POINTEURS SUR DESCRIPTEURS REELS     !
CARGU !  (NTGRRE)    !    !     !  DU MAILLAGE GROSSIER                !
CARGU ! IRESPR       ! TE ! <-  ! DESCRIPTION DES CELLULES GROSSIERES  !
CARGU !  (NCELF)     !    !     !                                      !
CARGU ! IFACEG       ! TE ! <-  ! DESCRIPTION DES FACES GROSSIERES     !
CARGU !  (2,NFACF)   !    !     !                                      !
CARGU ! DAG(NCELF)   ! TR ! <-  ! DIAGONALE MATRICE MAILLAGE GROSSIER  !
CARGU ! XAG          ! TR ! <-  ! EXTRADIAGONALE MATRICE MAILLAGE      !
CARGU !  (NFACF,ISYM)!    !     !  GROSSIER                            !
CARGU ! IVOIS(NCELF) ! TE !  -  ! INDICATEUR DE VOISINAGE CELLULES     !
CARGU ! IP(NCELF)    ! TE !  -  ! POINTEURS SUR VOISINS POUR CONNECTI- !
CARGU !              !    !     ! VITE INVERSE                         !
CARGU ! IW(NCELF)    ! TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! INDICF(NFACF)! TE !  -  ! INDICATEUR DE REGROUPEMENT DES FACES !
CARGU ! ICELFA       ! TE !  -  ! CONNCETIVITE CELLULES->FACES MAILLA- !
CARGU !  (2*NFACF)   !    !     ! GE FIN                               !
CARGU ! ICELCE       ! TE !  -  ! CONNECTIVITE CELLULES->CELLULES      !
CARGU !  (2*NFACF)   !    !     ! VOISINES DU MAILLAGE FIN             !
CARGU ! RW(NCELF)    ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          IDBIGR , IDBRGR
      INTEGER          NFECRA
      INTEGER          ISYM , IGR , NTGREN , NTGRRE
      INTEGER          JNCEL , JNFAC , JIFACE , JIRESP, JDA , JXA
      INTEGER          NCELF , NFACF
      INTEGER          NCELG , NFACG
C
      INTEGER          JPGREN(NTGREN) , JPGRRE(NTGRRE)
      INTEGER          IFACEF(2,NFACF) , IFACEG(2,NFACF)
      INTEGER          IRESPR(NCELF)
      INTEGER          IVOIS(NCELF) , IP(NCELF) , IW(NCELF)
      INTEGER          INDICF(NFACF)
      INTEGER          IRESFA(NFACF),ICELFA(2*NFACF) , ICELCE(2*NFACF)
C

      DOUBLE PRECISION DAF(NCELF), XAF(NFACF,ISYM)
      DOUBLE PRECISION DAG(NCELF), XAG(NFACF,ISYM)
      DOUBLE PRECISION RW(NCELF)
C
C
C VARIABLES LOCALES
C
      LOGICAL          CRITER
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL , IEL1 , IEL2 , IEL1V ,  IFAC , IFAC2
      INTEGER          IELG , IEL1G , IEL2G, IFACG , ICELG , NCELT
      INTEGER          ICOMP , IPASS , NPASS ,IDEBI , IFINI , IVOI , JJ
      DOUBLE PRECISION COCRIA , COCRIB , RTEST1 , RTEST2

C
C***********************************************************************
C
C ATTENTION :
C 1- le cas des faces periodiques n'est pas encore
C             pris en compte
C 2- le cas des faces redecoupees n'est pas encore
C             pris en compte
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1.  CREATION DES CELLULES GROSSIERES : IRESPR
C=======================================================================
C
C INITIALISATION
C
      ICELG = 0
      NPASS = 1
C CRITERE DEBUT
      COCRIA = 0.08D0*(0.5D0)**(IGR-1)
      COCRIB = 0.95D0
C CRITERE FIN
      DO IEL = 1, NCELF
       IW    (IEL) = 0
       IRESPR(IEL) = 0
       IVOIS (IEL) = 0
       RW    (IEL) = -1D12
      ENDDO
C CALCUL DES MAX PAR LIGNE : MAX (|Aij|,|Aji|) ,i<>j
      IF (ISYM.EQ.1) THEN
        DO IFAC = 1, NFACF
          IEL1 = IFACEF(1,IFAC)
          IEL2 = IFACEF(2,IFAC)
          RW(IEL1) = MAX(RW(IEL1),ABS(XAF(IFAC,1)))
          RW(IEL2) = MAX(RW(IEL2),ABS(XAF(IFAC,1)))
        ENDDO
      ELSE
        DO IFAC = 1, NFACF
          IEL1 = IFACEF(1,IFAC)
          IEL2 = IFACEF(2,IFAC)
          RW(IEL1) = MAX(RW(IEL1),ABS(XAF(IFAC,1)))
          RW(IEL2) = MAX(RW(IEL2),ABS(XAF(IFAC,2)))
        ENDDO
      ENDIF
      IF (IRANGP.GE.0) THEN
        CALL PARMAX (RW(IEL1))
C       ===========
        CALL PARMAX (RW(IEL2))
C       ===========
      ENDIF
C
C BOUCLE SUR LES FACES INTERNES
C
C CRITERE A
C ---------
C
      DO IPASS = 1, NPASS
C
      DO IFAC = 1, NFACF
        IEL1 = IFACEF(1,IFAC)
        IEL2 = IFACEF(2,IFAC)
        IF (IW(IEL1).EQ.0 .AND. IW(IEL2).EQ.0) THEN
C CRITERE DEBUT
CMO          RTEST1 = COCRIB*ABS(RW(IEL1))
CMO          RTEST2 = COCRIB*ABS(RW(IEL2))
          RTEST1 = COCRIA*SQRT(ABS(DAF(IEL1))*ABS(DAF(IEL2)))
          RTEST2 = RTEST1
C
          IF (ISYM.EQ.1) THEN
            CRITER = ABS(XAF(IFAC,1)).GE.RTEST1
     &               .OR. ABS(XAF(IFAC,1)).GE.RTEST2
          ELSE
            CRITER = ABS(XAF(IFAC,1)).GE.RTEST1
     &               .OR. ABS(XAF(IFAC,2)).GE.RTEST2
          ENDIF
C CRITERE FIN
          IF (CRITER) THEN
            IW(IEL1) = 1
            IW(IEL2) = 1
            ICELG = ICELG + 1
            IRESPR(IEL1) = ICELG
            IRESPR(IEL2) = ICELG
            IVOIS(IEL1) = IEL2
            IVOIS(IEL2) = IEL1
          ENDIF
        ENDIF
      ENDDO
C
C ACTUALISATION COEF DE COMPARAISON
      COCRIB = COCRIB*COCRIB
C
      ENDDO
C
C CAS DES CELLULES NON REGROUPEES
      DO IEL = 1, NCELF
        IF (IW(IEL).EQ.0) THEN
          ICELG = ICELG + 1
          IRESPR(IEL) = ICELG
          IVOIS(IEL) = IEL
        ENDIF
      ENDDO
C
      NCELG = ICELG
C
C VERIFICATION DEBUT : a enlever par la suite
      ICOMP = 0
      DO IEL = 1, NCELF
        IW(IEL) = 0
        IF (IRESPR(IEL).EQ.0) THEN
          ICOMP = ICOMP + 1
        ENDIF
      ENDDO
      IF (ICOMP.NE.0) THEN
        WRITE(NFECRA,1000) ICOMP
        CALL CSEXIT (1)
      ENDIF
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS calgrg (multigrille)             ',/,
     &'@    *********                                               ',/,
     &'@     IL EXISTE ',I10   ,' ELEMENTS NULS DANS IRESPR         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      DO IEL = 1, NCELF
        IELG = IRESPR(IEL)
        IW(IELG) = IW(IELG) + 1
      ENDDO
      NCELT = 0
      DO IEL = 1, NCELG
        NCELT = NCELT + IW(IEL)
      ENDDO
      IF (NCELT.NE.NCELF) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS calgrg (multigrille)             ',/,
     &'@    *********                                               ',/,
     &'@    ERREUR DANS IRESPR (NCELT.NE.NCELF)                     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C VERIFICATION FIN : a enlever par la suite
C
C=======================================================================
C 2.  CONSTRUCTION DE LA CONNECTIVITE        INVERSE : IP , ICELFA , ICELCE
C=======================================================================
C
C INITIALISATION
C
      DO IEL = 1, NCELF
        IW(IEL) = 0
      ENDDO
C
C CREATION DU TABLEAU IP DE POINTEURS SUR LES VOISINS DE CELLULES
C IP(1)                   : NOMBRE DE VOISINS DE LA CELLULE 1
C IP(K) - IP(K-1), K >= 2 : NOMBRE DE VOISINS DE LA CELLULE K
C IP(NCELF)                : = 2*NFACF
      DO IFAC = 1, NFACF
        IEL1 = IFACEF(1,IFAC)
        IEL2 = IFACEF(2,IFAC)
        IW(IEL1) = IW(IEL1) + 1
        IW(IEL2) = IW(IEL2) + 1
      ENDDO
      DO IEL = 2, NCELF
        IW(IEL) = IW(IEL) + IW(IEL-1)
      ENDDO
C RQ : CI-DESSOUS, IP EST DECREMENTE DU NOMBRE DE VOISINS
C      CAR IL SERA REINCREMENTE DANS LA BOUCLE SUIVANTE
      IP(1) = 0
      DO IEL = 2, NCELF
        IP(IEL) = IW(IEL-1)
      ENDDO
C
      DO IFAC = 1, NFACF
        IEL1 = IFACEF(1,IFAC)
        IEL2 = IFACEF(2,IFAC)
        IP(IEL1) = IP(IEL1) + 1
        IP(IEL2) = IP(IEL2) + 1
        ICELCE(IP(IEL1)) = IEL2
        ICELCE(IP(IEL2)) = IEL1
        ICELFA(IP(IEL1)) = IFAC
        ICELFA(IP(IEL2)) = IFAC
      ENDDO
C
C=======================================================================
C 3.  CREATION DES FACES GROSSIERES : IRESFA
C=======================================================================
C
C INITIALISATION
C
      IFACG = 0
      DO IFAC = 1, NFACF
        INDICF(IFAC)  = 0
        IRESFA(IFAC)  = 0
        IFACEG(1,IFAC) = 0
        IFACEG(2,IFAC) = 0
      ENDDO

C BOUCLE SUR LES FACES
C
      DO IFAC = 1, NFACF
        IF (INDICF(IFAC).EQ.0) THEN
C     PREMIERE FACE
          INDICF(IFAC) = 1
          IEL1  = IFACEF(1,IFAC)
          IEL2  = IFACEF(2,IFAC)
          IEL1G = IRESPR(IEL1)
          IEL2G = IRESPR(IEL2)
C  ON NE PREND PAS LES FACES INTERNES AUX CELLULES GROSSIERES
          IF (IEL1G.NE.IEL2G) THEN
            IEL1V = IVOIS(IEL1)
            IFACG = IFACG + 1
            IRESFA(IFAC) = IFACG
            IFACEG(1,IFACG) = IEL1G
            IFACEG(2,IFACG) = IEL2G
C  DEUXIEME FACE, SI NECESSAIRE
            IF (IEL1V.NE.IEL1) THEN
              IDEBI = 1
              IFINI = IP(1)
              IF (IEL1V.NE.1) THEN
                IDEBI = IP(IEL1V-1) + 1
                IFINI = IP(IEL1V)
              ENDIF
              DO IVOI = IDEBI, IFINI
                JJ = IRESPR(ICELCE(IVOI))
                IF (JJ.EQ.IEL2G) THEN
                  IFAC2 = ICELFA(IVOI)
                  IRESFA(IFAC2) = IFACG
                  INDICF(IFAC2) = 1
                  GOTO 110
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDIF
 110    CONTINUE
      ENDDO
C
      NFACG = IFACG
C
C=======================================================================
C 4.  RESTRICTION DE LA MATRICE
C=======================================================================
C
C INITIALISATIONS
C
      DO IEL = 1, NCELF
        DAG(IEL) = 0.D0
      ENDDO
      DO IFAC = 1, NFACF
        XAG(IFAC,1) = 0.D0
        IF (ISYM.EQ.2) THEN
          XAG(IFAC,2) = 0.D0
        ENDIF
      ENDDO
C
C TERMES EXTRA-DIAGONAUX
      DO IFAC = 1, NFACF
        IFACG = IRESFA(IFAC)
        IF (IFACG.NE.0) THEN
          IF (ISYM.EQ.1) THEN
            XAG(IFACG,1) = XAG(IFACG,1) + XAF(IFAC,1)
          ELSE IF (ISYM.EQ.2) THEN
            IEL1  = IFACEF(1,IFAC)
            IEL1G = IRESPR(IEL1)
            IF (IEL1G.EQ.IFACEG(1,IFACG)) THEN
              XAG(IFACG,1) = XAG(IFACG,1) + XAF(IFAC,1)
              XAG(IFACG,2) = XAG(IFACG,2) + XAF(IFAC,2)
            ELSE
              XAG(IFACG,1) = XAG(IFACG,1) + XAF(IFAC,2)
              XAG(IFACG,2) = XAG(IFACG,2) + XAF(IFAC,1)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C TERMES DIAGONAUX
C   ON MET DANS RW LA CONTRIBUTION NON DIFFERENTIELLE DE DAF
      DO IEL = 1, NCELF
        RW(IEL) = DAF(IEL)
      ENDDO
      DO IFAC = 1, NFACF
        IEL1 = IFACEF(1,IFAC)
        IEL2 = IFACEF(2,IFAC)
        RW(IEL1) = RW(IEL1) + XAF(IFAC,1)
        IF (ISYM.EQ.2) THEN
          RW(IEL2) = RW(IEL2) + XAF(IFAC,2)
        ELSE
          RW(IEL2) = RW(IEL2) + XAF(IFAC,1)
        ENDIF
      ENDDO
C   ON SOMME LES CONTRIBUTIONS NON DIFFERENTIELLES
      DO IEL = 1, NCELF
        IELG = IRESPR(IEL)
        DAG(IELG) = DAG(IELG) + RW(IEL)
      ENDDO
C   ON AJOUTE (-) LES CONTRIBUTIONS EXTRA-DIAG DIFFERENTIELLES
      DO IFAC = 1, NFACG
        IEL1 = IFACEG(1,IFAC)
        IEL2 = IFACEG(2,IFAC)
        DAG(IEL1) = DAG(IEL1) - XAG(IFAC,1)
        IF (ISYM.EQ.1) THEN
          DAG(IEL2) = DAG(IEL2) - XAG(IFAC,1)
        ELSE IF (ISYM.EQ.2) THEN
          DAG(IEL2) = DAG(IEL2) - XAG(IFAC,2)
        ENDIF
      ENDDO
C
C=======================================================================
C 5.  AFFECTATIONS DES POINTEURS DANS IA ET RA
C=======================================================================
C
C POINTEUR SUR NCELG
      JPGREN(JNCEL)  = IDBIGR
C POINTEUR SUR NFACG
      JPGREN(JNFAC)  = JPGREN(JNCEL)  + 1
C POINTEUR SUR IFACEG
      JPGREN(JIFACE) = JPGREN(JNFAC)  + 1
C POINTEUR SUR IRESPR
      JPGREN(JIRESP) = JPGREN(JIFACE) + 2*NFACG
C POINTEUR SUR CASE LIBRE DE IA
      JPGREN(NTGREN) = JPGREN(JIRESP) + NCELF
C POINTEUR SUR DAG
      JPGRRE(JDA)    = IDBRGR
C POINTEUR SUR XAG
      JPGRRE(JXA)    = JPGRRE(JDA) + NCELG
C POINTEUR SUR CASE LIBRE DE RA
      JPGRRE(NTGRRE) = JPGRRE(JXA) + NFACG*ISYM
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
      END
c@z
