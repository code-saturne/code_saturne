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
                        SUBROUTINE AUTMGR
C                       *****************
C
C     -------------------------------------------------------------
     &    ( IGR    , ISYM   , IAGMAX , NAGMAX ,
     &      NCELF  , NCELFE , NFACF  , IWARNP ,
     &      IFACEF ,
     &      DAF    , XAF    , SURFAF , VOLUMF , XYZFIN ,
     &      IRSCEL ,
     &      INDIC  , INOMBR , IRSFAC , INDICF , W1     , W2 )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC   MULTIGRILLE ALGEBRIQUE :
CFONC   CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER A PARTIR
CFONC   DU NIVEAU SUPERIEUR SUIVANT CRITERE AUTOMATIQUE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NFECRA       ! E  !  -> ! NUMERO DU FICHIER D'IMPRESSIONS      !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! IGR          ! E  !  -> ! NIVEAU DU MAILLAGE GROSSIER          !
CARGU ! NCELF        ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE FIN       !
CARGU ! NCELFE       ! E  !  -> ! NOMBRE D'ELEMENTS ETENDUS FIN        !
CARGU ! NFACF        ! E  !  -> ! NOMBRE DE FACES INTERNES MAILL. FIN  !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! IFACEF       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFACF)   !    !     !  DU MAILLAGE FIN                     !
CARGU ! DAF(NCELFE)  ! TR !  -> ! DIAGONALE MATRICE MAILLAGE FIN       !
CARGU ! XAF          ! TR !  -> ! EXTRADIAGONALE MATRICE MAILLAGE FIN  !
CARGU ! (NFACF, ISYM)!    !     !                                      !
CARGU ! SURFAF       ! TR !  -> ! SURFACES FACES INTERNES MAILLAGE FIN !
CARGU ! (3, NFACF)   !    !     !                                      !
CARGU ! VOLUMF       ! TR !  -> ! VOLUMES DES CELLULES DU MAILLAGE FIN !
CARGU ! (NCELFE)     !    !     !                                      !
CARGU ! XYZFIN       ! TR !  -> ! CENTRES DES CELLULES DU MAILLAGE FIN !
CARGU ! (3, NCELFE)  !    !     !                                      !
CARGU ! IRSCEL       ! TE ! <-  ! CELLULE FINE -> CELLULE GROSSIERE    !
CARGU !  (NCELFE)    !    !     !                                      !
CARGU ! INDIC(NCELFE)! TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! INOMBR       ! TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !  (NCELFE)    !    !     !                                      !
CARGU ! IRSFAC       ! TE !  -  ! FACE FINE -> FACE GROSSIERE          !
CARGU !  (NFACF)     !    !     !  (TABLEAU DE TRAVAIL)                !
CARGU ! INDICF(NFACF)! TE !  -  ! INDICATEUR DE REGROUPEMENT DES FACES !
CARGU ! ICELFA       ! TE !  -  ! CONNECTIVITE CELLULES->FACES MAILLA- !
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
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IGR, ISYM, IAGMAX, NAGMAX
      INTEGER          NCELF, NCELFE, NFACF
      INTEGER          IWARNP
C
      INTEGER          IFACEF(2, NFACF)
      INTEGER          IRSCEL(NCELFE)
      INTEGER          INDIC(NCELFE), INOMBR(NCELFE)
      INTEGER          INDICF(NFACF), IRSFAC(NFACF)
C
      DOUBLE PRECISION DAF(NCELFE), XAF(NFACF,ISYM)
      DOUBLE PRECISION SURFAF(3,NFACF), VOLUMF(NCELFE)
      DOUBLE PRECISION XYZFIN(3,NCELFE)
      DOUBLE PRECISION W1(NCELFE), W2(NCELFE)
C
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION CRITR, EPSLON
C
      INTEGER          NCELG, ICEL, IFAC , IFAC1, ICELG
      INTEGER          NFACN,NFACNR,NPASS,NPASMX
      INTEGER          INDITT, NOAGLO, NGROS, INCVOI
      INTEGER          IHIST(10)
      INTEGER          I, J, IMIN, IMAX
C
C
C***********************************************************************
C
C     PARAMETRES PAR DEFAUT
C
      EPSLON = +1.D-6
      NGROS  = 8
      NPASMX = 10
      INCVOI = 1
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      DO ICEL = 1, NCELFE
        INDIC(ICEL) = -1
        IRSCEL(ICEL) = 0
        INOMBR(ICEL) = 1
      ENDDO
      DO IFAC = 1, NFACF
        INDICF(IFAC) = IFAC
        IRSFAC(IFAC) = IFAC
      ENDDO
C
C     CALCUL DU CARDINAL (NOMBRE DE VOISINS DE CHAQUE CELLULE -1)
      DO IFAC = 1, NFACF
        I = IFACEF(1,IFAC)
        J = IFACEF(2,IFAC)
        INDIC(I) = INDIC(I) + 1
        INDIC(J) = INDIC(J) + 1
      ENDDO
C
      NCELG  = 0
      NFACNR = NFACF
      NPASS  = 0
      NOAGLO = NCELF
C
 100  CONTINUE
C
C     LES PASSES
C
      NPASS = NPASS+1
      NFACN = NFACNR
      IAGMAX= IAGMAX +1
      IAGMAX= MIN(IAGMAX, NAGMAX)
c
      DO IFAC=1,NFACN
        IRSFAC(IFAC) = INDICF(IFAC)
        INDICF(IFAC) = 0
      ENDDO
      IF (NFACN .LT. NFACF) THEN
        DO IFAC = NFACN+1, NFACF
          INDICF(IFAC) = 0
          IRSFAC(IFAC) = 0
        ENDDO
      ENDIF
C
      IF (IWARNP .GT. 3) THEN
        WRITE(NFECRA,*) '    autmgr.F : passe ', NPASS,
     &                  'NFACNR = ', NFACNR, ' NOAGLO = ',NOAGLO
      ENDIF
C
C     INCREMENTATION DU NOMBRE DE VOISINS
C
      DO ICEL = 1, NCELF
        INDIC(ICEL) = INDIC(ICEL) + INCVOI
      ENDDO
C
C     INITIALISATION DE FACES NON ELIMINEES
C
      NFACNR = 0
C
C     BOUCLE SUR LES FACES NON ELIMINEES
C
      DO IFAC1 = 1, NFACN
C
        IFAC = IRSFAC(IFAC1)
        I = IFACEF(1,IFAC)
        J = IFACEF(2,IFAC)
C
C       ON NE CONSIDERE PAS LES FACES EN FRONTIERE PARALLELE OU
C       PERIODIQUE, POUR NE PAS AGGLOMERER LA GRILLE A TRAVERS
C       CE TYPE DE FRONTIERES (CE QUI DEMANDERAIT UNE CONSTRUCTION
C       PLUS COMPLEXE PUIS CHANGERAIT LE SCEMA DE COMMUNICATION).
C
        IF (I.LE.NCELF .AND. J.LE.NCELF) THEN
C
          INDITT = 0
          CRITR  = (DAF(I)/INDIC(I))*(DAF(J)/INDIC(J))
     &             /( XAF(IFAC,1)*XAF(IFAC,ISYM))
C
          IF (       CRITR.LT.(1.D0-EPSLON)
     &         .AND. IRSCEL(I)*IRSCEL(J).LE.0) THEN
C
            IF (IRSCEL(I).GT.0 .AND. IRSCEL(J).LE.0) THEN
              IF(INOMBR(IRSCEL(I)) .LE. IAGMAX) THEN
                IRSCEL(J) = IRSCEL(I)
                INOMBR(IRSCEL(I)) = INOMBR(IRSCEL(I)) +1
                INDITT = INDITT +1
              ENDIF
            ELSE IF (IRSCEL(I).LE.0 .AND. IRSCEL(J).GT.0) THEN
              IF (INOMBR(IRSCEL(J)).LE.IAGMAX) THEN
                IRSCEL(I) = IRSCEL(J)
                INOMBR(IRSCEL(J)) = INOMBR(IRSCEL(J)) + 1
                INDITT = INDITT +1
              ENDIF
            ELSE IF (IRSCEL(I).LE.0 .AND. IRSCEL(J).LE.0) THEN
              NCELG = NCELG+1
              IRSCEL(I) = NCELG
              IRSCEL(J) = NCELG
              INOMBR(NCELG) = INOMBR(NCELG) +1
              INDITT = INDITT +1
            ENDIF
C
          ENDIF
C
          IF (INDITT.NE.0 .AND.INDITT.NE.1) THEN
            WRITE(NFECRA,*)' Bug dans autmgr.F, arret '
            CALL CSEXIT(1)
          ENDIF
C
          IF (INDITT.EQ.0 .AND. IRSCEL(I)*IRSCEL(J).LE.0) THEN
            NFACNR = NFACNR +1
            INDICF(NFACNR) = IFAC
          ENDIF
C
        ENDIF
C
      ENDDO
C
C     CONTROLE DU NOMBRE DE CELLULES GROSSIERES CREE
C
      NOAGLO = 0
      DO ICEL=1,NCELF
        IF (IRSCEL(ICEL).LE.0) NOAGLO = NOAGLO+1
      ENDDO
C
C     PASSES SUIVANTES SI AGGLOMERATION INSUFFISANTE
C
      IF (NOAGLO.GT.0) THEN
        IF ((NCELG+NOAGLO)*NGROS .GE. NCELF) THEN
          IF (NPASS.LT.NPASMX .AND. NFACNR.GT.0) THEN
            GOTO 100
          ENDIF
        ENDIF
      ENDIF
C
C     ON TERMINE L'ASSEMBLAGE
C
      DO ICEL = 1, NCELF
        IF (IRSCEL(ICEL).LE.0) THEN
          NCELG = NCELG+1
          IRSCEL(ICEL) = NCELG
        ENDIF
      ENDDO
C
C     CONTROLE DIVERS ET VARIES
C
      IMAX = 0
      IMIN = 2*NCELF
      DO ICELG =1,NCELG
        IMAX = MAX(IMAX, INOMBR(ICELG))
        IMIN = MIN(IMIN, INOMBR(ICELG))
      ENDDO
C
      IF (IRANGP .GE. 0) THEN
        CALL PARCMN(IMIN)
        CALL PARCMX(IMAX)
      ENDIF
C
      IF (IWARNP.GT.3) THEN
C
        WRITE(NFECRA,*) '    autmgr.F : INOMBR MIN = ', IMIN,
     &                  ' MAX = ', IMAX, ' CIBLE = ', NAGMAX
        WRITE(NFECRA,*) '      histogramme '
        NOAGLO=IMAX-IMIN+1
        IF (NOAGLO.GT.0) THEN
          IF (NOAGLO.GT.10) THEN
            WRITE(NFECRA,*) ' IHIST MAL DIMENSIONNE DANS autmgr.F'
            CALL CSEXIT(1)
          ENDIF
          DO I = 1, NOAGLO
            IHIST(I) = 0
          ENDDO
          DO ICELG = 1, NCELG
            DO I = 1, NOAGLO
              IF (INOMBR(ICELG).EQ.(IMIN+I-1))THEN
                IHIST(I)=IHIST(I)+1
              ENDIF
            ENDDO
          ENDDO
          IF (IRANGP .GE. 0) THEN
            CALL PARISM(NOAGLO, IHIST)
          ENDIF
          DO I = 1, NOAGLO
            EPSLON = 100.D0*IHIST(I)/NCELG
            WRITE(NFECRA,*) '        regroupement ',IMIN+I-1,'  =  ',
     &                      EPSLON,' % '
          ENDDO
        ENDIF
C
      ENDIF
C
      DO ICEL = 1, NCELF
        INDIC(ICEL) = 0
      ENDDO
      DO ICEL = 1, NCELF
        ICELG = IRSCEL(ICEL)
        INDIC(ICELG) = INDIC(ICELG) +1
      ENDDO
C
      I=0
      J=2*NCELF
      NOAGLO = 0
      DO ICELG = 1, NCELG
        I = MAX(I, INDIC(ICELG))
        J = MIN(J, INDIC(ICELG))
        NOAGLO = NOAGLO + INDIC(ICELG)
      ENDDO
C
      IF (IRANGP .GE. 0) THEN
        CALL PARCMN(J)
        CALL PARCMX(J)
      ENDIF
C
      IF (IWARNP.GT.3) THEN
        WRITE(NFECRA,*) '    autmgr.F : agglomeration MIN = ', J,
     &                  ' MAX= ',I
      ENDIF
C
      IF (NOAGLO .NE. NCELF) THEN
        WRITE(NFECRA,*) 'BUG DANS autmgr.f, Contacter l''assistance.'
        CALL CSEXIT(1)
      ENDIF
C
C**********************************************************************
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
      END
c@z

