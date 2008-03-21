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
                        SUBROUTINE PPTBHT
C                       *****************
C     -------------------------------------------------------------
     & ( NCOEL  ,
     &   NOMCOE , EHCOEL , CPCOEL , WMOLCE )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC          PHYSIQUES PARTICULIERES
CFONC
CFONC            CALCUL DE L'ENTHALPIE ET DU CP
CFONC                     A PARTIR DE LA BANDE DE JANAF
CFONC
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU ._______________.____._____.______________________________________.
CARGU !    NOM        !TYPE!MODE !                   ROLE               !
CARGU !_______________!____!_____!______________________________________!
CARGU ! NCOEL         ! E  !  -> ! NOMBRE DE CONST. ELEM.               !
CARGU ! NOMCOE(NGAZEM)! A  !  -> ! NOM DES CONSTITUANTS ELEMENTAIRES    !
CARGU ! EHCOEL        ! TR !  <- ! ENTHALPIE POUR CHAQUE CONSTITUANT    !
CARGU ! (NGAZEM,NPOT) !    !     !                ELEMENTAIRE           !
CARGU ! CPCOEL(NGAZEM)! TR !  <- ! CP POUR CHAQUE CONSTITUANT           !
CARGU ! WMOLCE(NGAZEM)! TR !  <- ! MASSE MOLAIRE DE CHAQUE CONSTITUANT  !
CARGU !_______________!____!_____!______________________________________!
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
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"

C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCOEL
C
      CHARACTER*12     NOMCOE(NGAZEM)
C
      DOUBLE PRECISION EHCOEL(NGAZEM,NPOT) , CPCOEL(NGAZEM)
      DOUBLE PRECISION WMOLCE (NGAZEM)
C
C VARIABLES LOCALES
C
      CHARACTER*40     DUMMY
      CHARACTER*12     NOMESP
C
      INTEGER          IND , ICHES , INDTP , INICFF , INJCFF
      INTEGER          NE   , NT  , NC , IOK
      INTEGER          ICOEFF(NGAZEM)
C
      DOUBLE PRECISION CTH  , CTC
C
      DOUBLE PRECISION TLIM(3) , WCOEFF(2,7) , COEFF(NGAZEM,2,7)
C
C***********************************************************************
C
C
C=======================================================================
C 2. LECTURE DU FICHIER DE DONNEES THERMODYNAMIQUES TABLE DE JANAF
C=======================================================================
C
C Initialisation
C
C
      DO ICHES= 1, 12
        NOMESP(ICHES:ICHES)=' '
      ENDDO
C
      DO NE = 1 , NGAZEM
        ICOEFF(NE) = 0
        DO INICFF = 1, 2
          DO INJCFF = 1, 7
            COEFF(NE,INICFF,INJCFF) = 0.D0
          ENDDO
        ENDDO
      ENDDO
C
      DO NE = 1 , NCOEL
        CPCOEL(NE)= 0.D0
        DO NT = 1, NPO
          EHCOEL(NE,NT)= 0.D0
        ENDDO
      ENDDO
C
      OPEN(UNIT=IMPJNF, FILE='JANAF', STATUS='OLD' , FORM='FORMATTED')
C
      READ (IMPJNF,'(A)') DUMMY
C
C Lecture des domaines de temperature
C
      READ (IMPJNF,*) (TLIM(INDTP) , INDTP=1,3)
C
C Boucle de lecture des especes chimiques avec stockage partiel
C
 5    CONTINUE
C
      READ (IMPJNF,'(A12,6X,A6)') NOMESP,DUMMY
C
      IF (NOMESP(1:3).EQ.'END') GOTO 100
C
      READ (IMPJNF,*) (WCOEFF(1,INJCFF), INJCFF=1,5)
      READ (IMPJNF,*) (WCOEFF(1,INJCFF), INJCFF=6,7),
     &                (WCOEFF(2,INJCFF), INJCFF=1,3)
      READ (IMPJNF,*) (WCOEFF(2,INJCFF), INJCFF=4,7)
C
C On ne stocke les coefficients que si
C  l'espece consideree fait partie de l'exemple
C
      DO NE = 1, NCOEL
        IF ( NOMCOE(NE).EQ.NOMESP ) THEN
          ICOEFF(NE) = 1
           DO INICFF = 1, 2
             DO INJCFF = 1, 7
               COEFF(NE,INICFF,INJCFF) = WCOEFF(INICFF,INJCFF)
             ENDDO
           ENDDO
        ENDIF
      ENDDO
C
      GOTO 5
C
 100  CONTINUE
C
C Arret de la lecture si tous les renseignements necessaires ont ete
C enregistres
C
      CLOSE(IMPJNF)
C
C Test et stop eventuel
C
      IOK = 0
      DO NE = 1, NCOEL
        IF(ICOEFF(NE).EQ.0) THEN
          IOK = IOK + 1
          WRITE(NFECRA,1000) NOMCOE(NE)
        ENDIF
      ENDDO
      IF(IOK.NE.0) THEN
        WRITE(NFECRA,1100) IOK
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,
     &'@    *********                                               ',/,
     &'@    L''ESPECE ',A12     ,' EST INCONNUE DANS JANAF          ',/,
     &'@                                                            ',/,
     &'@                                                            '  )
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE                                   ',/,
     &'@                                                            ',/,
     &'@    LE FICHIER PARAMETRIQUE FAIT REFERENCE A ',I10           ,/,
     &'@       ESPECE(S) INCONNUE(S) DANS JANAF (data/thch).        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C==============================================================================
C 3. CALCUL DES ENTHALPIES ET DES CP
C==============================================================================
C
C Constante des gaz parfaits en J/mol/K
C
      DO NT = 1,NPO
C
C Determination du jeu de coefficients utilises
C
        IF (TH(NT) .GT. TLIM(2)) THEN
          IND = 1
        ELSE
          IND = 2
        ENDIF
C
        DO NE = 1, NCOEL
          EHCOEL(NE,NT)  = COEFF(NE,IND,6) + COEFF(NE,IND,1) * TH(NT)
          CPCOEL(NE)     =                   COEFF(NE,IND,1)
          CTH = TH(NT)
          CTC = 1.D0
C
C Dans la table de Janaf, les COEFF sont adimensionnels (CP/R,H/R)
C
          DO NC = 2, 5
            CTH = CTH * TH(NT)
            CTC = CTC * TH(NT)
            EHCOEL(NE,NT) = EHCOEL(NE,NT)
     &                    + COEFF(NE,IND,NC) * CTH / DBLE(NC)
            CPCOEL(NE) = CPCOEL(NE) + COEFF(NE,IND,NC) * CTC
          ENDDO
C
C      Calcul du CP et du H pour chaque espece
C
          EHCOEL(NE,NT) = EHCOEL(NE,NT) * RR / WMOLCE(NE)
          CPCOEL(NE)    = CPCOEL(NE)    * RR / WMOLCE(NE)
        ENDDO
C
      ENDDO
C
      RETURN
C
      END
c@z
