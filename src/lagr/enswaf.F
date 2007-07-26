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
                        SUBROUTINE ENSWAF
C                       *****************
C     ------------------------------------------------------------------
     & ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NFIN   ,
     &   ITEPA  ,
     &   ETTP   , TEPA , TRAV   )
C     ------------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC    Ecriture des fichiers pour Ensight7 au format CASE pour la
CFONC    visualisation des deplacements des particules et de variables
CFONC    associees.
CFONC
CFONC    La visualisation des deplacement et le choix des variables
CFONC    associees est realise dans le sous-programme USLAG1.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NFIN         ! E  !  -> ! NFIN = 1 SI DERNIER PAS DE TEMPS     !
CARGU !              !    !     ! NFIN = 0 SINON                       !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES                     !
CARGU !              !    !     !   ETAPE COURANTE OU PRECEDENTE       !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! TRAV(NBPMAX,3! TR !  -  ! TABLEAUX DE TRAVAIL                  !
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
C
C***********************************************************************
C
      IMPLICIT NONE
C
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C**********************************************************************
C
C ARGUMENTS
C
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NFIN
      INTEGER          ITEPA(NBPMAX,NIVEP)
C
      DOUBLE PRECISION ETTP(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TRAV(NBPMAX,3)
C
C VARIABLES LOCALES
C
      INTEGER          NPT , IPT
      INTEGER          NP , NL
      INTEGER          II1 , II2 , LPOS , N1 , N2
C
      CHARACTER        FICH*80 , NAME*80 , ENTET*80
C
      INTEGER          IPWAF
      DATA             IPWAF /0/
      SAVE             IPWAF
C
C***********************************************************************
C
C=======================================================================
C 0. GESTION MEMOIRE
C=======================================================================
C
C=======================================================================
C 1. Initialisations
C=======================================================================
C
      IF (NFIN.EQ.0) IPWAF = IPWAF + 1
C
      IF (IPWAF.EQ.1) ITLAG = 0
C
      FICH = ' '
      FICH = 'deplacement'
      CALL VERLON (FICH,II1,II2,LPOS)
      ENTET = FICH(II1:II2)
C
C=======================================================================
C 2. ENREGISTREMENTS des deplacement.geom****
C=======================================================================
C
C-->Faut-il enregistrer ?
C
      IF ( (MOD(IPWAF-1,NVISLA).EQ.0 .AND. NFIN.EQ.0)   .OR.
     &     (NFIN.EQ.1 .AND. MOD(IPWAF-1,NVISLA).NE.0) ) THEN
C
C-->Nombre de particules a visualisees encore presentent dans le domaine
C
        NPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) NPT = NPT + 1
        ENDDO
C
C-->Y a t-il encore des particules a visualiser ?
C
        IF (NPT.EQ.0) GOTO 100
C
C-->Nombre d'enregistrements et incrementation du temps physique
C
        IF (ITLAG.LE.9999) THEN
          ITLAG = ITLAG + 1
          TIMLAG(ITLAG) = TTCLAG
        ELSE
          WRITE(NFECRA,9000) ITLAG
          GOTO 100
        ENDIF
C
C-->Ouverture des fichiers type deplacement.geo0001
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+4) = '.geo'
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        II2 = II2 + 5
        FICH(II2:II2+LPOS) = NAME(N1:N2)
C
        II2 = II2 + LPOS
        OPEN ( IMPLA1, FILE=FICH(II1:II2),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture de l'entete
C
        WRITE(IMPLA1,'(A)') 'geometrie deplacement'
        WRITE(IMPLA1,'(A)') 'au format ensight6'
        WRITE(IMPLA1,'(A)') 'node id given'
        WRITE(IMPLA1,'(A)') 'element id given'
        WRITE(IMPLA1,'(A)') 'coordinates'
        WRITE(IMPLA1,'(I8)') NPT
C
C-->Ecriture des points
C
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            WRITE(IMPLA1,'(I8,3E12.5)') NP,
     &                                ETTP(NP,JXP),
     &                                ETTP(NP,JYP),
     &                                ETTP(NP,JZP)
          ENDIF
        ENDDO
C
C-->Ecriture de la geometrie Ensight
C
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'deplacements'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(I8)') NPT
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            WRITE(IMPLA1,'(2I8)') IPT , NP
          ENDIF
        ENDDO
        CLOSE(IMPLA1)
C
      ELSE
C
        IF (NFIN.EQ.0) RETURN
        GOTO 100
C
      ENDIF
C
C=======================================================================
C 3. Ecriture de deplacement.tpssej0001
C=======================================================================
C
      IF (IVISTP.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = TEPA(NP,JRTSP)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.tpssej'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 4. Ecriture de deplacement.temper0001
C=======================================================================
C
      IF (IVISTE.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JTP)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.temper'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 5. Ecriture de deplacement.diamet0001
C=======================================================================
C
      IF (IVISDM.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JDP)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.diamet'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 6. Ecriture de deplacement.massep0001
C=======================================================================
C
      IF (IVISMP.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JMP)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.massep'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 7. Charbon : Ecriture de deplacement.temp_ch0001
C=======================================================================
C
      IF (IVISHP.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JHP)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.tempch'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 8. Charbon : Ecriture de deplacement.dck0001
C=======================================================================
C
      IF (IVISDK.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = TEPA(NP,JRDCK)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+4) = '.dck'
        II2 = II2 + 4
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 9. Charbon : Ecriture de deplacement.mch0001
C=======================================================================
C
      IF (IVISCH.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JMCH)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+4) = '.mch'
        II2 = II2 + 4
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 10. Charbon : Ecriture de deplacement.mck0001
C=======================================================================
C
      IF (IVISCK.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JMCK)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.mck'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 11. Ecriture de deplacement.vitflu0001
C=======================================================================
C
      IF (IVISV1.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JUF)
            TRAV(IPT,2) = ETTP(NP,JVF)
            TRAV(IPT,3) = ETTP(NP,JWF)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.vitflu'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( (REAL(TRAV(NP,NL)),NL=1,3),NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 12. Ecriture de deplacement.vitpar0001
C=======================================================================
C
      IF (IVISV2.EQ.1) THEN
C
        IPT = 0
        DO NL = 1,NBVIS
          NP = LISTE(NL)
          IF (NP.GE.1 .AND. ITEPA(NP,JISOR).NE.0) THEN
            IPT = IPT + 1
            TRAV(IPT,1) = ETTP(NP,JUP)
            TRAV(IPT,2) = ETTP(NP,JVP)
            TRAV(IPT,3) = ETTP(NP,JWP)
          ENDIF
        ENDDO
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+7) = '.vitpar'
        II2 = II2 + 7
        WRITE (NAME,'(I4.4)') ITLAG
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
C
        OPEN ( IMPLA1, FILE=FICH(II1:II2+LPOS),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C
C-->Ecriture
C
        WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
        WRITE(IMPLA1,'(A)') 'part   1'
        WRITE(IMPLA1,'(A)') 'point'
        WRITE(IMPLA1,'(6E12.5)') ( (REAL(TRAV(NP,NL)),NL=1,3),NP=1,NPT )
        CLOSE(IMPLA1)
C
      ENDIF
C
C=======================================================================
C 13. Ecriture du deplacement.case au dernier passage
C=======================================================================
C
 100  CONTINUE
C
      IF (NFIN.EQ.1) THEN
C
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        NAME = ' '
        NAME = '.CASE'
        CALL VERLON (NAME,N1,N2,LPOS)
        FICH(II2+1:II2+LPOS) = NAME(N1:N2)
        II2 = II2 + LPOS
        OPEN ( UNIT=IMPLA1, FILE=FICH (II1:II2),
     &         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
        REWIND ( UNIT=IMPLA1 )
C
        WRITE(IMPLA1,'(A)') 'FORMAT'
        WRITE(IMPLA1,'(A)') 'type:     ensight'
        WRITE(IMPLA1,'(A)') 'GEOMETRY'
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,II1,II2,LPOS)
        FICH(II2+1:II2+8) = '.geo****'
        CALL VERLON (FICH,II1,II2,LPOS)
        NAME = ' '
        NAME = 'model:                 1    '
        NAME(29:29+II2-II1+1) = FICH(II1:II2)
        CALL VERLON (NAME,II1,II2,LPOS)
        WRITE(IMPLA1,'(A)') NAME(II1:II2)
C
        WRITE(IMPLA1,'(A)') 'VARIABLE'
        FICH = ' '
        FICH = ENTET
        CALL VERLON (FICH,N1,N2,LPOS)
C
C Rem : les trois lignes suivantes sont pour eviter une erreur
C       de lecture fichier .CASE lors de sa lecture par ensight
C       s'il n'y a aucune VARIABLE a voir.
C
C       NAME = 'constant per case :    1    constant     1.0'
C       CALL VERLON (NAME,II1,II2,LPOS)
C       WRITE(IMPLA1,'(A)') NAME(II1:II2)
C
        IF (IVISTP.EQ.1) THEN
          NAME = 'scalar per node :      1    temps_de_sejour  '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 2
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.tpssej****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISTE.EQ.1) THEN
          NAME = 'scalar per node :      1    temperature      '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 6
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.temper****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISDM.EQ.1) THEN
          NAME = 'scalar per node :      1    diametre         '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 9
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.diamet****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISMP.EQ.1) THEN
          NAME = 'scalar per node :      1    masse            '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 12
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.massep****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISHP.EQ.1) THEN
          NAME = 'scalar per node :      1    tempch           '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 11
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.tempch****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISDK.EQ.1) THEN
          NAME = 'scalar per node :      1    dck              '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 14
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.dck****'
          II2 = II2 + 8
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISCH.EQ.1) THEN
          NAME = 'scalar per node :      1    mch              '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 14
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.mch****'
          II2 = II2 + 8
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISCK.EQ.1) THEN
          NAME = 'scalar per node :      1    mck              '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 14
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.mck****'
          II2 = II2 + 8
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISV1.EQ.1) THEN
          NAME = 'vector per node :      1    vitesse_fluide   '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 3
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.vitflu****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        IF (IVISV2.EQ.1) THEN
          NAME = 'vector per node :      1    vitesse_partic   '
          CALL VERLON (NAME,II1,II2,LPOS)
          II2 = II2 + 3
          NAME(II2+1:II2+N2)=FICH(N1:N2)
          II2 = II2 + N2
          CALL VERLON (NAME,II1,II2,LPOS)
          NAME(II2+1:II2+11) = '.vitpar****'
          II2 = II2 + 11
          WRITE(IMPLA1,'(A)') NAME(II1:II2)
        ENDIF
C
        WRITE(IMPLA1,'(A)') 'TIME'
        WRITE(IMPLA1,'(A)') 'time set:'
        FICH = ' '
        FICH = 'number of steps:'
        WRITE(NAME,'(I4)') ITLAG
        CALL VERLON(NAME,N1,N2,LPOS)
C       ===========
        FICH(25+1:25+LPOS) = NAME(N1:N2)
        WRITE(IMPLA1,'(A)') FICH(1:25+LPOS)
        WRITE(IMPLA1,'(A)') 'filename start number:  1'
        WRITE(IMPLA1,'(A)') 'filename increment:     1'
        WRITE(IMPLA1,'(A)') 'time values:'
        WRITE(IMPLA1,'(6E12.5)') (TIMLAG(NL),NL=1,ITLAG)
C
        CLOSE(IMPLA1)
C
      ENDIF
C
      RETURN
C
C-------
C FORMAT
C-------
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A L''EXECUTION DU MODULE LAGRANGIEN  ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    LE NOMBRE D''ENREGISTREMENTS TEMPORELS DEMANDES POUR    ',/,
     &'@      LE POST-PROCESSING EN MODE DEPLACEMENT DEPASSE        ',/,
     &'@      LE MAXIMUM ADMISSIBLE.                                ',/,
     &'@                                                            ',/,
     &'@    LE NOMBRE DE PAS DE TEMPS DEMANDE EST DE : ',I10         ,/,
     &'@      LE MAXIMUM ADMISSIBLE EST 9999                        ',/,
     &'@                                                            ',/,
     &'@  Le calcul continue, mais les enregistrements sont arretes.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
