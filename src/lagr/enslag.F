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
                        SUBROUTINE ENSLAG
C                       *****************
C     ------------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NFIN   , IFORCE ,
     &   ITEPA  ,
     &   ETTP   , TEPA   , RA)
C     ------------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     SAUVEGARDE DES TRAJECTOIRES AU FORMAT ENSIGHT GOLD
CFONC
CFONC     1. AU PREMIER PASSAGE : OUVERTURE  D'UN FICHIER
CFONC        SCRATCH ET DEBUT D'ECRITURE SEQUENTIELLE
CFONC
CFONC     2. A CHAQUE PASSAGE : ECRITURE DANS LE TEMPORAIRE
CFONC
CFONC     3. AU DERNIER PASSAGE : ECRITURE SUR LE FICHIER FINAL
CFONC        AU BON FORMAT POUR ENSIGHT
CFONC
CFONC     REMARQUE : on dispose de 15 fichiers de sortie,
CFONC                les unites logiques sont defini dans INIINI.F
CFONC                et on y accede par IMPLA5(1) a IMPLA5(15)
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
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NFIN         ! E  !  -> ! NFIN = 1 SI DERNIER PAS DE TEMPS     !
CARGU !              !    !     ! NFIN = 0 SINON                       !
CARGU ! IFORCE       ! E  !  -> ! FORCE L'ECRITURE SI = NUMERO DE LA   !
CARGU !              !    !     !   PARTICULE COURANTE                 !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES                     !
CARGU !              !    !     !   ETAPE COURANTE OU PRECEDENTE       !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NFIN   , IFORCE
      INTEGER          ITEPA(NBPMAX,NIVEP)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          NL , NP
C
      INTEGER          I , J
      INTEGER          NUML
      INTEGER          NUME(NLISTE)
      INTEGER          II1 , II2 , LPOS , NLMAX , IOS
      INTEGER          NPT , IPT , LMAX
      INTEGER          IX,IY,IZ,IUP,IVP,IWP,ITSJ,II
      INTEGER          IXPL , IYPL  , IZPL
      INTEGER          IU1L , IV1L  , IW1L
      INTEGER          IU2L , IV2L  , IW2L
      INTEGER          ITPL , IDML  , ITEL  , IMPL
      INTEGER          IHPL , IDCKL , IMCHL , IMCKL
      INTEGER          IFINIA, IFINRA
C
      DOUBLE PRECISION XPL , YPL , ZPL
      DOUBLE PRECISION U1L , V1L , W1L
      DOUBLE PRECISION U2L , V2L , W2L
      DOUBLE PRECISION TPL , DML , TEL , MPL
      DOUBLE PRECISION HPL , DCKL , MCHL , MCKL
      CHARACTER        FICH*80
      CHARACTER        NAME*80
C
      INTEGER IPASS
      DATA    IPASS /0/
      SAVE    IPASS
C
C**********************************************************************
C
C=======================================================================
C -1.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
C--> Compteur de passages
C
      IPASS = IPASS+1
C
      NLMAX = 0
C
C=======================================================================
C 1. OUVERTURE DU FICHIER TAMPON AU PREMIER PASSAGE
C=======================================================================
C
      IF (IPASS.EQ.1) THEN
C
        OPEN (IMPLA3,FILE='SCRATCH3.lag',
     &        STATUS='UNKNOWN',FORM='UNFORMATTED',
     &        ACCESS='SEQUENTIAL')
C
        DO NL = 1,NBVIS
          NPLIST(NL) = 0
          LIST0(NL) = LISTE(NL)
        ENDDO
C
      ENDIF
C
C=======================================================================
C 2. ECRITURE SEQUENTIELLE DES INFOS TRAJECTOIRES
C    POUR CHAQUE PARTICULE A VISUALISER SUIVANT LA FREQUENCE
C=======================================================================
C
      IF ((MOD(IPASS-1,NVISLA).EQ.0 .OR. IFORCE.GT.0)
     &   .AND. NFIN.EQ.0) THEN
C
C
        DO NL = 1,NBVIS
C
          NP = LISTE(NL)
C
          IF ( (NP.GE.1 .AND. IFORCE.EQ.0 ) .OR.
     &                       (IFORCE.EQ.NP)      ) THEN
C
C sortie du domaine ?
            IF (ITEPA(NP,JISOR).GT.0) THEN
C
C--->incrementation du nombre d'enregistrement pour la particule NP :
              NPLIST(NL) = NPLIST(NL)+1
C
              IF (NPLIST(NL).GT.NLMAX) NLMAX = NPLIST(NL)
C
C--->numero de liste :
              WRITE(IMPLA3) NL
C
C--->coordonnees de la particule NP :
              WRITE(IMPLA3) ETTP(NP,JXP), ETTP(NP,JYP), ETTP(NP,JZP)
C
C--->vitesse du fluide vu :
              IF (IVISV1.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JUF), ETTP(NP,JVF), ETTP(NP,JWF)
              ENDIF
C
C--->vitesse de la particule :
              IF (IVISV2.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JUP), ETTP(NP,JVP), ETTP(NP,JWP)
              ENDIF
C
C--->temps de sejour :
              IF (IVISTP.EQ.1) THEN
                WRITE(IMPLA3) TEPA(NP,JRTSP)
              ENDIF
C
C--->diametre :
              IF (IVISDM.EQ.1) THEN
                  WRITE(IMPLA3) ETTP(NP,JDP)
              ENDIF
C
C--->masse :
              IF (IVISMP.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JMP)
              ENDIF
C
C--->temperature :
              IF (IVISTE.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JTP)
              ENDIF
C
C--->Specifique charbon :
C        Temperature
              IF (IVISHP.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JHP)
              ENDIF
C        Diametre du coeur retrecisant
              IF (IVISDK.EQ.1) THEN
                WRITE(IMPLA3) TEPA(NP,JRDCK)
              ENDIF
C        Masse charbon reactif
              IF (IVISCH.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JMCH)
              ENDIF
C        Masse de coke
              IF (IVISCK.EQ.1) THEN
                WRITE(IMPLA3) ETTP(NP,JMCK)
              ENDIF
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
C
      ENDIF
C
C=======================================================================
C 3.1 FIN DU CALCUL (NFIN=1) : OUVERTURE DES FICHIERS RESU
C=======================================================================
C
      IF (NFIN.EQ.1) THEN
C
        NAME = ' '
        NAME = 'trajectoire'
C
C  0) ouverture du fichier .ensight.CASE :
C
        FICH = NAME
        CALL VERLON ( FICH, II1, II2, LPOS )
        FICH(II2+1:II2+14) = '.ensight.CASE'
        II2 = II2 + 14
        OPEN ( UNIT=IMPLA2, FILE=FICH (II1:II2),
     &         STATUS='UNKNOWN', FORM='FORMATTED',
     &         ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
        REWIND ( UNIT=IMPLA2,ERR=99 )
C
        REWIND(IMPLA2)
        WRITE(IMPLA2,5010)
        WRITE(IMPLA2,5011)
        WRITE(IMPLA2,5012)
C
C  1) ouverture du fichier .ensight.geom + entete fichier case(suite)
C
        FICH = NAME
        CALL VERLON ( FICH, II1, II2, LPOS )
        FICH(II2+1:II2+13) = '.ensight.geom'
        II2 = II2 + 13
C
        WRITE(IMPLA2,5013) FICH (II1:II2)
        WRITE(IMPLA2,5014)
C
        OPEN ( UNIT=IMPLA1, FILE=FICH (II1:II2),
     &         STATUS='UNKNOWN', FORM='FORMATTED',
     &         ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
C
        REWIND ( UNIT=IMPLA1,ERR=99 )
C
C  2) ouverture du fichier .vitflu + entete fichier case(suite)
C
        IF (IVISV1.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.vitflu'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(1), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(1),ERR=99 )
C
          WRITE(IMPLA2,5015) FICH (II1:II2)
C
        ENDIF
C
C  3) ouverture du fichier .vitpar + entete fichier case(suite)
C
        IF (IVISV2.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.vitpar'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(2), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(2),ERR=99 )
C
          WRITE(IMPLA2,5016) FICH (II1:II2)
C
        ENDIF
C
C  4) ouverture du fichier .tpssej + entete fichier case(suite)
C
        IF (IVISTP.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.tpssej'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(3), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(3),ERR=99 )
C
          WRITE(IMPLA2,5017) FICH (II1:II2)
C
        ENDIF
C
C  5) ouverture du fichier .diamet + entete fichier case(suite)
C
        IF (IVISDM.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.diamet'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(4), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(4),ERR=99 )
C
          WRITE(IMPLA2,5018) FICH (II1:II2)
C
        ENDIF
C
C  6) ouverture du fichier .masse + entete fichier case(suite)
C
        IF (IVISMP.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.masse'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(5), FILE=FICH (II1:II2),
     *           STATUS='UNKNOWN', FORM='FORMATTED',
     *           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(5),ERR=99 )
C
          WRITE(IMPLA2,5019) FICH (II1:II2)
C
        ENDIF
C
C  7) ouverture du fichier .temper + entete fichier case(suite)
C
        IF (IVISTE.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.temper'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(6), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(6),ERR=99 )
C
          WRITE(IMPLA2,5020) FICH (II1:II2)
C
        ENDIF
C
C  8) ouverture du fichier .tempch + entete fichier case(suite)
C
        IF (IVISHP.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.tempch'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(7), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(7),ERR=99 )
C
          WRITE(IMPLA2,5021) FICH (II1:II2)
C
        ENDIF
C
C  9) ouverture du fichier .dck + entete fichier case(suite)
C
        IF (IVISDK.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.dck'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(8), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(8),ERR=99 )
C
          WRITE(IMPLA2,5022) FICH (II1:II2)
C
        ENDIF
C
C  10) ouverture du fichier .mch + entete fichier case(suite)
C
        IF (IVISCH.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.mch'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(9), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(9),ERR=99 )
C
          WRITE(IMPLA2,5023) FICH (II1:II2)
C
        ENDIF
C
C  11) ouverture du fichier .mch + entete fichier case(suite)
C
        IF (IVISCK.EQ.1) THEN
          FICH = NAME
          CALL VERLON ( FICH, II1, II2, LPOS )
          FICH(II2+1:II2+7) = '.mck'
          II2 = II2 + 7
          OPEN ( UNIT=IMPLA5(10), FILE=FICH (II1:II2),
     &           STATUS='UNKNOWN', FORM='FORMATTED',
     &           ACCESS='SEQUENTIAL', IOSTAT=IOS, ERR=99 )
          REWIND ( UNIT=IMPLA5(10),ERR=99 )
C
          WRITE(IMPLA2,5024) FICH (II1:II2)
C
        ENDIF
C
C=======================================================================
C 3.2 FIN DU CALCUL (NFIN=1): ECRITURE SUR FICHIERS FINAUX
C                             ET FERMETURE DU TEMPORAIRE
C                             OU SAUVEGARDE INTERMEDIAIRE
C=======================================================================
C
C 1) ON COMPTE LE NOMBRE TOTAL D'ENREGISTREMENTS
C
        NPT = 0
        NLMAX = 0
        LMAX = 0
        DO NL = 1,NBVIS
          NPT = NPT + NPLIST(NL)
          LMAX = MAX(LMAX,NPLIST(NL))
          IF (NPLIST(NL).GT.NLMAX) NLMAX = NPLIST(NL)
          NUME(NL) = 0
        ENDDO
C
C 2) Allocation Memoire
C
        IFINIA = IDEBIA
        IX     = IDEBRA
        IY     = IX+LMAX
        IZ     = IY+LMAX
        IFINRA = IZ +LMAX
C
        IF (IVISV1.EQ.1) THEN
          IU1L   = IFINRA
          IV1L   = IU1L  +LMAX
          IW1L   = IV1L  +LMAX
          IFINRA = IW1L  +LMAX
        ENDIF
        IF (IVISV2.EQ.1) THEN
          IU2L   = IFINRA
          IV2L   = IU2L  +LMAX
          IW2L   = IV2L  +LMAX
          IFINRA = IW2L  +LMAX
        ENDIF
        IF (IVISTP.EQ.1) THEN
          ITPL = IFINRA
          IFINRA = ITPL +LMAX
        ENDIF
        IF (IVISDM.EQ.1) THEN
          IDML = IFINRA
          IFINRA = IDML +LMAX
        ENDIF
        IF (IVISMP.EQ.1) THEN
          IMPL = IFINRA
          IFINRA = IMPL +LMAX
        ENDIF
        IF (IVISTE.EQ.1) THEN
          ITEL = IFINRA
          IFINRA = ITEL +LMAX
        ENDIF
        IF (IVISHP.EQ.1) THEN
          IHPL = IFINRA
          IFINRA = IHPL +LMAX
        ENDIF
        IF (IVISDK.EQ.1) THEN
          IDCKL= IFINRA
          IFINRA = IDCKL +LMAX
        ENDIF
        IF (IVISCH.EQ.1) THEN
          IMCHL = IFINRA
          IFINRA = IMCHL +LMAX
        ENDIF
        IF (IVISCK.EQ.1) THEN
          IMCKL = IFINRA
          IFINRA = IMCKL +LMAX
        ENDIF
C
        CALL RASIZE('ENSLAG',IFINRA)
C       ===========
C
C 3) ON REMPLIT LES ENTETES DES FICHIERS : geo + variable
C
        REWIND(IMPLA1)
C
        WRITE(IMPLA1,3000)
        WRITE(IMPLA1,3001)
        WRITE(IMPLA1,3002)
        WRITE(IMPLA1,3003)
C
        IF (IVISV1.EQ.1) THEN
          REWIND(IMPLA5(1))
          WRITE(IMPLA5(1),4000)
        ENDIF
C
        IF (IVISV2.EQ.1) THEN
          REWIND(IMPLA5(2))
          WRITE(IMPLA5(2),4001)
        ENDIF
C
        IF (IVISTP.EQ.1) THEN
          REWIND(IMPLA5(3))
          WRITE(IMPLA5(3),4002)
        ENDIF
C
        IF (IVISDM.EQ.1) THEN
          REWIND(IMPLA5(4))
          WRITE(IMPLA5(4),4003)
        ENDIF
C
        IF (IVISMP.EQ.1) THEN
          REWIND(IMPLA5(5))
          WRITE(IMPLA5(5),4004)
        ENDIF
C
        IF (IVISTE.EQ.1) THEN
          REWIND(IMPLA5(6))
          WRITE(IMPLA5(6),4005)
        ENDIF
C
        IF (IVISHP.EQ.1) THEN
          REWIND(IMPLA5(7))
          WRITE(IMPLA5(7),4006)
        ENDIF
C
        IF (IVISDK.EQ.1) THEN
          REWIND(IMPLA5(8))
          WRITE(IMPLA5(8),4007)
        ENDIF
C
        IF (IVISCH.EQ.1) THEN
          REWIND(IMPLA5(9))
          WRITE(IMPLA5(9),4008)
       ENDIF
C
        IF (IVISCK.EQ.1) THEN
          REWIND(IMPLA5(10))
          WRITE(IMPLA5(10),4009)
        ENDIF
C
        DO NL = 1,NBVIS
C
          NP = LISTE(NL)
C
          REWIND(IMPLA3)
C
          IPT = 0
          DO I=1,NPT
C
            READ(IMPLA3) NUML
            READ(IMPLA3) XPL,YPL,ZPL
            IF (IVISV1.EQ.1) READ(IMPLA3) U1L,V1L,W1L
            IF (IVISV2.EQ.1) READ(IMPLA3) U2L,V2L,W2L
            IF (IVISTP.EQ.1) READ(IMPLA3) TPL
            IF (IVISDM.EQ.1) READ(IMPLA3) DML
            IF (IVISMP.EQ.1) READ(IMPLA3) MPL
            IF (IVISTE.EQ.1) READ(IMPLA3) TEL
            IF (IVISHP.EQ.1) READ(IMPLA3) HPL
            IF (IVISDK.EQ.1) READ(IMPLA3) DCKL
            IF (IVISCH.EQ.1) READ(IMPLA3) MCHL
            IF (IVISCK.EQ.1) READ(IMPLA3) MCKL
C
            IF (NUML .EQ. NL ) THEN
C
              IPT = IPT+1

              RA(IX+IPT-1) = XPL
              RA(IY+IPT-1) = YPL
              RA(IZ+IPT-1) = ZPL
              IF (IVISV1.EQ.1) THEN
                RA(IU1L+IPT-1) = U1L
                RA(IV1L+IPT-1) = V1L
                RA(IW1L+IPT-1) = W1L
              ENDIF
              IF (IVISV2.EQ.1) THEN
                RA(IU2L+IPT-1) = U2L
                RA(IV2L+IPT-1) = V2L
                RA(IW2L+IPT-1) = W2L
              ENDIF
              IF (IVISTP.EQ.1) THEN
                RA(ITPL+IPT-1) = TPL
              ENDIF
              IF (IVISDM.EQ.1) THEN
                RA(IDML+IPT-1) = DML
              ENDIF
              IF (IVISMP.EQ.1) THEN
                RA(IMPL+IPT-1) = MPL
              ENDIF
              IF (IVISTE.EQ.1) THEN
                RA(ITEL+IPT-1) = TEL
              ENDIF
              IF (IVISHP.EQ.1) THEN
                RA(IHPL+IPT-1) = HPL
              ENDIF
              IF (IVISDK.EQ.1) THEN
                RA(IDCKL+IPT-1) = DCKL
              ENDIF
              IF (IVISCH.EQ.1) THEN
                RA(IMCHL+IPT-1) = MCHL
              ENDIF
              IF (IVISCK.EQ.1) THEN
                RA(IMCKL+IPT-1) = MCKL
              ENDIF
C
            ENDIF
          ENDDO
C
C  Ecriture Fichier Geometrie
C
          WRITE(IMPLA1,3010)
          WRITE(IMPLA1,1010) NL
          WRITE(IMPLA1,3004) LIST0(NL)
          WRITE(IMPLA1,3005)
          WRITE(IMPLA1,1010) IPT
          DO II=1,IPT
            WRITE(IMPLA1,1030) RA(IX+II-1)
          ENDDO
          DO II=1,IPT
            WRITE(IMPLA1,1030) RA(IY+II-1)
          ENDDO
          DO II=1,IPT
            WRITE(IMPLA1,1030) RA(IZ+II-1)
          ENDDO
          WRITE(IMPLA1,3006)
          WRITE(IMPLA1,1010) IPT-1
          DO II=1,IPT-1
            WRITE(IMPLA1,1020) II,II+1
          ENDDO
C
C  Ecriture Fichiers Variables
C
          IF (IVISV1.EQ.1) THEN
C
            WRITE(IMPLA5(1),3010)
            WRITE(IMPLA5(1),1010) NL
            WRITE(IMPLA5(1),3005)
            DO II=1,IPT
              WRITE(IMPLA5(1),1030) RA(IU1L+II-1)
            ENDDO
            DO II=1,IPT
              WRITE(IMPLA5(1),1030) RA(IV1L+II-1)
            ENDDO
            DO II=1,IPT
              WRITE(IMPLA5(1),1030) RA(IW1L+II-1)
            ENDDO
          ENDIF
C
          IF (IVISV2.EQ.1) THEN
            WRITE(IMPLA5(2),3010)
            WRITE(IMPLA5(2),1010) NL
            WRITE(IMPLA5(2),3005)
            DO II=1,IPT
              WRITE(IMPLA5(2),1030) RA(IU2L+II-1)
            ENDDO
            DO II=1,IPT
              WRITE(IMPLA5(2),1030) RA(IV2L+II-1)
            ENDDO
            DO II=1,IPT
              WRITE(IMPLA5(2),1030) RA(IW2L+II-1)
            ENDDO
          ENDIF
C
          IF (IVISTP.EQ.1) THEN
            WRITE(IMPLA5(3),3010)
            WRITE(IMPLA5(3),1010) NL
            WRITE(IMPLA5(3),3005)
            DO II=1,IPT
              WRITE(IMPLA5(3),1030) RA(ITPL+II-1)
            ENDDO
          ENDIF
          IF (IVISDM.EQ.1) THEN
            WRITE(IMPLA5(4),3010)
            WRITE(IMPLA5(4),1010) NL
            WRITE(IMPLA5(4),3005)
            DO II=1,IPT
              WRITE(IMPLA5(4) ,1030) RA(IDML+II-1)
            ENDDO
          ENDIF
          IF (IVISMP.EQ.1) THEN
            WRITE(IMPLA5(5),3010)
            WRITE(IMPLA5(5),1010) NL
            WRITE(IMPLA5(5),3005)
            DO II=1,IPT
              WRITE(IMPLA5(5),1030) RA(IMPL+II-1)
            ENDDO
          ENDIF
          IF (IVISTE.EQ.1) THEN
            WRITE(IMPLA5(6),3010)
            WRITE(IMPLA5(6),1010) NL
            WRITE(IMPLA5(6),3005)
            DO II=1,IPT
              WRITE(IMPLA5(6),1030) RA(ITEL+II-1)
            ENDDO
          ENDIF
          IF (IVISHP.EQ.1) THEN
            WRITE(IMPLA5(7),3010)
            WRITE(IMPLA5(7),1010) NL
            WRITE(IMPLA5(7),3005)
            DO II=1,IPT
              WRITE(IMPLA5(7),1030) RA(IHPL+II-1)
            ENDDO
          ENDIF
          IF (IVISDK.EQ.1) THEN
            WRITE(IMPLA5(8),3010)
            WRITE(IMPLA5(8),1010) NL
            WRITE(IMPLA5(8),3005)
            DO II=1,IPT
              WRITE(IMPLA5(8),1030) RA(IDCKL+II-1)
            ENDDO
          ENDIF
          IF (IVISCH.EQ.1) THEN
            WRITE(IMPLA5(9),3010)
            WRITE(IMPLA5(9),1010) NL
            WRITE(IMPLA5(9),3005)
            DO II=1,IPT
             WRITE(IMPLA5(9),1030) RA(IMCHL+II-1)
            ENDDO
          ENDIF
          IF (IVISCK.EQ.1) THEN
            WRITE(IMPLA5(10),3010)
            WRITE(IMPLA5(10),1010) NL
            WRITE(IMPLA5(10),3005)
            DO II=1,IPT
              WRITE(IMPLA5(10),1030) RA(IMCKL+II-1)
            ENDDO
          ENDIF

C
        ENDDO
C
        CLOSE(IMPLA1)
        CLOSE(IMPLA2)
        CLOSE(IMPLA3)
        IF (IVISV1.EQ.1) CLOSE(IMPLA5(1))
        IF (IVISV2.EQ.1) CLOSE(IMPLA5(2))
        IF (IVISTP.EQ.1) CLOSE(IMPLA5(3))
        IF (IVISDM.EQ.1) CLOSE(IMPLA5(4))
        IF (IVISMP.EQ.1) CLOSE(IMPLA5(5))
        IF (IVISTE.EQ.1) CLOSE(IMPLA5(6))
        IF (IVISHP.EQ.1) CLOSE(IMPLA5(7))
        IF (IVISDK.EQ.1) CLOSE(IMPLA5(8))
        IF (IVISCH.EQ.1) CLOSE(IMPLA5(9))
        IF (IVISCK.EQ.1) CLOSE(IMPLA5(10))
C
      ENDIF
C
C
      RETURN
C
   99 CONTINUE
      WRITE (NFECRA,9999) FICH (II1:II2), IOS
      CALL CSEXIT(1)
C
C--------
C FORMATS
C--------
C
 1000 FORMAT (A)
 1010 FORMAT (I10)
 1020 FORMAT (I10,I10)
 1030 FORMAT (e12.5)
C
 3000 FORMAT('geometrie trajectoire')
 3001 FORMAT('au format ensight6 : .case')
 3002 FORMAT('node id assign')
 3003 FORMAT('element id assign')
 3004 FORMAT('trajectoire',I10)
 3005 FORMAT('coordinates')
 3006 FORMAT('bar2')
C
 3010 FORMAT('part')
C
 4000 FORMAT('vitesse fluide vu')
 4001 FORMAT('vitesse particules')
 4002 FORMAT('temps de sejour')
 4003 FORMAT('diametre')
 4004 FORMAT('masse')
 4005 FORMAT('temerature')
 4006 FORMAT('temperature')
 4007 FORMAT('dck')
 4008 FORMAT('mch')
 4009 FORMAT('mck')
C
 5010 FORMAT('FORMAT')
 5011 FORMAT('type: ensight gold')
 5012 FORMAT('GEOMETRY')
 5013 FORMAT('model: ',A)
 5014 FORMAT('VARIABLE')
 5015 FORMAT('vector per node: vitesse_fluide_vu       ',A )
 5016 FORMAT('vector per node: vitesse_particules      ',A )
 5017 FORMAT('scalar per node: temps_de_sejour         ',A )
 5018 FORMAT('scalar per node: diametre                ',A )
 5019 FORMAT('scalar per node: masse                   ',A )
 5020 FORMAT('scalar per node: temperature             ',A )
 5021 FORMAT('scalar per node: temperature             ',A )
 5022 FORMAT('scalar per node: dck                     ',A )
 5023 FORMAT('scalar per node: mch                     ',A )
 5024 FORMAT('scalar per node: mck                     ',A )
C
 6000 FORMAT('geometrie trajectoire')
 6001 FORMAT('au format ensight6 (case)')
 6002 FORMAT('node id assign')
 6003 FORMAT('element id assign')
C

C
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    ERREUR D''OUVERTURE SUR LE FICHIER : ',A                 ,/,
     &'@    AVEC UN IOSTAT EGAL A : ',I6                             ,/,
     &'@    (ENSLAG)                                                ',/,
     &'@                                                            ',/,
     &'@  Verifier les numero de fichiers utilises par le Lagrangien',/,
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
