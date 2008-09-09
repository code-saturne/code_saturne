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
                        SUBROUTINE COLECD
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
CFONC             RELATIF A LA COMBUSTION GAZ
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
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
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
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
C VARIABLES LOCALES
C
      CHARACTER*150    CHAIN1,CHAIN2
      CHARACTER*12     NOMCOE(NGAZEM)
C
      INTEGER          IT, IGG, IR, IGE, IAT, IOS
      INTEGER          NCGM, NRGM
      INTEGER          INICOE , INICHA
      INTEGER          IDEBCH , IFINCH , LONCH
      INTEGER          ICHAI , ICHCOE
      INTEGER          ATGAZE(NGAZEM,NATOM)
      INTEGER          IGOXY(NRGAZM), IGFUEL(NRGAZM)
      INTEGER          NCOEL
      INTEGER          IICO2 , IIH2O
C
      DOUBLE PRECISION TMIN , TMAX
      DOUBLE PRECISION KABSE(NGAZEM)
      DOUBLE PRECISION COMPOG(NGAZEM,NGAZGM)
      DOUBLE PRECISION EHCOEL (NGAZEM,NPOT) , WMOLCE (NGAZEM)
      DOUBLE PRECISION CPCOEL(NGAZEM)
      DOUBLE PRECISION MFUEL, MREAC, EPSI, NMOLG, BILAN
      DOUBLE PRECISION MOXYD
C
C***********************************************************************
C=======================================================================
C -0. INITIALISATION et VERIFICATIONS
C=======================================================================
C
C
      EPSI = 1.D-9
C
C
      IF (INDJON.NE.1.AND.INDJON.NE.0) THEN
        WRITE(NFECRA,9900) INDJON
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C -1. UTILISATION DE JANAF
C=======================================================================
C
      IF (INDJON.EQ.1) THEN
C
C
C 1.1 LECTURE DU FICHIER DONNEES SPECIFIQUES
C===========================================
C
C --> Ouverture du fichier
C
        OPEN ( UNIT=IMPFPP, FILE=FICFPP,
     &          STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &                                          IOSTAT=IOS, ERR=99 )
        REWIND ( UNIT=IMPFPP,ERR=99 )
C
C --> On se limite pour l'instant a de la combustion gaz
C     Plus exactement a une flamme de diffusion chimie 3 corps
C                  et a une flamme de premelange modele EBU
C     --------------------------------------------------------
C
C --> Lecture des donnees
C
C ---- Nb de constituants gazeux elementaires
C
        READ ( IMPFPP,*,ERR=999,END=999 ) NGAZE
        IF ( NGAZE.GT.NGAZEM ) THEN
          WRITE(NFECRA,9990) NGAZEM,NGAZE
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C ---- Nb de points de tabulation ENTH-TEMP
C
        READ ( IMPFPP,*,ERR=999,END=999 ) NPO
        IF ( NPO.GT.NPOT ) THEN
          WRITE(NFECRA,9991) NPOT, NPO
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C --- Temperature Min et Max
C
        READ (IMPFPP,*,ERR=999,END=999) TMIN
        READ (IMPFPP,*,ERR=999,END=999) TMAX
C
C ---- Lecture des noms des constituants elementaires gazeux
C
        DO IGE=1 , NGAZE
          DO INICOE=1,LEN(NOMCOE(IGE))
            NOMCOE(IGE)(INICOE:INICOE)=' '
          ENDDO
        ENDDO
C
        DO INICHA=1,LEN(CHAIN1)
          CHAIN1(INICHA:INICHA)=' '
        ENDDO
C
        DO INICHA=1,LEN(CHAIN2)
          CHAIN2(INICHA:INICHA)=' '
        ENDDO
C
        READ ( IMPFPP,*,ERR=999,END=999 )
        READ ( IMPFPP,1010,ERR=999,END=999 ) CHAIN1
        CALL VERLON (CHAIN1 , IDEBCH , IFINCH , LONCH)
        CHAIN2(1:LONCH)=CHAIN1(IDEBCH:IFINCH)
C
        IGE   =1
        ICHCOE=0
        DO ICHAI = 1, LONCH
          IF (CHAIN2(ICHAI:ICHAI).NE.' ') THEN
            ICHCOE=ICHCOE+1
            NOMCOE(IGE)(ICHCOE:ICHCOE)=CHAIN2(ICHAI:ICHAI)
          ELSE
            IF (ICHCOE.NE.0) THEN
              IGE=IGE+1
              ICHCOE=0
            ENDIF
          ENDIF
        ENDDO
C
 1010   FORMAT(A150)
C
C
C ---- Rayonnement
C      IRAYPP =  0 : pas de rayonnement
C             =  1 : DOM + calcul a partir des CABS des especes
C             =  2 : DOM + CABS par Modak
C             =  3 : P-1 + calcul a partir des CABS des especes
C             =  4 : P-1 + CABS par Modak
C Le rayonnement n'est autorise qu'avec des modeles permeatiques,
C car en adiabatique l'enthalpie est une grandeur algebrique qui
C ne prend pas en compte les pertes par rayonnement.
C
        READ ( IMPFPP,*,ERR=999,END=999 ) IRAYPP
        IF ( IRAYPP.GT.4 .OR. IRAYPP.LT.0 ) THEN
          WRITE(NFECRA,9992) IRAYPP
          CALL CSEXIT (1)
C         ===========
        ENDIF
        IF ( IRAYPP.GT.0 .AND. IPPMOD(ICOD3P).NE.1
     &       .AND. IPPMOD(ICOEBU).NE.1 .AND. IPPMOD(ICOEBU).NE.3
     &       .AND. IPPMOD(ICOLWC).NE.1 .AND. IPPMOD(ICOLWC).NE.3
     &       .AND. IPPMOD(ICOLWC).NE.5 ) THEN
          WRITE(NFECRA,9993)
     &         IRAYPP,IPPMOD(ICOD3P),IPPMOD(ICOEBU),IPPMOD(ICOLWC)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C ---- Coefficient d'absorption des especes courantes
C
        READ (IMPFPP,*,ERR=999,END=999 )
     &                          ( KABSE(IGE),IGE=1,NGAZE )
C
C ---- Nb especes atomiques (C, H, O, N, ...)
C
        READ (IMPFPP,*,ERR=999,END=999 ) NATO
        IF ( NATO.GT.NATOM ) THEN
          WRITE(NFECRA,9994) NATOM,NATO
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C ---- Masse molaire especes atomiques
C      Composition des especes courantes en fonction des especes
C        elementaires
C
        DO IAT = 1, NATO
          READ (IMPFPP,*,ERR=999,END=999 ) WMOLAT(IAT),
     &                      ( ATGAZE(IGE,IAT),IGE=1,NGAZE )
        ENDDO
C
C ---- Nb especes globales
C
        READ (IMPFPP,*,ERR=999,END=999 ) NGAZG
C       On ne considere qu'UNE SEULE REACTION GLOBALE
C         NGAZG = NCGM = 3 par consequent (F, O, P)
        NCGM = 3
        IF ( NGAZG.NE.NCGM ) THEN
          WRITE(NFECRA,9995) NCGM,NGAZG
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C ---- Composition des especes globales en fonction des especes
C        courantes
C
        DO IGG = 1, NGAZG
          READ (IMPFPP,*,ERR=999,END=999 )
     &                       ( COMPOG(IGE,IGG),IGE=1,NGAZE )
        ENDDO
C
C ----- Nb de reactions globales
C
        READ (IMPFPP,*,ERR=999,END=999 ) NRGAZ
C       On ne considere qu'UNE SEULE REACTION GLOBALE
        NRGM = 1
        IF ( NRGAZ.NE.NRGM ) THEN
          WRITE(NFECRA,9996) NRGM,NRGAZ
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C ---- No des especes concernees par la rapport stoechio
C       Stoechio en especes globales des reactions
C
        DO IR = 1, NRGAZ
          READ (IMPFPP,*,ERR=999,END=999 )
     &        IGFUEL(IR),IGOXY(IR),( STOEG(IGG,IR),IGG=1,NGAZG )
        ENDDO
C
C --> Fermeture du fichier
C
        CLOSE(IMPFPP)
C
C 1.2 CALCULS DE DONNEES COMPLEMENTAIRES
C=======================================
C
C ---- Calcul des masses molaires des especes courantes
C
        DO IGE = 1, NGAZE
          WMOLE(IGE) = 0.D0
          DO IAT = 1, NATO
            WMOLE(IGE)= WMOLE(IGE) + ATGAZE(IGE,IAT)*WMOLAT(IAT)
          ENDDO
        ENDDO
C
C --- Discretisation de la temperature
C
        DO IT = 1, NPO
          TH(IT) = DBLE(IT-1)*(TMAX-TMIN)/DBLE(NPO-1)+TMIN
        ENDDO
C
C ---Calcul des enthalpies par appel a la subroutine PPTBHT
C
        NCOEL  = NGAZE
        DO IGE = 1, NGAZE
          WMOLCE(IGE) = WMOLE (IGE)
        ENDDO
C
        CALL PPTBHT
C       ===========
     & ( NCOEL ,
     &  NOMCOE , EHCOEL , CPCOEL , WMOLCE )
C
        DO IGE = 1 ,NGAZE
          DO IT = 1 , NPO
            EHGAZE(IGE,IT) = EHCOEL (IGE,IT)
          ENDDO
        ENDDO
C
C ---- Calcul des masses molaires des especes globales
C          de la tabulation temperature - enthalpie massique
C          et des coefficients d'absorption des especes globales
c             si RAYONNEMENT
C
        DO IGG = 1 , NGAZG
          WMOLG(IGG) = 0.D0
          NMOLG      = 0.D0
          DO IGE = 1 , NGAZE
            WMOLG(IGG) = WMOLG(IGG)+COMPOG(IGE,IGG)*WMOLE(IGE)
            NMOLG      = NMOLG    +  COMPOG(IGE,IGG)
          ENDDO
          DO IT = 1,NPO
            EHGAZG(IGG,IT) = 0.D0
            IF ( ( IPPMOD(ICOD3P).EQ.1 .OR.
     &             IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 .OR.
     &             IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3 .OR.
     &             IPPMOD(ICOLWC).EQ.5 )
     &           .AND. (IRAYPP.GE.1) ) THEN
              CKABSG(IGG) = 0.D0
            ENDIF
            DO IGE = 1 , NGAZE
              EHGAZG(IGG,IT) = EHGAZG(IGG,IT)
     &             + COMPOG(IGE,IGG)*WMOLE(IGE)*EHGAZE(IGE,IT)
              IF ( ( IPPMOD(ICOD3P).EQ.1 .OR.
     &            IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 .OR.
     &            IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3 .OR.
     &            IPPMOD(ICOLWC).EQ.5 )
     &             .AND. (IRAYPP.GE.1)  ) THEN
                CKABSG(IGG) = CKABSG(IGG)
     &               + COMPOG(IGE,IGG)*KABSE(IGE)*WMOLE(IGE)
              ENDIF
            ENDDO
            EHGAZG(IGG,IT) = EHGAZG(IGG,IT)/WMOLG(IGG)
            IF ( ( IPPMOD(ICOD3P).EQ.1 .OR.
     &          IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 .OR.
     &          IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3 .OR.
     &          IPPMOD(ICOLWC).EQ.5 )
     &           .AND. (IRAYPP.GE.1) ) THEN
              CKABSG(IGG) = CKABSG(IGG)/WMOLG(IGG)
            ENDIF
          ENDDO
          WMOLG(IGG) = WMOLG(IGG) / NMOLG
          DO IGE = 1 , NGAZE
            COMPOG(IGE,IGG) = COMPOG(IGE,IGG) / NMOLG
          ENDDO
        ENDDO
C
C ---- Calcul des coefficients molaires XCO2 , XH2O
C
        DO IGE = 1 , NGAZE
          IF (NOMCOE(IGE)(1:3).EQ.'CO2') IICO2=IGE
          IF (NOMCOE(IGE)(1:3).EQ.'H2O') IIH2O=IGE
        ENDDO
C
        XCO2 =  COMPOG(IICO2,3)
        XH2O =  COMPOG(IIH2O,3)
C
C ---- Calcul bilan pour verification
C        et taux de melange a la stochio pour chaque reaction
C
        DO IR = 1, NRGAZ
          DO IAT = 1, NATO
            BILAN = 0.D0
            DO IGG = 1, NGAZG
              DO IGE = 1, NGAZE
                BILAN = BILAN
     &              + STOEG(IGG,IR)*COMPOG(IGE,IGG)*ATGAZE(IGE,IAT)
              ENDDO
            ENDDO
            IF ( ABS(BILAN) .GT. EPSI ) THEN
              WRITE(NFECRA,9997) IR, IAT
              CALL CSEXIT (1)
C             ===========
            ENDIF
          ENDDO
          MFUEL = STOEG(IGFUEL(IR),IR)*WMOLG(IGFUEL(IR))
          MOXYD = STOEG(IGOXY(IR),IR)*WMOLG(IGOXY(IR))
          MREAC = MFUEL + MOXYD
          FS(IR) = MFUEL/MREAC
        ENDDO
C
C ---- Calcul des coefficients de la fraction massique d'oxydant
C      du programme pdflwc
C
        COEFF1 = COMPOG(2,2)*WMOLE(2)/WMOLG(2)
C
        COEFF3 =   MOXYD / MFUEL
C
        COEFF2 = COEFF3*COEFF1
C
C=======================================================================
C -2. UTILISATION D'UNE TABULATION ENTHALPIE-TEMPERATURE
C=======================================================================
C
      ELSE
C
        OPEN ( UNIT=IMPFPP, FILE=FICFPP,
     &          STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &                                          IOSTAT=IOS, ERR=99 )
        REWIND ( UNIT=IMPFPP,ERR=99 )
C
C ---- Nb de points de tabulation ENTH-TEMP
C
        READ ( IMPFPP,*,ERR=999,END=999 ) NPO
        IF ( NPO.GT.NPOT ) THEN
          WRITE(NFECRA,9991) NPOT, NPO
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C --- Tabulation ENTH-TEMP pour les especes globales
C
        DO IT = 1, NPO
          READ (IMPFPP,*,ERR=999,END=999) TH(IT),
     &                 EHGAZG(1,IT),EHGAZG(2,IT),EHGAZG(3,IT)
        ENDDO
C
C       On ne considere qu'UNE SEULE REACTION GLOBALE
C         NGAZG = NCGM = 3 par consequent (F, O, P)
        NGAZG = 3
C
C       On ne considere qu'UNE SEULE REACTION GLOBALE
        NRGAZ = 1
C
C --- Masses molaires pour les especes globales
C
        READ (IMPFPP,*,ERR=999,END=999) WMOLG(1),WMOLG(2),WMOLG(3)
C
C --- Fraction de melange a la stoechiometrie
C
        READ (IMPFPP,*,ERR=999,END=999) FS(1)
C
C --- Rayonnement
C
        READ (IMPFPP,*,ERR=999,END=999) IRAYPP
        IF ( IRAYPP.GT.4 .OR. IRAYPP.LT.0 ) THEN
          WRITE(NFECRA,9992) IRAYPP
          CALL CSEXIT (1)
C         ===========
        ENDIF
        IF ( IRAYPP.GT.0 .AND. IPPMOD(ICOD3P).NE.1
     &       .AND. IPPMOD(ICOEBU).NE.1 .AND. IPPMOD(ICOEBU).NE.3
     &       .AND. IPPMOD(ICOLWC).NE.1 .AND. IPPMOD(ICOLWC).NE.3
     &       .AND. IPPMOD(ICOLWC).NE.5 ) THEN
          WRITE(NFECRA,9993)
     &         IRAYPP,IPPMOD(ICOD3P),IPPMOD(ICOEBU),IPPMOD(ICOLWC)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C --- Coefficients d'absorption des especes globales
C
        READ (IMPFPP,*,ERR=999,END=999) CKABSG(1),CKABSG(2),CKABSG(3)
C
C --- Coefficients molaires de CO2 et H2O dans les produits
C      (pour le rayonnement)
C
        READ (IMPFPP,*,ERR=999,END=999) XCO2, XH2O
C

C ---> Fermeture du fichier
C
        CLOSE (IMPFPP)
C
C
C ---- Calcul des coefficients de la fraction massique d oxydant
C      on considère que l'oxydant est un mélange d'O2 et N2
C
        COEFF1 = ((WMOLG(2)-0.028)/(0.032-0.028))* (0.032/WMOLG(2))
C
        COEFF3 = (1-FS(1))/FS(1)
C
        COEFF2 = COEFF3*COEFF1
C
      ENDIF
C
C
      RETURN
C
C====================
C 3. SORTIE EN ERREUR
C====================
C
  99  CONTINUE
      WRITE ( NFECRA,9998 )
      CALL CSEXIT (1)
C     ===========
C
  999 CONTINUE
      WRITE ( NFECRA,9999 )
      CALL CSEXIT (1)
C     ===========
C
C
C--------
C FORMATS
C--------
C
C
 9900 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  L''indicateur INDJON doit etre un entier de valeur        ',/,
     &'@    1 (utilisation de Janaf) ou 0 (tabulation utilisateur). ',/,
     &'@                                                            ',/,
     &'@  Il vaut ici ',I10                                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier INDJON dans usppmo.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9990 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes courantes doit etre ',I10             ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9991 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Le nombre de points de tabulation est limite a ',I10       ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9992 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  IRAYPP doit etre un entier egal a 0, 1,  2, 3 ou 4        ',/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9993 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  LE RAYONNEMENT NE PEUT ETRE ACTIVE QU''AVEC UN MODELE DE  ',/,
     &'@   COMBUSTION EN CONDITIONS PERMEATIQUES.                   ',/,
     &'@                                                            ',/,
     &'@  Le fichier parametrique specifie un mode de rayonnement   ',/,
     &'@   IRAYPP = ',I10                                            ,/,
     &'@  Or dans usppmo on a :                                     ',/,
     &'@   IPPMOD(ICOD3P) = ',I10                                    ,/,
     &'@   IPPMOD(ICOEBU) = ',I10                                    ,/,
     &'@   IPPMOD(ICOLWC) = ',I10                                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique et usppmo.               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9994 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes elementaires est limite a ',I10       ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9995 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes globales doit etre ',I10              ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9996 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Le nombre de reactions globales doit etre ',I10            ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9997 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Probleme de conservation rencontre dans la                ',/,
     &'@   reaction ',I10   ,' pour l''element ',I10                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9998 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Erreur a l''ouverture du fichier parametrique.            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,
     &'@                                                            ',/,
     &'@  Erreur a la lecture du fichier parametrique.              ',/,
     &'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,
     &'@    ou son format inadapte.                                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      END
C
C
c@z
