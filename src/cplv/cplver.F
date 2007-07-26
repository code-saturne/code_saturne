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
                        SUBROUTINE CPLVER
C                       *****************
C     -------------------------------------------------------------
     & ( IOK    )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
CFONC    --------------------------------------------------------------
CFONC
CFONC       COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
CFONC       TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON
CFONC
CFONC       VERIFICATION DES PARAMETRES DE CALCUL
CFONC         APRES INTERVENTION UTILISATEUR
CFONC         (COMMONS)
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
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
      INTEGER          IOK
C
C VARIABLES LOCALES
C
      CHARACTER        CHAINE*80
      INTEGER          IPP, ISC, JJ, IPHAS
      INTEGER          IGE
C
C***********************************************************************
C
C=======================================================================
C 1. ENTREES SORTIES entsor.h : formats 1000
C    Variables transportees
C=======================================================================
C
C --> Chrono
C
      DO ISC = 1, NSCAPP
        IPP = IPPRTP(ISCA(ISCAPP(ISC)))
        IF(ICHRVR(IPP).NE.1.AND.ICHRVR(IPP).NE.0) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1000)CHAINE(1:8),IPP,ICHRVR(IPP)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C --> Historiques
C
      DO ISC = 1, NSCAPP
        IPP = IPPRTP(ISCA(ISCAPP(ISC)))
        IF( IHISVR(IPP,1).GT.NCAPT.OR.
     &     (IHISVR(IPP,1).LT.0.AND.IHISVR(IPP,1).NE.-1) ) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1010)CHAINE(1:8),IPP,NCAPT,IHISVR(IPP,1)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
      DO ISC = 1, NSCAPP
        IPP = IPPRTP(ISCA(ISCAPP(ISC)))
        IF( (IHISVR(IPP,1).GT.0.AND.IHISVR(IPP,1).LT.NCAPT)) THEN
          DO JJ = 1, IHISVR(IPP,1)
            IF(IHISVR(IPP,JJ).LE.0.OR.IHISVR(IPP,JJ).GT.NCAPT) THEN
              CHAINE=NOMVAR(IPP)
              WRITE(NFECRA,1020)CHAINE(1:8),IPP,JJ,NCAPT,IHISVR(IPP,JJ)
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C --> Listing
C
      DO ISC = 1, NSCAPP
        IPP = IPPRTP(ISCA(ISCAPP(ISC)))
        IF( ILISVR(IPP).NE.0.AND.ILISVR(IPP).NE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1030)CHAINE(1:8),IPP,ILISVR(IPP)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C
C=======================================================================
C 2. ENTREES SORTIES entsor.h : formats 1000
C    Variables transportees
C=======================================================================
C
C --> Chrono
C
      IPP = IPPPRO(IPPROC(IMMEL))
      IF(ICHRVR(IPP).NE.1.AND.ICHRVR(IPP).NE.0) THEN
        CHAINE=NOMVAR(IPP)
        WRITE(NFECRA,1000)CHAINE(1:8),IPP,ICHRVR(IPP)
        IOK = IOK + 1
      ENDIF
C
      IPP = IPPPRO(IPPROC(ITEMP1))
      IF(ICHRVR(IPP).NE.1.AND.ICHRVR(IPP).NE.0) THEN
        CHAINE=NOMVAR(IPP)
        WRITE(NFECRA,1000)CHAINE(1:8),IPP,ICHRVR(IPP)
        IOK = IOK + 1
      ENDIF
      DO IGE = 1, (NGAZE-2*NCHARB)
        IPP = IPPPRO(IPPROC(IYM1(IGE)))
        IF(ICHRVR(IPP).NE.1.AND.ICHRVR(IPP).NE.0) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1000)CHAINE(1:8),IPP,ICHRVR(IPP)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C --> Historiques
C
      IPP = IPPPRO(IPPROC(IMMEL))
      IF( IHISVR(IPP,1).GT.NCAPT.OR.
     &   (IHISVR(IPP,1).LT.0.AND.IHISVR(IPP,1).NE.-1) ) THEN
        CHAINE=NOMVAR(IPP)
        WRITE(NFECRA,1010)CHAINE(1:8),IPP,NCAPT,IHISVR(IPP,1)
        IOK = IOK + 1
      ENDIF
C
      IPP = IPPPRO(IPPROC(ITEMP1))
      IF( IHISVR(IPP,1).GT.NCAPT.OR.
     &   (IHISVR(IPP,1).LT.0.AND.IHISVR(IPP,1).NE.-1) ) THEN
        CHAINE=NOMVAR(IPP)
        WRITE(NFECRA,1010)CHAINE(1:8),IPP,NCAPT,IHISVR(IPP,1)
        IOK = IOK + 1
      ENDIF
      DO IGE = 1, (NGAZE-2*NCHARB)
        IPP = IPPPRO(IPPROC(IYM1(IGE)))
        IF( IHISVR(IPP,1).GT.NCAPT.OR.
     &     (IHISVR(IPP,1).LT.0.AND.IHISVR(IPP,1).NE.-1) ) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1010)CHAINE(1:8),IPP,NCAPT,IHISVR(IPP,1)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C
      IPP = IPPPRO(IPPROC(IMMEL))
      IF( (IHISVR(IPP,1).GT.0.AND.IHISVR(IPP,1).LT.NCAPT)) THEN
        DO JJ = 1, IHISVR(IPP,1)
          IF(IHISVR(IPP,JJ).LE.0.OR.IHISVR(IPP,JJ).GT.NCAPT) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,1020)CHAINE(1:8),IPP,JJ,NCAPT,IHISVR(IPP,JJ)
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
      IPP = IPPPRO(IPPROC(ITEMP1))
      IF( (IHISVR(IPP,1).GT.0.AND.IHISVR(IPP,1).LT.NCAPT)) THEN
        DO JJ = 1, IHISVR(IPP,1)
          IF(IHISVR(IPP,JJ).LE.0.OR.IHISVR(IPP,JJ).GT.NCAPT) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,1020)CHAINE(1:8),IPP,JJ,NCAPT,IHISVR(IPP,JJ)
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
      DO IGE = 1, (NGAZE-2*NCHARB)
        IPP = IPPPRO(IPPROC(IYM1(IGE)))
        IF( (IHISVR(IPP,1).GT.0.AND.IHISVR(IPP,1).LT.NCAPT)) THEN
          DO JJ = 1, IHISVR(IPP,1)
            IF(IHISVR(IPP,JJ).LE.0.OR.IHISVR(IPP,JJ).GT.NCAPT) THEN
              CHAINE=NOMVAR(IPP)
              WRITE(NFECRA,1020)CHAINE(1:8),IPP,JJ,NCAPT,IHISVR(IPP,JJ)
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C
C --> Listing
C
      IPP = IPPPRO(IPPROC(IMMEL))
      IF( ILISVR(IPP).NE.0.AND.ILISVR(IPP).NE.1) THEN
        CHAINE=NOMVAR(IPP)
        WRITE(NFECRA,1030)CHAINE(1:8),IPP,ILISVR(IPP)
        IOK = IOK + 1
      ENDIF
      IPP = IPPPRO(IPPROC(ITEMP1))
      IF( ILISVR(IPP).NE.0.AND.ILISVR(IPP).NE.1) THEN
        CHAINE=NOMVAR(IPP)
        WRITE(NFECRA,1030)CHAINE(1:8),IPP,ILISVR(IPP)
        IOK = IOK + 1
      ENDIF
      DO IGE = 1, (NGAZE-2*NCHARB)
        IPP = IPPPRO(IPPROC(IYM1(IGE)))
        IF( ILISVR(IPP).NE.0.AND.ILISVR(IPP).NE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1030)CHAINE(1:8),IPP,ILISVR(IPP)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C
C=======================================================================
C 3. OPTIONS DU CALCUL : TABLEAUX DE ppincl.h : formats 2000
C=======================================================================
C
C --> Coefficient de relaxation de la masse volumique
C
      IF( SRROM.LT.0D0 .OR. SRROM.GT.1D0) THEN
        WRITE(NFECRA,2000)'SSROM ', SRROM
        IOK = IOK + 1
      ENDIF
C
C=======================================================================
C 4. TABLEAUX DE cstphy.h et ppthch.F : formats 3000
C=======================================================================
C
      IPHAS = 1
C
C --> Masse volumique
C
      IF( RO0(IPHAS).LT.0D0) THEN
          WRITE(NFECRA,3000)IPHAS,'RO0   ', RO0(IPHAS)
          IOK = IOK + 1
        ENDIF
C
C --> Diffusivite dynamique en kg/(m s) : DIFTL0
C
      IF( DIFTL0.LT.0D0) THEN
        WRITE(NFECRA,3010)'DIFTL0', DIFTL0
        IOK = IOK + 1
      ELSE
        VISLS0(IHM) = DIFTL0
      ENDIF
C
C=======================================================================
C 5. FORMATS VERIFICATION
C=======================================================================
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ICHRVR(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  ICHRVR indique si la variable doit etre incluse dans les  ',/,
     &'@    fichiers de post-traitement                             ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IHISVR(',I10   ,',1) DOIT ETRE UN ENTIER EGAL A -1      ',/,
     &'@      POSITIF OU NUL ET INFERIEUR OU EGAL A NCAPT =',I10     ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IHISVR(I,1) indique les sondes a utiliser pour la variable',/,
     &'@    I (-1 signifiant que toutes sont utilisees)             ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IHISVR(',I10   ,',',I10   ,') DOIT ETRE UN ENTIER       ',/,
     &'@      STRICTEMENT POSITIF ET                                ',/,
     &'@      INFERIEUR OU EGAL A NCAPT = ',I10                      ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  IHISVR(I,j) indique le numero de la jieme sonde a utiliser',/,
     &'@    pour la variable I                                      ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ILISVR(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  ILISVR(I) indique si la variable I sera suivie lors des   ',/,
     &'@    impressions dans le listing                             ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)

 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    ',A6,                            ' DOIT ETRE UN REEL    ',/,
     &'@    COMPRIS ENTRE 0 ET 1                                    ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN)          ',/,
     &'@                                                            ',/,
     &'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uscpi1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C 6. SORTIE
C=======================================================================
C
      RETURN
      END
c@z
