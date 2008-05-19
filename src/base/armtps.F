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
                        SUBROUTINE ARMTPS
C                       *****************
C     -------------------------------------------------------------
     & ( NTCABS , NTMABS )
C     -------------------------------------------------------------
C***********************************************************************
c@foncb
CFONC
CFONC FONCTION :
CFONC --------
CFONC           ROUTINE PERMETTANT DE STOPPER LE CALCUL PROPREMENT SI
CFONC             LE TEMPS RESTANT ALLOUE AU PROCESS EST INSUFFISANT
CFONC             UNIQUEMENT POUR VPP ET CLUSTER LINUX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                       ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NTCABS       ! E  !  -> ! NUMERO ABSOLU DU PAS DE TEMPS COURANT!
CARGU ! NTMABS       ! E  !  -> ! NUMERO ABSOLU DU PAS DE TEMPS FINAL  !
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
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NTCABS , NTMABS
C
C VARIABLES LOCALES
C
      INTEGER IRANGS, LNG, ITMP(1), ITLIM
      INTEGER IMETHO
      DATA    IMETHO /-1/
      SAVE    IMETHO
C
      INTEGER NTCAB0
      SAVE    NTCAB0
      DOUBLE PRECISION TREST0, TRESTP, TCPUPR
      SAVE             TREST0, TRESTP, TCPUPR
C
      DOUBLE PRECISION TMOY00, TITPRE, TRESTC, TMOYIT, ALPHA, TITSUP
      DOUBLE PRECISION TMARGE, AA, BB, CC, TMAMIN
      DOUBLE PRECISION TCPUCO
      DOUBLE PRECISION TRESMN, TITSMX
C
C***********************************************************************
C
C     La variable IMETHO indique si l'on a utilise avec succes
C     TREMAIN (1), TCPUMX (2) pour determiner le temps CPU limite,
C     ou s'il ne semble  avoir aucune limite ou que les deux methodes
C     on echoué (0). Avant le premier passage, elle vaut -1.
C
      IF (IMETHO.NE.0) THEN
C
C=======================================================================
C 1. AU PREMIER PASSAGE : INITIALISATIONS
C=======================================================================
C
        IF (IMETHO.EQ.-1) THEN
C
C ---     Premier passage : on essaie d'abord TREMAI
C           Si l'on obtient pas de temps CPU limite, on essaie
C           TCPUMX, qui se base sur la presence de la variable
C           d'environnement CS_MAXTIME.
C
          CALL TREMAI(TREST0, ITLIM)
C         ===========
          IF (ITLIM.EQ.1) THEN
            IMETHO = 1
          ELSE
            CALL TCPUMX(TREST0, ITLIM)
C           ===========
            IF (ITLIM.EQ.1) THEN
              IMETHO = 2
            ENDIF
          ENDIF
C
C         Si une des methodes fonctionne et indique une limite :
C
          IF (IMETHO.GE.0) THEN
C
            NTCAB0 = NTCABS
C
C ---       Temps restant et temps CPU a l'iteration courante
C           (qui sera ensuite la precedente)
C
            TRESTP = TREST0
C
            CALL DMTMPS(TCPUPR)
C           ===========
C
          ENDIF
C
        ELSE
C
C=======================================================================
C 2. TEMPS MOYEN PAR ITERATION
C=======================================================================
C
C --- Temps de l'iteration precedente
C
          CALL DMTMPS(TCPUCO)
C         ===========
          TITPRE = TCPUCO-TCPUPR
C
C
C --- Temps restant courant (temps restant pour le process en cours,
C      le temps des autres processes, compilation etc. non compris)
C --- Temps moyen par iteration depuis le debut
C
          IF (IMETHO.EQ.1) THEN
            CALL TREMAI(TRESTC, ITLIM)
C           ===========
            TMOY00 = (TREST0-TRESTC)/DBLE(NTCABS-NTCAB0)
          ELSE IF (IMETHO.EQ.2) THEN
C           Ici on utilise le temps alloue initialement
            TRESTC = MAX(TREST0 - TCPUCO,0.D0)
            TMOY00 = TCPUCO/DBLE(NTCABS-NTCAB0)
          ENDIF
C
C --- Estimation du temps par iteration
C       ALPHA -> 0 EST PLUS SUR
C
          ALPHA = 0.25D0
          TMOYIT = ALPHA*TITPRE + (1.D0-ALPHA)*TMOY00
C
C --- Temps restant iteration courante (qui sera ensuite la precedente)
C
          TRESTP = TRESTC
C
C --- Temps CPU iteration courante (qui sera ensuite la precedente)
C
          TCPUPR = TCPUCO
C
C=======================================================================
C 3. TEMPS NECESSAIRE POUR UNE ITERATION DE PLUS
C=======================================================================
C
C --- Marge pour les sorties ...
C      100 fois une iteration ou 10% du temps alloue au process (-lt)
C        et au moins 50s ou 1% du temps alloue alloue au process (-lt)
C
C      Soit pour des jobs de
C        moins de    1000 iter     :  10% du temps alloue
C        plus  de    1000 iter et
C          moins de 10000 iter     : 100 fois une iter
C        plus  de   10000 iter     :   1% du temps alloue
C
C
          IF (TMARUS.LT.0.D0) THEN
C
            AA = 100.D0
            BB = 0.10D0
            CC = 0.01D0
            TMAMIN = 50.D0
C
            TMARGE = MIN(TMOYIT*AA,TREST0*BB)
            TMARGE = MAX(TMARGE,TMAMIN)
            TMARGE = MAX(TMARGE,TREST0*CC)
C
          ELSE
C
            TMARGE = TMARUS
C
          ENDIF
C
C
C --- Temps necessaire pour une iteration de plus
C
          TITSUP = TMOYIT + TMARGE
C
C=======================================================================
C 4. TEST (en parallele, le processeur 0 commande)
C=======================================================================
C
C On s'arrete si le temps restant minimum est inferieur au
C   temps maximum estime pour une iter de plus.
C Les notions de min et de max sont relatives aux calculs en parallele.
C   En effet, les temps sont estimes separement par les differents
C   processeurs et il faut qu'ils aient tous le temps de finir.
C   Avec cette methode la marge est la meme pour tous.
C Bien noter que UN SEUL processeur doit decider d'arreter le calcul.
C
          TRESMN = TRESTC
          TITSMX = TITSUP
          IF(IRANGP.GE.0) THEN
            CALL PARMIN(TRESMN)
            CALL PARMAX(TITSMX)
                ENDIF
          IF(IRANGP.LT.0.OR.IRANGP.EQ.0) THEN
            IF (TRESMN.LT.TITSMX) THEN
              NTMABS = NTCABS
              WRITE(NFECRA,1000) NTMABS
            ENDIF
          ELSE
            NTMABS = 0
          ENDIF
C Broadcast
          IF(IRANGP.GE.0) THEN
            IRANGS  = 0
            LNG     = 1
            ITMP(1) = NTMABS
            CALL PARBCI(IRANGS,LNG,ITMP)
            NTMABS = ITMP(1)
          ENDIF
C
          IF (NTCABS.EQ.NTMABS) THEN
            WRITE(NFECRA,1100) TRESTC, TITSUP, TMOY00, TITPRE, TMARGE
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C 5. FORMATS
C=======================================================================
C
 1000 FORMAT(/,
     &'===============================================================',
     &/,'   ** ARRET PAR MANQUE DE TEMPS ',
     &/,'      ------------------------- ',
     &/,'      NOMBRE DE PAS DE TEMPS MAX IMPOSE A NTCABS : ',I10,    /,
     &'===============================================================',
     &                                                                /)
C
 1100 FORMAT(/,
     &'===============================================================',
     &/,'   ** GESTION DU TEMPS RESTANT ',
     &/,'      ------------------------ ',
     &/,'      TEMPS RESTANT ALLOUE AU PROCESS          : ',E14.5,
     &/,'      TEMPS ESTIME POUR UNE ITERATION DE PLUS  : ',E14.5,
     &/,'        DUREE MOYENNE D''UNE ITERATION EN TEMPS : ',E14.5,
     &/,'        DUREE DE L''ITERATION PRECEDENTE        : ',E14.5,
     &/,'        MARGE DE SECURITE                      : ',E14.5,   /,
     &'===============================================================',
     &                                                                /)
C
      END
c@z
