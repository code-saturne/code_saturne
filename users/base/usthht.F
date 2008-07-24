c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USTHHT
C                       *****************
C     --------------------------------------------------------------
     & ( MODE   , ENTHAL , TEMPER  )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     ROUTINE UTILISATEUR
CFONC     LOI ENTHALPIE   -> TEMPERATURE (MODE =  1)
CFONC     LOI TEMPERATURE -> ENTHALPIE   (MODE = -1)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! MODE         ! E  !  -> !  -1 : T -> H  ;   1 : H -> T         !
CARGU ! ENTHAL       ! R  ! <-> ! ENTHALPIE                            !
CARGU ! TEMPER       ! R  ! <-> ! TEMPERATURE                          !
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
      INCLUDE "period.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          MODE
C
      DOUBLE PRECISION ENTHAL, TEMPER
C
C VARIABLES LOCALES
C
      INTEGER         IT
C
C     Pour le second exemple ci-dessous,
C        donnee de NTAB > 1 valeurs (fictives) de H (=HT) tabulees
C        en fonction de NTAB valeurs de T (=TH) (attention a l'unite K ou C)
C
      INTEGER          NTAB
      PARAMETER       (NTAB=5)
      DOUBLE PRECISION HT(NTAB), TH(NTAB)
      DATA             HT /100000.D0,200000.D0,300000.D0,
     &                               400000.D0,500000.D0 /
      DATA             TH /   100.D0,   200.D0,   300.D0,
     &                                  400.D0,   500.D0 /
C
C***********************************************************************
C
C
C     ATTENTION, CE SOUS-PROGRAMME EST APPELE DANS DES BOUCLES :
C     =========                               ================
C
C       EVITER LES IMPRESSIONS
C       ======
C
C       EVITER LES OPERATIONS FAISANT APPEL AU PARALLELISME
C       ======
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C
      WRITE(NFECRA,9000)
      CALL CSEXIT (1)
C     ===========
C
C----
C FORMATS
C----
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DE LA TEMPERATURE      ',/,
     &'@    *********                                               ',/,
     &'@    LES TABLES ENTHALPIE TEMPERATURE NE SONT PAS DISPONIBLES',/,
     &'@                                                            ',/,
     &'@  Le sous-programme utilisateur usthht doit etre complete.  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Le couplage avec SYRTHES necessite la donne d''une        ',/,
     &'@    temperature de paroi.                                   ',/,
     &'@  Le scalaire choisi pour le couplage SYRTHES est ici une   ',/,
     &'@    enthalpie.                                              ',/,
     &'@  La loi donnant la temperature en fonction de l''enthalpie ',/,
     &'@    doit etre fournie par l''utilisateur dans le            ',/,
     &'@    sous-programme usthht.                                  ',/
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C EXEMPLES
C=======================================================================
C
C Premier exemple, correspondant a H=CpT avec Cp = 4000
C =====================================================
C
C --- Mode H -> T
      IF (MODE .EQ.  1) THEN
        TEMPER = ENTHAL / 4000.D0
C
C --- Mode T -> H
      ELSE
        ENTHAL = TEMPER * 4000.D0
C
      ENDIF
C
      RETURN
C
C
C
C
C
C
C
C
C
C
C
C
C Second exemple, correspondant a une interpolation simple
C === a partir d'une tabulation H=f(T) entree en DATA ====
C     ================================================
C
C --- Mode H -> T
      IF (MODE .EQ.  1) THEN
C
C       Initialisation par defaut
        TEMPER = 0.D0
C
C       Si H plus petit que la plus petite enthalpie tabulee
C         on limite arbitrairement T a la plus petite temperature
        IF(ENTHAL.LE.HT(1)) THEN
          TEMPER = TH(1)
C
C       Si H plus grand que la plus grande enthalpie tabulee
C         on limite arbitrairement T a la plus grande temperature
        ELSEIF(ENTHAL.GE.HT(NTAB)) THEN
          TEMPER = TH(NTAB)
C
C       Sinon, on interpole lineairement
        ELSE
          DO IT = 2, NTAB
            IF(ENTHAL.LE.HT(IT)) THEN
              TEMPER = TH(IT-1)
     &          +(ENTHAL-HT(IT-1))*(TH(IT)-TH(IT-1))/(HT(IT)-HT(IT-1))
            ENDIF
          ENDDO
        ENDIF
C
C --- Mode T -> H
      ELSE
C
C       Initialisation par defaut
        ENTHAL = 0.D0
C
C       Si T plus petit que la plus petite temperature tabulee
C         on limite arbitrairement H a la plus petite enthalpie
        IF(TEMPER.LE.TH(1)) THEN
          ENTHAL = HT(1)
C
C       Si T plus grand que la plus grande temperature tabulee
C         on limite arbitrairement H a la plus grande enthalpie
        ELSEIF(TEMPER.GE.TH(NTAB)) THEN
          ENTHAL = HT(NTAB)
C
C       Sinon, on interpole lineairement
        ELSE
          DO IT = 2, NTAB
            IF(TEMPER.LE.TH(IT)) THEN
              ENTHAL = HT(IT-1)
     &          +(TEMPER-TH(IT-1))*(HT(IT)-HT(IT-1))/(TH(IT)-TH(IT-1))
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
      RETURN
C
C----
C FIN
C----
C
      END
c@z
