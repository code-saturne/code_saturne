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
                        SUBROUTINE VORIMP
C                       *****************
C     -------------------------------------------------------------
     & ( IENT   )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  IMPRESSION DES PARAMETRES DE LA METHODE DES VORTEX
CFONC  APRES INTERVENTION UTILISATEUR DANS usvort
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IENT         ! E  !  -> ! NUMERO DE L'ENTREE                   !
CARGU ! NNENT        ! E  !  -> ! NOMBRE D'ENTREES UTILISANT DES VORTEX!
CARGU ! NVORT(NENTMX)! TE ! ->  ! NOMBRE DE VORTEX A CHAQUE ENTREE     !
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
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IENT
C
C VARIABLES LOCALES
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C On ne peut pas faire ces impressions dans impini car la geometrie
C  des entrees n'est pas encore connue.
C
C=======================================================================
C 1. Methode vortex
C=======================================================================
C
      IPASS = IPASS + 1
C
      WRITE(NFECRA,9010) IENT, NVORT(IENT), ICVOR(IENT), ICAS(IENT)
      WRITE(NFECRA,9020) DIR1(1,IENT),DIR1(2,IENT),DIR1(3,IENT),
     &                   DIR2(1,IENT),DIR2(2,IENT),DIR2(3,IENT),
     &                   DIR3(1,IENT),DIR3(2,IENT),DIR3(3,IENT),
     &                   CEN(1,IENT) ,CEN(2,IENT) ,CEN(3,IENT)
      IF(ICAS(IENT).EQ.1) THEN
        WRITE(NFECRA,9030) ICLVOR(1,IENT),ICLVOR(2,IENT),
     &                     ICLVOR(3,IENT),ICLVOR(4,IENT)
      ENDIF
      WRITE(NFECRA,9040) YMIN(IENT),YMAX(IENT),ZMIN(IENT),ZMAX(IENT)
      IF(ICAS(IENT).EQ.1) THEN
        WRITE(NFECRA,9050) LLY(IENT),LLZ(IENT)
      ENDIF
      IF(ICAS(IENT).EQ.2) THEN
        WRITE(NFECRA,9060) LLD(IENT)
      ENDIF
      WRITE(NFECRA,9070) ITLIVO(IENT)
      IF(ITLIVO(IENT).EQ.1) THEN
        WRITE(NFECRA,9080) TLIMVO(IENT)
      ENDIF
      WRITE(NFECRA,9090) ISGMVO(IENT)
      IF(ISGMVO(IENT).EQ.1) THEN
        WRITE(NFECRA,9100) XSGMVO(IENT)
      ENDIF
      WRITE(NFECRA,9110) IDEPVO(IENT)
      IF(IDEPVO(IENT).EQ.1) THEN
        WRITE(NFECRA,9120) UD(IENT)
      ENDIF
      IF(ICAS(IENT).EQ.1.OR.ICAS(IENT).EQ.2.OR.ICAS(IENT).EQ.3) THEN
        WRITE(NFECRA,9130) NDAT(IENT)
      ELSEIF(ICAS(IENT).EQ.4) THEN
        WRITE(NFECRA,9140) UDEBIT(IENT),KDEBIT(IENT),EDEBIT(IENT)
      ENDIF
C
 9010 FORMAT(
     &'-------------                                                ',/,
     &'  -- Entree : ',I10                                            /,
     &'-------------                                                ',/,
     &'       NVORT  = ',4X,I10,    ' (Nombre de vortex            )',/,
     &'       ICVOR  = ',4X,I10,   ' (Nombre de faces a l''entree  )',/,
     &'       ICAS   = ',4X,I10,    ' (1 : conduite rectangulaire   ',/,
     &'                ',14X,       '  2 : conduite circulaire      ',/,
     &'                ',14X,       '  3 : sans CL mais avec fichier',/,
     &'                ',14X,       '  4 : sans CL ni fichier      )',/)
 9020 FORMAT(
     &' --- Directions principales du repere local                  ',/,
     &'                  ---- X ----    ---- Y ----    ---- Z ----  ',/,
     &'       DIR1  = ', E14.5,' ',E14.5,' ',E14.5                    /,
     &'       DIR2  = ', E14.5,' ',E14.5,' ',E14.5                    /,
     &'       DIR3  = ', E14.5,' ',E14.5,' ',E14.5                    /,
     &'                                                             ',/,
     &' --- Coordonnees du centre de l''entree                      ',/,
     &'       CEN   = ', E14.5,' ',E14.5,' ',E14.5,                   /)
 9030 FORMAT(
     &' --- Conditions aux limites dans le repere local             ',/,
     &' Plan y = -LLY/2 ',4X,I10,    ' (1 : paroi                   ',/,
     &' Plan z =  LLZ/2 ',4X,I10,    '  2 : symetrie                ',/,
     &' Plan y =  LLY/2 ',4X,I10,    '  3 : periodicite            )',/,
     &' Plan z = -LLZ/2 ',4X,I10,    '                              ',/)
 9040 FORMAT(
     &' --- Dimentions de l''entree dans le repere local            ',/,
     &'                  ---- min ----    ---- max ----             ',/,
     &'       Y       = ',E14.5,'  ',E14.5,'                        ',/,
     &'       Z       = ',E14.5,'  ',E14.5,'                        ',/)
 9050 FORMAT(
     &'       LLY     = ',E14.5,    ' (longueur de la conduite dans ',/,
     &'       LLZ     = ',E14.5,    '  les directions DIR1 et DIR2) ',/)
 9060 FORMAT(
     &'       LLD     = ',E14.5,    ' (diametre de la conduite    ) ',/)
 9070 FORMAT(
     &' --- Duree de vie des vortex                                 ',/,
     &'       ITLIVO  = ',4X,I10,   ' (1 : constante                ',/,
     &'                 ',14X,      '  2 : en k^(3/2).U/epsilon   ) ',/)
 9080 FORMAT(
     &'       TLIMVO  = ',E14.5,    ' (1 : duree de vie imposee   ) ',/)
 9090 FORMAT(
     &' --- Taille des vortex                                       ',/,
     &'       ISGMVO  = ',4X,I10,   ' (1 : taille constante         ',/,
     &'                 ',14X,      '  2 : en k^(3/2)/epsilon       ',/,
     &'                 ',14X,      '  2 : en max[nu.k/eps,200.Lk]) ',/)
 9100 FORMAT(
     &'       XSGMVO  = ',E14.5,    ' (1 : taille imposee         ) ',/)
 9110 FORMAT(
     &' --- Marche en temps                                         ',/,
     &'       IDEPVO  = ',4X,I10,   ' (1 : deplacement aleatoire    ',/,
     &'                 ',14X,      '  2 : convection des vortex    ',/)
 9120 FORMAT(
     &'       UD      = ',E14.5,    ' (1 : vit. de deplacement max) ',/)
 9130 FORMAT(
     &' --- Fichier de donnees                                      ',/,
     &'       NDAT    = ',4X,I10,   ' (Nombre de lignes du fichier )',/)
 9140 FORMAT(
     &' --- Donnees a l''entree                                     ',/,
     &'       UDEBIT  = ',E14.5,    ' (vitesse debitante imposee)   ',/,
     &'       KDEBIT  = ',E14.5,    ' (energie cinetique imposee)   ',/,
     &'       EDEBIT  = ',E14.5,    ' (dissipation imposee)         ',/)
C
C=======================================================================
C 2. FIN
C=======================================================================
C
      RETURN
      END
c@z
