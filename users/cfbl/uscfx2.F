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
                        SUBROUTINE USCFX2
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DES OPTIONS DES VARIABLES POUR LE COMPRESSIBLE SANS CHOC
CFONC    EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1
CFONC
CFONC
CFONC     CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
CFONC     =============================================
CFONC
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
      INTEGER          IPHAS
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C=======================================================================
C
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@                      MODULE COMPRESSIBLE                   ',/,
     &'@                                                            ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uscfx2 DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@     Ce sous-programme utilisateur permet de definir les    ',/,
     &'@       options generales. Il est indispensable.             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1. PROPRIETES PHYSIQUES
C=======================================================================
C
C --> Pour chaque phase
C
      IPHAS = 1
C
C --> Conductivite thermique laminaire
C
C       Conductivite thermique constante : IVISLS = 0
C       Conductivite thermique variable  : IVISLS = 1
C
      IVISLS(ITEMPK(IPHAS)) = 0
C
C       Conductivite thermique de reference :
C
C       VISLS0 = LAMBDA0  (conductivite thermique en W/(m K))
C
C
C       ATTENTION ! IL FAUT QUE VISLS0 SOIT STRICTEMENT POSITIF
C         (donner une valeur meme si la conductivite est variable)
C
      VISLS0(ITEMPK(IPHAS)) = 3.D-2
C
C       Si la conductivite thermique est variable, il faut donner
C       sa loi de variation dans uscfpv.F
C
C
C --> Viscosite en volume
C
C       Viscosite en volume de reference :
C
C       VISCV0 = KAPPA0  (viscosite en volume en kg/(m s))
C       IVISCV = 0 : uniforme en espace et constant en temps
C              = 1 : variable en espace et  en temps
C
      IVISCV(IPHAS) = 0
      VISCV0(IPHAS) = 0.D0
C
C       Si la Viscosite en volume est variable, il faut donner
C       sa loi de variation dans uscfpv.F
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
