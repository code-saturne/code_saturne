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
                        SUBROUTINE USALIN
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DES PARAMETRES DE DE L'ALE
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
C
      INCLUDE "paramx.h"
      INCLUDE "optcal.h"
      INCLUDE "albase.h"
C
C
C***********************************************************************
C
C ARGUMENTS
C
C
C VARIABLES LOCALES
C
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
      IF(1.EQ.1) RETURN
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     CE SOUS-PROGRAMME PERMET DE RENSEIGNER LES PARAMETRES
C
C       SPECIFIQUES AU MODULE ALE
C
C
C     IL EST POSSIBLE D'AJOUTER OU DE RETRANCHER DES PARAMETRES
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera ci-dessous des exemples commentes.
C
C       L'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

C --- Activation de la methode ALE
      IALE = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C --- Nombre de sous-iterations d'initialisation du fluide. Dans le cas
C     d'un calcul suite, il s'agit du nombre d'iterations a partir du
C     debut de l'iteration en cours (i.e. pas un nombre absolu).
C     Dans le cas general, NALINF = 0 pour une suite de calcul.
C
        NALINF = 75
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C --- Nombre maximal d'iterations d'implicitation du deplacement des
C     structures (=1 pour le couplage explicite)
      NALIMX = 15
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C --- Precision relative d'implicitation du deplacement des structures
      EPALIM = 1.D-5
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C --- Type de viscosite de maillage (cf. usvima)
C     0 : isotrope
C     1 : orthotrope
      IORTVM = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z

