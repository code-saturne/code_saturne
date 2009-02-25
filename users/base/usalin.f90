!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

                  subroutine usalin
!================
!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DES PARAMETRES DE DE L'ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================


include "paramx.h"
include "optcal.h"
include "albase.h"


!===============================================================================

! Arguments


! VARIABLES LOCALES


!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) return

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     CE SOUS-PROGRAMME PERMET DE RENSEIGNER LES PARAMETRES

!       SPECIFIQUES AU MODULE ALE


!     IL EST POSSIBLE D'AJOUTER OU DE RETRANCHER DES PARAMETRES


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera ci-dessous des exemples commentes.

!       L'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

! --- Activation de la methode ALE
iale = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

! --- Nombre de sous-iterations d'initialisation du fluide. Dans le cas
!     d'un calcul suite, il s'agit du nombre d'iterations a partir du
!     debut de l'iteration en cours (i.e. pas un nombre absolu).
!     Dans le cas general, NALINF = 0 pour une suite de calcul.

  nalinf = 75

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

! --- Nombre maximal d'iterations d'implicitation du deplacement des
!     structures (=1 pour le couplage explicite)
nalimx = 15

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

! --- Precision relative d'implicitation du deplacement des structures
epalim = 1.d-5

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

! --- Type de viscosite de maillage (cf. usvima)
!     0 : isotrope
!     1 : orthotrope
iortvm = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!----
! FORMATS
!----



return
end

