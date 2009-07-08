!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine usppmo
!================


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR
! UTILISATION OU NON D'UNE PHYSIQUE PARTICULIERE


!       UNE SEULE PHYSIQUE PARTICULIERE A LA FOIS.


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
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "pointe.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"
include "ppcpfu.h"
include "atincl.h"

!===============================================================================



!===============================================================================


! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if(1.eq.1) then
  return
endif

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================
! 1.  DECLENCHEMENT DE L UTILISATION D'UNE PHYSIQUE PARTICULIERE
!===============================================================================


! ---- COD3P Flamme de diffusion en chimie complete rapide (3 points)
!        si = -1   modele non utilise
!        si =  0   modele utilise dans les conditions adiabatiques
!        si =  1   modele utilise dans les conditions permeatiques

ippmod(icod3p) = -1


!----- CODEQ Flamme de diffusion en chimie rapide vers l'equilibre
!      ATTENTION : la version CODEQ n'EST PAS OPERATIONNELLE
!      ==========
!        si = -1   modele non utilise
ippmod(icodeq) = -1


!----- COEBU Flamme premelangee en Eddy Break Up
!        si = -1   modele non utilise
!        si =  0   modele utilise dans les conditions adiabatiques
!        si =  1   modele utilise dans les conditions permeatiques (H)
!        si =  2   conditions adiabatiques    avec transport de f
!        si =  3   conditions permeatique (H) avec transport de f

ippmod(icoebu) = -1


!----- COBML premelange avec le modele Bray - Moss - Libby
!      ATTENTION : la version COBML n'EST PAS OPERATIONNELLE
!      ==========
!        si = -1   modele non utilise
ippmod(icobml) = -1


!----- COLWC non parfaitement premelange Libby Williams
!        si = -1   modele non utilise
!        si =  0   modele a 2 pics dans les conditions adiabatiques
!        si =  1   modele a 2 pics dans les conditions permeatiques
!                    (suppose rayonnement)
!        si =  2   modele a 3 pics dans les conditions adiabatiques
!        si =  3   modele a 3 pics dans les conditions permeatiques
!                    (suppose rayonnement)
!        si =  4   modele a 4 pics dans les conditions adiabatiques
!        si =  5   modele a 4 pics dans les conditions permeatiques
!                    (suppose rayonnement)

ippmod(icolwc) = -1


!----- Charbon pulverise avec trois combustibles gazeux et
!        granulometrie.
!        CP3PL Combustible moyen local
!        IPPMOD(ICP...) = 0 : Transport d'H2
!        IPPMOD(ICP...) = 1 : Transport d'H2 + sechage
ippmod(icp3pl) = -1

! Prise en compte de la comb. Heterog par le CO2 : attention il
! faut activer l'option "Equation sur le CO2"

ihtco2 = 0


!----- Charbon pulverise couple lagrangien avec trois combustibles
!      gazeux et granulometrie.
!        IPPMOD(ICPL3C) =-1 : Modele non utilise
!        IPPMOD(ICPL3C) = 0 : Transport d'H2
!        IPPMOD(ICPL3C) = 1 : Transport d'H2 + sechage (non operationnel)

ippmod(icpl3c) = -1


!----- Combustion Fuel
!        IPPMOD(ICFUEL) =-1 : modele non utilise
!        IPPMOD(ICFUEL) = 0 : modele active

ippmod(icfuel) = -1

!----- MODEL NOx : pour l'instant on le met ici mais il
!                          faudrait le deplacer, mais ou?
!      Valable uniquement pour le Fuel

!        IEQNOX = 1 ----> Model NOx

ieqnox = 0

!----- Equation sur YCO2 : pour l'instant on le met ici mais il
!                          faudrait le deplacer, mais ou?
!      Valable pour le charbon et pour le Fuel

!         IEQCO2 = 1 ----> Transport de CO2

ieqco2 = 0

!----- COMPF compressible sans choc
!      ==========
!        si = -1   modele non utilise
!        si = 0    modele active
ippmod(icompf) = -1

!----- VERSIONS ELECTRIQUES
!        Equation de l'energie obligatoire --> |IPPMOD(IEL...)| >= 1
!        + Possibilite de constituants

!       ELJOU : Effet Joule
!        IPPMOD(IELJOU) = 1 : Potentiel reel
!        IPPMOD(IELJOU) = 2 : Potentiel complexe
!        IPPMOD(IELJOU) = 3 : Potentiel reel     + CDL Transfo
!        IPPMOD(IELJOU) = 4 : Potentiel complexe + CDL Transfo

ippmod(ieljou) = -1

!       ELARC : Arc electrique
!        IPPMOD(IELARC) = 1 : Potentiel electrique
!        IPPMOD(IELARC) = 2 : Potentiel electrique +
!                             Potentiel vecteur (=>3D)

ippmod(ielarc) = -1

!       ELION : Mobilite ionique
!        IPPMOD(IELION) = 1 : Potentiel electrique

!       ATTENTION : la version ELION n'EST PAS OPERATIONNELLE
!       ==========

ippmod(ielion) = -1

!----- ATMOS ecoulements atmospheriques
!        si = -1   modele non utilise
!        si = 0    modele active
!        si = 1    atmosphere seche
!        si = 2    atmosphere humide (non operationnelle)
ippmod(iatmos) = -1

!----- Aerorefrigerants (cooling tower)
!        si = -1   non utilise
!        si = 0    active sans modele
!        si = 1    active avec modele de Poppe
!        si = 2    active avec modele de Merkel
ippmod(iaeros) = -1

!===============================================================================
! 2.  CHOIX DU FICHIER THERMOCHIMIE DANS LE CAS DE LA COMBUSTION GAZ
!===============================================================================


!-----Si INDJON=1 on utilise une tabulation ENTH-TEMP calculee par JANAF
!     sinon, l'utilisateur doit lui même fournir sa propre tabulation

indjon = 1

!----
! FORMATS
!----
return
end
