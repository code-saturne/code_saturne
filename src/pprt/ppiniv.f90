!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine ppiniv &
!================

 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL POUR
!      LA PHYSIQUE PARTICULIERE

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables


!===============================================================================

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================


!===============================================================================
! 2. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================


! ---> Combustion gaz
!      Flamme de diffusion : chimie 3 points

 if ( ippmod(icod3p).ge.0 ) then
  call d3pini                                                     &
 ( nvar   , nscal  ,                                              &
   dt     )
  endif

! ---> Combustion gaz
!      Flamme de premelange : modele EBU

 if (ippmod(icoebu).ge.0) then
  call ebuini(nvar, nscal, dt)
endif

! ---> Combustion gaz
!      Flamme de premelange : modele LWC

if (ippmod(icolwc).ge.0) then
  call lwcini(nvar, nscal, dt)
endif

! ---> Combustion charbon pulverise

if (ippmod(iccoal).ge.0) then
  call cs_coal_varini(nvar, nscal, dt)
endif

! ---> Combustion charbon pulverise couples Lagrangien

if (ippmod(icpl3c).ge.0) then
  call cplini
endif

! ---> Combustion fuel

if  (ippmod(icfuel).ge.0) then
  call cs_fuel_varini(nvar, nscal, dt)
endif

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1       ) then

  call eliniv(isuite)

  ! give back the hand to the user... in fortran.
  call cs_user_f_initialization(nvar, nscal, dt)

endif

! ---> Ecoulements atmospheriques

if (ippmod(iatmos).ge.0) then

  call atiniv                                                     &
 ( nvar   , nscal  ,                                              &
   dt     )

endif

! ---> Cooling towers

if (ippmod(iaeros).ge.0) then

  call ctiniv(nvar, nscal, dt)

endif

! Gas mixture modelling in presence of noncondensable gases and
! condensable gas as stream.
if (ippmod(igmix).ge.0) then

  call cs_gas_mix_initialization &
  ( nvar   , nscal  ,                                            &
    dt     )

endif

! Compressible
! Has to be called AFTER the gas mix initialization because the
! mixture composition is taken into account in the thermodynamic
! law, if gas mix specific physics is enabled.
if (ippmod(icompf).ge.0) then

  call cfiniv(nvar, nscal, dt)

endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine
