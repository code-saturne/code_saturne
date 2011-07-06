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
! Purpose:
! -------

!    User subroutine.

!    Define the use of a specific physics amongst the following:
!      - combustion with gaz / coal / heavy fuel oil
!      - compressible flows
!      - electric arcs
!      - atmospheric modelling
!      - cooling towers modelling

!    Only one specific physics module can be activated at once.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu

!===============================================================================

implicit none


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!===============================================================================

if(1.eq.1) then
  return
endif

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Choice for a specific physics
!===============================================================================

! --- cod3p: Diffusion flame with complete fast chemistry (3 points)
! ==========

!        if = -1   module not activated
!        if =  0   adiabatic model
!        if =  1   extended model with enthalpy source term

ippmod(icod3p) = -1


! --- coebu: Eddy-Break Up pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference Spalding model
!                   (adiabatic, homogeneous mixture fraction)
!        if =  1   extended model with enthalpy source term
!                   (homogeneous mixture fraction : perfect premix)
!        if =  2   extended model with mixture fraction transport
!                   (adiabatic, no variance of mixture fraction)
!        if =  3   extended model with enthalpy and mixture fraction transport
!                   (dilution, thermal losses, etc.)

ippmod(icoebu) = -1

! --- colwc: Libby-Williams pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference two-peak model with adiabatic condition
!        if =  1   extended two-peak model with enthapy source terms
!        if =  2   extended three-peak model, adiabatic
!        if =  3   extended three-peak model with enthalpy source terms
!        if =  4   extended four-peak model, adiabatic
!        if =  5   extended four-peak model with enthalpy source terms

ippmod(icolwc) = -1

! --- cp3pl: Pulverized coal combustion
! ==========

!        Description of granulometry
!        Assumption of diffusion flame around particles
!         (extension of 3-point fast chemistry "D3P")
!        Between a mixture of gaseous fuels (volatiles matters, CO from char
!                                            oxydation)
!            and a mixture of oxidisers (air and water vapor)
!        Enthalpy for both mix and solid phase are solved

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying

ippmod(icp3pl) = -1

! --- cpl3c: Pulverized coal with Lagrangian reciprocal approach
! ==========

!        Not recently tested... at least outdated, may be obsolete

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying (NOT functional)

ippmod(icpl3c) = -1

! --- cfuel: Heavy fuel oil combustion
! ==========

!        Progressive evaporation (temperature gap)
!        Char residue
!        Sulphur tracking

!        if = -1   module not activated
!        if = 0    module activated

ippmod(icfuel) = -1

! --- coal :
! ==========
!
!     Pulverized coal combustion
!        Description of granulometry
!        Assumption of diffusion flame around particles
!         (extension of 3-point fast chemistry "D3P")
!        Between a mixture of gaseous fuels (volatiles matters, CO from char
!                                            oxydation)
!            and a mixture of oxidisers (air and water vapor)
!        Enthalpy for both mix and solid phase are solved
!
!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying

ippmod(iccoal) = -1

! --- compf: Compressible flows
! ==========

!        if = -1   module not activated
!        if = 0    module activated

ippmod(icompf) = -1

! --- eljou: Joule effect
! ==========

!        if = -1   module not activated
!        if = 1    Potentiel reel
!        if = 2    Potentiel complexe
!        if = 3    Potentiel reel     + CDL Transfo
!        if = 4    Potentiel complexe + CDL Transfo

ippmod(ieljou) = -1

! --- elarc: Electric arcs
! ==========

!        if = -1   module not activated
!        if = 1    electric potential
!        if = 2    electric potential and vector potential (hence 3D modelling)

ippmod(ielarc) = -1

! --- atmos: Atmospheric flows
! ==========

!        if = -1   module not activated
!        if = 0    standard modelling
!        if = 1    dry atmosphere
!        if = 2    humid atmosphere (NOT functional)

ippmod(iatmos) = -1

! --- aeros: Cooling towers
! ==========

!        if = -1   module not activated
!        if = 0    no model (NOT functional)
!        if = 1    Poppe's model
!        if = 2    Merkel's model

ippmod(iaeros) = -1


!===============================================================================
! 2.  Specific physics module not available at the moment
!===============================================================================

! WARNING: The following modules ARE NOT functional!
! =======

! --- cobml: Premix model of Bray - Moss - Libby
! ==========

!        if = -1   module not activated

ippmod(icobml) = -1

! --- codeq: Diffusion flame with fast equilibrium chemistry
! ==========

!        if = -1   module not activated

ippmod(icodeq) = -1

! --- elion: Ionic mobility
! ==========

!        if = -1   module not activated

ippmod(ielion) = -1


!===============================================================================
! 3.  Specific options related to herebefore modules
!===============================================================================

! These options are defined here at the moment, this might change in the future

! --- Enthalpy-Temperature conversion law (for gas combustion modelling)

!       if = 0   user-specified
!       if = 1   tabulated by JANAF (default)

indjon = 1

! --- Kinetic model for NOx formation

!         Only compatible with heavy fuel oil combustion

!         if = 0  unused
!         if = 1  activated

ieqnox = 0

! --- Kinetic model for CO <=> CO2

!         Compatible with coal and heavy fuel oil combustion

!         if = 0  unused (maximal conversion in turbulent model)
!         if = 1  transport of CO2 mass fraction
!         if = 2  transport of CO mass fraction

ieqco2 = 0

! --- Heteregoneous combustion by CO2

!         Needs the activation of the CO2 transport equation
!         Account for the reaction between char and CO2: C(s) + CO2 => 2 CO

!         if = 0  unused
!         if = 1  activated

ihtco2 = 0

!----
! Formats
!----

!----
! End
!----

return
end subroutine
