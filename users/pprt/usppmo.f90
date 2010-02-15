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
!      - combustion with gaz / coal / heavy fioul oil
!      - compressible flows
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

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "entsor.h"
include "cstphy.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "ppcpfu.h"

!===============================================================================

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

! --- cod3p: Diffusion flame in the framework of "3 points" rapid complete
! ========== chemistry

!     if = -1   module not activated
!     if =  0   adiabatic conditions
!     if =  1   permeatic conditions

ippmod(icod3p) = -1


! --- coebu: Eddy-Break Up pre-mixed flame
! ==========

!     if = -1   module not activated
!     if =  0   adiabatic conditions at constant richness
!     if =  1   permeatic conditions at constant richness
!     if =  2   adiabatic conditions at variable richness
!     if =  3   permeatic conditions at variable richness

ippmod(icoebu) = -1

! --- colwc: Libby-Williams pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   two-peak model with adiabatic condition
!        if =  1   two-peak model with permeatic condition
!        if =  2   three-peak model with adiabatic condition
!        if =  3   three-peak model with permeatic condition
!        if =  4   four-peak model with adiabatic condition
!        if =  5   four-peak model with permeatic condition

ippmod(icolwc) = -1

! --- cp3pl: Pulverized coal, with three gaseous fuels and granulometry
! ========== CP3PL Combustible moyen local

!        if = -1   module not activated
!        if = 0    H2 transport
!        if = 1    H2 transport + drying

ippmod(icp3pl) = -1

! --- cpl3c: Pulverized coal with Lagrangian coupling, with three gaseous fuels
! ========== and granulometry

!        if = -1   module not activated
!        if = 0    H2 transport
!        if = 1    H2 transport + drying (NOT functional)

ippmod(icpl3c) = -1

! --- cfuel: Heavy fuel oil combustion
! ==========

!        if = -1   module not activated
!        if = 0    module activated

ippmod(icfuel) = -1

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
!        if = 0    electric potential
!        if = 1    electric potential and vector potential (hence 3D modelling)

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

! --- cobml: premelange avec le modele Bray - Moss - Libby
! ==========
!        if = -1   module not activated
ippmod(icobml) = -1

! --- codeq: Diffusion flame  en chimie rapide vers l'equilibre
! ==========
!        if = -1   module not activated
ippmod(icodeq) = -1

! --- elion: Ionic mobility
! ==========
!        if = -1   module not activated
!        if = 1    eletric potential
ippmod(ielion) = -1


!===============================================================================
! 3.  Specific options related to herebefore modules
!===============================================================================

! These options are defined here at the moment, this might change in the future

! --- Enthalpy-Temperature conversion law (for gas combustion modelling)
!       indjon = 0   user-specified
!       indjon = 1   tabulated by JANAF (default)

indjon = 1

! --- NOx modelling (ieqnox = 1)
!       Only compatible with heavy fuel oil combustion

ieqnox = 0

! --- CO2 transport equation (ieqco2 = 1)
!       Compatible with coal and heavy fuel oil combustion

ieqco2 = 0

! --- Heteregoneous combustion by CO2 (ihtco2 = 1)
!       Needs the activation of the CO2 transport equation

ihtco2 = 0

!----
! Formats
!----

!----
! End
!----

return
end subroutine
