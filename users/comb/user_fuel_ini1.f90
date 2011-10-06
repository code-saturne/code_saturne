!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine user_fuel_ini1
!========================

!===============================================================================
!  PURPOSE   :
!  ---------

!  USER ROUTINE FOR ALLOCATE COMPUTATION PARAMETERS DEALING WITH FUEL
!    (COMMONS)

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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu

!===============================================================================

implicit none

integer          ipp , icla

!===============================================================================

!===============================================================================
! 1. TRANSPORTED VARIABLES
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every probes defined in usini1


! --> Variables for the mix (carrying gas and coal particles)

!      - Enthalpy

ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! --> Variables for droplets

do icla = 1, nclafu
!       - Fuel mass fraction
  ipp = ipprtp(isca(iyfol(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Number of droplets in mix (1/kg)
  ipp = ipprtp(isca(ing(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Fuel enthalpy (J/kg)
  ipp = ipprtp(isca(ih2(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo


! --> Variables for carrying gas

!       - Mean of 1 mixture fraction (fuel vapor)
ipp = ipprtp(isca(ifvap))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Mean of 3 mixture fraction
!       (carbon from heterogeneous oxidation of char)
ipp = ipprtp(isca(if7m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Variance of 4 mixture fraction (air)
ipp = ipprtp(isca(ifvp2m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - YCO2

if ( ieqco2 .ge. 1 ) then
  ipp = ipprtp(isca(iyco2))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!     - HCN and NO

if ( ieqnox .eq. 1 ) then
  ipp = ipprtp(isca(iyhcn))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(iyno))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ihox))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!===============================================================================
! 2. State variables
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every probes defined in usini1


! --> Variables for the mix (carrying gas and coal particles)

!     - Mean Molar Mass of gases in kg
ipp = ipppro(ipproc(immel))
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! --> Variables for droplets

do icla = 1, nclafu
!       - Droplets' Temperature in K
  ipp = ipppro(ipproc(itemp2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Droplet's Density in kg/m3
  ipp = ipppro(ipproc(irom2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Droplet's Diameter
  ipp = ipppro(ipproc(idiam2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Heat flux (between gases and ICLA class droplets)
  ipp = ipppro(ipproc(ih1hlf(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Evaporation mass flow rate (s-1) < 0
  ipp = ipppro(ipproc(igmeva(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Char combsution mass flow rate
  ipp = ipppro(ipproc(igmhtf(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1
enddo

! --> State variables for carrying gas

!     - Temperature for gases only (not mixed with droplets)
ipp = ipppro(ipproc(itemp1))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!      - Nothing
ipp = ipppro(ipproc(iym1(1)))
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

!     -  Mass fraction of fuel vapor
!          (relative to pure gases : not mixed with droplets ..)
ipp = ipppro(ipproc(iym1(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO
ipp = ipppro(ipproc(iym1(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2S
ipp = ipppro(ipproc(iym1(4)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2
ipp = ipppro(ipproc(iym1(5)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of HCN
ipp = ipppro(ipproc(iym1(6)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of NH3
ipp = ipppro(ipproc(iym1(7)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of O2
ipp = ipppro(ipproc(iym1(8)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO2
ipp = ipppro(ipproc(iym1(9)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2O
ipp = ipppro(ipproc(iym1(10)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of SO2
ipp = ipppro(ipproc(iym1(11)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of N2
ipp = ipppro(ipproc(iym1(12)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Carbone Bilan
ipp = ipppro(ipproc(ibcarbone))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Oxygen bilan
ipp = ipppro(ipproc(iboxygen))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Hydrogen bilan
ipp = ipppro(ipproc(ibhydrogen))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!===============================================================================
! 3. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.7d0


!===============================================================================
! 4. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5
!----
! END
!----

return

end subroutine
