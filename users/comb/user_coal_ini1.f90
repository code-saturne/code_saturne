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

subroutine user_coal_ini1
!==================================


!===============================================================================
!  PURPOSE   :
!  ---------
!  User's routine to control outing of variables for pulverised coal combustion
!  (these parameters are in COMMON)
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
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl

!===============================================================================

implicit none

integer          ipp , icla , icha

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

! --> Variables for coal particles

do icla = 1, nclacp

!       - Char mass fraction (in class ICLA)
  ipp = ipprtp(isca(ixck(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Coal mass fraction (in class ICLA)
  ipp = ipprtp(isca(ixch(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Number of particles for 1 kg mix (from class ICLA)
  ipp = ipprtp(isca(inp(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Enthalpy J/kg (for class ICLA)
  ipp = ipprtp(isca(ih2(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Water mass fraction (in class ICLA)
  if ( ippmod(icp3pl) .eq. 1 ) then
    ipp = ipprtp(isca(ixwt(icla)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
enddo

! --> Variables for the carrier phase

do icha = 1, ncharb

!       - Mean of 1 mixture fraction
!         (from light volatiles of char ICHA)
  ipp = ipprtp(isca(if1m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Mean of 2 mixture fraction
!         (from heavy volatiles of char ICHA)
  ipp = ipprtp(isca(if2m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

enddo

! ---- Variables propres a la phase continue
  if ( noxyd .ge. 2 ) then
    ipp = ipprtp(isca(if4m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( noxyd .eq. 3 ) then
    ipp = ipprtp(isca(if5m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( ippmod(iccoal) .ge. 1 ) then
    ipp = ipprtp(isca(if6m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  ipp = ipprtp(isca(if7m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  if ( ihtco2 .eq. 1 ) then
    ipp = ipprtp(isca(if8m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( ihth2o .eq. 1 ) then
    ipp = ipprtp(isca(if9m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
!

  ipp = ipprtp(isca(ifvp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!
  if ( ieqco2 .ge. 1 ) then
    ipp = ipprtp(isca(iyco2))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( ieqnox .ge. 1 ) then
    ipp = ipprtp(isca(iyhcn))
    NOMVAR(IPP)  = 'FR_HCN'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(iyno))
    NOMVAR(IPP)  = 'FR_NO'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(ihox))
    NOMVAR(IPP)  = 'Enth_Ox'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

!===============================================================================
! 2. Sate variables
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every probes defined in usini1

! --> State varables for the mix

!     - Mean Molar Mass
ipp = ipppro(ipproc(immel))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! --> State variables for coal particles

do icla = 1, nclacp

!       - Particles' Temperature K (of class ICLA)
  ipp = ipppro(ipproc(itemp2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Particles' Density kg/m3 (of class ICLA)
  ipp = ipppro(ipproc(irom2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Particles' Diameter m (of class ICLA)
  ipp = ipppro(ipproc(idiam2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Rate of coal consumption  (s-1) < 0
!         (for class ICLA)
  ipp = ipppro(ipproc(igmdch(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of light volatiles exhaust (s-1) < 0
!         (for class ICLA)
  ipp = ipppro(ipproc(igmdv1(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of heavy volatile exhaust (s-1) < 0
!         (de la classe ICLA)
  ipp = ipppro(ipproc(igmdv2(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of char oxidation by O2 (s-1) < 0
!         (from class ICLA)
  ipp = ipppro(ipproc(igmhet(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of char gazeification by CO2 (s-1) < 0
!         (from class ICLA)
  if ( ihtco2 .eq. 1 ) then
    ipp = ipppro(ipproc(ighco2(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Rate of char gazeification by H2O (s-1) < 0
!         (from class ICLA)
  if ( ihth2o .eq. 1 ) then
    ipp = ipppro(ipproc(ighh2o(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Rate of drying (s-1) < 0
!         (from class ICLA)
  if ( ippmod(icp3pl) .eq. 1 ) then
    ipp = ipppro(ipproc(igmsec(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Mass fraction (of class ICLA) in mix
  ipp = ipppro(ipproc(ix2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

enddo

! --> State variables for carrier gas phase

!     - Temperature of gas mixture
ipp = ipppro(ipproc(itemp1))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Mass fraction (among gases) of  CHx1m
ipp = ipppro(ipproc(iym1(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CHx2m
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

!===============================================================================
! 3. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.95d0

!===============================================================================
! 4. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5

!----
! END
!----


end subroutine
