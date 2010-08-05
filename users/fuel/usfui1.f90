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

subroutine usfui1
!================

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
implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.f90"
include "dimens.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "entsor.f90"
include "cstnum.f90"
include "parall.f90"
include "period.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "fuincl.f90"
include "ppincl.f90"
include "ppcpfu.f90"

!===============================================================================

integer          jpp , icla

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  THIS TEST CERTIFY THIS VERY ROUINE IS USED
!     IN PLACE OF LIBRARY'S ONE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ Beware : Stop during DATA Inlet                         ',/,&
'@    =========                                               ',/,&
'@     Heavy Fuel Oil Combustion                              ',/,&
'@     user subroutine USFUI1 must be completed               ',/, &
'@                                                            ',/,&
'@  Computation will be stopped                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

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

jpp = ipprtp(isca(ihm))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

! --> Variables for droplets

do icla = 1, nclafu
!       - Fuel mass fraction
  jpp = ipprtp(isca(iyfol(icla)))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1

!       - Number of droplets in mix (1/kg)
  jpp = ipprtp(isca(ing(icla)))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1

!       - Fuel enthalpy (J/kg)
  jpp = ipprtp(isca(ihlf(icla)))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
enddo


! --> Variables for carrying gas

!       - Mean of 1 mixture fraction (fuel vapor)
jpp = ipprtp(isca(ifvap))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

!     - Mean of 3 mixture fraction
!       (carbon from heterogeneous oxidation of char)
jpp = ipprtp(isca(ifhtf))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

!     - Variance of 4 mixture fraction (air)
jpp = ipprtp(isca(if4p2m))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

!     - YCO2

if ( ieqco2 .ge. 1 ) then
  jpp = ipprtp(isca(iyco2))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
endif

!     - HCN and NO

if ( ieqnox .eq. 1 ) then
  jpp = ipprtp(isca(iyhcn))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
  jpp = ipprtp(isca(iyno))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
  jpp = ipprtp(isca(itaire))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
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
jpp = ipppro(ipproc(immel))
ichrvr(jpp)   = 0
ilisvr(jpp)   = 0
ihisvr(jpp,1) = -1

! --> Variables for droplets

do icla = 1, nclafu
!       - Droplets' Temperature in K
  jpp = ipppro(ipproc(itemp3(icla)))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1

!       - Droplet's Density in kg/m3
  jpp = ipppro(ipproc(irom3(icla)))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1

!       - Droplet's Diameter
  jpp = ipppro(ipproc(idiam3(icla)))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1

!       - Heat flux (between gases and ICLA class droplets)
  jpp = ipppro(ipproc(ih1hlf(icla)))
  ichrvr(jpp)   = 0
  ilisvr(jpp)   = 0
  ihisvr(jpp,1) = -1

!       - Evaporation mass flow rate (s-1) < 0
  jpp = ipppro(ipproc(igmeva(icla)))
  ichrvr(jpp)   = 0
  ilisvr(jpp)   = 0
  ihisvr(jpp,1) = -1

!       - Char combsution mass flow rate
  jpp = ipppro(ipproc(igmhtf(icla)))
  ichrvr(jpp)   = 0
  ilisvr(jpp)   = 0
  ihisvr(jpp,1) = -1
enddo

! --> State variables for carrying gas

!     - Temperature for gases only (not mixed with droplets)
jpp = ipppro(ipproc(itemp1))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!      - Mass fraction of fuel vapor
!          (relative to pure gases : not mixed with droplets ..)
jpp = ipppro(ipproc(iym1(ifov)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Mass fraction of CO ( pure gases)
jpp = ipppro(ipproc(iym1(ico)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Mass fraction of O2 ( same)
jpp = ipppro(ipproc(iym1(io2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Mass fraction of CO2
jpp = ipppro(ipproc(iym1(ico2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Mass fraction of H2O
jpp = ipppro(ipproc(iym1(ih2o)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Mass fraction of N2
jpp = ipppro(ipproc(iym1(in2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Mass fraction of H2S (exhaust form of sulphur during evaporation)
jpp = ipppro(ipproc(iym1(ih2s)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - mass fraction of SO2 (final form of sulphur)
jpp = ipppro(ipproc(iym1(iso2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - MODEL NOX :
if ( ieqnox .eq. 1 ) then
  jpp = ipppro(ipproc(ighcn1))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1
  jpp = ipppro(ipproc(ighcn2))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1
  jpp = ipppro(ipproc(ignoth))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1
endif


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

return

end subroutine
