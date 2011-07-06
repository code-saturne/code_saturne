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

subroutine uscpi1
!================


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

!===============================================================================

implicit none

integer          ipp , icla , icha

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  THIS TEST CERTIFY THIS VERY ROUTINE IS USED
!     IN PLACE OF LIBRARY'S ONE
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ BEWARE : STOP during data inlet for pulverised coal     ',/,&
'@    =========                                               ',/,&
'@     THE USER SUBROUTINE uscpi1 have to be modified         ',/,&
'@                                                            ',/,&
'@  The computation will not start                            ',/,&
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

!     - Mean of 3 mixture fraction
!       (C from heterogeneoux oxidation, of char, by O2)
ipp = ipprtp(isca(if3m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Meam of (6 ?) mixture fraction
!       (C from heterogeneous reaction between char and CO2)
if ( ihtco2 .eq. 1) then
  ipp = ipprtp(isca(if3mc2))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!     - Variance of 4 mixture fraction
!       (oxidisers)
ipp = ipprtp(isca(if4p2m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Mean of 5 mixture fraction
!       (water vapor from drying)
if ( ippmod(icp3pl) .eq. 1 ) then
  ipp = ipprtp(isca(if5m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!     - Mass fraction of CO2 or CO (relaxation to equilibrium)

if ( ieqco2 .ge. 1 ) then
  ipp = ipprtp(isca(iyco2))
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
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
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

!     - mass fraction (among gases) of O2
ipp = ipppro(ipproc(iym1(4)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO2
ipp = ipppro(ipproc(iym1(5)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2O
ipp = ipppro(ipproc(iym1(6)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of N2
ipp = ipppro(ipproc(iym1(7)))
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

return

end subroutine
