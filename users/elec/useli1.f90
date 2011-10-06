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

subroutine useli1
!================


!===============================================================================
!  Purpose  :
!  -------
!          User subroutines for input of calculation parameters,
!       and to initialize variables used for specific electric models,
!
!
!                 by addition of what is done in  USINI1
!
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
!

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
use ppppar
use ppthch
use ppincl
use elincl

!===============================================================================

implicit none

integer          ipp, iesp , idimve

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0. This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif
!
 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''useli1'' must be completed',/, &
'@     for electric module',/,                                    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. SOLVED VARIABLES
!===============================================================================

!  Chronological output, logging in listing, history output
!       if we do not assign the following array values,
!       default values will be used!
!
!     ichrvr( ) = chonological output (yes 1/no 0)
!     ilisvr( ) = logging in listing (yes 1/no 0)
!     ihisvr( ) = history output (number of probes and their numbers)
!     if ihisvr(.,1)  = -1, output for all probes
!
! --> Current variables for electric modules

! ---- Enthalpy
ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Real component of the electrical potential
ipp = ipprtp(isca(ipotr))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!---- Mass fraction of the different constituants of the phase
if ( ngazg .gt. 1 ) then
  do iesp = 1, ngazg-1
    ipp = ipprtp(isca(iycoel(iesp)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

! --> Specific variables for Joule effect for direct conduction
!     Imaginary component of electrical potential
if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ipp = ipprtp(isca(ipoti))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! --> Specific variables for electric arc in 3D
!     vector potential components
if ( ippmod(ielarc).ge.2 ) then
  do idimve = 1, ndimve
    ipp = ipprtp(isca(ipotva(idimve)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

! --> Ionic conduction module
!     Not available in the present version of the code

!===============================================================================
! 2. Algebric or state variables
!===============================================================================

! ---- Temperature
ipp = ipppro(ipproc(itemp) )
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Electric conductivity
ipp = ipppro(ipproc(ivisls(ipotr)))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Joule effect Power
ipp = ipppro(ipproc(iefjou) )
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Real component of the current density
do idimve = 1, ndimve
  ipp = ipppro(ipproc(idjr(idimve)) )
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo

! ---- Imaginary component of the current density
if ( ippmod(ieljou).eq.4 ) then
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(idji(idimve)) )
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

if ( ippmod(ielarc).ge.1 ) then

! ---- Electromagnetic Forces (Laplace forces)
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(ilapla(idimve)) )
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo

! ---- Absorption oefficient  or Radiative sources term
  if ( ixkabe.gt.0 ) then
    ipp = ipppro(ipproc(idrad) )
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
endif

! ---- Electric charge (volumic)
if ( ippmod(ielion).ge.1 ) then
  ipp = ipppro(ipproc(iqelec) )
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! 3. Calculation options
!===============================================================================

! --> Relaxation coefficient for mass density
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = 0.d0

! --> "Electric variables" scaling (Joule effect or electric arc version)
!      IELCOR = 0 : NO Correction
!      IELCOR = 1 : CORRECTION
ielcor = 0

!     Imposed current intensity (electric arc ) in Amp
!        and Imposed Power (Joule effect for glass melting applications) in Watt
!       These values have to be positive
!
couimp = 0.d0
puisim = 0.d0

!     Initial Potential Difference (positive value)
dpot = 0.d0


!----
! FIN
!----

return
end subroutine
