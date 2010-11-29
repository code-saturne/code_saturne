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

subroutine usebu1

!===============================================================================
!  PURPOSE:
!  --------
!  1. Variable Output
!     a. Transported Variables
!     b. Variables of State; User definied Variables
!
!  2. Additional Calculation Options
!     a. Density Relaxation
!
!  3. Physical Constants
!     a.Dynamic Diffusion Coefficient
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
use ppincl
use radiat

!===============================================================================

implicit none

integer          ipp

!===============================================================================
!===============================================================================
! 1. Variable Output
!===============================================================================
!    Function                             |  Key Word |   Indicator
!    ---------------------------------------------------------------
!    Variable Output in the result file   | ICHRVR()  | yes= 1  ; no=0
!    Variable Output in the listing file  | ILISVR()  | yes= 1  ; no=0
!    Output of the temporal evolution of  | IHISVR()  | yes=-1* ; no=0
!    the variable at monitoring points    |           |
!    -----------------------------------------------------------------
!    *: Output for all monitoring points defined in subroutine usini1.f90
!
!===============================================================================
! a. Transported Variables
!===============================================================================
! ---- Mass fraction of unburned (or fresh)  gas
if ( ippmod(icoebu).ge.0 ) then
  ipp = ipprtp(isca(iygfm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! ---- Mean Mixture Fraction
if ( ippmod(icoebu).ge.2 ) then
  ipp = ipprtp(isca(ifm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


! ---- Enthalpy
if ( ippmod(icoebu).eq.1 .or.                                     &
     ippmod(icoebu).eq.3      ) then
  ipp = ipprtp(isca(ihm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! b. Variables of State; User definied Variables
!===============================================================================

! ---- Temperature
ipp = ipppro(ipproc(itemp))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Mean mass fraction of Fuel:    YM_Fuel
ipp = ipppro(ipproc(iym(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Mean mass fraction of Oxidizer : YM_Oxy
ipp = ipppro(ipproc(iym(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Mean mass fraction of Product: YM_Prod
ipp = ipppro(ipproc(iym(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Premixed flame including gas radiation

if ( iirayo.gt.0 ) then

! ---- Absorption Coefficient
  ipp = ipppro(ipproc(ickabs))
  NOMVAR(IPP)   = 'KABS'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^4
  ipp = ipppro(ipproc(it4m))
  NOMVAR(IPP)   = 'TEMP4'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^3
  ipp = ipppro(ipproc(it3m))
  NOMVAR(IPP)   = 'TEMP3'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 2. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0


!===============================================================================
! 3. Physical Constants
!===============================================================================

!       DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

!       cebu: EBU-model constant

 cebu   = 2.5d0


!----
! END
!----

return
end subroutine
