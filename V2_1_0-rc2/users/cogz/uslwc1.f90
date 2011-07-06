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

subroutine uslwc1


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
!     b.Constants of the Libby-Williams Model
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

integer          ipp, idirac

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

! ---- Mean Mixture Fraction
if ( ippmod(icolwc).ge.0 ) then
  ipp = ipprtp(isca(ifm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance of Mixture Fraction
  ipp = ipprtp(isca(ifp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Fuel Mass fraction
  ipp = ipprtp(isca(iyfm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance of Fuel Mass fraction
  ipp = ipprtp(isca(iyfp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if (ippmod(icolwc).ge.2) then
    ipp = ipprtp(isca(icoyfp))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
endif

! ---- Enthalpy
if ( ippmod(icolwc).eq.1 .or.                                     &
     ippmod(icolwc).eq.3 .or.                                     &
     ippmod(icolwc).eq.5    ) then
  ipp = ipprtp(isca(ihm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! b. Variables of State; User definied Variables
!===============================================================================

! --- Source term
  ipp = ipppro(ipproc(itsc))
  NOMVAR(IPP)   = 'T.SOURCE'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Temperature in K
  ipp = ipppro(ipproc(itemp))
  NOMVAR(IPP)   = 'Temperature'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Fuel Mass fraction
  ipp = ipppro(ipproc(iym(1)))
  NOMVAR(IPP)   = 'YM_Fuel'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Oxidizer Mass fraction
  ipp = ipppro(ipproc(iym(2)))
  NOMVAR(IPP)   = 'YM_Oxyd'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Products Mass fraction
  ipp = ipppro(ipproc(iym(3)))
  NOMVAR(IPP)   = 'YM_Prod'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  do idirac = 1, ndirac
    ipp = ipppro(ipproc(irhol(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'RHOL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iteml(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'TEML',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmel(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'FMEL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmal(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'FMAL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iampl(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'AMPL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(itscl(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'TSCL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(imaml(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'MAML',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  enddo

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
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1))

srrom = 0.95d0


!===============================================================================
! 3. Physical Constants
!===============================================================================

! --> DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

! --> Constants of the Libby-Williams Model

! --- Reference velocity
 vref = 60.d0
! --- Reference length scale
 lref = 0.1d0
! --- Activation Temperature
 ta   = 0.2d5
! --- Cross-over Temperature (combustion of propane)
 tstar= 0.12d4

!----
! END
!----

return
end subroutine
