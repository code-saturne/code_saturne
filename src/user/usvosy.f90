!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine usvosy &
!================

 ( inbcou , ncecpl ,                                              &
   iscal  ,                                                       &
   dt     ,                                                       &
   lcecpl , hvol )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Compute a volume exchange coefficient for SYRTHES coupling

!
! Usage
! -----
! The routine is called in cpvosy() for each volume coupling
! therefore it is necessary to test the value of coupling number to separate
! the treatments of the different couplings
!
! Up to now temperature is the only scalar managed for volume couplings.
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! inbcou           ! i  ! <-- ! SYRTHES coupling number                        !
! ncecpl           ! i  ! <-- ! number of cells implied for this coupling      !
! iscal            ! i  ! <-- ! index number of the temperature scalar         !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! lcecpl(ncecpl)   ! ri ! <-- ! list of coupled cells                          !
! hvol(ncecpl)     ! ra ! --> ! volume exchange coefficient to compute         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use dimens, only: nvar
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          ncecpl
integer          iscal  , inbcou

integer          lcecpl(ncecpl)

double precision dt(ncelet)
double precision hvol(ncecpl)

! Local variables

integer          iiscvr, iel, iloc

double precision cp, mu, lambda, rho, uloc, L, sexcvo
double precision nu, re, pr
double precision hcorr, hvol_cst, lambda_over_cp
double precision, dimension(:), pointer ::  cpro_rom
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cpro_viscl, cpro_vscal, cpro_cp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Initialization
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), cvar_vel)

! Cell properties
call field_get_val_s(icrom, cpro_rom)
call field_get_val_s(iprpfl(iviscl), cpro_viscl)
if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

if (ivisls(iscal).gt.0) then
   call field_get_val_s(iprpfl(ivisls(iscal)), cpro_vscal)
else
   cpro_vscal => NULL()
endif

!===============================================================================
! 2. Example 1 of the computation of a volumic exchange coefficient
!
!    hvol(iel) = cst
!
!===============================================================================

! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return  ! (replace .true. with .false. or remove test to activate)

hvol_cst = 1.0d6

do iloc = 1, ncecpl  ! Loop on coupled cells
   hvol(iloc) = hvol_cst
enddo

!===============================================================================
! 2. Example 2 of the computation of a volumic exchange coefficient
!
!    hvol(iel) =  hsurf(iel) * exchange_surface_by_unit_vol
!
!    with: hsurf = Nusselt * lambda / L
!
!    lambda is the thermal conductivity coefficient
!    L is a characteristic length
!
!    Nusselt is computed by means of the Colburn correlation
!
!    Nu = 0.023 * Re^(0.8) * Pr^(1/3)
!
!    Re is the Reynolds number and Pr is the Prandtl number
!
!===============================================================================

! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return  ! (replace .true. with .false. or remove test to activate)

sexcvo = 36.18d0  ! Surface area where exchanges take place by unit of volume
L = 0.03d0        ! Characteristic length

! No test on the coupling number (inbcou). We assume that the same
! treatment is applied to all volume couplings

do iloc = 1, ncecpl  ! Loop on coupled cells

   iel = lcecpl(iloc)

   ! Get cell properties of the current element

   rho = cpro_rom(iel)
   mu = cpro_viscl(iel)

   if (icp.gt.0) then
      cp = cpro_cp(iel)
   else
      cp = cp0
   endif

   if (ivisls(iscal).gt.0) then ! lambda/Cp is variable
     if (iscacp(iscal).eq.1) then
       lambda =  cpro_vscal(iel)
       lambda_over_cp = lambda/cp
     else
       lambda_over_cp = cpro_vscal(iel)
       lambda =  lambda_over_cp * cp
     endif
   else
     if (iscacp(iscal).eq.1) then
       lambda =  visls0(iscal)
       lambda_over_cp = lambda/cp
     else
       lambda_over_cp = visls0(iscal)
       lambda = lambda_over_cp * cp
     endif
   endif

   ! Compute a local molecular Prandtl **(1/3)

   pr = mu / lambda_over_cp

   ! Compute a local Reynolds number

   uloc = sqrt(cvar_vel(1,iel)**2 + cvar_vel(2,iel)**2 + cvar_vel(3,iel)**2)
   re = max(uloc*rho*L/mu, 1.d0) ! To avoid division by zero

   ! Compute Nusselt number thanks to Colburn correlation

   nu = 0.023d0 * re**0.8d0 * pr**(1.d0/3.d0)
   hcorr = nu * lambda / L

   ! Compute hvol
   hvol(iloc) = hcorr * sexcvo

enddo

! ----------------------------------------------

!----
! End
!----

return
end subroutine usvosy
