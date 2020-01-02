!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!===============================================================================
! Purpose:
! -------

!> \file usvosy.f90
!>
!> \brief Compute a volume exchange coefficient for SYRTHES coupling
!>
!> See \subpage us_vosy for examples.


subroutine usvosy &
!================

 ( inbcou , ncecpl ,                                              &
   iscal  ,                                                       &
   dt     ,                                                       &
   lcecpl , hvol )

!===============================================================================
!> \brief Compute a volume exchange coefficient for SYRTHES coupling
!>
!> The routine is called in \ref cpvosy for each volume coupling
!> therefore it is necessary to test the value of coupling number to separate
!> the treatments of the different couplings
!>
!> Up to now temperature is the only scalar managed for volume couplings.
!>

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode          name          role                                           !
!______________________________________________________________________________!
!> \param[in]    inbcou        SYRTHES coupling number
!> \param[in]    ncecpl        number of cells implied for this coupling
!> \param[in]    iscal         index number of the temperature scalar
!> \param[in]    dt            time step (per cell)
!> \param[in]    lcecpl        list of coupled cells
!> \param[out]   hvol          volume exchange coefficient to compute
!______________________________________________________________________________!


!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

!< [arg]

! Arguments

integer          ncecpl
integer          iscal  , inbcou

integer          lcecpl(ncecpl)

double precision dt(ncelet)
double precision hvol(ncecpl)

!< [arg]

!< [loc_var_dec]

! Local variables

integer          iiscvr, iel, iloc, ifcvsl

double precision cp, mu, lambda, rho, uloc, L, sexcvo
double precision nu, re, pr
double precision hcorr, hvol_cst, lambda_over_cp
double precision, dimension(:), pointer ::  cpro_rom
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cpro_viscl, cpro_viscls, cpro_cp

!< [loc_var_dec]

!===============================================================================

!< [init]

!===============================================================================
! 1. Initialization
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), cvar_vel)

! Cell properties
call field_get_val_s(icrom, cpro_rom)
call field_get_val_s(iviscl, cpro_viscl)
if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
else
  cpro_viscls => NULL()
endif

!< [init]

!< [example_1]

!===============================================================================
! 2. Example 1 of the computation of a volumic exchange coefficient
!
!    hvol(iel) = cst
!
!===============================================================================

hvol_cst = 1.0d6

do iloc = 1, ncecpl  ! Loop on coupled cells
  hvol(iloc) = hvol_cst
enddo

!< [example_1]

!< [example_2]

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

sexcvo = 36.18d0  ! Surface area where exchanges take place by unit of volume
L = 0.03d0        ! Characteristic length

! No test on the coupling number (inbcou). We assume that the same
! treatment is applied to all volume couplings

do iloc = 1, ncecpl  ! Loop on coupled cells

  iel = lcecpl(iloc)

  ! Get cell properties of the current element

  rho = cpro_rom(iel)
  mu = cpro_viscl(iel)

  if (icp.ge.0) then
    cp = cpro_cp(iel)
  else
    cp = cp0
  endif

  if (ifcvsl.ge.0) then ! lambda/Cp is variable
    if (iscacp(iscal).eq.1) then
      lambda =  cpro_viscls(iel)
      lambda_over_cp = lambda/cp
    else
      lambda_over_cp = cpro_viscls(iel)
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

!< [example_2]

! ----------------------------------------------

!----
! End
!----

return
end subroutine usvosy