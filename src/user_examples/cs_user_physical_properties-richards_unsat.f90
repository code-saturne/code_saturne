!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_user_physical_properties-richards_unsat.f90
!>
!> \brief Definition of physical variable laws.
!>
!> See \subpage physical_properties for examples.
!>
!-------------------------------------------------------------------------------

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     )

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use field
use mesh
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

double precision darcy_h
double precision thetar_param, thetas_param, ks_param
double precision kr_param, m_param, n_param, se_param
double precision alpha_param, l_param, tmp_1, tmp_2
double precision darcy_anisotropic_dispersion_l, darcy_anisotropic_dispersion_t
double precision darcy_isotropic_dispersion, molecular_diffusion
double precision velocity_norm, rho
double precision ki, ki_xx, ki_yy, ki_zz, tmp_lt
integer          iel, ii, fid
integer          ncelt, icelt, ifcvsl
character*80     fname

integer, allocatable, dimension(:) :: lstcel
double precision, dimension(:), pointer :: capacity, permeability
double precision, dimension(:), pointer :: saturation, soil_density
double precision, dimension(:), pointer :: cpro_kd, cpro_kplus, cpro_kminus
double precision, dimension(:), pointer :: cpro_mxsol
double precision, dimension(:,:), pointer :: tensor_permeability, visten
double precision, dimension(:), pointer :: cpro_vscalt
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr

type(gwf_soilwater_partition) :: sorption_scal

!===============================================================================

! Index for hydraulic head (H=h+z) and darcian velocity
call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_v(ivarfl(iu), vel)

! Index field for the capacity table (C = grad(theta)/grad(h))
call field_get_val_s_by_name('capacity', capacity)

! Index field for the saturation table (theta)
call field_get_val_s_by_name('saturation', saturation)

! Index field for the soil density table (It is bulk density)
call field_get_val_s_by_name('soil_density', soil_density)

! Index field for tensorial dipersion
call field_get_val_v(ivsten, visten)

! Index field for soilwater_partition structure
call field_get_key_struct_gwf_soilwater_partition(ivarfl(isca(1)), sorption_scal)

! Index field for diffusion (dipersion and molecular diffusion) for the transport part
call field_get_key_int(ivarfl(isca(1)), kivisl, ifcvsl)

! Index field for permeability for the flow part
if (darcy_anisotropic_permeability.eq.0) then
  call field_get_val_s_by_name('permeability', permeability)
else
  call field_get_id('permeability',fid)
  call field_get_val_v(fid, tensor_permeability)
endif

!===============================================================================
! Example of the definition of the physical properties of a single variably saturated soil
!===============================================================================
! Flow part
!==========

!< [richards_set_genuch]
! Set intrinsic permeability (only depends on soil)
if (darcy_anisotropic_permeability.eq.0) then
  ki = 1.d0
else
  ki_xx = 1.d0
  ki_yy = 1.d0
  ki_zz = 1.d-1
endif

! Set values of the Van Genuchten model parameters
ks_param = 0.3d0
thetar_param = 0.078d0
thetas_param = 0.3d0
n_param = 1.56d0
m_param = 1-1/n_param !(Mualem condition)
l_param = 0.5d0
alpha_param = 0.036d0
!< [richards_set_genuch]

! Loop on all cell
do iel = 1, ncel

  !< [richards_set_press]
  ! Switch from hydraulic head (H=h+z) to pressure head (h)
  darcy_h = cvar_pr(iel)
  if (darcy_gravity.eq.1) then
    darcy_h = cvar_pr(iel) - xyzcen(3,iel)
  endif
  !< [richards_set_press]

  !< [richards_sat_part]
  ! Saturated part (h<=0)
  if (darcy_h.ge.0) then

    capacity(iel) = 0.d0
    saturation(iel) = thetas_param

    if (darcy_anisotropic_permeability.eq.0) then
      permeability(iel) = ks_param * ki
    else
      tensor_permeability(1,iel) = ks_param*ki_xx
      tensor_permeability(2,iel) = ks_param*ki_yy
      tensor_permeability(3,iel) = ks_param*ki_zz
      tensor_permeability(4:6,iel) = 0.d0
    endif
  !< [richards_sat_part]

  !< [richards_unsat_part]
  ! Unsaturated part (h<0)
  else

    tmp_1 = abs(alpha_param*darcy_h)**n_param
    tmp_2 = 1.d0 / (1.d0 + tmp_1)
    se_param = tmp_2**m_param

    capacity(iel) = -m_param*n_param*(tmp_1)* &
      (thetas_param-thetar_param)*se_param*tmp_2/darcy_h
    saturation(iel) = thetar_param + se_param*(thetas_param-thetar_param)

    if (darcy_anisotropic_permeability.eq.0) then
      permeability(iel) = ks_param*ki*se_param**l_param*(1.d0-(1.d0-tmp_2)**m_param)**2
    else
      tensor_permeability(1,iel) = ks_param*ki_xx*se_param**l_param*(1.d0-(1.d0-tmp_2)**m_param)**2
      tensor_permeability(2,iel) = ks_param*ki_yy*se_param**l_param*(1.d0-(1.d0-tmp_2)**m_param)**2
      tensor_permeability(3,iel) = ks_param*ki_zz*se_param**l_param*(1.d0-(1.d0-tmp_2)**m_param)**2
      tensor_permeability(4:6,iel) = 0.d0
    endif

  endif
  !< [richards_unsat_part]

enddo

!===============================================================================
! Transport part for one solute with anisotropic dispersion and sorption
!===============================================================================

!< [richards_unsat_trpt_init]
! Set values of the longitudinal and transversal dirpersivity
darcy_anisotropic_dispersion_l = 2.d0
darcy_anisotropic_dispersion_t = 1.d-1

! Set value of the molecular diffusion
molecular_diffusion = 1.d-3
!< [richards_unsat_trpt_init]

!< [richards_unsat_mol_diff]
! Computation of molecular diffusion of the diffusion term
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_vscalt)
  do iel = 1, ncel
    cpro_vscalt(iel) = saturation(iel)*molecular_diffusion
  enddo
else
  cpro_vscalt => NULL()
endif
!< [richards_unsat_mol_diff]

!< [richards_unsat_aniso_disp]
! Computation of the isotropic dispersivity
do iel = 1, ncel

  ! Computation of the norm of the velocity
  velocity_norm = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)

  ! Tensorial dispersion is stored in visten
  tmp_lt = darcy_anisotropic_dispersion_l-darcy_anisotropic_dispersion_t
  visten(1,iel) = darcy_anisotropic_dispersion_t*velocity_norm + tmp_lt*vel(1,iel)**2/(velocity_norm+1.d-15)
  visten(2,iel) = darcy_anisotropic_dispersion_t*velocity_norm + tmp_lt*vel(2,iel)**2/(velocity_norm+1.d-15)
  visten(3,iel) = darcy_anisotropic_dispersion_t*velocity_norm + tmp_lt*vel(3,iel)**2/(velocity_norm+1.d-15)
  visten(4,iel) = tmp_lt*vel(2,iel)*vel(1,iel)/(velocity_norm+1.d-15)
  visten(5,iel) = tmp_lt*vel(2,iel)*vel(3,iel)/(velocity_norm+1.d-15)
  visten(6,iel) = tmp_lt*vel(3,iel)*vel(1,iel)/(velocity_norm+1.d-15)

enddo
!< [richards_unsat_aniso_disp]

!< [richards_unsat_soilwater_partition]
! Set soil density (bulk density!) for delay computation (delay = 1 + soil_density * K_d / saturation)
do iel = 1, ncel
  soil_density(iel) = 1.5d0
enddo

! Get soil-water partition structure
call field_get_key_struct_gwf_soilwater_partition(ivarfl(isca(1)), &
                                                  sorption_scal)

! Index field for kd
call field_get_val_s(sorption_scal%ikd, cpro_kd)

! Index field for EK model parameters (kplus and kminus)
call field_get_val_s(sorption_scal%ikp, cpro_kplus)
call field_get_val_s(sorption_scal%ikm, cpro_kminus)

! Set sorption parameters
do iel=1, ncel
  cpro_kd(iel) = 5.d0
  ! if EK model is chosen, set specific parameters
  cpro_kplus(iel) =  1.d-3
  cpro_kminus(iel) = 1.d-4
enddo

!Index field for cpro_mxsol index (if precipitation option is activated)
call field_get_val_s(sorption_scal%imxsol, cpro_mxsol)

do iel=1, ncel
  cpro_mxsol(iel) = 10.d0
enddo

!< [richards_unsat_soilwater_partition]

!===============================================================================

return
end subroutine usphyv
