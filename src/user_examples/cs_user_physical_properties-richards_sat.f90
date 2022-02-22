!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_user_physical_properties-richards_sat.f90
!>
!> \brief Definition of physical variable laws example.
!>
!> See \ref physical_properties for examples.
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

double precision darcy_isotropic_dispersion, molecular_diffusion
double precision velocity_norm
integer          iel, ii, fid
integer          ncelt, icelt, ifcvsl
character*80     fname

integer, allocatable, dimension(:) :: lstcel
double precision, dimension(:), pointer :: capacity, permeability
double precision, dimension(:,:), pointer :: tensor_permeability
double precision, dimension(:), pointer :: cpro_vscalt, saturation
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr

type(gwf_soilwater_partition) :: sorption_scal

!===============================================================================
! Initialisation
!===============================================================================

! Index field for the capacity table (C = grad(theta)/grad(h))
call field_get_val_s_by_name('capacity', capacity)

! Index field for the saturation table (theta)
call field_get_val_s_by_name('saturation', saturation)

! Index field for permeability for the flow part
if (darcy_anisotropic_permeability.eq.0) then
  call field_get_val_s_by_name('permeability', permeability)
else
  call field_get_id('permeability',fid)
  call field_get_val_v(fid, tensor_permeability)
endif

! Index for hydraulic head (H=h+z) and darcian velocity
call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! Example of the physical properties definition of two fully saturated soils
!===============================================================================
! Flow part
!==========

allocate(lstcel(ncel))

!< [richards_flow_soils]
! Set parameters for soil 1
call getcel('SOIL1', ncelt, lstcel)

! Loop on cell of the list
do icelt = 1, ncelt

  ! Cell number
  iel = lstcel(icelt)

  ! Set saturation and capacity (0 if saturated)
  saturation(iel) = 0.6d0
  capacity(iel) = 0.d0

  ! Set permeability
  if (darcy_anisotropic_permeability.eq.0) then
    permeability(iel) = 1.d-1
  else
    tensor_permeability(1,iel) = 1.d-1
    tensor_permeability(2,iel) = 1.d-1
    tensor_permeability(3,iel) = 1.d-2
    tensor_permeability(4:6,iel) = 0.d0
  endif

enddo

! Set parameters for soil 2
call getcel('SOIL2', ncelt, lstcel)

! Loop on cell of the list
do icelt = 1, ncelt

  ! Cell number
  iel = lstcel(icelt)

  ! Set saturation and capacity (0 if saturated)
  saturation(iel) = 0.4d0
  capacity(iel) = 0.d0

  ! Set permeability
  if (darcy_anisotropic_permeability.eq.0) then
    permeability(iel) = 5.d-1
  else
    tensor_permeability(1,iel) = 5.d-1
    tensor_permeability(2,iel) = 5.d-1
    tensor_permeability(3,iel) = 5.d-2
    tensor_permeability(4:6,iel) = 0.d0
  endif

enddo
!< [richards_flow_soils]

!===============
! Transport part
!===============

! Loop on solute
do ii = 1, nscal

  ! Index field for diffusion (dispersion and molecular diffusion) for the transport part
  call field_get_key_int(ivarfl(isca(ii)), kivisl, ifcvsl)
  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_vscalt)
  else
    cpro_vscalt => NULL()
  endif

  ! Index field for soilwater_partition structure
  call field_get_key_struct_gwf_soilwater_partition(ivarfl(isca(ii)), &
                                                  sorption_scal)

  !< [richards_flow_solut]
  ! Definition of the isotropic diffusion (dispersion and moleculer diffusion)
  call getcel ('SOIL1', ncelt, lstcel)
  darcy_isotropic_dispersion = 1.d0
  molecular_diffusion = 1.d-6
  do icelt = 1, ncelt
    iel = lstcel(icelt)
    if (ifcvsl.ge.0) then
      velocity_norm = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
      cpro_vscalt(iel) = darcy_isotropic_dispersion * velocity_norm + &
                         saturation(iel) * molecular_diffusion
    endif
  enddo

  ! Definition of the isotropic diffusion (dispersion and moleculer diffusion)
  call getcel ('SOIL2', ncelt, lstcel)
  darcy_isotropic_dispersion = 0.2d0
  molecular_diffusion = 1.d-8
  do icelt = 1, ncelt
    iel = lstcel(icelt)
    if (ifcvsl.ge.0) then
      velocity_norm = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
      cpro_vscalt(iel) = darcy_isotropic_dispersion * velocity_norm + &
                         saturation(iel) * molecular_diffusion
    endif
  enddo

  !< [richards_flow_solut]

enddo

!===============================================================================

return
end subroutine usphyv
