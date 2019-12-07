!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
! Function:
! ---------

!> \file cs_user_boundary_conditions-mapped_inlet.f90
!>
!> Example of cs_f_user_boundary_conditions subroutine.f90 for inlet
!> with inlet profile mapped to profile inside the domain.

!> This example assumes the mesh is orthogonal at the inlet.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use atsoil
use ctincl
use cs_fuel_incl
use mesh
use field
use turbomachinery
use iso_c_binding
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

!< [loc_var_dec]
integer          ifac, iel, ii, ivar, iscal, ilelt, nlfac

integer          keyvar, keysca
integer          n_fields, f_id, normalize, interpolate
double precision xdh, rhomoy
double precision fmprsc, uref2

integer, allocatable, dimension(:) :: lstfac, lstcel
double precision, dimension(:), pointer :: brom

double precision, dimension(:,:), pointer :: vel
double precision, dimension(3,1) :: coord_shift
double precision, dimension(1) :: rvoid

type(c_ptr), save :: inlet_l = c_null_ptr
!< [loc_var_dec]

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstfac(nfabor))  ! temporary array for boundary faces selection

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

call field_get_val_s(ibrom, brom)

fmprsc = 1.d0 ! mean prescribed velocity
!< [init]

!===============================================================================
! Assign a pseudo-periodic channel type inlet to a set of boundary faces,
! using boundary condition mapping.

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - use boundary_conditions_map_set and boundary_conditions_mapped_set
!   to apply a profile from inside the domain to the inlet, renormalizing
!   for some variables.
!
! The impled feedback loop allows progressively reaching a state similar
! to that of a periodic channel at the inlet.
!===============================================================================

!< [example_1_base]
call getfbr('INLET', nlfac, lstfac)
!==========

do ilelt = 1, nlfac

  ifac = lstfac(ilelt)
  iel = ifabor(ifac)

  itypfb(ifac) = ientre

  rcodcl(ifac,iu,1) = fmprsc
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 0.d0

  uref2 =   rcodcl(ifac,iu,1)**2  &
          + rcodcl(ifac,iv,1)**2  &
          + rcodcl(ifac,iw,1)**2
  uref2 = max(uref2,1.d-12)

  !   Turbulence example computed using equations valid for a pipe.

  !   We will be careful to specify a hydraulic diameter adapted
  !     to the current inlet.

  !   We will also be careful if necessary to use a more precise
  !     formula for the dynamic viscosity use in the calculation of
  !     the Reynolds number (especially if it is variable, it may be
  !     useful to take the law from 'cs_user_physical_properties'. Here, we use by default
  !     the 'viscl0" value.
  !   Regarding the density, we have access to its value at boundary
  !     faces (romb) so this value is the one used here (specifically,
  !     it is consistent with the processing in 'cs_user_physical_properties', in case of
  !     variable density)

  !     Hydraulic diameter
  xdh     = 1.d0

  !   Calculation of turbulent inlet conditions using
  !     the turbulence intensity and standard laws for a circular pipe
  !     (their initialization is not needed here but is good practice)

  rhomoy  = brom(ifac)

  call turbulence_bc_inlet_hyd_diam(ifac, uref2, xdh, rhomoy, viscl0,  &
                                   rcodcl)

  ! Handle scalars
  if (nscal.gt.0) then
    do ii = 1, nscal
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1_base]

! Create locator at initialization

!< [example_1_map_init]
if (ntcabs.eq.ntpabs+1) then

  allocate(lstcel(ncel))

  do iel = 1, ncel
    lstcel(iel) = iel
  enddo

  coord_shift(1,1) = 5.95d0
  coord_shift(2,1) = 0.d0
  coord_shift(3,1) = 0.d0

  inlet_l = boundary_conditions_map(MESH_LOCATION_CELLS, ncel,                &
                                    nlfac, lstcel, lstfac,                    &
                                    coord_shift, 0,                           &
                                    0.1d0)

  deallocate(lstcel)

endif
!< [example_1_map_init]

! Subsequent time steps
!----------------------

!< [example_1_map_apply]
if (ntcabs.gt.1) then

  call field_get_n_fields(n_fields)
  call field_get_key_id("variable_id", keyvar)
  call field_get_key_id("scalar_id", keysca)

  interpolate = 0

  do f_id = 0, n_fields-1
    call field_get_key_int(f_id, keyvar, ivar)
    if (ivar.ge.1) then
      call field_get_key_int(f_id, keysca, iscal)
      if (ivar.eq.iu .or. iscal.gt.0) then
        normalize = 1
      else
        normalize = 0
      endif
      call boundary_conditions_mapped_set(f_id, inlet_l, MESH_LOCATION_CELLS, &
                                          normalize, interpolate,             &
                                          nlfac, lstfac, rvoid,               &
                                          nvar, rcodcl)
    endif
  enddo

endif
!< [example_1_map_apply]

! Destroy locator at end

!< [example_1_map_free]
if (ntcabs.eq.ntmabs) then
  call locator_destroy(inlet_l)
endif
!< [example_1_map_free]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstfac)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
