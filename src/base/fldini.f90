!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2019 EDF S.A., France

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

subroutine fldini
!================

!===============================================================================
! Purpose:
! --------

! Define main fields

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________.____._____.________________________________________________.

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use optcal
use cstphy
use numvar
use entsor
use pointe
use albase
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use lagran
use ihmpre
use radiat
use cplsat
use mesh
use post
use field
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, ivar
integer          iflid, kcvlim, ifctsl, clip_id
integer          kturt, kfturt, kislts, keyvar, kclipp, kfturt_alpha
integer          itycat, ityloc, idim1, idim3, idim6
logical          iprev, inoprv
integer          f_id, kscavr, f_type
integer          iopchr, ilog
integer          iscdri, icla, iclap
integer          keyccl, keydri
integer          idfm, iggafm, nfld
integer          iflidp, idimf, n_fans
integer          f_dim

character(len=80) :: name, f_name

type(gas_mix_species_prop) sasp
type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! Initialisation
!===============================================================================

! The itycat variable is used to define field categories. It is used in Fortran
! code with hard-coded values, but in the C API, those values are based on
! (much clearer) category mask definitions in cs_field.h.

itycat = FIELD_INTENSIVE + FIELD_VARIABLE  ! for most variables
ityloc = 1 ! variables defined on cells
idim1  = 1
idim3  = 3
idim6  = 6
iprev = .true.     ! variables have previous value
inoprv = .false.   ! variables have no previous value

call field_get_key_id('log', keylog)
call field_get_key_id('post_vis', keyvis)
call field_get_key_id('label', keylbl)
call field_get_key_id("variable_id", keyvar)

! If a scalar is a variance, store the id of the parent scalar
call field_get_key_id("first_moment_id", kscavr)

! Keys not stored globally
call field_get_key_id('turbulent_flux_model', kturt)
call field_get_key_id('turbulent_flux_id', kfturt)
call field_get_key_id('alpha_turbulent_flux_id', kfturt_alpha)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Number of fields
call field_get_n_fields(nfld)

!===============================================================================
! Set keywords and add some additional fields
!===============================================================================

! User variables
!---------------

idfm = 0
iggafm = 0

do ii = 1, nscal
  if (isca(ii) .gt. 0) then

    ivar = isca(ii)
    call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

    if (ityturt(ii).gt.0) then
      if (ityturt(ii).eq.3) then
        idfm = 1
      endif

      ! GGDH or AFM on current scalar
      ! and if DFM, GGDH on the scalar variance
      iggafm = 1

    ! If the user has chosen a tensorial diffusivity
    else if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
      idfm = 1
    endif

    ! Reference diffusivity value
    call field_set_key_double(ivarfl(ivar), kvisl0, visls0(ii))

    ! Additional fields for Drift scalars is done in addfld

  endif
enddo

! Reserved fields whose ids are not saved (may be queried by name)
!-----------------------------------------------------------------

itycat = FIELD_INTENSIVE
ityloc = 1 ! cells

if (iphydr.eq.1) then
  call field_find_or_create('volume_forces', &
                            itycat, ityloc, idim3, f_id)

  call field_set_key_int(f_id, keylog, 1)
  call field_set_key_int(f_id, keyvis, POST_ON_LOCATION)
endif

if (iphydr.eq.2) then
  call field_find_or_create('hydrostatic_pressure_prd', &
                            itycat, ityloc, idim1, f_id)
endif

! friction velocity at the wall, in the case of a LES calculation
! with van Driest-wall damping (delayed here rather than placed in
! addfld, as idries may be set in modini).

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces

if (itytur.eq.4 .and. idries.eq.1) then
  call field_find_or_create('ustar', itycat, ityloc, idim1, f_id)
endif

! Interior mass flux field
!-------------------------

itycat = FIELD_EXTENSIVE + FIELD_PROPERTY
ityloc = 2 ! inner faces

! Mass flux for the class on interior faces
f_name = 'inner_mass_flux'

if (istmpf.ne.1) then
  call field_create(f_name, itycat, ityloc, idim1, iprev, f_id)
else
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
endif

! The same mass flux for every variable, an other mass flux
! might be defined hereafterwards
do ivar = 1, nvar
  call field_set_key_int(ivarfl(ivar), kimasf, f_id)
enddo

! Boundary Mass flux field
!-------------------------

itycat = FIELD_EXTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces

! Mass flux for the class on interior faces
f_name = 'boundary_mass_flux'
if (istmpf.ne.1) then
  call field_create(f_name, itycat, ityloc, idim1, iprev, f_id)
else
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
endif

! The same mass flux for every variable, an other mass flux
! might be defined hereafterwards
do ivar = 1, nvar
  call field_set_key_int(ivarfl(ivar), kbmasf, f_id)
enddo

! Add mass flux for scalar with a drift (one mass flux per class)
!----------------------------------------------------------------

do iflid = 0, nfld-1

  call field_get_key_int(iflid, keydri, iscdri)

  if (btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

    call field_get_name(iflid, name)
    ! Index of the class, all member of the class share the same mass flux
    call field_get_key_int(iflid, keyccl, icla)

    ! The following fields are not solved, they are properties
    itycat = FIELD_PROPERTY
    ityloc = 2 ! variables defined on interior faces

    ! Mass flux for the class on interior faces
    f_name = 'inner_mass_flux_'//trim(name)
    call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)

    ! Set the inner mass flux index
    call field_set_key_int(iflid, kimasf, f_id)

    ! If the scalar is the representant of a class, then
    ! set the mass flux index to all members of the class
    if (icla.ne.0) then
      do ii = 0, nfld-1
        call field_get_key_int(ii, keyccl, iclap)
        if (icla.eq.iclap) then
          call field_set_key_int(ii, kimasf, f_id)
        endif
      enddo
    endif

    ityloc = 3 ! variables defined on boundary faces

    ! Mass flux for the class on boundary faces
    f_name = 'boundary_mass_flux_'//trim(name)
    call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)

    ! Set the boundary mass flux index
    call field_set_key_int(iflid, kbmasf, f_id)

    ! If the scalar is the representant of a class, then
    ! set the mass flux index to all members of the class
    if (icla.ne.0) then
      do ii = 0, nfld-1
        call field_get_key_int(ii, keyccl, iclap)
        if (icla.eq.iclap) then
          call field_set_key_int(ii, kbmasf, f_id)
        endif
      enddo
    endif

    ityloc = 1 ! variables defined on cells

    ! Get the scalar's output options
    ! (except non-reconstructed boundary output)
    call field_get_key_int(iflid, keyvis, iopchr)
    call field_get_key_int(iflid, keylog, ilog)
    if (iand(iopchr, POST_BOUNDARY_NR) .ne. 0) iopchr = iopchr - POST_BOUNDARY_NR

    ! If the mass flux is imposed, no need of drift_tau nor drift_vel
    if (.not.(btest(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX))) then

      ! Relaxation time
      f_name = 'drift_tau_'//trim(name)
      call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
      call field_set_key_str(f_id, keylbl, f_name)

      ! Set the same visualization options as the scalar,
      call field_set_key_int(f_id, keyvis, iopchr)
      call field_set_key_int(f_id, keylog, ilog)

      ! Store the drift velocity
      f_name = 'drift_vel_'//trim(name)
      call field_create(f_name, itycat, ityloc, idim3, inoprv, f_id)
      call field_set_key_str(f_id, keylbl, f_name)

      ! Set the same visualization options as the scalar
      call field_set_key_int(f_id, keyvis, iopchr)
      call field_set_key_int(f_id, keyvis, ilog)

    endif

    ! Interaction time particle--eddies
    if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS)) then
      f_name = 'drift_turb_tau_'//trim(name)
      call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
      call field_set_key_str(f_id, keylbl, f_name)

      ! Set the same visualization options as the scalar
      call field_set_key_int(f_id, keyvis, iopchr)
      call field_set_key_int(f_id, keyvis, ilog)
    endif

  endif
enddo

! Add weight field for variable to compute gradient
iflidp = -1
itycat = FIELD_PROPERTY
ityloc = 1         ! variables defined on cells
idimf  = -1        ! Field dimension


do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then

    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    if (vcopt%iwgrec.eq.1 .and. vcopt%idiff .ge. 1) then

      call field_get_name(f_id, name)
      f_name = 'gradient_weighting_'//trim(name)
      if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then
        idimf = 1
      else if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
        idimf = 6
      endif
      call field_create(f_name, itycat, ityloc, idimf, inoprv, iflid)
      call field_set_key_int(f_id, kwgrec, iflid)

    endif
  endif
enddo

! Postprocessing of slope tests

call field_get_key_id("slope_test_upwind_id", kislts)

itycat = FIELD_POSTPROCESS
ityloc = 1 ! cells

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_int(f_id, kislts, ifctsl)
    if (ifctsl.ge.0) then
      call field_get_key_struct_var_cal_opt(f_id, vcopt)

      ! Now create matching field
      if (vcopt%iconv.gt.0 .and. vcopt%blencv.gt.0 .and. vcopt%isstpc.eq.0) then
        ! Build name and label
        call field_get_name(f_id, f_name)
        name  = trim(f_name) // '_slope_upwind'
        call field_create(name, itycat, ityloc, idim1, inoprv, ifctsl)
        call field_set_key_int(ifctsl, keyvis, POST_ON_LOCATION)
        call field_set_key_int(f_id, kislts, ifctsl)
      endif
    endif
  endif
enddo

! Postprocessing of clippings

call field_get_key_id("clipping_id", kclipp)

itycat = FIELD_POSTPROCESS
ityloc = 1 ! cells

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_int(f_id, kclipp, clip_id)
    if (clip_id.ge.0) then

      ! Now create matching field
      ! Build name and label
      call field_get_name(f_id, f_name)

      call field_get_dim(f_id, f_dim)
      name  = trim(f_name) // '_clipped'
      call field_create(name, itycat, ityloc, f_dim, inoprv, clip_id)
      call field_set_key_int(clip_id, keyvis, POST_ON_LOCATION)
      call field_set_key_int(f_id, kclipp, clip_id)
    endif
  endif
enddo

! Fans output

itycat = FIELD_PROPERTY

n_fans = cs_fan_n_fans()
if (n_fans .gt. 0) then

  name = 'fan_id'
  ityloc = 1 ! cells

  call field_create(name, itycat, ityloc, idim1, inoprv, ifctsl)
  call field_set_key_int(ifctsl, keyvis, POST_ON_LOCATION)
  call field_set_key_int(ifctsl, keylog, 1)

endif

! Convection limiter

call field_get_key_id("convection_limiter_id", kcvlim)

itycat = FIELD_PROPERTY

do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then

    call field_get_key_struct_var_cal_opt(f_id, vcopt)

    ! Beta limiter or Roe-Sweby limiter
    if (vcopt%isstpc.eq.2 .or. vcopt%isstpc.eq.3) then
      ! Now create matching field
      ! Build name and label
      call field_get_name(f_id, f_name)

      call field_get_dim(f_id, f_dim)
      name = trim(f_name) // '_conv_lim'

      if (vcopt%isstpc.eq.2) then
        ityloc = 1 ! cells
      else
        ityloc = 2 ! Interior faces
      endif

      call field_create(name, itycat, ityloc, f_dim, inoprv, ifctsl)
      call field_set_key_int(ifctsl, keyvis, POST_ON_LOCATION)
      call field_set_key_int(ifctsl, keylog, 1)

      call field_set_key_int(f_id, kcvlim, ifctsl)
    endif
  endif
enddo

!===============================================================================

! VOF algorithm: the void fraction has its spectific convective flux
!-------------------------------------------------------------------

if (ivofmt.ge.0) then

  itycat = FIELD_EXTENSIVE + FIELD_PROPERTY

  ityloc = 2  ! inner faces
  f_name = 'inner_void_fraction_flux'
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call field_set_key_int(ivarfl(ivolf2), kimasf, f_id)

  ityloc = 3 ! boundary faces
  f_name = 'boundary_void_fraction_flux'
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call field_set_key_int(ivarfl(ivolf2), kbmasf, f_id)

endif

!===============================================================================

! Combustion
!-----------

if (     ippmod(icod3p).ne.-1         &
    .or. ippmod(iccoal).ne.-1         &
    .or. ippmod(icoebu).ne.-1         &
    .or. ippmod(icolwc).ne.-1 ) then

  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  ityloc = 3 ! boundary faces

  call field_create('boundary_ym_fuel',  &
                    itycat, ityloc, idim1, inoprv, ibym(1))
  call field_create('boundary_ym_oxydizer',  &
                    itycat, ityloc, idim1, inoprv, ibym(2))
  call field_create('boundary_ym_product',  &
                    itycat, ityloc, idim1, inoprv, ibym(3))
endif

!===============================================================================

! Turbulent anisotropic viscosity or user defined tensorial diffusivity for a
! scalar (exclusive or).
!----------------------------------------------------------------------

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! cells

if (idfm.eq.1 .or. itytur.eq.3 .and. idirsm.eq.1 &
    .or.darcy_anisotropic_dispersion.eq.1) then
  call field_create('anisotropic_turbulent_viscosity', itycat, ityloc, idim6, &
                    inoprv, ivsten)
  if (iturb.eq.32.and.iggafm.eq.1) then
    call field_create('anisotropic_turbulent_viscosity_scalar', itycat, &
                      ityloc, idim6, inoprv, ivstes)
  endif
endif

!===============================================================================
! Set some field keys
!===============================================================================

! Copy imrgra into the field structure
do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    vcopt%imrgra= imrgra
    call field_set_key_struct_var_cal_opt(f_id, vcopt)
  endif
enddo

! Copy field physical properties of species into the field structure

if (ippmod(igmix).ge.0) then

  call cs_parameters_define_field_key_gas_mix

  do f_id = 0, nfld - 1

    call field_get_name (f_id, name)

    if (name.eq.'y_o2') then

      sasp%mol_mas= 0.032d0
      sasp%cp = 930.d0
      sasp%vol_dif = 19.70d0
      sasp%mu_a = 5.086522d-8
      sasp%mu_b = 5.512391d-6
      sasp%lambda_a = 6.2d-5
      sasp%lambda_b = 8.1d-3
      sasp%muref = 1.919d-5
      sasp%lamref = 0.0244d0
      sasp%trefmu = 273.d0
      sasp%treflam = 273.d0
      sasp%smu = 139.d0
      sasp%slam = 240.d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_n2') then

      sasp%mol_mas = 0.028d0
      sasp%cp = 1042.d0
      sasp%vol_dif = 19.70d0
      sasp%mu_a = 4.210130d-8
      sasp%mu_b = 5.494348d-6
      sasp%lambda_a = 6.784141d-5
      sasp%lambda_b = 5.564317d-3
      sasp%muref = 1.663d-5
      sasp%lamref = 0.0242d0
      sasp%trefmu = 273.d0
      sasp%treflam = 273.d0
      sasp%smu = 107.d0
      sasp%slam = 150.d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_he') then

      sasp%mol_mas = 0.004d0
      sasp%cp = 5194.d0
      sasp%vol_dif = 2.67d0
      sasp%mu_a = 18.5752d-6
      sasp%mu_b = 0.0d0
      sasp%lambda_a = 0.144d0
      sasp%lambda_b = 0.0d0
      sasp%muref = 1.874d-5
      sasp%lamref = 0.647d0
      sasp%trefmu = 273.d0
      sasp%treflam = 273.d0
      sasp%smu = 78.d0
      sasp%slam = 78.d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_h2') then

      sasp%mol_mas = 0.002d0
      sasp%cp = 14560.d0
      sasp%vol_dif = 6.12d0
      sasp%mu_a = 1.93d-9
      sasp%mu_b = 8.40d-6
      sasp%lambda_a = 4.431d-4
      sasp%lambda_b = 5.334d-2
      sasp%muref = 8.411d-6
      sasp%lamref = 0.0168d0
      sasp%trefmu = 273.d0
      sasp%treflam = 273.d0
      sasp%smu = 97.d0
      sasp%slam = 120.d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_h2o_g') then

      sasp%mol_mas = 0.018d0
      sasp%cp = 2060.d0
      sasp%vol_dif = 13.10d0
      sasp%mu_a = 3.8496d-8
      sasp%mu_b = 8.2997d-6
      sasp%lambda_a = 7.6209d-5
      sasp%lambda_b = 0.016949d0
      sasp%muref = 1.12d-5
      sasp%lamref = 0.0181d0
      sasp%trefmu = 350.d0
      sasp%treflam = 300.d0
      sasp%smu = 1064.d0
      sasp%slam = 2200.d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    endif

  enddo

endif

! Check if scalars are buoyant and set n_buoyant_scal accordingly.
! It is then used in tridim to update buoyant scalars and density in U-P loop
call cs_parameters_set_n_buoyant_scalars

! If density is variable, a particular care must be taken when dealing with
! density in the unsteady term in the velocity pressure loop
if (irovar.eq.1) then
  ! EOS density, imposed after the correction step, so we need
  ! to keep the previous one, which is in balance with the mass
  f_name = 'density_mass'
  itycat = FIELD_PROPERTY
  ityloc = 1 ! cells
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  f_name = 'boundary_density_mass'
  itycat = FIELD_PROPERTY
  ityloc = 3 ! boundary faces
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
endif

! Some C mappings
call cs_field_pointer_map_boundary

return

end subroutine
