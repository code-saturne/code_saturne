!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
use radiat
use cplsat
use mesh
use post
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, ivar
integer          iflid, kcvlim, ifctsl, clip_id
integer          kdflim
integer          kturt, kfturt, kislts, keyvar, kclipp, kfturt_alpha
integer          turb_flux_model, turb_flux_model_type
integer          itycat, ityloc, idim1, idim3, idim6
logical          iprev, inoprv
integer          f_id, kscavr, f_type
integer          iopchr, ilog, ischcp
integer          iscdri, icla, iclap
integer          keyccl, keydri
integer          key_restart_id
integer          idfm, iggafm, nfld
integer          iflidp, idimf, idimc, n_fans
integer          f_dim, kiflux, kbflux

character(len=80) :: name, f_name, f_namec

type(var_cal_opt) :: vcopt

procedure() :: hide_property

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

! Key id for restart file
call field_get_key_id("restart_file", key_restart_id)
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

    call field_get_key_int(ivarfl(ivar), kturt, turb_flux_model)
    turb_flux_model_type = turb_flux_model / 10

    if (turb_flux_model_type.gt.0) then
      if (turb_flux_model_type.eq.3) then
        idfm = 1
      endif

      ! GGDH or AFM on current scalar
      ! and if DFM, GGDH on the scalar variance
      iggafm = 1

    ! If the user has chosen a tensorial diffusivity
    else if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
      idfm = 1
    endif

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
  call field_set_key_int(f_id, keyvis, 0)
  call field_set_key_int(f_id, key_restart_id, RESTART_AUXILIARY)
else if (iphydr.eq.2) then
  call field_find_or_create('hydrostatic_pressure_prd', &
                            itycat, ityloc, idim1, f_id)
  call field_set_key_int(f_id, key_restart_id, RESTART_AUXILIARY)
endif

! Hybrid blending field

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
ischcp = vcopt%ischcv
if (ischcp.eq.3) then
  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  ityloc = 1 ! cells
  call field_find_or_create('hybrid_blend', itycat, ityloc, idim1, f_id)
end if

! friction velocity at the wall, in the case of a LES calculation
! with van Driest-wall damping (delayed here rather than placed in
! addfld, as idries may be set in modini).

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces

if (itytur.eq.4 .and. idries.eq.1) then
  call field_find_or_create("boundary_ustar", itycat, ityloc, idim1, f_id)
endif

if (staggered.eq.1) then

  ! Head loss on interior faces

  itycat = FIELD_PROPERTY
  ityloc = 2 ! inner faces
  f_name = 'inner_face_head_loss'

  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)

  ! Head loss on boundary faces

  itycat = FIELD_PROPERTY
  ityloc = 3 ! boundary faces
  f_name = 'boundary_face_head_loss'

  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)

  ! Source term on interior faces

  itycat = FIELD_PROPERTY
  ityloc = 2 ! inner faces
  f_name = 'inner_face_source_term'

  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)

endif

! Interior mass flux field
!-------------------------

itycat = FIELD_EXTENSIVE
ityloc = 2 ! inner faces

! Mass flux for the class on interior faces
f_name = 'inner_mass_flux'

if (istmpf.ne.1 .or. staggered.eq.1) then
  call field_create(f_name, itycat, ityloc, idim1, iprev, f_id)
else
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
endif
call hide_property(f_id)

! The same mass flux for every variable, an other mass flux
! might be defined hereafterwards
do ivar = 1, nvar
  call field_set_key_int(ivarfl(ivar), kimasf, f_id)
enddo

! Rusanov flux
if (irijnu.eq.2) then
  ityloc = 2 ! inner faces
  call field_create('i_rusanov_diff', itycat, ityloc, idim1, inoprv, f_id)

  ityloc = 3 ! boundary faces
  call field_create('b_rusanov_diff', itycat, ityloc, idim1, inoprv, f_id)
endif

! Boundary Mass flux field
!-------------------------

itycat = FIELD_EXTENSIVE
ityloc = 3 ! boundary faces

! Mass flux for the class on interior faces
f_name = 'boundary_mass_flux'

if (istmpf.ne.1 .or. staggered.eq.1) then
  call field_create(f_name, itycat, ityloc, idim1, iprev, f_id)
else
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
endif
call hide_property(f_id)

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
    call field_set_key_int(f_id, keylog, 0)

    ! Set the inner mass flux index
    call field_set_key_int(iflid, kimasf, f_id)

    ! If the scalar is the representant of a class, then
    ! set the mass flux index to all members of the class
    if (icla.ne.0) then
      do ii = 0, nfld-1
        call field_get_key_int(ii, keyccl, iclap)
        call field_get_type(ii, f_type)
        if (icla.eq.iclap .and. &
          iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
          call field_set_key_int(ii, kimasf, f_id)
        endif
      enddo
    endif

    ityloc = 3 ! variables defined on boundary faces

    ! Mass flux for the class on boundary faces
    f_name = 'boundary_mass_flux_'//trim(name)
    call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    call field_set_key_int(f_id, keylog, 0)

    ! Set the boundary mass flux index
    call field_set_key_int(iflid, kbmasf, f_id)

    ! If the scalar is the representant of a class, then
    ! set the mass flux index to all members of the class
    if (icla.ne.0) then
      do ii = 0, nfld-1
        call field_get_key_int(ii, keyccl, iclap)
        call field_get_type(ii, f_type)
        if (icla.eq.iclap .and. &
          iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
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
      call field_set_key_int(f_id, keylog, ilog)

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
    ! Is this field not managed by CDO? Not useful for CDO
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

      call field_get_key_struct_var_cal_opt(f_id, vcopt)
      if (vcopt%iwgrec.eq.1 .and. vcopt%idiff .ge. 1) then

        call field_get_name(f_id, name)
        f_name = 'gradient_weighting_'//trim(name)
        if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then
          idimf = 1
        else if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
          idimf = 6
        endif

        ! Check if weighting has been assigned already
        ! (in case of head losses or tensorial porosity).
        ! If present; check compatibility.
        ! If not present yet, create it.
        call field_get_key_int(f_id, kwgrec, iflid)
        if (iflid .ge. 0) then
          call field_get_dim(iflid, idimc)
          if (idimc .ne. idimf) then
            call field_get_name(iflid, f_namec)
            write(nfecra,8100) name, idimf, f_namec, idimc
            call csexit(1)
          endif
        else
          call field_create(f_name, 0, ityloc, idimf, inoprv, iflid)
          call field_set_key_int(f_id, kwgrec, iflid)
        endif

      endif

    endif
  endif
enddo

! Postprocessing of slope tests

call field_get_key_id("slope_test_upwind_id", kislts)

itycat = FIELD_POSTPROCESS
ityloc = 1 ! cells

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE ?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    ! Is this field not managed by CDO ?
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

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

    endif ! CDO ?
  endif ! VARIABLE ?
enddo

! Postprocessing of clippings

call field_get_key_id("clipping_id", kclipp)

itycat = FIELD_POSTPROCESS
ityloc = 1 ! cells

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE ?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    ! Is this field not managed by CDO ?
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

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

    endif ! CDO ?
  endif ! VARIABLE ?
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

    call field_get_key_int(f_id, kcvlim, ifctsl)

    ! Beta limiter or Roe-Sweby limiter
    if (vcopt%isstpc.eq.2.or.ifctsl.ne.-1) then
      ! Now create matching field
      ! Build name and label
      call field_get_name(f_id, f_name)

      call field_get_dim(f_id, f_dim)
      name = trim(f_name) // '_conv_lim'

      ityloc = 1 ! cells

      call field_create(name, itycat, ityloc, f_dim, inoprv, ifctsl)
      call field_set_key_int(ifctsl, keylog, 1)

      call field_set_key_int(f_id, kcvlim, ifctsl)
    endif
  endif
enddo

! Diffusion limiter

call field_get_key_id("diffusion_limiter_id", kdflim)

itycat = FIELD_PROPERTY

do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    ! Is this field not managed by CDO ? Not useful with CDO schemes
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

      call field_get_key_int(f_id, kdflim, ifctsl)

      if (ifctsl.ne.-1) then
        ! Now create matching field
        ! Build name and label
        call field_get_name(f_id, f_name)

        call field_get_dim(f_id, f_dim)
        name = trim(f_name) // '_diff_lim'

        ityloc = 1 ! cells

        call field_create(name, itycat, ityloc, f_dim, inoprv, ifctsl)
        call field_set_key_int(ifctsl, keyvis, POST_ON_LOCATION)
        call field_set_key_int(ifctsl, keylog, 1)

        call field_set_key_int(f_id, kdflim, ifctsl)
      endif

    endif ! CDO ?
  endif ! VARIABLE ?
enddo

!===============================================================================

! VOF algorithm: the void fraction has its specific convective flux
!                and void fraction flux needs to be stored
!-------------------------------------------------------------------

if (ivofmt.gt.0) then

  itycat = FIELD_EXTENSIVE + FIELD_PROPERTY

  call field_get_key_id("inner_flux_id", kiflux)
  call field_get_key_id("boundary_flux_id", kbflux)

  ityloc = 2  ! inner faces
  f_name = 'inner_volume_flux'
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call field_set_key_int(ivarfl(ivolf2), kimasf, f_id)

  f_name = 'inner_void_fraction_flux'
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call field_set_key_int(ivarfl(ivolf2), kiflux, f_id)

  ityloc = 3 ! boundary faces
  f_name = 'boundary_volume_flux'
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call field_set_key_int(ivarfl(ivolf2), kbmasf, f_id)

  f_name = 'boundary_void_fraction_flux'
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call field_set_key_int(ivarfl(ivolf2), kbflux, f_id)

  if (idrift.gt.0) then
    ityloc = 2  ! inner faces
    f_name = 'inner_drift_velocity_flux'
    call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)

    ityloc = 3 ! boundary faces
    f_name = 'boundary_drift_velocity_flux'
    call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  endif

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

if (idfm.eq.1.or.iggafm.eq.1.or. itytur.eq.3 .and. idirsm.eq.1) then
  call field_create('anisotropic_turbulent_viscosity', itycat, ityloc, idim6, &
                    inoprv, ivsten)
  ! By default set log printing to false
  call field_set_key_int(ivsten, keylog, 0)
  if (iturb.eq.32.and.iggafm.eq.1) then
    call field_create('anisotropic_turbulent_viscosity_scalar', itycat, &
                      ityloc, idim6, inoprv, ivstes)
    ! By default set log printing to false
    call field_set_key_int(ivstes, keylog, 0)
  endif
endif

!===============================================================================
! Change some field settings
!===============================================================================

! If ALE, for fluid dtructure interaction, mass fluxes may be needed at the
! previous iteration
if (iale.ge.1) then
  call field_get_key_int(ivarfl(ipr), kimasf, f_id)
  call field_set_n_previous(f_id, 1)
  call field_get_key_int(ivarfl(ipr), kbmasf, f_id)
  call field_set_n_previous(f_id, 1)
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
    if (vcopt%imrgra .lt. 0) vcopt%imrgra= imrgra
    call field_set_key_struct_var_cal_opt(f_id, vcopt)
  endif
enddo

! Check if scalars are buoyant and set n_buoyant_scal accordingly.
! It is then used in tridim to update buoyant scalars and density in U-P loop
call cs_velocity_pressure_set_n_buoyant_scalars

! For Low Mach and compressible (increment) algorithms, a particular care
! must be taken when dealing with density in the unsteady term in the velocity
! pressure loop
if (irovar.eq.1.and.(idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3)) then
  ! EOS density, imposed after the correction step, so we need
  ! to keep the previous one, which is in balance with the mass
  f_name = 'density_mass'
  itycat = FIELD_PROPERTY
  ityloc = 1 ! cells
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call hide_property(f_id)
  f_name = 'boundary_density_mass'
  itycat = FIELD_PROPERTY
  ityloc = 3 ! boundary faces
  call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
  call hide_property(f_id)
endif

! Some C mappings
call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

!===============================================================================
! Formats
!===============================================================================

 8100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR: STOP IN THE INITIAL DATA SETUP',                   /,&
'@    =======',                                                 /,&
'@',                                                            /,&
'@    Variable ', a, ' should be assigned a',                   /,&
'@    gradient_weighting_field, of dimension ', i2,             /,&
'@    but has already been assigned field ', a,                 /,&
'@    ', a, ', of dimension ', i2,                              /,&
'@',                                                            /,&
'@  The calculation cannot be executed',                        /,&
'@',                                                            /,&
'@  Check parameters.',                                         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

return

end subroutine
