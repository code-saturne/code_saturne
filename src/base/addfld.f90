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

!===============================================================================
! Purpose:
! --------

!> \file addfld.f90
!>
!> \brief Add additional fields based on user options.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________


subroutine addfld

!===============================================================================
! Module files
!===============================================================================

use atincl, only: compute_z_ground, imeteo
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
use cplsat
use mesh
use post
use field
use turbomachinery
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, iscal
integer          iscacp, ifcvsl, kbfid
integer          iflid, iopchr
integer          itycat, ityloc, idim1, idim3
integer          f_id, potr, poti
integer          f_vis, f_log, ivtmp
integer          kturt, kfturt, turb_flux_model, turb_flux_model_type
integer          kfturt_alpha
integer          keycpl, keydri
integer          ivar, iscdri
logical          iprev, inoprv, is_set
integer          key_t_ext_id
integer          nfld
integer          n_prev
integer          t_ext
integer          kclipp
integer          key_turb_schmidt, kscavr, key_turb_diff, key_sgs_sca_coef
integer          key_restart_id
integer          var_f_id

character(len=80) :: name, f_name, f_label, s_label, s_name
type(var_cal_opt) :: vcopt_dfm, vcopt_alpha, vcopt

procedure() :: add_variable_field, add_property_field, hide_property

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! The itycat variable is used to define field categories. It is used in Fortran
! code with hard-coded values, but in the C API, those values are based on
! (much clearer) category mask definitions in cs_field.h.

itycat = FIELD_INTENSIVE + FIELD_VARIABLE  ! for most variables
ityloc = 1 ! variables defined on cells
idim1  = 1
idim3  = 3
iprev  = .true.    ! variables have previous value
inoprv = .false.   ! variables have no previous value
iopchr = 1         ! Postprocessing level for variables

! Keys not stored globally
call field_get_key_id('turbulent_flux_model', kturt)
call field_get_key_id('turbulent_flux_id', kfturt)
call field_get_key_id('alpha_turbulent_flux_id', kfturt_alpha)
call field_get_key_id('coupled', keycpl)
call field_get_key_id("first_moment_id", kscavr)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

! Restart file key
call field_get_key_id("restart_file", key_restart_id)

! Number of fields
call field_get_n_fields(nfld)

!===============================================================================
! 0. Initialization
!===============================================================================

call field_get_key_id("boundary_value_id", kbfid)

call field_get_key_id('log', keylog)
call field_get_key_id('label', keylbl)

!===============================================================================
! 1. Additional variable fields
!===============================================================================

! User variables
!---------------

do ii = 1, nscal

  if (isca(ii) .gt. 0) then

    ivar = isca(ii)
    f_id = ivarfl(ivar)

    call field_get_key_int(ivarfl(ivar), keyvis, f_vis)
    call field_get_key_int(ivarfl(ivar), keylog, f_log)
    call field_get_key_int(ivarfl(ivar), kturt, turb_flux_model)
    turb_flux_model_type = turb_flux_model / 10

    if (turb_flux_model_type.gt.0) then
      call field_get_name (f_id, name)
      f_name = trim(name)//'_turbulent_flux'

      if (turb_flux_model_type.eq.3) then
        call add_variable_field(f_name, f_name, 3, ivtmp)
        iflid = ivarfl(ivtmp)

        call field_set_key_int(iflid, keycpl, 1)
        ! Tensorial diffusivity
        call field_get_key_struct_var_cal_opt(iflid, vcopt_dfm)
        vcopt_dfm%idften = ANISOTROPIC_RIGHT_DIFFUSION
        call field_set_key_struct_var_cal_opt(iflid, vcopt_dfm)

        call field_get_key_id("clipping_id",kclipp)
        call field_set_key_int(iflid, kclipp, 1)

      else
        itycat = FIELD_INTENSIVE + FIELD_PROPERTY  ! for properties

        call field_create(f_name, itycat, ityloc, idim3, iprev, iflid)

        call field_set_key_int(iflid, keyvis, f_vis)
        call field_set_key_int(iflid, keylog, f_log)
      endif

      call field_set_key_int(ivarfl(ivar), kfturt, iflid)

      ! Elliptic Blending (AFM or DFM)
      if (     turb_flux_model.eq.11 .or. turb_flux_model.eq.21  &
          .or. turb_flux_model.eq.31) then
        f_name = trim(name)//'_alpha'

        call add_variable_field(f_name, f_name, 1, ivtmp)
        iflid = ivarfl(ivtmp)

        ! Elliptic equation (no convection, no time term)
        call field_get_key_struct_var_cal_opt(iflid, vcopt_alpha)
        vcopt_alpha%iconv = 0
        vcopt_alpha%istat = 0
        call field_set_key_struct_var_cal_opt(iflid, vcopt_alpha)

        call field_set_key_int(ivarfl(ivar), kfturt_alpha, iflid)
      endif

    endif

  endif

enddo

! Hydrostatic pressure used to update pressure BCs
if (icalhy.eq.1) then
  f_name  = 'hydrostatic_pressure'
  f_label = 'Hydrostatic Pressure'
  call add_variable_field(f_name, f_label, 1, ivar)
  f_id = ivarfl(ivar)

  call field_set_key_int(f_id, keyvis, 0)

  ! Elliptic equation (no convection, no time term)
  call field_get_key_struct_var_cal_opt(f_id, vcopt)
  vcopt%iconv = 0
  vcopt%istat = 0
  vcopt%nswrsm = 2
  vcopt%idifft = 0
  vcopt%relaxv = 1.d0 ! No relaxation, even for steady algorithm.
  call field_set_key_struct_var_cal_opt(f_id, vcopt)
endif

! Head losses weighting field in case of Lagrangian deposition and
! reentrainment model (general case in varpos, but the Lagrangian
! options are not know yet at the call site, so we have a similar
! code block here for this special case.

if (iflow .gt. 0 .and. idtten .lt. 0) then
  call field_create('dttens', itycat, ityloc, 6, .false., idtten)
  call field_set_key_int(idtten, keyvis, POST_ON_LOCATION)
  call field_set_key_int(idtten, keylog, 1)
  call field_set_key_int(ivarfl(ipr), kwgrec, idtten)
endif

!===============================================================================
! 2. Additional property fields
!===============================================================================

! Add a scalar diffusivity when defined as variable.
! The kivisl key should be equal to -1 for constant diffusivity,
! and f_id for a variable diffusivity defined by field f_id
! Assuming the first field created is not a diffusivity property
! (we define variables first), f_id > 0, so we use 0 to indicate
! the diffusivity is variable but its field has not been created yet.

do ii = 1, nscal
  f_id = ivarfl(isca(ii))
  call field_get_key_int(f_id, kivisl, ifcvsl)
  call field_get_key_int(f_id, kscavr, var_f_id)
  if (ifcvsl.eq.0 .and. var_f_id.lt.0) then
    ! Build name and label, using a general rule, with a
    ! fixed name for temperature or enthalpy
    call field_get_name(f_id, s_name)
    call field_get_label(f_id, s_label)
    if (ii.eq.iscalt) then
      s_name = 'thermal'
      s_label = 'Th'
    endif
    call field_get_key_int(f_id, kscacp, iscacp)
    if (iscacp.gt.0) then
      f_name  = trim(s_name) // '_conductivity'
      f_label = trim(s_label) // ' Cond'
    else
      f_name  = trim(s_name) // '_diffusivity'
      f_label = trim(s_label) // ' Diff'
    endif
    ! Special case for electric arcs: real and imaginary electric
    ! conductivity is the same (and ipotr < ipoti)
    if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
      call field_get_id('elec_pot_r', potr)
      call field_get_id('elec_pot_i', poti)
      if (f_id.eq.potr) then
        f_name = 'elec_sigma'
        f_label = 'Sigma'
      else if (f_id.eq.poti) then
        call field_get_key_int(potr, kivisl, ifcvsl)
        call field_set_key_int(poti, kivisl, ifcvsl)
        cycle ! go to next scalar in loop, avoid creating property
      endif
    endif
    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., ifcvsl)
    call field_set_key_int(ivarfl(isca(ii)), kivisl, ifcvsl)
  endif
enddo

! For variances, the diffusivity is that of the associated scalar,
! and must not be initialized first.

do ii = 1, nscal
  if (iscavr(ii).gt.0) then
    f_id = ivarfl(isca(ii))
    call field_get_key_int(ivarfl(isca(iscavr(ii))), kivisl, ifcvsl)
    call field_is_key_set(f_id, kivisl, is_set)
    if (is_set.eqv..true.) then
      write(nfecra,7040) f_id, ivarfl(isca(iscavr(ii))), ifcvsl
    else
      call field_set_key_int(f_id, kivisl, ifcvsl)
    endif
  endif
enddo

! Add a scalar turbulent diffusivity field
call field_get_key_id("turbulent_diffusivity_id", key_turb_diff)

do ii = 1, nvar
  f_id = ivarfl(ii)
  call field_get_key_int(f_id, key_turb_diff, ifcvsl)
  call field_get_key_int(f_id, kscavr, var_f_id)
  if (ifcvsl.ge.0 .and. var_f_id.lt.0) then
    ! Build name and label, using a general rule
    call field_get_name(f_id, s_name)
    call field_get_label(f_id, s_label)
    f_name  = trim(s_name) // '_turb_diffusivity'
    f_label = trim(s_label) // ' Turb Diff'

    if (ippmod(islfm).ge.0) then
      call field_get_key_int(ivarfl(isca(ifm)), key_turb_diff, ifcvsl)
      if (ii .ne. isca(ifm)) then
        call field_set_key_int(f_id,  key_turb_diff, ifcvsl)
        cycle
      endif
    endif
    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., ifcvsl)
    call field_set_key_int(ivarfl(ii), key_turb_diff, ifcvsl)
  endif
enddo

! For variances, the turbulent diffusivity is that of the associated scalar,
! and must not be initialized first.
do ii = 1, nscal
  if (iscavr(ii).gt.0) then
    f_id = ivarfl(isca(ii))
    call field_get_key_int(ivarfl(isca(iscavr(ii))), key_turb_diff, ifcvsl)
    call field_is_key_set(f_id, key_turb_diff, is_set)
    if (is_set.eqv..true.) then
      write(nfecra,7041) f_id, ivarfl(isca(iscavr(ii))), ifcvsl
    else
      call field_set_key_int(f_id, key_turb_diff, ifcvsl)
    endif
  endif
enddo

if (iturb.eq.41) then
  ! Add a subgrid-scale scalar flux coefficient field
  call field_get_key_id("sgs_scalar_flux_coef_id", key_sgs_sca_coef)

  do ii = 1, nvar
    f_id = ivarfl(ii)
    call field_get_key_int(f_id, key_sgs_sca_coef, ifcvsl)
    call field_get_key_int(f_id, kscavr, var_f_id)
    if (ifcvsl.ge.0 .and. var_f_id.lt.0) then
      ! Build name and label, using a general rule
      call field_get_name(f_id, s_name)
      call field_get_label(f_id, s_label)
      f_name  = trim(s_name) // '_sgs_flux_coef'
      f_label = trim(s_label) // ' SGS Flux Coef'

      if (ippmod(islfm).ge.0) then
        call field_get_key_int(ivarfl(isca(ifm)), key_sgs_sca_coef, ifcvsl)
        if (ii .ne. isca(ifm)) then
          call field_set_key_int(f_id,  key_sgs_sca_coef, ifcvsl)
          cycle
        endif
      endif
      ! Now create matching property
      call add_property_field(f_name, f_label, 1, .false., ifcvsl)
      call field_set_key_int(ivarfl(ii), key_sgs_sca_coef, ifcvsl)
    endif
  enddo

  ! For variances, the subgrid-scale flux is that of the associated scalar,
  ! and must not be initialized first.
  do ii = 1, nscal
    if (iscavr(ii).gt.0) then
      f_id = ivarfl(isca(ii))
      call field_get_key_int(ivarfl(isca(iscavr(ii))), key_sgs_sca_coef, ifcvsl)
      call field_is_key_set(f_id, key_sgs_sca_coef, is_set)
      if (is_set.eqv..true.) then
        write(nfecra,7042) f_id, ivarfl(isca(iscavr(ii))), ifcvsl
      else
        call field_set_key_int(f_id, key_sgs_sca_coef, ifcvsl)
      endif
    endif
  enddo
endif

! Add a scalar density when defined as variable and different from the bulk.
! WARNING: it must be consitent with continuity equation, this is used
! for fluid solid computation with passive scalars with different density in the solid.
! The kromsl key should be equal to -1 for constant density
! and f_id for a variable density defined by field f_id
! Assuming the first field created is not a density property
! (we define variables first), f_id > 0, so we use 0 to indicate
! the density is variable but its field has not been created yet.

do ii = 1, nscal
  f_id = ivarfl(isca(ii))
  call field_get_key_int(f_id, kromsl, ifcvsl)
  call field_get_key_int(f_id, kscavr, var_f_id)
  if (ifcvsl.eq.0 .and.var_f_id.lt.0) then
    ! Build name and label, using a general rule, with a
    ! fixed name for temperature or enthalpy
    call field_get_name(f_id, s_name)
    call field_get_label(f_id, s_label)
    f_name  = trim(s_name) // '_density'
    f_label = trim(s_label) // ' Rho'

    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., ifcvsl)
    call field_set_key_int(ivarfl(isca(ii)), kromsl, ifcvsl)
  endif
enddo

! For variances, the density is that of the associated scalar,
! and must not be initialized first.

do ii = 1, nscal
  if (iscavr(ii).gt.0) then
    f_id = ivarfl(isca(ii))
    call field_get_key_int(ivarfl(isca(iscavr(ii))), kromsl, ifcvsl)
    call field_is_key_set(f_id, kromsl, is_set)
    if (is_set.eqv..true.) then
      write(nfecra,7040) f_id, ivarfl(isca(iscavr(ii))), ifcvsl
    else
      call field_set_key_int(f_id, kromsl, ifcvsl)
    endif
  endif
enddo

! Add a scalar turbulent Schmidt field
call field_get_key_id("turbulent_schmidt_id", key_turb_schmidt)

do ii = 1, nvar
  f_id = ivarfl(ii)
  call field_get_key_int(f_id, key_turb_schmidt, ifcvsl)
  call field_get_key_int(f_id, kscavr, var_f_id)
  if (ifcvsl.ge.0 .and. var_f_id.lt.0) then
    ! Build name and label, using a general rule, with a
    ! fixed name for temperature or enthalpy
    call field_get_name(f_id, s_name)
    call field_get_label(f_id, s_label)
    f_name  = trim(s_name) // '_turb_schmidt'
    f_label = trim(s_label) // ' ScT'

    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., ifcvsl)
    call field_set_key_int(ivarfl(ii), key_turb_schmidt, ifcvsl)
  endif
enddo

! For variances, the Schmidt is that of the associated scalar,
! and must not be initialized first.
do ii = 1, nscal
  if (iscavr(ii).gt.0) then
    f_id = ivarfl(isca(ii))
    call field_get_key_int(ivarfl(isca(iscavr(ii))), key_turb_schmidt, ifcvsl)
    call field_is_key_set(f_id, key_turb_schmidt, is_set)
    if (is_set.eqv..true.) then
      write(nfecra,7040) f_id, ivarfl(isca(iscavr(ii))), ifcvsl
    else
      call field_set_key_int(f_id, key_turb_schmidt, ifcvsl)
    endif
  endif
enddo

! Boundary roughness (may be already created by the atmospheric module)
if (iwallf.eq.5.or.iwallf.eq.6) then
  idim1  = 1
  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  ityloc = 3 ! boundary faces

  call field_find_or_create('boundary_roughness', itycat, ityloc, idim1, iflid)
endif

! Van Driest damping
if (idries.eq.-1) then
  if (iturb.eq.40) then
    idries = 1
  elseif (iturb.eq.41.or.iturb.eq.42) then
    idries = 0
  endif
endif

! Wall distance for some turbulence models
! and for Lagrangian multilayer deposition
! for DRSM models, needed for inlets

if ( iturb.eq.23.or.                                &
    (itytur.eq.3).or.                               &
     (iturb.eq.30.and.irijec.eq.1).or.              &
     (itytur.eq.4.and.idries.eq.1).or.              &
     iturb.eq.60.or.iturb.eq.70.or.                 &
     iflow.eq.1.or.                                 &
     hybrid_turb.eq.4) then
  ineedy = 1
endif

! User may set ineedy to 1
if (ineedy.eq.1) then
  ! Working variable
  f_name  = 'wall_distance'
  f_label = 'Wall distance'
  call add_variable_field(f_name, f_label, 1, ivar)
  iflid = ivarfl(ivar)
  call field_set_key_int(iflid, key_restart_id, RESTART_AUXILIARY)

  ! Elliptic equation (no convection, no time term)
  call field_get_key_struct_var_cal_opt(iflid, vcopt)
  vcopt%iconv = 0
  vcopt%istat = 0
  vcopt%nswrsm = 2
  vcopt%idifft = 0
  vcopt%relaxv = 1.d0 ! No relaxation, even for steady algorithm.
  call field_set_key_struct_var_cal_opt(iflid, vcopt)

  ! Working field to store value of the solved variable at the previous
  ! time step if needed (ALE)
  if (iale.ne.0.or.iturbo.ne.0) then
    f_name  = 'wall_distance_aux_pre'
    call add_property_field(f_name, f_label, 1, .false., iflid)
    call hide_property(iflid)
  endif

  ! Dimensionless wall distance "y+"
  !> non-dimensional distance \f$y^+\f$ between a given volume and the closest
  !> wall, when it is necessary (LES with van Driest-wall damping).
  if (itytur.eq.4.and.idries.eq.1) then
    f_name  = 'wall_yplus'
    f_label = 'Wall Y+'
    call add_variable_field(f_name, f_label, 1, ivar)
    iflid = ivarfl(ivar)

    call field_set_key_int(iflid, keyvis, 1)
    call field_set_key_int(iflid, keylog, 1)

    ! Pure convection (no time term)
    call field_get_key_struct_var_cal_opt(iflid, vcopt)
    vcopt%iconv = 1 ! default
    vcopt%istat = 0
    vcopt%idiff = 0
    vcopt%idifft = 0
    vcopt%relaxv = 1.d0 ! No relaxation, even for steady algorithm.
    vcopt%blencv = 0.d0 ! Pure upwind
    vcopt%epsilo = 1.d-5 ! by default, not an high precision
    call field_set_key_struct_var_cal_opt(iflid, vcopt)

    ! Activate the drift for all scalars with key "drift" > 0
    iscdri = 1

    ! GNU function to return the value of iscdri
    ! with the bit value of iscdri at position
    ! 'DRIFT_SCALAR_ADD_DRIFT_FLUX' set to one
    iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

    iscdri = ibset(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX)

    call field_set_key_int(iflid, keydri, iscdri)

  endif

endif

if (ippmod(iatmos).ge.0.and.compute_z_ground) then
  f_name  = 'z_ground'
  f_label = 'Z ground'
  call add_variable_field(f_name, f_label, 1, ivar)
  iflid = ivarfl(ivar)

  ! Pure convection equation (no convection, no time term)
  call field_get_key_struct_var_cal_opt(iflid, vcopt)
  vcopt%iconv = 1
  vcopt%blencv= 0.d0 ! Pure upwind
  vcopt%istat = 0
  vcopt%nswrsm = 100
  vcopt%epsrsm = 1.d-3
  vcopt%idiff  = 0
  vcopt%idifft = 0
  vcopt%relaxv = 1.d0 ! No relaxation, even for steady algorithm.
  call field_set_key_struct_var_cal_opt(iflid, vcopt)

  ! Activate the drift for all scalars with key "drift" > 0
  iscdri = 1

  ! GNU function to return the value of iscdri
  ! with the bit value of iscdri at position
  ! 'DRIFT_SCALAR_ADD_DRIFT_FLUX' set to one
  iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

  iscdri = ibset(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX)

  call field_set_key_int(iflid, keydri, iscdri)
endif

if (imeteo.ge.2) then
  f_name  = 'meteo_pressure'
  f_label = 'Meteo pressure'
  call add_variable_field(f_name, f_label, 1, ivar)
  iflid = ivarfl(ivar)

  ! Pure convection equation (no convection, no time term)
  call field_get_key_struct_var_cal_opt(iflid, vcopt)
  vcopt%iconv = 0
  vcopt%blencv= 0.d0 ! Pure upwind
  vcopt%istat = 0
  vcopt%idircl = 1
  vcopt%nswrsm = 100
  vcopt%nswrgr = 100
  vcopt%imrgra = 0
  vcopt%imligr = -1
  vcopt%epsilo = 0.000001
  vcopt%epsrsm = 1.d-3
  vcopt%epsrgr = 0.0001
  vcopt%climgr = 1.5
  vcopt%idiff  = 1
  vcopt%idifft = 0
  vcopt%idften = 1
  vcopt%relaxv = 1.d0 ! No relaxation, even for steady algorithm.
  vcopt%thetav = 1
  call field_set_key_struct_var_cal_opt(iflid, vcopt)

  f_name  = 'meteo_density'
  f_label = 'Meteo density'
  ! Now create matching property
  call add_property_field(f_name, f_label, 1, .false., iflid)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'meteo_temperature'
  f_label = 'Meteo Temperature'
  ! Now create matching property
  call add_property_field(f_name, f_label, 1, .false., iflid)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'meteo_pot_temperature'
  f_label = 'Meteo pot Temperature'
  ! Now create matching property
  call add_property_field(f_name, f_label, 1, .false., iflid)
  call field_set_key_int(iflid, keylog, 1)

  if (ippmod(iatmos).eq.2) then
    f_name  = 'meteo_humidity'
    f_label = 'Meteo humidity'
    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., iflid)

    f_name  = 'meteo_drop_nb'
    f_label = 'Meteo drop. nb'
    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., iflid)
  endif
  f_name  = 'meteo_velocity'
  f_label = 'Meteo velocity'
  ! Now create matching property
  call add_property_field(f_name, f_label, 3, .false., iflid)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'meteo_tke'
  f_label = 'Meteo TKE'
  ! Now create matching property
  call add_property_field(f_name, f_label, 1, .false., iflid)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'meteo_eps'
  f_label = 'Meteo epsilon'
  ! Now create matching property
  call add_property_field(f_name, f_label, 1, .false., iflid)
  call field_set_key_int(iflid, keylog, 1)

  ! DRSM models, store Rxz/k
  if (itytur.eq.3) then
    f_name  = 'meteo_shear_anisotropy'
    f_label = 'meteo_shear_anisotropy'
    ! Now create matching property
    call add_property_field(f_name, f_label, 1, .false., iflid)
    call field_set_key_int(iflid, keylog, 1)
  endif

endif

if (compute_porosity_from_scan) then
  f_name  = 'nb_scan_points'
  f_label = 'Scan points number'
  call add_property_field(f_name, f_label, 1, .false., iflid)

  call field_set_key_int(iflid, keyvis, POST_ON_LOCATION)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'solid_roughness'
  f_label = 'Solid roughness'
  call add_property_field(f_name, f_label, 1, .false., iflid)

  call field_set_key_int(iflid, keyvis, POST_ON_LOCATION)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'cell_scan_points_cog'
  f_label = 'Point centers'
  call add_property_field(f_name, f_label, 3, .false., iflid)

  call field_set_key_int(iflid, keyvis, POST_ON_LOCATION)
  call field_set_key_int(iflid, keylog, 1)

  f_name  = 'cell_scan_points_color'
  f_label = 'Cell color'
  call add_property_field(f_name, f_label, 3, .false., iflid)

  call field_set_key_int(iflid, keyvis, POST_ON_LOCATION)
  call field_set_key_int(iflid, keylog, 1)
endif

!===============================================================================
! 3. Additional postprocessing fields
!===============================================================================

! Fields used to save postprocessing data

ityloc = 3 ! boundary faces

itycat = FIELD_INTENSIVE + FIELD_PROPERTY

! In case of ALE or postprocessing, ensure boundary forces are tracked

call field_get_id_try('boundary_forces', iforbr)  ! may already be present

if (iale.ge.1) then
  itycat = FIELD_EXTENSIVE + FIELD_POSTPROCESS
  call field_find_or_create('boundary_forces', itycat, ityloc, idim3, &
                            iforbr)
endif

! Boundary efforts postprocessing for immersed boundaries, create field

if (iporos.eq.3) then
  itycat = FIELD_EXTENSIVE + FIELD_POSTPROCESS
  call field_create('immersed_pressure_force', itycat, ityloc, idim3, inoprv, &
                    iflid )

  if (iturb.ne.0) then
    call field_create('immersed_boundary_uk', itycat, ityloc, idim1, inoprv, &
                      iflid)
    call field_create('immersed_boundary_yplus', itycat, ityloc, idim1, inoprv, &
                      iflid)
    call field_create('immersed_boundary_dplus', itycat, ityloc, idim1, inoprv, &
                      iflid)
  endif
endif

itycat = FIELD_INTENSIVE + FIELD_PROPERTY

! In case of condensation or postprocessing, ensure y+ is tracked

call field_get_id_try('yplus', iyplbr)  ! may already be present

if (icondb.ge.0.or.icondv.ge.0) then
  f_id = iyplbr ! Test if pre-existing
  call field_find_or_create('yplus', itycat, ityloc, idim1, iyplbr)
  if (f_id .lt. 0) then                ! Set some properties if new
    call field_set_key_str(iyplbr, keylbl, 'Yplus')
    call field_set_key_int(iyplbr, keylog, 1)
  endif
endif

! Some mappings

call cs_field_pointer_map_boundary

! Cooling towers mappings
if (ippmod(iaeros).ge.0) then
  call cs_ctwr_field_pointer_map
endif

call field_get_id_try('boundary_temperature', itempb)

if (itempb .ge. 0) then
  call field_get_id_try('temperature', f_id)
  if (f_id .ge. 0) then
    call field_get_label(f_id, f_label)
    call field_set_key_str(itempb, keylbl, f_label)
  endif
endif

!===============================================================================
! 4. Set some field keys and number of previous values if needed
!===============================================================================

! Density at the second previous time step for VOF algorithm
! or dilatable algorithm
if (ivofmt.gt.0.or.(idilat.gt.1.and.ipredfl.eq.0).or.irovar.eq.1) then
  call field_set_n_previous(icrom, 2)
  call field_set_n_previous(ibrom, 2)
! The density at the previous time step is required if
! we perform a hydrostatic pressure correction (icalhy=1)
else if (icalhy.eq.1.or.ipthrm.eq.1.or.ippmod(icompf).ge.0.or.idilat.gt.1) then
  call field_set_n_previous(icrom, 1)
  call field_set_n_previous(ibrom, 1)
endif

! Time extrapolation
!-------------------

! Density
call field_get_key_int(icrom, key_t_ext_id, t_ext)
if (t_ext.eq.-1) then
  if (ischtp.eq.1) then
    t_ext = 0
  else if (ischtp.eq.2) then
    ! not extrapolated by default
    t_ext = 0
  endif
  call field_set_key_int(icrom, key_t_ext_id, t_ext)
endif

! Molecular Viscosity
call field_get_key_int(iviscl, key_t_ext_id, t_ext)
if (t_ext.eq.-1) then
  if (ischtp.eq.1) then
    t_ext = 0
  else if (ischtp.eq.2) then
    ! not extrapolated by default
    t_ext = 0
  endif
  call field_set_key_int(iviscl, key_t_ext_id, t_ext)
endif

! Turbulent Viscosity
call field_get_key_int(ivisct, key_t_ext_id, t_ext)
if (t_ext.eq.-1) then
  if (ischtp.eq.1) then
    t_ext = 0
  else if (ischtp.eq.2) then
    ! not extrapolated by default
    t_ext = 0
  endif
  call field_set_key_int(ivisct, key_t_ext_id, t_ext)
endif

! Specific heat
if (icp.ge.0) then
  call field_get_key_int(icp, key_t_ext_id, t_ext)
  if (t_ext.eq.-1) then
    if (ischtp.eq.1) then
      t_ext = 0
    else if (ischtp.eq.2) then
      ! not extrapolated by default
      t_ext = 0
    endif
  endif
  call field_set_key_int(icp, key_t_ext_id, t_ext)
endif

! Scalar diffusivity time extrapolation
do iscal = 1, nscal
  ! Diffusivity of scalars
  call field_get_key_int(ivarfl(isca(iscal)), kivisl, f_id)
  if (f_id.ge.0) then
    call field_get_key_int(f_id, key_t_ext_id, t_ext)
    if (t_ext.eq.-1) then
      if (ischtp.eq.1) then
        t_ext = 0
      else if (ischtp.eq.2) then
        ! Pour le moment par defaut on ne prend pas l'ordre 2
        t_ext = 0
      endif
    endif
    call field_set_key_int(f_id, key_t_ext_id, t_ext)
  endif
enddo

! If time extrapolation, set previous values
do f_id = 0, nfld - 1
  call field_get_key_int(f_id, key_t_ext_id, t_ext)
  if (t_ext.gt.0) then
    call field_get_n_previous(iviscl, n_prev)
    if (n_prev .lt. 1) then
      call field_set_n_previous(iviscl, 1)
    endif
  endif
enddo

! Map pointers

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

!---
! Formats
!---

 7040 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP AT THE INITIAL DATA VERIFICATION'          ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The field ', i10, ' represents the variance'               ,/,&
'@    of fluctuations of the field ', i10                      ,/,&
'@    according to value of keyword first_moment_id'           ,/,&
'@'                                                            ,/,&
'@  The diffusivity_id keyword must not be set'                ,/,&
'@  It will be automatically set equal to that of the'         ,/,&
'@    associated scalar ',i10                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be executed.'                       ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

 7041 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP AT THE INITIAL DATA VERIFICATION'          ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The field ', i10, ' represents the variance'               ,/,&
'@    of fluctuations of the field ', i10                      ,/,&
'@    according to value of keyword first_moment_id'           ,/,&
'@'                                                            ,/,&
'@  The turbulent_diffusivity_id keyword must not be set'      ,/,&
'@  It will be automatically set equal to that of the'         ,/,&
'@    associated scalar ',i10                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be executed.'                       ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

 7042 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP AT THE INITIAL DATA VERIFICATION'          ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The field ', i10, ' represents the variance'               ,/,&
'@    of fluctuations of the field ', i10                      ,/,&
'@    according to value of keyword first_moment_id'           ,/,&
'@'                                                            ,/,&
'@  The sgs_scalar_flux_coef_id keyword must not be set'       ,/,&
'@  It will be automatically set equal to that of the'         ,/,&
'@    associated scalar ',i10                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be executed.'                       ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

return
end subroutine addfld
