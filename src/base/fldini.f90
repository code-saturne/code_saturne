!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2015 EDF S.A., France

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
use lagpar
use lagdim
use lagran
use ihmpre
use radiat
use cplsat
use mesh
use field
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, ivar
integer          keycpl, iflid
integer          kdiftn, kturt, kfturt, keyvar
integer          itycat, ityloc, idim1, idim3, idim6
logical          ilved, iprev, inoprv
integer          f_id, kscavr, f_vis, f_log, f_dften, f_type, f_loc
integer          kislts, ifctsl
integer          iopchr
integer          iscdri, icla, iclap
integer          keyccl, keydri
integer          idfm, iggafm, nfld

character(len=80) :: name, f_name

type(var_cal_opt) vcopt
type(gas_mix_species_prop) sasp

!===============================================================================

!===============================================================================
! 1. Initialisation
!===============================================================================

! The itycat variable is used to define field categories. It is used in Fortran
! code with hard-coded values, but in the C API, those values are based on
! (much clearer) category mask definitions in cs_field.h.

itycat = FIELD_INTENSIVE + FIELD_VARIABLE  ! for most variables
ityloc = 1 ! variables defined on cells
idim1  = 1
idim3  = 3
idim6  = 6
ilved  = .true.    ! interleaved by default
iprev = .true.     ! variables have previous value
inoprv = .false.   ! variables have no previous value

call field_get_key_id('log', keylog)
call field_get_key_id('post_vis', keyvis)
call field_get_key_id('label', keylbl)
call field_get_key_id('coupled', keycpl)
call field_get_key_id("variable_id", keyvar)

! If a scalar is a variance, store the id of the parent scalar
call field_get_key_id("first_moment_id", kscavr)

! Key id for diffusivity tensor
call field_get_key_id("diffusivity_tensor", kdiftn)

! Keys not stored globally

call field_get_key_id('turbulent_flux_model', kturt)
call field_get_key_id('turbulent_flux_id', kfturt)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Number of fields
call field_get_n_fields(nfld)

!===============================================================================
! 2. Mapping for post-processing
!===============================================================================

! User variables
!---------------

do ivar = 1, nvar
  ! Init key word: tensorial diffusivity
  call field_set_key_int(ivarfl(ivar), kdiftn, idften(ivar))
enddo

idfm = 0
iggafm = 0

do ii = 1, nscal

  if (isca(ii) .gt. 0) then

    ivar = isca(ii)
    f_id = ivarfl(ivar)

    call field_get_key_int(ivarfl(ivar), keyvis, f_vis)
    call field_get_key_int(ivarfl(ivar), keylog, f_log)

    call field_get_key_int(ivarfl(ivar), kdiftn, f_dften)

    if (ityturt(ii).gt.0) then
      call field_get_name (f_id, name)
      f_name = trim(name)//'_turbulent_flux'

      if (ityturt(ii).eq.3) then
        itycat = FIELD_INTENSIVE + FIELD_VARIABLE  ! for variables
      else
        itycat = FIELD_INTENSIVE + FIELD_PROPERTY  ! for properties
      endif

      ! GGDH or AFM on current scalar
      ! and if DFM, GGDH on the scalar variance
      iggafm = 1

      call field_create(f_name, itycat, ityloc, idim3, .true., iprev, iflid)
      if (ityturt(ii).eq.3) then
        call field_set_key_int(iflid, keycpl, 1)
        ! Tensorial diffusivity
        call field_set_key_int(iflid, kdiftn, 6)
        idfm = 1
      endif
      call field_set_key_int(iflid, keyvis, f_vis)
      call field_set_key_int(iflid, keylog, f_log)

      call field_set_key_int(ivarfl(ivar), kturt, iturt(ii))
      call field_set_key_int(ivarfl(ivar), kfturt, iflid)

    ! If the user has chosen a tensorial diffusivity
    else if (f_dften.eq.6) then
      idfm = 1
    endif

    ! Reference diffusivity value
    call field_set_key_double(ivarfl(ivar), kvisl0, visls0(ii))

    ! Additional fields for Drift scalars is done in addfld

  endif

enddo

! Reserved fields whose ids are not saved (may be queried by name)
!-----------------------------------------------------------------

! Interior mass flux field
!-------------------------

itycat = FIELD_EXTENSIVE + FIELD_PROPERTY
ityloc = 2 ! inner faces

! Mass flux for the class on interior faces
f_name = 'inner_mass_flux'

if (istmpf.ne.1) then
  call field_create(f_name, itycat, ityloc, idim1, ilved, iprev, f_id)
else
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
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
  call field_create(f_name, itycat, ityloc, idim1, ilved, iprev, f_id)
else
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
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

    ! The comming field are not solved, they are properties
    itycat = FIELD_PROPERTY
    ityloc = 2 ! variables defined on interior faces

    ! Mass flux for the class on interior faces
    f_name = 'inner_mass_flux_'//trim(name)
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
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
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
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

    ! Relaxation time
    f_name = 'drift_tau_'//trim(name)
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)

    ! Set the same visualization options as the scalar
    call field_get_key_int(iflid, keyvis, iopchr)
    if (iopchr.eq.1) then
      call field_set_key_int(f_id, keyvis, iopchr)
    endif

    ! Store the drift velocity
    f_name = 'drift_vel_'//trim(name)
    call field_create(f_name, itycat, ityloc, idim3, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)

    ! Set the same visualization options as the scalar
    call field_get_key_int(iflid, keyvis, iopchr)
    if (iopchr.eq.1) then
      call field_set_key_int(f_id, keyvis, iopchr)
    endif

    ! Interaction time particle--eddies
    if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS)) then
      f_name = 'drift_turb_tau_'//trim(name)
      call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
      call field_set_key_str(f_id, keylbl, f_name)
    endif

    ! Set the same visualization options as the scalar
    call field_get_key_int(iflid, keyvis, iopchr)
    if (iopchr.eq.1) then
      call field_set_key_int(f_id, keyvis, iopchr)
    endif

  endif
enddo


!===============================================================================

! Cavitation: the void fraction has its spectific convective flux
!----------------------------------------------------------------

if (icavit.ge.0) then

  itycat = FIELD_EXTENSIVE + FIELD_PROPERTY

  ityloc = 2  ! inner faces
  f_name = 'inner_void_fraction_flux'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_int(ivarfl(ivoidf), kimasf, f_id)

  ityloc = 3 ! boundary faces
  f_name = 'boundary_void_fraction_flux'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_int(ivarfl(ivoidf), kbmasf, f_id)

endif

!===============================================================================

! Combustion
!-----------

if (iirayo .gt. 0) then

  if (     ippmod(icod3p).eq.1                                    &
      .or. ippmod(iccoal).ge.0                                    &
      .or. (ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3)         &
      .or. (     ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3     &
            .or. ippmod(icolwc).eq.5)) then

    itycat = FIELD_INTENSIVE + FIELD_PROPERTY
    ityloc = 3 ! boundary faces

    call field_create('boundary_ym_fuel',  &
                      itycat, ityloc, idim1, ilved, inoprv, ibym(1))
    call field_create('boundary_ym_oxydizer',  &
                      itycat, ityloc, idim1, ilved, inoprv, ibym(2))
    call field_create('boundary_ym_product',  &
                      itycat, ityloc, idim1, ilved, inoprv, ibym(3))
  endif

endif

!===============================================================================

! Turbulent anisotropic viscosity or user defined tensorial diffusivity for a
! scalar (exclusive or).
!----------------------------------------------------------------------

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! cells
ilved = .true.

if (idfm.eq.1 .or. itytur.eq.3 .and. idirsm.eq.1 &
    .or.darcy_anisotropic_diffusion.eq.1) then
  call field_create('anisotropic_turbulent_viscosity', itycat, ityloc, idim6, &
                    ilved, inoprv, ivsten)
  if (iturb.eq.32.and.iggafm.eq.1) then
    call field_create('anisotropic_turbulent_viscosity_scalar', itycat, &
                      ityloc, idim6, ilved, inoprv, ivstes)
 endif
endif

! Additional fields
!------------------

! Fields used to save postprocessing data

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces

! If postprocessing of boundary temperature or boundary layer Nusselt required
if (ipstdv(ipsttb).gt.0 .or. ipstdv(ipstnu).gt.0) then
  call field_create('tplus', itycat, ityloc, idim1, ilved, inoprv, iflid)
  call field_create('tstar', itycat, ityloc, idim1, ilved, inoprv, iflid)
endif

ilved = .true.

if (ineedf.eq.1) then
  call field_create('boundary_forces', itycat, ityloc, idim3, ilved, inoprv, &
                    iforbr)
endif

if (ipstdv(ipstyp).ne.0) then
  call field_get_id_try('yplus', f_id)

  ! If it already exists with a different location, EXIT
  if (f_id.ge.0) then
    call field_get_location(f_id,f_loc)

    if (ityloc.ne.f_loc) call csexit(1)
    iyplbr = f_id
  else
    call field_create('yplus', itycat, ityloc, idim1, ilved, inoprv, iyplbr)
    call field_set_key_str(iyplbr, keylbl,'Yplus')
  endif
  ! yplus postreated and in the log
  call field_set_key_int(iyplbr, keyvis, 1)
  call field_set_key_int(iyplbr, keylog, 1)
endif

! Postprocessing of slope tests

call field_get_key_id("slope_test_upwind_id", kislts)

itycat = FIELD_POSTPROCESS
ityloc = 1 ! cells
ilved = .true.

do ii = 1, nvar
  f_id = ivarfl(ii)
  call field_get_key_int(f_id, kislts, ifctsl)
  if (ifctsl.eq.0) then
    ! Now create matching field
    if (iconv(ii).gt.0 .and. blencv(ii).gt.0 .and. isstpc(ii).eq.0) then
      ! Build name and label
      call field_get_name(f_id, f_name)
      name  = trim(f_name) // '_slope_upwind'
      call field_create(name, itycat, ityloc, idim1, ilved, inoprv, ifctsl)
      call field_set_key_int(ifctsl, keyvis, 1)
    else
      ifctsl = -1
    endif
    call field_set_key_int(f_id, kislts, ifctsl)
  endif
enddo

!===============================================================================
! 3. Set some field keys
!===============================================================================

! Copy field calculation options into the field structure
do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)

    call field_get_key_int(f_id, keyvar, ivar)

    vcopt%iwarni= iwarni(ivar)
    vcopt%iconv = iconv (ivar)
    vcopt%istat = istat (ivar)
    vcopt%idiff = idiff (ivar)
    vcopt%idifft= idifft(ivar)
    vcopt%idften= idften(ivar)
    vcopt%iswdyn= iswdyn(ivar)
    vcopt%ischcv= ischcv(ivar)
    vcopt%ibdtso= ibdtso(ivar)
    vcopt%isstpc= isstpc(ivar)
    vcopt%nswrgr= nswrgr(ivar)
    vcopt%nswrsm= nswrsm(ivar)
    vcopt%imrgra= imrgra
    vcopt%imligr= imligr(ivar)
    vcopt%ircflu= ircflu(ivar)
    vcopt%iwgrec= iwgrec(ivar)
    vcopt%thetav= thetav(ivar)
    vcopt%blencv= blencv(ivar)
    vcopt%epsilo= epsilo(ivar)
    vcopt%epsrsm= epsrsm(ivar)
    vcopt%epsrgr= epsrgr(ivar)
    vcopt%climgr= climgr(ivar)
    vcopt%extrag= extrag(ivar)
    vcopt%relaxv= relaxv(ivar)

    call field_set_key_struct_var_cal_opt(f_id, vcopt)
  endif
enddo

! Copy field physical properties of species into the field structure

if (ippmod(igmix).ge.0) then

  call cs_parameters_define_field_key_gas_mix

  do f_id = 1, nfld

    call field_get_name (f_id, name)

    if (name.eq.'y_o2') then

      sasp%mol_mas= 0.032d0
      sasp%cp = 930.d0
      sasp%vol_dif = 19.70d0
      sasp%mu_a = 5.086522d-8
      sasp%mu_b = 5.512391d-6
      sasp%lambda_a = 6.2d-5
      sasp%lambda_b = 8.1d-3

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_n2') then

      sasp%mol_mas = 0.028d0
      sasp%cp = 1042.d0
      sasp%vol_dif = 19.70d0
      sasp%mu_a = 4.210130d-8
      sasp%mu_b = 5.494348d-6
      sasp%lambda_a = 6.784141d-5
      sasp%lambda_b = 5.564317d-3

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_he') then

      sasp%mol_mas = 0.004d0
      sasp%cp = 5194.d0
      sasp%vol_dif = 2.67d0
      sasp%mu_a = 18.5752d-6
      sasp%mu_b = 0.0d0
      sasp%lambda_a = 0.144d0
      sasp%lambda_b = 0.0d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_h2') then

      sasp%mol_mas = 0.002d0
      sasp%cp = 14560.d0
      sasp%vol_dif = 6.12d0
      sasp%mu_a = 1.93d-9
      sasp%mu_b = 8.40d-6
      sasp%lambda_a = 4.431d-4
      sasp%lambda_b = 5.334d-2

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    else if (name.eq.'y_h2o_g') then

      sasp%mol_mas = 0.018d0
      sasp%cp = 2060.d0
      sasp%vol_dif = 13.10d0
      sasp%mu_a = 3.8496d-8
      sasp%mu_b = 8.2997d-6
      sasp%lambda_a = 7.6209d-5
      sasp%lambda_b = 0.016949d0

      call field_set_key_struct_gas_mix_species_prop(f_id, sasp)

    endif

  enddo

endif

return

end subroutine
