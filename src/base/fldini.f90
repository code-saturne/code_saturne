!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2014 EDF S.A., France

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

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, ivar
integer          keycpl, iflid, ikeyvl
integer          kdiftn
integer          itycat, ityloc, idim1, idim3, idim6
logical          ilved, iprev, inoprv, lprev
integer          f_id, kscavr, f_vis, f_log

character*80     name
character*80     f_name

type(var_cal_opt) vcopt

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
ilved  = .false.   ! not interleaved by default
iprev = .true.     ! variables have previous value
inoprv = .false.   ! variables have no previous value

name = 'log'
call field_get_key_id(name, keylog)

name = 'post_vis'
call field_get_key_id(name, keyvis)

name = 'label'
call field_get_key_id(name, keylbl)

name = 'coupled'
call field_get_key_id(name, keycpl)

! If a scalar is a variance, store the id of the parent scalar
call field_get_key_id("first_moment_id", kscavr)

! Key id for diffusivity tensor
call field_get_key_id("diffusivity_tensor", kdiftn)

!===============================================================================
! 2. Mapping for post-processing
!===============================================================================

! User variables
!---------------

do ii = 1, nscal

  if (isca(ii) .gt. 0) then

    ivar = isca(ii)
    f_id = ivarfl(ivar)

    call field_get_key_int(ivarfl(ivar), keyvis, f_vis)
    call field_get_key_int(ivarfl(ivar), keylog, f_log)

    if (ityturt(ii).gt.0) then
      call field_get_name (f_id, name)
      f_name = trim(name)//'_turbulent_flux'

      if (ityturt(ii).eq.3) then
        itycat = FIELD_INTENSIVE + FIELD_VARIABLE  ! for variables
      else
        itycat = FIELD_INTENSIVE + FIELD_PROPERTY  ! for properties
      endif
      call field_create(f_name, itycat, ityloc, idim3, .true., iprev, iflid)
      if (ityturt(ii).eq.3) then
        call field_set_key_int(iflid, keycpl, 1)
        ! Tensorial diffusivity
        call field_set_key_int(iflid, kdiftn, 6)
      endif
      call field_set_key_int(iflid, keyvis, f_vis)
      call field_set_key_int(iflid, keylog, f_log)
    endif

    ! Additional fields for Drift scalars is done in addfld

  endif

enddo

do ii = 1, nscal
  ! If it is a variance, store the id of the parent scalar
  if (iscavr(ii).gt.0) then
    iflid = ivarfl(isca(iscavr(ii)))
    call field_set_key_int(ivarfl(ivar), kscavr, iflid)
  endif
enddo

do ivar = 1, nvar
  ! Key word: tensorial diffusivity
  call field_set_key_int(ivarfl(ivar), kdiftn, idften(ivar))
enddo

! Copy field calculation options into the field structure
do ivar = 1, nvar
  if (ivar.ne.iv.and.ivar.ne.iw) then
    vcopt%iwarni= iwarni(ivar)
    vcopt%iconv = iconv (ivar)
    vcopt%istat = istat (ivar)
    vcopt%idiff = idiff (ivar)
    vcopt%idifft= idifft(ivar)
    vcopt%idften= idften(ivar)
    vcopt%iswdyn= iswdyn(ivar)
    vcopt%ischcv= ischcv(ivar)
    vcopt%isstpc= isstpc(ivar)
    vcopt%nswrgr= nswrgr(ivar)
    vcopt%nswrsm= nswrsm(ivar)
    vcopt%imrgra= imrgra
    vcopt%imligr= imligr(ivar)
    vcopt%ircflu= ircflu(ivar)
    vcopt%thetav= thetav(ivar)
    vcopt%blencv= blencv(ivar)
    vcopt%epsilo= epsilo(ivar)
    vcopt%epsrsm= epsrsm(ivar)
    vcopt%epsrgr= epsrgr(ivar)
    vcopt%climgr= climgr(ivar)
    vcopt%extrag= extrag(ivar)
    vcopt%relaxv= relaxv(ivar)

    call field_set_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  endif
enddo

! Density field
!--------------

itycat = FIELD_INTENSIVE + FIELD_PROPERTY

! The boundary density at the previous time step is not required
! if we perform a hydrostatic pressure correction (icalhy=1)
lprev = .false.
if (iroext.gt.0.or.idilat.gt.1) then
  lprev = .true.
endif

f_name = 'boundary_density'
ityloc = 3 ! boundary faces
call field_create(f_name, itycat, ityloc, idim1, ilved, lprev, ibrom)

! For now, base logging on that of cell density, and do not postprocess
! boundary density by default

call field_get_key_int(icrom, keylog, ikeyvl)
call field_set_key_int(ibrom, keylog, ikeyvl)

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
! might be defined afterwards in addfld.f90
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
! might be defined afterwards in addfld.f90
do ivar = 1, nvar
  call field_set_key_int(ivarfl(ivar), kbmasf, f_id)
enddo

!====================================================================

! Cavitation: the void fraction has its spectific convective flux

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

!====================================================================

! Combustion

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
  call field_create('yplus', itycat, ityloc, idim1, ilved, inoprv, iyplbr)
endif

! Porosity fields

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! cells
ilved = .true.

if (iporos.ge.1) then
  name = 'porosity'
  call field_create(name, itycat, ityloc, idim1, ilved, inoprv, ipori)
  if (iporos.eq.2) then
    name = 'tensorial_porosity'
    call field_create(name, itycat, ityloc, idim6, ilved, inoprv, iporf)
  endif
endif

return

end subroutine
