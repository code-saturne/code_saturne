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

integer          ii, ivar, iprop
integer          imom, idtnm
integer          keycpl, iflid, ikeyid, ikeyvl
integer          kdiftn
integer          nfld, itycat, ityloc, idim1, idim3, idim6
        integer          ipcrom, ipcroa
logical          ilved, iprev, inoprv, lprev
integer          ifvar(nvppmx), iapro(npromx)
integer          f_id, kscavr

character*80     name
character*80     f_name
character*80     fname(nvppmx)

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

! Velocity and pressure
!----------------------

ivar = ipr
call field_set_key_int(ivarfl(ivar), keyvis, ichrvr(ipprtp(ivar)))
call field_set_key_int(ivarfl(ivar), keylog, ilisvr(ipprtp(ivar)))

ivar = iu

call field_set_key_int(ivarfl(ivar), keyvis, ichrvr(ipprtp(ivar)))
call field_set_key_int(ivarfl(ivar), keylog, ilisvr(ipprtp(ivar)))
call field_set_key_int(ivarfl(ivar), keycpl, 1)

! Turbulence
!-----------

nfld = 0

if (itytur.eq.2) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iep
elseif (itytur.eq.3) then
  nfld = nfld + 1
  ifvar(nfld) = ir11
  nfld = nfld + 1
  ifvar(nfld) = ir22
  nfld = nfld + 1
  ifvar(nfld) = ir33
  nfld = nfld + 1
  ifvar(nfld) = ir12
  nfld = nfld + 1
  ifvar(nfld) = ir13
  nfld = nfld + 1
  ifvar(nfld) = ir23
  nfld = nfld + 1
  ifvar(nfld) = iep
  if (iturb.eq.32) then
    nfld = nfld + 1
    ifvar(nfld) = ial
  endif
elseif (itytur.eq.5) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iep
  nfld = nfld + 1
  ifvar(nfld) = iphi
  if (iturb.eq.50) then
    nfld = nfld + 1
    ifvar(nfld) = ifb
  elseif (iturb.eq.51) then
    nfld = nfld + 1
    ifvar(nfld) = ial
  endif
elseif (iturb.eq.60) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iomg
elseif (iturb.eq.70) then
  nfld = nfld + 1
  ifvar(nfld) = inusa
endif

! Set turbulence field options

do ii = 1, nfld
  ivar = ifvar(ii)
  name = fname(ii)
  call field_set_key_int(ivarfl(ivar), keyvis, ichrvr(ipprtp(ivar)))
  call field_set_key_int(ivarfl(ivar), keylog, ilisvr(ipprtp(ivar)))
enddo

nfld = 0

! Mesh velocity
!--------------

if (iale.eq.1) then
  ivar = iuma
  call field_set_key_int(ivarfl(ivar), keyvis, ichrvr(ipprtp(ivar)))
  call field_set_key_int(ivarfl(ivar), keylog, ilisvr(ipprtp(ivar)))
  call field_set_key_int(ivarfl(ivar), keycpl, 1)
endif

! User variables
!---------------

do ii = 1, nscal

  if (isca(ii) .gt. 0) then

    ivar = isca(ii)
    f_id = ivarfl(ivar)

    call field_set_key_int(ivarfl(ivar), keyvis, ichrvr(ipprtp(ivar)))
    call field_set_key_int(ivarfl(ivar), keylog, ilisvr(ipprtp(ivar)))

    if (ityturt(ii).gt.0) then
      call field_get_name (f_id, name)
      f_name = trim(name)//'_turbulent_flux'
      call field_create(f_name, itycat, ityloc, idim3, .true., iprev, iflid)
      call field_set_key_int(iflid, keycpl, 1)
      ! Tensorial diffusivity
      call field_set_key_int(iflid, kdiftn, 6)
      call field_set_key_int(iflid, keyvis, ichrvr(ipprtp(ivar)))
      call field_set_key_int(iflid, keylog, ilisvr(ipprtp(ivar)))
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
enddo




! Density field
!--------------

itycat = FIELD_INTENSIVE + FIELD_PROPERTY

! The density at the previous time step is required if idilat>1 or if
! we perform a hydrostatic pressure correction (icalhy=1)

f_name = 'density'
ityloc = 1 ! cells
ipcrom = ipproc(irom)
if (     iroext.gt.0 .or. icalhy.eq.1 .or. idilat.gt.1              &
    .or. ippmod(icompf).ge.0) then
  lprev = .true.
  ipcroa = ipproc(iroma)
else
  lprev = .false.
  ipcroa = -1
endif

call field_create(f_name, itycat, ityloc, idim1, ilved, lprev, icrom)

iprpfl(ipcrom) = icrom
if (lprev) iprpfl(ipcroa) = -1 ! could set icrom, but avoid this access mode

call field_set_key_str(icrom, keylbl, nomprp(ipcrom))
call field_set_key_int(icrom, keyvis, ichrvr(ipppro(ipcrom)))
call field_set_key_int(icrom, keylog, ilisvr(ipppro(ipcrom)))
call field_set_key_int(icrom, keyipp, ipppro(ipcrom))

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

call field_set_key_int(ibrom, keylog, ilisvr(ipppro(ipcrom)))

! Cell properties
!----------------

ityloc = 1 ! cells

! Flag moments

do ii = 1, npromx
  iapro(ii ) = 0
enddo

! For moments, this key defines the division by time mode
!  = 0: no division
!  > 0: property number for cumulative dt (property)
!  < 0: position in dtcmom of cumulative dt (uniform)

do imom = 1, nbmomt
  ! property id matching moment
  iprop = ipproc(icmome(imom))
  if (idtmom(imom).ne.0) then
    iapro(iprop) = 1
  endif
enddo

! Mark moment accumulators

do imom = 1, nbmomt
  idtnm = idtmom(imom)
  if (idtnm.gt.0) then
    iprop = ipproc(icdtmo(idtnm))
    iapro(iprop) = -idtnm
  endif
enddo

! The choice made in VARPOS specifies that we will only be interested in
! properties at cell centers (no mass flux, nor density at the boundary).

imom = 0
do iprop = 1, nproce
  if (iprop.eq.ipcrom .or. iprop.eq.ipcroa) cycle
  name = nomprp(iprop)
  if (iapro(iprop).eq.0) then
    if (name(1:4) .eq. '    ') then
      write(name, '(a, i3.3)') 'property_', iprop
    endif
    itycat = FIELD_PROPERTY
  else
    if (iapro(iprop).gt.0) then
      imom = imom + 1
      if (name(1:4) .eq. '    ') then
        write(name, '(a, i3.3)') 'moment_', imom
      endif
    else if (iapro(iprop).lt.0) then
      imom = imom + 1
      if (name(1:4) .eq. '    ') then
        write(name, '(a, i3.3)') 'accumulator_', -iapro(iprop)
      endif
    endif
    itycat = FIELD_PROPERTY + FIELD_ACCUMULATOR
  endif
  call field_create(name, itycat, ityloc, idim1, ilved, inoprv, iprpfl(iprop))
  call field_set_key_str(iprpfl(iprop), keylbl, name)
  if (ipppro(iprop).gt.1) then
    call field_set_key_int(iprpfl(iprop), keyvis, ichrvr(ipppro(iprop)))
    call field_set_key_int(iprpfl(iprop), keylog, ilisvr(ipppro(iprop)))
    call field_set_key_int(iprpfl(iprop), keyipp, ipppro(iprop))
  endif
enddo

! Add moment accumulators metadata
!---------------------------------

name = 'moment_dt'
call field_get_key_id(name, ikeyid)

do imom = 1, nbmomt
  ! property id matching moment
  iprop = ipproc(icmome(imom))
  ! dt type and number
  idtnm = idtmom(imom)
  ikeyvl = -1
  if (idtnm.gt.0) then
    ikeyvl = iprpfl(ipproc(icdtmo(idtnm)))
  elseif(idtnm.lt.0) then
    ikeyvl = idtnm - 1
  endif
  call field_set_key_int(iprpfl(iprop), ikeyid, ikeyvl)
enddo

! Reserved fields whose ids are not saved (may be queried by name)
!-----------------------------------------------------------------

itycat = FIELD_INTENSIVE

! Local time step

name = 'dt'
call field_create(name, itycat, ityloc, idim1, ilved, inoprv, iflid)
call field_set_key_str(iflid, keylbl, 'Local Time Step')
if (idtvar.eq.2.and.ichrvr(ippdt).gt.0) then
  call field_set_key_int(iflid, keyvis, ichrvr(ippdt))
endif
if (idtvar.eq.2.and.ilisvr(ippdt).gt.0) then
  call field_set_key_int(iflid, keylog, ilisvr(ippdt))
endif
call field_set_key_int(iflid, keyipp, ippdt)

! Transient velocity/pressure coupling, postprocessing field
! (variant used for computation is a tensorial field, not this one)

if (ipucou.ne.0 .or. ncpdct.gt.0) then
  name = 'dttens'
  call field_create(name, itycat, ityloc, 6, .true., inoprv, idtten)
  call field_set_key_int(idtten, keyipp, ipptx)
endif
if (ichrvr(ipptx).gt.0) then
  call field_set_key_int(idtten, keyvis, ichrvr(ipptx))
endif
if (ilisvr(ipptx).gt.0) then
  call field_set_key_int(idtten, keylog, ilisvr(ipptx))
endif

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

!====================================================================

! Radiative_model

if (iirayo.gt.0) then

  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  ityloc = 3 ! boundary faces

  call field_create('wall_temperature',  &
                    itycat, ityloc, idim1, ilved, inoprv, itparo)
  call field_create('incident_radiative_flux_density',  &
                    itycat, ityloc, idim1, ilved, inoprv, iqinci)
  call field_create('wall_thermal_conductivity',  &
                    itycat, ityloc, idim1, ilved, inoprv, ixlam)
  call field_create('wall_thickness',  &
                    itycat, ityloc, idim1, ilved, inoprv, iepa)
  call field_create('wall_emissivity',  &
                    itycat, ityloc, idim1, ilved, inoprv, ieps)
  call field_create('net_radiative_flux',  &
                    itycat, ityloc, idim1, ilved, inoprv, ifnet)
  call field_create('radiation_convective_flux',  &
                    itycat, ityloc, idim1, ilved, inoprv, ifconv)
  call field_create('radiation_exchange_coefficient',  &
                    itycat, ityloc, idim1, ilved, inoprv, ihconv)

endif

!=====================================================================

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
