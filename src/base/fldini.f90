!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2013 EDF S.A., France

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
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, ivar, iprop
integer          imom, idtnm
integer          keylog
integer          keyvis, keylbl, keycpl, iflid, ikeyid, ikeyvl, iopchr
integer          keysca, keyvar, kscmin, kscmax, kdiftn
integer          nfld, itycat, ityloc, idim1, idim3
logical          ilved, iprev, inoprv
integer          ifvar(nvppmx), iapro(npromx)
integer          f_id, kscavr

character*80     name
character*32     name1, name2, name3
character*80     f_name
character*80     fname(nvppmx)

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

! Key id for scalar id
call field_get_key_id("scalar_id", keysca)

! Key id for variable id
call field_get_key_id("variable_id", keyvar)

! Key id for scamin and scamax
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

! If a scalar is a variance, store the id of the parent scalar
call field_get_key_id("first_moment_id", kscavr)

! Key id for diffusivity tensor
call field_get_key_id("diffusivity_tensor", kdiftn)

! Postprocessing level for variables
iopchr = 1

!===============================================================================
! 2. Mapping for post-processing
!===============================================================================

! Velocity and pressure
!----------------------

ivar = ipr
name = 'pressure'
call field_create(name, itycat, ityloc, idim1, ilved, iprev, ivarfl(ivar))
call field_set_key_str(ivarfl(ivar), keylbl, nomvar(ipprtp(ivar)))
if (ichrvr(ipprtp(ivar)) .eq. 1) then
  call field_set_key_int(ivarfl(ivar), keyvis, iopchr)
endif
if (ilisvr(ipprtp(ivar)) .eq. 1) then
  call field_set_key_int(ivarfl(ivar), keylog, 1)
endif

ivar = iu
name = 'velocity'
call field_create(name, itycat, ityloc, idim3, ilved, iprev, ivarfl(iu))

! Change label for velocity to remove trailing coordinate name
name = nomvar(ipprtp(iu))
name1 = name(1:32)
name = nomvar(ipprtp(iv))
name2 = name(1:32)
name = nomvar(ipprtp(iw))
name3 = name(1:32)
call fldsnv (name1, name2, name3)
!==========
call field_set_key_str(ivarfl(ivar), keylbl, name1)
if (ichrvr(ipprtp(ivar)) .eq. 1) then
  call field_set_key_int(ivarfl(ivar), keyvis, iopchr)
endif
if (ilisvr(ipprtp(ivar)) .eq. 1) then
  call field_set_key_int(ivarfl(ivar), keylog, 1)
endif
call field_set_key_int(ivarfl(ivar), keycpl, 1)

! All components point to same field
ivarfl(iv) = ivarfl(iu)
ivarfl(iw) = ivarfl(iu)

! Turbulence
!-----------

nfld = 0

if (itytur.eq.2) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  fname(nfld) = 'k'
  nfld = nfld + 1
  ifvar(nfld) = iep
  fname(nfld) = 'epsilon'
elseif (itytur.eq.3) then
  nfld = nfld + 1
  ifvar(nfld) = ir11
  fname(nfld) = 'r11'
  nfld = nfld + 1
  ifvar(nfld) = ir22
  fname(nfld) = 'r22'
  nfld = nfld + 1
  ifvar(nfld) = ir33
  fname(nfld) = 'r33'
  nfld = nfld + 1
  ifvar(nfld) = ir12
  fname(nfld) = 'r12'
  nfld = nfld + 1
  ifvar(nfld) = ir13
  fname(nfld) = 'r13'
  nfld = nfld + 1
  ifvar(nfld) = ir23
  fname(nfld) = 'r23'
  nfld = nfld + 1
  ifvar(nfld) = iep
  fname(nfld) = 'epsilon'
  if (iturb.eq.32) then
    nfld = nfld + 1
    ifvar(nfld) = ial
    fname(nfld) = 'alpha'
  endif
elseif (itytur.eq.5) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  fname(nfld) = 'k'
  nfld = nfld + 1
  ifvar(nfld) = iep
  fname(nfld) = 'epsilon'
  nfld = nfld + 1
  ifvar(nfld) = iphi
  fname(nfld) = 'phi'
  if (iturb.eq.50) then
    nfld = nfld + 1
    ifvar(nfld) = ifb
    fname(nfld) = 'f_bar'
  elseif (iturb.eq.51) then
    nfld = nfld + 1
    ifvar(nfld) = ial
    fname(nfld) = 'alpha'
  endif
elseif (iturb.eq.60) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  fname(nfld) = 'k'
  nfld = nfld + 1
  ifvar(nfld) = iomg
  fname(nfld) = 'omega'
elseif (iturb.eq.70) then
  nfld = nfld + 1
  ifvar(nfld) = inusa
  fname(nfld) = 'nu_tilda'
endif

! Map fields

do ii = 1, nfld
  ivar = ifvar(ii)
  name = fname(ii)
  call field_create(name, itycat, ityloc, idim1, ilved, iprev, ivarfl(ivar))
  call field_set_key_str(ivarfl(ivar), keylbl, nomvar(ipprtp(ivar)))
  if (ichrvr(ipprtp(ivar)) .eq. 1) then
    call field_set_key_int(ivarfl(ivar), keyvis, iopchr)
  endif
  if (ilisvr(ipprtp(ivar)) .eq. 1) then
    call field_set_key_int(ivarfl(ivar), keylog, 1)
  endif
enddo

nfld = 0

! Mesh velocity
!--------------

if (iale.eq.1) then
  ivar = iuma
  name = 'mesh_velocity'
  call field_create(name, itycat, ityloc, idim3, ilved, iprev, ivarfl(ivar))
  name = nomvar(ipprtp(iuma))
  name1 = name(1:32)
  name = nomvar(ipprtp(ivma))
  name2 = name(1:32)
  name = nomvar(ipprtp(iwma))
  name3 = name(1:32)
  call fldsnv (name1, name2, name3)
  !==========
  call field_set_key_str(ivarfl(ivar), keylbl, name1)
  if (ichrvr(ipprtp(ivar)) .eq. 1) then
    call field_set_key_int(ivarfl(ivar), keyvis, iopchr)
  endif
  if (ilisvr(ipprtp(ivar)) .eq. 1) then
    call field_set_key_int(ivarfl(ivar), keylog, 1)
  endif
  call field_set_key_int(ivarfl(ivar), keycpl, 1)
  ivarfl(ivma) = ivarfl(iuma)
  ivarfl(iwma) = ivarfl(iuma)
endif

! User variables
!---------------

do ii = 1, nscal

  if (isca(ii) .gt. 0) then
    ivar = isca(ii)
    if (ii .eq. iscalt) then
      if (iscsth(iscalt) .eq. 2) then
        name = 'enthalpy'
      else
        if (iscalt.eq.ienerg) then
          name = 'total energy'
        else
          name = 'temperature'
        endif
      endif
    else
      name = nomvar(ipprtp(ivar))
    endif

    ! Test if the field has already been defined
    call field_get_id_try(trim(name), f_id)

    ! If not already created
    if (f_id.eq.-1) then
      call field_create(name, itycat, ityloc, idim1, ilved, iprev, ivarfl(ivar))
      call field_set_key_str(ivarfl(ivar), keylbl, nomvar(ipprtp(ivar)))
      if (ichrvr(ipprtp(ivar)).eq.1) then
        call field_set_key_int(ivarfl(ivar), keyvis, iopchr)
      endif
      if (ilisvr(ipprtp(ivar)) .eq. 1) then
        call field_set_key_int(ivarfl(ivar), keylog, 1)
      endif
      ! Set min and max clipping
      call field_set_key_double(ivarfl(ivar), kscmin, scamin(ii))
      call field_set_key_double(ivarfl(ivar), kscmax, scamax(ii))

    ! It already exists
    else
      ivarfl(ivar) = f_id
    endif

    ! Set the "scalar_id" key word (inverse of isca(ii))
    call field_set_key_int(ivarfl(ivar), keysca, ii)

    if (ityturt(ii).gt.0) then
      f_name = trim(name)//'_turbulent_flux'
      call field_create(f_name, itycat, ityloc, idim3, .true., iprev, iflid)
      call field_set_key_int(iflid, keycpl, 1)
      ! Tensorial diffusivity
      call field_set_key_int(iflid, kdiftn, 6)
      if (ichrvr(ipprtp(ivar)) .eq. 1) then
        call field_set_key_int(iflid, keyvis, iopchr)
      endif
      if (ilisvr(ipprtp(ivar)) .eq. 1) then
        call field_set_key_int(iflid, keylog, 1)
      endif
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
  ! Set the "variable_id" key word (inverse of ivarfl(ivar))
  call field_set_key_int(ivarfl(ivar), keyvar, ivar)
  ! Key word: tensorial diffusivity
  call field_set_key_int(ivarfl(ivar), kdiftn, idften(ivar))
enddo

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
  name = nomvar(ipppro(iprop))
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
  if (ichrvr(ipppro(iprop)) .eq. 1) then
    call field_set_key_int(iprpfl(iprop), keyvis, ichrvr(ipppro(iprop)))
  endif
  if (ilisvr(ipppro(iprop)) .eq. 1) then
    call field_set_key_int(iprpfl(iprop), keylog, ilisvr(ipppro(iprop)))
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
call field_set_key_str(iflid, keylbl, nomvar(ippdt))
if (idtvar.eq.2.and.ichrvr(ippdt).gt.0) then
  call field_set_key_int(iflid, keyvis, ichrvr(ippdt))
endif
if (idtvar.eq.2.and.ilisvr(ippdt).gt.0) then
  call field_set_key_int(iflid, keylog, ilisvr(ippdt))
endif

! Transient velocity/pressure coupling

if (ipucou.ne.0) then
  name = 'tpucou'
  call field_create(name, itycat, ityloc, idim3, ilved, inoprv, iflid)
  ! Change label to remove trailing coordinate name
  name = nomvar(ipptx)
  name1 = name(1:32)
  name = nomvar(ippty)
  name2 = name(1:32)
  name = nomvar(ipptz)
  name3 = name(1:32)
  call fldsnv (name1, name2, name3)
  !==========
  call field_set_key_str(iflid, keylbl, name1)
endif
if (ichrvr(ipptx).gt.0) then
  call field_set_key_int(iflid, keyvis, ichrvr(ipptx))
endif
if (ilisvr(ipptx).gt.0) then
  call field_set_key_int(iflid, keylog, ilisvr(ipptx))
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

return

end subroutine
