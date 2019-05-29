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

subroutine fldtri

!===============================================================================
! Purpose:
! --------

! Map fields

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________.____._____.________________________________________________.

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use cstphy
use numvar
use dimens, only: nvar
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
use cplsat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nscal

! Local variables

integer          ii, ivar
integer          nfld
integer          f_id

integer          ifvar(nvarmx)

character(len=80) :: fname

integer, save :: ipass = 0

logical :: has_exch_bc

type(var_cal_opt) :: vcopt

!===============================================================================


!===============================================================================
! 1. Initialisation
!===============================================================================

nfld = 0

ipass = ipass + 1

!===============================================================================
! 2. Mapping for post-processing
!===============================================================================

! Velocity and pressure
!----------------------

ivar = ipr

if (ipass .eq. 1) then
  call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
  call field_init_bc_coeffs(ivarfl(ivar))
endif

ivar = iu

if (ipass.eq.1) then
  if (ippmod(icompf).ge.0) then
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .true., .false.)
  else
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
  endif
  call field_init_bc_coeffs(ivarfl(ivar))
endif

! Void fraction for VOF algo.
!----------------------------

if (ivofmt.ge.0) then
  ivar = ivolf2
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
    call field_init_bc_coeffs(ivarfl(ivar))
  endif
endif

! Turbulence
!-----------

if (itytur.eq.2) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iep
elseif (itytur.eq.3) then
  if (irijco.eq.1) then
    nfld = nfld + 1
    ifvar(nfld) = ir11
    nfld = nfld + 1
    ifvar(nfld) = ir22
    nfld = nfld + 1
    ifvar(nfld) = ir33
    nfld = nfld + 1
    ifvar(nfld) = ir12
    nfld = nfld + 1
    ifvar(nfld) = ir23
    nfld = nfld + 1
    ifvar(nfld) = ir13
    nfld = nfld + 1
    ifvar(nfld) = iep
  else
    nfld = nfld + 1
    ifvar(nfld) = ir11
    nfld = nfld + 1
    ifvar(nfld) = ir22
    nfld = nfld + 1
    ifvar(nfld) = ir33
    nfld = nfld + 1
    ifvar(nfld) = ir12
    nfld = nfld + 1
    ifvar(nfld) = ir23
    nfld = nfld + 1
    ifvar(nfld) = ir13
    nfld = nfld + 1
    ifvar(nfld) = iep
  endif
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

! Map fields

do ii = 1, nfld
  ivar = ifvar(ii)
  if (ipass .eq. 1) then
    if (itytur.eq.3 ) then
      if (irijco.eq.1) then
        if(ivar.eq.irij) then
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .true., .false., .false.)
        else if (ivar.gt.ir13) then
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
        endif
      else
        if (ivar.ge.ir11 .and. ivar.le.ir13) then
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .true., .false., .false.)
        else
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
        endif
      endif
    else
      call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
    endif
    call field_init_bc_coeffs(ivarfl(ivar))
  endif
enddo

nfld = 0

! Mesh velocity
!--------------

! ALE legacy solver
if (iale.eq.1) then
  ivar = iuma
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., .false.)
    call field_init_bc_coeffs(ivarfl(ivar))
  endif
endif

! Wall distance
!--------------

call field_get_id_try("wall_distance", f_id)

if (f_id.ne.-1) then
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(f_id, .true., .false., .false., .false.)
    call field_init_bc_coeffs(f_id)
  endif
endif

call field_get_id_try("wall_yplus", f_id)

if (f_id.ne.-1) then
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(f_id, .true., .false., .false., .false.)
    call field_init_bc_coeffs(f_id)
  endif
endif

call field_get_id_try("z_ground", f_id)

if (f_id.ne.-1) then
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(f_id, .true., .false., .false., .false.)
    call field_init_bc_coeffs(f_id)
  endif
endif

call field_get_id_try("porosity_w_field", f_id)

if (f_id.ne.-1) then
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(f_id, .true., .false., .false., .false.)
    call field_init_bc_coeffs(f_id)
  endif
endif


! User variables
!---------------

nscal = nscaus + nscapp

do ii = 1, nscal
  if (isca(ii) .gt. 0) then
    ivar = isca(ii)
    has_exch_bc = .false.
    call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
    if (vcopt%icoupl.gt.0) has_exch_bc = .true.

    if (ipass .eq. 1) then
      if (ippmod(icompf).ge.0 .and. ii.eq.ienerg) then
        call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .true., has_exch_bc)
      else
        call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false., has_exch_bc)
      endif
      call field_init_bc_coeffs(ivarfl(ivar))
      ! Boundary conditions of the turbulent fluxes T'u'
      if (ityturt(ii).eq.3) then
        call field_get_name(ivarfl(ivar), fname)
        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)
        call field_allocate_bc_coeffs(f_id, .true., .true., .false., .false.)
        call field_init_bc_coeffs(f_id)
      endif
      ! Elliptic Blending (AFM or DFM)
      if (iturt(ii).eq.11 .or. iturt(ii).eq.21 .or. iturt(ii).eq.31) then
        call field_get_name(ivarfl(ivar), fname)
        call field_get_id(trim(fname)//'_alpha', f_id)
        call field_allocate_bc_coeffs(f_id, .true., .false., .false., .false.)
        call field_init_bc_coeffs(f_id)
      endif
    endif
  endif
enddo

! Reserved fields whose ids are not saved (may be queried by name)
!-----------------------------------------------------------------

! Local time step

! Set previous values for backward n order in time
do ivar = 1, nvar
  ! Here there is no problem if there are multiple
  ! set on non scalar fields.
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  if (vcopt%ibdtso.gt.1) then
    call field_set_n_previous(ivarfl(ivar), vcopt%ibdtso)
  endif
enddo

return
end subroutine fldtri
