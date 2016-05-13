!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2016 EDF S.A., France

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

subroutine fldtri &
!================

 ( nproce ,                                                            &
   dt     , propce )

!===============================================================================
! Purpose:
! --------

! Map fields

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nproce           ! i  ! <-- ! nombre de prop phy aux centres                 !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________.____._____.________________________________________________.

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use numvar

!===============================================================================

implicit none

! Arguments

integer          nproce, nscal
double precision dt(ncelet)
double precision propce(ncelet,*)

! Local variables

logical          lprev
integer          ii, ivar, iprop
integer          iflid, nfld
integer          ipcrom, ipcroa
integer          f_id

integer          ifvar(nvarmx)

character(len=80) :: fname

integer, save :: ipass = 0

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
  call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
  call field_init_bc_coeffs(ivarfl(ivar))
endif

ivar = iu

if (ipass.eq.1) then
  if (ippmod(icompf).ge.0) then
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .true.)
  else
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
  endif
  call field_init_bc_coeffs(ivarfl(ivar))
endif

! Void fraction for cavitation
!-----------------------------

if (icavit.ge.0) then
  ivar = ivoidf
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
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
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .true., .false.)
        else if (ivar.gt.ir13) then
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
        endif
      else
        if (ivar.ge.ir11 .and. ivar.le.ir13) then
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .true., .false.)
        else
          call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
        endif
      endif
    else
      call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
    endif
    call field_init_bc_coeffs(ivarfl(ivar))
  endif
enddo

nfld = 0

! Mesh velocity
!--------------

if (iale.eq.1) then
  ivar = iuma
  if (ipass .eq. 1) then
    call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
    call field_init_bc_coeffs(ivarfl(ivar))
  endif
endif

! User variables
!---------------

nscal = nscaus + nscapp

do ii = 1, nscal
  if (isca(ii) .gt. 0) then
    ivar = isca(ii)
    if (ipass .eq. 1) then
      if (ippmod(icompf).ge.0 .and. ii.eq.ienerg) then
        call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .true.)
      else
        call field_allocate_bc_coeffs(ivarfl(ivar), .true., .false., .false.)
      endif
      call field_init_bc_coeffs(ivarfl(ivar))
      ! Boundary conditions of the turbulent fluxes T'u'
      if (ityturt(ii).eq.3) then
        call field_get_name(ivarfl(ivar), fname)
        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)
        call field_allocate_bc_coeffs(f_id, .true., .true., .false.)
        call field_init_bc_coeffs(f_id)
      endif
    endif
  endif
enddo

! Density field

ipcrom = ipproc(irom)
call field_have_previous(iprpfl(ipcrom), lprev)
if (lprev) then
  ipcroa = ipproc(iroma)
  call field_map_values(iprpfl(ipcrom), propce(1, ipcrom), propce(1, ipcroa))
else
  ipcroa = -1
  call field_map_values(iprpfl(ipcrom), propce(1, ipcrom), propce(1, ipcrom))
endif

! The choice made in VARPOS specifies that we will only be interested in
! properties at cell centers (no mass flux, nor density at the boundary).

f_id = -1
do iprop = 1, nproce
  if (iprop.eq.ipcrom .or. iprop.eq.ipcroa) cycle
  if (iprpfl(iprop).ne.f_id) then
    f_id = iprpfl(iprop)
    call field_map_values(iprpfl(iprop), propce(1, iprop), propce(1, iprop))
  endif
enddo

! Reserved fields whose ids are not saved (may be queried by name)
!-----------------------------------------------------------------

! Local time step

call field_get_id('dt', iflid)
call field_map_values(iflid, dt, dt)

! Set previous values for backward n order in time
do ivar = 1, nvar
  ! Here there is no problem if there are multiple
  ! set on non scalar fields.
  if (ibdtso(ivar).gt.1) then
    call field_set_n_previous(ivarfl(ivar), ibdtso(ivar))
  endif
enddo

return
end subroutine fldtri
