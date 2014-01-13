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

!===============================================================================
! Purpose:
! --------

!> \file addfld.f90
!>
!> \brief Add additional fields based on user options.
!>
!> If the user has activated a drift for a scalar for instance,
!> additional fields are created, such an additional mass flux.
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

integer          ii
integer          iscdri, icla, iclap
integer          iflid, iopchr
integer          nfld, itycat, ityloc, idim1, idim3
integer          keyccl, keydri, kdiftn
logical          ilved, iprev, inoprv
integer          f_id

character*80     name
character*80     f_name

!===============================================================================

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
ilved  = .false.   ! not interleaved by default
iprev  = .true.    ! variables have previous value
inoprv = .false.   ! variables have no previous value
iopchr = 1         ! Postprocessing level for variables

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key id for diffusivity tensor
call field_get_key_id("diffusivity_tensor", kdiftn)

! Number of fields
call field_get_n_fields(nfld)

!===============================================================================
! 1. Initialisation
!===============================================================================

! Add mass flux for scalar with a drift (one mass flux per class)

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

return

end subroutine addfld
