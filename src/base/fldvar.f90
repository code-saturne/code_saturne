!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2025 EDF S.A.
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

!> \file fldvar.f90
!> \brief Variables definition initialization, according to calculation type
!> selected by the user.
!
!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine fldvar() &
 bind(C, name='cs_f_fldvar')

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use parall
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================
! Interfaces
!===============================================================================

procedure() :: fldvar_check_nvar

interface

  function cs_parameters_n_added_variables() result(n) &
    bind(C, name='cs_parameters_n_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int) :: n
  end function cs_parameters_n_added_variables

end interface

!===============================================================================
! Calcul de nscapp
! Verification du nombre de scalaires
! Calcul de nscal

!  A la sortie de cette section, NSCAL, NSCAUS et NSCAPP sont connus.
!  On en profite aussi pour remplir ITYTUR et ITYTURT puisque ITURB et ITURT
!    viennent d'etre definis.
!===============================================================================

! Define main variables
!======================

! Number of user variables

nscaus = cs_parameters_n_added_variables()

if (nscaus.gt.nscamx) then
  ! Since this error should be rare and Fortran code is deprecated,
  ! We do not try to ensure that the message here is ordered
  ! relative to C (as this would add a wrapper), so we prefer to send
  ! it to the error output (unit 0) in this case.
  if (irangp .lt. 1) then
    write(0, 6011) nscaus, nscamx, nscamx, nscaus
  endif
  call csexit (1)
endif

return

!===============================================================================
! 5. FORMATS
!===============================================================================

 6011 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR: STOP AT THE INITIAL DATA VERIFICATION',            /,&
'@    =====',                                                   /,&
'@     NUMBER OF SCALARS TOO LARGE',                            /,&
'@',                                                            /,&
'@  The number of users scalars',                               /,&
'@  requested                          is   NSCAUS = ', i10,    /,&
'@  The total number of scalars',                               /,&
'@    allowed    in   paramx           is   NSCAMX = ', i10,    /,&
'@',                                                            /,&
'@  The maximmum value allowed of   NSCAUS',                    /,&
'@                          is in   NSCAMX        = ', i10,     /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Verify   NSCAUS.',                                          /,&
'@',                                                            /,&
'@  NSCAMX must be at least     ', i10,                         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

!===============================================================================
! 5. FIN
!===============================================================================

return
end subroutine

!===============================================================================
! Local functions
!===============================================================================

!===============================================================================

!> \brief add fields defining user solved scalar variables,
!>        with default options
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine add_user_scalar_fields() &
 bind(C, name='cs_f_add_user_scalar_fields')

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar

! Arguments

! Local variables

integer  iscal, nfld1, nfld2
integer  dim, id, ii, ivar

integer :: keyvar, keysca

!===============================================================================
! Interfaces
!===============================================================================

interface

  ! Interface to C function building user-defined variables

  subroutine cs_parameters_create_added_variables() &
    bind(C, name='cs_parameters_create_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_create_added_variables

end interface

!===============================================================================

! Create fields

call field_get_n_fields(nfld1)

call cs_parameters_create_added_variables

call field_get_n_fields(nfld2)

! Now map those fields

iscal = 0

call field_get_key_id("scalar_id", keysca)
call field_get_key_id("variable_id", keyvar)

do id = nfld1, nfld2 - 1

  call field_get_dim(id, dim)

  if (dim.ne.1 .and. dim.ne.3) cycle

  iscal = iscal + 1

  ivar = nvar + 1
  nvar = nvar + dim
  nscal = nscal + 1

  ! Check we have enough slots
  call fldvar_check_nvar

  isca(iscal) = ivar
  ivarfl(ivar) = id

  call field_set_key_int(id, keyvar, ivar)
  call field_set_key_int(id, keysca, iscal)

  if (dim .gt. 1) then
    do ii = 2, dim
      ivarfl(ivar + ii - 1) = id
    enddo
  endif

enddo

return

!---
! Formats
!---

end subroutine add_user_scalar_fields

!===============================================================================
!
!> \brief add field indexes associated with a new non-user solved
!>        scalar variable, with default options
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  f_id           field id
!> \param[out] iscal          scalar id for defined field
!_______________________________________________________________________________

subroutine add_model_field_indexes &
 ( f_id, iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use field
use cs_c_bindings, only: csexit

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, fldvar_check_nscapp

! Arguments

integer, intent(in)  :: f_id
integer, intent(out) :: iscal

! Local variables

integer  dim, ivar, ii

integer, save :: keyvar = -1
integer, save :: keysca = -1

! Get field dimension

call field_get_dim(f_id, dim)

if (keysca.lt.0) then
  call field_get_key_id("scalar_id", keysca)
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim
nscal = nscal + 1
iscal = nscaus + nscapp + 1
nscapp = nscapp + 1

! Check we have enough slots
call fldvar_check_nvar
call fldvar_check_nscapp

isca(iscal) = ivar

do ii = 1, dim
  ivarfl(ivar + ii - 1) = f_id
enddo

call field_set_key_int(f_id, keyvar, ivar)
call field_set_key_int(f_id, keysca, iscal)

return

end subroutine add_model_field_indexes

!===============================================================================

!> \brief Check nvarmx is sufficient for the required number of variables.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldvar_check_nvar

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use parall
use cs_c_bindings, only: csexit

!===============================================================================

implicit none

! Arguments

! Local variables

if (nvar .gt. nvarmx) then
  ! Since this error should be rare and Fortran code is deprecated,
  ! We do not try to ensure that the message here is ordered
  ! relative to C (as this would add a wrapper), so we prefer to send
  ! it to the error output (unit 0) in this case.
  if (irangp .lt. 1) then
    write(0, 1000) nvar, nvarmx
  endif
  call csexit (1)
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP',              /,&
'@    ======',                                                  /,&
'@     NUMBER OF VARIABLES TOO LARGE',                          /,&
'@',                                                            /,&
'@  The type of calculation defined',                           /,&
'@    corresponds to a number of variables NVAR  >= ', i10,     /,&
'@  The maximum number of variables allowed',                   /,&
'@                      in   paramx   is  NVARMX  = ', i10,     /,&
'@',                                                            /,&
'@  The calculation cannot be executed',                        /,&
'@',                                                            /,&
'@  Verify   parameters.',                                      /,&
'@',                                                            /,&
'@  If NVARMX is increased, the code must be reinstalled.',     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

end subroutine fldvar_check_nvar

!===============================================================================

!> \brief Check nscamx is sufficient for the required number of model scalars.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldvar_check_nscapp

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use parall
use cs_c_bindings, only: csexit

!===============================================================================

implicit none

! Arguments

! Local variables

if ((nscaus+nscapp).gt.nscamx) then
  ! Since this error should be rare and Fortran code is deprecated,
  ! We do not try to ensure that the message here is ordered
  ! relative to C (as this would add a wrapper), so we prefer to send
  ! it to the error output (unit 0) in this case.
  if (irangp .lt. 1) then
    write(0, 1000) nscaus,nscamx,nscamx-nscaus
  endif
  call csexit (1)
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP              ',/,&
'@    ======                                                  ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested                       is   NSCAUS = ', i10    ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.f90       is   NSCAMX = ', i10    ,/,&
'@                                                            ',/,&
'@  The maximum value possible for NSCAPP                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAUS = ', i10    ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  If NSCAMX is increased, the code must be reinstalled.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine fldvar_check_nscapp

!===============================================================================
! C bindings (reverse)
!===============================================================================

!-------------------------------------------------------------------------------
!> \brief add field indexes associated with a new non-user solved
!>        variable, with default options
!
!> \param[in]  f_id    field id

!> \result             variable number for defined field
!-------------------------------------------------------------------------------

function cs_add_variable_field_indexes(f_id) result(ivar) &
  bind(C, name='cs_add_variable_field_indexes')

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use field
use cs_c_bindings, only: csexit

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar

! Arguments

integer(c_int), value :: f_id
integer(c_int) :: ivar

! Local variables

integer  f_id0
integer  dim, ii, keyvar

f_id0 = f_id

! Get field dimension

call field_get_dim(f_id0, dim)

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

do ii = 1, dim
  ivarfl(ivar + ii - 1) = f_id
enddo

call field_get_key_id("variable_id", keyvar)
call field_set_key_int(f_id, keyvar, ivar)

end function cs_add_variable_field_indexes

!-------------------------------------------------------------------------------
!> \brief add field indexes associated with a new non-user solved
!>        variable, with default options
!
!> \param[in]  f_id    field id

!> \result             scalar number for defined field
!-------------------------------------------------------------------------------

function cs_c_add_model_field_indexes(f_id) result(iscal) &
  bind(C, name='cs_add_model_field_indexes')

  use, intrinsic :: iso_c_binding
  use cs_c_bindings

  implicit none

  procedure() :: add_model_field_indexes

  ! Arguments

  integer(c_int), value :: f_id
  integer(c_int) :: iscal

  ! Local variables

  integer f_id0, iscal0

  f_id0 = f_id

  call add_model_field_indexes(f_id0, iscal0)

  iscal = iscal0

end function cs_c_add_model_field_indexes

!---------------------------------------------------------------------------
