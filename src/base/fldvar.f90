!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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
!> \param[out]    nmodpp        number of activated particle physic models
!______________________________________________________________________________

subroutine fldvar ( nmodpp ) &
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
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer(c_int), value ::  nmodpp

! Local variables

integer       ipp
integer       iok

!===============================================================================
! Interfaces
!===============================================================================

procedure() :: add_variable_field, add_cdo_variable_field
procedure() :: ppvarp, add_model_scalar_field
procedure() :: add_user_scalar_fields, fldvar_check_nvar
procedure() :: init_var_cal_opt

interface

  ! Interface to C function returning number of user-defined variables

  function cs_parameters_n_added_variables() result(n) &
    bind(C, name='cs_parameters_n_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int) :: n
  end function cs_parameters_n_added_variables

end interface

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
  endif
enddo

!===============================================================================
! Calcul de nscapp
! Verification du nombre de scalaires
! Construction de iscapp
! Calcul de nscal

!  A la sortie de cette section, NSCAL, NSCAUS et NSCAPP sont connus.
!  On en profite aussi pour remplir ITYTUR et ITYTURT puisque ITURB et ITURT
!    viennent d'etre definis.
!===============================================================================

! Define main variables
!======================

! Velocity (components point to same field)
call map_variable_field_try('velocity', iu)
iv = iu + 1
iw = iv + 1

! Pressure
call map_variable_field_try('pressure', ipr)

! void fraction (VoF algorithm)

call map_variable_field_try('void_fraction', ivolf2)

! Turbulence

if (itytur.eq.2) then
  call map_variable_field_try('k', ik)
  call map_variable_field_try('epsilon', iep)
else if (itytur.eq.3) then
  call map_variable_field_try('rij', irij)
  ! All rij components point to same field
  ir11 = irij
  ir22 = ir11 + 1
  ir33 = ir22 + 1
  ir12 = ir33 + 1
  ir23 = ir12 + 1
  ir13 = ir23 + 1

  call map_variable_field_try('epsilon', iep)
  if (iturb.eq.32) then
    call map_variable_field_try('alpha', ial)
  endif
else if (itytur.eq.5) then
  call map_variable_field_try('k', ik)
  call map_variable_field_try('epsilon', iep)
  call map_variable_field_try('phi', iphi)
  if (iturb.eq.50) then
    call map_variable_field_try('f_bar', ifb)
  else if (iturb.eq.51) then
    call map_variable_field_try('alpha', ial)
  endif
else if (iturb.eq.60) then
  call map_variable_field_try('k', ik)
  call map_variable_field_try('omega', iomg)
else if (iturb.eq.70) then
  call map_variable_field_try('nu_tilda', inusa)
endif

! Mesh velocity with ALE
call map_variable_field_try('mesh_velocity', iuma)
if (iuma.ge.1) then
  ivma = iuma + 1
  iwma = ivma + 1
endif

! Number of user variables

nscaus = cs_parameters_n_added_variables()

! Specific physics variables
call ppvarp

! Thermal model with no specific physics

if (nmodpp.eq.0) then

  if (itherm .eq. 1) then
    call add_model_scalar_field('temperature', 'Temperature', iscalt)
  else if (itherm .eq. 2) then
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
  else if (itherm .eq. 4) then
    call add_model_scalar_field('internal_energy', 'Eint', iscalt)
  endif

endif

call add_user_scalar_fields

! ---> Verifications

iok = 0

if (nscaus.gt.0 .or. nscapp.gt.0) then
  if ((nscaus+nscapp).gt.nscamx) then
    if (nscapp.le.0) then
      write(nfecra,6011) nscaus, nscamx, nscamx, nscaus
    else
      write(nfecra,6012) nscaus,nscapp,nscamx,nscamx-nscapp,nscaus+nscapp
    endif
    iok = iok + 1
  endif
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

return

!===============================================================================
! 5. FORMATS
!===============================================================================

 6011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@  requested                          is   NSCAUS = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx           is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximmum value allowed of   NSCAUS                    ',/,&
'@                          is in   NSCAMX        = ',I10      ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6012 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested                       is   NSCAUS = ',I10     ,/,&
'@  The number of scalars necessary for the specific physics'  ,/,&
'@    with the chosen model is              NSCAPP = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximum value allowed for  NSCAUS                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAPP = ',I10     ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 5. FIN
!===============================================================================

return
end subroutine

!===============================================================================
! Local functions
!===============================================================================

!===============================================================================

!> \brief add field defining a general solved variable, with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name          field name
!> \param[in]  label         field default label, or empty
!> \param[in]  dim           field dimension
!> \param[out] ivar          variable number for defined field
!_______________________________________________________________________________

subroutine add_variable_field &
 ( name, label, dim, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
integer, intent(out)         :: ivar

! Local variables

integer  id, ii

integer, save :: keyvar = -1

! Create field

call variable_field_create(name, label, MESH_LOCATION_CELLS, dim, id)

if (keyvar.lt.0) then
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

ivarfl(ivar) = id

call field_set_key_int(id, keyvar, ivar)

call init_var_cal_opt(id)

if (dim .gt. 1) then
  do ii = 2, dim
    ivarfl(ivar + ii - 1) = id
  enddo
endif

return

end subroutine add_variable_field

!===============================================================================

!> \brief Add a field defining a general solved variable, with default options
!>        This variable is solved with a CDO scheme.
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name          field name
!> \param[in]  label         field default label, or empty
!> \param[in]  dim           field dimension
!> \param[in]  location_id   id of the mesh location where the field is defined
!> \param[in]  has_previous  if greater than 0 then stores previous state
!> \param[out] ivar          variable number for defined field
!_______________________________________________________________________________

subroutine add_cdo_variable_field &
 ( name, label, dim, location_id, has_previous, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim, location_id, has_previous
integer, intent(out)         :: ivar

! Local variables

integer  id, ii

integer, save :: keyvar = -1

! Create field

call variable_cdo_field_create(name, label, location_id, dim, has_previous, id)

if (keyvar.lt.0) then
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

ivarfl(ivar) = id

call field_set_key_int(id, keyvar, ivar)

call init_var_cal_opt(id)

if (dim .gt. 1) then
  do ii = 2, dim
    ivarfl(ivar + ii - 1) = id
  enddo
endif

return

end subroutine add_cdo_variable_field

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

subroutine add_user_scalar_fields

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt

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
  call init_var_cal_opt(id)

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

!> \brief add field defining a non-user solved scalar variable,
!>        with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[out] iscal          scalar number for defined field
!_______________________________________________________________________________

subroutine add_model_scalar_field &
 ( name, label, iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: add_model_field

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(out)         :: iscal

! Local variables

integer  dim

dim = 1

call add_model_field(name, label, dim, iscal)

return

end subroutine add_model_scalar_field

!===============================================================================
!
!> \brief add field defining a non-user solved variable,
!>        with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[in]  dim            field dimension
!> \param[out] iscal          variable number for defined field
!_______________________________________________________________________________

subroutine add_model_field &
 ( name, label, dim, iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

procedure() :: add_model_field_indexes

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
integer, intent(out)         :: iscal

! Local variables

integer  id
integer  location_id

location_id = 1 ! variables defined on cells

! Create field

call variable_field_create(name, label, location_id, dim, id)

call add_model_field_indexes(id, iscal)

return

end subroutine add_model_field

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
!> \param[out] ivar           variable number for defined field
!_______________________________________________________________________________

subroutine add_variable_field_indexes &
 ( f_id, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt
procedure() :: csexit

! Arguments

integer, intent(in)  :: f_id
integer, intent(out) :: ivar

! Local variables

integer  dim, ii

integer, save :: keyvar = -1

! Get field dimension

call field_get_dim(f_id, dim)

if (keyvar.lt.0) then
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

do ii = 1, dim
  ivarfl(ivar + ii - 1) = f_id
enddo

call field_set_key_int(f_id, keyvar, ivar)
call init_var_cal_opt(f_id)

return

end subroutine add_variable_field_indexes

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
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, fldvar_check_nscapp, init_var_cal_opt
procedure() :: csexit

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
iscapp(nscapp) = iscal

do ii = 1, dim
  ivarfl(ivar + ii - 1) = f_id
enddo

call field_set_key_int(f_id, keyvar, ivar)
call field_set_key_int(f_id, keysca, iscal)
call init_var_cal_opt(f_id)

return

end subroutine add_model_field_indexes

!===============================================================================

!> \brief Map field defining a general solved variable to Fortran id

!> If field is not available, value is at 0

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   name          field name
!> \param[out]  ivar          variable id, or 0
!_______________________________________________________________________________

subroutine map_variable_field_try &
 ( name, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

character(len=*), intent(in) :: name
integer, intent(out)         :: ivar

! Local variables

integer f_id
integer keyvar
character(len=len_trim(name)+1, kind=c_char) :: c_name

! Find field id

ivar = 0

c_name = trim(name)//c_null_char

f_id = cs_f_field_id_by_name_try(c_name)

if (f_id .ge. 0) then
  call field_get_key_id("variable_id", keyvar)
  call field_get_key_int(f_id, keyvar, ivar)
endif

return

end subroutine map_variable_field_try

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
use entsor
use numvar

!===============================================================================

implicit none

procedure() :: csexit

! Arguments

! Local variables

if (nvar .gt. nvarmx) then
  write(nfecra,1000) nvar, nvarmx
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
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined                           ',/,&
'@    corresponds to a number of variables NVAR  >= ', i10     ,/,&
'@  The maximum number of variables allowed                   ',/,&
'@                      in   paramx   is  NVARMX  = ', i10     ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@  If NVARMX is increased, the code must be reinstalled.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine fldvar_check_nvar

!===============================================================================

!> \brief Initialize the given variable calculation option structure with
!>        legacy values (iniini) allowing to later test user modification.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine init_var_cal_opt &
 ( id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use cs_c_bindings
use field
use optcal

!===============================================================================

implicit none

! Arguments
integer id

! Local variables
type(var_cal_opt) :: vcopt

! Most values set by default at in _var_cal_opt default;
! see cs_parameters.c

call field_get_key_struct_var_cal_opt(id, vcopt)

! Undefined values, may be modified by cs_parameters_*_complete
vcopt%isstpc = -999
vcopt%nswrsm = -1
vcopt%thetav = -1.d0
vcopt%blencv = -1.d0
vcopt%epsilo = -1.d0
vcopt%epsrsm = -1.d0
vcopt%relaxv = -1.d0

call field_set_key_struct_var_cal_opt(id, vcopt)

return

end subroutine init_var_cal_opt

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
use entsor
use numvar

!===============================================================================

implicit none

procedure() :: csexit

! Arguments

! Local variables

if ((nscaus+nscapp).gt.nscamx) then
  write(nfecra,1000) nscaus,nscamx,nscamx-nscaus
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
'@    allowed    in   paramx.h         est  NSCAMX = ', i10    ,/,&
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

  use, intrinsic :: iso_c_binding
  use cs_c_bindings

  implicit none

  procedure() :: add_variable_field_indexes

  ! Arguments

  integer(c_int), value :: f_id
  integer(c_int) :: ivar

  ! Local variables

  integer f_id0, ivar0

  f_id0 = f_id

  call add_variable_field_indexes(f_id0, ivar0)

  ivar = ivar0

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

!-------------------------------------------------------------------------------
!> \brief add field indexes associated with a new solved thermal variable,
!>        with default options
!
!> \param[in]  f_id    field id
!-------------------------------------------------------------------------------

subroutine cs_c_add_model_thermal_field_indexes(f_id) &
  bind(C, name='cs_add_model_thermal_field_indexes')

  use, intrinsic :: iso_c_binding
  use optcal
  use cs_c_bindings

  implicit none

  procedure() :: add_model_field_indexes

  ! Arguments

  integer(c_int), value :: f_id

  ! Local variables

  integer f_id0, iscal0

  f_id0 = f_id

  call add_model_field_indexes(f_id0, iscal0)

  iscalt = iscal0

end subroutine cs_c_add_model_thermal_field_indexes

!---------------------------------------------------------------------------
