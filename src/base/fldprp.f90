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

!> \file fldprp.f90
!> \brief Properties definition initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine fldprp () &
 bind(C, name='cs_f_fldprp')

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
use atincl
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

! Local variables

integer :: iflid

!===============================================================================
! 1. PROPRIETES PRINCIPALES
!===============================================================================

! Base properties, always present

call field_get_id_try('density', icrom)
call field_get_id_try('boundary_density', ibrom)

call field_get_id_try('molecular_viscosity', iviscl)
call field_get_id_try('turbulent_viscosity', ivisct)

! Set itemp if temperature is present as a property

if (itemp .eq. 0) then
  if (itherm.eq.2 .or. itherm.eq.4) then
    call field_get_id_try('temperature', iflid)
    if (iflid.ge.0) itemp = iflid
  endif
endif

! ---> Atmospheric modules:
if (ippmod(iatmos).ge.0) then
  if (iatmst.ge.1) then
    call field_get_id('momentum_source_terms', imomst)
  endif
  call field_get_id_try("real_temperature", itempc)
  if (ippmod(iatmos).eq.2) then
    call field_get_id("liquid_water", iliqwt)
  end if
endif

!====
! End
!====

return
end subroutine fldprp

!===============================================================================

!> \brief add field defining a one-dimensional property field defined on cells,
!>        with no previous time values and with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     label         field default label, or empty
!> \param[out]    f_id          field id
!_______________________________________________________________________________

subroutine add_property_field_1d &
 ( name, label, f_id )

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

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(out)         :: f_id

! Local variables

character(len=len_trim(name)+1, kind=c_char) :: c_name
integer(c_int) :: c_type_flag
integer(c_int) :: c_location_id
integer(c_int) :: c_dim
logical(c_bool) :: c_has_previous

!===============================================================================

! Test if the field has already been defined
call field_get_id_try(trim(name), f_id)
if (f_id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit(1)
endif

! Create field
c_name = trim(name)//c_null_char
c_type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
c_location_id = 1
c_dim = 1
c_has_previous = .false.

call cs_physical_property_define_from_field(c_name, c_type_flag, &
  c_location_id, c_dim, c_has_previous)

f_id = cs_physical_property_field_id_by_name(c_name)

call field_set_key_int(f_id, keylog, 1)
call field_set_key_int(f_id, keyvis, 0)

if (len(trim(label)).gt.0) then
  call field_set_key_str(f_id, keylbl, trim(label))
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP'              ,/,&
'@    ======'                                                  ,/,&
'@     FIELD: ', a, ' HAS ALREADY BEEN DEFINED.'               ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

end subroutine add_property_field_1d

!===============================================================================

!> \brief disable logging and postprocessing for a property field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          field id
!_______________________________________________________________________________

subroutine hide_property &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use entsor
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

!===============================================================================

call field_set_key_int(f_id, keyvis, 0)
call field_set_key_int(f_id, keylog, 0)

return

end subroutine hide_property

!===============================================================================
