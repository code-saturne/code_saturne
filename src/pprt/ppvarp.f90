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

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine ppvarp() &
 bind(C, name='cs_f_ppvarp')

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use atchem
use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use coincl
use ppincl
use atincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          f_id

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_atmo_add_variable_fields() &
    bind(C, name='cs_atmo_add_variable_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_atmo_add_variable_fields

end interface

!===============================================================================

! 6. Atmospheric model
!---------------------

if (ippmod(iatmos).ge.0) then
  call cs_atmo_add_variable_fields()

  call init_chemistry_pointers()
  ! Update scalar id in Fortran no need in c version
  if (ippmod(iatmos).eq.2) then
    call field_get_id('ym_water', f_id)
    call field_get_key_int_by_name(f_id, "scalar_id", iymw)

    call field_get_id('number_of_droplets', f_id)
    call field_get_key_int_by_name(f_id, "scalar_id", intdrp)
  end if
endif

return
end subroutine
