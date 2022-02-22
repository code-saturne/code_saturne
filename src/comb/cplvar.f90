!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

subroutine cplvar
!================


!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      INIT DES POSITIONS DES VARIABLES TRANSPORTEES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use field

!===============================================================================

implicit none

integer           :: icha, isc, f_id
integer           :: kscmin, kscmax
character(len=80) :: f_label, f_name

!===============================================================================

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. Define variable fields
!===============================================================================

! Thermal model

itherm = 2
call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
iscalt = ihm

! Set min and max clipping
f_id = ivarfl(isca(iscalt))
call field_set_key_double(f_id, kscmin, -grand)
call field_set_key_double(f_id, kscmax, grand)

! Variables specific to continuous phase

do icha = 1, ncharb

  write(f_name,'(a7,i2.2)') 'mv1_fraction_', icha
  write(f_label,'(a6,i2.2)') 'Fr_mv1_', icha
  call add_model_scalar_field(f_name, f_label, if1m(icha))
  f_id = ivarfl(isca(if1m(icha)))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

enddo

do icha = 1, ncharb

  write(f_name,'(a7,i2.2)') 'mv2_fraction_', icha
  write(f_label,'(a6,i2.2)') 'Fr_mv2_', icha
  call add_model_scalar_field(f_name, f_label, if2m(icha))
  f_id = ivarfl(isca(if2m(icha)))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

enddo

call add_model_scalar_field('het_fraction', 'Fr_HET', if3m)
f_id = ivarfl(isca(if3m))

call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 1.d0)

call add_model_scalar_field('air_variance', 'Var_AIR', if4p2m)

f_id = ivarfl(isca(if4p2m))

call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 0.25d0)

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!    - PROPRES AUX SCALAIRES   : diffusivity_id
!      Rq : pas de variance associee a un scalaire dans notre cas
!    - PROPRES A LA SUSPENSION : icp
!===============================================================================

do isc = 1, nscapp

  if (iscavr(iscapp(isc)).le.0) then
    ! Reference dynamic viscosity relative to this scalar
    call field_set_key_int(ivarfl(isca(iscapp(isc))), kivisl, -1)
  endif

enddo

! Although we are in enthalpy formulation, we keep Cp constant

icp = -1

return
end subroutine
