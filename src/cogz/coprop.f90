!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine coprop
!================

!===============================================================================
! Purpose:
! --------

! Define state variables for gas combustion,
! diffusion and premixed flame.

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
use radiat

!===============================================================================

implicit none

! Local variables

integer       idirac

character(len=80) :: f_name, f_label

!===============================================================================

! Flamme de diffusion chimie 3 points
!====================================

if (ippmod(icod3p).ge.0) then

  call add_property_field_1d('temperature', 'Temperature', itemp)

  call add_property_field_1d('ym_fuel', 'Ym_Fuel', iym(1))
  call add_property_field_1d('ym_oxyd', 'Ym_Oxyd', iym(2))
  call add_property_field_1d('ym_prod', 'Ym_Prod', iym(3))

endif

! Premixed flame - EBU model
!===========================

if (ippmod(icoebu).ge.0) then

  call add_property_field_1d('temperature', 'Temperature', itemp)

  call add_property_field_1d('ym_fuel', 'Ym_Fuel', iym(1))
  call add_property_field_1d('ym_oxyd', 'Ym_Oxyd', iym(2))
  call add_property_field_1d('ym_prod', 'Ym_Prod', iym(3))

endif

! Premixed flame - LWC model
!===========================

if (ippmod(icolwc).ge.0) then

  call add_property_field_1d('temperature', 'Temperature', itemp)
  call add_property_field_1d('molar_mass',  'Molar_Mass',  imam)
  call add_property_field_1d('source_term', 'Source_Term', itsc)

  call add_property_field_1d('ym_fuel', 'Ym_Fuel', iym(1))
  call add_property_field_1d('ym_oxyd', 'Ym_Oxyd', iym(2))
  call add_property_field_1d('ym_prod', 'Ym_Prod', iym(3))

  do idirac = 1, ndirac

    write(f_name,  '(a,i1)') 'rho_local_', idirac
    write(f_label, '(a,i1)') 'Rho_Local_', idirac
    call add_property_field_1d(f_name, f_label, irhol(idirac))

    write(f_name,  '(a,i1)') 'temperature_local_', idirac
    write(f_label, '(a,i1)') 'Temperature_Local_', idirac
    call add_property_field_1d(f_name, f_label, iteml(idirac))

    write(f_name,  '(a,i1)') 'ym_local_', idirac
    write(f_label, '(a,i1)') 'Ym_Local_', idirac
    call add_property_field_1d(f_name, f_label, ifmel(idirac))

    write(f_name,  '(a,i1)') 'w_local_', idirac
    write(f_label, '(a,i1)') 'w_Local_', idirac
    call add_property_field_1d(f_name, f_label, ifmal(idirac))

    write(f_name,  '(a,i1)') 'amplitude_local_', idirac
    write(f_label, '(a,i1)') 'Amplitude_Local_', idirac
    call add_property_field_1d(f_name, f_label, iampl(idirac))

    write(f_name,  '(a,i1)') 'chemical_st_local_', idirac
    write(f_label, '(a,i1)') 'Chemical_ST_Local_', idirac
    call add_property_field_1d(f_name, f_label, itscl(idirac))

    write(f_name,  '(a,i1)') 'molar_mass_local_', idirac
    write(f_label, '(a,i1)') 'M_Local_', idirac
    call add_property_field_1d(f_name, f_label, imaml(idirac))

  enddo

endif

! Additional fields for radiation
!================================

if (iirayo.ge.1) then
  if (ippmod(icod3p).eq.1 .or.                                  &
      ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 .or.         &
      ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.         &
      ippmod(icolwc).eq.5) then
    call add_property_field_1d('kabs',          'KABS',  ickabs)
    call add_property_field_1d('temperature_4', 'Temp4', it4m)
    call add_property_field_1d('temperature_3', 'Temp3', it3m)
  endif
endif

return

end subroutine
