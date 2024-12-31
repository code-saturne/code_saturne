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

!===============================================================================
! Function :
! --------

!> \brief Convert enthalpy to temperature at cells for gas combustion.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     xespec       mass fraction of constituents
!> \param[in]     enthal       enthalpy at cells
!> \param[out]    temper       temperature at cells
!_______________________________________________________________________________


function cs_gas_combustion_h_to_t(xespec, enthal) result(temper)  &
   bind(C, name='cs_gas_combustion_h_to_t')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding
use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppthch
use coincl

!===============================================================================

implicit none

! Arguments

real(kind=c_double), dimension(*) :: xespec
real(c_double), value :: enthal
real(c_double) :: temper

! Local variables

integer          it , iesp

double precision eh1 , eh0

it  = npo-1
eh1 = zero
do iesp = 1, ngazg
  eh1 = eh1 + xespec(iesp)*ehgazg(iesp,it+1)
enddo
if (enthal.ge.eh1) temper = th(it+1)

it  = 1
eh0 = zero
do iesp = 1, ngazg
  eh0 = eh0 + xespec(iesp)*ehgazg(iesp,it  )
enddo
if (enthal.le.eh0) temper = th(it)

do it = 1, npo-1
  eh0 = zero
  eh1 = zero
  do iesp = 1, ngazg
    eh0 = eh0 + xespec(iesp)*ehgazg(iesp,it  )
    eh1 = eh1 + xespec(iesp)*ehgazg(iesp,it+1)
  enddo
  if (enthal.ge.eh0 .and. enthal.le.eh1)                        &
    temper = th(it) + (enthal-eh0)*(th(it+1)-th(it))/(eh1-eh0)
enddo

!----
! End
!----

return
end function cs_gas_combustion_h_to_t

!===============================================================================
! Function :
! --------

!> \brief Convert a temperature value to enthalpy for gas combustion.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     xespec       mass fraction of constituants
!> \param[in]     temper       temperature at cells
!> \param[out]    enthal       enthalpy at cells
!_______________________________________________________________________________


function cs_gas_combustion_t_to_h(xespec, temper) result(enthal)  &
   bind(C, name='cs_gas_combustion_t_to_h')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding
use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppthch
use coincl

!===============================================================================

implicit none

! Arguments

real(kind=c_double), dimension(*) :: xespec
real(c_double), value :: temper
real(c_double) :: enthal

! Local variables

integer          it , iesp

double precision eh1 , eh0

!===============================================================================
! 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
!===============================================================================

it = npo
if (temper.ge.th(it)) then
  enthal = zero
  do iesp = 1, ngazg
    enthal = enthal + xespec(iesp)*ehgazg(iesp,it)
  enddo
  go to 11
endif

it = 1
if (temper.le.th(it)) then
  enthal = zero
  do iesp = 1, ngazg
    enthal = enthal + xespec(iesp)*ehgazg(iesp,it)
  enddo
  go to 11
endif
10 continue

it = it + 1
if (temper.le.th(it)) then
  eh0 = zero
  eh1 = zero
  do iesp = 1, ngazg
    eh0 = eh0 + xespec(iesp)*ehgazg(iesp,it-1)
    eh1 = eh1 + xespec(iesp)*ehgazg(iesp,it  )
  enddo
  enthal = eh0 + (eh1-eh0)*(temper-th(it-1))/(th(it)-th(it-1))
  goto 11
endif
goto 10
11 continue

!----
! End
!----

return
end function cs_gas_combustion_t_to_h
