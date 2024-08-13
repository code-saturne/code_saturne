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
! Purpose:
! --------

!> \file ppprop.f90
!>
!> \brief Add additional property fields for dedicated modules
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________


subroutine ppprop () &
 bind(C, name='cs_f_ppprop')

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
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

!===============================================================================

! ---> Physique particuliere : Combustion Gaz

if (ippmod(icod3p).ge.0 .or. ippmod(islfm ).ge.0 .or.            &
    ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0) then
  call coprop
endif

! ---> Physique particuliere :  Combustion Charbon Pulverise

if (ippmod(iccoal).ge.0) then
  call cs_coal_prop
endif

end subroutine ppprop
