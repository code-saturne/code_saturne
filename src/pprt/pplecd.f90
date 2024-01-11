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
!  Function:
!  ---------

!> file pplecd.f90
!>
!> \brief Read specific physical model data file

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!

subroutine pplecd

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_gui_coal_model()  &
    bind(C, name='cs_gui_coal_model')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_coal_model

  subroutine cs_coal_read_data()  &
    bind(C, name='cs_coal_read_data')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_coal_read_data

end interface

!===============================================================================

! ---> Diffusion flame - 3-point chemistry
!      Premix flame    - EBU model
!      Premix flame    - LWC model

if (ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                      &
                        .or. ippmod(icolwc).ge.0) then
  call colecd
endif

! ---> Diffusion flame - Steady laminar flamelet approach

if (ippmod(islfm).ge.0) then
  call cs_steady_laminar_flamelet_read_base
endif

! ---> Pulverized coal combustion

if (ippmod(iccoal).ge.0) then
  call cs_gui_coal_model
  call cs_coal_read_data
endif

! ---> Joule effect, electric arc, or ionic conduction

if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
  call ellecd(ippmod(ieljou), ippmod(ielarc))
endif

!----
! End
!----

return
end subroutine

