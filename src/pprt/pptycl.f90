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

!> \file pptycl.f90
!>
!> \brief Boundary conditions for specific physics modules.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     init          partial treatment (before time loop) if true
!> \param[in,out] itypfb        boundary face types
!_______________________________________________________________________________

subroutine pptycl           &
 ( init , itypfb )  &
 bind(C, name='cs_f_pptycl')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! Arguments

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use coincl
use ppincl
use atincl
use mesh
use field
use cs_c_bindings
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

logical(c_bool), value :: init

integer(c_int) ::  itypfb(nfabor)

! Local variables

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_coal_boundary_conditions(bc_type)  &
    bind(C, name='cs_coal_boundary_conditions')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), dimension(*) :: bc_type
  end subroutine cs_coal_boundary_conditions

  subroutine cs_combustion_boundary_conditions(bc_type)  &
    bind(C, name='cs_combustion_boundary_conditions')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), dimension(*) :: bc_type
  end subroutine cs_combustion_boundary_conditions

  subroutine cs_combustion_boundary_conditions_ebu(bc_type)  &
    bind(C, name='cs_combustion_boundary_conditions_ebu')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), dimension(*) :: bc_type
  end subroutine cs_combustion_boundary_conditions_ebu

  subroutine cs_combustion_boundary_conditions_lw(bc_type)  &
    bind(C, name='cs_combustion_boundary_conditions_lw')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), dimension(*) :: bc_type
  end subroutine cs_combustion_boundary_conditions_lw

  subroutine cs_ctwr_bcond()  &
    bind(C, name='cs_ctwr_bcond')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_ctwr_bcond

  subroutine cs_atmo_bcond()  &
    bind(C, name='cs_atmo_bcond')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_atmo_bcond

end interface

!===============================================================================
! Call to boundary conditions computations, model by model.
! Computations should not be called under initialization for most
! models; those for which it should be called are tested first.
!===============================================================================

! Atmospheric flows
if (ippmod(iatmos).ge.0) then
  call cs_atmo_bcond()
endif

! Cooling towers
if (ippmod(iaeros).ge.0) then
  call cs_ctwr_bcond()
endif

if (init .eqv. .true.) return

! 3-point chemistry or steady laminar flamelet

if (ippmod(icod3p).ge.0 .or. ippmod(islfm).ge.0) then
  call cs_combustion_boundary_conditions(itypfb)

! ---> Combustion gaz USEBUC
!      Flamme de premelange modele EBU

elseif (ippmod(icoebu).ge.0) then
  call cs_combustion_boundary_conditions_ebu(itypfb)

! ---> Combustion gaz USLWCC
!      Flamme de premelange modele LWC

elseif (ippmod(icolwc).ge.0) then
  call cs_combustion_boundary_conditions_lw(itypfb)

! ---> Combustion charbon pulverise

elseif (ippmod(iccoal).ge.0) then
  call cs_coal_boundary_conditions(itypfb)

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
