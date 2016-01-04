!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!> \file cfxtcl.f90
!> \brief Initialisation of the variables if the compressible flow model is
!> enabled.
!>
!> This subroutine is called at the beginning of a computation (or when a
!> computation is resumed) before the start of the time loop.
!>
!> It allows to initialise or modify (for resumed computations) the variables
!> and the time step values.
!>
!> Before this subroutine call, the density and the molecular viscosity have
!> been initialised at ro0 and viscl0 respectively or they have been read in
!> a checkpoint file in the case of a resumed computation.
!> If the scalar diffusivities (visls) and the isobaric specific heat (cp) were
!> defined (i.e. variable), their values are here at hand only if a computation
!> is resumed.
!>
!> Any modification of a physical property (density, molecular viscosity,
!> scalar diffusivity, isobaric specific heat) shall be performed in the ppphyv
!> subroutine and never here.
!>
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!______________________________________________________________________________

subroutine cfiniv &
 ( nvar   , nscal  , dt     )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use cs_cf_bindings

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!===============================================================================

!===============================================================================
! Initialisation of the variables only if this not a resumed computation
!===============================================================================

if (isuite.eq.0) then
  call cs_user_initialization(nvar, nscal, dt)
else
  ! Initialisation of the isochoric specific heat
  call cs_cf_thermo_default_init(isuite)
endif

!----
! FORMATS
!----


!----
! END
!----

return
end subroutine
