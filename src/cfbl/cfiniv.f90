!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file cfiniv.f90
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
use field
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

double precision, dimension(:), pointer :: cpro_cp, cpro_cv, mix_mol_mas

!===============================================================================

!===============================================================================
! User initialisation of the variables only if this not a resumed computation
! and if it has not already been done in gas mix initialization
!===============================================================================

if (isuite.eq.0.and.ippmod(igmix).lt.0) then
  call cs_user_f_initialization(nvar, nscal, dt)
endif

!===============================================================================
! Computation of variable Cv in order to have a correct initialization
! of the total energy (computed in inivar by a call to a thermodynamic
! function), now that initial gas mixture composition is known.
! Note that the only eos with a variable Cv is the ideal gas mix (ieos=3).
!===============================================================================

if (icv.ge.0) then
  call field_get_val_s(icp, cpro_cp)
  call field_get_val_s(icv, cpro_cv)
  call field_get_val_s(igmxml, mix_mol_mas)

  call cs_cf_thermo_cv(cpro_cp, mix_mol_mas, cpro_cv, ncel)
endif ! Constant cv computed in iniusi

!----
! FORMATS
!----


!----
! END
!----

return
end subroutine
