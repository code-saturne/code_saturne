!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
! Function:
! ---------

!> \file cs_user_boundary_conditions_ale.f90
!>
!> \brief User subroutine dedicated the use of ALE (Arbitrary Lagrangian
!> Eulerian) Method:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     itrale        number of iterations for ALE method
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!> \param[in,out] itypfb        boundary face types
!> \param[out]    ialtyb        boundary face types for mesh velocity
!> \param[in]     impale        indicator for fixed node displacement
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!> \param[in]     xyzno0        vertex coordinates of initial mesh
!> \param[in,out] disale        nodes total displacement (from time 0 to n)
!_______________________________________________________________________________

subroutine usalcl &
 ( itrale , itypfb , ialtyb , impale ,                                         &
   dt     , xyzno0 , disale )

!===============================================================================

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
use mesh
use field
use dimens, only: nvar, nscal
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itrale

integer          itypfb(nfabor), ialtyb(nfabor)
integer          impale(nnod)

double precision dt(ncelet)
double precision disale(3,nnod), xyzno0(3,nnod)

! Local variables

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================

call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

!--------
! Formats
!--------

!----
! End
!----


return
end subroutine usalcl
