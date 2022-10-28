!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Function:
! ---------

!> \file cs_user_boundary_conditions-atmospheric.f90
!>
!> Atmospheric example of cs_user_boundary_conditions.f90 subroutine
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[in,out] izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

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
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atchem
use atincl
use atsoil
use ctincl
use cs_fuel_incl
use mesh
use field
use turbomachinery
use iso_c_binding
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

!< [loc_var_dec]
integer          ifac, ii
integer          izone
integer          ilelt, nlelt
integer          f_id_rough, f_id_t_rough
double precision d2s3
double precision zref, xuref
double precision ustar, rugd, rugt
double precision zent, xuent, xvent
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: bpro_roughness
double precision, dimension(:), pointer :: bpro_roughness_t
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0

! Paremeters for the analytical rough wall law (neutral)
zref = 10.d0
xuref = 10.d0
rugd = 0.1d0
rugt = rugd

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! --- For boundary faces of color 11,
!       assign an inlet boundary condition for all phases prescribed from the
!       meteo profile with automatic choice between inlet/ outlet according to
!       the meteo profile

!< [example_1]
call getfbr('11',nlelt,lstelt)

! Get a new zone number (1 <= izone <= nozppm)
izone = maxval(izfppp) + 1

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the face belongs
  izfppp(ifac) = izone

  ! - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

  ! - Chemical boundary conditions are prescribed from the chemistry profile
  iprofc(izone) = 1

  ! - boundary condition type can be set to ientre or i_convective_inlet

  itypfb(ifac) = ientre

  ! - automatic determination of type (inlet/outlet) according to sign of
  !   mass flux

  iautom(ifac) = 1

enddo
!< [example_1]


! ---For boundary faces of color 21,
!     assign an inlet boundary condition for all phases prescribed from the
!     meteo profile

!< [example_2]
call getfbr('21',nlelt,lstelt)

! Get a new zone number (1 <= izone <= nozppm)
izone = maxval(izfppp) + 1

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the face belongs
  izfppp(ifac) = izone

  ! - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

  ! - Chemical boundary conditions are prescribed from the chemistry profile
  iprofc(izone) = 1

  ! - Assign inlet boundary conditions
  itypfb(ifac) = ientre

enddo
!< [example_2]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
