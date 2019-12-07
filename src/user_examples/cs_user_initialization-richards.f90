!-------------------------------------------------------------------------------
!
!VERS
!
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_initialization-richards.f90
!>
!> \brief Basic example
!>
!> See \subpage cs_user_initialization for examples.
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
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
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments
integer          :: nvar, nscal
double precision :: dt(ncelet)

! Local variables

!< [loc_var_dec]
integer :: iel, icelt, ncelt, isorb, keysrb, igwfpr, keypre
integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cvar_scal_1, cvar_scal_2
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cpro_sorb, cpro_precip
!< [loc_var_dec]

!===============================================================================
! Initialization
!===============================================================================

!< [alloc]
allocate(lstelt(ncel)) ! temporary array for cells selection
!< [alloc]

!===============================================================================
! Variables initialization for groundwater flow module:
!
! * isca(i) is the number related to the solute number i
! * cvar_scal_i(iel) is the value of this variable in cell number iel.
! * cvar_pr(iel) is the value of the hydraulic pressure H in cell number iel
! * cvar_vel(j,iel) is the value of the component j of the velocity in cell iel
!
!   ONLY done if there is no restart computation
!===============================================================================

if (isuite.eq.0) then
  call field_get_val_s(ivarfl(isca(1)), cvar_scal_1)
  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_v(ivarfl(iu), cvar_vel)

!< [richards_init_cell]
! Global initialisation
  do iel = 1, ncel

!   Initialisation of the hydraulic pressure H with a contant gradient of 1
!   among z axis and -0.01 among x axis
    cvar_pr(iel) = 1.d0*xyzcen(3,iel) - 1.d-2*xyzcen(1,iel)

!   Null concentration by default
    cvar_scal_1(iel) = 0.d0

!   Velocity initialisation for security reason
    cvar_vel(1,iel) = 0.d0
    cvar_vel(2,iel) = 0.d0
    cvar_vel(3,iel) = 0.d0
  enddo
!< [richards_init_cell]

!< [richards_init_grp]
! Initialisation per group of volume
  call getcel('SOURCE', ncelt, lstelt)
  do icelt = 1, ncelt
    iel = lstelt(icelt)
    cvar_scal_1(iel) = 1.d0
  enddo
!< [richards_init_grp]

!< [richards_init_sorb]
  ! Index field of sorbed concentration
  call field_get_key_id("gwf_sorbed_concentration_id", keysrb)
  call field_get_key_int(ivarfl(isca(1)), keysrb, isorb)
  call field_get_val_s(isorb, cpro_sorb)

  do iel = 1, ncel
    ! no initial contamination of sorbed phase
    cpro_sorb(iel) = 0.d0
  enddo
!< [richards_init_sorb]

!< [richards_init_precip]
  ! Index field of precipitated concentration
  call field_get_key_id("gwf_precip_concentration_id", keypre)
  call field_get_key_int(ivarfl(isca(1)), keypre, igwfpr)
  call field_get_val_s(igwfpr, cpro_precip)

  do iel = 1, ncel
    ! no initial precipitation phase
    cpro_precip(iel) = 0.d0
  enddo
!< [richards_init_precip]

endif

!< [finalize]
deallocate(lstelt)  ! temporary array for cells selection
!< [finalize]

return
end subroutine cs_user_f_initialization
