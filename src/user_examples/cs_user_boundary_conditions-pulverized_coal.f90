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
!> \file  cs_user_boundary_conditions-pulverized_coal.f90
!> \brief Example of cs_f_user_boundary_conditions subroutine.f90
!>        for pulverized coal
!
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
use atincl
use ctincl
use cs_fuel_incl
use mesh
use field

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
integer          icha, iclapc
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! Air inlet with pulverized coal, e.g. secondary or tertiary air.

!< [example_1]
call getfbr('inlet', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Kind of boundary conditions for standard variables
  itypfb(ifac) = ientre

  ! zone's number (from 1 to n)
  izone = 1

  ! Assign zone number to the face
  izfppp(ifac) = izone

  ! For theses inlet faces, mass flux is fixed
  !  The speed direction is given
  !  (speed vector norm is irrelevant)

  ientat(izone) = 1
  iqimp(izone)  = 1

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

  ! Oxidizer's number (1 to 3)
  inmoxy(izone) = 1
  ! Oxidizer mass flow rate in kg/s
  qimpat(izone) = 1.46d-03
  ! Oxidizer Temperature in K
  timpat(izone) = 400.d0 + tkelvi

  ! Turbulence variables BC's calculated based on icalke(izone):
  ! - If 1: hydraulic diameter and reference velocity (similar to
  !         turbulence_bc_inlet_hyd_diam)
  ! - If 2: hydraulic diameter, reference velocity and turbulence intensity
  !         (similar to turbulence_bc_inlet_turb_intensity)
  ! - If 0: must be defined here for all turbulence variables
  !        (k, epsilon, omega, Rij, ... depending on active model)

  icalke(izone) = 1

  dh(izone)     = 0.032d0
  xintur(izone) = 0.d0

  ! Automatic treatment of non-user scalars

  ! Treatment of user-defined scalars

  if ((nscal-nscapp).gt.0) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1]

! Primary air inlet with pulverized coal.

!< [example_2]
call getfbr('primary_inlet', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Kind of boundary conditions for standard variables
  itypfb(ifac) = ientre

  ! zone's number (from 1 to n)
  izone = 2

  ! Assign zone number to the face
  izfppp(ifac) = izone

  ! For this inlet, mass flux is fixed
  !  The speed direction is given
  !  (speed vector norm is irrelevant)

  ientcp(izone) = 1
  iqimp(izone)  = 1

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

  ! Oxidizer's number (1 to 3)
  inmoxy(izone) = 1
  ! Oxidizer's mass flow rate in kg/s
  qimpat(izone) = 1.46d-03
  ! Oxidizer's Temperature in K
  timpat(izone) = 800.d0  + tkelvi

  !  Coal inlet, initialization
  do icha = 1, ncharm
    qimpcp(izone,icha) = zero
    timpcp(izone,icha) = zero
    do iclapc = 1, ncpcmx
      distch(izone,icha,iclapc) = zero
    enddo
  enddo

  ! code_saturne deals with NCHA different coals (component of blend)
  !       every coal is described by NCLPCH(icha) class of particles
  !       (each of them described by an inlet diameter)

  ! Treatment for the first coal
  icha = 1
  ! Coal mass flow rate in kg/s
  qimpcp(izone,icha) = 1.46d-4
  ! Percentage mass fraction of each granulometric class
  do iclapc = 1, nclpch(icha)
    distch(izone,icha,iclapc) = 100.d0/dble(nclpch(icha))
  enddo
  ! Inlet temperature for coal & primary air
  timpcp(izone,icha) = 800.d0 + tkelvi

  ! Turbulence variables BC's calculated based on icalke(izone):
  ! - If 1: hydraulic diameter and reference velocity (similar to
  !         turbulence_bc_inlet_hyd_diam)
  ! - If 2: hydraulic diameter, reference velocity and turbulence intensity
  !         (similar to turbulence_bc_inlet_turb_intensity)
  ! - If 0: must be defined here for all turbulence variables
  !        (k, epsilon, omega, Rij, ... depending on active model)

  icalke(izone) = 1

  dh(izone)     = 0.1d0
  xintur(izone) = 0.1d0

enddo
!< [example_2]

! Air inlet with pulverized coal with Lagragian model,
! e.g. secondary or tertiary air.

!< [example_3]
call getfbr('inlet', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Kind of boundary conditions for standard variables
  itypfb(ifac) = ientre

  ! zone's number (from 1 to n)
  izone = 1

  ! Assign zone number to the face
  izfppp(ifac) = izone

  ! For this inlet, mass flux is fixed
  !  The speed direction is given
  !  (speed vector norm is irrelevant)

  ientat(izone) = 1
  iqimp(izone)  = 1

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

  ! Oxidizer mass flow rate in kg/s
  qimpat(izone) = 1.46d-03
  ! Oxidizer Temperature in K
  timpat(izone) = 400.d0 + tkelvi

  ! Turbulence variables BC's calculated based on icalke(izone):
  ! - If 1: hydraulic diameter and reference velocity (similar to
  !         turbulence_bc_inlet_hyd_diam)
  ! - If 2: hydraulic diameter, reference velocity and turbulence intensity
  !         (similar to turbulence_bc_inlet_turb_intensity)
  ! - If 0: must be defined here for all turbulence variables
  !        (k, epsilon, omega, Rij, ... depending on active model)

  icalke(izone) = 1

  dh(izone)     = 0.1d0
  xintur(izone) = 0.1d0

  ! Treatment of user-defined scalars

  if ((nscal-nscapp).gt.0) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_3]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
