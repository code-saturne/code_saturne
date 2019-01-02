!-------------------------------------------------------------------------------

!VERS

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
! Function:
! ---------

! Basic example of cs_user_boundary_conditions subroutine.f90
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
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
use ihmpre
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
integer          ifac, iel, ii
integer          ilelt, nlelt
double precision uref2
double precision rhomoy, xdh
double precision xitur

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: bfpro_rom
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

call field_get_val_s(ibrom, bfpro_rom)
!< [init]

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! Assign an inlet to boundary faces of group '2' and x < 0.01,

!< [example_1]
call getfbr('2 and x < 0.01', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel = ifabor(ifac)

  itypfb(ifac) = ientre

  rcodcl(ifac,iu,1) = 1.1d0
  rcodcl(ifac,iv,1) = 1.1d0
  rcodcl(ifac,iw,1) = 1.1d0

  uref2 = rcodcl(ifac,iu,1)**2  &
        + rcodcl(ifac,iv,1)**2  &
        + rcodcl(ifac,iw,1)**2
  uref2 = max(uref2,1.d-12)

  !   Turbulence example computed using equations valid for a pipe.

  !   We will be careful to specify a hydraulic diameter adapted
  !     to the current inlet.

  !   We will also be careful if necessary to use a more precise
  !     formula for the dynamic viscosity use in the calculation of
  !     the Reynolds number (especially if it is variable, it may be
  !     useful to take the law from 'usphyv'. Here, we use by default
  !     the 'viscl0" value.
  !   Regarding the density, we have access to its value at boundary
  !     faces (romb) so this value is the one used here (specifically,
  !     it is consistent with the processing in 'usphyv', in case of
  !     variable density)

  !     Hydraulic diameter
  xdh     = 0.075d0

  !   Calculation of turbulent inlet conditions using
  !     standard laws for a circular pipe
  !     (their initialization is not needed here but is good practice).
  rhomoy  = bfpro_rom(ifac)

  call turbulence_bc_inlet_hyd_diam(ifac, uref2, xdh, rhomoy, viscl0,  &
                                    rcodcl)

  ! Handle scalars
  if (nscal.gt.0) then
    do ii = 1, nscal
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1]

! Assign an inlet to boundary faces of group '3'

!< [example_2]
call getfbr('3', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel  = ifabor(ifac)

  itypfb(ifac) = ientre

  rcodcl(ifac,iu,1) = 1.1d0
  rcodcl(ifac,iv,1) = 1.1d0
  rcodcl(ifac,iw,1) = 1.1d0

  uref2 = rcodcl(ifac,iu,1)**2   &
        + rcodcl(ifac,iv,1)**2   &
        + rcodcl(ifac,iw,1)**2
  uref2 = max(uref2,1.d-12)

  ! Turbulence example computed using turbulence intensity data.

  ! We will be careful to specify a hydraulic diameter adapted
  !   to the current inlet.

  ! Hydraulic diameter

  xdh   = 0.075d0
  ! Turbulence intensity
  xitur = 0.02d0

  ! Calculation of turbulent inlet conditions using
  !   the turbulence intensity and standard laws for a circular pipe
  !   (their initialization is not needed here but is good practice)

  call turbulence_bc_inlet_turb_intensity(ifac, uref2, xitur, xdh,  &
                                          rcodcl)

  ! --- Handle scalars
  if (nscal.gt.0) then
    do ii = 1, nscal
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_2]

! Assign an outlet to boundary faces of group 'outlet'

!< [example_3]
call getfbr('outlet', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Outlet: zero flux for velocity and temperature, prescribed pressure
  !         Note that the pressure will be set to P0 at the first
  !         free outlet face (isolib)

  itypfb(ifac) = isolib

enddo
!< [example_3]

! Assign a wall to boundary faces of group '5'

!< [example_4]
call getfbr('5', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  itypfb(ifac) = iparoi

  ! If sliding wall with velocity u = 1:
  ! rcodcl(ifac, iu, 1) = 1.d0

  ! If sliding wall with velocity u = 0: nothing to do

  if (nscal.gt.0) then

    ! If temperature prescribed to 20 with wall law (scalar ii=1):
    ii = 1
    icodcl(ifac, isca(ii))    = 5
    rcodcl(ifac, isca(ii), 1) = 20.d0

    ! If temperature prescribed to 50 with no wall law (simple Dirichlet)
    !   with exchange coefficient 8 (scalar ii=2):
    ii = 2
    icodcl(ifac, isca(ii))    = 1
    rcodcl(ifac, isca(ii),1)  = 50.d0
    rcodcl(ifac, isca(ii), 2) = 8.d0

    ! If flux prescribed to 4.d0 (scalar ii=3):
    ii = 3
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 4.d0

  endif
enddo
!< [example_4]

! Assign a rough wall to boundary faces of group '7'

!< [example_5]
call getfbr('7', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       rough friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  itypfb(ifac) = iparug

  ! Roughness for velocity: 1cm
  rcodcl(ifac,iu,3) = 0.01d0

  ! Roughness for scalar (if required): 1cm
  rcodcl(ifac,iv,3) = 0.01d0

  ! If sliding wall with velocity u = 1:
  rcodcl(ifac, iu, 1) = 1.d0

  ! If sliding wall with velocity u = 0: nothing to do
  if (nscal.gt.0) then

    ! If temperature prescribed to 20 (scalar ii=1)
    ! (with thermal roughness specified in rcodcl(ifac,iv,3)) :
    ii = 1
    icodcl(ifac, isca(ii))    = 6
    rcodcl(ifac, isca(ii), 1) = 20.d0

    ! If flux prescribed to 4.d0 (scalar ii=3):
    ii = 3
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 4.d0

  endif

enddo
!< [example_5]

! Assign a symmetry to boundary faces of group '4'

!< [example_6]
call getfbr('4', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Symmetries

  itypfb(ifac) = isymet

enddo
!< [example_6]

!--------
! Formats
!--------

!----
! End
!----

!< [finalize]
deallocate(lstelt)  ! temporary array for boundary faces selection
!< [finalize]

return
end subroutine cs_f_user_boundary_conditions
