!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
!>                               - 6 rought wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!> \param[in]                    (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughtness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________


subroutine cs_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce ,                            &
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
use elincl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

!< [loc_var_dec]
integer          ifac, iel, ii, ivar
integer          izone
integer          ilelt, nlelt
double precision uref2, d2s3
double precision rhomoy, xdh, xustar2
double precision xitur
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: brom
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0

call field_get_val_s(ibrom, brom)
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

  !   Calculation of friction velocity squared (ustar2)
  !     and of k and epsilon at the inlet (xkent and xeent) using
  !     standard laws for a circular pipe
  !     (their initialization is not needed here but is good practice).
  rhomoy  = brom(ifac)
  xustar2 = 0.d0
  xkent   = epzero
  xeent   = epzero

  call keendb &
  !==========
( uref2, xdh, rhomoy, viscl0, cmu, xkappa,   &
  xustar2, xkent, xeent )

  ! itytur is a flag equal to iturb/10
  if    (itytur.eq.2) then

    rcodcl(ifac,ik,1)  = xkent
    rcodcl(ifac,iep,1) = xeent

  elseif (itytur.eq.3) then

    rcodcl(ifac,ir11,1) = d2s3*xkent
    rcodcl(ifac,ir22,1) = d2s3*xkent
    rcodcl(ifac,ir33,1) = d2s3*xkent
    rcodcl(ifac,ir12,1) = 0.d0
    rcodcl(ifac,ir13,1) = 0.d0
    rcodcl(ifac,ir23,1) = 0.d0
    rcodcl(ifac,iep,1)  = xeent
    if (iturb.eq.32) then
      rcodcl(ifac,ial,1)  = 1.d0
    endif

  elseif (itytur.eq.5) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iep,1)  = xeent
    rcodcl(ifac,iphi,1) = d2s3
    if (iturb.eq.50) then
      rcodcl(ifac,ifb,1)  = 0.d0
    elseif (iturb.eq.51) then
      rcodcl(ifac,ial,1)  = 0.d0
    endif

  elseif (iturb.eq.60) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iomg,1) = xeent/cmu/xkent

  elseif (iturb.eq.70) then

    rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

  endif

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

  ! Calculation of k and epsilon at the inlet (xkent and xeent) using
  !   the turbulence intensity and standard laws for a circular pipe
  !   (their initialization is not needed here but is good practice)
  xkent  = epzero
  xeent  = epzero

  call keenin &
  !==========
( uref2, xitur, xdh, cmu, xkappa, xkent, xeent )

  ! itytur is a flag equal to iturb/10
  if    (itytur.eq.2) then

    rcodcl(ifac,ik,1)  = xkent
    rcodcl(ifac,iep,1) = xeent

  elseif (itytur.eq.3) then

    rcodcl(ifac,ir11,1) = d2s3*xkent
    rcodcl(ifac,ir22,1) = d2s3*xkent
    rcodcl(ifac,ir33,1) = d2s3*xkent
    rcodcl(ifac,ir12,1) = 0.d0
    rcodcl(ifac,ir13,1) = 0.d0
    rcodcl(ifac,ir23,1) = 0.d0
    rcodcl(ifac,iep,1)  = xeent
    if (iturb.eq.32) then
      rcodcl(ifac,ial,1)  = 1.d0
    endif

  elseif (itytur.eq.5) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iep,1)  = xeent
    rcodcl(ifac,iphi,1) = d2s3
    if (iturb.eq.50) then
      rcodcl(ifac,ifb,1)  = 0.d0
    elseif (iturb.eq.51) then
      rcodcl(ifac,ial,1)  = 0.d0
    endif

  elseif (iturb.eq.60) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iomg,1) = xeent/cmu/xkent

  elseif (iturb.eq.70) then

    rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

  endif

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
end subroutine cs_user_boundary_conditions
