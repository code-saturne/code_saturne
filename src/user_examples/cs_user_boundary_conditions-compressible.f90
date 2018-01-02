!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \file  cs_user_boundary_conditions-compressible.f90
!> \brief Example of cs_f_user_boundary_conditions.f90 for compressible
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
integer          ifac  , iel   , ii
integer          izone , iutile
integer          ilelt, nlelt

double precision uref2 , dhy  , rhomoy
double precision d2s3

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0

call field_get_val_s(icrom, cpro_rom)

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! --- Example of inlet/outlet for which everything is known

!       Without assuming the subsonic or supersonic nature of the inlet
!       the user wishes to impose all the characteristics of the flow,
!       a supersonic inlet is a particular case.

!       The turbulence and the user scalars take a zero flux if the
!       velocity is outward.

!< [example_1]
call getfbr('1 and X <= 1.0 ', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Cell adjacent to boundary face

  iel = ifabor(ifac)

  ! Number of zones from 1 to n...
  izone = 1
  izfppp(ifac) = izone

  itypfb(ifac) = iesicf

  ! Velocity
  rcodcl(ifac,iu,1) = 5.0d0
  rcodcl(ifac,iv,1) = 0.0d0
  rcodcl(ifac,iw,1) = 0.0d0

  ! Pressure, Density, Temperature, Total Specific Energy

  !     Only 2 of the 4 variables are independant,
  !     hence one can impose values for any couple of variables
  !     (except Temperature-Energy) and the two other variables
  !     will be computed automatically

  !  ** Choose a couple of variables which values are to be imposed
  !     and delete the others (that will be computed with the help of
  !     the thermodynamic laws in cfther.f90).

  ! Pressure (in Pa)
  rcodcl(ifac,ipr,1) = 5.d5

  ! Temperature (in K)
  rcodcl(ifac,isca(itempk),1) = 300.d0

  ! Specific Total Energy (in J/kg)
  ! rcodcl(ifac,isca(ienerg),1) = 355.d3

  ! Turbulence

  uref2 = rcodcl(ifac,iu,1)**2 + rcodcl(ifac,iv,1)**2 + rcodcl(ifac,iw,1)**2
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
  dhy   = 0.075d0

  rhomoy = cpro_rom(iel)

  call turbulence_bc_inlet_hyd_diam &
  ( ifac, uref2, dhy, rhomoy, viscl0, rcodcl )

  ! Handle scalars
  ! (do not loop on nscal to avoid modifying rho and energy)
  if(nscaus.gt.0) then
    do ii = 1, nscaus
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1]

! --- Supersonic outlet example

!     All the characteristics are outward,
!     nothing needs to be imposed (only internal values are used
!     to compute the boundary flux).

!     for the turbulence and the scalar, if values of rcodcl are
!     provided here, we impose them as Dirichlet if the mass flux is
!     inward ; otherwise a zero flux is imposed (outward mass flux or
!     RCODCL values given here).
!     Note that for turbulence, RCODCL has to be filled in for all
!     turbulent variable (otherwise a zero flux is imposed).

!< [example_2]
call getfbr('2', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Number zones from 1 to n...
  izone = 2
  izfppp(ifac) = izone

  itypfb(ifac) = isspcf

enddo
!< [example_2]

! --- Example of subsonic inlet (flow rate, entalpy flow rate)

!     2 characteristics out of 3 are inward : 2 informations have
!     to be given, the third is deduced by a 2-contact and
!     3-rarefaction scenario in the domain
!     here it is chosen to give (rho*(U.n), rho*(U.n)*H)
!     with H = 1/2 U*U + P/rho + e
!            n being the unit inward normal

!     WARNING, flux DENSITIES have to be given (per area unit)

!< [example_3]
call getfbr('3', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Number of zones from 1 to n...
  izone = 3
  izfppp(ifac) = izone

  itypfb(ifac) = ieqhcf

  ! - flow rate density (kg/(m2 s))
  rcodcl(ifac,irun,1) = 5.d5

  ! - enthalpy flow rate density (J/(m2 s))
  rcodcl(ifac,irunh,1) = 5.d5

  !   Unavailable B.C. in current version
  call csexit (1)
  !==========

enddo
!< [example_3]

! --- Example of subsonic inlet (total pressure, total enthalpy)

!< [example_4]
call getfbr('4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Number zones from 1 to n...
  izone = 4
  izfppp(ifac) = izone

  itypfb(ifac) = iephcf

  ! Total pressure (in Pa)
  rcodcl(ifac,ipr,1) = 1.d5

  ! Total enthalpy
  rcodcl(ifac,isca(ienerg),1) = 294465.d0

  ! Direction of the velocity: normal to inlet faces
  rcodcl(ifac,iu,1) = -surfbo(1,ifac)
  rcodcl(ifac,iv,1) = -surfbo(2,ifac)
  rcodcl(ifac,iw,1) = -surfbo(3,ifac)

  ! Turbulence (no turbulence)

  ! Handle scalars
  ! (do not loop on nscal to avoid modifying rho and energy)
  if(nscaus.gt.0) then
    do ii = 1, nscaus
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_4]

! --- Subsonic outlet example

! 1 characteristic out of 3 exits: 1 information must be given
! the 2 others are deduced by a 2-contact and 3-relaxation in the domain.
! Here we choose to definer P.

! Turbulence and user scalars take a zero flux.

!< [example_5]
call getfbr('5', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Number zones from 1 to n...
  izone = 5
  izfppp(ifac) = izone

  itypfb(ifac) = isopcf

  ! Pressure (in Pa)
  rcodcl(ifac,ipr,1) = 5.d5

enddo
!< [example_5]

! --- Wall example

!< [example_6]
call getfbr('7', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Number zones from 1 to n...
  izone = 7
  izfppp(ifac) = izone

  itypfb(ifac) = iparoi

  ! --- Sliding wall

  ! By default, the wall does not slide.
  ! If the wall slides, define the nonzero components of its velocity.
  ! The velocity will be projected in a plane tangent to the wall.
  ! In the following example, we prescribe Ux = 1.
  ! (example activated if iutile=1)

  iutile = 0
  if(iutile.eq.1) then
    rcodcl(ifac,iu,1) = 1.d0
  endif

  ! --- Prescribed temperature

  ! By default, the wall is adiabatic.
  ! If the wall has a prescribed temperature, indicate it by setting
  ! icodcl = 5 and define a value in Kelvin in rcodcl(., ., 1)
  ! In the following example, we prescribe T = 293.15 K
  ! (example activated if iutile=1)

  iutile = 0
  if(iutile.eq.1) then
    icodcl(ifac,isca(itempk))   = 5
    rcodcl(ifac,isca(itempk),1) = 20.d0 + 273.15d0
  endif

  ! --- Prescribed flux

  ! By default, the wall is adiabatic.
  ! If the wall has a prescribed flux, indicate it by setting
  ! icodcl = 3 and define the value in Watt/m2 in rcodcl(., ., 3)
  ! In the following example, we prescribe a flux of 1000 W/m2
  ! - a midday in the summer - (example is activated if iutile=1)

  iutile = 0
  if(iutile.eq.1) then
    icodcl(ifac,isca(itempk))   = 3
    rcodcl(ifac,isca(itempk),3) = 1000.d0
  endif

enddo
!< [example_6]

! --- Symmetry example

!< [example_7]
call getfbr('8', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Number zones from 1 to n...
  izone = 8
  izfppp(ifac) = izone

  itypfb(ifac) = isymet

enddo
!< [example_7]

! It is not recommended to use other boundary condition types than
! the ones provided above.

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
