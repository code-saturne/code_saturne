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
!> \file  cs_user_boundary_conditions-richards.f90
!> \brief Example of cs_f_user_boundary_conditions subroutine.f90 for
!>        groundwater flow module
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
integer          ifac, ilelt, nlelt, ii, iflmab
double precision Vel, Cref, R_act, S_act, T_act

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: bmasfl
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria on group (face name) to filter boundary faces of a
!   given subset
! - loop on faces from a subset
! - set the boundary condition for each face
!===============================================================================

! For groundwater flow module, the undefined type face iindef is always used.
! BC for the flow part are on the pressure head H=h+z

! Example 1: Dirichlet on hydraulic head and solutes

!< [richards_BC_ex1]
call getfbr('FACE1', nlelt, lstelt)

do ilelt = 1, nlelt

! Face number
  ifac = lstelt(ilelt)

! Undefined type face
  itypfb(ifac) = iindef

! Velocity is not solved but deduced from darcy law, then no BC are required.
! However, no flux BC are imposed for safety.
  icodcl(ifac,iu) = 3
  icodcl(ifac,iv) = 3
  icodcl(ifac,iw) = 3
  rcodcl(ifac,iu,3) = 0.d0
  rcodcl(ifac,iv,3) = 0.d0
  rcodcl(ifac,iw,3) = 0.d0

! Dirichlet BC on hydraulic head (H = h + z) to impose a constant value
  icodcl(ifac,ipr) = 1
  rcodcl(ifac,ipr,1) = 10.d0

! Dirichlet BC on centration C (g/m^3)
  do ii = 1, nscal
    icodcl(ifac,isca(ii)) = 1
    rcodcl(ifac,isca(ii),1) = 1.d0
  enddo
enddo
!< [richards_BC_ex1]

!< [richards_BC_ex2]
call getfbr('FACE2', nlelt, lstelt)

do ilelt = 1, nlelt

! Face number
  ifac = lstelt(ilelt)

! Undefined type face
  itypfb(ifac) = iindef

! Velocity is not solved but deduced from darcy law, then no BC are required.
! However, no flux BC are imposed for safety.
  icodcl(ifac,iu) = 3
  icodcl(ifac,iv) = 3
  icodcl(ifac,iw) = 3
  rcodcl(ifac,iu,3) = 0.d0
  rcodcl(ifac,iv,3) = 0.d0
  rcodcl(ifac,iw,3) = 0.d0

! Dirichlet BC on hydraulic head (H = h + z) to impose a constant gradient over z axis
! here \f$ \grad H \cdot \vect{z} = -0.1 \f$
  icodcl(ifac,ipr) = 1
  rcodcl(ifac,ipr,1) = -1.d0 * cdgfbo(3,ifac)
enddo
!< [richards_BC_ex2]

! ---- Neumann on hydraulic head and solutes

!< [richards_BC_ex3]
call getfbr('FACE3', nlelt, lstelt)

do ilelt = 1, nlelt

! Face number
  ifac = lstelt(ilelt)

! Undefined type face
  itypfb(ifac) = iindef

! Velocity is not solved but deduced from darcy law, then no BC are required.
! However, no flux BC are imposed for safety.
  icodcl(ifac,iu) = 3
  icodcl(ifac,iv) = 3
  icodcl(ifac,iw) = 3
  rcodcl(ifac,iu,3) = 0.d0
  rcodcl(ifac,iv,3) = 0.d0
  rcodcl(ifac,iw,3) = 0.d0

! Neumann BC on hydraulic head (H = h + z) to impose a gradient among the surface normal
! Here \f$ \grad H \cdot \vect{n} = 0 \f$
  icodcl(ifac,ipr) = 3
  rcodcl(ifac,ipr,3) = 0.d0

! Neumann BC on centration C for boundary surface with outward or null normal flow $V_out$:
! It allows to impose a gradient among the surface normal for a diffusive flux as the
! convective flux is defined by the $C*V_out$.
! Here \f$ \grad C \cdot \vect{n} = 0 \f$
  do ii = 1, nscal
    icodcl(ifac,isca(ii)) = 3
    rcodcl(ifac,isca(ii),3) = 0.d0
  enddo
enddo
!< [richards_BC_ex3]

! ---- Mixed (or Robin) BC to impose a total flux (diffusive and convective flux)

!< [richards_BC_ex4]
call getfbr('FACE4', nlelt, lstelt)

! Pointers to the mass flux
call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

do ilelt = 1, nlelt

! Face number
  ifac = lstelt(ilelt)

! Undefined type face
  itypfb(ifac) = iindef

! Velocity is not solved but deduced from darcy law, then no BC are required.
! However, no flux BC are imposed for safety.
  icodcl(ifac,iu) = 3
  icodcl(ifac,iv) = 3
  icodcl(ifac,iw) = 3
  rcodcl(ifac,iu,3) = 0.d0
  rcodcl(ifac,iv,3) = 0.d0
  rcodcl(ifac,iw,3) = 0.d0

! Dirichlet BC on hydraulic head (H = h + z) to impose a constant value
  icodcl(ifac,ipr) = 1
  rcodcl(ifac,ipr,1) = 10.d0

! In order to impose a radioactive activity R_act (Bq) at a surface S_act (m^2)
! during a period of time T_act, a mixed (or Robin) BC is used.
! To do so, two quantities are required:
! * the velocity at an entrance is Vel = - mass_flux / surf
!   or Vel = mass_flux / surf at an exit
! * the reference concentration Cref = R_act / (S_act * dt * Vel)
  R_act = 0.5d0
  S_act = 50.d0
  T_act = 10.d0

! The total flux is imposed from 0 to T_act
  Vel = - bmasfl(ifac)/surfbn(ifac)
  Cref = R_act / (S_act * dt(ntcabs) * Vel)
  if (ttcabs.gt.T_act) Cref = 0.d0
  do ii = 1, nscal
    icodcl(ifac,isca(ii)) = 1
    rcodcl(ifac,isca(ii),1) = Cref
    rcodcl(ifac,isca(ii),2) = Vel
  enddo
enddo
!< [richards_BC_ex4]

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
