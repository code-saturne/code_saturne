!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
integer          ifac, ii
integer          izone
integer          ilelt, nlelt
double precision d2s3
double precision zref, xuref
double precision ustar, rugd, rugt
double precision zent, xuent, xvent
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt
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
rugt = 0.1d0

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
!==========
!   - Zone number (from 1 to n)
izone = 1

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the face belongs
  izfppp(ifac) = izone

  ! - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

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
!==========
!   -  Zone number (from 1 to n)
izone = 2

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the face belongs
  izfppp(ifac) = izone

  ! - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

  ! - Assign inlet boundary conditions
  itypfb(ifac) = ientre

enddo
!< [example_2]

! --- For boundary faces of color 31,
!       assign an inlet boundary condition for all phases prescribed from the
!       meteo profile except for dynamical variables which are prescribed with
!       a rough log law.

!< [example_3]
call getfbr('31',nlelt,lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 3

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the face belongs
  izfppp(ifac) = izone

  ! - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

  ! - Dynamical variables are prescribed with a rough log law
  zent=cdgfbo(3,ifac)

  ustar=xkappa*xuref/log((zref+rugd)/rugd)
  xuent=ustar/xkappa*log((zent+rugd)/rugd)
  xvent = 0.d0
  xkent=ustar**2/sqrt(cmu)
  xeent=ustar**3/xkappa/(zent+rugd)

  itypfb(ifac) = ientre

  rcodcl(ifac,iu,1) = xuent
  rcodcl(ifac,iv,1) = xvent
  rcodcl(ifac,iw,1) = 0.d0

  ! itytur is a flag equal to iturb/10
  if    (itytur.eq.2) then

    rcodcl(ifac,ik,1)  = xkent
    rcodcl(ifac,iep,1) = xeent

  elseif(itytur.eq.3) then

    rcodcl(ifac,ir11,1) = d2s3*xkent
    rcodcl(ifac,ir22,1) = d2s3*xkent
    rcodcl(ifac,ir33,1) = d2s3*xkent
    rcodcl(ifac,ir12,1) = 0.d0
    rcodcl(ifac,ir13,1) = 0.d0
    rcodcl(ifac,ir23,1) = 0.d0
    rcodcl(ifac,iep,1)  = xeent

  elseif(iturb.eq.50) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iep,1)  = xeent
    rcodcl(ifac,iphi,1) = d2s3
    rcodcl(ifac,ifb,1)  = 0.d0

  elseif(iturb.eq.60) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iomg,1) = xeent/cmu/xkent

  elseif(iturb.eq.70) then

    rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

  endif

enddo
!< [example_3]

! --- Prescribe at boundary faces of color '12' an outlet for all phases
call getfbr('12', nlelt, lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 4

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the zone belongs
  izfppp(ifac) = izone

  ! Outlet: zero flux for velocity and temperature, prescribed pressure
  !         Note that the pressure will be set to P0 at the first
  !         free outlet face (isolib)

  itypfb(ifac)   = isolib

enddo

! --- Prescribe at boundary faces of color 15 a rough wall for all phases
!< [example_4]
call getfbr('15', nlelt, lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 5

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       rough friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  ! - Zone to which the zone belongs
  izfppp(ifac) = izone

  itypfb(ifac)   = iparug

  ! Roughness for velocity: rugd
  rcodcl(ifac,iu,3) = rugd

  ! Roughness for scalars (if required):
  rcodcl(ifac,iv,3) = rugt


  if(iscalt.ne.-1) then

    ! If temperature prescribed to 20 with a rough wall law (scalar ii=1)
    ! (with thermal roughness specified in rcodcl(ifac,iv,3)) :
    ii = iscalt
    icodcl(ifac, isca(ii))    = 6
    rcodcl(ifac, isca(ii),1)  = 293.15d0

    ! If flux prescribed to 4.d0 (scalar ii=2):
    ii = 2
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 4.d0

  endif
enddo
!< [example_4]

! --- Prescribe at boundary faces of color 4 a symmetry for all phases
!< [example_5]
call getfbr('4', nlelt, lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 6

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! - Zone to which the zone belongs
  izfppp(ifac) = izone

  itypfb(ifac)   = isymet

enddo
!< [example_5]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
