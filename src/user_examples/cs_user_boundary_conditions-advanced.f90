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

!> \file cs_user_boundary_conditions-advanced.f90
!>
!> \brief Advanced example of cs_user_boundary_conditions.f90 subroutine
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
integer          ifac, iel, ii, ivar
integer          izone
integer          ilelt, nlelt
double precision uref2, d2s3
double precision rhomoy, xdh, xustar2
double precision xitur
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: boundary_roughness
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0
!< [init]

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! Example of specific boundary conditions fully defined by the user,
! on the basis of wall conditions.
! selection (mass flow computation, specific logging, ...)
! We prescribe for group '1234' a wall, with in addition:
!   - a Dirichlet condition on velocity (sliding wall with no-slip condition)
!   - a Dirichlet condition on the first scalar.

!< [example_1]
call getfbr('1234', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac) = iparoi

  icodcl(ifac,iu )  = 1
  rcodcl(ifac,iu,1) = 1.d0
  rcodcl(ifac,iu,2) = rinfin
  rcodcl(ifac,iu,3) = 0.d0
  icodcl(ifac,iv )  = 1
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iv,2) = rinfin
  rcodcl(ifac,iv,3) = 0.d0
  icodcl(ifac,iw )  = 1
  rcodcl(ifac,iw,1) = 0.d0
  rcodcl(ifac,iw,2) = rinfin
  rcodcl(ifac,iw,3) = 0.d0

  ivar = isca(1)
  icodcl(ifac,ivar )  = 1
  rcodcl(ifac,ivar,1) = 10.d0
  rcodcl(ifac,ivar,2) = rinfin
  rcodcl(ifac,ivar,3) = 0.d0

enddo
!< [example_1]

! Example of specific boundary conditions fully defined by the user,
! with no definition of a specific type.
! We prescribe at group '5678' a homogeneous Neumann condition for
! all variables.

!< [example_2]
call getfbr('5678', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! CAUTION: the value of itypfb must be assigned to iindef

  itypfb(ifac) = iindef

  do ii = 1, nvar
    icodcl(ifac,ii )  = 3
    rcodcl(ifac,ii,1) = 0.d0
    rcodcl(ifac,ii,2) = rinfin
    rcodcl(ifac,ii,3) = 0.d0
  enddo

enddo
!< [example_2]

! Example of specific boundary conditions fully defined by the user,
! with the definition of a specific type, for example for future
! selection (mass flow computation, specific logging, ...)
! We prescribe for group '6789' a homogeneous Neumann condition for
! all variables, except for the first
! scalar, for which we select a homogeneous Dirichlet.

!< [example_3]
call getfbr('6789', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! CAUTION: the value of itypfb must be different from
  !          iparoi, ientre, isymet, isolib, iindef,
  !          greater than or equal to 1, and
  !          less than or equal to ntypmx;
  !          these integers are defined in paramx.h

  itypfb(ifac) = 89

  do ii = 1, nvar
    icodcl(ifac,ii )  = 3
    rcodcl(ifac,ii,1) = 0.d0
    rcodcl(ifac,ii,2) = rinfin
    rcodcl(ifac,ii,3) = 0.d0
  enddo

  icodcl(ifac,isca(1) )  = 1
  rcodcl(ifac,isca(1),1) = 0.d0
  rcodcl(ifac,isca(1),2) = rinfin
  rcodcl(ifac,isca(1),3) = 0.d0

enddo
!< [example_3]

! Example of wall boundary condition with automatic continuous switch
! between rough and smooth.

!< [example_4]
call field_get_val_s_by_name("boundary_roughness", boundary_roughness)

call getfbr('6789', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! CAUTION: wall function keyword "iwallf" must be set to 6
  !          in cs_user_parameters.

  itypfb(ifac) = iparoi

  ! Boundary roughtness (in meter)
  boundary_roughness(ifac) = 0.05

enddo
!< [example_4]

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
