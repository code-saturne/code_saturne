!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> \file  cs_user_boundary_conditions-fuel.f90
!> \brief Example of cs_f_user_boundary_conditions subroutine.f90 for fuel
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
integer          iclafu
integer          ilelt, nlelt

double precision uref2, d2s3
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection


d2s3 = 2.d0/3.d0

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

!  Each kind of condition for extended physic is allocated with a number
!  IZONE ( 0<IZONE<= NOZPPM ; NOZPPM allocated in ppppar.h)

! ---- The 12 color is a pure air inlet

!< [example_1]
call getfbr('12', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   kind of boundary condition for standard variables
  itypfb(ifac) = ientre

!   Zone number allocation
  izone = 1

!   Zone number storage
  izfppp(ifac) = izone

! ------ This inlet have a fixed mass flux

  ientat(izone) = 1
  iqimp(izone)  = 1
  inmoxy(izone) = 1
!     - Air mass flow rate kg/s
  qimpat(izone) = 1.46d-03
!     - Air inlet temperature K
  timpat(izone) = 400.d0 + tkelvi
!     - Fuel mass flow rate kg/s
  qimpfl(izone) = 0.d0

! ----- The 12 color is a fixed mass flow rate inlet
!        user gives only the speed vector direction
!        (spedd vector norm is irrelevent)
!
  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

!   Boundary conditions of turbulence
  icalke(izone) = 1
!
!    - If ICALKE = 0 the boundary conditions of turbulence at
!      the inlet are calculated as follows:

  if (icalke(izone).eq.0) then

    uref2 = rcodcl(ifac,iu,1)**2                           &
           +rcodcl(ifac,iv,1)**2                           &
           +rcodcl(ifac,iw,1)**2
    uref2 = max(uref2,1.d-12)
    xkent  = epzero
    xeent  = epzero

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

    elseif (iturb.eq.50) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iep,1)  = xeent
      rcodcl(ifac,iphi,1) = d2s3
      rcodcl(ifac,ifb,1)  = 0.d0

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    endif

  endif
!      the inlet refer to both, a hydraulic diameter and a
!      reference velocity.
!
  dh(izone)     = 0.032d0
!
!    - If ICALKE = 2 the boundary conditions of turbulence at
!      the inlet refer to a turbulence intensity.
!
  xintur(izone) = 0.d0

! ------ Treatment of user's scalars

  if ( (nscal-nscapp).gt.0 ) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1]

! ---- Inlet of both primary Air and Fuel

!< [example_2]
call getfbr('11', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   kind of boundary condition for standard variables
  itypfb(ifac) = ientre

!   Zone number allocation
  izone = 2

!   Zone number storage
  izfppp(ifac) = izone

! ------ This inlet have a fixed mass flux

  ientfl(izone) = 1
  iqimp(izone)  = 1
  inmoxy(izone) = 1
!     - Air mass flow rate in kg/s
  qimpat(izone) = 1.46d-03
!     - Air Temperature at inlet in K
  timpat(izone) = 800.d0  + tkelvi

!     - Fuel mass flow rate in kg/s
  qimpfl(izone) = 1.46d-04/360.d0

!     - PERCENTAGE mass fraction of each granulometric class
!       ICLAFU (1 < ICLAFU < NCLAFU )
  iclafu = 1
  distfu(izone,iclafu) = 100.d0

!     - Fuel Temperature at inlet in K
  timpfl(izone) = 100.d0  + tkelvi

! ----- The 11 color is a fixed mass flow rate inlet
!        user gives only the speed vector direction
!        (spedd vector norm is irrelevent)
!

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

!   Boundary conditions of turbulence
  icalke(izone) = 1
!
!    - If ICALKE = 0 the boundary conditions of turbulence at
!      the inlet are calculated as follows:

  if(icalke(izone).eq.0) then

    uref2 = rcodcl(ifac,iu,1)**2                           &
           +rcodcl(ifac,iv,1)**2                           &
           +rcodcl(ifac,iw,1)**2
    uref2 = max(uref2,1.d-12)
    xkent  = epzero
    xeent  = epzero

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

    elseif (iturb.eq.50) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iep,1)  = xeent
      rcodcl(ifac,iphi,1) = d2s3
      rcodcl(ifac,ifb,1)  = 0.d0

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    endif

  endif
!
!    - If ICALKE = 1 the boundary conditions of turbulence at
!      the inlet refer to both, a hydraulic diameter and a
!      reference velocity.
!
  dh(izone)     = 0.032d0
!
!    - If ICALKE = 2 the boundary conditions of turbulence at
!      the inlet refer to a turbulence intensity.
!
  xintur(izone) = 0.d0
enddo
!< [example_2]

! --- Color 15 is a wall

!< [example_3]
call getfbr('15', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          WALL  : nul mass flow rate (nul pressure flux)
!                  rubbing for speed (and turbulence)
!                  nul scalar fluxes

!   kind of boundary condition for standard variables
  itypfb(ifac)   = iparoi


!   Zone number allocation
  izone = 3

!   Zone number storage
  izfppp(ifac) = izone

enddo
!< [example_3]

! --- Color 19 is an outlet

!< [example_4]
call getfbr('19', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          OUTLET : nul fluxes for speed and scalar
!                   pressure fixed

!   kind of boundary condition for standard variables
    itypfb(ifac)   = isolib

!   Zone number allocation
    izone = 4

!   Zone number storage
    izfppp(ifac) = izone

  enddo
!< [example_4]

! --- 14 and 4 are symetry

!< [example_5]
call getfbr('14 or 4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRY

!   kind of boundary condition for standard variables
  itypfb(ifac)   = isymet

!   Zone number allocation
  izone = 5

!   Zone number storage
  izfppp(ifac) = izone

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
