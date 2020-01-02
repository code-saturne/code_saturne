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
!> \file  cs_user_boundary_conditions-gas_3pthem.f90
!> \brief Example of cs_user_boundary_conditions subroutine.f90 for 3 PTHEM gas
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
integer          iel, ifac, izone, ii
integer          ilelt, nlelt

double precision uref2, d2s3
double precision xkent, xeent
double precision dp, rho, S0, Kadm, radm, madm, Qadm, Padm, Ploc
double precision, dimension(:), pointer :: crom

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

!   Definition of a fuel flow inlet for each face of color 11

!< [example_1]
call getfbr('11', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type of pre-defined boundary condidition (see above)
  itypfb(ifac) = ientre

!   Zone number (arbitrary number between 1 and n)
  izone = 1

!   Allocation of the actual face to the zone
  izfppp(ifac) = izone

!   Indicating the inlet as a fuel flow inlet
  ientfu(izone) = 1

!   Inlet Temperature in K
  tinfue = 436.d0

!   The incoming fuel flow refers to:
!   a) a massflow rate   -> iqimp()  = 1
  iqimp(izone)  = 1
  qimp(izone)   = 2.62609d-4 / 72.d0
!
!   b) an inlet velocity -> iqimp()  = 0
  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 21.47d0
! ATTENTION: If iqimp()  = 1 the direction vector of the massfow has
!            to begiven here.
!
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

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

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
!< [example_1]

!   Definition of an air flow inlet for each face of color 21

!< [example_2]
call getfbr('21', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type of pre-defined boundary condidition (see above)
  itypfb(ifac) = ientre

!   Zone number (arbitrary number between 1 and n)
  izone = 2

!   Allocation of the actual face to the zone
  izfppp(ifac) = izone

!   Indicating the inlet as a air flow inlet
  ientox(izone) = 1

!   Inlet Temperature in K
  tinoxy = 353.d0

!   The inflowing fuel flow refers to:
!   a) a massflow rate   -> iqimp()  = 1
  iqimp(izone)  = 1
  qimp(izone)   = 4.282d-3 / 72.d0
!
!   b) an inlet velocity -> iqimp()  = 0
  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 0.097d0
! ATTENTION: If iqimp()  = 1 the direction vector of the massfow has
!            to be given here.

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

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

    endif

  endif
!
!    - If ICALKE = 1 the boundary conditions of turbulence at
!      the inlet refer to both, a hydraulic diameter and a
!      reference velocity.
!
  dh(izone)     = 0.218d0
!
!    - If ICALKE = 2 the boundary conditions of turbulence at
!      the inlet refer to a turbulence intensity.
!
  xintur(izone) = 0.d0
!
enddo
!< [example_2]

!  Definition of a wall for each face of color 51 up to 59

!< [example_3]
call getfbr('51 to 59', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type of pre-defined boundary condidition (see above)
  itypfb(ifac)   = iparoi

!   Zone number (arbitrary number between 1 and n)
  izone = 3

!   Allocation of the actual face to the zone
  izfppp(ifac) = izone

enddo
!< [example_3]


!  Definition of an exit for each face of color 91

!< [example_4]
call getfbr('91', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type of pre-defined boundary condidition (see above)
  itypfb(ifac)   = isolib

!   Zone number (arbitrary number between 1 and n)
  izone = 4

!   Allocation of the actual face to the zone
  izfppp(ifac) = izone

enddo
!< [example_4]

!  Definition of symmetric boundary conditions for each
!  face of color 41 and 4.

!< [example_5]
call getfbr('41 or 4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type of pre-defined boundary condidition (see above)
  itypfb(ifac)   = isymet

!   Zone number (arbitrary number between 1 and n)
  izone = 5

!   Allocation of the actual face to the zone
  izfppp(ifac) = izone

enddo
!< [example_5]

! Definition of an air supply for each face of color 61

!< [example_6]
call getfbr('61', nlelt, lstelt)
!==========

call field_get_val_s(icrom, crom)

! Head losses in the ventilation pipe (often from exp data)

! nominal room pressure
Ploc= 51.d0 ! Pa
! vent surface
S0 = acos(-1.d0) * (75.d-3)**2 ! m2
! density
radm = (P0+ploc)/(cs_physical_constants_r*tinoxy/wmolg(2))
! nominal vent pressure
Padm = 119.d0 ! Pa
! nominal flow rate
Qadm = 504.d0 ! m3/h
madm = Qadm*radm/3600.d0
! head loss
Kadm = 2.d0*radm*(S0/madm)**2 * abs(Ploc-Padm)

! pressure gradient (Pther computed if ipthrm = 1)
dp   = (Pther - P0) - Padm

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel  = ifabor(ifac)

  ! Density
  if (dp.gt.0.d0) then
    rho = crom(iel)
  else
    rho = radm
  endif

  ! Mass flow rate (opposed to the boundary face normal vector)
  madm = - sign(1.d0,dp) * s0 * sqrt(2.d0*rho/kadm*abs(dp))

  ! Type of pre-defined boundary condidition (see above)
  itypfb(ifac) = i_convective_inlet

  ! Zone number (arbitrary number between 1 and n)
  izone = 6

  ! Allocation of the actual face to the zone
  izfppp(ifac) = izone

  ! Indicating the inlet as a air flow inlet
  ientox(izone) = 1

  ! Inlet Temperature in K
  tinoxy = t0

  ! The inflowing fuel flow refers to:
  ! a) a massflow rate   -> iqimp()  = 1
  iqimp(izone)  = 1
  qimp(izone)   = madm

  ! b) an inlet velocity -> iqimp()  = 0
  rcodcl(ifac,iu,1) = 1.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 0.d0

  ! Boundary conditions of turbulence
  icalke(izone) = 1

  ! - If ICALKE = 0 the boundary conditions of turbulence at
  !   the inlet are calculated as follows:
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

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

    endif

  endif

  ! - If ICALKE = 1 the boundary conditions of turbulence at
  !   the inlet refer to both, a hydraulic diameter and a
  !   reference velocity.
  dh(izone)     = 0.218d0

  ! - If ICALKE = 2 the boundary conditions of turbulence at
  !   the inlet refer to a turbulence intensity.
  xintur(izone) = 0.d0

enddo
!< [example_6]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
