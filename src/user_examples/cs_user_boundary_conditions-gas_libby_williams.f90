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
!> \file  cs_user_boundary_conditions-gas_libby_williams.f90
!> \brief Example of cs_f_user_boundary_conditions subroutine.f90 for
!>        Libby-Williams gas
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
integer          ifac, izone, ii
integer          ilelt, nlelt

double precision uref2, xkent, xeent, d2s3

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

! Definition of a burned gas inlet (pilot flame) for each face of color 11

!< [example_1]
call getfbr('11', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  !  Type of pre-defined boundary condidition (see above)
  itypfb(ifac) = ientre

  !  Zone number (arbitrary number between 1 and n)
  izone = 1

  !  Allocation of the actual face to the zone
  izfppp(ifac) = izone

  !  Indicating the inlet as a burned gas inlet
  ientgb(izone) = 1
  !  The incoming burned gas flow refers to:
  !  a) a massflow rate   -> iqimp()  = 1
  iqimp(izone) = 0
  qimp (izone) = zero

  !  b) an inlet velocity -> iqimp()  = 0
  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 21.47d0
  ! ATTENTION: If iqimp()  = 1 the direction vector of the massfow has
  !           to be given here.

  !  Mean Mixture Fraction at Inlet
  fment(izone) = 1.d0*fs(1)
  !  Inlet Temperature in K
  tkent(izone) = 2000.d0

  !  Boundary Conditions of Turbulence
  icalke(izone) = 1

  !   - If ICALKE = 0 the boundary conditions of turbulence at
  !     the inlet are calculated as follows:

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

  !   - If ICALKE = 1 the boundary conditions of turbulence at
  !     the inlet refer to both, a hydraulic diameter and a
  !     reference velocity.
  !
  dh(izone)     = 0.032d0

  !   - If ICALKE = 2 the boundary conditions of turbulence at
  !     the inlet refer to a turbulence intensity.

  xintur(izone) = 0.d0

enddo
!< [example_1]

! Definition of an unburned gas inlet for each face of color 12

!< [example_2]
call getfbr('12', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  !  Type of pre-defined boundary condidition (see above)
  itypfb(ifac) = ientre

  !  Zone number (arbitrary number between 1 and n)
  izone = 2

  !  Allocation of the actual face to the zone
  izfppp(ifac) = izone

  !  Indicating the inlet as an unburned gas inlet
  ientgf(izone) = 1
  !  The incoming unburned gas flow refers to:
  !  a) a massflow rate   -> iqimp()  = 1
  iqimp(izone) = 0
  qimp(izone)  = zero

  !  b) an inlet velocity -> iqimp()  = 0
  rcodcl(ifac,iu,1) = 60.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 0.d0
  ! ATTENTION: If iqimp()  = 1 the direction vector of the massfow has
  !           to be given here.

  !  Mean Mixture Fraction at Inlet
  fment(izone) = 8.d-1*fs(1)

  !  Inlet Temperature in K
  tkent(izone) = 600.d0

  !  Boundary Conditions of Turbulence
  icalke(izone) = 1

  !   - If ICALKE = 1 the boundary conditions of turbulence at
  !    the inlet refer to both, a hydraulic diameter and a
  !    reference velocity.

  dh(izone)     = 0.08d0

  !  - If ICALKE = 2 the boundary conditions of turbulence at
  !    the inlet refer to a turbulence intensity.

  xintur(izone) = 0.1d0

enddo
!< [example_2]

! Definition of a wall for each face of color 51 and 5

!< [example_3]
call getfbr('51 or 5', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Type de condition aux limites pour les variables standard
  itypfb(ifac)   = iparoi

  ! Zone number (arbitrary number between 1 and n)
  izone = 4

  ! Allocation of the actual face to the zone
  izfppp(ifac) = izone

enddo
!< [example_3]

! Definition of an exit for each face of color 91 and 9

!< [example_4]
call getfbr('91 or 9', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Type de condition aux limites pour les variables standard
  itypfb(ifac)   = isolib

  ! Zone number (arbitrary number between 1 and n)
  izone = 5

  ! Allocation of the actual face to the zone 9
  izfppp(ifac) = izone

enddo
!< [example_4]

! Definition of symmetric boundary conditions for each
! face of color 41 and 4.

!< [example_5]
call getfbr('41 or 4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Type de condition aux limites pour les variables standard
  itypfb(ifac)   = isymet

  ! Zone number (arbitrary number between 1 and n)
  izone = 6

  ! Allocation of the actual face to the zonec
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
