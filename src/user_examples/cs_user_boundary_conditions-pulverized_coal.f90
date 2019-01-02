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
integer          icha, iclapc
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

! ---- BOUNDARY FACE corresponding to AIR INLET
!      e.g. : secondary or tertiary air

!< [example_1]
call getfbr('12', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Kind of boundary conditions for standard variables
  itypfb(ifac) = ientre

!   zone's number (from 1 to n)
  izone = 1

!      - Allocation of the zone's number to the face
  izfppp(ifac) = izone

!      - For theses inlet faces, mass flux is fixed

  ientat(izone) = 1
  iqimp(izone)  = 1
!      - Oxidizer's number (1 to 3)
  inmoxy(izone) = 1
!      - Oxidizer mass flow rate  in kg/s
  qimpat(izone) = 1.46d-03
!      - Oxidizer's Temperature in K
  timpat(izone) = 400.d0 + tkelvi

!      - The color 12 becomes an fixed flow rate inlet
!        The user gives speed vector direction
!        (speed vector norm is irrelevent)

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

! ------ Turbulence treatment
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

! ------ Automatic treatment of scalars for extended physic


! ------ treatment of user's scalars

  if ( (nscal-nscapp).gt.0 ) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1]

! ---- BOUNDARY FACE for pulverised COAL & primary air INLET

!< [example_2]
call getfbr('11', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Kind of boundary conditions for standard variables
  itypfb(ifac) = ientre

!   zone's number (from 1 to n)
  izone = 2

!      - Allocation of the zone's number to the face
  izfppp(ifac) = izone

!      - For theses inlet faces, mass flux is fixed

  ientcp(izone) = 1
  iqimp(izone)  = 1
!      - Oxidizer's number (1 to 3)
  inmoxy(izone) = 1
!      - Oxidizer's mass flow rate in kg/s
  qimpat(izone) = 1.46d-03
!      - Oxidizer's Temperature in K
  timpat(izone) = 800.d0  + tkelvi

!        Coal inlet, initialization
  do icha = 1, ncharm
    qimpcp(izone,icha) = zero
    timpcp(izone,icha) = zero
    do iclapc = 1, ncpcmx
      distch(izone,icha,iclapc) = zero
    enddo
  enddo

! Code_Saturne deals with NCHA different coals (component of blend)
!       every coal is described by NCLPCH(icha) class of particles
!       (each of them described by an inlet diameter)
!
!      - Treatment for the first coal
  icha = 1
!      - Coal mass flow rate in kg/s
   qimpcp(izone,icha) = 1.46d-4
!      - PERCENTAGE mass fraction of each granulometric class
  do iclapc = 1, nclpch(icha)
    distch(izone,icha,iclapc) = 100.d0/dble(nclpch(icha))
  enddo
!      - Inlet temperature for coal & primary air
  timpcp(izone,icha) = 800.d0 + tkelvi

!      - The color 11 becomes an fixed flow rate inlet
!        The user gives speed vector direction
!        (speed vector norm is irrelevent)

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

! PPl
! ------ Traitement de la turbulence

!        La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

!      Choix pour le calcul automatique ICALKE = 1 ou 2
  icalke(izone) = 1
!      Saisie des donnees
  dh(izone)     = 0.1d0
  xintur(izone) = 0.1d0

! PPl
!
enddo
!< [example_2]

!     The color 15 become a WALL

!< [example_3]
call getfbr('15', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          WALL : NUL MASS FLUX (PRESSURE FLUX is zero valued)
!                 FRICTION FOR SPEED (& TURBULENCE)
!                 NUL SCALAR FLUX

!   Kind of boundary conditions for standard variables
  itypfb(ifac)   = iparoi


!   zone's number (from 1 to n)
  izone = 3

!      - Allocation of the zone's number to the face
  izfppp(ifac) = izone

enddo
!< [example_3]


!     The color 19 becomes an OUTLET

!< [example_4]
call getfbr('19', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          OUTLET : NUL FLUX for SPEED & SCALARS, FIXED PRESSURE

!   Kind of boundary conditions for standard variables
  itypfb(ifac)   = isolib

!   zone's number (from 1 to n)
  izone = 4

!      - Allocation of the zone's number to the face
  izfppp(ifac) = izone

enddo
!< [example_4]

!     The color 14 becomes a symetry plane

!< [example_5]
call getfbr('14 or 4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRIES

!   Kind of boundary conditions for standard variables
  itypfb(ifac)   = isymet

!   zone's number (from 1 to n)
  izone = 5

!      - Allocation of the zone's number to the face
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
