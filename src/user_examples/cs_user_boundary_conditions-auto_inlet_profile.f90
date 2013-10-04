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

! Example of cs_user_boundary_conditions subroutine.f90 for inlet
! with automatic inlet profile.

! This example assumes the mesh is orthogonal at the inlet.
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
integer          ifac, iel, ii, ivar, irangv, iun, ilelt, nlelt
double precision xkent, xeent, xphi, xfb

integer          ibc, ix, ny, ios
double precision uent, vent, went, xustar2, xdh, d2s3, rhomoy
double precision acc(2), fmprsc, fmul, uref2, vnrm

integer, allocatable, dimension(:) :: lstelt, mrkcel
double precision, dimension(:), pointer :: brom
!< [loc_var_dec]

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0
!< [init]

!===============================================================================
! Assign a pseudo-periodic channel type inlet to a set of boundary faces.

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face

! A feedback loop is used so as to progressively reach a state similar
! to that of a periodic channel at the inlet.
!===============================================================================

!< [example_1]
call getfbr('INLET', nlelt, lstelt)
!==========

fmprsc = 1.d0 ! mean prescribed velocity

if (ntcabs.eq.1) then

  ! For the Rij-EBRSM model (and possibly V2f), we need a non-flat profile,
  ! so as to ensure turbulent production, and avoid laminarization;
  ! here, we simply divide the initial velocity by 10 for inlet
  ! faces adjacent to the wall.

  ! The loop below assumes wall conditions have been defined first
  ! (in the GUI, or in this file, before the current test).

  if (iturb.eq.32 .or. itytur.eq.5) then

    allocate(mrkcel(ncelet))
    do iel = 1, ncelet
      mrkcel(iel) = 0
    enddo

    do ifac = 1, nfabor
      if (itypfb(ifac) .eq. iparoi) then
        iel = ifabor(ifac)
        mrkcel(iel) = 1
      endif
    enddo

  endif

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel = ifabor(ifac)

    itypfb(ifac) = ientre

    rcodcl(ifac,iu,1) = fmprsc * surfbo(1,ifac) / surfbn(ifac)
    rcodcl(ifac,iv,1) = fmprsc * surfbo(2,ifac) / surfbn(ifac)
    rcodcl(ifac,iw,1) = fmprsc * surfbo(3,ifac) / surfbn(ifac)

    if (iturb.eq.32 .or. itytur.eq.5) then
      rcodcl(ifac,iu,1) = fmprsc/10.d0
    endif

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
    xdh     = 1.d0

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
    if (itytur.eq.2) then

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

  if (iturb.eq.32 .or. itytur.eq.5) then
    deallocate(mrkcel)
  endif

else

! Subsequent time steps
!----------------------

  acc(1) = 0.d0
  acc(2) = 0.d0

  ! Estimate multiplier

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel = ifabor(ifac)

    vnrm = sqrt(rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2)
    acc(1) = acc(1) + vnrm*surfbn(ifac)
    acc(2) = acc(2) + surfbn(ifac)

  enddo

  if (irangp.ge.0) then
    call parrsm(2, acc)
  endif

  fmul = fmprsc/(acc(1)/acc(2)) ! 1 / estimate flow multiplier

  ! Apply BC

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel = ifabor(ifac)

    itypfb(ifac) = ientre

    vnrm = sqrt(rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2)

    rcodcl(ifac,iu,1) = fmul * vnrm * surfbo(1,ifac) / surfbn(ifac)
    rcodcl(ifac,iv,1) = fmul * vnrm * surfbo(2,ifac) / surfbn(ifac)
    rcodcl(ifac,iw,1) = fmul * vnrm * surfbo(3,ifac) / surfbn(ifac)

    if (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = rtp(iel,ik)
      rcodcl(ifac,iep,1) = rtp(iel,iep)

    elseif (itytur.eq.3) then

      rcodcl(ifac,ir11,1) = rtp(iel,ir11)
      rcodcl(ifac,ir22,1) = rtp(iel,ir22)
      rcodcl(ifac,ir33,1) = rtp(iel,ir33)
      rcodcl(ifac,ir12,1) = rtp(iel,ir12)
      rcodcl(ifac,ir13,1) = rtp(iel,ir13)
      rcodcl(ifac,ir23,1) = rtp(iel,ir23)
      rcodcl(ifac,iep,1)  = rtp(iel,iep)

      if (iturb.eq.32) then
        rcodcl(ifac,ial,1)  = rtp(iel,ial)
      endif

    elseif (itytur.eq.5) then

      rcodcl(ifac,ik,1)  = rtp(iel,ik)
      rcodcl(ifac,iep,1) = rtp(iel,iep)
      rcodcl(ifac,iphi,1) = rtp(iel,iphi)

      if (iturb.eq.50) then
        rcodcl(ifac,ifb,1)  = rtp(iel,ifb)
      elseif (iturb.eq.51) then
        rcodcl(ifac,ial,1)  = rtp(iel,ial)
      endif

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)  = rtp(iel,ik)
      rcodcl(ifac,iomg,1) = rtp(iel,iomg)

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = rtp(iel,inusa)

    endif

    ! Handle scalars (a correction similar to that of velocity is suggested
    !                 rather than the simpler code below)
    if (nscal.gt.0) then
      do ii = 1, nscal
        rcodcl(ifac,isca(ii),1) = rtp(iel,isca(ii))
      enddo
    endif

  enddo

endif
!< [example_1]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_user_boundary_conditions
