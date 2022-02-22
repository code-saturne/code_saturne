!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

! Example of cs_f_user_boundary_conditions subroutine.f90 for inlet
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
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[in,out] izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2
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
use atsoil
use ctincl
use cs_fuel_incl
use mesh
use field
use turbomachinery
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
integer          ifac, iel, ii, ilelt, nlelt

double precision xustar2, xdh, d2s3, rhomoy
double precision acc(2), fmprsc, fmul, uref2, vnrm

integer, allocatable, dimension(:) :: lstelt, mrkcel

double precision, dimension(:), pointer :: bfpro_rom
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_omg
double precision, dimension(:), pointer :: cvar_al, cvar_fb
double precision, dimension(:), pointer :: cvar_nusa
double precision, dimension(:), pointer :: cvar_scal
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:,:), pointer :: cvar_rij
!< [loc_var_dec]

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

! Map field arrays
call field_get_val_v(ivarfl(iu), cvar_vel)
call field_get_val_s(ibrom, bfpro_rom)

if (itytur.eq.2) then

  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)

elseif (itytur.eq.3) then
  call field_get_val_v(ivarfl(irij), cvar_rij)
  call field_get_val_s(ivarfl(iep), cvar_ep)

  if (iturb.eq.32) then
    call field_get_val_s(ivarfl(ial), cvar_al)
  endif

elseif (itytur.eq.5) then

  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(iphi), cvar_phi)

  if (iturb.eq.50) then
    call field_get_val_s(ivarfl(ifb), cvar_fb)
  elseif (iturb.eq.51) then
    call field_get_val_s(ivarfl(ial), cvar_al)
  endif

elseif (iturb.eq.60) then

  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)

elseif (iturb.eq.70) then

  call field_get_val_s(ivarfl(inusa), cvar_nusa)

endif

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

    rcodcl(ifac,iu,1) = - fmprsc * surfbo(1,ifac) / surfbn(ifac)
    rcodcl(ifac,iv,1) = - fmprsc * surfbo(2,ifac) / surfbn(ifac)
    rcodcl(ifac,iw,1) = - fmprsc * surfbo(3,ifac) / surfbn(ifac)

    if (iturb.eq.32 .or. itytur.eq.5) then
      if (mrkcel(iel) .eq. 1) then
        rcodcl(ifac,iu,1) = fmprsc/10.d0
      endif
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

    ! Calculation of turbulent inlet conditions using
    !   the turbulence intensity and standard laws for a circular pipe
    !   (their initialization is not needed here but is good practice)

    rhomoy  = bfpro_rom(ifac)

    call turbulence_bc_inlet_hyd_diam(ifac, uref2, xdh, rhomoy, viscl0,  &
                                      rcodcl)

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

    vnrm = sqrt(cvar_vel(1,iel)**2 + cvar_vel(2,iel)**2 + cvar_vel(3,iel)**2)
    acc(1) = acc(1) + vnrm*surfbn(ifac)
    acc(2) = acc(2) + surfbn(ifac)

  enddo

  if (irangp.ge.0) then
    call parrsm(2, acc)
  endif

  if (acc(1).le.epzero) then
     fmul = 0.d0 ! zero velocity in bulk domain
  else
     fmul = fmprsc/(acc(1)/acc(2)) ! 1 / estimate flow multiplier
  endif

  ! Apply BC

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel = ifabor(ifac)

    itypfb(ifac) = ientre

    vnrm = sqrt(cvar_vel(1,iel)**2 + cvar_vel(2,iel)**2 + cvar_vel(3,iel)**2)

    rcodcl(ifac,iu,1) = - fmul * vnrm * surfbo(1,ifac) / surfbn(ifac)
    rcodcl(ifac,iv,1) = - fmul * vnrm * surfbo(2,ifac) / surfbn(ifac)
    rcodcl(ifac,iw,1) = - fmul * vnrm * surfbo(3,ifac) / surfbn(ifac)

    if (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = cvar_k(iel)
      rcodcl(ifac,iep,1) = cvar_ep(iel)

    elseif (itytur.eq.3) then

        rcodcl(ifac,ir11,1) = cvar_rij(1,iel)
        rcodcl(ifac,ir22,1) = cvar_rij(2,iel)
        rcodcl(ifac,ir33,1) = cvar_rij(3,iel)
        rcodcl(ifac,ir12,1) = cvar_rij(4,iel)
        rcodcl(ifac,ir13,1) = cvar_rij(6,iel)
        rcodcl(ifac,ir23,1) = cvar_rij(5,iel)

      rcodcl(ifac,iep,1)  = cvar_ep(iel)

      if (iturb.eq.32) then
        rcodcl(ifac,ial,1)  = cvar_al(iel)
      endif

    elseif (itytur.eq.5) then

      rcodcl(ifac,ik,1)  = cvar_k(iel)
      rcodcl(ifac,iep,1) = cvar_ep(iel)
      rcodcl(ifac,iphi,1) = cvar_phi(iel)

      if (iturb.eq.50) then
        rcodcl(ifac,ifb,1)  = cvar_fb(iel)
      elseif (iturb.eq.51) then
        rcodcl(ifac,ial,1)  = cvar_al(iel)
      endif

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)  = cvar_k(iel)
      rcodcl(ifac,iomg,1) = cvar_omg(iel)

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = cvar_nusa(iel)

    endif

    ! Handle scalars (a correction similar to that of velocity is suggested
    !                 rather than the simpler code below)
    if (nscal.gt.0) then
      do ii = 1, nscal
        call field_get_val_s(ivarfl(isca(ii)), cvar_scal)
        rcodcl(ifac,isca(ii),1) = cvar_scal(iel)
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
end subroutine cs_f_user_boundary_conditions
