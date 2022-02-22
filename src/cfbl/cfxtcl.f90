!-------------------------------------------------------------------------------

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

!> \file cfxtcl.f90
!> \brief Handle boundary condition type code (\ref itypfb) when the
!> compressible model is enabled.
!>
!> Please refer to the
!> <a href="../../theory.pdf#cfxtcl"><b>cfxtcl</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in,out] icodcl        face boundary condition code:
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
!> \param[in]     itypfb        boundary face types
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!______________________________________________________________________________

subroutine cfxtcl &
 ( nvar   ,                                                       &
   icodcl , itypfb , dt     , rcodcl )

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
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar

integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac  , iel
integer          ii    , iccfth
integer          icalep
integer          ien   , itk, niv
integer          nvarcf

integer          nvcfmx
parameter       (nvcfmx=6)
integer          ivarcf(nvcfmx)

double precision hint, bmasfl, drom

double precision, allocatable, dimension(:) :: w5
double precision, allocatable, dimension(:) :: w7
double precision, allocatable, dimension(:) :: wbfb, wbfa
double precision, allocatable, dimension(:) :: bc_en, bc_pr, bc_tk
double precision, allocatable, dimension(:) :: bc_fracv, bc_fracm, bc_frace
double precision, allocatable, dimension(:,:) :: bc_vel

double precision, dimension(:), pointer :: coefbp
double precision, dimension(:), pointer :: crom, brom, cpro_cv, cvar_en
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_fracv, cvar_fracm, cvar_frace

!===============================================================================

! Map field arrays
 call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 1. Initializations
!===============================================================================

! Allocate temporary arrays
allocate(w5(ncelet))
allocate(w7(nfabor), wbfa(nfabor), wbfb(nfabor))
allocate(bc_en(nfabor))
allocate(bc_pr(nfabor))
allocate(bc_tk(nfabor))
allocate(bc_fracv(nfabor))
allocate(bc_fracm(nfabor))
allocate(bc_frace(nfabor))
allocate(bc_vel(3,nfabor))

ien = isca(ienerg)
itk = isca(itempk)

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_s(ivarfl(ien), cvar_en)

if (icv.ge.0) call field_get_val_s(icv, cpro_cv)

! mixture fractions for the homogeneous two-phase flows
if (ippmod(icompf).eq.2) then
  call field_get_val_s(ivarfl(isca(ifracv)), cvar_fracv)
  call field_get_val_s(ivarfl(isca(ifracm)), cvar_fracm)
  call field_get_val_s(ivarfl(isca(ifrace)), cvar_frace)
endif

! list of the variables of the compressible model
ivarcf(1) = ipr
ivarcf(2) = iu
ivarcf(3) = iv
ivarcf(4) = iw
ivarcf(5) = ien
ivarcf(6) = itk
nvarcf    = 6

call field_get_coefb_s(ivarfl(ipr), coefbp)
do ifac = 1, nfabor
  wbfb(ifac) = coefbp(ifac)
enddo

! Computation of epsilon_sup = e - CvT
! Needed if walls with imposed temperature are set.

icalep = 0
do ifac = 1, nfabor
  if(icodcl(ifac,itk).eq.5) then
    icalep = 1
  endif
enddo
if(icalep.ne.0) then
  ! At cell centers
  call cs_cf_thermo_eps_sup(crom, w5, ncel)

  ! At boundary faces centers
  call cs_cf_thermo_eps_sup(brom, w7, nfabor)
endif

! Loop on all boundary faces and treatment of types of BCs given by itypfb

do ifac = 1, nfabor
  iel = ifabor(ifac)

!===============================================================================
! 2. Treatment of all wall boundary faces
!===============================================================================

  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

    ! pressure :
    ! if the gravity is prevailing: hydrostatic pressure
    ! (warning: the density is here explicit and the term is an approximation)

    if(icfgrp.eq.1) then

      icodcl(ifac,ipr) = 15
      hint = dt(iel)/distb(ifac)
      rcodcl(ifac,ipr,3) = -hint                                  &
           * ( gx*(cdgfbo(1,ifac)-xyzcen(1,iel))                  &
           + gy*(cdgfbo(2,ifac)-xyzcen(2,iel))                    &
           + gz*(cdgfbo(3,ifac)-xyzcen(3,iel)) )                  &
           * crom(iel)

    else

      ! generally proportional to the bulk value
      ! (Pboundary = COEFB*Pi)
      ! The part deriving from pinf in stiffened gas is set in explicit for now
      call cs_cf_thermo_wall_bc(wbfa, wbfb, ifac-1)

      if (wbfb(ifac).lt.rinfin*0.5d0.and.wbfb(ifac).gt.0.d0) then
        icodcl(ifac,ipr) = 12
        rcodcl(ifac,ipr,1) = wbfa(ifac)
        rcodcl(ifac,ipr,2) = wbfb(ifac)
      else
        ! If rarefaction is too strong : homogeneous Dirichlet
        icodcl(ifac,ipr) = 13
        rcodcl(ifac,ipr,1) = 0.d0
      endif

    endif

    ! Velocity and turbulence are treated in a standard manner in condli.

    ! For thermal B.C., a pre-treatment has be done here since the solved
    ! variable is the total energy
    ! (internal energy + epsilon_sup + cinetic energy).
    ! Especially, when a temperature is imposed on a wall, clptur treatment
    ! has to be prepared. Except for the solved energy all the variables rho
    ! and s will take arbitrarily a zero flux B.C. (their B.C. are only used
    ! for the gradient reconstruction and imposing something else than zero
    ! flux could bring out spurious values near the boundary layer).

    ! adiabatic by default
    if(  icodcl(ifac,itk).eq.0.and.                          &
         icodcl(ifac,ien).eq.0) then
      icodcl(ifac,itk) = 3
      rcodcl(ifac,itk,3) = 0.d0
    endif

    ! imposed temperature
    if(icodcl(ifac,itk).eq.5) then

      ! The value of the energy that leads to the right flux is imposed.
      ! However it should be noted that it is the B.C. for the diffusion
      ! flux. For the gradient reconstruction, something else will be
      ! needed. For example, a zero flux or an other B.C. respecting a
      ! profile: it may be possible to treat the total energy as the
      ! temperature, keeping in mind that the total energy contains
      ! the cinetic energy, which could make the choice of the profile more
      ! difficult.

      icodcl(ifac,ien) = 5
      if(icv.eq.-1) then
        rcodcl(ifac,ien,1) = cv0*rcodcl(ifac,itk,1)
      else
        rcodcl(ifac,ien,1) = cpro_cv(iel)*rcodcl(ifac,itk,1)
      endif
      rcodcl(ifac,ien,1) = rcodcl(ifac,ien,1)             &
           + 0.5d0*(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)          &
           + w5(iel)
      ! w5 contains epsilon_sup

      ! fluxes in grad(epsilon_sup and cinetic energy) have to be zero
      ! since they are already accounted for in the energy diffusion term
      ifbet(ifac) = 1

      ! Dirichlet condition on the temperature for gradient reconstruction
      ! used only in post-processing (typically Nusselt computation)
      icodcl(ifac,itk) = 1

    ! imposed flux
    elseif(icodcl(ifac,itk).eq.3) then

      ! zero flux on energy
      icodcl(ifac,ien) = 3
      rcodcl(ifac,ien,3) = rcodcl(ifac,itk,3)

      ! fluxes in grad(epsilon_sup and cinetic energy) have to be zero
      ! since they are already accounted for in the energy diffusion term
      ifbet(ifac) = 1

      ! zero flux for the possible temperature reconstruction
      icodcl(ifac,itk) = 3
      rcodcl(ifac,itk,3) = 0.d0

    endif

!===============================================================================
! 3. Treatment of all inlet/outlet boundary faces and thermo step
!===============================================================================

!===============================================================================
! 3.1 Imposed Inlet/outlet (for example: supersonic inlet)
!===============================================================================

  elseif ( itypfb(ifac).eq.iesicf ) then

    ! we have
    !   - velocity,
    !   - 2 variables among P, rho, T, E (but not the couple (T,E)),
    !   - turbulence variables
    !   - scalars

    ! we look for the variable to be initialized
    ! (if a zero value has been given, it is not adapted, so it will
    ! be considered as not initialized and the computation will stop
    ! displaying an error message. The boundary density may
    ! be pre-initialized to the cell density also, so is tested last.

    iel = ifabor(ifac)
    drom = abs(crom(iel) - brom(ifac))

    niv = 0
    iccfth = 10000
    if (rcodcl(ifac,ipr,1).lt.rinfin*0.5d0) then
      iccfth = 2*iccfth
      niv = niv + 1
    endif
    if (rcodcl(ifac,itk,1).lt.rinfin*0.5d0) then
      iccfth = 5*iccfth
      niv = niv + 1
    endif
    if (rcodcl(ifac,ien,1).lt.rinfin*0.5d0) then
      iccfth = 7*iccfth
      niv = niv + 1
    endif

    if (brom(ifac).gt.0.d0 .and. (niv.lt.2 .or. drom.gt.epzero)) then
      iccfth = 3*iccfth
      niv = niv + 1
    endif

    if (niv .ne. 2) then
      write(nfecra,1000) iccfth
      call csexit (1)
    endif
    iccfth = iccfth + 900

    ! missing thermo variables among P, rho, T, E are computed
    bc_en(ifac) = rcodcl(ifac,ien,1)
    bc_pr(ifac) = rcodcl(ifac,ipr,1)
    bc_tk(ifac) = rcodcl(ifac,itk,1)
    bc_vel(1,ifac) = rcodcl(ifac,iu,1)
    bc_vel(2,ifac) = rcodcl(ifac,iv,1)
    bc_vel(3,ifac) = rcodcl(ifac,iw,1)

    call cs_cf_thermo(iccfth, ifac-1, bc_en, bc_pr, bc_tk, bc_vel)

!===============================================================================
! 3.2 Outlet with imposed pressure
!===============================================================================

  elseif ( itypfb(ifac).eq.isopcf ) then

    ! If no value was given for P or if its value is negative, the computation
    ! stops (a negative value could be possible, but in most cases it would be
    ! an error).
    if(rcodcl(ifac,ipr,1).lt.-rinfin*0.5d0) then
      write(nfecra,1100)
      call csexit (1)
    endif

    bc_en(ifac) = rcodcl(ifac,ien,1)
    bc_pr(ifac) = rcodcl(ifac,ipr,1)
    bc_tk(ifac) = rcodcl(ifac,itk,1)
    bc_vel(1,ifac) = rcodcl(ifac,iu,1)
    bc_vel(2,ifac) = rcodcl(ifac,iv,1)
    bc_vel(3,ifac) = rcodcl(ifac,iw,1)

    call cs_cf_thermo_subsonic_outlet_bc(bc_en, bc_pr, bc_vel, ifac-1)

!===============================================================================
! 3.3 Inlet with Ptot, Htot imposed (reservoir boundary conditions)
!===============================================================================

  elseif ( itypfb(ifac).eq.iephcf ) then

    ! If values for Ptot and Htot were not given, the computation stops.

    ! rcodcl(ifac,isca(ienerg),1) contains the boundary total enthalpy values
    ! prescribed by the user

    if(rcodcl(ifac,ipr ,1).lt.-rinfin*0.5d0.or.               &
         rcodcl(ifac,isca(ienerg) ,1).lt.-rinfin*0.5d0) then
      write(nfecra,1200)
      call csexit (1)
    endif

    bc_en(ifac) = rcodcl(ifac,ien,1)
    bc_pr(ifac) = rcodcl(ifac,ipr,1)
    bc_tk(ifac) = rcodcl(ifac,itk,1)
    bc_vel(1,ifac) = rcodcl(ifac,iu,1)
    bc_vel(2,ifac) = rcodcl(ifac,iv,1)
    bc_vel(3,ifac) = rcodcl(ifac,iw,1)

    call cs_cf_thermo_ph_inlet_bc(bc_en, bc_pr, bc_vel, ifac-1)

!===============================================================================
! 3.4 Inlet with imposed rho*U and rho*U*H
!===============================================================================

  elseif ( itypfb(ifac).eq.ieqhcf ) then

    ! TODO to be implemented
    write(nfecra,1301)
    call csexit (1)

    !     On utilise un scenario dans lequel on a un 2-contact et une
    !       3-détente entrant dans le domaine. On détermine les conditions
    !       sur l'interface selon la thermo et on passe dans Rusanov
    !       ensuite pour lisser.

    !     Si rho et u ne sont pas donnés, erreur
    if(rcodcl(ifac,irunh,1).lt.-rinfin*0.5d0) then
      write(nfecra,1300)
      call csexit (1)
    endif

  endif ! end of test on boundary condition types

!===============================================================================
! 4. Complete the treatment for inlets and outlets:
!    - boundary convective fluxes computation (analytical or Rusanov) if needed
!    - B.C. code (Dirichlet or Neumann)
!    - Dirichlet values
!===============================================================================

  if (itypfb(ifac).eq.iesicf.or.                    &
      itypfb(ifac).eq.isspcf.or.                    &
      itypfb(ifac).eq.iephcf.or.                    &
      itypfb(ifac).eq.isopcf.or.                    &
      itypfb(ifac).eq.ieqhcf) then

!===============================================================================
! 4.1 Boundary convective fluxes computation (analytical or Rusanov) if needed
!     (gamma should already have been computed if Rusanov fluxes are computed)
!===============================================================================

    ! Rusanov fluxes are computed only for the imposed inlet for stability
    ! reasons.
    if (itypfb(ifac).eq.iesicf) then

    ! Dirichlet for velocity and pressure are computed in order to
    ! impose the Rusanov fluxes in mass, momentum and energy balance.
      call cfrusb(ifac, bc_en, bc_pr, bc_vel)

    ! For the other types of inlets/outlets (subsonic outlet, QH inlet,
    ! PH inlet), analytical fluxes are computed
    elseif (itypfb(ifac).ne.isspcf) then

      ! the pressure part of the boundary analytical flux is not added here,
      ! but set through the pressure gradient boundary conditions (Dirichlet)
      call cffana(ifac, bc_en, bc_pr, bc_vel)

    endif

!===============================================================================
! 4.2 Copy of boundary values into the Dirichlet values array
!===============================================================================

    if (itypfb(ifac).ne.isspcf) then
      rcodcl(ifac,ien,1) = bc_en(ifac)
      rcodcl(ifac,ipr,1) = bc_pr(ifac)
      rcodcl(ifac,itk,1) = bc_tk(ifac)
      rcodcl(ifac,iu,1)  = bc_vel(1,ifac)
      rcodcl(ifac,iv,1)  = bc_vel(2,ifac)
      rcodcl(ifac,iw,1)  = bc_vel(3,ifac)
      if (ippmod(icompf).eq.2) then ! FIXME fill bc_frac...
        rcodcl(ifac,isca(ifracv),1) = bc_fracv(ifac)
        rcodcl(ifac,isca(ifracm),1) = bc_fracm(ifac)
        rcodcl(ifac,isca(ifrace),1) = bc_frace(ifac)
      endif
    else ! supersonic outlet
      rcodcl(ifac,ien,3) = 0.d0
      rcodcl(ifac,ipr,3) = 0.d0
      rcodcl(ifac,itk,3) = 0.d0
      rcodcl(ifac,iu,3)  = 0.d0
      rcodcl(ifac,iv,3)  = 0.d0
      rcodcl(ifac,iw,3)  = 0.d0
      if (ippmod(icompf).eq.2) then
        rcodcl(ifac,isca(ifracv),3) = 0.d0
        rcodcl(ifac,isca(ifracm),3) = 0.d0
        rcodcl(ifac,isca(ifrace),3) = 0.d0
      endif
    endif

!===============================================================================
! 4.3 Boundary conditions codes (Dirichlet or Neumann)
!===============================================================================

!     P               : Neumann but pressure part of momentum flux is imposed
!                       as a Dirichlet BC for the pressure gradient (code 13)
!     rho, U, E, T    : Dirichlet
!     k, R, eps, scal : Dirichlet/Neumann depending on the flux mass value

    if (itypfb(ifac).ne.isspcf) then
      ! Pressure : - Dirichlet for the gradient computation, allowing to have the
      !              pressure part of the convective flux at the boundary
      !            - Homogeneous Neumann for the diffusion
      icodcl(ifac,ipr)   = 13
      ! velocity
      icodcl(ifac,iu)    = 1
      icodcl(ifac,iv)    = 1
      icodcl(ifac,iw)    = 1
      ! total energy
      icodcl(ifac,ien)   = 1
      ! temperature
      icodcl(ifac,itk)   = 1
      ! mixture fractions
      if (ippmod(icompf).eq.2) then
        icodcl(ifac,isca(ifracv))   = 1
        icodcl(ifac,isca(ifracm))   = 1
        icodcl(ifac,isca(ifrace))   = 1
      endif
    else ! supersonic outlet
      icodcl(ifac,ipr)   = 3
      icodcl(ifac,iu)    = 3
      icodcl(ifac,iv)    = 3
      icodcl(ifac,iw)    = 3
      icodcl(ifac,ien)   = 3
      icodcl(ifac,itk)   = 3
      ! mixture fractions
      if (ippmod(icompf).eq.2) then
        icodcl(ifac,isca(ifracv))   = 3
        icodcl(ifac,isca(ifracm))   = 3
        icodcl(ifac,isca(ifrace))   = 3
      endif
    endif

!-------------------------------------------------------------------------------
! Turbulence and passive scalars: Dirichlet / Neumann depending on the mass flux
!-------------------------------------------------------------------------------

    ! Dirichlet or homogeneous Neumann
    ! A Dirichlet is chosen if the mass flux is ingoing and if the user provided
    ! a value in rcodcl(ifac,ivar,1)

    bmasfl =  brom(ifac)                                                &
             *(  bc_vel(1,ifac)*suffbo(1,ifac)                          &
               + bc_vel(2,ifac)*suffbo(2,ifac)                          &
               + bc_vel(3,ifac)*suffbo(3,ifac) )

    if (itypfb(ifac).ne.isspcf.and.bmasfl.ge.0.d0) then
      if(itytur.eq.2) then
        icodcl(ifac,ik ) = 3
        icodcl(ifac,iep) = 3
      elseif(itytur.eq.3) then
        icodcl(ifac,ir11) = 3
        icodcl(ifac,ir22) = 3
        icodcl(ifac,ir33) = 3
        icodcl(ifac,ir12) = 3
        icodcl(ifac,ir13) = 3
        icodcl(ifac,ir23) = 3
        icodcl(ifac,iep ) = 3
      elseif(iturb.eq.50) then
        icodcl(ifac,ik  ) = 3
        icodcl(ifac,iep ) = 3
        icodcl(ifac,iphi) = 3
        icodcl(ifac,ifb ) = 3
      elseif(iturb.eq.60) then
        icodcl(ifac,ik  ) = 3
        icodcl(ifac,iomg) = 3
      elseif(iturb.eq.70) then
        icodcl(ifac,inusa) = 3
      endif
      if (nscaus.gt.0) then
        do ii = 1, nscaus
          icodcl(ifac,isca(ii)) = 3
        enddo
      endif
      if (nscasp.gt.0) then
        do ii = 1, nscasp
          icodcl(ifac,iscasp(ii)) = 3
        enddo
      endif
    else
      if(itytur.eq.2) then
        if(rcodcl(ifac,ik ,1).lt.rinfin*0.5d0  .and.    &
           rcodcl(ifac,iep,1).lt.rinfin*0.5d0) then
          icodcl(ifac,ik ) = 1
          icodcl(ifac,iep) = 1
        else
          icodcl(ifac,ik ) = 3
          icodcl(ifac,iep) = 3
        endif
      elseif(itytur.eq.3) then
        if(rcodcl(ifac,ir11,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,ir22,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,ir33,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,ir12,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,ir13,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,ir23,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,iep ,1).lt.rinfin*0.5d0) then
          icodcl(ifac,ir11) = 1
          icodcl(ifac,ir22) = 1
          icodcl(ifac,ir33) = 1
          icodcl(ifac,ir12) = 1
          icodcl(ifac,ir13) = 1
          icodcl(ifac,ir23) = 1
          icodcl(ifac,iep ) = 1
        else
          icodcl(ifac,ir11) = 3
          icodcl(ifac,ir22) = 3
          icodcl(ifac,ir33) = 3
          icodcl(ifac,ir12) = 3
          icodcl(ifac,ir13) = 3
          icodcl(ifac,ir23) = 3
          icodcl(ifac,iep ) = 3
        endif
      elseif(iturb.eq.50) then
        if(rcodcl(ifac,ik  ,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,iep ,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,iphi,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,ifb ,1).lt.rinfin*0.5d0) then
          icodcl(ifac,ik  ) = 1
          icodcl(ifac,iep ) = 1
          icodcl(ifac,iphi) = 1
          icodcl(ifac,ifb ) = 1
        else
          icodcl(ifac,ik  ) = 3
          icodcl(ifac,iep ) = 3
          icodcl(ifac,iphi) = 3
          icodcl(ifac,ifb ) = 3
        endif
      elseif(iturb.eq.60) then
        if(rcodcl(ifac,ik  ,1).lt.rinfin*0.5d0  .and.      &
           rcodcl(ifac,iomg,1).lt.rinfin*0.5d0) then
          icodcl(ifac,ik  ) = 1
          icodcl(ifac,iomg) = 1
        else
          icodcl(ifac,ik  ) = 3
          icodcl(ifac,iomg) = 3
        endif
      elseif(iturb.eq.70) then
        if(rcodcl(ifac,inusa,1).gt.0.d0) then
          icodcl(ifac,inusa) = 1
        else
          icodcl(ifac,inusa) = 3
        endif
      endif
      if (nscaus.gt.0) then
        do ii = 1, nscaus
          if(rcodcl(ifac,isca(ii),1).lt.rinfin*0.5d0) then
            icodcl(ifac,isca(ii)) = 1
          else
            icodcl(ifac,isca(ii)) = 3
          endif
        enddo
      endif
      if (nscasp.gt.0) then
        do ii = 1, nscasp
          if(rcodcl(ifac,iscasp(ii),1).lt.rinfin*0.5d0) then
            icodcl(ifac,iscasp(ii)) = 1
          else
            icodcl(ifac,iscasp(ii)) = 3
          endif
        enddo
      endif
    endif

  endif ! end of test on inlet/outlet faces

enddo ! end of loop on boundary faces

! Free memory
deallocate(w5)
deallocate(w7, wbfb, wbfa)
deallocate(bc_en, bc_pr, bc_tk, bc_fracv, bc_fracm, bc_frace, bc_vel)

!----
! FORMATS
!----

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : Error during execution,                       ',/,&
'@    =========                                               ',/,&
'@    two and only two independant variables among            ',/,&
'@    P, rho, T and E have to be imposed at boundaries of type',/,&
'@    iesicf in uscfcl (iccfth = ',i10,').                  ',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Check the boundary condition definitions.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : Error during execution,                       ',/,&
'@    =========                                               ',/,&
'@    The pressure was not provided at outlet with pressure   ',/,&
'@    imposed.                                                ',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Check the boundary conditions in                        ',/,&
'@    cs_user_boundary_conditions                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : Error during execution,                       ',/,&
'@    =========                                               ',/,&
'@    The total pressure or total enthalpy were not provided  ',/,&
'@    at inlet with total pressure and total enthalpy imposed.',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Check the boundary conditions in                        ',/,&
'@    cs_user_boundary_conditions                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : Error during execution,                       ',/,&
'@    =========                                               ',/,&
'@    The mass or enthalpy flow rate were not provided        ',/,&
'@    at inlet with mass and enthalpy flow rate imposed.      ',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Check the boundary conditions in                        ',/,&
'@    cs_user_boundary_conditions                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1301 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : Error during execution,                       ',/,&
'@    =========                                               ',/,&
'@    Inlet with mass and enthalpy flow rate not provided.    ',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Check the boundary conditions in                        ',/,&
'@    cs_user_boundary_conditions                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! END
!----

return
end subroutine
