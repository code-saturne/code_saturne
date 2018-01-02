!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
! --------
!> \file cs_ctwr_bcond.f90
!>
!> \brief Automatic boundary condition for cooling towers
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[in]     izfppp        zone number for the boundary face for
!>                                      the specific physic module
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
!> \param[in,out] rcodcl        value of the boundary conditions to edge faces
!>
!>                              boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                               -  coefficient (infinite if no exchange)
!>                               -  rcodcl(3) value flux density
!>                               -  (negative if gain) \f$w.m^{-2} \f$ or
!>                               -  roughness in \f$m\f$ if  icodcl=6
!>                                -# for velocity:
!>                                           \f$(\mu+\mu_T)\gradv \vect{u}\f$
!>                                -# for pressure: \f$ \Delta \grad P
!>                                                 \cdot \vect{n} \f$
!>                                -# for scalar:   \f$ C_p \left ( K +
!>                                                 \dfrac{K_T}{\sigma_T} \right)
!>                                                 \grad T \cdot \vect{n} \f$
!______________________________________________________________________________!

subroutine cs_ctwr_bcond &
 ( itypfb , izfppp ,                                              &
   icodcl , rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only : nvar
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use ctincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)
integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, izone, iel
integer          icke

double precision viscla, uref2, rhomoy, dhy, xiturb
double precision humidity, t_h
double precision t_l, h_l
double precision, dimension(:), pointer :: brom

!===============================================================================
! 0. Initializations
!===============================================================================

call field_get_val_s(ibrom, brom)

!===============================================================================
! 1.  Parallel exchanges for the user data
!===============================================================================

!===============================================================================
! 2.  Filling the table of the boundary conditions
!       Loop on all input faces
!      =========================
!         Determining the family and its properties
!         Imposing boundary conditions for the turbulence

!===============================================================================

do ifac = 1, nfabor

  izone = izfppp(ifac)

  if (itypfb(ifac).eq.ientre .or.itypfb(ifac).eq.ifrent) then

    !       The turbulence is calculated by default if icalke different from 0
    !          - or from hydraulic diameter and a reference velocity adapted
    !            for the current input if icalke = 1
    !          - either from the hydraulic diameter, a reference velocity and
    !            a turbulence intensity adapted to the current input if icalke = 2
    if (icalke(izone).ne.0) then

      uref2 = rcodcl(ifac,iu,1)**2                         &
            + rcodcl(ifac,iv,1)**2                         &
            + rcodcl(ifac,iw,1)**2
      uref2 = max(uref2,1.d-12)
      rhomoy = brom(ifac)
      iel    = ifabor(ifac)
      viscla = viscl0
      icke   = icalke(izone)
      dhy    = dh(izone)
      xiturb = xintur(izone)

      if (icke.eq.1) then
        !   Calculation of turbulent inlet conditions using
        !     standard laws for a circular pipe
        !     (their initialization is not needed here but is good practice).
        call turbulence_bc_inlet_hyd_diam(ifac, uref2, dhy, rhomoy, viscla,  &
                                          rcodcl)
      else if (icke.eq.2) then

        ! Calculation of turbulent inlet conditions using
        !   the turbulence intensity and standard laws for a circular pipe
        !   (their initialization is not needed here but is good practice)

        call turbulence_bc_inlet_turb_intensity(ifac, uref2, xiturb, dhy,  &
                                                rcodcl)

      endif

    endif

    ! -- Boundary conditions for the transported temperature
    ! -- of the humid air and of the liquid water injected in the packing zones
    ! Bulk values if not set by the user
    ! Assume that the humid air is at conditions '0'
    humidity = humidity0
    rhomoy = brom(ifac)
    if (icodcl(ifac, isca(iscalt)).eq.0) then
       t_h = t0-tkelvi
       icodcl(ifac, isca(iscalt)) = 1
       rcodcl(ifac, isca(iscalt), 1) = t_h
    endif
    if (icodcl(ifac, isca(iymw)).eq.0) then
      icodcl(ifac, isca(iymw)) = 1
      rcodcl(ifac, isca(iymw), 1) = humidity / (1.d0 + humidity)
    endif

    ! For injected liquid
    if (icodcl(ifac, isca(iyml)).eq.0) then
      icodcl(ifac, isca(iyml)) = 1
      rcodcl(ifac, isca(iyml), 1) = 0.d0
    endif
    if (icodcl(ifac, isca(ihml)).eq.0) then
       t_l = t0-tkelvi
       call h_liqwater(t_l, h_l)

       ! Y_l . h_l is transported (and not h_l)
       h_l = h_l*rcodcl(ifac, isca(iyml), 1)

       icodcl(ifac, isca(ihml)) = 1
       rcodcl(ifac, isca(ihml), 1) = h_l
    endif

    ! Wall BCs: zero Flux
  else if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then

    icodcl(ifac, isca(iscalt)) = 3
    rcodcl(ifac, isca(iscalt), 3) = 0.d0
    icodcl(ifac, isca(iymw)) = 3
    rcodcl(ifac, isca(iymw), 3) = 0.d0

    icodcl(ifac, isca(ihml)) = 3
    rcodcl(ifac, isca(ihml), 3) = 0.d0
    icodcl(ifac, isca(iyml)) = 3
    rcodcl(ifac, isca(iyml), 3) = 0.d0

    icodcl(ifac, isca(iy_p_l)) = 1
    rcodcl(ifac, isca(iy_p_l), 1) = 0.d0

  endif

enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_ctwr_bcond
