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
! Purpose:
! -------

!> \file cs_user_initialization-atmospheric.f90
!>
!> \brief Atmospheric example
!>
!> See \ref cs_user_initialization for examples.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_f_initialization &
!================================

 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
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
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel
double precision d2s3
double precision zent,xuent,xvent,xkent,xeent,tpent

double precision, dimension(:,:), pointer :: cvar_vel

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi, cvar_fb
double precision, dimension(:), pointer :: cvar_omg, cvar_nusa
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: cvar_scalt

!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

!< [init]
! Map field arrays
call field_get_val_v(ivarfl(iu), cvar_vel)

d2s3 = 2.d0/3.d0

!===============================================================================
! Initialize variables using an input meteo profile
!   (only if we are not doing a restart)
!===============================================================================

if (isuite.eq.0) then

  if (itytur.eq.2) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (itytur.eq.3) then
    call field_get_val_v(ivarfl(irij), cvar_rij)

    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (iturb.eq.50) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)
    call field_get_val_s(ivarfl(iphi), cvar_phi)
    call field_get_val_s(ivarfl(ifb), cvar_fb)
  elseif (iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)
  elseif (iturb.eq.70) then
    call field_get_val_s(ivarfl(inusa), cvar_nusa)
  endif

  do iel = 1, ncel

    zent = xyzcen(3,iel)

    call intprf(nbmetd, nbmetm, zdmet, tmmet, umet, zent, ttcabs, xuent)

    call intprf(nbmetd, nbmetm, zdmet, tmmet, vmet, zent, ttcabs, xvent)

    call intprf(nbmetd, nbmetm, zdmet, tmmet, ekmet, zent, ttcabs, xkent)

    call intprf(nbmetd, nbmetm, zdmet, tmmet, epmet, zent, ttcabs, xeent)

    cvar_vel(1,iel) = xuent
    cvar_vel(2,iel) = xvent
    cvar_vel(3,iel) = 0.d0

    ! Initiliation of turbulence variables

    if (itytur.eq.2) then

      cvar_k(iel)  = xkent
      cvar_ep(iel) = xeent

    elseif (itytur.eq.3) then

      cvar_rij(1,iel) = d2s3*xkent
      cvar_rij(2,iel) = d2s3*xkent
      cvar_rij(3,iel) = d2s3*xkent
      cvar_rij(4,iel) = 0.d0
      cvar_rij(5,iel) = 0.d0
      cvar_rij(6,iel) = 0.d0
      cvar_ep(iel)  = xeent

    elseif (iturb.eq.50) then

      cvar_k(iel)   = xkent
      cvar_ep(iel)  = xeent
      cvar_phi(iel) = d2s3
      cvar_fb(iel)  = 0.d0

    elseif (iturb.eq.60) then

      cvar_k(iel)   = xkent
      cvar_omg(iel) = xeent/cmu/xkent

    elseif (iturb.eq.70) then

      cvar_nusa(iel) = cmu*xkent**2/xeent

    endif

    if (iscalt.ge.0) then

      ! Assume the scalar is a potential temperature
      call intprf(nbmett, nbmetm, ztmet, tmmet, tpmet, zent, ttcabs, tpent)

      call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
      cvar_scalt(iel) = tpent

    endif
  enddo

endif
!< [init]

!----
! End
!----

return
end subroutine cs_user_f_initialization
