!-------------------------------------------------------------------------------

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
!> \file cs_at_source_term.f90
!
!> \brief Additional right-hand side source terms for momentum equation in case
!>        of free inlet
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   crvexp          explicit part of the momentum source term
!-------------------------------------------------------------------------------

subroutine cs_at_source_term_for_inlet ( crvexp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use mesh
use atincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision crvexp(3, ncelet)

! Local variables

integer              iel, vel_id
integer              level_id, n_level, id
integer              iz1, iz2
double precision  :: xent, yent
double precision  :: xuent, xvent, zent, dpdtx_ent, dpdty_ent
double precision  :: dist_ent, dist_min
double precision  :: mom_met_norm, mom_norm, mom_met_norm_a, mom_norm_a
double precision, dimension(3) :: dir_met, dir_met_a

double precision, dimension(:), pointer :: dt
double precision, allocatable, dimension (:,:) :: mom_met_a, mom_a
double precision, allocatable, dimension (:) :: tot_vol, dpdtx, dpdty
double precision, dimension(:), pointer :: crom
double precision, dimension(:,:), pointer :: vel, cpro_momst, cpro_vel_target
double precision, allocatable, dimension (:,:), target :: wvel_target

!===============================================================================
! 1. Initialisation
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_id_try('meteo_velocity', vel_id)

if (vel_id.gt.0) then
  call field_get_val_v_by_name('meteo_velocity', cpro_vel_target)
else
  allocate(wvel_target(3, ncelet))
  cpro_vel_target => wvel_target(1:3, 1:ncelet)
endif

call field_get_val_s_by_name('dt', dt)

call field_get_val_v(imomst, cpro_momst)

! --- Density
call field_get_val_s(icrom, crom)

! Bulk momentum
if (iatmst.eq.1) then
  n_level = 1

! Variable in z
else
  n_level = nbmetd
endif

allocate(tot_vol(n_level))
allocate(dpdtx(n_level))
allocate(dpdty(n_level))
allocate(mom_a(3, n_level))
allocate(mom_met_a(3, n_level))

! Save previous values
mom_a = mom
mom_met_a = mom_met ! target meteo value

do level_id = 1, n_level
  mom(1, level_id) = 0.d0
  mom(2, level_id) = 0.d0
  mom(3, level_id) = 0.d0
  mom_met(1, level_id) = 0.d0
  mom_met(2, level_id) = 0.d0
  mom_met(3, level_id) = 0.d0
  tot_vol(level_id) = 0.d0
enddo

!===============================================================================
! 2. Computation of the target momentum bulk using the interpolated
!    mean velocity field
!===============================================================================

do iel = 1, ncel

  ! Get level id (nearest)
  zent = xyzcen(3,iel)
  dist_min = abs(zent - zdmet(1))
  level_id = 1
  do id = 2, n_level
    dist_ent = abs(zent -zdmet(id))
    if (dist_ent.lt.dist_min) then
      level_id = id
      dist_min = dist_ent
    endif
  enddo

  if (theo_interp.eq.1) then

    xuent = cpro_vel_target(1,iel)
    xvent = cpro_vel_target(2,iel)

  else

    call intprf &
    (nbmetd, nbmetm,                                               &
      zdmet, tmmet, umet , zent  , ttcabs, xuent )

    call intprf &
    (nbmetd, nbmetm,                                               &
      zdmet, tmmet, vmet , zent  , ttcabs, xvent )

    cpro_vel_target(1,iel) = xuent
    cpro_vel_target(2,iel) = xvent
    cpro_vel_target(3,iel) = 0.d0

  endif

  mom_met(1, level_id) = mom_met(1, level_id) + crom(iel) * cell_f_vol(iel) * xuent
  mom_met(2, level_id) = mom_met(2, level_id) + crom(iel) * cell_f_vol(iel) * xvent
  tot_vol(level_id) = tot_vol(level_id) + cell_f_vol(iel)

enddo

if (irangp.ge.0) then
  call parrsm(3 * n_level,mom_met)
  call parrsm(n_level, tot_vol)
endif

do level_id = 1, n_level
  if (tot_vol(level_id).gt.0.d0) then
    mom_met(1, level_id) = mom_met(1, level_id) / tot_vol(level_id)
    mom_met(2, level_id) = mom_met(2, level_id) / tot_vol(level_id)
    mom_met(3, level_id) = mom_met(3, level_id) / tot_vol(level_id)
  endif
enddo

!===============================================================================
! 3. Computation of the momentum using the computed velocity
!===============================================================================

do iel = 1, ncel

  ! Get level id (nearest)
  zent = xyzcen(3,iel)
  dist_min = abs(zent - zdmet(1))
  level_id = 1
  do id = 2, n_level
    dist_ent = abs(zent -zdmet(id))
    if (dist_ent.lt.dist_min) then
      level_id = id
      dist_min = dist_ent
    endif
  enddo

  mom(1, level_id) = mom(1, level_id) + crom(iel) * cell_f_vol(iel)*vel(1,iel)
  mom(2, level_id) = mom(2, level_id) + crom(iel) * cell_f_vol(iel)*vel(2,iel)
  mom(3, level_id) = mom(3, level_id) + crom(iel) * cell_f_vol(iel)*vel(3,iel)
enddo

if (irangp.ge.0) then
  call parrsm(3 * n_level, mom)
endif

do level_id = 1, n_level
  if (tot_vol(level_id).gt.0.d0) then
    mom(1, level_id) = mom(1, level_id) / tot_vol(level_id)
    mom(2, level_id) = mom(2, level_id) / tot_vol(level_id)
    mom(3, level_id) = mom(3, level_id) / tot_vol(level_id)
  endif
enddo

!===============================================================================
! 4. Computation of the momentum source term
!===============================================================================

! First pass, reset previous values
if (ntcabs.le.1.or.ntcabs.eq.ntpabs+1) then
  do level_id = 1, n_level
    mom_met_a(1, level_id) = mom_met(1, level_id)
    mom_met_a(2, level_id) = mom_met(2, level_id)
    mom_met_a(3, level_id) = mom_met(3, level_id)
    mom_a(1, level_id) = mom(1, level_id)
    mom_a(2, level_id) = mom(2, level_id)
    mom_a(3, level_id) = mom(3, level_id)
    dpdt_met(level_id) = 0.d0
  enddo

  ! Initialisation of the momentum source term with an initial momentum balance
!  call cs_balance_vector(idtvar, ivarfl(iu), imasac, inc, ivisse, )

endif

! Delta of pressure integrated over a time step for each level
do level_id = 1, n_level

  ! Momentum of CS and of the target
  mom_norm = sqrt(mom(1, level_id)**2 + mom(2, level_id)**2 + mom(3, level_id)**2)

  mom_norm_a = &
    sqrt(mom_a(1, level_id)**2 + mom_a(2, level_id)**2 + mom_a(3, level_id)**2)

  mom_met_norm = &
    sqrt(mom_met(1, level_id)**2 &
    + mom_met(2, level_id)**2 &
    + mom_met(3, level_id)**2)

  mom_met_norm_a = &
    sqrt(mom_met_a(1, level_id)**2 &
    + mom_met_a(2, level_id)**2 &
    + mom_met_a(3, level_id)**2)

  dpdt_met(level_id) = dpdt_met(level_id) &
    + 0.5d0*(2.d0*(mom_norm - mom_met_norm) &
    - (mom_norm_a - mom_met_norm_a))

  ! target meteo directions (current and previous)
  dir_met(1) = 0.d0
  dir_met(2) = 0.d0
  dir_met(3) = 0.d0
  if (mom_met_norm.gt. epzero * uref * ro0) then
    dir_met(1) = mom_met(1, level_id) / mom_met_norm
    dir_met(2) = mom_met(2, level_id) / mom_met_norm
    dir_met(3) = mom_met(3, level_id) / mom_met_norm
  endif
  dir_met_a(1) = 0.d0
  dir_met_a(2) = 0.d0
  dir_met_a(3) = 0.d0
  if (mom_met_norm_a.gt. epzero * uref * ro0) then
    dir_met_a(1) = mom_met_a(1, level_id) / mom_met_norm_a
    dir_met_a(2) = mom_met_a(2, level_id) / mom_met_norm_a
    dir_met_a(3) = mom_met_a(3, level_id) / mom_met_norm_a
  endif

  ! Delta of pressure in the target direction
  ! Steady state DP and Rotating term due to transient meteo
  dpdtx(level_id) = dpdt_met(level_id) * dir_met(1) &
    - mom_met_norm*(dir_met(1)-dir_met_a(1))!FIXME use directly umet?
  dpdty(level_id) = dpdt_met(level_id) * dir_met(2) &
    - mom_met_norm*(dir_met(2)-dir_met_a(2))

enddo

do iel = 1, ncel

  xent = xyzcen(1,iel)
  yent = xyzcen(2,iel)
  ! Z Interpolation of dpdtx and dpdty
  zent = xyzcen(3,iel)
  if (n_level.gt.1) then
    call intprz &
      (nbmetd, zdmet, &
      dpdtx, zent  , &
      iz1, iz2, &
      dpdtx_ent)
    call intprz &
      (nbmetd, zdmet, &
      dpdty, zent  , &
      iz1, iz2, &
      dpdty_ent)
  else
    dpdtx_ent = dpdtx(1)
    dpdty_ent = dpdty(1)
  endif

  cpro_momst(1,iel) = - dpdtx_ent / dt(iel)
  cpro_momst(2,iel) = - dpdty_ent / dt(iel)
  !FIXME not ok
  if (iz1.ne.iz2.and..false.) then
    cpro_momst(3,iel) = &
      + (dpdtx(iz2) - dpdtx(iz1)) / (zdmet(iz2) - zdmet(iz1)) * xent / dt(iel) &
      + (dpdty(iz2) - dpdty(iz1)) / (zdmet(iz2) - zdmet(iz1)) * yent / dt(iel)
  else
    cpro_momst(3,iel) = 0.d0
  endif

  crvexp(1,iel) = crvexp(1,iel) + cpro_momst(1,iel) * cell_f_vol(iel)
  crvexp(2,iel) = crvexp(2,iel) + cpro_momst(2,iel) * cell_f_vol(iel)
  crvexp(3,iel) = crvexp(3,iel) + cpro_momst(3,iel) * cell_f_vol(iel)
enddo

if (allocated(wvel_target)) deallocate(wvel_target)
deallocate(mom_met_a, mom_a, tot_vol)
deallocate(dpdtx, dpdty)

return

end subroutine cs_at_source_term_for_inlet
