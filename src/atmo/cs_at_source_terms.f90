!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
double precision  :: xuent, xvent, zent, tot_vol
double precision  :: norm_bulk, norm_mom, norm_bulk_a, norm_mom_a
double precision, save :: dp
double precision, dimension(3), save :: mom_a, mom_bulk_a, dir_a
double precision, dimension(3) :: mom_bulk, mom, dir, dir_qdm, dir_var
double precision, dimension(:), pointer :: crom
double precision, dimension(:,:), pointer :: vel, cpro_momst, met_vel

!===============================================================================
! 1. Initialisation
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_id_try('meteo_velocity', vel_id)
if (vel_id.gt.0) call field_get_val_v_by_name('meteo_velocity', met_vel)
call field_get_val_v(imomst, cpro_momst)

! --- Density
call field_get_val_s(icrom, crom)

! Computation of the total fluid volume
tot_vol = 0.d0

do iel = 1, ncel
  tot_vol = tot_vol + volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(tot_vol)
endif

!===============================================================================
! 2. Computation of the momentum bulk using the interpolated mean velocity field
!===============================================================================

mom_bulk = 0.d0

do iel = 1, ncel

  zent = xyzcen(3,iel)

  if (vel_id.gt.0.and.theo_interp.eq.1) then

    xuent = met_vel(1,iel)
    xvent = met_vel(2,iel)

  else

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdmet, tmmet, umet , zent  , ttcabs, xuent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdmet, tmmet, vmet , zent  , ttcabs, xvent )

    if (vel_id.gt.0) then

      met_vel(1, iel) = xuent
      met_vel(2, iel) = xvent

    endif
  endif

  mom_bulk(1) = mom_bulk(1) + crom(iel)*volume(iel)*xuent/tot_vol
  mom_bulk(2) = mom_bulk(2) + crom(iel)*volume(iel)*xvent/tot_vol

  if (vel_id.gt.0) then
    met_vel(1,iel) = xuent
    met_vel(2,iel) = xvent
  endif

enddo

if (irangp.ge.0) then
  call parrsm(3,mom_bulk)
endif

norm_bulk = sqrt(mom_bulk(1)**2 + mom_bulk(2)**2 + mom_bulk(3)**2)

!===============================================================================
! 3. Computation of the momentum using the computed velocity
!===============================================================================

mom = 0.d0

do iel = 1, ncel
  mom(1) = mom(1) + crom(iel)*volume(iel)*vel(1,iel)/tot_vol
  mom(2) = mom(2) + crom(iel)*volume(iel)*vel(2,iel)/tot_vol
  mom(3) = mom(3) + crom(iel)*volume(iel)*vel(3,iel)/tot_vol
enddo

if (irangp.ge.0) then
  call parrsm(3,mom)
endif

norm_mom = sqrt(mom(1)**2 + mom(2)**2 + mom(3)**2)

!===============================================================================
! 4. Computation of the direction vectors
!===============================================================================

dir = 0.d0

if (norm_bulk.gt.1.d-12*uref) then
  dir(1) = mom_bulk(1)/norm_bulk
  dir(2) = mom_bulk(2)/norm_bulk
  dir(3) = mom_bulk(3)/norm_bulk
endif

dir_qdm = 0.d0

if (norm_mom.gt.1.d-12*uref) then
  dir_qdm(1) = mom(1)/norm_mom
  dir_qdm(2) = mom(2)/norm_mom
  dir_qdm(3) = mom(3)/norm_mom
endif

!===============================================================================
! 5. Computation of the momentum source term
!===============================================================================

if (ntcabs.eq.1.or.ntcabs.eq.ntpabs+1) then
  mom_bulk_a = mom_bulk
  mom_a = mom
  dp = 0.d0
  dir_a = dir
endif

norm_mom_a = sqrt(mom_a(1)**2 + mom_a(2)**2 + mom_a(3)**2)
norm_bulk_a = sqrt(mom_bulk_a(1)**2 + mom_bulk_a(2)**2 + mom_bulk_a(3)**2)

dp = dp + 0.5d0*(2.d0*(norm_mom - norm_bulk) &
                     -(norm_mom_a - norm_bulk_a))/dtref

dir_var(1) = norm_bulk*(dir(1)-dir_a(1))/dtref
dir_var(2) = norm_bulk*(dir(2)-dir_a(2))/dtref
dir_var(3) = 0.d0

do iel = 1, ncel
  crvexp(1,iel) = -(dp*dir_qdm(1)-dir_var(1))*volume(iel)
  crvexp(2,iel) = -(dp*dir_qdm(2)-dir_var(2))*volume(iel)
  crvexp(3,iel) = -(dp*dir_qdm(3)-dir_var(3))*volume(iel)
  cpro_momst(1,iel) = -(dp*dir_qdm(1)-dir_var(1))
  cpro_momst(2,iel) = -(dp*dir_qdm(2)-dir_var(2))
  cpro_momst(3,iel) = -(dp*dir_qdm(3)-dir_var(3))
enddo

mom_a = mom
mom_bulk_a = mom_bulk
dir_a = dir

return

end subroutine cs_at_source_term_for_inlet
