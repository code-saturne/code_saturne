!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file compute_siream.f90
!> \brief Computation of atmospheric aerosol chemistry
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     dt             time step (per cell)
!> \param[in,out] rtp            calculated variables at cell centers at current time step
!> \param[in]     propce         physical properties at cell centers
!______________________________________________________________________________

subroutine compute_siream &
!================

 ( dt     , rtp,   propce)

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
use ppppar
use ppthch
use ppincl
use elincl
use mesh
use field
use dimens
use atincl
use atchem
use siream

implicit none

!===============================================================================

! Arguments

double precision dt(ncelet), rtp(ncelet,*)
double precision propce(ncelet,*)

! Local Variables

integer ii, jb, jsp, iel
double precision rvoid(1)
double precision dlconc(nespg_siream)
double precision dlconc_aer(nbin_aer,nesp_aer)
double precision dlnum_aer(nbin_aer)
double precision dltemp, dlpress, dlhumid, dldens
double precision zent                               ! Z coordinate of a cell

double precision, dimension(:), pointer :: crom

! Initialisation

dltemp = t0
dldens = ro0
dlpress = dldens*rair*dltemp ! ideal gas law
dlhumid = 0.0d0

! Get Physical quantities
! Densisty
call field_get_val_s(icrom, crom)

! Clipping before calling siream

do ii = 1, nespg_siream+nesp_aer*nbin_aer
  call clpsca(ncelet, ncel, ii, rvoid, rtp)
enddo

! Loop on cells
do iel = 1, ncel

  zent = xyzcen(3,iel)
  ! Verifier l'ordre dans lequel on rentre rtp, dlconc et dlconc_aer

  ! Filling working arrays
  do jsp = 1, nespg_siream
    dlconc(jsp) = rtp(iel,isca(jsp))
  enddo

  do jb = 1, nbin_aer
    do jsp = 1, nesp_aer
      dlconc_aer(jb,jsp) = rtp(iel,isca(nespg_siream+jb+(jsp-1)*nbin_aer))
    enddo
  enddo

  do jb = 1, nbin_aer
    dlnum_aer(jb) = rtp(iel,isca(nespg_siream+nesp_aer*nbin_aer+jb))
  enddo

  ! Temperature and density
  ! Dry or humid atmosphere
  if (ippmod(iatmos).ge.1) then
    dltemp = propce(iel,ipproc(itempc)) + tkelvi
    dldens = crom(iel)
    dlpress = dldens*rair*dltemp

  ! Constant density with a meteo file
  else if (imeteo.eq.1) then

    ! Hydrostatic pressure
    call intprf                                                   &
    !==========
   (nbmett, nbmetm,                                               &
    ztmet , tmmet, phmet, zent, ttcabs, dlpress)

    ! Temperature
    call intprf                                                   &
    !==========
   (nbmett, nbmetm,                                               &
    ztmet , tmmet, ttmet, zent, ttcabs, dltemp )
    dltemp = dltemp + tkelvi
  endif

  ! Specific hymidity
  ! Humid atmosphere
  if (ippmod(iatmos).ge.2) then
    dlhumid = (rtp(iel,isca(itotwt))-propce(iel,ipproc(iliqwt)))   &
              /(1.d0-propce(iel,ipproc(iliqwt)))

  ! Constant density or dry atmosphere with a meteo file
  else if (imeteo.eq.1) then

    call intprf                                                     &
    !==========
         (nbmett, nbmetm,                                           &
         ztmet , tmmet, qvmet, zent, ttcabs, dlhumid )

  endif

  ! Siream options
  options_aer(1) = icoag_siream
  options_aer(2) = icond_siream
  options_aer(3) = inucl_siream
  options_aer(4) = 0 ! Useful for aqueous chemistry only. Not yet implemented
  options_aer(5) = 0 ! Useful for aqueous chemistry only. Not yet implemented
  options_aer(6) = ikelv_siream
  options_aer(7) = 0 ! Fixed density
  options_aer(8) = icut_siream
  options_aer(9) = isulfcond_siream
  options_aer(10) = kdslv_siream
  options_aer(11) = iredist_siream
  options_aer(12) = itern_siream
  options_aer(13) = ithrm_siream
  options_aer(14) = 0 ! Forced use of isorropia as the thermodynamics module

  ! Calling Siream
  call plug_aerosol(1, 1, 1, nespg_siream, 0.0d0, dlhumid,              &
       dltemp, dlpress, dt(iel), dlconc, noptions_aer,                  &
       options_aer, nesp_aer, nbin_aer, ncycle_aer,                     &
       bin_bound_aer, fixed_density_aer, density_aer,                   &
       couples_coag, first_index_coag, second_index_coag,               &
       coefficient_coag, dlconc_aer, dlnum_aer)

  ! Updating rtp
  do jsp = 1, nespg_siream
    rtp(iel,isca(jsp)) = dlconc(jsp)
  enddo

  do jb = 1, nbin_aer
    do jsp = 1, nesp_aer
      rtp(iel,isca(nespg_siream+jb+(jsp-1)*nbin_aer)) = dlconc_aer(jb,jsp)
    enddo
  enddo

  do jb = 1, nbin_aer
    rtp(iel,isca(nespg_siream+nesp_aer*nbin_aer+jb)) = dlnum_aer(jb)
  enddo

enddo ! Loop on cells

! Clipping after calling Siream

do ii = 1, nespg_siream+nesp_aer*nbin_aer
  call clpsca(ncelet, ncel, ii, rvoid, rtp)
enddo

!--------
! FORMATS
!--------
!----
! END
!----

return
end subroutine compute_siream
