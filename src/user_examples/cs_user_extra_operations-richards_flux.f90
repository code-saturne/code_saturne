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

!> \file cs_user_extra_operations-richards_flux.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \ref cs_user_extra_operations_examples and
!> \ref cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

subroutine cs_f_user_extra_operations &
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
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
double precision dt(ncelet)

! Local variables

integer iel
integer impout
integer ifac
integer nlelt, ilelt, iflmas, iflmab
integer iel11, iel22
integer iprev, inc, iccocg
integer ifcvsl
integer pas_iter

double precision tra_face, tra_diff, flux_tmp, eps_geom
double precision flux_in, flux_out

integer, allocatable, dimension(:) :: lstelt

double precision, allocatable, dimension(:) :: viscf, viscb, scalar_diffusion

double precision, dimension(:), pointer :: cvar_scal_1
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cpro_vscalt
double precision, dimension(:), pointer :: imasfl, bmasfl

type(var_cal_opt), target   :: vcopt

!===============================================================================

allocate(lstelt(nfac))
allocate(scalar_diffusion(ncelet))
allocate(viscf(nfac))
allocate(viscb(nfabor))

call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_v(ivarfl(iu), cvar_vel)
call field_get_val_s(ivarfl(isca(1)), cvar_scal_1)

!< [richards_flux_comp]
pas_iter = 10 ! Number of iterations separating two flux calculations.

if (mod(ntcabs,pas_iter).eq.0) then

  iprev = 0
  inc = 1
  iccocg = 1
  eps_geom = 1.d-10

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)
  call field_get_key_int(ivarfl(isca(1)), kivisl, ifcvsl)
  call field_get_val_s(ifcvsl, cpro_vscalt)
  scalar_diffusion(1:ncel) = cpro_vscalt(1:ncel)
  call field_get_key_struct_var_cal_opt(ivarfl(isca(1)), vcopt)
  call viscfa (vcopt%imvisf, scalar_diffusion, viscf, viscb)

! Example of tracer flux computation through an internal surface.
! Fluxes are calculated whithout reconstruction, and with a simple definition of the concentration at faces
! (i.e. mean of the two adjacent concentrations).
! Thus, the fluxes are not coherent with the real fluxes calculated in bilsc2 (variable "flux")
! They are approximated and should not be used to verify the exact conservation of mass

  call getfac('PLANE1', nlelt, lstelt)

  flux_in = 0.d0
  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    iel11 = ifacel(1,ifac)
    iel22 = ifacel(2,ifac)

!   Test locality of the cell (only for mpi computation)
    if (iel22.le.ncel) then
      tra_face = 5.d-1*(cvar_scal_1(iel11)+cvar_scal_1(iel22))
      tra_diff = cvar_scal_1(iel11) - cvar_scal_1(iel22)
    else
!     Remove duplicated boundary boundary faces for parallel computation
      tra_face = 0.d0
      tra_diff = 0.d0
    endif

    flux_tmp = imasfl(ifac)*tra_face + viscf(ifac)*tra_diff

!   We impose a norm on the direction of the flux, based on the direction
!   of main coordinates, to be sure that fluxes of all faces are all summed
!   in the same direction
    if ( (surfac(1,ifac).le.(-eps_geom) ) .or.   &
       ( (abs(surfac(1,ifac)).le.eps_geom).and.  &
         (surfac(2,ifac).le.(-eps_geom)) ) .or.  &
       ( (abs(surfac(1,ifac)).le.eps_geom).and.  &
         (abs(surfac(2,ifac)).le.eps_geom).and.  &
         (surfac(3,ifac).le.(-eps_geom)) ) ) then
      flux_tmp = -flux_tmp
    endif
    flux_in = flux_in + flux_tmp
  enddo

! Example of boundary flux computation
! Fluxes are calculated whithout reconstruction, and with a simple definition of the concentration at faces
! (i.e. mean of the two adjacent concentrations).

  call getfac('TOP_SOIL1', nlelt, lstelt)

  flux_out = 0.d0
  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    iel11 = ifacel(1,ifac)
    iel22 = ifacel(2,ifac)

!   Test locality of the cell
    if (iel22.le.ncel) then
      tra_face = 5.d-1*(cvar_scal_1(iel11)+cvar_scal_1(iel22))
      tra_diff = cvar_scal_1(iel11) - cvar_scal_1(iel22)
    else
!     Remove duplicated boundary boundary faces for parallel computation
      tra_face = 0.d0
      tra_diff = 0.d0
    endif
!
    flux_tmp = imasfl(ifac)*tra_face + viscf(ifac)*tra_diff

!   We impose a norm on the direction of the flux, based on the direction
!   of main coordinates, to be sure that fluxes of all faces are all summed
!   in the same direction
    if ( (surfac(1,ifac).le.(-eps_geom) ) .or.   &
       ( (abs(surfac(1,ifac)).le.eps_geom).and.  &
         (surfac(2,ifac).le.(-eps_geom)) ) .or.  &
       ( (abs(surfac(1,ifac)).le.eps_geom).and.  &
         (abs(surfac(2,ifac)).le.eps_geom).and.  &
         (surfac(3,ifac).le.(-eps_geom)) ) ) then
      flux_tmp = -flux_tmp
    endif
    flux_out = flux_out + flux_tmp
  enddo

  ! Compute sum for parallel computations
  if (irangp.ge.0) then
    call parsom(flux_in)
    call parsom(flux_out)
  endif

  ! Write fluxes in file
  impout = impusr(2)
  if (irangp.le.0) then
    open(impout,file='flux_tracer1.dat',position="append")
    write(impout,97) ttcabs, flux_in, flux_out
97      format(3g17.9)
    close(impout)
  endif

endif
!< [richards_flux_comp]

!< [richards_time_modif]
! Example of time modification
do iel = 1, ncel
  dt(iel) = (ttcabs**0.95d0) * 5.d-2
  if (dt(iel).lt.5.d-2) dt(iel) = 5.d-2
  if (dt(iel).gt.1.d3) dt(iel) = 1.d3
enddo
!< [richards_time_modif]

deallocate(lstelt)
deallocate(scalar_diffusion)
deallocate(viscf)
deallocate(viscb)

return
end subroutine cs_f_user_extra_operations
