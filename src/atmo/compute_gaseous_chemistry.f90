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

!===============================================================================
!  Purpose:
!  --------

!> \file compute_gaseous_chemistry.f90
!> \brief Calls the rosenbrock resolution for atmospheric chemistry
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine compute_gaseous_chemistry ( dt )

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
use pointe
use ppppar
use ppthch
use ppincl
use mesh
use field
use dimens
use atchem
use siream

implicit none

!===============================================================================

! Arguments

double precision dt(ncelet)

! Local Variables

integer iel,ii
double precision rom

!  Variables used for Rosenbrock resolution
double precision  dlconc(nespg)
double precision  source(nespg)
double precision  dchema(nespg)
double precision  conv_factor(nespg) ! conversion factors for reaction rates
double precision rk(nrg)
double precision dtc
integer ncycle
double precision dtrest

double precision, dimension(:), pointer :: crom
type(pmapper_double_r1), dimension(:), allocatable :: cvar_espg, cvara_espg

!===============================================================================

allocate(cvar_espg(nespg), cvara_espg(nespg))

call field_get_val_s(icrom, crom)

! Arrays of pointers containing the fields values for each species
! (loop on cells outside loop on species)
do ii = 1, nespg
  call field_get_val_s(ivarfl(isca(isca_chem(ii))), cvar_espg(ii)%p)
  call field_get_val_prev_s(ivarfl(isca(isca_chem(ii))), cvara_espg(ii)%p)
enddo

do iel = 1, ncel

  ! time step
  dtc = dt(iel)

  ! density
  rom = crom(iel)

  ! Filling working array rk
  do ii = 1, nrg
    rk(ii) = reacnum((ii-1)*ncel+iel)
  enddo

  do ii = 1, nespg
    conv_factor(chempoint(ii)) = rom*navo*(1.0d-9)/dmmk(ii)
    source(ii) = 0.0d0
  enddo

  if ((isepchemistry.eq.1).or.(ntcabs.eq.1)) then
    ! -----------------------------
    ! -- splitted Rosenbrock solver
    ! -----------------------------

    ! Filling working array dlconc with values at current time step
    do ii = 1, nespg
      dlconc(chempoint(ii)) = cvar_espg(ii)%p(iel)
    enddo

  else
    ! -----------------------------
    ! -- semi-coupled Rosenbrock solver
    ! -----------------------------

    ! Filling working array dlconc with values at previous time step
    do ii = 1, nespg
      dlconc(chempoint(ii)) = cvara_espg(ii)%p(iel)
    enddo

    ! Computation of C(Xn)
    if (ichemistry.eq.1) then
      call fexchem_1 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
    else if (ichemistry.eq.2) then
      call fexchem_2 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
    else if (ichemistry.eq.3) then
      if (iaerosol.eq.1) then
        call fexchem_siream (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
      else
        call fexchem_3 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
      endif
    else if (ichemistry.eq.4) then
      call fexchem_4 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
    endif

    ! Explicit contribution from dynamics as a source term:
    ! (X*-Xn)/dt(dynamics) - C(Xn). See usatch.f90
    ! The first nespg user scalars are supposed to be chemical species
    do ii = 1, nespg
      source(chempoint(ii)) = (cvar_espg(ii)%p(iel)-cvara_espg(ii)%p(iel))/dtc  &
                            - dchema(chempoint(ii))
    enddo

  endif ! End test isepchemistry

  ! Rosenbrock resoluion

  ! The maximum time step used for chemistry resolution is dtchemmax
  if (dtc.le.dtchemmax) then
    call roschem (dlconc,source,source,conv_factor,dtc,rk,rk)
  else
    ncycle = int(dtc/dtchemmax)
    dtrest = mod(dtc,dtchemmax)
    do ii = 1, ncycle
      call roschem (dlconc,source,source,conv_factor,dtchemmax,rk,rk)
    enddo
    call roschem (dlconc,source,source,conv_factor,dtrest,rk,rk)
  endif

  ! Update of values at current time step
  do ii = 1, nespg
    cvar_espg(ii)%p(iel) = dlconc(chempoint(ii))
  enddo

enddo

deallocate(cvar_espg, cvara_espg)

! Clipping
do ii = 1, nespg
  call clpsca(isca_chem(ii))
enddo

!--------
! FORMATS
!--------
!----
! END
!----

return
end subroutine compute_gaseous_chemistry
