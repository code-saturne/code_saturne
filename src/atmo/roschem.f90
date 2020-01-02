!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file roschem.f90
!> \brief Rosenbrock solver for atmospheric chemistry
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in,out] dlconc        concentrations vector
!> \param[in]     zcsourc       source term for first iteration
!> \param[in]     zcsourcf      source term for second iteration
!> \param[in]     conv_factor   conversion factor
!> \param[in]     dlstep        time step
!> \param[in]     dlrki         kinetic rates for first iteration
!> \param[in]     dlrkf         kinetic rates for second iteration
!______________________________________________________________________________

subroutine roschem (dlconc,zcsourc,zcsourcf,conv_factor,                      &
                    dlstep,dlrki,dlrkf)

!===============================================================================
! Module files
!===============================================================================

use atchem
use siream

implicit none

! Arguments

double precision dlconc(nespg)
double precision zcsourc(nespg),zcsourcf(nespg)
double precision conv_factor(nespg)
double precision dlstep
double precision dlrki(nrg),dlrkf(nrg)


! Local variables

integer ji, jj
double precision dlconcbis(nespg)
! chemical production terms for every species
double precision dlr(nespg)
double precision igamma
double precision dlk1(nespg), dlk2(nespg)
double precision dlb1(nespg), dlb2(nespg)
double precision dlmat(nespg,nespg)
double precision dlmatlu(nespg,nespg)
! jacobian matrix
double precision dldrdc(nespg,nespg)

!------------------------------------------------------------------------
!*    0. Setup

igamma = 1.d0 + 1.d0/dsqrt(2.d0)

!------------------------------------------------------------------------
!*    1. Computes the chemistry

if (ichemistry.eq.1) then
  call fexchem_1 (nespg,nrg,dlconc,dlrki,zcsourc,conv_factor,dlr)
else if (ichemistry.eq.2) then
  call fexchem_2 (nespg,nrg,dlconc,dlrki,zcsourc,conv_factor,dlr)
else if (ichemistry.eq.3) then
  if (iaerosol.eq.1) then
    call fexchem_siream (nespg,nrg,dlconc,dlrki,zcsourc,conv_factor,dlr)
  else
    call fexchem_3 (nespg,nrg,dlconc,dlrki,zcsourc,conv_factor,dlr)
  endif
else if (ichemistry.eq.4) then
  call fexchem_4 (nespg,nrg,dlconc,dlrki,zcsourc,conv_factor,dlr)
endif

!------------------------------------------------------------------------
!*    2. Compute the jacobian

if (ichemistry.eq.1) then
  call jacdchemdc_1 (nespg,nrg,dlconc,conv_factor,conv_factor_jac,dlrki,dldrdc)
else if (ichemistry.eq.2) then
  call jacdchemdc_2 (nespg,nrg,dlconc,conv_factor,conv_factor_jac,dlrki,dldrdc)
else if (ichemistry.eq.3) then
  if (iaerosol.eq.1) then
    call jacdchemdc_siream (nespg,nrg,dlconc,conv_factor,conv_factor_jac,       &
                            dlrki,dldrdc)
  else
    call jacdchemdc_3 (nespg,nrg,dlconc,conv_factor,conv_factor_jac,dlrki,dldrdc)
  endif
else if (ichemistry.eq.4) then
  call jacdchemdc (nespg,nrg,dlconc,conv_factor,conv_factor_jac,dlrki,dldrdc)
endif

!------------------------------------------------------------------------
!*    4. Computes K1    system: DLmat * K1 = DLb1

do ji = 1, nespg
  dlb1(ji) = dlr(ji)
  do jj = 1, nespg
    dlmat(ji,jj) = -igamma*dlstep*dldrdc(ji,jj)
  enddo
  dlmat(ji,ji) = 1.d0 + dlmat(ji,ji)
enddo

call solvlin (0,dlmat,dlmatlu,dlk1,dlb1)

!------------------------------------------------------------------------
!*    5. Computes K2    system: DLmat * K2 = DLb2

do ji = 1, nespg
  dlconcbis(ji) = dlconc(ji) + dlstep * dlk1(ji)
  if (dlconcbis(ji) .lt. 0.d0) then
    dlconcbis(ji) = 0.d0
    dlk1(ji) = (dlconcbis(ji) - dlconc(ji)) / dlstep
  endif
enddo

if (ichemistry.eq.1) then
  call fexchem_1 (nespg,nrg,dlconcbis,dlrkf,zcsourcf,conv_factor,dlr)
else if (ichemistry.eq.2) then
  call fexchem_2 (nespg,nrg,dlconcbis,dlrkf,zcsourcf,conv_factor,dlr)
else if (ichemistry.eq.3) then
  if (iaerosol.eq.1) then
    call fexchem_siream (nespg,nrg,dlconcbis,dlrkf,zcsourcf,conv_factor,dlr)
  else
    call fexchem_3 (nespg,nrg,dlconcbis,dlrkf,zcsourcf,conv_factor,dlr)
  endif
else if (ichemistry.eq.4) then
  call fexchem_4 (nespg,nrg,dlconcbis,dlrkf,zcsourcf,conv_factor,dlr)
endif

do ji = 1, nespg
  dlb2(ji) =  dlr(ji) - 2.d0*dlk1(ji)
enddo

call solvlin (1,dlmat,dlmatlu,dlk2,dlb2)

!------------------------------------------------------------------------
!*    6. Outputs - Compute DLconc - Advance the time

do ji = 1, nespg
  dlconc(ji) = dlconc(ji) + 1.5d0 * dlstep * dlk1(ji)         &
             + 0.5d0 * dlstep * dlk2(ji)

  if (dlconc(ji) .lt. 0.0d0) then
    dlconc(ji) = 0.d0
  endif
enddo

return
end subroutine roschem
