!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file chem_source_terms.f90
!> \brief Computes the explicit chemical source term for atmospheric chemistry in
!>        case of a semi-coupled resolution
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar number
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine chem_source_terms &
 ( iscal  ,                  &
   crvexp , crvimp)          &
   bind(C, name='cs_atmo_chem_source_terms')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding
use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use atchem
use sshaerosol

implicit none

!===============================================================================

! Arguments

integer(c_int), value :: iscal
real(c_double) ::  crvexp(ncelet), crvimp(ncelet)

! Local variables

integer iel, ii
double precision dlconc(nespg), source(nespg), dchema(nespg)
double precision rk(nrg)
double precision rom
double precision conv_factor(nespg) ! conversion factors for reaction rates

double precision, dimension(:), pointer :: crom
type(pmapper_double_r1), dimension(:), allocatable :: cvara_espg

!===============================================================================

! If gaseous chemistry is computed with the external library SSH-aerosol
! The source term was estimated at the beginning of scalai
! We update the source term and return directly
if (iaerosol.ne.CS_ATMO_AEROSOL_OFF) then
  ! This is not ready yet
  write(nfecra,*) "Partially coupled chemistry combined with external aerosol library not implemented yet"
  call csexit(1)
  ! c_iscal = iscal - isca_chem(1)
  ! call fexchem_sshaerosol(c_iscal, c_crvexp)
  ! ! FIXME? The negative part could be put in crvimp
  ! do iel = 1, ncel
  !   crvexp(iel) = crvexp(iel) + c_crvexp(iel)
  ! enddo
  ! return
endif

allocate(cvara_espg(nespg))

call field_get_val_s(icrom, crom)

! Arrays of pointers containing the fields values for each species
! (loop on cells outside loop on species)
do ii = 1, nespg
  call field_get_val_prev_s(ivarfl(isca(isca_chem(ii))), cvara_espg(ii)%p)
enddo

do iel = 1, ncel

  ! density
  rom = crom(iel)

  ! Filling working array rk
  do ii = 1, nrg
    rk(ii) = reacnum((ii-1)*ncel+iel)
  enddo

  ! Filling working arrays
  do ii = 1, nespg
    dlconc(chempoint(ii)) = cvara_espg(ii)%p(iel)
    conv_factor(chempoint(ii)) = rom*navo*(1.0d-9)/dmmk(ii)
    source(ii) = 0.0d0
  enddo

  ! Computation of C(Xn)
  if (ichemistry.eq.1) then
    call fexchem_1 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
  else if (ichemistry.eq.2) then
    call fexchem_2 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
  else if (ichemistry.eq.3) then
    call fexchem_3 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
  else if (ichemistry.eq.4) then
    call fexchem_4 (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
  endif

  ! Adding source term to crvexp
  ! The first nespg user scalars are supposed to be chemical species
  ! TODO: try to implicit the ST
  crvexp(iel) = crvexp(iel) + rom * cell_f_vol(iel) &
              * dchema(chempoint(iscal-isca_chem(1)+1))

enddo

deallocate(cvara_espg)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine chem_source_terms
