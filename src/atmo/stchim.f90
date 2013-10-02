!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file stchim.f90
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
!> \param[in]     rtpa          calculated variables at cell centers
!>                               (preceding time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine stchim &
!================

 ( iscal  , rtpa   ,  propce ,                            &
   crvexp , crvimp)


!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use atchem

implicit none

!===============================================================================

! Arguments

integer          iscal
double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

!------------------------------------------------------------------------------------
!  Variables used for computation of the explicit chemical source term
integer iel, ii
double precision dlconc(nespg), source(nespg), dchema(nespg)
double precision rk(nrg)
double precision rom
double precision conv_factor(nespg) ! conversion factors for reaction rates

!===============================================================================

do iel = 1, ncel

  ! density
  rom = propce(iel,ipproc(irom))

  ! Filling working array rk
  do ii = 1, nrg
    rk(ii) = reacnum((ii-1)*ncel+iel)
  enddo

  ! Filling working arrays
  do ii = 1, nespg
    dlconc(chempoint(ii)) = rtpa(iel,isca(ii))  ! rtpa
    conv_factor(chempoint(ii)) = rom*navo*(1.0d-12)/dmmk(ii)
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
    call fexchem (nespg,nrg,dlconc,rk,source,conv_factor,dchema)
  endif

  ! Adding source term to crvexp
  ! The first nespg user scalars are supposed to be chemical species
  ! TODO: try to implicit the ST
  crvexp(iel) = crvexp(iel)+dchema(chempoint(iscal))*rom*volume(iel)

enddo

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine stchim
