!-------------------------------------------------------------------------------

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
! Function:
! ---------
!> \file rijthe2.f90
!> \brief Gravity terms
!>        For \f$R_{ij}\f$
!>        \f[ var = R_{11} \: R_{22} \: R_{33} \:R_{12} \:R_{23} \:R_{13}\f]

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     gradro        work array for \f$ \grad{\rho} \f$
!> \param[in,out] buoyancy      Buoyancy term for the Reynolds stress model
!______________________________________________________________________________!

subroutine rijthe2 &
 ( gradro , buoyancy )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision gradro(3,ncelet)
double precision buoyancy(6,ncelet)

! Local variables

integer          iel, isou, i, j
double precision dij

double precision uns3, const, kseps
double precision rit(3), gij(3, 3), grav(3)
double precision gkks3
double precision turb_schmidt

double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:,:), pointer :: cvara_rij

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

if (iscalt.gt.0) then
  call field_get_key_double(ivarfl(isca(iscalt)), ksigmas, turb_schmidt)
  const = -1.5d0 * cmu / turb_schmidt
else
  const = -1.5d0 * cmu
endif

uns3  = 1.d0/3.d0

grav(1) = gx
grav(2) = gy
grav(3) = gz

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_prev_v(ivarfl(irij), cvara_rij)

!===============================================================================
! 2. Terms for Rij:
!      rom*volume*dRij/dt =
!                     ... + (Gij - crij3*(Gij-Delta ij Gkk/3))*volume
!            With Gij = -(1.5 cmu/PrT) (k/eps) (Rit Gj + Rjt Gi)
!                 Rit = Rik drom/dxk (sum on k)
!            Note that in tensorial notation Gij is:
!                 G = a [R.( Grho x g) + (g x Grho).R]
!===============================================================================
do iel = 1, ncel
  rit(1) = cvara_rij(1,iel)*gradro(1,iel)                            &
         + cvara_rij(4,iel)*gradro(2,iel)                            &
         + cvara_rij(6,iel)*gradro(3,iel)
  rit(2) = cvara_rij(4,iel)*gradro(1,iel)                            &
         + cvara_rij(2,iel)*gradro(2,iel)                            &
         + cvara_rij(5,iel)*gradro(3,iel)
  rit(3) = cvara_rij(6,iel)*gradro(1,iel)                            &
         + cvara_rij(5,iel)*gradro(2,iel)                            &
         + cvara_rij(3,iel)*gradro(3,iel)

  kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

  do i = 1, 3
    do j = 1, 3
      gij(i,j) = const*kseps* (rit(i) * grav(j) + rit(j) * grav(i))
    enddo
  enddo

  gkks3 = uns3*(gij(1,1) + gij(2,2) + gij(3,3))

  do isou = 1, 6

    if     (isou.eq.1) then
      i = 1
      j = 1
      dij = 1.0d0
    elseif (isou.eq.2) then
      i = 2
      j = 2
      dij = 1.0d0
    elseif (isou.eq.3) then
      i = 3
      j = 3
      dij = 1.0d0
    elseif (isou.eq.4) then
      i = 1
      j = 2
      dij = 0.0d0
    elseif (isou.eq.5) then
      i = 2
      j = 3
      dij = 0.0d0
    elseif (isou.eq.6) then
      i = 1
      j = 3
      dij = 0.0d0
    endif

    buoyancy(isou,iel) = gij(i,j) * (1.d0 - crij3) + crij3 * dij * gkks3

  enddo
enddo

return

end subroutine rijthe2
