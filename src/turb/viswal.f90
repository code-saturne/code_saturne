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
!===============================================================================
! Function:
! ---------

!> \file viswal.f90
!>
!> \brief Compute the turbulent viscosity for the WALE LES model
!>
!> The turbulent viscosity is:
!> \f$ \mu_T = \rho (C_{wale} L)^2 * \dfrac{(\tens{S}:\tens{Sd})^{3/2}}
!>                                         {(\tens{S} :\tens{S})^(5/2)
!>                                         +(\tens{Sd}:\tens{Sd})^(5/4)} \f$
!> with \f$ \tens{S}  = 1/2(\gradt \vect{u} + \transpose{\gradt \vect{u}})\f$
!> and  \f$ \tens{Sd} = \deviator{(\symmetric{(\tens{S}^2)})}\f$
!-------------------------------------------------------------------------------

subroutine viswal

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use dimens, only: nvar
use cstphy
use entsor
use parall
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iel, inc
integer          ipcvst, ipcvis, iprev
integer          i, j, k

double precision coef, delta, third
double precision sij, sijd, s, sd, sinv
double precision con
double precision dudx(ndim,ndim), kdelta(ndim,ndim)

double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: visct

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_val_s(iprpfl(ivisct), visct)
call field_get_val_s(icrom, crom)

third  = 1.d0/3.d0

!===============================================================================
! 2. Computation of the velocity gradient
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

! Kronecker delta Dij

kdelta(1,1) = 1
kdelta(1,2) = 0
kdelta(1,3) = 0
kdelta(2,1) = 0
kdelta(2,2) = 1
kdelta(2,3) = 0
kdelta(3,1) = 0
kdelta(3,2) = 0
kdelta(3,3) = 1

coef = sqrt(2.d0) * cwale**2

do iel = 1, ncel

  ! Dudx is interleaved, but not gradv...
  ! gradv(iel, xyz, uvw)
  dudx(1,1) = gradv(1,1,iel)
  dudx(1,2) = gradv(2,1,iel)
  dudx(1,3) = gradv(3,1,iel)
  dudx(2,1) = gradv(1,2,iel)
  dudx(2,2) = gradv(2,2,iel)
  dudx(2,3) = gradv(3,2,iel)
  dudx(3,1) = gradv(1,3,iel)
  dudx(3,2) = gradv(2,3,iel)
  dudx(3,3) = gradv(3,3,iel)

  s  = 0.d0
  sd = 0.d0

  do i = 1, ndim
    do j = 1, ndim

      ! Sij = 0.5 * (dUi/dXj + dUj/dXi)

      sij = 0.5d0*(dudx(i,j)+dudx(j,i))

      s = s + sij**2

      do k = 1, ndim

        ! traceless symmetric part of the square of the velocity gradient tensor
        !   Sijd = 0.5*( dUi/dXk dUk/dXj + dUj/dXk dUk/dXi)
        !        - 1/3 Dij dUk/dXk dUk/dXk

        sijd = 0.5d0*(dudx(i,k)*dudx(k,j)+ dudx(j,k)*dudx(k,i)) &
              -third*kdelta(i,j)*dudx(k,k)**2

        sd = sd + sijd**2

      enddo
    enddo
  enddo

!===============================================================================
! 3. Computation of turbulent viscosity
!===============================================================================

  ! Turbulent inverse time scale =
  !   (Sijd Sijd)^3/2 / [ (Sij Sij)^5/2 + (Sijd Sijd)^5/4 ]

  sinv = (s**2.5d0 + sd**1.25d0)
  if (sinv.gt.0.d0) then
    con = sd**1.5d0 / sinv
  else
    con = 0.d0
  endif

  delta = xlesfl* (ales*volume(iel))**bles
  delta = coef * delta**2

  visct(iel) = crom(iel) * delta * con

enddo

! Free memory
deallocate(gradv)

!-------
! Format
!-------

!----
! End
!----

return

end subroutine
