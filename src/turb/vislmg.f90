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

!===============================================================================
! Function:
! --------
!> \file vislmg.f90
!> \brief Calculation of turbulent viscosity for
!>        a model of length of simple mixture
!>
!> \f[ \mu_T = \rho (\kappa L)^2 \cdot \sqrt{2 S_{ij} S_{ij}} \f]
!> \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
!>
!> Edge face types are available at previous time step (except at the first time
!> step, when the itypfb and itrifb tables have not been filled).
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!______________________________________________________________________________!

subroutine vislmg

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iel, inc
integer          iprev

double precision coef, deux

double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: visct

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

call field_get_val_s(ivisct, visct)
call field_get_val_s(icrom, crom)

!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iu), iprev, 0, inc, gradv)

do iel = 1, ncel
  visct(iel) = &
      gradv(1, 1, iel)**2 + gradv(2, 2, iel)**2 + gradv(3, 3, iel)**2  &
    + 0.5d0*( (gradv(2, 1, iel) + gradv(1, 2, iel))**2                 &
            + (gradv(3, 1, iel) + gradv(1, 3, iel))**2                 &
            + (gradv(3, 2, iel) + gradv(2, 3, iel))**2 )
enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  Calculation of (dynamic) velocity
!===============================================================================

deux = 2.d0
coef = (xkappa*xlomlg)**2 * sqrt(deux)

do iel = 1, ncel
  visct(iel) = crom(iel) * coef * sqrt(visct(iel))
enddo

!----
! Format
!----


!----
! End
!----

return
end subroutine
