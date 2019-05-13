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
! Function:
! --------
!> \file vissma.f90
!> \brief Calculation of turbulent viscosity for
!>        a Smagorinsky LES model
!>
!> \f[ \mu_T = \rho (C_{S} l)^2  \sqrt{2 S_{ij}S_{ij}} \f]
!> \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
!>
!> Edge faces types are available at the previous time step
!> (except at the first time step, when the itypfb and itrifb
!> have not been filled).
!>
!> Please refer to the
!> <a href="../../theory.pdf#smago"><b>standard Smagorinsky model</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!______________________________________________________________________________!

subroutine vissma

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

double precision coef, delta
double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xfil, xa  , xb

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

! --- For the calculation of viscosity on the sub-mesh
xfil   = xlesfl
xa     = ales
xb     = bles

!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

do iel = 1, ncel

  s11  = gradv(1, 1, iel)
  s22  = gradv(2, 2, iel)
  s33  = gradv(3, 3, iel)
  dudy = gradv(2, 1, iel)
  dvdx = gradv(1, 2, iel)
  dudz = gradv(3, 1, iel)
  dwdx = gradv(1, 3, iel)
  dvdz = gradv(3, 2, iel)
  dwdy = gradv(2, 3, iel)

  visct(iel) = s11**2 + s22**2 + s33**2               &
                     + 0.5d0*((dudy+dvdx)**2          &
                     +        (dudz+dwdx)**2          &
                     +        (dvdz+dwdy)**2)
enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  Calcualtion of (dynamic) velocity
!===============================================================================

coef = csmago**2 * sqrt(2.d0)

do iel = 1, ncel
  delta  = xfil* (xa*volume(iel))**xb
  delta  = coef * delta**2
  visct(iel) = crom(iel) * delta * sqrt(visct(iel))
enddo

!----
! Format
!----


!----
! End
!----

return
end subroutine
