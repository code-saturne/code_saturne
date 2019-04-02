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
!> \file visqke.f90
!> \brief Calculation of turbulent viscosity for
!>         the non-linear quadratic K-epsilon from
!>         Baglietto et al. (2005)  
!>
!> Edge faces types are vailable at the previous time step
!>  (except at the first time step, when the itypfb anf itrifb tables
!>  have not been filled)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!______________________________________________________________________________!

subroutine visqke

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

integer          iel, inc, f_id
integer          iprev

double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xk, xe, xrom, xmu, xmut, xcmu, xfmu, xss, xrey, xdist
double precision xttke

double precision, allocatable, dimension(:) :: s2
double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:), pointer :: w_dist
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

! --- Memory
allocate(s2(ncelet))

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
call field_get_val_s(icrom, crom)

call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iep), cvar_ep)

call field_get_id("wall_distance", f_id)
call field_get_val_s(f_id, w_dist)

!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S2 = S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet))

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

do iel = 1, ncel

  s11  = gradv(1,1,iel)
  s22  = gradv(2,2,iel)
  s33  = gradv(3,3,iel)
  dudy = gradv(2,1,iel)
  dudz = gradv(3,1,iel)
  dvdx = gradv(1,2,iel)
  dvdz = gradv(3,2,iel)
  dwdx = gradv(1,3,iel)
  dwdy = gradv(2,3,iel)

  s2(iel) = s11**2 + s22**2 + s33**2  &
          + 0.5d0*(dudy+dvdx)**2      &
          + 0.5d0*(dudz+dwdx)**2      &
          + 0.5d0*(dvdz+dwdy)**2

  s2(iel) = max(s2(iel),1.d-10)

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  Calculation of viscosity
!===============================================================================


do iel = 1, ncel
  xk   = cvar_k(iel)
  xe   = cvar_ep(iel)
  xrom = crom(iel)
  xmu  = viscl(iel)
  xdist= max(w_dist(iel),1.d-10)

  xmut = xrom*xk**2/xe
  xrey = xdist*sqrt(xk)*xrom/xmu
  xttke = xk/xe
  xss = xttke*sqrt(.5d0*s2(iel))

  xfmu = 1.d0-exp(-2.9d-2*xrey**.5d0-1.1d-4*xrey**2.d0)
  xcmu = 2.d0/3.d0/(3.9d0 + xss)

  visct(iel) = xcmu*xfmu*xmut
enddo

! Free memory
deallocate(s2)

!----
! Format
!----


!----
! End
!----

return
end subroutine
