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
!> \file visv2f.f90
!> \brief Calculation of turbulent viscosity for
!>         the V2F-BL model
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

subroutine visv2f

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

double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xk, xe, xrom, xnu
double precision ttke, ttmin, ttlim, tt

double precision, allocatable, dimension(:) :: s2
double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! --- Memory
allocate(s2(ncelet))

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
call field_get_val_s(icrom, crom)

call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iep), cvar_ep)
call field_get_val_s(ivarfl(iphi), cvar_phi)

!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S2 = 2* (S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
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

  s2(iel) = 2.d0*(s11**2 + s22**2 + s33**2)                   &
       + (dudy+dvdx)**2 + (dudz+dwdx)**2 + (dvdz+dwdy)**2
  s2(iel) = sqrt(max(s2(iel),1.d-10))

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  Calculation of viscosity
!===============================================================================

do iel = 1, ncel

  xk = cvar_k(iel)
  xe = cvar_ep(iel)
  xrom = crom(iel)
  xnu = viscl(iel)/xrom

  ttke = xk / xe
  ttmin = cpalct*sqrt(xnu/xe)
  ttlim = 0.6d0/cvar_phi(iel)/sqrt(3.d0)/cpalmu/s2(iel)
  tt = min(ttlim,sqrt(ttke**2 + ttmin**2))

  visct(iel) = cpalmu*xrom*tt*cvar_phi(iel)*cvar_k(iel)

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
