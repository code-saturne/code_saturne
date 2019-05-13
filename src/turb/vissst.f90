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
!> \file vissst.f90
!> \brief Calculation of turbulent viscosity for
!>        the \f$ k - \omega \f$ SST model
!>
!> \f[ \mu_T = \rho A1 \dfrac{k}{\max(A1 \omega; \; S f_2)} \f]
!> with
!> \f[ S = \sqrt{  2 S_{ij} S_{ij}} \f]
!> \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
!>
!> and \f$ f_2 = \tanh(arg2^2) \f$
!> \f[ arg2^2 = \max(2 \dfrac{\sqrt{k}}{C_\mu \omega y}; \;
!>                   500 \dfrac{\nu}{\omega y^2}) \f]
!> where \f$ y \f$ is the distance to the wall.
!>
!> \f$ \divs{\vect{u}} \f$ is calculated at the same time than \f$ S \f$
!> for being reused in turbkm
!>
!> Edge faces types are available at the previous time step
!> (except at the first time step, when the itypfb and itrifb tables
!>  have not been filled).
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode            name         role
!______________________________________________________________________________!
!______________________________________________________________________________!

subroutine vissst

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use numvar
use optcal
use cstphy
use entsor
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iel, inc
integer          iprev
integer          f_id

double precision d1s3, d2s3
double precision xk, xw, rom, xmu, xdist, xarg2, xf2

double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:), pointer :: cvar_k, cvar_omg
double precision, dimension(:), pointer :: cpro_s2kw, cpro_divukw
double precision, dimension(:), pointer :: w_dist

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
call field_get_val_s(icrom, crom)
call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iomg), cvar_omg)

call field_get_id("wall_distance", f_id)
call field_get_val_s(f_id, w_dist)

d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

!===============================================================================
! 2. Compute the scalar s2kw rate SijSij and the trace of the velocity
!    gradient

!      (Sij^D) (Sij^D)  is stored in    s2kw (deviatoric s2kw tensor rate)
!      tr(Grad u)       is stored in    divukw
!===============================================================================


! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet))

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

! s2kw = Stain rate of the deviatoric part of the s2kw tensor
!      = 2 (Sij^D).(Sij^D)
! divukw   = trace of the velocity gradient
!          = dudx + dvdy + dwdz

call field_get_val_s(is2kw, cpro_s2kw)
call field_get_val_s(idivukw, cpro_divukw)

do iel = 1, ncel

  cpro_s2kw(iel) = 2.d0                                                        &
    *( ( d2s3*gradv(1,1,iel) - d1s3*gradv(2,2,iel) - d1s3*gradv(3,3,iel))**2   &
     + (-d1s3*gradv(1,1,iel) + d2s3*gradv(2,2,iel) - d1s3*gradv(3,3,iel))**2   &
     + (-d1s3*gradv(1,1,iel) - d1s3*gradv(2,2,iel) + d2s3*gradv(3,3,iel))**2   &
     )                                                                         &
    + (gradv(2,1,iel) + gradv(1,2,iel))**2                                     &
    + (gradv(3,1,iel) + gradv(1,3,iel))**2                                     &
    + (gradv(3,2,iel) + gradv(2,3,iel))**2

  cpro_divukw(iel) = gradv(1,1,iel) + gradv(2,2,iel) + gradv(3,3,iel)

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  Calculation of viscosity
!===============================================================================

do iel = 1, ncel

  xk = cvar_k(iel)
  xw = cvar_omg(iel)
  rom = crom(iel)
  xmu = viscl(iel)

  ! Wall distance
  xdist = max(w_dist(iel), epzero)

  ! FIXME should be a check on xw...
  if (xk > 0.d0) then
    ! Wall distance has no value at the first time step, we consider it as infinite
    if (ntcabs.eq.1) then
      xf2 = 0.d0
    else
      xarg2 = max (2.d0*sqrt(xk)/cmu/xw/xdist,                  &
                   500.d0*xmu/rom/xw/xdist**2)
      xf2 = tanh(xarg2**2)
    endif
    visct(iel) =   rom*ckwa1*xk                               &
                 / max(ckwa1*xw, sqrt(cpro_s2kw(iel))*xf2)
  else
    visct(iel) = 1.d-30
  endif

enddo

!-------
! Format
!-------

!----
! End
!----

return
end subroutine
