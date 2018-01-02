!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
! Purpose:
! -------

!> \file cs_user_extra_operations-extract_1d_profile.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!> This is an example of cs_user_extra_operations.f90 which
!> performs 1D profile.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel
integer          iel1
integer          impout
integer          ii     , irangv , irang1 , npoint
integer          iun

double precision xyz(3), xabs, xu, xv, xw, xk, xeps

double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_omg
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33

!< [loc_var_dec]

!===============================================================================

!< [example_1]
if (ntcabs.eq.ntmabs) then

  ! Map field arrays
  call field_get_val_v(ivarfl(iu), cvar_vel)

  if (itytur.eq.2 .or. iturb.eq.50    &
       .or. iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
  elseif (itytur.eq.3) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
  endif
  if (     itytur.eq.2 .or. itytur.eq.3    &
       .or. iturb.eq.50) then
    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)
  endif

  ! Only process of rank 0 (parallel) or -1 (scalar) writes to this file.
  ! We use 'user' Fortran units.
  impout = impusr(1)
  if (irangp.le.0) then
    open(impout,file='profile.dat')
    write(impout,*)  &
         '# z(m) U(m/s) V(m/s) W(m/s) k(m2/s2) eps(m2/s3)'
  endif

  npoint = 200
  iel1   = -999
  irang1 = -999
  do ii = 1, npoint

    xyz(1) = 0.d0
    xyz(2) = float(ii-1)/float(npoint-1)*0.1d0
    xyz(3) = 0.d0

    call findpt(ncelet, ncel, xyzcen, xyz(1), xyz(2), xyz(3), iel, irangv)
    !==========

    if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
      iel1   = iel
      irang1 = irangv

      ! Set temporary variables xu, xv, ... for the process containing
      ! the point and then send it to other processes.
      if (irangp.eq.irangv) then
        xabs = xyzcen(2,iel)
        xu   = cvar_vel(1,iel)
        xv   = cvar_vel(2,iel)
        xw   = cvar_vel(3,iel)
        xk   = 0.d0
        xeps = 0.d0
        if (     itytur.eq.2 .or. iturb.eq.50    &
            .or. iturb.eq.60) then
          xk = cvar_k(iel)
        elseif (itytur.eq.3) then
          xk = (  cvar_r11(iel) + cvar_r22(iel)          &
                + cvar_r33(iel)) / 2.d0
        endif
        if (     itytur.eq.2 .or. itytur.eq.3    &
            .or. iturb.eq.50) then
          xeps = cvar_ep(iel)
        elseif (iturb.eq.60) then
          xeps = cmu*cvar_k(iel)*cvar_omg(iel)
        endif
      else
        xabs = 0.d0
        xu   = 0.d0
        xv   = 0.d0
        xw   = 0.d0
        xk   = 0.d0
        xeps = 0.d0
      endif

      ! Broadcast to other ranks in parallel
      if (irangp.ge.0) then
        call parall_bcast_r(irangv, xabs)
        call parall_bcast_r(irangv, xu)
        call parall_bcast_r(irangv, xv)
        call parall_bcast_r(irangv, xw)
        call parall_bcast_r(irangv, xk)
        call parall_bcast_r(irangv, xeps)
      endif

      if (irangp.le.0) write(impout,99) xabs, xu, xv, xw, xk, xeps

99    format(6g17.9)

    endif

  enddo

  if (irangp.le.0) close(impout)

endif
!< [example_1]

return
end subroutine cs_f_user_extra_operations
