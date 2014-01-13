!-------------------------------------------------------------------------------

!VERS

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
! Purpose:
! -------

! This is an example of cs_user_extra_operations.f90 which
! performs 1D profile.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     , rtpa   , rtp    , propce )

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
use lagpar
use lagran
use lagdim
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)

! Local variables

!< [loc_var_dec]
integer          iel
integer          iel1
integer          impout
integer          ii     , irangv , irang1 , npoint
integer          iun

double precision xyz(3), xabs, xu, xv, xw, xk, xeps
!< [loc_var_dec]

!===============================================================================

!< [example_1]
if (ntcabs.eq.ntmabs) then

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
        xu   = rtp(iel,iu)
        xv   = rtp(iel,iv)
        xw   = rtp(iel,iw)
        xk   = 0.d0
        xeps = 0.d0
        if (     itytur.eq.2 .or. iturb.eq.50    &
            .or. iturb.eq.60) then
          xk = rtp(iel,ik)
        elseif (itytur.eq.3) then
          xk = (  rtp(iel,ir11) + rtp(iel,ir22)  &
                + rtp(iel,ir33)) / 2.d0
        endif
        if (     itytur.eq.2 .or. itytur.eq.3    &
            .or. iturb.eq.50) then
          xeps = rtp(iel,iep)
        elseif (iturb.eq.60) then
          xeps = cmu*rtp(iel,ik)*rtp(iel,iomg)
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
        iun = 1
        call parbcr(irangv, iun, xabs)
        call parbcr(irangv, iun, xu)
        call parbcr(irangv, iun, xv)
        call parbcr(irangv, iun, xw)
        call parbcr(irangv, iun, xk)
        call parbcr(irangv, iun, xeps)
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
