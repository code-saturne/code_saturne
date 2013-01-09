!-------------------------------------------------------------------------------

!VERS

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

subroutine cs_user_extra_operations &
!==================================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   itepa  ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Called at end of each time step, very general purpose
!    (i.e. anything that does not have another dedicated user subroutine)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! max. number of particles allowed               !
! nvp              ! i  ! <-- ! number of particle-defined variables           !
! nvep             ! i  ! <-- ! number of real particle properties             !
! nivep            ! i  ! <-- ! number of integer particle properties          !
! ntersl           ! i  ! <-- ! number of return coupling source terms         !
! nvlsta           ! i  ! <-- ! number of Lagrangian statistical variables     !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! itepa            ! ia ! <-- ! integer particle attributes                    !
!  (nbpmax, nivep) !    !     !   (containing cell, ...)                       !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! ettp, ettpa      ! ra ! <-- ! particle-defined variables                     !
!  (nbpmax, nvp)   !    !     !  (at current and previous time steps)          !
! tepa             ! ra ! <-- ! real particle properties                       !
!  (nbpmax, nvep)  !    !     !  (statistical weight, ...                      !
! statis           ! ra ! <-- ! statistic means                                !
!  (ncelet, nvlsta)!    !     !                                                !
! stativ(ncelet,   ! ra ! <-- ! accumulator for variance of volume statisitics !
!        nvlsta -1)!    !     !                                                !
! tslagr           ! ra ! <-- ! Lagrangian return coupling term                !
!  (ncelet, ntersl)!    !     !  on carrier phase                              !
! parbor           ! ra ! <-- ! particle interaction properties                !
!  (nfabor, nvisbr)!    !     !  on boundary faces                             !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta), stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)


! Local variables

integer          iel
integer          iel1
integer          impout
integer          ii     , irangv , irang1 , npoint
integer          iun

double precision xyz(3), xabs, xu, xv, xw, xk, xeps

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!===============================================================================
! Example: extraction of a 1D profile
! -----------------------------------

! We seek here to extract the profile of U, V, W, k and epsilon on an
! arbitrary 1D curve based on a curvilear abscissa.
! The profile is described in the 'profile.dat' file (do not forget to
! define it as user data in the run script).

! - the curve used here is the segment: [(0;0;0),(0;0.1;0)], but the
!   generalization to an arbitrary curve is simple.
! - the routine handles parallelism an periodicity, as well as the different
!   turbulence models.
! - the 1D curve is discretized into 'npoint' points. For each of these
!   points, we search for the closest cell center and we output the variable
!   values at this cell center. For better consistency, the coordinate
!   which is output is that of the cell center (instead of the initial point).
! - we avoid using the same cell multiple times (in case several points
!   an the curve are associated with the same cell).
!===============================================================================

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

return
end subroutine cs_user_extra_operations
