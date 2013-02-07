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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_extra_operations-extract_1d_profile.f90
!>
!> \brief This is an example of cs_user_extra_operations.f90 which
!>  performs 1D profile.
!>
!>
!> \section loc_var Local variables to be added
!>
!> \snippet cs_user_extra_operations-extract_1d_profile.f90 loc_var_dec
!>
!>
!> \subsection ex_1 Example 1
!>
!>  We seek here to extract the profile of U, V, W, k and epsilon on an
!>  arbitrary 1D curve based on a curvilear abscissa.
!>  The profile is described in the 'profile.dat' file (do not forget to
!>  define it as user data in the run script).
!>
!>  - the curve used here is the segment: [(0;0;0),(0;0.1;0)], but the
!>    generalization to an arbitrary curve is simple.
!>  - the routine handles parallelism an periodicity, as well as the different
!>    turbulence models.
!>  - the 1D curve is discretized into 'npoint' points. For each of these
!>    points, we search for the closest cell center and we output the variable
!>    values at this cell center. For better consistency, the coordinate
!>    which is output is that of the cell center (instead of the initial point).
!>  - we avoid using the same cell multiple times (in case several points
!>    an the curve are associated with the same cell).
!>
!> Here is the corresponding code:
!>
!> \snippet cs_user_extra_operations-extract_1d_profile.f90 example_1
!>
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     nbpmax        max. number of particles allowed
!> \param[in]     nvp           number of particle-defined variables
!> \param[in]     nvep          number of real particle properties
!> \param[in]     nivep         number of integer particle properties
!> \param[in]     ntersl        number of return coupling source terms
!> \param[in]     nvlsta        number of Lagrangian statistical variables
!> \param[in]     nvisbr        number of boundary statistics
!> \param[in]     itepa         integer particle attributes
!>                                (containing cell, ...)
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     ettp, ettpa   particle-defined variables
!> \param[in]                    (at current and previous time steps)
!> \param[in]     tepa          real particle properties
!> \param[in]                    (statistical weight, ...
!> \param[in]     statis        statistic means
!> \param[in]     stativ        accumulator for variance of volume statisitics
!> \param[in]     tslagr        Lagrangian return coupling term
!> \param[in]                    on carrier phase
!> \param[in]     parbor        particle interaction properties
!> \param[in]                    on boundary faces
!_______________________________________________________________________________


subroutine cs_user_extra_operations &
!==================================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   itepa  ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor )

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
end subroutine cs_user_extra_operations
