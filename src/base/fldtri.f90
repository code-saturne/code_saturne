!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2012 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine fldtri &
!================

 ( nproce ,                                                            &
   dt     , tpucou , rtpa   , rtp    , propce , propfa , propfb , ra )

!===============================================================================
! Purpose:
! --------

! Map fields

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nproce           ! i  ! <-- ! nombre de prop phy aux centres                 !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! tpucou(ncelet,3) ! ra ! <-- ! velocity-pressure coupling                     !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! ra(*)            ! ra ! <-- ! main work array                                !
!__________________.____._____.________________________________________________.

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use optcal
use cstphy
use numvar
use entsor
use pointe
use albase
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use lagpar
use lagdim
use lagran
use ihmpre
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nproce
double precision dt(ncelet), tpucou(ncelet,3), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision frcxt(ncelet,3)
double precision ra(*)

! Local variables

integer          ii, ippu, ippv, ippw, ivar, iprop
integer          imom, idtnm
integer          iflid, nfld

integer          ifvar(nvppmx)

character*80     name
character*32     name1, name2, name3

!===============================================================================


!===============================================================================
! 1. Initialisation
!===============================================================================

nfld = 0

!===============================================================================
! 2. Mapping for post-processing
!===============================================================================

! Velocity and pressure
!----------------------

ivar = ipr
call fldmap(ivarfl(ivar), rtp(1,ivar), rtpa(1,ivar))
!==========

ivar = iu
call fldmap(ivarfl(ivar), rtp(1,ivar), rtpa(1,ivar))
!==========

! Turbulence
!-----------

if (itytur.eq.2) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iep
elseif (itytur.eq.3) then
  nfld = nfld + 1
  ifvar(nfld) = ir11
  nfld = nfld + 1
  ifvar(nfld) = ir22
  nfld = nfld + 1
  ifvar(nfld) = ir33
  nfld = nfld + 1
  ifvar(nfld) = ir12
  nfld = nfld + 1
  ifvar(nfld) = ir13
  nfld = nfld + 1
  ifvar(nfld) = ir23
  nfld = nfld + 1
  ifvar(nfld) = iep
elseif (itytur.eq.5) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iep
  nfld = nfld + 1
  ifvar(nfld) = iphi
  if (iturb.eq.50) then
    nfld = nfld + 1
    ifvar(nfld) = ifb
  elseif (iturb.eq.51) then
    nfld = nfld + 1
    ifvar(nfld) = ial
  endif
elseif (iturb.eq.60) then
  nfld = nfld + 1
  ifvar(nfld) = ik
  nfld = nfld + 1
  ifvar(nfld) = iomg
elseif (iturb.eq.70) then
  nfld = nfld + 1
  ifvar(nfld) = inusa
endif

! Map fields

do ii = 1, nfld
  ivar = ifvar(ii)
  call fldmap(ivarfl(ivar), rtp(1,ivar), rtpa(1,ivar))
  !==========
enddo

nfld = 0

! Mesh velocity
!--------------

if (iale.eq.1) then
  ivar = iuma
  call fldmap(ivarfl(ivar), rtp(1,ivar), rtpa(1,ivar))
  !==========
endif

! User variables
!---------------

do ii = 1, nscaus
  if (isca(ii) .gt. 0) then
    ivar = isca(ii)
    call fldmap(ivarfl(ivar), rtp(1,ivar), rtpa(1,ivar))
    !==========
  endif
enddo

! The choice made in VARPOS specifies that we will only be interested in
! properties at cell centers (no mass flux, nor density at the boundary).

do iprop = 1, nproce
  call fldmap(iprpfl(iprop), propce(1, iprop), propce(1, iprop))
  !==========
enddo

! Reserved fields whose ids are not saved (may be queried by name)
!-----------------------------------------------------------------

! Local time step

name = 'dt'
call fldfid(name, iflid)
!==========
call fldmap(iflid, dt, dt)
!==========

! Transient velocity/pressure coupling

if (ipucou.ne.0) then
  name = 'tpucou'
  call fldfid(name, iflid)
  !==========
  call fldmap(iflid, tpucou, tpucou)
  !==========
endif

return
end subroutine
