!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine usatcl &
!================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )

!===============================================================================
! Purpose:
! -------

!    User subroutine for the atmospheric module.

!    Fill boundary conditions arrays (icodcl, rcodcl) for unknown variables.

!    (similar to usclim.f90)


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Boundary condition types
! ========================

! Boundary conditions setup for standard variables (pressure, velocity,
! turbulence, scalars) is described precisely in the 'usclim' subroutine.

! Detailed explanation will be found in the theory guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb           ! ia ! --> ! boundary face types                            !
! izfppp(nfabor)   ! te ! --> ! boundary face zone number                      !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!  (nfabor,nvar,3) !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ihmpre
use ppppar
use atincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iel, ii, ivar
integer          izone
integer          ilelt, nlelt
double precision uref2, d2s3
double precision rhomoy, dh, ustar2
double precision xintur
double precision zref,xuref
double precision ustar,rugd, rugt
double precision zent,xuent,xvent
double precision xkent, xeent
double precision tpent

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in definition of boundary conditions',/,   &
'@    =======',/,                                                 &
'@      for the atmospheric module                          ',/,  &
'@     The user subroutine ''usatcl'' must be completed.',/,      &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END



!===============================================================================
! 1.  Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


d2s3 = 2.d0/3.d0

! Paremeters for the analytical rough wall law (neutral)
zref=10.d0
xuref=10.d0
rugd=0.1d0
rugt=0.1d0

!===============================================================================
! 2.  Assign boundary conditions to boundary faces here

!     One may use selection criteria to filter boundary case subsets
!       Loop on faces from a subset
!         Set the boundary condition for each face
!===============================================================================

! --- For boundary faces of color 11,
!       assign an inlet boundary condition for all phases prescribed from the meteo profile
!       with automatic choice between inlet/ outlet according to the meteo profile

call getfbr('11',nlelt,lstelt)
!==========
!   - Zone number (from 1 to n)
izone = 1

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     - Zone to which the face belongs
  izfppp(ifac) = izone

!     - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

enddo


! ---For boundary faces of color 21,
!     assign an inlet boundary condition for all phases prescribed from the meteo profile

call getfbr('21',nlelt,lstelt)
!==========
!   -  Zone number (from 1 to n)
izone = 2

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     - Zone to which the face belongs
  izfppp(ifac) = izone

!     - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

!     - Assign inlet boundary conditions
  itypfb(ifac) = ientre

enddo

! --- For boundary faces of color 31,
!       assign an inlet boundary condition for all phases prescribed from the meteo profile
!       except for dynamical variables which are prescribed with a rough log law.

call getfbr('31',nlelt,lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 3

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     - Zone to which the face belongs
  izfppp(ifac) = izone

!     - Boundary conditions are prescribed from the meteo profile
  iprofm(izone) = 1

!     - Dynamical variables are prescribed with a rough log law
  zent=cdgfbo(3,ifac)

  ustar=xkappa*xuref/log((zref+rugd)/rugd)
  xuent=ustar/xkappa*log((zent+rugd)/rugd)
  xvent = 0.d0
  xkent=ustar**2/sqrt(cmu)
  xeent=ustar**3/xkappa/(zent+rugd)

  itypfb(ifac) = ientre

  rcodcl(ifac,iu,1) = xuent
  rcodcl(ifac,iv,1) = xvent
  rcodcl(ifac,iw,1) = 0.d0

  ! itytur is a flag equal to iturb/10
  if    (itytur.eq.2) then

    rcodcl(ifac,ik,1)  = xkent
    rcodcl(ifac,iep,1) = xeent

  elseif(itytur.eq.3) then

    rcodcl(ifac,ir11,1) = d2s3*xkent
    rcodcl(ifac,ir22,1) = d2s3*xkent
    rcodcl(ifac,ir33,1) = d2s3*xkent
    rcodcl(ifac,ir12,1) = 0.d0
    rcodcl(ifac,ir13,1) = 0.d0
    rcodcl(ifac,ir23,1) = 0.d0
    rcodcl(ifac,iep,1)  = xeent

  elseif(iturb.eq.50) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iep,1)  = xeent
    rcodcl(ifac,iphi,1) = d2s3
    rcodcl(ifac,ifb,1)  = 0.d0

  elseif(iturb.eq.60) then

    rcodcl(ifac,ik,1)   = xkent
    rcodcl(ifac,iomg,1) = xeent/cmu/xkent

  elseif(iturb.eq.70) then

    rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

  endif

enddo

! --- Prescribe at boundary faces of color '12' an outlet for all phases
call getfbr('12', nlelt, lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 4

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     - Zone to which the zone belongs
  izfppp(ifac) = izone

  ! Outlet: zero flux for velocity and temperature, prescribed pressure
  !         Note that the pressure will be set to P0 at the first
  !         free outlet face (isolib)

  itypfb(ifac)   = isolib

enddo

! --- Prescribe at boundary faces of color 15 a rough wall for all phases
call getfbr('15', nlelt, lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 5

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       rough friction for velocities (+ turbulent variables)
  !       zero flux for scalars

!     - Zone to which the zone belongs
  izfppp(ifac) = izone

  itypfb(ifac)   = iparug

  !     Roughness for velocity: rugd
  rcodcl(ifac,iu,3) = rugd

  !     Roughness for scalars (if required):
  !   rcodcl(ifac,iv,3) = rugd


  if(iscalt.ne.-1) then

    ! If temperature prescribed to 20 with a rough wall law (scalar ii=1)
    ! (with thermal roughness specified in rcodcl(ifac,iv,3)) :
    ! ii = 1
    ! icodcl(ifac, isca(ii))    = 6
    ! rcodcl(ifac, isca(ii),1)  = 293.15d0

    ! If flux prescribed to 4.d0 (scalar ii=2):
    ! ii = 2
    ! icodcl(ifac, isca(ii))    = 3
    ! rcodcl(ifac, isca(ii), 3) = 4.D0

  endif
enddo

! --- Prescribe at boundary faces of color 4 a symmetry for all phases
call getfbr('4', nlelt, lstelt)
!==========

!   - Zone number (from 1 to n)
izone = 6

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     - Zone to which the zone belongs
  izfppp(ifac) = izone

  itypfb(ifac)   = isymet

enddo

!----
! Formats
!----

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
