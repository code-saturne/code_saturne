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

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   maxelt , lstelt ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
!  (nfabor, nphas) !    !     !                                                !
! itypfb           ! ia ! --> ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! izfppp(nfabor)   ! te ! --> ! boundary face zone number                      !
! ia(*)            ! ia ! --- ! main integer work array                        !
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
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! coefu            ! ra ! --- ! work array                                     !
!  (nfabor, 3)     !    !     !  (computation of pressure gradient)            !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
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

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          maxelt, lstelt(maxelt)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, ii, ivar, iphas
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

idebia = idbia0
idebra = idbra0

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
  do iphas = 1, nphas
    itypfb(ifac,iphas) = ientre
  enddo

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

  do iphas = 1, nphas

    itypfb(ifac,iphas) = ientre

    rcodcl(ifac,iu(iphas),1) = xuent
    rcodcl(ifac,iv(iphas),1) = xvent
    rcodcl(ifac,iw(iphas),1) = 0.d0

    ! itytur is a flag equal to iturb/10
    if    (itytur(iphas).eq.2) then

      rcodcl(ifac,ik(iphas),1)  = xkent
      rcodcl(ifac,iep(iphas),1) = xeent

    elseif(itytur(iphas).eq.3) then

      rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir12(iphas),1) = 0.d0
      rcodcl(ifac,ir13(iphas),1) = 0.d0
      rcodcl(ifac,ir23(iphas),1) = 0.d0
      rcodcl(ifac,iep(iphas),1)  = xeent

    elseif(iturb(iphas).eq.50) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iep(iphas),1)  = xeent
      rcodcl(ifac,iphi(iphas),1) = d2s3
      rcodcl(ifac,ifb(iphas),1)  = 0.d0

    elseif(iturb(iphas).eq.60) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

    endif

  enddo

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

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = isolib
  enddo

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

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = iparug

!     Roughness for velocity: rugd
    rcodcl(ifac,iu(iphas),3) = rugd

!     Roughness for scalars (if required):
!   rcodcl(ifac,iv(iphas),3) = rugd


    if(iscalt(iphas).ne.-1) then

    ! If temperature prescribed to 20 with a rough wall law (scalar ii=1)
    ! (with thermal roughness specified in rcodcl(ifac,iv(iphas),3)) :
    ! ii = 1
    ! icodcl(ifac, isca(ii))    = 6
    ! rcodcl(ifac, isca(ii),1)  = 293.15d0

    ! If flux prescribed to 4.d0 (scalar ii=2):
    ! ii = 2
    ! icodcl(ifac, isca(ii))    = 3
    ! rcodcl(ifac, isca(ii), 3) = 4.D0

    endif
  enddo
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

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = isymet
  enddo


enddo

!----
! Formats
!----

!----
! End
!----

return
end subroutine
