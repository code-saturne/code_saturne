!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine uslain &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nptnew ,                                                       &
   itypfb , itrifb , itepa  , ifrlag , injfac ,                   &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)

! User subroutine for the boundary conditions for the particles
! (inlet and treatment for the other boundaries)
!
! This routine is called after the initialization of the
! ettp, tepa and itepa arrays for the new particles in order to modify them
! to inject new particle profiles.
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! nptnew           ! i  ! <-- ! total number of new particles for all the      !
!                  !    !     ! injection zones                                !
! itrifb(nfabor)   ! ia ! <-- ! indirection for the sorting of the boundary    !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! ifrlag(nfabor)   ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! injfac(npbmax)   ! ia ! <-- ! number of the injection boundary face          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous timestep !
! (ncelet,*)       !    !     !                                                !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! vagaus           ! ra ! --> ! Gaussian random variables                      !
!(nbpmax,nvgaus    !    !     !                                                !
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
use cstnum
use cstphy
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nptnew

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)
integer          injfac(nbpmax)

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*)

! Local variables

integer          iclas , izone , ifac
integer          ii , ip , npt , npar1 , npar2, ipnorm

! User-defined local variables
! (the dimension of vgauss is 3, but 2 would be sufficient here)

double precision vgauss(3)

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

!     By default, we do not modify them

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

if (nbpnew.eq.0) return

!===============================================================================
! 1. Memory management
!===============================================================================


!===============================================================================
! 2. Initializations
!===============================================================================



!===============================================================================
! 3. Modification of properties of the new particles (injection profiles,
!    position of the injection point, statistical
!    weights, correction of the diameter if the standard-deviation option
!    is activated.)
!===============================================================================

!    These modifications occur after all the initializations related to
!    the particle injection, but before the treatment of the continuous
!    injection: it is thus possible to impose an injection profile with
!    the continous-injection option.
!

!   reinitialization of the counter of the new particles
npt = nbpart

!   for each boundary zone:
do ii = 1,nfrlag
  izone = ilflag(ii)

!       for each class:
  do iclas = 1, iusncl(izone)

!         if new particles must enter the domain:
    if (mod(ntcabs,iuslag(iclas,izone,ijfre)).eq.0) then

      do ip = npt+1 , npt+iuslag(iclas,izone,ijnbp)

!         number of the original boundary face of injection

      ifac = injfac(ip)
!
!-----------------------------------------------------------
!        EXAMPLE OF MODIFICATION OF THE INJECTION VELOCITY
!        WITH RESPECT TO THE INJECTION POSITION
!-----------------------------------------------------------
!    For instance, the user can call his own subroutine that provides
!    the three components of the instantaneous velocities ettp(ip,jup)
!    ettp(ip,jvp) and  ettp(ip,jwp) with respect to  ettp(ip,jzp)
!    (through interpolation for instance). More simply, the user can provide
!    the three components of the instantaneous velocities, under the form
!    of a mean value (taken arbitrarily here equal to (2,0,0) m/s) added
!    to a fluctuating value (equal here to 0,2 m/s for the 1st and 3rd components)
!
!
        ipnorm = 2
        call normalen(ipnorm,vgauss)
        ettp(ip,jup) = 2.d0 + vgauss(1) * 0.2d0
        ettp(ip,jvp) = 0.d0
        ettp(ip,jwp) = 0.d0 + vgauss(2) * 0.2d0

      enddo

      npt = npt + iuslag(iclas,izone,ijnbp)

    endif

  enddo
enddo

!===============================================================================
! 4. SIMULATION OF THE INSTANTANEOUS TURBULENT FLUID FLOW VELOCITIES SEEN
!    BY THE SOLID PARTICLES ALONG THEIR TRAJECTORIES.
!===============================================================================
!
! Entering this subroutine, the ettp(ip,juf) ettp(ip,jvf) and ettp(ip,jwf) arrays
! are filled with the components of the instantaneous velocity (fluctuation + mean value)
! seen by the particles
!
! When the velocity of the flow is modified just above, most of the time
! the user knows only the mean value. In some flow configurations and some
! injection conditions, it may be necessary to reconstruct the fluctuating part.
! That is why the following routine is called.
!
! Caution: this turbulent component must be reconstructed only on the modified
! velocities of the flow seen.
!
! The reconstruction is unactivated here and must be adapted to the case.
!

if ( 1.eq.0 ) then

  npar1 = nbpart+1
  npar2 = nbpart+nbpnew

  call lagipn                                                     &
  !==========
  ( ncelet , ncel   ,                                             &
    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
    npar1  , npar2  ,                                             &
    itepa  ,                                                      &
    rtpa   ,                                                      &
    ettp   , tepa   , vagaus )

endif

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
