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

subroutine uslain &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nptnew ,                                                       &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , itepa  , ifrlag , injfac ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , w1     , w2     , w3     ,          &
   rdevel , rtuser , ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
!                  !    !     !                                                !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
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
! nideve nrdeve    ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse nrtuse    ! i  ! <-- ! sizes of ituser and rtuser arrays              !
!                  !    !     !                                                !
! ifacel           ! ia ! <-- ! interior faces -> cells connectivity           !
! (2, nfac)        !    !     !                                                !
! ifabor           ! ia ! <-- ! boundary faces -> cells connectivity           !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! ia ! <-- ! boundary face family numbers                   !
! (nfabor)         !    !     !                                                !
! ifmcel           ! ia ! <-- ! cell family numbers                            !
! (ncelet)         !    !     !                                                !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! ia ! <-- ! interior faces -> vertices index (optional)    !
!   (lndfac)       !    !     !                                                !
! nodfac           ! ia ! <-- ! interior faces -> vertices list (optional)     !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! ia ! <-- ! boundary faces -> vertices index (optional)    !
!   (lndfbr)       !    !     !                                                !
! nodfbr           ! ia ! <-- ! boundary faces -> vertices list  (optional)    !
!   (nfabor+1)     !    !     !                                                !
! itrifb(nfabor    ! ia ! <-- ! indirection for the sorting of the boundary    !
!  nphas      )    !    !     ! faces                                          !
! itypfb(nfabor    ! ia ! <-- ! type of the boundary faces                     !
!  nphas      )    !    !     !                                                !
! ifrlag(nfabor    ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! injfac(nptnew    ! ia ! <-- ! number of the injection boundary face          !
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! <-- ! complementary user array of integers           !
! ia(*)            ! ia ! --- ! macro array of integers                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
! (ndim,nnod)      !    !     !                                                !
! volume           ! ra ! <-- ! cell volumes                                   !
! (ncelet          !    !     !                                                !
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
! w1..w3(ncelet    ! ra ! --- ! work arrays                                    !
! rdevel(nrdeve    ! ra ! <-- ! dev. complementary array of reals              !
! rtuser(nrtuse    ! ra ! <-- ! user complementary array of reals              !
! ra(*)            ! ra ! --- ! macro array of reals                           !
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

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nptnew
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)
integer          injfac(nbpnew)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra
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

idebia = idbia0
idebra = idbra0

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
  ( idebia , idebra ,                                             &
    ncelet , ncel   ,                                             &
    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
    npar1  , npar2  ,                                             &
    nideve , nrdeve , nituse , nrtuse ,                           &
    itepa  ,                                                      &
    idevel , ituser , ia     ,                                    &
    rtpa   ,                                                      &
    ettp   , tepa   , vagaus ,                                    &
    w1     , w2     , w3     ,                                    &
    rdevel , rtuser , ra     )

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
