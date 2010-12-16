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

subroutine uslatp &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   numpt  , itepa  ,                                              &
   ia     ,                                                       &
   rep    , uvwr   , romf   , romp   , xnul   , taup   ,          &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   ,                                     &
   ra     )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)
!
! Modification of the calculation of the particle relaxation time
! with respect to the chosen formulation for the drag coefficient

! This subroutine being called in a loop on the particle number,
! be careful not to "load" it too heavily..
!
!            rho             4 d
!               p               p
!      Tau = ---- --------------------------------
!         p
!            rho   3 C     | U [X (t),t] - V (t) |
!               f     drag    f  p          p

!     Tau  : Particle relaxation time
!        p

!     rho  : Particle density
!        p

!     rho  : Fluid density
!        f

!     C    : Drag coefficient
!      drag

!     d    : Particle diameter
!      p

!     U [X (t),t] : Instantaneous velocity of the flow seen
!      f  p

!     V (t) : Particle velocity
!      p

!
!

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
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! numpt            ! i  ! <-- !                                                !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! ia(*)            ! ia ! --- ! macro array of integers                        !
! rep              ! r  ! <-- ! particle Reynolds number                       !
!                  !    !     ! rep = uvwr * ettp(numpt,jdp) / xnul            !
! uvwr             ! r  ! <-- ! particle relative velocity                     !
!                  !    !     ! uvwr= |flow-seen velocity - part. velocity|    !
! romf             ! r  ! <-- ! fluid density at  particle position            !
!                  !    !     !                                                !
! romp             ! r  ! <-- ! particle density                               !
! xnul             ! r  ! <-- ! kinematic viscosity of the fluid at            !
!                  !    !     ! particle position                              !
! taup             ! r  ! --> ! particle relaxation time                       !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! transported variables at cells centers         !
! (ncelet,*)       !    !     ! at the previous time step                      !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use cstnum
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          numpt

integer          itepa(nbpmax,nivep)
integer          ia(*)

double precision rep    , uvwr   , romf   , romp   , xnul  , taup

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ip
double precision fdr

! User-defined local variables

double precision cd1 , cd2 , d2
double precision rec1, rec2, rec3, rec4

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0.  Memory management
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Initializations
!===============================================================================

ip = numpt

!===============================================================================
! 2. Relaxation time with the standard (Wen-Yu) formulation of the drag coefficient
!===============================================================================

! This example is unactivated, it gives the standard relaxation time
! as an indication:

if (1.eq.0) then

  cd1  = 0.15d0
  cd2  = 0.687d0

  if (rep.le.1000) then
      d2 = ettp(ip,jdp) * ettp(ip,jdp)
      fdr = 18.d0 * xnul * (1.d0 + cd1 * rep**cd2) / d2
  else
      fdr = (0.44d0 * 3.d0 / 4.d0) * uvwr / ettp(ip,jdp)
  endif

  taup = romp / romf / fdr

endif

!===============================================================================
! 3. Computation of the relaxation time with the drag coefficient of
!    S.A. Morsi and A.J. Alexander, J. of Fluid Mech., Vol.55, pp 193-208 (1972)
!===============================================================================

rec1 =  0.1d0
rec2 =  1.0d0
rec3 =  10.d0
rec4 = 200.d0

d2 = ettp(ip,jdp) * ettp(ip,jdp)

if ( rep.le.rec1 ) then
  fdr = 18.d0 * xnul / d2

else if ( rep.le.rec2 ) then
  fdr = 3.d0/4.d0 * xnul / d2                                     &
      * (22.73d0 + 0.0903d0/rep + 3.69d0*rep )

else if ( rep.le.rec3 ) then
  fdr = 3.d0/4.d0 * xnul / d2                                     &
      * (29.1667d0 - 3.8889d0/rep + 1.222d0*rep)

else if ( rep.le.rec4 ) then
    fdr = 18.d0*xnul/d2 *(1.d0 + 0.15d0*rep**0.687d0)

else
   fdr = (0.44d0 * 3.d0 / 4.d0) * uvwr / ettp(ip,jdp)
endif

taup = romp / romf / fdr


!==============================================================================

!--------
! Formats
!--------


!----
! End
!----

end subroutine
