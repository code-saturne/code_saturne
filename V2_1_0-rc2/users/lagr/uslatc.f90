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

subroutine uslatc &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   numpt  , itepa  ,                                              &
   rep    , uvwr   , romf   , romp   , xnul   ,                   &
   xcp    , xrkl   , tauc   ,                                     &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)
!
! Modification of the computation of the thermal relaxation time
! of the particles with respect to the chosen formulation of the
! Nusselt number.

! This subroutine being called in a loop on the particle number,
! be careful not to "load" it to heavily..
!
!

!               m   Cp
!                p    p
!      Tau = ---------------
!         c          2
!               PI d    h
!                   p    e

!     Tau  : Thermal relaxation time (value to be computed)
!        c

!     m    : Particle mass
!      p

!     Cp   : Particle specific heat
!       p

!     d    : Particle diameter
!      p

!     h    : Coefficient of thermal exchange
!      e

!  The coefficient of thermal exchange is calculated from a Nusselt number,
!  itself evaluated by a correlation (Ranz-Marshall by default)
!
!

!            h  d
!             e  p
!     Nu = --------  = 2 + 0.55 Re **(0.5) Prt**(0.33)
!           Lambda                p

!     Lambda : Thermal conductivity of the carrier field

!     Re     : Particle Reynolds number
!       p

!     Prt    : Prandtl number

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
! numpt            ! i  ! <-- !                                                !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! rep              ! r  ! <-- ! particle Reynolds number                       !
!                  !    !     ! rep = uvwr * ettp(numpt,jdp) / xnul            !
! uvwr             ! r  ! <-- ! relative velocity of the particle              !
!                  !    !     ! uvwr = |flow-seen velocity - part. velocity |  !
! romf             ! r  ! <-- ! fluid density at  particle position            !
!                  !    !     !                                                !
! romp             ! r  ! <-- ! particle density                               !
! xnul             ! r  ! <-- ! kinematic viscosity of the fluid at            !
!                  !    !     ! particle position                              !
! xcp              ! r  ! <-- ! specific heat of the fluid at particle         !
!                  !    !     ! position                                       !
! xrkl             ! r  ! <-- ! diffusion coefficient of the fluid at particle !
!                  !    !     ! position                                       !
! tauc             ! r  ! --> ! thermal relaxation time                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! transported variables at cell centers at       !
! (ncelet,*)       !    !     ! the current time step                          !
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

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          numpt

integer          itepa(nbpmax,nivep)

double precision rep    , uvwr   , romf   , romp   , xnul
double precision xcp    , xrkl   , tauc

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)

! Local variables

integer          ip

! User-defined local variables

double precision prt, fnus

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. Memory management
!===============================================================================


!===============================================================================
! 1. Initializations
!===============================================================================

ip = numpt

!===============================================================================
! 2. Standard thermal relaxation time
!===============================================================================

!   This example is unactivated, it gives the standard thermal relaxation time
!   as an indication.


if (1.eq.0) then

  prt  = xnul / xrkl

  fnus = 2.d0 + 0.55d0 * rep**0.5d0 * prt**(1.d0/3.d0)

  tauc = ettp(ip,jdp) *ettp(ip,jdp) * romp * ettp(ip,jcp)         &
           / ( fnus * 6.d0 * romf * xcp * xrkl )

endif

!==============================================================================

!--------
! Formats
!--------


!----
! End
!----

end subroutine
