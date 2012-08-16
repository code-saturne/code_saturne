!===============================================================================
! User synthetic turbulence inlet definition.
!
! 1) Global caracteristics of synthetic turbulence inlets
! 2) Caracteristics of one specific inlet
! 3) Accurate specification of target statistics at inlet
!===============================================================================

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

subroutine cs_user_les_inflow_init (nent)
!=================================

!===============================================================================
! Purpose :
! --------

! Generation of synthetic turbulence at LES inlets

! Definition of global caracteristics of synthetic turbulence inlets

! nent and isuisy might be defined.

! nent = Number of inlets
! isuisy = 1: Reading of the LES inflow module restart file
!        = 0: not activated (synthetic turbulence reinitialized)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nent             ! i  ! --> ! number of synthetic turbulence inlets          !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================


!===============================================================================
! Module files
!===============================================================================

use optcal

!===============================================================================

implicit none

! Arguments

integer nent

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! INSERT_MAIN_CODE_HERE

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine

!===============================================================================

subroutine cs_user_les_inflow_define &
!===================================
( nument, typent, nelent, iverbo,                                             &
  nfbent, lfbent,                                                             &
  vitent, enrent, dspent                                                      &
)

!===============================================================================
! Purpose :
! --------

! Generation of synthetic turbulence at LES inlets

! Definition of the caracteristics of the synthetic turbulence inlet 'nument'

! For each LES inlet, the following parameters might be defined:

!  1. Data relatve to the method employed

!     * typent indicates the synthetic turbulence method:

!         0 : laminar, no turbulent fluctations
!         1 : random gaussian noise
!         2 : Batten method, based on Fourier mode decomposition
!         3 : Synthetic Eddy Method (SEM)

!     * nelent indicates the number of "entities" relative to the method
!       (usefull only for the Batten method and the SEM):

!         for Batten : number of Fourier modes of the turbulent fluctuations
!         for SEM    : number of synthetic eddies building the fluctuations

!     * iverbo indicates the verbosity level (listing)

!          0  no specific output
!         > 0 additionnal output (only for SEM)


!  2. Data relative to the LES inflow boundary faces

!       nfbent: number of boundary faces of the LES inflow
!       lfbent: list of boundary faces of the LES inflow


!  3. Data relative to the flow

!       vitent(3): reference mean velocity vector
!       enrent   : reference turbulent kinetic energy
!       dspent   : reference dissipation rate

!       Note :
!       ----
!       - dspent usefull only for typent = 2 (Batten) or typent = 3 (SEM).
!       - Strictly positive values are required for enrent and dspent.
!       - Accurate specification of the statistics of the flow at LES inlet
!         can be made via the user subroutine cs_user_les_inflow_advanced.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nument           ! i  ! <-- ! id of the inlet                                !
! typent           ! i  ! --> ! type of inflow method at the inlet             !
! nelent           ! i  ! --> ! numb. of entities of the inflow meth           !
! iverbo           ! i  ! --> ! verbosity level                                !
! nfbent           ! i  ! --> ! numb. of bound. faces of the inlet             !
! lfbent           ! ra ! --> ! list of bound. faces of the inlet              !
! vitent           ! ra ! --> ! ref. mean velocity at the inlet                !
! enrent           ! r  ! --> ! ref. turb. kin. ener. at the inlet             !
! dspent           ! r  ! --> ! ref. turb. dissipation at the inlet            !
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
use entsor
use mesh
use cstnum
use cstphy
use parall

!===============================================================================

implicit none

! Arguments

integer          nument, iverbo
integer          nfbent, typent, nelent
integer          lfbent(nfabor)

double precision vitent(3), enrent, dspent

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! INSERT_MAIN_CODE_HERE

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine

!===============================================================================

subroutine cs_user_les_inflow_advanced &
!=====================================

 ( nument , nfbent ,                                              &
   nvar   , nscal ,                                               &
   lfbent ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   uvwent , rijent , epsent )

!===============================================================================
! Purpose :
! --------

!    User subroutine.

!    Generation of synthetic turbulence at LES inlets

!    Accurate definition of mean velocity, Reynolds stresses and dissipation
!    rate for each boundary face of the synthetic turbulence inlet 'nument'

! Usage
! -----
!     uvwent(ndim,nfbent) : mean velocity vector
!     rijent(   6,nfbent) : Reynolds stresses!
!     epsent(     nfbent) : dissipation rate

!    rijent components are ordonned as follow : 11, 22, 33, 12, 13, 23

!    Arrays are initialized before this subroutine is called by
!    (see the user subroutine cs_user_les_inflow_define):

!       uvwent(idim,ifac) = vitent(idim)

!       rijent(1,ifac)    = 2.d0/3.d0*enrent
!       rijent(2,ifac)    = 2.d0/3.d0*enrent
!       rijent(3,ifac)    = 2.d0/3.d0*enrent
!       rijent(4,ifac)    = 0.d0
!       rijent(5,ifac)    = 0.d0
!       rijent(6,ifac)    = 0.d0

!       epsent(ifac)      = dspent

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nument           ! i  ! --> ! id of the inlet                                !
! nfbent           ! i  ! --> ! numb. of bound. faces of the inlet             !
! nvar             ! i  ! --> ! number of variables                            !
! nscal            ! i  ! --> ! number of scalars                              !
! lfbent           ! i  ! --> ! list of bound. faces of the inlet              !
! dt               ! r  ! --> ! time step                                      !
! rtpa             ! ra ! --> ! variables at cells (previous)                  !
! rtp              ! ra ! --> ! variables at cells                             !
! propce           ! ra ! --> ! physical properties at cells                   !
! propfa           ! ra ! --> ! physical properties at faces                   !
! propfb           ! ra ! --> ! physical properties at bound. faces            !
! coefa            ! ra ! --> ! boundary conditions array                      !
! coefb            ! ra ! --> ! boundary conditions array                      !
! uent             ! ra ! <-- ! mean velocity at the inlet faces               !
! rijent           ! ra ! <-- ! turb. kin. ener. at the inlet faces            !
! epsent           ! ra ! <-- ! turb. dissipation at the inlet faces           !
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
use entsor
use cstphy
use cstnum
use mesh

!===============================================================================

implicit none

! Arguments

integer          nument , nfbent
integer          nvar   , nscal
integer          lfbent(nfbent)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision uvwent(ndim,nfbent), rijent(6,nfbent)
double precision epsent(nfbent)

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! INSERT_MAIN_CODE_HERE

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
