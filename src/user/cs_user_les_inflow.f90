!===============================================================================
! User synthetic turbulence inlet definition.
!
! 1) Global characteristics of synthetic turbulence inlets
! 2) Caracteristics of one specific inlet
! 3) Accurate specification of target statistics at inlet
!===============================================================================

!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file cs_user_les_inflow.f90
!>
!> \brief Generation of synthetic turbulence at LES inlets initialization
!>
!> See \subpage les_inflow for examples.
!>
!> \c nent and \c isuisy might be defined.
!>
!> \c nent = Number of inlets
!> \c isuisy = 1: Reading of the LES inflow module restart file
!>           = 0: not activated (synthetic turbulence reinitialized)
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________!
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    nent          number of synthetic turbulence inlets
!______________________________________________________________________________!

subroutine cs_user_les_inflow_init (nent)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use optcal

!===============================================================================

implicit none

! Arguments

integer nent

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_user_les_inflow_init

!===============================================================================

!===============================================================================
! Purpose :
! --------

!> \brief Definition of the characteristics of the synthetic turbulence inlet
!>  \c nument
!>
!> For each LES inlet, the following parameters might be defined:
!>
!>  1. Data relatve to the method employed
!>
!>    - typent indicates the synthetic turbulence method:
!>
!>       - 0: laminar, no turbulent fluctations
!>       - 1: random gaussian noise
!>       - 2: Batten method, based on Fourier mode decomposition
!>       - 3: Synthetic Eddy Method (SEM)
!>       .
!>
!>    - nelent indicates the number of "entities" relative to the method
!>       (useful only for the Batten method and the SEM):
!>
!>       - for Batten : number of Fourier modes of the turbulent fluctuations
!>       - for SEM    : number of synthetic eddies building the fluctuations
!>       .
!>
!>    - iverbo indicates the verbosity level (log)
!>
!>       -   0  no specific output
!>       - > 0 additionnal output (only for SEM)
!>       .
!>
!>
!>  2. Data relative to the LES inflow boundary faces
!>
!>     - nfbent: number of boundary faces of the LES inflow
!>     - lfbent: list of boundary faces of the LES inflow
!>
!>
!>  3. Data relative to the flow
!>
!>     - vitent(3): reference mean velocity vector
!>     - enrent   : reference turbulent kinetic energy
!>     - dspent   : reference dissipation rate
!>
!>       \remarks
!>       - dspent useful only for typent = 2 (Batten) or typent = 3 (SEM).
!>       - Strictly positive values are required for enrent and dspent.
!>       - Accurate specification of the statistics of the flow at LES inlet
!>         can be made via the user subroutine cs_user_les_inflow_advanced.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nument        id of the inlet
!> \param[out]    typent        type of inflow method at the inlet
!> \param[out]    nelent        numb. of entities of the inflow meth
!> \param[out]    iverbo        verbosity level
!> \param[out]    nfbent        numb. of bound. faces of the inlet
!> \param[out]    lfbent        list of bound. faces of the inlet
!> \param[out]    vitent        ref. mean velocity at the inlet
!> \param[out]    enrent        ref. turb. kin. ener. at the inlet
!> \param[out]    dspent        ref. turb. dissipation at the inlet
!_______________________________________________________________________________

subroutine cs_user_les_inflow_define &
( nument, typent, nelent, iverbo,                                             &
  nfbent, lfbent,                                                             &
  vitent, enrent, dspent                                                      &
)

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

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_user_les_inflow_define

!===============================================================================

!===============================================================================
! Purpose :
! --------

!> \brief Generation of synthetic turbulence at LES inlets advanced mode
!>
!> Accurate definition of mean velocity, Reynolds stresses and dissipation
!> rate for each boundary face of the synthetic turbulence inlet \c nument
!>
!> \section Usage
!> \code
!>   uvwent(ndim,nfbent) ! mean velocity vector
!>   rijent(   6,nfbent) ! Reynolds stresses!
!>   epsent(     nfbent) ! dissipation rate
!> \endcode
!>
!> \c rijent components are ordonned as follows: 11, 22, 33, 12, 13, 23
!>
!> Arrays are initialized before this subroutine is called by
!> (see the user subroutine \ref cs_user_les_inflow_define):
!>   \code
!>     uvwent(idim,ifac) = vitent(idim)
!>
!>     rijent(1,ifac)    = 2.d0/3.d0*enrent
!>     rijent(2,ifac)    = 2.d0/3.d0*enrent
!>     rijent(3,ifac)    = 2.d0/3.d0*enrent
!>     rijent(4,ifac)    = 0.d0
!>     rijent(5,ifac)    = 0.d0
!>     rijent(6,ifac)    = 0.d0
!>
!>     epsent(ifac)      = dspent
!>   \endcode
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nument        id of the inlet
!> \param[in]     nfbent        numb. of bound. faces of the inlet
!> \param[in]     nvar          number of variables
!> \param[in]     nscal         number of scalars
!> \param[in]     lfbent        list of bound. faces of the inlet
!> \param[in]     dt            time step
!> \param[out]    uvwent        mean velocity at the inlet faces
!> \param[out]    rijent        turb. kin. ener. at the inlet faces
!> \param[out]    epsent        turb. dissipation at the inlet faces
!_______________________________________________________________________________

subroutine cs_user_les_inflow_advanced &
 ( nument , nfbent ,                                              &
   nvar   , nscal ,                                               &
   lfbent ,                                                       &
   dt     ,                                                       &
   uvwent , rijent , epsent )

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

double precision dt(ncelet)
double precision uvwent(ndim,nfbent), rijent(6,nfbent)
double precision epsent(nfbent)

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_user_les_inflow_advanced
