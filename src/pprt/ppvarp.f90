!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[out]    nmodpp        number of activated paricle physic models
!______________________________________________________________________________

subroutine ppvarp ( nmodpp )
!================

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use elincl
use cs_coal_incl
use cs_fuel_incl
use ppcpfu
use atincl
use ihmpre

!===============================================================================

implicit none

! Arguments

integer       nmodpp

!===============================================================================

! 1. Gas combustion
!------------------

if (     ippmod(icod3p).ge.0   &
    .or. ippmod(icoebu).ge.0   &
    .or. ippmod(icolwc).ge.0) then
  call covarp
  !==========
endif

! Number of Diracs for LWC model

if (     ippmod(icolwc).eq.0  &
    .or. ippmod(icolwc).eq.1) then
  ndirac = 2

else if (      ippmod(icolwc).eq.2  &
         .or.  ippmod(icolwc).eq.3) then
  ndirac = 3

else if (     ippmod(icolwc).eq.4  &
         .or. ippmod(icolwc).eq.5) then
  ndirac = 4
endif

! 2. Coal combustion
!-------------------

if (ippmod(iccoal).ge.0) then
  call cs_coal_varpos
  !==================
endif

! Pulverized coal combustion coupled with Lagrangian model

if (ippmod(icpl3c).ge.0) then
  call cplvar
  !==========
endif

! 3. Compressible model
!----------------------

if (ippmod(icompf).ge.0) then
  call cfvarp
  !==========
endif

! 4. Electric arcs model
!-----------------------

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then
  call elvarp
  !==========
endif

! 5. Fuel combustion model
!-------------------------

if (ippmod(icfuel).ge.0) then
  call cs_fuel_varpos
  !==================
endif

! 6. Atmospheric model
!---------------------

if (ippmod(iatmos).ge.1) then
  call atvarp
  !==========
endif

! 7. Cooling towers model
!------------------------

if (ippmod(iaeros).ge.0) then
  call ctvarp
  !==========
endif

return
end subroutine
