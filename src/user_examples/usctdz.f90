!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!>
!> \file usctdz.f90
!>
!> \brief Define cooling tower parameters.
!>
!> See \subpage us_ctdz for examples.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine usctdz
!================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use entsor
use cstphy
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl
use mesh

!===============================================================================

implicit none

! Local variables

integer          imzech,ntypze,idimze,neleze

double precision teaueze,qeaueze,deltat
double precision xap,xnp,surface,dgout

!===============================================================================

!< [usctdz]

! The chosen model must be consistent with ippmod

! IMZECH = 0 - no model
!          1 - Merkels model
!          2 - Poppe's model

imzech  = ippmod(iaeros)

! Define exchange zone

idimze  = 2
ntypze  = 2
neleze  = 20
deltat  = 10.0d0
teaueze = 36d0
qeaueze = 33737.d0
xap     = 0.2d0
xnp     = 0.5d0
surface = 5183.6d0
dgout   = 0.005d0

call defct &
!=========
 (idimze, '2 or 3', imzech, ntypze, neleze,           &
  deltat, teaueze, qeaueze, xap, xnp, surface, dgout)

!< [usctdz]

!----
! End
!----

return
end subroutine usctdz