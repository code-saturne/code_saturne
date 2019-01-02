!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file usvosy.f90
!>
!> \brief Compute a volume exchange coefficient for SYRTHES coupling
!>
!> See \subpage us_vosy for examples.


subroutine usvosy &
!================

 ( inbcou , ncecpl ,                                              &
   iscal  ,                                                       &
   dt     ,                                                       &
   lcecpl , hvol )

!===============================================================================
!> \brief Compute a volume exchange coefficient for SYRTHES coupling
!>
!> The routine is called in \ref cpvosy for each volume coupling
!> therefore it is necessary to test the value of coupling number to separate
!> the treatments of the different couplings
!>
!> Up to now temperature is the only scalar managed for volume couplings.
!>

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode          name          role                                           !
!______________________________________________________________________________!
!> \param[in]    inbcou        SYRTHES coupling number
!> \param[in]    ncecpl        number of cells implied for this coupling
!> \param[in]    iscal         index number of the temperature scalar
!> \param[in]    dt            time step (per cell)
!> \param[in]    lcecpl        list of coupled cells
!> \param[out]   hvol          volume exchange coefficient to compute
!______________________________________________________________________________!


!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          ncecpl
integer          iscal  , inbcou

integer          lcecpl(ncecpl)

double precision dt(ncelet)
double precision hvol(ncecpl)

!===============================================================================

!----
! End
!----

return
end subroutine usvosy