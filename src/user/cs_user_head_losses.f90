!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
! --------

!> \file cs_user_head_losses.f90
!>
!> \brief Define Head losses
!>
!> The subroutine cs_user_head_losses is called at three different stages
!> in the code (iappel = 1, 2 or 3)
!>
!> iappel = 1:
!>    Calculation of the number of cells where a head loss term is
!>    imposed: ncepdp
!>    Called once at the beginning of the calculation
!>
!> iappel = 2
!>    Identification of the cells where a head loss term is imposed:
!>    array icepdc(ncepdc)
!>    Called once at the beginning of the calculation
!>
!> iappel = 3
!>    Calculation of the values of the head loss term
!>    Called at each time step
!>
!> Note that calling this subroutine completely overwrites head losses
!> defined using the GUI.
!>
!> ckupdc is the local head loss term
!>
!> It appears on the momentum as follows:
!>    rho du/dt = - grad p + head_loss        (+ other terms)
!>                      with head_loss = - rho ckupdc u (in kg/(m2 s2))
!>
!> For a distributed head loss,
!>
!>    let ksil = dhl/(0.5 rho u**2) given by the litterature
!>    (dhl is the head loss per unit length)
!>
!>    the source term tspdc is equal to dhl = - ksil *(0.5 rho u**2)
!>
!>    we have ckupdc = 0.5 ksil abs(U)
!>
!>
!> For a singular head loss,
!>
!>    let ksil = dhs/(0.5 rho u**2) given by the litterature
!>    (dhs is the singular head loss)
!>
!>    the source term tspdc is equal to dhs/L = - ksil/L *(0.5 rho u**2)
!>
!>    we have ckupdc = 0.5 ksil/L abs(u)
!>
!>    where L is the length over which we have chosen to represent the
!>    singular head loss
!>
!>
!> Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     iappel        stage in the code
!> \param[in]     icepdc        numbers of ncepdp cells with head loss
!> \param[in]     izcpdc        cells zone for head loss definition
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head loss
!_______________________________________________________________________________


subroutine cs_f_user_head_losses &
 ( ncepdp , iappel ,                                              &
   icepdc , izcpdc ,                                              &
   dt     ,                                                       &
   ckupdc )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          ncepdp
integer          iappel

integer          icepdc(ncepdp)
integer          izcpdc(ncel)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6)

!===============================================================================

return
end subroutine cs_f_user_head_losses
