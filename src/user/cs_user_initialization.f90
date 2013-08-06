!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file cs_user_initialization.f90
!>
!> \brief Initialize variables
!>
!> This subroutine is called at beginning of the computation
!> (restart or not) before the loop time step.
!>
!> This subroutine enables to initialize or modify (for restart)
!> unkown variables and time step values.
!>
!> \c rom and \c viscl values are equal to \c ro0 and \c viscl0 or initialize
!> by reading the restart file.
!> viscls and cp variables (when there are defined) have no value
!> excepted if they are read from a restart file.
!>
!> Physical quantities are defined in the following arrays:
!> \code
!>  propce ! physical quantities defined at cell center
!>  propfa ! physical quantities defined at border face center
!> \endcode
!>
!> Examples:
!> \code
!>  propce(iel, ipproc(irom  )) ! means rom  (iel)
!>  propce(iel, ipproc(iviscl)) ! means viscl(iel)
!>  propce(iel, ipproc(icp   )) ! means cp   (iel)
!>  propce(iel, ipproc(ivisls(iscal))) ! means visls(iel, iscal)
!>  propfb(ifac, ipprob(irom )) ! means romb  (ifac)
!> \endcode
!>
!> Modification of the behaviour law of physical quantities (rom, viscl,
!> viscls, cp) is not done here. It is the purpose of the user subroutine
!> \ref cs_user_physical_properties.
!>
!> \section cell_id Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time step)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!_______________________________________________________________________________

subroutine cs_user_initialization &
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfa , propfb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)

! Local variables

! INSERT_VARIABLE_DEFINITIONS_HERE

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

if (1.eq.1) then
!       Tag to know if a call to this subroutine has already been done
  iusini = 0
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(ncel)) ! temporary array for cells selection

! INSERT_MAIN_CODE_HERE

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_initialization
