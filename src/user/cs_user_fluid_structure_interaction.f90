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
! --------

!> \file cs_user_fluid_structure_interaction.f90
!>
!> \brief User subroutines dedicated to Fluid-Structure interaction modeling
!>
!>   \par Management of internal Fluid-Structure coupled calculations:
!>        a simplified solid model is used
!>        (linear "mass, friction and spring" modeling).
!>        Here are 2 differents subroutines that need to be filled :
!>
!>     - \ref usstr1 : Called at the beginning of the calculation. It enables
!>                     one to define internal structures and corresponding
!>                     initial conditions (initial displacement and velocity).
!>
!>     - \ref usstr2 : Called at each time step of the calculation. Here one
!>                     structural parameters
!>                     (considered to be potentially time dependent),
!>                     i.e. Mass, Friction, Stiffness and Fluid Stresses.
!>
!>   \par Fluid-Structure coupling with code_aster:
!>        the user subroutine \ref usaste has to be used.
!>
!>   \par Examples of data settings for fluid-structure interaction (FSI):
!>        Several examples are available
!>        \ref cs_user_fluid_structure_interaction "here".
!
!-------------------------------------------------------------------------------

!===============================================================================
! Function :
! ----------

!> \brief Definition of internal structures and corresponding initial conditions
!>       (initial displacement and velocity ).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    idfstr         boundary faces -> structure definition
!> \param[in]    aexxst,bexxst  prediction coefficients of structural data
!> \param[in]    cfopre         prediction coefficients of fluid forces
!> \param[in]    xstr0          initial displacement of internal structures
!> \param[in]    vstr0          initial velocity of internal structures
!> \param[in]    xstreq         displacement of initial mesh compared to
!>                              structures position at equilibrium
!______________________________________________________________________________!

subroutine usstr1 &
 ( idfstr ,                                                       &
   aexxst , bexxst , cfopre ,                                     &
   xstr0  , vstr0  , xstreq )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use optcal
use entsor
use albase
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idfstr(nfabor)

double precision aexxst, bexxst, cfopre
double precision xstr0(3,nstrmx), xstreq(3,nstrmx)
double precision vstr0(3,nstrmx)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!===============================================================================
! 2.  Definition of internal structures
!===============================================================================


! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine usstr1


!===============================================================================
! Purpose :
! ---------

!> \brief Definition of structural parameters in case of Fluid Structure
!>        internal coupling : Mass, Friction, Stiffness anf Fluid Stresses.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    nbstru         nombre de structures definies
!> \param[in]    idfstr         definition des structures
!> \param[in]    dtcel          time step (per cell
!> \param[in]    xmstru         matrix of structural mass
!> \param[out]   xcstru         matrix of structural friction
!> \param[out]   xkstru         matrix of structural stiffness
!> \param[in]    xstreq         displacement of initial mesh compared to
!> \param[in]    xstr           structural displacement
!> \param[in]    vstr           structural velocity
!> \param[in]    forstr         forces acting on structures (take forces)
!> \param[out]   dtstr          structural time step
!______________________________________________________________________________!

subroutine usstr2 &
 ( nbstru ,                                                       &
   idfstr ,                                                       &
   dtcel  ,                                                       &
   xmstru , xcstru , xkstru , xstreq , xstr   , vstr   , forstr , &
   dtstr  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use optcal
use albase
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbstru

integer          idfstr(nfabor)

double precision dtcel(ncelet)
double precision xmstru(3,3,nstrmx)
double precision xcstru(3,3,nstrmx)
double precision xkstru(3,3,nstrmx)
double precision xstreq(3,nstrmx)
double precision xstr(3,nstrmx)
double precision vstr(3,nstrmx)
double precision forstr(3,nstrmx)
double precision dtstr(nstrmx)

! Local variables

!===============================================================================

!===============================================================================
! Definition of structural parameters: mass, friction, stiffness and stresses.
!===============================================================================

return

end subroutine usstr2


!===============================================================================
! Purpose:
! -------

!> \brief User subroutine dedicated the Fluid-Structure external coupling
!>        with code_aster :
!>          Here one defines the boundary faces coupled
!>          with code_aster and the fluid forces components
!>          which are given to structural calculations.
!

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     idfstr        boundary faces -> structure definition
!______________________________________________________________________________!

subroutine usaste &
 ( idfstr )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use optcal
use entsor
use albase
use parall
use period
use alaste
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbstru

integer          idfstr(nfabor)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!===============================================================================
! 1.  initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!===============================================================================
! 2.  Definition of external structures
!===============================================================================

!----
! Formats
!----

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine usaste
