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
!>   \par Fluid-Structure coupling with Code_Aster:
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

integer          ifac
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!< [usstr1_init]

!===============================================================================
! 1.  Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!< [usstr1_init]

!< [usstr1_example_a]

!===============================================================================
! 2.  Definition of internal structures
!===============================================================================

call getfbr('4', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = 1

enddo

call getfbr('6', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = 2

enddo

!< [usstr1_example_a]

!< [usstr1_example_b]

if (ntpabs.le.1) then  ! Avoid resetting in case of restart
  xstr0(2,1)  = 2.d0
  xstreq(2,1) = 1.d0
  vstr0(3,2)  =-0.5d0
endif

!< [usstr1_example_b]

!< [usstr1_example_c]

aexxst =  0.5d0
bexxst =  0.0d0
cfopre =  2.d0

!< [usstr1_example_c]

!< [usstr1_example_d]

ihistr = 1

!< [usstr1_example_d]

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
!> \param[in]    nvar           total number of variables
!> \param[in]    nscal          total number of scalars
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

integer          ii, jj, istr
double precision theta, sint, cost, xm, xc, xk, fx, fy

!===============================================================================

!< [usstr2_init]

!===============================================================================
! Definition of structural parameters: mass, friction, stiffness and stresses.
!===============================================================================

! --- Matrices xmstru, xcstru and xkstru are initialized to the value of 0.
do istr = 1, nbstru

  do ii = 1, 3
    do jj = 1, 3
      xmstru(ii,jj,istr) = 0.d0
      xcstru(ii,jj,istr) = 0.d0
      xkstru(ii,jj,istr) = 0.d0
    enddo
  enddo

enddo

!< [usstr2_init]

!< [usstr2_example_1]

! --- Example 1): In following example structure '1' is defined as an isotropic
!     system (i.e. matrices M, C and K are diagonal) : mass equals 5 kg,
!     stiffness equals 2 N/m and friction
!     =  =     =
!     coefficient equals 3 kg.s .

do ii = 1, 3
  xmstru(ii,ii,1) = 5.d0
  xcstru(ii,ii,1) = 2.d0
  xkstru(ii,ii,1) = 3.d0
enddo

!< [usstr2_example_1]

!< [usstr2_example_2]

! --- Example 2): In this example structure '2' is subjected to the following
!                 movement :
!               - In plane xOy the movement is locally defined along an axis
!                 (OX). Structural parameters in X direction are called
!                 xm, xc and xk. The angle of inclination between global (Ox)
!                 axis and local (OX) axis is called theta. Movement in local (OY)
!                 direction is imposed to be rigid.
!               - In 'z' direction the movement is modeled to be oscillating and
!                 harmonic (meaning that there is no friction). Mass equals 1. kg
!                 and stiffness equals 1. N/m. Fluid stresses in that direction
!                 are taken into account. Moreover the structure is also
!                 subjected to an external oscillating
!                 force Fz_ext = 3 * cos(4*t).


!                 This leads to the following local equations :
!                 xm.X'' + xc.X' + xk.X = FX
!                                     Y = 0
!                    Z''         +    Z = FZ + 3.cos(4.t)

theta = pi/6.d0
cost = cos(theta)
sint = sin(theta)

! FX, FY, and FZ stand for the local fluid forces components.
! They are defined as follows, using gobal components of
! fluid forces Fx, Fy and Fz.
!   fx =  cost*Fx + sint*Fy
!   fy = -sint*Fx + cost*Fy
!   fz = Fz

! After changing of basis, the problem can be described as follows,
! using global coordinates:

xm = 1.d0
xc = 3.d-1
xk = 2.d0
fx = forstr(1,2)
fy = forstr(2,2)

xmstru(1,1,2) = xm*cost**2
xmstru(1,2,2) = xm*cost*sint
xmstru(2,1,2) = xm*cost*sint
xmstru(2,2,2) = xm*sint**2
xmstru(3,3,2) = 1.d0

xcstru(1,1,2) = xc*cost**2
xcstru(1,2,2) = xc*cost*sint
xcstru(2,1,2) = xc*cost*sint
xcstru(2,2,2) = xc*sint**2

xkstru(1,1,2) = (xk-1.d0)*cost**2 + 1.d0
xkstru(1,2,2) = (xk-1.d0)*cost*sint
xkstru(2,1,2) = (xk-1.d0)*cost*sint
xkstru(2,2,2) = (xk-1.d0)*sint**2 + 1.d0
xkstru(3,3,2) = 1.d0

forstr(1,2) = fx*cost**2   + fy*sint*cost
forstr(2,2) = fx*sint*cost + fy*sint**2
forstr(3,2) = forstr(3,2) + 3.d0*cos(4.d0*ttcabs)

do istr = 1, nbstru
  dtstr(istr) = dtcel(1)
enddo

!< [usstr2_example_2]

return

end subroutine usstr2


!===============================================================================
! Purpose:
! -------

!> \brief User subroutine dedicated the Fluid-Structure external coupling
!>        with Code_Aster :
!>          Here one defines the boundary faces coupled
!>          with Code_Aster and the fluid forces components
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

integer          ifac
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!< [usaste_init]

!===============================================================================
! 1.  initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!< [usaste_init]

!< [usaste_example]

!===============================================================================
! 2.  Definition of external structures
!===============================================================================

call getfbr('2 and X < 2.0', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = -1

enddo


call getfbr('2 and X > 2.0', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = -2

enddo

! The movement of external structure called '-1' is blocked in Z direction.

asddlf(3,1) = 0

! The movement of external structure called '-2' is blocked in Z direction.

asddlf(3,2) = 0

!< [usaste_example]

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
