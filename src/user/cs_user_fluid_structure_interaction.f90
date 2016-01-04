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

! Purpose :
! ---------
!     User subroutine dedicated to Fluid - Structure internal coupling

! --- Definition and management of internal Fluid Structure coupled calculations,
!     using a simplified solid modeling (linear "mass, friction and spring" modeling).

!     Here are 2 differents subroutines that need to be filled :

!  -   USSTR1 : Called at the beginning of the calculation. It enables one to define
!               internal structures and corresponding initial conditions (initial
!               displacement and velocity).

!  -   USSTR2 : Called at each time step of the calculation. Here one defines
!               structural parameters (considered to be potentially time dependent),
!               i.e. Mass, Friction, Stiffness and Fluid Stresses.

! --- Boundary faces identification
!     =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
subroutine usstr1 &
!================

 ( idfstr ,                                                       &
   aexxst , bexxst , cfopre ,                                     &
   xstr0  , vstr0  , xstreq )

!===============================================================================
! Purpose :
! ---------

! --- Definition of internal structures and corresponding initial conditions
!     (initial displacement and velocity )

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idfstr(nfabor)   ! ia ! <-- ! boundary faces -> structure definition         !
! aexxst,bexxst    ! r  ! <-- ! prediction coefficients of structural data     !
! cfopre           ! r  ! <-- ! prediction coefficients of fluid forces        !
! xstr0(ndim,      ! ra ! <-- ! initial displacement of internal structures    !
!       nbstru)    !    !     !                                                !
! vstr0(ndim,      ! ra ! <-- ! initial velocity of internal structures        !
!       nbstru)    !    !     !                                                !
! xstreq(ndim,     ! ra ! <-- ! displacement of initial mesh compared to       !
!       nbstru)    !    !     ! the structures position at equilibrium         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START


if (1.eq.1) then
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 2.  Definition of internal structures
!===============================================================================

!    Here one fills array IDFSTR(NFABOR)
!    For each boundary face IFAC, IDFSTR(IFAC) is the number of the structure
!    the face belongs to (if IDFSTR(IFAC) = 0, the face IFAC doesn't
!    belong to any structure.)
!    When using internal coupling, structure number necessarily
!    needs to be positive (as shown in following examples).

!    The number of "internal" structures is automatically defined with the
!    maximum value of IDFSTR table, meaning that
!    internal structure numbers must be defined sequentially with positive values,
!    beginning with integer value '1'.


!    In following example, boundary faces with color 4 belong
!    to internal structure '1'.
!    Boundary faces with color 2 belong to internal structure '2'.
!    The total number of internal structures equals 2.

!    Boundary faces identification
!    =============================

!    Boundary faces may be identified using the 'getfbr' subroutine.
!    The syntax of this subroutine is described in the
!    'cs_user_boundary_conditions' subroutine,
!    but a more thorough description can be found in the user guide.


call getfbr('4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = 1

enddo

call getfbr('6', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = 2

enddo

! --- For each internal structure one can here define :
!     - an initial velocity VSTR0
!     - an initial displacement XSTR0 (i.e. XSTR0 is the value of the
!       displacement XSTR compared to the initial mesh at time t = 0)
!     - a displacement compared to equilibrium XSTREQ (i.e. XSTREQ is the
!       initial displacement of the internal structure compared to its position
!       at equilibrium; at each time step t and for a displacement XSTR(t)
!       associated internal structure will be subjected to a force
!       -k*(XSTR(t)+XSTREQ) due to the spring).

! --- Note that XSTR0, XSTREQ and VSTR0 arrays are initialized at the beginning
!     of the calculations to the value of 0.

! --- When starting a calculation using ALE, or re-starting a calculation with
!     ALE basing on a first calculation without ALE, an initial iteration 0 is
!     automatically calculated in order to take initial arrays XSTR0, VSTR0 and
!     XSTREQ into account. In another case add the following expression 'italin=1'
!     in subroutine 'usipsu', so that the code can deal with arrays
!     XSTR0, VSTR0 or XSTREQ.

! --- In the following example :
!     - internal structure '1' has got an initial displacement
!       XSTR0 = 2 (m) in 'y' direction and
!     a displacement compared to equilibrium XSTREQ = 1 (m) in 'y' direction, too.
!     - Initial velocity in 'z' direction of structure '2' equals VSTR0=-0.5 (m/s).

if (ntpabs.le.1) then  ! Avoid resetting in case of restart
  xstr0(2,1)  = 2.d0
  xstreq(2,1) = 1.d0
  vstr0(3,2)  =-0.5d0
endif

! --- Here one can modify the values of the prediction coefficients for
!     displacements anf fluid forces in internal FSI coupled algorithm.
!
!     the use of these coefficients is reminded here :
!       - predicted displacement = X(n) + aexxst * DT * X'(n)
!                                       + bexxst * DT * (X'(n)-X'(n-1))
!             X(n) stands for the displacement at iteration 'n'
!             X'(n) and X'(n-1) represent internal structure velocity respectively
!              at iteration 'n' and 'n-1'.
!       - fluid force sent to structure internal solver = CFOPRE * F(n)
!                                                        + (1.D0-CFOPRE) * F(n-1)
!             F(n) and F(n-1) stand for fluid force acting on the structure
!              respectively at iteration 'n' and 'n-1'.

aexxst =  0.5d0
bexxst =  0.0d0
cfopre =  2.d0

! --- Activation of structural history output (i.e. displacement, structural
!     velocity, structural acceleration anf fluid force)
!     (ihistr=0, not activated ; ihistr=1, activated)
!     The value of structural history output step is the same as the one for
!     standard variables ('NTHIST').
ihistr = 1

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine usstr1


!===============================================================================


subroutine usstr2 &
!================

 ( nbstru ,                                                       &
   idfstr ,                                                       &
   dtcel  ,                                                       &
   xmstru , xcstru , xkstru , xstreq , xstr   , vstr   , forstr , &
   dtstr  )

!===============================================================================
! Purpose :
! ---------

! --- Definition of structural parameters in case of Fluid Structure internal
!     coupling :
!                 Mass, Friction, Stiffness anf Fluid Stresses.
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbstru           ! e  ! <-- ! nombre de structures definies                  !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! idfstr(nfabor    ! te ! <-- ! definition des structures                      !
! dtcel(ncelet)    ! ra ! <-- ! time step (per cell)                           !
! xmstru(ndim,     ! ra ! --> ! matrix of structural mass                      !
!  ndim,nbstru)    !    !     !                                                !
! xcstru(ndim,     ! ra ! --> ! matrix of structural friction                   !
!  ndim,nbstru)    !    !     !                                                !
! xkstru(ndim,     ! ra ! --> ! matrix of structural stiffness                 !
!  ndim,nbstru)    !    !     !                                                !
! xstreq(ndim,     ! ra ! <-- ! displacement of initial mesh compared to       !
!       nbstru)    !    !     ! the structures position at equilibrium         !
! xstr(ndim,       ! ra ! <-- ! structural displacement                        !
!       nbstru)    !    !     !                                                !
! vstr(ndim,       ! ra ! <-- ! structural velocity                            !
!       nbstru)    !    !     !                                                !
! forstr(ndim      ! ra ! <-- ! forces acting on structures (take forces       !
!       nbstru)    !    !     !         due to fluid effects into account   )  !
! dtstr(nbstru)    ! ra ! --> ! structural time step                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

! VARIABLES LOCALES

integer          ii, jj, istr
double precision theta, sint, cost, xm, xc, xk, fx, fy

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START


if (1.eq.1) return

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALIZATION

!===============================================================================


!===============================================================================
! 2.  Structural parameters (subroutine usstr2 is called at each time step
!                            of the calculation)
!===============================================================================

! --- For each internal structure one defines here :
!     - its Mass                  (XMSTRU)
!     - its Friction coefficient C (XCSTRU)
!     - its Stiffness K           (XKSTRU)

!     FORSTR array gives fluid stresses acting on each internal structure.
!     Moreover it's possible to take external forces (gravity for example) into
!     account, too.

!     XSTR array indicates the displacement of the structure compared to its
!     position in the initial mesh.

!     XSTR0 array gives the displacement of the structures in initial mesh
!     compared to structural equilibrium.

!     VSTR array stands for structural velocity.

!     XSTR, XSTR0, and VSTR arrays are DATA tables that can be used to define
!     arrays Mass, Friction and Stiffness. THOSE ARE NOT TO BE MODIFIED.

!     The 3D structural equation that is solved is the following one :

!       M.X'' + C.X' + K.(X+X0) = F   (1)
!       = -     = -    =  - --    -

!       X stands for the structural displacement compared to initil mesh
!         postition (XSTR)
!       -
!       X0 represents the displacement of the structure in initial mesh
!          compared to equilibrium
!       --

!       Note that M, C and K are 3x3 matrices.
!                 =  =     =
!       Equation (1) is solved using a Newmark HHT algorithm.
!       Note that the time step used to solve this equation (DTSTR) can be
!       different from the one of fluid calculations. USER is free to define
!       DTSTR array. At the beginning of the calculation DTSTR is
!       initialized to the value of DTCEL (Fluid time step).
!

! --- Matrices XMSTRU, XCSTRU and XKSTRU are initialized to the value of 0.
do istr = 1, nbstru

  do ii = 1, 3
    do jj = 1, 3
      xmstru(ii,jj,istr) = 0.d0
      xcstru(ii,jj,istr) = 0.d0
      xkstru(ii,jj,istr) = 0.d0
    enddo
  enddo

enddo

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

! --- Example 2): In this example structure '2' is subjected to the following
!                 movement :
!               - In plane xOy the movement is locally defined along an axis
!                 (OX). Structural parameters in X direction are called
!                 xm, xc and xk. The angle of inclination between global (Ox)
!                 axis and local (OX) axis is called THETA. Movement in local (OY)
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


return

end subroutine usstr2


!===============================================================================


subroutine usaste &
!================

 ( idfstr )

!===============================================================================
! Purpose:
! -------

!    User subroutine dedicated the Fluid-Structure external coupling
!    with Code_Aster :
!          Here one defines the boundary faces coupled
!          with Code_Aster and the fluid forces components
!          which are given to structural calculations

!    Boundary faces identification
!    =============================

!    Boundary faces may be identified using the 'getfbr' subroutine.
!    The syntax of this subroutine is described in the
!    'cs_user_boundary_conditions' subroutine,
!    but a more thorough description can be found in the user guide.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idfstr(nfabor)   ! ia ! <-- ! boundary faces -> structure definition         !
! ia(*)            ! ia ! --- ! main integer work array                        !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START


if (1.eq.1) then
  nbaste = 0
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALIZATION
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 2.  Definition of external structures
!===============================================================================

!    Here one fills array IDFSTR(NFABOR)
!    For each boundary face IFAC, IDFSTR(IFAC) is the number of the structure
!    the face belongs to (if IDFSTR(IFAC) = 0, the face IFAC doesn't
!    belong to any structure.)
!    When using external coupling with Code_Aster, structure number necessarily
!    needs to be negative (as shown in following examples).

!    The number of "external" structures with Code Aster is automatically
!    defined with the maximum absolute value of IDFSTR table, meaning that
!    external structure numbers must be defined sequentially with negative values
!    beginning with integer value '-1'.


!    In following example, boundary faces with color 2 and which abscissa X < 2.0
!    belong to external structure '-1'.
!    Boundary faces with color 2 and which abscissa X > 2.0 belong to external
!    structure '-2'. The total number of external structures coupled with Code_Aster
!    equals 2.

!    Boundary faces identification
!    =============================

!    Boundary faces may be identified using the 'getfbr' subroutine.
!    The syntax of this subroutine is described in the
!    'cs_user_boundary_conditions' subroutine,
!    but a more thorough description can be found in the user guide.

!================================================================================

call getfbr('2 and X < 2.0', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = -1

enddo


call getfbr('2 and X > 2.0', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = -2

enddo

! --- The movement of external structure called '-1' is blocked in Z direction.

asddlf(3,1) = 0

! --- The movement of external structure called '-2' is blocked in Z direction.

asddlf(3,2) = 0

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
