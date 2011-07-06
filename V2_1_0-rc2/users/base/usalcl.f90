!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine usalcl &
!================

 ( itrale ,                                                       &
   nvar   , nscal  ,                                              &
   icodcl , itypfb , ialtyb , impale ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl , xyzno0 , depale )

!===============================================================================
! Purpose:
! -------

! --- User subroutine dedicated the use of ALE (Arbitrary Lagrangian Eulerian Method) :
!                 Fills boundary conditions (ialtyb, icodcl, rcodcl) for mesh velocity.
!                 This subroutine also enables one to fix displacement on nodes.

!
! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in 'usclim' subroutine,
! but a more thorough description can be found in the user guide.

! Boundary conditions setup for standard variables (pressure, velocity,
! turbulence, scalars) is described precisely in 'usclim' subroutine.

! Detailed explanation will be found in the theory guide.

! Boundary condition types
! ========================

! Boundary conditions may be assigned in two ways.

!
!    For "standard" boundary conditions:
!    -----------------------------------

!     (fixed boundary, sliding mesh boundary, fixed velocity), one defines a code in the 'ialtyb'
!     array (of dimensions number of boundary faces, number of phases).

! * ialtyb(ifac) = ibfixe : the face IFAC is considered to be motionless. A zero Dirichlet
!         boundary condition is automatically imposed on mesh velocity. Moreover the displacement
!         of corresponding nodes will automatically be set to 0 (for further information please
!         read the paragraph dedicated to the description of IMPALE array in 'usalcl' subroutine),
!         unless the USER has modified the condition of at least one  component of mesh velocity
!         (modification of ICOCL array, please read the following paragraph 'For "non-standard"
!         conditions')

! * ialtyb(ifac) = igliss : The mesh slides on corresponding face IFAC. The normal component of mesh
!          viscosity is automatically set to 0. A homogeneous Neumann condition is automatically
!          prescribed for the other components, as it's the case for 'Symmetry' fluid condition (Please
!          note that homogeneous Neumann condition is only partially implicit in case of boudary face
!          that is not aligned with axis).

! * ialtyb(ifac) = ivimpo : the mesh velocity is imposed on face IFAC. Thus, the users needs to
!          specify the mesh velocity values filling RCODCL arrays as follows :
!          rcodcl(ifac,iuma,1) = mesh velocity in 'x' direction
!          rcodcl(ifac,ivma,1) = mesh velocity in 'y' direction
!          rcodcl(ifac,iwma,1) = mesh velocity in 'z' direction
!          Components of rcodcl(.,i.ma,1) arrays that are not specified by user
!          will automatically be set to 0, meaning that user only needs to specify
!          non zero mesh velocity components.

!    For "non-standard" conditions:
!    ------------------------------

!     Other than (fixed boundary, sliding mesh boundary, fixed velocity), one defines
!     for each face and each component IVAR = IUMA, IVMA, IWMA :
!        -> a code             icodcl(ifac, ivar)
!        -> three real values  rcodcl(ifac, ivar, 1)
!                              rcodcl(ifac, ivar, 2)
!                              rcodcl(ifac, ivar, 3)
!     The value of 'icodcl' is taken from the following:
!       1: Dirichlet
!       3: Neumann
!       4: Symmetry
!     The values of the 3 'rcodcl' components are:
!      rcodcl(ifac, ivar, 1):
!         Dirichlet for the variable if icodcl(ifac, ivar) =  1
!         The dimension of rcodcl(ifac, ivar, 1) is in m/s
!      rcodcl(ifac, ivar, 2):
!         "exterior" exchange coefficient (between the prescribed value
!                          and the value at the domain boundary)
!                          rinfin = infinite by default
!           rcodcl(ifac,ivar,2) =  (VISCMA) / d
!              (D has the dimension of a distance in m, VISCMA stands for
!              the mesh viscosity)
!         NB : the definition of rcodcl(.,.,2) is based on the manner
!              other standard variables are managed in the same case.
!              This type of boundary condition appears nonsense
!              concerning mesh in that context.

!      rcodcl(ifac,ivar,3) :
!        Flux density (in kg/m s2) = J if icodcl(ifac, ivar) = 3
!                     (<0 if gain, n outwards-facing normal)
!             rcodcl(ifac,ivar,3) = -(VISCMA)* (grad Um).n
!                  (Um represents mesh velocity)
!         NB : note that the definition of condition rcodcl(ifac,ivar,3)
!              is based on the manner other standard variables are
!              managed in the same case.
!              rcodcl(.,.,3) = 0.d0 enables one to specify a homogeneous
!              Neuman condition on mesh velocity. Any other value will be
!              physically nonsense in that context.

!      Note that if the user assigns a value to ialtyb equal to ibfixe, igliss,
!      or ivimpo and does not modify icodcl (zero value by
!       default), ialtyb will define the boundary condition type.

!      To the contrary, if the user prescribes icodcl(ifac, ivar) (nonzero),
!        the values assigned to rcodcl will be used for the considered face
!        and variable (if rcodcl values are not set, the default values will
!        be used for the face and variable, so:
!                                 rcodcl(ifac, ivar, 1) = 0.d0
!                                 rcodcl(ifac, ivar, 2) = rinfin
!                                 rcodcl(ifac, ivar, 3) = 0.d0)


!      If the user decides to prescribe his own non-standard boundary conditions
!      it will be necessary to assign values to icodcl AND to rcodcl for ALL
!      mesh viscosity components. Thus, the user does not need to assign values
!      to IALTYB for each associated face, as it will not be taken into account
!      in the code.



! Consistency rules
! =================

!       A consistency rules between 'icodcl' codes for variables with
!       non-standard boundary conditions:
!            If a symmetry code (ICODCL=4) is imposed for one mesh velocity
!            component, one must have the same condition for all other mesh
!            velocity components.


! Fixed displacement on nodes
! ============================
!  For a better precision concerning mesh displacement, one can also assign values
!    of displacement to certain internal and/or boundary nodes. Thus, one
!    need to fill DEPALE and IMPALE arrays :
!    depale(inod,1) = displacement of node inod in 'x' direction
!    depale(inod,2) = displacement of node inod in 'y' direction
!    depale(inod,3) = displacement of node inod in 'z' direction
!    This array is defined as the total displacement of the node compared
!    its initial position in initial mesh.
!    impale(inod) = 1 indicates that the displacement of node inod is imposed
!    (Note that IMPALE array is initialized to the value of 0; if its value
!    is not modified, corresponding value in DEPALE array will not be
!    taken into account)

!  During mesh's geometry re-calculation at each time step, the position of the nodes, which
!    displacement is fixed ( i.e. IMPALE=1), is not calculated using the value of mesh viscosity
!    at the center of corresponding cell, but directly filled using the values of DEPALE.
!  If the displacement is fixed for all nodes of a boundary face it's not necessary to
!    prescribe boundary conditions at this face on mesh viscosity. ICODCL and RCODCL values will
!    be overwritten :
!    -> ICODCL is automatically set to 1 (Dirichlet)
!    -> RCODCL value will be automatically set to face's mean mesh velocity value, that is
!       calculated using DEPALE array.

!  If a fixed boundary condition (ialtyb(ifac)=ibfixe) is imposed to the face ifac,
!    the displacement of each node inod belonging to ifac is considered to be fixed,
!    meaning that impale(inod) = 1 and depale(inod,.) = 0.d0.


! Description of nodes
! ====================
! NNOD gives the total (internal and boundary) number of nodes.
! Vertices coordinates are given by XYZNOD(3, NNOD) array. This table is
! updated at each time step of the calculation.
! XYZNO0(3,NNOD) gives the coordinates of initial mesh at the beginning
! of the calculation.

! The faces - nodes connectivity is stored by means of four integer arrays :
! IPNFAC, NODFAC, IPNFBR, NODFBR.

! NODFAC (NODFBR) stores sequentially the index-numbers of the nodes of each
! internal (boundary) face.

! IPNFAC (IPNFBR) gives the position of the first node of each internal
! (boundary) face in the array NODFAC (NODFBR).

! For example, in order to get all nodes of internal face IFAC, one can
! use the following loop :
!   DO II = IPNFAC(IFAC), IPNFAC(IFAC+1)-1 <- index number of NODFAC array
!                                             corresponding to IFAC
!
!     INOD = NODFAC(II)                    <- index-number IIth node of face IFAC.
!
!
!     ...
!   ENDDO

! Influence on boundary conditions related to fluid velocity
! ==========================================================
! The effect of fluid velocity and ALE modeling on boundary faces that
! are declared as walls (ITYPFB = IPAROI or IPARUG) really depends on
! the physical nature of this interface.
! Indeed when studying an immersed structure the motion of corresponding
! boundary faces is the one of the structure, meaning that it leads to
! fluid motion. On the other hand when studying a piston the motion of vertices
! belonging to lateral boundaries has no physical meaning therefore it has
! no influence on fluid motion.
! Whatever the case, mesh velocity component that is normal to the boundary
! face is always taken into account (Ufluid.n = Wmesh.n). The modeling
!                                    -      -   -     -
! of tangential mesh velocity component differs from one case to another.

! The influence of mesh velocity on boundary conditions for fluid modeling is
! managed and modeled in Code_Saturne as follows :
!  - If ialtyb(ifac) = ibfixe : mesh velocity equals 0. (In case of 'fluid sliding
!  wall' modeling corresponding condition will be specified in Code_Saturne
!  Interface or in 'usclim' subroutine.)
!  - If ialtyb(ifac) = ivimpo : tangential mesh velocity is modeled as a sliding
!  wall velocity in fluid boundary conditions unless a value for fluid sliding
!  wall velocity has been specified by USER in Code_Saturne Interface
!  or in 'usclim' subroutine.
!  - If ialtyb(ifac) = igliss : tangential mesh velocity is not taken into account
!  in fluid boundary conditions (In case of 'fluid sliding wall' modeling
!  corresponding condition will be specified in Code_Saturne Interface
!  or in 'usclim' subroutine.)
!  - If impale(inod) = 1 for all vertices of a boundary face : tangential mesh
!  velocity value that has been derived from nodes displacement is modeled as a
!  sliding wall velocity in fluid boundary conditions unless a value for fluid
!  sliding wall velocity has been specified by USER in Code_Saturne Interface or
!  in 'usclim' subroutine.

! Note that mesh velocity has no influence on modeling of
! boundary faces with 'inlet' or 'free outlet' fluid boundary condition.

! For "non standard" conditions USER has to manage the influence of boundary
! conditions for ALE method (i.e. mesh velocity) on the ones for Navier Stokes
! equations(i.e. fluid velocity). (Note that fluid boundary conditions can be
! specified in this subroutine.)

! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.

! Faces identification
! ====================

! Faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
! itrale           ! i  ! <-- ! number of iterations for ALE method            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itypfb           ! ia ! --> ! boundary face types                            !
! ialtyb (nfabor)  ! ia ! --> ! boundary face types for mesh velocity          !
! impale(nnod)     ! ia ! <-- ! indicator for fixed node displacement          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!  (nfabor,nvar,3) !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
! depale(nnod,3)   ! ra ! <-- ! nodes displacement                             !
! xyzno0           ! ra ! <-- ! vertex coordinates of initial mesh             !
!  (3, nnod)       !    !     !                                                !
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
use cstphy
use cstnum
use entsor
use parall
use period
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          itrale
integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor), ialtyb(nfabor)
integer          impale(nnod)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision depale(nnod,3), xyzno0(3,nnod)

! Local variables

integer          ifac, iel, ii
integer          inod
integer          ilelt, nlelt

double precision delta, deltaa

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================
if(iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : stop in definition of boundary conditions   ',/,&
'@    =========                                               ',/,&
'@     ALE Method has been activated                          ',/,&
'@     User subroutine ''usalcl'' must be completed           ',/, &
'@                                                            ',/,&
'@  The calculation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))



!===============================================================================
! 2.  Assign boundary conditions to boundary faces here

!     One may use selection criteria to filter boundary case subsets
!       Loop on faces from a subset
!         Set the boundary condition for each face
!===============================================================================

!     Calculation of displacement at current time step
deltaa = sin(3.141596d0*(ntcabs-1)/50.d0)
delta  = sin(3.141596d0*ntcabs/50.d0)

! --- For boundary faces of color 4 assign a fixed velocity

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     ELEMENT ADJACENT A LA FACE DE BORD
  iel = ifabor(ifac)

  ialtyb(ifac) = ivimpo
  rcodcl(ifac,iuma,1) = 0.d0
  rcodcl(ifac,ivma,1) = 0.d0
  rcodcl(ifac,iwma,1) = (delta-deltaa)/dt(iel)

enddo

! --- For boundary faces of color 5 assign a fixed displacement on nodes

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
    inod = nodfbr(ii)
    if (impale(inod).eq.0) then
      depale(inod,1) = 0.d0
      depale(inod,2) = 0.d0
      depale(inod,3) = delta
      impale(inod) = 1
    endif
  enddo

enddo

! --- For boundary faces of color 6 assign a sliding boundary

CALL GETFBR('6',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ialtyb(ifac) = igliss

enddo

! --- prescribe elsewhere a fixed boundary

CALL GETFBR( 'not (4 or 5 or 6)',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ialtyb(ifac) = ibfixe

enddo

!----
! FORMATS
!----

!----
! FIN
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
