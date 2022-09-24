/*============================================================================
 * code_saturne documentation page
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*-----------------------------------------------------------------------------*/



/*!

  \page cs_user_boundary_conditions_ale Examples of data settings for boundary conditions ale (cs_user_boundary_conditions_ale.f90)


 \brief User subroutine dedicated the use of ALE (Arbitrary Lagrangian
 Eulerian) Method:
  - Fills boundary conditions (\c ialtyb, \c icodcl, \c rcodcl) for mesh velocity.
  - This subroutine also enables one to fix displacement on nodes.

 \section intro_ale Introduction

 Here one defines boundary conditions on a per-face basis.

 Boundary faces may be identified using the \ref getfbr subroutine.
 The syntax of this subroutine is described in
 \c cs_user_boundary_conditions.f90 subroutine,
 but a more thorough description can be found in the user guide.

 Boundary conditions setup for standard variables (pressure, velocity,
 turbulence, scalars) is described precisely in
 \c cs_user_boundary_conditions.f90 subroutine.

 Detailed explanation will be found in the theory guide.

 \section bc_types_ale Boundary condition types

 Boundary conditions may be assigned in two ways.


 \subsection std_bcs_ale For "standard" boundary conditions

 One defines a code in the \c ialtyb array (of dimensions number of
 boundary faces). The available codes are:

  - \c ialtyb(ifac) = \c ibfixe: the face \c ifac is considered to be motionless.
           A zero Dirichlet boundary condition is automatically imposed on mesh
           velocity. Moreover the displacement of corresponding nodes will
           automatically be set to 0 (for further information please
           read the paragraph dedicated to the description of \c impale array in the
           cs_user_boundary_conditions_ale.f90 subroutine), unless the USER has modified the condition of
           at least one mesh velocity component (modification of \c icodcl array,
           please read the following paragraph \ref non_std_bc_ale)

  - \c ialtyb(ifac) = \c igliss: The mesh slides on corresponding face \c ifac.
           The normal component of mesh velocity is automatically set to 0.
           A homogeneous Neumann condition is automatically prescribed for the
           other components, as it's the case for 'Symmetry' fluid condition.

  - \c ialtyb(ifac) = \c ivimpo: the mesh velocity is imposed on face \c ifac. Thus,
           the users needs to specify the mesh velocity values filling \c rcodcl
           arrays as follows:
            - \c rcodcl(ifac,iuma,1) = mesh velocity in 'x' direction
            - \c rcodcl(ifac,ivma,1) = mesh velocity in 'y' direction
            - \c rcodcl(ifac,iwma,1) = mesh velocity in 'z' direction
            .
           Components of \c rcodcl(.,i.ma,1) arrays that are not specified by user
           will automatically be set to 0, meaning that user only needs to specify
           non zero mesh velocity components.


 \subsection non_std_bc_ale For "non-standard" conditions

 Other than (fixed boundary, sliding mesh boundary, fixed velocity), one
 defines for each face and each component \c IVAR = IUMA, IVMA, IWMA:
  - a code
    - \c icodcl(ifac, ivar)
  - three real values:
    - \c rcodcl(ifac, ivar, 1)
    - \c rcodcl(ifac, ivar, 2)
    - \c rcodcl(ifac, ivar, 3)

 The value of \c icodcl is taken from the following:
  - 1: Dirichlet
  - 3: Neumann
  - 4: Symmetry

 The values of the 3 \c rcodcl components are:
  - \c rcodcl(ifac, ivar, 1):
     Dirichlet for the variable if \c icodcl(ifac, ivar) = 1
     The dimension of \c rcodcl(ifac, ivar, 1) is in \f$m \cdot s^{-1}\f$
  - \c rcodcl(ifac, ivar, 2):
    "exterior" exchange coefficient (between the prescribed value
                      and the value at the domain boundary),\n
                      rinfin = infinite by default
    \f$  rcodcl(ifac,ivar,2) =  \dfrac{VISCMA}{d} \f$
          (d has the dimension of a distance in \f$m\f$, \f$VISCMA\f$ stands for
          the mesh viscosity)


 \remark
  - The definition of \c rcodcl(.,.,2) is based on the manner
            other standard variables are managed in the same case.
            This type of boundary condition appears nonsense
            concerning mesh in that context.

      - \c rcodcl(ifac,ivar,3) :
    Flux density (in \f$kg \cdot m \cdot s^2\f$) = J if icodcl(ifac, ivar) = 3
                 (<0 if gain, n outwards-facing normal)
    \f$ rcodcl(ifac,ivar,3) = -(VISCMA) \grad {Um}.\vect{n} \f$
              \f$(Um\f$ represents mesh velocity)


  - The definition of condition \c rcodcl(ifac,ivar,3)
            is based on the manner other standard variables are
            managed in the same case.
            \c rcodcl(.,.,3) = 0.d0 enables one to specify a homogeneous
            Neuman condition on mesh velocity. Any other value will be
            physically nonsense in that context.

  - If the user assigns a value to \c ialtyb equal to \c ibfixe, \c igliss,
 or \c ivimpo and does not modify \c icodcl (zero value by
 default), \c ialtyb will define the boundary condition type.

 To the contrary, if the user prescribes \c icodcl(ifac, ivar) (nonzero),
 the values assigned to rcodcl will be used for the considered face
 and variable (if rcodcl values are not set, the default values will
 be used for the face and variable, so:
                         - \c rcodcl(ifac, ivar, 1) = 0.d0
                         - \c rcodcl(ifac, ivar, 2) = rinfin
                         - \c rcodcl(ifac, ivar, 3) = 0.d0)

 If the user decides to prescribe his own non-standard boundary conditions,
 it will be necessary to assign values to \c icodcl AND to \c rcodcl for ALL
 mesh velocity components. Thus, the user does not need to assign values
 to \c ialtyb for each associated face, as it will not be taken into account
 in the code.


 \subsection cons_rul_ale Consistency rules

 A consistency rules between \c icodcl codes for variables with
 non-standard boundary conditions:
  - If a symetry code (\c icodcl = 4) is imposed for one mesh velocity
    component, one must have the same condition for all other mesh
    velocity components.


 \subsection fix_nod_ale Fixed displacement on nodes

 For a better precision concerning mesh displacement, one can also assign values
 of displacement to certain internal and/or boundary nodes. Thus, one
 need to fill \c disale and \c impale arrays :
  - \c disale(1,inod) = displacement of node inod in 'x' direction
  - \c disale(2,inod) = displacement of node inod in 'y' direction
  - \c disale(3,inod) = displacement of node inod in 'z' direction
 This array is defined as the total displacement of the node compared
 its initial position in initial mesh.
 \c impale(inod) = 1 indicates that the displacement of node inod is imposed.


 \note \c impale array is initialized to the value of 0; if its value
       is not modified, corresponding value in \c DEPALE array will not be
       taken into account.

 During mesh's geometry re-calculation at each time step, the position of the
 nodes, which displacement is fixed (\c i.e. \c impale=1), is not calculated
 using the value of mesh velocity at the center of corresponding cell, but
 directly filled using the values of \c disale.

 If the displacement is fixed for all nodes of a boundary face it's not
 necessary to prescribe boundary conditions at this face on mesh velocity.
 \c icodcl and \c rcodcl values will be overwritten:
  - \c icodcl is automatically set to 1 (Dirichlet)
  - \c rcodcl value will be automatically set to face's mean mesh velocity
    value, that is calculated using \c DEPALE array.

 If a fixed boundary condition (\c ialtyb(ifac)=ibfixe) is imposed to the face
 \c ifac, the displacement of each node \c inod belonging to \c ifac is considered
 to be fixed, meaning that \c impale(inod) = 1 and \c disale(.,inod) = 0.d0.


 \subsubsection nod_des_ale Description of nodes

 \c nnod gives the total (internal and boundary) number of nodes.
 Vertices coordinates are given by \c xyznod(3, nnod) array. This table is
 updated at each time step of the calculation.
 \c xyzno0(3,nnod) gives the coordinates of initial mesh at the beginning
 of the calculation.

 The faces - nodes connectivity is stored by means of four integer arrays :
 \c ipnfac, \c nodfac, \c ipnfbr, \c nodfbr.

 \c nodfac(nodfbr) stores sequentially the index-numbers of the nodes of each
 internal (boundary) face.
 \c ipnfac(ipnfbr) gives the position of the first node of each internal
 (boundary) face in the array \c nodfac(nodfbr).

 For example, in order to get all nodes of internal face \c ifac, one can
 use the following loop:

 \code
 do ii = ipnfac(ifac), ipnfac(ifac+1)-1 !! index number of nodfac array
                                        !! corresponding to ifac

   inod = nodfac(ii)                    !! index-number iith node of face ifac.
   !! ...
 enddo
 \endcode


 \subsection flui_bc_ale Influence on boundary conditions related to fluid velocity

 The effect of fluid velocity and ALE modeling on boundary faces that
 are declared as walls (\c itypfb = \c iparoi or \c iparug) really depends on
 the physical nature of this interface.

 Indeed when studying an immersed structure the motion of corresponding
 boundary faces is the one of the structure, meaning that it leads to
 fluid motion. On the other hand when studying a piston the motion of vertices
 belonging to lateral boundaries has no physical meaning therefore it has
 no influence on fluid motion.

 Whatever the case, mesh velocity component that is normal to the boundary
 face is always taken into account
 (\f$ \vect{u}_{fluid} \cdot \vect{n} = \vect{w}_{mesh} \cdot \vect{n} \f$).

 The modeling of tangential mesh velocity component differs from one case
 to another.

 The influence of mesh velocity on boundary conditions for fluid modeling is
 managed and modeled in \c code_saturne as follows:
  - If \c ialtyb(ifac) = ibfixe: mesh velocity equals 0. (In case of 'fluid sliding
  wall' modeling corresponding condition will be specified in code_saturne
  Interface or in cs_user_boundary_conditions.f90 subroutine.)
  - If \c ialtyb(ifac) = ivimpo: tangential mesh velocity is modeled as a sliding
  wall velocity in fluid boundary conditions unless a value for fluid sliding
  wall velocity has been specified by USER in code_saturne Interface
  or in cs_user_boundary_conditions.f90 subroutine.
  - If \c ialtyb(ifac) = igliss: tangential mesh velocity is not taken into account
  in fluid boundary conditions (In case of 'fluid sliding wall' modeling
  corresponding condition will be specified in code_saturne Interface
  or in cs_user_boundary_conditions.f90 subroutine.)
  - If \c impale(inod) = 1 for all vertices of a boundary face: tangential mesh
  velocity value that has been derived from nodes displacement is modeled as a
  sliding wall velocity in fluid boundary conditions unless a value for fluid
  sliding wall velocity has been specified by USER in code_saturne Interface or
  in cs_user_boundary_conditions.f90 subroutine.

 Note that mesh velocity has no influence on modeling of
 boundary faces with 'inlet' or 'free outlet' fluid boundary condition.

 For "non standard" conditions USER has to manage the influence of boundary
 conditions for ALE method (i.e. mesh velocity) on the ones for Navier Stokes
 equations(i.e. fluid velocity). (Note that fluid boundary conditions can be
 specified in this subroutine.)


\subsubsection cell_id_ale Cells identification

 Cells may be identified using the getcel subroutine.
 The syntax of this subroutine is described in the
 cs_user_boundary_conditions.f90 subroutine,
 but a more thorough description can be found in the user guide.

 \subsubsection fac_id_ale Faces identification

 Faces may be identified using the \ref getfbr subroutine.
 The syntax of this subroutine is described in the
 cs_user_boundary_conditions.f90 subroutine,
 but a more thorough description can be found in the user guide.

  \section example_ale Example of boundary conditions ale

Here is the list of examples:

  \subpage example_ale2
*/
//______________________________________________________________________________________
/*!

  \page example_ale2 Examples of boundary conditions ale

  \section loc_var_ale Local variables

  \snippet  cs_user_boundary_conditions_ale-base.f90 loc_var

  \section init_fin_ale Initialization and finalization

  The following initialization block needs to be added for the following examples:
   \snippet cs_user_boundary_conditions_ale-base.f90 allocate_ale

At the end of the subroutine, it is recommended to deallocate the work array:
   \snippet cs_user_boundary_conditions_ale-base.f90 deallocate_ale

In theory Fortran 95 deallocates locally-allocated arrays automatically, but deallocating arrays in a symetric manner to their allocation is good pratice, and avoids using a different logic C and Fortran.

  \section assign_ale Assign boundary conditions to boundary faces

One may use selection criteria to filter boundary case subsets.\n Loop on faces from a subset. \n Set the boundary condition for each face.

   \subsection calcualtion_ale Calculation of displacement at current time step

   \snippet  cs_user_boundary_conditions_ale-base.f90 calcul

   \subsection example1_ale Example 1
 Example : For boundary faces of color 4 assign a fixed velocity

   \snippet  cs_user_boundary_conditions_ale-base.f90 example_1

   \subsection example2_ale Example 2
Example: For boundary face of color 5 assign a fixed displacement on nodes

   \snippet  cs_user_boundary_conditions_ale-base.f90 example_2

   \subsection  example3_ale Example 3
Example : For boundary faces of color 6 assign a sliding boundary

   \snippet  cs_user_boundary_conditions_ale-base.f90  example_3

   \subsection  example4_ale Example 4
Example : Prescribe elsewhere a fixed boundary

   \snippet  cs_user_boundary_conditions_ale-base.f90 example_4


*/
