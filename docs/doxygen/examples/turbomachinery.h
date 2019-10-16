/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
  \page turbomachinery Turbomachinery module usage

  \section tbm_intro Introduction

  Two classical models are available in \c Code_Saturne for rotor/stator
  interactions modelling in turbomachinery computations: the steady approach
  which is based on the so-called <I>Frozen Rotor</I> modelling and the
  <I>transient rotor/stator</I> approach which is based on a sliding mesh
  technique.

  \warning This page describes these functionalities based on a single
  \c Code_Saturne computation. An alternative rotor/stator coupling based on the
  code/code coupling of \c Code_Saturne with itself is still possible
  (and briefly described in this page) but it is not recommended
  (will be highly modified soon).

  \section tbm_mesh Mesh requirements

  \subsection tbm_mesh_perio Periodicity

  The rotational periodicity treatment is possible only in <I>Frozen Rotor</I>.
  However, the interface plane between rotor and stator must match in the
  azimutal \f$\theta\f$ direction:
  \f$\theta_{\text{min}}^{\text{rotor}}(z)=\theta_{\text{min}}^{\text{stator}}(z),\quad\theta_{\text{max}}^{\text{rotor}}(z)=\theta_{\text{max}}^{\text{stator}}(z)\f$
  for all <I>z</I> through the rotation axis direction.

  \subsection tbm_mesh_interface Rotor/stator interface

   - <I>Transient rotor/stator</I>: in the input mesh(es), the interface
   between rotor and stator domains has to be composed of
   <I>boundary faces</I>. Then the inteface boundary faces are joined during
   the computation and transformed into internal faces, as the preprocessing
   step can do in usual calculations. The classical way to fulfill this
   requirement is to provide <I>separate meshes</I> for each rotor or stator
   domains, but placing both sections in the same mesh is allowed, as long
   as they are topologically disjoint.

   - <I>Frozen Rotor</I>: the interface can be composed of boundary
   faces (then the inteface boundary faces are joint at the beginning of
   the computation and turn into internal faces) or internal faces.

  \section tbm_data_setting Data setting

  The data setting is made through the GUI or through the
  cs_user_turbomachinery.c source file. Through the latter, the type of
  turbomachinery model can be first set in \ref cs_user_turbomachinery function:

  \snippet cs_user_turbomachinery.c user_tbm_set_model

  Then, the rotor region and the rotor/stator interface boundaries has to be
  specified in \ref cs_user_turbomachinery_rotor function. The rotor
  region is first specified, as follows:

  \snippet cs_user_turbomachinery.c user_tbm_set_rotor

  In this example, rotation_velocity refers to the rotation angular velocity
  of the rotor, in rad/s. The rotation axis can be normalized (it is eventually
  normalized afterwards by the code). Note that if a rotor is added,
  the code appends an instance of \ref cs_rotation_t
  structure in the global list of rotors (see \ref tbm_user_basic_op
  section in the following).

  Then the rotor/stator interface boundary is specified. This step is mandatory
  only for the CS_TURBOMACHINERY_TRANSIENT model. For the
  CS_TURBOMACHINERY_FROZEN model, this step is required only if the interface
  is made of boundary faces (see \ref tbm_mesh_interface subsection above).

  \snippet cs_user_turbomachinery.c user_tbm_set_interface

  The rotation velocity can be modified during the calculation. The following example
  shows how to set a linearly increasing rotation velocity in
  \ref cs_user_turbomachinery_set_rotation_velocity function:

  \snippet cs_user_turbomachinery.c user_tbm_set_linear_rotation_velocity

  \section tbm_user_bpg User recomandations

  \subsection tbm_user_bpg_mesh Meshing recomandations at interface

  As mentioned above, when a rotor/stator inteface boundary exists (in particular
  for the CS_TURBOMACHINERY_TRANSIENT model), it is then joined by the solver
  during the computation. It is thus important to be aware that the success of
  a joining operation is strongly dependant on the <I>quality of the mesh at the
  interface</I>. More precisely, the refinement must be as similar as possible
  on both sides of the interface. Moreover, it is reminded that the tolerance
  parameter of a joining is a fraction of the shortest edge linked with a vertex
  to join. Consequently, cells where the refinement in the azimutal \f$\theta\f$
  direction is much coarser than those in one of the two others can also lead
  to a joining failure. In particular, the user should be wary of boundary
  layers.

  \remark
  - The tolerance parameter of a joining should always be lower than 0.5.
  - If the meshes at both sides of the interface are very different such that
  the joining fails, advanced joining parameters are available. However,
  modifying the mesh is more likely to succeed. The introduction of some kind of
  buffer cells layer on both sides of the interface is strongly recommended.
  Ideally, each of the two layers should have the same refinement and a constant
  azimutal step (this latter recommandation is relevant only for
  CS_TURBOMACHINERY_TRANSIENT model).

  \subsection tbm_user_bpg_cpl Alternative rotor/stator coupling

  If the meshes on both sides of the interface are very different and can not
  be modified, the user should turn to a rotor/stator coupling based on the
  code/code coupling of \c Code_Saturne with itself. This simply requires
  replacing the definition of a face joining for an interface (using the
  GUI or a call to \ref cs_turbomachinery_join_add by a call to
  \ref cs_turbomachinery_coupling_add, as in the following example:

  \snippet cs_user_turbomachinery.c user_tbm_set_coupling

  If a coupling is defined and the rotation is not null in at least one of the
  cases, \ref cstphy::icorio "icorio" determines the type of rotor/stator
  coupling: icorio = 1 for a <I>Frozen Rotor</I> computation and icorio = 0 for
  a sliding mesh approach.

  \warning Contrary to the mesh joining approach, the code/code coupling
  approach is not fully conservative.

  \section tbm_user_op User operations

  \subsection tbm_user_basic_op Basic operations

  A specific rotor is identified by an integer value. The "default rotor" is
  actually the stator that is the fixed part of the domain, with identifier 0.
  Then the identifiers of the rotors are incremented as the user add rotors
  (\ref cs_turbomachinery_add_rotor, see above).

  Once the rotor identifier is known, one can acces to its parameters:

  \snippet cs_user_extra_operations-turbomachinery.c extra_tbm_get_rotor_info

  Importantly, one access to the rotor associated to each cell of the domain
  thanks to a rotor identifier list, like in this example:

  \snippet cs_user_extra_operations-turbomachinery.c extra_tbm_get_rotor

  \subsection tbm_user_post Turbomachinery dedicated postprocessing functions

  Useful postprocessing functions relative to the machinery characteristics are
  available: postprocessing of the couple on the rotor walls
  (\ref cs_post_moment_of_force) and postprocessing of the head generated by the
  machinery (\ref cs_post_turbomachinery_head). In the latter one, the mesh
  locations (CS_MESH_LOCATION_CELLS, CS_MESH_LOCATION_BOUNDARY_FACES or
  CS_MESH_LOCATION_INTERNAL_FACES) associated with the given selection
  criteria must be specified.

  \snippet cs_user_extra_operations-turbomachinery.c extra_tbm_post_util

  \subsection tbm_user_fortran Fortran naming

  Useful fortran variables and/or functions specific to turbomachinery computations
  are contained int the fortran module \ref at_turbomachinery "turbomachinery"
  and \ref at_rotation "rotation".

*/
