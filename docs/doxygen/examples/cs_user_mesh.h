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
  \page cs_user_mesh Examples of mesh modifications

  \section cs_user_mesh_h_intro Introduction

  C user functions for optional modification of the mesh.
    These subroutines are called in all cases.

  Several functions are present in the file, each specific to different
    modification types.

  \section cs_user_mesh_h_cs_user_mesh_modifiy  General mesh modifications

  Mesh modifications not available through specialized functions
  should be defined in \ref cs_user_mesh_modify.

  \subsection cs_user_mesh_h_cs_user_mesh_modifiy_coords Coordinates modification

  For example, to modify coordinates, the
  following code can be added:

  \snippet cs_user_mesh-modify.c mesh_modify_coords

  \subsection cs_user_mesh_h_cs_user_mesh_modifiy_extrude_1 Boundary mesh patch extrusion

  It is possible to add cells by extruding selected boundary faces.
  A simplified usage is available through the \ref cs_mesh_extrude_constant,
  while extruding a mesh with vertex-local extrusion parameters is available
  through \ref cs_mesh_extrude.

  The example which follows illustrates the use of the simplified function.

  \snippet cs_user_mesh-modify.c mesh_modify_extrude_1

  \subsection cs_user_mesh_h_cs_user_mesh_modifiy_boundary_layer_1 Boundary layer insertion

  Boundary faces extrusion can also be used to insert boundary layer cells,
  by first shrinking the mesh around the selected zones so as to accomodate for
  added cells. The following example shows how this can be done for 2 different
  zones, using different parameters for each zone. Note that adjacent zones with
  a different number of inserted layers are possible. Also a positive zone
  thickness implies a fixed thicknesss, while a negative thickness is interpreted
  as a ratio relative to the mean adjacent cell size, allowing for automatic
  and adapted local thickness definition.

  \snippet cs_user_mesh-modify.c mesh_modify_boundary_layer

  Groups of cells, interior, and boundary faces may be created or
  modified, using the \ref cs_mesh_group_cells_set, \ref cs_mesh_group_i_faces_set,
  \ref cs_mesh_group_b_faces_set functions to assign a group to selected
  elements, removing all previous group attributions for those elements,
  and \ref cs_mesh_group_cells_add, \ref cs_mesh_group_i_faces_add,
  \ref cs_mesh_group_b_faces_add may be used to add those elements to a group
  while also keeping previous group information for those elements.

  The mesh is not marked as modified by default for this "light" modification,
  so the user may force this using a modification flag, as in the example here:

  \snippet cs_user_mesh-modify.c mesh_modify_groups_1

  The user can also add groups for boundary faces:

  \snippet cs_user_mesh-modify.c mesh_modify_groups_2

  \subsection cs_user_mesh_h_cs_user_mesh_modifiy_refine_1 Mesh refinement

  Cells may be refined automatically, using templates adapted to each cell
  type (see \ref cs_mesh_refine.c").

  This simply required selecting the cells to refine and whether the
  refinement should be partially propagated to neighboring cells
  (conforming) or not.

  The following code shows an example of mesh refinement for a given region.

  \snippet cs_user_mesh-modify.c mesh_modify_refine_1

  \subsection  cs_user_mesh_h_cs_user_mesh_input Mesh reading and modification

  The user function \ref cs_user_mesh_input allows a detailed selection of imported
  meshes read, reading files multiple times, applying geometric transformations,
  and renaming groups.

  The following code shows an example of mesh reading with no transformation.

  \snippet cs_user_mesh-input-save.c mesh_input_1

  A mesh can also be read while its groups are renamed, and its geometry
  transformed.

  \snippet cs_user_mesh-input-save.c mesh_input_2

  \subsection  cs_user_mesh_h_cs_user_mesh_save Mesh saving

  The user function \ref cs_user_mesh_save can enable or disable mesh saving.
  By default, mesh is saved when modified. The following code shows an example
  of disabled saving.

  \snippet cs_user_mesh-input-save.c mesh_save

  \section cs_user_mesh_h_cs_user_mesh_quality  Mesh quality modifications

  \subsection cs_user_mesh_h_cs_user_mesh_warping Mesh warping

  The \ref cs_user_mesh_warping function allows the user to cut the warped
  faces of his mesh using the \ref cs_mesh_warping_set_defaults function to
  define the maximum warped angle.

  \snippet cs_user_mesh-quality.c mesh_warping

  \subsection cs_user_mesh_h_cs_user_mesh_smoothing Mesh smoothing

  The smoothing utilities may be useful when the calculation mesh has local
  defects. The principle of smoothers is to mitigate the local defects by
  averaging the mesh quality. This procedure can help for calculation
  robustness or/and results quality. The user function \ref cs_user_mesh_smoothe
  allows to use different smoothing functions detailed below.

  The following code shows an example of use of the cs_mesh_smoother functions,
  \ref cs_mesh_smoother_fix_by_feature which fixes all boundary vertices that have
  one of their feature angles less than the maximum feature angle defined by the
  user and \ref cs_mesh_smoother_unwarp which reduces face warping in the calculation
  mesh.

  \snippet cs_user_mesh-quality.c mesh_smoothing

  \subsection cs_user_mesh_h_cs_user_mesh_tag_bad_cells Bad cells tagging

  Bad cells of a mesh can be tagged based on user-defined geometric criteria.
  The following example shows how to tag cells that have a volume below a
  certain value and then post-process the tagged cells. This is done using the
  \ref cs_user_mesh_bad_cells_tag function.

  \snippet cs_user_mesh-quality.c mesh_tag_bad_cells

  \section cs_user_mesh_h_cs_user_mesh_joining  Mesh joining

  \subsection cs_user_mesh_h_cs_user_mesh_add_simple_joining Simple mesh joining

  Conforming joining of possibly non-conforming meshes may be done by the
  \ref cs_user_join user function. For a simple mesh joining, the
  \ref cs_join_add subroutine is sufficient.

  \snippet cs_user_mesh-joining.c mesh_add_simple_joining

  \subsection cs_user_mesh_h_cs_user_mesh_add_advanced_joining Advanced mesh joining

  For a more complex mesh, or a mesh with thin walls which we want to avoid
  transforming into interior faces, the user can use the \ref cs_join_set_advanced_param
  function to define a specific mesh joining.

  \snippet cs_user_mesh-joining.c mesh_add_advanced_joining

  \section cs_user_mesh_h_cs_user_mesh_periodicity  Mesh periodicity

  Handling of periodicity can be performed with the \ref cs_user_periodicity function.

  \subsection cs_user_mesh_h_cs_user_mesh_translation_perio Periodicity of translation

  The following example illustrates the periodicity of translation case using the
  \ref cs_join_perio_add_translation subroutine.

  \snippet cs_user_mesh-periodicity.c mesh_periodicity_1

  \subsection cs_user_mesh_h_cs_user_mesh_rotation_perio Periodicity of rotation

  The following example illustrates the periodicity of rotation case using the
  \ref cs_join_perio_add_rotation subroutine.

  \snippet cs_user_mesh-periodicity.c mesh_periodicity_2

  \subsection cs_user_mesh_h_cs_user_mesh_mixed_perio General periodicity

  The following example illustrates a more general case of periodicity which
  combines different kinds of transformation. The function \ref cs_join_perio_add_mixed
  is used.

  \snippet cs_user_mesh-periodicity.c mesh_periodicity_3

  \subsection cs_user_mesh_h_cs_user_mesh_advanced_perio Periodicity advanced parameters

  As with the \ref cs_user_mesh_h_cs_user_mesh_add_advanced_joining subsection,
  a more complex periodicity can be defined using the \ref cs_join_set_advanced_param
  subroutine.

  \snippet cs_user_mesh-periodicity.c mesh_periodicity_4

  \section cs_user_mesh_h_cs_user_mesh_boundary  Mesh boundary insertion

  The user function \ref cs_user_mesh_boundary allows insertion of boundaries in
  the calculation mesh. This function transforms the selected interior faces
  into boundary faces, on which boundary conditions can (and must) be applied.
  Vertices are also split so as to be distinct on each side of the boundary.

  Boundaries can be directly inserted based on a selection of interior faces,
  such as shown here:

  \snippet cs_user_mesh-boundary.c mesh_thinwall

  Boundaries can also be inserted between a set of selected cells and the
  rest of the mesh. In this case, a mesh group name can be assigned to
  the added boundary faces.

  \snippet cs_user_mesh-boundary.c mesh_boundary_cells

*/
