/*============================================================================
 * Doxygen documentation for specific keywords
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

/*----------------------------------------------------------------------------*/

#include "cs_function_default.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file function_defaults.h
        Predefined function objects
*/

/*----------------------------------------------------------------------------*/

/*!
 * \defgroup function_object_defaults Predefined function objects

  Note that a function object of a given name may always be accessed
  using \ref cs_function_by_name.

 */
/*!@{*/

/*!
  \var cells_mpi_rank_id

  MPI rank id to which a cell is assigned.

  Accessed or activated using:
  \code{.c}
  cs_function_define_mpi_rank_id(CS_MESH_LOCATION_CELLS);
  \endcode
*/
char *cells_mpi_rank_id;

/*!
  \var interior_faces_mpi_rank_id

  MPI rank id to which interior faces are assigned. Interior faces
  at parallel boundaries are seen from 2 domains, but may be assigned
  specifically to one domain when range sets (\ref cs_range_set_t) are used.
  If no range set has already been assigned to interior faces, a default one
  (ignoring periodicity) is used. Otherwise (such as when using CDO face-based
  schemes), the one actually used for computation is used.

  Accessed or activated using:
  \code{.c}
  cs_function_define_mpi_rank_id(CS_MESH_LOCATION_INTERIOR_FACES);
  \endcode
*/
char *interior_faces_mpi_rank_id;

/*!
  \var boundary_faces_mpi_rank_id

  MPI rank id to which boundary faces are assigned. This rank is always the
  same as the adjacent cell.

  Accessed or activated using:
  \code{.c}
  cs_function_define_mpi_rank_id(CS_MESH_LOCATION_BOUNDARY_FACES);
  \endcode
*/
char *boundary_faces_mpi_rank_id;

/*!
  \var vertices_mpi_rank_id

  MPI rank id to which vertices are assigned. Vertices at parallel boundaries
  are seen from multiple, but may be assigned specifically to one domain when
  range sets (\ref cs_range_set_t) are used.
  If no range set has already been assigned to vertices, a default one
  (ignoring periodicity) is used. Otherwise (such as when using CDO vertex-based
  schemes), the one actually used for computation is used.

  Accessed or activated using:
  \code{.c}
  cs_function_define_mpi_rank_id(CS_MESH_LOCATION_VERTICES);
  \endcode
*/
char *vertices_mpi_rank_id;

/*!
  \var cells_r_gen

  Refinement generation of mesh cells. This is considered to be the highest
  refinement level of adjacent interior faces.

  Accessed or activated using:
  \code{.c}
  cs_function_define_refinement_generation(CS_MESH_LOCATION_CELLS);
  \endcode
*/
char *cells_r_gen;

/*!
  \var interior_faces_r_gen

  Refinement generation of interior faces. This is determined and stored
  when refining or coarsening a mesh. A refined face's sub-faces have the
  same refinement generation as their parent. Faces added in the course
  of refinement (i.e. faces separating sub-cells) have a refinement level
  one higher than the parent cell.

  Accessed or activated using:
  \code{.c}
  cs_function_define_refinement_generation(CS_MESH_LOCATION_INTERIOR_FACES);
  \endcode
*/
char *interior_faces_r_gen;

/*!
  \var boundary_faces_r_gen

  Refinement generation of boundary faces. This should always be that of
  the adjacent cell.

  Accessed or activated using:
  \code{.c}
  cs_function_define_refinement_generation(CS_MESH_LOCATION_BOUNDARY_FACES);
  \endcode
*/
char *boundary_faces_r_gen;

/*!
  \var boundary_zone_class_id

  Optional boundary face class or zone ids. If no face classes have been
  defined by \ref cs_boundary_zone_face_class_id, the boundary face zone
  id is used instead.

  Activated by default.
*/
char *boundary_zone_class_id;

/*!
  \var vertices_r_gen

  Refinement generation of vertices. This is determined and stored
  when refining or coarsening a mesh.

  Accessed or activated using:
  \code{.c}
  cs_function_define_refinement_generation(CS_MESH_LOCATION_VERTICES);
  \endcode
*/
char *vertices_r_gen;

/*!
  \var "boundary_stress"

  Stress exerted by fluid forces over boundary. Activating this function
  also leads to the creation of the "boundary_forces" field, used to store
  and update the (extensive) boundary forces.

  Accessed or activated using:
  \code{.c}
  cs_function_define_boundary_stress();
  \endcode
*/
char *boundary_stress;

/*!
  \var "boundary_stress_normal"

  Normal component of the stress exerted by fluid forces over boundary.
  Activating this function also leads to the creation of the "boundary_forces"
  field, used to store and update the (extensive) boundary forces.

  Accessed or activated using:
  \code{.c}
  cs_function_define_boundary_stress_normal();
  \endcode
*/
char *boundary_stress_normal;

/*!
  \var "boundary_stress_tangential"

  Tangential component of the stress exerted by fluid forces over boundary.
  Activating this function also leads to the creation of the "boundary_forces"
  field, used to store and update the (extensive) boundary forces.

  Accessed or activated using:
  \code{.c}
  cs_function_define_boundary_stress_tangential();
  \endcode
*/
char *boundary_stress_tangential;

/*!
  \var "boundary_thermal_flux"

  Thermal flux density over the boundary. Incoming thermal flux leads to a
  positive sign; outgoing thermal flux to a negative sign.

  Accessed or activated using:
  \code{.c}
  cs_function_define_boundary_thermal_flux();
  \endcode
*/
char *boundary_thermal_flux;

/*!
  \var "boundary_nusselt"

  Boundary layer Nusselt number. This is baed only on the near-boundary
  temperature value and heat flux, so does not require a bulk temperature value
  but is not a "true" Nusselt number computation, as the fluid value used may
  be sensitive to mesh refinement.

  Accessed or activated using:
  \code{.c}
  cs_function_define_boundary_nusselt();
  \endcode
*/
char *boundary_thermal_flux;

/*!
  \var "relative_pressure"

  Relative pressure (at cells) for turbomachinery computations.

  Automatically activated for turbomachinery computations.
*/
char *relative_pressure;

/*!
  \var "relative_velocity"

  Relative velocity (at cells) for turbomachinery computations.

  Automatically activated for turbomachinery computations.
*/
char *relative_velocity;

/*!
  \var "absolute_pressure"

  Absolute pressure (at cells) with Corilolis forces.

  Automatically activated when Coriolis forces are present.
*/
char *absolute_pressure;

/*!
  \var "absolute_velocity"

  Absolute velocity (at cells) with Coriolis forces.

  Automatically activated when Coriolis forces are present.
*/
char *absolute_velocity;

/*!
  \var "elec_pot_gradient_im"

  For Joule Heating by direct conduction, gradient of the imaginary component
  of the potential.

  Automatically activated when Joule Heating by direct conduction is active
  cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 2 or 4.
*/
char *elec_pot_gradient_im;

/*!
  \var "elec_current_gradient_im"

  For Joule Heating by direct conduction, imaginary component of the
  current density

  Automatically activated when Joule Heating by direct conduction is active
  cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 2 or 4.
*/
char *elec_current_gradient_im;

/*!
  \var "elec_pot_module"

  For Joule Heating by direct conduction, module of the complexe potential.

  Automatically activated when Joule Heating by direct conduction is active
  cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 4.
*/
char *elec_pot_module;

/*!
  \var "elec_pot_arg"

  For Joule Heating by direct conduction, argument of the complexe potential.

  Automatically activated when Joule Heating by direct conduction is active
  cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 4.
*/
char *elec_pot_arg;

/*!
  \var "q_criterion"

  Compute the Q criterion from Hunt et. al.

  \f[
    Q = \tens{\Omega}:\tens{\Omega} -
    \deviator{ \left(\tens{S} \right)}:\deviator{ \left(\tens{S} \right)}
  \f]
  where \f$\tens{\Omega}\f$ is the vorticity tensor and
  \f$\deviator{ \left(\tens{S} \right)}\f$ the deviatoric of the rate of strain
  tensor.

  Accessed or activated using:
  \code{.c}
  cs_function_define_q_criterion();
  \endcode
*/
char *q_criterion;

/*!@}*/

/*----------------------------------------------------------------------------*/

