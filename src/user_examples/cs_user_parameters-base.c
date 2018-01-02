/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_convection_diffusion.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_util.h"
#include "cs_grid.h"
#include "cs_internal_coupling.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_selector.h"
#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-base.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /* If the GUI is used, user fields should preferably be defined with the GUI,
     so that associated numerical options, boundary conditions, initializations,
     and such may also be defined using the GUI. */

  if (cs_gui_file_is_loaded())
    return;

  /*--------------------------------------------------------------------------*/

  /* Example: add 2 scalar variables ("species" in the GUI nomenclature).
   *
   * Note that at this (very early) stage of the data setup, fields are
   * not defined yet. Associated fields will be defined later (after
   * model-defined fields) in the same order as that used here, and
   * after user-defined variables defined throught the GUI, if used.
   *
   * Currently, only 1d (scalar) fields are handled.
   *
   * parameters for cs_parameters_add_variable():
   *   name             <-- name of variable and associated field
   *   dim              <-- variable dimension
   */

  cs_parameters_add_variable("species_1", 1);
  cs_parameters_add_variable("tracer", 1);

  /*--------------------------------------------------------------------------*/

  /* Example: add the variance of a user variable.
   *
   * parameters for cs_parameters_add_variable_variance():
   *   name          <-- name of variance and associated field
   *   variable_name <-- name of associated variable
   */

  cs_parameters_add_variable_variance("variance_1",
                                      "species_1");

  /*--------------------------------------------------------------------------*/

  /* Example: add a user property defined on boundary faces.
   *
   * parameters for cs_parameters_add_property():
   *   name        <-- name of property and associated field
   *   dim         <-- property dimension
   *   location_id <-- id of associated mesh location, which must be one of:
   *                     CS_MESH_LOCATION_CELLS
   *                     CS_MESH_LOCATION_INTERIOR_FACES
   *                     CS_MESH_LOCATION_BOUNDARY_FACES
   *                     CS_MESH_LOCATION_VERTICES
   */

  cs_parameters_add_property("user_b_property_1",
                             1,
                             CS_MESH_LOCATION_BOUNDARY_FACES);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(void)
{
  /* Example: choose a limiter for a given scalar */
  /*----------------------------------------------*/

  /*! [param_var_limiter_choice] */

  /* retrieve scalar field by its name */
  cs_field_t *sca1 = cs_field_by_name("scalar1");

  /* isstpc:
     0: swich on the slope test
     1: swich off the slope test (default)
     2: continuous limiter ensuring boundedness (beta limiter)
     3: NVD/TVD Scheme */

  cs_var_cal_opt_t vcopt;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  cs_field_get_key_struct(sca1, key_cal_opt_id, &vcopt);
  vcopt.isstpc = 3;
  cs_field_set_key_struct(sca1, key_cal_opt_id, &vcopt);

  /* Min/Max limiter or NVD/TVD limiters
     then "limiter_choice" keyword must be set:
     0: Gamma
     1: SMART
     2: CUBISTA
     3: SUPERBEE
     4: MUSCL
     5: MINMOD
     6: CLAM
     7: STOIC
     8: OSHER
     9: WASEB
     --- VOF scheme ---
     10: M-HRIC
     11: M-CICSAM       */

  int key_lim_id = cs_field_key_id("limiter_choice");
  cs_field_set_key_int(sca1, key_lim_id, CS_NVD_SUPERBEE);

  /*! [param_var_limiter_choice] */

  /* Example: add boundary values for all scalars */
  /*----------------------------------------------*/

  /*! [param_var_boundary_vals_1] */

  int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t  *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE)
      cs_parameters_add_boundary_values(f);

  }

  /*! [param_var_boundary_vals_1] */

  /* Example: post-process clippings for Rij tensor */
  /*------------------------------------------------*/

  /*! [param_var_rij_clipping] */

  cs_field_set_key_int(CS_F_(rij), cs_field_key_id("clipping_id"), 1);
  cs_field_set_key_int(CS_F_(eps), cs_field_key_id("clipping_id"), 1);

  /*! [param_var_rij_clipping */

  /* Example: post-process the Q-critertion on the whole domain mesh */
  /*-----------------------------------------------------------------*/

  /*! [param_var_q_criterion] */

  cs_glob_post_util_flag[CS_POST_UTIL_Q_CRITERION] = 1;

  /*! [param_var_q_criterion] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volumes as internal coupling zones.
 *
 * These zones will be separated from the rest of the domain using automatically
 * defined thin walls.
 *
 * \param[in, out]  mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_add_volumes(cs_mesh_t *mesh)
{
  /* Example: Internal coupling */
  /*--------------------------- */

  /*! [param_internal_coupling] */

  cs_internal_coupling_add_volume(mesh,
                                  "x<.5"); /* Solid volume criterion */

  /*! [param_internal_coupling] */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volumesi from separated meshes as internal coupling zones.
 *
 * These zones must be disjoint and the face selection criterion must be
 * specified.
 *
 * \param[in, out]  mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_from_disjoint_meshes(cs_mesh_t *mesh)
{
  /* Example: Internal coupling */
  /*--------------------------- */

  /*! [param_internal_coupling] */

  cs_internal_coupling_add(mesh,
                           "solid_volume_criterion",
                           "interface_criterion");

  /*! [param_internal_coupling] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling options.
 *
 * Options are usually defined using cs_internal_coupling_add_entity.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling(void)
{
  /* Example: Internal coupling */
  /*--------------------------- */

  /*! [param_internal_coupling] */

  int f_id = cs_field_id_by_name("scalar1");

  cs_internal_coupling_add_entity(f_id);  /* Field to be coupled */


  /*! [param_internal_coupling] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
