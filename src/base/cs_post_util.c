/*============================================================================
 * Postprocessing utility functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_interface.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_join.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_math.h"
#include "cs_matrix_default.h"
#include "cs_mesh.h"
#include "cs_mesh_coherency.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_renumber.h"
#include "cs_rotation.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_turbomachinery.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post_util.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the head of a turbomachinery (total pressure increase)
 *
 * \param[in]   criteria_in   selection criteria of turbomachinery suction
 * \param[in]   location_in   mesh location of turbomachinery suction
 * \param[in]   criteria_out  selection criteria of turbomachinery discharge
 * \param[in]   location_out  mesh location of turbomachinery discharge
 *
 * \return turbomachinery head
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_turbomachinery_head(const char               *criteria_in,
                            cs_mesh_location_type_t   location_in,
                            const char               *criteria_out,
                            cs_mesh_location_type_t   location_out)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;

  cs_real_t *total_pressure = cs_field_by_name("total_pressure")->val;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;
  cs_real_t *density = CS_F_(rho)->val;

  cs_real_t pabs_in = 0.;
  cs_real_t sum_in = 0.;
  cs_real_t pabs_out = 0.;
  cs_real_t sum_out = 0.;

  for (int _n = 0; _n < 2; _n++) {

    cs_lnum_t n_elts = 0;
    cs_lnum_t *elt_list = NULL;
    cs_real_t pabs = 0.;
    cs_real_t sum = 0.;

    cs_mesh_location_type_t location;
    const char *criteria = NULL;

    if (_n == 0) {
      location = location_in;
      criteria = criteria_in;
    } else {
      location = location_out;
      criteria = criteria_out;
    }

    switch(location) {
    case CS_MESH_LOCATION_CELLS:

      BFT_MALLOC(elt_list, mesh->n_cells, cs_lnum_t);
      cs_selector_get_cell_num_list(criteria, &n_elts, elt_list);

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t cell_id = elt_list[i];
        cs_real_t weight = mesh_quantities->cell_vol[cell_id];
        pabs += weight*(total_pressure[cell_id] + 0.5*density[cell_id]*
                        cs_math_3_square_norm(vel[cell_id]));
        sum += weight;
      }
      BFT_FREE(elt_list);
      break;

    case CS_MESH_LOCATION_BOUNDARY_FACES:

      BFT_MALLOC(elt_list, mesh->n_b_faces, cs_lnum_t);
      cs_selector_get_b_face_num_list(criteria, &n_elts, elt_list);

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t face_id = elt_list[i];
        cs_lnum_t cell_id = mesh->b_face_cells[face_id];
        cs_real_t surf = mesh_quantities->b_face_surf[face_id];
        pabs += surf*(total_pressure[cell_id] + 0.5*density[cell_id]
                      *cs_math_3_square_norm(vel[cell_id]));
        sum += surf;
      }
      BFT_FREE(elt_list);
      break;

    case CS_MESH_LOCATION_INTERIOR_FACES:

      BFT_MALLOC(elt_list, mesh->n_i_faces, cs_lnum_t);
      cs_selector_get_i_face_num_list(criteria, &n_elts, elt_list);

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t face_id = elt_list[i];
        cs_lnum_t c_i = mesh->i_face_cells[face_id][0];
        cs_lnum_t c_j = mesh->i_face_cells[face_id][1];
        cs_real_t w = mesh_quantities->i_face_surf[face_id];

        cs_real_t pt = w*total_pressure[c_i] + (1.-w)*total_pressure[c_j];
        cs_real_t r = w*density[c_i] + (1.-w)*density[c_j];
        cs_real_3_t v = {w*vel[c_i][0] + (1.-w)*vel[c_j][0],
                         w*vel[c_i][1] + (1.-w)*vel[c_j][1],
                         w*vel[c_i][2] + (1.-w)*vel[c_j][2]};
        pabs += w*(pt + 0.5*r*cs_math_3_square_norm(v));
        sum += w;
      }
      BFT_FREE(elt_list);
      break;

    default:
      pabs = 0.;
      sum = 1.;
      bft_printf
        (_("Warning: while post-processing the turbomachinery head.\n"
           "         Mesh location %d is not supported, so the computed head\n"
           "         is erroneous.\n"
           "         The %s parameters should be checked.\n"),
           location, __func__);
      break;
    }

    if (_n == 0) {
      pabs_in = pabs;
      sum_in = sum;
    } else {
      pabs_out = pabs;
      sum_out = sum;
    }

  }

  double _s[4] = {pabs_in, pabs_out, sum_in, sum_out};
  cs_parall_sum(4, CS_DOUBLE, _s);

  pabs_in  = _s[0] / _s[2];
  pabs_out = _s[1] / _s[3];

  return pabs_out - pabs_in;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the magnitude of a moment of force torque) given an
 *         axis and the stress on a specific boundary.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_list  face list
 * \param[in]   axis         axis
 *
 * \return couple about the axis
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_moment_of_force(cs_lnum_t     n_b_faces,
                        cs_lnum_t    *b_face_list,
                        cs_real_3_t   axis)
{
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  const cs_real_3_t *b_forces
    = (const cs_real_3_t *)cs_field_by_name("boundary_forces")->val;

  cs_real_3_t moment = {0., 0., 0.};

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_real_3_t m;
    cs_lnum_t face_id = b_face_list[i];
    cs_math_3_cross_product(b_face_cog[face_id], b_forces[face_id], m);

    /* b_forces is the stress on the solid boundary,
       thus it comes with a '-' sign here */
    for (int j = 0; j < 3; j++)
      moment[j] -= m[j];
  }
  cs_parall_sum(3, CS_DOUBLE, moment);

  return cs_math_3_dot_product(moment, axis);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
