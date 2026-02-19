/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User operations called at the end of each time step.
 *
 * This function has a very general purpose, although it is recommended to
 * handle mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations([[maybe_unused]] cs_domain_t  *domain)
{
  FILE *file = nullptr;

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)domain->mesh_quantities->b_face_cog;
  const cs_nreal_3_t *b_face_u_normal
    = domain->mesh_quantities->b_face_u_normal;

  const cs_real_t    *f_b_temp
    = (const cs_real_t *)cs_field("boundary_temperature")->val;
  const cs_real_3_t  *b_stress
    = (const cs_real_3_t *)cs_field("boundary_stress")->val;

  cs_field_t *f = cs_thermal_model_field();
  const double visls_0 = f->get_key_double("diffusivity_ref");

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max) {

    /* Compute the thermal fluxes at selected boundary faces */
    const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
    const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];
    cs_array<cs_real_t> boundary_flux(n_elts);

    cs_post_boundary_flux(f->name, n_elts, nullptr, boundary_flux.data());

    /* Some declarations -> Declared like this, the arrays are of size 0*/
    cs_array<cs_real_t> glo_nusselt;
    cs_array<cs_real_t> glo_friction;
    cs_array<cs_real_t> glo_coords;

    /* Print Nusselt and friction coeff file header*/
    if (cs_glob_rank_id <= 0) {
       file = fopen("Surface_values.dat","w");
       fprintf(file, "# This routine writes values at walls\n");
       fprintf(file, "# 1:Coords, 2:Cf, 3:Nu \n");
    }

    /* Select boundary faces to print */

    /* -------- TO MODIFY ---------- */
    const char criteria[] = "Wall";
    /* -------- TO MODIFY ---------- */

    cs_lnum_t   n_selected_faces   = 0;
    cs_gnum_t   n_selected_faces_g = 0;
    cs_lnum_t  *selected_faces = nullptr;

    CS_MALLOC(selected_faces, n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(criteria,
                                &n_selected_faces,
                                selected_faces);

    /* Get the total number of faces shared on all ranks */
    n_selected_faces_g = n_selected_faces;
    cs_parall_sum(1 , CS_GNUM_TYPE, &n_selected_faces_g);

    /*Allocate local and global arrays to store the desired data */
    cs_array<cs_real_t> loc_nusselt(n_selected_faces);
    cs_array<cs_real_t> loc_friction(n_selected_faces);
    cs_array<cs_real_t> loc_coords(n_selected_faces);

    if (cs_glob_n_ranks > 1) {
      glo_nusselt.reshape(n_selected_faces_g);
      glo_friction.reshape(n_selected_faces_g);
      glo_coords.reshape(n_selected_faces_g);
    }

    /* Add some reference values */

    /* -------- TO MODIFY ---------- */
    cs_real_t length_ref = 1.0;
    cs_real_t temp_ref = 1.0;
    /* -------- TO MODIFY ---------- */

    for (cs_lnum_t ielt = 0; ielt < n_selected_faces; ielt++) {
      cs_lnum_t f_id = selected_faces[ielt];

      /* Compute the Nusselt number */
      cs_real_t tfac = f_b_temp[f_id];
      loc_nusselt[ielt] =   boundary_flux[f_id] * length_ref
                          / (visls_0 * (tfac-temp_ref));

      /* Compute the friction coefficient */
      const cs_nreal_t *u_n = b_face_u_normal[f_id];
      cs_real_t s_nor = cs_math_3_dot_product(b_stress[f_id], u_n);

      cs_real_t stresses[3];
      for (cs_lnum_t j = 0; j < 3; j++)
        stresses[j] = b_stress[f_id][j] - s_nor*u_n[j];

      loc_friction[ielt] = cs_math_3_norm(stresses);

      /* Here we plot the results with respect to X */
      loc_coords[ielt]   = b_face_cog[f_id][0];
    }

    /* Gather the data of all ranks */
    if (cs_glob_n_ranks > 1) {
      cs_parall_allgather_r(n_selected_faces, n_selected_faces_g,
                            loc_nusselt.data(),      glo_nusselt.data());
      cs_parall_allgather_r(n_selected_faces, n_selected_faces_g,
                            loc_friction.data(),     glo_friction.data());
      cs_parall_allgather_r(n_selected_faces, n_selected_faces_g,
                            loc_coords.data(),       glo_coords.data());
    }

    /* Print in file */
    for (cs_gnum_t ielt = 0; ielt < n_selected_faces_g; ielt++) {
      if (cs_glob_rank_id == -1)
        fprintf(file,"%17.9e %17.9e %17.9e\n",
                loc_coords[ielt], loc_friction[ielt], loc_nusselt[ielt]);
      if (cs_glob_rank_id == 0)
        fprintf(file,"%17.9e %17.9e %17.9e\n",
                glo_coords[ielt], glo_friction[ielt], glo_nusselt[ielt]);
    }

    /* Free allocated memory */
    CS_FREE(selected_faces);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
