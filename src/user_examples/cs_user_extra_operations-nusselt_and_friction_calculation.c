/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user subroutine)
 *============================================================================*/

/* VERS */

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  FILE *file = NULL;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_real_3_t *surfbo = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
  const cs_real_t *surfbn = cs_glob_mesh_quantities->b_face_surf;

  const cs_real_t *f_b_temp = (const cs_real_t *)cs_field_by_name("boundary_temperature")->val;
  const cs_real_3_t  *forbr = (const cs_real_3_t *)cs_field_by_name("boundary_forces")->val;

  cs_field_t *f = cs_thermal_model_field();
  const double visls0 = cs_field_get_key_double(f, cs_field_key_id("diffusivity_ref"));

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max) {

    /* We compute the thermal fluxes at boundaries */
    const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
    const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];
    cs_real_t *boundary_flux ;
    BFT_MALLOC(boundary_flux, n_elts, cs_real_t);

    cs_post_boundary_flux(f->name, n_elts, NULL, boundary_flux);

    /* Some declarations */
    cs_real_t srfbn, srfnor[3], fornor, stresses[3] ;
    cs_real_t *loc_nusselt , *glo_nusselt  ;
    cs_real_t *loc_friction, *glo_friction ;
    cs_real_t *loc_coords  , *glo_coords   ;

    /* Print Nusselt and friction coeff file header*/
    if (cs_glob_rank_id <= 0) {
       file = fopen("Surface_values.dat","w");
       fprintf(file, "# This routine writes values at walls\n");
       fprintf(file, "# 1:Coords, 2:Cf, 3:Nu \n");
    }

    /* Select boudary faces to print */

    /* -------- TO MODIFY ---------- */
    const char criteria[] = "Wall";
    /* -------- TO MODIFY ---------- */

    cs_lnum_t   n_selected_faces   = 0;
    cs_lnum_t   n_selected_faces_g = 0;
    cs_lnum_t  *selected_faces = NULL;

    BFT_MALLOC(selected_faces, n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(criteria,
                                &n_selected_faces,
                                selected_faces);

    /* Get the total number of faces shared on each procs */
    n_selected_faces_g = n_selected_faces ;
    if (cs_glob_rank_id >= 0) {
      cs_parall_sum(1 , CS_LNUM_TYPE, &n_selected_faces_g);
    }

    /*Allocate local and global arrays to store the desired data */
    BFT_MALLOC(loc_nusselt , n_selected_faces, cs_real_t);
    BFT_MALLOC(loc_friction, n_selected_faces, cs_real_t);
    BFT_MALLOC(loc_coords  , n_selected_faces, cs_real_t);

    if (cs_glob_rank_id >= 0) {
      BFT_MALLOC(glo_nusselt , n_selected_faces_g, cs_real_t);
      BFT_MALLOC(glo_friction, n_selected_faces_g, cs_real_t);
      BFT_MALLOC(glo_coords  , n_selected_faces_g, cs_real_t);
    }

    /* Add some reference values */

    /* -------- TO MODIFY ---------- */
    cs_real_t longueur_ref = 1.0 ;
    cs_real_t temp_ref     = 1.0 ;
    cs_real_t vel_ref      = 1.0 ;
    /* -------- TO MODIFY ---------- */

    cs_real_t tfac ;

    for (cs_lnum_t ielt = 0; ielt < n_selected_faces; ielt++) {
      cs_lnum_t f_id = selected_faces[ielt];

      // Compute the Nusselt number
      tfac = f_b_temp[f_id];
      loc_nusselt[ielt]  = boundary_flux[f_id] * longueur_ref
                         / ( visls0 * (tfac-temp_ref) );

      // Compute the friction coefficient
      srfbn = surfbn[f_id];
      srfnor[0] = surfbo[f_id][0] / srfbn;
      srfnor[1] = surfbo[f_id][1] / srfbn;
      srfnor[2] = surfbo[f_id][2] / srfbn;
      fornor = forbr[f_id][0]*srfnor[0]
             + forbr[f_id][1]*srfnor[1]
             + forbr[f_id][2]*srfnor[2];

      stresses[0] = (forbr[f_id][0] - fornor*srfnor[0]) / srfbn ;
      stresses[1] = (forbr[f_id][1] - fornor*srfnor[1]) / srfbn ;
      stresses[2] = (forbr[f_id][2] - fornor*srfnor[2]) / srfbn ;

      loc_friction[ielt] = sqrt( stresses[0]*stresses[0] +
                                 stresses[1]*stresses[1] +
                                 stresses[2]*stresses[2] );

      // Here we plot the results with respects to X
      loc_coords[ielt]   = b_face_cog[f_id][0];
    }

    /* Gather the data of each procs */
    if (cs_glob_rank_id >= 0) {
      cs_parall_allgather_r(n_selected_faces , n_selected_faces_g,
                            loc_nusselt      , glo_nusselt       );
      cs_parall_allgather_r(n_selected_faces , n_selected_faces_g,
                            loc_friction     , glo_friction      );
      cs_parall_allgather_r(n_selected_faces , n_selected_faces_g,
                            loc_coords       , glo_coords        );
    }

    /* Print in file */
    for (cs_lnum_t ielt = 0; ielt < n_selected_faces_g; ielt++) {
      if (cs_glob_rank_id == -1)
        fprintf(file,"%17.9e %17.9e %17.9e\n",
                loc_coords[ielt], loc_friction[ielt], loc_nusselt[ielt]);
      if (cs_glob_rank_id == 0)
        fprintf(file,"%17.9e %17.9e %17.9e\n",
                glo_coords[ielt], glo_friction[ielt], glo_nusselt[ielt]);
    }

    /* Free allocated memmory */
    BFT_FREE(boundary_flux);
    BFT_FREE(selected_faces);

    BFT_FREE(loc_nusselt);
    BFT_FREE(loc_friction);
    BFT_FREE(loc_coords);

    if (cs_glob_rank_id >= 0) {
      BFT_FREE(glo_nusselt);
      BFT_FREE(glo_friction);
      BFT_FREE(glo_coords);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
