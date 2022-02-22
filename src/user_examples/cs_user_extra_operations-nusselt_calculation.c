/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations-nusselt_calculation.c
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
  /* Calculation of the Nusselt number */

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max) {

    /*! [loc_var_f_user] */
    const cs_lnum_t n_cells     = domain->mesh->n_cells;
    const cs_lnum_t n_b_faces    = domain->mesh->n_b_faces;
    const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

    const cs_real_3_t *cell_cen
      = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;
    const cs_real_3_t *cdgfbo
      = (const cs_real_3_t *)domain->mesh_quantities->b_face_cog;
    const cs_real_t *volume  = domain->mesh_quantities->cell_vol;

    FILE *file = NULL;

    cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
    const cs_fluid_properties_t *phys_pro = cs_glob_fluid_properties;

    /* Parameters, case dependant */
    const cs_real_t prandtl = 0.71, height = 1., qwall = 1.;

    const cs_field_t *f = CS_F_(t);
    const cs_field_t *mu = CS_F_(mu);

    cs_lnum_t nlelt;
    cs_lnum_t *lstelt;
    BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list
      ("normal[0, -1, 0, 0.1] and box[-1000, -1000, -1000, 1000, 0.01, 1000]",
       &nlelt, lstelt);
    /*! [loc_var_f_user] */

    /* Compute value reconstructed at I' for boundary faces */

    /*! [compute_nusselt] */
    cs_real_t *treloc;
    BFT_MALLOC(treloc, nlelt, cs_real_t);

    int iortho = 0;
    /*! [compute_nusselt] */

    /* General case (for non-orthogonal meshes) */

    /*! [gen_nusselt] */
    if (iortho == 0) {
      cs_post_field_cell_to_b_face_values(f, nlelt, lstelt, treloc);
    }
    /*! [gen_nusselt] */

    /* Case of orthogonal meshes (no reconstruction needed) */

    else {
      /*! [else_nusselt] */
      const cs_real_t *coefap = CS_F_(t)->bc_coeffs->a;
      const cs_real_t *coefbp = CS_F_(t)->bc_coeffs->b;

      for (cs_lnum_t ielt = 0; ielt < nlelt; ielt++) {
        cs_lnum_t face_id = lstelt[ielt];
        cs_lnum_t c_id = b_face_cells[face_id];
        treloc[ielt] = coefap[face_id] + coefbp[face_id] * f->val[c_id];
      }
      /*! [else_nusselt] */
    }

    /*! [value_ortho_nusselt] */
    if  (cs_glob_rank_id < 1) {
      file = fopen("Nusselt.dat", "w");
    }

    cs_gnum_t neltg = nlelt;
    cs_parall_counter(&neltg, 1);

    cs_real_t *xabs = NULL, *xabsg = NULL, *xnusselt = NULL;
    cs_real_t *treglo = NULL;

    BFT_MALLOC(xabs, nlelt, cs_real_t);

    BFT_MALLOC(xnusselt, neltg, cs_real_t);
    BFT_MALLOC(xabsg, neltg, cs_real_t);
    BFT_MALLOC(treglo, neltg, cs_real_t);

    for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt ++) {
      cs_lnum_t face_id = lstelt[ilelt];
      xabs[ilelt] = cdgfbo[face_id][0];
    }

    cs_parall_allgather_ordered_r(nlelt, neltg, 1, xabs, xabs, xabsg);
    cs_parall_allgather_ordered_r(nlelt, neltg, 1, xabs, treloc, treglo);

    BFT_FREE(xabs);
    BFT_FREE(treloc);
    BFT_FREE(lstelt);
    /*! [value_ortho_nusselt] */

    /* Calculation of the bulk temperature and compute Nusselt number */
    /*! [bulk_nusselt] */
    for (cs_gnum_t ilelt = 0; ilelt < neltg; ilelt ++) {

      cs_real_t xtbulk = 0;
      cs_real_t xubulk = 0;
      cs_real_t lambda = 0;

      int npoint = 200;
      cs_lnum_t c_id_prev = -1;
      int irang1 = -1;
      int irangv;

      for (int ii = 0; ii < npoint; ii++) {

        cs_lnum_t c_id;

        cs_real_t xyz[3] = {xabsg[ilelt],
                            (cs_real_t)ii/(cs_real_t)(npoint-1),
                            0.};

        cs_geom_closest_point(n_cells,
                              cell_cen,
                              xyz,
                              &c_id, &irangv);

        if (c_id != c_id_prev || irangv != irang1) {
          c_id_prev = c_id;
          irang1 = irangv;
          if (cs_glob_rank_id == irangv) {
            cs_real_t xtb = volume[c_id]*f->val[c_id]*vel[c_id][0];
            cs_real_t xub = volume[c_id]*vel[c_id][0];
            xtbulk += xtb;
            xubulk += xub;
            lambda = phys_pro->cp0*mu->val[c_id] / prandtl;
          }
        }

      }

      /* Note: only last location used for lambda.
         Would be more appropriate (but still not general at all)
         with visls0 in definition of lambda and from the derivation
         of Nusselt number; use of lambda at the wall (ii = 0)
         would seem better */
      cs_parall_bcast(irangv, 1, CS_REAL_TYPE, &lambda);

      cs_real_t xbulk[2] = {xtbulk, xubulk};
      cs_parall_sum(2, CS_REAL_TYPE, xbulk);

      xtbulk = xbulk[0] / xbulk[1];
      xubulk = xbulk[1];

      cs_real_t tfac = treglo[ilelt];

      xnusselt[ilelt] = qwall * 2. * height / lambda / (tfac - xtbulk);
    }

    for (cs_gnum_t ii = 0; ii < neltg; ii++)
      fprintf(file,
              "%17.9e %17.9e\n",
              xabsg[ii]*10.,
              xnusselt[ii]/(0.023*pow(30000., 0.8)*pow(0.71, 0.4)));

    if (file != NULL)
      fclose(file);

    BFT_FREE(xnusselt);
    BFT_FREE(xabsg);
    BFT_FREE(treglo);
    /*! [bulk_nusselt] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
