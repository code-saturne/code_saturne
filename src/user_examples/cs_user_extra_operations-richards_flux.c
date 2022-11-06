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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations-richards_flux.c
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
  const cs_lnum_t n_cells = domain->mesh->n_cells;

  const int nt_cur = domain->time_step->nt_cur;
  const cs_real_t t_cur =  domain->time_step->t_cur;

  /*! [richards_flux_comp] */

  const int iter_step = 10; /* Number of iterations between flux calculations */

  if (nt_cur%iter_step == 0) {

    const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
    const cs_lnum_t n_i_faces = domain->mesh->n_i_faces;
    const cs_lnum_t n_cells_ext  = domain->mesh->n_cells_with_ghosts;

    const cs_lnum_2_t *restrict i_face_cells
      = (const cs_lnum_2_t *restrict)domain->mesh->i_face_cells;
    const cs_real_3_t *restrict surfac
      = (const cs_real_3_t *restrict)domain->mesh_quantities->i_face_normal;

    const cs_field_t *scalar1 = cs_field_by_name("scalar1");
    const cs_real_t *cvar_scal_1 = scalar1->val;

    cs_real_t *viscf, *viscb;
    BFT_MALLOC(viscb, n_b_faces, cs_real_t);
    BFT_MALLOC(viscf, n_i_faces, cs_real_t);

    const cs_real_t eps_geom = 1e-10;

    const int kimasf = cs_field_key_id("inner_mass_flux_id");
    const int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
    const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;

    const int kivisl = cs_field_key_id("diffusivity_id");
    int ifcvsl = cs_field_get_key_int(scalar1, kivisl);
    const cs_real_t *cpro_vscalt = cs_field_by_id(ifcvsl)->val;

    cs_real_t *scalar_diffusion;
    BFT_MALLOC(scalar_diffusion, n_cells_ext, cs_real_t);
    cs_array_copy_real(n_cells, 1, cpro_vscalt, scalar_diffusion);

    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(scalar1);

    cs_face_viscosity(domain->mesh,
                      domain->mesh_quantities,
                      eqp->imvisf,
                      scalar_diffusion,
                      viscf,
                      viscb);

    BFT_FREE(scalar_diffusion);

    /* Example of tracer flux computation through an internal surface.
     * Fluxes are calculated whithout reconstruction, and with a simple definition
     * of the concentration at faces */

    cs_lnum_t nlelt = 0;
    cs_lnum_t *lstelt;
    BFT_MALLOC(lstelt, n_i_faces, cs_lnum_t);

    cs_real_t flux_surf[2] = {0, 0};

    const char *criteria[] = {"PLANE1",        /* flux in */
                              "TOP_SOIL1"};    /* flux out */

    for (int surf_id = 0; surf_id < 2; surf_id++) {

      cs_selector_get_i_face_list(criteria[surf_id], &nlelt, lstelt);

      for (cs_lnum_t ilelt = 0 ; ilelt < nlelt ; ilelt++) {

        const cs_lnum_t face_id = lstelt[ilelt];
        const cs_lnum_t c_id1 = i_face_cells[face_id][0];
        const cs_lnum_t c_id2 = i_face_cells[face_id][1];

        cs_real_t tra_face, tra_diff;
        if (c_id1 < n_cells) { /* Avoids duplicate parallel contribution */
          tra_face = 0.5 * (cvar_scal_1[c_id1] + cvar_scal_1[c_id2]);
          tra_diff = cvar_scal_1[c_id1] - cvar_scal_1[c_id2];
        }
        else {
          tra_face = 0;
          tra_diff = 0;
        }

        cs_real_t flux_tmp = imasfl[face_id]*tra_face + viscf[face_id]*tra_diff;

        /* We impose a norm on the direction of the flux, based on the direction
         * of main coordinates, to be sure that fluxes of all faces are all summed
         * in the same direction */
        if (   (surfac[face_id][0] <= (-eps_geom))
            || (   (fabs(surfac[face_id][0]) <= eps_geom)
                && (surfac[face_id][1] <= (-eps_geom)))
            || (   (fabs(surfac[face_id][0]) <= eps_geom)
                && (fabs(surfac[face_id][1]) <= eps_geom)
                && (surfac[face_id][2] <= (-eps_geom))))
          flux_tmp = -flux_tmp;

        flux_surf[surf_id] += flux_tmp;

      }

      cs_parall_sum(2, CS_REAL_TYPE, flux_surf);

    }

    /* Write fluxes in file */
    if (cs_glob_rank_id < 1) {

      FILE *file = NULL;
      file = fopen("flux_tracer1.dat", "a");

      fprintf(file, "#      Time         Flow in        Flow out\n"
                    "#      (s)          (kg/s)          (kg/s)\n");
      fprintf(file, "%14.5e %14.5e %14.5e\n",
              t_cur, flux_surf[0], flux_surf[1]);

      fclose(file);
    }

    BFT_FREE(lstelt);
    BFT_FREE(viscb);
    BFT_FREE(viscf);
  }
  /*![richards_flux_comp]*/

  /*![richards_time_modif]*/

  /* Example of time modification */

  cs_real_t *dt = CS_F_(dt)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    dt[c_id] = pow(t_cur, 0.95)*5.e-2;

    if (dt[c_id] < 5e-2)
      dt[c_id] = 5e-2;

    if (dt[c_id] < 1e3)
      dt[c_id] = 1e3;
  }

  /*![richards_time_modif]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
