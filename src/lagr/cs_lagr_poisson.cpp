/*============================================================================
 * Methods for particle localization
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_equation_iterative_solve.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_gradient.h"
#include "cs_face_viscosity.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_stat.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_poisson.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* -------------------------------------------------------
 *        COMPUTE THE DIVERGENCE OF A VECTEUR
 *       (Simply call the gradient three times)
 * ------------------------------------------------------- */

static void
diverv (cs_real_t                  *diverg,
        cs_real_3_t                *u,
        const cs_field_bc_coeffs_t *bc_coeffs_v)
{

  /* ====================================================================
   * 1. INITIALISATIONS
   * ====================================================================*/
  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells   = cs_glob_mesh->n_cells;

  /* Allocate work arrays */
  cs_real_33_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_33_t);

  /* ====================================================================
   * Calcul du gradient de U
   * ====================================================================*/

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_vector("Work array",
                     gradient_type,
                     halo_type,
                     1,                      /* inc */
                     100,                    /* n_r_sweeps, */
                     2,                      /* iwarnp */
                     CS_GRADIENT_LIMIT_NONE, /* imligp */
                     1e-8,                   /* epsrgp */
                     1.5,                    /* climgp */
                     bc_coeffs_v,
                     u,
                     NULL,                   /* weighted gradient */
                     NULL,                   /* cpl */
                     grad);

  /* ====================================================================
   * Calcul de la divergence du vecteur
   * ====================================================================*/

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    diverg[c_id]  = grad[c_id][0][0] + grad[c_id][1][1] + grad[c_id][2][2];

  /* Free memory     */
  BFT_FREE(grad);
}

/*-------------------------------------------------------------------
 *          RESOLUTION D'UNE EQUATION DE POISSON
 *            div[ALPHA grad(PHI)] = div(ALPHA <Up>)
 *-------------------------------------------------------------------*/

static void
_lageqp(cs_real_t   *vitessel,
        cs_real_t   *alphal,
        cs_real_t   *phi,
        const int    itypfb[])
{
  /* ====================================================================
   * 1. INITIALISATION
   * ====================================================================*/

  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  cs_lnum_t n_cells = m->n_cells;
  cs_lnum_t n_i_faces = m->n_i_faces;
  cs_lnum_t n_b_faces = m->n_b_faces;

  cs_real_t *viscf, *viscb;
  cs_real_t *smbrs;
  cs_real_t *rovsdt;
  cs_real_t *fmala, *fmalb;
  cs_real_t *phia;
  cs_real_t *dpvar;

  /* Allocate temporary arrays */

  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);
  BFT_MALLOC(smbrs, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);
  BFT_MALLOC(fmala, n_i_faces, cs_real_t);
  BFT_MALLOC(fmalb, n_b_faces, cs_real_t);
  BFT_MALLOC(phia, n_cells_ext, cs_real_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  /* Allocate work arrays */
  cs_real_3_t *w;
  BFT_MALLOC(w, n_cells_ext, cs_real_3_t);

  bft_printf(_("   ** RESOLUTION POUR LA VARIABLE Pressure correction"));

  /* ====================================================================
   * 2. TERMES SOURCES
   * ==================================================================== */

  /* --> Initialization   */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    smbrs[c_id]  = 0.0;
    rovsdt[c_id] = 0.0;
    phi[c_id]    = 0.0;
    phia[c_id]   = 0.0;
  }

  /*     "VITESSE" DE DIFFUSION FACE     */
  cs_face_viscosity(m,
                    fvq,
                    cs_glob_space_disc->imvisf,
                    alphal,
                    viscf,
                    viscb);

  /* CALCUL  de div(Alpha Up) avant correction     */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    for (cs_lnum_t isou = 0; isou < 3; isou++)
      w[c_id][isou] = -vitessel[isou + c_id * 3] * alphal[c_id];

  }

  /* --> Calcul du gradient de W1   */
  /*     ========================   */
  /* Allocate temporary arrays */

  cs_field_bc_coeffs_t bc_coeffs_v_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_v_loc);
  BFT_MALLOC(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t);
  BFT_MALLOC(bc_coeffs_v_loc.b, 9*n_b_faces, cs_real_t);

  cs_real_3_t  *coefaw = (cs_real_3_t  *)bc_coeffs_v_loc.a;
  cs_real_33_t *coefbw = (cs_real_33_t *)bc_coeffs_v_loc.b;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t c_id = m->b_face_cells[f_id];

    for (cs_lnum_t isou = 0; isou < 3; isou++)
      coefaw[f_id][isou]  = w[c_id][isou];

  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    for (cs_lnum_t isou = 0; isou < 3; isou++) {

      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        coefbw[jsou][isou][f_id] = 0.0;

    }

  }

  diverv(smbrs, w, &bc_coeffs_v_loc);

  /* Free memory */
  BFT_FREE(bc_coeffs_v_loc.a);
  BFT_FREE(bc_coeffs_v_loc.b);

  /* --> Boundary condition for PHI  */
  /*     ==============================  */
  /* Allocate temporary arrays */

  cs_field_bc_coeffs_t bc_coeffs_phi_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_phi_loc);

  BFT_MALLOC(bc_coeffs_phi_loc.a,  n_b_faces, cs_real_t);
  BFT_MALLOC(bc_coeffs_phi_loc.b,  n_b_faces, cs_real_t);
  BFT_MALLOC(bc_coeffs_phi_loc.af, n_b_faces, cs_real_t);
  BFT_MALLOC(bc_coeffs_phi_loc.bf, n_b_faces, cs_real_t);

  cs_real_t *coefap = bc_coeffs_phi_loc.a;
  cs_real_t *coefbp = bc_coeffs_phi_loc.b;
  cs_real_t *cofafp = bc_coeffs_phi_loc.af;
  cs_real_t *cofbfp = bc_coeffs_phi_loc.bf;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t c_id  = m->b_face_cells[f_id];
    cs_real_t hint = alphal[c_id] / fvq->b_dist[f_id];

    if (   itypfb[f_id] == CS_INLET
        || itypfb[f_id] == CS_SMOOTHWALL
        || itypfb[f_id] == CS_ROUGHWALL
        || itypfb[f_id] == CS_SYMMETRY) {

      /* Neumann Boundary Conditions    */

      cs_boundary_conditions_set_neumann_scalar(f_id,
                                                &bc_coeffs_phi_loc,
                                                0.0,
                                                hint);
      coefap[f_id] = 0.0;
      coefbp[f_id] = 1.0;

    }
    else if (itypfb[f_id] == CS_OUTLET) {

      /* Dirichlet Boundary Condition   */

      cs_boundary_conditions_set_dirichlet_scalar(f_id,
                                                  &bc_coeffs_phi_loc,
                                                  phia[c_id],
                                                  hint,
                                                  -1);

    }
    else
      bft_error
        (__FILE__, __LINE__, 0,
         _("\n%s (Lagrangian module):\n"
           " unexpected boundary conditions for Phi."), __func__);

  }

  /* ====================================================================
   * 3. RESOLUTION
   * ====================================================================   */

  int idtva0 = 0;  /* No steady option here */

  /* Cancel mass fluxes */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    fmala[f_id] = 0.0;
    fmalb[f_id] = 0.0;
  }

  /* In the theta-scheme case, set theta to 1 (order 1) */


  cs_equation_param_t eqp_loc = cs_parameters_equation_param_default();

  eqp_loc.verbosity = 2;  /* quasi-debug at this stage, TODO clean */
  eqp_loc.iconv  = 0;     /* no convection, pure diffusion here */
  eqp_loc.istat  = -1;
  eqp_loc.ndircl = 1;
  eqp_loc.idifft = -1;
  eqp_loc.isstpc = 0;
  eqp_loc.nswrgr = 10000;
  eqp_loc.nswrsm = 2;
  eqp_loc.imrgra = cs_glob_space_disc->imrgra;
  eqp_loc.imligr = 1;

  cs_equation_iterative_solve_scalar(idtva0,
                                     1,            /* external sub-iteration? */
                                     -1,           /* field_id (not a field) */
                                     "PoissonL",   /* name */
                                     0,            /* iescap */
                                     0,            /* imucpp */
                                     -1,           /* normp */
                                     &eqp_loc,
                                     phia, phia,
                                     &bc_coeffs_phi_loc,
                                     fmala, fmalb,
                                     viscf, viscb,
                                     viscf, viscb,
                                     NULL,         /* viscel */
                                     NULL,         /* weighf */
                                     NULL,         /* weighb */
                                     0,            /* icvflb (all upwind) */
                                     NULL,         /* icvfli */
                                     rovsdt,
                                     smbrs,
                                     phi,
                                     dpvar,
                                     NULL,         /* xcpp */
                                     NULL);        /* eswork */

  /* Free memory */

  BFT_FREE(viscf);
  BFT_FREE(viscb);
  BFT_FREE(smbrs);
  BFT_FREE(rovsdt);
  BFT_FREE(fmala);
  BFT_FREE(fmalb);
  BFT_FREE(phia);
  BFT_FREE(w);
  BFT_FREE(dpvar);

  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(cofafp);
  BFT_FREE(cofbfp);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------*/
/*! \brief Solve Poisson equation for mean particle velocities
 * and correct particle instantaneous velocities
 *
 * \param[in]  itypfb  boundary face type
 */
/*-----------------------------------------------------------------*/

void
cs_lagr_poisson(const int  itypfb[])
{
  cs_lnum_t n_cells   = cs_glob_mesh->n_cells;
  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* Allocate a temporary array     */

  cs_real_t *phil;
  BFT_MALLOC(phil, n_cells_ext, cs_real_t);

  /* Initialization */

  cs_lagr_particle_set_t *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  /* Means of global class */

  int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

  cs_field_t *mean_vel
    = cs_lagr_stat_get_moment(stat_type,
                              CS_LAGR_STAT_GROUP_PARTICLE,
                              CS_LAGR_MOMENT_MEAN,
                              0,
                              -1);

  cs_field_t *mean_fv
    = cs_lagr_stat_get_moment(CS_LAGR_STAT_VOLUME_FRACTION,
                              CS_LAGR_STAT_GROUP_PARTICLE,
                              CS_LAGR_MOMENT_MEAN,
                              0,
                              -1);

  cs_field_t  *stat_weight = cs_lagr_stat_get_stat_weight(0);

  _lageqp(mean_vel->val, mean_fv->val, phil, itypfb);

  /* Compute gradient of phi corrector */

  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_field_bc_coeffs_t bc_coeffs_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_loc);

  BFT_MALLOC(bc_coeffs_loc.a,  n_b_faces, cs_real_t);
  BFT_MALLOC(bc_coeffs_loc.b,  n_b_faces, cs_real_t);
  cs_real_t *coefap = bc_coeffs_loc.a;
  cs_real_t *coefbp = bc_coeffs_loc.b;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = cs_glob_mesh->b_face_cells[f_id];
    coefap[f_id] = phil[c_id];
    coefbp[f_id] = 0.0;
  }

  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_scalar("Work array",
                     gradient_type,
                     halo_type,
                     1,                      /* inc */
                     100,                    /* n_r_sweeps */
                     0,                      /* hyd_p_flag */
                     1,                      /* w_stride */
                     2,                      /* iwarnp */
                     CS_GRADIENT_LIMIT_NONE, /* imligp */
                     1e-8,                   /* epsrgp */
                     1.5,                    /* climgp */
                     NULL,                   /* f_ext */
                     &bc_coeffs_loc,
                     phil,
                     NULL,            /* c_weight */
                     NULL, /* internal coupling */
                     grad);

  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(phil);

  /* Correct mean velocities */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (stat_weight->val[c_id] > cs_glob_lagr_stat_options->threshold) {
      for (cs_lnum_t id = 0; id < 3; id++)
        mean_vel->val[c_id * 3 + id] += - grad[c_id][id];
    }

  }

  /* Correct instant velocities */

  for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

    unsigned char *part = p_set->p_buffer + p_am->extents * npt;
    cs_lnum_t      c_id  = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_CELL_ID);

    if (c_id >= 0) {

      cs_real_t *part_vel =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(part, p_am, CS_LAGR_VELOCITY);

      for (cs_lnum_t id = 0; id < 3; id++)
        part_vel[id] += -grad[c_id][id];

    }

  }

  /* Free memory */

  BFT_FREE(grad);

}

END_C_DECLS
