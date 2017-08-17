/*============================================================================
 * Methods for particle localization
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_boundary_conditions.h"
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
 *        CALCUL DE LA DIVERGENCE D'UN VECTEUR
 *    (On ne s'embete pas, on appelle 3 fois le gradient)
 * ------------------------------------------------------- */

static void
diverv (cs_real_t    *diverg,
        cs_real_3_t  *u,
        cs_real_3_t  *coefa,
        cs_real_33_t *coefb)
{

  /* ====================================================================
   * 1. INITIALISATIONS
   * ====================================================================*/
  cs_lnum_t ncelet = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t ncel   = cs_glob_mesh->n_cells;

  /* Allocate work arrays */
  cs_real_33_t *grad;
  BFT_MALLOC(grad, ncelet, cs_real_33_t);

  /* ====================================================================
   * Calcul du gradient de U
   * ====================================================================*/

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_vector("Work array",
                     gradient_type,
                     halo_type,
                     1,      /* inc */
                     100,    /* n_r_sweeps, */
                     2,      /* iwarnp */
                     -1,     /* imligp */
                     1e-8,   /* epsrgp */
                     1.5,    /* climgp */
                     (const cs_real_3_t *)coefa,
                     (const cs_real_33_t *)coefb,
                     u,
                     NULL, /* weighted gradient */
                     NULL, /* cpl */
                     grad);

  /* ====================================================================
   * Calcul de la divergence du vecteur
   * ====================================================================*/

  for (cs_lnum_t iel = 0; iel < ncel; iel++)
    diverg[iel]  = grad[iel][0][0] + grad[iel][1][1] + grad[iel][2][2];

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

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  cs_lnum_t ncelet = m->n_cells_with_ghosts;
  cs_lnum_t ncel   = m->n_cells;
  cs_lnum_t nfac   = m->n_i_faces;
  cs_lnum_t nfabor = m->n_b_faces;

  cs_real_t *viscf, *viscb;
  cs_real_t *smbrs;
  cs_real_t *rovsdt;
  cs_real_t *fmala, *fmalb;
  cs_real_t *phia;
  cs_real_t *dpvar;

  /* Allocate temporary arrays */

  BFT_MALLOC(viscf , nfac  , cs_real_t);
  BFT_MALLOC(viscb , nfabor, cs_real_t);
  BFT_MALLOC(smbrs , ncelet, cs_real_t);
  BFT_MALLOC(rovsdt, ncelet, cs_real_t);
  BFT_MALLOC(fmala , nfac  , cs_real_t);
  BFT_MALLOC(fmalb , nfabor, cs_real_t);
  BFT_MALLOC(phia  , ncelet, cs_real_t);
  BFT_MALLOC(dpvar , ncelet, cs_real_t);

  /* Allocate work arrays */
  cs_real_3_t *w;
  BFT_MALLOC(w, ncelet, cs_real_3_t);

  bft_printf(_("   ** RESOLUTION POUR LA VARIABLE Pressure correction"));

  /* ====================================================================
   * 2. TERMES SOURCES
   * ==================================================================== */

  /* --> Initialisation   */

  for (cs_lnum_t iel = 0; iel < ncel; iel++) {

    smbrs[iel]  = 0.0;
    rovsdt[iel] = 0.0;
    phi[iel]    = 0.0;
    phia[iel]   = 0.0;

  }

  /*     "VITESSE" DE DIFFUSION FACE     */
  cs_face_viscosity(m,
                    fvq,
                    cs_glob_space_disc->imvisf,
                    alphal,
                    viscf,
                    viscb);

  /* CALCUL  de div(Alpha Up) avant correction     */
  for (cs_lnum_t iel = 0; iel < ncel; iel++) {

    for (cs_lnum_t isou = 0; isou < 3; isou++)
      w[isou][iel] = -vitessel[isou + iel * 3] * alphal[iel];

  }

  /* --> Calcul du gradient de W1   */
  /*     ========================   */
  /* Allocate temporary arrays */
  cs_real_3_t *coefaw;
  cs_real_33_t *coefbw;

  BFT_MALLOC(coefaw, nfabor, cs_real_3_t);
  BFT_MALLOC(coefbw, nfabor, cs_real_33_t);

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    cs_lnum_t iel = m->b_face_cells[ifac];

    for (cs_lnum_t isou = 0; isou < 3; isou++)
      coefaw[isou][ifac]  = w[isou][iel];

  }

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    for (cs_lnum_t isou = 0; isou < 3; isou++) {

      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        coefbw[jsou][isou][ifac]  = 0.0;

    }

  }

  diverv(smbrs, w, coefaw, coefbw);

  /* Free memory     */
  BFT_FREE (coefaw);
  BFT_FREE (coefbw);

  /* --> Conditions aux limites sur PHI  */
  /*     ==============================  */
  /* Allocate temporary arrays */
  cs_real_t *coefap, *coefbp;
  cs_real_t *cofafp, *cofbfp;
  BFT_MALLOC(coefap, nfabor, cs_real_t);
  BFT_MALLOC(coefbp, nfabor, cs_real_t);
  BFT_MALLOC(cofafp, nfabor, cs_real_t);
  BFT_MALLOC(cofbfp, nfabor, cs_real_t);

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {

    cs_lnum_t iel  = m->b_face_cells[ifac];
    cs_real_t hint = alphal[iel] / fvq->b_dist[ifac];

    if (   itypfb[ifac] == CS_INLET
        || itypfb[ifac] == CS_SMOOTHWALL
        || itypfb[ifac] == CS_ROUGHWALL
        || itypfb[ifac] == CS_SYMMETRY) {

      /* Neumann Boundary Conditions    */

      cs_boundary_conditions_set_neumann_scalar(&coefap[ifac],
                                                &cofafp[ifac],
                                                &coefbp[ifac],
                                                &cofbfp[ifac],
                                                0.0,
                                                hint);
      coefap[ifac]      = 0.0;
      coefbp[ifac]      = 1.0;

    }
    else if (itypfb[ifac] == CS_OUTLET) {

      /* Dirichlet Boundary Condition   */

      cs_boundary_conditions_set_dirichlet_scalar(&coefap[ifac],
                                                  &cofafp[ifac],
                                                  &coefbp[ifac],
                                                  &cofbfp[ifac],
                                                  phia[iel],
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

  for (cs_lnum_t ifac = 0; ifac < nfac; ifac++) {
    fmala[ifac] = 0.0;
    fmalb[ifac] = 0.0;
  }

  /* In the theta-scheme case, set theta to 1 (order 1) */

  cs_var_cal_opt_t  var_cal_opt = cs_parameters_var_cal_opt_default();

  var_cal_opt.iwarni = 2;  /* quasi-debug at this stage, TODO clean */
  var_cal_opt.iconv  = 0;  /* no convection, pure diffusion here */
  var_cal_opt.istat  = -1;
  var_cal_opt.idifft = -1;
  var_cal_opt.isstpc = 0;
  var_cal_opt.nswrgr = 10000;
  var_cal_opt.nswrsm = 2;
  var_cal_opt.imrgra = cs_glob_space_disc->imrgra;
  var_cal_opt.imligr = 1;

  cs_equation_iterative_solve_scalar(idtva0,
                                     -1,           /* field_id (not a field) */
                                     "PoissonL",   /* name */
                                     1,            /* ndircp */
                                     0,            /* iescap */
                                     0,            /* imucpp */
                                     &var_cal_opt,
                                     phia, phia,
                                     coefap, coefbp,
                                     cofafp, cofbfp,
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
  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(cofafp);
  BFT_FREE(cofbfp);
  BFT_FREE(phia);
  BFT_FREE(w);
  BFT_FREE(dpvar);
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
  cs_lnum_t ncel   = cs_glob_mesh->n_cells;
  cs_lnum_t ncelet = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;

  /* Allocate a temporary array     */

  cs_real_t *phil;
  BFT_MALLOC(phil, ncelet, cs_real_t);

  /* Initialization */

  cs_lagr_particle_set_t *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  /* Means of global class */

  int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

  cs_field_t *mean_vel
    = cs_lagr_stat_get_moment(stat_type,
                              CS_LAGR_MOMENT_MEAN,
                              0,
                              -1);

  cs_field_t *mean_fv
    = cs_lagr_stat_get_moment(CS_LAGR_STAT_VOLUME_FRACTION,
                              CS_LAGR_MOMENT_MEAN,
                              0,
                              -1);

  cs_field_t  *stat_weight = cs_lagr_stat_get_stat_weight(0);

  _lageqp(mean_vel->val, mean_fv->val, phil, itypfb);

  /* Compute gradient of phi corrector */

  cs_real_3_t *grad;

  BFT_MALLOC(grad, ncelet, cs_real_3_t);

  cs_real_t *coefap, *coefbp;

  BFT_MALLOC(coefap, nfabor, cs_real_t);
  BFT_MALLOC(coefbp, nfabor, cs_real_t);

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {
    cs_lnum_t iel = cs_glob_mesh->b_face_cells[ifac];
    coefap[ifac] = phil[iel];
    coefbp[ifac] = 0.0;
  }

  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_scalar("Work array",
                     gradient_type,
                     halo_type,
                     1,               /* inc */
                     1,               /* recompute_cocg */
                     100,             /* n_r_sweeps */
                     0,               /* idimtr */
                     0,               /* hyd_p_flag */
                     1,               /* w_stride */
                     2,               /* iwarnp */
                     -1,              /* imligp */
                     1e-8,            /* epsrgp */
                     0.,              /* extrap */
                     1.5,             /* climgp */
                     NULL,            /* f_ext */
                     coefap,
                     coefbp,
                     phil,
                     NULL,            /* c_weight */
                     NULL, /* internal coupling */
                     grad);

  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(phil);

  /* Correct mean velocities */

  for (cs_lnum_t iel = 0; iel < ncel; iel++) {

    if (stat_weight->val[iel] > cs_glob_lagr_stat_options->threshold) {
      for (cs_lnum_t id = 0; id < 3; id++)
        mean_vel->val[iel * 3 + id] += - grad[iel][id];
    }

  }

  /* Correct instant velocities */

  for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

    unsigned char *part = p_set->p_buffer + p_am->extents * npt;
    cs_lnum_t      iel  = cs_lagr_particle_get_cell_id(part, p_am);

    if (iel >= 0) {

      cs_real_t *part_vel = cs_lagr_particle_attr(part, p_am, CS_LAGR_VELOCITY);

      for (cs_lnum_t id = 0; id < 3; id++)
        part_vel[id] += -grad[id][iel];

    }

  }

  /* Free memory */

  BFT_FREE(grad);

}

END_C_DECLS
