/*============================================================================
 * Rij-epsilon turbulence model.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_base.h"
#include "cs_equation_iterative_solve.h"
#include "cs_equation_param.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_log_iteration.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_time_step.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_rij.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of alpha in the framwork of the Rij-EBRSM model.
 *
 *  \param[in]  f_id          field id of alpha variable
 *  \param[in]  n_cells       number of cells
 *  \param[in]  alpha_min     minimum acceptable value for alpha
 */
/*----------------------------------------------------------------------------*/

static void
_clip_alpha(const int          f_id,
            const cs_lnum_t    n_cells,
            const cs_real_t    alpha_min[])
{
  cs_real_t *cvar_al = cs_field_by_id(f_id)->val;

  int kclipp = cs_field_key_id("clipping_id");
  cs_gnum_t nclp[2] =  {0, 0};  /* Min and max clipping values respectively */

  /* Postprocess clippings ? */
  cs_real_t *cpro_a_clipped = NULL;
  int clip_a_id = cs_field_get_key_int(cs_field_by_id(f_id), kclipp);
  if (clip_a_id > -1) {
    cpro_a_clipped = cs_field_by_id(clip_a_id)->val;
    cs_array_set_value_real(n_cells, 1, 0, cpro_a_clipped);
  }

  /* Store Min and Max for logging */
  cs_real_t vmin[1] = {cs_math_big_r};
  cs_real_t vmax[1] = {-cs_math_big_r};

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t var = cvar_al[c_id];
    vmin[0] = cs_math_fmin(vmin[0], var);
    vmax[0] = cs_math_fmax(vmax[0], var);
  }

  cs_parall_min(1, CS_REAL_TYPE, vmin);
  cs_parall_max(1, CS_REAL_TYPE, vmax);

  /* Clipping (edit to avoid exactly zero values) */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (cvar_al[c_id] < alpha_min[c_id]) {
      if (clip_a_id > -1)
        cpro_a_clipped[c_id] = alpha_min[c_id]-cvar_al[c_id];
      nclp[0] += 1;
      cvar_al[c_id] = alpha_min[c_id];
    }
    else if (cvar_al[c_id] >= 1) {
      if (clip_a_id > -1)
        cpro_a_clipped[c_id] = cvar_al[c_id]-1.0;
      nclp[1] += 1;
      cvar_al[c_id] = 1.0;
    }
  }

  cs_parall_counter(nclp, 2);

  cs_lnum_t iclpmn[1] = {nclp[0]}, iclpmx[1] = {nclp[1]};
  cs_log_iteration_clipping_field(f_id,
                                  iclpmn[0],
                                  iclpmx[0],
                                  vmin, vmax,
                                  iclpmn, iclpmx);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Solve the equation on alpha in the framework of the Rij-EBRSM model.
 *
 * Also called for alpha of scalars for EB-DFM.
 *
 * \param[in]  f_id          field id of alpha variable
 * \param[in]  c_durbin_l    constant for the Durbin length
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_alpha(int        f_id,
                        cs_real_t  c_durbin_l)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *distb = fvq->b_dist;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;

  cs_real_t *cvar_al = cs_field_by_id(f_id)->val;

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;
  const cs_real_t *cvara_al = cs_field_by_id(f_id)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab =  cs_field_get_key_int(CS_F_(vel), kbmasf);

  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  const cs_real_t d1s2 = 0.50;
  const cs_real_t d1s4 = 0.25;
  const cs_real_t d3s2 = 1.50;
  const cs_real_t uref = cs_glob_turb_ref_values->uref;

  /* Resolving the equation of alpha
     =============================== */

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param(cs_field_by_id(f_id));

  if (eqp->iwarni == 1) {
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" Solving the variable %s\n"),
                  cs_field_get_label(cs_field_by_id(f_id)));
  }

  cs_real_t *rovsdt, *smbr;

  /* Allocate temporary arrays */
  BFT_MALLOC(smbr, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  cs_array_set_value_real(n_cells, 1, 0, smbr);
  cs_array_set_value_real(n_cells, 1, 0, rovsdt);

  /* Source term of alpha
   *  \f$ smbr = \dfrac{1}{L^2 (\alpha)} - \dfrac{1}{L^2}\f$
   * In fact there is a mark "-" because the solved equation is
   *   \f$-\div{\grad {alpha}} = smbr \f$
   *================================================================*/

  /* Matrix */

  const cs_real_t thetv = eqp->thetav;
  cs_real_t thetap = (cs_glob_time_scheme->isto2t > 0) ? thetv : 1.0;

  // FIXME the source term extrapolation is not well done!!!!

  /* For automatic initialization, the length scale is fixed at L^+ =50 */

  if (   cs_glob_time_step->nt_cur == 1
      && cs_glob_turb_rans_model->reinit_turb == 1) {

    const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
    const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;
    const cs_real_t xlldrb = 50.0 * viscl0 / ro0 / (0.050*uref);
    const cs_real_t l2 = cs_math_pow2(xlldrb);

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Explicit term */
      smbr[c_id] = cell_f_vol[c_id]*(1.0-cvara_al[c_id]) / l2;

      /* Implicit term */
      rovsdt[c_id] = (rovsdt[c_id]+cell_f_vol[c_id]*thetap) / l2;
    }

  }
  else {

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xk
        = d1s2 * (cvara_rij[c_id][0] + cvara_rij[c_id][1] + cvara_rij[c_id][2]);
      const cs_real_t xnu = viscl[c_id] / crom[c_id];

      /* Integral length scale */
      const cs_real_t xllke = pow(xk, d3s2) / cvara_ep[c_id];

      /* Kolmogorov length scale */
      const cs_real_t xllkmg =   cs_turb_xceta
                               * pow(cs_math_pow3(xnu)/cvara_ep[c_id], d1s4);

      /* Durbin length scale */
      const cs_real_t xlldrb = c_durbin_l*cs_math_fmax(xllke, xllkmg);

      const cs_real_t l2 = cs_math_pow2(xlldrb);

      /* Explicit term */
      smbr[c_id] = cell_f_vol[c_id]*(1.0-cvara_al[c_id]) / l2;

      /* Implicit term */
      rovsdt[c_id] = (rovsdt[c_id]+cell_f_vol[c_id]*thetap) / l2;
    }

  }

  /* Calculation of viscf and viscb for cs_equation_iterative_solve_scalar. */
  cs_real_t *w1, *viscf, *viscb;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  cs_array_set_value_real(n_cells, 1, 1, w1);

  cs_face_viscosity(m,
                    fvq,
                    eqp->imvisf,
                    w1,
                    viscf,
                    viscb);

  BFT_FREE(w1);

  /* Effective resolution of the equation of alpha
     ============================================= */

  const cs_real_t *coefap = cs_field_by_id(f_id)->bc_coeffs->a;
  const cs_real_t *coefbp = cs_field_by_id(f_id)->bc_coeffs->b;
  const cs_real_t *cofafp = cs_field_by_id(f_id)->bc_coeffs->af;
  const cs_real_t *cofbfp = cs_field_by_id(f_id)->bc_coeffs->bf;

  cs_equation_param_t eqp_loc = *eqp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.thetav = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_real_t *dpvar;
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1, /* init */
                                     f_id,
                                     NULL,
                                     0, /* iescap */
                                     0, /* imucpp */
                                     -1, /* normp */
                                     &eqp_loc,
                                     cvara_al,
                                     cvara_al,
                                     coefap,
                                     coefbp,
                                     cofafp,
                                     cofbfp,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     NULL,
                                     NULL,
                                     NULL,
                                     0, /* boundary convective upwind flux */
                                     NULL,
                                     rovsdt,
                                     smbr,
                                     cvar_al,
                                     dpvar,
                                     NULL,
                                     NULL);

  BFT_FREE(dpvar);
  BFT_FREE(smbr);

  /* Clipping
     ======== */

  cs_real_t *alpha_min;
  BFT_MALLOC(alpha_min, n_cells_ext, cs_real_t);

  /* Compute a first estimator of the minimal value of alpha per cell.
   * This is deduced from "alpha/L^2 - div(grad alpha) = 1/L^2" and assuming that
   * boundary cell values are 0. This value is thefore non zero but
   * much smaller than the wanted value. */

  cs_array_set_value_real(n_cells_ext, 1, 0, alpha_min);
  cs_array_copy_real(n_cells, 1, rovsdt, alpha_min);

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    const cs_lnum_t ii = i_face_cells[face_id][0];
    const cs_lnum_t jj = i_face_cells[face_id][1];
    alpha_min[ii] += viscf[face_id];
    alpha_min[jj] += viscf[face_id];
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    const cs_lnum_t ii = b_face_cells[face_id];
    alpha_min[ii] += viscb[face_id]/distb[face_id];
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    alpha_min[c_id] = rovsdt[c_id]/alpha_min[c_id];

  _clip_alpha(f_id, n_cells, alpha_min);

  BFT_FREE(alpha_min);
  BFT_FREE(rovsdt);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
}

/*----------------------------------------------------------------------------*/

 END_C_DECLS
