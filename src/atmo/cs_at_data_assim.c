/*============================================================================
 * Data assimilation
 *============================================================================*/

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_at_opt_interp.h"
#include "cs_ext_neighborhood.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_measures_util.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_math.h"
#include "cs_time_step.h"
#include "cs_restart_default.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_at_data_assim.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_at_data_assim.c
        Data assimilation functions.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/* debugging switch */
#define _DA_DEBUG_ 1

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool _data_assim_activated = false;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize data assimilation structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_initialize(void)
{
  const int key_ms = cs_field_key_id("measures_set_id");
  const int key_oi = cs_field_key_id("opt_interp_id");
  const int key_oia = cs_field_key_id("opt_interp_analysis_id");
  const int key_vis = cs_field_key_id("post_vis");
  const int key_log = cs_field_key_id("log");
  const int key_rst = cs_field_key_id("restart_file");

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    /* optimal interpolation can be defined only on variables
       and can not be defined on pressure */

    if (!(f->type & CS_FIELD_VARIABLE) || f->id == CS_F_(p)->id)
      continue;

    /* Is an optimal interpolation defined for current variable? */

    int oi_id = cs_field_get_key_int(f, key_oi);
    if (oi_id == -1)
      continue;

    _data_assim_activated = true;

    /* create measures set for current variable */

    size_t fn_size = strlen(f->name);
    char *name_buf;
    size_t suf_size = 3;
    BFT_MALLOC(name_buf, fn_size + suf_size + 1, char);
    snprintf(name_buf, fn_size + suf_size + 1, "%s_ms", f->name);

    int type_flag = 0;
    bool ilved = true;
    cs_measures_set_t *ms = cs_measures_set_create(name_buf,
                                                   type_flag,
                                                   f->dim,
                                                   ilved);
    cs_field_set_key_int(f, key_ms, ms->id);

    /* create interpol grid for current variable */

    snprintf(name_buf, fn_size + suf_size + 1, "%s_ig", f->name);

    cs_interpol_grid_t *ig = cs_interpol_grid_create(name_buf);

    /* create optimal interpolation for current variable */

    snprintf(name_buf, fn_size + suf_size + 1, "%s_oi", f->name);

    cs_at_opt_interp_t *oi = cs_at_opt_interp_create(name_buf);
    BFT_FREE(name_buf);
    cs_field_set_key_int(f, key_oi, oi->id);

    oi->ig_id = ig->id;

    char filename[50];
    sprintf(filename, "%s_%s", "measures", f->name);

    cs_at_opt_interp_read_file(filename, ms, oi, f->dim);

    cs_at_opt_interp_map_values(oi, ms);

    cs_at_data_assim_log(ms, oi, f);

    /* create analysis field for current variable */

    suf_size = 9;
    BFT_MALLOC(name_buf, fn_size + suf_size + 1, char);
    snprintf(name_buf, fn_size + suf_size + 1, "%s_analysis", f->name);

    cs_field_t *oia_f = cs_field_create(name_buf,
                                        CS_FIELD_PROPERTY,
                                        CS_MESH_LOCATION_CELLS,
                                        f->dim,
                                        false);
    BFT_FREE(name_buf);
    cs_field_set_key_int(f, key_oia, oia_f->id);

    cs_field_set_key_int(oia_f, key_vis, CS_POST_ON_LOCATION);
    cs_field_set_key_int(oia_f, key_log, 1);

    /* back up analysis in auxiliary restart file */
    cs_field_set_key_int(oia_f, key_rst, CS_RESTART_AUXILIARY);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build operators.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_build_ops(void)
{
  const int key_ms = cs_field_key_id("measures_set_id");
  const int key_oi = cs_field_key_id("opt_interp_id");

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    /* optimal interpolation can be defined only on variables
       and can not be defined on pressure */

    if (!(f->type & CS_FIELD_VARIABLE) || f->id == CS_F_(p)->id)
      continue;

    /* Is an optimal interpolation defined for current variable? */

    int oi_id = cs_field_get_key_int(f, key_oi);
    if (oi_id == -1)
      continue;

    cs_at_opt_interp_t *oi = cs_at_opt_interp_by_id(oi_id);

    int ms_id = cs_field_get_key_int(f, key_ms);
    cs_measures_set_t *ms = cs_measures_set_by_id(ms_id);
    cs_lnum_t n_obs = ms->nb_measures;

    int ig_id = oi->ig_id;
    cs_interpol_grid_t *ig = cs_interpol_grid_by_id(ig_id);
    cs_interpol_grid_init(ig, n_obs, ms->coords);

#if _DA_DEBUG_
    bft_printf("\n *Start processing variable %s\n\n", f->name);
#endif

    /* Computing observation operator (H) */

    cs_at_opt_interp_obs_operator(ms, oi, ig);

#if _DA_DEBUG_
    if (cs_glob_rank_id <= 0) {
      cs_real_t *proj = oi->model_to_obs_proj;
      cs_lnum_t *proj_idx = oi->model_to_obs_proj_idx;

      for (int ii = 0; ii < n_obs; ii++) {
        bft_printf("    Obs %i\n",ii);
        for (cs_lnum_t jj = proj_idx[ii]; jj < proj_idx[ii+1]; jj++)
          bft_printf("    Point %i x %.2f y %.2f z %.2f coef %.2f\n",
                     jj, (proj + jj*4)[1], (proj + jj*4)[2],
                     (proj + jj*4)[3], (proj + jj*4)[0]);
        bft_printf("\n");
      }
      bft_printf("    Sum of interpolation coefficients\n");
      for (int ii = 0; ii < n_obs; ii++) {
        bft_printf("    ");
        cs_real_t sum = 0.;
        for (cs_lnum_t jj = proj_idx[ii]; jj < proj_idx[ii+1]; jj++)
          sum += (proj + jj*4)[0];
        bft_printf("Obs %i Sum %.5f\n", ii, sum);
      }
      bft_printf("\n");
    }
#endif

    /* pass model covariance matrix in observations space */

    cs_at_opt_interp_project_model_covariance(ms, oi);

#if _DA_DEBUG_
    if (cs_glob_rank_id <= 0) {
      bft_printf("   *Building HBHT\n");
      for (int ii = 0; ii < n_obs; ii++) {
        bft_printf("    ");
        for (int jj = 0; jj < n_obs; jj++)
          bft_printf("%.8f ", oi->b_proj[ii*n_obs + jj]);
        bft_printf("\n");
      }
      bft_printf("\n");

      bft_printf("   *Building R\n");
      for (int kk = 0; kk < ms->dim; kk++) {
        bft_printf("   Comp. %i\n", kk);
        for (int ii = 0; ii < n_obs; ii++) {
          bft_printf("    ");
          for (int jj = 0; jj < n_obs; jj++)
            if (!oi->obs_cov_is_diag)
              bft_printf("%.2f ", oi->obs_cov[ms->dim*(ii*n_obs + jj) + kk]);
            else if (oi->obs_cov_is_diag && ii == jj)
              bft_printf("%.2f ", oi->obs_cov[ms->dim*ii + kk]);
            else if (oi->obs_cov_is_diag && ii != jj)
              bft_printf("%.2f ",0.);
          bft_printf("\n");
        }
        bft_printf("\n");
      }

      bft_printf(" *End of processing variable %s\n\n\n", f->name);
    }
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all structures linked to data assimilation.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_finalize(void)
{
  if (_data_assim_activated) {
      cs_measures_sets_destroy();
      cs_interpol_grids_destroy();
      cs_at_opt_interps_destroy();
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log data assimilation.
 *
 * \param[in]  ms         pointer to measures set structure
 * \param[in]  oi         pointer to optimal interpolation structure
 * \param[in]  fname      field name
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_log(cs_measures_set_t     *ms,
                     cs_at_opt_interp_t    *oi,
                     cs_field_t            *f)
{
  int n_log_data = oi->n_log_data;
  int n_obs = ms->nb_measures;

  bft_printf("\n   * Variable %s\n", f->name);
  bft_printf("\n\n");
  bft_printf("  Number of observations : %i\n", n_obs);
  bft_printf("  Coordinates (x,y,z) :\n");

  for (int ii = 0; ii < n_obs; ii++) {
    bft_printf("   ");
    if (ii == n_log_data) {
      bft_printf("... ... ...\n");
      break;
    }
    for (int jj = 0; jj < 3; jj++)
      bft_printf("%.2f ", ms->coords[3*ii + jj]);
    bft_printf("\n");
  }
  if (oi->steady <= 0) {
    bft_printf("  Number of times : %i\n", oi->nb_times);
    bft_printf("  Times (s) :\n");
    bft_printf("   ");
    for (cs_lnum_t ii = 0; ii < oi->nb_times; ii++) {
      if (ii == n_log_data) {
        bft_printf("...");
        break;
      }
      bft_printf("%.2f ", oi->times_read[ii]);
    }
    bft_printf("\n");
  } else {
    bft_printf("  3DVar mode activated for iteration : %i\n", oi->steady);
  }

  bft_printf("  Measures :\n");

  for (int ii = 0; ii < n_obs; ii++) {
    if (ii == n_log_data) {
      for (cs_lnum_t jj = 0; jj < oi->nb_times; jj++) {
        bft_printf("... ");
        if (jj == n_log_data)
          break;
      }
      bft_printf("\n");
      break;
    }
    for (int kk = 0; kk < ms->dim; kk++) {
      bft_printf("   Comp. %i\n   ", kk);
      int ic = oi->measures_idx[ii*ms->dim+kk];
      for (cs_lnum_t jj = 0; jj < oi->nb_times; jj++) {
        if (jj == n_log_data) {
          bft_printf("...");
          break;
        }
        if (ic < oi->measures_idx[(ii+1)*ms->dim+kk]
            && (oi->steady > 0 || CS_ABS(oi->times[ic] - oi->times_read[jj]) < 1.e-12)) {
          bft_printf("%.2f ", ms->measures[ic]);
          ic++;
        }
        else
          bft_printf("NaN ");
      }
      bft_printf("\n  ");
    }
  }

  bft_printf("\n  Observation covariance error matrix :\n");
  for (int kk = 0; kk < ms->dim; kk++) {
    bft_printf("   Comp. %i\n", kk);
    for (int ii = 0; ii < n_obs; ii++) {
      bft_printf("    ");

      if (ii == n_log_data) {
        for (int jj = 0; jj < n_obs; jj++) {
          bft_printf("... ");
          if (jj == n_log_data)
            break;
        }
        bft_printf("\n");
        break;
      }

      for (int jj = 0; jj < n_obs; jj++) {
        if (jj == n_log_data) {
          bft_printf("...");
          break;
        }

        if (!oi->obs_cov_is_diag)
          bft_printf("%.2f ", oi->obs_cov[ms->dim*(ii*n_obs + jj) + kk]);
        else if (ii == jj)
          bft_printf("%.2f ", oi->obs_cov[ms->dim*ii + kk]);
        else
          bft_printf("%.2f ", 0.);
      }
      bft_printf("\n");
    }
  }

  bft_printf("  Interpolation type for Observation Operator : ");
  if (oi->interp_type == CS_AT_OPT_INTERP_P0)
    bft_printf("P0 (from cells containing the observations)\n");
  else if (oi->interp_type == CS_AT_OPT_INTERP_P1)
    bft_printf("P1 (from neighborhood used for gradient calculation)\n");
  if (oi->steady <= 0) {
    bft_printf("  Temporal window times of observations (s) : ");
    for (int ii = 0; ii < 4; ii++)
      bft_printf("%.2f ", oi->time_window[ii]);
    bft_printf("\n");
  }

  if (oi->steady <= 0)
    bft_printf("  Frequency of analysis calculation : each %i iteration(s)\n",oi->frequency);

  bft_printf("  Influence radii of observations (m, used for Model covariance error matrix) : %.2f %.2f\n", oi->ir[0], oi->ir[1]);
  for (int kk = 0; kk < f->dim; kk++) {
    bft_printf("  Relaxation factor (1/s) for comp. %i: %.1e\n",
               kk, oi->relax[kk]);
  }
  bft_printf("  Nudging type : ");

  if (oi->type_nudging == 1)
    bft_printf("explicit\n");
  else
    bft_printf("implicit\n");

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute source terms for a given variable and add them up to source
 *        term arrays.
 *
 * \param[in]  f_id     field id of variable
 * \param[in]  exp_st   array containing the explicit part of the source term
 * \param[in]  imp_st   array containing the implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_source_term(int        f_id,
                             cs_real_t *exp_st,
                             cs_real_t *imp_st)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells = m->n_cells;

  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  /* time step */
  const cs_time_step_t *ts = cs_glob_time_step;

  const int key_ms = cs_field_key_id("measures_set_id");
  const int key_oi = cs_field_key_id("opt_interp_id");
  const int key_oia = cs_field_key_id("opt_interp_analysis_id");

  cs_field_t *f = cs_field_by_id(f_id);

  int ms_id = cs_field_get_key_int(f, key_ms);
  cs_measures_set_t *ms = cs_measures_set_by_id(ms_id);

  int oi_id = cs_field_get_key_int(f, key_oi);
  cs_at_opt_interp_t *oi = cs_at_opt_interp_by_id(oi_id);

  int f_oia_id = cs_field_get_key_int(f, key_oia);
  cs_field_t *f_oia = cs_field_by_id(f_oia_id);

  bool nudging = false;

  if (oi->steady > 0 && ts->nt_cur > oi->steady) {//steady before start of OI
    /* 3DVar mode */
    nudging = true;
  } else if (   (oi->steady > 0 && ts->nt_cur == oi->steady) //steady, start of OI
             || (oi->steady <= 0 && (ts->nt_cur-1)%(oi->frequency) == 0)) {//unst.
#if _OI_DEBUG_
    bft_printf(" * Var %s", f->name);
    bft_printf(" - Iteration %i - Time %.2f\n\n", ts->nt_cur, ts->t_cur);
#endif

    /* if steady, first step of 3DVar mode */

    int **ao_idx = NULL; /* active obs. index */
    bool *inverse = NULL;
    BFT_MALLOC(inverse, ms->dim, bool);
    int *n_active_obs = cs_at_opt_interp_get_active_obs(ms,
                                                        oi,
                                                        f_oia,
                                                        &inverse,
                                                        &ao_idx);

    for (int kk = 0; kk < ms->dim; kk++) {
      if (n_active_obs[kk] > 0) {
        cs_at_opt_interp_compute_analysis(f,
                                          oi,
                                          f_oia,
                                          n_active_obs[kk],
                                          ao_idx[kk],
                                          inverse[kk],
                                          kk);
        nudging = true;
      }
      BFT_FREE(ao_idx[kk]);
    }
    BFT_FREE(inverse);
    BFT_FREE(ao_idx);
  }

  if (nudging) {
    int dim = f->dim;

    /* explicit */
    if (oi->type_nudging == 1) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t rovol = cell_f_vol[c_id]*CS_F_(rho)->val[c_id];
        for (int ii = 0; ii < dim; ii++) {
          exp_st[dim*c_id+ii] = exp_st[dim*c_id+ii]
                              +  rovol*oi->relax[ii]
                               * (f_oia->val[dim*c_id+ii]-f->val[dim*c_id+ii]);
        }
      }
    } else { /* implicit */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t rovol = cell_f_vol[c_id]*CS_F_(rho)->val[c_id];
        for (int ii = 0; ii < dim; ii++) {
          exp_st[dim*c_id+ii] = exp_st[dim*c_id+ii]
                              + rovol*oi->relax[ii]*f_oia->val[dim*c_id+ii];

          imp_st[dim*dim*c_id+(dim+1)*ii] = imp_st[dim*dim*c_id+(dim+1)*ii]
                                          - rovol*oi->relax[ii];
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
