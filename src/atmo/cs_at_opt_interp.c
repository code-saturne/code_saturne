/*============================================================================
 * Optimal interpolation
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

#include "cs_ext_neighborhood.h"
#include "cs_field.h"
#include "cs_map.h"
#include "cs_measures_util.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_math.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_at_opt_interp.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_at_opt_interp.c
        Optimal interpolation functions.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/* debugging switch */
#define _OI_DEBUG_ 1

/* maximum size of line in measures files */
#define MAX_LINE_SIZE 1000

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* optimal interpolation definitions */

static int  _n_opt_interps = 0;
static int  _n_opt_interps_max = 0;
static cs_at_opt_interp_t  *_opt_interps = NULL;
static cs_map_name_to_id_t *_opt_interps_map = NULL;
static int _p1_projection_needed = 0;

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

/*----------------------------------------------------------------------------
 * Compute observation operator (H) as a P0 interpolation operator.
 *
 * parameters:
 *   ms        <--   pointer to a measures set
 *   oi        <--   pointer to an optimal interpolation
 *   ig        <--   pointer to an interpol grid
 *----------------------------------------------------------------------------*/

static void
_p0_projection(cs_measures_set_t  *ms,
               cs_at_opt_interp_t *oi,
               cs_interpol_grid_t *ig)
{
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_3_t  *restrict cell_cen =
    (const cs_real_3_t *restrict)fvq->cell_cen;

  const int dim = ms->dim;
  const int stride = dim + 3; /* dimension of field + dimension of space */

  cs_lnum_t n_obs = ms->nb_measures;

  BFT_MALLOC(oi->model_to_obs_proj_idx, n_obs + 1, cs_lnum_t);
  cs_lnum_t *proj_idx = oi->model_to_obs_proj_idx;

  for (cs_lnum_t ii = 0; ii < n_obs + 1; ii++)
    proj_idx[ii] = 0;

  /* P0 => neighborhood of 1 cell per obs. */

  BFT_MALLOC(oi->model_to_obs_proj, n_obs * stride, cs_real_t);
  BFT_MALLOC(oi->model_to_obs_proj_c_ids, n_obs, cs_lnum_t);

  cs_real_t *proj = oi->model_to_obs_proj;
  cs_lnum_t *proj_c_ids = oi->model_to_obs_proj_c_ids;

  for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
    cs_lnum_t c_id0 = ig->cell_connect[ii];

    int r_id0 = 0;
    if (cs_glob_rank_id > -1) r_id0 = ig->rank_connect[ii];

    if (cs_glob_rank_id < 0 || cs_glob_rank_id == r_id0) {
      proj_idx[ii+1] += 1;
      for (cs_lnum_t kk = 0; kk < dim; kk++)
        (proj + stride*ii)[kk] = 1.;

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        (proj + stride*ii)[kk+dim] = cell_cen[c_id0][kk];
      proj_c_ids[ii] = c_id0;
    }
  }

  /* exchange size of neighborhood of each observation */

  if (cs_glob_rank_id > -1) {
    for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
      int r_id0 = ig->rank_connect[ii];
      cs_parall_bcast(r_id0, 1, CS_INT_TYPE, &(proj_idx[ii+1]));
    }
  }

  /* build observations neighborhood index */

  for (cs_lnum_t ii = 0; ii < n_obs; ii++)
    proj_idx[ii+1] += proj_idx[ii];

  /* exchange projection values on neighborhood of each observation */

  if (cs_glob_rank_id > -1) {
    for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
      cs_lnum_t s_idx = proj_idx[ii];
      cs_lnum_t l_size = stride*(proj_idx[ii+1] - s_idx);
      int r_id0 = ig->rank_connect[ii];
      cs_parall_bcast(r_id0, l_size, CS_DOUBLE, proj + stride*s_idx);
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute observation operator (H) as a P1 interpolation operator.
 *
 * parameters:
 *   ms        <--   pointer to a measures set
 *   oi        <--   pointer to a optimal interpolation
 *   ig        <--   pointer to an interpol grid
 *----------------------------------------------------------------------------*/

static void
_p1_projection(cs_measures_set_t  *ms,
               cs_at_opt_interp_t *oi,
               cs_interpol_grid_t *ig)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_halo_type_t halo_type = mesh->halo_type;

  /* cells -> cells connectivity (standard) */

  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)ma->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict)ma->cell_cells;

  /* cells -> cells connectivity (extended) */

  const cs_lnum_t *restrict cell_cells_e_idx
    = (const cs_lnum_t *restrict)ma->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_cells_e
    = (const cs_lnum_t *restrict)ma->cell_cells_e;

  const cs_real_3_t  *restrict cell_cen =
    (const cs_real_3_t *restrict)fvq->cell_cen;

  cs_lnum_t n_obs = ms->nb_measures;

  BFT_MALLOC(oi->model_to_obs_proj_idx, n_obs + 1, cs_lnum_t);
  cs_lnum_t *proj_idx = oi->model_to_obs_proj_idx;

  for (cs_lnum_t ii = 0; ii < n_obs + 1; ii++)
    proj_idx[ii] = 0;

  /* P1 => neighborhood =  standard + eventually (partially) extended */

  /* count size of neighborhood of each observation */

  for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
    int r_id0 = 0;
    if (cs_glob_rank_id > -1) r_id0 = ig->rank_connect[ii];

    if (cs_glob_rank_id < 0 || cs_glob_rank_id == r_id0) {
      proj_idx[ii+1] += 1;
      cs_lnum_t c_id0 = ig->cell_connect[ii];

      /* standard neighborhood */

      for (cs_lnum_t jj = cell_cells_idx[c_id0];
           jj < cell_cells_idx[c_id0+1];
           jj++)
        proj_idx[ii+1] += 1;

      /* (partially) extended neighborhood */

      if (halo_type == CS_HALO_EXTENDED)
        for (cs_lnum_t jj = cell_cells_e_idx[c_id0];
             jj < cell_cells_e_idx[c_id0+1];
             jj++)
          proj_idx[ii+1] += 1;
    }
  }

  /* exchange size of neighborhood of each observation */

  if (cs_glob_rank_id > -1) {
    for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
      int r_id0 = ig->rank_connect[ii];
      cs_parall_bcast(r_id0, 1, CS_INT_TYPE, &(proj_idx[ii+1]));
    }
  }

  /* now allocate projection matrix values and cell ids */

  const int dim = ms->dim;
  const int stride = dim + 3; /* dimension of field + dimension of space */

  BFT_MALLOC(oi->model_to_obs_proj, proj_idx[n_obs] * stride, cs_real_t);
  BFT_MALLOC(oi->model_to_obs_proj_c_ids, proj_idx[n_obs], cs_lnum_t);
  cs_real_t *proj = oi->model_to_obs_proj;
  cs_lnum_t *proj_c_ids = oi->model_to_obs_proj_c_ids;

  /* compute max. size and build observations neighborhood index */

  cs_lnum_t n_max_size = 0;
  for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
    n_max_size = CS_MAX(proj_idx[ii+1], n_max_size);
    proj_idx[ii+1] += proj_idx[ii];
  }

  /* compute projection matrix coefficients */

  cs_real_t *dist = NULL;
  BFT_MALLOC(dist, n_max_size, cs_real_t);

  for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
    int r_id0 = 0;
    if (cs_glob_rank_id > -1) r_id0 = ig->rank_connect[ii];

    if (cs_glob_rank_id < 0 || cs_glob_rank_id == r_id0) {
      cs_lnum_t c_id0 = ig->cell_connect[ii];

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        (proj + stride*proj_idx[ii])[kk+dim] = cell_cen[c_id0][kk];

      proj_c_ids[proj_idx[ii]] = c_id0;
      cs_lnum_t ccount = 1;

      /* standard neighborhood */

      for (cs_lnum_t jj = cell_cells_idx[c_id0];
           jj < cell_cells_idx[c_id0+1];
           jj++) {
        cs_lnum_t tot_ccount = proj_idx[ii] + ccount;
        cs_lnum_t c_id1 = cell_cells[jj];

        for (cs_lnum_t kk = 0; kk < 3; kk++)
          (proj + stride*tot_ccount)[kk+dim] = cell_cen[c_id1][kk];

        proj_c_ids[tot_ccount] = c_id1;
        ccount += 1;
      }

      /* (partially) extended neighborhood */

      if (halo_type == CS_HALO_EXTENDED) {
        for (cs_lnum_t jj = cell_cells_idx[c_id0];
             jj < cell_cells_idx[c_id0+1];
             jj++) {
          cs_lnum_t tot_ccount = proj_idx[ii] + ccount;
          cs_lnum_t c_id1 = cell_cells_e[jj];

          for (cs_lnum_t kk = 0; kk < 3; kk++)
            (proj + stride*tot_ccount)[kk+dim] = cell_cen[c_id1][kk];

          proj_c_ids[tot_ccount] = c_id1;
          ccount += 1;
        }
      }

      /* interpolation coefficients */

      cs_real_t xobs = ms->coords[ii*3    ];
      cs_real_t yobs = ms->coords[ii*3 + 1];
      cs_real_t zobs = ms->coords[ii*3 + 2];

      for (cs_lnum_t jj = 0; jj < ccount; jj++) {
        cs_lnum_t tot_ccount = proj_idx[ii] + jj;
        cs_real_t x = (proj + stride*tot_ccount)[dim  ];
        cs_real_t y = (proj + stride*tot_ccount)[dim+1];
        cs_real_t z = (proj + stride*tot_ccount)[dim+2];
        dist[jj] = sqrt( cs_math_sq(xobs - x)
                        +cs_math_sq(yobs - y)
                        +cs_math_sq(zobs - z));
      }

      /* handle case of null distance between obs and one of the
         cell center in its neighborhood */

      bool obs_on_cell_center = false;

      for (cs_lnum_t jj = 0; jj < ccount; jj++) {
        if (dist[jj] < cs_math_epzero) {
          for (int ll = 0; ll < 3; ll++)
            (proj + stride*(proj_idx[ii] + jj))[ll] = 1.;

          for (cs_lnum_t kk = 0; kk < ccount; kk++) {
            if (kk != jj)
              for (int ll = 0; ll < dim; ll++)
                (proj + stride*(proj_idx[ii] + kk))[ll] = 0.;
          }

          obs_on_cell_center = true;
          break;
        }
      }

      /* general case: harmonic mean of interpolation coefficients for 1 obs
         must be equal to 1 */

      if (!obs_on_cell_center) {
        for (cs_lnum_t jj = 0; jj < ccount; jj++) {
          cs_real_t v = 0.;
          for (cs_lnum_t kk = 0; kk < ccount; kk++)
            if (kk != jj)
              v += 1 / dist[kk];
          v = v*dist[jj] + 1.;
          for (int ll = 0; ll < dim; ll++)
            (proj + stride*(proj_idx[ii] + jj))[ll] = 1/v;
        }
      }
    }
  }

  BFT_FREE(dist);

  /* exchange projection values on neighborhood of each observation */

  if (cs_glob_rank_id > -1) {
    for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
      cs_lnum_t s_idx = proj_idx[ii];
      cs_lnum_t l_size = stride*(proj_idx[ii+1] - s_idx);
      int r_id0 = ig->rank_connect[ii];
      cs_parall_bcast(r_id0, l_size, CS_DOUBLE, proj + stride*s_idx);
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute coefficient bij of model covariance matrix (B)
 * of size n_cells*n_cells.
 *
 * parameters:
 *   xi     <-- x coordinate of point I
 *   yi     <-- y coordinate of point I
 *   zi     <-- z coordinate of point I
 *   xj     <-- x coordinate of point J
 *   yj     <-- y coordinate of point J
 *   zj     <-- z coordinate of point J
 *   ir_xy2 <-- square of influence radius with respect to x and y
 *   ir_z2  <-- square of influence radius with respect to z
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_b_matrix(cs_real_t  xi,
          cs_real_t  yi,
          cs_real_t  zi,
          cs_real_t  xj,
          cs_real_t  yj,
          cs_real_t  zj,
          cs_real_t  ir_xy2,
          cs_real_t  ir_z2)
{
  cs_real_t dist = sqrt( ( cs_math_sq(xi - xj)
                         + cs_math_sq(yi - yj) )/ir_xy2
                       + cs_math_sq(zi - zj)/ir_z2);

  return (1. + dist) * exp(-dist);
}

/*----------------------------------------------------------------------------
 * Assemble full matrix HB(H)t+R.
 *----------------------------------------------------------------------------*/

static cs_real_t *
_assembly_adding_obs_covariance(cs_measures_set_t  *ms,
                                cs_at_opt_interp_t *oi,
                                int                *ao_idx,
                                int                 n_active_obs,
                                int                 mc_id)
{
  cs_lnum_t n_obs = ms->nb_measures;
  cs_real_t *obs_cov = oi->obs_cov;
  cs_real_t *b_proj = oi->b_proj;
  cs_real_t *a, r;

  int m_dim = ms->dim;
  int a_l_size = n_active_obs;
  int a_size = cs_math_sq(a_l_size);
  BFT_MALLOC(a, a_size, cs_real_t);

  /* filling in the full matrix */

  for (cs_lnum_t ii = 0; ii < a_size; ii++)
    a[ii] = 0.;

  for (cs_lnum_t ii = 0; ii < n_active_obs; ii++) {
    for (cs_lnum_t jj = 0; jj < n_active_obs; jj++) {

      int id = ii*a_l_size + jj;
      a[id] = b_proj[m_dim*(ao_idx[ii] * n_obs + ao_idx[jj])+mc_id];

      /* time weighting of variances */
      if (ii == jj) {
        if (!oi->obs_cov_is_diag) {
          r = obs_cov[m_dim*(ao_idx[ii] * n_obs + ao_idx[jj])+mc_id];
        } else if (oi->obs_cov_is_diag) {
          r = obs_cov[m_dim*ao_idx[ii]+mc_id];
        }

        if (oi->steady <= 0) {
          a[id] += (r + 1.) / oi->time_weights[m_dim*ao_idx[ii]+mc_id] - 1.;
        } else {
          a[id] += r;
        }

      } else if (!oi->obs_cov_is_diag) {
        a[id] += obs_cov[m_dim*(ao_idx[ii] * n_obs + ao_idx[jj])+mc_id];
      }
    }
  }

  return a;
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an optimal interpolation descriptor.
 *
 * \param[in]  name   optimal interpolation name
 *
 */
/*----------------------------------------------------------------------------*/

cs_at_opt_interp_t *
cs_at_opt_interp_create(const char   *name)
{
  const char *addr_0 = NULL, *addr_1 = NULL;

  cs_at_opt_interp_t *oi =  NULL;

  /* Initialize if necessary */

  if (_opt_interps_map == NULL)
    _opt_interps_map = cs_map_name_to_id_create();
  else
    addr_0 = cs_map_name_to_id_reverse(_opt_interps_map, 0);

  if (strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Defining an optimal interpolation requires a name."));

  /* Find or insert entry in map */

  int opt_interp_id = cs_map_name_to_id(_opt_interps_map, name);

  /* Move name pointers of previous measure set if necessary
     (i.e. reallocation of map names array) */

  addr_1 = cs_map_name_to_id_reverse(_opt_interps_map, 0);

  if (addr_1 != addr_0) {
    int i;
    ptrdiff_t addr_shift = addr_1 - addr_0;
    for (i = 0; i < opt_interp_id; i++)
      (_opt_interps + i)->name += addr_shift;
  }

  bool reall = true;
  if (opt_interp_id == _n_opt_interps) {
    _n_opt_interps = opt_interp_id + 1;
    reall = false;
  }

  /* Reallocate optimal interpolations if necessary */

  if (_n_opt_interps > _n_opt_interps_max) {
    if (_n_opt_interps_max == 0)
      _n_opt_interps_max = 8;
    else
      _n_opt_interps_max *= 2;
    BFT_REALLOC(_opt_interps, _n_opt_interps_max, cs_at_opt_interp_t);
  }

  /* Assign optimal interpolation */

  oi = _opt_interps + opt_interp_id;

  oi->name = cs_map_name_to_id_reverse(_opt_interps_map, opt_interp_id);

  oi->id = opt_interp_id;

  oi->ig_id = -1;

  if (!reall) {
    oi->b_proj = NULL;
    oi->relax = NULL;
    oi->times = NULL;
    oi->times_read = NULL;
    oi->obs_cov = NULL;
    oi->measures_idx = NULL;
    oi->model_to_obs_proj = NULL;
    oi->model_to_obs_proj_idx = NULL;
    oi->model_to_obs_proj_c_ids = NULL;
    oi->active_time = NULL;
    oi->time_weights = NULL;
    oi->time_window = NULL;
  }
  else {
    BFT_FREE(oi->b_proj);
    BFT_FREE(oi->relax);
    BFT_FREE(oi->times);
    BFT_FREE(oi->times_read);
    BFT_FREE(oi->obs_cov);
    BFT_FREE(oi->measures_idx);
    BFT_FREE(oi->model_to_obs_proj);
    BFT_FREE(oi->model_to_obs_proj_idx);
    BFT_FREE(oi->model_to_obs_proj_c_ids);
    BFT_FREE(oi->active_time);
    BFT_FREE(oi->time_weights);
    BFT_FREE(oi->time_window);
  }

  return oi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to an optimal interpolation based on its id.
 *
 * This function requires that an optimal interpolation of the given id is
 * defined.
 *
 * \param[in]  id   optimal interpolation id
 *
 * \return  pointer to the optimal interpolation structure.
 */
/*----------------------------------------------------------------------------*/

cs_at_opt_interp_t  *
cs_at_opt_interp_by_id(int  id)
{
  if (id > -1 && id < _n_opt_interps)
    return _opt_interps + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Optimal interpolation with id %d is not defined."), id);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to an optimal interpolation based on its name.
 *
 * This function requires that an optimal interpolation of the given name is
 * defined.
 *
 * \param[in]  name   optimal interpolation name
 *
 * \return  pointer to the optimal interpolation structure.
 */
/*----------------------------------------------------------------------------*/

cs_at_opt_interp_t  *
cs_at_opt_interp_by_name(const char  *name)
{
  int id = cs_map_name_to_id_try(_opt_interps_map, name);

  if (id > -1)
    return _opt_interps + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Optimal interpolation \"%s\" is not defined."), name);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all defined optimal interpolations.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interps_destroy(void)
{
  for (int i = 0; i < _n_opt_interps; i++) {
    cs_at_opt_interp_t  *oi = _opt_interps + i;
    BFT_FREE(oi->b_proj);
    BFT_FREE(oi->relax);
    BFT_FREE(oi->obs_cov);
    BFT_FREE(oi->times);
    BFT_FREE(oi->times_read);
    BFT_FREE(oi->measures_idx);
    BFT_FREE(oi->model_to_obs_proj);
    BFT_FREE(oi->model_to_obs_proj_idx);
    BFT_FREE(oi->model_to_obs_proj_c_ids);
    BFT_FREE(oi->active_time);
    BFT_FREE(oi->time_weights);
    BFT_FREE(oi->time_window);
  }

  BFT_FREE(_opt_interps);

  cs_map_name_to_id_destroy(&_opt_interps_map);

  _n_opt_interps = 0;
  _n_opt_interps_max = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read an optimal interpolation file for a given variable
 *        and fill in the matching measures set and optimal interpolation
 *        structures.
 *
 * \param[in]  filename   name of interpolation file
 * \param[in]  ms         measures set structure
 * \param[in]  oi         optimal interpolation structure
 * \param[in]  f_dim      dimension of field
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_read_file(char const           filename[50],
                           cs_measures_set_t   *ms,
                           cs_at_opt_interp_t  *oi,
                           const int            f_dim)
{
  cs_real_t val;
  char line[MAX_LINE_SIZE];
  FILE* fichier = NULL;

  BFT_MALLOC(oi->relax, f_dim, cs_real_t);
  for (int kk = 0; kk < f_dim; kk++)
    oi->relax[kk] = 0.; /* nudging disabled by default */

  /* Defaults */
  ms->nb_measures = 0;
  oi->nb_times = 0;
  oi->ir[0] = 100.;
  oi->ir[1] = 100.;
  oi->n_log_data = 10;
  oi->interp_type = CS_AT_OPT_INTERP_P0;
  oi->steady = -1;
  oi->frequency = 1;
  oi->type_nudging = 0;

  /* First reading - test if file exists */
  fichier = fopen(filename, "r");
  if (fichier == NULL) {
    bft_printf(" Error when opening file %s - exit\n",filename);
    cs_exit(EXIT_FAILURE);
  }
  fclose(fichier);

  /* Second reading */
  fichier = fopen(filename, "r");
  while (fgets(line, MAX_LINE_SIZE, fichier)) {
    /* Reading Number of observations */
    if (strncmp(line, "_nobs_", 6) == 0) {
      fscanf(fichier, "%i", &(ms->nb_measures));

#if _OI_DEBUG_
      bft_printf("   *Reading _nobs_ : %i\n", ms->nb_measures);
#endif

      if (ms->nb_measures <= 0) {
        bft_printf(" File %s. The number of observations (_nobs_) must be "
                   "strictly positive.\n"
                   " Otherwise turn off nudging for the corresponding variable"
                   "- exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }
    }

    /* Reading dimension of measures */
    if (strncmp(line, "_dim_", 5) == 0) {
      fscanf(fichier, "%i", &(ms->dim));

#if _OI_DEBUG_
      bft_printf("   *Reading _dim_ : %i\n", ms->dim);
#endif

      if (ms->dim <= 0) {
        bft_printf(" File %s. The observation dimension (_dim_) must be "
                   "strictly positive.\n"
                   " Otherwise turn off nudging for the corresponding variable"
                   "- exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }
      BFT_MALLOC(ms->comp_ids, ms->dim, int);
    }


    /* Reading observation times */
    if (strncmp(line, "_times_", 7) == 0) {
      fscanf(fichier, "%i", &(oi->nb_times));

#if _OI_DEBUG_
      bft_printf("   * Reading _times_ : %i\n    ", oi->nb_times);
#endif

      if (oi->nb_times < 0) {
        bft_printf(" File %s. The number of times (_times_) must be positive.\n"
                   " Otherwise turn off nudging for the corresponding variable"
                   "- exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }

      /* 3DVar mode */
      if (oi->nb_times == 0) {
        oi->nb_times = 1;
        fscanf(fichier, "%i", &(oi->steady));

#if _OI_DEBUG_
        bft_printf("3DVar mode active for iteration %i\n",oi->steady);
#endif

        if (oi->steady <= 0) {
          bft_printf(" File %s. The iteration identifier (_times_ with 0 time) "
                     "must be strictly positive for 3DVar mode.\n"
                     " Otherwise turn off nudging for the corresponding variable"
                     " - exit\n", filename);
          cs_exit(EXIT_FAILURE);
        }
      }

      /* Time management */
      else {
        BFT_MALLOC(oi->times_read, oi->nb_times, cs_real_t);

        for (cs_lnum_t i = 0; i < oi->nb_times; i++)
          fscanf(fichier, "%lf", oi->times_read+i);

#if _OI_DEBUG_
        for (cs_lnum_t i = 0; i < oi->nb_times; i++)
          bft_printf("%.2f ", oi->times_read[i]);
        bft_printf("\n");
#endif

      }
    }

  }
  fclose(fichier);

  int n_obs = ms->nb_measures;

  /* Third reading */
  fichier = fopen(filename, "r");

  /* initialize count */
  int *_n_readable_measures;
  BFT_MALLOC(_n_readable_measures, ms->dim+1, int);
  for (int kk = 0; kk < ms->dim+1; kk++)
    _n_readable_measures[kk] = 0;

  while (fgets(line, MAX_LINE_SIZE, fichier)) {

    /* counting non nan measures */

    if (strncmp(line, "_measures_", 10) == 0) {
      for (int kk = 0; kk < ms->dim-1; kk++) {
        fscanf(fichier, "%i ", &ms->comp_ids[kk]);
      }
      fscanf(fichier, "%i", &ms->comp_ids[ms->dim-1]);

      for (int ii = 0; ii < n_obs; ii++)
        for (int jj = 0; jj < oi->nb_times; jj++) {
          for (int kk = 0; kk < ms->dim-1; kk++) {
            fscanf(fichier, "%lf ", &val);
            if (!isnan(val))
              _n_readable_measures[kk+1] += 1;
          }
          fscanf(fichier, "%lf", &val); /* read after because of ending space */
          if (!isnan(val))
            _n_readable_measures[ms->dim] += 1;
        }
    }
  }
  fclose(fichier);

  int _tot_n_readable_measures = 0;
  for (int kk = 0; kk < ms->dim + 1; kk++)
    _tot_n_readable_measures += _n_readable_measures[kk];

  /* 4th reading */
  fichier = fopen(filename, "r");
  while (fgets(line, MAX_LINE_SIZE, fichier)) {

    /*Reading influence radius*/
    if (strncmp(line, "_L_", 3) == 0) {
      fscanf(fichier, "%lf", &(oi->ir[0]));
      fscanf(fichier, "%lf", &(oi->ir[1]));

#if _OI_DEBUG_
      bft_printf("   * Reading _L_ : %.2f %.2f\n", oi->ir[0], oi->ir[1]);
#endif

      if (oi->ir[0] <= 0. || oi->ir[1] <= 0.) {
        bft_printf(" File %s. The observation influence radii (_L_ in m)"
                   " must be strictly positive - exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }
    }

    /* Reading relaxation time */
    if (strncmp(line, "_t_", 3) == 0) {
      cs_real_t *tau = NULL;
      BFT_MALLOC(tau, ms->dim, cs_real_t);
      for (int kk = 0; kk < ms->dim-1; kk++)
        fscanf(fichier, "%lf ", tau+kk);
      fscanf(fichier, "%lf", tau+ms->dim-1);

#if _OI_DEBUG_
      bft_printf("   * Reading _t_ :");
      for (int kk = 0; kk < ms->dim; kk++)
        bft_printf(" %.1e", tau[kk]);
      bft_printf("\n");
#endif

      for (int kk = 0; kk < ms->dim; kk++) {
        if (tau[kk] <= 0.) {
          bft_printf(" File %s. Each relaxation time (_t_ in s) must be"
                     " strictly positive - exit\n", filename);
          cs_exit(EXIT_FAILURE);
        }
        oi->relax[ms->comp_ids[kk]] = 1/tau[kk];
      }
    }

    /* Reading the number of obs/times to print in the listing */
    if (strncmp(line, "_n_log_data_", 12) == 0) {
      fscanf(fichier, "%i", &(oi->n_log_data));
      if (oi->n_log_data < 0)
        oi->n_log_data = 0;

#if _OI_DEBUG_
      bft_printf("   * Reading _n_log_data_ : %i\n", oi->n_log_data);
#endif

    }

    /* Reading the frequency of analysis computation */
    if (strncmp(line, "_frequency_", 11) == 0) {
      fscanf(fichier, "%i", &(oi->frequency));
      if (oi->frequency < 1)
        oi->frequency = 1;

#if _OI_DEBUG_
      bft_printf("   * Reading _frequency_ : %i\n", oi->frequency);
#endif

    }

    /* Reading the type of nudging applied */
    if (strncmp(line, "_nudging_", 9) == 0) {
      fgets(line, MAX_LINE_SIZE, fichier);
      if (strncmp(line, "explicit", 8) == 0) {
        oi->type_nudging = 1;

#if _OI_DEBUG_
        bft_printf("   * Reading _nudging_ : explicit\n");
#endif

      }
      else if (strncmp(line, "implicit", 8) == 0) {
        oi->type_nudging = 0;

#if _OI_DEBUG_
        bft_printf("   * Reading _nudging_ : implicit\n");
#endif

      }
      else {
        bft_printf(" File %s. Section _nudging_: wrong key"
                    " (admissible keys: explicit, implicit) - exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }
    }

    /* Reading interpolation type */
    if (strncmp(line, "_interp_", 8) == 0) {
      fgets(line, MAX_LINE_SIZE, fichier);
      if (strncmp(line, "P0", 2) == 0) {
        oi->interp_type = CS_AT_OPT_INTERP_P0;

#if _OI_DEBUG_
        bft_printf("   * Reading _interp_ : P0\n");
#endif

      }
      else if (strncmp(line, "P1", 2) == 0) {
        oi->interp_type = CS_AT_OPT_INTERP_P1;
        _p1_projection_needed = 1;
#if _OI_DEBUG_
        bft_printf("   * Reading _interp_ : P1\n");
#endif

      }
      else {
        bft_printf(" File %s. Section _interp_: wrong key"
                   " (admissible keys: P0, P1) - exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }
    }

    /* Reading coordinates */
    if (strncmp(line,"_coordinates_", 13) == 0) {
      BFT_MALLOC(ms->coords, n_obs*3, cs_real_t);
      for (cs_lnum_t i = 0; i < n_obs; i++)
        fscanf(fichier, "%lf %lf %lf", ms->coords+3*i,
                                       ms->coords+3*i+1,
                                       ms->coords+3*i+2);

#if _OI_DEBUG_
      bft_printf("   * Reading _coordinates_\n");
      for (cs_lnum_t i = 0; i < n_obs; i++)
        bft_printf("    %.2f %.2f %.2f\n", ms->coords[3*i],
                                           ms->coords[3*i+1],
                                           ms->coords[3*i+2]);
#endif

    }

    /* reading measures */
    if (strncmp(line, "_measures_", 10) == 0) {

      BFT_MALLOC(oi->measures_idx, ms->dim*(n_obs + 1), int);
      BFT_MALLOC(ms->measures, _tot_n_readable_measures, cs_real_t);

      if (oi->steady <= 0)
        BFT_MALLOC(oi->times, _tot_n_readable_measures, cs_real_t);

      /* initialising index list */
      for (int ii = 0; ii < ms->dim*(n_obs + 1); ii++)
        oi->measures_idx[ii] = 0;

      /* build index */
      for (int kk = 0; kk < ms->dim; kk++)
        _n_readable_measures[kk+1] += _n_readable_measures[kk];

#if _OI_DEBUG_
      bft_printf("   * Reading _measures_\n");
#endif

      /* try to read measures value at each defined time and
         count successful reads by measure (measures by time index) */

      /* read again components to skip the line */

      for (int kk = 0; kk < ms->dim-1; kk++) {
        fscanf(fichier, "%i ", &ms->comp_ids[kk]);
      }
      fscanf(fichier, "%i", &ms->comp_ids[ms->dim-1]);

      /* read measures */

      int *_n_readable_measures_c;
      BFT_MALLOC(_n_readable_measures_c, ms->dim, int);
      for (int kk = 0; kk < ms->dim; kk++)
        _n_readable_measures_c[kk] = 0;

      for (int ii = 0; ii < n_obs; ii++) {
#if _OI_DEBUG_
        bft_printf("    ");
#endif
        for (int jj = 0; jj < oi->nb_times; jj++) {
          val = 0;
          for (int kk = 0; kk < ms->dim-1; kk++) {
            fscanf(fichier, "%lf ", &val);
            if (!isnan(val)) {
              int cc = _n_readable_measures[kk] + _n_readable_measures_c[kk];

              ms->measures[cc] = val;
              if (oi->steady <= 0)
                oi->times[cc] = oi->times_read[jj];
              oi->measures_idx[(ii+1)*ms->dim+kk] += 1;
              _n_readable_measures_c[kk] += 1;
#if _OI_DEBUG_
              bft_printf("%.2f ", val);
#endif
            }
          }
          val = 0;
          fscanf(fichier, "%lf", &val); /* read after because of ending space */
          if (!isnan(val)) {
            int cc = _n_readable_measures[ms->dim-1]
                   + _n_readable_measures_c[ms->dim-1];

            ms->measures[cc] = val;
            if (oi->steady <= 0)
              oi->times[cc] = oi->times_read[jj];
            oi->measures_idx[(ii+1)*ms->dim+ms->dim-1] += 1;
            _n_readable_measures_c[ms->dim-1] += 1;
#if _OI_DEBUG_
            bft_printf("%.2f", val);
#endif
          }
        }
#if _OI_DEBUG_
        bft_printf("\n");
#endif
      }

      for (int kk = 0; kk < ms->dim; kk++) {
        for (int ii = 0; ii < n_obs; ii++) {
          oi->measures_idx[(ii+1)*ms->dim+kk] +=
            oi->measures_idx[ii*ms->dim+kk];
        }
        oi->measures_idx[kk] = 0;
      }

      for (int kk = 0; kk < ms->dim; kk++) {
        for (int ii = 0; ii < n_obs + 1; ii++) {
          oi->measures_idx[ii*ms->dim+kk] += _n_readable_measures[kk];
        }
      }

      BFT_FREE(_n_readable_measures_c);
      BFT_FREE(_n_readable_measures);

#if _OI_DEBUG_
      bft_printf("   * Building index list\n    ");
      for (int kk = 0; kk < ms->dim; kk++) {
        bft_printf("    Comp. %i\n    ", kk);
        for (int ii = 0; ii < n_obs + 1; ii++) {
          bft_printf("%i ", oi->measures_idx[ms->dim*ii+kk]);
        }
        bft_printf("\n\n");
      }

      bft_printf("\n");
      bft_printf("   * Building measures list\n");
      for (int kk = 0; kk < ms->dim; kk++) {
        bft_printf("    Comp. %i\n    ", kk);
        for (int ii = 0; ii < n_obs; ii++) {
          for (int jj = oi->measures_idx[ms->dim*ii+kk];
               jj < oi->measures_idx[ms->dim*(ii+1)+kk];
               jj++) {
            bft_printf("%.2f ", ms->measures[jj]);
          }
          bft_printf("\n    ");
        }
        bft_printf("\n\n");
      }

      if (oi->steady <= 0) {
        bft_printf("   * Building times list\n");
        for (int kk = 0; kk < ms->dim; kk++) {
          bft_printf("    Comp. %i\n    ", kk);
          for (int ii = 0; ii < n_obs; ii++) {
            for (int jj = oi->measures_idx[ms->dim*ii+kk];
                 jj < oi->measures_idx[ms->dim*(ii+1)+kk];
                 jj++) {
              bft_printf("%.2f ", oi->times[jj]);
            }
            bft_printf("\n    ");
          }
          bft_printf("\n\n");
        }
      }
#endif
    }

    /* Reading observations error covariances */
    if (strncmp(line, "_errors_", 8) == 0) {
      fgets(line, MAX_LINE_SIZE, fichier);

      if (strncmp(line, "diagonal", 8) == 0) {
        BFT_MALLOC(oi->obs_cov, ms->dim*n_obs, cs_real_t);
        oi->obs_cov_is_diag = true;

        for (int ii = 0; ii < ms->dim*n_obs; ii++)
          fscanf(fichier, "%lf", oi->obs_cov + ii);

#if _OI_DEBUG_
        bft_printf("   * Reading _errors_ : diagonal\n    ");
        for (int ii = 0; ii < ms->dim*n_obs; ii++)
          bft_printf("%.2f ", oi->obs_cov[ii]);
        bft_printf("\n");
#endif

      } else if (strncmp(line,"full", 4) == 0) {
        BFT_MALLOC(oi->obs_cov, ms->dim*n_obs*ms->dim*n_obs, cs_real_t);
        oi->obs_cov_is_diag = false;

        for (int ii = 0; ii < n_obs; ii++)
          for (int jj = 0; jj < n_obs; jj++)
            for (int kk = 0; kk < ms->dim; kk++)
              fscanf(fichier, "%lf", oi->obs_cov + ms->dim*(ii*n_obs + jj) + kk);

#if _OI_DEBUG_
        bft_printf("   * Reading _errors_ : full\n");
        for (int ii = 0; ii < n_obs; ii++) {
          bft_printf("    ");
          for (int jj = 0; jj < n_obs; jj++)
            for (int kk = 0; kk < ms->dim; kk++)
              bft_printf("%.2f ", oi->obs_cov[ms->dim*(ii*n_obs + jj) + kk]);
          bft_printf("\n");
        }
#endif

      } else {
        bft_printf(" File %s. Section _errors_: wrong key (admissible"
                   " keys: diagonal, full) - exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }
    }

    /* Reading observations time window */
    if (strncmp(line, "_temps_", 7) == 0) {
      BFT_MALLOC(oi->time_window, 4, cs_real_t);
      fgets(line, MAX_LINE_SIZE, fichier);
      if (strncmp(line, "sym", 3) == 0) {
        fscanf(fichier, "%lf", oi->time_window + 2);
        fscanf(fichier, "%lf", oi->time_window + 3);
        oi->time_window[1] = - oi->time_window[2];
        oi->time_window[0] = - oi->time_window[3];
      } else if (strncmp(line,"asym", 4) == 0) {
        for (cs_lnum_t i = 0; i < 4; i++)
          fscanf(fichier, "%lf", oi->time_window + i);
      } else {
        bft_printf(" File %s. Section _temps_: wrong key (admissible"
                   " keys: sym, asym) - exit\n", filename);
        cs_exit(EXIT_FAILURE);
      }

      /* Automatic corrections */
      for (cs_lnum_t i = 0; i < 2; i++) {
        if (oi->time_window[i] >= 0)
          oi->time_window[i] = - oi->time_window[i];
        if (oi->time_window[2 + i] <= 0)
          oi->time_window[2 + i] = - oi->time_window[2+i];
      }

      if (oi->time_window[0] >= oi->time_window[1])
        oi->time_window[0] = oi->time_window[1];
      if (oi->time_window[3] <= oi->time_window[2])
        oi->time_window[3] = oi->time_window[2];

#if _OI_DEBUG_
      if (oi->steady <= 0) {
        bft_printf("   * Reading _temps_\n    ");
        for (cs_lnum_t i = 0; i < 4; i++)
          bft_printf("%.2f ", oi->time_window[i]);
        bft_printf("\n");
      }
#endif

    }
  }
  fclose(fichier);

#if _OI_DEBUG_
  bft_printf("\n   * File %s closed\n\n", filename);
#endif

  /* Verifications */
  if (ms->coords == NULL) {
    bft_printf(" File %s. The observation coordinates are missing"
               " (_coordinates_) - exit\n", filename);
    cs_exit(EXIT_FAILURE);
  }
  if (oi->steady <= 0 && oi->times == NULL) {
    bft_printf(" File %s. The observation times are missing"
               " (_times_) - exit\n", filename);
    cs_exit(EXIT_FAILURE);
  }
  if (ms->measures == NULL) {
    bft_printf(" File %s. The measures are missing (_measures_) - exit\n", filename);
    cs_exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 1 if a p1 projection has been enabled for at least one
 *        optimal interpolation. This function is used to determine if
 *        extended neighborhood is needed.
 *
 * \return  1 if a p1 proj. is needed, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_at_opt_interp_is_p1_proj_needed(void)
{
  return _p1_projection_needed;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief (re)Allocate and fill in an optimal interpolation structure from an
 *        optimal interpolation file.
 *
 * \param[in]  oi  pointer to the optimal interpolation
 * \param[in]  ms  pointer to the associated measures set
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_map_values(cs_at_opt_interp_t *oi,
                            cs_measures_set_t  *ms)
{
  int n_obs = ms->nb_measures;

  /* Initializing remaining variables */
  /* Observations errors */
  if (oi->obs_cov == NULL) {
    oi->obs_cov_is_diag = true;
    BFT_MALLOC(oi->obs_cov, ms->dim*n_obs, cs_real_t);
    for (int ii = 0; ii < ms->dim*n_obs; ii++)
      oi->obs_cov[ii] = 1.;
  }

  /* Time window */
  if (oi->time_window == NULL) {
    BFT_MALLOC(oi->time_window, 4, cs_real_t);
    oi->time_window[2] = 300.;
    oi->time_window[3] = 360.;
    oi->time_window[1] = -oi->time_window[2];
    oi->time_window[0] = -oi->time_window[3];
  }

  /* Active Times */
  BFT_MALLOC(oi->active_time, ms->dim*n_obs, int);
  for (int ii = 0; ii < n_obs; ii++)
    for (int kk = 0; kk < ms->dim; kk++)
      oi->active_time[ms->dim*ii+kk] = oi->measures_idx[ms->dim*ii+kk];

  /* Initialising time weighting coefficients */
  if (oi->steady <= 0) {
    BFT_MALLOC(oi->time_weights, ms->dim*n_obs, cs_real_t);
    for (int ii = 0; ii < n_obs; ii++)
      for (int kk = 0; kk < ms->dim; kk++)
        oi->time_weights[ms->dim*ii+kk] = -999.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute observation operator (H).
 *
 * \param[in]  ms  pointer to measures set
 * \param[in]  oi  pointer to an optimal interpolation
 * \param[in]  ig  pointer to interpol grid
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_obs_operator(cs_measures_set_t  *ms,
                              cs_at_opt_interp_t *oi,
                              cs_interpol_grid_t *ig)
{
  switch(oi->interp_type) {
  case CS_AT_OPT_INTERP_P0:
#if _OI_DEBUG_
    bft_printf("   *Computing P0 interpolator\n");
#endif
    _p0_projection(ms, oi, ig);
    break;
  case CS_AT_OPT_INTERP_P1:
#if _OI_DEBUG_
    bft_printf("   *Computing P1 interpolator\n");
#endif
    _p1_projection(ms, oi, ig);
    break;
  default:
    assert(0);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute $\tens{H}\tens{B}\transpose{\tens{H}}$.
 *
 * \param[in]  ms  pointer to measures set
 * \param[in]  oi  pointer to an optimal interpolation
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_project_model_covariance(cs_measures_set_t  *ms,
                                          cs_at_opt_interp_t *oi)
{
  int n_obs = ms->nb_measures;
  cs_real_t *proj = oi->model_to_obs_proj;
  cs_lnum_t *proj_idx = oi->model_to_obs_proj_idx;

  const int dim = ms->dim;
  const int stride = dim + 3; /* dimension of field + dimension of space */

  BFT_MALLOC(oi->b_proj, n_obs*n_obs*dim, cs_real_t);
  cs_real_t *b_proj = oi->b_proj;

  const cs_real_t ir_xy2 = cs_math_sq(oi->ir[0]);
  const cs_real_t ir_z2 = cs_math_sq(oi->ir[1]);

  for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
    for (cs_lnum_t jj = 0; jj < n_obs; jj++) {
      for (int pp = 0; pp < dim; pp++)
        b_proj[dim*(ii*n_obs + jj) + pp] = 0;

      for (cs_lnum_t kk = proj_idx[ii]; kk < proj_idx[ii+1]; kk++) {
        cs_real_t x1 = (proj + kk*stride)[dim  ];
        cs_real_t y1 = (proj + kk*stride)[dim+1];
        cs_real_t z1 = (proj + kk*stride)[dim+2];

        for (cs_lnum_t ll = proj_idx[jj]; ll < proj_idx[jj+1]; ll++) {
          cs_real_t x2 = (proj + ll*stride)[dim  ];
          cs_real_t y2 = (proj + ll*stride)[dim+1];
          cs_real_t z2 = (proj + ll*stride)[dim+2];

          cs_real_t influ = _b_matrix(x1, y1, z1, x2, y2, z2, ir_xy2, ir_z2);

          for (int pp = 0; pp < dim; pp++)
            b_proj[dim*(ii*n_obs + jj) + pp] += (proj + kk*stride)[pp] * (proj + ll*stride)[pp]
                                              * influ;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count active observations and compute time weights in case of
 *        unsteady.
 *
 * \param[in]  ms           pointer to measures set
 * \param[in]  oi           optimal interpolation for field variable
 * \param[in]  f_oia        analysis field of field variable
 * \param[in]  inverse      boolean, true if it necessary to recompute the
 *                          inverse of HB(H)
 * \param[in]  ao_idx       index of active observations
 *
 * \return  number of active observations for each measures component.
 */
/*----------------------------------------------------------------------------*/

int *
cs_at_opt_interp_get_active_obs(cs_measures_set_t   *ms,
                                cs_at_opt_interp_t  *oi,
                                cs_field_t          *f_oia,
                                bool               **inverse,
                                int               ***ao_idx)
{
  const cs_mesh_t *mesh = cs_glob_mesh;

  const int n_obs = ms->nb_measures;
  const int m_dim = ms->dim;

  int *n_active_obs = NULL;
  BFT_MALLOC(n_active_obs, m_dim, int);
  for (int kk = 0; kk < ms->dim; kk++)
    n_active_obs[kk] = 0;

  for (int kk = 0; kk < ms->dim; kk++)
    (*inverse)[kk] = true;

  /* steady mode */

  if (oi->steady > 0) {

    /* count active observations */

    for (int ii = 0; ii < n_obs; ii++)
      for (int kk = 0; kk < ms->dim; kk++)
        if (oi->measures_idx[m_dim*ii+kk] < oi->measures_idx[m_dim*(ii+1)+kk])
          n_active_obs[kk] += 1;

    /* build indirection from active observations to observations */

    BFT_MALLOC(*ao_idx, m_dim, int *);
    for (int kk = 0; kk < ms->dim; kk++)
      BFT_MALLOC((*ao_idx)[kk], n_active_obs[kk], int);

    int *ao_count = NULL;
    BFT_MALLOC(ao_count, m_dim, int);
    for (int kk = 0; kk < ms->dim; kk++)
      ao_count[kk] = 0;

    for (int ii = 0; ii < n_obs; ii++) {
      for (int kk = 0; kk < ms->dim; kk++) {
        if (oi->measures_idx[m_dim*ii+kk] < oi->measures_idx[m_dim*(ii+1)+kk]) {
          (*ao_idx)[kk][ao_count[kk]] = ii;
          ao_count[kk] += 1;
        }
      }
    }

    BFT_FREE(ao_count);

  /* unsteady mode - time management */

  } else {

    const cs_time_step_t *ts = cs_glob_time_step;

    BFT_MALLOC(*ao_idx, m_dim, int *);
    int *ao_count = NULL;
    BFT_MALLOC(ao_count, m_dim, int);
    for (int kk = 0; kk < ms->dim; kk++)
      ao_count[kk] = 0;

    /* storing previous time coefficients */
    cs_real_t *temp = NULL;
    BFT_MALLOC(temp, m_dim*n_obs, cs_real_t);
    for (int ii = 0; ii < m_dim*n_obs; ii++)
      temp[ii] = oi->time_weights[ii];

    for (int kk = 0; kk < ms->dim; kk++) {
      for (int ii = 0; ii < n_obs; ii++) {
        /* finding closest observation time */
        cs_real_t deltat = cs_math_infinite_r;
        for (int jj = oi->active_time[m_dim*ii+kk];
             jj < oi->measures_idx[m_dim*(ii+1)+kk];
             jj++) {
          cs_real_t dt = ts->t_cur - oi->times[jj];
          if (CS_ABS(dt) > CS_ABS(deltat))
            break;
          deltat = dt;
          oi->active_time[m_dim*ii+kk] = jj;
        }

        /* computing time coefficient */

        if (deltat < oi->time_window[0] || deltat > oi->time_window[3]) {
          oi->time_weights[m_dim*ii+kk] = 0.;
        } else if (   deltat >= oi->time_window[1]
                   && deltat <= oi->time_window[2]) {
          oi->time_weights[m_dim*ii+kk] = 1.;
        } else if (deltat < oi->time_window[1]) {
          oi->time_weights[m_dim*ii+kk] = (deltat - oi->time_window[0])
            / (oi->time_window[1] - oi->time_window[0]);
        } else if (deltat > oi->time_window[2]) {
          oi->time_weights[m_dim*ii+kk] = (oi->time_window[3] - deltat)
            / (oi->time_window[3] - oi->time_window[2]);
        }

        /* incrementing number of active observations */

        if (oi->time_weights[m_dim*ii+kk] > cs_math_epzero)
          n_active_obs[kk] += 1;
      }

      /* stop if no active obs */

      if (n_active_obs[kk] == 0) {
        int ll = ms->comp_ids[kk];
        for (cs_lnum_t c_id = 0; c_id < mesh->n_cells ; c_id++) {
          f_oia->val[f_oia->dim*c_id+ll] = 0.;
        }
#if _OI_DEBUG_
        bft_printf("   * No active observation\n\n");
#endif

      } else {

        /* build indirection from active observations to observations */

        BFT_MALLOC((*ao_idx)[kk], n_active_obs[kk], int);

        for (int ii = 0; ii < n_obs; ii++) {
          if (oi->time_weights[m_dim*ii+kk] > 1.e-300) { // FIXME too small?
            (*ao_idx)[kk][ao_count[kk]] = ii;
            ao_count[kk] += 1;
          }
        }

        /* L1 norm ||weights(n+1)-weights(n)|| */

        cs_real_t sum = 0.;
        for (cs_lnum_t ii = 0; ii < n_obs; ii++) {
          sum += CS_ABS(oi->time_weights[m_dim*ii+kk] - temp[m_dim*ii+kk]);
        }
        BFT_FREE(temp);
        (*inverse)[kk] = (sum > 1.e-6); // FIXME define tolerance

      } /* end if on n_active_obs */

    } /* end for on measures dimension */

    BFT_FREE(ao_count);

  } /* end if steady */

  return n_active_obs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute analysis for a given variable.
 *
 * \param[in]  f            field variable of which analysis will be computed
 * \param[in]  oi           optimal interpolation for field variable
 * \param[in]  f_oia        analysis field of field variable
 * \param[in]  n_active_obs number of active observations.
 * \param[in]  ao_idx       index of active observations
 * \param[in]  inverse      boolean, true if it necessary to recompute the
 *                          inverse of HB(H)
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_compute_analysis(cs_field_t         *f,
                                  cs_at_opt_interp_t *oi,
                                  cs_field_t         *f_oia,
                                  int                 n_active_obs,
                                  int                *ao_idx,
                                  bool                inverse,
                                  int                 mc_id)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_3_t  *restrict cell_cen =
    (const cs_real_3_t *restrict)fvq->cell_cen;

  /* measures set */
  const int key_ms = cs_field_key_id("measures_set_id");
  const int ims = cs_field_get_key_int(f, key_ms);
  cs_measures_set_t *ms = cs_measures_set_by_id(ims);

  cs_real_t *proj = oi->model_to_obs_proj;
  cs_lnum_t *proj_idx = oi->model_to_obs_proj_idx;
  cs_lnum_t *proj_c_ids = oi->model_to_obs_proj_c_ids;

  cs_interpol_grid_t *ig = cs_interpol_grid_by_id(oi->ig_id);

#if _OI_DEBUG_
  bft_printf("   * %i active observations\n    ", n_active_obs);
  for (int ii = 0; ii < n_active_obs; ii++)
    bft_printf("%i ", ao_idx[ii]);
  bft_printf("\n");
  if (oi->steady <= 0) {
    bft_printf("\n   * Time coefficients\n    ");
    for (int ii = 0; ii < n_active_obs; ii++)
      bft_printf("%.2f ", oi->time_weights[ao_idx[ii]]);
    bft_printf("\n");
  }
#endif

  /* building innovation/increment (Y0 - HX0) */

  int f_dim = f->dim;
  int m_dim = ms->dim;
  const int stride = m_dim + 3; /* dimension of field + dimension of space */

  int a_l_size = n_active_obs;
  int a_size = cs_math_sq(a_l_size);

  cs_real_t *inc = NULL;
  BFT_MALLOC(inc, a_l_size, cs_real_t);

  for (int ii = 0; ii < n_active_obs; ii++) {
    int obs_id = ao_idx[ii]; /* retrieve obs id */

    int r_id0 = 0;
    if (cs_glob_rank_id > -1) r_id0 = ig->rank_connect[obs_id];

    if (cs_glob_rank_id < 0 || cs_glob_rank_id == r_id0) {
      inc[ii] = ms->measures[oi->active_time[m_dim*obs_id+mc_id]];

      for (cs_lnum_t jj = proj_idx[obs_id]; jj < proj_idx[obs_id+1]; jj++) {
        cs_lnum_t c_id1 = proj_c_ids[jj];
        inc[ii] -= ((proj + stride*jj)[mc_id])
                  *f->val_pre[c_id1*f_dim + ms->comp_ids[mc_id]];
      }
    }
  }

  /* exchange innovation */

  if (cs_glob_rank_id > -1) {
    for (int ii = 0; ii < n_active_obs; ii++) {
      int obs_id = ao_idx[ii];
      int r_id0 = ig->rank_connect[obs_id];
      cs_parall_bcast(r_id0, 1, CS_DOUBLE, &(inc[ii]));
    }
  }

#if _OI_DEBUG_
  bft_printf("\n   * Observation increments\n    ");
  for (int ii = 0; ii < n_active_obs; ii++) {
    bft_printf("\n");
    bft_printf("%.2f ", inc[ii]);
  }
  bft_printf("\n");
#endif

  cs_real_t *alu = NULL;
  if (inverse) {
    cs_real_t *a = _assembly_adding_obs_covariance(ms,
                                                   oi,
                                                   ao_idx,
                                                   n_active_obs,
                                                   mc_id);
    BFT_MALLOC(alu, a_size, cs_real_t);

    cs_math_fact_lu(1, a_l_size, a, alu);
    BFT_FREE(a);

#if _OI_DEBUG_
    bft_printf("\n   * LU Matrix\n");
    for (int ii = 0; ii < n_active_obs; ii++) {
      bft_printf("    ");
      for (int jj = 0; jj < n_active_obs; jj++) {
        int id = ii*a_l_size + jj;
        bft_printf("%.8f ", alu[id]);
      }
      bft_printf("\n");
    }
#endif
  }

#if _OI_DEBUG_
  bft_printf("\n   * Computing (HBHT + R)^-1*I\n");
#endif

  cs_real_t *zeros = NULL;
  BFT_MALLOC(zeros, a_l_size, cs_real_t);
  for (cs_lnum_t ii = 0; ii < a_l_size; ii++)
    zeros[ii] = 0.;

  cs_real_t *vect = NULL;
  BFT_MALLOC(vect, a_l_size, cs_real_t);

  /* Forward and backward */

  cs_math_fw_and_bw_lu(alu, a_l_size, vect, zeros, inc);

  BFT_FREE(alu);
  BFT_FREE(zeros);
  BFT_FREE(inc);

  const cs_real_t ir_xy2 = cs_math_sq(oi->ir[0]);
  const cs_real_t ir_z2 = cs_math_sq(oi->ir[1]);

  for (cs_lnum_t ii = 0; ii < mesh->n_cells; ii++) {
    f_oia->val[ii*f_dim+ms->comp_ids[mc_id]] =
      f->val_pre[ii*f_dim+ms->comp_ids[mc_id]];

    for (int ll = 0; ll < n_active_obs; ll++) {
      for (int mm = proj_idx[ao_idx[ll]];
           mm < proj_idx[ao_idx[ll]+1];
           mm++) {
        cs_real_t x = (proj + mm*stride)[m_dim  ];
        cs_real_t y = (proj + mm*stride)[m_dim+1];
        cs_real_t z = (proj + mm*stride)[m_dim+2];
        f_oia->val[ii*f_dim+ms->comp_ids[mc_id]] +=  (proj + mm*stride)[mc_id] * vect[ll]
                                                   * _b_matrix(cell_cen[ii][0],
                                                               cell_cen[ii][1],
                                                               cell_cen[ii][2],
                                                               x, y, z, ir_xy2, ir_z2);
      }
    }
  }

  BFT_FREE(vect);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
