/*============================================================================
 * Log field and other array statistics at relevant time steps.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array_reduce.h"
#include "cs_base.h"
#include "cs_blas.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_function.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_notebook.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_range_set.h"
#include "cs_time_moment.h"
#include "cs_time_plot.h"
#include "cs_time_step.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_log.h"
#include "fvm_convert_array.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_log_iteration.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_log_iteration.c

  \brief Log field and other array statistics at relevant time steps.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/* Simple statistics */
/*-------------------*/

typedef struct {

  int     name_id;   /* Associated name id */
  int     cat_id;    /* Associated category id */
  int     loc_id;    /* Associated mesh location id */
  bool    intensive; /* Are associated values intensive ? */
  int     dim;       /* Associated dimension */
  int     v_idx;     /* Start index of values */

} cs_log_sstats_t;

/* Clipping info */
/*---------------*/

typedef struct {

  int     f_id;     /* Associated field id, or -1 */
  int     name_id;  /* Associated name id if not a field, -1 for a field */
  int     dim;      /* Associated dimension */
  int     c_idx;    /* Start index of counts */
  int     v_idx;    /* Start index of values */

} cs_log_clip_t;

/*----------------------------------------------------------------------------
 * Prototypes for local static functions
 *----------------------------------------------------------------------------*/

static bool
_log_time_control_automatic(const cs_time_step_t  *ts,
                          void                  *input);

static bool
_log_time_control_interval(const cs_time_step_t  *ts,
                           void                  *input);

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_map_name_to_id_t  *_name_map = nullptr;

static cs_map_name_to_id_t  *_category_map = nullptr;
static int _sstats_val_size = 0;
static int _sstats_val_size_max = 0;
static int _n_sstats = 0;
static int _n_sstats_max = 0;
static double *_sstats_vmin = nullptr;
static double *_sstats_vmax = nullptr;
static double *_sstats_vsum = nullptr;
static double *_sstats_wsum = nullptr;
static cs_log_sstats_t  *_sstats = nullptr;

static int _clips_val_size = 0;
static int _clips_val_size_max = 0;
static int _n_clips = 0;
static int _n_clips_max = 0;
static cs_gnum_t *_clips_count = nullptr;
static double *_clips_vmin = nullptr;
static double *_clips_vmax = nullptr;
static cs_log_clip_t  *_clips = nullptr;

static cs_time_plot_t  *_l2_residual_plot = nullptr;

static int _log_interval_base = 10;

static cs_time_control_t  _log_time_control
  = {.type = CS_TIME_CONTROL_FUNCTION,
     .at_start = true,
     .at_first = true,
     .at_end = true,
     .start_nt = 0,
     .end_nt = 0,
     .interval_nt = 0,
     .control_func = _log_time_control_automatic,
     .control_input = nullptr,
     .current_state = true,
     .current_time_step = -1,
     .last_nt = -1,
     .last_t = -1};

cs_time_control_t  *cs_glob_log_iteration_time_control = &_log_time_control;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Default function for main log activation.
 *
 * The activation occurs:
 * - At each time step for the first n absolute or restarted time steps.
 * - Every 5 time steps between n and 5.n time steps.
 * - Every 10 time steps between 5.n and 10.n time steps.
 * - Every 50 time steps between 10.n and 50.n time steps.
 * - Every 100 time steps between 50.n and 100.n time steps.
 * - And so on...
 *
 * \param[in]  ts     current time step structure
 * \param[in]  input  pointer to integer n, or nullptr
 *
 * \return  true if log is active at given time step, false otherwise
 */
/*----------------------------------------------------------------------------*/

static bool
_log_time_control_automatic(const cs_time_step_t  *ts,
                            void                  *input)
{
  int n = (input != nullptr) ? *((int *)input) : 10;

  bool active = false;

  int nt_abs = ts->nt_cur;
  int nt_rel = ts->nt_cur - ts->nt_prev;

  if (nt_abs < n || nt_rel < n)
    active = true;

  else if (ts->nt_max > -1 && nt_abs >= ts->nt_max)
    active = true;

  else if (ts->t_max > 0 && ts->t_cur >= ts->t_max)
    active = true;

  else {
    int n10 = 10*n, n5 = 5*n;
    while (nt_abs > n10 && nt_abs%n5 == 0) {
      nt_abs /= n;
    }
    if (nt_abs <= n10) {
      int interval = 10;
      if (nt_abs <= n5)
        interval = 5;
      if (nt_abs%interval == 0)
        active = true;
    }
  }

  return active;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Time interval function for main log activation.
 *
 * The activation occurs:
 * - Each of the first 10 time steps
 * - At the end of the computation
 * - Every "nt_interval" time step
 *
 * \param[in]  ts     current time step structure
 * \param[in]  input  pointer to integer nt_interval, or nullptr
 *
 * \return  true if log is active at given time step, false otherwise
 */
/*----------------------------------------------------------------------------*/

static bool
_log_time_control_interval(const cs_time_step_t  *ts,
                           void                  *input)
{
  int nt_interval = (input != nullptr) ? *((int *)input) : 10;

  bool active = false;

  int nt_abs = ts->nt_cur;
  int nt_rel = ts->nt_cur - ts->nt_prev;

  if (ts->nt_max > -1 && nt_abs >= ts->nt_max)
    active = true;

  else if (ts->t_max > 0 && ts->t_cur >= ts->t_max)
    active = true;

  else if (nt_interval > 0) {
    if (nt_abs <= 10 || nt_rel <= 10)
      active = true;
    else if (nt_abs%nt_interval == 0)
      active = true;
  }

  return active;
}

/*----------------------------------------------------------------------------
 * Compare simple stats elements (qsort function).
 *
 * parameters:
 *   x <-> pointer to first element
 *   y <-> pointer to second element
 *
 * returns:
 *   -1 if x < y, 0 if x = y, or 1 if x > y
 *----------------------------------------------------------------------------*/

static int
_compare_sstats(const void *x,
                const void *y)
{
  int retval = 1;

  const cs_log_sstats_t *s0 = static_cast<const cs_log_sstats_t *>(x);
  const cs_log_sstats_t *s1 = static_cast<const cs_log_sstats_t *>(y);

  if (s0->cat_id < s1->cat_id)
    retval = -1;

  else if (s0->cat_id == s1->cat_id) {
    if (s0->name_id < s1->name_id)
      retval = -1;
    else if (s0->name_id == s1->name_id)
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Find simple stats by id
 *
 * parameters:
 *   cat_id     <-- category id
 *   name_id    <-- name id
 *
 * returns:
 *   id of simple stats, or -1 if not found
 *----------------------------------------------------------------------------*/

static int
_find_sstats(int  cat_id,
             int  name_id)
{
  /* Use binary search to find entry */

  int sstat_id = -1;

  if (_n_sstats > 0) {

    int start_id = 0;
    int end_id = _n_sstats;

    while (start_id < end_id) {
      int mid_id = start_id + ((end_id -start_id) / 2);
      int cmp_ret = 0;
      int cmp_cat = _sstats[mid_id].cat_id;
      if (cmp_cat < cat_id)
        cmp_ret = -1;
      else if (cmp_cat > cat_id)
        cmp_ret = 1;
      else {
        int cmp_name = _sstats[mid_id].name_id;
        if (cmp_name < name_id)
          cmp_ret = -1;
        else if (cmp_name > name_id)
          cmp_ret = 1;
        else {
          sstat_id = mid_id;
          break;
        }
      }
      if (cmp_ret < 0)
        start_id = mid_id + 1;
      else if (cmp_ret > 0)
        end_id = mid_id;
    }
  }

  return sstat_id;
}

/*----------------------------------------------------------------------------
 * Compare clipping elements (qsort function).
 *
 * parameters:
 *   x <-> pointer to first element
 *   y <-> pointer to second element
 *
 * returns:
 *   -1 if x < y, 0 if x = y, or 1 if x > y
 *----------------------------------------------------------------------------*/

static int _compare_clips(const void *x, const void *y)
{
  int retval = 1;

  const cs_log_clip_t *c0 = static_cast<const cs_log_clip_t *>(x);
  const cs_log_clip_t *c1 = static_cast<const cs_log_clip_t *>(y);

  if (c0->name_id < c1->name_id)
    retval = -1;

  else if (c0->name_id == c1->name_id) {
    if (c0->f_id < c1->f_id)
      retval = -1;
    else if (c0->f_id == c1->f_id)
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Find clippings by id
 *
 * parameters:
 *   f_id    <-- field id
 *   name_id <-- name id
 *
 * returns:
 *   id of clipping, or -1 if not found
 *----------------------------------------------------------------------------*/

static int
_find_clip(int  f_id,
           int  name_id)
{
  /* Use binary search to find entry */

  int clip_id = -1;

  if (_n_clips > 0) {

    int start_id = 0;
    int end_id = _n_clips;

    while (start_id < end_id) {
      int mid_id = start_id + ((end_id -start_id) / 2);
      int cmp_ret = 0;
      int cmp_name = _clips[mid_id].name_id;
      if (cmp_name < name_id)
        cmp_ret = -1;
      else if (cmp_name > name_id)
        cmp_ret = 1;
      else {
        int cmp_f = _clips[mid_id].f_id;
        if (cmp_f < f_id)
          cmp_ret = -1;
        else if (cmp_f > f_id)
          cmp_ret = 1;
        else {
          clip_id = mid_id;
          break;
        }
      }
      if (cmp_ret < 0)
        start_id = mid_id + 1;
      else if (cmp_ret > 0)
        end_id = mid_id;
    }
  }

  return clip_id;
}

/*----------------------------------------------------------------------------
 * Log information for a given array.
 *
 * parameters:
 *   prefix       <-- string inserted before name
 *   name         <-- array name or label
 *   name_width   <-- width of "name" column
 *   dim          <-- array dimension
 *   n_g_elts     <-- global number of associated elements,
 *   total_weight <-- if > 0, weight (volume or surface) of array location
 *   vmin         <-- minimum values of each component or norm
 *   vmax         <-- maximum values of each component or norm
 *   vsum         <-- sum of each component or norm
 *   wsum         <-- weighted sum of each component or norm, or nullptr
 *   fpe_flag     <-- was a "not a number" or floating-point error detected ?
 *----------------------------------------------------------------------------*/

static void
_log_array_info(const char        *prefix,
                const char        *name,
                size_t             name_width,
                int                dim,
                cs_gnum_t          n_g_elts,
                double             total_weight,
                double             vmin[],
                const double       vmax[],
                const double       vsum[],
                const double      *wsum,
                int               *fpe_flag)
{
  const int _dim = (dim == 3) ? 4 : dim;

  char tmp_s[2][64] =  {"", ""};

  for (int c_id = 0; c_id < _dim; c_id++) {

    if (dim == 3) {
      if (c_id < 3)
        snprintf(tmp_s[1], 63, "%s%s", name, cs_glob_field_comp_name_3[c_id]);
      else
        snprintf(tmp_s[1], 63, "ǁ%sǁ", name);
      tmp_s[1][63] = '\0';
      cs_log_strpad(tmp_s[0], tmp_s[1], name_width, 64);
    }
    else if (dim == 6) {
      snprintf(tmp_s[1], 63, "%s%s", name, cs_glob_field_comp_name_6[c_id]);
      tmp_s[1][63] = '\0';
      cs_log_strpad(tmp_s[0], tmp_s[1], name_width, 64);
    }
    else if (dim == 9) {
      snprintf(tmp_s[1], 63, "%s%s", name, cs_glob_field_comp_name_9[c_id]);
      tmp_s[1][63] = '\0';
      cs_log_strpad(tmp_s[0], tmp_s[1], name_width, 64);
    }
    else
      cs_log_strpad(tmp_s[0], name, name_width, 64);

    if (total_weight >= 0)
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s%s  %14.5g  %14.5g  %14.5g  %14.5g\n",
                    prefix,
                    tmp_s[0],
                    vmin[c_id],
                    vmax[c_id],
                    vsum[c_id] / n_g_elts,
                    wsum[c_id] / total_weight);
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s%s  %14.5g  %14.5g  %14.5g\n",
                    prefix,
                    tmp_s[0],
                    vmin[c_id],
                    vmax[c_id],
                    vsum[c_id] / n_g_elts);

    /* Check NAN  and exit */
    if (isnan(vsum[c_id]))
      *fpe_flag = 1;
  }

}

/*----------------------------------------------------------------------------
 * Log information for a given clipping.
 *
 * parameters:
 *   prefix       <-- string inserted before name
 *   name         <-- array name or label
 *   name_width   <-- width of "name" column
 *   dim          <-- array dimension
 *   count_min    <-- global number of clips to minimum
 *   count_max    <-- global number of clips to maximum
 *   vmin         <-- minimum values of each component
 *   vmax         <-- maximum values of each component
 *----------------------------------------------------------------------------*/

static void
_log_clip_info(const char        *prefix,
               const char        *name,
               size_t             name_width,
               int                dim,
               cs_gnum_t          count_min[],
               cs_gnum_t          count_max[],
               const double       vmin[],
               const double       vmax[])
{
  int c_id;
  const int _dim = (dim == 1) ? dim : dim + 1;

  char tmp_s[2][64] =  {"", ""};

  for (c_id = 0; c_id < _dim; c_id++) {

    int    _count_max, _count_min;
    double _vmin, _vmax;
    const char *_name = (c_id == 0) ? name : " ";

    if (c_id < _dim) {
      _count_min = count_min[2*c_id];
      _count_max = count_max[2*c_id];
      _vmin = vmin[c_id];
      _vmax = vmax[c_id];
    }

    if (dim > 1) {
      if (c_id == 0) {
        snprintf(tmp_s[1], 63, "%s", _name);
        tmp_s[1][63] = '\0';
        cs_log_strpad(tmp_s[0], tmp_s[1], name_width, 64);
        _vmin = 0.;
        _vmax = 0.;
        for (int i = 0; i< dim; i++) {
          _vmin += vmin[i]*vmin[i];
          _vmax += vmax[i]*vmax[i];
        }
        _vmin = sqrt(_vmin);
        _vmax = sqrt(_vmax);
      }
      else {
        if (dim == 3)
          snprintf(tmp_s[1], 63, "%s%s", name, cs_glob_field_comp_name_3[c_id - 1]);
        else if (dim == 6)
          snprintf(tmp_s[1], 63, "%s%s", name, cs_glob_field_comp_name_6[c_id - 1]);
        tmp_s[1][63] = '\0';
        cs_log_strpad(tmp_s[0], tmp_s[1], name_width, 64);
      }
    }
    else
      cs_log_strpad(tmp_s[0], name, name_width, 64);

    if (_count_min > 0 && _count_max > 0)
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s%s  %14.5g  %14.5g  %12llu  %12llu\n",
                    prefix,
                    tmp_s[0],
                    _vmin,
                    _vmax,
                    (unsigned long long)_count_min,
                    (unsigned long long)_count_max);
    else if (_count_min > 0)
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s%s  %14.5g                  %12llu\n",
                    prefix,
                    tmp_s[0],
                    _vmin,
                    (unsigned long long)_count_min);
    else if (_count_max > 0)
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s%s                  %14.5g                %12llu\n",
                    prefix,
                    tmp_s[0],
                    _vmax,
                    (unsigned long long)_count_max);
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s%s\n",
                    prefix,
                    tmp_s[0]);
  }
}

/*----------------------------------------------------------------------------
 * Main logging output of variables.
 *----------------------------------------------------------------------------*/

static void
_log_fields_and_functions(void)
{
  int log_count;

  int log_count_max = 0;
  int fpe_flag = 0;
  int *log_id = nullptr, *moment_id = nullptr, *f_location_id = nullptr;
  bool *location_log = nullptr;
  double  *vmin = nullptr, *vmax = nullptr, *vsum = nullptr, *wsum = nullptr;

  char tmp_s[5][64] =  {"", "", "", "", ""};

  const char _underline[] = "---------------------------------";
  const int n_locations = cs_mesh_location_n_locations();
  const int n_fields = cs_field_n_fields();
  const int n_functions = cs_function_n_functions();
  const int n_moments = cs_time_moment_n_moments();
  const int log_key_id = cs_field_key_id("log");
  const int label_key_id = cs_field_key_id("label");

  const int n_ff = n_fields + n_functions;

  cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Allocate working arrays */

  log_count_max = n_fields + n_functions;

  BFT_MALLOC(log_id, n_ff, int);
  BFT_MALLOC(vmin, log_count_max, double);
  BFT_MALLOC(vmax, log_count_max, double);
  BFT_MALLOC(vsum, log_count_max, double);
  BFT_MALLOC(wsum, log_count_max, double);

  BFT_MALLOC(location_log, n_locations, bool);
  for (int i = 0; i < n_locations; i++)
    location_log[i] = false;

  BFT_MALLOC(f_location_id, n_ff, int);
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, log_key_id)) {
      f_location_id[f_id] = f->location_id;
      location_log[f->location_id] = true;
    }
    else
      f_location_id[f_id] = -1;
  }
  for (int f_id = 0; f_id < n_functions; f_id++) {
    const cs_function_t  *f = cs_function_by_id(f_id);
    if (f->log) {
      f_location_id[n_fields + f_id] = f->location_id;
      location_log[f->location_id] = true;
    }
    else
      f_location_id[n_fields + f_id] = -1;
  }

  if (n_moments > 0) {
    BFT_MALLOC(moment_id, n_fields, int);
    for (int f_id = 0; f_id < n_fields; f_id++)
      moment_id[f_id] = -1;
    for (int m_id = 0; m_id < n_moments; m_id++) {
      const cs_field_t *f = cs_time_moment_get_field(m_id);
      if (f != nullptr) {
        moment_id[f->id] = m_id;

        /* Only log active moments */
        if (!cs_time_moment_is_active(m_id))
          f_location_id[f->id] = -1;
      }
    }
  }

  /* Loop on locations */

  for (int loc_id = 0; loc_id < n_locations; loc_id++) {

    if (location_log[loc_id] == false)
      continue;

    size_t max_name_width = cs_log_strlen(_("field"));
    cs_lnum_t have_weight = 0;
    double total_weight = -1;
    cs_gnum_t n_g_elts = 0;
    cs_real_t *gather_array = nullptr; /* only if CS_MESH_LOCATION_VERTICES */
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(loc_id);
    const cs_lnum_t _n_elts = n_elts[0];
    const cs_lnum_t *elt_ids = nullptr;
    const cs_real_t *weight = nullptr;

    if (mq != nullptr) {

      switch(loc_id) {

      case CS_MESH_LOCATION_CELLS:
        n_g_elts = m->n_g_cells;
        weight = mq->cell_vol;
        have_weight = 1;
        total_weight = mq->tot_vol;
        break;

      case CS_MESH_LOCATION_INTERIOR_FACES:
        n_g_elts = m->n_g_i_faces;
        weight = mq->i_face_surf;
        cs_array_reduce_sum_l(_n_elts, 1, nullptr, weight, &total_weight);
        cs_parall_sum(1, CS_DOUBLE, &total_weight);
        have_weight = 1;
        break;

      case CS_MESH_LOCATION_BOUNDARY_FACES:
        n_g_elts = m->n_g_b_faces;
        weight = mq->b_face_surf;
        cs_array_reduce_sum_l(_n_elts, 1, nullptr, weight, &total_weight);
        cs_parall_sum(1, CS_DOUBLE, &total_weight);
        have_weight = 1;
        break;

      case CS_MESH_LOCATION_VERTICES:
        n_g_elts = m->n_g_vertices;
        have_weight = 0;
        BFT_MALLOC(gather_array, m->n_vertices, cs_real_t);
        break;

      default:
        {
          cs_mesh_location_type_t
            loc_type = cs_mesh_location_get_type(loc_id);

          elt_ids = cs_mesh_location_get_elt_ids_try(loc_id);

          /* FIXME: using sum is correct for cells and boundary faces,
           *        would need range set for interior faces and vertices. */
          n_g_elts = _n_elts;
          cs_parall_counter(&n_g_elts, 1);

          switch(loc_type) {
          case CS_MESH_LOCATION_CELLS:
            weight = mq->cell_vol;
            have_weight = 1;
            break;
          case CS_MESH_LOCATION_INTERIOR_FACES:
            weight = mq->i_face_surf;
            have_weight = 1;
            break;
          case CS_MESH_LOCATION_BOUNDARY_FACES:
            weight = mq->b_face_surf;
            have_weight = 1;
            break;
          default:
            break;
          }

          if (have_weight) {
            cs_array_reduce_sum_l(_n_elts, 1, elt_ids, weight, &total_weight);
            cs_parall_sum(1, CS_DOUBLE, &total_weight);
          }
        }
        break;
      }
    }

    if (n_g_elts == 0)
      continue;

    log_count = 0;

    /* First loop on fields and functions */

    for (int f_id = 0; f_id < n_ff; f_id++) {

      if (f_location_id[f_id] != loc_id) {
        log_id[f_id] = -1;
        continue;
      }

      bool use_weight = false;

      const char *name;
      cs_lnum_t f_dim, _dim, c_id;
      cs_real_t *f_val, *_f_val;

      if (f_id < n_fields) { /* Field */
        const cs_field_t  *f = cs_field_by_id(f_id);

        name = cs_field_get_key_str(f, label_key_id);
        if (name == nullptr)
          name = f->name;

        f_dim = f->dim;

        if (have_weight && (f->type & CS_FIELD_INTENSIVE))
          use_weight = true;

        _f_val = nullptr;
        f_val = f->val;
      }
      else {  /* Function */
        const cs_function_t  *f = cs_function_by_id(f_id - n_fields);

        name = f->label;
        if (name == nullptr)
          name = f->name;

        f_dim = f->dim;

        if (have_weight && (f->type & CS_FUNCTION_INTENSIVE))
          use_weight = true;

        BFT_MALLOC(_f_val, f_dim*_n_elts, cs_real_t);
        f_val = _f_val;

        const cs_time_step_t *ts = cs_glob_time_step;

        if (f->datatype == CS_REAL_TYPE) {
          cs_function_evaluate(f,
                               ts,
                               f->location_id,
                               _n_elts,
                               nullptr,  /* elt_ids */
                               _f_val);
        }
        else {
          const cs_lnum_t  parent_id_shift[] = {0};
          char *buffer;
          const void *src_data[1];
          BFT_MALLOC(buffer, f_dim*_n_elts*cs_datatype_size[f->datatype], char);
          src_data[0] = buffer;
          cs_function_evaluate(f,
                               ts,
                               f->location_id,
                               _n_elts,
                               nullptr,  /* elt_ids */
                               (void *)buffer);
          fvm_convert_array(f->dim, 0, f->dim,
                            0, _n_elts,
                            CS_INTERLACE,
                            f->datatype,
                            CS_REAL_TYPE,
                            0, /* n_parent_lists */
                            parent_id_shift,
                            nullptr, /* parent_id */
                            src_data,
                            _f_val);
          BFT_FREE(buffer);
        }
      }

      size_t l_name_width = cs_log_strlen(name);
      if (f_dim == 3)
        l_name_width += 3;
      else if (f_dim > 3)
        l_name_width += 4;

      max_name_width = CS_MAX(max_name_width, l_name_width);

      /* Position in log */

      log_id[f_id] = log_count;

      _dim = (f_dim == 3) ? 4 : f_dim;

      while (log_count + _dim > log_count_max) {
        log_count_max *= 2;
        BFT_REALLOC(vmin, log_count_max, double);
        BFT_REALLOC(vmax, log_count_max, double);
        BFT_REALLOC(vsum, log_count_max, double);
        BFT_REALLOC(wsum, log_count_max, double);
      }

      if (use_weight) {
        cs_array_reduce_simple_stats_l_w(_n_elts,
                                         f_dim,
                                         nullptr,
                                         elt_ids,
                                         f_val,
                                         weight,
                                         vmin + log_count,
                                         vmax + log_count,
                                         vsum + log_count,
                                         wsum + log_count);

      }
      else {

        cs_real_t  *field_val = f_val;
        cs_lnum_t  _n_cur_elts = _n_elts;

        /* Eliminate shared values whose local rank is not owner and compact */
        if (loc_id == CS_MESH_LOCATION_VERTICES) {
          if (m->vtx_range_set == nullptr)
            m->vtx_range_set = cs_range_set_create(m->vtx_interfaces,
                                                   nullptr,
                                                   m->n_vertices,
                                                   false, /* balance */
                                                   2,  /* tr_ignore */
                                                   0); /* g_id_base */

          if (f_dim > 1)
            BFT_REALLOC(gather_array, (f_dim * m->n_vertices), cs_real_t);

          cs_range_set_gather(m->vtx_range_set,
                              CS_REAL_TYPE,
                              f_dim,
                              f_val,
                              gather_array);
          field_val = gather_array;
          _n_cur_elts = m->vtx_range_set->n_elts[0];
        }

        cs_array_reduce_simple_stats_l(_n_cur_elts,
                                       f_dim,
                                       nullptr,
                                       field_val,
                                       vmin + log_count,
                                       vmax + log_count,
                                       vsum + log_count);

        if (have_weight) {
          for (c_id = 0; c_id < _dim; c_id++)
            wsum[log_count + c_id] = 0.;
        }
      }

      BFT_FREE(_f_val);

      log_count += _dim;

    } /* End of first loop on fields */

    if (gather_array != nullptr)
      BFT_FREE(gather_array);

    if (log_count < 1)
      continue;

    /* Group MPI operations if required */

    cs_parall_min(log_count, CS_DOUBLE, vmin);
    cs_parall_max(log_count, CS_DOUBLE, vmax);
    cs_parall_sum(log_count, CS_DOUBLE, vsum);

#   if defined(HAVE_MPI)
    cs_parall_counter_max(&have_weight, 1);
    if (have_weight)
      cs_parall_sum(log_count, CS_DOUBLE, wsum);
#   endif

    /* Print headers */

    max_name_width = CS_MIN(max_name_width, 63);

    const char *loc_name = _(cs_mesh_location_get_name(loc_id));
    size_t loc_name_w = cs_log_strlen(loc_name);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("\n"
                    "  ** Field values on %s\n"
                    "     ----------------%.*s\n"),
                  loc_name, (int)loc_name_w, _underline);

    cs_log_strpad(tmp_s[0], _("field"), max_name_width, 64);
    cs_log_strpadl(tmp_s[1], _("minimum"), 14, 64);
    cs_log_strpadl(tmp_s[2], _("maximum"), 14, 64);
    cs_log_strpadl(tmp_s[3], _("set mean"), 14, 64);
    if (have_weight) {
      cs_log_strpadl(tmp_s[4], _("spatial mean"), 14, 64);
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n   %s  %s  %s  %s  %s\n",
                    tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);
    }
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n   %s  %s  %s  %s\n",
                    tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

    /* Underline */

    size_t n_cols = (have_weight) ? 5 : 4;

    for (size_t col = 0; col < n_cols; col++) {
      size_t i;
      size_t w0 = (col == 0) ? max_name_width : 14;
      for (i = 0; i < w0; i++)
        tmp_s[col][i] = '-';
      tmp_s[col][w0] = '\0';
    }
    if (have_weight) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "-  %s  %s  %s  %s  %s\n",
                    tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);
    }
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "-  %s  %s  %s  %s\n",
                    tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

    /* Second loop on fields and functions */

    log_count = 0;

    for (int f_id = 0; f_id < n_ff; f_id++) {

      if (log_id[f_id] < 0)
        continue;

      const char *name;
      int f_dim;
      double t_weight = -1;

      char prefix[] = "v  ";

      if (f_id < n_fields) { /* Field */
        const cs_field_t  *f = cs_field_by_id(f_id);
        name = cs_field_get_key_str(f, label_key_id);
        if (name == nullptr)
          name = f->name;
        f_dim = f->dim;
        if (total_weight > 0 && (f->type & CS_FIELD_INTENSIVE))
          t_weight = total_weight;
        if (moment_id != nullptr) {
          if (moment_id[f_id] > -1)
            prefix[0] = 'm';
        }
        if (f->type & CS_FIELD_ACCUMULATOR)
          prefix[0] = 'm';
      }
      else {  /* Function */
        const cs_function_t  *f = cs_function_by_id(f_id - n_fields);
        name = f->label;
        if (name == nullptr)
          name = f->name;
        f_dim = f->dim;
        if (total_weight > 0 && (f->type & CS_FUNCTION_INTENSIVE))
          t_weight = total_weight;
        prefix[0] = 'f';
      }

      int _dim = (f_dim == 3) ? 4 : f_dim;

      /* Position in log */

      _log_array_info(prefix,
                      name,
                      max_name_width,
                      f_dim,
                      n_g_elts,
                      t_weight,
                      vmin + log_count,
                      vmax + log_count,
                      vsum + log_count,
                      wsum + log_count,
                      &fpe_flag);

      log_count += _dim;

    } /* End of loop on fields */

  } /* End of loop on mesh locations */

  /* Check NaN and exit */
  if (fpe_flag == 1)
    bft_error(__FILE__, __LINE__, 0,
                _("Invalid (not-a-number) values detected for a field."));

  BFT_FREE(moment_id);
  BFT_FREE(wsum);
  BFT_FREE(vsum);
  BFT_FREE(vmax);
  BFT_FREE(vmin);
  BFT_FREE(log_id);
  BFT_FREE(f_location_id);
  BFT_FREE(location_log);

  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*----------------------------------------------------------------------------
 * Main logging output of additional simple statistics
 *----------------------------------------------------------------------------*/

static void
_log_sstats(void)
{
  int     stat_id;
  int     fpe_flag = 0;
  double _boundary_surf = -1;
  double _interior_surf = -1;
  double  *vmin = nullptr, *vmax = nullptr, *vsum = nullptr, *wsum = nullptr;

  char tmp_s[5][64] =  {"", "", "", "", ""};

  const char _underline[] = "---------------------------------";

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Allocate working arrays */

  BFT_MALLOC(vmin, _sstats_val_size, double);
  BFT_MALLOC(vmax, _sstats_val_size, double);
  BFT_MALLOC(vsum, _sstats_val_size, double);
  BFT_MALLOC(wsum, _sstats_val_size, double);

  memcpy(vmin, _sstats_vmin, _sstats_val_size*sizeof(double));
  memcpy(vmax, _sstats_vmax, _sstats_val_size*sizeof(double));
  memcpy(vsum, _sstats_vsum, _sstats_val_size*sizeof(double));
  memcpy(wsum, _sstats_wsum, _sstats_val_size*sizeof(double));

  /* Group MPI operations if required */

  cs_parall_min(_sstats_val_size, CS_DOUBLE, vmin);
  cs_parall_max(_sstats_val_size, CS_DOUBLE, vmax);
  cs_parall_sum(_sstats_val_size, CS_DOUBLE, vsum);
  cs_parall_sum(_sstats_val_size, CS_DOUBLE, wsum);

  /* Loop on statistics */

  int sstat_cat_start = 0;

  while (sstat_cat_start < _n_sstats) {

    int sstat_cat_end = sstat_cat_start;

    while (sstat_cat_end < _n_sstats) {
      if (_sstats[sstat_cat_end].cat_id != _sstats[sstat_cat_start].cat_id)
        break;
      else
        sstat_cat_end ++;
    }

    const char *cat_name
      = cs_map_name_to_id_reverse(_category_map,
                                  _sstats[sstat_cat_start].cat_id);

    const size_t cat_name_w = cs_log_strlen(_(cat_name));
    size_t max_name_width = cat_name_w;

    /* Now separate by mesh location */

    int loc_min = cs_mesh_location_n_locations() + 1;
    int loc_max = -1;

    for (stat_id = sstat_cat_start; stat_id < sstat_cat_end; stat_id++) {
      int loc_id = _sstats[stat_id].loc_id;
      loc_min = CS_MIN(loc_min, loc_id);
      loc_max = CS_MAX(loc_max, loc_id);
    }

    for (int loc_id = loc_min; loc_id <= loc_max; loc_id++) {

      int n_loc_stats = 0;
      for (stat_id = sstat_cat_start; stat_id < sstat_cat_end; stat_id++) {
        if (_sstats[stat_id].loc_id == loc_id)
          n_loc_stats++;
      }
      if (n_loc_stats == 0)
        continue;

      cs_gnum_t n_g_elts = 0;
      int have_weight = 0;
      double total_weight = -1;
      const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(loc_id);;
      const cs_lnum_t _n_elts = n_elts[0];
      const cs_real_t *weight = nullptr;
      const char *loc_name = _(cs_mesh_location_get_name(loc_id));
      size_t loc_name_w = cs_log_strlen(loc_name);

      if (mq != nullptr) {
        switch(loc_id) {
        case CS_MESH_LOCATION_CELLS:
          n_g_elts = m->n_g_cells;
          weight = mq->cell_vol;
          have_weight = 1;
          total_weight = mq->tot_vol;
          break;
        case CS_MESH_LOCATION_INTERIOR_FACES:
          n_g_elts = m->n_g_i_faces;
          weight = mq->i_face_surf;
          if (_interior_surf < 0) {
            cs_array_reduce_sum_l(_n_elts, 1, nullptr, weight, &_interior_surf);
            cs_parall_sum(1, CS_DOUBLE, &_interior_surf);
            if (_interior_surf < 0) _interior_surf = 0; /* just to be safe */
          }
          total_weight = _interior_surf;
          have_weight = 1;
          break;
        case CS_MESH_LOCATION_BOUNDARY_FACES:
          n_g_elts = m->n_g_b_faces;
          weight = mq->b_face_surf;
          if (_boundary_surf < 0) {
            cs_array_reduce_sum_l(_n_elts, 1, nullptr, weight, &_boundary_surf);
            cs_parall_sum(1, CS_DOUBLE, &_boundary_surf);
            if (_boundary_surf < 0) _boundary_surf = 0; /* just to be safe */
          }
          total_weight = _boundary_surf;
          have_weight = 1;
          break;
        case CS_MESH_LOCATION_VERTICES:
          n_g_elts = m->n_g_vertices;
          have_weight = 0;
          break;
        default:
          n_g_elts = _n_elts;
          cs_parall_counter(&n_g_elts, 1);
          break;
        }
      }

      for (stat_id = sstat_cat_start; stat_id < sstat_cat_end; stat_id++) {
        if (_sstats[stat_id].loc_id == loc_id) {
          const char *stat_name
            = cs_map_name_to_id_reverse(_name_map, _sstats[stat_id].name_id);
          size_t name_width = strlen(stat_name);
          max_name_width = CS_MAX(max_name_width, name_width);
        }
      }
      max_name_width = CS_MIN(max_name_width, 63);

      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n"
                      "  ** Computed values on %s\n"
                      "     -------------------%.*s\n"),
                    loc_name, (int)loc_name_w, _underline);

      cs_log_strpad(tmp_s[0], _(cat_name), max_name_width, 64);
      cs_log_strpadl(tmp_s[1], _("minimum"), 14, 64);
      cs_log_strpadl(tmp_s[2], _("maximum"), 14, 64);
      cs_log_strpadl(tmp_s[3], _("set mean"), 14, 64);
      if (have_weight) {
        cs_log_strpadl(tmp_s[4], _("spatial mean"), 14, 64);
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n   %s  %s  %s  %s  %s\n",
                      tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);
      }
      else
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n   %s  %s  %s  %s\n",
                      tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

      /* Underline */

      size_t n_cols = (have_weight) ? 5 : 4;

      for (size_t col = 0; col < n_cols; col++) {
        size_t i;
        size_t w0 = (col == 0) ? max_name_width : 14;
        for (i = 0; i < w0; i++)
          tmp_s[col][i] = '-';
        tmp_s[col][w0] = '\0';
      }
      if (have_weight) {
        cs_log_printf(CS_LOG_DEFAULT,
                      "   %s  %s  %s  %s  %s\n",
                      tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);
      }
      else
        cs_log_printf(CS_LOG_DEFAULT,
                      "   %s  %s  %s  %s\n",
                      tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

      /* Print values */

      for (stat_id = sstat_cat_start; stat_id < sstat_cat_end; stat_id++) {

        if (_sstats[stat_id].loc_id != loc_id)
          continue;

        const char *name
          = cs_map_name_to_id_reverse(_name_map, _sstats[stat_id].name_id);

        double t_weight = -1;
        if (total_weight > 0 && _sstats[stat_id].intensive)
          t_weight = total_weight;

        _log_array_info("   ",
                        name,
                        max_name_width,
                        _sstats[stat_id].dim,
                        n_g_elts,
                        t_weight,
                        vmin + stat_id,
                        vmax + stat_id,
                        vsum + stat_id,
                        wsum + stat_id,
                        &fpe_flag);

      } /* End of loop on stats */

    } /* End of loop on mesh locations */

    sstat_cat_start = sstat_cat_end;

  } /* End of loop on mesh categories */

  /* Check NaN and exit */
  if (fpe_flag == 1)
    bft_error(__FILE__, __LINE__, 0,
                _("Invalid (not-a-number) values detected for a statistic."));

  BFT_FREE(wsum);
  BFT_FREE(vsum);
  BFT_FREE(vmax);
  BFT_FREE(vmin);

  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*----------------------------------------------------------------------------
 * Add or update clipping info for a given array
 *
 * parameters:
 *   name_id       Associated name id if not a field, -1 for a field
 *   f_id          associated field id, or -1
 *   dim           associated dimension
 *   n_clip_min    number of local clippings to minimum value
 *   n_clip_max    number of local clippings to maximum value
 *   min_pre_clip  minimum values prior to clipping
 *   max_pre_clip  maximum values prior to clipping
 *----------------------------------------------------------------------------*/

static void
_add_clipping(int               name_id,
              int               f_id,
              int               dim,
              cs_lnum_t         n_clip_min,
              cs_lnum_t         n_clip_max,
              const cs_real_t   min_pre_clip[],
              const cs_real_t   max_pre_clip[],
              cs_lnum_t         n_clip_min_comp[],
              cs_lnum_t         n_clip_max_comp[])
{
  bool need_sort = false;

  int clip_id = _find_clip(f_id, name_id);

  /* If not found, insert statistic */

  if (clip_id < 0) {

    _n_clips += 1;
    _clips_val_size += (dim == 1) ? 1 : dim + 1;

    /* Reallocate key definitions if necessary */

    if (_n_clips > _n_clips_max) {
      if (_n_clips_max == 0)
        _n_clips_max = 1;
      else
        _n_clips_max *= 2;
      BFT_REALLOC(_clips, _n_clips_max, cs_log_clip_t);
    }

    if (_clips_val_size > _clips_val_size_max) {
      if (_clips_val_size_max == 0)
        _clips_val_size_max = dim;
      while (_clips_val_size > _clips_val_size_max)
        _clips_val_size_max *= 2;
      BFT_REALLOC(_clips_vmin, _clips_val_size_max, double);
      BFT_REALLOC(_clips_vmax, _clips_val_size_max, double);
      BFT_REALLOC(_clips_count, _clips_val_size_max*2, cs_gnum_t);
    }

    need_sort = true;  /* allow for binary search */

    clip_id = _n_clips - 1;
    _clips[clip_id].f_id = f_id;
    _clips[clip_id].name_id = name_id;
    _clips[clip_id].dim = dim;
    _clips[clip_id].v_idx = (dim == 1) ?
      _clips_val_size - dim : _clips_val_size - dim - 1;

  }

  if (_clips[clip_id].dim != dim) {
    if (f_id > -1)
      bft_error(__FILE__, __LINE__, 0,
                "Clipping of field id %d previously defined in %s\n"
                "with dimension %d, redefined with dimension %d.",
                f_id, __func__,
                _clips[clip_id].dim, dim);
    else
      bft_error(__FILE__, __LINE__, 0,
                "Clipping of name %s previously defined in %s\n"
                "with dimension %d, redefined with dimension %d.",
                cs_map_name_to_id_reverse(_name_map, name_id),
                __func__,
                _clips[clip_id].dim, dim);
  }

  /* Update clips */

  _clips[clip_id].dim = dim;

  int v_idx = _clips[clip_id].v_idx;

  /* Prepare for future binary search */

  if (need_sort)
    qsort(_clips, _n_clips, sizeof(cs_log_clip_t), &_compare_clips);

  /* Update values */
  if (dim > 1) {
    _clips_count[(v_idx)*2] = n_clip_min;
    _clips_count[(v_idx)*2 + 1] = n_clip_max;
    _clips_vmin[v_idx] = min_pre_clip[0];
    _clips_vmax[v_idx] = max_pre_clip[0];
    for (int i = 0; i < dim; i++) {
      _clips_vmin[v_idx + i + 1] = min_pre_clip[i];
      _clips_vmax[v_idx + i + 1] = max_pre_clip[i];
      _clips_count[(v_idx + i + 1)*2] = n_clip_min_comp[i];
      _clips_count[(v_idx + i + 1)*2 + 1] = n_clip_max_comp[i];
    }
  }
  else {
    _clips_vmin[v_idx] = min_pre_clip[0];
    _clips_vmax[v_idx] = max_pre_clip[0];
    _clips_count[(v_idx)*2] = n_clip_min;
    _clips_count[(v_idx)*2 + 1] = n_clip_max;
  }
}

/*----------------------------------------------------------------------------
 * Main logging output of additional clippings
 *----------------------------------------------------------------------------*/

static void
_log_clips(void)
{
  int     clip_id;
  int     type_idx[] = {0, 0, 0};
  double  *vmin = nullptr, *vmax = nullptr;
  cs_gnum_t  *vcount = nullptr;
  size_t max_name_width = cs_log_strlen(_("field"));
  const int label_key_id = cs_field_key_id("label");

  char tmp_s[5][64] =  {"", "", "", "", ""};

  const char *_cat_name[] = {N_("field"), N_("value")};
  const char *_cat_prefix[] = {"a  ", "a   "};

  /* Allocate working arrays */

  BFT_MALLOC(vmin, _clips_val_size, double);
  BFT_MALLOC(vmax, _clips_val_size, double);
  BFT_MALLOC(vcount, _clips_val_size*2, cs_gnum_t);

  memcpy(vmin, _clips_vmin, _clips_val_size*sizeof(double));
  memcpy(vmax, _clips_vmax, _clips_val_size*sizeof(double));
  memcpy(vcount, _clips_count, _clips_val_size*sizeof(cs_gnum_t)*2);

  /* Group MPI operations if required */

  cs_parall_min(_clips_val_size, CS_DOUBLE, vmin);
  cs_parall_max(_clips_val_size, CS_DOUBLE, vmax);
  cs_parall_sum(_clips_val_size*2, CS_GNUM_TYPE, vcount);

  /* Fist loop on clippings for counting */

  for (clip_id = 0; clip_id < _n_clips; clip_id++) {

    const char *name = nullptr;
    int f_id = _clips[clip_id].f_id;
    int f_dim = 0;
    if (f_id > -1) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      name = cs_field_get_key_str(f, label_key_id);
      if (name == nullptr)
        name = f->name;
      type_idx[1] = clip_id + 1;
      f_dim = f->dim;
    }
    else {
      name = cs_map_name_to_id_reverse(_name_map, _clips[clip_id].name_id);
      type_idx[2] = clip_id + 1;
    }

    assert(name != nullptr);

    size_t l_name_width = cs_log_strlen(name);
    if (f_dim == 3)
      l_name_width += 3;
    else if (f_dim > 3)
      l_name_width += 4;
    max_name_width = CS_MAX(max_name_width, l_name_width);

  }

  if (type_idx[2] - type_idx[1] > 0) {
    size_t v_name_w = cs_log_strlen(_("value"));
    max_name_width = CS_MAX(max_name_width, v_name_w);
  }

  max_name_width = CS_MIN(max_name_width, 63);

  /* Loop on types */

  for (int cat_id = 0; cat_id < 2; cat_id++) {

    int start_id = type_idx[cat_id];
    int end_id = type_idx[cat_id+1];

    if (end_id - start_id < 1)
      continue;

    /* Print headers */

    if (cat_id == 0)
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n"
                      "  ** Clippings for computed fields\n"
                      "     -----------------------------\n"));
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n"
                      "  ** Clippings for auxiliary values\n"
                      "     ------------------------------\n"));

    cs_log_strpad(tmp_s[0], _(_cat_name[cat_id]), max_name_width, 64);
    cs_log_strpadl(tmp_s[1], _("initial min"), 14, 64);
    cs_log_strpadl(tmp_s[2], _("initial max"), 14, 64);
    cs_log_strpadl(tmp_s[3], _("clips to min"), 12, 64);
    cs_log_strpadl(tmp_s[4], _("clips to max"), 12, 64);
    cs_log_printf(CS_LOG_DEFAULT,
                    "\n   %s  %s  %s  %s  %s\n",
                    tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);

    /* Underline */

    for (size_t col = 0; col < 5; col++) {
      size_t i;
      size_t w0;
      if (col == 0)
        w0 = max_name_width;
      else if (col < 3)
        w0 = 14;
      else
        w0 = 12;
      for (i = 0; i < w0; i++)
        tmp_s[col][i] = '-';
      tmp_s[col][w0] = '\0';
    }
    cs_log_printf(CS_LOG_DEFAULT,
                  "-  %s  %s  %s  %s  %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);

    /* Second loop on clippings */

    for (clip_id = start_id; clip_id < end_id; clip_id++) {

      int v_idx = _clips[clip_id].v_idx;

      const char *name = nullptr;
      int f_id = _clips[clip_id].f_id;
      if (f_id > -1) {
        const cs_field_t  *f = cs_field_by_id(f_id);
        name = cs_field_get_key_str(f, label_key_id);
        if (name == nullptr)
          name = f->name;
      }
      else
        name = cs_map_name_to_id_reverse(_name_map, _clips[clip_id].name_id);

      assert(name != nullptr);

      _log_clip_info(_cat_prefix[cat_id],
                      name,
                      max_name_width,
                     _clips[clip_id].dim,
                     vcount + v_idx*2,
                     vcount + v_idx*2 + 1,
                     vmin + v_idx,
                     vmax + v_idx);
    }

  }

  BFT_FREE(vcount);
  BFT_FREE(vmax);
  BFT_FREE(vmin);

  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free arrays possible used by logging of array statistics.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_destroy_all(void)
{
  if (_category_map != nullptr) {
    _sstats_val_size = 0;
    _sstats_val_size_max = 0;
    _n_sstats = 0;
    _n_sstats_max = 0;
    BFT_FREE(_sstats_vmin);
    BFT_FREE(_sstats_vmax);
    BFT_FREE(_sstats_vsum);
    BFT_FREE(_sstats_wsum);
    BFT_FREE(_sstats);
    cs_map_name_to_id_destroy(&_category_map);
  }

  if (_n_clips_max > 0) {
    _clips_val_size = 0;
    _clips_val_size_max = 0;
    _n_clips = 0;
    _n_clips_max = 0;
    BFT_FREE(_clips_count);
    BFT_FREE(_clips_vmin);
    BFT_FREE(_clips_vmax);
    BFT_FREE(_clips);
  }

  if (_name_map != nullptr)
    cs_map_name_to_id_destroy(&_name_map);

  if (_l2_residual_plot != nullptr)
    cs_time_plot_finalize(&_l2_residual_plot);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log field and other array statistics for the current time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration(void)
{

  if (_n_clips > 0)
    _log_clips();

  _log_fields_and_functions();

  if (_n_sstats > 0)
    _log_sstats();

  cs_time_moment_log_iteration();
  cs_lagr_stat_log_iteration();
  cs_lagr_log_iteration();

  cs_fan_log_iteration();
  cs_ctwr_log_balance();

  cs_notebook_log();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Default function for equation convergence log info.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_equation_convergence_info_write(void)
{
  if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_ONLY)
    return;

  const int n_fields = cs_field_n_fields();
  const int keylog = cs_field_key_id("log");
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_t *w1, *w2;
  BFT_MALLOC(w1, n_cells, cs_real_t);
  BFT_MALLOC(w2, n_cells, cs_real_t);

  char title[128] = "   Variable    ";
  char line[128] = "";
  int ic = 15;

  /* Compute largest name width to adjust log output */
  int max_name_width = 12;
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      const char *f_label = cs_field_get_label(f);
      const int f_dim = f->dim;
      int length = strlen(f_label);

      if (f_dim == 3)
        length += 3;
      else if (f_dim == 6 || f_dim == 9)
        length += 4;

      if (length > max_name_width)
        max_name_width = length;
    }
  }

  for (int i = ic; i < ic + fmax(0, max_name_width - 11); i++)
    strcat(title, " ");
  strcat(title, "Rhs norm      N_iter  Norm. residual   Drift   Time residual");

  for (int i = 0; i < (int) strlen(title); i++)
    strcat(line, "-");

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       "  ** Information on convergence\n"
       "     --------------------------\n\n"));

  cs_log_printf(CS_LOG_DEFAULT, _("%s\n%s\n%s\n"), line, title, line);

  /* Print convergence information for each solved field */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    cs_real_t dervar[9], varres[9], varnrm[9];
    assert(f->dim <= 9);

    int log_flag = cs_field_get_key_int(f, keylog);

    if (!(f->type & CS_FIELD_VARIABLE))
      continue;

    char chain[128] = "c  ";

    const char *f_label = cs_field_get_label(f);
    char var_name_trunc[max_name_width + 1];
    strncpy(var_name_trunc, f_label, max_name_width);
    var_name_trunc[max_name_width] = '\0';
    strcat(chain, var_name_trunc);

    for (int i = strlen(chain); i < 3 + max_name_width; i++)
      strcat(chain, " ");

    /* Check if the variable was solved in the current time step */
    cs_solving_info_t sinfo;
    cs_field_get_key_struct(f, cs_field_key_id("solving_info"), &sinfo);
    if (sinfo.n_it < 0)
      continue;

    /* Compute the time drift */
    /* Cell based variables */
    if (f->location_id == CS_MESH_LOCATION_CELLS) {
      const int dim = f->dim;
      cs_real_t *dt = CS_F_(dt)->val;

      /* Pressure time drift (computed in cs_pressure_correction.c) */
      dervar[0] = sinfo.derive;

        /* Time drift for cell based variables (except pressure) */
      if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1
          || strcmp(f->name, "pressure") != 0) {
        for (int isou = 0; isou < dim; isou++) {
          for (int c_id = 0; c_id < n_cells; c_id++)
            w1[c_id] =   (f->val[dim*c_id + isou] - f->val_pre[dim*c_id + isou])
                       / sqrt(dt[c_id]);

          dervar[isou] = cs_gres(n_cells, cell_vol, w1, w1);
        }

        for (int isou = 1; isou < dim; isou++)
          dervar[0] += dervar[isou];
        /* We don't update the sinfo attribute since it may not be
         * updated at each time step (only when logging)
         * NOTE: it should be added in the bindings again
         * if needed */
        // sinfo.derive = dervar[0];
      }

      /* L2 time normalized residual */
      for (int isou = 0; isou < dim; isou++) {
        for (int c_id = 0; c_id < n_cells; c_id++) {
          w1[c_id] =   (f->val[dim*c_id + isou] - f->val_pre[dim*c_id + isou])
                     / dt[c_id];
          w2[c_id] = f->val[dim * c_id + isou];
        }

        varres[isou] = cs_gres(n_cells, cell_vol, w1, w1);
        varnrm[isou] = cs_gres(n_cells, cell_vol, w2, w2);

        if (isou > 0) {
          varres[0] += varres[isou];
          varnrm[0] += varnrm[isou];
        }
      }

      if (varnrm[0] > 0.)
        varres[0] = varres[0] / varnrm[0];
      /* We don't update the sinfo attribute since it may not be
       * updated at each time step (only when logging)
       * NOTE: it should be added in the bindings again
       * if needed */
      // sinfo.l2residual = sqrt(cs_math_fabs(varres[0]));
    }

    char var_log[128];
    snprintf(var_log, 127, "%12.5e %7d   %12.5e %12.5e %12.5e",
             sinfo.rhs_norm, sinfo.n_it, sinfo.res_norm,
             dervar[0], sqrt(cs_math_fabs(varres[0])));

    strcat(chain, var_log);

    if (log_flag > 0)
      cs_log_printf(CS_LOG_DEFAULT, _("%s\n"), chain);

  } /* End loop on fields */

  cs_log_printf(CS_LOG_DEFAULT, _("%s\n"), line);

  BFT_FREE(w1);
  BFT_FREE(w2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set adaptive interval for "per time step" logging information.
 *
 * Logging will also occur:
 * - Each time step for the first n absolute or restarted time steps.
 * - Every 5 time steps between n and 5.n time steps.
 * - Every 10 time steps between 5.n and 10.n time steps.
 * - Every 50 time steps between 10.n and 50.n time steps.
 * - Every 100 time steps between 50.n and 100.n time steps.
 * - ...
 * - At the last time step\n\n"),
 *
 * \param[in]  n  base interval for output.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_set_automatic(int  n)
{
  _log_interval_base = n;

  cs_time_control_init_by_func(cs_glob_log_iteration_time_control,
                               _log_time_control_automatic,
                               &_log_interval_base,
                               true,
                               true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set interval for "per time step" logging information.
 *
 * Logging will also occur for the 10 first time steps, as well as the last one.
 *
 * \param[in]  n  interval between 2 output time steps.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_set_interval(int  n)
{
  _log_interval_base = n;

  cs_time_control_init_by_func(cs_glob_log_iteration_time_control,
                               _log_time_control_interval,
                               &_log_interval_base,
                               true,
                               true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate or deactivate default log for current iteration.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_set_active(void)
{
  bool is_active
    = cs_time_control_is_active(cs_glob_log_iteration_time_control,
                                cs_glob_time_step);

  cs_log_default_activate(is_active);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add or update array not saved as permanent field to iteration log.
 *
 * \param[in]  name         array name
 * \param[in]  category     category name
 * \param[in]  loc_id       associated mesh location id
 * \param[in]  is_intensive are the matching values intensive ?
 * \param[in]  dim          associated dimension (interleaved)
 * \param[in]  val          associated values
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_add_array(const char                     *name,
                           const char                     *category,
                           const cs_mesh_location_type_t   loc_id,
                           bool                            is_intensive,
                           int                             dim,
                           const cs_real_t                 val[])
{
  /* Initialize if necessary */

  if (_name_map == nullptr)
    _name_map = cs_map_name_to_id_create();

  if (_category_map == nullptr)
    _category_map = cs_map_name_to_id_create();

  /* Find or insert entries in map */

  bool need_sort = false;

  int cat_id = cs_map_name_to_id(_category_map, category);
  int name_id = cs_map_name_to_id(_name_map, name);

  int sstat_id = _find_sstats(cat_id, name_id);

  /* If not found, insert statistic */

  if (sstat_id < 0) {

    int _dim = (dim == 3) ? 4 : dim;

    _n_sstats += 1;
    _sstats_val_size += _dim;

    /* Reallocate key definitions if necessary */

    if (_n_sstats > _n_sstats_max) {
      if (_n_sstats_max == 0)
        _n_sstats_max = 1;
      else
        _n_sstats_max *= 2;
      BFT_REALLOC(_sstats, _n_sstats_max, cs_log_sstats_t);
    }

    if (_sstats_val_size > _sstats_val_size_max) {
      if (_sstats_val_size_max == 0)
        _sstats_val_size_max = dim;
      while (_sstats_val_size > _sstats_val_size_max)
        _sstats_val_size_max *= 2;
      BFT_REALLOC(_sstats_vmin, _sstats_val_size_max, double);
      BFT_REALLOC(_sstats_vmax, _sstats_val_size_max, double);
      BFT_REALLOC(_sstats_vsum, _sstats_val_size_max, double);
      BFT_REALLOC(_sstats_wsum, _sstats_val_size_max, double);
    }

    need_sort = true;  /* allow for binary search */

    sstat_id = _n_sstats - 1;
    _sstats[sstat_id].name_id = name_id;
    _sstats[sstat_id].cat_id = cat_id;
    _sstats[sstat_id].dim = dim;
    _sstats[sstat_id].v_idx = _sstats_val_size - _dim;

  }

  if (_sstats[sstat_id].dim != dim)
    bft_error(__FILE__, __LINE__, 0,
              "Array of name %s and category %s previously defined in %s\n"
              "with dimension %d, redefined with dimension %d.",
              cs_map_name_to_id_reverse(_name_map, name_id),
              cs_map_name_to_id_reverse(_category_map, cat_id),
              __func__,
              _sstats[sstat_id].dim, dim);

  /* Update stats */

  _sstats[sstat_id].loc_id = loc_id;
  _sstats[sstat_id].intensive = is_intensive;
  _sstats[sstat_id].dim = dim;

  int v_idx = _sstats[sstat_id].v_idx;

  /* Prepare for future binary search */

  if (need_sort)
    qsort(_sstats, _n_sstats, sizeof(cs_log_sstats_t), &_compare_sstats);

  /* Compute sstats */

  bool have_weight = false;
  const cs_real_t *weight = nullptr;

  if (is_intensive != false) {
    cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
    switch(loc_id) {
    case CS_MESH_LOCATION_CELLS:
      weight = mq->cell_vol;
      have_weight = true;
      break;
    case CS_MESH_LOCATION_INTERIOR_FACES:
      weight = mq->i_face_surf;
      have_weight = true;
      break;
    case CS_MESH_LOCATION_BOUNDARY_FACES:
      weight = mq->b_face_surf;
      have_weight = true;
      break;
    case CS_MESH_LOCATION_VERTICES:
      have_weight = false;
      break;
    default:
      break;
    }
  }

  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(loc_id);
  const cs_lnum_t *elt_list = cs_mesh_location_get_elt_ids_try(loc_id);

  if (have_weight)
    cs_array_reduce_simple_stats_l_w(n_elts[0],
                                     dim,
                                     elt_list,
                                     elt_list,
                                     val,
                                     weight,
                                     _sstats_vmin + v_idx,
                                     _sstats_vmax + v_idx,
                                     _sstats_vsum + v_idx,
                                     _sstats_wsum + v_idx);
  else {
    cs_array_reduce_simple_stats_l(n_elts[0],
                                   dim,
                                   elt_list,
                                   val,
                                   _sstats_vmin + v_idx,
                                   _sstats_vmax + v_idx,
                                   _sstats_vsum + v_idx);
    for (int i = 0; i < dim; i++)
      _sstats_wsum[v_idx + i] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add or update clipping info for a given array.
 *
 * \param[in]  name          field or array name
 * \param[in]  dim           associated dimension
 * \param[in]  n_clip_min    number of local clippings to minimum value
 * \param[in]  n_clip_max    number of local clippings to maximum value
 * \param[in]  min_pre_clip  minimum values prior to clipping
 * \param[in]  max_pre_clip  maximum values prior to clipping
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_clipping(const char       *name,
                          int               dim,
                          cs_lnum_t         n_clip_min,
                          cs_lnum_t         n_clip_max,
                          const cs_real_t   min_pre_clip[],
                          const cs_real_t   max_pre_clip[])
{
  /* Initialize if necessary */

  if (_name_map == nullptr)
    _name_map = cs_map_name_to_id_create();

  int name_id = cs_map_name_to_id(_name_map, name);

  _add_clipping(name_id, -1, dim,
                n_clip_min, n_clip_max,
                min_pre_clip, max_pre_clip, nullptr, nullptr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add or update clipping info for a field
 *
 * \param[in]  f_id            associated field id
 * \param[in]  n_clip_min      number of local clippings to minimum value
 * \param[in]  n_clip_max      number of local clippings to maximum value
 * \param[in]  min_pre_clip    minimum values prior to clipping
 * \param[in]  max_pre_clip    maximum values prior to clipping
 * \param[in]  n_clip_min_comp number of clip min by component
 * \param[in]  n_clip_max_comp number of clip max by component
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_clipping_field(int               f_id,
                                cs_lnum_t         n_clip_min,
                                cs_lnum_t         n_clip_max,
                                const cs_real_t   min_pre_clip[],
                                const cs_real_t   max_pre_clip[],
                                cs_lnum_t         n_clip_min_comp[],
                                cs_lnum_t         n_clip_max_comp[])
{
  const cs_field_t  *f = cs_field_by_id(f_id);

  _add_clipping(-1, f_id, f->dim,
                n_clip_min, n_clip_max,
                min_pre_clip, max_pre_clip,n_clip_min_comp,n_clip_max_comp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize structures used for logging for new iteration.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_prepare(void)
{
  const int n_fields = cs_field_n_fields();

  int si_k_id = cs_field_key_id("solving_info");

  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      auto *sinfo = static_cast<cs_solving_info_t *>(
        cs_field_get_key_struct_ptr(f, si_k_id));
      sinfo->n_it = -1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log L2 time residual for variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_l2residual(void)
{
  if (cs_glob_rank_id > 0)
    return;

  const cs_time_step_t *ts = cs_glob_time_step;
  const int n_fields = cs_field_n_fields();

  /* write header */

  if (_l2_residual_plot == nullptr) {

    int                    _plot_buffer_steps = -1;
    double                 _plot_flush_wtime = 3600;
    cs_time_plot_format_t  _plot_format = CS_TIME_PLOT_CSV;
    bool                   use_iteration = (ts->is_local) ? true : false;

    const char **labels;
    BFT_MALLOC(labels, n_fields + 1, const char *);

    int n_variables = 0;
    for (int f_id = 0 ; f_id < n_fields ; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        labels[n_variables] = f->name;
        n_variables++;
      }
    }

    _l2_residual_plot = cs_time_plot_init_probe("residuals",
                                                "",
                                                _plot_format,
                                                use_iteration,
                                                _plot_flush_wtime,
                                                _plot_buffer_steps,
                                                n_variables,
                                                nullptr,
                                                nullptr,
                                                labels);

    BFT_FREE(labels);
  }

  {
    cs_real_t *vals;
    BFT_MALLOC(vals, n_fields, cs_real_t);

    int si_k_id = cs_field_key_id("solving_info");

    int n_variables = 0;
    for (int f_id = 0 ; f_id < n_fields ; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        auto sinfo = static_cast<const cs_solving_info_t *>(
          cs_field_get_key_struct_const_ptr(f, si_k_id));
        vals[n_variables] = sinfo->l2residual;
        n_variables += 1;
      }
    }

    cs_time_plot_vals_write(_l2_residual_plot,
                            ts->nt_cur,
                            ts->t_cur,
                            n_variables,
                            vals);

    BFT_FREE(vals);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print default log per iteration options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_log_setup(void)
{
  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "Logging options\n"
       "---------------\n\n"));

  cs_log_printf
    (CS_LOG_SETUP,
     _("  run_solver.log output interval:\n\n"));

  if (   cs_glob_log_iteration_time_control->control_func
      == _log_time_control_automatic) {
    int n = _log_interval_base;
    if (n == 0)
      n = 10;
    cs_log_printf
      (CS_LOG_SETUP,
       _("    Automatic:\n"
         "     - Each time step for the first %d absolute "
         "or restarted time steps.\n"
         "     - Every 5 time steps between %d and %d time steps.\n"
         "     - Every 10 time steps between %d and %d time steps.\n"
         "     - Every 50 time steps between %d and %d time steps.\n"
         "     - Every 100 time steps between %d and %d time steps.\n"
         "     - ...\n"
         "     - At the last time step.\n"),
       n, n, 5*n, 5*n, 10*n, 10*n, 50*n, 50*n, 100*n);
  }
  else if (   cs_glob_log_iteration_time_control->control_func
           == _log_time_control_interval) {
    int n = _log_interval_base;
    if (n > 0)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Interval-based:\n"
           "     - Each time step for the first 10 absolute "
           "or restarted time steps.\n"
           "     - Every %d time step(s).\n"
           "     - At the last time step.\n\n"), n);
    else
      cs_log_printf
        (CS_LOG_SETUP,
         _("    At the last time step only.\n"));
  }
  else {
    char desc[256];
    cs_time_control_get_description(cs_glob_log_iteration_time_control,
                                    desc,
                                    256);
    cs_log_printf(CS_LOG_SETUP,
                  _("    %s\n"), desc);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
