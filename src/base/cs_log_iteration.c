/*============================================================================
 * Log field and other array statistics at relevant time steps.
 *============================================================================*/

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
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_range_set.h"
#include "cs_time_moment.h"
#include "cs_time_plot.h"
#include "cs_time_step.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_log.h"


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

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_map_name_to_id_t  *_name_map = NULL;

static cs_map_name_to_id_t  *_category_map = NULL;
static int _sstats_val_size = 0;
static int _sstats_val_size_max = 0;
static int _n_sstats = 0;
static int _n_sstats_max = 0;
static double *_sstats_vmin = NULL;
static double *_sstats_vmax = NULL;
static double *_sstats_vsum = NULL;
static double *_sstats_wsum = NULL;
static cs_log_sstats_t  *_sstats = NULL;

static int _clips_val_size = 0;
static int _clips_val_size_max = 0;
static int _n_clips = 0;
static int _n_clips_max = 0;
static cs_gnum_t *_clips_count = NULL;
static double *_clips_vmin = NULL;
static double *_clips_vmax = NULL;
static cs_log_clip_t  *_clips = NULL;

static cs_time_plot_t  *_l2_residual_plot = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran function prototypes for subroutines from field.f90.
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

static int _compare_sstats(const void *x, const void *y)
{
  int retval = 1;

  const cs_log_sstats_t *s0 = x;
  const cs_log_sstats_t *s1 = y;

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

  const cs_log_clip_t *c0 = x;
  const cs_log_clip_t *c1 = y;

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
 *   wsum         <-- weighted sum of each component or norm, or NULL
 *   fpe_flag     <-- was a "not a number" or floating-point error detected ?
 *----------------------------------------------------------------------------*/

static void
_log_array_info(const char        *prefix,
                const char        *name,
                size_t             name_width,
                int                dim,
                int                n_g_elts,
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
_log_fields(void)
{
  int f_id, li, log_count;

  int log_count_max = 0;
  int fpe_flag = 0;
  int     *log_id = NULL, *moment_id = NULL;
  double  *vmin = NULL, *vmax = NULL, *vsum = NULL, *wsum = NULL;

  char tmp_s[5][64] =  {"", "", "", "", ""};

  const char _underline[] = "---------------------------------";
  const int n_fields = cs_field_n_fields();
  const int n_moments = cs_time_moment_n_moments();
  const int log_key_id = cs_field_key_id("log");
  const int label_key_id = cs_field_key_id("label");

  const cs_mesh_location_type_t m_l[] = {CS_MESH_LOCATION_CELLS,
                                         CS_MESH_LOCATION_INTERIOR_FACES,
                                         CS_MESH_LOCATION_BOUNDARY_FACES,
                                         CS_MESH_LOCATION_VERTICES
                                         };

  cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Allocate working arrays */

  log_count_max = n_fields;

  BFT_MALLOC(log_id, n_fields, int);
  BFT_MALLOC(vmin, log_count_max, double);
  BFT_MALLOC(vmax, log_count_max, double);
  BFT_MALLOC(vsum, log_count_max, double);
  BFT_MALLOC(wsum, log_count_max, double);

  if (n_moments > 0) {
    BFT_MALLOC(moment_id, n_fields, int);
    for (f_id = 0; f_id < n_fields; f_id++)
      moment_id[f_id] = -1;
    for (int m_id = 0; m_id < n_moments; m_id++) {
      const cs_field_t *f = cs_time_moment_get_field(m_id);
      if (f != NULL)
        moment_id[f->id] = m_id;
    }
  }

  /* Loop on locations */

  for (li = 0; li < 4; li++) {

    size_t max_name_width = cs_log_strlen(_("field"));
    int loc_id = m_l[li];
    int have_weight = 0;
    double total_weight = -1;
    cs_gnum_t n_g_elts = 0;
    cs_real_t *gather_array = NULL; /* only if CS_MESH_LOCATION_VERTICES */
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(loc_id);
    const cs_lnum_t _n_elts = n_elts[0];
    const cs_real_t *weight = NULL;

    if (mq != NULL) {
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
        cs_array_reduce_sum_l(_n_elts, 1, NULL, weight, &total_weight);
        cs_parall_sum(1, CS_DOUBLE, &total_weight);
        have_weight = 1;
        break;
      case CS_MESH_LOCATION_BOUNDARY_FACES:
        n_g_elts = m->n_g_b_faces;
        weight = mq->b_face_surf;
        cs_array_reduce_sum_l(_n_elts, 1, NULL, weight, &total_weight);
        cs_parall_sum(1, CS_DOUBLE, &total_weight);
        have_weight = 1;
        break;
      case CS_MESH_LOCATION_VERTICES:
        n_g_elts = m->n_g_vertices;
        have_weight = 0;
        BFT_MALLOC(gather_array, m->n_vertices, cs_real_t);
        break;
      default:
        n_g_elts = _n_elts;
        cs_parall_counter(&n_g_elts, 1);
        break;
      }
    }

    if (n_g_elts == 0)
      continue;

    /* First loop on fields */

    log_count = 0;

    for (f_id = 0; f_id < n_fields; f_id++) {

      int _dim, c_id;

      const cs_field_t  *f = cs_field_by_id(f_id);

      if (f->location_id != loc_id || ! (cs_field_get_key_int(f, log_key_id))) {
        log_id[f_id] = -1;
        continue;
      }

      /* Only log active moments */

      if (moment_id != NULL) {
        if (moment_id[f_id] > -1) {
          if (!cs_time_moment_is_active(moment_id[f_id])) {
            log_id[f_id] = -1;
            continue;
          }
        }
      }

      /* Position in log */

      log_id[f_id] = log_count;

      _dim = (f->dim == 3) ? 4 : f->dim;

      while (log_count + _dim > log_count_max) {
        log_count_max *= 2;
        BFT_REALLOC(vmin, log_count_max, double);
        BFT_REALLOC(vmax, log_count_max, double);
        BFT_REALLOC(vsum, log_count_max, double);
        BFT_REALLOC(wsum, log_count_max, double);
      }

      if (have_weight && (f->type & CS_FIELD_INTENSIVE)) {
        cs_array_reduce_simple_stats_l_w(_n_elts,
                                         f->dim,
                                         NULL,
                                         NULL,
                                         f->val,
                                         weight,
                                         vmin + log_count,
                                         vmax + log_count,
                                         vsum + log_count,
                                         wsum + log_count);

      }
      else {

        cs_real_t  *field_val = f->val;
        cs_lnum_t  _n_cur_elts = _n_elts;
        if (gather_array != NULL) { /* Eliminate shared values whose local
                                       rank is not owner and compact */

          if (m->vtx_range_set == NULL)
            m->vtx_range_set = cs_range_set_create(m->vtx_interfaces,
                                                   NULL,
                                                   m->n_vertices,
                                                   false, /* balance */
                                                   0); /* g_id_base */

          if (f->dim > 1)
            BFT_REALLOC(gather_array, (f->dim * m->n_vertices), cs_real_t);

          cs_range_set_gather(m->vtx_range_set,
                              CS_REAL_TYPE,
                              f->dim,
                              f->val,
                              gather_array);
          field_val = gather_array;
          _n_cur_elts = m->vtx_range_set->n_elts[0];
        }
        cs_array_reduce_simple_stats_l(_n_cur_elts,
                                       f->dim,
                                       NULL,
                                       field_val,
                                       vmin + log_count,
                                       vmax + log_count,
                                       vsum + log_count);

        if (have_weight) {
          for (c_id = 0; c_id < _dim; c_id++)
            wsum[log_count + c_id] = 0.;
        }
      }

      log_count += _dim;

      const char *name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;

      size_t l_name_width = cs_log_strlen(name);
      if (f->dim == 3)
        l_name_width += 3;
      else if (f->dim > 3)
        l_name_width += 4;

      max_name_width = CS_MAX(max_name_width, l_name_width);

    } /* End of first loop on fields */

    if (gather_array != NULL)
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

    /* Second loop on fields */

    log_count = 0;

    for (f_id = 0; f_id < n_fields; f_id++) {

      int _dim;
      if (log_id[f_id] < 0)
        continue;

      const cs_field_t  *f = cs_field_by_id(f_id);

      /* Position in log */

      _dim = (f->dim == 3) ? 4 : f->dim;

      const char *name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;

      double t_weight = -1;
      if (total_weight > 0 && (f->type & CS_FIELD_INTENSIVE))
        t_weight = total_weight;

      char prefix[] = "v  ";
      if (moment_id != NULL) {
        if (moment_id[f_id] > -1)
          prefix[0] = 'm';
      }
      if (f->type & CS_FIELD_ACCUMULATOR)
        prefix[0] = 'm';

      _log_array_info(prefix,
                      name,
                      max_name_width,
                      f->dim,
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
  double  *vmin = NULL, *vmax = NULL, *vsum = NULL, *wsum = NULL;

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
      const cs_real_t *weight = NULL;
      const char *loc_name = _(cs_mesh_location_get_name(loc_id));
      size_t loc_name_w = cs_log_strlen(loc_name);

      if (mq != NULL) {
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
            cs_array_reduce_sum_l(_n_elts, 1, NULL, weight, &_interior_surf);
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
            cs_array_reduce_sum_l(_n_elts, 1, NULL, weight, &_boundary_surf);
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
    _clips[clip_id].v_idx = (dim == 1) ? _clips_val_size - dim : _clips_val_size - dim - 1;

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
  double  *vmin = NULL, *vmax = NULL;
  cs_gnum_t  *vcount = NULL;
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

    const char *name = NULL;
    int f_id = _clips[clip_id].f_id;
    int f_dim = 0;
    if (f_id > -1) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;
      type_idx[1] = clip_id + 1;
      f_dim = f->dim;
    }
    else {
      name = cs_map_name_to_id_reverse(_name_map, _clips[clip_id].name_id);
      type_idx[2] = clip_id + 1;
    }

    assert(name != NULL);

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

      const char *name = NULL;
      int f_id = _clips[clip_id].f_id;
      if (f_id > -1) {
        const cs_field_t  *f = cs_field_by_id(f_id);
        name = cs_field_get_key_str(f, label_key_id);
        if (name == NULL)
          name = f->name;
      }
      else
        name = cs_map_name_to_id_reverse(_name_map, _clips[clip_id].name_id);

      assert(name != NULL);

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

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

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
  if (_category_map != NULL) {
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

  if (_name_map != NULL)
    cs_map_name_to_id_destroy(&_name_map);

  if (_l2_residual_plot != NULL)
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

  _log_fields();

  if (_n_sstats > 0)
    _log_sstats();

  cs_time_moment_log_iteration();
  cs_lagr_stat_log_iteration();
  cs_lagr_log_iteration();

  cs_fan_log_iteration();
  cs_ctwr_log_balance();
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

  if (_name_map == NULL)
    _name_map = cs_map_name_to_id_create();

  if (_category_map == NULL)
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
  const cs_real_t *weight = NULL;

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
  const cs_lnum_t *elt_list = cs_mesh_location_get_elt_list(loc_id);

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

  if (_name_map == NULL)
    _name_map = cs_map_name_to_id_create();

  int name_id = cs_map_name_to_id(_name_map, name);

  _add_clipping(name_id, -1, dim,
                n_clip_min, n_clip_max,
                min_pre_clip, max_pre_clip,0,0);
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
 * \brief Log L2 time residual for every variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_l2residual(void)
{
  if (cs_glob_rank_id > 0)
    return;

  const cs_time_step_t *ts = cs_glob_time_step;
  const int n_fields = cs_field_n_fields();

  /* write header */

  if (_l2_residual_plot == NULL) {

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
                                                NULL,
                                                NULL,
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
        const cs_solving_info_t *sinfo
          = cs_field_get_key_struct_const_ptr(f, si_k_id);
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

END_C_DECLS
