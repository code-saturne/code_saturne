/*============================================================================
 * Log field and other array statistics at relevant time steps.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
#include "cs_field.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"

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

/*============================================================================
 * Static global variables
 *============================================================================*/

static int _sstats_val_size = 0;
static int _sstats_val_size_max = 0;
static int _n_sstats = 0;
static int _n_sstats_max = 0;
static double *_sstats_vmin = NULL;
static double *_sstats_vmax = NULL;
static double *_sstats_vsum = NULL;
static double *_sstats_wsum = NULL;
static cs_log_sstats_t  *_sstats = NULL;
static cs_map_name_to_id_t  *_category_map = NULL;
static cs_map_name_to_id_t  *_name_map = NULL;

/* Additional variable for moments */

static const cs_real_t  *_cumulative_mom_time = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_log_iteration(void);

/*! \endcond (end ignore by Doxygen) */

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
    if (s0->cat_id < s1->cat_id)
      retval = -1;
    else if (s0->cat_id == s1->cat_id)
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
 * Build moments array for post-processing.
 *
 * If of dimension > 1, the moments array is always interleaved, whether
 * the accumulator field is interleaved or not.
 *
 * parameters:
 *   f          <-- pointer to field structure
 *   moment_id  <-- id of associated moment divisor:
 *                  - if moment_id == -1, the field is not a moment;
 *                  - if moment_id >= 0, it is the field id for the divisor;
 *                  - if moment_id < -1, (-1 -moment_id) is the moment id
 *                    in the Fortran "dtcmom" array of the optcal module
 *   n_elts     <-- local number of elements
 *   moment     --> resulting moment array (size: n_elts)
 *----------------------------------------------------------------------------*/

static void
_build_moment(const cs_field_t  *f,
              int                moment_id,
              cs_lnum_t          n_elts,
              cs_real_t          moment[])
{
  cs_lnum_t i, j;
  cs_lnum_t  d_mult = 0;
  const cs_real_t ep_zero = 1.e-12;
  const cs_real_t *denom = NULL;

  assert(moment_id != -1);

  if (moment_id > -1) {
    const cs_field_t  *fd = cs_field_by_id(moment_id);
    assert(fd->dim == 1);
    denom = fd->val;
    d_mult = 1;
  }
  else if (moment_id < -1) {
    denom = &(_cumulative_mom_time[(- moment_id - 1) - 1]);
    /* d_mult = 0 is set above */
  }

  if (f->dim == 1) {
    for (i = 0; i < n_elts; i++)
      moment[i] = f->val[i] / CS_MAX(denom[i*d_mult], ep_zero);
  }
  else {
    cs_lnum_t i_mult = 1, j_mult = 1;
    if (f->interleaved)
      i_mult = f->dim;
    else {
      const cs_lnum_t *n_loc_elts
        = cs_mesh_location_get_n_elts(f->location_id);
      j_mult = n_loc_elts[2];
    }
    for (i = 0; i < n_elts; i++) {
      for (j = 0; j < f->dim; j++)
        moment[i*f->dim + j] =   f->val[i*i_mult + j*j_mult]
                                 / CS_MAX(denom[i*d_mult], ep_zero);
    }
  }
}

/*----------------------------------------------------------------------------
 * Log information for a given array.
 *
 * parameters:
 *   prefix       <-- string inserted before name
 *   name         <-- field name or label
 *   name_width   <-- width of "name" column
 *   interleaved  <-- is field interleaved ?
 *   dim          <-- field dimension
 *   n_g_elts     <-- global number of associated elements,
 *   total_weight <-- if > 0, weight (volume or surface) of field location
 *   vmin         <-- minimum values of each component or norm
 *   vmax         <-- maximum values of each component or norm
 *   vsum         <-- sum of each component or norm
 *   wsum         <-- weighted sum of each component or norm, or NULL
 *----------------------------------------------------------------------------*/

static void
_log_array_info(const char        *prefix,
                const char        *name,
                size_t             name_width,
                bool               interleaved,
                int                dim,
                int                n_g_elts,
                double             total_weight,
                double             vmin[],
                const double       vmax[],
                const double       vsum[],
                const double      *wsum)
{
  int c_id;
  const int _dim = (interleaved && dim == 3) ? 4 : dim;

  char tmp_s[2][64] =  {"", ""};
  const char *ext_3[] = {"[x]", "[y]", "[z]", ""};
  const char *ext_6[] = {"[11]", "[22]", "[33]", "[12]", "[13]", "[23]"};

  for (c_id = 0; c_id < _dim; c_id++) {

    if (dim == 3) {
      snprintf(tmp_s[1], 63, "%s%s", name, ext_3[c_id]);
      tmp_s[1][63] = '\0';
      cs_log_strpad(tmp_s[0], tmp_s[1], name_width, 64);
    }
    else if (dim == 6) {
      snprintf(tmp_s[1], 63, "%s%s", name, ext_6[c_id]);
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
                    tmp_s[c_id],
                    vmin[c_id],
                    vmax[c_id],
                    vsum[c_id] / n_g_elts);

  }
}

/*----------------------------------------------------------------------------
 * Main post-processing output of variables.
 *----------------------------------------------------------------------------*/

static void
_cs_log_fields(void)
{
  int f_id, li, log_count;

  int log_count_max = 0;
  int     *log_id = NULL;
  double  *vmin = NULL, *vmax = NULL, *vsum = NULL, *wsum = NULL;

  char tmp_s[5][64] =  {"", "", "", "", ""};

  const char _underline[] = "---------------------------------";
  const int n_fields = cs_field_n_fields();
  const int log_key_id = cs_field_key_id("log");
  const int label_key_id = cs_field_key_id("label");
  const int moment_key_id = cs_field_key_id("moment_dt");

  const cs_mesh_location_type_t m_l[] = {CS_MESH_LOCATION_CELLS,
                                         CS_MESH_LOCATION_BOUNDARY_FACES};

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Allocate working arrays */

  log_count_max = n_fields;

  BFT_MALLOC(log_id, n_fields, int);
  BFT_MALLOC(vmin, log_count_max, double);
  BFT_MALLOC(vmax, log_count_max, double);
  BFT_MALLOC(vsum, log_count_max, double);
  BFT_MALLOC(wsum, log_count_max, double);

  /* Loop on locations */

  for (li = 0; li < 2; li++) {

    size_t max_name_width = cs_log_strlen(_("field"));
    int loc_id = m_l[li];
    int have_weight = 0;
    double total_weight = -1;
    cs_gnum_t n_g_elts = 0;
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(loc_id);;
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
      case CS_MESH_LOCATION_BOUNDARY_FACES:
        n_g_elts = m->n_g_b_faces;
        weight = mq->b_face_surf;
        cs_array_reduce_sum_l(_n_elts, 1, NULL, weight, &total_weight);
        cs_parall_sum(1, CS_DOUBLE, &total_weight);
        have_weight = 1;
        break;
      default:
        n_g_elts = _n_elts;
        cs_parall_counter(&n_g_elts, 1);
        break;
      }
    }

    if (have_weight)

    /* First loop on fields */

    log_count = 0;

    for (f_id = 0; f_id < n_fields; f_id++) {

      int _dim, c_id;

      cs_real_t *vmom = NULL;

      const cs_field_t  *f = cs_field_by_id(f_id);
      const cs_real_t *val = f->val;

      if (f->location_id != loc_id || ! (cs_field_get_key_int(f, log_key_id))) {
        log_id[f_id] = -1;
        continue;
      }

      /* Position in log */

      log_id[f_id] = log_count;

      _dim = (f->interleaved && f->dim == 3) ? 4 : f->dim;

      while (log_count + _dim > log_count_max) {
        log_count_max *= 2;
        BFT_REALLOC(vmin, log_count_max, double);
        BFT_REALLOC(vmax, log_count_max, double);
        BFT_REALLOC(vsum, log_count_max, double);
        BFT_REALLOC(wsum, log_count_max, double);
      }

      /* A property field might be a moment */

      if (f->type & CS_FIELD_ACCUMULATOR) {

        int moment_id = cs_field_get_key_int(f, moment_key_id);

        /* if moment_id == -1, the field is not a moment;
           if moment_id > 0, it is the field id for the divisor;
           if moment_id < -1, (-1 -moment_id) is the moment id in "dtcmom" */

        if (moment_id != -1) {
          BFT_MALLOC(vmom, _n_elts*f->dim, cs_real_t);
          _build_moment(f, moment_id, _n_elts, vmom);
        }

        val = vmom;

      } /* End of test on properties/moments */

      if (have_weight && (f->type | CS_FIELD_INTENSIVE)) {
        if (f->interleaved)
          cs_array_reduce_simple_stats_l_w(_n_elts,
                                           f->dim,
                                           NULL,
                                           NULL,
                                           val,
                                           weight,
                                           vmin + log_count,
                                           vmax + log_count,
                                           vsum + log_count,
                                           wsum + log_count);

        else {
          for (c_id = 0; c_id < f->dim; c_id++)
            cs_array_reduce_simple_stats_l_w(_n_elts,
                                             1,
                                             NULL,
                                             NULL,
                                             val + n_elts[2]*c_id,
                                             weight,
                                             vmin + log_count + c_id,
                                             vmax + log_count + c_id,
                                             vsum + log_count + c_id,
                                             wsum + log_count + c_id);
        }
      }
      else {
        if (f->interleaved)
          cs_array_reduce_simple_stats_l(_n_elts,
                                         f->dim,
                                         NULL,
                                         val,
                                         vmin + log_count,
                                         vmax + log_count,
                                         vsum + log_count);

        else {
          for (c_id = 0; c_id < f->dim; c_id++)
            cs_array_reduce_simple_stats_l(_n_elts,
                                           1,
                                           NULL,
                                           val + n_elts[2]*c_id,
                                           vmin + log_count + c_id,
                                           vmax + log_count + c_id,
                                           vsum + log_count + c_id);
        }
        if (have_weight) {
          for (c_id = 0; c_id < _dim; c_id++)
            wsum[log_count + c_id] = 0.;
        }
      }

      log_count += _dim;

      val = NULL;
      BFT_FREE(vmom);

      const char *name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;

      size_t l_name_width = cs_log_strlen(name);

      max_name_width = CS_MAX(max_name_width, l_name_width);

    } /* End of first loop on fields */

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
                    "  ** Computed fields on %s\n"
                    "     -------------------%.*s\n"),
                  loc_name, loc_name_w, _underline);

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
      size_t w1 = (col == 0) ? max_name_width : 14;
      for (i = 0; i < w0; i++)
        tmp_s[col][i] = '-';
      for (i = w0; i < w1; i++)
        tmp_s[col][i] = ' ';
      tmp_s[col][w1] = '\0';
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

      _dim = (f->interleaved && f->dim == 3) ? 4 : f->dim;

      const char *name = cs_field_get_key_str(f, label_key_id);
      if (name == NULL)
        name = f->name;

      double t_weight = -1;
      if (total_weight > 0 && (f->type | CS_FIELD_INTENSIVE))
        t_weight = total_weight;

      _log_array_info("v  ",
                      name,
                      max_name_width,
                      f->interleaved,
                      f->dim,
                      n_g_elts,
                      t_weight,
                      vmin + log_count,
                      vmax + log_count,
                      vsum + log_count,
                      wsum + log_count);

      log_count += _dim;

    } /* End of loop on fields */

  } /* End of loop on mesh locations */

  BFT_FREE(wsum);
  BFT_FREE(vsum);
  BFT_FREE(vmax);
  BFT_FREE(vmin);
  BFT_FREE(log_id);

  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*----------------------------------------------------------------------------
 * Main post-processing output of additional simple statistics
 *----------------------------------------------------------------------------*/

static void
_cs_log_sstats(void)
{
  int     stat_id;
  double _boundary_surf = -1;
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
                    loc_name, loc_name_w, _underline);

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
        size_t w1 = (col == 0) ? max_name_width : 14;
        for (i = 0; i < w0; i++)
          tmp_s[col][i] = '-';
        for (i = w0; i < w1; i++)
          tmp_s[col][i] = ' ';
        tmp_s[col][w1] = '\0';
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
                        true,
                        _sstats[stat_id].dim,
                        n_g_elts,
                        t_weight,
                        vmin + stat_id,
                        vmax + stat_id,
                        vsum + stat_id,
                        wsum + stat_id);

      } /* End of loop on stats */

    } /* End of loop on mesh locations */

    sstat_cat_start = sstat_cat_end;

  } /* End of loop on mesh categories */

  BFT_FREE(wsum);
  BFT_FREE(vsum);
  BFT_FREE(vmax);
  BFT_FREE(vmin);

  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize logging of moments
 *
 * Currently, an external cumulative time array is simply mapped to
 * the post-processing API.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_init_moments(const cs_real_t  *cumulative_time)
{
  _cumulative_mom_time = cumulative_time;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free arrays possible used by logging of array statistics.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_destroy_all(void)
{
  if (_name_map != NULL) {
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
    cs_map_name_to_id_destroy(&_name_map);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log field and other array statistics for the current time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration(void)
{
  _cs_log_fields();

  if (_n_sstats > 0)
    _cs_log_sstats();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add array not saved as permanent field to logging of fields.
 *
 * \param[in]  name         array name
 * \param[in]  category     category name
 * \param[in]  loc_id       associated mesh location id
 * \param[in]  is_intensive are the matching values intensive ?
 * \param[in]  dimension    associated dimension (interleaved)
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

  if (_name_map == NULL) {
    _name_map = cs_map_name_to_id_create();
    _category_map = cs_map_name_to_id_create();
  }

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
    case CS_MESH_LOCATION_BOUNDARY_FACES:
      weight = mq->b_face_surf;
      have_weight = true;
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

END_C_DECLS
