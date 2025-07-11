/*============================================================================
 * Management of temporal moments.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_parall.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"
#include "base/cs_prototypes.h"
#include "base/cs_time_step.h"
#include "base/cs_function.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_time_moment.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_time_moment.cpp

  \brief Temporal moments management.

  \enum cs_time_moment_type_t

  \brief Moment type

  \var CS_TIME_MOMENT_MEAN
       Moment is a mean
  \var CS_TIME_MOMENT_VARIANCE
       Moment is a variance

  \enum cs_time_moment_restart_t

  \brief Moment restart behavior.

  \var CS_TIME_MOMENT_RESTART_RESET
       Moment is reset in case of a restart file:
       starting time step will be no older than the restart time step.
  \var CS_TIME_MOMENT_RESTART_AUTO
       Moment uses restart information if available:
       if the requested time step is older than the restart time step,
       restart information will be used, though the starting time step
       is not guaranteed.
  \var CS_TIME_MOMENT_RESTART_EXACT
       If the requested time step is older than the restart time step,
       restart information will be used, and if the start time step
       or time does not match the one available, an error is thrown.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Moment time accumulator definition */
/*------------------------------------*/

typedef struct {

  int                     restart_id;   /* Matching id in restart info */

  int                     nt_start;     /* Associated starting time step;
                                           if < 0 (and f_id < 0), t_start is
                                           used directly, but nt_start is
                                           always set to a non-negative value
                                           once accumulation starts) */

  double                  t_start;      /* Associated starting time value
                                           (may be initialized to -1 if
                                           accumulation starts at nt_start,
                                           but is always set to a non-negative
                                           value once accumulation starts) */

  int                     location_id;  /* Associated mesh location id */

  cs_time_moment_data_t  *data_func;    /* Associated data value computation
                                           function (1 assumed if null) */
  const void             *data_input;   /* pointer to optional (untyped)
                                           value or structure */

  cs_real_t               val0;         /* Associated value if location_id
                                           is CS_MESH_LOCATION_NONE */
  cs_real_t              *val;          /* Pointer to associated values
                                           otherwise */

} cs_time_moment_wa_t;

/* Moment definitions */
/*--------------------*/

typedef struct {

  cs_time_moment_type_t   type;         /* Moment type */

  int                     restart_id;   /* Matching id in restart info */

  int                     wa_id;        /* Associated weight accumulator id */

  int                     f_id;         /* Associated field id, or -1 */

  int                     dim;          /* Associated field dimensions */
  int                     data_dim;     /* Associated data field dimensions */
  int                     location_id;  /* Associated mesh location id */

  cs_function_t          *eval_func;    /* Associated evaluated elements computation
                                           function, or null */
  cs_time_moment_data_t  *data_func;    /* Associated data elements computation
                                           function, or null */
  const void             *data_input;   /* pointer to optional (untyped)
                                           value or structure */

  int                     l_id;         /* Associated id of lower order moment
                                           (mean for variance), or -1 */

  char                   *name;         /* Associated name, if f_id < 0 */
  double                 *val;          /* Associated value, if f_id < 0 */

  int                     nt_cur;       /* Time step number of last update */

} cs_time_moment_t;

/* Moment restart metadata */
/*-------------------------*/

typedef struct {

  int                     nt_prev;        /* Restart time step */
  double                  t_prev;         /* Restart time */

  int                     n_wa;           /* Number of weight accumulators */
  int                     n_moments;      /* Number of moments */

  const char            **name;           /* Moment name */
  char                   *name_buf;       /* Buffer for names */

  int                    *wa_location_id; /* Weight accumulator location ids */
  int                    *wa_nt_start;    /* Weight accumulator start iters. */
  cs_real_t              *wa_t_start;     /* Weight accumulator start times */
  cs_real_t              *wa_val0;        /* Weight accumulator values for
                                             global (loc 0) accumulators */

  int                    *m_type;         /* Moment types */
  int                    *location_id;    /* Moment location */
  int                    *dimension;      /* Moment dimension */
  int                    *wa_id;          /* Associated accumulator ids */
  int                    *l_id;           /* Associated lower order ids */

} cs_time_moment_restart_info_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static int  _n_moment_sd_defs = 0;
static int  _n_moment_sd_defs_max = 0;

static int  _n_moment_wa = 0;
static int  _n_moment_wa_max = 0;

static int  _n_moments = 0;
static int  _n_moments_max = 0;

static int **_moment_sd_defs = nullptr;

static cs_time_moment_wa_t *_moment_wa = nullptr;
static cs_time_moment_t *_moment = nullptr;

static  bool _restart_info_checked = false;
static  bool _restart_uses_main = false;
static  cs_time_moment_restart_info_t *_restart_info = nullptr;

static double _t_prev_iter = 0.;

static const cs_real_t *_p_dt = nullptr; /* Mapped reference time step */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Names associated with moment types */

const char  *cs_time_moment_type_name[] = {N_("mean"),
                                           N_("variance")};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

int
cs_f_time_moment_field_id(int  m_num);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Abort in case expected restart read failed.
 *
 * parameters:
 *   retcode <-- previous return code for restart read operation
 *----------------------------------------------------------------------------*/

static void
_assert_restart_success(int retcode)
{
  if (retcode != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0,
              _("Error reading expected section in restart file."));
}

/*----------------------------------------------------------------------------
 * Read restart metadata.
 *
 * parameters:
 *   r <-- pointer to restart file
 *----------------------------------------------------------------------------*/

static void
_restart_info_read_auxiliary(cs_restart_t  *r)
{
  cs_lnum_t sizes[3];
  int retcode;

  const cs_time_step_t  *ts = cs_glob_time_step;

  retcode = cs_restart_read_section(r,
                                    "time_moments:sizes",
                                    CS_MESH_LOCATION_NONE,
                                    3,
                                    CS_TYPE_int,
                                    sizes);

  if (retcode == CS_RESTART_ERR_EXISTS)
    return;

  /* Now read main metadata */

  CS_MALLOC(_restart_info, 1, cs_time_moment_restart_info_t);

  cs_time_moment_restart_info_t  *ri = _restart_info;

  ri->nt_prev = ts->nt_prev;
  ri->t_prev = ts->t_prev;

  ri->n_wa = sizes[0];
  ri->n_moments = sizes[1];

  CS_MALLOC(ri->name, ri->n_moments, const char*);
  CS_MALLOC(ri->name_buf, sizes[2] + 1, char);

  retcode = cs_restart_read_section(r,
                                    "time_moments:names",
                                    CS_MESH_LOCATION_NONE,
                                    sizes[2],
                                    CS_TYPE_char,
                                    ri->name_buf);
  _assert_restart_success(retcode);

  ri->name[0] = ri->name_buf;
  for (int i = 0, j = 1; j < ri->n_moments; i++) {
    if (ri->name_buf[i] == '\0') {
      ri->name[j] = ri->name_buf + i + 1;
      j++;
    }
  }

  CS_MALLOC(ri->wa_location_id, ri->n_wa, int);
  CS_MALLOC(ri->wa_nt_start, ri->n_wa, int);
  CS_MALLOC(ri->wa_t_start, ri->n_wa, cs_real_t);
  ri->wa_val0 = nullptr;

  cs_restart_read_section(r,
                          "time_moments:wa:location_id",
                          CS_MESH_LOCATION_NONE,
                          ri->n_wa,
                          CS_TYPE_int,
                          ri->wa_location_id);
  _assert_restart_success(retcode);

  int n_val0 = 0;
  for (int i = 0; i < ri->n_wa; i++) {
    if (ri->wa_location_id[i] == CS_MESH_LOCATION_NONE)
      n_val0 += 1;
  }

  cs_restart_read_section(r,
                          "time_moments:wa:nt_start",
                          CS_MESH_LOCATION_NONE,
                          ri->n_wa,
                          CS_TYPE_int,
                          ri->wa_nt_start);
  _assert_restart_success(retcode);

  cs_restart_read_section(r,
                          "time_moments:wa:t_start",
                          CS_MESH_LOCATION_NONE,
                          ri->n_wa,
                          CS_TYPE_cs_real_t,
                          ri->wa_t_start);
  _assert_restart_success(retcode);

  if (n_val0 > 0) {
    CS_MALLOC_HD(ri->wa_val0, ri->n_wa, cs_real_t, cs_alloc_mode);
    cs_restart_read_section(r,
                            "time_moments:wa:val_g",
                            CS_MESH_LOCATION_NONE,
                            ri->n_wa,
                            CS_TYPE_cs_real_t,
                            ri->wa_val0);
    _assert_restart_success(retcode);
  }

  /* Information on moments proper */

  CS_MALLOC(ri->m_type, ri->n_moments, int);
  CS_MALLOC(ri->location_id, ri->n_moments, int);
  CS_MALLOC(ri->dimension, ri->n_moments, int);
  CS_MALLOC(ri->wa_id, ri->n_moments, int);
  CS_MALLOC(ri->l_id, ri->n_moments, int);

  retcode = cs_restart_read_section(r,
                                    "time_moments:type",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_int,
                                    ri->m_type);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "time_moments:location_id",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_int,
                                    ri->location_id);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "time_moments:dimension",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_int,
                                    ri->dimension);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "time_moments:wa_id",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_int,
                                    ri->wa_id);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "time_moments:lower_order_id",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_int,
                                    ri->l_id);
  _assert_restart_success(retcode);
}

/*----------------------------------------------------------------------------
 * Read restart metadata.
 *----------------------------------------------------------------------------*/

static void
_restart_info_read(void)
{
  const cs_time_step_t  *ts = cs_glob_time_step;

  if (ts->nt_prev < 1 && !cs_restart_present())
    return;

  cs_restart_t *r = nullptr;

  /* Read previous time step if not already done */

  if (ts->nt_prev < 1) {
    r = cs_restart_create("main.csc", "restart", CS_RESTART_MODE_READ);
    cs_restart_read_time_step_info(r);
    if (_restart_uses_main == false)
      cs_restart_destroy(&r);
  }

  /* Now read time-moment specific data */

  if (r == nullptr) {
    if (_restart_uses_main)
      r = cs_restart_create("main.csc", nullptr, CS_RESTART_MODE_READ);
    else
      r = cs_restart_create("auxiliary.csc", nullptr, CS_RESTART_MODE_READ);
  }

  _restart_info_read_auxiliary(r);

  cs_restart_destroy(&r);

  /* Now change checked status */

  _restart_info_checked = true;
}

/*----------------------------------------------------------------------------
 * Free restart metadata.
 *----------------------------------------------------------------------------*/

static void
_restart_info_free(void)
{
  cs_time_moment_restart_info_t  *ri = _restart_info;

  if (ri != nullptr) {

    CS_FREE(ri->l_id);
    CS_FREE(ri->wa_id);
    CS_FREE(ri->dimension);
    CS_FREE(ri->location_id);
    CS_FREE(ri->m_type);

    CS_FREE(ri->wa_val0);
    CS_FREE(ri->wa_t_start);
    CS_FREE(ri->wa_nt_start);
    CS_FREE(ri->wa_location_id);

    CS_FREE(ri->name_buf);
    CS_FREE(ri->name);

    CS_FREE(ri);

    _restart_info = ri;
  }
}

/*----------------------------------------------------------------------------
 * Check if a moment can use previous data.
 *
 * Depending on the restart mode, restart time and time step may also
 * be updated.
 *
 * parameters:
 *   name           <-- moment name
 *   ts             <-- time step status
 *   ri             <-> resource info
 *   location_id    <-- associated mesh location id
 *   wa_location_id <-> associated weigh accumulator mesh location id
 *   dim            <-- dimension associated with moment
 *   type           <-- moment type
 *   nt_start       <-> starting time step
 *   t_start        <-> starting time
 *   restart_mode   <-- behavior in case of restart (reset, automatic, strict)
 *   restart_name   <-- if non-null, previous name in case of restart
 *
 * returns:
 *   id of matching restart moment id, or -1 if none matches
 *----------------------------------------------------------------------------*/

static int
_check_restart(const char                     *name,
               const cs_time_step_t           *ts,
               cs_time_moment_restart_info_t  *ri,
               int                             location_id,
               int                            &wa_location_id,
               int                             dim,
               cs_time_moment_type_t           type,
               int                            *nt_start,
               double                         *t_start,
               cs_time_moment_restart_t        restart_mode,
               const char                     *restart_name)
{
  int i;
  int prev_id = -1;
  int prev_wa_id = -1;

  if (   (*nt_start > -1 && *nt_start >= ri->nt_prev)
      || (*t_start >= 0 && *t_start >= ri->t_prev))
    return prev_id;

  /* Adjust accumulator info if moment should be restarted */

  if (restart_mode == CS_TIME_MOMENT_RESTART_RESET) {
    *nt_start = ri->nt_prev + 1;
    *t_start = ri->t_prev;
    return prev_id;
  }

  /* If we reach here, restart info is required */

  /* Find matching restart data if required, and adjust accumulator
     info if moment should be restarted, or if available data does
     not match and we do not require exact mode. */

  const char *_r_name = (restart_name != nullptr) ? restart_name : name;
  for (i = 0; i < ri->n_moments; i++) {
    if (strcmp(ri->name[i], _r_name) == 0) {
      bool matching_restart = true;
      prev_id = i;
      prev_wa_id = ri->wa_id[i];
      if (   (   ri->wa_location_id[prev_wa_id] != wa_location_id
              && wa_location_id != 0)
          || ri->m_type[i] != (int)type
          || ri->location_id[i] != location_id
          || ri->dimension[i] != dim)
        matching_restart = false;
      if (   restart_mode == CS_TIME_MOMENT_RESTART_EXACT
           && (   ri->wa_nt_start[prev_wa_id] != *nt_start
               || (   !ts->is_local
                   && fabs(ri->wa_t_start[prev_wa_id] - *t_start) > 1.e-18)))
        matching_restart = false;
      if (matching_restart == false) {
        bft_printf(_("\nRestart data for time moment \"%s\"\n"
                     " (previously \"%s\") does not match.\n"
                     "  previous values:\n"
                     "    weight accumulator location_id: %d\n"
                     "    type:                           %d\n"
                     "    location_id:                    %d\n"
                     "    dimension:                      %d\n"
                     "    start time step:                %d\n"
                     "    start time:                     %12.5e\n"),
                   name, _r_name, ri->wa_location_id[prev_wa_id],
                   ri->m_type[i], ri->location_id[i], ri->dimension[i],
                   ri->wa_nt_start[prev_wa_id],
                   ri->wa_t_start[prev_wa_id]);
        if (restart_mode == CS_TIME_MOMENT_RESTART_AUTO) {
          bft_printf
            (_("\nWarning: computation of time moment \"%s\""
               " will be reset,\n"
               "         as restart data for \"%s\" does not match.\n"),
             name, _r_name);
          *nt_start = ri->nt_prev + 1;
          *t_start = ri->t_prev;
        }
        else if (restart_mode == CS_TIME_MOMENT_RESTART_EXACT)
          bft_error(__FILE__, __LINE__, 0,
                    _("Restart data for time moment \"%s\"\n"
                      " (previously \"%s\") does not match."),
                    name, _r_name);
      }
      if (matching_restart == false)
        prev_id = -1;
      else {
        // wa_location_id > 0 only for compatibility with older checkpoints.
        wa_location_id = ri->wa_location_id[prev_wa_id];
        *nt_start = ri->wa_nt_start[prev_wa_id];
        *t_start = ri->wa_t_start[prev_wa_id];
      }
      break;
    }
  }

  if (i >= ri->n_moments) {
    if (restart_mode == CS_TIME_MOMENT_RESTART_AUTO) {
      bft_printf
        (_("\nWarning: computation of time moment \"%s\""
           "will be reset,\n"
           "           as restart data for \"%s\" is not available.\n"),
         name, _r_name);
      prev_id = -1;
      *nt_start = ri->nt_prev + 1;
      *t_start = ri->t_prev;
    }
    else if (restart_mode == CS_TIME_MOMENT_RESTART_EXACT)
      bft_error(__FILE__, __LINE__, 0,
                _("Restart data for time moment \"%s\"\n"
                  " (previously \"%s\") not available."),
                name, _r_name);

  }

  /* Also check for presence of sub-moment restart info in case of
     higer order restart */

  if (prev_id > -1) {

    for (int m_type = (cs_time_moment_type_t)type;
         m_type > CS_TIME_MOMENT_MEAN;
         m_type--) {

      cs_time_moment_type_t s_type = (cs_time_moment_type_t)(m_type -1);
      int l_dim = (dim == 6 && m_type == CS_TIME_MOMENT_VARIANCE) ? 3 : dim;

      int l_id = ri->l_id[prev_id];

      if (   ri->wa_id[l_id] != prev_wa_id
             || ri->m_type[l_id] != (int)s_type
          || ri->location_id[l_id] != location_id
          || ri->dimension[l_id] != l_dim)
        bft_error(__FILE__, __LINE__, 0,
                  _("Restart data for time moment \"%s\"\n"
                    " (previously \"%s\") seems inconsistent:\n"
                    "   lower order moment of type %s was \"%s\",\n"
                    "   but has non-matching attributes:\n"
                    "    weight accumulator id: %d (expected %d)\n"
                    "    type:                  %d\n"
                    "    location_id:           %d\n"
                    "    dimension:             %d\n"),
                  name, _r_name, cs_time_moment_type_name[s_type],
                  ri->name[l_id], ri->wa_id[l_id], prev_wa_id,
                  ri->m_type[l_id], ri->location_id[l_id], ri->dimension[l_id]);

    }

  }

  /* Return previous moment id */

  return prev_id;
}

/*----------------------------------------------------------------------------
 * Build simple data description string (for error messages)
 *
 * parameters:
 *   n_fields  <-- number of multiplying fields
 *   field_id  <-- array of ids of multiplying fields
 *   comp_id   <-- array of ids of multiplying components
 *   desc_size <-- string size
 *   desc      --> descritpion string
 *----------------------------------------------------------------------------*/

static void
_build_sd_desc(int        n_fields,
               const int  f_id[],
               const int  c_id[],
               size_t     desc_size,
               char       desc[])
{
  char  *s = desc;
  size_t c_size = 0;

  assert(desc_size > 4);

  for (int i = 0; i < n_fields; i++) {
    size_t r_size = desc_size - c_size;
    if (r_size > 4) {
      s = desc + c_size;
      snprintf(s, r_size - 4, "(%d, %d)", f_id[i], c_id[i]);
    }
    else {
      s = desc + desc_size - 4;
      snprintf(s, 4, "...");
    }
    desc[desc_size - 1] = '\0';
    c_size = strlen(desc);
  }
}

/*----------------------------------------------------------------------------
 * Add or find simple data definition.
 *
 * Simple data is defined by an array of integers of size:
 *   3 + (2+dim)*n_fields
 * It contains, in succession: location_id, field_dimension, n_fields,
 * then (field_id, component_id, component_id_0, component_id_dim-1) tuples.
 * Negative component ids in the second position mean all components are used.
 *
 * parameters:
 *   name         <-- name of associated decription
 *   n_fields     <-- number of multiplying fields
 *   field_id     <-- array of ids of multiplying fields
 *   comp_id      <-- array of ids of multiplying components
 *   is_intensive <-> is resulting data intensive ?
 *
 * returns:
 *   id of matching simple data definition
 *----------------------------------------------------------------------------*/

static int
_find_or_add_sd(const char  *name,
                int          n_fields,
                const int    f_id[],
                const int    c_id[],
                bool        *is_intensive)
{
  char sd_desc[256];

  int sd_id = -1;

  if (n_fields < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Definition of simple data requires at least one field id."));

  /* Check if all associated fields are intensive */

  for (int i = 0; i < n_fields; i++) {
    const cs_field_t *f = cs_field_by_id(f_id[i]);
    *is_intensive = (*is_intensive) && (f->type & CS_FIELD_INTENSIVE);
  }

  /* Check if this definition has already been provided (assume field and
     component ids are given in same order; at worse, if this is not the case
     some data which could be shared will be duplicated, leading to slightly
     higher memory usage and computational cost)

     Also check if the time moment is a combination of intensive fields
     and therefore is an intensive field
     */

  for (sd_id = 0; sd_id < _n_moment_sd_defs; sd_id++) {
    bool is_different = false;
    const int *msd = _moment_sd_defs[sd_id];
    const int stride = 2 + msd[1];
    if (n_fields != msd[2])
      is_different = true;
    else {
      for (int i = 0; i < n_fields; i++) {
        const cs_field_t *f = cs_field_by_id(f_id[i]);
        const int _c_id = (f->dim > 1) ? c_id[i] : 0;
        if (   msd[3 + i*stride] != f_id[i]
            || msd[3 + i*stride+1] != _c_id)
          is_different = true;
      }
    }
    if (! is_different)
      return sd_id;
  }

  /* If we did not return yet, a new structure must be added */

  /* Reallocate if necessary */

  if (_n_moment_sd_defs + 1 > _n_moment_sd_defs_max) {
    if (_n_moment_sd_defs_max < 1)
      _n_moment_sd_defs_max = 2;
    else
      _n_moment_sd_defs_max *= 2;
    CS_REALLOC(_moment_sd_defs,
               _n_moment_sd_defs_max,
               int *);
  }

  sd_id = _n_moment_sd_defs;
  _n_moment_sd_defs += 1;

  /* Determine location and dimension */

  int location_id = CS_MESH_LOCATION_NONE;
  int dim = 1;

  for (int i = 0; i < n_fields; i++) {
    const cs_field_t *f = cs_field_by_id(f_id[i]);
    if (location_id != f->location_id) {
      if (location_id != CS_MESH_LOCATION_NONE) {
        _build_sd_desc(n_fields, f_id, c_id, 256, sd_desc);
        bft_error
          (__FILE__, __LINE__, 0,
           _("Definition of simple data used for %s:\n"
             "%s\n"
             "mixes fields with location id %d and location id %d."),
           name, sd_desc, location_id, f->location_id);
      }
      else
        location_id = f->location_id;
    }
    if (c_id[i] < 0) { /* All components */
      if (f->dim != 1  && f->dim != 3 && f->dim != 6 && f->dim != 9) {
        _build_sd_desc(n_fields, f_id, c_id, 256, sd_desc);
        bft_error
          (__FILE__, __LINE__, 0,
           _("Definition of simple data used for %s:\n"
             "%s\n"
             "includes field of dimension different from 1, 3, 6, or 9.\n"
             "The definition must be split."),
           name, sd_desc);
      }
      if (dim == 3 && f->dim == 3)
        dim = 6;
      else
        dim *= f->dim;
      if (dim > 9) {
        _build_sd_desc(n_fields, f_id, c_id, 256, sd_desc);
        bft_error
          (__FILE__, __LINE__, 0,
           _("Definition of simple data used for %s:\n"
             "%s\n"
             "leads to a field of dimension > 9.\n"
             "The definition must be split."),
           name, sd_desc);
      }
    }
    else if (c_id[i] >= f->dim) {
      _build_sd_desc(n_fields, f_id, c_id, 256, sd_desc);
      bft_error
        (__FILE__, __LINE__, 0,
         _("Definition of simple data used for %s:\n"
           "%s\n"
           "includes a component id incompatible with field dimension."),
         name, sd_desc);
    }
  }

  /* Now initialize members */

  int stride = 2 + dim;
  int cur_dim = 1;
  int *msd;

  CS_MALLOC(msd, 3 + n_fields*stride, int);

  _moment_sd_defs[sd_id] = msd;

  msd[0] = location_id;
  msd[1] = dim;
  msd[2] = n_fields;

  for (int i = 0; i < n_fields; i++) {

    const cs_field_t *f = cs_field_by_id(f_id[i]);
    const int _c_id = (f->dim > 1) ? c_id[i] : 0;

    msd[3 + i*stride]     = f_id[i];
    msd[3 + i*stride + 1] = _c_id;

    if (_c_id > -1) {
      for (int j = 0; j < dim; j++)
        msd[3 + i*stride + 2 + j] = _c_id;
    }
    else if (f->dim == dim) {
      assert(cur_dim == 1);
      for (int j = 0; j < dim; j++)
        msd[3 + i*stride + 2 + j] = j;
      cur_dim = dim;
    }
    else {
      assert(dim == 6);
      assert(f->dim == 3);
      msd[3 + i*stride + 2 + 0] = 0;
      msd[3 + i*stride + 2 + 1] = 1;
      msd[3 + i*stride + 2 + 2] = 2;
      if (cur_dim == 1) {
        msd[3 + i*stride + 2 + 3] = 0;
        msd[3 + i*stride + 2 + 4] = 1;
        msd[3 + i*stride + 2 + 5] = 0;
        cur_dim = 3;
      }
      else {
        msd[3 + i*stride + 2 + 3] = 1;
        msd[3 + i*stride + 2 + 4] = 2;
        msd[3 + i*stride + 2 + 5] = 2;
        cur_dim = 6;
      }
    }

  }

  /* Structure is now initialized */

  return sd_id;
}

/*----------------------------------------------------------------------------
 * Free all moment weight and time accumulators
 *----------------------------------------------------------------------------*/

static void
_free_all_sd_defs(void)
{
  int i;

  for (i = 0; i < _n_moment_sd_defs; i++)
    CS_FREE(_moment_sd_defs[i]);

  CS_FREE(_moment_sd_defs);

  _n_moment_sd_defs = 0;
  _n_moment_sd_defs_max = 0;
}

/*----------------------------------------------------------------------------
 * Function pointer for computation of data values for moments computation
 * using simple data (field component product) definitions
 *
 * Simple data is defined by an array of integers of size:
 *   3 + 2*n_fields
 * It contains, in succession: location_id, field_dimension, n_fields,
 * the (field_id, component_id) couples. Negative component ids
 * mean all components are used.
 *
 * parameters:
 *   input <-- pointer to simple data array
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_sd_moment_data(const void  *input,
                cs_real_t   *vals)
{
  const int *msd = reinterpret_cast<const int *>(input);

  const cs_lnum_t dim = msd[1];
  const int stride = 2 + dim;
  const int location_id = msd[0];
  const int n_fields = msd[2];

  const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];

  int _f_dim[16*2];
  int *f_dim;

  const cs_real_t * _f_val[16];
  const cs_real_t **f_val;

  if (n_fields*2 > 16*3)
    CS_MALLOC(f_dim, n_fields*2, int);
  else
    f_dim = _f_dim;
  if (n_fields > 16)
    CS_MALLOC(f_val, n_fields, const cs_real_t *);
  else
    f_val = _f_val;

  for (int i = 0; i < n_fields; i++) {

    const cs_field_t *f = cs_field_by_id(msd[3 + stride*i]);

    f_val[i] = (const cs_real_t *)f->val;

    /* Field access strides for consistent access method */

    if (f->location_id != 0) {
      f_dim[i*2]     = f->dim;
      f_dim[i*2 + 1] = 1;
    }
    else {
      f_dim[i*2]     = 0;
      f_dim[i*2 + 1] = 1;
    }

  }

  /* Now compute values */

  for (cs_lnum_t  i = 0; i < n_elts; i++) {
    const cs_real_t *restrict v = f_val[0];
    cs_lnum_t m0 = f_dim[0];
    cs_lnum_t m1 = f_dim[1];
    for (cs_lnum_t k = 0; k < dim; k++) {
      cs_lnum_t c_id = msd[3 + 2 + k]; /* as below, with j = 0 */
      vals[i*dim + k] = v[m0*i + m1*c_id];
    }
    for (int j = 1; j < n_fields; j++) {
      v = f_val[j];
      m0 = f_dim[j*2];
      m1 = f_dim[j*2 + 1];
      for (cs_lnum_t k = 0; k < dim; k++) {
        cs_lnum_t c_id = msd[3 + j*stride + 2 + k];
        vals[i*dim + k] *= v[m0*i + m1*c_id];
      }
    }
  }

  /* Free temporary memory */

  if (f_dim != _f_dim)
    CS_FREE(f_dim);
  if (f_val != _f_val)
    CS_FREE(f_val);
}

/*----------------------------------------------------------------------------
 * Add or find moment weight and time accumulator.
 *
 * If no data function is provided, a constant weight of 1 is assumed
 * (this weight will be multiplied by the time step).
 *
 * Note that if the data_input associated with a data_func pointer is not
 * null, the lifecycle of the data pointed to must be handled separately
 * (and the pointer must remain valid throughout the time moment updates).
 *
 * parameters:
 *   data_func   <-- function used to define data values
 *   data_input  <-- pointer to optional (untyped) value or structure
 *   location_id <-- associated mesh location id
 *   nt_start    <-- starting time step
 *   t_start     <-- starting time
 *   prev_wa_id  <-- previous weight accumulator id, or -1
 *
 * returns:
 *   id of matching time accumulator
 *----------------------------------------------------------------------------*/

static int
_find_or_add_wa(cs_time_moment_data_t  *data_func,
                const void             *data_input,
                int                     location_id,
                int                     nt_start,
                double                  t_start,
                int                     prev_wa_id)
{
  int wa_id = -1;
  int _nt_start = nt_start;
  double _t_start = t_start;

  cs_time_moment_wa_t *mwa = nullptr;

  /* Reduce number of possible options */

  if (_nt_start < 0)
    _nt_start = -1;

  if (_t_start < 0. && _nt_start < 0)
    _nt_start = 0;

  if (nt_start >= 0)
    _t_start = -1.;

  /* Check if this accumulator is already defined */

  for (int i = 0; i < _n_moment_wa; i++) {
    mwa = _moment_wa + i;
    if (   nt_start == mwa->nt_start && fabs(mwa->t_start - _t_start) < 1e-18
        && data_func == mwa->data_func && data_input == mwa->data_input
        && prev_wa_id == mwa->restart_id)
      return i;
  }

  /* If we did not return yet, a new structure must be added */

  /* Reallocate if necessary */

  if (_n_moment_wa + 1 > _n_moment_wa_max) {
    if (_n_moment_wa_max < 1)
      _n_moment_wa_max = 2;
    else
      _n_moment_wa_max *= 2;
    CS_REALLOC_HD(_moment_wa, _n_moment_wa_max, cs_time_moment_wa_t, cs_alloc_mode);
  }

  /* Now initialize members */

  wa_id = _n_moment_wa;

  mwa = _moment_wa + _n_moment_wa;

  _n_moment_wa += 1;

  mwa->restart_id = prev_wa_id;

  mwa->nt_start = _nt_start;
  mwa->t_start = _t_start;

  mwa->location_id = location_id;

  mwa->data_func = data_func;
  mwa->data_input = data_input;

  mwa->val0 = 0;
  if (prev_wa_id > -1 && mwa->location_id == CS_MESH_LOCATION_NONE)
    mwa->val0 = _restart_info->wa_val0[prev_wa_id];

  mwa->val = nullptr;

  /* Structure is now initialized */

  return wa_id;
}

/*----------------------------------------------------------------------------
 * Free all moment weight and time accumulators
 *----------------------------------------------------------------------------*/

static void
_free_all_wa(void)
{
  int i;

  for (i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    CS_FREE(mwa->val);
  }

  CS_FREE(_moment_wa);

  _n_moment_wa = 0;
  _n_moment_wa_max = 0;
}

/*----------------------------------------------------------------------------
 * Add or find moment structure.                 .
 *
 * parameters:
 *   name        <-- name of associated moment, or null
 *   location_id <-- id of associated mesh location
 *   dim         <-- dimension associated with element data
 *   eval_func   <-- function used to define data values
 *   data_func   <-- function used to define data values
 *   data_input  <-- pointer to optional (untyped) value or structure
 *                   to be used by data_func
 *   type        <-- moment type
 *   wa_id       <-- weight accumulator id
 *   prev_id     <-- restart moment id
 *
 * return:
 *   id of matching moment
 *----------------------------------------------------------------------------*/

static int
_find_or_add_moment(int                     location_id,
                    int                     dim,
                    cs_function_t          *eval_func,
                    cs_time_moment_data_t  *data_func,
                    const void             *data_input,
                    cs_time_moment_type_t   type,
                    int                     wa_id,
                    int                     prev_id)
{
  cs_time_moment_t *mt = nullptr;
  int moment_id = -1;

  /* Check if this moment is already defined;
     ignore associated field at this stage, as a moment defined automatically
     to satisfy a dependency (i.e. a mean for a variance) may not have an
     associated field. */

  for (int i = 0; i < _n_moments; i++) {
    mt = _moment + i;

    if (   location_id == mt->location_id && dim == mt->data_dim
        && eval_func == mt->eval_func
        && data_func == mt->data_func && data_input == mt->data_input
        && type == mt->type && wa_id == mt->wa_id
        && prev_id == mt->restart_id)
      return i;
  }

  /* If we did not return yet, a new structure must be added */

  /* Reallocate if necessary */

  if (_n_moments + 1 > _n_moments_max) {
    if (_n_moments_max < 1)
      _n_moments_max = 2;
    else
      _n_moments_max *= 2;
    CS_REALLOC(_moment, _n_moments_max, cs_time_moment_t);
  }

  /* Now define moment */

  moment_id = _n_moments;
  _n_moments += 1;

  mt = _moment + moment_id;

  mt->type = type;
  mt->restart_id = prev_id;
  mt->wa_id = wa_id;
  mt->f_id = -1;

  mt->dim = (dim == 3 && type == CS_TIME_MOMENT_VARIANCE) ? 6 : dim;
  mt->data_dim = dim;
  mt->location_id = location_id;

  mt->eval_func = eval_func;
  mt->data_func = data_func;
  mt->data_input = data_input;

  mt->l_id = -1;

  mt->name = nullptr;
  mt->val = nullptr;

  mt->nt_cur = -1;

  return moment_id;
}

/*----------------------------------------------------------------------------
 * Free all moments
 *----------------------------------------------------------------------------*/

static void
_free_all_moments(void)
{
  int i;

  for (i = 0; i < _n_moments; i++) {
    cs_time_moment_t *mt = _moment + i;
    CS_FREE(mt->name);
    CS_FREE(mt->val);
  }

  CS_FREE(_moment);

  _n_moments = 0;
  _n_moments_max = 0;
}

/*----------------------------------------------------------------------------
 * Compute current weight.
 *
 * This function either returns a pointer to an allocated array,
 * or to w0. If the returned value is different from w0 (i.e. allocated),
 * the caller is responsible for freeing it.
 *
 * parameters:
 *   mwa  <-- moment weight accumulator
 *   dt   <-- cell time step values
 *   w0   <-- pointer to buffer in case weight values is of size 1
 *
 * returns:
 *   pointer to weight array (w0 or allocated array)
 *----------------------------------------------------------------------------*/

static cs_real_t *
_compute_current_weight(cs_time_moment_wa_t  *mwa,
                        cs_real_t             dt,
                        cs_real_t             w0[1])
{
  cs_lnum_t  n_w_elts;
  cs_real_t *w;

  const cs_time_step_t  *ts = cs_glob_time_step;

  assert(mwa->nt_start <= ts->nt_cur);

  if (mwa->location_id == CS_MESH_LOCATION_NONE) {
    n_w_elts = 1;
    w = w0;
  }
  else {
    n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
    CS_MALLOC_HD(w, n_w_elts, cs_real_t, cs_alloc_mode);
  }

  /* global time step multiplier */

  double _dt = dt;;
  if (mwa->nt_start == ts->nt_cur)
    _dt = ts->t_cur - mwa->t_start;

  /* Base weight */

  if (mwa->data_func != nullptr) {
    mwa->data_func(mwa->data_input, w);
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      w[i] *= _dt;
  }
  else {
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      w[i] = _dt;
  }

  return w;
}

/*----------------------------------------------------------------------------
 * Initialize weight accumulator if required
 *
 * parameters:
 *   mwa <-- moment weight accumulator
 *----------------------------------------------------------------------------*/

static void
_ensure_init_weight_accumulator(cs_time_moment_wa_t  *mwa)
{
  if (mwa->location_id != CS_MESH_LOCATION_NONE && mwa->val == nullptr) {
    cs_lnum_t n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
    CS_MALLOC_HD(mwa->val, n_w_elts, cs_real_t, cs_alloc_mode);
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      mwa->val[i] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Reset weight accumulator
 *
 * parameters:
 *   mwa <-- moment weight accumulator
 *----------------------------------------------------------------------------*/

static void
_reset_weight_accumulator(cs_time_moment_wa_t  *mwa)
{
  if (mwa->location_id == CS_MESH_LOCATION_NONE)
    mwa->val0 = 0.;
  else if (mwa->val != nullptr) { /* null possible if not active yet */
    cs_lnum_t n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      mwa->val[i] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Update weight accumulator
 *
 * parameters:
 *   mwa <-- moment weight accumulator
 *   w   <-- pointer to current weight values
 *----------------------------------------------------------------------------*/

static void
_update_weight_accumulator(cs_time_moment_wa_t  *mwa,
                           cs_real_t            *restrict w)
{
  if (mwa->location_id == CS_MESH_LOCATION_NONE)
    mwa->val0 += w[0];
  else {
    cs_lnum_t n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      mwa->val[i] += w[i];
  }
}

/*----------------------------------------------------------------------------
 * Initialize moment value if required
 *
 * parameters:
 *   mwa <-- moment
 *----------------------------------------------------------------------------*/

static void
_ensure_init_moment(cs_time_moment_t  *mt)
{
  if (mt->f_id < 0 && mt->val == nullptr) {
    cs_lnum_t n_elts = cs_mesh_location_get_n_elts(mt->location_id)[2];
    cs_lnum_t n_d_elts = n_elts*(cs_lnum_t)(mt->dim);
    CS_MALLOC_HD(mt->val, n_d_elts, cs_real_t, cs_alloc_mode);
    for (cs_lnum_t i = 0; i < n_d_elts; i++)
      mt->val[i] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a moment whose data values will be computed using a
 * specified function.
 *
 * \param[in]  name           name of associated moment
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  dim            dimension associated with element data
 * \param[in]  is_intensive   is the time moment intensive?
 * \param[in]  data_func      function used to define data values
 * \param[in]  data_input     pointer to optional (untyped) value or structure
 *                            to be used by data_func
 * \param[in]  w_data_func    function used to define weight values
 * \param[in]  w_data_input   pointer to optional (untyped) value or structure
 *                            to be used by w_data_func
 * \param[in]  type           moment type
 * \param[in]  nt_start       starting time step (or -1 to use t_start)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 * \param[in]  restart_name   if non-null, previous name in case of restart
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

static int
_time_moment_define_by_func(const char                *name,
                            int                        location_id,
                            int                        dim,
                            bool                       is_intensive,
                            cs_function_t             *eval_func,
                            cs_time_moment_data_t     *data_func,
                            const void                *data_input,
                            cs_time_moment_data_t     *w_data_func,
                            void                      *w_data_input,
                            cs_time_moment_type_t      type,
                            int                        nt_start,
                            double                     t_start,
                            cs_time_moment_restart_t   restart_mode,
                            const char                *restart_name)
{
  int i;
  cs_field_t  *f;

  int wa_location_id = 0; // > 0 only if restarted from older checkpoint.

  cs_time_moment_t *mt = nullptr;

  int moment_dim = (   dim == 3
                    && type == CS_TIME_MOMENT_VARIANCE) ? 6 : dim;
  int moment_id = -1;
  int prev_id = -1, prev_wa_id = -1;
  int _nt_start = nt_start;
  double _t_start = t_start;

  const cs_time_step_t  *ts = cs_glob_time_step;

  if (w_data_func == nullptr && ts->is_local == 0)
    wa_location_id = 0;

  /* If this is the first moment to be defined, ensure
     restart data is read if available */

  if (_restart_info_checked == false)
    _restart_info_read();

  /* Find matching restart data if required, and adjust accumulator
     info if moment should be restarted, or if available data does
     not match and we do not require exact mode. */

  if (_restart_info != nullptr) {
    cs_time_moment_restart_info_t  *ri = _restart_info;
    prev_id = _check_restart(name,
                             ts,
                             ri,
                             location_id,
                             wa_location_id,
                             moment_dim,
                             type,
                             &_nt_start,
                             &_t_start,
                             restart_mode,
                             restart_name);
    if (prev_id > -1)
      prev_wa_id = ri->wa_id[prev_id];
  }

  if (_nt_start < 0 && _t_start < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Time moment definition for \"%s\" is inconsistent:\n"
                " either starting time step or physical time must be >= 0."),
              name);

  /* Find or define matching weight accumulator info
     (if the weight value is constant,
     do not assign a mesh location to it as a single value is enough) */

  const int wa_id = _find_or_add_wa(w_data_func,
                                    w_data_input,
                                    wa_location_id,
                                    _nt_start,
                                    _t_start,
                                    prev_wa_id);

  /* Check for possible previous definition */

  f = cs_field_by_name_try(name);
  if (f != nullptr) {
    for (i = 0; i < _n_moments; i++) {
      mt = _moment + i;
      if (mt->f_id == f->id) {
        moment_id = i;
        break;
      }
    }
  }
  else { /* Build field matching moment */
    int type_flag = CS_FIELD_POSTPROCESS | CS_FIELD_ACCUMULATOR;
    if (is_intensive)
      type_flag |= CS_FIELD_INTENSIVE;

    f = cs_field_create(name,
                        type_flag,
                        location_id,
                        moment_dim,
                        false);  /* no previous values */
  }

  moment_id = _find_or_add_moment(location_id,
                                  dim,
                                  eval_func,
                                  data_func,
                                  data_input,
                                  type,
                                  wa_id,
                                  prev_id);

  /* Now define moment */

  mt = _moment + moment_id;

  mt->f_id = f->id;
  CS_FREE(mt->name); /* in case previously defined as sub-moment */

  /* Define sub moments */

  for (int m_type = type;
       m_type > CS_TIME_MOMENT_MEAN;
       m_type--) {

    const cs_time_moment_restart_info_t  *ri = _restart_info;

    int s_prev_id = (ri != nullptr && prev_id != -1) ? ri->l_id[prev_id] : prev_id;

    cs_time_moment_type_t s_type = (cs_time_moment_type_t)(m_type -1);

    int l_id = _find_or_add_moment(location_id,
                                   dim,
                                   eval_func,
                                   data_func,
                                   data_input,
                                   s_type,
                                   wa_id,
                                   s_prev_id);

    /* Redefine the parent moment as pointer _moment might have
     * been reallocated during the previous _find_or_add_moment call */
    mt = _moment + moment_id;
    mt->l_id = l_id;

    /* Now defining sub moment */
    mt = _moment + l_id;

    if (mt->f_id < 0) {
      char s[64];
      snprintf(s, 64, "<auto_%s_moment_%d>",
               cs_time_moment_type_name[mt->type], l_id);
      s[63] = '\0';
      CS_MALLOC(mt->name, strlen(s)+1, char);
      strcpy(mt->name, s);
    }

  }

  return moment_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset selected time step for starting time step of selected moment.
 *
 * All other time moments sharing the same start time will also start
 * at the same time step.
 *
 * \param[in]   moment_id  id of associated moment, or -1 for all
 * \param[in]   nt_start   starting time step (if < 0, t_start is used)
 * \param[in]   t_start    starting time value
 */
/*----------------------------------------------------------------------------*/

static void
_time_moment_set_start_time(int     moment_id,
                            int     nt_start,
                            double  t_start)
{
  const cs_time_step_t  *ts = cs_glob_time_step;

  int *nt_start_prev;
  CS_MALLOC(nt_start_prev, _n_moment_wa, int);

  for (int i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    nt_start_prev[i] = mwa->nt_start;
  }

  int s_id = moment_id;
  int e_id = moment_id+1;
  if (moment_id == -1) {
    s_id = 0;
    e_id = _n_moments;
  }

  for (int i = s_id; i < e_id; i++) {
    cs_time_moment_t *mt = _moment + i;
    cs_time_moment_wa_t *mwa = _moment_wa + mt->wa_id;
    if (nt_start > -1)
      mwa->nt_start = nt_start;
    if (t_start > -1)
      mwa->t_start = t_start;
  }

  for (int i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    if (   nt_start_prev[i] != mwa->nt_start
        && mwa->nt_start < ts->nt_cur) {
      _reset_weight_accumulator(mwa);
      mwa->t_start = ts->t_cur;
    }
    else if (t_start > -1 && mwa->nt_start > ts->nt_cur)
      mwa->nt_start = -1;
  }

  for (int i = 0; i < _n_moments; i++) {
    cs_time_moment_t *mt = _moment + i;
    cs_time_moment_wa_t *mwa = _moment_wa + mt->wa_id;
    if (   nt_start_prev[i] != mwa->nt_start
        && mwa->nt_start < ts->nt_cur) {
      cs_lnum_t n_elts = cs_mesh_location_get_n_elts(mt->location_id)[0];
      cs_lnum_t nd = n_elts * mt->dim;

      /* Reset values */
      cs_real_t *restrict val = mt->val;
      if (mt->f_id > -1) {
        cs_field_t *f = cs_field_by_id(mt->f_id);
        val = f->val;
      }
      cs_array_real_fill_zero(nd, val);
    }
  }

  CS_FREE(nt_start_prev);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return field id associated with a given moment
 *
 * parameters:
 *   m_num <-- moment number (1 to n)
 *
 * returns:
 *   field id, or -1
 *----------------------------------------------------------------------------*/

int
cs_f_time_moment_field_id(int m_num)

{
  int retval = -1;

  const cs_field_t *f = cs_time_moment_get_field(m_num - 1);
  if (f != nullptr)
    retval = f->id;

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all moments management metadata.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_destroy_all(void)
{
  _free_all_moments();
  _free_all_wa();
  _free_all_sd_defs();

  _p_dt = nullptr;
  _restart_info_checked = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a moment of a product of existing fields components.
 *
 * Moments will involve the tensor products of their component fields,
 * and only scalar, vector, or rank-2 tensors are handled (for
 * post-processing output reasons), so a moment may not involve more than
 * 2 vectors or 1 tensor, unless single components are specified.
 *
 * \param[in]  name           name of associated moment
 * \param[in]  n_fields       number of associated fields
 * \param[in]  field_id       ids of associated fields
 * \param[in]  component_id   ids of matching field components (-1 for all)
 * \param[in]  type           moment type
 * \param[in]  nt_start       starting time step (or -1 to use t_start)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 * \param[in]  restart_name   if non-null, previous name in case of restart
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_define_by_field_ids(const char                *name,
                                   int                        n_fields,
                                   const int                  field_id[],
                                   const int                  component_id[],
                                   cs_time_moment_type_t      type,
                                   int                        nt_start,
                                   double                     t_start,
                                   cs_time_moment_restart_t   restart_mode,
                                   const char                *restart_name)
{
  int m_id = -1;
  bool is_intensive = true;
  int sd_id =_find_or_add_sd(name, n_fields, field_id, component_id,
                             &is_intensive);

  const int *msd = _moment_sd_defs[sd_id];

  m_id = _time_moment_define_by_func(name,
                                     msd[0],
                                     msd[1],
                                     is_intensive,
                                     nullptr,
                                     _sd_moment_data,
                                     msd,
                                     nullptr,
                                     nullptr,
                                     type,
                                     nt_start,
                                     t_start,
                                     restart_mode,
                                     restart_name);

  return m_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time moment of an existing field.
 *
 * Moments will involve the tensor products of their component fields,
 * and only scalar, vector, or rank-2 tensors are handled (for
 * post-processing output reasons), so a 1st-order moment (i.e. mean) may
 * involve a scalar, vector, or tensor, while a second-order moment
 * (i.e. variance) may only involve a scalar or vector.
 *
 * \param[in]  name           name of associated moment
 * \param[in]  f              pointer to associated field
 * \param[in]  type           moment type
 * \param[in]  nt_start       starting time step (or -1 to use t_start)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 * \param[in]  restart_name   if non-null, previous name in case of restart
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_define_by_field(const char                *name,
                               const cs_field_t          *f,
                               cs_time_moment_type_t      type,
                               int                        nt_start,
                               double                     t_start,
                               cs_time_moment_restart_t   restart_mode,
                               const char                *restart_name)
{
  int m_id = -1;
  bool is_intensive = true;
  int component_ids[1] = {-1};
  int field_ids[1] = {f->id};
  int sd_id = _find_or_add_sd(name, 1, field_ids, component_ids,
                              &is_intensive);

  const int *msd = _moment_sd_defs[sd_id];

  m_id = _time_moment_define_by_func(name,
                                     msd[0],
                                     msd[1],
                                     is_intensive,
                                     nullptr,
                                     _sd_moment_data,
                                     msd,
                                     nullptr,
                                     nullptr,
                                     type,
                                     nt_start,
                                     t_start,
                                     restart_mode,
                                     restart_name);

  return m_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time moment based on an evaluation function.
 *
 * Moments will involve the tensor products of their component fields,
 * and only scalar, vector, or rank-2 tensors are handled (for
 * post-processing output reasons), so a 1st-order moment (i.e. mean) may
 * involve a scalar, vector, or tensor, while a second-order moment
 * (i.e. variance) may only involve a scalar or vector.
 *
 * \param[in]  name           name of associated moment
 * \param[in]  f              pointer to function object
 * \param[in]  type           moment type
 * \param[in]  nt_start       starting time step (or -1 to use t_start)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 * \param[in]  restart_name   if non-null, previous name in case of restart
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_define_by_function(const char                *name,
                                  cs_function_t             *f,
                                  cs_time_moment_type_t      type,
                                  int                        nt_start,
                                  double                     t_start,
                                  cs_time_moment_restart_t   restart_mode,
                                  const char                *restart_name)
{
  bool is_intensive = f->type & CS_FUNCTION_INTENSIVE;

  int m_id = _time_moment_define_by_func(name,
                                         f->location_id,
                                         f->dim,
                                         is_intensive,
                                         f,
                                         nullptr,
                                         f->func_input,
                                         nullptr,
                                         nullptr,
                                         type,
                                         nt_start,
                                         t_start,
                                         restart_mode,
                                         restart_name);

  return m_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a moment whose data values will be computed using a
 * specified function.
 *
 * \param[in]  name           name of associated moment
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  dim            dimension associated with element data
 * \param[in]  is_intensive   is the time moment intensive?
 * \param[in]  data_func      function used to define data values
 * \param[in]  data_input     pointer to optional (untyped) value or structure
 *                            to be used by data_func
 * \param[in]  w_data_func    function used to define weight values
 * \param[in]  w_data_input   pointer to optional (untyped) value or structure
 *                            to be used by w_data_func
 * \param[in]  type           moment type
 * \param[in]  nt_start       starting time step (or -1 to use t_start)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 * \param[in]  restart_name   if non-null, previous name in case of restart
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_define_by_func(const char                *name,
                              int                        location_id,
                              int                        dim,
                              bool                       is_intensive,
                              cs_time_moment_data_t     *data_func,
                              const void                *data_input,
                              cs_time_moment_data_t     *w_data_func,
                              void                      *w_data_input,
                              cs_time_moment_type_t      type,
                              int                        nt_start,
                              double                     t_start,
                              cs_time_moment_restart_t   restart_mode,
                              const char                *restart_name)
{
  int moment_id = _time_moment_define_by_func(name,
                                              location_id,
                                              dim,
                                              is_intensive,
                                              nullptr,
                                              data_func,
                                              data_input,
                                              w_data_func,
                                              w_data_input,
                                              type,
                                              nt_start,
                                              t_start,
                                              restart_mode,
                                              restart_name);

  return moment_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined time moments.
 *
 * \return  number of defined time moments
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_n_moments(void)
{
  return _n_moments;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of time moments in the restart file, if any
 *
 * \return  number of defined moments in restart file, or 0
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_n_moments_restart(void)
{
  int n_restart_moments = 0;

  if (_restart_info_checked == false)
    _restart_info_read();

  if (_restart_info != nullptr)
    n_restart_moments = _restart_info->n_moments;

  return n_restart_moments;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a moment restart mode and name by an id.
 *
 * This is a utility function, to allow simplification of automatic setups.
 * It must be called just before defining a moment to work properly if
 * restart_id < -1 (automatic mode).
 *
 * \param[in]   restart_id    -2: automatic, -1: reset, >= 0: id of
 *                            matching moment in restart data
 * \param[out]  restart_mode  matching restart mode
 * \param[out]  restart_name  matching restart name
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_restart_options_by_id(int                         restart_id,
                                     cs_time_moment_restart_t   *restart_mode,
                                     const char                **restart_name)
{
  *restart_name = nullptr;
  if (restart_id < -1) {
    *restart_mode = CS_TIME_MOMENT_RESTART_AUTO;
    if (_restart_info_checked == false)
      _restart_info_read();
  }
  else if (restart_id == -1)
    *restart_mode = CS_TIME_MOMENT_RESTART_RESET;
  else {
    *restart_name = cs_time_moment_restart_name(restart_id);
    *restart_mode = CS_TIME_MOMENT_RESTART_AUTO;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return name of a given time moments in the restart file, if any
 *        (check also \ref cs_time_moment_n_moments_restart).
 *
 * \param[in]  restart_id  id of time moment in restart data
 *
 * \return  name of defined moment in restart file, or null
 */
/*----------------------------------------------------------------------------*/

const char *
cs_time_moment_restart_name(int  restart_id)
{
  const char *retval = nullptr;

  if (_restart_info_checked == false)
    _restart_info_read();

  if (_restart_info != nullptr) {
    if (restart_id >= 0 && restart_id < _restart_info->n_moments)
      retval = _restart_info->name[restart_id];
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to field associated with a given moment.
 *
 * For moments defined automatically to assist computation of higher order
 * moments, which do not have an associated field, a null pointer is returned.
 *
 * \param[in]  moment_id  id of associated moment
 *
 * \return  pointer to field associated with given moment, or null
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_time_moment_get_field(int  moment_id)
{
  assert(moment_id >= 0 && moment_id < _n_moments);

  const cs_time_moment_t *mt = _moment + moment_id;
  if (mt->f_id > -1) {
    cs_field_t *f = cs_field_by_id(mt->f_id);
    return f;
  }

  return nullptr;
}

END_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset selected time step for starting time step of selected moment.
 *
 * All other time moments sharing the same start time will also start
 * at the same time step.
 *
 * \param[in]   moment_id  id of associated moment, or -1 for all
 * \param[in]   nt_start   starting time step
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_set_start_time(int   moment_id,
                              int   nt_start)
{
  _time_moment_set_start_time(moment_id, nt_start, -1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset selected time step for starting time step of selected moment.
 *
 * All other time moments sharing the same start time will also start
 * at the same time step.
 *
 * \param[in]   moment_id  id of associated moment, or -1 for all
 * \param[in]   nt_start   starting time value
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_set_start_time(int     moment_id,
                              double  t_start)
{
  _time_moment_set_start_time(moment_id, -1, t_start);
}

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 1 if moment is active, 0 if it is not active yet.
 *
 * \param[in]   moment_id  id of associated moment
 *
 * \return 1 if moment is active, 0 if it is not active yet
 */
/*----------------------------------------------------------------------------*/

int
cs_time_moment_is_active(int  moment_id)
{
  int retval = 1;

  assert(moment_id >= 0 && moment_id < _n_moments);

  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_time_moment_t *mt = _moment + moment_id;
  const cs_time_moment_wa_t *mwa = _moment_wa + mt->wa_id;

  if (mwa->nt_start < 0 || mwa->nt_start > ts->nt_cur)
    retval = 0;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current iteration as starting time step of selected moment.
 *
 * All other time moments sharing the same start time should also be reset.
 *
 * \param[in]   moment_id  id of associated moment
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_reset(int   moment_id)
{
  const cs_time_step_t  *ts = cs_glob_time_step;

  int s_id = moment_id;
  int e_id = moment_id+1;
  if (moment_id == -1) {
    s_id = 0;
    e_id = _n_moments;
  }

  for (int m_id = s_id; m_id < e_id; m_id++) {

    cs_time_moment_t *mt = _moment + m_id;
    cs_time_moment_wa_t *mwa = _moment_wa + mt->wa_id;

    /* Set current iteration as starting time step for current moment */
    mt->nt_cur = -1;
    mwa->nt_start = ts->nt_cur;
    mwa->t_start = -1.;

    cs_lnum_t n_elts = cs_mesh_location_get_n_elts(mt->location_id)[0];
    cs_lnum_t nd = n_elts * mt->dim;

    cs_real_t *restrict val = mt->val;
    if (mt->f_id > -1) {
      cs_field_t *f = cs_field_by_id(mt->f_id);
      val = f->val;
    }

    for (cs_lnum_t i = 0; i < nd; i++) {
      /* reset moment values */
      val[i] = 0.;
    }

    _reset_weight_accumulator(mwa);

    /* sub-moments (means for variance) */

    if (mt->l_id > -1) {
      int l_id = mt->l_id;

      mt = _moment + l_id;
      mwa = _moment_wa + mt->wa_id;

      mt->nt_cur = -1;
      mwa->nt_start = ts->nt_cur;
      mwa->t_start = -1.;

      n_elts = cs_mesh_location_get_n_elts(mt->location_id)[0];
      nd = n_elts * mt->dim;

      val = mt->val;
      if (mt->f_id > -1) {
        cs_field_t *f = cs_field_by_id(mt->f_id);
        val = f->val;
      }

      for (cs_lnum_t i = 0; i < nd; i++) {
        /* reset moment values */
        val[i] = 0.;
      }

      _reset_weight_accumulator(mwa);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map time step values for temporal moments.
 *
 * If this function is not called, the value from cs_glob_time_step->dt_ref
 * will be used instead.
 *
 * \param[in]   dt   pointer to time step values array
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_map_cell_dt(const cs_real_t  *dt)
{
  _p_dt = dt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update all moment accumulators.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_update_all(void)
{
  bool active_moments = false;

  const cs_time_step_t  *ts = cs_glob_time_step;

  const cs_real_t dt_val = (_p_dt != nullptr) ? *_p_dt : ts->dt_ref;

  /* Prepare accumulators */

  double **wa_cur_data;
  double *wa_cur_data0;

  for (int i = 0; i < _n_moment_wa; i++) {

    cs_time_moment_wa_t *mwa = _moment_wa + i;

    /* start time step equal to current iteration */
    if (mwa->t_start < 0. && mwa->nt_start <= ts->nt_cur)
      mwa->t_start = _t_prev_iter;
    /* start time value in interval [t_prev-0.01*dt^(n-1), t_cur-0.01*dt^n[ */
    else if (mwa->nt_start < 0 && mwa->t_start < ts->t_cur - 0.01*dt_val)
      mwa->nt_start = ts->nt_cur;

    if (mwa->nt_start > -1 && mwa->nt_start <= ts->nt_cur)
      active_moments = true;

  }

  _t_prev_iter = ts->t_cur;

  if (!active_moments)
    return;

  cs_dispatch_context ctx;
  bool on_device = ctx.use_gpu();

  CS_MALLOC(wa_cur_data, _n_moment_wa, cs_real_t *);
  CS_MALLOC_HD(wa_cur_data0, _n_moment_wa, cs_real_t, cs_alloc_mode);

  /* Compute current weight data */

  for (int i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    if (mwa->nt_start > -1 && mwa->nt_start <= ts->nt_cur) {
      _ensure_init_weight_accumulator(mwa);
      wa_cur_data[i] = _compute_current_weight(mwa,
                                               dt_val,
                                               wa_cur_data0 + i);
    }
    else
      wa_cur_data[i] = nullptr;
  }

  /* Loop on variances first */

  for (int m_type = CS_TIME_MOMENT_VARIANCE;
       m_type >= (int)CS_TIME_MOMENT_MEAN;
       m_type --) {

    for (int i = 0; i < _n_moments; i++) {

      cs_time_moment_t *mt = _moment + i;
      cs_time_moment_wa_t *mwa = _moment_wa + mt->wa_id;

      if (   mt->nt_cur < ts->nt_cur
          && (int)(mt->type) == m_type
          && (mwa->nt_start > -1 && mwa->nt_start <= ts->nt_cur)) {

        /* Current and accumulated weight */

        cs_lnum_t  wa_stride;
        cs_real_t *restrict wa_sum, *restrict x;

        cs_real_t *const restrict w = wa_cur_data[mt->wa_id];

        if (mwa->location_id == CS_MESH_LOCATION_NONE) {
          wa_sum = &(mwa->val0);
          wa_stride = 0;
        }
        else {
          wa_sum = mwa->val;
          wa_stride = 1;
        }

        /* Current value */

        const cs_lnum_t n_elts
          = cs_mesh_location_get_n_elts(mt->location_id)[0];
        const cs_lnum_t nd = n_elts * mt->dim;

        CS_MALLOC_HD(x, nd, cs_real_t, cs_alloc_mode);

        assert(mt->data_func == nullptr || mt->eval_func == nullptr);

        if (mt->data_func != nullptr) {
          mt->data_func(mt->data_input, x);
        }
        else {
          const int location_id = mt->location_id;
          const cs_lnum_t* elt_ids = cs_mesh_location_get_elt_ids(mt->location_id);
          cs_function_evaluate(mt->eval_func,
                               ts,
                               location_id,
                               n_elts,
                               elt_ids,
                               x);
        }

        _ensure_init_moment(mt);

        cs_real_t *restrict val = mt->val;
        if (mt->f_id > -1) {
          cs_field_t *f = cs_field_by_id(mt->f_id);
          val = f->val;
        }

        if (mt->type == CS_TIME_MOMENT_VARIANCE) {

          assert(mt->l_id > -1);

          cs_time_moment_t *mt_mean = _moment + mt->l_id;

          _ensure_init_moment(mt_mean);
          cs_real_t *restrict m = mt_mean->val;
          if (mt_mean->f_id > -1) {
            cs_field_t *f_mean = cs_field_by_id(mt_mean->f_id);
            m = f_mean->val;
          }

          if (mt->dim == 6) { /* variance-covariance matrix */
            assert(mt->data_dim == 3);
            ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t je) {
              double delta[3], delta_n[3], r[3], m_n[3];
              const cs_lnum_t k = je*wa_stride;
              const double wa_sum_n = w[k] + wa_sum[k];
              for (cs_lnum_t l = 0; l < 3; l++) {
                cs_lnum_t jl = je*6 + l, jml = je*3 + l;
                delta[l]   = x[jml] - m[jml];
                r[l] = delta[l] * (w[k] / wa_sum_n);
                m_n[l] = m[jml] + r[l];
                delta_n[l] = x[jml] - m_n[l];
                val[jl] =   (val[jl]*wa_sum[k] + (w[k]*delta[l]*delta_n[l]))
                          / wa_sum_n;
              }
              /* Covariance terms.
                 Note we could have a symmetric formula using
                   0.5*(delta[i]*delta_n[j] + delta[j]*delta_n[i])
                 instead of
                   delta[i]*delta_n[j]
                 but unit tests in cs_moment_test.c do not seem to favor
                 one variant over the other; we use the simplest one.
              */
              cs_lnum_t j3 = je*6 + 3, j4 = je*6 + 4, j5 = je*6 + 5;
              val[j3] =   (val[j3]*wa_sum[k] + (w[k]*delta[0]*delta_n[1]))
                        / wa_sum_n;
              val[j4] =   (val[j4]*wa_sum[k] + (w[k]*delta[1]*delta_n[2]))
                        / wa_sum_n;
              val[j5] =   (val[j5]*wa_sum[k] + (w[k]*delta[0]*delta_n[2]))
                        / wa_sum_n;
              for (cs_lnum_t l = 0; l < 3; l++)
                m[je*3 + l] += r[l];
            });
          }

          else { /* simple variance */
            cs_lnum_t k_stride = wa_stride / mt->dim;
            ctx.parallel_for(nd, [=] CS_F_HOST_DEVICE (cs_lnum_t j) {
              const cs_lnum_t k = j*k_stride;
              double wa_sum_n = w[k] + wa_sum[k];
              double delta = x[j] - m[j];
              double r = delta * (w[k] / wa_sum_n);
              double m_n = m[j] + r;
              val[j] = (val[j]*wa_sum[k] + (w[k]*delta*(x[j]-m_n))) / wa_sum_n;
              m[j] += r;
            });
          }

          mt_mean->nt_cur = ts->nt_cur;
        }

        else if (mt->type == CS_TIME_MOMENT_MEAN) {

          cs_lnum_t k_stride = wa_stride / mt->dim;
          ctx.parallel_for(nd, [=] CS_F_HOST_DEVICE (cs_lnum_t j) {
            const cs_lnum_t k = j*k_stride;
            val[j] += (x[j] - val[j]) * (w[k] / (w[k] + wa_sum[k]));
          });

        }

        mt->nt_cur = ts->nt_cur;

        ctx.wait();

        CS_FREE(x);

        /* Sync ghost cells so downstream use is safe */

        if (mt->location_id == CS_MESH_LOCATION_CELLS) {
          const cs_halo_t *halo = cs_glob_mesh->halo;
          if (halo != nullptr) {
            switch (mt->dim) {
            case 1:
              cs_halo_sync(halo, CS_HALO_EXTENDED, on_device, val);
              break;
            case 3:
              cs_halo_sync_r(halo, CS_HALO_EXTENDED, on_device,
                             reinterpret_cast<cs_real_3_t *>(val));
              break;
            case 6:
              cs_halo_sync_r(halo, CS_HALO_EXTENDED, on_device,
                             reinterpret_cast<cs_real_6_t *>(val));
              break;
            default:
              cs_halo_sync_var_strided(halo, CS_HALO_EXTENDED, val, mt->dim);
            }
          }
        }

      } /* End of test if moment is active */

    } /* End of loop on moments */

  } /* End of loop on moment types */

  /* Update and free weight data */

  for (int i = 0; i < _n_moment_wa; i++) {
    if (wa_cur_data[i] != nullptr) {
      _update_weight_accumulator(_moment_wa + i, wa_cur_data[i]);
      if (wa_cur_data[i] != wa_cur_data0 + i)
        CS_FREE(wa_cur_data[i]);
    }
  }

  CS_FREE(wa_cur_data0);
  CS_FREE(wa_cur_data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log moment definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_log_setup(void)
{
  if (_n_moments < 1)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Temporal moments\n"
                  "----------------\n"));

  /* Accumulator info */

  {
    char s[64];
    char tmp_s[4][64] =  {"", "", "", ""};

    /* Header */

    cs_log_strpad(tmp_s[0], _("Accumulator"), 16, 64);
    cs_log_strpad(tmp_s[1], _("Location"), 20, 64);
    cs_log_strpad(tmp_s[2], _("Start"), 16, 64);
    cs_log_strpad(tmp_s[3], _("Weight"), 16, 64);

    cs_log_printf(CS_LOG_SETUP, "\n");

    cs_log_printf(CS_LOG_SETUP, "  %s %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

    for (int j = 0; j < 4; j++)
      memset(tmp_s[j], '-', 64);

    tmp_s[0][16] = '\0';
    tmp_s[1][20] = '\0';
    tmp_s[2][16] = '\0';
    tmp_s[3][16] = '\0';

    cs_log_printf(CS_LOG_SETUP, "  %s %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

    /* Accumulators */

    for (int i = 0; i < _n_moment_wa; i++) {

      cs_time_moment_wa_t *mwa = _moment_wa + i;

      cs_log_strpad(tmp_s[1],
                    _(cs_mesh_location_get_name(mwa->location_id)),
                    20,
                    64);

      if (mwa->t_start >= 0)
        snprintf(s, 17, "%10.5g s", mwa->t_start);
      else
        snprintf(s, 17, "nt %d", mwa->nt_start);
      cs_log_strpad(tmp_s[2], s, 16, 64);

      if (mwa->data_func != nullptr)
        cs_log_strpad(tmp_s[3], _("user"), 16, 64);
      else
        cs_log_strpad(tmp_s[3], "-", 16, 64);

      cs_log_printf(CS_LOG_SETUP,
                    "  %-16d %s %s %s\n",
                    i, tmp_s[1], tmp_s[2], tmp_s[3]);

    }
  }

  /* Base moments info */

  {
    char s[64];
    char tmp_s[8][64] =  {"", "", "", "", "", "", "", ""};

    size_t name_width = 16;

    /* First loop to determine name width */

    for (int i = 0; i < _n_moments; i++) {
      cs_time_moment_t *mt = _moment + i;
      if (mt->f_id > -1) {
        const cs_field_t *f = cs_field_by_id(mt->f_id);
        size_t l = strlen(f->name);
        if (l > name_width)
          name_width = l;
      }
    }
    if (name_width > 63)
      name_width = 63;

    /* Header */

    cs_log_strpad(tmp_s[0], _("Moment"), name_width, 64);
    cs_log_strpad(tmp_s[1], _("Dim."), 4, 64);
    cs_log_strpad(tmp_s[2], _("Location"), 20, 64);
    cs_log_strpad(tmp_s[3], _("Type"), 8, 64);
    cs_log_strpad(tmp_s[4], _("Id"), 4, 64);
    cs_log_strpad(tmp_s[5], _("Acc."), 4, 64);
    cs_log_strpad(tmp_s[6], _("Lower"), 6, 64);
    cs_log_strpad(tmp_s[7], _("Field"), 6, 64);

    cs_log_printf(CS_LOG_SETUP, "\n");

    cs_log_printf(CS_LOG_SETUP, "  %s %s %s %s %s %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3],
                  tmp_s[4], tmp_s[5], tmp_s[6], tmp_s[7]);

    for (int j = 0; j < 8; j++)
      memset(tmp_s[j], '-', 64);

    tmp_s[0][name_width] = '\0';
    tmp_s[1][4] = '\0';
    tmp_s[2][20] = '\0';
    tmp_s[3][8] = '\0';
    tmp_s[4][4] = '\0';
    tmp_s[5][4] = '\0';
    tmp_s[6][6] = '\0';
    tmp_s[7][6] = '\0';

    cs_log_printf(CS_LOG_SETUP, "  %s %s %s %s %s %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3],
                  tmp_s[4], tmp_s[5], tmp_s[6], tmp_s[7]);

    /* Main loop on moments */

    for (int i = 0; i < _n_moments; i++) {

      cs_time_moment_t *mt = _moment + i;

      if (mt->f_id > -1) {
        const cs_field_t *f = cs_field_by_id(mt->f_id);
        cs_log_strpad(tmp_s[0], f->name, name_width, 64);
      }
      else
        cs_log_strpad(tmp_s[0], mt->name, name_width, 64);

      cs_log_strpad(tmp_s[2],
                    _(cs_mesh_location_get_name(mt->location_id)),
                    20,
                    64);

      cs_log_strpad(tmp_s[3],
                    _(cs_time_moment_type_name[mt->type]),
                    8,
                    64);

      if (mt->l_id > -1)
        snprintf(s, 64, "%d", mt->l_id);
      else
        snprintf(s, 64, "-");
      cs_log_strpad(tmp_s[6], s, 6, 64);

      if (mt->f_id > -1)
        snprintf(tmp_s[7], 64, "%d", mt->f_id);
      else
        snprintf(tmp_s[7], 64, "-");

      cs_log_printf(CS_LOG_SETUP, "  %s %-4d %s %s %-4d %-4d %s %s\n",
                    tmp_s[0], mt->dim, tmp_s[2], tmp_s[3],
                    i, mt->wa_id, tmp_s[6], tmp_s[7]);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log moment definition information for a given iteration.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_log_iteration(void)
{
  if (_n_moment_wa < 1)
    return;

  int n_active_wa[] = {0, 0};

  const cs_time_step_t  *ts = cs_glob_time_step;

  for (int i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    if (mwa->nt_start <= ts->nt_cur) {
      if (mwa->location_id == 0)
        n_active_wa[0] += 1;
      else
        n_active_wa[1] += 1;
    }
  }

  if (n_active_wa[0] + n_active_wa[1] < 1)
    return;

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  "  ** Temporal moment accumulated weights\n"
                  "     -----------------------------------\n"));

  cs_dispatch_context ctx;
  /* Info for accumulators on global locations */

  if (n_active_wa[0] > 0) {

    char tmp_s[3][64] =  {"", "", ""};

    /* Header */

    cs_log_strpad(tmp_s[0], _("id"), 4, 64);
    cs_log_strpad(tmp_s[1], _("n it."), 8, 64);
    cs_log_strpadl(tmp_s[2], _("value"), 14, 64);

    cs_log_printf(CS_LOG_DEFAULT, "\n");

    cs_log_printf(CS_LOG_DEFAULT, "   %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2]);

    for (int j = 0; j < 3; j++)
      memset(tmp_s[j], '-', 64);

    tmp_s[0][4] = '\0';
    tmp_s[1][8] = '\0';
    tmp_s[2][14] = '\0';

    cs_log_printf(CS_LOG_DEFAULT, "   %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2]);

    /* Now log values */

    for (int i = 0; i < _n_moment_wa; i++) {
      cs_time_moment_wa_t *mwa = _moment_wa + i;
      if (mwa->nt_start <= ts->nt_cur && mwa->location_id == 0) {
        int nt_acc = ts->nt_cur - mwa->nt_start + 1;
        cs_log_printf(CS_LOG_DEFAULT, "   %-4d %-8d %14.5g\n",
                      i, nt_acc, mwa->val0);
      }
    }

  }

  /* Info for accumulators on non-global locations */

  if (n_active_wa[1] > 0) {

    char tmp_s[6][64] =  {"", "", "", "", "", ""};

    /* Header */

    cs_log_strpad(tmp_s[0], _("id"), 4, 64);
    cs_log_strpad(tmp_s[1], _("location"), 20, 64);
    cs_log_strpad(tmp_s[2], _("n it."), 8, 64);
    cs_log_strpadl(tmp_s[3], _("minimum"), 14, 64);
    cs_log_strpadl(tmp_s[4], _("maximum"), 14, 64);
    cs_log_strpadl(tmp_s[5], _("set mean"), 14, 64);

    cs_log_printf(CS_LOG_DEFAULT, "\n");

    cs_log_printf(CS_LOG_DEFAULT, "   %s %s %s %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4], tmp_s[5]);

    for (int j = 0; j < 6; j++)
      memset(tmp_s[j], '-', 64);

    tmp_s[0][4] = '\0';
    tmp_s[1][20] = '\0';
    tmp_s[2][8] = '\0';
    tmp_s[3][14] = '\0';
    tmp_s[4][14] = '\0';
    tmp_s[5][14] = '\0';

    cs_log_printf(CS_LOG_DEFAULT, "   %s %s %s %s %s %s\n",
                  tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4], tmp_s[5]);

    cs_gnum_t *n_g_elts;
    double *vmin, *vmax, *vsum;

    CS_MALLOC(n_g_elts, n_active_wa[1], cs_gnum_t);
    CS_MALLOC(vmin, n_active_wa[1], double);
    CS_MALLOC(vmax, n_active_wa[1], double);
    CS_MALLOC(vsum, n_active_wa[1], double);

    n_active_wa[1] = 0;

    /* Determine min, max, sum */

    for (int i = 0; i < _n_moment_wa; i++) {
      cs_time_moment_wa_t *mwa = _moment_wa + i;
      if (mwa->nt_start <= ts->nt_cur && mwa->location_id > 0) {
        const cs_lnum_t n_elts
          = cs_mesh_location_get_n_elts(mwa->location_id)[0];
        const cs_mesh_location_type_t loc_type
          = cs_mesh_location_get_type(mwa->location_id);
        if (   loc_type == CS_MESH_LOCATION_CELLS
            || loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
          n_g_elts[n_active_wa[1]] = n_elts;
        else
          n_g_elts[n_active_wa[1]] = 0;
        _ensure_init_weight_accumulator(mwa);
        cs_array_reduce_simple_stats_l(ctx,
                                       n_elts,
                                       1,
                                       nullptr,
                                       mwa->val,
                                       vmin + n_active_wa[1],
                                       vmax + n_active_wa[1],
                                       vsum + n_active_wa[1]);
        n_active_wa[1] += 1;
      }
    }

    /* Group MPI operations if required */

    cs_parall_counter(n_g_elts, n_active_wa[1]);
    cs_parall_min(n_active_wa[1], CS_DOUBLE, vmin);
    cs_parall_max(n_active_wa[1], CS_DOUBLE, vmax);
    cs_parall_sum(n_active_wa[1], CS_DOUBLE, vsum);

    /* Now log values */

    n_active_wa[1] = 0;

    for (int i = 0; i < _n_moment_wa; i++) {
      cs_time_moment_wa_t *mwa = _moment_wa + i;
      if (mwa->nt_start <= ts->nt_cur && mwa->location_id > 0) {

        cs_log_strpad(tmp_s[1],
                      _(cs_mesh_location_get_name(mwa->location_id)),
                      20,
                      64);

        int nt_acc = ts->nt_cur - mwa->nt_start + 1;

        if (n_g_elts[n_active_wa[1]] > 0) {
          double v_mean = vsum[n_active_wa[1]] / n_g_elts[n_active_wa[1]];
          snprintf(tmp_s[5], 63, " %14.5g", v_mean);
          tmp_s[5][63] = '\0';
        }
        else
          tmp_s[5][0] = '\0';

        cs_log_printf(CS_LOG_DEFAULT, "   %-4d %s %-8d %14.5g %14.5g%s\n",
                      i, tmp_s[1], nt_acc,
                      vmin[n_active_wa[1]], vmax[n_active_wa[1]], tmp_s[5]);

        n_active_wa[1] += 1;

      }
    }

    CS_FREE(vsum);
    CS_FREE(vmax);
    CS_FREE(vmin);
    CS_FREE(n_g_elts);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if restart API should use "main" instead of "auxiliary"
 *        file.
 *
 * \param[in]  use_main  use "main" restart if nonzero, "auxiliary" otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_restart_use_main(int  use_main)
{
  if (use_main)
    _restart_uses_main = true;
  else
    _restart_uses_main = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read restart moment data
 *
 * \param[in]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_restart_read(cs_restart_t  *restart)
{
  if (_n_moment_wa < 1)
    return;

  int retcode;

  /* Initialize */

  const cs_time_step_t  *ts = cs_glob_time_step;
  _t_prev_iter = ts->t_prev;

  if (_restart_info == nullptr)
    _restart_info_read_auxiliary(restart);

  if (_restart_info == nullptr)
    return;

  cs_time_moment_restart_info_t  *ri = _restart_info;

  /* Read information proper */

  for (int i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    if (mwa->restart_id > -1 && mwa->location_id > CS_MESH_LOCATION_NONE) {
      char s[64];
      snprintf(s, 64, "time_moments:wa:%02d:val", mwa->restart_id);
      _ensure_init_weight_accumulator(mwa);
      retcode = cs_restart_read_section(restart,
                                        s,
                                        mwa->location_id,
                                        1,
                                        CS_TYPE_cs_real_t,
                                        mwa->val);
      _assert_restart_success(retcode);
    }
  }

  for (int i = 0; i < _n_moments; i++) {
    cs_time_moment_t *mt = _moment + i;
    if (mt->restart_id > -1) {
      _ensure_init_moment(mt);
      cs_real_t *val = mt->val;
      if (mt->f_id > -1) {
        cs_field_t *f = cs_field_by_id(mt->f_id);
        val = f->val;
      }
      retcode = cs_restart_read_section(restart,
                                        ri->name[mt->restart_id],
                                        mt->location_id,
                                        mt->dim,
                                        CS_TYPE_cs_real_t,
                                        val);
      _assert_restart_success(retcode);
    }
  }

  /* Free info */

  _restart_info_free();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Checkpoint moment data
 *
 * \param[in]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_time_moment_restart_write(cs_restart_t  *restart)
{
  int *nt_start, *location_id, *dimension, *m_type, *wa_id, *l_id;
  cs_real_t *t_start, *val0;

  int n_active_wa = 0, n_active_moments = 0;
  int *active_wa_id = nullptr, *active_moment_id = nullptr;

  if (_n_moments < 1)
    return;

  const cs_time_step_t  *ts = cs_glob_time_step;

  /* General information */
  /* ------------------- */

  CS_MALLOC(active_wa_id, _n_moment_wa, int);
  CS_MALLOC(active_moment_id, _n_moments, int);

  /* Check for active moments */

  for (int i = 0; i < _n_moment_wa; i++) {
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    if (mwa->nt_start > -1 && mwa->nt_start <= ts->nt_cur) {
      active_wa_id[i] = n_active_wa;
      n_active_wa += 1;
    }
    else
      active_wa_id[i] = -1;
  }

  for (int i = 0; i < _n_moments; i++) {
    cs_time_moment_t *mt = _moment + i;
    if (active_wa_id[mt->wa_id] > -1) {
      active_moment_id[i] = n_active_moments;
      n_active_moments += 1;
    }
    else
      active_moment_id[i] = -1;
  }

  if (n_active_moments < 1) {
    CS_FREE(active_wa_id);
    CS_FREE(active_moment_id);
    return;
  }

  /* Build global names array */

  size_t names_max_size = 32;
  int *names_idx;
  char *names;

  CS_MALLOC(names_idx, n_active_moments + 1, int);
  CS_MALLOC(names, names_max_size, char);

  names_idx[0] = 0;

  for (int i = 0; i < _n_moments; i++) {

    const int j = active_moment_id[i];
    if (j > -1) {

      cs_time_moment_t *mt = _moment + i;
      const char *name = nullptr;
      if (mt->f_id > -1) {
        const cs_field_t *f = cs_field_by_id(mt->f_id);
        name = f->name;
      }
      else
        name = mt->name;
      const size_t l = strlen(name) + 1;
      if (names_idx[j] + l > names_max_size) {
        while (names_idx[j] + l > names_max_size)
          names_max_size *= 2;
        CS_REALLOC(names, names_max_size, char);
      }
      strcpy(names + names_idx[j], name);
      names[names_idx[j] + l - 1] = '\0';
      names_idx[j+1] = names_idx[j] + l;

    }

  }

  cs_lnum_t sizes[3] = {n_active_wa,
                        n_active_moments,
                        names_idx[n_active_moments]};

  cs_restart_write_section(restart,
                           "time_moments:sizes",
                           CS_MESH_LOCATION_NONE,
                           3,
                           CS_TYPE_int,
                           sizes);

  cs_restart_write_section(restart,
                           "time_moments:names",
                           CS_MESH_LOCATION_NONE,
                           names_idx[n_active_moments],
                           CS_TYPE_char,
                           names);

  CS_FREE(names_idx);
  CS_FREE(names);

  /* Information on moment weight accumulators */

  CS_MALLOC(location_id, n_active_wa, int);
  CS_MALLOC(nt_start, n_active_wa, int);
  CS_MALLOC(t_start, n_active_wa, cs_real_t);
  CS_MALLOC_HD(val0, n_active_wa, cs_real_t, cs_alloc_mode);

  int n_val0 = 0;
  for (int i = 0; i < _n_moment_wa; i++) {
    int j = active_wa_id[i];
    if (j > -1) {
      cs_time_moment_wa_t *mwa = _moment_wa + i;
      location_id[j] = mwa->location_id;
      nt_start[j] = mwa->nt_start;
      t_start[j] = mwa->t_start;
      val0[j] = mwa->val0;
      if (mwa->location_id == CS_MESH_LOCATION_NONE)
        n_val0 += 1;
    }
  }

  cs_restart_write_section(restart,
                           "time_moments:wa:location_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_wa,
                           CS_TYPE_int,
                           location_id);

  cs_restart_write_section(restart,
                           "time_moments:wa:nt_start",
                           CS_MESH_LOCATION_NONE,
                           n_active_wa,
                           CS_TYPE_int,
                           nt_start);

  cs_restart_write_section(restart,
                           "time_moments:wa:t_start",
                           CS_MESH_LOCATION_NONE,
                           n_active_wa,
                           CS_TYPE_cs_real_t,
                           t_start);

  if (n_val0 > 0)
    cs_restart_write_section(restart,
                             "time_moments:wa:val_g",
                             CS_MESH_LOCATION_NONE,
                             n_active_wa,
                             CS_TYPE_cs_real_t,
                             val0);

  CS_FREE(val0);
  CS_FREE(t_start);
  CS_FREE(nt_start);
  CS_FREE(location_id);

  for (int i = 0; i < _n_moment_wa; i++) {
    int j = active_wa_id[i];
    cs_time_moment_wa_t *mwa = _moment_wa + i;
    if (j > -1 && mwa->location_id > CS_MESH_LOCATION_NONE) {
      char s[64];
      snprintf(s, 64, "time_moments:wa:%02d:val", i);
      cs_restart_write_section(restart,
                               s,
                               mwa->location_id,
                               1,
                               CS_TYPE_cs_real_t,
                               mwa->val);
    }
  }

  /* Information on moments proper */

  CS_MALLOC(m_type, n_active_moments, int);
  CS_MALLOC(location_id, n_active_moments, int);
  CS_MALLOC(dimension, n_active_moments, int);
  CS_MALLOC(wa_id, n_active_moments, int);
  CS_MALLOC(l_id, n_active_moments, int);

  for (int i = 0; i < _n_moments; i++) {
    int j = active_moment_id[i];
    if (j > -1) {
      cs_time_moment_t *mt = _moment + i;
      m_type[j] = mt->type;
      location_id[j] = mt->location_id;
      dimension[j] = mt->dim;
      wa_id[j] = active_wa_id[mt->wa_id];
      if (mt->l_id > -1)
        l_id[j] = active_moment_id[mt->l_id];
      else
        l_id[j] = -1;
    }
  }

  cs_restart_write_section(restart,
                           "time_moments:type",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_int,
                           m_type);

  cs_restart_write_section(restart,
                           "time_moments:location_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_int,
                           location_id);

  cs_restart_write_section(restart,
                           "time_moments:dimension",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_int,
                           dimension);

  cs_restart_write_section(restart,
                           "time_moments:wa_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_int,
                           wa_id);

  cs_restart_write_section(restart,
                           "time_moments:lower_order_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_int,
                           l_id);

  CS_FREE(l_id);
  CS_FREE(wa_id);
  CS_FREE(dimension);
  CS_FREE(location_id);
  CS_FREE(m_type);

  for (int i = 0; i < _n_moments; i++) {
    int j = active_moment_id[i];
    if (j > -1) {
      cs_time_moment_t *mt = _moment + i;
      if (mt->f_id > -1) {
        const cs_field_t *f = cs_field_by_id(mt->f_id);
        cs_restart_write_section(restart,
                                 f->name,
                                 f->location_id,
                                 f->dim,
                                 CS_TYPE_cs_real_t,
                                 f->val);
      }
      else
        cs_restart_write_section(restart,
                                 mt->name,
                                 mt->location_id,
                                 mt->dim,
                                 CS_TYPE_cs_real_t,
                                 mt->val);
    }
  }

  CS_FREE(active_moment_id);
  CS_FREE(active_wa_id);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
