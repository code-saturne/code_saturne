/*============================================================================
 * Methods for particle statistics
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

/*============================================================================
 * Functions dealing with particle statistics
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"
#include "cs_file.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_restart_default.h"
#include "cs_timer_stats.h"

#include "cs_log.h"
#include "cs_array_reduce.h"
#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr_tracking.h"
#include "cs_lagr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_stat.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Moment time accumulator definition */
/*------------------------------------*/

typedef struct {

  int                       class;        /* Matching statistical class number */

  int                       restart_id;   /* Matching id in restart info */

  int                       f_id;         /* Associated field id, or -1 */

  int                       nt_start;     /* Associated starting time step;
                                             if < 0 (and f_id < 0), t_start is
                                             used directly, but nt_start is
                                             always set to a non-negative value
                                             once accumulation starts) */

  double                    t_start;      /* Associated starting time value
                                             (may be initialized to -1 if
                                             accumulation starts at nt_start,
                                             but is always set to a non-negative
                                             value once accumulation starts) */

  int                       location_id;  /* Associated mesh location id */

  cs_lagr_moment_p_data_t  *p_data_func;  /* Associated particle data value
                                             computation function (statistical
                                             weight assumed if NULL) */
  cs_lagr_moment_m_data_t  *m_data_func;  /* Associated mesh data value
                                             computation function, or NULL */
  const void               *data_input;   /* pointer to optional (untyped)
                                             value or structure */
  cs_real_t                 val0;         /* Associated value if location_id
                                             is CS_MESH_LOCATION_NONE */
  cs_real_t                *val;          /* Pointer to associated values,
                                             if f_id < 0 */

} cs_lagr_moment_wa_t;

/* Moment definitions */
/*--------------------*/

typedef struct {

  cs_lagr_stat_moment_t     m_type;       /* Moment type */

  int                       restart_id;   /* Matching id in restart info */

  int                       wa_id;        /* Associated weight accumulator id */

  int                       f_id;         /* Associated field id, or -1 */

  int                       dim;          /* Associated field dimensions */
  int                       data_dim;     /* Associated data field dimensions */
  int                       location_id;  /* Associated mesh location id */


  cs_lagr_moment_p_data_t  *p_data_func;  /* Associated particle data elements
                                             computation function, or NULL */
  cs_lagr_moment_m_data_t  *m_data_func;  /* Associated mesh data elements
                                             computation function, or NULL */
  const void               *data_input;   /* pointer to optional (untyped)
                                             value or structure */

  int                       l_id;         /* Associated id of lower order moment
                                             (mean for variance), or -1 */

  int                       stat_type;    /* Associated type id, or -1 */
  int                       component_id; /* Associated component id, or -1 */

  int                       class;        /* Class number */

  char                     *name;         /* Associated name, if f_id < 0 */
  double                   *val;          /* Associated value, if f_id < 0 */

  int                       nt_cur;      /* Time step number of last update */

} cs_lagr_moment_t;

/* Moment restart metadata */
/*-------------------------*/

typedef struct {

  int                     nt_prev;        /* Restart time step */
  cs_real_t               t_prev;         /* Restart time */

  int                     n_wa;           /* Number of weight accumulators */
  int                     n_moments;      /* Number of moments */

  const char            **name;           /* Moment name */
  char                   *name_buf;       /* Buffer for names */

  int                    *wa_location_id; /* Weight accumulator location ids */
  int                    *wa_nt_start;    /* Weight accumulator start iters. */
  cs_real_t              *wa_t_start;     /* Weight accumulator start times */
  int                    *wa_class;       /* Weight accumulator class */

  int                    *m_type;         /* Moment types */
  int                    *class;          /* Moment class */
  int                    *location_id;    /* Moment location */
  int                    *dimension;      /* Moment dimension */
  int                    *stat_type;      /* Moment pointer number */
  int                    *wa_id;          /* Associated accumulator ids */
  int                    *l_id;           /* Associated lower order ids */

} cs_lagr_moment_restart_info_t;

/* Input helper structure */

typedef struct {

  int  class;        /* Matching statistical class number */
  int  location_id;  /* Matching location id */

} cs_lagr_moment_input_t;

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static  bool *_vol_stat_activate = NULL;

static  bool _restart_info_checked = false;
static  cs_lagr_moment_restart_info_t *_restart_info = NULL;

/* Global Lagragian statistics parameters */
static cs_lagr_moment_wa_t  *_lagr_stats_wa = NULL;
static cs_lagr_moment_t     *_lagr_stats = NULL;

static int  _n_lagr_stats_wa = 0;
static int  _n_lagr_stats_wa_max = 0;

static int  _n_lagr_stats = 0;
static int  _n_lagr_stats_max = 0;

static double _t_prev_iter = 0.;

static const cs_real_t *_p_dt = NULL; /* Mapped cell time step */

/* Names associated with moment types */

const char  *cs_lagr_moment_type_name[] = {N_("MEAN"),
                                           N_("VARIANCE")};

/* lagr statistics structure and associated pointer */

static cs_lagr_stat_options_t _lagr_stat_options
  = {.isuist = 0,
     .idstnt = 0,
     .nstist = 0,
     .threshold = 1e-12};

cs_lagr_stat_options_t *cs_glob_lagr_stat_options = &_lagr_stat_options;

/*============================================================================
 * Private functions definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Particle data function computing unit value for mesh-based data
 *
 * parameters:
 *   input       <-- pointer to optional (untyped) value or structure.
 *   location_id <-- associated mesh location id
 *   class_id    <-- associated particle class id (0 for all)
 *   vals        --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_unit_value_m_cells(const void  *input,
                    int          location_id,
                    int          class_id,
                    cs_real_t    vals[])
{
  CS_UNUSED(input);
  CS_UNUSED(class_id);

  const cs_lnum_t n_elts
    = (location_id == CS_MESH_LOCATION_NONE) ?
      1 : cs_mesh_location_get_n_elts(location_id)[0];

  for (cs_lnum_t i = 0; i < n_elts; i++)
    vals[i] = 1.;
}

/*----------------------------------------------------------------------------
 * Particle data function returning volume fraction
 *
 * parameters:
 *   input       <-- pointer to optional (untyped) value or structure.
 *   location_id <-- associated mesh location id
 *   class_id    <-- associated particle class id (0 for all)
 *   vals        --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_vol_fraction(const void  *input,
              int          location_id,
              int          class_id,
              cs_real_t    vals[])
{
  CS_UNUSED(input);

  cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];
  cs_lagr_particle_set_t *p_set = cs_lagr_get_particle_set();

  for (cs_lnum_t i = 0; i < n_elts; i++)
    vals[i] = 0.;

  if (class_id == 0) {

    for (cs_lnum_t part = 0; part < p_set->n_particles; part++) {

      unsigned char *particle
        = p_set->p_buffer + p_set->p_am->extents * part;

      cs_real_t diam = cs_lagr_particle_get_real(particle, p_set->p_am,
                                                 CS_LAGR_DIAMETER);

      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_set->p_am);

      cs_real_t p_weight = cs_lagr_particle_get_real(particle, p_set->p_am,
                                                     CS_LAGR_STAT_WEIGHT);

      cs_real_t vol = cs_glob_mesh_quantities->cell_vol[cell_id];

      vals[cell_id] += diam*diam*diam * cs_math_pi / 6.0 * p_weight / vol;
    }

  }
  else {

    assert(p_set->p_am->displ[0][CS_LAGR_STAT_CLASS] > 0);

    for (cs_lnum_t part = 0; part < p_set->n_particles; part++) {

      unsigned char *particle
        = p_set->p_buffer + p_set->p_am->extents * part;

      int p_class = cs_lagr_particle_get_lnum(particle,
                                              p_set->p_am,
                                              CS_LAGR_STAT_CLASS);

      if (p_class == class_id) {
        cs_real_t diam = cs_lagr_particle_get_real(particle, p_set->p_am,
                                                   CS_LAGR_DIAMETER);

        cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_set->p_am);

        cs_real_t p_weight = cs_lagr_particle_get_real(particle, p_set->p_am,
                                                       CS_LAGR_STAT_WEIGHT);

        cs_real_t vol = cs_glob_mesh_quantities->cell_vol[cell_id];

        vals[cell_id] += diam*diam*diam * cs_math_pi / 6.0 * p_weight / vol;
      }
    }

  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Check statistics type is in possible range
 *
 * \param[in]   type          moment type
 */
/*---------------------------------------------------------------------------*/

static void
_check_moment_type(int  type)
{
  if (type < CS_LAGR_MOMENT_MEAN || type > CS_LAGR_MOMENT_VARIANCE)
    bft_error(__FILE__, __LINE__,0,
              _("Out-of range statistics type: %d"),
              (int)type);
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Return number of possible statistics types.
 *
 * \return number of possible statistics types
 */
/*---------------------------------------------------------------------------*/

inline static int
_n_stat_types(void)
{
  return CS_LAGR_STAT_PARTICLE_ATTR + CS_LAGR_N_ATTRIBUTES;
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Create statistical weight
 *
 * \param[in]   class  statistical class id, or 0
 * \param[out]  name   resulting name
 */
/*---------------------------------------------------------------------------*/

static void
_statistical_weight_name(int   class,
                         char  name[64])
{
  char _class_name[12];

  _class_name[0] = '\0';

  if (class > 0)
    snprintf(_class_name, 12, "_c%d", class);

  size_t l0 =  strlen(_class_name);

  snprintf(name,
           63 - l0,
           "particle_cumulative_weight");

  name[63] = '\0';
  strcat(name, _class_name);
  name[63] = '\0';
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Create statistics name
 *
 * \param[in]   stat_type     particle statistics type
 * \param[in]   component_id  component id, or -1
 * \param[in]   class_id      statistical class id, or 0
 * \param[in]   moment_type   moment type
 * \param[out]  name          resulting name
 */
/*---------------------------------------------------------------------------*/

static void
_stat_name(int                    stat_type,
           int                    component_id,
           int                    class_id,
           cs_lagr_stat_moment_t  moment_type,
           char                   name[64])
{
  _check_moment_type(moment_type);

  char _class_name[12];
  char _comp_name[12];

  int attr_id = cs_lagr_stat_type_to_attr_id(stat_type);

  const char *type_name[2] = {"mean", "var"};

  _comp_name[0] = '\0';
  _class_name[0] = '\0';

  if (component_id > -1)
    snprintf(_comp_name, 12, "_l%d", component_id);

  if (class_id > 0)
    snprintf(_class_name, 12, "_c%d", class_id);

  size_t l0 =   strlen(_comp_name) + strlen(_class_name)
              + strlen(type_name[moment_type]);

  if (attr_id > -1) {
    snprintf(name,
             63 - l0,
             "%s_particle_%s",
             type_name[moment_type],
             cs_lagr_attribute_name[attr_id]);
    name[63] = '\0';
  }
  else {
    switch(stat_type) {
    case CS_LAGR_STAT_VOLUME_FRACTION:
      snprintf(name,
               63 - l0,
               "%s_particle_%s",
               type_name[moment_type],
               "volume_fraction");
      break;
    default:
      assert(0);
      break;
    }
  }

  name[63] = '\0';
  strcat(name, _comp_name);
  strcat(name, _class_name);
  name[63] = '\0';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to weight accumulator values.
 *
 * \param[in]   mwa  pointer to weight accumulator structure
 *
 * \return pointer to weight accumulator values
 */
/*----------------------------------------------------------------------------*/

static cs_real_t *
_mwa_val(cs_lagr_moment_wa_t  *mwa)
{
  assert(mwa != NULL);

  cs_real_t *val = mwa->val;

  if (mwa->f_id >= 0)
    val = cs_field_by_id(mwa->f_id)->val;
  else if (mwa->location_id == CS_MESH_LOCATION_NONE)
    val = &(mwa->val0);

  return val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to weight accumulator values.
 *
 * \param[in]   mwa  pointer to weight accumulator structure
 *
 * \return pointer to weight accumulator values
 */
/*----------------------------------------------------------------------------*/

static const cs_real_t *
_mwa_const_val(const cs_lagr_moment_wa_t  *mwa)
{
  assert(mwa != NULL);

  const cs_real_t *val = mwa->val;

  if (mwa->f_id >= 0)
    val = cs_field_by_id(mwa->f_id)->val;
  else if (mwa->location_id == CS_MESH_LOCATION_NONE)
    val = &(mwa->val0);

  return val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of elements associated with a weight accumulator.
 *
 * \return  number of elements associated with a weight accumulator.
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_n_w_elts(const cs_lagr_moment_wa_t  *mwa)
{
  cs_lnum_t n_w_elts = 1;

  if (mwa->location_id != CS_MESH_LOCATION_NONE)
    n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];

  return n_w_elts;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Statistics initialization
 *
 * Reset moment to zero
 * Called once for stationnary moments,
 *        at every step for non-stationnary statistics
 */
/*----------------------------------------------------------------------------*/

static void
_cs_lagr_moment_reset(cs_lagr_moment_t *mt)
{
  cs_field_t *f = cs_field_by_id(mt->f_id);

  cs_field_set_values(f, 0.);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize weight accumulator if required and reset to 0
 *
 * \param[in, out]  mwa  moment weight accumulator
 */
/*----------------------------------------------------------------------------*/

static void
_ensure_init_wa(cs_lagr_moment_wa_t  *mwa)
{
  if (   mwa->location_id != CS_MESH_LOCATION_NONE
      && mwa->val == NULL
      && mwa->f_id < 0) {

    cs_lnum_t n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
    BFT_MALLOC(mwa->val, n_w_elts, cs_real_t);
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      mwa->val[i] = 0.;

  }

  else if (cs_glob_lagr_time_scheme->isttio == 0) {

    cs_lnum_t n_w_elts = _n_w_elts(mwa);
    cs_real_t *val = _mwa_val(mwa);
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      val[i] = 0.;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize moment value if required
 *
 * \param[in, out]  mt  moment
 */
/*----------------------------------------------------------------------------*/

static void
_ensure_init_moment(cs_lagr_moment_t  *mt)
{
  assert(mt->f_id > 0);

  cs_field_t *f = cs_field_by_id(mt->f_id);

  if (f->vals[0] == NULL)
    cs_field_allocate_values(f);

  else if (cs_glob_lagr_time_scheme->isttio == 0)
    _cs_lagr_moment_reset(mt);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to time step values.
 *
 * \return pointer to time step values
 */
/*----------------------------------------------------------------------------*/

static const cs_real_t *
_dt_val(void)
{
  const cs_real_t *dt_val;
  const cs_field_t *f = cs_field_by_name_try("dt");

  if (_p_dt != NULL)
    dt_val = _p_dt;
  else if (f != NULL)
    dt_val = f->val;
  else
    dt_val = &(cs_glob_time_step_options->dtref);

  return dt_val;
}

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
 *   wa_location_id <-- associated weigh accumulator mesh location id
 *   dim            <-- dimension associated with moment
 *   moment_type    <-- moment type
 *   stat_type      <-- particle statistics type, or -1
 *   nt_start       <-> starting time step
 *   t_start        <-> starting time
 *   restart_mode   <-- behavior in case of restart (reset, automatic, strict)
 *   restart_name   <-- if not NULL, previous name in case of restart
 *
 *----------------------------------------------------------------------------*/

static int
_check_restart(const char                     *name,
               const cs_time_step_t           *ts,
               cs_lagr_moment_restart_info_t  *ri,
               int                             location_id,
               int                             wa_location_id,
               int                             dim,
               cs_lagr_stat_moment_t           moment_type,
               int                             stat_type,
               int                             class_id,
               int                            *nt_start,
               double                         *t_start,
               cs_lagr_stat_restart_t          restart_mode)
{
  int i;
  int prev_id = -1;
  int prev_wa_id = -1;

  if (   (*nt_start > -1 && *nt_start >= ri->nt_prev)
      || (*t_start >= 0  && *t_start >= ri->t_prev))
    return prev_id;

  /* Adjust accumulator info if moment should be restarted */

  if (restart_mode == CS_LAGR_MOMENT_RESTART_RESET) {
    *nt_start = ri->nt_prev + 1;
    *t_start = ri->t_prev;
    return prev_id;
  }

  /* If we reach here, restart info is required */

  /* Find matching restart data if required, and adjust accumulator
     info if moment should be restarted, or if available data does
     not match and we do not require exact mode. */

  char _r_name[64];
  strncpy(_r_name, name, 63);
  _r_name[63] = '\0';

  for (i = 0; i < ri->n_moments; i++) {

    if (strcmp(ri->name[i], _r_name) == 0) {

      bool matching_restart = true;
      prev_id = i;
      prev_wa_id = ri->wa_id[i];

      if (   ri->wa_location_id[prev_wa_id] != wa_location_id
          || ri->m_type[i] != (int)moment_type
          || ri->location_id[i] != location_id
          || ri->stat_type[i]  != stat_type
          || ri->class[i] != class_id
          || ri->dimension[i] != dim)
        matching_restart = false;

      if (   restart_mode == CS_LAGR_MOMENT_RESTART_EXACT
          && (   ri->wa_nt_start[prev_wa_id] != *nt_start
              || (   !ts->is_local
                  && fabs(ri->wa_t_start[prev_wa_id] - *t_start) > 1.e-18)))
        matching_restart = false;

      if (matching_restart == false) {

        bft_printf(_("\nRestart data for particle statistics \"%s\"\n"
                     " (previously \"%s\") does not match.\n"
                     "  previous values:\n"
                     "    weight accumulator location_id: %d\n"
                     "    moment_type:                           %d\n"
                     "    location_id:                    %d\n"
                     "    dimension:                      %d\n"
                     "    start time step:                %d\n"
                     "    start time:                     %12.5e\n"),
                   name, _r_name, ri->wa_location_id[prev_wa_id],
                   ri->m_type[i], ri->location_id[i], ri->dimension[i],
                   ri->wa_nt_start[prev_wa_id],
                   ri->wa_t_start[prev_wa_id]);

        if (restart_mode == CS_LAGR_MOMENT_RESTART_AUTO) {

          bft_printf
            (_("\nWarning: computation of time moment \"%s\""
               " will be reset,\n"
               "         as restart data for \"%s\" does not match.\n"),
             name, _r_name);
          *nt_start = ri->nt_prev + 1;
          *t_start = ri->t_prev;

        }
        else if (restart_mode == CS_LAGR_MOMENT_RESTART_EXACT)

          bft_error(__FILE__, __LINE__, 0,
                    _("Restart data for time moment \"%s\"\n"
                      " (previously \"%s\") does not match."),
                    name, _r_name);

      }

      if (matching_restart == false)
        prev_id = -1;
      /* else { */
      /*   *nt_start = ri->wa_nt_start[prev_wa_id]; */
      /*   *t_start = ri->wa_t_start[prev_wa_id]; */
      /* } */

      break;
    }

  }

  if (i >= ri->n_moments) {
    if (restart_mode == CS_LAGR_MOMENT_RESTART_AUTO) {
      bft_printf
        (_("\nWarning: computation of time moment \"%s\""
           "will be reset,\n"
           "           as restart data for \"%s\" is not available.\n"),
         name, _r_name);
      prev_id = -1;
      *nt_start = ri->nt_prev + 1;
      *t_start = ri->t_prev;
    }
    else if (restart_mode == CS_LAGR_MOMENT_RESTART_EXACT)
      bft_error(__FILE__, __LINE__, 0,
                _("Restart data for time moment \"%s\"\n"
                  " (previously \"%s\") not available."),
                name, _r_name);

  }

  /* Also check for presence of sub-moment restart info in case of
     higer order restart */

  if (prev_id > -1 && moment_type == CS_LAGR_MOMENT_VARIANCE) {

      cs_lagr_stat_moment_t s_m_type = CS_LAGR_MOMENT_MEAN;
      int l_dim = (dim == 6) ? 3 : dim;

      int l_id = ri->l_id[prev_id];

      if (   ri->wa_id[l_id] != prev_wa_id
          || ri->m_type[l_id] != (int)s_m_type
          || ri->location_id[l_id] != location_id
          || ri->dimension[l_id] != l_dim)
        bft_error(__FILE__, __LINE__, 0,
                  _("Restart data for time moment \"%s\"\n"
                    " (previously \"%s\") seems inconsistent:\n"
                    "   lower order moment of moment_type %s was \"%s\",\n"
                    "   but has non-matching attributes:\n"
                    "    weight accumulator id: %d (expected %d)\n"
                    "    moment_type:           %d\n"
                    "    location_id:           %d\n"
                    "    dimension:             %d\n"),
                  name, _r_name, cs_lagr_moment_type_name[s_m_type],
                  ri->name[l_id], ri->wa_id[l_id], prev_wa_id,
                  ri->m_type[l_id], ri->location_id[l_id], ri->dimension[l_id]);

  }

  /* Return previous moment id */

  return prev_id;
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
  int sizes[3];
  int retcode;

  const cs_time_step_t  *ts = cs_glob_time_step;
  int itysup;

  /* Step 1: check version */
  {
    itysup = 0;

    cs_int_t ivers;

    retcode = cs_restart_read_section
                (r, "version_fichier_suite_Lagrangien_statistiques",
                 itysup, 1, CS_TYPE_cs_int_t, &ivers);
  }

  /* Now read main metadata */

  BFT_MALLOC(_restart_info, 1, cs_lagr_moment_restart_info_t);

  cs_lagr_moment_restart_info_t  *ri = _restart_info;

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:sizes",
                                    CS_MESH_LOCATION_NONE,
                                    3,
                                    CS_TYPE_cs_int_t,
                                    sizes);
  ri->nt_prev = ts->nt_prev;
  ri->t_prev = ts->t_prev;

  ri->n_wa = sizes[0];
  ri->n_moments = sizes[1];

  BFT_MALLOC(ri->name, ri->n_moments, const char*);
  BFT_MALLOC(ri->name_buf, sizes[2] + 1, char);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:names",
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

  BFT_MALLOC(ri->wa_location_id, ri->n_wa, int);
  BFT_MALLOC(ri->wa_nt_start, ri->n_wa, int);
  BFT_MALLOC(ri->wa_t_start, ri->n_wa, cs_real_t);

  cs_restart_read_section(r,
                          "lagr_stats:wa:location_id",
                          CS_MESH_LOCATION_NONE,
                          ri->n_wa,
                          CS_TYPE_cs_int_t,
                          ri->wa_location_id);
  _assert_restart_success(retcode);

  cs_restart_read_section(r,
                          "lagr_stats:wa:nt_start",
                          CS_MESH_LOCATION_NONE,
                          ri->n_wa,
                          CS_TYPE_cs_int_t,
                          ri->wa_nt_start);
  _assert_restart_success(retcode);

  cs_restart_read_section(r,
                          "lagr_stats:wa:t_start",
                          CS_MESH_LOCATION_NONE,
                          ri->n_wa,
                          CS_TYPE_cs_real_t,
                          ri->wa_t_start);
  _assert_restart_success(retcode);

  /* Information on moments proper */

  BFT_MALLOC(ri->m_type, ri->n_moments, int);
  BFT_MALLOC(ri->class, ri->n_moments, int);
  BFT_MALLOC(ri->location_id, ri->n_moments, int);
  BFT_MALLOC(ri->dimension, ri->n_moments, int);
  BFT_MALLOC(ri->wa_id, ri->n_moments, int);
  BFT_MALLOC(ri->l_id, ri->n_moments, int);
  BFT_MALLOC(ri->stat_type, ri->n_moments, int);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:type",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->m_type);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:class",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->class);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:location_id",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->location_id);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:dimension",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->dimension);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:wa_id",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->wa_id);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:lower_order_id",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->l_id);
  _assert_restart_success(retcode);

  retcode = cs_restart_read_section(r,
                                    "lagr_stats:stat_type",
                                    CS_MESH_LOCATION_NONE,
                                    ri->n_moments,
                                    CS_TYPE_cs_int_t,
                                    ri->stat_type);
  _assert_restart_success(retcode);
}

/*----------------------------------------------------------------------------
 * Read restart metadata.
 *----------------------------------------------------------------------------*/

static void
_restart_info_read(void)
{
  const cs_time_step_t  *ts = cs_glob_time_step;

  if (   ts->nt_prev < 1
      || !cs_file_isreg("restart/lagrangian_stats")
      || cs_glob_lagr_stat_options->isuist < 1) {
    _restart_info_checked = true;
    return;
  }

  cs_restart_t *r = NULL;

  /* Read previous time step if not already done */

  if (ts->nt_prev < 1) {
    r = cs_restart_create("main", "restart", CS_RESTART_MODE_READ);
    cs_restart_read_time_step_info(r);
    cs_restart_destroy(&r);
  }

  /* Now read lagr-moment specific data */

  r = cs_restart_create("lagrangian_stats", NULL, CS_RESTART_MODE_READ);

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
  cs_lagr_moment_restart_info_t  *ri = _restart_info;

  if (ri != NULL) {

    BFT_FREE(ri->l_id);
    BFT_FREE(ri->wa_id);
    BFT_FREE(ri->dimension);
    BFT_FREE(ri->location_id);
    BFT_FREE(ri->m_type);
    BFT_FREE(ri->class);

    BFT_FREE(ri->wa_t_start);
    BFT_FREE(ri->wa_nt_start);
    BFT_FREE(ri->wa_location_id);

    BFT_FREE(ri->name_buf);
    BFT_FREE(ri->name);

    BFT_FREE(ri);

    _restart_info = ri;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read restart moment data
 */
/*----------------------------------------------------------------------------*/

static void
_cs_lagr_moment_restart_read(void)
{
  int retcode;

  /* Initialize */

  const cs_time_step_t  *ts = cs_glob_time_step;
  _t_prev_iter = ts->t_prev;

  static cs_restart_t  *cs_lag_stat_restart = NULL;

  char const *ficsui = "lagrangian_stats";
  cs_lag_stat_restart = cs_restart_create(ficsui, NULL, CS_RESTART_MODE_READ);
  if (cs_lag_stat_restart == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error opening Lagrangian statistics restart file.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              ficsui);

  if (_restart_info == NULL)
    _restart_info_read_auxiliary(cs_lag_stat_restart);

  if (_restart_info == NULL) {
    cs_restart_destroy(&cs_lag_stat_restart);
    return;
  }

  cs_lagr_moment_restart_info_t  *ri = _restart_info;

  /* Read information proper */

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (mwa->restart_id > -1 && mwa->location_id > CS_MESH_LOCATION_NONE) {
      char s[64];
      snprintf(s, 64, "lagr_stats:wa:%02d:val", mwa->restart_id);
      _ensure_init_wa(mwa);
      retcode = cs_restart_read_section(cs_lag_stat_restart,
                                        s,
                                        mwa->location_id,
                                        1,
                                        CS_TYPE_cs_real_t,
                                        _mwa_val(mwa));
      _assert_restart_success(retcode);
    }
  }

  for (int i = 0; i < _n_lagr_stats; i++) {
    cs_lagr_moment_t *mt = _lagr_stats + i;
    if (mt->restart_id > -1) {
      _ensure_init_moment(mt);
      cs_real_t *val = mt->val;
      if (mt->f_id > -1) {
        cs_field_t *f = cs_field_by_id(mt->f_id);
        val = f->val;
      }
      retcode = cs_restart_read_section(cs_lag_stat_restart,
                                        ri->name[mt->restart_id],
                                        mt->location_id,
                                        mt->dim,
                                        CS_TYPE_cs_real_t,
                                        val);
      _assert_restart_success(retcode);
    }
  }

  _restart_info_checked = true;

  cs_restart_destroy(&cs_lag_stat_restart);
}

/*----------------------------------------------------------------------------
 * Initializa particle attribute mapping to statistics variable number
 *----------------------------------------------------------------------------*/

static void
_init_vars_attribute(void)
{
  if (cs_glob_lagr_model->physical_model == 1) {
    if (cs_glob_lagr_specific_physics->idpvar)
      cs_lagr_stat_activate_attr(CS_LAGR_DIAMETER);
    if (cs_glob_lagr_specific_physics->impvar)
      cs_lagr_stat_activate_attr(CS_LAGR_MASS);
    if (cs_glob_lagr_specific_physics->itpvar)
      cs_lagr_stat_activate_attr(CS_LAGR_TEMPERATURE);
  }

  else if (cs_glob_lagr_model->physical_model == 2) {
    cs_lagr_stat_activate_attr(CS_LAGR_MASS);
    cs_lagr_stat_activate_attr(CS_LAGR_WATER_MASS);
    cs_lagr_stat_activate_attr(CS_LAGR_COAL_MASS);
    cs_lagr_stat_activate_attr(CS_LAGR_COKE_MASS);
    cs_lagr_stat_activate_attr(CS_LAGR_TEMPERATURE);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add or find moment weight and time accumulator.
 *
 * If no data function is provided, a constant weight of 1 is assumed
 * (this weight will be multiplied by the time step).
 *
 * Note that if the data_input associated with a data_func pointer is not
 * NULL, the lifecycle of the data pointed to must be handled separately
 * (and the pointer must remain valid throughout the time moment updates).
 *
 * \param[in]  m_type       associated moment type
 * \param[in]  p_data_func  particle-based function used to define data
 *                          values, or NULL
 * \param[in]  m_data_func  mesh-based function used to define data
 *                          values, or NULL
 * \param[in]  data_input   pointer to optional (untyped) value or structure
 * \param[in]  class_id     statistical class number
 * \param[in]  location_id  associated mesh location id
 * \param[in]  nt_start     starting time step
 * \param[in]  t_start      starting time
 * \param[in]  prev_wa_id   previous weight accumulator id, or -1
 *
 * \return id of matching time accumulator
 */
/*----------------------------------------------------------------------------*/

static int
_find_or_add_wa(cs_lagr_moment_p_data_t  *p_data_func,
                cs_lagr_moment_m_data_t  *m_data_func,
                const void               *data_input,
                int                       class_id,
                int                       location_id,
                int                       nt_start,
                double                    t_start,
                int                       prev_wa_id)
{
  int wa_id = -1;
  int _nt_start = nt_start;
  double _t_start = t_start;

  cs_lagr_moment_wa_t *mwa = NULL;

  /* Reduce number of possible options */

  if (_nt_start < 0)
    _nt_start = -1;

  if (_t_start < 0. && _nt_start < 0)
    _nt_start = 0;

  if (nt_start >= 0)
    _t_start = -1.;

  /* Check if this accumulator is already defined */

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    mwa = _lagr_stats_wa + i;
    if (   nt_start == mwa->nt_start && fabs(mwa->t_start - _t_start) < 1e-18
        && prev_wa_id == mwa->restart_id && class_id == mwa->class
        && p_data_func == mwa->p_data_func && m_data_func == mwa->m_data_func
        && data_input == mwa->data_input)
      return i;
  }

  /* If we did not return yet, a new structure must be added */

  /* Reallocate if necessary */

  if (_n_lagr_stats_wa + 1 > _n_lagr_stats_wa_max) {
    if (_n_lagr_stats_wa_max < 1)
      _n_lagr_stats_wa_max = 2;
    else
      _n_lagr_stats_wa_max *= 2;
    BFT_REALLOC(_lagr_stats_wa, _n_lagr_stats_wa_max, cs_lagr_moment_wa_t);
  }

  /* Now initialize members */

  wa_id = _n_lagr_stats_wa;

  mwa = _lagr_stats_wa + _n_lagr_stats_wa;

  _n_lagr_stats_wa += 1;

  mwa->restart_id = prev_wa_id;

  mwa->f_id = -1;

  mwa->nt_start = _nt_start;
  mwa->t_start = _t_start;

  mwa->location_id = location_id;

  mwa->class = class_id;

  mwa->p_data_func = p_data_func;
  mwa->m_data_func = m_data_func;
  mwa->data_input = data_input;

  /* Create field in specific case of statistical weight */

  if (   location_id == CS_MESH_LOCATION_CELLS
      && p_data_func == NULL
      && m_data_func == NULL
      && nt_start == 0) {

    char name[64];
    _statistical_weight_name(class_id, name);

    if (cs_field_by_name_try(name) == NULL) {
      cs_field_t *f
        = cs_field_create(name,
                          CS_FIELD_POSTPROCESS | CS_FIELD_ACCUMULATOR,
                          location_id,
                          1,
                          false);
      mwa->f_id = f->id;
    }

  }

  mwa->val0 = 0.;
  mwa->val = NULL;

  /* Structure is now initialized */

  return wa_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and associate a field to a moment
 *
 * \param[in]  name         field name
 * \param[in]  location_id  mesh location id
 * \param[in]  dim          field dimension
 *
 * \return associated field
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_cs_lagr_moment_associate_field(const char  *name,
                                int          location_id,
                                int          dim)
{
  cs_field_t *f
    = cs_field_find_or_create(name,
                              CS_FIELD_POSTPROCESS | CS_FIELD_ACCUMULATOR,
                              location_id,
                              dim);

  /* cs_field_allocate_values(f); */
  const int log_key_id = cs_field_key_id("log");
  cs_field_set_key_int(f, log_key_id, 1);
  const int vis_key_id = cs_field_key_id("post_vis");
  cs_field_set_key_int(f, vis_key_id, 1);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add or find moment structure.                 .
 *
 * \param[in]  location_id id  of associated mesh location
 * \param[in]  component_id    attribute component id, or < 0 for all
 * \param[in]  class_id        statistical class
 * \param[in]  stat_type       statistics type id, or -1
 * \param[in]  dim             dimension associated with element data
 * \param[in]  p_data_func     particle-based function used to define data
 *                             values, or NULL
 * \param[in]  m_data_func     mesh-based function used to define data
 *                             values, or NULL
 * \param[in]  data_input      pointer to optional value or structure
 *                             to be used by data_func
 * \param[in]  m_type          moment type, mean or variance
 * \param[in]  wa_id           weight accumulator id
 * \param[in]  prev_id         restart moment id
 *
 * \return
 *   id of matching moment
 */
/*----------------------------------------------------------------------------*/

static int
_find_or_add_moment(int                       location_id,
                    int                       component_id,
                    int                       class_id,
                    int                       stat_type,
                    int                       dim,
                    cs_lagr_moment_p_data_t  *p_data_func,
                    cs_lagr_moment_m_data_t  *m_data_func,
                    const void               *data_input,
                    cs_lagr_stat_moment_t     m_type,
                    int                       wa_id,
                    int                       prev_id)
{
  cs_lagr_moment_t *mt = NULL;
  int _dim = (dim == 3 && m_type == CS_LAGR_MOMENT_VARIANCE) ? 6 : dim;

  int moment_id = -1;

  /* Check if this moment is already defined;
     ignore associated field at this stage, as a moment defined automatically
     to satisfy a dependency (i.e. a mean for a variance) may not have an
     associated field. */

  for (int i = 0; i < _n_lagr_stats; i++) {

    mt = _lagr_stats + i;

    if (   location_id  == mt->location_id
        && component_id == mt->component_id
        && stat_type    == mt->stat_type
        && _dim         == mt->dim
        && dim          == mt->data_dim
        && p_data_func  == mt->p_data_func
        && m_data_func  == mt->m_data_func
        && data_input   == mt->data_input
        && m_type       == mt->m_type
        && wa_id        == mt->wa_id
        && class_id     == mt->class
        && prev_id      == mt->restart_id)
      return i;

  }

  /* If we did not return yet, a new structure must be added */

  /* Reallocate if necessary */

  if (_n_lagr_stats + 1 > _n_lagr_stats_max) {

    if (_n_lagr_stats_max < 1)
      _n_lagr_stats_max = 2;

    else
      _n_lagr_stats_max *= 2;

    BFT_REALLOC(_lagr_stats, _n_lagr_stats_max, cs_lagr_moment_t);

  }

  /* Now define moment */

  moment_id = _n_lagr_stats;
  _n_lagr_stats += 1;

  mt = _lagr_stats + moment_id;

  mt->m_type = m_type;
  mt->restart_id = prev_id;
  mt->wa_id = wa_id;
  mt->f_id = -1;

  mt->dim = _dim;
  mt->data_dim = dim;
  mt->location_id = location_id;

  mt->p_data_func = p_data_func;
  mt->m_data_func = m_data_func;
  mt->data_input = data_input;

  mt->l_id = -1;
  mt->stat_type = stat_type;
  mt->class = class_id;
  mt->component_id = component_id;

  mt->name = NULL;
  mt->val = NULL;

  mt->nt_cur = -1;

  return moment_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update weight accumulator for a mesh-based weight array
 *
 * \param[in, out]   mwa  moment weight accumulator
 * \param[in]        w    pointer to current weight values
 */
/*----------------------------------------------------------------------------*/

static void
_update_wa_m(cs_lagr_moment_wa_t  *mwa,
             cs_real_t            *restrict w)
{
  assert(w != NULL);

  if (mwa->location_id == CS_MESH_LOCATION_NONE)
    mwa->val0 += w[0];
  else {
    cs_lnum_t n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      mwa->val[i] += w[i];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Vector-vector multiplication, where the multiplying vector is
 *        defined on cells, and the multiplied vector defined on a given
 *        mesh location.
 *
 * \param[in]      location_id_y  location id for y
 * \param[in]      x              multiplier array
 * \param[in, out] y              multiplied array
 */
/*----------------------------------------------------------------------------*/

static void
_vv_mesh_location_cells(int              location_id_y,
                        const cs_real_t  x[],
                        cs_real_t        y[])
{
  const cs_mesh_location_type_t loc_type
    = cs_mesh_location_get_type(location_id_y);
  const cs_lnum_t *elt_list
    = cs_mesh_location_get_elt_list(location_id_y);
  const cs_mesh_t *mesh = cs_glob_mesh;

  cs_lnum_t n_y_elts = cs_mesh_location_get_n_elts(location_id_y)[0];
  assert(location_id_y != CS_MESH_LOCATION_NONE);

  switch(loc_type) {
  case CS_MESH_LOCATION_CELLS:
    {
      if (elt_list == NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_y_elts; c_id++)
          y[c_id] *= x[c_id];
      }
      else {
        for (cs_lnum_t i = 0; i < n_y_elts; i++) {
          cs_lnum_t c_id = elt_list[i];
          y[i] *= x[c_id];
        }
      }
    }
    break;
  case CS_MESH_LOCATION_INTERIOR_FACES:
    {
      const cs_lnum_2_t *i_face_cells
        = (const cs_lnum_2_t *)mesh->i_face_cells;
      if (elt_list == NULL) {
        for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
          cs_lnum_t c_id_0 = i_face_cells[f_id][0];
            cs_lnum_t c_id_1 = i_face_cells[f_id][1];
            y[f_id] *= (x[c_id_0] + x[c_id_1]) * 0.5;
        }
      }
      else {
        for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++) {
          cs_lnum_t f_id = elt_list[i];
          cs_lnum_t c_id_0 = i_face_cells[f_id][0];
          cs_lnum_t c_id_1 = i_face_cells[f_id][1];
          y[i] *= (x[c_id_0] + x[c_id_1]) * 0.5;
        }
      }
    }
    break;
  case CS_MESH_LOCATION_BOUNDARY_FACES:
    {
      const cs_lnum_t *b_face_cells = (const cs_lnum_t *)mesh->b_face_cells;
      if (elt_list == NULL) {
        for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
          cs_lnum_t c_id = b_face_cells[f_id];
          y[f_id] *= x[c_id];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
          cs_lnum_t f_id = elt_list[i];
          cs_lnum_t c_id = b_face_cells[f_id];
          y[i] *= x[c_id];
        }
      }
    }
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Multiplication for mesh locations of type:\n"
                "%s is not currently supported."),
              cs_mesh_location_type_name[loc_type]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute current weight if weight is based on a mesh (not particle)
 *        based function
 *
 * When applicable, this function either returns a pointer to an allocated
 * array, or to w0. If the returned value is different from w0
 * (i.e. allocated), the caller is responsible for freeing it.
 *
 * When not applicable (i.e. when no mesh-based weight computation function
 * is defined), NULL is returned.
 *
 * \parameter[in, out]  mwa       moment weight accumulator
 * \parameter[in]       dt        cell time step values
 * \parameter[in, out]  w0        pointer to buffer in case weight values
 *                                is of size 1
 *
 * \return  pointer to weight array (w0 or allocated array), or NULL
 */
/*----------------------------------------------------------------------------*/

static cs_real_t *
_compute_current_weight_m(cs_lagr_moment_wa_t  *mwa,
                          const cs_real_t      *restrict dt,
                          cs_real_t             w0[1])
{
  cs_real_t *w = NULL;

  if (mwa->m_data_func == NULL)
    return w;

  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_lnum_t        n_w_elts = _n_w_elts(mwa);

  assert(mwa->nt_start <= ts->nt_cur);

  if (n_w_elts == 1)
    w = w0;
  else {
    BFT_MALLOC(w, n_w_elts, cs_real_t);
  }

  /* Base weight */

  mwa->m_data_func(mwa->data_input, mwa->location_id, mwa->class, w);

  /* Multiply time step */

  if (! ts->is_local) {
    double _dt;
    if (mwa->nt_start == ts->nt_cur)
      _dt = ts->t_cur - mwa->t_start;
    else
      _dt = dt[0];
    for (cs_lnum_t i = 0; i < n_w_elts; i++)
      w[i] *= _dt;
  }

  else
    _vv_mesh_location_cells(mwa->location_id, dt, w);

  return w;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset unsteady stats
 */
/*----------------------------------------------------------------------------*/

static void
_cs_lagr_moment_reset_unsteady_stats(void)
{
  const cs_time_step_t  *ts = cs_glob_time_step;

  /* reset moment values*/

  for (int i = 0; i < _n_lagr_stats; i++) {

    cs_lagr_moment_t *mt = _lagr_stats + i;
    _cs_lagr_moment_reset(mt);

  }

  /* reset weight accumulator values*/

  for (int i = 0; i < _n_lagr_stats_wa; i++) {

    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;

    if (mwa->nt_start > 0 && mwa->nt_start <= ts->nt_cur) {

      mwa->nt_start = ts->nt_cur;
      mwa->t_start = ts->t_prev;

      mwa->val0 = 0;

      cs_real_t *val = _mwa_val(mwa);

      if (val != NULL) {
        cs_lnum_t n_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];
        for (cs_lnum_t j = 0; j < n_elts; j++)
          val[j] = 0.;
      }

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update a given mesh based date function moment.
 *
 * \param[in, out]  mt      pointer to associated moment
 * \param[in]       mwa     pointer to associated weight accumulator
 * \param[in]       w       weight values for current time step
 * \param[in]       nt_cur
 */
/*----------------------------------------------------------------------------*/

static void
_cs_lagr_stat_update_mesh_moment(cs_lagr_moment_t           *mt,
                                 const cs_lagr_moment_wa_t  *mwa,
                                 const cs_real_t            *restrict w,
                                 int                         nt_cur)
{
  /* Return if moment already updated */

  if (mt->nt_cur >= nt_cur)
    return;

  /* Dimensions and previous accumulation */

  const cs_lnum_t n_w_elts = _n_w_elts(mwa);
  const cs_lnum_t n_elts   = cs_mesh_location_get_n_elts(mt->location_id)[0];

  const cs_lnum_t nd = n_elts * mt->dim;

  const cs_real_t  *restrict wa_sum = _mwa_const_val(mwa);

  /* Current and accumulated weight */

  cs_lnum_t  wa_stride = (n_w_elts > 1) ? 1 : 0;

  cs_real_t *restrict x;
  BFT_MALLOC(x, nd, cs_real_t);

  mt->m_data_func(mt->data_input, mt->location_id, mt->class, x);

  cs_real_t *restrict val = mt->val;
  if (mt->f_id > -1) {
    cs_field_t *f = cs_field_by_id(mt->f_id);
    val = f->val;
  }

  if (mt->m_type == CS_LAGR_MOMENT_VARIANCE) {

    assert(mt->l_id > -1);

    cs_lagr_moment_t *mt_mean = _lagr_stats + mt->l_id;

    _ensure_init_moment(mt_mean);
    cs_real_t *restrict m = mt_mean->val;
    if (mt_mean->f_id > -1) {
      cs_field_t *f_mean = cs_field_by_id(mt_mean->f_id);
      m = f_mean->val;
    }

    if (mt->dim == 6) { /* variance-covariance matrix */
      assert(mt->data_dim == 3);
      for (cs_lnum_t je = 0; je < n_elts; je++) {
        double delta[3], delta_n[3], r[3], m_n[3];
        const cs_lnum_t k = je*wa_stride;
        const double wa_sum_n = w[k] + wa_sum[k];
        for (cs_lnum_t l = 0; l < 3; l++) {
          cs_lnum_t jl = je*6 + l, jml = je*3 + l;
          delta[l]   = x[jml] - m[jml];
          r[l] = delta[l] * (w[k] / (fmax(wa_sum_n, 1e-100)));
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
      }
    }

    else { /* simple variance */
      for (cs_lnum_t j = 0; j < nd; j++) {
        const cs_lnum_t k = (j*wa_stride) / mt->dim;
        double wa_sum_n = w[k] + wa_sum[k];
        double delta = x[j] - m[j];
        double r = delta * (w[k] / (fmax(wa_sum_n, 1e-100)));
        double m_n = m[j] + r;
        val[j] = (val[j]*wa_sum[k] + (w[k]*delta*(x[j]-m_n))) / wa_sum_n;
        m[j] += r;
      }
    }

    mt_mean->nt_cur = nt_cur;
  }

  else if (mt->m_type == CS_LAGR_MOMENT_MEAN) {

    for (cs_lnum_t j = 0; j < nd; j++) {
      const cs_lnum_t k = (j*wa_stride) / mt->dim;
      val[j] += (x[j] - val[j]) * (w[k] / (w[k] + wa_sum[k]));
    }

  }

  mt->nt_cur = nt_cur;

  BFT_FREE(x);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update all moment accumulators.
 */
/*----------------------------------------------------------------------------*/

static void
_cs_lagr_stat_update_all(void)
{
  const cs_time_step_t  *ts = cs_glob_time_step;
  cs_lagr_particle_set_t *p_set = cs_lagr_get_particle_set();
  const cs_real_t *dt_val = _dt_val();

  _t_prev_iter = ts->t_prev;

  /* Outer loop in weight accumulators, to avoid recomputing weights
     too many times */

  for (int wa_id = 0; wa_id < _n_lagr_stats_wa; wa_id++) {

    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + wa_id;

    /* Check if accumulator and associated moments are active */

    if (mwa->nt_start == 0 && cs_glob_lagr_stat_options->idstnt <= ts->nt_cur) {
      mwa->nt_start = ts->nt_cur;
      mwa->t_start = _t_prev_iter;
    }
    else if (mwa->t_start < 0. && mwa->nt_start <= ts->nt_cur)
      mwa->t_start = _t_prev_iter;
    else if (mwa->nt_start < 0 && mwa->t_start <= ts->t_cur)
      mwa->nt_start = ts->nt_cur;

    if (mwa->nt_start > ts->nt_cur)
      continue;

    /* Here, only active accumulators are considered */

    _ensure_init_wa(mwa);
    cs_real_t *g_wa_sum = _mwa_val(mwa);

    const cs_lnum_t n_w_elts = cs_mesh_location_get_n_elts(mwa->location_id)[0];

    /* Local weight array allocation */

    cs_real_t *restrict l_wa_sum = NULL;

    /* Compute mesh-based weight now if applicable
       (possibly sharing it across moments) */

    cs_real_t m_w0[1];
    cs_real_t *restrict m_weight = _compute_current_weight_m(mwa, dt_val, m_w0);

    if (m_weight == NULL)
      BFT_MALLOC(l_wa_sum, n_w_elts, cs_real_t);

    /* Loop on variances first, then means */

    for (int m_type = CS_LAGR_MOMENT_VARIANCE;
         m_type >= (int)CS_LAGR_MOMENT_MEAN;
         m_type--) {

      for (int i = 0; i < _n_lagr_stats; i++) {

        cs_lagr_moment_t *mt = _lagr_stats + i;

        if (   (int)mt->m_type == m_type
            && mt->wa_id == wa_id
            && mwa->nt_start > -1
            && mwa->nt_start <= ts->nt_cur
            && mt->nt_cur < ts->nt_cur) {

          int attr_id = cs_lagr_stat_type_to_attr_id(mt->stat_type);

          _ensure_init_moment(mt);

          /* Copy weight sum content to a local array
             for every new moment inside the current class */

          for (cs_lnum_t j = 0; j < n_w_elts; j++)
            l_wa_sum[j]= g_wa_sum[j];

          /* Case where data is particle-based */
          /*-----------------------------------*/

          if (mt->m_data_func == NULL) {

            cs_real_t *restrict val = mt->val;

            if (mt->f_id > -1) {
              cs_field_t *f = cs_field_by_id(mt->f_id);
              val = f->val;
            }

            /* prepare submoment definition */

            cs_lagr_moment_t *mt_mean = NULL;
            cs_real_t *restrict mean_val = NULL;

            /* Check if lower moment is defined and attached */

            if (mt->m_type == CS_LAGR_MOMENT_VARIANCE) {

              assert(mt->l_id > -1);
              mt_mean = _lagr_stats + mt->l_id;
              _ensure_init_moment(mt_mean);

              mean_val = mt_mean->val;

              if (mt_mean->f_id > -1) {
                cs_field_t *f_mean = cs_field_by_id(mt_mean->f_id);
                mean_val = f_mean->val;
              }

            }

            cs_real_t *pval = NULL;
            if (mt->p_data_func != NULL)
              BFT_MALLOC(pval, mt->data_dim, cs_real_t);

            for (cs_lnum_t part = 0; part < p_set->n_particles; part++) {

              unsigned char *particle
                = p_set->p_buffer + p_set->p_am->extents * part;

              cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle,
                                                               p_set->p_am);

              int p_class = 0;
              if (p_set->p_am->displ[0][CS_LAGR_STAT_CLASS] > 0)
                p_class = cs_lagr_particle_get_lnum(particle,
                                                    p_set->p_am,
                                                    CS_LAGR_STAT_CLASS);

              if (cell_id >= 0 && (p_class == mt->class || mt->class == 0)) {

                /* weight associated to current particle */

                cs_real_t p_weight;

                if (mwa->p_data_func == NULL)
                  p_weight = cs_lagr_particle_get_real(particle,
                                                       p_set->p_am,
                                                       CS_LAGR_STAT_WEIGHT);
                else
                  mwa->p_data_func(mwa->data_input,
                                   particle,
                                   p_set->p_am,
                                   &p_weight);
                p_weight *= dt_val[cell_id];

                if (mt->p_data_func == NULL)
                  pval = cs_lagr_particle_attr(particle, p_set->p_am, attr_id);
                else
                  mt->p_data_func(mt->data_input, particle, p_set->p_am, pval);

                /* update weight sum with new particle weight */
                const cs_real_t wa_sum_n = p_weight + l_wa_sum[cell_id];

                if (mt->m_type == CS_LAGR_MOMENT_VARIANCE) {

                  if (mt->dim == 6) { /* variance-covariance matrix */

                    assert(mt->data_dim == 3);

                    double delta[3], delta_n[3], r[3], m_n[3];

                    for (int l = 0; l < 3; l++) {

                      cs_lnum_t jl = cell_id*6 + l;
                      cs_lnum_t jml = cell_id*3 + l;
                      delta[l]   = pval[l] - mean_val[jml];
                      r[l] = delta[l] * (p_weight / (fmax(wa_sum_n, 1e-100)));
                      m_n[l] = mean_val[jml] + r[l];
                      delta_n[l] = pval[l] - m_n[l];
                      val[jl] = (  val[jl]*l_wa_sum[cell_id]
                                 + p_weight*delta[l]*delta_n[l]) / wa_sum_n;

                    }

                    /* Covariance terms.
                       Note we could have a symmetric formula using
                       0.5*(delta[i]*delta_n[j] + delta[j]*delta_n[i])
                       instead of
                       delta[i]*delta_n[j]
                       but unit tests in cs_moment_test.c do not seem to favor
                       one variant over the other; we use the simplest one.  */

                    cs_lnum_t j3 = cell_id*6 + 3,
                              j4 = cell_id*6 + 4,
                              j5 = cell_id*6 + 5;

                    val[j3] = (  val[j3]*l_wa_sum[cell_id]
                               + p_weight*delta[0]*delta_n[1]) / wa_sum_n;
                    val[j4] = (  val[j4]*l_wa_sum[cell_id]
                               + p_weight*delta[1]*delta_n[2]) / wa_sum_n;
                    val[j5] = (  val[j5]*l_wa_sum[cell_id]
                               + p_weight*delta[0]*delta_n[2]) / wa_sum_n;

                    /* update mean value */

                    for (cs_lnum_t l = 0; l < 3; l++)
                      mean_val[cell_id*3 + l] += r[l];

                  }

                  else { /* simple variance */

                    /* new weight for the cell: weight attached to
                       current particle (=dt*weight) plus old weight */

                    const cs_lnum_t dim = mt->dim;

                    for (cs_lnum_t l = 0; l < dim; l++) {

                      double delta = pval[l] - mean_val[cell_id*dim+l];
                      double r = delta * (p_weight / (fmax(wa_sum_n, 1e-100)));
                      double m_n = mean_val[cell_id*dim+l] + r;

                      val[cell_id*dim+l]
                        = (  val[cell_id*dim+l]*l_wa_sum[cell_id]
                           + (p_weight*delta*(pval[l]-m_n))) / wa_sum_n;

                      /* update mean value */

                      mean_val[cell_id*dim+l] += r;

                    }

                  }

                }

                else if (mt->m_type == CS_LAGR_MOMENT_MEAN) {

                  const cs_lnum_t dim = mt->dim;

                  for (cs_lnum_t l = 0; l < dim; l++)
                    val[cell_id*dim+l] +=   (pval[l] - val[cell_id*dim+l])
                                          * p_weight / (fmax(wa_sum_n, 1e-100));

                } /* End of test if moment is a variance or a mean */

                /* update local weight associated to current moment and class */

                l_wa_sum[cell_id] += p_weight;

              } /* End of test if particle is in a cell
                   and if particle class corresponds to moment class */

            } /* end of loop on particles */

            if (mt->p_data_func != NULL)
              BFT_FREE(pval);

            mt->nt_cur = ts->nt_cur;
            if (mt->m_type == CS_LAGR_MOMENT_VARIANCE)
              mt_mean->nt_cur = ts->nt_cur;
          }

          /* Case where data is mesh-based */
          /*-------------------------------*/

          else
            _cs_lagr_stat_update_mesh_moment(mt,
                                             mwa,
                                             m_weight,
                                             ts->nt_cur);

        } /* end of test if moment is for the current class */

      } /* End of loop on moment types */

    } /* End of loop on moments */

    /* At end of loop on moments inside a class, update
       global class weight array */

    if (l_wa_sum != NULL) {
      for (cs_lnum_t i = 0; i < n_w_elts; i++)
        g_wa_sum[i] = l_wa_sum[i];
      BFT_FREE(l_wa_sum);
    }
    else if (m_weight != NULL) {
      _update_wa_m(mwa, m_weight);
      if (m_weight != m_w0)
        BFT_FREE(m_weight);
    }

  } /* End of loop on active weigh accumulators */
}


/*----------------------------------------------------------------------------
 * Free all moments
 *----------------------------------------------------------------------------*/

static void
_free_all_moments(void)
{
  int i;

  for (i = 0; i < _n_lagr_stats; i++) {

    cs_lagr_moment_t *mt = _lagr_stats + i;
    BFT_FREE(mt->name);
    BFT_FREE(mt->val);
    //TODO
    /* field associated with moments are destroyed in cs_field_destroy_all */
  }

  BFT_FREE(_lagr_stats);

  _n_lagr_stats = 0;
  _n_lagr_stats_max = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all moment weight and time accumulators
 */
/*----------------------------------------------------------------------------*/

static void
_free_all_wa(void)
{
  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    BFT_FREE(mwa->val);
  }

  BFT_FREE(_lagr_stats_wa);

  _n_lagr_stats_wa = 0;
  _n_lagr_stats_wa_max = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particle statistic.
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name             statistics base name
 * \param[in]  location_id      id of associated mesh location
 * \param[in]  stat_type        predefined statistics type, or -1
 * \param[in]  m_type           moment type
 * \param[in]  class_id         particle class id, or 0 for all
 * \param[in]  dim              dimension associated with element data
 * \param[in]  component_id     attribute component id, or < 0 for all
 * \param[in]  p_data_func      pointer to particle function to compute
 *                              statistics (if stat_type < 0)
 * \param[in]  m_data_func      pointer to mesh location function to compute
 *                              statistics (if stat_type < 0)
 * \param[in]  data_input       associated input
 * \param[in]  w_p_data_func    pointer to particle function to compute weight
 *                              (if NULL, statistic weight assumed)
 * \param[in]  w_m_data_func    pointer to mesh location function to compute
 *                              weight (if stat_type < 0)
 * \param[in]  w_data_input     associated input for w_data_func
 * \param[in]  nt_start         starting time step (or -1 to use t_start,
 *                              0 to use idstnt)
 * \param[in]  t_start          starting time
 * \param[in]  restart_mode     behavior in case of restart (reset,
 *                              automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

static int
_stat_define(const char                *name,
             int                        location_id,
             int                        stat_type,
             cs_lagr_stat_moment_t      m_type,
             int                        class_id,
             int                        dim,
             int                        component_id,
             cs_lagr_moment_p_data_t   *p_data_func,
             cs_lagr_moment_m_data_t   *m_data_func,
             void                      *data_input,
             cs_lagr_moment_p_data_t   *w_p_data_func,
             cs_lagr_moment_m_data_t   *w_m_data_func,
             void                      *w_data_input,
             int                        nt_start,
             double                     t_start,
             cs_lagr_stat_restart_t     restart_mode)
{
  char name_buf[64];
  const char *_name = name;

  const int attr_id = cs_lagr_stat_type_to_attr_id(stat_type);

  if (attr_id > 0) {
    _stat_name(stat_type,
               component_id,
               class_id,
               m_type,
               name_buf);
    _name = name_buf;
  }

  int wa_location_id = location_id;

  cs_lagr_moment_t *mt = NULL;

  int moment_dim = (dim == 3 && m_type == CS_LAGR_MOMENT_VARIANCE) ? 6 : dim;
  int moment_id = -1;
  int prev_id = -1, prev_wa_id = -1;
  int _nt_start = nt_start;
  double _t_start = t_start;

  const cs_time_step_t  *ts = cs_glob_time_step;

  /* Optimization for constant mesh-based weights */

  if (w_m_data_func == _unit_value_m_cells && ts->is_local == 0)
    wa_location_id = 0;

  /* If this is the first moment to be defined, ensure
     restart data is read if available */

  if (_restart_info_checked == false)
    _restart_info_read();

  /* Find matching restart data if required, and adjust accumulator
     info if moment should be restarted, or if available data does
     not match and we do not require exact mode. */

  if (_restart_info != NULL) {

    prev_id = _check_restart(_name,
                             ts,
                             _restart_info,
                             location_id,
                             wa_location_id,
                             moment_dim,
                             m_type,
                             stat_type,
                             class_id,
                             &_nt_start,
                             &_t_start,
                             restart_mode);

    if (prev_id > -1)
      prev_wa_id = _restart_info->wa_id[prev_id];

  }

  if (_nt_start < 0 && _t_start < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Lagrangian statistics definition for \"%s\" is inconsistent:\n"
                " either starting time step or physical time must be >= 0."),
              name);

  /* Find or define matching weight accumulator info */

  const int wa_id = _find_or_add_wa(w_p_data_func,
                                    w_m_data_func,
                                    w_data_input,
                                    class_id,
                                    wa_location_id,
                                    _nt_start,
                                    _t_start,
                                    prev_wa_id);

  /* Check for possible previous definition */

  cs_field_t *f = cs_field_by_name_try(name);

  if (f != NULL) {

    for (int i = 0; i < _n_lagr_stats; i++) {

      mt = _lagr_stats + i;

      if (mt->f_id == f->id) {

        moment_id = i;
        return moment_id;

      }

    }

  }

  /* if not returned yet, build new moment */

  moment_id = _find_or_add_moment(location_id,
                                  component_id,
                                  class_id,
                                  stat_type,
                                  dim,
                                  p_data_func,
                                  m_data_func,
                                  data_input,
                                  m_type,
                                  wa_id,
                                  prev_id);

  mt = _lagr_stats + moment_id;
  BFT_FREE(mt->name); /* in case previously defined as sub-moment */

  /* matching field */

  f = _cs_lagr_moment_associate_field(_name, location_id, mt->dim);

  mt->f_id = f->id;

  /* Define sub moments */

  if (mt->m_type == CS_LAGR_MOMENT_VARIANCE) {

    prev_id = -1;
    if (_restart_info != NULL) {
      char m_name[128];
      snprintf(m_name, 127, "mean%s", _name+3);
      prev_id = _check_restart(m_name,
                               ts,
                               _restart_info,
                               location_id,
                               wa_location_id,
                               dim,
                               CS_LAGR_MOMENT_MEAN,
                               stat_type,
                               class_id,
                               &_nt_start,
                               &_t_start,
                               restart_mode);
    }

    int l_id = _find_or_add_moment(location_id,
                                   component_id,
                                   class_id,
                                   stat_type,
                                   dim,
                                   p_data_func,
                                   m_data_func,
                                   data_input,
                                   CS_LAGR_MOMENT_MEAN,
                                   wa_id,
                                   prev_id);

    mt->l_id = l_id;
    mt = _lagr_stats + l_id;

    if (mt->f_id < 0) {
      char s[64];
      snprintf(s, 64, "<auto_mean_particle_stat_%d>",
               l_id);
      s[63] = '\0';
      BFT_MALLOC(mt->name, strlen(s)+1, char);
      strcpy(mt->name, s);
    }

  }

  return moment_id;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particle statistic.
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_type      predefined statistics type, or -1
 * \param[in]  m_type         moment type
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  dim            dimension associated with element data
 * \param[in]  component_id   attribute component id, or < 0 for all
 * \param[in]  data_func      pointer to function to compute statistics
 *                            (if stat_type < 0)
 * \param[in]  data_input     associated input
 * \param[in]  w_data_func    pointer to function to compute weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  w_data_input   associated input for w_data_func
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_define(const char                *name,
                    int                        location_id,
                    int                        stat_type,
                    cs_lagr_stat_moment_t      m_type,
                    int                        class_id,
                    int                        dim,
                    int                        component_id,
                    cs_lagr_moment_p_data_t   *data_func,
                    void                      *data_input,
                    cs_lagr_moment_p_data_t   *w_data_func,
                    void                      *w_data_input,
                    int                        nt_start,
                    double                     t_start,
                    cs_lagr_stat_restart_t     restart_mode)
{
  return  _stat_define(name,
                       location_id,
                       stat_type,
                       m_type,
                       class_id,
                       dim,
                       component_id,
                       data_func,
                       NULL,
                       data_input,
                       w_data_func,
                       NULL,
                       w_data_input,
                       nt_start,
                       t_start,
                       restart_mode);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particle weight type statistic.
 *
 * Weights are automatically associated to general statitistics, but defining
 * them explicitely allows activation of standard logging and postprocessing
 * for those weights, as well as defining specific weights.
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  w_data_func    pointer to function to compute particle weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  w_data_input   associated input for w_data_func
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_accumulator_define(const char                *name,
                                int                        location_id,
                                int                        class_id,
                                cs_lagr_moment_p_data_t   *w_data_func,
                                void                      *w_data_input,
                                int                        nt_start,
                                double                     t_start,
                                cs_lagr_stat_restart_t     restart_mode)
{
  int wa_location_id = location_id;

  int prev_wa_id = -1;
  int _nt_start = nt_start;
  double _t_start = t_start;

  const cs_time_step_t  *ts = cs_glob_time_step;

  /* If this is the first moment to be defined, ensure
     restart data is read if available */

  if (_restart_info_checked == false)
    _restart_info_read();

  /* Find matching restart data if required, and adjust accumulator
     info if moment should be restarted, or if available data does
     not match and we do not require exact mode. */

  if (_restart_info != NULL) {

    /* TODO adjust this */

    prev_wa_id = _check_restart(name,
                                ts,
                                _restart_info,
                                location_id,
                                wa_location_id,
                                1,
                                0,
                                -1,
                                class_id,
                                &_nt_start,
                                &_t_start,
                                restart_mode);

  }

  if (_nt_start < 0 && _t_start < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Lagrangian statistics definition for \"%s\" is inconsistent:\n"
                " either starting time step or physical time must be >= 0."),
              name);

  /* Find or define matching weight accumulator info */

  const int wa_id = _find_or_add_wa(w_data_func,
                                    NULL,
                                    w_data_input,
                                    class_id,
                                    wa_location_id,
                                    _nt_start,
                                    _t_start,
                                    prev_wa_id);

  /* matching field */

  _cs_lagr_moment_associate_field(name, location_id, 1);

  return wa_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time moment associated to particle statistics.
 *
 * This is similar to general time moments (see \ref cs_time_moment.c),
 * with restart, logging, and unsteady reinitialization behavior
 * similar to other particle statistics.
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_type      predefined statistics type, or -1
 * \param[in]  m_type         moment type
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  dim            dimension associated with element data
 * \param[in]  component_id   attribute component id, or < 0 for all
 * \param[in]  data_func      pointer to function to compute statistics
 *                            (if stat_type < 0)
 * \param[in]  data_input     associated input
 * \param[in]  w_data_func    pointer to function to compute weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  w_data_input   associated input for w_data_func
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_time_moment_define(const char                *name,
                                int                        location_id,
                                int                        stat_type,
                                cs_lagr_stat_moment_t      m_type,
                                int                        class_id,
                                int                        dim,
                                int                        component_id,
                                cs_lagr_moment_m_data_t   *data_func,
                                void                      *data_input,
                                cs_lagr_moment_m_data_t   *w_data_func,
                                void                      *w_data_input,
                                int                        nt_start,
                                double                     t_start,
                                cs_lagr_stat_restart_t     restart_mode)
{
  return _stat_define(name,
                      location_id,
                      stat_type,
                      m_type,
                      class_id,
                      dim,
                      component_id,
                      NULL,
                      data_func,
                      data_input,
                      NULL,
                      w_data_func,
                      w_data_input,
                      nt_start,
                      t_start,
                      restart_mode);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate Lagrangian statistics for a given statistics type.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  stat_type   particle statistics type
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_activate(int  stat_type)
{
  const int n_stat_types = _n_stat_types();

  const int attr_id = cs_lagr_stat_type_to_attr_id(stat_type);

  if (attr_id > -1)
    cs_lagr_particle_attr_in_range(attr_id);
  else if (stat_type < 0)
    return;

  /* Setup flag if not already done */

  if (_vol_stat_activate == NULL) {
    BFT_MALLOC(_vol_stat_activate, n_stat_types, bool);
    for (int i = 0; i < n_stat_types; i++)
      _vol_stat_activate[i] = false;
  }

  _vol_stat_activate[stat_type] = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deactivate Lagrangian statistics for a given statistics type.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  stat_type   particle statistics type
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_deactivate(int  stat_type)
{
  const int attr_id = cs_lagr_stat_type_to_attr_id(stat_type);

  if (attr_id > -1)
    cs_lagr_particle_attr_in_range(attr_id);
  else if (stat_type < 0 || stat_type >= _n_stat_types())
    return;

  if (_vol_stat_activate != NULL)
    _vol_stat_activate[stat_type] = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate Lagrangian statistics for a given particle attribute.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  attr_id   particle attribute id
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_activate_attr(int  attr_id)
{
  cs_lagr_stat_activate(cs_lagr_stat_type_from_attr_id(attr_id));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deactivate Lagrangian statistics for a given particle attribute.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  attr_id   particle attribute id
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_deactivate_attr(int  attr_id)
{
  cs_lagr_stat_deactivate(cs_lagr_stat_type_from_attr_id(attr_id));
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Return statistics type associated with a given particle
 *        attribute id.
 *
 * \param[in]  attr_id  particle  attribute id
 *
 * \return  associated  particle statistics type id
 */
/*---------------------------------------------------------------------------*/

int
cs_lagr_stat_type_from_attr_id(int attr_id)
{
  cs_lagr_particle_attr_in_range(attr_id);

  return (attr_id + CS_LAGR_STAT_PARTICLE_ATTR);
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Return attribute id associated with a given statistics type.
 *
 * \param[in]   stat_type     particle statistics type
 *
 * \return attribute id, or -1 if not applicable
 */
/*---------------------------------------------------------------------------*/

int
cs_lagr_stat_type_to_attr_id(int  stat_type)
{
  int attr_id = -1;

  if (stat_type >= CS_LAGR_STAT_PARTICLE_ATTR)
    attr_id = stat_type - CS_LAGR_STAT_PARTICLE_ATTR;

  return attr_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map time step values array for Lagrangian statistics.
 *
 * If this function is not called, the field referenced by field pointer
 * CS_F_(dt) will be used instead.
 *
 * \param[in]   dt   pointer to time step values array
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_map_cell_dt(const cs_real_t  *dt)
{
  _p_dt = dt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lagrangian statistics initialization.
 *
 * Statistics activated or deactivated by previous calls to
 * \ref cs_lagr_stat_activate, \ref cs_lagr_stat_deactivate,
 * \ref cs_lagr_stat_activate_attr, and \ref cs_lagr_stat_deactivate_attr
 * will be initialized here.
 *
 * Restart info will be used after to fill in the moments structure
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_initialize(void)
{
  const cs_lagr_attribute_map_t *p_am = cs_lagr_particle_get_attr_map();
  const cs_lagr_stat_options_t *stat_options = cs_glob_lagr_stat_options;
  cs_lagr_stat_restart_t restart_mode = (stat_options->isuist) ?
    CS_LAGR_MOMENT_RESTART_AUTO : CS_LAGR_MOMENT_RESTART_RESET;

  /* Automatic initializations based on physical options */

  _init_vars_attribute();

  /* init moments */
  char name[64];

  if (_vol_stat_activate != NULL) {

    for (int class = 0;
         class < cs_glob_lagr_model->n_stat_classes + 1;
         class++) {

      for (int stat_type = 0; stat_type < _n_stat_types(); stat_type++) {

        if (_vol_stat_activate[stat_type] == false)
          continue;

        if (stat_type == CS_LAGR_STAT_CUMULATIVE_WEIGHT) {
          _statistical_weight_name(class, name);
          cs_lagr_stat_accumulator_define(name,
                                          CS_MESH_LOCATION_CELLS,
                                          class,
                                          NULL,
                                          NULL,
                                          0,
                                          -1,
                                          restart_mode);
          continue;
        }

        int dim = 1;

        int attr_id = cs_lagr_stat_type_to_attr_id(stat_type);

        if (   attr_id == CS_LAGR_VELOCITY
            || attr_id == CS_LAGR_COORDS
            || attr_id == CS_LAGR_VELOCITY_SEEN
            || attr_id == CS_LAGR_TURB_STATE_1
            || attr_id == CS_LAGR_PRED_VELOCITY
            || attr_id == CS_LAGR_PRED_VELOCITY_SEEN)
          dim = 3;

        for (cs_lagr_stat_moment_t m_type = CS_LAGR_MOMENT_MEAN;
             m_type <= CS_LAGR_MOMENT_VARIANCE;
             m_type++) {

          _stat_name(stat_type, -1, class, m_type, name);

          if (stat_type == CS_LAGR_STAT_VOLUME_FRACTION) {

            cs_lagr_stat_time_moment_define
              (name,
               CS_MESH_LOCATION_CELLS,
               stat_type,
               m_type,
               class,
               dim,
               -1,                   /* component_id, */
               _vol_fraction,        /* data_func */
               NULL,                 /* data_input */
               _unit_value_m_cells,  /* w_data_func */
               NULL,                 /* w_data_input */
               0,
               -1,
               restart_mode);

          }

          else if (attr_id > -1) {

            int n_comp = p_am->count[0][attr_id];
            if (n_comp == dim)
              cs_lagr_stat_define(name,
                                  CS_MESH_LOCATION_CELLS,
                                  stat_type,
                                  m_type,
                                  class,
                                  dim,
                                  -1,      /* component_id, */
                                  NULL,    /* data_func */
                                  NULL,    /* data_input */
                                  NULL,    /* w_data_func */
                                  NULL,    /* w_data_input */
                                  0,
                                  -1,
                                  restart_mode);

            else {
              for (int c_id = 0; c_id < n_comp; c_id++)
                cs_lagr_stat_define(name,
                                    CS_MESH_LOCATION_CELLS,
                                    stat_type,
                                    m_type,
                                    class,
                                    1,       /* dim */
                                    c_id,
                                    NULL,    /* data_func */
                                    NULL,    /* data_input */
                                    NULL,    /* w_data_func */
                                    NULL,    /* w_data_input */
                                    0,
                                    -1,
                                    restart_mode);
            }

          }

        }

      }

    }

  }

  /* Activation status not needed after this stage */

  BFT_FREE(_vol_stat_activate);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update particle statistics for a given time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_update(void)
{
  if (_restart_info != NULL) {
    /* ensure correct particle attribute is associated with statistics number */
    if (cs_glob_lagr_stat_options->isuist == 1)
      _cs_lagr_moment_restart_read();
    _restart_info_free();
  }

  /* if unsteady statistics, reset stats, wa, and durations */
  if (   cs_glob_lagr_time_scheme->isttio == 0
      || (   cs_glob_lagr_time_scheme->isttio == 1
          && cs_glob_time_step->nt_cur <= cs_glob_lagr_stat_options->nstist))
    _cs_lagr_moment_reset_unsteady_stats();

  _cs_lagr_stat_update_all();

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all moments management metadata.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_finalize(void)
{
  _free_all_moments();
  _free_all_wa();

 _restart_info_checked = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log moment definition information for a given iteration.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_log_iteration(void)
{
  if (_n_lagr_stats_wa < 1)
    return;

  int n_active_wa = 0;

  const cs_time_step_t  *ts = cs_glob_time_step;

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (mwa->nt_start > 0 && mwa->nt_start <= ts->nt_cur)
      n_active_wa += 1;
  }

  if (n_active_wa < 1)
    return;

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  "  ** Particle moment accumulated weights\n"
                  "     -----------------------------------\n"));

  /* Info for accumulators on non-global locations */

  char tmp_s[5][64] =  {"", "", "", "", ""};

  /* Header */

  cs_log_strpad(tmp_s[0], _("id"), 4, 64);
  cs_log_strpad(tmp_s[1], _("n it."), 8, 64);
  cs_log_strpadl(tmp_s[2], _("minimum"), 14, 64);
  cs_log_strpadl(tmp_s[3], _("maximum"), 14, 64);
  cs_log_strpadl(tmp_s[4], _("set mean"), 14, 64);

  cs_log_printf(CS_LOG_DEFAULT, "\n");

  cs_log_printf(CS_LOG_DEFAULT, "   %s %s %s %s %s\n",
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);

  for (int j = 0; j < 5; j++)
    memset(tmp_s[j], '-', 64);

  tmp_s[0][4] = '\0';
  tmp_s[1][8] = '\0';
  tmp_s[2][14] = '\0';
  tmp_s[3][14] = '\0';
  tmp_s[4][14] = '\0';

  cs_log_printf(CS_LOG_DEFAULT, "   %s %s %s %s %s\n",
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);

  cs_gnum_t *n_g_elts;
  double *vmin, *vmax, *vsum;

  BFT_MALLOC(n_g_elts, n_active_wa, cs_gnum_t);
  BFT_MALLOC(vmin, n_active_wa, double);
  BFT_MALLOC(vmax, n_active_wa, double);
  BFT_MALLOC(vsum, n_active_wa, double);

  n_active_wa = 0;

  /* Determine min, max, sum */

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (   mwa->nt_start > 0 && mwa->nt_start <= ts->nt_cur
        && mwa->location_id > 0) {
      const cs_lnum_t n_elts
        = cs_mesh_location_get_n_elts(mwa->location_id)[0];
      const cs_mesh_location_type_t loc_type
        = cs_mesh_location_get_type(mwa->location_id);
      if (   loc_type == CS_MESH_LOCATION_CELLS
          || loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
        n_g_elts[n_active_wa] = n_elts;
      else
        n_g_elts[n_active_wa] = 0;

      const cs_real_t *val = _mwa_val(mwa);
      cs_array_reduce_simple_stats_l(n_elts,
                                     1,
                                     NULL,
                                     val,
                                     vmin + n_active_wa,
                                     vmax + n_active_wa,
                                     vsum + n_active_wa);
      n_active_wa += 1;
    }
  }

  /* Group MPI operations if required */

  cs_parall_counter(n_g_elts, n_active_wa);
  cs_parall_min(n_active_wa, CS_DOUBLE, vmin);
  cs_parall_max(n_active_wa, CS_DOUBLE, vmax);
  cs_parall_sum(n_active_wa, CS_DOUBLE, vsum);

  /* Now log values */

  n_active_wa = 0;

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (   mwa->nt_start > 0 && mwa->nt_start <= ts->nt_cur
        && mwa->location_id > 0) {

      int nt_acc = ts->nt_cur - mwa->nt_start + 1;

      if (n_g_elts[n_active_wa] > 0) {
        double v_mean = vsum[n_active_wa] / n_g_elts[n_active_wa];
        snprintf(tmp_s[4], 63, " %14.5g", v_mean);
        tmp_s[4][63] = '\0';
      }
      else
        tmp_s[4][0] = '\0';

      cs_log_printf(CS_LOG_DEFAULT, "   %-4d %-8d %14.5g %14.5g%s\n",
                    i, nt_acc,
                    vmin[n_active_wa], vmax[n_active_wa], tmp_s[4]);

      n_active_wa += 1;

    }
  }

  BFT_FREE(vsum);
  BFT_FREE(vmax);
  BFT_FREE(vmin);
  BFT_FREE(n_g_elts);

  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Checkpoint moment data
 *
 * \param[in]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_restart_write(cs_restart_t  *restart)
{
  int *nt_start, *location_id, *dimension;
  int *m_type, *wa_id, *l_id, *class, *stat_type;
  cs_real_t *t_start;

  int n_active_wa = 0, n_active_moments = 0;
  int *active_wa_id = NULL, *active_moment_id = NULL;

  if (_n_lagr_stats < 1)
    return;

  const cs_time_step_t  *ts = cs_glob_time_step;

  /* General information */
  /* ------------------- */

  BFT_MALLOC(active_wa_id, _n_lagr_stats_wa, int);
  BFT_MALLOC(active_moment_id, _n_lagr_stats, int);

  /* Check for active moments */

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (mwa->nt_start > 0 && mwa->nt_start <= ts->nt_cur) {
      active_wa_id[i] = n_active_wa;
      n_active_wa += 1;
    }
    else
      active_wa_id[i] = -1;
  }

  for (int i = 0; i < _n_lagr_stats; i++) {
    cs_lagr_moment_t *mt = _lagr_stats + i;
    if (active_wa_id[mt->wa_id] > -1) {
      active_moment_id[i] = n_active_moments;
      n_active_moments += 1;
    }
    else
      active_moment_id[i] = -1;
  }

  if (n_active_moments < 1) {
    BFT_FREE(active_wa_id);
    BFT_FREE(active_moment_id);
    return;
  }

  /* Build global names array */

  size_t names_max_size = 32;
  int *names_idx;
  char *names;

  BFT_MALLOC(names_idx, n_active_moments + 1, int);
  BFT_MALLOC(names, names_max_size, char);

  names_idx[0] = 0;

  for (int i = 0; i < _n_lagr_stats; i++) {

    const int j = active_moment_id[i];
    if (j > -1) {

      cs_lagr_moment_t *mt = _lagr_stats + i;
      const char *name = NULL;
      if (mt->f_id > -1) {
        const cs_field_t *f = cs_field_by_id(mt->f_id);
        name = f->name;
      }
      else
        name = mt->name;
      const size_t l = strlen(name) + 1;
      if (names_idx[i] + l > names_max_size) {
        while (names_idx[i] + l > names_max_size)
          names_max_size *= 2;
        BFT_REALLOC(names, names_max_size, char);
      }
      strcpy(names + names_idx[i], name);
      names[names_idx[i] + l - 1] = '\0';
      names_idx[i+1] = names_idx[i] + l;

    }

  }

  int sizes[3] = {n_active_wa,
                  n_active_moments,
                  names_idx[n_active_moments]};

  cs_restart_write_section(restart,
                           "lagr_stats:sizes",
                           CS_MESH_LOCATION_NONE,
                           3,
                           CS_TYPE_cs_int_t,
                           sizes);

  cs_restart_write_section(restart,
                           "lagr_stats:names",
                           CS_MESH_LOCATION_NONE,
                           names_idx[n_active_moments],
                           CS_TYPE_char,
                           names);

  BFT_FREE(names_idx);
  BFT_FREE(names);

  /* Information on moment weight accumulators */

  BFT_MALLOC(location_id, n_active_wa, int);
  BFT_MALLOC(nt_start, n_active_wa, int);
  BFT_MALLOC(t_start, n_active_wa, cs_real_t);

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    int j = active_wa_id[i];
    if (j > -1) {
      cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
      location_id[j] = mwa->location_id;
      nt_start[j] = mwa->nt_start;
      t_start[j] = mwa->t_start;
    }
  }

  cs_restart_write_section(restart,
                           "lagr_stats:wa:location_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_wa,
                           CS_TYPE_cs_int_t,
                           location_id);

  cs_restart_write_section(restart,
                           "lagr_stats:wa:nt_start",
                           CS_MESH_LOCATION_NONE,
                           n_active_wa,
                           CS_TYPE_cs_int_t,
                           nt_start);

  cs_restart_write_section(restart,
                           "lagr_stats:wa:t_start",
                           CS_MESH_LOCATION_NONE,
                           n_active_wa,
                           CS_TYPE_cs_real_t,
                           t_start);

  BFT_FREE(t_start);
  BFT_FREE(nt_start);
  BFT_FREE(location_id);

  /* To be decided for wa save, already in lagout*/
  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    int j = active_wa_id[i];
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (j > -1 && mwa->location_id > CS_MESH_LOCATION_NONE) {
      char s[64];
      snprintf(s, 64, "lagr_stats:wa:%02d:val", i);
      cs_restart_write_section(restart,
                               s,
                               mwa->location_id,
                               1,
                               CS_TYPE_cs_real_t,
                               _mwa_val(mwa));
    }
  }

  /* Information on moments proper */

  BFT_MALLOC(m_type, n_active_moments, int);
  BFT_MALLOC(class, n_active_moments, int);
  BFT_MALLOC(location_id, n_active_moments, int);
  BFT_MALLOC(dimension, n_active_moments, int);
  BFT_MALLOC(wa_id, n_active_moments, int);
  BFT_MALLOC(l_id, n_active_moments, int);
  BFT_MALLOC(stat_type, n_active_moments, int);

  for (int i = 0; i < _n_lagr_stats; i++) {
    int j = active_moment_id[i];
    if (j > -1) {
      cs_lagr_moment_t *mt = _lagr_stats + i;
      m_type[j] = mt->m_type;
      class[j] = mt->class;
      location_id[j] = mt->location_id;
      dimension[j] = mt->dim;
      wa_id[j] = active_wa_id[mt->wa_id];
      stat_type[j] = mt->stat_type;
      if (mt->l_id > -1)
        l_id[j] = active_moment_id[mt->l_id];
      else
        l_id[j] = -1;
    }
  }

  cs_restart_write_section(restart,
                           "lagr_stats:type",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           m_type);

  cs_restart_write_section(restart,
                           "lagr_stats:class",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           class);

  cs_restart_write_section(restart,
                           "lagr_stats:location_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           location_id);

  cs_restart_write_section(restart,
                           "lagr_stats:dimension",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           dimension);

  cs_restart_write_section(restart,
                           "lagr_stats:wa_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           wa_id);

  cs_restart_write_section(restart,
                           "lagr_stats:lower_order_id",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           l_id);

  cs_restart_write_section(restart,
                           "lagr_stats:stat_type",
                           CS_MESH_LOCATION_NONE,
                           n_active_moments,
                           CS_TYPE_cs_int_t,
                           stat_type);

  BFT_FREE(l_id);
  BFT_FREE(wa_id);
  BFT_FREE(dimension);
  BFT_FREE(location_id);
  BFT_FREE(m_type);
  BFT_FREE(class);
  BFT_FREE(stat_type);

  /* Write of moments value */
  for (int i = 0; i < _n_lagr_stats; i++) {

    int j = active_moment_id[i];

    if (j > -1) {

      cs_lagr_moment_t *mt = _lagr_stats + i;
      const cs_field_t *f = cs_field_by_id(mt->f_id);
      cs_restart_write_section(restart,
                               f->name,
                               f->location_id,
                               f->dim,
                               CS_TYPE_cs_real_t,
                               f->val);

    }

  }

  BFT_FREE(active_moment_id);
  BFT_FREE(active_wa_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return field associated with a given Lagrangian statistic,
 *        given a statistics type (i.e. variable), moment order, and
 *        given a statistics type (i.e. variable), moment order,
 *        statistical class, and component id.
 *
 * \param[in]  stat_type     statistics type
 * \param[in]  m_type        moment type (mean or variance)
 * \param[in]  class_id      particle statistical class
 * \param[in]  component_id  component id, or -1 for all
 *
 * \returns pointer to the field associated to the corresponding moment
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_lagr_stat_get_moment(int                    stat_type,
                        cs_lagr_stat_moment_t  m_type,
                        int                    class_id,
                        int                    component_id)
{
  assert(class_id >= 0);
  assert(m_type == CS_LAGR_MOMENT_MEAN || m_type == CS_LAGR_MOMENT_VARIANCE);

  for (int m_id = 0; m_id < _n_lagr_stats; m_id++) {

    cs_lagr_moment_t *mt = _lagr_stats + m_id;

    if (   mt->m_type       == m_type
        && mt->stat_type    == stat_type
        && mt->class        == class_id
        && mt->component_id == component_id)
      return cs_field_by_id(mt->f_id);

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return statistical weight
 *
 * \param[in]  class_id    particle statistical class
 *
 * \returns pointer to the field associated to the corresponding weight
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_lagr_stat_get_stat_weight(int  class_id)
{
  assert(class_id >= 0);

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (   mwa->f_id > -1
        && mwa->class == class_id)
      return cs_field_by_id(mwa->f_id);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return global volume statistics age
 *
 * \returns age of volume statistics, or -1 if statistics not active yet
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_stat_get_age(void)
{
  cs_real_t retval = -1.;

  for (int i = 0; i < _n_lagr_stats_wa; i++) {
    cs_lagr_moment_wa_t *mwa = _lagr_stats_wa + i;
    if (mwa->f_id > -1 && mwa->class == 0) {
      const cs_time_step_t  *ts = cs_glob_time_step;
      if (mwa->nt_start >= ts->nt_cur)
        retval = ts->t_cur - mwa->t_start;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return statistics age for a given moment
 *
 * \param[in]  f  field associated with given statistic
 *
 * \returns age of give statistic, or -1 if not active yet
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_stat_get_moment_age(cs_field_t  *f)
{
  cs_real_t retval = -1.;

  for (int m_id = 0; m_id < _n_lagr_stats; m_id++) {
    cs_lagr_moment_t *mt = _lagr_stats + m_id;
    if (mt->f_id == f->id) {
      cs_lagr_moment_wa_t  *mwa = _lagr_stats_wa + mt->wa_id;
      const cs_time_step_t  *ts = cs_glob_time_step;
      if (mwa->nt_start >= ts->nt_cur)
        retval = ts->t_cur - mwa->t_start;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
