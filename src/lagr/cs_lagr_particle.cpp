/*============================================================================
 * Lagrangian particle model
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

#include "bft/bft_printf.h"
#include "bft/bft_error.h"

#include "fvm/fvm_periodicity.h"

#include "base/cs_base.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_random.h"
#include "base/cs_timer_stats.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_post.h"
#include "lagr/cs_lagr_clogging.h"
#include "lagr/cs_lagr_roughness.h"
#include "lagr/cs_lagr_dlvo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_particle.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_particle.cpp
        Particle structure.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define  N_GEOL 13
#define  CS_LAGR_MIN_COMM_BUF_SIZE  8

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/* keys to sort attributes by type.

   rvar real values at current and previous time steps
   ivar integer values at current and previous time steps
   rprp real properties at current and previous time steps
   iprp integer properties at current and previous time steps
   rkid values are for rank ids, useful and valid only for previous
   time steps */

typedef enum {
  CS_LAGR_P_RVAR_TS = 1, /* CS_LAGR_P_RVAR with possible source terms */
  CS_LAGR_P_RVAR,
  CS_LAGR_P_IVAR,
  CS_LAGR_P_RPRP,
  CS_LAGR_P_IPRP,
  CS_LAGR_P_RKID,
} _array_map_id_t;

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/* Private tracking data associated to each particle */
/* --------------------------------------------------*/

/* This structure is a copy of the one defined in cs_lagr_tracking.c,
   which is currently mapped to the beginning of each
 * particle's data, and contains values which are used during the
 * tracking algorithm only.
 * It could be separated in the future, but this would require
 * keeping track of a particle's local id in most functions. */

typedef struct {

  cs_real_t  start_coords[3];       /* starting coordinates for
                                       next displacement */

  cs_lnum_t  last_face_num;         /* last face number encountered */

  int        state;                 /* current state (actually an enum) */

  int        integ_tracked_loc;     /* useful when using cell_wise_integ
                                       0 : determnistic virtual partner tracked
                                       1 : stochastic particle tracked
                                       2 : particle tracking finished*/
} cs_lagr_tracking_info_t;

/* Particle data value */
/*---------------------*/

union cs_lagr_value_t {
  cs_lnum_t      l; /* v_lnum_t */
  cs_gnum_t      g; /* v_gnum_t */
  cs_real_t      f; /* v_real_t */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Enumerator names */

const char *cs_lagr_attribute_name[] = {
  "state_flag",
  "cell_id",
  "rank_id",
  "rebound_id",
  "random_value",
  "stat_weight",
  "residence_time",
  "mass",
  "diameter",
  "taup_aux",
  "coords",
  "velocity",
  "velocity_seen",
  "velocity_seen_velocity_cov",
  "tr_truncate",
  "tr_reposition",

  /* Arrays for 2nd order scheme */
  "turb_state_1",
  "brown_state_1",
  "pred_velocity",
  "pred_velocity_seen",
  "v_gauss",
  "br_gauss",

  /* Non-spherical particles submodel additoinal parameters */
  "shape",
  "orientation",
  "quaternion",
  "radii",
  "angular_vel",
  "euler",
  "shape_param",

  /* Deposition submodel additional parameters */
  "yplus",
  "interf",
  "neighbor_face_id",
  "marko_value",
  "fouling_index",

  "n_large_asperities",
  "n_small_asperities",
  "adhesion_force",
  "adhesion_torque",
  "displacement_norm",

  "height",
  "cluster_nb_part",
  "depo_time",
  "consol_height",

  "temperature",
  "fluid_temperature",
  "cp",

  "water_mass",
  "coal_mass",
  "coke_mass",
  "shrinking_diameter",
  "initial_diameter",
  "coal_id",
  "coal_density",

  "emissivity",

  "stat_class",

  "agglo_class_id",
  "agglo_fractal_dim",

  "remaining_integration_time",
  "user",
  "<none>"};

/* Global particle attributes map */

static cs_lagr_attribute_map_t  *_p_attr_map = nullptr;

/* Particle set reallocation parameters */

static  double              _reallocation_factor = 2.0;
static  unsigned long long  _n_g_max_particles = ULLONG_MAX;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointer to the main particle set */

cs_lagr_particle_set_t *cs_glob_lagr_particle_set = nullptr;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * Compute new extents to ensure alignment of data
 *
 * returns:
 *   padded extents ensuring alignement
 *----------------------------------------------------------------------------*/

static size_t
_align_extents(size_t  size)
{
  size_t retval = size;

  size_t align_size = sizeof(union cs_lagr_value_t);

  size_t r = size % align_size;
  if (r > 0)
    retval += (align_size - r);

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump a particle structure
 *
 * parameter
 *   p_set   <-- cs_lagr_particle_set_t structure to dump
 *   p_id <-- id of particle to dump
 *----------------------------------------------------------------------------*/

static void
_dump_particle(const cs_lagr_particle_set_t  *p_set,
               cs_lnum_t                      p_id)
{
  const cs_lagr_attribute_map_t *am = p_set->p_am;

  bft_printf("  particle: %lu\n", (unsigned long)p_id);

  for (int time_id = 0; time_id < p_set->p_am->n_time_vals; time_id++) {

    if (time_id == 0)
      bft_printf("    values at time n:\n");
    else
      bft_printf("    values at time: n-%d\n", time_id);

    for (int i_attr = 0;
         i_attr < CS_LAGR_N_ATTRIBUTES;
         i_attr++) {
      auto attr = static_cast<cs_lagr_attribute_t>(i_attr);
      if (am->count[time_id][attr] > 0) {
        const char *attr_name = cs_lagr_attribute_name[attr];
        switch (am->datatype[attr]) {
        case CS_LNUM_TYPE:
          {
            const cs_lnum_t *v
              = cs_lagr_particles_attr_n_get_const_ptr<cs_lnum_t>(p_set,
                                                                  p_id,
                                                                  time_id,
                                                                  attr);
            bft_printf("      %24s: %10ld\n", attr_name, (long)v[0]);
            for (int i = 1; i < am->count[time_id][attr]; i++)
              bft_printf("      %24s: %10ld\n", " ", (long)v[i]);
          }
          break;
        case CS_GNUM_TYPE:
          {
            const cs_gnum_t *v
              = cs_lagr_particles_attr_n_get_const_ptr<cs_gnum_t>(p_set,
                                                                  p_id,
                                                                  time_id,
                                                                  attr);
            bft_printf("      %24s: %10lu\n", attr_name, (unsigned long)v[0]);
            for (int i = 1; i < am->count[time_id][attr]; i++)
              bft_printf("      %24s: %10lu\n", " ", (unsigned long)v[i]);
        }
          break;
        case CS_REAL_TYPE:
          {
            const cs_real_t *v
              = cs_lagr_particles_attr_n_get_const_ptr<cs_real_t>(p_set,
                                                                  p_id,
                                                                  time_id,
                                                                  attr);
            bft_printf("      %24s: %10.3g\n", attr_name, v[0]);
            for (int i = 1; i < am->count[time_id][attr]; i++)
              bft_printf("      %24s: %10.3g\n", " ", v[i]);
          }
          break;
        default:
          break;
        }
      }
    }
  }
  bft_printf("\n");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle map based on defined options.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_attr_initialize(void)
{
  cs_lagr_model_t  *lagr_model = cs_glob_lagr_model;
  cs_lagr_time_scheme_t  *lagr_time_scheme = cs_glob_lagr_time_scheme;
  const cs_lagr_extra_module_t *extra_i = cs_glob_lagr_extra_module;
  const cs_lagr_extra_module_t *extra = extra_i;

  int  i;
  int n_phases = extra->n_phases;

  int loc_count = 0;

  int pepa_loc_add = 1000; /* should be above any j* pointer if used */

  /* attr_keys contains:
   * 0: type (0 if unused attr)
   * 1: attr_id (0 if unused attr)
   * 2: dim (0 if unused attr)
   * */
  cs_lnum_t attr_keys[CS_LAGR_N_ATTRIBUTES][3];

  /* Initialize global parameter relative to the Lagrangian module */

  /*  cs_glob_lagr_brownian->lamvbr = *lamvbr; */

  if (lagr_model->physical_model == CS_LAGR_PHYS_COAL)
    lagr_model->n_temperature_layers = cs_glob_lagr_const_dim->nlayer;
  else
    lagr_model->n_temperature_layers = 1;

  /* Set indexes */

  for (i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
    attr_keys[i][0] = CS_LAGR_P_RVAR; /* default */
    attr_keys[i][1] = 0;
    attr_keys[i][2] = 0;
  }

  /* Special case:
     cell id is also needed in for previous time step.
  */

  attr_keys[CS_LAGR_CELL_ID][0] = CS_LAGR_P_IVAR;
  attr_keys[CS_LAGR_CELL_ID][1] = ++loc_count;

  attr_keys[CS_LAGR_RANK_ID][0] = CS_LAGR_P_RKID;
  attr_keys[CS_LAGR_RANK_ID][1] = ++loc_count;

  /* Other attributes */

  attr_keys[CS_LAGR_P_FLAG][0] = CS_LAGR_P_IPRP;
  attr_keys[CS_LAGR_P_FLAG][1] = ++loc_count;

  attr_keys[CS_LAGR_REBOUND_ID][0] = CS_LAGR_P_IPRP;
  attr_keys[CS_LAGR_REBOUND_ID][1] = ++loc_count;

  attr_keys[CS_LAGR_RANDOM_VALUE][0] = CS_LAGR_P_RPRP;
  attr_keys[CS_LAGR_RANDOM_VALUE][1] = ++loc_count;

  attr_keys[CS_LAGR_STAT_WEIGHT][0] = CS_LAGR_P_RPRP;
  attr_keys[CS_LAGR_STAT_WEIGHT][1] = ++loc_count;

  attr_keys[CS_LAGR_RESIDENCE_TIME][0] = CS_LAGR_P_RPRP;
  attr_keys[CS_LAGR_RESIDENCE_TIME][1] = ++loc_count;

  if (lagr_model->clogging == 1)
    attr_keys[CS_LAGR_HEIGHT][1] = ++loc_count;

  attr_keys[CS_LAGR_MASS][0] = CS_LAGR_P_RVAR_TS;
  attr_keys[CS_LAGR_MASS][1] = ++loc_count;

  attr_keys[CS_LAGR_DIAMETER][0] = CS_LAGR_P_RVAR_TS;
  attr_keys[CS_LAGR_DIAMETER][1] = ++loc_count;

  /* Activate Euler angles for spheroids / generic GLM
   * in order to compute a local frame of reference */
  // FIXME: change to activate all the time ?
  if (lagr_model->shape != CS_LAGR_SHAPE_SPHERE_MODEL
      || lagr_model->transport_GLM_rotated ) {
    attr_keys[CS_LAGR_EULER][1] = ++loc_count;
    attr_keys[CS_LAGR_EULER][2] = 4;
  }

  /* Non-sphere model
   * TODO activate only required arrays */
  if (lagr_model->shape != CS_LAGR_SHAPE_SPHERE_MODEL) {
    attr_keys[CS_LAGR_SHAPE][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_SHAPE][1] = ++loc_count;

    attr_keys[CS_LAGR_ORIENTATION][1] = ++loc_count;
    attr_keys[CS_LAGR_ORIENTATION][2] = 3;

    attr_keys[CS_LAGR_QUATERNION][1] = ++loc_count;
    attr_keys[CS_LAGR_QUATERNION][2] = 4;

    attr_keys[CS_LAGR_RADII][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_RADII][1] = ++loc_count;
    attr_keys[CS_LAGR_RADII][2] = 3;

    attr_keys[CS_LAGR_ANGULAR_VEL][1] = ++loc_count;
    attr_keys[CS_LAGR_ANGULAR_VEL][2] = 3;

    attr_keys[CS_LAGR_SHAPE_PARAM][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_SHAPE_PARAM][1] = ++loc_count;
    attr_keys[CS_LAGR_SHAPE_PARAM][2] = 4;
  }

  attr_keys[CS_LAGR_COORDS][1] = ++loc_count;
  attr_keys[CS_LAGR_COORDS][2] = 3;

  attr_keys[CS_LAGR_VELOCITY][1] = ++loc_count;
  attr_keys[CS_LAGR_VELOCITY][2] = 3;

  attr_keys[CS_LAGR_VELOCITY_SEEN][1] = ++loc_count;
  attr_keys[CS_LAGR_VELOCITY_SEEN][2] = 3 * n_phases;

  attr_keys[CS_LAGR_TR_TRUNCATE][0] = CS_LAGR_P_RPRP;
  attr_keys[CS_LAGR_TR_TRUNCATE][1] = ++loc_count;

  attr_keys[CS_LAGR_TR_REPOSITION][0] = CS_LAGR_P_IPRP;
  attr_keys[CS_LAGR_TR_REPOSITION][1] = ++loc_count;

  attr_keys[CS_LAGR_VELOCITY_SEEN_VELOCITY_COV][1] = ++loc_count;
  attr_keys[CS_LAGR_VELOCITY_SEEN_VELOCITY_COV][2] = 9 * n_phases;

  if (lagr_time_scheme->t_order > 1) {
    attr_keys[CS_LAGR_TAUP_AUX][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_TAUP_AUX][1] = ++pepa_loc_add;

    attr_keys[CS_LAGR_TURB_STATE_1][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_TURB_STATE_1][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_TURB_STATE_1][2] = 3;

    attr_keys[CS_LAGR_BROWN_STATE_1][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_BROWN_STATE_1][1] = ++pepa_loc_add;

    attr_keys[CS_LAGR_PRED_VELOCITY][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_PRED_VELOCITY][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_PRED_VELOCITY][2] = 3;

    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][2] = 3;

    if (cs_glob_lagr_model->idistu == 1) {
      attr_keys[CS_LAGR_V_GAUSS][0] = CS_LAGR_P_RPRP;
      attr_keys[CS_LAGR_V_GAUSS][1] = ++pepa_loc_add;
      attr_keys[CS_LAGR_V_GAUSS][2] = 3 * (n_phases + 2);
    }

    if (cs_glob_lagr_brownian->lamvbr == 1) {
      attr_keys[CS_LAGR_BR_GAUSS][0] = CS_LAGR_P_RPRP;
      attr_keys[CS_LAGR_BR_GAUSS][1] = ++pepa_loc_add;
      attr_keys[CS_LAGR_BR_GAUSS][2] = 6;
    }
  }

  if (lagr_model->deposition == 1) {

    attr_keys[CS_LAGR_YPLUS][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_YPLUS][1] = ++loc_count;

    attr_keys[CS_LAGR_INTERF][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_INTERF][1] = ++loc_count;

    attr_keys[CS_LAGR_NEIGHBOR_FACE_ID][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_NEIGHBOR_FACE_ID][1] = ++loc_count;

    attr_keys[CS_LAGR_MARKO_VALUE][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_MARKO_VALUE][1] = ++loc_count;

  }

  attr_keys[CS_LAGR_FOULING_INDEX][0] = CS_LAGR_P_RPRP;
  attr_keys[CS_LAGR_FOULING_INDEX][1] = ++loc_count;

  if (lagr_model->resuspension == 1) {

    attr_keys[CS_LAGR_N_LARGE_ASPERITIES][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_N_LARGE_ASPERITIES][1] = ++loc_count;

    attr_keys[CS_LAGR_N_SMALL_ASPERITIES][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_N_SMALL_ASPERITIES][1] = ++loc_count;

    attr_keys[CS_LAGR_ADHESION_FORCE][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_ADHESION_FORCE][1] = ++loc_count;

    attr_keys[CS_LAGR_ADHESION_TORQUE][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_ADHESION_TORQUE][1] = ++loc_count;

    attr_keys[CS_LAGR_DISPLACEMENT_NORM][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_DISPLACEMENT_NORM][1] = ++loc_count;

  }

  if (lagr_model->clogging == 1) {

    attr_keys[CS_LAGR_CLUSTER_NB_PART][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_CLUSTER_NB_PART][1] = ++loc_count;

    attr_keys[CS_LAGR_DEPO_TIME][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_DEPO_TIME][1] = ++loc_count;

    attr_keys[CS_LAGR_CONSOL_HEIGHT][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_CONSOL_HEIGHT][1] = ++loc_count;

  }

  if (lagr_model->physical_model == CS_LAGR_PHYS_HEAT) {

    if (cs_glob_lagr_specific_physics->solve_temperature_seen == 1) {
      attr_keys[CS_LAGR_TEMPERATURE_SEEN][0] = CS_LAGR_P_RVAR_TS;
      attr_keys[CS_LAGR_TEMPERATURE_SEEN][1] = ++loc_count;
    }

    if (cs_glob_lagr_specific_physics->solve_temperature == 1) {

      attr_keys[CS_LAGR_CP][1] = ++loc_count;

      attr_keys[CS_LAGR_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
      attr_keys[CS_LAGR_TEMPERATURE][1] = ++loc_count;

      if (extra->radiative_model > 0)
        attr_keys[CS_LAGR_EMISSIVITY][1] = ++loc_count;

    }

  }

  if (lagr_model->physical_model == CS_LAGR_PHYS_COAL) {

    attr_keys[CS_LAGR_CP][1] = ++loc_count;

    attr_keys[CS_LAGR_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
    attr_keys[CS_LAGR_TEMPERATURE][1] = ++loc_count;
    attr_keys[CS_LAGR_TEMPERATURE][2]
      = lagr_model->n_temperature_layers;

    attr_keys[CS_LAGR_TEMPERATURE_SEEN][0] = CS_LAGR_P_RVAR_TS;
    attr_keys[CS_LAGR_TEMPERATURE_SEEN][1] = ++loc_count;

    attr_keys[CS_LAGR_WATER_MASS][1] = ++loc_count;

    attr_keys[CS_LAGR_COAL_MASS][1] = ++loc_count;
    attr_keys[CS_LAGR_COAL_MASS][2]
      = lagr_model->n_temperature_layers;

    attr_keys[CS_LAGR_COKE_MASS][1] = ++loc_count;
    attr_keys[CS_LAGR_COKE_MASS][2]
      = lagr_model->n_temperature_layers;

    attr_keys[CS_LAGR_SHRINKING_DIAMETER][1] = ++loc_count;

    attr_keys[CS_LAGR_INITIAL_DIAMETER][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_INITIAL_DIAMETER][1] = ++loc_count;

    attr_keys[CS_LAGR_COAL_DENSITY][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_COAL_DENSITY][1] = ++loc_count;
    attr_keys[CS_LAGR_COAL_DENSITY][2]
      = lagr_model->n_temperature_layers;

    attr_keys[CS_LAGR_COAL_ID][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_COAL_ID][1] = ++loc_count;

  }

  if (lagr_model->physical_model == CS_LAGR_PHYS_CTWR) {

      attr_keys[CS_LAGR_CP][0] = CS_LAGR_P_RVAR;
      attr_keys[CS_LAGR_CP][1] = ++loc_count;

      attr_keys[CS_LAGR_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
      attr_keys[CS_LAGR_TEMPERATURE][1] = ++loc_count;

      attr_keys[CS_LAGR_TEMPERATURE_SEEN][0] = CS_LAGR_P_RVAR_TS;
      attr_keys[CS_LAGR_TEMPERATURE_SEEN][1] = ++loc_count;

// FIXME: accounting for droplet radiation will be done later
/*      if (extra->radiative_model > 0)
        attr_keys[CS_LAGR_EMISSIVITY][1] = ++loc_count;*/

  }

  if (lagr_model->n_stat_classes > 0) {
    attr_keys[CS_LAGR_STAT_CLASS][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_STAT_CLASS][1] = ++loc_count;
  }

  if (lagr_model->agglomeration == 1 ||
      lagr_model->fragmentation == 1 ) {
    attr_keys[CS_LAGR_AGGLO_CLASS_ID][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_AGGLO_CLASS_ID][1] = ++loc_count;

    attr_keys[CS_LAGR_AGGLO_FRACTAL_DIM][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_AGGLO_FRACTAL_DIM][1] = ++loc_count;
  }

  if (lagr_time_scheme->cell_wise_integ == 1)
    attr_keys[CS_LAGR_REMAINING_INTEG_TIME][1] = ++loc_count;

  if (lagr_model->n_user_variables > 0) {
    attr_keys[CS_LAGR_USER][1] = ++loc_count;
    attr_keys[CS_LAGR_USER][2] = lagr_model->n_user_variables;
  }

  /* Default count of 1 */

  for (i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
    if (attr_keys[i][1] > 0 && attr_keys[i][2] == 0)
      attr_keys[i][2] = 1;
    else if (attr_keys[i][1] < 1)
      attr_keys[i][0] = 0;
  }

  /* Build mappings
     (in the future, they should be created first, then marked,
     then built) */

  _p_attr_map = new cs_lagr_attribute_map_t(attr_keys);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return const pointer to the main particle attribute map structure.
 *
 * \return pointer to current particle attribute map, or nullptr
 */
/*----------------------------------------------------------------------------*/

const cs_lagr_attribute_map_t *
cs_lagr_particle_get_attr_map(void)
{
  const cs_lagr_attribute_map_t *p_am = _p_attr_map;
  return p_am;
}

/*----------------------------------------------------------------------------*/
/*!
 * Allocate main cs_lagr_particle_set_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_create(void)
{
  cs_glob_lagr_particle_set = new cs_lagr_particle_set_t(128, _p_attr_map);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER CREATION\n");
  cs_lagr_particle_set_dump(cs_glob_lagr_particle_set);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy main particle set and map if they exist.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_finalize(void)
{
  delete cs_glob_lagr_particle_set;
  delete _p_attr_map;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy attributes from one particle to another
 *
 * The random value associated with the particle is modified.
 *
 * \param  dest  id (0-based) of destination particle
 * \param  src   id (0-based) of source particle
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_part_copy(cs_lnum_t  dest,
                  cs_lnum_t  src)
{
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  memcpy(p_set->p_buffer + p_set->p_am->extents*(dest),
         p_set->p_buffer + p_set->p_am->extents*(src),
         p_set->p_am->extents);
  cs_real_t random = -1;
  cs_random_uniform(1, &random);
  cs_lagr_particles_set_real(p_set, dest, CS_LAGR_RANDOM_VALUE,
                             random);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data extents for a given particle attribute.
 *
 * For attributes not currently present, the displacement and data
 * size should be -1 and 0 respectively.
 *
 * \param[in]   p_set      pointer to particle set
 * \param[in]   time_id    associated time id (0: current, 1: previous)
 * \param[in]   attr       particle attribute
 * \param[out]  extents    size (in bytes) of particle structure, or nullptr
 * \param[out]  size       size (in bytes) of attribute in particle structure,
 *                         or nullptr
 * \param[out]  displ      displacement (in bytes) in particle structure,
 *                         or nullptr
 * \param[out]  datatype   datatype of associated attribute, or nullptr
 * \param[out]  count      number of type values associated with attribute,
 *                         or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_get_attr_info(const cs_lagr_particle_set_t  *p_set,
                      int                            time_id,
                      cs_lagr_attribute_t            attr,
                      size_t                        *extents,
                      size_t                        *size,
                      ptrdiff_t                     *displ,
                      cs_datatype_t                 *datatype,
                      int                           *count)
{
  if (extents)
    *extents = p_set->p_am->extents;
  if (size)
    *size = p_set->p_am->size[attr];
  if (displ)
    *displ = p_set->p_am->displ[time_id][attr];
  if (datatype)
    *datatype = p_set->p_am->datatype[attr];
  if (count)
    *count = p_set->p_am->count[time_id][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that query for a given particle attribute is valid.
 *
 * \param[in]   p_set             associated particle set
 * \param[in]   attr              attribute whose values are required
 * \param[in]   datatype          associated value type
 * \param[in]   stride            number of values per particle
 * \param[in]   component_id      if -1 : extract the whole attribute
 *                                if >0 : id of the component to extract
 *
 * \return 0 in case of success, 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_check_attr_query(const cs_lagr_particle_set_t  *p_set,
                         cs_lagr_attribute_t            attr,
                         cs_datatype_t                  datatype,
                         int                            stride,
                         int                            component_id)
{
  int retval = 0;

  assert(p_set != nullptr);

  int _count;
  cs_datatype_t _datatype;

  cs_lagr_get_attr_info(p_set, 0, attr,
                        nullptr, nullptr, nullptr, &_datatype, &_count);

  if (   datatype != _datatype || stride != _count
      || component_id < -1 || component_id >= stride) {

    char attr_name[128];
    attr_name[127] = '\0';
    const char *_attr_name = attr_name;
    if (attr < CS_LAGR_N_ATTRIBUTES) {
      snprintf(attr_name, 127, "CS_LAGR_%s", cs_lagr_attribute_name[attr]);
      size_t l = strlen(attr_name);
      for (size_t i = 0; i < l; i++)
        attr_name[i] = toupper(attr_name[i]);
    }
    else {
      snprintf(attr_name, 127, "%d", (int)attr);
    }

    if (datatype != _datatype || stride != _count)
      bft_error(__FILE__, __LINE__, 0,
                _("Attribute %s is of datatype %s and stride %d\n"
                  "but %s and %d were requested."),
                _attr_name,
                cs_datatype_name[_datatype], _count,
                cs_datatype_name[datatype], stride);

    else if (component_id < -1 || component_id >= stride)
      bft_error(__FILE__, __LINE__, 0,
                _("Attribute %s has a number of components equal to %d\n"
                  "but component %d is requested."),
                _attr_name,
                stride,
                component_id);

    retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a particle attribute is in a valid range.
 *
 * If this is not the case, a fatal error is provoked.

 * \param[in]   attr       particle attribute
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_attr_in_range(int  attr)
{
  if (attr < 0 || attr >= CS_LAGR_N_ATTRIBUTES)
    bft_error(__FILE__, __LINE__,0,
              _("Out-of range attribute type: %d"),
              (int)attr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main cs_lagr_particle_set_t structure.
 *
 * \return  pointer to current particle set, or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_set_t  *
cs_lagr_get_particle_set(void)
{
  return cs_glob_lagr_particle_set;
}

cs_lagr_particle_set_t&
cs_lagr_get_particle_set_ref(void)
{
  return *cs_glob_lagr_particle_set;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set reallocation factor for particle sets.
 *
 * This factor determines the multiplier used for reallocations when
 * the particle set's buffers are too small to handle the new number of
 * particles.
 *
 * \param[in]  f  reallocation size multiplier
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_reallocation_factor(double  f)
{
  if (f > 1)
    _reallocation_factor = f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global maximum number of particles.
 *
 * By default, the number is limited only by local \ref cs_lnum_t and global
 * \ref cs_gnum_t data representation limits.
 *
 * \return  global maximum number of particles
 */
/*----------------------------------------------------------------------------*/

unsigned long long
cs_lagr_get_n_g_particles_max(void)
{
  return _n_g_max_particles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global maximum number of particles.
 *
 * By default, the number is limited only by local \ref cs_lnum_t and global
 * \ref cs_gnum_t data representation limits.
 *
 * \param[in]  n_g_particles_max  global maximum number of particles
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_n_g_particles_max(unsigned long long  n_g_particles_max)
{
  _n_g_max_particles = n_g_particles_max;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy current attributes to previous attributes.
 *
 * \param[in, out]  p_set        reference to particle set
 * \param[in]       p_id         id of particle
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particles_current_to_previous(cs_lagr_particle_set_t  &p_set,
                                      cs_lnum_t                p_id)
{
  const cs_lagr_attribute_map_t  *p_am = p_set.p_am;
  unsigned char *p_buf = p_set.p_buffer + p_am->extents*(p_id);

  for (int i_attr = 0;
       i_attr < CS_LAGR_N_ATTRIBUTES;
       i_attr++) {
    auto attr = static_cast<cs_lagr_attribute_t>(i_attr);
    if (p_am->count[1][attr] > 0 && p_am->count[0][attr] > 0) {
      memcpy(p_buf + p_am->displ[1][attr],
             p_buf + p_am->displ[0][attr],
             p_am->size[attr]);
    }
  }
  *((cs_lnum_t *)(p_buf + p_am->displ[1][CS_LAGR_RANK_ID])) = cs_glob_rank_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_lagr_particle_set_t structure
 *
 * \param[in]  p_set  cs_lagr_particle_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_dump(const cs_lagr_particle_set_t  *p_set)
{
  if (p_set != nullptr) {

    bft_printf("Particle set\n");
    bft_printf("------------\n");
    bft_printf("  n_particles:      %10ld\n", (long)p_set->n_particles);
    bft_printf("  n_particles_max:  %10ld\n", (long)p_set->n_particles_max);

    bft_printf_flush();

    for (cs_lnum_t i = 0; i < p_set->n_particles; i++) {
      _dump_particle(p_set, i);
    }

  }
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set number of user particle variables.
 *
 * \param[in]  n_user_variables  number of user variables
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_n_user_variables(int  n_user_variables)
{
  cs_glob_lagr_model->n_user_variables = n_user_variables;
}

/*----------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*!
 * \brief Default constructor
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_attribute_map_t::cs_lagr_attribute_map_t()
{
  /* Private members */
  extents = 0;
  lb = 0;
  n_time_vals = 0;

  for (int i_attr = 0; i_attr < CS_LAGR_N_ATTRIBUTES; i_attr++) {
    cs_lagr_attribute_t attr = static_cast<cs_lagr_attribute_t>(i_attr);
    size[attr] = 0;
    datatype[attr] = CS_REAL_TYPE;
  }

  /* Public members (to be removed when possible) */
  private_to_public_();

}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Constructor based on a list of attribute keys
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST
cs_lagr_attribute_map_t::cs_lagr_attribute_map_t
(
  cs_lnum_t attr_keys[CS_LAGR_N_ATTRIBUTES][3]
)
{
  /* Start of buffer is used for private tracking state info */

  lb = _align_extents(sizeof(cs_lagr_tracking_info_t));
  extents =  lb;

  /* Currently, current and previous time values are managed */

  n_time_vals = 2;

  _displ.reshape(2, CS_LAGR_N_ATTRIBUTES);
  _count.reshape(2, CS_LAGR_N_ATTRIBUTES);

  for (int i_attr = 0; i_attr < CS_LAGR_N_ATTRIBUTES; i_attr++) {
    cs_lagr_attribute_t attr = static_cast<cs_lagr_attribute_t>(i_attr);
    size[attr] = 0;
    datatype[attr] = CS_REAL_TYPE;

    for (int time_id = 0; time_id < n_time_vals; time_id++) {
      _displ(time_id, attr) = -1;
      _count(time_id, attr) = 1;
    }
  }

  cs_array<cs_lnum_t> order(CS_LAGR_N_ATTRIBUTES);

  cs_order_lnum_allocated_s(nullptr,
                            (const cs_lnum_t *)attr_keys,
                            3,
                            order.data(),
                            CS_LAGR_N_ATTRIBUTES);

  /* Loop on available times */

  for (int time_id = 0; time_id < n_time_vals; time_id++) {

    int array_prev = 0;

    /* Now loop on ordered attributes */

    for (int i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {

      cs_datatype_t datatype_ = CS_REAL_TYPE;
      int min_time_id = 0;
      int max_time_id = 0;

      cs_lagr_attribute_t attr = static_cast<cs_lagr_attribute_t>(order[i]);

      if (time_id == 0)
        datatype[attr] = CS_DATATYPE_NULL;

      _displ(time_id, attr) = -1;
      _count(time_id, attr) = 0;

      if (attr_keys[attr][0] < 1) continue;

      /*
        integer values at current and previous time steps
        real values at current time step
        integer values at current time step */

      /* Behavior depending on array */

      switch(attr_keys[attr][0]) {
      case CS_LAGR_P_RVAR_TS:
      case CS_LAGR_P_RVAR:
        max_time_id = 1;
        break;
      case CS_LAGR_P_IVAR:
        datatype_ = CS_LNUM_TYPE;
        max_time_id = 1;
        break;
      case CS_LAGR_P_RPRP:
        break;
      case CS_LAGR_P_IPRP:
        datatype_ = CS_LNUM_TYPE;
        break;
      case CS_LAGR_P_RKID:
        datatype_ = CS_LNUM_TYPE;
        min_time_id = 1;
        max_time_id = 1;
        break;
      default:
        continue;
      }

      if (time_id < min_time_id || time_id > max_time_id)
        continue;

      /* Add padding for alignment when changing array */

      if (attr_keys[attr][0] != array_prev) {
        extents = _align_extents(extents);
        array_prev = attr_keys[attr][0];
      }

      /* Add attribute to map */

      _displ(time_id, attr) = extents;
      _count(time_id, attr) = attr_keys[attr][2];
      if (time_id == min_time_id) {
        datatype[attr] = datatype_;
        size[attr] = _count(time_id, attr)
                    * cs_datatype_size[datatype[attr]];
      }

      extents += size[attr];
    }

    extents = _align_extents(extents);
  }

  /* Add source terms for 2nd order */

  if (cs_glob_lagr_time_scheme->t_order > 1) {
    _source_term_displ.reshape(CS_LAGR_N_ATTRIBUTES);

    /* loop again on ordered attributes */

    for (int i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {

      cs_lagr_attribute_t attr = static_cast<cs_lagr_attribute_t>(order[i]);

      /* Add attribute to map if in CS_LAGR_P_RVARS array */

      if (   attr_keys[attr][0] == CS_LAGR_P_RVAR_TS
          && _count(0, attr) > 0) {
        _source_term_displ[attr] = extents;
        extents += size[attr];
      }
      else
        _source_term_displ[attr] = -1;
    }

    extents = _align_extents(extents);
  }

  /* Public members (to be removed when possible) */
  private_to_public_();

}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Copy constructor using a shallow copy (no allocation)
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_attribute_map_t::cs_lagr_attribute_map_t
(
  const cs_lagr_attribute_map_t& other
)
{
  extents = other.extents;
  lb = other.lb;
  n_time_vals = other.n_time_vals;

  for (int i_attr = 0; i_attr < CS_LAGR_N_ATTRIBUTES; i_attr++) {
    cs_lagr_attribute_t attr = static_cast<cs_lagr_attribute_t>(i_attr);
    size[attr] = other.size[attr];
    datatype[attr] = other.datatype[attr];
  }
  /* Non owner copies of arrays */
  _count = cs_array_2d<int>(other._count);
  _displ = cs_array_2d<ptrdiff_t>(other._displ);
  _source_term_displ = cs_array<ptrdiff_t>(other._source_term_displ);

  private_to_public_();
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Move constructor
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_attribute_map_t::cs_lagr_attribute_map_t
(
  cs_lagr_attribute_map_t&& other /*!<[in] Original reference to move */
)
{
  swap_members_(*this, other);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Default destructor
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_attribute_map_t::~cs_lagr_attribute_map_t()
{
  count = nullptr;
  displ = nullptr;
  source_term_displ = nullptr;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Swap method
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
void
cs_lagr_attribute_map_t::swap_members_
(
  cs_lagr_attribute_map_t& first,
  cs_lagr_attribute_map_t& second
)
{
  /* Swap the different members */
  cs::swap_objects(first.extents, second.extents);
  cs::swap_objects(first.lb, second.lb);
  cs::swap_objects(first.n_time_vals, second.n_time_vals);
  cs::swap_objects(first.size, second.size);
  cs::swap_objects(first.datatype, second.datatype);
  cs::swap_objects(first._count, second._count);
  cs::swap_objects(first.count, second.count);
  cs::swap_objects(first._displ, second._displ);
  cs::swap_objects(first.displ, second.displ);
  cs::swap_objects(first._source_term_displ, second._source_term_displ);
  cs::swap_objects(first.source_term_displ, second.source_term_displ);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Assignment operator
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_attribute_map_t&
cs_lagr_attribute_map_t::operator=(cs_lagr_attribute_map_t other)
{
  swap_members_(*this, other);

  return *this;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Update public pointers based on private structures.
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
void
cs_lagr_attribute_map_t::private_to_public_()
{
  using  lagr_attr_ptrdiff_t = ptrdiff_t[CS_LAGR_N_ATTRIBUTES];
  using  lagr_attr_int_t = int[CS_LAGR_N_ATTRIBUTES];

  count = _count.data<lagr_attr_int_t>();
  displ = _displ.data<lagr_attr_ptrdiff_t>();
  source_term_displ = _source_term_displ.data();
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Default constructor
 */
/*--------------------------------------------------------------------------*/

cs_lagr_particle_set_t::cs_lagr_particle_set_t()
{
  n_particles = 0;
  n_part_new  = 0;
  n_part_merged = 0;
  n_part_out    = 0;
  n_part_dep    = 0;
  n_part_fou    = 0;
  n_part_resusp = 0;
  n_failed_part = 0;

  weight = 0.0;
  weight_new = 0.0;
  weight_merged = 0.0;
  weight_out = 0.0;
  weight_dep = 0.0;
  weight_fou = 0.0;
  weight_resusp = 0.0;
  weight_failed = 0.0;
  n_particles_max = 0;

  private_to_public_();
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Constructor
 */
/*--------------------------------------------------------------------------*/

cs_lagr_particle_set_t::cs_lagr_particle_set_t
(
  cs_lnum_t                      n_part_max, /*!<[in] Maximal number of particles */
  const cs_lagr_attribute_map_t *map         /*!<[in] Attribute map */
)
{
  assert(n_part_max >= 1);

  n_particles = 0;
  n_part_new  = 0;
  n_part_merged = 0;
  n_part_out    = 0;
  n_part_dep    = 0;
  n_part_fou    = 0;
  n_part_resusp = 0;
  n_failed_part = 0;

  weight = 0.0;
  weight_new = 0.0;
  weight_merged = 0.0;
  weight_out = 0.0;
  weight_dep = 0.0;
  weight_fou = 0.0;
  weight_resusp = 0.0;
  weight_failed = 0.0;
  n_particles_max = n_part_max;

  _p_am = cs_lagr_attribute_map_t(*map);
  _p_buffer.reshape(n_part_max * map->extents);
  private_to_public_();
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Copy constructor (shallow copy)
 */
/*--------------------------------------------------------------------------*/

cs_lagr_particle_set_t::cs_lagr_particle_set_t
(
  const cs_lagr_particle_set_t& other
)
{
  n_particles = other.n_particles;
  n_part_new = other.n_part_new;
  n_part_merged = other.n_part_merged;
  n_part_out = other.n_part_out;
  n_part_dep = other.n_part_dep;
  n_part_fou = other.n_part_fou;
  n_part_resusp = other.n_part_resusp;
  n_failed_part = other.n_failed_part;

  weight = other.weight;
  weight_new = other.weight_new;
  weight_merged = other.weight_merged;
  weight_out = other.weight_out;
  weight_dep = other.weight_dep;
  weight_fou = other.weight_fou;
  weight_resusp = other.weight_resusp;
  weight_failed = other.weight_failed;
  n_particles_max = other.n_particles_max;

  _p_buffer = cs_array<unsigned char>(other._p_buffer);
  _p_am = cs_lagr_attribute_map_t(other._p_am);

  private_to_public_();
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Move constructor
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_particle_set_t::cs_lagr_particle_set_t
(
  cs_lagr_particle_set_t&& other
)
{
  swap_members_(*this, other);
}
/*--------------------------------------------------------------------------*/
/*!
 * \brief Default destructor
 */
/*--------------------------------------------------------------------------*/

cs_lagr_particle_set_t::~cs_lagr_particle_set_t()
{
  p_buffer = nullptr;
  p_am = nullptr;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Swap method
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
void
cs_lagr_particle_set_t::swap_members_
(
  cs_lagr_particle_set_t& first,
  cs_lagr_particle_set_t& second
)
{
  cs::swap_objects(first.n_particles, second.n_particles);
  cs::swap_objects(first.n_part_new, second.n_part_new);
  cs::swap_objects(first.n_part_merged, second.n_part_merged);
  cs::swap_objects(first.n_part_out, second.n_part_out);
  cs::swap_objects(first.n_part_dep, second.n_part_dep);
  cs::swap_objects(first.n_part_fou, second.n_part_fou);
  cs::swap_objects(first.n_part_resusp, second.n_part_resusp);
  cs::swap_objects(first.n_failed_part, second.n_failed_part);

  cs::swap_objects(first.weight, second.weight);
  cs::swap_objects(first.weight_new, second.weight_new);
  cs::swap_objects(first.weight_merged, second.weight_merged);
  cs::swap_objects(first.weight_out, second.weight_out);
  cs::swap_objects(first.weight_dep, second.weight_dep);
  cs::swap_objects(first.weight_fou, second.weight_fou);
  cs::swap_objects(first.weight_resusp, second.weight_resusp);
  cs::swap_objects(first.weight_failed, second.weight_failed);
  cs::swap_objects(first.n_particles_max, second.n_particles_max);

  cs::swap_objects(first._p_buffer, second._p_buffer);
  cs::swap_objects(first.p_buffer, second.p_buffer);
  cs::swap_objects(first._p_am, second._p_am);
  cs::swap_objects(first.p_am, second.p_am);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Resize the particle set
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST
int
cs_lagr_particle_set_t::resize
(
  const cs_lnum_t n_particles_max_min /*!<[in] minimum local max number of particles */
)
{
  int retval = 0;
  assert(n_particles_max_min >= 0);

  /* Do we have a limit ? */

  if (_n_g_max_particles < ULLONG_MAX) {
    cs_gnum_t _n_g_min_particles = n_particles_max_min;
    cs_parall_counter(&_n_g_min_particles, 1);
    if (_n_g_min_particles > _n_g_max_particles)
      retval = -1;
  }
  else if (n_particles_max < n_particles_max_min) {

    if (n_particles_max == 0)
      n_particles_max = 1;

    while (n_particles_max < n_particles_max_min)
      n_particles_max *= _reallocation_factor;

    _p_buffer.reshape(n_particles_max * _p_am.extents);

    private_to_public_();

    retval = 1;
  }

  return retval;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Assignment operator
 */
/*--------------------------------------------------------------------------*/

CS_F_HOST_DEVICE
cs_lagr_particle_set_t&
cs_lagr_particle_set_t::operator=(cs_lagr_particle_set_t other)
{
  swap_members_(*this, other);

  return *this;
}

CS_F_HOST_DEVICE
void
cs_lagr_particle_set_t::private_to_public_()
{
  p_buffer = _p_buffer.data();
  p_am = &(_p_am);
}

