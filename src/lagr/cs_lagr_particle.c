/*============================================================================
 * Lagrangian particle model
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

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_math.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_random.h"
#include "cs_timer_stats.h"

#include "cs_lagr.h"
#include "cs_lagr_post.h"
#include "cs_lagr_clogging.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_dlvo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_particle.c
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

   rvar real values at current an previous time steps
   ivar integer values at current and previous time steps
   rprp real properties at current an previous time steps
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
  "tr_truncate",
  "tr_reposition",
  "turb_state_1",
  "pred_velocity",
  "pred_velocity_seen",
  "v_gauss",
  "br_gauss",
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
  "user",
  "<none>"};

/* Global particle attributes map */

static cs_lagr_attribute_map_t  *_p_attr_map = NULL;

/* Particle set reallocation parameters */

static  double              _reallocation_factor = 2.0;
static  unsigned long long  _n_g_max_particles = ULLONG_MAX;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointer to the main particle set */

cs_lagr_particle_set_t *cs_glob_lagr_particle_set = NULL;

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

/*----------------------------------------------------------------------------*
 * Map particle attributes for a given configuration.
 *
 * parameters:
 *   attr_keys   <-> keys to sort attributes by Fortran array and index
 *                   for each attribute: array, index in array, count
 *
 * returns:
 *   pointer to structure mapping particle attributes
 *----------------------------------------------------------------------------*/

static cs_lagr_attribute_map_t *
_create_attr_map(cs_lnum_t attr_keys[CS_LAGR_N_ATTRIBUTES][3])
{
  cs_lagr_attribute_t attr;
  cs_lnum_t *order;

  cs_lagr_attribute_map_t  *p_am;

  BFT_MALLOC(p_am, 1, cs_lagr_attribute_map_t);

  /* Start of buffer is used for private tracking state info */

  p_am->lb = _align_extents(sizeof(cs_lagr_tracking_info_t));

  p_am->extents = p_am->lb;

  /* Currently, current and previous time values are managed */

  p_am->n_time_vals = 2;

  if (true) { /* Allocation requires typedef due to cast */

    typedef ptrdiff_t lagr_attr_ptrdiff_t[CS_LAGR_N_ATTRIBUTES];
    typedef int       lagr_attr_int_t[CS_LAGR_N_ATTRIBUTES];

    BFT_MALLOC(p_am->displ, 2, lagr_attr_ptrdiff_t);
    BFT_MALLOC(p_am->count, 2, lagr_attr_int_t);

  }

  else { /* Variant:
            to avoid issue with cast and no typdef, use lower level function */

    p_am->displ = bft_mem_malloc(2, sizeof(p_am->displ[0]),
                                 "p_am->displ", __FILE__, __LINE__);
    p_am->count = bft_mem_malloc(2, sizeof(p_am->count[0]),
                                 "p_am->count", __FILE__, __LINE__);

  }

  p_am->source_term_displ = NULL;

  for (attr = 0; attr < CS_LAGR_N_ATTRIBUTES; attr++) {
    p_am->size[attr] = 0;
    p_am->datatype[attr] = CS_REAL_TYPE;
    for (int time_id = 0; time_id < p_am->n_time_vals; time_id++) {
      p_am->displ[time_id][attr] = -1;
      p_am->count[time_id][attr] = 1;
    }
  }

  BFT_MALLOC(order, CS_LAGR_N_ATTRIBUTES, cs_lnum_t);

  cs_order_lnum_allocated_s(NULL,
                            (const cs_lnum_t *)attr_keys,
                            3,
                            order,
                            CS_LAGR_N_ATTRIBUTES);

  /* Loop on available times */

  for (int time_id = 0; time_id < p_am->n_time_vals; time_id++) {

    int array_prev = 0;

    /* Now loop on ordered attributes */

    for (int i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {

      cs_datatype_t datatype = CS_REAL_TYPE;
      int min_time_id = 0;
      int max_time_id = 0;

      attr = order[i];

      if (time_id == 0)
        p_am->datatype[attr] = CS_DATATYPE_NULL;
      p_am->displ[time_id][attr] =-1;
      p_am->count[time_id][attr] = 0;

      if (attr_keys[attr][0] < 1) continue;

      /*
        ieptp/ieptpa integer values at current and previous time steps
        pepa real values at current time step
        ipepa integer values at current time step */

      /* Behavior depending on array */

      switch(attr_keys[attr][0]) {
      case CS_LAGR_P_RVAR_TS:
      case CS_LAGR_P_RVAR:
        max_time_id = 1;
        break;
      case CS_LAGR_P_IVAR:
        datatype = CS_LNUM_TYPE;
        max_time_id = 1;
        break;
      case CS_LAGR_P_RPRP:
        break;
      case CS_LAGR_P_IPRP:
        datatype = CS_LNUM_TYPE;
        break;
      case CS_LAGR_P_RKID:
        datatype = CS_LNUM_TYPE;
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
        p_am->extents = _align_extents(p_am->extents);
        array_prev = attr_keys[attr][0];
      }

      /* Add attribute to map */

      p_am->displ[time_id][attr] = p_am->extents;
      p_am->count[time_id][attr] = attr_keys[attr][2];
      if (time_id == min_time_id) {
        p_am->datatype[attr] = datatype;
        p_am->size[attr] =   p_am->count[time_id][attr]
                           * cs_datatype_size[p_am->datatype[attr]];
      }

      p_am->extents += p_am->size[attr];

    }

    p_am->extents = _align_extents(p_am->extents);

  }

  /* Add source terms for 2nd order */

  if (cs_glob_lagr_time_scheme->t_order > 1) {

    BFT_MALLOC(p_am->source_term_displ, CS_LAGR_N_ATTRIBUTES, ptrdiff_t);

    /* loop again on ordered attributes */

    for (int i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {

      attr = order[i];

      /* Add attribute to map if in CS_LAGR_P_RVARS array */

      if (   attr_keys[attr][0] == CS_LAGR_P_RVAR_TS
          && p_am->count[0][attr] > 0) {
        p_am->source_term_displ[attr] = p_am->extents;
        p_am->extents += p_am->size[attr];
      }
      else
        p_am->source_term_displ[attr] = -1;

    }

    p_am->extents = _align_extents(p_am->extents);

  }

  BFT_FREE(order);

  return p_am;
}

/*----------------------------------------------------------------------------*
 * Free particle attributes for a given configuration.
 *----------------------------------------------------------------------------*/

static void
_destroy_attr_map(cs_lagr_attribute_map_t  **p_am)
{
  if (*p_am != NULL) {
    cs_lagr_attribute_map_t  *_p_am = *p_am;

    BFT_FREE(_p_am->source_term_displ);

    BFT_FREE(_p_am->displ);
    BFT_FREE(_p_am->count);

    BFT_FREE(*p_am);
  }
}

/*----------------------------------------------------------------------------
 * Allocate a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *   n_particles_max <-- local max. number of particles
 *   p_am            <-- particle attributes map
 *
 * returns:
 *   a new allocated cs_lagr_particle_set_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_particle_set_t *
_create_particle_set(cs_lnum_t                       n_particles_max,
                     const cs_lagr_attribute_map_t  *p_am)

{
  cs_lagr_particle_set_t  *new_set = NULL;

  if (n_particles_max == 0)
    return NULL;

  BFT_MALLOC(new_set, 1, cs_lagr_particle_set_t);

  BFT_MALLOC(new_set->p_buffer, n_particles_max * p_am->extents, unsigned char);

  new_set->n_particles = 0;
  new_set->n_part_new = 0;
  new_set->n_part_merged = 0;
  new_set->n_part_out = 0;
  new_set->n_part_dep = 0;
  new_set->n_part_fou = 0;
  new_set->n_part_resusp = 0;
  new_set->n_failed_part = 0;

  new_set->weight = 0.0;
  new_set->weight_new = 0.0;
  new_set->weight_merged = 0.0;
  new_set->weight_out = 0.0;
  new_set->weight_dep = 0.0;
  new_set->weight_fou = 0.0;
  new_set->weight_resusp = 0.0;
  new_set->weight_failed = 0.0;

  new_set->n_particles_max = n_particles_max;

  assert(n_particles_max >= 1);

  new_set->p_am = p_am;

  return new_set;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *   particle_set  <->  a cs_lagr_particle_set_t structure
 *----------------------------------------------------------------------------*/

static void
_destroy_particle_set(cs_lagr_particle_set_t **set)
{
  if (set != NULL) {

    cs_lagr_particle_set_t *_set = *set;
    BFT_FREE(_set->p_buffer);

    BFT_FREE(*set);
  }
}

/*----------------------------------------------------------------------------
 * Dump a particle structure
 *
 * parameter
 *   particles   <-- cs_lagr_particle_set_t structure to dump
 *   particle_id <-- id of particle to dump
 *----------------------------------------------------------------------------*/

static void
_dump_particle(const cs_lagr_particle_set_t  *particles,
               cs_lnum_t                      particle_id)
{
  const unsigned char *p =   particles->p_buffer
                           + particles->p_am->extents*particle_id;
  const cs_lagr_attribute_map_t *am = particles->p_am;

  bft_printf("  particle: %lu\n", (unsigned long)particle_id);

  for (int time_id = 0; time_id < particles->p_am->n_time_vals; time_id++) {

    if (time_id == 0)
      bft_printf("    values at time n:\n");
    else
      bft_printf("    values at time: n-%d\n", time_id);

    for (cs_lagr_attribute_t attr = 0;
         attr < CS_LAGR_N_ATTRIBUTES;
         attr++) {
      if (am->count[time_id][attr] > 0) {
        const char *attr_name = cs_lagr_attribute_name[attr];
        switch (am->datatype[attr]) {
        case CS_LNUM_TYPE:
          {
            const cs_lnum_t *v
              = cs_lagr_particle_attr_n_const(p, particles->p_am, time_id, attr);
            bft_printf("      %24s: %10ld\n", attr_name, (long)v[0]);
            for (int i = 1; i < am->count[time_id][attr]; i++)
              bft_printf("      %24s: %10ld\n", " ", (long)v[i]);
          }
          break;
        case CS_GNUM_TYPE:
          {
            const cs_gnum_t *v
              = cs_lagr_particle_attr_n_const(p, particles->p_am, time_id, attr);
            bft_printf("      %24s: %10lu\n", attr_name, (unsigned long)v[0]);
            for (int i = 1; i < am->count[time_id][attr]; i++)
              bft_printf("      %24s: %10lu\n", " ", (unsigned long)v[i]);
        }
          break;
        case CS_REAL_TYPE:
          {
            const cs_real_t *v
              = cs_lagr_particle_attr_n_const(p, particles->p_am, time_id, attr);
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

/*----------------------------------------------------------------------------
 * Resize a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *   particle_set        <-> pointer to a cs_lagr_particle_set_t structure
 *   n_particles_max_min <-- minimum local max. number of particles
 *
 * returns:
 *   1 if resizing was required, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_particle_set_resize(cs_lagr_particle_set_t   *particle_set,
                     const cs_lnum_t           n_particles_max_min)
{
  int retval = 0;

  assert(n_particles_max_min >= 0);

  if (particle_set->n_particles_max < n_particles_max_min) {

    if (particle_set->n_particles_max == 0)
      particle_set->n_particles_max = 1;

    while (particle_set->n_particles_max < n_particles_max_min)
      particle_set->n_particles_max *= _reallocation_factor;

    BFT_REALLOC(particle_set->p_buffer,
                particle_set->n_particles_max * particle_set->p_am->extents,
                unsigned char);

    retval = 1;
  }

  return retval;
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
  const  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  int  i;

  int loc_count = 0;

  int pepa_loc_add = 1000; /* should be abore any j* pointer if used */

  cs_lnum_t attr_keys[CS_LAGR_N_ATTRIBUTES][3];

  /* Initialize global parameter relative to the lagrangian module */

  /*  cs_glob_lagr_brownian->lamvbr = *lamvbr; */

  if (lagr_model->physical_model == 2)
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
  attr_keys[CS_LAGR_RANK_ID][1] = 1;

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

  attr_keys[CS_LAGR_COORDS][1] = ++loc_count;
  attr_keys[CS_LAGR_COORDS][2] = 3;

  attr_keys[CS_LAGR_VELOCITY][1] = ++loc_count;
  attr_keys[CS_LAGR_VELOCITY][2] = 3;

  attr_keys[CS_LAGR_VELOCITY_SEEN][1] = ++loc_count;
  attr_keys[CS_LAGR_VELOCITY_SEEN][2] = 3;

  attr_keys[CS_LAGR_TR_TRUNCATE][0] = CS_LAGR_P_RPRP;
  attr_keys[CS_LAGR_TR_TRUNCATE][1] = ++loc_count;

  attr_keys[CS_LAGR_TR_REPOSITION][0] = CS_LAGR_P_IPRP;
  attr_keys[CS_LAGR_TR_REPOSITION][1] = ++loc_count;

  if (lagr_time_scheme->t_order > 1) {
    attr_keys[CS_LAGR_TAUP_AUX][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_TAUP_AUX][1] = ++pepa_loc_add;

    attr_keys[CS_LAGR_TURB_STATE_1][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_TURB_STATE_1][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_TURB_STATE_1][2] = 3;

    attr_keys[CS_LAGR_PRED_VELOCITY][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_PRED_VELOCITY][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_PRED_VELOCITY][2] = 3;

    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][0] = CS_LAGR_P_RPRP;
    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][2] = 3;

    if (cs_glob_lagr_time_scheme->idistu == 1) {
      attr_keys[CS_LAGR_V_GAUSS][0] = CS_LAGR_P_RPRP;
      attr_keys[CS_LAGR_V_GAUSS][1] = ++pepa_loc_add;
      attr_keys[CS_LAGR_V_GAUSS][2] = 9;
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

  if (lagr_model->physical_model == 1) {

    if (cs_glob_lagr_specific_physics->itpvar == 1) {

      attr_keys[CS_LAGR_CP][1] = ++loc_count;

      attr_keys[CS_LAGR_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
      attr_keys[CS_LAGR_TEMPERATURE][1] = ++loc_count;

      attr_keys[CS_LAGR_FLUID_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
      attr_keys[CS_LAGR_FLUID_TEMPERATURE][1] = ++loc_count;

      if (extra->radiative_model > 0)
        attr_keys[CS_LAGR_EMISSIVITY][1] = ++loc_count;

    }

  }

  if (lagr_model->physical_model == 2) {

    attr_keys[CS_LAGR_CP][1] = ++loc_count;

    attr_keys[CS_LAGR_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
    attr_keys[CS_LAGR_TEMPERATURE][1] = ++loc_count;
    attr_keys[CS_LAGR_TEMPERATURE][2]
      = lagr_model->n_temperature_layers;

    attr_keys[CS_LAGR_FLUID_TEMPERATURE][0] = CS_LAGR_P_RVAR_TS;
    attr_keys[CS_LAGR_FLUID_TEMPERATURE][1] = ++loc_count;

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

  if (lagr_model->n_stat_classes > 0) {
    attr_keys[CS_LAGR_STAT_CLASS][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_STAT_CLASS][1] = ++loc_count;
  }

  if (lagr_model->n_particle_aggregates > 0) {
    attr_keys[CS_LAGR_PARTICLE_AGGREGATE][0] = CS_LAGR_P_IPRP;
    attr_keys[CS_LAGR_PARTICLE_AGGREGATE][1] = ++loc_count;
  }

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

  _p_attr_map = _create_attr_map(attr_keys);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return const pointer to the main particle attribute map structure.
 *
 * \return pointer to current particle attrbute map, or NULL
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
  cs_glob_lagr_particle_set = _create_particle_set(128, _p_attr_map);

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
  _destroy_particle_set(&cs_glob_lagr_particle_set);

  _destroy_attr_map(&_p_attr_map);
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
  cs_lagr_particle_set_t  *particles = cs_glob_lagr_particle_set;
  memcpy(particles->p_buffer + particles->p_am->extents*(dest),
         particles->p_buffer + particles->p_am->extents*(src),
         particles->p_am->extents);
  cs_real_t random = -1;
  cs_random_uniform(1, &random);
  cs_lagr_particles_set_real(particles, (dest-1), CS_LAGR_RANDOM_VALUE,
                             random);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data extents for a given particle attribute.
 *
 * For attributes not currently present, the displacement and data
 * size should be -1 and 0 respectively.
 *
 * \param[in]   particles  associated particle set
 * \param[in]   time_id    associated time id (0: current, 1: previous)
 * \param[in]   attr       particle attribute
 * \param[out]  extents    size (in bytes) of particle structure, or NULL
 * \param[out]  size       size (in bytes) of attribute in particle structure,
 *                         or NULL
 * \param[out]  displ      displacement (in bytes) in particle structure,
 *                         or NULL
 * \param[out]  datatype   datatype of associated attribute, or NULL
 * \param[out]  count      number of type values associated with attribute,
 *                         or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_get_attr_info(const cs_lagr_particle_set_t  *particles,
                      int                            time_id,
                      cs_lagr_attribute_t            attr,
                      size_t                        *extents,
                      size_t                        *size,
                      ptrdiff_t                     *displ,
                      cs_datatype_t                 *datatype,
                      int                           *count)
{
  if (extents)
    *extents = particles->p_am->extents;
  if (size)
    *size = particles->p_am->size[attr];
  if (displ)
    *displ = particles->p_am->displ[time_id][attr];
  if (datatype)
    *datatype = particles->p_am->datatype[attr];
  if (count)
    *count = particles->p_am->count[time_id][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that query for a given particle attribute is valid.
 *
 * \param[in]   particles             associated particle set
 * \param[in]   attr                  attribute whose values are required
 * \param[in]   datatype              associated value type
 * \param[in]   stride                number of values per particle
 * \param[in]   component_id          if -1 : extract the whole attribute
 *                                    if >0 : id of the component to extract
 *
 * \return 0 in case of success, 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_check_attr_query(const cs_lagr_particle_set_t  *particles,
                         cs_lagr_attribute_t            attr,
                         cs_datatype_t                  datatype,
                         int                            stride,
                         int                            component_id)
{
  int retval = 0;

  assert(particles != NULL);

  int _count;
  cs_datatype_t _datatype;

  cs_lagr_get_attr_info(particles, 0, attr,
                        NULL, NULL, NULL, &_datatype, &_count);

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
 * \return  pointer to current particle set, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_set_t  *
cs_lagr_get_particle_set(void)
{
  return cs_glob_lagr_particle_set;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Resize particle set buffers if needed.
 *
 * By default, the total number of particles is not limited. A global limit
 * may be set using \ref cs_lagr_set_n_g_particles_max.
 *
 * \param[in]  n_min_particles  minimum number of particles required
 *
 * \return  1 if resizing was required, -1 if the global minimum number
 *          of particles would exceed the global limit, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_particle_set_resize(cs_lnum_t  n_min_particles)
{
  int retval = 0;

  /* Do we have a limit ? */

  if (_n_g_max_particles < ULLONG_MAX) {
    cs_gnum_t _n_g_min_particles = n_min_particles;
    cs_parall_counter(&_n_g_min_particles, 1);
    if (_n_g_min_particles > _n_g_max_particles)
      retval = -1;
  }
  else
    retval = _particle_set_resize(cs_glob_lagr_particle_set, n_min_particles);

  return retval;
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
 * \param[in, out]  particles     associated particle set
 * \param[in]       particle_id  id of particle
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particles_current_to_previous(cs_lagr_particle_set_t  *particles,
                                      cs_lnum_t                particle_id)
{
  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  unsigned char *p_buf = particles->p_buffer + p_am->extents*(particle_id);

  for (cs_lagr_attribute_t attr = 0;
       attr < CS_LAGR_N_ATTRIBUTES;
       attr++) {
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
 * \param[in]  particles  cs_lagr_particle_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_dump(const cs_lagr_particle_set_t  *particles)
{
  if (particles != NULL) {

    bft_printf("Particle set\n");
    bft_printf("------------\n");
    bft_printf("  n_particles:      %10d\n", particles->n_particles);
    bft_printf("  n_particles_max:  %10d\n", particles->n_particles_max);

    bft_printf_flush();

    for (cs_lnum_t i = 0; i < particles->n_particles; i++) {
      _dump_particle(particles, i);
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

END_C_DECLS
