/*============================================================================
 * Methods for particle localization
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Functions dealing with the particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

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
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_search.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr_utils.h"
#include "cs_lagr_clogging.h"
#include "cs_lagr_roughness.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define  N_GEOL 13
#define  CS_LAGR_MIN_COMM_BUF_SIZE  8
#define  CS_LAGR_MAX_PROPAGATION_LOOPS  30

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/* Tracking error types */

typedef enum {

  CS_LAGR_TRACKING_OK,

  CS_LAGR_TRACKING_ERR_INSIDE_B_FACES,
  CS_LAGR_TRACKING_ERR_MAX_LOOPS,
  CS_LAGR_TRACKING_ERR_DISPLACEMENT_I_FACE,
  CS_LAGR_TRACKING_ERR_DISPLACEMENT_B_FACE

} cs_lagr_tracking_error_t;

/* Linked list */
/* ----------- */

struct _cs_lagr_tracking_list_t {

  cs_lnum_t   prev_id;  /* id in particle set of the previous particle */
  cs_lnum_t   next_id;  /* id in particle set of the next particle */

};

/* Base particle description */
/* ------------------------- */

struct _cs_lagr_particle_t {

  cs_lnum_t   cur_cell_num;    /* current local cell number */
  cs_lnum_t   last_face_num;

  cs_lnum_t   switch_order_1;
  cs_lnum_t   state;         /* < 0 : - number of the boundary face where
                                      the particle is kept
                                0   : particle has to be destroyed
                                1   : particle has to be synchronized
                                2   : particle treated. End of displacement */

  cs_real_t   random_value;   /* random value associated with the particle */

  cs_real_t   stat_weight;
  cs_real_t   residence_time;
  cs_real_t   mass;
  cs_real_t   diameter;
  cs_real_t   taup_aux;
  cs_real_t   coord[3];
  cs_real_t   velocity[3];
  cs_real_t   velocity_seen[3];

  /* Deposition submodel additional parameters */

  cs_real_t   yplus;
  cs_real_t   interf;
  cs_lnum_t   close_face_id;
  cs_lnum_t   marko_val;
  cs_lnum_t   depo;                  /* jdepo   */
  cs_lnum_t   rank_flag;

  /* Resuspension model additional parameters */

  cs_lnum_t   nb_large_asperities;   /* jnbasg  */
  cs_lnum_t   nb_small_asperities;   /* jnbasg  */
  cs_real_t   adhesion_force;        /* jfadh   */
  cs_real_t   adhesion_torque;       /* jmfadh  */
  cs_real_t   displacement_norm;     /* jndisp  */

  /* Thermal model additional parameters */

  cs_real_t   temp[CS_LAGR_N_LAYERS]; /* jhp */
  cs_real_t   fluid_temp;             /* jtf */
  cs_real_t   cp;                     /* jcp */

  /* Coal combustion additional parameters */

  cs_real_t   water_mass;                  /* jmwat */
  cs_real_t   coal_mass[CS_LAGR_N_LAYERS]; /* jmch  */
  cs_real_t   coke_mass[CS_LAGR_N_LAYERS]; /* jmck  */

  cs_real_t   shrinking_diam;  /* jrdck */
  cs_real_t   initial_diam;    /* jrd0p */

  cs_lnum_t   coal_number;                    /* jinch  */
  cs_real_t   coal_density[CS_LAGR_N_LAYERS]; /* jrhock */

  /* Radiative model additional parameters */

  cs_real_t   emissivity;      /* jreps */

};

/* Particle data value */
/*---------------------*/

union cs_lagr_value_t {
  cs_lnum_t      l; /* v_lnum_t */
  cs_gnum_t      g; /* v_gnum_t */
  cs_real_t      f; /* v_real_t */
};

/* face_yplus auxiliary type */
/* ------------------------- */

typedef struct {

  cs_real_t  yplus;
  cs_lnum_t  face_id;

} face_yplus_t;

/* Structures useful to manage the exchange of particles
   between communicating ranks */

typedef struct {

  cs_lnum_t   n_cells;        /* Number of cells in the halo */

  cs_lnum_t  *rank;           /* value between [0, n_c_domains-1]
                                 (cf. cs_halo.h) */
  cs_lnum_t  *dist_cell_num;  /* local cell num. on distant ranks */
  cs_lnum_t  *transform_id;   /* In case of periodicity, transformation
                                 associated to a given halo cell */

  /* Buffer used to exchange particle between communicating ranks */

  size_t      send_buf_size;  /* Current maximum send buffer size */
  size_t      recv_buf_size;  /* Current maximum send buffer size */
  size_t      cur_extents;    /* Extents for current particle set */
  size_t      prev_extents;   /* Extents for previous particle set */

  cs_lnum_t  *send_count;     /* number of particles to send to
                                 each communicating rank */
  cs_lnum_t  *recv_count;     /* number of particles to receive from
                                 each communicating rank */

  cs_lnum_t  *send_shift;
  cs_lnum_t  *recv_shift;

  unsigned char  *send_buf;
  unsigned char  *recv_buf;

#if defined(HAVE_MPI)
  MPI_Request  *request;
  MPI_Status   *status;
#endif

} cs_lagr_halo_t;

/* Structures useful to build and manage the Lagrangian computation:
   - exchanging of particles between communicating ranks
   - finding the next cells where the particle moves on to
   - controlling the flow of particles coming in/out
*/

typedef struct {

  cs_lnum_t  max_face_connect_size;
  cs_lnum_t *face_connect_buffer;

  /* Cell -> Face connectivity */

  cs_lnum_t  *cell_face_idx;
  cs_lnum_t  *cell_face_lst;

  cs_lagr_halo_t    *halo;   /* Lagrangian halo structure */

  cs_interface_set_t  *face_ifs;

} cs_lagr_track_builder_t;

/* Structures useful to deal with boundary conditions
   For USLABO => _boundary_track_treatment */

typedef struct {

  cs_lnum_t   n_b_zones;  /* NFRLAG */
  cs_lnum_t   n_b_max_zones;

  cs_lnum_t  *b_zone_lst; /* ILFLAG */
  cs_lnum_t  *b_zone_classes; /* IUSNCL */
  cs_lnum_t  *b_zone_natures; /* IUSCLB */

  cs_lnum_t  *b_face_zone_num; /* IFRLAG */

  cs_lnum_t   continuous_injection;  /* INJCON */
  bool        steady_bndy_conditions;

  cs_real_t  *particle_flow_rate; /* DEBLAG -> post-processing use */

} cs_lagr_bdy_condition_t;

typedef struct {

  int  physic_mode;  /* FIXME: => enum: CS_LAGR_PHYS_STD,
                        CS_LAGR_PHYS_COAL,
                        CS_LAGR_PHYS_HEAT... */

  int  n_temperature_layers;

  int  deposition;
  int  rough;
  int  resuspension;
  int  clogging;

  int  n_stat_classes;
  int  n_user_variables;

} cs_lagr_param_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Kind of boundary type we can encounter */

enum {
  CS_LAGR_IENTRL = 1,
  CS_LAGR_ISORTL = 2,
  CS_LAGR_IREBOL = 3,
  CS_LAGR_IDEPO1 = 4,
  CS_LAGR_IDEPO2 = 5,
  CS_LAGR_IENCRL = 7,
  CS_LAGR_IDEPFA = 13,
  CS_LAGR_ISYMTL = 14,
};

/* State where a particle can be. */

enum {
  CS_LAGR_PART_TO_DELETE = 0,
  CS_LAGR_PART_TO_SYNC   = 1,
  CS_LAGR_PART_TREATED   = 2,
  CS_LAGR_PART_STICKED   = 3,
  CS_LAGR_PART_OUT       = 4,
  CS_LAGR_PART_ERR       = 5
};

/* Physical state where a particle can be. */

enum {
  CS_LAGR_PART_IN_FLOW        = 0,
  CS_LAGR_PART_DEPOSITED      = 1,
  CS_LAGR_PART_ROLLING        = 2,
  CS_LAGR_PART_NO_MOTION      = 10,
};

enum {
  CS_LAGR_PART_MOVE_OFF = 0,
  CS_LAGR_PART_MOVE_ON  = 1
};

/* According to the scheme order is degenerated to order 1 */

enum {
  CS_LAGR_SWITCH_OFF = 0,
  CS_LAGR_SWITCH_ON = 1
};

/* Enumerator names */

const char *cs_lagr_attribute_name[] = {
  "CS_LAGR_CUR_CELL_NUM",
  "CS_LAGR_LAST_FACE_NUM",
  "CS_LAGR_SWITCH_ORDER_1",
  "CS_LAGR_STATE",
  "CS_LAGR_RANDOM_VALUE",
  "CS_LAGR_STAT_WEIGHT",
  "CS_LAGR_RESIDENCE_TIME",
  "CS_LAGR_MASS",
  "CS_LAGR_DIAMETER",
  "CS_LAGR_TAUP_AUX",
  "CS_LAGR_COORDS",
  "CS_LAGR_VELOCITY",
  "CS_LAGR_VELOCITY_SEEN",
  "CS_LAGR_YPLUS",
  "CS_LAGR_INTERF",
  "CS_LAGR_NEIGHBOR_FACE_ID",
  "CS_LAGR_MARKO_VALUE",
  "CS_LAGR_DEPOSITION_FLAG",
  "CS_LAGR_RANK_FLAG",
  "CS_LAGR_N_LARGE_ASPERITIES",
  "CS_LAGR_N_SMALL_ASPERITIES",
  "CS_LAGR_ADHESION_FORCE",
  "CS_LAGR_ADHESION_TORQUE",
  "CS_LAGR_DISPLACEMENT_NORM",
  "CS_LAGR_TEMPERATURE",
  "CS_LAGR_FLUID_TEMPERATURE",
  "CS_LAGR_CP",
  "CS_LAGR_WATER_MASS",
  "CS_LAGR_COAL_MASS",
  "CS_LAGR_COKE_MASS",
  "CS_LAGR_SHRINKING_DIAMETER",
  "CS_LAGR_INITIAL_DIAMETER",
  "CS_LAGR_COAL_NUM",
  "CS_LAGR_COAL_DENSITY",
  "CS_LAGR_EMISSIVITY",
  "CS_LAGR_N_ATTRIBUTES"};

/* Global variable for the current subroutines */

static  cs_lagr_particle_set_t  *_particle_set = NULL;
static  cs_lagr_particle_set_t  *_prev_particle_set = NULL;
static  cs_lagr_track_builder_t  *_particle_track_builder = NULL;
static  cs_lagr_bdy_condition_t  *_lagr_bdy_conditions = NULL;

static  cs_lagr_param_t  cs_glob_lagr_param; // Should move to cs_lagr.c

static cs_lagr_attribute_map_t  *_p_attr_map = NULL;
static cs_lagr_attribute_map_t  *_p_attr_map_prev = NULL;

enum {X, Y, Z};  /* Used for _get_norm() and _get_dot_prod() */

/* MPI datatype associated to each particle structures */

#if defined(HAVE_MPI)
static  MPI_Datatype  _CS_MPI_PARTICLE;
#endif

/* Indexes in Fortran arrays */

static int _jisor = - 1, _jrval = - 1, _jrpoi = -1, _jrtsp = - 1;
static int _jmp = -1, _jdp = -1, _jxp = -1, _jyp = -1, _jzp = -1;
static int _jup = -1, _jvp = -1, _jwp = -1, _juf = -1, _jvf = -1, _jwf = -1;
static int _jtaux = - 1, _jdepo = - 1, _jrank_flag = -1, _jryplu = - 1, _jrinpf = - 1;
static int _jdfac = - 1, _jimark = -1, _jnbasg = - 1, _jnbasp = - 1;
static int _jfadh = - 1, _jmfadh = - 1, _jndisp = - 1;
static int _jtf = - 1, _jmwat = - 1;
static int _jcp = -1, _jrdck = -1, _jrd0p = - 1;
static int _jinch = - 1, _jreps = -1, _jclst = -1;
static int *_jthp = NULL, *_jmch = NULL, *_jmck = NULL, *_jrhock = NULL;
static int *_jvls = NULL;

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
 * returns:
 *   pointer to structure mapping particle attributes
 *----------------------------------------------------------------------------*/

static cs_lagr_attribute_map_t *
_create_attr_map(void)
{
  cs_lagr_attribute_t attr;

  cs_lagr_attribute_map_t  *am;

  BFT_MALLOC(am, 1, cs_lagr_attribute_map_t);

  am->n_attributes = CS_LAGR_N_ATTRIBUTES;
  am->extents = sizeof(cs_lagr_particle_t);

  BFT_MALLOC(am->size, am->n_attributes, size_t);
  BFT_MALLOC(am->displ, am->n_attributes, ptrdiff_t);
  BFT_MALLOC(am->datatype, am->n_attributes, cs_datatype_t);
  BFT_MALLOC(am->count, am->n_attributes, int);
  BFT_MALLOC(am->have_previous, am->n_attributes, bool);

  for (attr = 0; attr < (cs_lagr_attribute_t)am->n_attributes; attr++) {
    am->size[attr] = 0;
    am->displ[attr] = -1;
    am->datatype[attr] = CS_REAL_TYPE;
    am->count[attr] = 1;
    am->have_previous[attr] = true;
  }

  attr = CS_LAGR_CUR_CELL_NUM;
  am->displ[attr] = offsetof(cs_lagr_particle_t, cur_cell_num);
  am->datatype[attr] = CS_LNUM_TYPE;

  attr = CS_LAGR_LAST_FACE_NUM;
  am->displ[attr] = offsetof(cs_lagr_particle_t, last_face_num);
  am->datatype[attr] = CS_LNUM_TYPE;

  attr = CS_LAGR_SWITCH_ORDER_1;
  am->displ[attr] = offsetof(cs_lagr_particle_t, switch_order_1);
  am->datatype[attr] = CS_INT_TYPE;

  attr = CS_LAGR_STATE;
  am->displ[attr] = offsetof(cs_lagr_particle_t, state);
  am->datatype[attr] = CS_LNUM_TYPE;

  attr = CS_LAGR_RANDOM_VALUE;
  am->displ[attr] = offsetof(cs_lagr_particle_t, random_value);

  attr = CS_LAGR_STAT_WEIGHT;
  am->displ[attr] = offsetof(cs_lagr_particle_t, stat_weight);

  attr = CS_LAGR_RESIDENCE_TIME;
  am->displ[attr] = offsetof(cs_lagr_particle_t, residence_time);

  attr = CS_LAGR_MASS;
  am->displ[attr] = offsetof(cs_lagr_particle_t, mass);

  attr = CS_LAGR_DIAMETER;
  am->displ[attr] = offsetof(cs_lagr_particle_t, diameter);

  attr = CS_LAGR_TAUP_AUX;
  am->displ[attr] = offsetof(cs_lagr_particle_t, taup_aux);

  attr = CS_LAGR_COORDS;
  am->displ[attr] = offsetof(cs_lagr_particle_t, coord);
  am->count[attr] = 3;

  attr = CS_LAGR_VELOCITY;
  am->displ[attr] = offsetof(cs_lagr_particle_t, velocity);
  am->count[attr] = 3;

  attr = CS_LAGR_VELOCITY_SEEN;
  am->displ[attr] = offsetof(cs_lagr_particle_t, velocity_seen);
  am->count[attr] = 3;

  /* Deposition submodel additional parameters */

  attr = CS_LAGR_YPLUS;
  if (_jryplu > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, yplus);

  attr = CS_LAGR_INTERF;
  if (_jrinpf > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, interf);

  attr = CS_LAGR_NEIGHBOR_FACE_ID;
  if (_jdfac > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, close_face_id);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  attr = CS_LAGR_MARKO_VALUE;
  if (_jimark > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, marko_val);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  attr = CS_LAGR_DEPOSITION_FLAG;
  if (_jdepo > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, depo);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  attr = CS_LAGR_RANK_FLAG;
  if (_jrank_flag > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, rank_flag);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  /* Resuspension model additional parameters */

  attr = CS_LAGR_N_LARGE_ASPERITIES;
  if (_jnbasg > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, nb_large_asperities);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  attr = CS_LAGR_N_SMALL_ASPERITIES;
  if (_jnbasp > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, nb_small_asperities);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  attr = CS_LAGR_ADHESION_FORCE;
  if (_jfadh > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, adhesion_force);

  attr = CS_LAGR_ADHESION_TORQUE;
  if (_jmfadh > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, adhesion_torque);

  attr = CS_LAGR_DISPLACEMENT_NORM;
  if (_jndisp > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, displacement_norm);

  /* Thermal model additional parameters */

  attr = CS_LAGR_TEMPERATURE;
  if (_jthp[0] > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, temp);
    am->count[attr] = cs_glob_lagr_param.n_temperature_layers;
  }

  attr = CS_LAGR_FLUID_TEMPERATURE;
  if (_jtf > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, fluid_temp);

  attr = CS_LAGR_CP;
  if (_jcp > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, cp);
  }

  /* Coal combustion additional parameters */

  attr = CS_LAGR_WATER_MASS;
  if (_jmwat > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, water_mass);

  attr = CS_LAGR_COAL_MASS;
  if (_jmch[0] > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, coal_mass);
    am->count[attr] = CS_LAGR_N_LAYERS;
  }

  attr = CS_LAGR_COKE_MASS;
  if (_jmck[0] > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, coke_mass);
    am->count[attr] = CS_LAGR_N_LAYERS;
  }

  attr = CS_LAGR_SHRINKING_DIAMETER;
  if (_jrdck > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, shrinking_diam);

  attr = CS_LAGR_INITIAL_DIAMETER;
  if (_jrd0p > -1)
    am->displ[attr] = offsetof(cs_lagr_particle_t, initial_diam);

  attr = CS_LAGR_COAL_NUM;
  if (_jinch > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, coal_number);
    am->datatype[attr] = CS_LNUM_TYPE;
  }

  attr = CS_LAGR_COAL_DENSITY;
  if (_jrhock[0] > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, coal_density);
    am->count[attr] = CS_LAGR_N_LAYERS;
  }

  attr = CS_LAGR_EMISSIVITY;
  if (_jreps > -1) {
    am->displ[attr] = offsetof(cs_lagr_particle_t, emissivity);
  }

  attr = CS_LAGR_STAT_CLASS;
  if (_jclst > -1) {
    am->displ[attr] = am->extents;
    am->datatype[attr] = CS_LNUM_TYPE;
    am->count[attr] = 1;
    am->extents += _align_extents(sizeof(cs_lnum_t));
  }

  attr = CS_LAGR_USER;
  if (cs_glob_lagr_param.n_user_variables > 0) {
    am->displ[attr] = am->extents;
    am->count[attr] = cs_glob_lagr_param.n_user_variables;
    am->extents += _align_extents(sizeof(cs_real_t)*am->count[attr]);
  }

  for (attr = 0; attr < (cs_lagr_attribute_t)am->n_attributes; attr++) {
    if (am->displ[attr] < 0) {
      am->datatype[attr] = CS_DATATYPE_NULL;
      am->count[attr] = 0;
      am->have_previous[attr] = false;
    }
    else
      am->size[attr] = am->count[attr] * cs_datatype_size[am->datatype[attr]];
  }

  return am;
}

/*----------------------------------------------------------------------------*
 * Free particle attributes for a given configuration.
 *----------------------------------------------------------------------------*/

static void
_destroy_attr_map(cs_lagr_attribute_map_t  **am)
{
  if (*am != NULL) {
    cs_lagr_attribute_map_t  *_am = *am;

    BFT_FREE(_am->size);
    BFT_FREE(_am->displ);
    BFT_FREE(_am->datatype);
    BFT_FREE(_am->count);
    BFT_FREE(_am->have_previous);

    BFT_FREE(*am);
  }
}

/*----------------------------------------------------------------------------
 * Compute norm squared of 3D vector difference
 *
 * parameters:
 *   v1 <-- first vector
 *   v2 <-- second vector
 *
 * return:
 *   (v1-v2).(v1-v2)
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_delta_norm_2_3d(const cs_real_t  v1[],
                 const cs_real_t  v2[])
{
  cs_real_t d[3] = {v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]};
  return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

/*----------------------------------------------------------------------------
 * Compute the norm of a 3D vector of double (cs_real_t)
 *
 * parameters:
 *   vector <-- 3D vector to treat
 *
 * returns:
 *   norm associated to the vector
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_get_norm(const cs_real_t  vect[])
{
  return sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z]);
}

/*----------------------------------------------------------------------------
 * Compute the norm of a 3D vector of double (cs_real_t)
 *
 * parameters:
 *   v1 <-- 3D vector to treat
 *   v2 <-- 3D vector to treat
 *
 * returns:
 *   norm associated to the vector
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_get_dot_prod(const cs_real_t  v1[],
              const cs_real_t  v2[])
{
  return (v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z]);
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a given
 * vector of coordinates.
 *
 * parameters:
 *   matrix[3][4] <-- matrix of the transformation in homogeneous coord
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   v            <-> vector
 *----------------------------------------------------------------------------*/

inline static void
_apply_vector_transfo(const cs_real_t  matrix[3][4],
                      cs_real_t        v[])
{
  cs_lnum_t  i, j;

  /* Define a vector in homogeneous coordinates before transformation */

  cs_real_t  t[4] = {v[0], v[1], v[2], 1};

  /* Initialize output */

  for (i = 0; i < 3; i++)
    v[i] = 0.;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++)
      v[i] += matrix[i][j]*t[j];
  }
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a rotation to a given vector.
 *
 * parameters:
 *   matrix[3][4] <-- matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   v            <-> vector to rotate
 *----------------------------------------------------------------------------*/

inline static void
_apply_vector_rotation(const cs_real_t   matrix[3][4],
                       cs_real_t         v[3])
{
  const cs_real_t v_in[3] = {v[0], v[1], v[2]};

  v[0] = matrix[0][0]*v_in[0] + matrix[0][1]*v_in[1] + matrix[0][2]*v_in[2];
  v[1] = matrix[1][0]*v_in[0] + matrix[1][1]*v_in[1] + matrix[1][2]*v_in[2];
  v[2] = matrix[2][0]*v_in[0] + matrix[2][1]*v_in[1] + matrix[2][2]*v_in[2];
}

/*----------------------------------------------------------------------------
 * Remove a particle from a particle set.
 *
 * parameters:
 *   set    <->  particles set
 *   cur_id <-- id of particle to remove
 *----------------------------------------------------------------------------*/

static void
_remove_particle(cs_lagr_particle_set_t   *set,
                 cs_lnum_t                 cur_id)
{
  /* Remove cur_part from "used particle list" */

  cs_lnum_t prev_id = set->used_id[cur_id].prev_id;
  cs_lnum_t next_id = set->used_id[cur_id].next_id;

  if (prev_id != -1)
    set->used_id[prev_id].next_id = next_id;
  else
    set->first_used_id = next_id;

  if (next_id != set->n_particles_max && next_id != -1)
    set->used_id[next_id].prev_id = prev_id;

  /* Add cur_part to "free particle list" */

  if (cur_id < set->first_free_id) {

    cs_lnum_t old_first_free = set->first_free_id;
    cs_lnum_t old_first_free_prev = set->used_id[old_first_free].prev_id;

    set->first_free_id = cur_id;

    set->used_id[set->first_free_id].next_id = old_first_free;
    set->used_id[set->first_free_id].prev_id = old_first_free_prev;

    set->used_id[old_first_free].prev_id = cur_id;

  }
  else { /* We place the cur_part just behind the first free particle. */

    cs_lnum_t first_free = set->first_free_id;
    cs_lnum_t old_next = set->used_id[first_free].next_id;

    set->used_id[first_free].next_id = cur_id;

    set->used_id[cur_id].next_id = old_next;
    set->used_id[cur_id].prev_id = first_free;

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype which maps main particle characteristics.
 *
 * parameters:
 *   am_cur  <-- attributes map for current set
 *   am_prev <-- attributes map for previous set
 *
 * returns:
 *   MPI_Datatype matching given attributes map
 *----------------------------------------------------------------------------*/

static MPI_Datatype
_define_particle_datatype(const cs_lagr_attribute_map_t  *am_cur,
                          const cs_lagr_attribute_map_t  *am_prev)
{
  size_t i;
  MPI_Datatype  new_type;
  int           count;
  cs_datatype_t *cs_type;
  int           *blocklengths;
  MPI_Datatype  *types;
  MPI_Aint      *displacements;

  size_t tot_extents = am_cur->extents + am_prev->extents;
  const cs_lagr_attribute_map_t *a_maps[2] = {am_cur, am_prev};

  /* Mark bytes with associated type */

  BFT_MALLOC(cs_type, tot_extents, cs_datatype_t);

  for (i = 0; i < tot_extents; i++)
    cs_type[i] = CS_CHAR;

  size_t shift = 0;

  for (int j = 0; j < 2; j++) {

    const cs_lagr_attribute_map_t *am = a_maps[j];

    for (size_t attr = 0; attr < am->n_attributes; attr++) {
      if (am->count[attr] > 0) {
        assert(am->displ[attr] > -1);
        size_t b_size = am->count[attr] * cs_datatype_size[am->datatype[attr]];
        for (i = 0; i < b_size; i++)
          cs_type[am->displ[attr] + i + shift] = am->datatype[attr];
      }
    }

    shift += am->extents;

  }

  /* Count type groups */

  count = 0;

  i = 0;
  while (i < tot_extents) {
    size_t j;
    for (j = i; j < tot_extents; j++) {
      if (cs_type[j] != cs_type[i])
        break;
    }
    count += 1;
    i = j;
  }

  /* Assign types */

  BFT_MALLOC(blocklengths, count, int);
  BFT_MALLOC(types, count, MPI_Datatype);
  BFT_MALLOC(displacements, count, MPI_Aint);

  count = 0;

  i = 0;
  while (i < tot_extents) {
    size_t j;
    types[count] = cs_datatype_to_mpi[cs_type[i]];
    displacements[count] = i;
    for (j = i; j < tot_extents; j++) {
      if (cs_type[j] != cs_type[i])
        break;
    }
    blocklengths[count] = (j-i) / cs_datatype_size[cs_type[i]];
    count += 1;
    i = j;
  }

  /* Create new datatype */

  MPI_Type_create_struct(count, blocklengths, displacements, types,
                         &new_type);

  MPI_Type_commit(&new_type);

  BFT_FREE(displacements);
  BFT_FREE(types);
  BFT_FREE(blocklengths);
  BFT_FREE(cs_type);

  return new_type;
}

/*----------------------------------------------------------------------------
 * Delete all the MPI_Datatypes related to particles.
 *----------------------------------------------------------------------------*/

static void
_delete_particle_datatypes(void)
{
  MPI_Type_free(&_CS_MPI_PARTICLE);
}
#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Allocate a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *   time_id         <-- 0 for current time, 1 for previous
 *   n_particles_max <-- local max. number of particles
 *   p_am            <-- particle attributes map
 *
 * returns:
 *   a new allocated cs_lagr_particle_set_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_particle_set_t *
_create_particle_set(int                             time_id,
                     cs_lnum_t                       n_particles_max,
                     const cs_lagr_attribute_map_t  *p_am)

{
  cs_lnum_t  i;

  cs_lagr_particle_set_t  *new_set = NULL;

  if (n_particles_max == 0)
    return NULL;

  BFT_MALLOC(new_set, 1, cs_lagr_particle_set_t);

  new_set->time_id = time_id;

  BFT_MALLOC(new_set->p_buffer, n_particles_max * p_am->extents, unsigned char);

  new_set->n_particles_max = n_particles_max;
  new_set->n_particles = 0;

  if (time_id == 0) {
    BFT_MALLOC(new_set->used_id, n_particles_max, cs_lagr_tracking_list_t);
    new_set->first_used_id = -1;
    new_set->first_free_id = 0;
    for (i = 0; i < n_particles_max; i++) {
      new_set->used_id[i].prev_id = i-1;
      new_set->used_id[i].next_id = i+1;
    }
  }
  else {
    new_set->used_id = NULL;
    new_set->first_used_id = -1;
    new_set->first_free_id = -1;
  }

  assert(n_particles_max >= 1);

  new_set->p_am = p_am;

  return new_set;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *   particle_set    <->  a cs_lagr_particle_set_t structure
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

static cs_lagr_particle_set_t *
_destroy_particle_set(cs_lagr_particle_set_t *set)
{
  if (set == NULL)
    return set;

  BFT_FREE(set->used_id);
  BFT_FREE(set->p_buffer);

  BFT_FREE(set);

  return NULL;
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

  for (cs_lagr_attribute_t attr = 0;
       attr < (cs_lagr_attribute_t)am->n_attributes;
       attr++) {
    if (am->count[attr] > 0) {
      char attr_name[64];
      strncpy(attr_name, cs_lagr_attribute_name[attr] + 8, 63);
      attr_name[63] = '\0';
      for (int i = 0; attr_name[i] != '\0'; i++)
        attr_name[i] = tolower(attr_name[i]);
      switch (am->datatype[attr]) {
      case CS_LNUM_TYPE:
        {
          const cs_lnum_t *v
            = cs_lagr_particle_attr_const(p, particles->p_am, attr);
          bft_printf("    %24s: %10ld\n", attr_name, (long)v[0]);
          for (int i = 1; i < am->count[attr]; i++)
            bft_printf("    %24s: %10ld\n", " ", (long)v[i]);
        }
        break;
      case CS_GNUM_TYPE:
        {
          const cs_gnum_t *v
            = cs_lagr_particle_attr_const(p, particles->p_am, attr);
          bft_printf("    %24s: %10lu\n", attr_name, (unsigned long)v[0]);
          for (int i = 1; i < am->count[attr]; i++)
            bft_printf("    %24s: %10lu\n", " ", (unsigned long)v[i]);
        }
        break;
      case CS_REAL_TYPE:
        {
          const cs_real_t *v
            = cs_lagr_particle_attr_const(p, particles->p_am, attr);
          bft_printf("    %24s: %10.3g\n", attr_name, v[0]);
          for (int i = 1; i < am->count[attr]; i++)
            bft_printf("    %24s: %10.3g\n", " ", v[i]);
        }
        break;
      default:
        break;
      }
    }
  }
  bft_printf("\n");
}

#if 0 /* TODO will be used later, when C and Fortran versions are mapped */

/*----------------------------------------------------------------------------
 * Resize a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *   p_particle_set    <-> pointer to a cs_lagr_particle_set_t structure
 *   p_am              <-- particle attributes map
 *   time_id           <-- 0 for current time, 1 for previous
 *   n_particles_max   <-- local max. number of particles
 *----------------------------------------------------------------------------*/

static void
_resize_particle_set(cs_lagr_particle_set_t        **p_particle_set,
                     const cs_lagr_attribute_map_t  *p_am,
                     int                             time_id,
                     const cs_lnum_t                 n_particles_max)
{
  cs_lagr_particle_set_t  *particle_set = *p_particle_set;

  assert(n_particles_max >= 0);

  if (n_particles_max == 0)
    particle_set = _destroy_particle_set(particle_set);

  else if (particle_set == NULL && n_particles_max > 0)
    particle_set = _create_particle_set(time_id, n_particles_max, p_am);

  else if (particle_set->n_particles_max < n_particles_max) {

    particle_set->n_particles_max = n_particles_max;

    BFT_REALLOC(particle_set->particles, n_particles_max, cs_lagr_particle_t);
    particle_set->p_buffer = (unsigned char* )particle_set->particles;

    if (particle_set->time_id == 0)
      BFT_REALLOC(particle_set->used_id, n_particles_max, cs_lagr_tracking_list_t);

  }
  else {

    /* FIX ME */
    /* bft_printf("nmax local = %d\n",particle_set->n_particles_max); */
    /* bft_printf("n demande : %d\n",n_particles_max); */

    /* bft_error(__FILE__, __LINE__, 0, */
    /*           _(" The current situation is not managed.\n")); */

  }

  /* Returns pointer */

  if (particle_set != NULL) {
    _particle_set->p_size = p_am->extents;
    _particle_set->p_am = p_am;
  }

  *p_particle_set = particle_set;
}

#endif

/*----------------------------------------------------------------------------
 * Create a cs_lagr_halo_t structure to deal with parallelism and
 * periodicity
 *
 * parameters:
 *   cur_extents  <-- extents for particles of current set
 *   prev_extents <-- extents for particles of previous set
 *
 * returns:
 *   a new allocated cs_lagr_halo_t structure.
 *----------------------------------------------------------------------------*/

static cs_lagr_halo_t *
_create_lagr_halo(size_t  cur_extents,
                  size_t  prev_extents)
{
  cs_lnum_t  i, rank, tr_id, shift, start, end, n;

  cs_lnum_t  halo_cell_id = 0;
  cs_lnum_t  *cell_num = NULL;
  cs_lagr_halo_t  *lagr_halo = NULL;

  const size_t tot_extents = cur_extents + prev_extents;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const cs_lnum_t  n_halo_cells = halo->n_elts[CS_HALO_EXTENDED];

  BFT_MALLOC(lagr_halo, 1, cs_lagr_halo_t);

  assert(n_halo_cells == halo->index[2*halo->n_c_domains]);
  assert(n_halo_cells == mesh->n_ghost_cells);

  lagr_halo->n_cells = n_halo_cells;

  /* Allocate buffers to enable the exchange between communicating ranks */

  BFT_MALLOC(lagr_halo->send_shift, halo->n_c_domains, cs_lnum_t);
  BFT_MALLOC(lagr_halo->send_count, halo->n_c_domains, cs_lnum_t);
  BFT_MALLOC(lagr_halo->recv_shift, halo->n_c_domains, cs_lnum_t);
  BFT_MALLOC(lagr_halo->recv_count, halo->n_c_domains, cs_lnum_t);

  lagr_halo->send_buf_size = CS_LAGR_MIN_COMM_BUF_SIZE;
  lagr_halo->recv_buf_size = CS_LAGR_MIN_COMM_BUF_SIZE;

  lagr_halo->cur_extents = cur_extents;
  lagr_halo->prev_extents = prev_extents;

  BFT_MALLOC(lagr_halo->send_buf,
             lagr_halo->send_buf_size * tot_extents,
             unsigned char);
  BFT_MALLOC(lagr_halo->recv_buf,
             lagr_halo->recv_buf_size * tot_extents,
             unsigned char);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    cs_lnum_t  request_size = 2 * halo->n_c_domains;

    BFT_MALLOC(lagr_halo->request, request_size, MPI_Request);
    BFT_MALLOC(lagr_halo->status,  request_size, MPI_Status);

  }
#endif

  /* Fill rank */

  BFT_MALLOC(lagr_halo->rank, n_halo_cells, cs_lnum_t);

  for (rank = 0; rank < halo->n_c_domains; rank++) {

    for (i = halo->index[2*rank]; i < halo->index[2*rank+2]; i++)
      lagr_halo->rank[halo_cell_id++] = rank;

  }

  assert(halo_cell_id == n_halo_cells);

  /* Fill transform_id */

  BFT_MALLOC(lagr_halo->transform_id, n_halo_cells, cs_lnum_t);

  for (i = 0; i < n_halo_cells; i++)
    lagr_halo->transform_id[i] = -1; /* Undefined transformation */

  if (mesh->n_init_perio > 0) { /* Periodicity */

    for (tr_id = 0; tr_id < mesh->n_transforms; tr_id++) {

      shift = 4 * halo->n_c_domains * tr_id;

      for (rank = 0; rank < halo->n_c_domains; rank++) {

        /* standard */
        start = halo->perio_lst[shift + 4*rank];
        n =  halo->perio_lst[shift + 4*rank + 1];
        end = start + n;

        for (i = start; i < end; i++)
          lagr_halo->transform_id[i] = tr_id;

        /* extended */
        start = halo->perio_lst[shift + 4*rank + 2];
        n =  halo->perio_lst[shift + 4*rank + 3];
        end = start + n;

        for (i = start; i < end; i++)
          lagr_halo->transform_id[i] = tr_id;

      }

    } /* End of loop on transforms */

  } /* End of periodicity handling */

  /* Fill dist_cell_num */

  BFT_MALLOC(lagr_halo->dist_cell_num, n_halo_cells, cs_lnum_t);

  BFT_MALLOC(cell_num, mesh->n_cells_with_ghosts, cs_lnum_t);

  for (i = 0; i < mesh->n_cells_with_ghosts; i++)
    cell_num[i] = i+1;

  cs_halo_sync_num(halo, CS_HALO_EXTENDED, cell_num);

  for (i = 0; i < n_halo_cells; i++)
    lagr_halo->dist_cell_num[i] = cell_num[mesh->n_cells + i];

  /* Free memory */

  BFT_FREE(cell_num);

  return lagr_halo;
}

/*----------------------------------------------------------------------------
 * Delete a cs_lagr_halo_t structure.
 *
 * parameters:
 *   halo <-- pointer to cs_lagr_halo_t structure to delete
 *----------------------------------------------------------------------------*/

static void
_delete_lagr_halo(cs_lagr_halo_t   **halo)
{
  if (*halo != NULL) {

    cs_lagr_halo_t *h = *halo;

    BFT_FREE(h->rank);
    BFT_FREE(h->transform_id);
    BFT_FREE(h->dist_cell_num);

    BFT_FREE(h->send_shift);
    BFT_FREE(h->send_count);
    BFT_FREE(h->recv_shift);
    BFT_FREE(h->recv_count);

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {
      BFT_FREE(h->request);
      BFT_FREE(h->status);
    }
#endif

    BFT_FREE(h->send_buf);
    BFT_FREE(h->recv_buf);

    BFT_FREE(*halo);
  }

}

/*----------------------------------------------------------------------------
 * Resize a halo's buffers.
 *
 * parameters:
 *   lag_halo         <->  pointer to a cs_lagr_halo_t structure
 *   n_send_particles <-- number of particles to send
 *   n_recv_particles <-- number of particles to receive
 *----------------------------------------------------------------------------*/

static void
_resize_lagr_halo(cs_lagr_halo_t  *lag_halo,
                  cs_lnum_t        n_send_particles,
                  cs_lnum_t        n_recv_particles)
{
  size_t n_target[] = {n_send_particles, n_recv_particles};
  size_t *n_halo[] = {&(lag_halo->send_buf_size),
                      &(lag_halo->recv_buf_size)};
  unsigned char **buf[] = {&(lag_halo->send_buf),
                           &(lag_halo->recv_buf)};

  size_t tot_extents = lag_halo->cur_extents + lag_halo->prev_extents;

  for (int i = 0; i < 2; i++) {

    /* If increase is required */

    if (*(n_halo[i]) < n_target[i]) {
      if (*(n_halo[i]) < CS_LAGR_MIN_COMM_BUF_SIZE)
        *(n_halo[i]) = CS_LAGR_MIN_COMM_BUF_SIZE;
      while (*(n_halo[i]) < n_target[i])
        *(n_halo[i]) *= 2;
      BFT_REALLOC(*(buf[i]), *(n_halo[i])*tot_extents, unsigned char);
    }

    /* If decrease is allowed, do it progressively, and with a wide
       margin, so as to avoid re-increasing later if possible */

    else if (*(n_halo[i]) > n_target[i]*16) {
      *(n_halo[i]) /= 8;
      BFT_REALLOC(*(buf[i]), *(n_halo[i])*tot_extents, unsigned char);
    }

    /* Otherwise, keep current size */

  }
}

/*----------------------------------------------------------------------------
 * Define a cell -> face connectivity. Index begins with 0.
 *
 * parameters:
 *   builder   <--  pointer to a cs_lagr_track_builder_t structure
 *----------------------------------------------------------------------------*/

static void
_define_cell_face_connect(cs_lagr_track_builder_t   *builder)
{
  cs_lnum_t  i, j;

  cs_lnum_t  *counter = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  BFT_MALLOC(counter, mesh->n_cells, cs_lnum_t);
  BFT_MALLOC(builder->cell_face_idx, mesh->n_cells + 1, cs_lnum_t);

  /* Initialize */

  builder->cell_face_idx[0] = 0;
  for (i = 0; i < mesh->n_cells; i++) {
    builder->cell_face_idx[i+1] = 0;
    counter[i] = 0;
  }

  /* Count of the number of faces per cell: loop on interior faces */

  for (i = 0; i < mesh->n_i_faces; i++)
    for (j = 0; j < 2; j++) {
      cs_lnum_t iel = mesh->i_face_cells[2*i+j];
      if (iel <= mesh->n_cells)
                builder->cell_face_idx[iel] += 1;
    }

  /* Count of the number of faces per cell: loop on border faces */

  for (i = 0; i < mesh->n_b_faces; i++)
    builder->cell_face_idx[mesh->b_face_cells[i]] += 1;

  /* Build index */

  for (i = 0; i < mesh->n_cells; i++)
    builder->cell_face_idx[i+1] += builder->cell_face_idx[i];

  BFT_MALLOC(builder->cell_face_lst,
             builder->cell_face_idx[mesh->n_cells], cs_lnum_t);

  /* Build list: border faces are < 0 and interior faces > 0 */

  for (i = 0; i < mesh->n_i_faces; i++) {
    for (j = 0; j < 2; j++) {

      cs_lnum_t iel = mesh->i_face_cells[2*i+j];

      if (iel <= mesh->n_cells) {

        cs_lnum_t  cell_id = iel - 1;
        cs_lnum_t  shift = builder->cell_face_idx[cell_id] + counter[cell_id];

        builder->cell_face_lst[shift] = i+1;
        counter[cell_id] += 1;
      }
    }
  }

  for (i = 0; i < mesh->n_b_faces; i++) {

    cs_lnum_t  cell_id = mesh->b_face_cells[i] - 1;
    cs_lnum_t  shift = builder->cell_face_idx[cell_id] + counter[cell_id];

    builder->cell_face_lst[shift] = -(i+1);
    counter[cell_id] += 1;

  }

  /* Free memory */

  BFT_FREE(counter);
}

/*----------------------------------------------------------------------------
 * Initialize a cs_lagr_track_builder_t structure.
 *
 * parameters:
 *   n_particles_max <-- local max number of particles
 *   cur_extents     <-- extents for particles of current set
 *   prev_extents    <-- extents for particles of previous set
 *
 * returns:
 *   a new defined cs_lagr_track_builder_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_track_builder_t *
_init_track_builder(cs_lnum_t  n_particles_max,
                    size_t     cur_extents,
                    size_t     prev_extents)
{
  cs_lnum_t  i;
  cs_mesh_t  *mesh = cs_glob_mesh;

  cs_lagr_track_builder_t  *builder = NULL;

  if (n_particles_max == 0)
    return NULL;

  BFT_MALLOC(builder, 1, cs_lagr_track_builder_t);

  /* Define _max_face_connect_size and _face_connect_buffer */

  builder->max_face_connect_size = 0;

  /* Loop on interior faces */

  for (i = 0; i < mesh->n_i_faces; i++)
    builder->max_face_connect_size =
      CS_MAX(builder->max_face_connect_size,
             mesh->i_face_vtx_idx[i+1] - mesh->i_face_vtx_idx[i]);

  /* Loop on border faces */

  for (i = 0; i < mesh->n_b_faces; i++)
    builder->max_face_connect_size =
      CS_MAX(builder->max_face_connect_size,
             mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i]);

  builder->max_face_connect_size += 1;

  BFT_MALLOC(builder->face_connect_buffer,
             builder->max_face_connect_size,
             cs_lnum_t);

  for (i = 0; i < builder->max_face_connect_size; i++)
    builder->face_connect_buffer[i] = -1;

  /* Define a cell->face connectivity */

  _define_cell_face_connect(builder);

  /* Define a cs_lagr_halo_t structure to deal with parallelism and
     periodicity */

  if (cs_glob_mesh->n_init_perio > 0 || cs_glob_n_ranks > 1)
    builder->halo = _create_lagr_halo(cur_extents, prev_extents);
  else
    builder->halo = NULL;

  /* Define an interface set on interior faces for keeping up-to-date
     the last_face_num value across ranks. Not used in serial mode */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    builder->face_ifs = cs_interface_set_create(mesh->n_i_faces,
                                                NULL,
                                                mesh->global_i_face_num,
                                                NULL,
                                                0,
                                                NULL,
                                                NULL,
                                                NULL);

    cs_interface_set_add_match_ids(builder->face_ifs);
  }

  else
    builder->face_ifs = NULL;
#endif

  return builder;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_lagr_track_builder_t structure.
 *
 * parameters:
 *   builder   <--  pointer to a cs_lagr_track_builder_t structure
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

static cs_lagr_track_builder_t *
_destroy_track_builder(cs_lagr_track_builder_t  *builder)
{
  if (builder == NULL)
    return builder;

  BFT_FREE(builder->face_connect_buffer);
  BFT_FREE(builder->cell_face_idx);
  BFT_FREE(builder->cell_face_lst);

  /* Destroy the cs_lagr_halo_t structure */

  _delete_lagr_halo(&(builder->halo));
  cs_interface_set_destroy(&(builder->face_ifs));

  /* Destroy the builder structure */

  BFT_FREE(builder);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Initialize a cs_lagr_bdy_condition_t structure.
 *
 * parameters:
 *   n_max_zones     <--  number max. of boundary zones
 *
 * returns:
 *   a new defined cs_lagr_bdy_condition_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_bdy_condition_t *
_create_bdy_cond_struct(cs_lnum_t  n_max_zones)
{
  cs_lnum_t  i;

  cs_lagr_bdy_condition_t *bdy_cond = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  BFT_MALLOC(bdy_cond, 1, cs_lagr_bdy_condition_t);

  bdy_cond->n_b_zones = 0;
  bdy_cond->n_b_max_zones = n_max_zones;

  BFT_MALLOC(bdy_cond->particle_flow_rate, n_max_zones, cs_real_t);
  BFT_MALLOC(bdy_cond->b_zone_lst, n_max_zones, cs_lnum_t);
  BFT_MALLOC(bdy_cond->b_zone_classes, n_max_zones, cs_lnum_t);
  BFT_MALLOC(bdy_cond->b_zone_natures, n_max_zones, cs_lnum_t);

  for (i = 0; i < n_max_zones; i++) {

    bdy_cond->particle_flow_rate[i] = 0.0;
    bdy_cond->b_zone_lst[i] = -1;
    bdy_cond->b_zone_classes[i] = -1;
    bdy_cond->b_zone_natures[i] = -1;

  }

  BFT_MALLOC(bdy_cond->b_face_zone_num, mesh->n_b_faces, cs_lnum_t);

  for (i = 0; i < cs_glob_mesh->n_b_faces; i++)
    bdy_cond->b_face_zone_num[i] = -1;

  bdy_cond->continuous_injection = 0.0;
  bdy_cond->steady_bndy_conditions = false;

  return bdy_cond;
}

/*----------------------------------------------------------------------------
 * Initialize a cs_lagr_bdy_condition_t structure.
 *
 * parameters:
 *   n_max_zones     <--  number max. of boundary zones
 *
 * returns:
 *   a new defined cs_lagr_bdy_condition_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_bdy_condition_t *
_resize_bdy_cond_struct(cs_lnum_t  n_max_zones)
{
  cs_lagr_bdy_condition_t *bdy_cond = _lagr_bdy_conditions;

  assert(bdy_cond != NULL);

  bdy_cond->n_b_zones = n_max_zones;
  bdy_cond->n_b_max_zones = n_max_zones;

  BFT_REALLOC(bdy_cond->particle_flow_rate, bdy_cond->n_b_zones, cs_real_t);
  BFT_REALLOC(bdy_cond->b_zone_lst, bdy_cond->n_b_zones, cs_lnum_t);
  BFT_REALLOC(bdy_cond->b_zone_classes, bdy_cond->n_b_zones, cs_lnum_t);
  BFT_REALLOC(bdy_cond->b_zone_natures, bdy_cond->n_b_zones, cs_lnum_t);

  return bdy_cond;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_lagr_bdy_condition_t structure.
 *
 * returns:
 *   a NULL pointer.
 *----------------------------------------------------------------------------*/

static cs_lagr_bdy_condition_t *
_destroy_bdy_cond_struct(cs_lagr_bdy_condition_t  *bdy_cond)
{
  if (bdy_cond != NULL) {

    BFT_FREE(bdy_cond->b_zone_lst);
    BFT_FREE(bdy_cond->b_zone_natures);
    BFT_FREE(bdy_cond->b_zone_classes);

    BFT_FREE(bdy_cond->b_face_zone_num);

    BFT_FREE(bdy_cond->particle_flow_rate);

    BFT_FREE(bdy_cond);

  }

  return NULL;
}

/*----------------------------------------------------------------------------
 * Manage detected errors
 *
 * parameters:
 *   failsafe_mode            <-- indicate if failsafe mode is used
 *   particle                 <-- pointer to particle data
 *   attr_map                 <-- pointer to attribute map
 *   error_type               <-- error code
 *   p_n_failed_particles     <->
 *   p_failed_particle_weight <->
 *   msg                      <-- error message
 *----------------------------------------------------------------------------*/

static void
_manage_error(cs_lnum_t                       failsafe_mode,
              void                           *particle,
              const cs_lagr_attribute_map_t  *attr_map,
              cs_lagr_tracking_error_t        error_type,
              cs_lnum_t                      *p_n_failed_particles,
              cs_real_t                      *p_failed_particle_weight)
{
  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight = *p_failed_particle_weight;

  cs_lagr_particle_set_lnum(particle, attr_map, CS_LAGR_CUR_CELL_NUM, 0);

  n_failed_particles++;
  failed_particle_weight
    += cs_lagr_particle_get_real(particle, attr_map, CS_LAGR_STAT_WEIGHT);

  if (failsafe_mode == 1) {
    switch (error_type) {
    case CS_LAGR_TRACKING_ERR_INSIDE_B_FACES:
      bft_error(__FILE__, __LINE__, 0,
                _(" Error during boundary treatment.\n"
                  " Part of trajectography inside the boundary faces."));
      break;
    case CS_LAGR_TRACKING_ERR_MAX_LOOPS:
      bft_error(__FILE__, __LINE__, 0,
                _("Max number of loops reached in particle displacement."));
      break;
    case CS_LAGR_TRACKING_ERR_DISPLACEMENT_I_FACE:
      bft_error(__FILE__, __LINE__, 0,
                _("Error during the particle displacement (Interior face)."));
      break;
    case CS_LAGR_TRACKING_ERR_DISPLACEMENT_B_FACE:
      bft_error(__FILE__, __LINE__, 0,
                _("Error during the particle displacement (Boundary face)."));
      break;
    default:
      break;
    }
  }

  /* Return pointers */

  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;
}

/*----------------------------------------------------------------------------
 * Test if all displacements are finished for all ranks.
 *
 * parameters:
 *   n_particles     <--  local number of particles
 *
 * returns:
 *   true if there is a need to move particles or false, otherwise
 *----------------------------------------------------------------------------*/

static bool
_continue_displacement(void)
{
  cs_lnum_t  i, j;
  cs_lnum_t  _test = 1, test = 1;

  const cs_lagr_particle_set_t  *set = _particle_set;
  const cs_lnum_t  n_particles = set->n_particles;

  for (i = 0, j = set->first_used_id; i < n_particles; i++) {
    if (   cs_lagr_particles_get_lnum(set, j, CS_LAGR_STATE)
        == CS_LAGR_PART_TO_SYNC) {
      _test = 0;
      break;
    }
    j = set->used_id[j].next_id;
  }

  if (cs_glob_n_ranks == 1)
    test = _test;

  else {

    assert(cs_glob_n_ranks > 1);

#if defined(HAVE_MPI)
    /*  MPI_Allreduce(&_test, &test, 1, CS_MPI_INT, MPI_MAX,
        cs_glob_mpi_comm); */
    MPI_Allreduce(&_test, &test, 1, CS_MPI_INT, MPI_MIN,
                  cs_glob_mpi_comm);
#endif /* HAVE_MPI */

  }

  if (test == 0)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------
 * Test if the current particle moves to the next cell through this face.
 *
 *                               |
 *      x------------------------|--------x Q: particle location
 *   P: prev. particle location  |
 *                               x Face (Center of Gravity)
 *             x current         |
 *               cell center     |
 *                               |
 *                           Face number
 *
 * parameters:
 *   face_num      <-- local number of the studied face
 *   n_vertices    <-- size of the face connectivity
 *   face_connect  <-- face -> vertex connectivity
 *   prev_particle <-- particle attributes for previous time step
 *   particle      <-- particle attributes for current time step
 *   prev_p_am     <-- pointer to attributes map for previous time step
 *   p_am          <-- pointer to attributes map for current time step
 *   p_error       <-> pointer to an error indicator
 *
 * returns:
 *   -1: if the particle remains inside the same cell
 *    0: if the trajectory doesn't go through the current face
 *    1: if the trajectory goes through the current face
 *----------------------------------------------------------------------------*/

static int
_where_are_you(cs_lnum_t                       face_num,
               cs_lnum_t                       n_vertices,
               const cs_lnum_t                 face_connect[],
               const void                     *prev_particle,
               const void                     *particle,
               const cs_lagr_attribute_map_t  *prev_p_am,
               const cs_lagr_attribute_map_t  *p_am,
               int                            *p_error)
{
  cs_lnum_t  i, j, vtx_id1, vtx_id2;
  cs_real_t  face_cog[3], cell_cen[3], vtx1[3], vtx2[3];

  cs_lnum_t  cur_cell_id
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM) - 1;
  const cs_real_t  *prev_location
    = cs_lagr_particle_attr_const(prev_particle, prev_p_am, CS_LAGR_COORDS);
  const cs_real_t  *next_location
    = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_COORDS);

  cs_mesh_t  *mesh = cs_glob_mesh;

  /* Initialize local parameters */

  cs_real_t  max_value = DBL_MIN;
  cs_lnum_t  orient_count = 0;
  cs_lnum_t  face_orient = 0;
  int        orient = 0;
  int        orient_test = 0;
  int        first_orient = 0;
  int        ijkl_ref = 0;
  int        colocalization = -1;
  int        error = 0;
  int        retcode = -999; /* initialize to an incoherent value */

  assert(sizeof(cs_real_t) == 8);

  /* Initialization */

  for (j = 0; j < 3; j++)
    cell_cen[j] = cs_glob_mesh_quantities->cell_cen[3*cur_cell_id+j];

  if (face_num > 0) { /* Interior  face */

    cs_lnum_t  face_id = face_num - 1;

    for (j = 0; j < 3; j++)
      face_cog[j] = cs_glob_mesh_quantities->i_face_cog[3*face_id+j];

  }
  else { /* Border face */

    cs_lnum_t  face_id = CS_ABS(face_num) - 1;

    for (j = 0; j < 3; j++)
      face_cog[j] = cs_glob_mesh_quantities->b_face_cog[3*face_id+j];

  }

  /* First: compute max_value to set a calculation grid */

  /* Vertex coordinates of the studied face */

  for (i = 0; i < n_vertices - 1; i++) {
    cs_lnum_t  vtx_id = face_connect[i] - 1;
    for (j = 0; j < 3; j++)
      max_value = CS_MAX(max_value, CS_ABS(mesh->vtx_coord[3*vtx_id+j]));
  }

  /* Center of the current cell */
  for (j = 0; j < 3; j++)
    max_value = CS_MAX(max_value, CS_ABS(cell_cen[j]));

  /* Center of gravity of the current face */
  for (j = 0; j < 3; j++)
    max_value = CS_MAX(max_value, CS_ABS(face_cog[j]));

  /* Starting/Ending location of the particle */
  for (j = 0; j < 3; j++) {
    max_value = CS_MAX(max_value, CS_ABS(prev_location[j]));
  }

  /* Starting/Ending location of the particle */
  for (j = 0; j < 3; j++) {
    max_value = CS_MAX(max_value, CS_ABS(next_location[j]));
  }

  /* Check if the two location are different */

  colocalization = cs_lagrang_check_colocalization(prev_location,
                                                   next_location);

  if (colocalization == 1) {
    retcode = -1;
    return retcode;
  }

  /* Check face orientation with the tetrahedron [P, */

  assert(colocalization == 0);

  vtx_id1 = face_connect[0] - 1; /* First vertex of the face */
  vtx_id2 = face_connect[1] - 1; /* Second vertex of the face */

  for (i = 0; i < 3; i++) {
    vtx1[i] = mesh->vtx_coord[3*vtx_id1+i];
    vtx2[i] = mesh->vtx_coord[3*vtx_id2+i];
  }

  face_orient = cs_lagrang_tetra_orientation(prev_location,
                                             face_cog,
                                             vtx1,
                                             vtx2);

  /* Special treatment in case of periodicity  */

  if (mesh->n_init_perio > 0)
    face_orient = 0;

  if (face_orient == 0)  /* => coplanar. Change prev_location
                            by cell center */
    face_orient = cs_lagrang_tetra_orientation(cell_cen,
                                               face_cog,
                                               vtx1,
                                               vtx2);

  if (face_orient == 0) { /* points are still coplanar */
#if 1 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf(_(" Lagrangian module warning.\n"
                 "  Failure during the particle tracking.\n"
                 "  Wrong face orientation"
                 "  => particle is lost.\n"
                 "  Local face num.:  %d\n"
                 "  Local cell num.:  %d\n"), face_num, cur_cell_id+1);
#endif
    error = 1;
    *p_error = error;

    return retcode;
  }

  /* Test first vertex of the face: [P, Q, Face_CoG, V1] */

  first_orient = cs_lagrang_tetra_orientation(prev_location,
                                              next_location,
                                              face_cog,
                                              vtx1);

  first_orient *= face_orient;

  if (first_orient == 0) { /* points are coplanar */
#if 1 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf(_(" Lagrangian module warning.\n"
                 "  Failure during the particle tracking.\n"
                 "  [P, Q, Face CoG, V(1)] orientation failed"
                 "  => Particle is lost.\n"
                 "  Local face num.:  %d\n"
                 "  Local cell num.:  %d\n"), face_num, cur_cell_id+1);
#endif
    error = 1;
    *p_error = error;

    return retcode;
  }

  /* Loop on all the vertices of the face and test orientation */

  for (i = 1; i < n_vertices; i++) {

    vtx_id1 = face_connect[i] - 1;
    for (j = 0; j < 3; j++)
      vtx1[j] = mesh->vtx_coord[3*vtx_id1+j];

    /* Test first vertex of the face: [P, Q, Face_CoG, V(i)] */

    orient = cs_lagrang_tetra_orientation(prev_location,
                                          next_location,
                                          face_cog,
                                          vtx1);

    orient *= face_orient;

    if (orient == 0) { /* points are coplanar */
#if 1 && defined(DEBUG) && !defined(NDEBUG)
      bft_printf(_(" Lagrangian module warning.\n"
                   "  Failure during the particle tracking.\n"
                   "  [P, Q, Face CoG, V(%d)] orientation failed"
                   "  => Particle is lost.\n"
                   "  Local face num.:  %d\n"
                   "  Local cell num.:  %d\n"), i+1, face_num, cur_cell_id+1);
#endif
      error = 1;
      *p_error = error;

      return retcode;
    }

    if (first_orient == -orient) {

      /* Inversed orientation between faces */

      if (first_orient == 1)
        ijkl_ref = i;

      first_orient = orient;

      /* Test orienation of [P, Q, V(i-1), V(i)] */

      vtx_id2 = face_connect[i-1] - 1;
      for (j = 0; j < 3; j++)
        vtx2[j] = mesh->vtx_coord[3*vtx_id2+j];

      orient_test = cs_lagrang_tetra_orientation(prev_location,
                                                 next_location,
                                                 vtx2,
                                                 vtx1);

      orient_test *= face_orient;
      orient_count += orient_test;

      if (orient_test == 0) { /* points are coplanar */
#if 1 && defined(DEBUG) && !defined(NDEBUG)
        bft_printf(_(" Lagrangian module warning.\n"
                     "  Failure during the particle tracking.\n"
                     "  [P, Q, V(%d), V(%d)] orientation failed"
                     "  => Particle is lost.\n"
                     "  Local face num.:  %d\n"
                     "  Local cell num.:  %d\n"),
                   i, i+1, face_num, cur_cell_id+1);
#endif
        error = 1;
        *p_error = error;

        return retcode;
      }

    } /* orient = -first_orient */

  } /* End of loop on face vertices */

  if (orient_count != -2 && orient_count != 0 && orient_count != 2) {
#if 1 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf(_(" Lagrangian module warning.\n"
                 "  Failure during the particle tracking.\n"
                 "  Local orientation counter must be -2, 0 or 2. Here is %d"
                 "  => Particle is lost.\n"
                 "  Local face num.:  %d\n"
                 "  Local cell num.:  %d\n"),
               orient_count, face_num, cur_cell_id+1);
#endif
    error = 1;
    *p_error = error;

    return retcode;
  }
  else if ( (orient_count == -2 || orient_count == 2) && ijkl_ref == 0 ) {
#if 1 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf(_(" Lagrangian module warning.\n"
                 "  Failure during the particle tracking.\n"
                 "  => Particle is lost.\n"
                 "  Local face num.:  %d\n"
                 "  Local cell num.:  %d\n"),
               face_num, cur_cell_id+1);
#endif
    error = 1;
    *p_error = error;

    return retcode;
  }

  if (orient_count == 0 || orient_count == -2) {
    retcode = 0;
    return retcode;
  }

  /* Relative position between the current face and the starting and ending
     particle locations */

  assert(orient_count == 2);

  vtx_id1 = face_connect[ijkl_ref - 1] - 1;
  vtx_id2 = face_connect[ijkl_ref] - 1;

  for (j = 0; j < 3; j++) {
    vtx1[j] = mesh->vtx_coord[3*vtx_id1+j];
    vtx2[j] = mesh->vtx_coord[3*vtx_id2+j];
  }

  orient = cs_lagrang_tetra_orientation(next_location,
                                        face_cog,
                                        vtx1,
                                        vtx2);

  if (orient == 0) { /* points are coplanar */
#if 1 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf(_(" Lagrangian module warning.\n"
                 "  Failure during the particle tracking.\n"
                 "  Do not find the relative position between:\n"
                 "   - P [%9.4f, %9.4f, %9.4f]\n"
                 "   - Q [%9.4f, %9.4f, %9.4f]\n"
                 "   - and the current face: CoG [%9.4f, %9.4f, %9.4f]\n"
                 "  => Particle is lost.\n"
                 "  Local face num.:  %d\n"
                 "  Local cell num.:  %d\n"),
               prev_location[0], prev_location[1], prev_location[2],
               next_location[0], next_location[1], next_location[2],
               face_cog[0], face_cog[1], face_cog[2],
               face_num, cur_cell_id+1);
#endif
    error = 1;
    *p_error = error;

    return retcode;
  }

  retcode = -orient * face_orient;

  /* Returns pointers */

  *p_error = error;

  return retcode;
}

/*----------------------------------------------------------------------------
 * Determine the number of the closest wall face from the particle
 * as well as the corresponding wall normal distance (y_p^+)
 *
 * Used for the deposition model.
 *
 * parameters:
 *   particle      <-> particle attributes for current time step
 *   p_am          <-- pointer to attributes map for current time step
 *   visc_length   <--
 *   dlgeo         <-- array with various geometry values for particles
 *----------------------------------------------------------------------------*/

static void
_test_wall_cell(void                           *particle,
                const cs_lagr_attribute_map_t  *p_am,
                const cs_real_t                 visc_length[],
                const cs_real_t                 dlgeo[])
{
  cs_lnum_t cell_num
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM);

  if (cell_num < 0) return;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;
  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;
  cs_lnum_t cell_id = cell_num - 1;

  cs_lagr_particle_set_real(particle, p_am, CS_LAGR_YPLUS, 10000);
  cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID, -1);

  cs_lnum_t  start = cell_face_idx[cell_id];
  cs_lnum_t  end =  cell_face_idx[cell_id + 1];

  for (cs_lnum_t i = start; i < end; i++) {
    cs_lnum_t  face_num = cell_face_lst[i];

    if (face_num < 0) {
      cs_lnum_t face_id = CS_ABS(face_num) - 1;
      cs_lnum_t b_zone_id = bdy_conditions->b_face_zone_num[face_id]-1;

      if (   (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IDEPO1)
          || (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IDEPO2)
          || (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IREBOL)
          || (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IDEPFA)) {

        cs_real_t x_face = dlgeo[face_id];
        cs_real_t y_face = dlgeo[face_id + (mesh->n_b_faces)];
        cs_real_t z_face = dlgeo[face_id + (mesh->n_b_faces) * 2];
        cs_real_t offset_face = dlgeo[face_id + (mesh->n_b_faces) * 3];
        const cs_real_t  *particle_coord
          = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_COORDS);

        cs_real_t dist_norm =   CS_ABS(  particle_coord[0] * x_face
                                       + particle_coord[1] * y_face
                                       + particle_coord[2] * z_face
                                       + offset_face) / visc_length[face_id];
        if (  dist_norm
            < cs_lagr_particle_get_real(particle, p_am, CS_LAGR_YPLUS)) {
          cs_lagr_particle_set_real(particle, p_am,
                                    CS_LAGR_YPLUS, dist_norm);
          cs_lagr_particle_set_lnum(particle, p_am,
                                    CS_LAGR_NEIGHBOR_FACE_ID, face_id);
        }
      }
    }

  }

}

/*----------------------------------------------------------------------------
 * Test if the current particle moves to the next cell through  this face
 *
 * parameters:
 *   set           <-- pointer to particle set
 *   prev_particle <-- particle attributes for previous time step
 *   particle      <-- particle attributes for current time step
 *   prev_p_am     <-- pointer to attributes map for previous time step
 *   p_am          <-- pointer to attributes map for current time step
 *   p_error       <-> pointer to an error indicator
 *
 * returns:
 *   particle state
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_bdy_treatment(cs_lagr_particle_set_t         *prev_set,
               cs_lagr_particle_set_t         *set,
               void                           *prev_particle,
               void                           *particle,
               cs_lnum_t                       face_num,
               cs_real_t                      *boundary_stat,
               cs_lnum_t                       boundary_zone,
               cs_lnum_t                       failsafe_mode,
               cs_lnum_t                      *p_move_particle,
               cs_lnum_t                      *p_n_failed_particles,
               cs_real_t                      *p_failed_particle_weight,
               const cs_lnum_t                *iensi3,
               const cs_lnum_t                *inbr,
               const cs_lnum_t                *inbrbd,
               const cs_lnum_t                *iflm,
               const cs_lnum_t                *iflmbd,
               const cs_lnum_t                *iang,
               const cs_lnum_t                *iangbd,
               const cs_lnum_t                *ivit,
               const cs_lnum_t                *ivitbd,
               const cs_lnum_t                *iencnb,
               const cs_lnum_t                *iencma,
               const cs_lnum_t                *iencdi,
               const cs_lnum_t                *iencck,
               const cs_lnum_t                *iencnbbd,
               const cs_lnum_t                *iencmabd,
               const cs_lnum_t                *iencdibd,
               const cs_lnum_t                *iencckbd,
               const cs_lnum_t                *inclg,
               const cs_lnum_t                *iscovc,
               const cs_lnum_t                *nusbor,
               cs_lnum_t                       iusb[],
               cs_real_t                       energt[],
               const cs_real_t                 tprenc[],
               const cs_real_t                 visref[],
               const cs_real_t                 enc1[],
               const cs_real_t                 enc2[],
               const cs_real_t                *tkelvi)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const double pi = 4 * atan(1);

  const cs_lagr_attribute_map_t  *prev_p_am = prev_set->p_am;
  const cs_lagr_attribute_map_t  *p_am = set->p_am;

  cs_lnum_t n_b_faces  = mesh->n_b_faces;

  cs_lnum_t  k;
  cs_real_t  tmp;
  cs_real_t  depl[3], face_normal[3], face_cog[3], intersect_pt[3];

  cs_real_t  compo_vel[3] = {0.0, 0.0, 0.0};
  cs_real_t  norm_vel = 0.0;

  cs_real_t  abs_curv = 0.0;
  cs_lnum_t  move_particle = *p_move_particle;
  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight = *p_failed_particle_weight;

  cs_lnum_t  face_id = face_num - 1;
  cs_lnum_t  particle_state = -999;
  cs_lnum_t  depch = 1; /* Indicator of mass flux calculation */

  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;

  cs_lnum_t  contact_number = 0;
  cs_real_t  *surface_coverage = NULL;

  cs_real_t  *prev_particle_coord
    = cs_lagr_particle_attr(prev_particle, prev_p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_coord
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_velocity
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);
  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

  const cs_real_t particle_stat_weight
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
  const cs_real_t particle_mass
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

  assert(bdy_conditions != NULL);

  for (k = 0; k < 3; k++)
    depl[k] = particle_coord[k] - prev_particle_coord[k];

  if (fabs(depl[0]) < 1e-15 && fabs(depl[1]) < 1e-15 && fabs(depl[2]) < 1e-15)
    return 0; /* move_particle = 0 */

  for (k = 0; k < 3; k++) {
    face_normal[k] = cs_glob_mesh_quantities->b_face_normal[3*face_id+k];
    face_cog[k] = cs_glob_mesh_quantities->b_face_cog[3*face_id+k];
  }

  cs_real_t face_area  = _get_norm(face_normal);

  cs_real_t  face_norm[3] = {face_normal[0]/face_area,
                             face_normal[1]/face_area,
                             face_normal[2]/face_area};

  /* Save particle impacting velocity */
  if (*iangbd > 0 || *ivitbd > 0) {
    norm_vel = _get_norm(particle_velocity);
    for (k = 0; k < 3; k++)
      compo_vel[k] = particle_velocity[k];
  }

  tmp = 0.0;
  for (k = 0; k < 3; k++)
    tmp += depl[k] * face_normal[k];

  if (fabs(tmp) < 1e-15)
    _manage_error(failsafe_mode,
                  particle,
                  p_am,
                  CS_LAGR_TRACKING_ERR_INSIDE_B_FACES,
                  &n_failed_particles,
                  &failed_particle_weight);

  /* 3D geometry refresher:
     ~~~~~~~~~~~~~~~~~~~~~~
     1)  equation of a plane with normal (a,b,c):
     a*x + b*y + c*z + d = 0
     2)  equation of a line passing through P and Q :
     x = XP + (XQ-XP) * AA
     y = YP + (YQ-YP) * AA
     z = ZP + (ZQ-ZP) * AA
     where AA is a real parameter
  */

  for (k = 0; k < 3; k++)
    abs_curv
      += face_normal[k]*face_cog[k] - face_normal[k]*prev_particle_coord[k];
  abs_curv /= tmp;

  for (k = 0; k < 3; k++)
    intersect_pt[k] = depl[k] * abs_curv + prev_particle_coord[k];

  if (   bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_ISORTL
      || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IENTRL
      || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_OUT;

    if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) {
      _particle_set->n_part_dep += 1;
      _particle_set->weight_dep += particle_stat_weight;
    }

    bdy_conditions->particle_flow_rate[boundary_zone]
      -= particle_stat_weight * particle_mass;

    /* FIXME: For post-processing by trajectory purpose */

    for (k = 0; k < 3; k++)
      particle_coord[k] = intersect_pt[k];
  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO2) {

    const cs_real_t particle_diameter
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

    move_particle = CS_LAGR_PART_MOVE_OFF;

    for (k = 0; k < 3; k++) {
      particle_velocity[k] = 0.0;
      particle_coord[k] =   intersect_pt[k]
                          - 0.5 * particle_diameter * face_norm[k];
    }

    _particle_set->n_part_dep += 1;
    _particle_set->weight_dep += particle_stat_weight;

    /* Specific treatment in case of particle resuspension modeling */

    cs_lnum_t *cur_cell_num = cs_lagr_particle_attr(particle,
                                                    p_am,
                                                    CS_LAGR_CUR_CELL_NUM);

    if (cs_glob_lagr_param.resuspension == 0) {

      *cur_cell_num = - *cur_cell_num;

      for (k = 0; k < 3; k++)
        particle_velocity_seen[k] = 0.0;

      particle_state = CS_LAGR_PART_STICKED;

    } else {

      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 1);
      *cur_cell_num = cs_glob_mesh->b_face_cells[face_id];

      particle_state = CS_LAGR_PART_TREATED;

    }

  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPFA) {

    const cs_real_t particle_diameter
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id]);

    cs_real_t uxn = particle_velocity[0] * face_norm[0];
    cs_real_t vyn = particle_velocity[1] * face_norm[1];
    cs_real_t wzn = particle_velocity[2] * face_norm[2];

    cs_real_t energ = 0.5 * particle_mass * (uxn+vyn+wzn) * (uxn+vyn+wzn);
    cs_real_t min_porosity ;
    cs_real_t limit;

    if (cs_glob_lagr_param.clogging) {

      /* If the clogging modeling is activated,                 */
      /* computation of the number of particles in contact with */
      /* the depositing particle                                */

      surface_coverage = &boundary_stat[(*iscovc -1) * n_b_faces + face_id];

      contact_number = cs_lagr_clogging_barrier(particle,
                                                p_am,
                                                face_id,
                                                face_area,
                                                &energt[face_id],
                                                surface_coverage,
                                                &limit,
                                                &min_porosity);

    }

    if (cs_glob_lagr_param.rough > 0)
      cs_lagr_roughness_barrier(particle,
                                p_am,
                                face_id,
                                &energt[face_id]);

    if (energ > energt[face_id] * 0.5 * particle_diameter) {

      /* The particle deposits*/
      if (!cs_glob_lagr_param.clogging && !cs_glob_lagr_param.resuspension) {
        move_particle = CS_LAGR_PART_MOVE_OFF;

        /* Set negative value for current cell number */
        cs_lagr_particle_set_lnum
          (particle, p_am, CS_LAGR_CUR_CELL_NUM,
           - cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM));

        _particle_set->n_part_dep += 1;
        _particle_set->weight_dep += particle_stat_weight;

        particle_state = CS_LAGR_PART_STICKED;
      }

      if (cs_glob_lagr_param.resuspension > 0) {

        move_particle = CS_LAGR_PART_MOVE_OFF;
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 1);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                  cs_glob_mesh->b_face_cells[face_id]);

        for (k = 0; k < 3; k++) {
          particle_velocity[k] = 0.0;
          particle_coord[k]
            = intersect_pt[k] - 0.5 * particle_diameter * face_norm[k];
        }
        _particle_set->n_part_dep += 1;
        _particle_set->weight_dep += particle_stat_weight;
        particle_state = CS_LAGR_PART_TREATED;

      }

      if (cs_glob_lagr_param.clogging) {

        cs_real_t depositing_radius = particle_diameter * 0.5;

        if (contact_number == 0) {

          /* The surface coverage increases if the particle
             has deposited on a naked surface */

          *surface_coverage += (pi * pow(depositing_radius,2))
                               *  particle_stat_weight / face_area;

          for (k = 0; k < 3; k++) {
            particle_coord[k]
              = intersect_pt[k] - (0.5 * particle_diameter * face_norm[k]);
            particle_velocity[k] = 0.0;
            particle_velocity_seen[k] = 0.0;
          }

        }
        else {

          cs_lnum_t i, j, one = 1, two = 2;

          cs_real_t nb_depo_part, spher_angle[2];
          cs_lnum_t contact;
          const void *cur_part = NULL;

          nb_depo_part = boundary_stat[(*inclg -1) * n_b_faces + face_id];

          double unit_vect[3];

          cs_lnum_t compt,compt2;
          cs_lnum_t compt_max = 100;
          cs_real_t dist;

          compt2 = 0;
          do {
            /* We choose randomly a deposited particle to interact
               with the depositing one */
            cs_real_t random;
            CS_PROCF(zufall, ZUFALL)(&one, &random);
            int part_num = (int) (random * nb_depo_part);

            k = 0;
            for (i = 0, j = _particle_set->first_used_id;
                 i < _particle_set->n_particles;
                 i++) {

              cur_part = (const void *)(set->p_buffer + p_am->extents * j);

              cs_lnum_t cur_part_depo
                = cs_lagr_particle_get_lnum(cur_part, p_am,
                                            CS_LAGR_DEPOSITION_FLAG);
              cs_lnum_t cur_part_close_face_id
                = cs_lagr_particle_get_lnum(cur_part, p_am,
                                            CS_LAGR_NEIGHBOR_FACE_ID);

              if ((cur_part_depo) && (cur_part_close_face_id == face_id)) {

                if (k == part_num)
                  break;
                else
                  k += 1;

              }
              j = _particle_set->used_id[j].next_id;
            }

            double norm_vect = 0.5 * (  cs_lagr_particle_get_real
                                          (cur_part, p_am, CS_LAGR_DIAMETER)
                                      + particle_diameter);

            compt2 += 1;
            compt = 0;
            do {
              do {
                /* The depositing particle is relocated in contact
                   with the randomly chosen deposited */

                CS_PROCF(zufall, ZUFALL)(&two, spher_angle);

                spher_angle[0] *= pi;
                spher_angle[1] *= 2 * pi;

                unit_vect[0] =  sin(spher_angle[0]) * cos(spher_angle[1]);
                unit_vect[1] =  sin(spher_angle[0]) * sin(spher_angle[1]);
                unit_vect[2] =  cos(spher_angle[0]);
              } while (_get_dot_prod(unit_vect,face_norm) > -0.3);

              const cs_real_t  *cur_part_coord
                = cs_lagr_particle_attr_const(cur_part, p_am,
                                              CS_LAGR_COORDS);
              for (k = 0; k < 3; k++) {
                particle_velocity[k] = 0.0;
                particle_coord[k]
                  = cur_part_coord[k] + unit_vect[k] * norm_vect;
              }

              k = 0;
              contact=0;
              for (i = 0, j = _particle_set->first_used_id;
                   i < _particle_set->n_particles;
                   i++) {

                const void * cur_part2 = set->p_buffer + p_am->extents * j;

                cs_lnum_t cur_part2_depo
                  = cs_lagr_particle_get_lnum(cur_part2, p_am,
                                              CS_LAGR_DEPOSITION_FLAG);
                cs_lnum_t cur_part2_close_face_id
                  = cs_lagr_particle_get_lnum(cur_part2, p_am,
                                              CS_LAGR_NEIGHBOR_FACE_ID);

                /* Calculation of the distance of two particles */
                if ((cur_part2_depo) && (cur_part2_close_face_id == face_id)) {

                  const cs_real_t  *cur_part2_coord
                    = cs_lagr_particle_attr_const(cur_part2, p_am,
                                                  CS_LAGR_COORDS);
                  cs_real_t cur_part2_diameter
                  = cs_lagr_particle_get_real(cur_part2, p_am,
                                              CS_LAGR_DIAMETER);

                  dist = sqrt(_delta_norm_2_3d(particle_coord,
                                               cur_part2_coord));

                  if (dist < (cur_part2_diameter/2 + particle_diameter/2))
                    contact = contact + 1;
                }
                j = _particle_set->used_id[j].next_id;
              }

              compt += 1;

            } while (contact != 0 && compt < compt_max);
            /* Test of an other angle if contact between particles */

          } while (contact != 0 && compt2 < compt_max);
          /* Test to prevent the covering of particles */
        }

        cs_lnum_t  ncel = cs_glob_mesh->n_cells;
        cs_lnum_t  node;
        cs_real_t volp[ncel];
        cs_real_t porosity = 0.;
        const cs_real_t *xyzcen = cs_glob_mesh_quantities->cell_cen;
        const cs_real_t *volume  = cs_glob_mesh_quantities->cell_vol;

        /* Determination of the cell number of the particle */
        /* FIXME for parallel cases */

        node = (ncel + 1) / 2;

        cs_real_t xyz1[3] = {xyzcen[3 * (node - 1)],
                             xyzcen[3 * (node - 1) + 1],
                             xyzcen[3 * (node - 1) + 2]};

        cs_real_t dis2mn = _delta_norm_2_3d(particle_coord, xyz1);

        for (cs_lnum_t ii = 0; ii < ncel ; ii++) {
          xyz1[0] = xyzcen[3*ii];
          xyz1[1] = xyzcen[3*ii + 1];
          xyz1[2] = xyzcen[3*ii + 2];

          cs_real_t dis2 = _delta_norm_2_3d(particle_coord, xyz1);

          if (dis2 < dis2mn) {
            node = ii + 1;
            dis2mn = dis2;
          }
        }

        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                  node);

        cs_lnum_t p_cur_cell_id = node - 1;

        /* Calculation of the cell porosity */
        porosity = (  volume[p_cur_cell_id]
                    - volp[p_cur_cell_id]) / volume[p_cur_cell_id];

        if (porosity > min_porosity) {

          move_particle = CS_LAGR_PART_MOVE_OFF;

          volp[p_cur_cell_id]
            +=  particle_stat_weight * 4./3. * pi * pow(particle_diameter/2,3);

          /* Set negative value for current cell number */
          cs_lagr_particle_set_lnum
            (particle, p_am, CS_LAGR_CUR_CELL_NUM,
             - cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM));

          _particle_set->n_part_dep += 1;
          _particle_set->weight_dep += particle_stat_weight;

          particle_state = CS_LAGR_PART_STICKED;
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 1);

          /* Update of the number of stat. weight of deposited particles */
          /* in case of clogging modeling                                */

          boundary_stat[(*inclg -1) * n_b_faces + face_id]
            += particle_stat_weight;

        }
        else {
          /*The particle does not deposit: It 'rebounds' */
          /* because the porosity threshold is reached   */

          /* The particle mass flux is not calculated */
          depch = 0;

          move_particle = CS_LAGR_PART_MOVE_ON;
          particle_state = CS_LAGR_PART_TO_SYNC;
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                    cs_glob_mesh->b_face_cells[face_id]);

          for (k = 0; k < 3; k++)
            prev_particle_coord[k] = intersect_pt[k];

          /* Modify the ending point. */

          for (k = 0; k < 3; k++)
            depl[k] = particle_coord[k] - intersect_pt[k];

          tmp = CS_ABS(_get_dot_prod(depl, face_norm));
          tmp *= 2.0;

          for (k = 0; k < 3; k++)
            particle_coord[k] -= tmp * face_norm[k];

          /* Modify particle velocity and velocity seen */

          tmp = CS_ABS(_get_dot_prod(particle_velocity, face_norm));
          tmp *= 2.0;

          for (k = 0; k < 3; k++)
            particle_velocity[k] -= tmp * face_norm[k];

          tmp = CS_ABS(_get_dot_prod(particle_velocity_seen, face_norm));
          tmp *= 2.0;

          for (k = 0; k < 3; k++)
            particle_velocity_seen[k] -= tmp * face_norm[k];
        }

      }

      /* For post-processing purpose */

      if (!cs_glob_lagr_param.clogging && !cs_glob_lagr_param.resuspension ) {
        for (k = 0; k < 3; k++) {
          particle_coord[k] = intersect_pt[k]
                              - (0.5 * particle_diameter * face_norm[k]);
          particle_velocity[k] = 0.0;
          particle_velocity_seen[k] = 0.0;
        }
      }
    }

    else  {
      /*The particle does not deposit:
        It 'rebounds' on the energy barrier*/

      /* The particle mass flux is not calculated */
      depch = 0;

      move_particle = CS_LAGR_PART_MOVE_ON;
      particle_state = CS_LAGR_PART_TO_SYNC;
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                cs_glob_mesh->b_face_cells[face_id]);

      for (k = 0; k < 3; k++)
        prev_particle_coord[k] = intersect_pt[k];

      /* Modify the ending point. */

      for (k = 0; k < 3; k++)
        depl[k] = particle_coord[k] - intersect_pt[k];

      tmp = CS_ABS(_get_dot_prod(depl, face_norm));
      tmp *= 2.0;

      for (k = 0; k < 3; k++)
        particle_coord[k] -= tmp * face_norm[k];

      /* Modify particle velocity and velocity seen */

      tmp = CS_ABS(_get_dot_prod(particle_velocity, face_norm));
      tmp *= 2.0;

      for (k = 0; k < 3; k++)
        particle_velocity[k] -= tmp * face_norm[k];

      tmp = CS_ABS(_get_dot_prod(particle_velocity_seen, face_norm));
      tmp *= 2.0;

      for (k = 0; k < 3; k++)
        particle_velocity_seen[k] -= tmp * face_norm[k];
    }
  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IREBOL) {

    move_particle = CS_LAGR_PART_MOVE_ON;
    particle_state = CS_LAGR_PART_TO_SYNC;
    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id]);

    for (k = 0; k < 3; k++)
      prev_particle_coord[k] = intersect_pt[k];

    /* Modify the ending point. */

    for (k = 0; k < 3; k++)
      depl[k] = particle_coord[k] - intersect_pt[k];

    tmp = CS_ABS(_get_dot_prod(depl, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle_coord[k] -= tmp * face_norm[k];

    /* Modify particle velocity and velocity seen */

    tmp = CS_ABS(_get_dot_prod(particle_velocity, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle_velocity[k] -= tmp * face_norm[k];

    tmp = CS_ABS(_get_dot_prod(particle_velocity_seen, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle_velocity_seen[k] -= tmp * face_norm[k];

  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_ISYMTL) {

    move_particle = CS_LAGR_PART_MOVE_ON;
    particle_state = CS_LAGR_PART_TO_SYNC;
    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id]);

    for (k = 0; k < 3; k++)
      prev_particle_coord[k] = intersect_pt[k];

    /* Modify the ending point. */

    for (k = 0; k < 3; k++)
      depl[k] = particle_coord[k] - intersect_pt[k];

    tmp = CS_ABS(_get_dot_prod(depl, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle_coord[k] -= tmp * face_norm[k];

    /* Modify particle velocity and velocity seen */

    tmp = CS_ABS(_get_dot_prod(particle_velocity, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle_velocity[k] -= tmp * face_norm[k];

    tmp = CS_ABS(_get_dot_prod(particle_velocity_seen, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle_velocity_seen[k] -= tmp * face_norm[k];

  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IENCRL) {

    /* Fouling of the particle, if its properties make it possible and
       with respect to a probability
       ICI if  Tp     > TPENC
           if viscp <= VISREF ==> Probability of fouling equal to 1
           if viscp  > VISREF ==> Probability equal to TRAP = 1-VISREF/viscp
                              ==> Fouling if VNORL is between TRAP et 1. */

    cs_real_t  random = -1, viscp = -1, trap = -1;

    /* Selection of the fouling coefficient*/

    const cs_lnum_t particle_coal_number
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_COAL_NUM);
    const cs_lnum_t n_layers = p_am->count[CS_LAGR_TEMPERATURE];
    const cs_real_t *particle_temp
      = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_TEMPERATURE);

    cs_real_t  temp_ext_part = particle_temp[n_layers - 1];
    cs_real_t  tprenc_icoal  = tprenc[particle_coal_number - 1];
    cs_real_t  visref_icoal  = visref[particle_coal_number - 1];
    cs_real_t  enc1_icoal    = enc1[particle_coal_number - 1];
    cs_real_t  enc2_icoal    = enc2[particle_coal_number - 1];

    if (temp_ext_part > tprenc_icoal+*tkelvi) {

      /* Coal viscosity*/
      tmp = (  (1.0e7*enc1_icoal)
             / ((temp_ext_part-150.e0-*tkelvi)*(temp_ext_part-150.e0-*tkelvi)))
            + enc2_icoal;
      if (tmp <= 0.0) {
        bft_error
          (__FILE__, __LINE__, 0,
           _("Coal viscosity calculation impossible, tmp = %e is < 0.\n"),
           tmp);
      }
      else
        viscp = 0.1e0 * exp(log(10.e0)*tmp);

      if (viscp >= visref_icoal) {
        int  one = 1;
        CS_PROCF(zufall, ZUFALL)(&one, &random);
        trap = 1.e0- (visref_icoal / viscp);
      }

      if (   (viscp <= visref_icoal)
          || (viscp >= visref_icoal  &&  random >= trap)) {

        move_particle = CS_LAGR_PART_MOVE_OFF;
        particle_state = CS_LAGR_PART_OUT;

        /* Recording for listing/listla*/
        _particle_set->n_part_fou += 1;
        _particle_set->weight_fou += particle_stat_weight;

        /* Recording for statistics*/
        if (*iencnbbd > 0) {
          boundary_stat[(*iencnb -1) * n_b_faces + face_id]
            += particle_stat_weight;
        }
        if (*iencmabd > 0) {
          boundary_stat[(*iencma -1) * n_b_faces + face_id]
            += particle_stat_weight * particle_mass / face_area;
        }
        if (*iencdibd > 0) {
          boundary_stat[(*iencdi -1) * n_b_faces + face_id]
            +=   particle_stat_weight
               * cs_lagr_particle_get_real(particle, p_am,
                                           CS_LAGR_SHRINKING_DIAMETER);
        }
        if (*iencckbd > 0) {
          if (particle_mass > 0) {
            const cs_real_t *particle_coal_mass
              = cs_lagr_particle_attr_const(particle, p_am,
                                            CS_LAGR_COAL_MASS);
            const cs_real_t *particle_coke_mass
              = cs_lagr_particle_attr_const(particle, p_am,
                                            CS_LAGR_COKE_MASS);
            for (k = 0; k < n_layers; k++) {
              boundary_stat[(*iencck -1) * n_b_faces + face_id]
                +=   particle_stat_weight
                   * (particle_coal_mass[k] + particle_coke_mass[k])
                   / particle_mass;
            }
          }
        }

        /* FIXME: For post-processing by trajectory purpose */

        for (k = 0; k < 3; k++) {
          particle_coord[k] = intersect_pt[k];
          particle_velocity[k] = 0.0;
          particle_velocity_seen[k] = 0.0;
        }
      }
    }

    /*--> if there is no fouling, then it is an elastic rebound*/
    if (move_particle != CS_LAGR_PART_MOVE_OFF) {

      move_particle = CS_LAGR_PART_MOVE_ON;
      particle_state = CS_LAGR_PART_TO_SYNC;
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                cs_glob_mesh->b_face_cells[face_id]);

      for (k = 0; k < 3; k++)
        prev_particle_coord[k] = intersect_pt[k];

      /* Modify the ending point. */

      for (k = 0; k < 3; k++)
        depl[k] = particle_coord[k] - intersect_pt[k];

      tmp = CS_ABS(_get_dot_prod(depl, face_norm));
      tmp *= 2.0;

      for (k = 0; k < 3; k++)
        particle_coord[k] -= tmp * face_norm[k];

      /* Modify particle velocity and velocity seen */

      tmp = CS_ABS(_get_dot_prod(particle_velocity, face_norm));
      tmp *= 2.0;

      for (k = 0; k < 3; k++)
        particle_velocity[k] -= tmp * face_norm[k];

      tmp = CS_ABS(_get_dot_prod(particle_velocity_seen, face_norm));
      tmp *= 2.0;

      for (k = 0; k < 3; k++) {
        particle_velocity_seen[k] -= tmp * face_norm[k];
        particle_velocity_seen[k] = 0.0;}

    }

  }

  /* FIXME: JBORD* (user-defined boundary condition) not yet implemented
     nor defined by a macro */
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Boundary condition %d not recognized.\n"),
              bdy_conditions->b_zone_natures[boundary_zone]);

  /* FIXME: Post-treatment not yet implemented... */

  /* Return pointer */

  *p_move_particle = move_particle;
  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;

  if (*iensi3 > 0) {
    if  (   bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO2
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPFA
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IREBOL
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IENCRL) {

      /* Number of particle-boundary interactions  */
      if (*inbrbd > 0)
        boundary_stat[(*inbr -1) * n_b_faces + face_id]
          += particle_stat_weight;

      /* Particulate boundary mass flux */
      if ((*iflmbd > 0) && (depch == 1))
        boundary_stat[(*iflm -1) * n_b_faces + face_id]
          += particle_stat_weight * particle_mass / face_area;

      /* Particle impact angle and velocity*/
      if (*iangbd > 0) {
        cs_real_t imp_ang = acos(_get_dot_prod(compo_vel, face_normal)
                                 / (face_area * norm_vel));
        boundary_stat[(*iang -1) * n_b_faces + face_id]
          += imp_ang * particle_stat_weight;
      }

      if (*ivitbd > 0)
        boundary_stat[(*ivit -1) * n_b_faces + face_id]
          += norm_vel * particle_stat_weight;

      /* User statistics management. By defaut, set to zero */
      if (*nusbor > 0)
        for (int n1 = 0; n1 < *nusbor; n1++)
            boundary_stat[(iusb[n1] -1) * n_b_faces + face_id] = 0.0;
    }
  }

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Move locally a particle as far as it is possible while keeping on the same
 * rank.
 *
 * parameters:
 *   prev_particle            <->  pointer to particle data
 *   particle                 <->  pointer to particle data
 *   scheme_order             <--  current order of the scheme used for
 *                                 the lagrangian algorithm (1 or 2)
 *   failsafe_mode            <--  with (0) / without (1) failure capability
 *   p_n_failed_particles     <-> number of failed particles
 *   p_failed_particle_weight <-> stat. weight of the failed particles
 *
 * returns:
 *   a state associated to the status of the particle (treated, to be deleted,
 *   to be synchonised)
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_local_propagation(void                           *prev_particle,
                   void                           *particle,
                   const cs_lagr_attribute_map_t  *prev_p_am,
                   const cs_lagr_attribute_map_t  *p_am,
                   int                             scheme_order,
                   int                             failsafe_mode,
                   cs_real_t                       boundary_stat[],
                   cs_lnum_t                      *p_n_failed_particles,
                   cs_real_t                      *p_failed_particle_weight,
                   const cs_lnum_t                *iensi3,
                   const cs_lnum_t                *inbr,
                   const cs_lnum_t                *inbrbd,
                   const cs_lnum_t                *iflm,
                   const cs_lnum_t                *iflmbd,
                   const cs_lnum_t                *iang,
                   const cs_lnum_t                *iangbd,
                   const cs_lnum_t                *ivit,
                   const cs_lnum_t                *ivitbd,
                   const cs_lnum_t                *iencnb,
                   const cs_lnum_t                *iencma,
                   const cs_lnum_t                *iencdi,
                   const cs_lnum_t                *iencck,
                   const cs_lnum_t                *iencnbbd,
                   const cs_lnum_t                *iencmabd,
                   const cs_lnum_t                *iencdibd,
                   const cs_lnum_t                *iencckbd,
                   const cs_lnum_t                *inclg,
                   const cs_lnum_t                *iscovc,
                   const cs_lnum_t                *nusbor,
                   cs_lnum_t                       iusb[],
                   cs_real_t                       visc_length[],
                   cs_real_t                       dlgeo[],
                   const cs_field_t               *u,
                   cs_real_t                       energt[],
                   const cs_real_t                 tprenc[],
                   const cs_real_t                 visref[],
                   const cs_real_t                 enc1[],
                   const cs_real_t                 enc2[],
                   const cs_real_t                *tkelvi,
                   cs_lnum_t                       ipass)
{
  cs_lnum_t  i, j, k;
  cs_real_t  depl[3];
  cs_real_t  null_yplus;

  cs_lnum_t  *neighbor_face_id;
  cs_real_t  *particle_yplus;

  cs_lnum_t  error = 0;
  cs_lnum_t  n_loops = 0;
  cs_lnum_t  move_particle = CS_LAGR_PART_MOVE_ON;
  cs_lnum_t  particle_state = CS_LAGR_PART_TO_SYNC;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_lagr_param_t *lagr_param = &cs_glob_lagr_param;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;
  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;
  cs_lnum_t  *face_connect = builder->face_connect_buffer;

  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight = *p_failed_particle_weight;
  cs_lnum_t  new_face_num
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_LAST_FACE_NUM);

  cs_real_t  *prev_particle_coord
    = cs_lagr_particle_attr(prev_particle, prev_p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_coord
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);

  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

  for (k = 0; k < 3; k++)
    depl[k] = particle_coord[k] - prev_particle_coord[k];

  if (fabs(depl[0]) < 1e-15 && fabs(depl[1]) < 1e-15 && fabs(depl[2]) < 1e-15) {
    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_TREATED;
  }

  if (lagr_param->deposition > 0) {
    neighbor_face_id
      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
    particle_yplus
      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_YPLUS);
  }
  else {
    neighbor_face_id = NULL;
    null_yplus = 0;
    particle_yplus = &null_yplus;  /* allow tests even without particle y+ */
  }

  /*  particle_state is defined at the top of this file */

  while (move_particle == CS_LAGR_PART_MOVE_ON) {

    cs_lnum_t  cur_cell_id
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM) - 1;
    cs_lnum_t  old_face_num = new_face_num;

    assert(cur_cell_id < mesh->n_cells);
    assert(cur_cell_id > -1);

    cs_lnum_t  start = cell_face_idx[cur_cell_id];
    cs_lnum_t  end =  cell_face_idx[cur_cell_id+1];

    n_loops++;

    if (n_loops > CS_LAGR_MAX_PROPAGATION_LOOPS) { /* Manage error */

      if (ipass == 1)
        _manage_error(failsafe_mode,
                      particle,
                      p_am,
                      CS_LAGR_TRACKING_ERR_MAX_LOOPS,
                      &n_failed_particles,
                      &failed_particle_weight);

      move_particle  = CS_LAGR_PART_MOVE_OFF;
      particle_state = CS_LAGR_PART_ERR;

    }

    /*Treatment for particles which change rank*/

    if (lagr_param->deposition > 0 && *particle_yplus < 0.) {

      _test_wall_cell(particle, p_am, visc_length, dlgeo);

      if (*particle_yplus < 100.) {

        cs_real_t flow_velo_x, flow_velo_y, flow_velo_z;
        if (u->interleaved) {
          flow_velo_x = u->val[cur_cell_id*3];
          flow_velo_y = u->val[cur_cell_id*3 + 1];
          flow_velo_z = u->val[cur_cell_id*3 + 2];
        }
        else {
          flow_velo_x = u->val[cur_cell_id];
          flow_velo_y = u->val[cur_cell_id + mesh->n_cells_with_ghosts];
          flow_velo_z = u->val[cur_cell_id + 2*mesh->n_cells_with_ghosts];
        }

        /* e1 (normal) vector coordinates */
        cs_real_t e1_x = dlgeo[*neighbor_face_id];
        cs_real_t e1_y = dlgeo[*neighbor_face_id + n_b_faces];
        cs_real_t e1_z = dlgeo[*neighbor_face_id + n_b_faces*2];

        /* e2 vector coordinates */
        cs_real_t e2_x = dlgeo[*neighbor_face_id + n_b_faces*7];
        cs_real_t e2_y = dlgeo[*neighbor_face_id + n_b_faces*8];
        cs_real_t e2_z = dlgeo[*neighbor_face_id + n_b_faces*9];

        /* e3 vector coordinates */
        cs_real_t e3_x = dlgeo[*neighbor_face_id + n_b_faces*10];
        cs_real_t e3_y = dlgeo[*neighbor_face_id + n_b_faces*11];
        cs_real_t e3_z = dlgeo[*neighbor_face_id + n_b_faces*12];

        /* V_n * e1 */

        cs_real_t v_n_e1[3] = {particle_velocity_seen[0] * e1_x,
                               particle_velocity_seen[0] * e1_y,
                               particle_velocity_seen[0] * e1_z};

        /* (U . e2) * e2 */

        cs_real_t flow_e2 =   flow_velo_x * e2_x
                            + flow_velo_y * e2_y
                            + flow_velo_z * e2_z;

        cs_real_t u_e2[3] = {flow_e2 * e2_x,
                             flow_e2 * e2_y,
                             flow_e2 * e2_z};

        /* (U . e3) * e3 */

        cs_real_t flow_e3 =   flow_velo_x * e3_x
                            + flow_velo_y * e3_y
                            + flow_velo_z * e3_z;

        cs_real_t u_e3[3] = {flow_e3 * e3_x,
                             flow_e3 * e3_y,
                             flow_e3 * e3_z};

        /* Update of the flow seen velocity */

        particle_velocity_seen[0] =  v_n_e1[0] + u_e2[0] + u_e3[0];
        particle_velocity_seen[1] =  v_n_e1[1] + u_e2[1] + u_e3[1];
        particle_velocity_seen[2] =  v_n_e1[2] + u_e2[2] + u_e3[2];
      }
    }

    /* Loop on faces connected to the current cell */

    for (i = start; i < end && move_particle == CS_LAGR_PART_MOVE_ON; i++) {

      cs_lnum_t  face_num = cell_face_lst[i];
      cs_lnum_t  indian = 0;

      if (face_num > 0 && face_num != old_face_num) {

        /* Interior face which is different from the incoming face */

        cs_lnum_t  face_id = face_num - 1;
        cs_lnum_t  vtx_start = mesh->i_face_vtx_idx[face_id] - 1;
        cs_lnum_t  vtx_end = mesh->i_face_vtx_idx[face_id+1] - 1;
        cs_lnum_t  n_vertices = vtx_end - vtx_start + 1;

        for (k = 0, j = vtx_start; j < vtx_end; j++, k++)
          face_connect[k] = mesh->i_face_vtx_lst[j];
        face_connect[n_vertices-1] = face_connect[0];

        /*
          indian = -1 : keep inside the same cell.
          indian = 0  : trajectory doesn't go through the current face.
          indian = 1  : trajectory goes through the current face.
        */

        indian = _where_are_you(face_num,
                                n_vertices,
                                face_connect,
                                prev_particle,
                                particle,
                                prev_p_am,
                                p_am,
                                &error);

        if (error == 1 && ipass == 1) {

          _manage_error(failsafe_mode,
                        particle,
                        p_am,
                        CS_LAGR_TRACKING_ERR_DISPLACEMENT_I_FACE,
                        &n_failed_particles,
                        &failed_particle_weight);

          move_particle = CS_LAGR_PART_MOVE_OFF;
          particle_state = CS_LAGR_PART_ERR;

        }
        else if (indian == 1) { /* Particle moves to the neighbor cell
                                   through the current face "face_num" */

          cs_lnum_t  cell_num1 = mesh->i_face_cells[2*face_id];
          cs_lnum_t  cell_num2 = mesh->i_face_cells[2*face_id+1];

          new_face_num = face_num;

          if (cur_cell_id + 1 == cell_num1) {
            cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                      cell_num2);
            cur_cell_id = cell_num2 - 1;
          }

          else {
            cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CUR_CELL_NUM,
                                      cell_num1);
            cur_cell_id = cell_num1 - 1;
          }

          if (cur_cell_id >= mesh->n_cells) {

            particle_state = CS_LAGR_PART_TO_SYNC;
            move_particle = CS_LAGR_PART_ERR;

            if (lagr_param->deposition > 0 && *particle_yplus < 100.) {

              cs_real_t x_p_q = particle_coord[0] - prev_particle_coord[0];
              cs_real_t y_p_q = particle_coord[1] - prev_particle_coord[1];
              cs_real_t z_p_q = particle_coord[2] - prev_particle_coord[2];

              cs_real_t face_normal[3], face_cog[3];

              for (k = 0; k < 3; k++) {
                face_normal[k] = mesh_quantities->i_face_normal[3*face_id+k];
                face_cog[k] = mesh_quantities->i_face_cog[3*face_id+k];
              }

              cs_real_t aa =   x_p_q * face_normal[0]
                             + y_p_q * face_normal[1]
                             + z_p_q * face_normal[2];


              cs_real_t bb = (  face_normal[0] * face_cog[0]
                              + face_normal[1] * face_cog[1]
                              + face_normal[2] * face_cog[2]
                              - face_normal[0] * prev_particle_coord[0]
                              - face_normal[1] * prev_particle_coord[1]
                              - face_normal[2] * prev_particle_coord[2]) / aa;

              cs_real_t xk =  prev_particle_coord[0] + bb * x_p_q;
              cs_real_t yk =  prev_particle_coord[1] + bb * y_p_q;
              cs_real_t zk =  prev_particle_coord[2] + bb * z_p_q;

              cs_real_t *xyzcen = mesh_quantities->cell_cen;

              particle_coord[0] = xk + 1e-8 * (xyzcen[3*cur_cell_id] - xk);
              particle_coord[1] = yk + 1e-8 * (xyzcen[3*cur_cell_id + 1] - yk);
              particle_coord[2] = zk + 1e-8 * (xyzcen[3*cur_cell_id + 2] - zk);

              /* Marking of particles */

              *particle_yplus = - *particle_yplus;

              /*Saving of dot product */

              cs_real_t bdy_normal[3];

              for (k = 0; k < 3; k++)
                bdy_normal[k]
                  = mesh_quantities->b_face_normal[(*neighbor_face_id)*3 + k];

              cs_real_t area  = _get_norm(bdy_normal);

              cs_real_t  face_norm[3] = {bdy_normal[0]/area,
                                         bdy_normal[1]/area,
                                         bdy_normal[2]/area};

              particle_velocity_seen[0]
                =   particle_velocity_seen[0] * face_norm[0]
                  + particle_velocity_seen[1] * face_norm[1]
                  + particle_velocity_seen[2] * face_norm[2];

            }

          }
          else {

            /* Specific treatment for the particle deposition model */

            if (lagr_param->deposition > 0) {

              cs_lnum_t save_close_face_id = *neighbor_face_id;
              cs_real_t save_yplus = *particle_yplus;

              /*Wall cell detection */

              _test_wall_cell(particle, p_am, visc_length, dlgeo);

              if (save_yplus < 100.) {

                cs_real_t x_p_q = particle_coord[0] - prev_particle_coord[0];
                cs_real_t y_p_q = particle_coord[1] - prev_particle_coord[1];
                cs_real_t z_p_q = particle_coord[2] - prev_particle_coord[2];

                cs_real_t face_normal[3], face_cog[3];

                for (k = 0; k < 3; k++) {
                  face_normal[k]
                    = mesh_quantities->i_face_normal[3*face_id+k];
                  face_cog[k] = mesh_quantities->i_face_cog[3*face_id+k];
                }

                cs_real_t aa =   x_p_q * face_normal[0]
                               + y_p_q * face_normal[1]
                               + z_p_q * face_normal[2];


                cs_real_t bb = (  face_normal[0] * face_cog[0]
                                + face_normal[1] * face_cog[1]
                                + face_normal[2] * face_cog[2]
                                - face_normal[0] *  prev_particle_coord[0]
                                - face_normal[1] *  prev_particle_coord[1]
                                - face_normal[2] *  prev_particle_coord[2]) / aa;

                cs_real_t xk =  prev_particle_coord[0] + bb * x_p_q;
                cs_real_t yk =  prev_particle_coord[1] + bb * y_p_q;
                cs_real_t zk =  prev_particle_coord[2] + bb * z_p_q;

                cs_real_t* xyzcen = mesh_quantities->cell_cen;

                particle_coord[0] = xk + 1e-8 * (xyzcen[3*cur_cell_id] - xk);
                particle_coord[1] = yk + 1e-8 * (xyzcen[3*cur_cell_id + 1] - yk);
                particle_coord[2] = zk + 1e-8 * (xyzcen[3*cur_cell_id + 2] - zk);

                /* Second test with the new particle position */

                _test_wall_cell(particle, p_am, visc_length, dlgeo);

                if (*particle_yplus < 100.e0) {

                  cs_real_t flow_velo_x, flow_velo_y, flow_velo_z;
                  if (u->interleaved) {
                    flow_velo_x = u->val[cur_cell_id*3];
                    flow_velo_y = u->val[cur_cell_id*3 + 1];
                    flow_velo_z = u->val[cur_cell_id*3 + 2];
                  }
                  else {
                    flow_velo_x = u->val[  cur_cell_id];
                    flow_velo_y = u->val[  cur_cell_id
                                         + mesh->n_cells_with_ghosts];
                    flow_velo_z = u->val[  cur_cell_id
                                         + 2*mesh->n_cells_with_ghosts];
                  }

                  /* The particle is still in the boundary layer */

                  cs_real_t old_bdy_normal[3];

                  for (k = 0; k < 3; k++)
                    old_bdy_normal[k]
                      = mesh_quantities->b_face_normal[3*save_close_face_id + k];

                  cs_real_t old_area  = _get_norm(old_bdy_normal);

                  cs_real_t  old_face_norm[3] = {old_bdy_normal[0]/old_area,
                                                 old_bdy_normal[1]/old_area,
                                                 old_bdy_normal[2]/old_area};

                  cs_real_t old_fl_seen_norm
                    =   particle_velocity_seen[0] * old_face_norm[0]
                      + particle_velocity_seen[1] * old_face_norm[1]
                      + particle_velocity_seen[2] * old_face_norm[2];

                  /* e1 (normal) vector coordinates */
                  cs_real_t e1_x = dlgeo[*neighbor_face_id];
                  cs_real_t e1_y = dlgeo[*neighbor_face_id + n_b_faces];
                  cs_real_t e1_z = dlgeo[*neighbor_face_id + n_b_faces*2];

                  /* e2 vector coordinates */
                  cs_real_t e2_x = dlgeo[*neighbor_face_id + n_b_faces*7];
                  cs_real_t e2_y = dlgeo[*neighbor_face_id + n_b_faces*8];
                  cs_real_t e2_z = dlgeo[*neighbor_face_id + n_b_faces*9];

                  /* e3 vector coordinates */
                  cs_real_t e3_x = dlgeo[*neighbor_face_id + n_b_faces*10];
                  cs_real_t e3_y = dlgeo[*neighbor_face_id + n_b_faces*11];
                  cs_real_t e3_z = dlgeo[*neighbor_face_id + n_b_faces*12];

                  /* V_n * e1 */

                  cs_real_t v_n_e1[3] = {old_fl_seen_norm * e1_x,
                                         old_fl_seen_norm * e1_y,
                                         old_fl_seen_norm * e1_z};

                  /* (U . e2) * e2 */

                  cs_real_t flow_e2 =   flow_velo_x * e2_x
                                      + flow_velo_y * e2_y
                                      + flow_velo_z * e2_z;

                  cs_real_t u_e2[3] = {flow_e2 * e2_x,
                                       flow_e2 * e2_y,
                                       flow_e2 * e2_z};

                  /* (U . e3) * e3 */

                  cs_real_t flow_e3 =   flow_velo_x * e3_x
                                      + flow_velo_y * e3_y
                                      + flow_velo_z * e3_z;

                  cs_real_t u_e3[3] = {flow_e3 * e3_x,
                                       flow_e3 * e3_y,
                                       flow_e3 * e3_z};

                  /* Update of the flow seen velocity */

                  particle_velocity_seen[0] = v_n_e1[0] + u_e2[0] + u_e3[0];
                  particle_velocity_seen[1] = v_n_e1[1] + u_e2[1] + u_e3[1];
                  particle_velocity_seen[2] = v_n_e1[2] + u_e2[2] + u_e3[2];
                }

                move_particle =  CS_LAGR_PART_MOVE_OFF;
                particle_state = CS_LAGR_PART_TREATED;

              }
            }
          }
        } else if (indian == -1) {

          move_particle =  CS_LAGR_PART_MOVE_OFF;
          particle_state = CS_LAGR_PART_TREATED;
        }

      } /* End if face_num > 0 (interior face) && face_num != old_face_num */

      else if (face_num < 0 &&  face_num != old_face_num) {

        /* Boundary faces */

        cs_lnum_t  face_id = CS_ABS(face_num) - 1;
        cs_lnum_t  vtx_start = mesh->b_face_vtx_idx[face_id] - 1;
        cs_lnum_t  vtx_end = mesh->b_face_vtx_idx[face_id+1] - 1;
        cs_lnum_t  n_vertices = vtx_end - vtx_start + 1;

        for (k = 0, j = vtx_start; j < vtx_end; j++, k++)
          face_connect[k] = mesh->b_face_vtx_lst[j];
        face_connect[n_vertices-1] = face_connect[0];

        /* Interior face which is different from the incoming face */

        indian = _where_are_you(face_num,
                                n_vertices,
                                face_connect,
                                prev_particle,
                                particle,
                                prev_p_am,
                                p_am,
                                &error);

        face_num = -face_num;

        if (error == 1 && ipass == 1) {

          _manage_error(failsafe_mode,
                        particle,
                        p_am,
                        CS_LAGR_TRACKING_ERR_DISPLACEMENT_I_FACE,
                        &n_failed_particles,
                        &failed_particle_weight);

          move_particle = CS_LAGR_PART_MOVE_OFF;
          particle_state = CS_LAGR_PART_ERR;

        }
        else if (indian == 1) { /* Particle moves to the neighbor cell
                                   through the current face "face_num" */

          /* particle / boundary condition interaction
             1 - modify particle cur_cell_num : 0 or boundary_cell_num
             2 -

             P -->  *         *  <-- Q
             \       /
             \     /
             \   /
             \ /
             ------------------      boundary condition
             K

             3 - move_particle = 0: end of particle tracking
             move_particle = 1: continue particle tracking

          */
          particle_state
            = _bdy_treatment(_prev_particle_set,
                             _particle_set,
                             prev_particle,
                             particle,
                             face_num,
                             boundary_stat,
                             bdy_conditions->b_face_zone_num[face_num-1]-1,
                             failsafe_mode,
                             &move_particle,
                             &n_failed_particles,
                             &failed_particle_weight,
                             iensi3,
                             inbr,
                             inbrbd,
                             iflm,
                             iflmbd,
                             iang,
                             iangbd,
                             ivit,
                             ivitbd,
                             iencnb,
                             iencma,
                             iencdi,
                             iencck,
                             iencnbbd,
                             iencmabd,
                             iencdibd,
                             iencckbd,
                             inclg,
                             iscovc,
                             nusbor,
                             iusb,
                             energt,
                             tprenc,
                             visref,
                             enc1,
                             enc2,
                             tkelvi);

          if (scheme_order == 2)
            cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_SWITCH_ORDER_1,
                                      CS_LAGR_SWITCH_ON);

          if (   move_particle != CS_LAGR_PART_MOVE_ON
              && move_particle != CS_LAGR_PART_MOVE_OFF)
            bft_error(__FILE__, __LINE__, 0,
                      _(" Incoherent value for move_particle = %d."
                        " Value must be 0 or 1.  \n"), move_particle);

          new_face_num = -face_num; /* To be sure that it's a boundary face */

          break;

        } /* End if indian == 1 */

        else if (indian == -1) {

          particle_state = CS_LAGR_PART_TREATED;
          move_particle = CS_LAGR_PART_MOVE_OFF;

        }

      } /* End if face_num < 0 (boundary face) && face_num != old_face_num */

    } /* End of loop on faces
         Cell -> face connect && error == 1 */

  } /* End of while : local displacement */

  cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_LAST_FACE_NUM,
                            new_face_num);

  /* Return pointers */

  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Exchange counters on the number of particles to send and to receive
 *
 * parameters:
 *  halo        <--  pointer to a cs_halo_t structure
 *  lag_halo    <--  pointer to a cs_lagr_halo_t structure
 *----------------------------------------------------------------------------*/

static void
_exchange_counter(const cs_halo_t  *halo,
                  cs_lagr_halo_t   *lag_halo)
{
  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    int  rank;
    int  request_count = 0;
    const int  local_rank = cs_glob_rank_id;

    /* Receive data from distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      if (halo->c_domain_rank[rank] != local_rank)
        MPI_Irecv(&(lag_halo->recv_count[rank]),
                  1,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank],
                  halo->c_domain_rank[rank],
                  cs_glob_mpi_comm,
                  &(lag_halo->request[request_count++]));
      else
        local_rank_id = rank;

    }

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(cs_glob_mpi_comm);

    /* Send data to distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank] != local_rank)
        MPI_Isend(&(lag_halo->send_count[rank]),
                  1,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank],
                  local_rank,
                  cs_glob_mpi_comm,
                  &(lag_halo->request[request_count++]));

    }

    /* Wait for all exchanges */

    MPI_Waitall(request_count, lag_halo->request, lag_halo->status);

  }
#endif /* defined(HAVE_MPI) */

  /* Copy local values in case of periodicity */

  if (halo->n_transforms > 0)
    if (local_rank_id > -1)
      lag_halo->recv_count[local_rank_id] = lag_halo->send_count[local_rank_id];

}

/*----------------------------------------------------------------------------
 * Exchange counters on the number of particles to send and to receive
 *
 * parameters:
 *  halo     <-- pointer to a cs_halo_t structure
 *  lag_halo <-> pointer to a cs_lagr_halo_t structure
 *----------------------------------------------------------------------------*/

static void
_exchange_particles(const cs_halo_t  *halo,
                    cs_lagr_halo_t   *lag_halo)
{
  cs_lnum_t  shift;

  void  *recv_buf = NULL, *send_buf = NULL;

  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

  const size_t tot_extents = lag_halo->cur_extents + lag_halo->prev_extents;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    int  rank;
    int  request_count = 0;
    const int  local_rank = cs_glob_rank_id;

    /* Receive data from distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      shift = lag_halo->recv_shift[rank];

      if (lag_halo->recv_count[rank] == 0)
        recv_buf = NULL;
      else
        recv_buf = lag_halo->recv_buf + tot_extents*shift;

      if (halo->c_domain_rank[rank] != local_rank)
        MPI_Irecv(recv_buf,
                  lag_halo->recv_count[rank],
                  _CS_MPI_PARTICLE,
                  halo->c_domain_rank[rank],
                  halo->c_domain_rank[rank],
                  cs_glob_mpi_comm,
                  &(lag_halo->request[request_count++]));
      else
        local_rank_id = rank;

    }

    /* We wait for posting all receives
       (often recommended in the past, apparently not anymore) */

#if 0
    MPI_Barrier(cs_glob_mpi_comm)
#endif

    /* Send data to distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank] != local_rank) {

        shift = lag_halo->send_shift[rank];
        if (lag_halo->send_count[rank] == 0)
          send_buf = NULL;
        else
          send_buf = lag_halo->send_buf + tot_extents*shift;

        MPI_Isend(send_buf,
                  lag_halo->send_count[rank],
                  _CS_MPI_PARTICLE,
                  halo->c_domain_rank[rank],
                  local_rank,
                  cs_glob_mpi_comm,
                  &(lag_halo->request[request_count++]));

      }

    }

    /* Wait for all exchanges */

    MPI_Waitall(request_count, lag_halo->request, lag_halo->status);

  }
#endif /* defined(HAVE_MPI) */

  /* Copy local values in case of periodicity */

  if (halo->n_transforms > 0) {
    if (local_rank_id > -1) {
      cs_lnum_t  recv_shift = lag_halo->recv_shift[local_rank_id];
      cs_lnum_t  send_shift = lag_halo->send_shift[local_rank_id];

      assert(   lag_halo->recv_count[local_rank_id]
             == lag_halo->send_count[local_rank_id]);

      for (cs_lnum_t i = 0; i < lag_halo->send_count[local_rank_id]; i++)
        memcpy(lag_halo->recv_buf + tot_extents*(recv_shift + i),
               lag_halo->send_buf + tot_extents*(send_shift + i),
               tot_extents);
    }
  }
}

/*----------------------------------------------------------------------------
 * Update cs_lagr_particle_set structures after a parallel and/or periodic
 * exchange between communicating ranks.
 *
 * parameters:
 *   n_recv_particles <-- number of particles received during the exchange
 *   lag_halo         <-- pointer to a cs_lagr_halo_t structure
 *   cur_set          <-- set of current particle data to update
 *   prev_set         <-- set of previous particle data to update
 *----------------------------------------------------------------------------*/

static void
_update_particle_sets(cs_lnum_t                n_recv_particles,
                      cs_lagr_halo_t          *lag_halo,
                      cs_lagr_particle_set_t  *cur_set,
                      cs_lagr_particle_set_t  *prev_set)
{
  cs_lnum_t  i, new_id;
  const size_t cur_extents = cur_set->p_am->extents;
  const size_t prev_extents = cur_set->p_am->extents;
  const size_t tot_extents = cur_extents + prev_extents;
  const unsigned char *recv_buf = lag_halo->recv_buf;

  assert(cur_extents == lag_halo->cur_extents);
  assert(prev_extents == lag_halo->prev_extents);

  for (i = 0; i < n_recv_particles; i++) {

    const unsigned char *new_part = recv_buf + tot_extents*i;

    new_id = cur_set->first_free_id;

    /* Note: used and free ids for current and previous data are the same */

    if (cur_set->used_id[new_id].next_id != -1) {
      cur_set->first_free_id = cur_set->used_id[new_id].next_id;
    }
    else {
      cur_set->first_free_id = new_id + 1 ;
      cur_set->used_id[cur_set->first_free_id].next_id = -1;
    }

    /* Add new_particle at the beginning of the "used list"
       Update first_used_id */

    if (cur_set->first_used_id != -1) {
      cur_set->used_id[cur_set->first_used_id].prev_id = new_id;
    }

    memcpy(cur_set->p_buffer + cur_extents*new_id,
           new_part,
           cur_extents);
    memcpy(prev_set->p_buffer + prev_extents*new_id,
           new_part + cur_extents,
           prev_extents);

    cur_set->used_id[new_id].prev_id = -1;
    cur_set->used_id[new_id].next_id = cur_set->first_used_id;

    cur_set->first_used_id = new_id;
  }
}

/*----------------------------------------------------------------------------
 * Synchronize particle displacement over the ranks.
 *
 * parameters:
 *  lag_halo  <--  pointer to a cs_lagr_halo_t structure
 *  set       <--  set of particles to update
 *----------------------------------------------------------------------------*/

static void
_sync_particle_sets(cs_lagr_halo_t           *lag_halo,
                    cs_interface_set_t       *face_ifs,
                    cs_lagr_particle_set_t   *prev_set,
                    cs_lagr_particle_set_t   *cur_set)
{
  cs_lnum_t  i, j, k, tr_id, rank, shift, ghost_id;
  cs_real_t matrix[3][4];

  cs_lnum_t  n_recv_particles = 0;

  const size_t cur_extents = cur_set->p_am->extents;
  const size_t prev_extents = cur_set->p_am->extents;
  const size_t tot_extents = cur_extents + prev_extents;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const fvm_periodicity_t *periodicity = mesh->periodicity;

  for (i = 0; i < halo->n_c_domains; i++){
    n_recv_particles += lag_halo->recv_count[i];
    lag_halo->send_count[i] = 0;
  }

  /* Fill send_buf for current and previous particle sets */

  for (i = 0, j = cur_set->first_used_id; i < cur_set->n_particles; i++) {

    cs_lnum_t next_id = cur_set->used_id[j].next_id;

    if (   cs_lagr_particles_get_lnum(cur_set, j, CS_LAGR_STATE)
        == CS_LAGR_PART_TO_SYNC) {

      ghost_id =   cs_lagr_particles_get_lnum(cur_set, j, CS_LAGR_CUR_CELL_NUM)
                 - halo->n_local_elts - 1;
      rank = lag_halo->rank[ghost_id];
      tr_id = lag_halo->transform_id[ghost_id];

      cs_lagr_particles_set_lnum(cur_set, j, CS_LAGR_CUR_CELL_NUM,
                                 lag_halo->dist_cell_num[ghost_id]);

      shift = lag_halo->send_shift[rank] + lag_halo->send_count[rank];

      /* Update if needed last_face_num */

      if (tr_id >= 0) { /* Same initialization as in previous algorithm */

        cs_lagr_particles_set_lnum(cur_set, j, CS_LAGR_LAST_FACE_NUM, 0);

      }

      else {

        if (cs_glob_n_ranks > 1) {

          assert(face_ifs != NULL);

          int  distant_rank, n_entities, id;
          const int* local_num, * dist_num;

          const int search_rank = halo->c_domain_rank[rank];
          const cs_interface_t  *interface = NULL;
          const int  n_interfaces = cs_interface_set_size(face_ifs);

          for (k = 0; k < n_interfaces; k++) {

            interface = cs_interface_set_get(face_ifs,k);

            distant_rank = cs_interface_rank(interface);

            if (distant_rank == search_rank)
              break;

          }

          if (k == n_interfaces) {
            bft_error(__FILE__, __LINE__, 0,
                      _(" Cannot find the relative distant rank.\n"));

          }
          else {

            n_entities = cs_interface_size(interface);
            local_num = cs_interface_get_elt_ids(interface);

            cs_lnum_t cur_part_last_face_num
              = cs_lagr_particles_get_lnum(cur_set, j, CS_LAGR_LAST_FACE_NUM);

            id = cs_search_binary(n_entities,
                                  cur_part_last_face_num - 1,
                                  local_num);

            if (id == -1)
              bft_error(__FILE__, __LINE__, 0,
                        _(" Cannot find the relative distant face num.\n"));

            dist_num = cs_interface_get_match_ids(interface);

            cs_lagr_particles_set_lnum(cur_set, j, CS_LAGR_LAST_FACE_NUM,
                                       dist_num[id] + 1);

          }

        }
      }

      /* Periodicity treatment */

      if (tr_id >= 0) {

        /* Transform coordinates */

        fvm_periodicity_type_t  perio_type
          = fvm_periodicity_get_type(periodicity, tr_id);

        int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, tr_id);

        fvm_periodicity_get_matrix(periodicity, rev_id, matrix);

        /* Apply transformation to the coordinates in any case */

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              cs_lagr_particles_attr(cur_set, j,
                                                     CS_LAGR_COORDS));

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              cs_lagr_particles_attr(prev_set, j,
                                                     CS_LAGR_COORDS));

        /* Apply rotation to velocity vectors in case of rotation */

        if (perio_type >= FVM_PERIODICITY_ROTATION) {

          /* Rotation of the velocity */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(cur_set, j,
                                                        CS_LAGR_VELOCITY));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(prev_set, j,
                                                        CS_LAGR_VELOCITY));

          /* Rotation of the velocity seen */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(cur_set, j,
                                                        CS_LAGR_VELOCITY_SEEN));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(prev_set, j,
                                                        CS_LAGR_VELOCITY_SEEN));

        } /* Specific treatment in case of rotation for the velocities */

      } /* End of periodicity treatment */

      memcpy(lag_halo->send_buf + tot_extents*shift,
             cur_set->p_buffer + cur_extents*j,
             cur_extents);

      memcpy(lag_halo->send_buf + tot_extents*shift + cur_extents,
             prev_set->p_buffer + prev_extents*j,
             prev_extents);

      lag_halo->send_count[rank] += 1;

      /* Pick out the particle from the particle_sets */

      _remove_particle(cur_set, j);

    } /* TO_SYNC */

    j = next_id;

  } /* End of loop on particles */

  /* Exchange particles, then update sets */

  _exchange_particles(halo, lag_halo);

  _update_particle_sets(n_recv_particles, lag_halo, cur_set, prev_set);
}

/*----------------------------------------------------------------------------
 * Move particles between communicating ranks and
 * update particle set structures
 *----------------------------------------------------------------------------*/

static void
_lagr_halo_sync(void)
{
  cs_lnum_t  i, j, ghost_id;
  cs_lnum_t  delta_particles;

  cs_lnum_t  n_recv_particles = 0, n_send_particles = 0;
  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_halo_t  *lag_halo = builder->halo;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;

  assert(set->n_particles == prev_set->n_particles);

  /* Check this because this assumption is done. */

  /* Initialization */

  for (i = 0; i < halo->n_c_domains; i++) {
    lag_halo->send_count[i] = 0;
    lag_halo->recv_count[i] = 0;
  }

  /* Loop on particles to count number of particles to send on each rank */

  for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

    cs_lnum_t next_id = set->used_id[j].next_id;

    if (   cs_lagr_particles_get_lnum(set, j, CS_LAGR_STATE)
        == CS_LAGR_PART_TO_SYNC) {

      ghost_id =   cs_lagr_particles_get_lnum(set, j, CS_LAGR_CUR_CELL_NUM)
                 - mesh->n_cells - 1;

      assert(ghost_id >= 0);
      lag_halo->send_count[lag_halo->rank[ghost_id]] += 1;

    }

    j = next_id;

  } /* End of loop on particles */

  /* Exchange counters */

  _exchange_counter(halo, lag_halo);

  for (i = 0; i < halo->n_c_domains; i++) {
    n_recv_particles += lag_halo->recv_count[i];
    n_send_particles += lag_halo->send_count[i];
  }

  delta_particles = n_recv_particles - n_send_particles;

  lag_halo->send_shift[0] = 0;
  lag_halo->recv_shift[0] = 0;

  for (i = 1; i < halo->n_c_domains; i++) {

    lag_halo->send_shift[i] =  lag_halo->send_shift[i-1]
                             + lag_halo->send_count[i-1];

    lag_halo->recv_shift[i] =  lag_halo->recv_shift[i-1]
                             + lag_halo->recv_count[i-1];

  }

  /* Resize halo only if needed */

  _resize_lagr_halo(lag_halo, n_send_particles, n_recv_particles);

  /* Get the updated particle set after synchronization */

  _sync_particle_sets(lag_halo, builder->face_ifs, prev_set, set);

  set->n_particles += delta_particles;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id == 1) {
    bft_printf(" delta part = %d\n", delta_particles);
    bft_printf(" set -> npart = %d\n", set->n_particles);
  }
#endif

  prev_set->n_particles += delta_particles;

  if (delta_particles > set->n_particles_max - set->n_particles)
    bft_error(__FILE__, __LINE__, 0,
              _(" Not enough memory to receive particles.\n"
                " We can still receive %d particles and"
                " we have to receive %d additional particles.\n"
                " Check n_particles_max (%d).\n"),
              set->n_particles_max - set->n_particles,
              delta_particles, set->n_particles_max);

  /* TODO: Do a resize to fit to the new size of the particle set */
}

/*----------------------------------------------------------------------------
 * Update C particle structures data with Fortran data.
 *
 * The associated metadata should have been reinitialized or updated
 * separately.
 *
 * parameters:
 *   nbpmax <-- n_particles max.
 *   ...
 *----------------------------------------------------------------------------*/

static void
_update_c_from_fortran(const cs_lnum_t   *nbpmax,
                       const cs_real_t    ettp[],
                       const cs_real_t    ettpa[],
                       const cs_lnum_t    itepa[],
                       const cs_real_t    tepa[],
                       const cs_lnum_t    ibord[],
                       const cs_lnum_t    indep[])
{
  cs_lnum_t  i, id;

  cs_lagr_particle_set_t  *cur = _particle_set;
  cs_lagr_particle_set_t  *prv = _prev_particle_set;

  cs_lagr_param_t *lagr_param = &cs_glob_lagr_param;

  const cs_lagr_attribute_map_t  *am = cur->p_am;

  assert(*nbpmax == cur->n_particles_max); /* Up to now, we don't manage
                                              a mofification of nbpmax */

  /* When we receive particles from FORTRAN, we keep a compact
     storage of particles */

  if (cur->n_particles > 0)
    cur->first_used_id = 0;
  else
    cur->first_used_id = -1;

  cur->first_free_id = cur->n_particles;
  cur->used_id[cur->first_free_id].next_id = -1;

  /* Fill set and prv structures */

  for (i = 0; i < cur->n_particles; i++) {

    if (i > 0)
      cur->used_id[i].prev_id = i-1;
    else  /* Not defined */
      cur->used_id[i].prev_id = -1;

    if (i < cur->n_particles - 1)
      cur->used_id[i].next_id = i+1;
    else  /* Not defined */
      cur->used_id[i].next_id = -1;

    cs_lagr_particles_set_lnum(cur, i, CS_LAGR_CUR_CELL_NUM,
                               itepa[i + _jisor * (*nbpmax)]);
    cs_lagr_particles_set_lnum(prv, i, CS_LAGR_CUR_CELL_NUM,
                               indep[i]);

    cs_lagr_particles_set_real(cur, i, CS_LAGR_RANDOM_VALUE,
                               tepa[i + _jrval * (*nbpmax)]);
    cs_lagr_particles_set_real(prv, i, CS_LAGR_RANDOM_VALUE,
                               tepa[i + _jrval * (*nbpmax)]);

    cs_lnum_t cur_part_cur_cell_num
      = cs_lagr_particles_get_lnum(cur, i, CS_LAGR_CUR_CELL_NUM);
    if (cur_part_cur_cell_num < 0)
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_STATE, CS_LAGR_PART_STICKED);
    else if (cur_part_cur_cell_num == 0)
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_STATE, CS_LAGR_PART_TO_DELETE);
    else
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_STATE, CS_LAGR_PART_TO_SYNC);

    cs_lagr_particles_set_lnum(cur, i, CS_LAGR_LAST_FACE_NUM, 0);
    cs_lagr_particles_set_lnum(prv, i, CS_LAGR_LAST_FACE_NUM, 0);

    cs_lagr_particles_set_lnum(cur, i, CS_LAGR_SWITCH_ORDER_1, ibord[i]);
    cs_lagr_particles_set_lnum(prv, i, CS_LAGR_SWITCH_ORDER_1, ibord[i]);

    id = _jrpoi * (*nbpmax) + i;
    cs_lagr_particles_set_real(cur, i, CS_LAGR_STAT_WEIGHT, tepa[id]);
    cs_lagr_particles_set_real(prv, i, CS_LAGR_STAT_WEIGHT, tepa[id]);

    id = _jrtsp * (*nbpmax) + i;
    cs_lagr_particles_set_real(cur, i, CS_LAGR_RESIDENCE_TIME, tepa[id]);
    cs_lagr_particles_set_real(prv, i, CS_LAGR_RESIDENCE_TIME, tepa[id]);

    id = _jmp * (*nbpmax) + i;
    cs_lagr_particles_set_real(cur, i, CS_LAGR_MASS, ettp[id]);
    cs_lagr_particles_set_real(prv, i, CS_LAGR_MASS, ettpa[id]);

    id = _jdp * (*nbpmax) + i;
    cs_lagr_particles_set_real(cur, i, CS_LAGR_DIAMETER, ettp[id]);
    cs_lagr_particles_set_real(prv, i, CS_LAGR_DIAMETER, ettpa[id]);

    /* Coordinates of the particle */

    cs_real_t *cur_part_coord = cs_lagr_particles_attr(cur, i, CS_LAGR_COORDS);
    cs_real_t *prv_part_coord = cs_lagr_particles_attr(prv, i, CS_LAGR_COORDS);

    id = _jxp * (*nbpmax) + i;
    cur_part_coord[0] = ettp[id];
    prv_part_coord[0] = ettpa[id];

    id = _jyp * (*nbpmax) + i;
    cur_part_coord[1] = ettp[id];
    prv_part_coord[1] = ettpa[id];

    id = _jzp * (*nbpmax) + i;
    cur_part_coord[2] = ettp[id];
    prv_part_coord[2] = ettpa[id];

    /* Velocity of the particle */

    cs_real_t *cur_part_velocity
      = cs_lagr_particles_attr(cur, i, CS_LAGR_VELOCITY);
    cs_real_t *prv_part_velocity
      = cs_lagr_particles_attr(prv, i, CS_LAGR_VELOCITY);

    id = _jup * (*nbpmax) + i;
    cur_part_velocity[0] = ettp[id];
    prv_part_velocity[0] = ettpa[id];

    id = _jvp * (*nbpmax) + i;
    cur_part_velocity[1] = ettp[id];
    prv_part_velocity[1] = ettpa[id];

    id = _jwp * (*nbpmax) + i;
    cur_part_velocity[2] = ettp[id];
    prv_part_velocity[2] = ettpa[id];

    /* Velocity seen by the fluid */

    cs_real_t *cur_part_velocity_seen
      = cs_lagr_particles_attr(cur, i, CS_LAGR_VELOCITY_SEEN);
    cs_real_t *prv_part_velocity_seen
      = cs_lagr_particles_attr(prv, i, CS_LAGR_VELOCITY_SEEN);

    id = _juf * (*nbpmax) + i;
    cur_part_velocity_seen[0] = ettp[id];
    prv_part_velocity_seen[0] = ettpa[id];

    id = _jvf * (*nbpmax) + i;
    cur_part_velocity_seen[1] = ettp[id];
    prv_part_velocity_seen[1] = ettpa[id];

    id = _jwf * (*nbpmax) + i;
    cur_part_velocity_seen[2] = ettp[id];
    prv_part_velocity_seen[2] = ettpa[id];

    id = _jtaux * (*nbpmax) + i;
    cs_lagr_particles_set_real(cur, i, CS_LAGR_TAUP_AUX, ettp[id]);

    /* Data needed if the deposition model is activated */
    if (lagr_param->deposition > 0) {

      id = _jdepo  * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_DEPOSITION_FLAG, itepa[id]);

      id = _jryplu * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_YPLUS, tepa[id]);

      id = _jrinpf * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_INTERF, tepa[id]);

      id = _jdfac  * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_NEIGHBOR_FACE_ID,
                                 itepa[id] - 1);

      id = _jimark  * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_MARKO_VALUE, itepa[id]);

    } else {

      if (am->size[CS_LAGR_DEPOSITION_FLAG] > 0)
        cs_lagr_particles_set_lnum(cur, i, CS_LAGR_DEPOSITION_FLAG, 0);

    }

    /* Data needed if the resuspension model is activated */

    if (am->size[CS_LAGR_N_LARGE_ASPERITIES] > 0) {
      id = _jnbasg  * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_N_LARGE_ASPERITIES,
                                 itepa[id]);
    }

    if (am->size[CS_LAGR_N_SMALL_ASPERITIES] > 0) {
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_N_SMALL_ASPERITIES, 0);
      id = _jnbasp  * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_N_SMALL_ASPERITIES,
                                 itepa[id]);
    }

    if (am->size[CS_LAGR_ADHESION_FORCE] > 0) {
      id = _jfadh  * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_ADHESION_FORCE, tepa[id]);
    }

    if (am->size[CS_LAGR_DISPLACEMENT_NORM] > 0) {
      id = _jndisp  * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_DISPLACEMENT_NORM, tepa[id]);
    }

    if (am->size[CS_LAGR_ADHESION_TORQUE] > 0) {
      id = _jmfadh  * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_ADHESION_TORQUE, tepa[id]);
    }

    if (lagr_param->resuspension > 0) {
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_RANK_FLAG, 0);
    }

    if (am->size[CS_LAGR_TEMPERATURE] > 0) {

      cs_real_t *cur_part_temp
        = cs_lagr_particles_attr(cur, i, CS_LAGR_TEMPERATURE);
      cs_real_t *prv_part_temp
        = cs_lagr_particles_attr(prv, i, CS_LAGR_TEMPERATURE);

      for (cs_lnum_t j = 0; j < am->count[CS_LAGR_TEMPERATURE]; j++) {
        if (_jthp[j] > -1) {
          id = _jthp[j] * (*nbpmax) + i;
          cur_part_temp[j] = ettp[id];
          prv_part_temp[j] = ettpa[id];
        }
        else {
          cur_part_temp[j] = 0.0;
          prv_part_temp[j] = 0.0;
        }
      }

    }

    if (am->size[CS_LAGR_FLUID_TEMPERATURE] > 0) {
      id = _jtf * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_FLUID_TEMPERATURE, ettp[id]);
      cs_lagr_particles_set_real(prv, i, CS_LAGR_FLUID_TEMPERATURE, ettpa[id]);
    }

    /* Data needed if the coal combustion is activated */

    if (lagr_param->physic_mode ==2) {

      /* ettp and ettpa arrays */

      id = _jmwat * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_WATER_MASS, ettp[id]);
      cs_lagr_particles_set_real(prv, i, CS_LAGR_WATER_MASS, ettpa[id]);

      if (am->count[CS_LAGR_TEMPERATURE] > 0) {

        cs_real_t *cur_part_coal_mass
          = cs_lagr_particles_attr(cur, i, CS_LAGR_COAL_MASS);
        cs_real_t *prv_part_coal_mass
          = cs_lagr_particles_attr(prv, i, CS_LAGR_COAL_MASS);

        cs_real_t *cur_part_coke_mass
          = cs_lagr_particles_attr(cur, i, CS_LAGR_COKE_MASS);
        cs_real_t *prv_part_coke_mass
          = cs_lagr_particles_attr(prv, i, CS_LAGR_COKE_MASS);

        for (cs_lnum_t j = 0; j < am->count[CS_LAGR_TEMPERATURE]; j++) {
          id = _jmch[j] * (*nbpmax) + i;
          cur_part_coal_mass[j] = ettp[id];
          prv_part_coal_mass[j] = ettpa[id];

          id = _jmck[j] * (*nbpmax) + i;
          cur_part_coke_mass[j] = ettp[id];
          prv_part_coke_mass[j] = ettpa[id];
        }

      }

      id = _jcp * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_CP, ettp[id]);
      cs_lagr_particles_set_real(prv, i, CS_LAGR_CP, ettpa[id]);

      /* tepa array */

      id = _jrdck * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_SHRINKING_DIAMETER, tepa[id]);
      /*  FIXME for post-processing by trajectory purpose. To have better
          results shrinking_diam should be transfered to ETTPA/ETTP arrays*/
      cs_lagr_particles_set_real(prv, i, CS_LAGR_SHRINKING_DIAMETER, tepa[id]);

      id = _jrd0p * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_INITIAL_DIAMETER, tepa[id]);

      if  (am->count[CS_LAGR_TEMPERATURE] > 0) {
        cs_real_t *cur_part_coal_density
          = cs_lagr_particles_attr(cur, i, CS_LAGR_COAL_DENSITY);
        for (cs_lnum_t j = 0; j < am->count[CS_LAGR_TEMPERATURE]; j++) {
          id = _jrhock[j] * (*nbpmax) + i;
          cur_part_coal_density[j] = tepa[id];
        }
      }

      /* itepa array */

      id = _jinch * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_COAL_NUM, itepa[id]);

    } /* iphyla == 2 */

    /* Values based on other physical model options */

    if (am->size[CS_LAGR_EMISSIVITY] > 0) {
      id = _jreps * (*nbpmax) + i;
      cs_lagr_particles_set_real(cur, i, CS_LAGR_EMISSIVITY, ettp[id]);
      cs_lagr_particles_set_real(cur, i, CS_LAGR_EMISSIVITY, ettpa[id]);
    }

    /* Statistical class */

    if (am->size[CS_LAGR_STAT_CLASS] > 0) {
      id = _jclst * (*nbpmax) + i;
      cs_lagr_particles_set_lnum(cur, i, CS_LAGR_STAT_CLASS, itepa[id]);
    }

    /* User variables */

    if (am->count[CS_LAGR_USER] > 0) {
      cs_real_t *usr = cs_lagr_particles_attr(cur, i, CS_LAGR_USER);
      cs_real_t *usra = cs_lagr_particles_attr(prv, i, CS_LAGR_USER);
      for (cs_lnum_t j = 0; j < cur->p_am->count[CS_LAGR_USER]; j++) {
        id = _jvls[j] * (*nbpmax) + i;
        usr[j] = ettp[id];
        usra[j] = ettpa[id];
      }
    }

  }

}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate cs_lagr_particle_set_t structure and initialize useful buffers
 * and indexes
 *
 * parameters:
 *   n_particles_max <--  local max. number of particles
 *   iphyla          <--  kind of physics used for the lagrangian approach
 *   nvls            <--  number of user-defined variables
 *   nbclst          <--  number of stat. class to study sub-set of particles
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagbeg, LAGBEG)(const cs_int_t    *n_particles_max,
                          const cs_int_t    *nlayer,
                          const cs_int_t    *iphyla,
                          const cs_int_t    *idepst,
                          const cs_int_t    *irough,
                          const cs_int_t    *ireent,
                          const cs_int_t    *iclogst,
                          const cs_int_t    *nvls,
                          const cs_int_t    *nbclst,
                          cs_lnum_t          icocel[],
                          cs_lnum_t          itycel[],
                          const cs_lnum_t   *jisor,
                          const cs_lnum_t   *jrval,
                          const cs_lnum_t   *jrpoi,
                          const cs_lnum_t   *jrtsp,
                          const cs_lnum_t   *jdp,
                          const cs_lnum_t   *jmp,
                          const cs_lnum_t   *jxp,
                          const cs_lnum_t   *jyp,
                          const cs_lnum_t   *jzp,
                          const cs_lnum_t   *jup,
                          const cs_lnum_t   *jvp,
                          const cs_lnum_t   *jwp,
                          const cs_lnum_t   *juf,
                          const cs_lnum_t   *jvf,
                          const cs_lnum_t   *jwf,
                          const cs_lnum_t   *jtaux,
                          const cs_lnum_t   *jryplu,
                          const cs_lnum_t   *jrinpf,
                          const cs_lnum_t   *jdfac,
                          const cs_lnum_t   *jimark,
                          const cs_lnum_t   *jtp,
                          const cs_lnum_t    jhp[],
                          const cs_lnum_t   *jtf,
                          const cs_lnum_t   *jmwat,
                          const cs_lnum_t    jmch[],
                          const cs_lnum_t    jmck[],
                          const cs_lnum_t   *jcp,
                          const cs_lnum_t   *jrdck,
                          const cs_lnum_t   *jrd0p,
                          const cs_lnum_t   *jinch,
                          const cs_lnum_t    jrhock[],
                          const cs_lnum_t   *jreps,
                          const cs_lnum_t   *jdepo,
                          const cs_lnum_t   *jnbasg,
                          const cs_lnum_t   *jnbasp,
                          const cs_lnum_t   *jfadh,
                          const cs_lnum_t   *jmfadh,
                          const cs_lnum_t   *jndisp,
                          const cs_lnum_t   *jclst,
                          const cs_lnum_t   *jvls)
{
  cs_lnum_t  i;
  cs_mesh_t  *mesh = cs_glob_mesh;

  /* Initialize global parameter relative to the lagrangian module */

  cs_glob_lagr_param.physic_mode = *iphyla;

  if (cs_glob_lagr_param.physic_mode == 2 && *nlayer != CS_LAGR_N_LAYERS)
    bft_error(__FILE__, __LINE__, 0,"%s %i \n %s %i \n %s \n %s %i\n %s\n",
              "Multi-layer computation with iphyla = ",
              cs_glob_lagr_param.physic_mode,
              "The number of layers (defined in lagpar) is nlayer = ",
              *nlayer,
              "For the particle structure defined in cs_lagr_tracking.h,",
              "the size of temp, coal_mass, coke_mass, coal_density arrays "
              "is CS_LAGR_N_LAYERS = ",
              CS_LAGR_N_LAYERS,
              "Please correct this difference.");

  if (cs_glob_lagr_param.physic_mode == 2)
    cs_glob_lagr_param.n_temperature_layers = CS_LAGR_N_LAYERS;
  else
    cs_glob_lagr_param.n_temperature_layers = 1;

  cs_glob_lagr_param.deposition = *idepst;
  cs_glob_lagr_param.rough = *irough;
  cs_glob_lagr_param.resuspension = *ireent;
  cs_glob_lagr_param.clogging = *iclogst;

  cs_glob_lagr_param.n_user_variables = *nvls;
  cs_glob_lagr_param.n_stat_classes = *nbclst;

  /* Allocate _jthp _jmch _jmck _jrhock */
  BFT_MALLOC(_jthp   , CS_LAGR_N_LAYERS, int);
  BFT_MALLOC(_jmch   , CS_LAGR_N_LAYERS, int);
  BFT_MALLOC(_jmck   , CS_LAGR_N_LAYERS, int);
  BFT_MALLOC(_jrhock , CS_LAGR_N_LAYERS, int);
  for (i = 0; i < CS_LAGR_N_LAYERS ; i++) {
    _jthp[i]=-1;
    _jmch[i]=-1;
    _jmck[i]=-1;
    _jrhock[i]=-1;
  }

  /* Set indexes */

  _jisor = *jisor - 1;
  _jrval = *jrval - 1;
  _jrpoi = *jrpoi - 1;
  _jrtsp = *jrtsp - 1;

  _jmp = *jmp - 1;
  _jdp = *jdp - 1;

  _jxp = *jxp - 1;
  _jyp = *jyp - 1;
  _jzp = *jzp - 1;

  _jup = *jup - 1;
  _jvp = *jvp - 1;
  _jwp = *jwp - 1;

  _juf = *juf - 1;
  _jvf = *jvf - 1;
  _jwf = *jwf - 1;

  _jtaux = *jtaux - 1;
  _jdepo = *jdepo - 1;

  _jryplu = *jryplu - 1;
  _jrinpf = *jrinpf - 1;

  _jdfac = *jdfac - 1;
  _jimark = *jimark - 1;

  _jnbasg = *jnbasg - 1;
  _jnbasp = *jnbasp - 1;

  _jfadh = *jfadh - 1;
  _jmfadh = *jmfadh - 1;
  _jndisp = *jndisp - 1;

  if (*jtp > 0) {
    _jthp[0] = *jtp - 1;
  }
  else{
    for (i = 0; i < CS_LAGR_N_LAYERS; i++)
      _jthp[i]=jhp[i]-1;
  }

  _jtf = *jtf - 1;
  _jmwat = *jmwat - 1;
  for (i = 0; i < CS_LAGR_N_LAYERS; i++){
    _jmch[i] = jmch[i] - 1;
    _jmck[i] = jmck[i] - 1;
  }
  _jcp = *jcp - 1;

  _jrdck = *jrdck - 1;
  _jrd0p = *jrd0p - 1;
  for (i = 0; i < CS_LAGR_N_LAYERS; i++){
    _jrhock[i] = jrhock[i] - 1;
  }
  _jinch = *jinch - 1;

  _jreps = *jreps - 1;

  if (*nbclst > 0)
    _jclst = *jclst - 1;

  if (*nvls > 0) {
    BFT_MALLOC(_jvls, *nvls, int);
    for (i = 0; i < *nvls; i++){
      _jvls[i] = jvls[i] - 1;
    }
  }

  /* Build mappings
     (in the future, they should be created first, then marked,
     then built) */

  _p_attr_map = _create_attr_map();
  _p_attr_map_prev = _create_attr_map();

  /* Initialize particle set : prev and current */

  _particle_set = _create_particle_set(0, *n_particles_max, _p_attr_map);
  _prev_particle_set = _create_particle_set(1, *n_particles_max, _p_attr_map_prev);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER CREATION\n");
  _dump_particle_set(_particle_set);
  bft_printf("\n PREV PARTICLE SET AFTER CREATION\n");
  _dump_particle_set(_prev_particle_set);
#endif

  /* Initialization */

  for (i = 0; i < _particle_set->n_particles_max; i++) {

    cs_lagr_particles_set_lnum(_particle_set, i, CS_LAGR_SWITCH_ORDER_1,
                               CS_LAGR_SWITCH_OFF);
    cs_lagr_particles_set_lnum(_prev_particle_set, i, CS_LAGR_SWITCH_ORDER_1,
                               CS_LAGR_SWITCH_OFF);

    cs_lagr_particles_set_lnum(_particle_set, i, CS_LAGR_STATE,
                               CS_LAGR_PART_TO_SYNC);
    cs_lagr_particles_set_lnum(_prev_particle_set, i, CS_LAGR_STATE,
                               CS_LAGR_PART_TO_SYNC);

  }

  /* Initialize builder */

  _particle_track_builder
    = _init_track_builder(*n_particles_max,
                          _particle_set->p_am->extents,
                          _prev_particle_set->p_am->extents);

  /* Create all useful MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _CS_MPI_PARTICLE = _define_particle_datatype(_particle_set->p_am,
                                                 _prev_particle_set->p_am);
    }
#endif

  /* Saving of itycel and icocel */
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  for (i = 0; i <= mesh->n_cells; i++)
    itycel[i] = builder->cell_face_idx[i];
  for (i = 0; i < builder->cell_face_idx[mesh->n_cells] ; i++)
    icocel[i] = builder->cell_face_lst[i];
}

/*----------------------------------------------------------------------------
 * Get variables and parameters associated to each particle and keep it in
 * a new structure
 *
 * parameters:
 *   nbpmax <-- n_particles max.
 *   nbpart --> number of current particles
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (prtget, PRTGET)(const cs_lnum_t   *nbpmax,
                          const cs_lnum_t   *nbpart,
                          const cs_real_t    ettp[],
                          const cs_real_t    ettpa[],
                          const cs_lnum_t    itepa[],
                          const cs_real_t    tepa[],
                          const cs_lnum_t    ibord[],
                          const cs_lnum_t    indep[])
{
  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prv = _prev_particle_set;

  assert(*nbpmax == set->n_particles_max); /* Up to now, we don't manage
                                              a mofification of nbpmax */
  set->n_particles = *nbpart;
  prv->n_particles = *nbpart;

  set->weight = 0.0;
  prv->weight =0.0;

  set->n_part_out = 0;
  prv->n_part_out = 0;

  set->n_part_dep = 0;
  prv->n_part_dep = 0;

  set->n_part_fou = 0;
  prv->n_part_fou = 0;

  set->weight_out = 0.0;
  prv->weight_out = 0.0;

  set->weight_dep = 0.0;
  prv->weight_dep = 0.0;

  set->weight_fou = 0.0;
  prv->weight_fou = 0.0;

  set->n_failed_part = 0;
  prv->n_failed_part = 0;

  set->weight_failed = 0.0;
  prv->weight_failed = 0.0;

  _update_c_from_fortran(nbpmax,
                         ettp,
                         ettpa,
                         itepa,
                         tepa,
                         ibord,
                         indep);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER PRTGET\n");
  _dump_particle_set(set);
  bft_printf("\n PREV PARTICLE SET AFTER PRTGET\n");
  _dump_particle_set(prv);
#endif
}

/*----------------------------------------------------------------------------
 * Put variables and parameters associated to each particles into FORTRAN
 * arrays.
 *
 * parameters:
 *   nbpmax <-- n_particles max.
 *   nbpart --> number of current particles
 *   dnbpar --> particle total weight
 *   nbpout --> number of outgoing particles
 *   dnbpou --> outgoing particle total weight
 *   nbperr --> number of failed particles
 *   dnbper --> failed particles total weight
 *   nbpdep --> number of depositing particles
 *   dnbdep --> depositing particles total weight
 *   npencr --> number of fouled particles (coal)
 *   dnpenc --> fouled particles (coal) total weight
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (prtput, PRTPUT)(const cs_int_t   *nbpmax,
                          cs_int_t         *nbpart,
                          cs_real_t        *dnbpar,
                          cs_int_t         *nbpout,
                          cs_real_t        *dnbpou,
                          cs_int_t         *nbperr,
                          cs_real_t        *dnbper,
                          cs_int_t         *nbpdep,
                          cs_real_t        *dnbdep,
                          cs_int_t         *npencr,
                          cs_real_t        *dnpenc,
                          cs_real_t         ettp[],
                          cs_real_t         ettpa[],
                          cs_int_t          itepa[],
                          cs_real_t         tepa[],
                          cs_int_t          ibord[])
{
  cs_lnum_t  i, j, id, nbp;

  cs_lagr_particle_set_t  *cur = _particle_set;
  cs_lagr_particle_set_t  *prv = _prev_particle_set;

  cs_lagr_param_t *lagr_param = &cs_glob_lagr_param;

  const cs_lnum_t n_temperature_layers = cur->p_am->count[CS_LAGR_TEMPERATURE];

  assert(*nbpmax == cur->n_particles_max);

  j = cur->first_used_id;

  nbp = 0;

  for (i = 0; (i < cur-> n_particles )&(j != -1); i++) {

    nbp++;

    cs_lnum_t cur_part_state
      = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_STATE);

    if  (   cur_part_state == CS_LAGR_PART_TREATED
         || cur_part_state == CS_LAGR_PART_STICKED) {
      itepa[_jisor * (*nbpmax) + i]
        = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_CUR_CELL_NUM);
    }
    else if (cur_part_state == CS_LAGR_PART_OUT) {
      itepa[_jisor * (*nbpmax) + i] = 0;
    }
    else {
      assert(cur_part_state = CS_LAGR_PART_ERR);

      itepa[_jisor * (*nbpmax) + i] = 0;
    }

    tepa[_jrval * (*nbpmax) + i]
      = cs_lagr_particles_get_real(cur, j, CS_LAGR_RANDOM_VALUE);

    /* Data needed if the deposition model is activated */

    if (lagr_param->deposition > 0) {

      tepa[_jryplu * (*nbpmax) + i]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_YPLUS);

      tepa[_jrinpf * (*nbpmax) + i]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_INTERF);

      itepa[_jdfac * (*nbpmax) + i]
        = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_NEIGHBOR_FACE_ID) + 1;

      itepa[_jimark * (*nbpmax) + i]
        = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_MARKO_VALUE);

      itepa[_jdepo * (*nbpmax) + i]
        = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_DEPOSITION_FLAG);

    }

    /* Data needed if the resuspension model is activated */

    if (lagr_param->resuspension > 0) {

      itepa[_jnbasg * (*nbpmax) + i]
        = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_N_LARGE_ASPERITIES);

      itepa[_jnbasp * (*nbpmax) + i]
        = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_N_SMALL_ASPERITIES);

      tepa[_jfadh * (*nbpmax) + i]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_ADHESION_FORCE);

      tepa[_jmfadh * (*nbpmax) + i]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_ADHESION_TORQUE);

      tepa[_jndisp * (*nbpmax) + i]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_DISPLACEMENT_NORM);

    }

    /* Data needed if coal combustion is activated */

    if (n_temperature_layers > 0) {

      const cs_real_t *cur_part_temp
        = cs_lagr_particles_attr_const(cur, j, CS_LAGR_TEMPERATURE);
      const cs_real_t *prev_part_temp
        = cs_lagr_particles_attr_const(prv, j, CS_LAGR_TEMPERATURE);

      for (cs_lnum_t k = 0; k < n_temperature_layers; k++) {
        if (_jthp[k] > -1) {
          id = _jthp[k] * (*nbpmax) + i;
          ettp[id] = cur_part_temp[k];
          ettpa[id] = prev_part_temp[k];
        }
      }

    }

    if (lagr_param->physic_mode == 2) {

      /* ettp and ettpa arrays */

      id = _jtf * (*nbpmax) + i;
      ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_FLUID_TEMPERATURE);
      ettpa[id] = cs_lagr_particles_get_real(prv, j, CS_LAGR_FLUID_TEMPERATURE);

      id = _jmwat * (*nbpmax) + i;
      ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_WATER_MASS);
      ettpa[id] = cs_lagr_particles_get_real(prv, j, CS_LAGR_WATER_MASS);

      if (n_temperature_layers > 0) {

        const cs_real_t *cur_part_coal_mass
          = cs_lagr_particles_attr_const(cur, j, CS_LAGR_COAL_MASS);
        const cs_real_t *prev_part_coal_mass
          = cs_lagr_particles_attr_const(prv, j, CS_LAGR_COAL_MASS);

        const cs_real_t *cur_part_coke_mass
          = cs_lagr_particles_attr_const(cur, j, CS_LAGR_COKE_MASS);
        const cs_real_t *prev_part_coke_mass
          = cs_lagr_particles_attr_const(prv, j, CS_LAGR_COKE_MASS);

        for (cs_lnum_t k = 0; k < n_temperature_layers; k++) {
          id = _jmch[k] * (*nbpmax) + i;
          ettp[id] = cur_part_coal_mass[k];
          ettpa[id] = prev_part_coal_mass[k];

          id = _jmck[k] * (*nbpmax) + i;
          ettp[id] = cur_part_coke_mass[k];
          ettpa[id] = prev_part_coke_mass[k];
        }

      }

      id = _jcp * (*nbpmax) + i;
      ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_CP);
      ettpa[id] = cs_lagr_particles_get_real(prv, j, CS_LAGR_CP);

      /* tepa array */

      id = _jrdck * (*nbpmax) + i;
      tepa[id]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_SHRINKING_DIAMETER);

      id = _jrd0p * (*nbpmax) + i;
      tepa[id]
        = cs_lagr_particles_get_real(cur, j, CS_LAGR_INITIAL_DIAMETER);

      if (n_temperature_layers > 0) {
        const cs_real_t *cur_part_coal_density
          = cs_lagr_particles_attr_const(cur, j, CS_LAGR_COAL_DENSITY);
        for (cs_lnum_t k = 0; k < n_temperature_layers; k++) {
          id = _jrhock[k] * (*nbpmax) + i;
          tepa[id] = cur_part_coal_density[k];
        }
      }

      /* itepa array */

      id = _jinch * (*nbpmax) + i;
      itepa[id] = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_COAL_NUM);

    }

    ibord[i] = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_SWITCH_ORDER_1);

    id = _jrpoi * (*nbpmax) + i;
    tepa[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_STAT_WEIGHT);

    id = _jrtsp * (*nbpmax) + i;
    tepa[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_RESIDENCE_TIME);

    id = _jmp * (*nbpmax) + i;
    ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_MASS);
    ettpa[id] = cs_lagr_particles_get_real(prv, j, CS_LAGR_MASS);

    id = _jdp * (*nbpmax) + i;
    ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_DIAMETER);
    ettpa[id] = cs_lagr_particles_get_real(prv, j, CS_LAGR_DIAMETER);

    /* Coordinates of the particle */

    const cs_real_t *cur_part_coord
      = cs_lagr_particles_attr_const(cur, j, CS_LAGR_COORDS);
    const cs_real_t *prv_part_coord
      = cs_lagr_particles_attr_const(prv, j, CS_LAGR_COORDS);

    id = _jxp * (*nbpmax) + i;
    ettp[id] = cur_part_coord[0];
    ettpa[id] = prv_part_coord[0];

    id = _jyp * (*nbpmax) + i;
    ettp[id] = cur_part_coord[1];
    ettpa[id] = prv_part_coord[1];

    id = _jzp * (*nbpmax) + i;
    ettp[id] = cur_part_coord[2];
    ettpa[id] = prv_part_coord[2];

    /* Velocity of the particle */

    const cs_real_t *cur_part_velocity
      = cs_lagr_particles_attr_const(cur, j, CS_LAGR_VELOCITY);
    const cs_real_t *prv_part_velocity
      = cs_lagr_particles_attr_const(prv, j, CS_LAGR_VELOCITY);

    id = _jup * (*nbpmax) + i;
    ettp[id] = cur_part_velocity[0];
    ettpa[id] = prv_part_velocity[0];

    id = _jvp * (*nbpmax) + i;
    ettp[id] = cur_part_velocity[1];
    ettpa[id] = prv_part_velocity[1];

    id = _jwp * (*nbpmax) + i;
    ettp[id] = cur_part_velocity[2];
    ettpa[id] = prv_part_velocity[2];

    /* Velocity seen by the fluid */

    const cs_real_t *cur_part_velocity_seen
      = cs_lagr_particles_attr_const(cur, j, CS_LAGR_VELOCITY_SEEN);
    const cs_real_t *prv_part_velocity_seen
      = cs_lagr_particles_attr_const(prv, j, CS_LAGR_VELOCITY_SEEN);

    id = _juf * (*nbpmax) + i;
    ettp[id] = cur_part_velocity_seen[0];
    ettpa[id] = prv_part_velocity_seen[0];

    id = _jvf * (*nbpmax) + i;
    ettp[id] = cur_part_velocity_seen[1];
    ettpa[id] = prv_part_velocity_seen[1];

    id = _jwf * (*nbpmax) + i;
    ettp[id] = cur_part_velocity_seen[2];
    ettpa[id] = prv_part_velocity_seen[2];

    id = _jtaux * (*nbpmax) + i;
    ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_TAUP_AUX);

    /* Other values based on physical model */

    if (_jreps > -1) {
      id = _jreps * (*nbpmax) + i;
      ettp[id] = cs_lagr_particles_get_real(cur, j, CS_LAGR_EMISSIVITY);
      ettpa[id] = cs_lagr_particles_get_real(prv, j, CS_LAGR_EMISSIVITY);
    }

    /* Statistical class */

    if (_jclst > -1) {
      id = _jclst * (*nbpmax) + i;
      itepa[id] = cs_lagr_particles_get_lnum(cur, j, CS_LAGR_STAT_CLASS);
    }

    /* User variables */

    if (cur->p_am->count[CS_LAGR_USER] > 0) {
      const cs_real_t
        *usr = cs_lagr_particles_attr_const(cur, j, CS_LAGR_USER);
      const cs_real_t
        *usra = cs_lagr_particles_attr_const(prv, j, CS_LAGR_USER);
      for (cs_lnum_t k = 0; k < cur->p_am->count[CS_LAGR_USER]; k++) {
        id = _jvls[k] * (*nbpmax) + i;
        ettp[id] = usr[k];
        ettpa[id] = usra[k];
      }
    }

    /* Next particle id to treat */

    j = cur->used_id[j].next_id;

  } /* End of loop on particles */

  /* New number of particles */
  *nbpart= cur->n_particles;

  /* New weight */
  *dnbpar= cur->weight;

  /* Number of exiting particles */
  *nbpout = cur->n_part_out;

  /* weight of exiting particles */
  *dnbpou = cur->weight_out;

  /* Number of depositing particles */
  *nbpdep = cur->n_part_dep;

  /* weight of depositing particles */
  *dnbdep = cur->weight_dep;

  /* Number of fouled particles (coal) */
  *npencr = cur->n_part_fou;

  /* weight of fouled particles (coal) */
  *dnpenc = cur->weight_fou;

  /* Number of failed particles */
  *nbperr = cur->n_failed_part;

  /* weight of failed particles */
  *dnbper = cur->weight_failed;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER PRTPUT\n");
  _dump_particle_set(set);
  bft_printf("\n PREV PARTICLE SET AFTER PRTPUT\n");
  _dump_particle_set(prv);
#endif
}

/*----------------------------------------------------------------------------
 * Get variables and parameters associated to each particles and keep it in
 * a new structure
 *
 * parameters:
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (getbdy, GETBDY)(const cs_int_t    *nflagm,
                          const cs_int_t    *nfrlag,
                          const cs_int_t    *injcon,
                          const cs_int_t     ilflag[],
                          const cs_int_t     iusncl[],
                          const cs_int_t     iusclb[],
                          const cs_real_t    deblag[],
                          const cs_int_t     ifrlag[])
{
  cs_lnum_t  i;
  bool  steady = false;

  assert(cs_glob_mesh != NULL);

  if (_lagr_bdy_conditions == NULL) {

    /* Define a structure with default parameters */

    _lagr_bdy_conditions = _create_bdy_cond_struct(*nflagm);

  }
  else { /* Update structure if needed */

    if (*nflagm > _lagr_bdy_conditions->n_b_max_zones) {

      _resize_bdy_cond_struct(*nflagm);

      assert(_lagr_bdy_conditions->steady_bndy_conditions == false);

    }

    if (_lagr_bdy_conditions->steady_bndy_conditions == true)
      steady = true;

  }

  if (steady == false) {

    _lagr_bdy_conditions->n_b_zones = *nfrlag;
    _lagr_bdy_conditions->continuous_injection = *injcon;

    for (i = 0; i < _lagr_bdy_conditions->n_b_zones; i++) {

      cs_lnum_t zone_id = ilflag[i] - 1;

      assert(zone_id > -1);

      _lagr_bdy_conditions->particle_flow_rate[zone_id] = deblag[zone_id];
      _lagr_bdy_conditions->b_zone_lst[zone_id] = ilflag[zone_id]; // To be deleted
      _lagr_bdy_conditions->b_zone_classes[zone_id] = iusncl[zone_id];
      _lagr_bdy_conditions->b_zone_natures[zone_id] = iusclb[zone_id];

    }

    for (i = 0; i < cs_glob_mesh->n_b_faces; i++)
      _lagr_bdy_conditions->b_face_zone_num[i] = ifrlag[i];

  } /* End if steady == false */

}

/*----------------------------------------------------------------------------
 * Displacement of particles.
 *
 * parameters:
 *   p_n_particles <-> pointer to the number of particles
 *   scheme_order  <-> current order of the scheme used for Lagrangian
 *----------------------------------------------------------------------------*/

void
CS_PROCF (dplprt, DPLPRT)(cs_lnum_t        *p_n_particles,
                          cs_lnum_t        *p_scheme_order,
                          cs_real_t         boundary_stat[],
                          const cs_lnum_t  *iensi3,
                          const cs_lnum_t  *inbr,
                          const cs_lnum_t  *inbrbd,
                          const cs_lnum_t  *iflm,
                          const cs_lnum_t  *iflmbd,
                          const cs_lnum_t  *iang,
                          const cs_lnum_t  *iangbd,
                          const cs_lnum_t  *ivit,
                          const cs_lnum_t  *ivitbd,
                          const cs_lnum_t  *iencnb,
                          const cs_lnum_t  *iencma,
                          const cs_lnum_t  *iencdi,
                          const cs_lnum_t  *iencck,
                          const cs_lnum_t  *iencnbbd,
                          const cs_lnum_t  *iencmabd,
                          const cs_lnum_t  *iencdibd,
                          const cs_lnum_t  *iencckbd,
                          const cs_lnum_t  *inclg,
                          const cs_lnum_t  *iscovc,
                          const cs_lnum_t  *nusbor,
                          cs_lnum_t         iusb[],
                          cs_real_t         visc_length[],
                          cs_real_t         dlgeo[],
                          cs_real_t         energt[],
                          const cs_real_t   tprenc[],
                          const cs_real_t   visref[],
                          const cs_real_t   enc1[],
                          const cs_real_t   enc2[],
                          const cs_real_t   *tkelvi)
{
  cs_lnum_t  i, j;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  cs_lnum_t  n_delete_particles = 0;
  cs_lnum_t  n_failed_particles = 0;

  cs_real_t  failed_particle_weight = 0.0;
  cs_real_t  r_weight = 0.0;
  cs_real_t  tot_weight = 0.0;

  cs_lnum_t  scheme_order = *p_scheme_order;
  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;

  const cs_lagr_attribute_map_t  *p_am = set->p_am;
  const cs_lagr_attribute_map_t  *prev_p_am = prev_set->p_am;

  const cs_lagr_param_t *lagr_param = &cs_glob_lagr_param;

  const cs_lnum_t  failsafe_mode = 0; /* If 1 : stop as soon as an error is
                                         detected */

  const cs_field_t *u = CS_F_(u);
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;

  cs_lnum_t n_b_faces  = mesh->n_b_faces;

  assert(set != NULL && prev_set != NULL);

  /* Main loop on  particles: global propagation */

  while (_continue_displacement()) {

    n_delete_particles = 0;
    n_failed_particles = 0;

    r_weight = 0.0;
    tot_weight = 0.0;
    failed_particle_weight = 0.0;

    assert(set->first_free_id != -1);

    /* Local propagation */

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

      unsigned char *cur_part = set->p_buffer + p_am->extents * j;
      unsigned char *prev_part = prev_set->p_buffer + prev_p_am->extents * j;

      cs_lnum_t next_id = set->used_id[j].next_id;

      /* Local copies of the current and previous particles state vectors
         to be used in case of the first pass of _local_propagation fails */

      cs_lnum_t cur_part_state
        = cs_lagr_particles_get_lnum(set, j, CS_LAGR_STATE);

      if (cur_part_state == CS_LAGR_PART_TO_SYNC) {
        cur_part_state = _local_propagation(prev_part,
                                             cur_part,
                                             prev_p_am,
                                             p_am,
                                             scheme_order,
                                             failsafe_mode,
                                             boundary_stat,
                                             &n_failed_particles,
                                             &failed_particle_weight,
                                             iensi3,
                                             inbr,
                                             inbrbd,
                                             iflm,
                                             iflmbd,
                                             iang,
                                             iangbd,
                                             ivit,
                                             ivitbd,
                                             iencnb,
                                             iencma,
                                             iencdi,
                                             iencck,
                                             iencnbbd,
                                             iencmabd,
                                             iencdibd,
                                             iencckbd,
                                             inclg,
                                             iscovc,
                                             nusbor,
                                             iusb,
                                             visc_length,
                                             dlgeo,
                                             u,
                                             energt,
                                             tprenc,
                                             visref,
                                             enc1,
                                             enc2,
                                             tkelvi,
                                             0);
        cs_lagr_particles_set_lnum(set, j, CS_LAGR_STATE, cur_part_state);
      }

      /* Treatment of particle mass flow at the boundary faces in case
        the particle rolls and changes rank */

      cs_lnum_t cur_part_depo = (lagr_param->deposition > 0) ?
        cs_lagr_particles_get_lnum(set, j, CS_LAGR_DEPOSITION_FLAG) : 0;

      if (   cur_part_state == CS_LAGR_PART_TO_SYNC
          && cur_part_depo == CS_LAGR_PART_ROLLING) {

        cs_lagr_particles_set_lnum(set, j, CS_LAGR_RANK_FLAG, 0);

        _test_wall_cell(prev_part, prev_p_am, visc_length, dlgeo);

        cs_lnum_t prev_neighbor_face_id
          = cs_lagr_particles_get_lnum(prev_set, j, CS_LAGR_NEIGHBOR_FACE_ID);
        cs_real_t face_area
          = _get_norm(b_face_normal[prev_neighbor_face_id]);

        cs_real_t prev_stat_weight
          = cs_lagr_particles_get_real(prev_set, j, CS_LAGR_STAT_WEIGHT);
        cs_real_t prev_mass
          = cs_lagr_particles_get_real(prev_set, j, CS_LAGR_MASS);

        boundary_stat[(*iflm -1) * n_b_faces + prev_neighbor_face_id]
          -=  prev_stat_weight * prev_mass / face_area;

        cs_lagr_particles_set_lnum(set, j, CS_LAGR_RANK_FLAG, 1);

      }

      if (set->p_am->count[CS_LAGR_DEPOSITION_FLAG] > 0)
        cs_lagr_particles_set_lnum
          (prev_set, j, CS_LAGR_DEPOSITION_FLAG,
           cs_lagr_particles_get_lnum(set, j, CS_LAGR_DEPOSITION_FLAG));

      cs_lagr_particles_set_lnum
        (prev_set, j, CS_LAGR_STATE,
         cs_lagr_particles_get_lnum(set, j, CS_LAGR_STATE));

      assert(next_id == set->used_id[j].next_id);

      j = next_id;

    } /* End of loop on particles */

    /* Update of the particle set structure. Delete particles. */

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

      cs_lnum_t next_id = set->used_id[j].next_id;

      // FIXME: assert(cur_part.state == prev_part.state);

      cs_lnum_t cur_part_state
        = cs_lagr_particles_get_lnum(set, j, CS_LAGR_STATE);

      cs_real_t cur_part_stat_weight
        = cs_lagr_particles_get_real(set, j, CS_LAGR_STAT_WEIGHT);

      if (   cur_part_state == CS_LAGR_PART_TO_DELETE
          || cur_part_state == CS_LAGR_PART_OUT
          || cur_part_state == CS_LAGR_PART_ERR) {

        _remove_particle(set, j);

        n_delete_particles++;

        r_weight += cur_part_stat_weight;

        cs_lnum_t cur_part_depo = (lagr_param->deposition > 0) ?
          cs_lagr_particles_get_lnum(set, j, CS_LAGR_DEPOSITION_FLAG) : 0;

        if (cur_part_depo ==  CS_LAGR_PART_ROLLING) {

          /* TODO check if we should use previous or current values */

          cs_lnum_t prev_neighbor_face_id
            = cs_lagr_particles_get_lnum(prev_set, j, CS_LAGR_NEIGHBOR_FACE_ID);
          cs_real_t face_area
            = _get_norm(b_face_normal[prev_neighbor_face_id]);

          unsigned char *prev_part = prev_set->p_buffer + prev_p_am->extents*j;

          cs_real_t prev_stat_weight
            = cs_lagr_particles_get_real(prev_set, j, CS_LAGR_STAT_WEIGHT);
          cs_real_t prev_mass
            = cs_lagr_particles_get_real(prev_set, j, CS_LAGR_MASS);

          _test_wall_cell(&prev_part, prev_p_am, visc_length, dlgeo);
          boundary_stat[(*iflm -1) * n_b_faces + prev_neighbor_face_id]
            -=  prev_stat_weight * prev_mass / face_area;
        }
      }
      else {

        tot_weight += cur_part_stat_weight;

      }

      /* next_id was saved before potentially calling _remove_particle() */

      j = next_id;

    }

    set->n_particles -= n_delete_particles;
    prev_set->n_particles -= n_delete_particles;

    set->weight = tot_weight;
    prev_set->weight = tot_weight;

    set->n_part_out += n_delete_particles;
    set->weight_out += r_weight;

    set->n_failed_part += n_failed_particles;
    set->weight_failed = failed_particle_weight;

    /*  assert(j == -1);  After a loop on particles, next_id of the last
        particle must not be defined */

    if (mesh->halo != NULL) {

      /* Synchronisation of a selection of particles for parallelism and
         periodicity. Number of particles on the local rank may change. */

      _lagr_halo_sync();
    }

  } /* End of while (global displacement) */

   /* Particle resuspension specific treatment */
   /* Update of particle mass flow at the boundary faces in
      case the particle rolls */

  if (lagr_param->resuspension > 0) {

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

      unsigned char *cur_part = set->p_buffer + p_am->extents * j;
      unsigned char *prev_part = prev_set->p_buffer + prev_p_am->extents * j;

      cs_lnum_t next_id = set->used_id[j].next_id;

      _test_wall_cell(prev_part, prev_p_am, visc_length, dlgeo);
      _test_wall_cell(cur_part, p_am, visc_length, dlgeo);

      cs_lnum_t cur_neighbor_face_id
        = cs_lagr_particles_get_lnum(set, j, CS_LAGR_NEIGHBOR_FACE_ID);
      cs_lnum_t prev_neighbor_face_id
        = cs_lagr_particles_get_lnum(prev_set, j, CS_LAGR_NEIGHBOR_FACE_ID);

      cs_lnum_t prev_part_depo
        = cs_lagr_particles_get_lnum(prev_set, j, CS_LAGR_DEPOSITION_FLAG);

      cs_real_t face_area
        = _get_norm(b_face_normal[cur_neighbor_face_id]);

      if (   cur_neighbor_face_id != prev_neighbor_face_id
          && prev_part_depo == CS_LAGR_PART_ROLLING) {

        cs_real_t cur_stat_weight
          = cs_lagr_particles_get_real(set, j, CS_LAGR_STAT_WEIGHT);
        cs_real_t cur_mass
          = cs_lagr_particles_get_real(set, j, CS_LAGR_MASS);

        boundary_stat[(*iflm -1) * n_b_faces + cur_neighbor_face_id]
          += cur_stat_weight * cur_mass / face_area;

        cs_real_t prev_stat_weight
          = cs_lagr_particles_get_real(prev_set, j, CS_LAGR_STAT_WEIGHT);
        cs_real_t prev_mass
          = cs_lagr_particles_get_real(prev_set, j, CS_LAGR_MASS);

        face_area = _get_norm(b_face_normal[prev_neighbor_face_id]);

        if (cs_lagr_particles_get_lnum(set, j, CS_LAGR_RANK_FLAG) == 0)
          boundary_stat[(*iflm -1) * n_b_faces + prev_neighbor_face_id]
            -= prev_stat_weight * prev_mass / face_area;

      }
      cs_lagr_particles_set_lnum(set, j, CS_LAGR_RANK_FLAG, 0);
      j = next_id;
    }
  }

  /* Deposition sub-model additional loop */

  if (lagr_param->deposition > 0) {

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

      unsigned char *cur_part = set->p_buffer + p_am->extents * j;

      cs_lnum_t next_id = set->used_id[j].next_id;

      _test_wall_cell(cur_part, p_am, visc_length, dlgeo);

      cs_real_t cur_part_yplus
        = cs_lagr_particles_get_real(set, j, CS_LAGR_YPLUS);

      if (cur_part_yplus < 100.e0) {

        /* TODO: specific treatment */

      }

      /* Particle resuspension specific treatment */

      if (lagr_param->resuspension > 0) {

        if (   cs_lagr_particles_get_lnum(set, j, CS_LAGR_DEPOSITION_FLAG)
            == 2) {

          const cs_real_t *cur_part_coord
            = cs_lagr_particles_attr_const(set, j, CS_LAGR_COORDS);
          const cs_real_t *prev_part_coord
            = cs_lagr_particles_attr_const(prev_set, j, CS_LAGR_COORDS);

          double traj[3] = {cur_part_coord[0] - prev_part_coord[0],
                            cur_part_coord[1] - prev_part_coord[1],
                            cur_part_coord[2] - prev_part_coord[2]};

          cs_real_t *cur_part_displacement_norm
            = cs_lagr_particles_attr(set, j, CS_LAGR_DISPLACEMENT_NORM);

          *cur_part_displacement_norm += _get_norm(traj);

        }

      }

      j = next_id;

    }
  }

  /* Returns pointers */

  *p_n_particles =  set->n_particles;
}

/*----------------------------------------------------------------------------
 * Update C structures metadata after particle computations.
 *
 * This metadata is overwritten and rebuilt at each time step, so
 * it is useful only for a possible postprocessing step.
 *
 * The matching data is copied separately, as it may not need to be
 * updated at each time step.
 *
 * parameters:
 *   nbpmax <-- n_particles max.
 *   nbpart <-- number of current particles
 *   dnbpar <-- particle total weight
 *   nbpout <-- number of outgoing particles
 *   dnbpou <-- outgoing particle total weight
 *   nbperr <-- number of failed particles
 *   dnbper <-- failed particles total weight
 *   nbpdep <-- number of depositing particles
 *   dnbdep <-- depositing particles total weight
 *   npencr <-- number of fouled particles (coal)
 *   dnpenc <-- fouled particles (coal) total weight
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (ucdprt, UCDPRT)(const cs_lnum_t   *nbpmax,
                          const cs_lnum_t   *nbpart,
                          const cs_real_t   *dnbpar,
                          const cs_int_t    *nbpout,
                          const cs_real_t   *dnbpou,
                          const cs_int_t    *nbperr,
                          const cs_real_t   *dnbper,
                          const cs_int_t    *nbpdep,
                          const cs_real_t   *dnbdep,
                          const cs_int_t    *npencr,
                          const cs_real_t   *dnpenc,
                          const cs_real_t    ettp[],
                          const cs_real_t    ettpa[],
                          const cs_int_t     itepa[],
                          const cs_real_t    tepa[],
                          const cs_int_t     ibord[],
                          const cs_lnum_t    indep[])
{
  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;

  assert(*nbpmax == set->n_particles_max); /* Up to now, we don't manage
                                              a mofification of nbpmax */

  set->n_particles = *nbpart;
  prev_set->n_particles = *nbpart;

  set->weight = *dnbpar;
  set->n_part_out = *nbpout;
  set->n_part_dep = *nbpdep;
  set->n_part_fou = *npencr;
  set->weight_out = *dnbpou;
  set->weight_dep = *dnbdep;
  set->weight_fou = *dnpenc;
  set->n_failed_part = *nbperr;
  set->weight_failed = *dnbper;

  _update_c_from_fortran(nbpmax,
                         ettp,
                         ettpa,
                         itepa,
                         tepa,
                         ibord,
                         indep);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data extents for a given particle attribute.
 *
 * For attributes not currently present, the displacement and data
 * size should be -1 and 0 respectively.
 *
 * \param[in]   particles  associated particle set
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
    *displ = particles->p_am->displ[attr];
  if (datatype)
    *datatype = particles->p_am->datatype[attr];
  if (count)
    *count = particles->p_am->count[attr];
}

/*----------------------------------------------------------------------------
 * Return pointers to the main cs_lagr_particle_set_t structures.
 *
 * parameters:
 *   current_set  --> pointer to current particle set, or NULL
 *   previous_set --> pointer to previous particle set, or NULL
 *----------------------------------------------------------------------------*/

void
cs_lagr_get_particle_sets(cs_lagr_particle_set_t  **current_set,
                          cs_lagr_particle_set_t  **previous_set)
{
  if (current_set != NULL)
    *current_set = _particle_set;
  if (previous_set != NULL)
    *previous_set = _prev_particle_set;
}

/*----------------------------------------------------------------------------
 * Delete cs_lagr_particle_set_t structure and delete other useful buffers.
 *----------------------------------------------------------------------------*/

void
cs_lagr_destroy(void)
{
  if (_particle_set == NULL)
    return;

  /* Destroy particle sets */

  _prev_particle_set = _destroy_particle_set(_prev_particle_set);
  _particle_set = _destroy_particle_set(_particle_set);

  _destroy_attr_map(&_p_attr_map);
  _destroy_attr_map(&_p_attr_map_prev);

  /* Destroy builder */
  _particle_track_builder = _destroy_track_builder(_particle_track_builder);

  /* Destroy boundary condition structure */

  _lagr_bdy_conditions = _destroy_bdy_cond_struct(_lagr_bdy_conditions);

  /* Destroy the structure dedicated to clogging modeling */

  if (cs_glob_lagr_param.clogging)
    cs_lagr_clogging_finalize();

  /* Destroy the structure dedicated to roughness surface modeling */

  if (cs_glob_lagr_param.rough)
    cs_lagr_roughness_finalize();

  /* Delete MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)  _delete_particle_datatypes();
#endif

  /* Free local indexes */

  BFT_FREE(_jthp);
  BFT_FREE(_jmch);
  BFT_FREE(_jmck);
  BFT_FREE(_jrhock);
  BFT_FREE(_jvls);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_lagr_particle_t structure
 *
 * \param[in]  particles  cs_lagr_particle_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_dump(const cs_lagr_particle_set_t  *particles)
{
  if (particles != NULL) {

    cs_lnum_t  i, j;

    bft_printf("Particle set\n");
    bft_printf("------------\n");
    bft_printf("  n_particles:      %10d\n", particles->n_particles);
    bft_printf("  n_particles_max:  %10d\n", particles->n_particles_max);

    if (particles->used_id != NULL) {

      bft_printf("  first_used_id:    %10d\n", particles->first_used_id);
      bft_printf("  first_free_id:    %10d\n", particles->first_free_id);
      bft_printf_flush();

      if (particles->n_particles > 0) {
        for (i = 0, j = particles->first_used_id;
             i < particles->n_particles && j > -1;
             i++) {
          bft_printf("  dump_particle_set i j = %d %d \n", i, j);
          _dump_particle(particles, j);
          j = particles->used_id[j].next_id;
        }
        assert(j == -1); /* The next_id is not defined for the last particle
                            of the particle set */
      }

    }
    else {

      bft_printf_flush();

      if (particles->n_particles > 0) {
        for (i = 0; i < particles->n_particles; i++) {
          bft_printf("  dump_particle_set i = %d \n", i);
          _dump_particle(particles, i);
          j = particles->used_id[j].next_id;
        }
        assert(j == -1); /* The next_id is not defined for the last particle
                            of the particle set */
      }

    }
  }
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

/* Delete local macro definitions */

END_C_DECLS
