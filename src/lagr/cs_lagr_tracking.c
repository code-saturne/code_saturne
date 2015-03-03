/*============================================================================
 * Methods for particle localization
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
 * Functions dealing with particle tracking
 *============================================================================*/

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
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
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

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define  N_GEOL 13
#define  CS_LAGR_MIN_COMM_BUF_SIZE  8

/*=============================================================================
 * Local Enumeration definitions
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
  CS_LAGR_PART_STUCK     = 3,
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

/* Tracking error types */

typedef enum {

  CS_LAGR_TRACKING_OK,
  CS_LAGR_TRACKING_ERR_MAX_LOOPS

} cs_lagr_tracking_error_t;

/* keys to sort attributes by Fortran array and index (for mapping)
   eptp/eptpa real values at current and previous time steps

   ieptp/ieptpa integer values at current and previous time steps
   pepa real values at current time step
   _int_loc local integer values at current time step
   ipepa integer values at current time step
   iprkid values are for rank ids, useful and valid only for previous
   time steps */

typedef enum {
  EPTP_TS = 1, /* EPTP with possible source terms */
  EPTP,
  IEPTP,
  PEPA,
  IPEPA,
  IPRKID,
} _array_map_id_t;

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/* Private tracking data associated to each particle */
/* --------------------------------------------------*/

/* This structure is currently mapped to the beginning of each
 * particle's data, and contains values which are used during the
 * tracking algorithm only.
 * It could be separated in the future, but this would require
 * keeping track of a particle's local id in most functions
 * of this file. */

typedef struct {

  cs_real_t  start_coords[3];  /* starting coordinates for next displacement */

  cs_lnum_t  last_face_num;    /* last face number encountered */

  int        state;            /* CS_LAGR_PART_TO_DELETE, CS_LAGR_PART_TO_SYNC,
                                  CS_LAGR_PART_TREATED, CS_LAGR_PART_STUCK,
                                  CS_LAGR_PART_OUT, or CS_LAGR_PART_ERR */

} cs_lagr_tracking_info_t;

/* Linked list */
/* ----------- */

struct _cs_lagr_tracking_list_t {

  cs_lnum_t   prev_id;  /* id in particle set of the previous particle */
  cs_lnum_t   next_id;  /* id in particle set of the next particle */

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

/* Manage the exchange of particles between communicating ranks */
/* -------------------------------------------------------------*/

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
  size_t      extents;        /* Extents for particle set */

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

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Enumerator names */

const char *cs_lagr_attribute_name[] = {
  "CS_LAGR_CELL_NUM",
  "CS_LAGR_RANK_ID",
  "CS_LAGR_SWITCH_ORDER_1",
  "CS_LAGR_RANDOM_VALUE",
  "CS_LAGR_STAT_WEIGHT",
  "CS_LAGR_RESIDENCE_TIME",
  "CS_LAGR_MASS",
  "CS_LAGR_DIAMETER",
  "CS_LAGR_TAUP_AUX",
  "CS_LAGR_COORDS",
  "CS_LAGR_VELOCITY",
  "CS_LAGR_VELOCITY_SEEN",
  "CS_LAGR_TURB_STATE_1",
  "CS_LAGR_PRED_VELOCITY",
  "CS_LAGR_PRED_VELOCITY_SEEN",
  "CS_LAGR_YPLUS",
  "CS_LAGR_INTERF",
  "CS_LAGR_NEIGHBOR_FACE_ID",
  "CS_LAGR_MARKO_VALUE",
  "CS_LAGR_DEPOSITION_FLAG",
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
  "CS_LAGR_STAT_CLASS",
  "CS_LAGR_USER",
  "CS_LAGR_N_ATTRIBUTES"};

/* Global variable for the current subroutines */

static  cs_lagr_particle_set_t  *_particle_set = NULL;
static  cs_lagr_track_builder_t  *_particle_track_builder = NULL;
static  cs_lagr_bdy_condition_t  *_lagr_bdy_conditions = NULL;

static cs_lagr_attribute_map_t  *_p_attr_map = NULL;

enum {X, Y, Z};  /* Used for _get_norm() and _get_dot_prod() */

static  double              _reallocation_factor = 2.0;
static  unsigned long long  _n_g_max_particles = ULLONG_MAX;

static  int                 _max_propagation_loops = 100;

/* MPI datatype associated to each particle "structure" */

#if defined(HAVE_MPI)
static  MPI_Datatype  _cs_mpi_particle_type;
#endif

/* Global Lagragian module parameters and associated pointer
   Should move to cs_lagr.c */

static cs_lagr_param_t  _lagr_param = {0, 1, 0, 0, 0, 0, 0, 0, 1};
const  cs_lagr_param_t  *cs_glob_lagr_params = &_lagr_param;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_lagr_pointers(int          dim_ipepa[2],
                   int          dim_pepa[2],
                   int          dim_eptp[2],
                   int          dim_eptpa[2],
                   int          dim_ptsvar[2],
                   cs_int_t   **p_ipeta,
                   cs_real_t  **p_peta,
                   cs_real_t  **p_eptp,
                   cs_real_t  **p_eptpa,
                   cs_real_t  **p_ptsvar,
                   cs_int_t   **p_nbpart,
                   cs_real_t  **p_dnbpar,
                   cs_int_t   **p_nbpout,
                   cs_real_t  **p_dnbpou,
                   cs_int_t   **p_nbperr,
                   cs_real_t  **p_dnbper,
                   cs_int_t   **p_nbpdep,
                   cs_real_t  **p_dnbdep,
                   cs_int_t   **p_npencr,
                   cs_real_t  **p_dnpenc);

void
cs_f_lagr_part_copy(cs_lnum_t  dest,
                    cs_lnum_t  src);

void
cs_f_lagr_current_to_previous(cs_lnum_t  pn);

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

  p_am->extents = _align_extents(sizeof(cs_lagr_tracking_info_t));

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
      case EPTP_TS:
      case EPTP:
        max_time_id = 1;
        break;
      case IEPTP:
        datatype = CS_LNUM_TYPE;
        max_time_id = 1;
        break;
      case PEPA:
        break;
      case IPEPA:
        datatype = CS_LNUM_TYPE;
        break;
      case IPRKID:
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

  if (cs_glob_lagr_params->t_order > 1) {

    BFT_MALLOC(p_am->source_term_displ, CS_LAGR_N_ATTRIBUTES, ptrdiff_t);

    /* loop again on ordered attributes */

    for (int i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {

      attr = order[i];

      /* Add attribute to map if in EPTPS array */

      if (   attr_keys[attr][0] == EPTP_TS
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

    BFT_FREE(_p_am->displ);
    BFT_FREE(_p_am->count);

    BFT_FREE(*p_am);
  }
}

/*----------------------------------------------------------------------------
 * Get pointer to a particle's tracking information.
 *
 * parameters:
 *   particle_set <-> pointer to particle set
 *   particle_id  <-- particle id
 *
 * returns:
 *   pointer to particle state structure
 *----------------------------------------------------------------------------*/

inline static cs_lagr_tracking_info_t *
_tracking_info(cs_lagr_particle_set_t  *particle_set,
               cs_lnum_t                particle_id)
{
  return (cs_lagr_tracking_info_t *)(  particle_set->p_buffer
                                     + particle_set->p_am->extents*particle_id);
}

/*----------------------------------------------------------------------------
 * Get const pointer to a particle's tracking information.
 *
 * parameters:
 *   particle_set <-> pointer to particle set
 *   particle_id  <-- particle id
 *
 * returns:
 *   pointer to particle state structure
 *----------------------------------------------------------------------------*/

inline static const cs_lagr_tracking_info_t *
_get_tracking_info(const cs_lagr_particle_set_t  *particle_set,
                   cs_lnum_t                      particle_id)
{
  return (const cs_lagr_tracking_info_t *)
    (particle_set->p_buffer + particle_set->p_am->extents*particle_id);
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
 *   set <->  particles set
 *   id  <-- id of particle to remove
 *----------------------------------------------------------------------------*/

static void
_remove_particle(cs_lagr_particle_set_t   *set,
                 cs_lnum_t                 id)
{
  /* Remove particle from "used particle list" */

  cs_lnum_t prev_id = set->used_id[id].prev_id;
  cs_lnum_t next_id = set->used_id[id].next_id;

  if (prev_id != -1)
    set->used_id[prev_id].next_id = next_id;
  else
    set->first_used_id = next_id;

  if (next_id != set->n_particles_max && next_id != -1)
    set->used_id[next_id].prev_id = prev_id;

  /* Add particle to "free particle list" */

  if (id < set->first_free_id) {

    cs_lnum_t old_first_free = set->first_free_id;
    cs_lnum_t old_first_free_prev = set->used_id[old_first_free].prev_id;

    set->first_free_id = id;

    set->used_id[set->first_free_id].next_id = old_first_free;
    set->used_id[set->first_free_id].prev_id = old_first_free_prev;

    set->used_id[old_first_free].prev_id = id;

  }
  else { /* We place the particle just behind the first free particle. */

    cs_lnum_t first_free = set->first_free_id;
    cs_lnum_t old_next = set->used_id[first_free].next_id;

    set->used_id[first_free].next_id = id;

    set->used_id[id].next_id = old_next;
    set->used_id[id].prev_id = first_free;

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype which maps main particle characteristics.
 *
 * parameters:
 *   am  <-- attributes map
 *
 * returns:
 *   MPI_Datatype matching given attributes map
 *----------------------------------------------------------------------------*/

static MPI_Datatype
_define_particle_datatype(const cs_lagr_attribute_map_t  *p_am)
{
  size_t i;
  MPI_Datatype  new_type;
  int           count;
  cs_datatype_t *cs_type;
  int           *blocklengths;
  MPI_Datatype  *types;
  MPI_Aint      *displacements;

  size_t tot_extents = p_am->extents;

  /* Mark bytes with associated type */

  BFT_MALLOC(cs_type, tot_extents, cs_datatype_t);

  for (i = 0; i < tot_extents; i++)
    cs_type[i] = CS_CHAR;

  /* Map tracking info */

  size_t attr_start, attr_end;

  attr_start = offsetof(cs_lagr_tracking_info_t, start_coords);
  attr_end = attr_start + 3*sizeof(cs_real_t);
  for (i = attr_start; i < attr_end; i++)
    cs_type[i] = CS_REAL_TYPE;

  attr_start = offsetof(cs_lagr_tracking_info_t, last_face_num);
  attr_end = attr_start + sizeof(cs_lnum_t);
  for (i = attr_start; i < attr_end; i++)
    cs_type[i] = CS_LNUM_TYPE;

  attr_start = offsetof(cs_lagr_tracking_info_t, last_face_num);
  attr_end = attr_start + sizeof(int);
  for (i = attr_start; i < attr_end; i++)
    cs_type[i] = CS_INT_TYPE;

  /* Map attributes */

  for (int j = 0; j < p_am->n_time_vals; j++) {

    for (size_t attr = 0; attr < CS_LAGR_N_ATTRIBUTES; attr++) {
      if (p_am->count[j][attr] > 0) {
        assert(p_am->displ[j][attr] > -1);
        size_t b_size
          = p_am->count[j][attr] * cs_datatype_size[p_am->datatype[attr]];
        for (i = 0; i < b_size; i++)
          cs_type[p_am->displ[j][attr] + i] = p_am->datatype[attr];
      }
    }

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

  MPI_Type_commit(&new_type);

  return new_type;
}

/*----------------------------------------------------------------------------
 * Delete all the MPI_Datatypes related to particles.
 *----------------------------------------------------------------------------*/

static void
_delete_particle_datatypes(void)
{
  MPI_Type_free(&_cs_mpi_particle_type);
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

  for (int time_id = 0; time_id < particles->p_am->n_time_vals; time_id++) {

    if (time_id == 0)
      bft_printf("    values at time n:\n");
    else
      bft_printf("    values at time: n-%d\n", time_id);

    for (cs_lagr_attribute_t attr = 0;
         attr < CS_LAGR_N_ATTRIBUTES;
         attr++) {
      if (am->count[time_id][attr] > 0) {
        char attr_name[64];
        strncpy(attr_name, cs_lagr_attribute_name[attr] + 8, 63);
        attr_name[63] = '\0';
        for (int i = 0; attr_name[i] != '\0'; i++)
          attr_name[i] = tolower(attr_name[i]);
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
_resize_particle_set(cs_lagr_particle_set_t   *particle_set,
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

    if (particle_set->time_id == 0)
      BFT_REALLOC(particle_set->used_id,
                  particle_set->n_particles_max,
                  cs_lagr_tracking_list_t);

    retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a cs_lagr_halo_t structure to deal with parallelism and
 * periodicity
 *
 * parameters:
 *   extents  <-- extents for particles of set
 *
 * returns:
 *   a new allocated cs_lagr_halo_t structure.
 *----------------------------------------------------------------------------*/

static cs_lagr_halo_t *
_create_lagr_halo(size_t  extents)
{
  cs_lnum_t  i, rank, tr_id, shift, start, end, n;

  cs_lnum_t  halo_cell_id = 0;
  cs_lnum_t  *cell_num = NULL;
  cs_lagr_halo_t  *lagr_halo = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const cs_lnum_t  n_halo_cells = halo->n_elts[CS_HALO_EXTENDED];

  BFT_MALLOC(lagr_halo, 1, cs_lagr_halo_t);

  assert(n_halo_cells == halo->index[2*halo->n_c_domains]);
  assert(n_halo_cells == mesh->n_ghost_cells);

  lagr_halo->extents = extents;
  lagr_halo->n_cells = n_halo_cells;

  /* Allocate buffers to enable the exchange between communicating ranks */

  BFT_MALLOC(lagr_halo->send_shift, halo->n_c_domains, cs_lnum_t);
  BFT_MALLOC(lagr_halo->send_count, halo->n_c_domains, cs_lnum_t);
  BFT_MALLOC(lagr_halo->recv_shift, halo->n_c_domains, cs_lnum_t);
  BFT_MALLOC(lagr_halo->recv_count, halo->n_c_domains, cs_lnum_t);

  lagr_halo->send_buf_size = CS_LAGR_MIN_COMM_BUF_SIZE;
  lagr_halo->recv_buf_size = CS_LAGR_MIN_COMM_BUF_SIZE;

  BFT_MALLOC(lagr_halo->send_buf,
             lagr_halo->send_buf_size * extents,
             unsigned char);
  BFT_MALLOC(lagr_halo->recv_buf,
             lagr_halo->recv_buf_size * extents,
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

  size_t tot_extents = lag_halo->extents;

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
      cs_lnum_t iel = mesh->i_face_cells[i][j] + 1;
      if (iel <= mesh->n_cells)
        builder->cell_face_idx[iel] += 1;
    }

  /* Count of the number of faces per cell: loop on border faces */

  for (i = 0; i < mesh->n_b_faces; i++)
    builder->cell_face_idx[mesh->b_face_cells[i] + 1] += 1;

  /* Build index */

  for (i = 0; i < mesh->n_cells; i++)
    builder->cell_face_idx[i+1] += builder->cell_face_idx[i];

  BFT_MALLOC(builder->cell_face_lst,
             builder->cell_face_idx[mesh->n_cells], cs_lnum_t);

  /* Build list: border faces are < 0 and interior faces > 0 */

  for (i = 0; i < mesh->n_i_faces; i++) {
    for (j = 0; j < 2; j++) {

      cs_lnum_t iel = mesh->i_face_cells[i][j] + 1;

      if (iel <= mesh->n_cells) {

        cs_lnum_t  cell_id = iel - 1;
        cs_lnum_t  shift = builder->cell_face_idx[cell_id] + counter[cell_id];

        builder->cell_face_lst[shift] = i+1;
        counter[cell_id] += 1;
      }
    }
  }

  for (i = 0; i < mesh->n_b_faces; i++) {

    cs_lnum_t  cell_id = mesh->b_face_cells[i];
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
 *   extents         <-- extents for particles of set
 *
 * returns:
 *   a new defined cs_lagr_track_builder_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_track_builder_t *
_init_track_builder(cs_lnum_t  n_particles_max,
                    size_t     extents)
{
  cs_mesh_t  *mesh = cs_glob_mesh;

  cs_lagr_track_builder_t  *builder = NULL;

  if (n_particles_max == 0)
    return NULL;

  BFT_MALLOC(builder, 1, cs_lagr_track_builder_t);

  /* Define a cell->face connectivity */

  _define_cell_face_connect(builder);

  /* Define a cs_lagr_halo_t structure to deal with parallelism and
     periodicity */

  if (cs_glob_mesh->n_init_perio > 0 || cs_glob_n_ranks > 1)
    builder->halo = _create_lagr_halo(extents);
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

  cs_lagr_particle_set_lnum(particle, attr_map, CS_LAGR_CELL_NUM, 0);

  n_failed_particles++;
  failed_particle_weight
    += cs_lagr_particle_get_real(particle, attr_map, CS_LAGR_STAT_WEIGHT);

  if (failsafe_mode == 1) {
    switch (error_type) {
    case CS_LAGR_TRACKING_ERR_MAX_LOOPS:
      bft_error(__FILE__, __LINE__, 0,
                _("Max number of loops reached in particle displacement."));
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
    if (_get_tracking_info(set, j)->state == CS_LAGR_PART_TO_SYNC) {
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
 *   particle      <-- particle attributes
 *   p_am          <-- pointer to attributes map
 *
 * returns:
 *   -1 if the cell-center -> particle segment does not go through the face's
 *   plane, minimum relative distance (in terms of barycentric coordinates)
 *   of intersection point to face.
 *----------------------------------------------------------------------------*/

static double
_intersect_face(cs_lnum_t                       face_num,
                cs_lnum_t                       n_vertices,
                const cs_lnum_t                 face_connect[],
                const void                     *particle,
                const cs_lagr_attribute_map_t  *p_am)
{
  cs_lnum_t  i, j;
  cs_real_t  face_surf;
  const cs_real_t  *face_norm, *face_cog, *cell_cen;

  cs_lnum_t  cur_cell_id
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM) - 1;
  const cs_real_t  *next_location
    = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_COORDS);
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const double epsilon = 1.e-15;
  const double bc_epsilon = 1.e-2;

  double retval = -1;

  cs_mesh_t  *mesh = cs_glob_mesh;

  assert(sizeof(cs_real_t) == 8);

  /* Initialization */

  cell_cen = fvq->cell_cen + (3*cur_cell_id);

  if (face_num > 0) { /* Interior  face */

    cs_lnum_t  face_id = face_num - 1;
    face_cog = fvq->i_face_cog + (3*face_id);
    face_norm = fvq->i_face_normal + (3*face_id);
    face_surf = fvq->i_face_surf[face_id];

  }
  else { /* Boundary face */

    cs_lnum_t  face_id = -face_num - 1;
    face_cog = fvq->b_face_cog + (3*face_id);
    face_norm = fvq->b_face_normal + (3*face_id);
    face_surf = fvq->b_face_surf[face_id];

    /* Avoid tangential trajectory in _boundary_treatment */
    const cs_real_t  *prev_location
      = ((const cs_lagr_tracking_info_t *)particle)->start_coords;
    cs_real_3_t depl = {next_location[0] - prev_location[0],
                        next_location[1] - prev_location[1],
                        next_location[2] - prev_location[2]};
    if ((cs_math_3_dot_product(depl, face_norm) / face_surf) < 1.1e-15)
      return retval;

  }

  cs_real_3_t s01 = {next_location[0] - cell_cen[0],
                     next_location[1] - cell_cen[1],
                     next_location[2] - cell_cen[2]};
  cs_real_3_t scf = {face_cog[0] - cell_cen[0],
                     face_cog[1] - cell_cen[1],
                     face_cog[2] - cell_cen[2]};

  /* Check basic orientation relative to plane */

  double rin = (  face_norm[0]*scf[0]
                + face_norm[1]*scf[1]
                + face_norm[2]*scf[2]) / face_surf;
  double rid = (  face_norm[0]*s01[0]
                + face_norm[1]*s01[1]
                + face_norm[2]*s01[2]) / face_surf;

  /* If direction is parallel or nearly parallel, we consider that
     exit is not done through this face */
  if (rid < epsilon && rid > -epsilon)
    return retval;

  double ri = rin/rid;

  if (ri < 0 || ri >= 1)
    return retval;

  int n_intersects = 0;

  retval = HUGE_VAL;

  /* If segment intersects plane, compute intersection barycentric
     coordinates of successive triangles:
     - loop on sub-triangles of the face and test for each, recomputing
       intersection on each triangle's plane, in case the face is not plane
     - use Moller-Trumbore algorithm (without backface culling)
       for each triangle */

  for (i = 0; i < n_vertices; i++) {

    cs_real_3_t  vtx_0, vtx_1, e0, e1;

    const cs_lnum_t vtx_id_0 = face_connect[i];
    const cs_lnum_t vtx_id_1 = face_connect[(i+1)%n_vertices];

    for (j = 0; j < 3; j++) {
      vtx_0[j] = mesh->vtx_coord[3*vtx_id_0+j];
      vtx_1[j] = mesh->vtx_coord[3*vtx_id_1+j];
      e0[j] = vtx_0[j] - face_cog[j];
      e1[j] = vtx_1[j] - face_cog[j];
    }

    /* determinant */

    const cs_real_3_t pvec = {s01[1]*e1[2] - s01[2]*e1[1],
                              s01[2]*e1[0] - s01[0]*e1[2],
                              s01[0]*e1[1] - s01[1]*e1[0]};

    const double det = cs_math_3_dot_product(e0, pvec);

    if (det/face_surf < epsilon && det/face_surf > -epsilon)
      continue;

    const double inv_det = 1./det;

    const double u = cs_math_3_dot_product(s01, pvec) * inv_det;

    const cs_real_3_t qvec = {scf[1]*e0[2] - scf[2]*e0[1],
                              scf[2]*e0[0] - scf[0]*e0[2],
                              scf[0]*e0[1] - scf[1]*e0[0]};

    const double v = cs_math_3_dot_product(s01, qvec) * inv_det;

    /* We have an intersection if u in [0,1] and u+v in [0,1] */

    const double bc1 = fabs(u - 0.5);
    const double bc2 = fabs(u + v - 0.5);
    const double max_bc = 2.0 * CS_MAX(bc1, bc2);

    /*
      intersection point (not currently used):
      const double t = cs_math_3_dot_product(e1, qvec) * inv_det;
    */

    /* Update return code */

    if (max_bc < retval)
      retval = max_bc;

    /* In case of "true" intersect */

    if (max_bc < 1.0 + bc_epsilon) {

      /* Compare orientation of triangle normal with face normal
         so as to handle non-convex cases; */

      const cs_real_3_t tn = {e0[1]*e1[2] - e0[2]*e1[1],
                              e0[2]*e1[0] - e0[0]*e1[2],
                              e0[0]*e1[1] - e0[1]*e1[0]};

      if (cs_math_3_dot_product(tn, face_norm) < 0)
        n_intersects -= 1;
      else if (max_bc < 1.0 - bc_epsilon)
        n_intersects += 1;

    }

  }

  /* In case intersections were removed due to non-convex cases,
     n_intersects < 0, but we may not be "far" from cutting this
     face. To trace this with a single "distance", we set the
     distance to 1 + nc_epsilon; this way, a better fitting face
     will intersect with a better distance, but if none is found,
     this face will be chosen. */

  if (n_intersects < 1 && retval < 1 + bc_epsilon)
    retval = 1 + bc_epsilon;

  return retval;
}

/*----------------------------------------------------------------------------
 * Determine the number of the closest wall face from the particle
 * as well as the corresponding wall normal distance (y_p^+)
 *
 * Used for the deposition model.
 *
 * parameters:
 *   particle      <-- particle attributes for current time step
 *   p_am          <-- pointer to attributes map for current time step
 *   visc_length   <--
 *   dlgeo         <-- array with various geometry values for particles
 *   yplus         --> associated yplus value
 *   face_id       --> associated neighbor wll face, or -1
 *----------------------------------------------------------------------------*/

static void
_test_wall_cell(const void                     *particle,
                const cs_lagr_attribute_map_t  *p_am,
                const cs_real_t                 visc_length[],
                const cs_real_t                 dlgeo[],
                cs_real_t                      *yplus,
                cs_lnum_t                      *face_id)
{
  cs_lnum_t cell_num
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM);

  if (cell_num < 0) return;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;
  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;
  cs_lnum_t cell_id = cell_num - 1;

  *yplus = 10000;
  *face_id = -1;

  cs_lnum_t  start = cell_face_idx[cell_id];
  cs_lnum_t  end =  cell_face_idx[cell_id + 1];

  for (cs_lnum_t i = start; i < end; i++) {
    cs_lnum_t  face_num = cell_face_lst[i];

    if (face_num < 0) {
      cs_lnum_t f_id = CS_ABS(face_num) - 1;
      cs_lnum_t b_zone_id = bdy_conditions->b_face_zone_num[f_id]-1;

      if (   (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IDEPO1)
          || (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IDEPO2)
          || (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IREBOL)
          || (bdy_conditions->b_zone_natures[b_zone_id] == CS_LAGR_IDEPFA)) {

        cs_real_t x_face = dlgeo[f_id];
        cs_real_t y_face = dlgeo[f_id + (mesh->n_b_faces)];
        cs_real_t z_face = dlgeo[f_id + (mesh->n_b_faces) * 2];
        cs_real_t offset_face = dlgeo[f_id + (mesh->n_b_faces) * 3];
        const cs_real_t  *particle_coord
          = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_COORDS);

        cs_real_t dist_norm =   CS_ABS(  particle_coord[0] * x_face
                                       + particle_coord[1] * y_face
                                       + particle_coord[2] * z_face
                                       + offset_face) / visc_length[f_id];
        if (dist_norm  < *yplus) {
          *yplus = dist_norm;
          *face_id = f_id;
        }
      }
    }

  }

}

/*----------------------------------------------------------------------------
 * Handle particles moving to boundary
 *
 * parameters:
 *   particles  <-- pointer to particle set
 *   particle   <-> particle data for current particle
 *   ...        <-> pointer to an error indicator
 *
 * returns:
 *   particle state
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_boundary_treatment(cs_lagr_particle_set_t    *particles,
                    void                      *particle,
                    cs_lnum_t                  face_num,
                    cs_real_t                 *boundary_stat,
                    cs_lnum_t                  boundary_zone,
                    cs_lnum_t                 *p_move_particle,
                    cs_lnum_t                 *p_n_failed_particles,
                    cs_real_t                 *p_failed_particle_weight,
                    const cs_lnum_t           *iensi3,
                    const cs_lnum_t           *inbr,
                    const cs_lnum_t           *inbrbd,
                    const cs_lnum_t           *iang,
                    const cs_lnum_t           *iangbd,
                    const cs_lnum_t           *ivit,
                    const cs_lnum_t           *ivitbd,
                    const cs_lnum_t           *iencnb,
                    const cs_lnum_t           *iencma,
                    const cs_lnum_t           *iencdi,
                    const cs_lnum_t           *iencck,
                    const cs_lnum_t           *iencnbbd,
                    const cs_lnum_t           *iencmabd,
                    const cs_lnum_t           *iencdibd,
                    const cs_lnum_t           *iencckbd,
                    const cs_lnum_t           *inclg,
                    const cs_lnum_t           *iscovc,
                    const cs_lnum_t           *nusbor,
                    cs_lnum_t                  iusb[],
                    cs_real_t                 *part_b_mass_flux,
                    cs_real_t                  energt[],
                    const cs_real_t            tprenc[],
                    const cs_real_t            visref[],
                    const cs_real_t            enc1[],
                    const cs_real_t            enc2[],
                    cs_real_t                  tkelvi)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const double pi = 4 * atan(1);

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_lnum_t  k;
  cs_real_t  tmp;
  cs_real_t  depl[3], face_normal[3], face_cog[3], intersect_pt[3];

  cs_real_t  compo_vel[3] = {0.0, 0.0, 0.0};
  cs_real_t  norm_vel = 0.0;

  cs_real_t  s = 0.0;
  cs_lnum_t  move_particle = *p_move_particle;
  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight = *p_failed_particle_weight;

  cs_lnum_t  face_id = face_num - 1;
  cs_lnum_t  particle_state = -999;

  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;

  cs_lnum_t  contact_number = 0;
  cs_real_t  *surface_coverage = NULL;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

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

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  assert(bdy_conditions != NULL);

  for (k = 0; k < 3; k++)
    depl[k] = particle_coord[k] - p_info->start_coords[k];

  assert(! (   fabs(depl[0]) < 1e-15
            && fabs(depl[1]) < 1e-15 && fabs(depl[2]) < 1e-15));

  for (k = 0; k < 3; k++) {
    face_normal[k] = fvq->b_face_normal[3*face_id+k];
    face_cog[k] = fvq->b_face_cog[3*face_id+k];
  }

  cs_real_t face_area  = fvq->b_face_surf[face_id];

  cs_real_t  face_norm[3] = {face_normal[0]/face_area,
                             face_normal[1]/face_area,
                             face_normal[2]/face_area};

  /* Save particle impacting velocity */
  if (*iangbd > 0 || *ivitbd > 0) {
    norm_vel = _get_norm(particle_velocity);
    for (k = 0; k < 3; k++)
      compo_vel[k] = particle_velocity[k];
  }

  tmp = cs_math_3_dot_product(depl, face_normal);

  /* 3D geometry refresher:
     ~~~~~~~~~~~~~~~~~~~~~~
     1)  equation of a plane with normal (a,b,c):
     a*x + b*y + c*z + d = 0
     2)  equation of a line passing through P and Q :
     x = XP + (XQ-XP) * s
     y = YP + (YQ-YP) * s
     z = ZP + (ZQ-ZP) * s
     where s is a real parameter

     If trajectory is in (or nearly in) the face's plane,
     the "arrival" position should be inside the face (otherwise,
     the cell center -> arrival position segment should not have
     intersected the face), so s = 1 in this case.
  */

  if (fabs(tmp)/face_area < 1e-15)
    s = 1.;

  else {
    for (k = 0; k < 3; k++)
      s += face_normal[k] * (face_cog[k] - p_info->start_coords[k]);
    s /= tmp;
  }

  for (k = 0; k < 3; k++)
    intersect_pt[k] = depl[k]*s + p_info->start_coords[k];

  if (   bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_ISORTL
      || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IENTRL
      || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_OUT;

    if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) {
      particles->n_part_dep += 1;
      particles->weight_dep += particle_stat_weight;
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

    particles->n_part_dep += 1;
    particles->weight_dep += particle_stat_weight;

    /* Specific treatment in case of particle resuspension modeling */

    cs_lnum_t *cell_num = cs_lagr_particle_attr(particle,
                                                    p_am,
                                                    CS_LAGR_CELL_NUM);

    if (cs_glob_lagr_params->resuspension == 0) {

      *cell_num = - *cell_num;

      for (k = 0; k < 3; k++)
        particle_velocity_seen[k] = 0.0;

      particle_state = CS_LAGR_PART_STUCK;

    } else {

      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_DEPOSITED);
      *cell_num = cs_glob_mesh->b_face_cells[face_id] + 1;

      particle_state = CS_LAGR_PART_TREATED;

    }

  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPFA) {

    const cs_real_t particle_diameter
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id] + 1);

    cs_real_t uxn = particle_velocity[0] * face_norm[0];
    cs_real_t vyn = particle_velocity[1] * face_norm[1];
    cs_real_t wzn = particle_velocity[2] * face_norm[2];

    cs_real_t energ = 0.5 * particle_mass * (uxn+vyn+wzn) * (uxn+vyn+wzn);
    cs_real_t min_porosity ;
    cs_real_t limit;

    if (cs_glob_lagr_params->clogging) {

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

    if (cs_glob_lagr_params->roughness > 0)
      cs_lagr_roughness_barrier(particle,
                                p_am,
                                face_id,
                                &energt[face_id]);

    if (energ > energt[face_id] * 0.5 * particle_diameter) {

      /* The particle deposits*/
      if (!cs_glob_lagr_params->clogging && !cs_glob_lagr_params->resuspension) {
        move_particle = CS_LAGR_PART_MOVE_OFF;

        /* Set negative value for current cell number */
        cs_lagr_particle_set_lnum
          (particle, p_am, CS_LAGR_CELL_NUM,
           - cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM));

        particles->n_part_dep += 1;
        particles->weight_dep += particle_stat_weight;

        particle_state = CS_LAGR_PART_STUCK;
      }

      if (cs_glob_lagr_params->resuspension > 0) {

        move_particle = CS_LAGR_PART_MOVE_OFF;
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                  CS_LAGR_PART_DEPOSITED);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                  cs_glob_mesh->b_face_cells[face_id] + 1);

        for (k = 0; k < 3; k++) {
          particle_velocity[k] = 0.0;
          particle_coord[k]
            = intersect_pt[k] - 0.5 * particle_diameter * face_norm[k];
        }
        particles->n_part_dep += 1;
        particles->weight_dep += particle_stat_weight;
        particle_state = CS_LAGR_PART_TREATED;

      }

      if (cs_glob_lagr_params->clogging) {

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
            for (i = 0, j = particles->first_used_id;
                 i < particles->n_particles;
                 i++) {

              cur_part = (const void *)(particles->p_buffer + p_am->extents * j);

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
              j = particles->used_id[j].next_id;
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
              for (i = 0, j = particles->first_used_id;
                   i < particles->n_particles;
                   i++) {

                const void * cur_part2 = particles->p_buffer + p_am->extents * j;

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
                j = particles->used_id[j].next_id;
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
        const cs_real_t *xyzcen = fvq->cell_cen;
        const cs_real_t *volume  = fvq->cell_vol;

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

        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
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
            (particle, p_am, CS_LAGR_CELL_NUM,
             - cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM));

          particles->n_part_dep += 1;
          particles->weight_dep += particle_stat_weight;

          particle_state = CS_LAGR_PART_STUCK;
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                    CS_LAGR_PART_DEPOSITED);

          /* Update of the number of stat. weight of deposited particles */
          /* in case of clogging modeling                                */

          boundary_stat[(*inclg -1) * n_b_faces + face_id]
            += particle_stat_weight;

        }
        else {

          /*The particle does not deposit: It 'rebounds' */
          /* because the porosity threshold is reached   */

          move_particle = CS_LAGR_PART_MOVE_ON;
          particle_state = CS_LAGR_PART_TO_SYNC;
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                    cs_glob_mesh->b_face_cells[face_id] + 1);
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                    CS_LAGR_PART_IN_FLOW);

          for (k = 0; k < 3; k++)
            p_info->start_coords[k] = intersect_pt[k];

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

      if (!cs_glob_lagr_params->clogging && !cs_glob_lagr_params->resuspension ) {
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

      move_particle = CS_LAGR_PART_MOVE_ON;
      particle_state = CS_LAGR_PART_TO_SYNC;
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                cs_glob_mesh->b_face_cells[face_id] + 1);
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_IN_FLOW);

      for (k = 0; k < 3; k++)
        p_info->start_coords[k] = intersect_pt[k];

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
    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id] + 1);

    for (k = 0; k < 3; k++)
      p_info->start_coords[k] = intersect_pt[k];

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
    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id] + 1);

    for (k = 0; k < 3; k++)
      p_info->start_coords[k] = intersect_pt[k];

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
    const cs_lnum_t n_layers = p_am->count[0][CS_LAGR_TEMPERATURE];
    const cs_real_t *particle_temp
      = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_TEMPERATURE);

    cs_real_t  temp_ext_part = particle_temp[n_layers - 1];
    cs_real_t  tprenc_icoal  = tprenc[particle_coal_number - 1];
    cs_real_t  visref_icoal  = visref[particle_coal_number - 1];
    cs_real_t  enc1_icoal    = enc1[particle_coal_number - 1];
    cs_real_t  enc2_icoal    = enc2[particle_coal_number - 1];

    if (temp_ext_part > tprenc_icoal+tkelvi) {

      /* Coal viscosity*/
      tmp = (  (1.0e7*enc1_icoal)
             / ((temp_ext_part-150.e0-tkelvi)*(temp_ext_part-150.e0-tkelvi)))
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
        particles->n_part_fou += 1;
        particles->weight_fou += particle_stat_weight;

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
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                cs_glob_mesh->b_face_cells[face_id] + 1);

      for (k = 0; k < 3; k++)
        p_info->start_coords[k] = intersect_pt[k];

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

  /* Return pointer */

  *p_move_particle = move_particle;
  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;

  /* Particulate boundary mass flux (contribution of rolling particles
     will be added  at the end of their movement) */

  if (p_am->size[CS_LAGR_DEPOSITION_FLAG] > 0) {
    cs_lnum_t deposition_flag
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG);
    if (deposition_flag == CS_LAGR_PART_DEPOSITED && part_b_mass_flux != NULL)
      part_b_mass_flux[face_id]
        += particle_stat_weight * particle_mass / face_area;
  }

  /* FIXME: Post-treatment not yet implemented... */

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

  assert(particle_state > -1);

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Move a particle as far as possible while remaining on a given rank.
 *
 * parameters:
 *   particle                 <-> pointer to particle data
 *   p_am                     <-- particle attribute map
 *   displacement_step_id     <-- id of displacement step
 *   scheme_order             <-- current order of the scheme used for
 *                                the lagrangian algorithm (1 or 2)
 *   failsafe_mode            <-- with (0) / without (1) failure capability
 *   p_n_failed_particles     <-> number of failed particles
 *   p_failed_particle_weight <-> stat. weight of the failed particles
 *
 * returns:
 *   a state associated to the status of the particle (treated, to be deleted,
 *   to be synchonised)
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_local_propagation(void                           *particle,
                   const cs_lagr_attribute_map_t  *p_am,
                   int                             displacement_step_id,
                   int                             scheme_order,
                   int                             failsafe_mode,
                   cs_real_t                       boundary_stat[],
                   cs_lnum_t                      *p_n_failed_particles,
                   cs_real_t                      *p_failed_particle_weight,
                   const cs_lnum_t                *iensi3,
                   const cs_lnum_t                *inbr,
                   const cs_lnum_t                *inbrbd,
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
                   cs_real_t                      *part_b_mass_flux,
                   const cs_field_t               *u,
                   cs_real_t                       energt[],
                   const cs_real_t                 tprenc[],
                   const cs_real_t                 visref[],
                   const cs_real_t                 enc1[],
                   const cs_real_t                 enc2[],
                   cs_real_t                       tkelvi)
{
  cs_lnum_t  i, k;
  cs_real_t  depl[3];
  cs_real_t  null_yplus;

  cs_lnum_t  *neighbor_face_id;
  cs_real_t  *particle_yplus;

  cs_lnum_t  n_loops = displacement_step_id;
  cs_lnum_t  move_particle = CS_LAGR_PART_MOVE_ON;
  cs_lnum_t  particle_state = CS_LAGR_PART_TO_SYNC;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  const cs_lagr_param_t *lagr_params = cs_glob_lagr_params;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;
  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;

  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight = *p_failed_particle_weight;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  cs_real_t  *particle_coord
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);

  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

  for (k = 0; k < 3; k++)
    depl[k] = particle_coord[k] - p_info->start_coords[k];

  if (fabs(depl[0]) < 1e-15 && fabs(depl[1]) < 1e-15 && fabs(depl[2]) < 1e-15) {
    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_TREATED;
  }

  if (lagr_params->deposition > 0) {

    neighbor_face_id
      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
    particle_yplus
      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_YPLUS);

    /* Remove contribution from rolling particles to boundary mass flux
       at the beginning of their movement, so as to later add it at the
       end of their movement */

    if (displacement_step_id == 0 && move_particle == CS_LAGR_PART_MOVE_ON) {

      if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
          == CS_LAGR_PART_ROLLING) {

        assert(*neighbor_face_id > -1);

        if (part_b_mass_flux != NULL) {
          cs_real_t cur_stat_weight
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
          cs_real_t cur_mass
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
          cs_real_t face_area
            = fvq->b_face_surf[*neighbor_face_id];

          part_b_mass_flux[*neighbor_face_id]
            -= cur_stat_weight * cur_mass / face_area;
        }

      }

    }

  }
  else {
    neighbor_face_id = NULL;
    null_yplus = 0;
    particle_yplus = &null_yplus;  /* allow tests even without particle y+ */
  }

  /*  particle_state is defined at the top of this file */

  while (move_particle == CS_LAGR_PART_MOVE_ON) {

    cs_lnum_t  cur_cell_id
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM) - 1;
    cs_lnum_t  old_face_num = p_info->last_face_num;

    assert(cur_cell_id < mesh->n_cells);
    assert(cur_cell_id > -1);

    cs_lnum_t  start = cell_face_idx[cur_cell_id];
    cs_lnum_t  end =  cell_face_idx[cur_cell_id+1];

    n_loops++;

    if (n_loops > _max_propagation_loops) { /* Manage error */

      _manage_error(failsafe_mode,
                    particle,
                    p_am,
                    CS_LAGR_TRACKING_ERR_MAX_LOOPS,
                    &n_failed_particles,
                    &failed_particle_weight);

      move_particle  = CS_LAGR_PART_MOVE_OFF;
      particle_state = CS_LAGR_PART_ERR;

    }

    /* Treatment for depositing particles */

    if (lagr_params->deposition > 0 && *particle_yplus < 0.) {

      _test_wall_cell(particle, p_am, visc_length, dlgeo,
                      particle_yplus, neighbor_face_id);

      if (*particle_yplus < 100.) {

        cs_real_t flow_velo_x, flow_velo_y, flow_velo_z;

        assert(u->interleaved);

        flow_velo_x = u->val[cur_cell_id*3];
        flow_velo_y = u->val[cur_cell_id*3 + 1];
        flow_velo_z = u->val[cur_cell_id*3 + 2];

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

    cs_lnum_t exit_face = 0; /* > 0 for interior faces,
                                < 0 for boundary faces */

    double adist_min = HUGE_VAL;

    for (i = start; i < end && move_particle == CS_LAGR_PART_MOVE_ON; i++) {

      cs_lnum_t  face_id, vtx_start, vtx_end, n_vertices;
      const cs_lnum_t  *face_connect;

      cs_lnum_t  face_num = cell_face_lst[i];

      if (face_num == old_face_num)
        continue;

      if (face_num > 0) {

        /* Interior face */

        face_id = face_num - 1;
        vtx_start = mesh->i_face_vtx_idx[face_id];
        vtx_end = mesh->i_face_vtx_idx[face_id+1];
        n_vertices = vtx_end - vtx_start;

        face_connect = mesh->i_face_vtx_lst + vtx_start;

      }
      else {

        assert(face_num < 0);

        /* Boundary faces */

        face_id = -face_num - 1;
        vtx_start = mesh->b_face_vtx_idx[face_id];
        vtx_end = mesh->b_face_vtx_idx[face_id+1];
        n_vertices = vtx_end - vtx_start;

        face_connect = mesh->b_face_vtx_lst + vtx_start;

      }

      /*
        adimensional distance estimation of face intersection
        (-1 if no chance of intersection)
      */

      double adist = _intersect_face(face_num,
                                     n_vertices,
                                     face_connect,
                                     particle,
                                     p_am);
      if (adist > -0.5 && adist < adist_min) {
        exit_face = face_num;
        adist_min = adist;
        if (adist < 1.0)
          break;
      }

    }

    if (exit_face == 0) {
      move_particle =  CS_LAGR_PART_MOVE_OFF;
      particle_state = CS_LAGR_PART_TREATED;
    }

    else if (exit_face > 0) { /* Particle moves to the neighbor cell
                                 through the current face "face_num" */

      cs_lnum_t face_id = exit_face - 1;

      cs_lnum_t  c_id1 = mesh->i_face_cells[face_id][0];
      cs_lnum_t  c_id2 = mesh->i_face_cells[face_id][1];

      p_info->last_face_num = exit_face;

      if (cur_cell_id == c_id1) {
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                  c_id2 + 1);
        cur_cell_id = c_id2;
      }

      else {
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                  c_id1 + 1);
        cur_cell_id = c_id1;
      }

      /* Specific treatment for the particle deposition model */

      if (lagr_params->deposition > 0 && *particle_yplus < 100.) {

        cs_real_t x_p_q = particle_coord[0] - p_info->start_coords[0];
        cs_real_t y_p_q = particle_coord[1] - p_info->start_coords[1];
        cs_real_t z_p_q = particle_coord[2] - p_info->start_coords[2];

        cs_real_t face_normal[3], face_cog[3];

        for (k = 0; k < 3; k++) {
          face_normal[k] = fvq->i_face_normal[3*face_id+k];
          face_cog[k] = fvq->i_face_cog[3*face_id+k];
        }

        cs_real_t aa =   x_p_q * face_normal[0]
                       + y_p_q * face_normal[1]
                       + z_p_q * face_normal[2];


        cs_real_t bb = (  face_normal[0] * face_cog[0]
                        + face_normal[1] * face_cog[1]
                        + face_normal[2] * face_cog[2]
                        - face_normal[0] * p_info->start_coords[0]
                        - face_normal[1] * p_info->start_coords[1]
                        - face_normal[2] * p_info->start_coords[2]) / aa;

        cs_real_t xk =  p_info->start_coords[0] + bb * x_p_q;
        cs_real_t yk =  p_info->start_coords[1] + bb * y_p_q;
        cs_real_t zk =  p_info->start_coords[2] + bb * z_p_q;

        cs_real_t *xyzcen = fvq->cell_cen;

        particle_coord[0] = xk + 1e-8 * (xyzcen[3*cur_cell_id] - xk);
        particle_coord[1] = yk + 1e-8 * (xyzcen[3*cur_cell_id + 1] - yk);
        particle_coord[2] = zk + 1e-8 * (xyzcen[3*cur_cell_id + 2] - zk);

      }

      /* Particle changes rank */

      if (cur_cell_id >= mesh->n_cells) {

        particle_state = CS_LAGR_PART_TO_SYNC;
        move_particle = CS_LAGR_PART_MOVE_OFF;

        /* Specific treatment for the particle deposition model */

        if (lagr_params->deposition > 0 && *particle_yplus < 100.) {

          /* Marking of particles */

          *particle_yplus = - *particle_yplus;

          /*Saving of dot product */

          cs_real_t bdy_normal[3];

          for (k = 0; k < 3; k++)
            bdy_normal[k] = fvq->b_face_normal[(*neighbor_face_id)*3 + k];

          cs_real_t area  = _get_norm(bdy_normal);

          cs_real_t  face_norm[3] = {bdy_normal[0]/area,
                                     bdy_normal[1]/area,
                                     bdy_normal[2]/area};

          particle_velocity_seen[0]
            =   particle_velocity_seen[0] * face_norm[0]
              + particle_velocity_seen[1] * face_norm[1]
              + particle_velocity_seen[2] * face_norm[2];

        }

      } /* end of case where particle changes rank */

      else if (lagr_params->deposition > 0) {

        /* Specific treatment for the particle deposition model */

        cs_lnum_t save_close_face_id = *neighbor_face_id;
        cs_real_t save_yplus = *particle_yplus;

        /* Wall cell detection */

        _test_wall_cell(particle, p_am, visc_length, dlgeo,
                        particle_yplus, neighbor_face_id);

        if (save_yplus < 100.) {

          /* Second test with the new particle position */

          _test_wall_cell(particle, p_am, visc_length, dlgeo,
                          particle_yplus, neighbor_face_id);

          if (*particle_yplus < 100.e0) {

            assert(u->interleaved);

            cs_real_t flow_velo_x = u->val[cur_cell_id*3];
            cs_real_t flow_velo_y = u->val[cur_cell_id*3 + 1];
            cs_real_t flow_velo_z = u->val[cur_cell_id*3 + 2];

            /* The particle is still in the boundary layer */

            const cs_real_t *old_bdy_normal
              = fvq->b_face_normal + 3*save_close_face_id;

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

      } /* End of case for deposition */

    }
    else if (exit_face < 0) { /* Particle moves to the boundary
                                 through the current face "face_num" */

      cs_lnum_t face_num = -exit_face;

      /* particle / boundary condition interaction
         1 - modify particle cell_num : 0 or boundary_cell_num
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
        = _boundary_treatment(_particle_set,
                              particle,
                              face_num,
                              boundary_stat,
                              bdy_conditions->b_face_zone_num[face_num-1]-1,
                              &move_particle,
                              &n_failed_particles,
                              &failed_particle_weight,
                              iensi3,
                              inbr,
                              inbrbd,
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
                              part_b_mass_flux,
                              energt,
                              tprenc,
                              visref,
                              enc1,
                              enc2,
                              tkelvi);

      if (scheme_order == 2)
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_SWITCH_ORDER_1,
                                  CS_LAGR_SWITCH_ON);

      assert(   move_particle == CS_LAGR_PART_MOVE_ON
             || move_particle == CS_LAGR_PART_MOVE_OFF);

      p_info->last_face_num = -face_num;

    } /* end if exit_face < 0 */

  } /* End of while : local displacement */

  /* Add contribution from rolling particles to boundary mass flux
     at the beginning of their movement. */

  if (lagr_params->deposition > 0) {

    if (   move_particle == CS_LAGR_PART_MOVE_OFF
        && particle_state != CS_LAGR_PART_OUT
        && particle_state != CS_LAGR_PART_ERR) {

      if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
           == CS_LAGR_PART_ROLLING) {

        assert(*neighbor_face_id > -1);

        if (part_b_mass_flux != NULL) {
          cs_real_t cur_stat_weight
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
          cs_real_t cur_mass
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
          cs_real_t face_area = fvq->b_face_surf[*neighbor_face_id];

          part_b_mass_flux[*neighbor_face_id]
            += cur_stat_weight * cur_mass / face_area;
        }

      }

    }

  }

  assert(   move_particle != CS_LAGR_PART_MOVE_ON
         || particle_state != CS_LAGR_PART_TO_SYNC);

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

    /* We wait for posting all receives
       (often recommended in the past, apparently not anymore) */

#if 0
    MPI_Barrier(cs_glob_mpi_comm)
#endif

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

  const size_t tot_extents = lag_halo->extents;

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
                  _cs_mpi_particle_type,
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
    MPI_Barrier(cs_glob_mpi_comm);
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
                  _cs_mpi_particle_type,
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
 *   particles        <-- set of particle data to update
 *----------------------------------------------------------------------------*/

static void
_update_particle_set(cs_lnum_t                n_recv_particles,
                     cs_lagr_halo_t          *lag_halo,
                     cs_lagr_particle_set_t  *particles)
{
  cs_lnum_t  i, new_id;
  const size_t extents = particles->p_am->extents;
  const unsigned char *recv_buf = lag_halo->recv_buf;

  assert(extents == lag_halo->extents);

  for (i = 0; i < n_recv_particles; i++) {

    const unsigned char *new_part = recv_buf + extents*i;

    new_id = particles->first_free_id;

    /* Note: used and free ids for current and previous data are the same */

    if (particles->used_id[new_id].next_id != -1) {
      particles->first_free_id = particles->used_id[new_id].next_id;
    }
    else {
      particles->first_free_id = new_id + 1 ;
      particles->used_id[particles->first_free_id].next_id = -1;
    }

    /* Add new_particle at the beginning of the "used list"
       Update first_used_id */

    if (particles->first_used_id != -1)
      particles->used_id[particles->first_used_id].prev_id = new_id;

    memcpy(particles->p_buffer + extents*new_id,
           new_part,
           extents);

    particles->used_id[new_id].prev_id = -1;
    particles->used_id[new_id].next_id = particles->first_used_id;

    particles->first_used_id = new_id;
  }
}

/*----------------------------------------------------------------------------
 * Synchronize particle displacement over the ranks.
 *
 * parameters:
 *  lag_halo  <--  pointer to a cs_lagr_halo_t structure
 *  particles <--  set of particles to update
 *----------------------------------------------------------------------------*/

static void
_sync_particle_set(cs_lagr_halo_t           *lag_halo,
                   cs_interface_set_t       *face_ifs,
                   cs_lagr_particle_set_t   *particles)
{
  cs_lnum_t  i, j, k, tr_id, rank, shift, ghost_id;
  cs_real_t matrix[3][4];

  cs_lnum_t  n_recv_particles = 0;

  const size_t extents = particles->p_am->extents;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const fvm_periodicity_t *periodicity = mesh->periodicity;

  for (i = 0; i < halo->n_c_domains; i++){
    n_recv_particles += lag_halo->recv_count[i];
    lag_halo->send_count[i] = 0;
  }

  /* Fill send_buf for particle set */

  for (i = 0, j = particles->first_used_id; i < particles->n_particles; i++) {

    cs_lnum_t next_id = particles->used_id[j].next_id;

    if (_get_tracking_info(particles, j)->state == CS_LAGR_PART_TO_SYNC) {

      ghost_id =   cs_lagr_particles_get_lnum(particles, j, CS_LAGR_CELL_NUM)
                 - halo->n_local_elts - 1;
      rank = lag_halo->rank[ghost_id];
      tr_id = lag_halo->transform_id[ghost_id];

      cs_lagr_particles_set_lnum(particles, j, CS_LAGR_CELL_NUM,
                                 lag_halo->dist_cell_num[ghost_id]);

      shift = lag_halo->send_shift[rank] + lag_halo->send_count[rank];

      /* Update if needed last_face_num */

      if (tr_id >= 0) { /* Same initialization as in previous algorithm */

        _tracking_info(particles, j)->last_face_num = 0;

      }

      else {

        if (cs_glob_n_ranks > 1) {

          assert(face_ifs != NULL);

          int  distant_rank;
          cs_lnum_t n_entities, id;
          const cs_lnum_t *local_num, *dist_num;

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

            id = cs_search_binary
                   (n_entities,
                    _get_tracking_info(particles, j)->last_face_num  - 1,
                    local_num);

            if (id == -1)
              bft_error(__FILE__, __LINE__, 0,
                        _(" Cannot find the relative distant face num.\n"));

            dist_num = cs_interface_get_match_ids(interface);

            _tracking_info(particles, j)->last_face_num = dist_num[id] + 1;

          }

        }
      }

      /* Periodicity treatment.
         Note that for purposes such as postprocessing of trajectories,
         we also apply periodicity transformations to values at the previous
         time step, so that previous/current data is consistent relative to the
         new position */

      if (tr_id >= 0) {

        /* Transform coordinates */

        fvm_periodicity_type_t  perio_type
          = fvm_periodicity_get_type(periodicity, tr_id);

        int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, tr_id);

        fvm_periodicity_get_matrix(periodicity, rev_id, matrix);

        /* Apply transformation to the coordinates in any case */

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              cs_lagr_particles_attr(particles, j,
                                                     CS_LAGR_COORDS));

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              _tracking_info(particles, j)->start_coords);

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              cs_lagr_particles_attr_n(particles, j, 1,
                                                       CS_LAGR_COORDS));

        /* Apply rotation to velocity vectors in case of rotation */

        if (perio_type >= FVM_PERIODICITY_ROTATION) {

          /* Rotation of the velocity */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(particles, j,
                                                        CS_LAGR_VELOCITY));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr_n(particles, j, 1,
                                                          CS_LAGR_VELOCITY));

          /* Rotation of the velocity seen */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(particles, j,
                                                        CS_LAGR_VELOCITY_SEEN));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr_n(particles, j, 1,
                                                          CS_LAGR_VELOCITY_SEEN));

        } /* Specific treatment in case of rotation for the velocities */

      } /* End of periodicity treatment */

      memcpy(lag_halo->send_buf + extents*shift,
             particles->p_buffer + extents*j,
             extents);

      lag_halo->send_count[rank] += 1;

      /* Pick out the particle from the particle_set */

      _remove_particle(particles, j);

    } /* TO_SYNC */

    j = next_id;

  } /* End of loop on particles */

  /* Exchange particles, then update set */

  _exchange_particles(halo, lag_halo);

  _update_particle_set(n_recv_particles, lag_halo, particles);
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
  cs_lagr_particle_set_t  *particles = _particle_set;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_halo_t  *lag_halo = builder->halo;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;

  /* Check this because this assumption is done. */

#if 0
  bft_printf("\n Particle set before sync\n");
  cs_lagr_particle_set_dump(particles);
#endif

  /* Initialization */

  for (i = 0; i < halo->n_c_domains; i++) {
    lag_halo->send_count[i] = 0;
    lag_halo->recv_count[i] = 0;
  }

  /* Loop on particles to count number of particles to send on each rank */

  for (i = 0, j = particles->first_used_id; i < particles->n_particles; i++) {

    cs_lnum_t next_id = particles->used_id[j].next_id;

    if (_get_tracking_info(particles, j)->state == CS_LAGR_PART_TO_SYNC) {

      ghost_id =   cs_lagr_particles_get_lnum(particles, j, CS_LAGR_CELL_NUM)
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

  _sync_particle_set(lag_halo, builder->face_ifs, particles);

  particles->n_particles += delta_particles;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id == 1) {
    bft_printf(" delta part = %d\n", delta_particles);
    bft_printf(" particles -> npart = %d\n", particles->n_particles);
  }
#endif

#if 0
  bft_printf("\n Particle set after sync\n");
  cs_lagr_particle_set_dump(particles);
#endif

  if (delta_particles > particles->n_particles_max - particles->n_particles)
    _resize_particle_set(particles,
                         particles->n_particles + delta_particles);
}

/*----------------------------------------------------------------------------
 * Prepare for particle movement phase
 *
 * parameters:
 *   particles <-> pointer to particles structure
 *----------------------------------------------------------------------------*/

static void
_initialize_displacement(cs_lagr_particle_set_t  *particles)
{
  cs_lnum_t  i;

  const cs_lagr_param_t *lagr_params = cs_glob_lagr_params;

  const cs_lagr_attribute_map_t  *am = particles->p_am;

  /* Prepare tracking info */

  for (i = 0; i < particles->n_particles; i++) {

    particles->used_id[i].prev_id = i-1;
    particles->used_id[i].next_id = i+1;

    cs_lnum_t cur_part_cell_num
      = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_NUM);
    if (cur_part_cell_num < 0)
      _tracking_info(particles, i)->state = CS_LAGR_PART_STUCK;
    else if (cur_part_cell_num == 0)
      _tracking_info(particles, i)->state = CS_LAGR_PART_TO_DELETE;
    else
      _tracking_info(particles, i)->state = CS_LAGR_PART_TO_SYNC;

    _tracking_info(particles, i)->last_face_num = 0;

    assert(   cs_lagr_particles_get_lnum(particles, i, CS_LAGR_SWITCH_ORDER_1)
           != 999);

    /* Coordinates of the particle */

    cs_real_t *prv_part_coord
      = cs_lagr_particles_attr_n(particles, i, 1, CS_LAGR_COORDS);

    _tracking_info(particles, i)->start_coords[0] = prv_part_coord[0];
    _tracking_info(particles, i)->start_coords[1] = prv_part_coord[1];
    _tracking_info(particles, i)->start_coords[2] = prv_part_coord[2];

    /* Data needed if the deposition model is activated */
    if (   lagr_params->deposition <= 0
        && am->size[CS_LAGR_DEPOSITION_FLAG] > 0)
      cs_lagr_particles_set_lnum(particles, i, CS_LAGR_DEPOSITION_FLAG,
                                 CS_LAGR_PART_IN_FLOW);

  }

  if (particles->n_particles > 0) {
    particles->first_used_id = 0;
    particles->used_id[0].prev_id = -1;
    particles->used_id[particles->n_particles - 1].next_id = -1;
  }
  else
    particles->first_used_id = -1;

  particles->first_free_id = particles->n_particles;
  particles->used_id[particles->first_free_id].next_id = -1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Particle set after %s\n", __func__);
  cs_lagr_particle_set_dump(particles);
#endif
}

/*----------------------------------------------------------------------------
 * Update particle set structures: transform list to array.
 *
 * parameters:
 *   particles <-> pointer to particle set
 *----------------------------------------------------------------------------*/

static void
_finalize_displacement(cs_lagr_particle_set_t  *particles)
{
  cs_lnum_t  i, j;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  const cs_lnum_t  n_particles = particles->n_particles;

  unsigned char *swap_buffer;
  size_t swap_buffer_size = p_am->extents * ((size_t)n_particles);

  BFT_MALLOC(swap_buffer, swap_buffer_size, unsigned char);

  /* For particle state */

  for (i = 0, j = particles->first_used_id; i < n_particles; i++) {

    cs_lnum_t cur_part_state = _get_tracking_info(particles, j)->state;

    if  (cur_part_state == CS_LAGR_PART_ERR)
      cs_lagr_particles_set_lnum(particles, j, CS_LAGR_CELL_NUM, 0);

#if defined(DEBUG) || !defined(NDEBUG)

    assert(cur_part_state != CS_LAGR_PART_TO_DELETE);
    assert(cur_part_state != CS_LAGR_PART_OUT);
    assert(cur_part_state != CS_LAGR_PART_ERR);
    assert(cur_part_state != CS_LAGR_PART_TO_SYNC);

#endif

    memcpy(swap_buffer + p_am->extents*i,
           particles->p_buffer + p_am->extents*j,
           p_am->extents);

    j = particles->used_id[j].next_id;

  }

  memcpy(particles->p_buffer, swap_buffer, swap_buffer_size);

  BFT_FREE(swap_buffer);

  /* Shift from 0-based to 1-based (temporary) */

  if (p_am->count[0][CS_LAGR_NEIGHBOR_FACE_ID] > 0) {
    for (i = 0; i < n_particles; i++) {
      const cs_lnum_t n_face_num
        = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_NEIGHBOR_FACE_ID) + 1;
        cs_lagr_particles_set_lnum(particles, i, CS_LAGR_NEIGHBOR_FACE_ID,
                                   n_face_num);
    }
  }

  /* Also mark list as not used */

  particles->first_used_id = -1;
  particles->first_free_id = -1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Particle set after %s\n", __func__);
  cs_lagr_particle_set_dump(particles);
#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointers to particle attributes
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   dim_ipepa   --> dimensions for ipepa pointer
 *   dim_peta    --> dimensions for pepa pointer
 *   dim_eptp    --> dimensions for eptp pointer
 *   dim_eptpa   --> dimensions for eptpa pointer
 *   dim_ptsvar  --> dimensions for ptsvar pointer
 *   p_ipeta     --> ipepa pointer
 *   p_pepa      --> pepa pointer
 *   p_eptp      --> eptp pointer
 *   p_eptpa     --> eptpa pointer
 *   p_ptsvar    --> ptsvar pointer
 *   p_nbpart    --> nbpart pointer
 *   p_dnbpar    --> dnbpar pointer
 *   p_nbpout    --> nbpout pointer
 *   p_dnbpou    --> dnbpou pointer
 *   p_nbperr    --> nbperr pointer
 *   p_dnbper    --> dnbper  pointer
 *   p_nbpdep    --> nbpdep pointer
 *   p_dnbdep    --> dnbdep  pointer
 *   p_npencr    --> npencr pointer
 *   p_dnpenc    --> dnpenc  pointer
 *
 * returns:
 *   pointer to the field structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_f_lagr_pointers(int          dim_ipepa[2],
                   int          dim_pepa[2],
                   int          dim_eptp[2],
                   int          dim_eptpa[2],
                   int          dim_ptsvar[2],
                   cs_int_t   **p_ipeta,
                   cs_real_t  **p_peta,
                   cs_real_t  **p_eptp,
                   cs_real_t  **p_eptpa,
                   cs_real_t  **p_ptsvar,
                   cs_int_t   **p_nbpart,
                   cs_real_t  **p_dnbpar,
                   cs_int_t   **p_nbpout,
                   cs_real_t  **p_dnbpou,
                   cs_int_t   **p_nbperr,
                   cs_real_t  **p_dnbper,
                   cs_int_t   **p_nbpdep,
                   cs_real_t  **p_dnbdep,
                   cs_int_t   **p_npencr,
                   cs_real_t  **p_dnpenc)
{
  cs_lagr_particle_set_t  *particles = _particle_set;
  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  const size_t extents = p_am->extents;
  const size_t max_data_size
    = p_am->extents * ((size_t)particles->n_particles_max);

  dim_ipepa[0] = extents / sizeof(cs_int_t);
  dim_ipepa[1] = max_data_size / sizeof(cs_int_t);

  dim_pepa[0] = extents / sizeof(cs_real_t);
  dim_pepa[1] = max_data_size / sizeof(cs_real_t);

  dim_eptp[0] = extents / sizeof(cs_real_t);
  dim_eptp[1] = max_data_size / sizeof(cs_real_t);

  dim_eptpa[0] = extents / sizeof(cs_real_t);
  dim_eptpa[1] = max_data_size / sizeof(cs_real_t);

  if (cs_glob_lagr_params->t_order > 1) {
    dim_ptsvar[0] = extents / sizeof(cs_real_t);
    dim_ptsvar[1] = max_data_size / sizeof(cs_real_t);
  }
  else {
    dim_ptsvar[0] = 0;
    dim_ptsvar[1] = 0;
  }

  *p_ipeta = (cs_int_t *)particles->p_buffer;
  *p_peta = (cs_real_t *)particles->p_buffer;

  *p_eptp = (cs_real_t *)particles->p_buffer;
  *p_eptpa = (cs_real_t *)(  particles->p_buffer
                           + p_am->displ[1][CS_LAGR_COORDS]
                           - p_am->displ[0][CS_LAGR_COORDS]);

  if (cs_glob_lagr_params->t_order > 1) {
    assert(p_am->displ[0][CS_LAGR_MASS] > -1);
    *p_ptsvar = (cs_real_t *)(  particles->p_buffer
                              + p_am->source_term_displ[CS_LAGR_MASS]
                              - p_am->displ[0][CS_LAGR_MASS]);
    assert(p_am->displ[0][CS_LAGR_TURB_STATE_1] > 0);
  }
  else
    *p_ptsvar = NULL;

  *p_nbpart = &(particles->n_particles);
  *p_dnbpar = &(particles->weight);
  *p_nbpout = &(particles->n_part_out);
  *p_dnbpou = &(particles->weight_out);
  *p_nbperr = &(particles->n_failed_part);
  *p_dnbper = &(particles->weight_failed);
  *p_nbpdep = &(particles->n_part_dep);
  *p_dnbdep = &(particles->weight_dep);
  *p_npencr = &(particles->n_part_fou);
  *p_dnpenc = &(particles->weight_fou);
}

/*----------------------------------------------------------------------------
 * Copy attributes from one particle to another
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   dest <-- number (1-based) of destination particle
 *   src  <-- number (1-based) of source particle
 *----------------------------------------------------------------------------*/

void
cs_f_lagr_part_copy(cs_lnum_t  dest,
                    cs_lnum_t  src)
{
  cs_lagr_particle_set_t  *particles = _particle_set;
  memcpy(particles->p_buffer + particles->p_am->extents*(dest-1),
         particles->p_buffer + particles->p_am->extents*(src-1),
         particles->p_am->extents);
  cs_lagr_particles_set_real(particles, (dest-1), CS_LAGR_RANDOM_VALUE,
                             rand()/RAND_MAX);
}

/*----------------------------------------------------------------------------
 * Copy current attributes to previous attributes.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   pn <-- number (1-based) of destination particle
 *   src  <-- number (1-based) of source particle
 *----------------------------------------------------------------------------*/

void
cs_f_lagr_current_to_previous(cs_lnum_t  pn)
{
  cs_lagr_particle_set_t  *particles = _particle_set;
  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  unsigned char *p_buf = particles->p_buffer + p_am->extents*(pn-1);

  for (cs_lagr_attribute_t attr = 0;
       attr < CS_LAGR_N_ATTRIBUTES;
       attr++) {
    if (p_am->count[1][attr] > 0 && p_am->count[0][attr] > 0) {
      memcpy(p_buf + p_am->displ[1][attr],
             p_buf + p_am->displ[0][attr],
             p_am->size[attr]);
    }
  }
  *((cs_lnum_t *)(  p_buf + p_am->displ[1][CS_LAGR_RANK_ID])) = cs_glob_rank_id;
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate cs_lagr_particle_set_t structure and initialize useful buffers
 * and indexes
 *
 * parameters:
 *   nordre          <--  time algorithm order (1 or 2)
 *   iphyla          <--  kind of physics used for the lagrangian approach
 *   nvls            <--  number of user-defined variables
 *   nbclst          <--  number of stat. class to study sub-set of particles
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagbeg, LAGBEG)(const cs_int_t    *nordre,
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
                          cs_int_t          *jisor,
                          cs_int_t          *jisora,
                          cs_int_t          *jirka,
                          cs_int_t          *jord1,
                          cs_int_t          *jrval,
                          cs_int_t          *jrpoi,
                          cs_int_t          *jrtsp,
                          cs_int_t          *jdp,
                          cs_int_t          *jmp,
                          cs_int_t          *jxp,
                          cs_int_t          *jyp,
                          cs_int_t          *jzp,
                          cs_int_t          *jup,
                          cs_int_t          *jvp,
                          cs_int_t          *jwp,
                          cs_int_t          *juf,
                          cs_int_t          *jvf,
                          cs_int_t          *jwf,
                          cs_int_t          *jtaux,
                          cs_int_t           jbx1[3],
                          cs_int_t           jtsup[3],
                          cs_int_t           jtsuf[3],
                          cs_int_t          *jryplu,
                          cs_int_t          *jrinpf,
                          cs_int_t          *jdfac,
                          cs_int_t          *jimark,
                          cs_int_t          *jtp,
                          cs_int_t           jhp[],
                          cs_int_t          *jtf,
                          cs_int_t          *jmwat,
                          cs_int_t           jmch[],
                          cs_int_t           jmck[],
                          cs_int_t          *jcp,
                          cs_int_t          *jrdck,
                          cs_int_t          *jrd0p,
                          cs_int_t          *jinch,
                          cs_int_t           jrhock[],
                          cs_int_t          *jreps,
                          cs_int_t          *jdepo,
                          cs_int_t          *jnbasg,
                          cs_int_t          *jnbasp,
                          cs_int_t          *jfadh,
                          cs_int_t          *jmfadh,
                          cs_int_t          *jndisp,
                          cs_int_t          *jclst,
                          cs_int_t          *jvls)
{
  cs_lnum_t  i;
  cs_mesh_t  *mesh = cs_glob_mesh;

  int n_layers = 0;
  int ieptp_loc_count = 0;
  int pepa_loc_add = 1000; /* should be abore any j* pointer if used */

  cs_lnum_t attr_keys[CS_LAGR_N_ATTRIBUTES][3];

  /* Initialize global parameter relative to the lagrangian module */

  _lagr_param.physical_model = *iphyla;

  if (_lagr_param.physical_model == 2)
    _lagr_param.n_temperature_layers = *nlayer;
  else
    _lagr_param.n_temperature_layers = 1;

  _lagr_param.deposition = *idepst;
  _lagr_param.roughness = *irough;
  _lagr_param.resuspension = *ireent;
  _lagr_param.clogging = *iclogst;

  _lagr_param.n_user_variables = *nvls;
  _lagr_param.n_stat_classes = *nbclst;

  _lagr_param.t_order = *nordre;

  /* Set indexes */

  for (i = 0; i < CS_LAGR_N_ATTRIBUTES; i++) {
    attr_keys[i][0] = EPTP; /* default */
    attr_keys[i][1] = 0;
    attr_keys[i][2] = 0;
  }

  /* Special case:
     jisor is first index in ipepa, but cell number is also
     needed in indep for previous time step;
     so we cheat here and assign this attribute to "ieptp".
  */

  assert(*jisor == 1);               /* for future mapping of Fortran to C */

  attr_keys[CS_LAGR_CELL_NUM][0] = IEPTP;
  attr_keys[CS_LAGR_CELL_NUM][1] = ++ieptp_loc_count;

  attr_keys[CS_LAGR_RANK_ID][0] = IPRKID;
  attr_keys[CS_LAGR_RANK_ID][1] = 1;

  /* Other attributes */

  attr_keys[CS_LAGR_SWITCH_ORDER_1][0] = IPEPA;
  attr_keys[CS_LAGR_SWITCH_ORDER_1][1] = *jord1;

  attr_keys[CS_LAGR_RANDOM_VALUE][0] = PEPA;
  attr_keys[CS_LAGR_RANDOM_VALUE][1] = *jrval;

  attr_keys[CS_LAGR_STAT_WEIGHT][0] = PEPA;
  attr_keys[CS_LAGR_STAT_WEIGHT][1] = *jrpoi;

  attr_keys[CS_LAGR_RESIDENCE_TIME][0] = PEPA;
  attr_keys[CS_LAGR_RESIDENCE_TIME][1] = *jrtsp;

  attr_keys[CS_LAGR_MASS][0] = EPTP_TS;
  attr_keys[CS_LAGR_MASS][1] = *jmp;

  attr_keys[CS_LAGR_DIAMETER][0] = EPTP_TS;
  attr_keys[CS_LAGR_DIAMETER][1] = *jdp;

  attr_keys[CS_LAGR_COORDS][1] = *jxp;
  attr_keys[CS_LAGR_COORDS][2] = 3;

  attr_keys[CS_LAGR_VELOCITY][1] = *jup;
  attr_keys[CS_LAGR_VELOCITY][2] = 3;

  attr_keys[CS_LAGR_VELOCITY_SEEN][1] = *juf;
  attr_keys[CS_LAGR_VELOCITY_SEEN][2] = 3;

  if (_lagr_param.t_order > 1) {
    attr_keys[CS_LAGR_TAUP_AUX][0] = PEPA;
    attr_keys[CS_LAGR_TAUP_AUX][1] = ++pepa_loc_add;

    attr_keys[CS_LAGR_TURB_STATE_1][0] = PEPA;
    attr_keys[CS_LAGR_TURB_STATE_1][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_TURB_STATE_1][2] = 3;

    attr_keys[CS_LAGR_PRED_VELOCITY][0] = PEPA;
    attr_keys[CS_LAGR_PRED_VELOCITY][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_PRED_VELOCITY][2] = 3;

    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][0] = PEPA;
    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][1] = ++pepa_loc_add;
    attr_keys[CS_LAGR_PRED_VELOCITY_SEEN][2] = 3;
  }

  attr_keys[CS_LAGR_DEPOSITION_FLAG][0] = IPEPA;
  attr_keys[CS_LAGR_DEPOSITION_FLAG][1] = *jdepo;

  attr_keys[CS_LAGR_YPLUS][0] = PEPA;
  attr_keys[CS_LAGR_YPLUS][1] = *jryplu;

  attr_keys[CS_LAGR_INTERF][0] = PEPA;
  attr_keys[CS_LAGR_INTERF][1] = *jrinpf;

  attr_keys[CS_LAGR_NEIGHBOR_FACE_ID][0] = IPEPA;
  attr_keys[CS_LAGR_NEIGHBOR_FACE_ID][1] = *jdfac;

  attr_keys[CS_LAGR_MARKO_VALUE][0] = IPEPA;
  attr_keys[CS_LAGR_MARKO_VALUE][1] = *jimark;

  attr_keys[CS_LAGR_N_LARGE_ASPERITIES][0] = IPEPA;
  attr_keys[CS_LAGR_N_LARGE_ASPERITIES][1] = *jnbasg;

  attr_keys[CS_LAGR_N_SMALL_ASPERITIES][0] = IPEPA;
  attr_keys[CS_LAGR_N_SMALL_ASPERITIES][1] = *jnbasp;

  attr_keys[CS_LAGR_ADHESION_FORCE][0] = PEPA;
  attr_keys[CS_LAGR_ADHESION_FORCE][1] = *jfadh;

  attr_keys[CS_LAGR_ADHESION_TORQUE][0] = PEPA;
  attr_keys[CS_LAGR_ADHESION_TORQUE][1] = *jmfadh;

  attr_keys[CS_LAGR_DISPLACEMENT_NORM][0] = PEPA;
  attr_keys[CS_LAGR_DISPLACEMENT_NORM][1] = *jndisp;

  attr_keys[CS_LAGR_TEMPERATURE][0] = EPTP_TS;
  if (*jtp > 0) {
    attr_keys[CS_LAGR_TEMPERATURE][1] = *jtp;
  }
  else {
    attr_keys[CS_LAGR_TEMPERATURE][1] = jhp[0];
    attr_keys[CS_LAGR_TEMPERATURE][2]
      = cs_glob_lagr_params->n_temperature_layers;
  }

  attr_keys[CS_LAGR_FLUID_TEMPERATURE][0] = EPTP_TS;
  attr_keys[CS_LAGR_FLUID_TEMPERATURE][1] = *jtf;

  attr_keys[CS_LAGR_WATER_MASS][1] = *jmwat;

  attr_keys[CS_LAGR_COAL_MASS][1] = jmch[0];
  attr_keys[CS_LAGR_COAL_MASS][2]
    = cs_glob_lagr_params->n_temperature_layers;

  attr_keys[CS_LAGR_COKE_MASS][1] = jmck[0];
  attr_keys[CS_LAGR_COKE_MASS][2]
    = cs_glob_lagr_params->n_temperature_layers;

  attr_keys[CS_LAGR_CP][1] = *jcp;

  /*  FIXME for post-processing by trajectory purpose. To have better
      results shrinking_diam should be transfered to EPTPA/EPTP arrays*/

  attr_keys[CS_LAGR_SHRINKING_DIAMETER][0] = PEPA;
  attr_keys[CS_LAGR_SHRINKING_DIAMETER][1] = *jrdck;

  attr_keys[CS_LAGR_INITIAL_DIAMETER][0] = PEPA;
  attr_keys[CS_LAGR_INITIAL_DIAMETER][1] = *jrd0p;

  attr_keys[CS_LAGR_COAL_DENSITY][0] = PEPA;
  attr_keys[CS_LAGR_COAL_DENSITY][1] = jrhock[0];
  attr_keys[CS_LAGR_COAL_DENSITY][2]
    = cs_glob_lagr_params->n_temperature_layers;

  attr_keys[CS_LAGR_COAL_NUM][0] = IPEPA;
  attr_keys[CS_LAGR_COAL_NUM][1] = *jinch;

  attr_keys[CS_LAGR_EMISSIVITY][1] = *jreps;

  if (*nbclst > 0) {
    attr_keys[CS_LAGR_STAT_CLASS][0] = IPEPA;
    attr_keys[CS_LAGR_STAT_CLASS][1] = *jclst;
  }

  if (*nvls > 0) {
    attr_keys[CS_LAGR_USER][1] = *jvls;
    attr_keys[CS_LAGR_USER][2] = *nvls;
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

  /* Initialize particle set : prev and current */

  _particle_set = _create_particle_set(0, 128, _p_attr_map);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER CREATION\n");
  cs_lagr_particle_set_dump(_particle_set);
#endif

  /* Initialization */

  for (i = 0; i < _particle_set->n_particles_max; i++) {

    cs_lagr_particles_set_lnum(_particle_set, i, CS_LAGR_SWITCH_ORDER_1,
                               CS_LAGR_SWITCH_OFF);

    _tracking_info(_particle_set, i)->state = CS_LAGR_PART_TO_SYNC;

  }

  /* Create all useful MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _cs_mpi_particle_type = _define_particle_datatype(_particle_set->p_am);
  }
#endif

  /* Initialize builder */

  _particle_track_builder
    = _init_track_builder(_particle_set->n_particles_max,
                          _particle_set->p_am->extents);

  /* Saving of itycel and icocel */
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  for (i = 0; i <= mesh->n_cells; i++)
    itycel[i] = builder->cell_face_idx[i];
  for (i = 0; i < builder->cell_face_idx[mesh->n_cells] ; i++)
    icocel[i] = builder->cell_face_lst[i];

  /* Now remap all pointers */
  /*------------------------*/

  const cs_lagr_attribute_map_t *p_am = _particle_set->p_am;

  *jisor = p_am->displ[0][CS_LAGR_CELL_NUM]/sizeof(cs_lnum_t) + 1;
  *jisora = p_am->displ[1][CS_LAGR_CELL_NUM]/sizeof(cs_lnum_t) + 1;

  *jirka = p_am->displ[1][CS_LAGR_RANK_ID]/sizeof(cs_lnum_t) + 1;

  *jord1 = p_am->displ[0][CS_LAGR_SWITCH_ORDER_1]/sizeof(cs_lnum_t) + 1;

  *jrval = p_am->displ[0][CS_LAGR_RANDOM_VALUE]/sizeof(cs_real_t) + 1;
  *jrpoi = p_am->displ[0][CS_LAGR_STAT_WEIGHT]/sizeof(cs_real_t) + 1;
  *jrtsp = p_am->displ[0][CS_LAGR_RESIDENCE_TIME]/sizeof(cs_real_t) + 1;
  *jdp = p_am->displ[0][CS_LAGR_DIAMETER]/sizeof(cs_real_t) + 1;
  *jmp = p_am->displ[0][CS_LAGR_MASS]/sizeof(cs_real_t) + 1;

  *jxp =  p_am->displ[0][CS_LAGR_COORDS]/sizeof(cs_real_t) + 1;
  *jyp = *jxp + 1;
  *jzp = *jyp + 1;

  *jup =  p_am->displ[0][CS_LAGR_VELOCITY]/sizeof(cs_real_t) + 1;
  *jvp = *jup + 1;
  *jwp = *jvp + 1;

  *juf =  p_am->displ[0][CS_LAGR_VELOCITY_SEEN]/sizeof(cs_real_t) + 1;
  *jvf = *juf + 1;
  *jwf = *jvf + 1;

  *jtaux = p_am->displ[0][CS_LAGR_TAUP_AUX]/sizeof(cs_real_t) + 1;

  if (p_am->count[0][CS_LAGR_TURB_STATE_1] > 0) {
    for (int id = 0; id < 3; id++)
      jbx1[id]
        = p_am->displ[0][CS_LAGR_TURB_STATE_1]/sizeof(cs_real_t) + 1 + id;
  }
  else {
    for (int id = 0; id < 3; id++)
      jbx1[id] = -1;
  }

  if (p_am->count[0][CS_LAGR_PRED_VELOCITY] > 0) {
    for (int id = 0; id < 3; id++)
      jtsup[id]
        = p_am->displ[0][CS_LAGR_PRED_VELOCITY]/sizeof(cs_real_t) + 1 + id;
  }
  else {
    for (int id = 0; id < 3; id++)
      jtsup[id] = -1;
  }

  if (p_am->count[0][CS_LAGR_PRED_VELOCITY_SEEN] > 0) {
    for (int id = 0; id < 3; id++)
      jtsuf[id]
        = p_am->displ[0][CS_LAGR_PRED_VELOCITY_SEEN]/sizeof(cs_real_t) + 1 + id;
  }
  else {
    for (int id = 0; id < 3; id++)
      jtsuf[id] = -1;
  }

  if (p_am->count[0][CS_LAGR_YPLUS] > 0)
    *jryplu = p_am->displ[0][CS_LAGR_YPLUS]/sizeof(cs_real_t) + 1;
  else
    *jryplu = -1;

  if (p_am->count[0][CS_LAGR_INTERF] > 0)
    *jrinpf = p_am->displ[0][CS_LAGR_INTERF]/sizeof(cs_real_t) + 1;
  else
    *jrinpf = -1;

  if (p_am->count[0][CS_LAGR_NEIGHBOR_FACE_ID] > 0)
    *jdfac = p_am->displ[0][CS_LAGR_NEIGHBOR_FACE_ID]/sizeof(cs_lnum_t) + 1;
  else
    *jdfac = -1;

  if (p_am->count[0][CS_LAGR_MARKO_VALUE] > 0)
    *jimark = p_am->displ[0][CS_LAGR_MARKO_VALUE]/sizeof(cs_lnum_t) + 1;
  else
    *jimark = -1;

  n_layers = p_am->count[0][CS_LAGR_TEMPERATURE];
  if (*jtp > 0) {
    assert(n_layers == 1);
    *jtp = p_am->displ[0][CS_LAGR_TEMPERATURE]/sizeof(cs_real_t) + 1;
  }
  else {
    *jtp = -1;
    for (int l_id = 0; l_id < n_layers; l_id++)
      jhp[l_id]
        = p_am->displ[0][CS_LAGR_TEMPERATURE]/sizeof(cs_real_t) + 1 + l_id;
  }

  if (p_am->count[0][CS_LAGR_FLUID_TEMPERATURE] > 0)
    *jtf = p_am->displ[0][CS_LAGR_FLUID_TEMPERATURE]/sizeof(cs_real_t) + 1;
  else
    *jtf = -1;

  if (p_am->count[0][CS_LAGR_WATER_MASS] > 0)
    *jmwat = p_am->displ[0][CS_LAGR_WATER_MASS]/sizeof(cs_real_t) + 1;
  else
    *jmwat = -1;

  n_layers = p_am->count[0][CS_LAGR_COAL_MASS];
  for (int l_id = 0; l_id < n_layers; l_id++)
    jmch[l_id] = p_am->displ[0][CS_LAGR_COAL_MASS]/sizeof(cs_real_t) + 1 + l_id;

  n_layers = p_am->count[0][CS_LAGR_COKE_MASS];
  for (int l_id = 0; l_id < n_layers; l_id++)
    jmck[l_id] = p_am->displ[0][CS_LAGR_COKE_MASS]/sizeof(cs_real_t) + 1 + l_id;

  if (p_am->count[0][CS_LAGR_CP] > 0)
    *jcp = p_am->displ[0][CS_LAGR_CP]/sizeof(cs_real_t) + 1;
  else
    *jcp = -1;

  if (p_am->count[0][CS_LAGR_SHRINKING_DIAMETER] > 0)
    *jrdck = p_am->displ[0][CS_LAGR_SHRINKING_DIAMETER]/sizeof(cs_real_t) + 1;
  else
    *jrdck = -1;

  if (p_am->count[0][CS_LAGR_INITIAL_DIAMETER] > 0)
    *jrd0p = p_am->displ[0][CS_LAGR_INITIAL_DIAMETER]/sizeof(cs_real_t) + 1;
  else
    *jrd0p = -1;

  if (p_am->count[0][CS_LAGR_COAL_NUM] > 0)
    *jinch = p_am->displ[0][CS_LAGR_COAL_NUM]/sizeof(cs_lnum_t) + 1;
  else
    *jinch = -1;

  n_layers = p_am->count[0][CS_LAGR_COAL_DENSITY];
  for (int l_id = 0; l_id < n_layers; l_id++) {
    jrhock[l_id]
      = p_am->displ[0][CS_LAGR_COAL_DENSITY]/sizeof(cs_real_t) + 1 + l_id;
  }

  if (p_am->count[0][CS_LAGR_EMISSIVITY] > 0)
    *jreps = p_am->displ[0][CS_LAGR_EMISSIVITY]/sizeof(cs_real_t) + 1;
  else
    *jreps = -1;

  if (p_am->count[0][CS_LAGR_DEPOSITION_FLAG] > 0)
    *jdepo = p_am->displ[0][CS_LAGR_DEPOSITION_FLAG]/sizeof(cs_lnum_t) + 1;
  else
    *jdepo = -1;

  if (p_am->count[0][CS_LAGR_N_LARGE_ASPERITIES] > 0)
    *jnbasg = p_am->displ[0][CS_LAGR_N_LARGE_ASPERITIES]/sizeof(cs_lnum_t) + 1;
  else
    *jnbasg = -1;

  if (p_am->count[0][CS_LAGR_N_SMALL_ASPERITIES] > 0)
    *jnbasp = p_am->displ[0][CS_LAGR_N_SMALL_ASPERITIES]/sizeof(cs_lnum_t) + 1;
  else
    *jnbasp = -1;

  if (p_am->count[0][CS_LAGR_ADHESION_FORCE] > 0)
    *jfadh = p_am->displ[0][CS_LAGR_ADHESION_FORCE]/sizeof(cs_real_t) + 1;
  else
    *jfadh = -1;

  if (p_am->count[0][CS_LAGR_ADHESION_TORQUE] > 0)
    *jmfadh = p_am->displ[0][CS_LAGR_ADHESION_TORQUE]/sizeof(cs_real_t) + 1;
  else
    *jmfadh = -1;

  if (p_am->count[0][CS_LAGR_DISPLACEMENT_NORM] > 0)
    *jndisp = p_am->displ[0][CS_LAGR_DISPLACEMENT_NORM]/sizeof(cs_real_t) + 1;
  else
    *jndisp = -1;

  if (p_am->count[0][CS_LAGR_STAT_CLASS] > 0)
    *jclst = p_am->displ[0][CS_LAGR_STAT_CLASS]/sizeof(cs_lnum_t) + 1;
  else
    *jclst = -1;

  int n_user = p_am->count[0][CS_LAGR_USER];
  for (int u_id = 0; u_id < n_user; u_id++) {
    jvls[u_id]
      = p_am->displ[0][CS_LAGR_USER]/sizeof(cs_real_t) + 1 + u_id;
  }
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
 *   scheme_order  <-> current order of the scheme used for Lagrangian
 *----------------------------------------------------------------------------*/

void
CS_PROCF (dplprt, DPLPRT)(cs_lnum_t        *p_scheme_order,
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
                          const cs_real_t  *tkelvi)
{
  cs_lnum_t  i, j;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  int        n_displacement_steps = 0;

  cs_lnum_t  n_delete_particles = 0;
  cs_lnum_t  n_failed_particles = 0;

  cs_real_t  failed_particle_weight = 0.0;
  cs_real_t  r_weight = 0.0;
  cs_real_t  tot_weight = 0.0;

  cs_lnum_t  scheme_order = *p_scheme_order;
  cs_lagr_particle_set_t  *particles = _particle_set;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  const cs_lagr_param_t *lagr_params = cs_glob_lagr_params;

  const cs_lnum_t  failsafe_mode = 0; /* If 1 : stop as soon as an error is
                                         detected */

  const cs_field_t *u = CS_F_(u);

  cs_real_t *part_b_mass_flux = NULL;

  if (*iflmbd) {
    assert(*iflm > 0);
    part_b_mass_flux = boundary_stat + ((*iflm -1) * mesh->n_b_faces);
  }

  assert(particles != NULL);

  particles->weight = 0.0;
  particles->n_part_out = 0;
  particles->n_part_dep = 0;
  particles->n_part_fou = 0;
  particles->weight_out = 0.0;
  particles->weight_dep = 0.0;
  particles->weight_fou = 0.0;
  particles->n_failed_part = 0;
  particles->weight_failed = 0.0;

  _initialize_displacement(particles);

  /* Shift from 1-based to 0-based (temporary) */

  if (p_am->count[0][CS_LAGR_NEIGHBOR_FACE_ID] > 0) {
    for (i = 0; i < particles->n_particles; i++) {
      const cs_lnum_t n_face_id
        = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_NEIGHBOR_FACE_ID) - 1;
        cs_lagr_particles_set_lnum(particles, i, CS_LAGR_NEIGHBOR_FACE_ID,
                                   n_face_id);
    }
  }

  /* Main loop on  particles: global propagation */

  while (_continue_displacement()) {

    n_delete_particles = 0;
    n_failed_particles = 0;

    r_weight = 0.0;
    tot_weight = 0.0;
    failed_particle_weight = 0.0;

    assert(particles->first_free_id != -1);

    /* Local propagation */

    for (i = 0, j = particles->first_used_id; i < particles->n_particles; i++) {

      unsigned char *particle = particles->p_buffer + p_am->extents * j;

      cs_lnum_t next_id = particles->used_id[j].next_id;

      /* Local copies of the current and previous particles state vectors
         to be used in case of the first pass of _local_propagation fails */

      cs_lnum_t cur_part_state = _get_tracking_info(particles, j)->state;

      if (cur_part_state == CS_LAGR_PART_TO_SYNC) {

        /* Main particle displacement stage */

        cur_part_state = _local_propagation(particle,
                                            p_am,
                                            n_displacement_steps,
                                            scheme_order,
                                            failsafe_mode,
                                            boundary_stat,
                                            &n_failed_particles,
                                            &failed_particle_weight,
                                            iensi3,
                                            inbr,
                                            inbrbd,
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
                                            part_b_mass_flux,
                                            u,
                                            energt,
                                            tprenc,
                                            visref,
                                            enc1,
                                            enc2,
                                            *tkelvi);

        _tracking_info(particles, j)->state = cur_part_state;

      }

      assert(next_id == particles->used_id[j].next_id);

      j = next_id;

    } /* End of loop on particles */

    /* Update of the particle set structure. Delete particles. */

    for (i = 0, j = particles->first_used_id; i < particles->n_particles; i++) {

      cs_lnum_t next_id = particles->used_id[j].next_id;

      cs_lnum_t cur_part_state = _get_tracking_info(particles, j)->state;

      cs_real_t cur_part_stat_weight
        = cs_lagr_particles_get_real(particles, j, CS_LAGR_STAT_WEIGHT);

      if (   cur_part_state == CS_LAGR_PART_TO_DELETE
          || cur_part_state == CS_LAGR_PART_OUT
          || cur_part_state == CS_LAGR_PART_ERR) {

        _remove_particle(particles, j);

        n_delete_particles++;

        r_weight += cur_part_stat_weight;

      }
      else {

        tot_weight += cur_part_stat_weight;

      }

      /* next_id was saved before potentially calling _remove_particle() */

      j = next_id;

    }

    particles->n_particles -= n_delete_particles;

    particles->weight = tot_weight;

    particles->n_part_out += n_delete_particles;
    particles->weight_out += r_weight;

    particles->n_failed_part += n_failed_particles;
    particles->weight_failed = failed_particle_weight;

    /*  assert(j == -1);  After a loop on particles, next_id of the last
        particle must not be defined */

    if (mesh->halo != NULL) {

      /* Synchronisation of a selection of particles for parallelism and
         periodicity. Number of particles on the local rank may change. */

      _lagr_halo_sync();
    }

    n_displacement_steps++;

  } /* End of while (global displacement) */

  /* Deposition sub-model additional loop */

  if (lagr_params->deposition > 0) {

    for (i = 0, j = particles->first_used_id; i < particles->n_particles; i++) {

      unsigned char *particle = particles->p_buffer + p_am->extents * j;

      cs_lnum_t next_id = particles->used_id[j].next_id;

      cs_lnum_t *cur_neighbor_face_id
        = cs_lagr_particle_attr(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
      cs_real_t *cur_part_yplus
        = cs_lagr_particle_attr(particle, p_am, CS_LAGR_YPLUS);

      _test_wall_cell(particle, p_am, visc_length, dlgeo,
                      cur_part_yplus, cur_neighbor_face_id);

      if (*cur_part_yplus < 100.e0) {

        /* TODO: specific treatment */

      }

      /* Particle resuspension specific treatment */

      if (lagr_params->resuspension > 0) {

        if (   cs_lagr_particles_get_lnum(particles, j, CS_LAGR_DEPOSITION_FLAG)
            == CS_LAGR_PART_ROLLING) {

          const cs_real_t *cur_part_coord
            = cs_lagr_particles_attr_n_const(particles, j, 0, CS_LAGR_COORDS);
          const cs_real_t *prev_part_coord
            = cs_lagr_particles_attr_n_const(particles, j, 1, CS_LAGR_COORDS);

          double traj[3] = {cur_part_coord[0] - prev_part_coord[0],
                            cur_part_coord[1] - prev_part_coord[1],
                            cur_part_coord[2] - prev_part_coord[2]};

          cs_real_t *cur_part_displacement_norm
            = cs_lagr_particles_attr(particles, j, CS_LAGR_DISPLACEMENT_NORM);

          *cur_part_displacement_norm += _get_norm(traj);

        }

      }

      j = next_id;

    }
  }

  _finalize_displacement(particles);
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
 * \brief Return pointer to the main cs_lagr_particle_set_t structure.
 *
 * \return
 *   pointer to current particle set, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_set_t  *
cs_lagr_get_particle_set(void)
{
  return _particle_set;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Delete particle set structure and other useful buffers.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_destroy(void)
{
  if (_particle_set == NULL)
    return;

  /* Destroy particle set */

  _particle_set = _destroy_particle_set(_particle_set);

  _destroy_attr_map(&_p_attr_map);

  /* Destroy builder */
  _particle_track_builder = _destroy_track_builder(_particle_track_builder);

  /* Destroy boundary condition structure */

  _lagr_bdy_conditions = _destroy_bdy_cond_struct(_lagr_bdy_conditions);

  /* Destroy the structure dedicated to clogging modeling */

  if (cs_glob_lagr_params->clogging)
    cs_lagr_clogging_finalize();

  /* Destroy the structure dedicated to roughness surface modeling */

  if (cs_glob_lagr_params->roughness)
    cs_lagr_roughness_finalize();

  /* Delete MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)  _delete_particle_datatypes();
#endif
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
cs_lagr_resize_particle_set(cs_lnum_t  n_min_particles)
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
    retval = _resize_particle_set(_particle_set, n_min_particles);

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

    if (particles->n_particles > 0 && particles->first_used_id > 0) {

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

      for (i = 0; i < particles->n_particles; i++) {
        bft_printf("  dump_particle_set i = %d \n", i);
        _dump_particle(particles, i);
      }

    }
  }
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

/* Delete local macro definitions */

END_C_DECLS
