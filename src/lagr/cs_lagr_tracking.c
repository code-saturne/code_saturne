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
 * Functions dealing with the particle tracking
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_printf.h>
#include <bft_error.h>
#include <bft_mem.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <cs_parall.h>
#include <cs_interface.h>
#include <fvm_periodicity.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_search.h"
#include "cs_lagr_utils.h"
#include "cs_halo.h"

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
#define  CS_LAGR_MIN_COMM_BUF_SIZE  10
#define  CS_LAGR_MAX_PROPAGATION_LOOPS  30
#define  N_VAR_PART_STRUCT  21
#define  N_VAR_PART_COAL     1
#define  N_VAR_PART_HEAT     1
#define  N_VAR_PART_AUX      1

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/* face_yplus auxiliary type */
/* ------------------------- */

typedef struct {

  cs_real_t  yplus;
  cs_lnum_t  face_id;

} face_yplus_t;

/* Main particle description */
/* ------------------------- */

typedef struct {

  cs_gnum_t   global_num;  /* global number of the particle */

  cs_lnum_t   cur_cell_num;
  cs_lnum_t   last_face_num;

#if defined(HAVE_MPI)
  int         cell_rank;
#endif

  int         switch_order_1;
  int         state;         /* < 0 : - number of the border face where
                                the particle is kept
                                0   : particle has to be destroyed
                                1   : particle has to be synchronized
                                2   : particle treated. End of displacement */
  int         visualized;    /* -1  : particle not visualized in displacement
                                      or trajectory mode
                                1   : particle visualized */

  cs_lnum_t   prev_id;  /* id in particle set of the previous particle */
  cs_lnum_t   next_id;  /* id in particle set of the next particle */

  cs_real_t   stat_weight;
  cs_real_t   residence_time;
  cs_real_t   mass;
  cs_real_t   diameter;
  cs_real_t   taup_aux;
  cs_real_t   coord[3];
  cs_real_t   velocity[3];
  cs_real_t   velocity_seen[3];

  /* Deposition submodel parameters */
  cs_real_t   yplus;
  cs_real_t   interf;
  cs_lnum_t   close_face_id;
  cs_lnum_t   marko_val;

} cs_lagr_particle_t;

/* Particle description for simulation with coal */
/* --------------------------------------------- */

typedef struct { /* Defined if IPHYLA == 2 */

  int         coal_number;

  cs_real_t   temp;
  cs_real_t   fluid_temp;
  cs_real_t   cp;
  cs_real_t   coal_mass;
  cs_real_t   coke_mass;
  cs_real_t   coke_diameter;
  cs_real_t   coke_fraction;

  /* 2 others left */

} cs_lagr_coal_particle_t;


/* Particle description for simulation with heat transfert */
/* ------------------------------------------------------- */

typedef struct { /* Defined if IPHYLA == 1 */

  cs_real_t   temp;
  cs_real_t   fluid_temp;
  cs_real_t   cp;
  cs_real_t   emissivity; /* Only useful if there is an equation on temperature
                             and if radiative transfert is active */

} cs_lagr_heat_particle_t;


typedef struct { /* User-defined variables. Max. 10 */

  cs_lnum_t   stat_class;  /* Only if NBCLST > 0 */
  cs_real_t   aux[10];

} cs_lagr_aux_particle_t;


typedef struct {

  cs_lnum_t  n_particles;
  cs_lnum_t  n_part_out;
  cs_lnum_t  n_part_dep;
  cs_lnum_t  n_failed_part;

  cs_real_t  weight;
  cs_real_t  weight_out;
  cs_real_t  weight_dep;
  cs_real_t  weight_failed;

  cs_lnum_t  n_particles_max;

  cs_lnum_t  first_used_id;
  cs_lnum_t  first_free_id;

  cs_lagr_particle_t       *particles;  /* Main  particle description */

  cs_lagr_coal_particle_t  *coal_desc;  /* Additional description for study
                                           with coal */
  cs_lagr_heat_particle_t  *heat_desc;  /* Additional description for study
                                           with heat transfert */
  cs_lagr_aux_particle_t   *aux_desc;   /* Additional description for study
                                           with user-defined variables */

} cs_lagr_particle_set_t;


/* Structures useful to manage the exchange of particles
   between communicating ranks */

typedef struct {

  cs_lnum_t  n_cells;        /* Number of cells in the halo */
  cs_lnum_t *rank;           /* value between [0, n_c_domains-1]
                                (cf. cs_halo.h) */
  cs_lnum_t *dist_cell_num;  /* local cell num. on distant ranks */
  cs_lnum_t *transform_id;   /* In case of periodicity, transformation
                                associated to a given halo cell */

  /* Buffer used to exchange particle between communicating ranks */

  cs_lnum_t *send_count;     /* To store the number of particles to send to
                                each communicating rank */
  cs_lnum_t *recv_count;     /* To store the number of particles to receive from
                                each communicating rank */

  cs_lnum_t *send_shift;
  cs_lnum_t *recv_shift;

  cs_lagr_particle_set_t  *send_buf;
  cs_lagr_particle_set_t  *recv_buf;

#if defined(HAVE_MPI)
  MPI_Request   *request;
  MPI_Status    *status;
#endif

} cs_lagr_halo_t;

#if 0
/* Structures useful to manage the dynamic flow of particles coming in/out */

typedef struct _cs_lagr_list_t cs_lagr_list_t;
struct _cs_lagr_list_t {

  cs_lnum_t        val;  /* id of the free space in particle set */
  cs_lagr_list_t  *next; /* pointer to the next item */

};

typedef struct {

  cs_lagr_list_t   *free_spaces;  /* List of free spaces in particle set */
  cs_lnum_t         size;         /* Current size of the list
                                     Max size available in particle set */

} cs_lagr_stack_t;
#endif

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
  int  n_stat_classes;
  int   n_user_variables;

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
  CS_LAGR_IDEPFA = 13
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

enum {
  CS_LAGR_PART_MOVE_OFF = 0,
  CS_LAGR_PART_MOVE_ON  = 1
};

/* According to the scheme order is degenerated to order 1 */

enum {
  CS_LAGR_SWITCH_OFF = 0,
  CS_LAGR_SWITCH_ON = 1
};

/* Global variable for the current subroutines */

static  cs_lagr_particle_set_t  *_particle_set = NULL;
static  cs_lagr_particle_set_t  *_prev_particle_set = NULL;
static  cs_lagr_track_builder_t  *_particle_track_builder = NULL;
static  cs_lagr_bdy_condition_t  *_lagr_bdy_conditions = NULL;

static  cs_lagr_param_t  cs_glob_lagr_param; // Should move to cs_lagr.c

enum {X, Y, Z};  /* Used for _get_norm() and _get_dot_prod() */

/* MPI datatype associated to each particle structures */

#if defined(HAVE_MPI)
static  MPI_Datatype  _CS_MPI_PARTICLE;
static  MPI_Datatype  _CS_MPI_COAL_PARTICLE;
static  MPI_Datatype  _CS_MPI_HEAT_PARTICLE;
static  MPI_Datatype  _CS_MPI_AUX_PARTICLE;
#endif

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the norm of a 3D vector of double (cs_real_t)
 *
 * parameters:
 *  vector  --> 3D vector to treat
 *
 * returns:
 *  norm associated to the vector
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_get_norm(cs_real_t  vect[])
{
  return sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z]);
}

/*----------------------------------------------------------------------------
 * Compute the norm of a 3D vector of double (cs_real_t)
 *
 * parameters:
 *  v1  --> 3D vector to treat
 *  v2  --> 3D vector to treat
 *
 * returns:
 *  norm associated to the vector
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_get_dot_prod(cs_real_t   v1[],
              cs_real_t   v2[])
{
  return (v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z]);
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a given
 * vector of coordinates.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   xyz_in       --> input vector
 *   xyz_out      --> output vector
 *----------------------------------------------------------------------------*/

inline static void
_apply_vector_transfo(cs_real_t    matrix[3][4],
                      cs_real_t    xyz_in[],
                      cs_real_t    xyz_out[])
{
  cs_lnum_t  i, j;
  cs_real_t  xyz_a[3 + 1], xyz_b[3];

  /* Define a vector in homogeneous coordinates before transformation */

  for (j = 0; j < 3; j++)
    xyz_a[j] = xyz_in[j];
  xyz_a[3] = 1;

  /* Initialize output */

  for (i = 0; i < 3; i++)
    xyz_b[i] = 0.;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      xyz_b[i] += matrix[i][j]*xyz_a[j];

  /* Store updated cell center */

  for (j = 0; j < 3; j++)
    xyz_out[j] = xyz_b[j];

}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a rotation to a given vector.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   x_in         --> X coord. of the incoming vector
 *   y_in         --> Y coord. of the incoming vector
 *   z_in         --> Z coord. of the incoming vector
 *   x_out        <-- pointer to the X coord. of the output
 *   y_out        <-- pointer to the Y coord. of the output
 *   z_out        <-- pointer to the Z coord. of the output
 *----------------------------------------------------------------------------*/

inline static void
_apply_vector_rotation(cs_real_t   matrix[3][4],
                       cs_real_t   x_in,
                       cs_real_t   y_in,
                       cs_real_t   z_in,
                       cs_real_t   *x_out,
                       cs_real_t   *y_out,
                       cs_real_t   *z_out)
{
  *x_out = matrix[0][0] * x_in + matrix[0][1] * y_in + matrix[0][2] * z_in;
  *y_out = matrix[1][0] * x_in + matrix[1][1] * y_in + matrix[1][2] * z_in;
  *z_out = matrix[2][0] * x_in + matrix[2][1] * y_in + matrix[2][2] * z_in;
}

/*----------------------------------------------------------------------------
 * Remove a particle from a particle set.
 *
 * parameters:
 *
 *
 * returns:
 *  NULL
 *----------------------------------------------------------------------------*/

static void
_remove_particle(cs_lagr_particle_set_t   *set,
                 cs_lagr_particle_t        cur_part,
                 cs_lnum_t                  cur_id)

{

  /* Remove cur_part from "used particle list" */
  if (cur_part.prev_id != -1)
    set->particles[cur_part.prev_id].next_id = cur_part.next_id;
  else
    set->first_used_id = cur_part.next_id;

  if (cur_part.next_id != set->n_particles_max)
  {
    if ( cur_part.next_id != -1)
      set->particles[cur_part.next_id].prev_id = cur_part.prev_id;
  }
  /* Add cur_part to "free particle list" */

  if (cur_id < set->first_free_id) {

    cs_lnum_t old_first_free = set->first_free_id;

    set->first_free_id = cur_id;

    set->particles[set->first_free_id].next_id = old_first_free;

    set->particles[set->first_free_id].prev_id
      = set->particles[old_first_free].prev_id;

    set->particles[old_first_free].prev_id = cur_id;

  }
  else { /* We place the cur_part just behind the first free particle " */

    cs_lnum_t first_free = set->first_free_id;

    cs_lnum_t old_next = set->particles[first_free].next_id;

    set->particles[first_free].next_id = cur_id;
    set->particles[cur_id].next_id = set->particles[first_free].next_id;
    set->particles[cur_id].prev_id = first_free;

    set->particles[set->particles[first_free].next_id].next_id = old_next;
  }
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------
 * Create a MPI_Datatype which maps main particle characteristics.
 *
 * returns:
 *  a MPI_Datatype
 *----------------------------------------------------------------------------*/

static MPI_Datatype
_define_particle_datatype(void)
{
  int  i, j;
  MPI_Datatype  new_type;
  cs_lagr_particle_t  part;

  int  count = 0;
  int  blocklengths[N_VAR_PART_STRUCT]
    = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1, 1, 1, 1};
  MPI_Aint  displacements[N_VAR_PART_STRUCT]
    = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Datatype  types[N_VAR_PART_STRUCT] = {CS_MPI_GNUM,
                                            CS_MPI_LNUM,
                                            CS_MPI_LNUM,
                                            MPI_INT,
                                            MPI_INT,
                                            MPI_INT,
                                            MPI_INT,
                                            CS_MPI_LNUM,
                                            CS_MPI_LNUM,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_REAL,
                                            CS_MPI_LNUM,
                                            CS_MPI_LNUM,
  };

  /* Initialize a default particle */

  part.global_num = 1;
  part.cur_cell_num = 1;
  part.last_face_num = 1;
  part.cell_rank = 0;
  part.visualized = -1;
  part.switch_order_1 = CS_LAGR_SWITCH_OFF;
  part.state = CS_LAGR_PART_TO_SYNC;
  part.prev_id = 0;
  part.next_id = 0;
  part.stat_weight = 0.0;
  part.residence_time = 0.0;
  part.mass = 0.0;
  part.diameter = 0.0;
  part.taup_aux = 0.0;

  for (j = 0; j < 3; j++) {
    part.coord[j] = 0.0;
    part.velocity[j] = 0.0;
    part.velocity_seen[j] = 0.0;
  }

  part.yplus = 0.0;
  part.interf = 0.0;
  part.close_face_id = 0;
  part.marko_val = 0;

  /* Define array of displacements */

  MPI_Get_address(&part, displacements + count++);
  MPI_Get_address(&part.cur_cell_num, displacements + count++);
  MPI_Get_address(&part.last_face_num, displacements + count++);
  MPI_Get_address(&part.cell_rank, displacements + count++);
  MPI_Get_address(&part.switch_order_1, displacements + count++);
  MPI_Get_address(&part.state, displacements + count++);
  MPI_Get_address(&part.visualized, displacements + count++);
  MPI_Get_address(&part.next_id, displacements + count++);
  MPI_Get_address(&part.prev_id, displacements + count++);
  MPI_Get_address(&part.stat_weight, displacements + count++);
  MPI_Get_address(&part.residence_time, displacements + count++);
  MPI_Get_address(&part.mass, displacements + count++);
  MPI_Get_address(&part.diameter, displacements + count++);
  MPI_Get_address(&part.taup_aux, displacements + count++);
  MPI_Get_address(&part.coord, displacements + count++);
  MPI_Get_address(&part.velocity, displacements + count++);
  MPI_Get_address(&part.velocity_seen, displacements + count++);
  MPI_Get_address(&part.yplus, displacements + count++);
  MPI_Get_address(&part.interf, displacements + count++);
  MPI_Get_address(&part.close_face_id, displacements + count++);
  MPI_Get_address(&part.marko_val, displacements + count++);

  assert(count == N_VAR_PART_STRUCT);

  for (i = N_VAR_PART_STRUCT - 1; i >= 0; i--)
    displacements[i] -= displacements[0];

  assert(fabs(displacements[0]) < 1e-15);


  /* Create new datatype */

  MPI_Type_create_struct(N_VAR_PART_STRUCT,
                         blocklengths, displacements, types, &new_type);

  MPI_Type_commit(&new_type);

  return new_type;
}

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype which maps coal particle characteristics.
 *
 * parameters:
 *
 * returns:
 *  a MPI_Datatype
 *----------------------------------------------------------------------------*/
static MPI_Datatype
_define_coal_particle_datatype(void)
{
  MPI_Datatype  new_type;

  int  blocklengths[N_VAR_PART_COAL] = {1};
  MPI_Aint  displacements[N_VAR_PART_COAL] = {0};
  MPI_Datatype  types[N_VAR_PART_COAL] = {CS_MPI_GNUM};

  MPI_Type_create_struct(N_VAR_PART_COAL,
                         blocklengths, displacements, types, &new_type);

  MPI_Type_commit(&new_type);

  return new_type;
}

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype which maps heat particle characteristics.
 *
 * parameters:
 *
 * returns:
 *  a MPI_Datatype
 *----------------------------------------------------------------------------*/
static MPI_Datatype
_define_heat_particle_datatype(void)
{
  MPI_Datatype  new_type;

  int  blocklengths[N_VAR_PART_HEAT] = {1};
  MPI_Aint  displacements[N_VAR_PART_HEAT] = {0};
  MPI_Datatype  types[N_VAR_PART_HEAT] = {CS_MPI_GNUM};

  MPI_Type_create_struct(N_VAR_PART_HEAT,
                         blocklengths, displacements, types, &new_type);

  MPI_Type_commit(&new_type);

  return new_type;
}

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype which maps aux. particle characteristics.
 *
 * returns:
 *  a MPI_Datatype
 *----------------------------------------------------------------------------*/
static MPI_Datatype
_define_aux_particle_datatype(void)
{
  MPI_Datatype  new_type;

  int  blocklengths[N_VAR_PART_AUX] = {1};
  MPI_Aint  displacements[N_VAR_PART_AUX] = {0};
  MPI_Datatype  types[N_VAR_PART_AUX] = {CS_MPI_GNUM};

  MPI_Type_create_struct(N_VAR_PART_AUX,
                         blocklengths, displacements, types, &new_type);

  MPI_Type_commit(&new_type);

  return new_type;
}

/*----------------------------------------------------------------------------
 * Delete all the MPI_Datatypes related to particles.
 *----------------------------------------------------------------------------*/

static void
_delete_particle_datatypes(void)
{
  MPI_Type_free(&_CS_MPI_PARTICLE);
  MPI_Type_free(&_CS_MPI_COAL_PARTICLE);
  MPI_Type_free(&_CS_MPI_HEAT_PARTICLE);
  MPI_Type_free(&_CS_MPI_AUX_PARTICLE);
}
#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Allocate a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *  n_particles_max   -->  local max. number of particles
 *
 * returns:
 *  a new allocated cs_lagr_particle_set_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_particle_set_t *
_create_particle_set(const cs_lnum_t n_particles_max)
{
  cs_lnum_t  i;

  cs_lagr_particle_set_t  *new_set = NULL;

  if (n_particles_max == 0)
    return NULL;

  BFT_MALLOC(new_set, 1, cs_lagr_particle_set_t);
  BFT_MALLOC(new_set->particles, n_particles_max, cs_lagr_particle_t);

  new_set->n_particles_max = n_particles_max;
  new_set->n_particles = 0;
  new_set->first_used_id = -1;
  new_set->first_free_id = 0;

  assert(n_particles_max >= 1);

  new_set->particles[0].prev_id = -1;
  new_set->particles[0].next_id = 1;

  for (i = 1; i < n_particles_max; i++) {
    new_set->particles[i].prev_id = i-1;
    new_set->particles[i].next_id = i+1;
  }

  new_set->coal_desc = NULL;
  new_set->heat_desc = NULL;
  new_set->aux_desc = NULL;

  if (cs_glob_lagr_param.physic_mode == 1)
    BFT_MALLOC(new_set->heat_desc, n_particles_max, cs_lagr_heat_particle_t);

  else if (cs_glob_lagr_param.physic_mode == 2)
    BFT_MALLOC(new_set->coal_desc, n_particles_max, cs_lagr_coal_particle_t);

  if (   cs_glob_lagr_param.n_user_variables > 0
      || cs_glob_lagr_param.n_stat_classes > 0)
    BFT_MALLOC(new_set->aux_desc, n_particles_max, cs_lagr_aux_particle_t);

  return new_set;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *  particle_set    <->  a cs_lagr_particle_set_t structure
 *
 * returns:
 *  NULL
 *----------------------------------------------------------------------------*/

static cs_lagr_particle_set_t *
_destroy_particle_set(cs_lagr_particle_set_t *set)
{
  if (set == NULL)
    return set;

  BFT_FREE(set->particles);

  if (set->coal_desc != NULL)
    BFT_FREE(set->coal_desc);

  if (set->heat_desc != NULL)
    BFT_FREE(set->heat_desc);

  if (set->aux_desc != NULL)
    BFT_FREE(set->aux_desc);

  BFT_FREE(set);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Dump a cs_lagr_particle_t structure
 *
 * parameter
 *  part   --> cs_lagr_particle_t structure to dump
 *----------------------------------------------------------------------------*/

static void
_dump_particle(cs_lagr_particle_t  part)
{
  bft_printf("  particle:\n"
             "\tglobal num    : %llu\n"
             "\tcur_cell_num  : %d\n"
             "\tstate         : %d\n"
             "\tprev_id       : %d\n"
             "\tnext_id       : %d\n"
             "\tcoord         : [%e, %e, %e]\n",
             (unsigned long long)part.global_num,
             (int)part.cur_cell_num,
             (int)part.state,
             (int)part.prev_id,
             (int)part.next_id,
             (double)part.coord[0],
             (double)part.coord[1],
             (double)part.coord[2]);
#if defined(HAVE_MPI)
  bft_printf("\tcell_rk       : %d\n", part.cell_rank);
#endif
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Dump a cs_lagr_particle_t structure
 *
 * parameter
 *  part   --> cs_lagr_particle_t structure to dump
 *----------------------------------------------------------------------------*/

static void
_dump_particle_set(cs_lagr_particle_set_t   *set)
{
  if (set != NULL) {

    cs_lnum_t  i, j;

    bft_printf("    n_particles      :  %9d\n", set->n_particles);
    bft_printf("    n_particles_max  :  %9d\n", set->n_particles_max);
    bft_printf("    first_used_id:  %9d\n", set->first_used_id);
    bft_printf("    first_free_id:  %9d\n", set->first_free_id);
    bft_printf_flush();

    if (set->n_particles > 0) {

      for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {
        bft_printf("dump_particle_set i j = %d %d \n",i,j);
        _dump_particle(set->particles[j]);
        j = set->particles[j].next_id;
      }

      assert(j == -1); /* The next_id is not defined for the last particle
                          of the particle set */

    }

  }
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Resize a cs_lagr_particle_set_t structure.
 *
 * parameters:
 *  p_particle_set    <->  pointer to a cs_lagr_particle_set_t structure
 *  n_particles_max   -->  local max. number of particles
 *----------------------------------------------------------------------------*/

static void
_resize_particle_set(cs_lagr_particle_set_t  **p_particle_set,
                     const cs_lnum_t           n_particles_max)
{
  cs_lagr_particle_set_t  *particle_set = *p_particle_set;

  assert(n_particles_max >= 0);

  if (n_particles_max == 0)
    particle_set = _destroy_particle_set(particle_set);

  else if (particle_set == NULL && n_particles_max > 0)
    particle_set = _create_particle_set(n_particles_max);

  else if (particle_set->n_particles_max < n_particles_max) {

    particle_set->n_particles_max = n_particles_max;

    BFT_REALLOC(particle_set->particles, n_particles_max, cs_lagr_particle_t);

    particle_set->coal_desc = NULL;
    particle_set->heat_desc = NULL;
    particle_set->aux_desc = NULL;

    if (cs_glob_lagr_param.physic_mode == 1)
      BFT_REALLOC(particle_set->heat_desc,
                  n_particles_max,
                  cs_lagr_heat_particle_t);

    else if (cs_glob_lagr_param.physic_mode == 2)
      BFT_REALLOC(particle_set->coal_desc,
                  n_particles_max,
                  cs_lagr_coal_particle_t);

    if (   cs_glob_lagr_param.n_user_variables > 0
        || cs_glob_lagr_param.n_stat_classes > 0)
      BFT_REALLOC(particle_set->aux_desc,
                  n_particles_max,
                  cs_lagr_aux_particle_t);

  }
  else {

    /* FIX ME */
    /* bft_printf("nmax local = %d\n",particle_set->n_particles_max); */
    /* bft_printf("n demande : %d\n",n_particles_max); */

    /* bft_error(__FILE__, __LINE__, 0, */
    /*           _(" The current situation is not managed.\n")); */

  }

  /* Returns pointer */

  *p_particle_set = particle_set;
}

/*----------------------------------------------------------------------------
 * Define a cs_lagr_halo_t structure to deal with parallelism and
 * periodicity
 *
 * parameters:
 *  n_particles_max   -->  local max number of particles
 *
 * returns:
 *  a new allocated cs_lagr_halo_t structure.
 *----------------------------------------------------------------------------*/

static cs_lagr_halo_t *
_create_lagr_halo(cs_lnum_t  n_particles_max)
{
  cs_lnum_t  i, rank, tr_id, shift, start, end, n;

  cs_lnum_t  buf_size = CS_LAGR_MIN_COMM_BUF_SIZE;
  cs_lnum_t  halo_cell_id = 0;
  cs_lnum_t  *cell_num = NULL;
  cs_lagr_halo_t  *lagr_halo = NULL;

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

  /* FIXME: check the rule for the size of buf_size */
  buf_size = CS_MAX(buf_size, n_particles_max / CS_LAGR_MIN_COMM_BUF_SIZE);

  lagr_halo->send_buf = _create_particle_set(buf_size);
  lagr_halo->recv_buf = _create_particle_set(buf_size);

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

  if (mesh->n_init_perio > 0) { /* Periodicity is activate */

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

    } /* End of loop on transformation */

  } /* End if periodicity is activate */

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
 *  halo    --> cs_lagr_halo_t structure to delete
 *
 * returns:
 *  NULL
 *----------------------------------------------------------------------------*/

static cs_lagr_halo_t *
_delete_lagr_halo(cs_lagr_halo_t   *halo)
{
  if (halo == NULL)
    return NULL;

  BFT_FREE(halo->rank);
  BFT_FREE(halo->transform_id);
  BFT_FREE(halo->dist_cell_num);

  BFT_FREE(halo->send_shift);
  BFT_FREE(halo->send_count);
  BFT_FREE(halo->recv_shift);
  BFT_FREE(halo->recv_count);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    BFT_FREE(halo->request);
    BFT_FREE(halo->status);
  }
#endif

  halo->send_buf = _destroy_particle_set(halo->send_buf);
  halo->recv_buf = _destroy_particle_set(halo->recv_buf);

  BFT_FREE(halo);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Define a cell -> face connectivity. Index begins with 0.
 *
 * parameters:
 *  builder   -->  pointer to a cs_lagr_track_builder_t structure
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
 *  n_particles_max   -->  local max number of particles
 *
 * returns:
 *  a new defined cs_lagr_track_builder_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_track_builder_t *
_init_track_builder(cs_lnum_t  n_particles_max)
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
    builder->halo = _create_lagr_halo(n_particles_max);
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
 *  builder   -->  pointer to a cs_lagr_track_builder_t structure
 *
 * returns:
 *  NULL
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

  builder->halo = _delete_lagr_halo(builder->halo);

  /* Destrou the builder structure */

  BFT_FREE(builder);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Initialize a cs_lagr_bdy_condition_t structure.
 *
 * parameters:
 *  n_max_zones     -->  number max. of boundary zones
 *
 * returns:
 *  a new defined cs_lagr_bdy_condition_t structure
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
 *  n_max_zones     -->  number max. of boundary zones
 *
 * returns:
 *  a new defined cs_lagr_bdy_condition_t structure
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
 *  a NULL pointer.
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
 *  failsafe_mode           -->
 *  particle                -->
 *  p_n_failed_particles     <->
 *  p_failed_particle_weight <->
 *
 *----------------------------------------------------------------------------*/

static void
_manage_error(cs_lnum_t              failsafe_mode,
              cs_lagr_particle_t     particle,
              cs_lnum_t             *p_error,
              cs_lnum_t             *p_n_failed_particles,
              cs_real_t             *p_failed_particle_weight,
              const char            *msg)

{
  cs_lnum_t  error = *p_error;
  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight =*p_failed_particle_weight;

  particle.cur_cell_num = 0;
  n_failed_particles++;
  failed_particle_weight += particle.stat_weight;

  error = 0;

  if (failsafe_mode == 1)
    bft_error(__FILE__, __LINE__, 0, _("%s\n"), msg);

  /* Return pointers */

  *p_error = error;
  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;
}

/*----------------------------------------------------------------------------
 * Test if all displacements are finished for all ranks.
 *
 * parameters:
 *  n_particles     -->  local number of particles
 *
 * returns:
 *  true if there is a need to move particles or false, otherwise
 *----------------------------------------------------------------------------*/

static bool
_continue_displacement(void)
{
  cs_lnum_t  i, j;
  cs_lnum_t  _test = 1, test = 1;

  const cs_lagr_particle_set_t  *set = _particle_set;
  const cs_lnum_t  n_particles = set->n_particles;

  for (i = 0, j = set->first_used_id; i < n_particles; i++) {
    if (set->particles[j].state == CS_LAGR_PART_TO_SYNC) {
      _test = 0;
      break;
    }
    j = set->particles[j].next_id;
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
 *  face_num        --> local number of the studied face
 *  n_vertices      --> size of the face connectivity
 *  face_connect    --> face -> vertex connectivity
 *  prev_particle   --> data relative to the particle for the previous time
 *                      step
 *  particle        --> data relative to the particle for the current time
 *                      step
 *  p_error         <-> pointer to an error indicator
 *
 * returns:
 *  -1: if the particle keeps inside the same cell
 *   0: if the trajectory doesn't go through the current face
 *   1: if the trajectory goes through the current face
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_where_are_you(cs_lnum_t            face_num,
               cs_lnum_t            n_vertices,
               cs_lnum_t           *face_connect,
               cs_lagr_particle_t   prev_particle,
               cs_lagr_particle_t   particle,
               cs_lnum_t           *p_error)
{
  cs_lnum_t  i, j, vtx_id1, vtx_id2;
  cs_real_t  face_cog[3], cell_cen[3], vtx1[3], vtx2[3];

  cs_lnum_t  cur_cell_id = particle.cur_cell_num - 1;
  cs_real_t  prev_location[3] = { prev_particle.coord[0],
                                  prev_particle.coord[1],
                                  prev_particle.coord[2]};
  cs_real_t  next_location[3] = { particle.coord[0],
                                  particle.coord[1],
                                  particle.coord[2]};

  cs_mesh_t  *mesh = cs_glob_mesh;

  /* Initialize local parameters */

  cs_real_t  max_value = DBL_MIN;
  cs_lnum_t  orient_count = 0;
  cs_lnum_t  face_orient = 0;
  cs_lnum_t  orient = 0;
  cs_lnum_t  orient_test = 0;
  cs_lnum_t  first_orient = 0;
  cs_lnum_t  ijkl_ref = 0;
  cs_lnum_t  colocalization = -1;
  cs_lnum_t  error = 0;
  cs_lnum_t  indian = -999; /* initialize to an incoherent value */
  cs_lnum_t  perturbation = 0; /* No perturbation algorithm used */

  assert(sizeof(cs_real_t) == 8*sizeof(char));

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
    indian = -1;
    return indian;
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
                                             vtx2,
                                             perturbation);

  /* Special treatment in case of periodicity  */

  if (mesh->n_init_perio > 0)
    face_orient = 0;

  if (face_orient == 0)  /* => coplanar. Change prev_location
                            by cell center */
    face_orient = cs_lagrang_tetra_orientation(cell_cen,
                                               face_cog,
                                               vtx1,
                                               vtx2,
                                               perturbation);

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

    return indian;
  }

  /* Test first vertex of the face: [P, Q, Face_CoG, V1] */

  first_orient = cs_lagrang_tetra_orientation(prev_location,
                                              next_location,
                                              face_cog,
                                              vtx1,
                                              perturbation);

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

    return indian;
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
                                          vtx1,
                                          perturbation);

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

      return indian;
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
                                                 vtx1,
                                                 perturbation);

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

        return indian;
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

    return indian;
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

    return indian;
  }

  if (orient_count == 0 || orient_count == -2) {
    indian = 0;
    return indian;
  }

  /* Relative position between the current face and the starting and ending
     particle locations */

  assert(orient_count == 2);

  perturbation = 1;

  vtx_id1 = face_connect[ijkl_ref - 1] - 1;
  vtx_id2 = face_connect[ijkl_ref] - 1;

  for (j = 0; j < 3; j++) {
    vtx1[j] = mesh->vtx_coord[3*vtx_id1+j];
    vtx2[j] = mesh->vtx_coord[3*vtx_id2+j];
  }

  orient = cs_lagrang_tetra_orientation(next_location,
                                        face_cog,
                                        vtx1,
                                        vtx2,
                                        perturbation);

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

    return indian;
  }

  indian = -orient * face_orient;

  /* Returns pointers */

  *p_error = error;

  return indian;
}

/*----------------------------------------------------------------------------
 * Deposition model
 * ----------------
 * Calculate the number of the closest face to the particle
 * as wall as the corresponding wall normal distance (y_p^+)
 *----------------------------------------------------------------------------*/

static void
_test_wall_cell(cs_lagr_particle_t *particle,
                cs_real_t           visc_length[],
                cs_real_t           dlgeo[])
{

  cs_lnum_t cell_num = particle->cur_cell_num;

  if (cell_num < 0) return;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;
  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;
  cs_lnum_t cell_id = cell_num - 1;

  particle->yplus = 10000;
  particle->close_face_id = -1;

  cs_lnum_t  start = cell_face_idx[cell_id];
  cs_lnum_t  end =  cell_face_idx[cell_id + 1];

  cs_lnum_t i;

  for (i = start; i < end; i++)
  {
    cs_lnum_t  face_num = cell_face_lst[i];

    if (face_num < 0)
    {
      cs_lnum_t face_id = CS_ABS(face_num) - 1;

      cs_lnum_t boundary_zone = bdy_conditions->b_face_zone_num[face_id]-1;

      if ( (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) ||
           (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO2) ||
           (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IREBOL) ||
           (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPFA) )
      {

        cs_real_t x_faceid = dlgeo[face_id];
        cs_real_t y_faceid = dlgeo[face_id + (mesh->n_b_faces)];
        cs_real_t z_faceid = dlgeo[face_id + (mesh->n_b_faces) * 2];
        cs_real_t offset_faceid = dlgeo[face_id + (mesh->n_b_faces) * 3];


        cs_real_t dist_norm =   CS_ABS( particle->coord[0] * x_faceid +
                                      particle->coord[1] * y_faceid +
                                      particle->coord[2] * z_faceid +
                                        offset_faceid ) / visc_length[face_id];


        if (dist_norm < particle->yplus)
        {
          particle->yplus = dist_norm;
          particle->close_face_id = face_id;
        }
      }
    }

  }

}


/*----------------------------------------------------------------------------
 * Test if the current particle moves to the next cell through  this face
 *
 * parameters:
 *
 *
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_bdy_treatment(cs_lagr_particle_t   *p_prev_particle,
               cs_lagr_particle_t   *p_particle,
               cs_lnum_t             face_num,
               cs_real_t            *boundary_stat,
               cs_lnum_t             boundary_zone,
               cs_lnum_t             failsafe_mode,
               cs_lnum_t            *p_move_particle,
               cs_lnum_t            *p_n_failed_particles,
               cs_real_t            *p_failed_particle_weight,
               const cs_lnum_t      *iensi3,
               const cs_lnum_t      *nvisbr,
               const cs_lnum_t      *inbr,
               const cs_lnum_t      *inbrbd,
               const cs_lnum_t      *iflm,
               const cs_lnum_t      *iflmbd,
               const cs_lnum_t      *iang,
               const cs_lnum_t      *iangbd,
               const cs_lnum_t      *ivit,
               const cs_lnum_t      *ivitbd,
               const cs_lnum_t      *nusbor,
               cs_lnum_t             iusb[])

{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  int nfabor  = mesh->n_b_faces;

  cs_lnum_t  k;
  cs_real_t  tmp;
  cs_real_t  depl[3], face_normal[3], face_cog[3], intersect_pt[3];

  cs_real_t compo_vit[3] = {0.0, 0.0, 0.0};
  cs_real_t norm_vit = 0.0;

  cs_real_t  abs_curv = 0.0;
  cs_lagr_particle_t  prev_particle = *p_prev_particle;
  cs_lagr_particle_t  particle = *p_particle;
  cs_lnum_t  move_particle = *p_move_particle;
  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight = *p_failed_particle_weight;

  cs_lnum_t  face_id = face_num - 1;
  cs_lnum_t  error = 0; /* FIXME: Not very useful -> _manage_error() */
  cs_lnum_t  particle_state = -999;


  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;

  assert(builder != NULL);
  assert(bdy_conditions != NULL);

  for (k = 0; k < 3; k++)
    depl[k] = particle.coord[k] - prev_particle.coord[k];

  if (fabs(depl[0]) < 1e-15 && fabs(depl[1]) < 1e-15 && fabs(depl[2]) < 1e-15)
    return 0; /* move_particle = 0 */

  for (k = 0; k < 3; k++) {
    face_normal[k] = cs_glob_mesh_quantities->b_face_normal[3*face_id+k];
    face_cog[k] = cs_glob_mesh_quantities->b_face_cog[3*face_id+k];
  }

  cs_real_t face_area  = _get_norm(face_normal);

  /* Saving of particle impacting velocity */
  if (*iangbd > 0 || *ivitbd > 0) {
    norm_vit = _get_norm(particle.velocity);
    for (k = 0; k < 3; k++)
      compo_vit[k] = particle.velocity[k];
  }

  tmp = 0.0;
  for (k = 0; k < 3; k++)
    tmp += depl[k] * face_normal[k];

  if (fabs(tmp) < 1e-15) {
    const char msg[] = " Error during boundary treatment.\nPiece of trajectography inside the boundary faces.\n";
    _manage_error(failsafe_mode,
                  particle,
                  &error,
                  &n_failed_particles,
                  &failed_particle_weight,
                  msg);
  }

  /* Petit rappel de geometrie 3D :
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     1)  equation d'un plan de normal (a,b,c) :
     a*x + b*y + c*z + d = 0
     2)  equation d'une droite qui passe par P et Q :
     x = XP + (XQ-XP) * AA
     y = YP + (YQ-YP) * AA
     z = ZP + (ZQ-ZP) * AA
     ou AA est un parametre qui varie dans l'ensemble des reels
  */

  for (k = 0; k < 3; k++)
    abs_curv +=
      face_normal[k]*face_cog[k] - face_normal[k]*prev_particle.coord[k];
  abs_curv /= tmp;

  for (k = 0; k < 3; k++)
    intersect_pt[k] = depl[k] * abs_curv + prev_particle.coord[k];

  if (   bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_ISORTL
      || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IENTRL
      || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_OUT;

    if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1) {
      _particle_set->n_part_dep += 1;
      _particle_set->weight_dep += particle.stat_weight;
    }

    bdy_conditions->particle_flow_rate[boundary_zone]
      -= particle.stat_weight * particle.mass;

    /* FIXME: For post-processing by trajectory purpose */

    for (k = 0; k < 3; k++)
      particle.coord[k] = intersect_pt[k];
  }

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO2) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle.cur_cell_num =  - particle.cur_cell_num; /* Store a negative value */

    _particle_set->n_part_dep += 1;
    _particle_set->weight_dep += particle.stat_weight;

    particle_state = CS_LAGR_PART_STICKED;

    /* For post-processing purpose */

    for (k = 0; k < 3; k++) {
      particle.coord[k] = intersect_pt[k];
      particle.velocity[k] = 0.0;
      particle.velocity_seen[k] = 0.0;
    }
  }

  /* FIXME: IDEPFA not yet implemented */

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPFA)
    bft_error(__FILE__, __LINE__, 0,
              " Boundary condition IDEPFA not yet implemented.\n");

  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IREBOL) {

    cs_real_t  face_norm[3] = {face_normal[0]/face_area,
                               face_normal[1]/face_area,
                               face_normal[2]/face_area};

    move_particle = 1;
    particle_state = CS_LAGR_PART_TO_SYNC;
    particle.cur_cell_num = cs_glob_mesh->b_face_cells[face_id];

    for (k = 0; k < 3; k++)
      prev_particle.coord[k] = intersect_pt[k];

    /* Modify the ending point. */

    for (k = 0; k < 3; k++)
      depl[k] = particle.coord[k] - intersect_pt[k];

    tmp = CS_ABS(_get_dot_prod(depl, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle.coord[k] -= tmp * face_norm[k];

    /* Modify particle velocity and velocity seen */

    tmp = CS_ABS(_get_dot_prod(particle.velocity, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle.velocity[k] -= tmp * face_norm[k];

    tmp = CS_ABS(_get_dot_prod(particle.velocity_seen, face_norm));
    tmp *= 2.0;

    for (k = 0; k < 3; k++)
      particle.velocity_seen[k] -= tmp * face_norm[k];

  }

  /* FIXME: Coal particles boundary treatment */
  else if (bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IENCRL)
    bft_error(__FILE__, __LINE__, 0,
              " Boundary condition IENCRL not yet implemented.\n");

  /* FIXME: JBORD* (user-defined boundary condition) not yet implemented
     nor defined by a macro */
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Boundary condition %d not recognized.\n"),
              bdy_conditions->b_zone_natures[boundary_zone]);

  /* FIXME: Post-treatment not yet implemented... */

  /* Return pointer */

  *p_prev_particle = prev_particle;
  *p_particle = particle;
  *p_move_particle = move_particle;
  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;

  if (*iensi3 > 0) {
    if  (   bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO1
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPO2
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IDEPFA
         || bdy_conditions->b_zone_natures[boundary_zone] == CS_LAGR_IREBOL) {

      /* Number of particle-boundary interactions  */
      if (*inbrbd > 0)
        boundary_stat[(*inbr -1) * nfabor + face_id] +=  particle.stat_weight;

      /* Particulate boundary mass flux */
      if (*iflmbd > 0)
        boundary_stat[(*iflm -1) * nfabor + face_id] +=
          particle.stat_weight * particle.mass / face_area;

      /* Particle impact angle and velocity*/
      if (*iangbd > 0) {
        cs_real_t imp_ang = acos(_get_dot_prod(compo_vit, face_normal)
                                 / (face_area * norm_vit));
        boundary_stat[(*iang -1) * nfabor + face_id]
          += imp_ang * particle.stat_weight;
      }

      if (*ivitbd > 0)
        boundary_stat[(*ivit -1) * nfabor + face_id]
          += norm_vit * particle.stat_weight;

      /* User statistics management. By defaut, set to zero */
      if (*nusbor > 0)
        for (int n1 = 0; n1 < *nusbor; n1++)
            boundary_stat[(iusb[n1] -1) * nfabor + face_id] = 0.0;
    }
  }

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Move locally a particle as far as it is possible while keeping on the same
 * rank.
 *
 * parameters:
 *  p_prev_particle <->  pointer on a cs_lagr_particle_t structure
 *  p_particle      <->  pointer on a cs_lagr_particle_t structure
 *  displacement_step_id     <-- id of displacement step
 *  scheme_order    -->  current order of the scheme used for the lagrangian
 *                       algorithm (1 or 2)
 *  failsafe_mode   -->  with (0) / without (1) failure capability
 *  p_n_failed_particles     <->  number of failed particles
 *  p_failed_particle_weight <->  stat. weight of the failed particles
 *
 * returns:
 *  a state associated to the status of the particle (treated, to be deleted,
 *  to be synchonised)
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_local_propagation(cs_lagr_particle_t     *p_prev_particle,
                   cs_lagr_particle_t     *p_particle,
                   int                     displacement_step_id,
                   cs_lnum_t               scheme_order,
                   cs_lnum_t               failsafe_mode,
                   cs_real_t               boundary_stat[],
                   cs_lnum_t              *p_n_failed_particles,
                   cs_real_t              *p_failed_particle_weight,
                   const cs_lnum_t         *iensi3,
                   const cs_lnum_t         *nvisbr,
                   const cs_lnum_t         *inbr,
                   const cs_lnum_t         *inbrbd,
                   const cs_lnum_t         *iflm,
                   const cs_lnum_t         *iflmbd,
                   const cs_lnum_t         *iang,
                   const cs_lnum_t         *iangbd,
                   const cs_lnum_t         *ivit,
                   const cs_lnum_t         *ivitbd,
                   const cs_lnum_t         *nusbor,
                   cs_lnum_t               iusb[],
                   cs_real_t               visc_length[],
                   cs_real_t               dlgeo[],
                   cs_real_t               rtp[],
                   const cs_lnum_t  *iu,
                   const cs_lnum_t  *iv,
                   const cs_lnum_t  *iw,
                   cs_lnum_t               *idepst)
{
  cs_lnum_t  i, j, k;
  cs_real_t  depl[3];

  cs_lnum_t  error = 0;
  cs_lnum_t  n_loops = displacement_step_id;
  cs_lnum_t  move_particle = CS_LAGR_PART_MOVE_ON;
  cs_lnum_t  particle_state = CS_LAGR_PART_TO_SYNC;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_bdy_condition_t  *bdy_conditions = _lagr_bdy_conditions;
  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;
  cs_lnum_t  *face_connect = builder->face_connect_buffer;

  cs_lagr_particle_t  particle = *p_particle;
  cs_lagr_particle_t  prev_particle = *p_prev_particle;
  cs_lnum_t  n_failed_particles = *p_n_failed_particles;
  cs_real_t  failed_particle_weight =*p_failed_particle_weight;
  cs_lnum_t  new_face_num = particle.last_face_num;

  for (k = 0; k < 3; k++)
    depl[k] = particle.coord[k] - prev_particle.coord[k];

  if (fabs(depl[0]) < 1e-15 && fabs(depl[1]) < 1e-15 && fabs(depl[2]) < 1e-15) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_TREATED;

  }

  /*  particle_state is defined at the top of this file */

  while (move_particle == CS_LAGR_PART_MOVE_ON) {

    cs_lnum_t  cur_cell_id = particle.cur_cell_num - 1;
    cs_lnum_t  old_face_num = new_face_num;
    cs_lnum_t  start = cell_face_idx[cur_cell_id];
    cs_lnum_t  end =  cell_face_idx[cur_cell_id+1];

    assert(cur_cell_id < mesh->n_cells);
    assert(cur_cell_id > -1);

    n_loops++;

    if (n_loops > CS_LAGR_MAX_PROPAGATION_LOOPS) { /* Manage error */

      const char msg[] = "Max number of loops reached in particle displacement";
      _manage_error(failsafe_mode,
                    particle,
                    &error,
                    &n_failed_particles,
                    &failed_particle_weight,
                    msg);

      move_particle  = CS_LAGR_PART_MOVE_OFF;
      particle_state = CS_LAGR_PART_ERR;

    }

/*Treatment for particles which change rank*/
    if (*idepst > 0 && particle.yplus < 0.) {
      _test_wall_cell(&particle,visc_length,dlgeo);

      if (particle.yplus < 100.) {

        cs_real_t flow_velo_x = rtp[ (particle.cur_cell_num - 1) + (*iu - 1) * mesh->n_cells_with_ghosts ];
        cs_real_t flow_velo_y = rtp[ (particle.cur_cell_num - 1) + (*iv - 1) * mesh->n_cells_with_ghosts ];
        cs_real_t flow_velo_z = rtp[ (particle.cur_cell_num - 1) + (*iw - 1) * mesh->n_cells_with_ghosts ];

        /* e1 (normal) vector coordinates */
        cs_real_t e1_x = dlgeo[particle.close_face_id];
        cs_real_t e1_y = dlgeo[particle.close_face_id + (mesh->n_b_faces)];
        cs_real_t e1_z = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 2];

        /* e2 vector coordinates */
        cs_real_t e2_x = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 7];
        cs_real_t e2_y = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 8];
        cs_real_t e2_z = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 9];

        /* e3 vector coordinates */
        cs_real_t e3_x = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 10];
        cs_real_t e3_y = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 11];
        cs_real_t e3_z = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 12];

        /* V_n * e1 */

        cs_real_t v_n_e1[3] = { particle.velocity_seen[0] * e1_x,
                                particle.velocity_seen[0] * e1_y,
                                particle.velocity_seen[0] * e1_z};


        /* (U . e2) * e2 */

        cs_real_t flow_e2 = flow_velo_x * e2_x +
          flow_velo_y * e2_y +
          flow_velo_z * e2_z;

        cs_real_t u_times_e2[3] = {flow_e2 * e2_x,
                                   flow_e2 * e2_y,
                                   flow_e2 * e2_z};

        /* (U . e3) * e3 */

        cs_real_t flow_e3 = flow_velo_x * e3_x +
          flow_velo_y * e3_y +
          flow_velo_z * e3_z;

        cs_real_t u_times_e3[3] = {flow_e3 * e3_x,
                                   flow_e3 * e3_y,
                                   flow_e3 * e3_z};
        /* Update of the flow seen velocity */

        particle.velocity_seen[0] =  v_n_e1[0] + u_times_e2[0] + u_times_e3[0];
        particle.velocity_seen[1] =  v_n_e1[1] + u_times_e2[1] + u_times_e3[1];
        particle.velocity_seen[2] =  v_n_e1[2] + u_times_e2[2] + u_times_e3[2];
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
                                &error);

        if (error == 1) {

          const char msg[]
            = "Error during the particle displacement (Interior face)";
          _manage_error(failsafe_mode,
                        particle,
                        &error,
                        &n_failed_particles,
                        &failed_particle_weight,
                        msg);

          move_particle = CS_LAGR_PART_MOVE_OFF;
          particle_state = CS_LAGR_PART_ERR;

        }
        else if (indian == 1) { /* Particle moves to the neighbor cell
                                   through the current face "face_num" */

          cs_lnum_t  cell_num1 = mesh->i_face_cells[2*face_id];
          cs_lnum_t  cell_num2 = mesh->i_face_cells[2*face_id+1];

          new_face_num = face_num;

          if (particle.cur_cell_num == cell_num1)
            particle.cur_cell_num = cell_num2;

          else
            particle.cur_cell_num = cell_num1;

          if (particle.cur_cell_num > mesh->n_cells) {

            particle_state = CS_LAGR_PART_TO_SYNC;
            move_particle = CS_LAGR_PART_ERR;

            if (*idepst > 0 && particle.yplus < 100.) {

              cs_real_t x_p_q = particle.coord[0] - prev_particle.coord[0];
              cs_real_t y_p_q = particle.coord[1] - prev_particle.coord[1];
              cs_real_t z_p_q = particle.coord[2] - prev_particle.coord[2];

              cs_real_t face_normal[3], face_cog[3];

              for (k = 0; k < 3; k++) {
                face_normal[k] = cs_glob_mesh_quantities->i_face_normal[3*face_id+k];
                face_cog[k] = cs_glob_mesh_quantities->i_face_cog[3*face_id+k];
              }

              cs_real_t aa = x_p_q * face_normal[0] +
                y_p_q * face_normal[1] +
                z_p_q * face_normal[2];


              cs_real_t bb = ( face_normal[0] * face_cog[0] +
                               face_normal[1] * face_cog[1] +
                               face_normal[2] * face_cog[2] -
                               face_normal[0] *  prev_particle.coord[0] -
                               face_normal[1] *  prev_particle.coord[1] -
                               face_normal[2] *  prev_particle.coord[2] ) / aa;

              cs_real_t xk =  prev_particle.coord[0] + bb * x_p_q;
              cs_real_t yk =  prev_particle.coord[1] + bb * y_p_q;
              cs_real_t zk =  prev_particle.coord[2] + bb * z_p_q;


              cs_real_t* xyzcen = cs_glob_mesh_quantities->cell_cen;

              particle.coord[0] =  xk + 1e-8 * (xyzcen[ 3* (particle.cur_cell_num - 1)] - xk);
              particle.coord[1] =  yk + 1e-8 * (xyzcen[ 3* (particle.cur_cell_num - 1) + 1] - yk);
              particle.coord[2] =  zk + 1e-8 * (xyzcen[ 3* (particle.cur_cell_num - 1) + 2] - zk);

              /*Marking of particles */

              particle.yplus = -particle.yplus;

              /*Saving of dot product */

              cs_real_t bdy_normal[3];


              for (k = 0; k < 3; k++)

                bdy_normal[k] = cs_glob_mesh_quantities->b_face_normal[3 * particle.close_face_id + k];

              cs_real_t area  = _get_norm(bdy_normal);

              cs_real_t  face_norm[3] = {bdy_normal[0]/area,
                                         bdy_normal[1]/area,
                                         bdy_normal[2]/area};



              particle.velocity_seen[0] = particle.velocity_seen[0] * face_norm[0] +
                particle.velocity_seen[1] * face_norm[1] +
                particle.velocity_seen[2] * face_norm[2];

            }

          }
          else {

            /* Specific treatment for the particle deposition model */

            if (*idepst > 0) {

              cs_lnum_t save_close_face_id = particle.close_face_id;
              cs_real_t save_yplus = particle.yplus;

              /* Wall cell detection */

              _test_wall_cell(&particle,visc_length,dlgeo);

              if (save_yplus < 100.) {

                cs_real_t x_p_q = particle.coord[0] - prev_particle.coord[0];
                cs_real_t y_p_q = particle.coord[1] - prev_particle.coord[1];
                cs_real_t z_p_q = particle.coord[2] - prev_particle.coord[2];

                cs_real_t face_normal[3], face_cog[3];

                for (k = 0; k < 3; k++) {
                  face_normal[k]
                    = cs_glob_mesh_quantities->i_face_normal[3*face_id+k];
                  face_cog[k] = cs_glob_mesh_quantities->i_face_cog[3*face_id+k];
                }

                cs_real_t aa =   x_p_q * face_normal[0]
                               + y_p_q * face_normal[1]
                               + z_p_q * face_normal[2];


                cs_real_t bb = (  face_normal[0] * face_cog[0]
                                + face_normal[1] * face_cog[1]
                                + face_normal[2] * face_cog[2]
                                - face_normal[0] *  prev_particle.coord[0]
                                - face_normal[1] *  prev_particle.coord[1]
                                - face_normal[2] *  prev_particle.coord[2]) / aa;

                cs_real_t xk =  prev_particle.coord[0] + bb * x_p_q;
                cs_real_t yk =  prev_particle.coord[1] + bb * y_p_q;
                cs_real_t zk =  prev_particle.coord[2] + bb * z_p_q;

                cs_real_t* xyzcen = cs_glob_mesh_quantities->cell_cen;

                particle.coord[0] = xk + 1e-8 * (xyzcen[3* (particle.cur_cell_num - 1)] - xk);
                particle.coord[1] = yk + 1e-8 * (xyzcen[3* (particle.cur_cell_num - 1) + 1] - yk);
                particle.coord[2] = zk + 1e-8 * (xyzcen[3* (particle.cur_cell_num - 1) + 2] - zk);

                /* Second test with the new particle position */

                _test_wall_cell(&particle,visc_length,dlgeo);


                if (particle.yplus < 100.e0 ) {

                  cs_real_t flow_velo_x = rtp[ (particle.cur_cell_num - 1) + (*iu - 1) * mesh->n_cells_with_ghosts ];
                  cs_real_t flow_velo_y = rtp[ (particle.cur_cell_num - 1) + (*iv - 1) * mesh->n_cells_with_ghosts ];
                  cs_real_t flow_velo_z = rtp[ (particle.cur_cell_num - 1) + (*iw - 1) * mesh->n_cells_with_ghosts ];

                  /* The particle is still in the boundary layer */

                  cs_real_t old_bdy_normal[3];

                  for (k = 0; k < 3; k++)
                    old_bdy_normal[k] = cs_glob_mesh_quantities->b_face_normal[3 * save_close_face_id + k];

                  cs_real_t old_area  = _get_norm(old_bdy_normal);

                  cs_real_t  old_face_norm[3] = {old_bdy_normal[0]/old_area,
                                                 old_bdy_normal[1]/old_area,
                                                 old_bdy_normal[2]/old_area};

                  cs_real_t old_fl_seen_norm
                    =   particle.velocity_seen[0] * old_face_norm[0]
                      + particle.velocity_seen[1] * old_face_norm[1]
                      + particle.velocity_seen[2] * old_face_norm[2];

                  /* e1 (normal) vector coordinates */
                  cs_real_t e1_x = dlgeo[particle.close_face_id];
                  cs_real_t e1_y = dlgeo[particle.close_face_id + (mesh->n_b_faces)];
                  cs_real_t e1_z = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 2];

                  /* e2 vector coordinates */
                  cs_real_t e2_x = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 7];
                  cs_real_t e2_y = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 8];
                  cs_real_t e2_z = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 9];

                  /* e3 vector coordinates */
                  cs_real_t e3_x = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 10];
                  cs_real_t e3_y = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 11];
                  cs_real_t e3_z = dlgeo[particle.close_face_id + (mesh->n_b_faces) * 12];


                  /* V_n * e1 */

                  cs_real_t v_n_e1[3] = {old_fl_seen_norm * e1_x,
                                         old_fl_seen_norm * e1_y,
                                         old_fl_seen_norm * e1_z};

                  /* (U . e2) * e2 */

                  cs_real_t flow_e2 =   flow_velo_x * e2_x
                                      + flow_velo_y * e2_y
                                      + flow_velo_z * e2_z;

                  cs_real_t u_times_e2[3] = {flow_e2 * e2_x,
                                             flow_e2 * e2_y,
                                             flow_e2 * e2_z};

                  /* (U . e3) * e3 */

                  cs_real_t flow_e3 =   flow_velo_x * e3_x
                                      + flow_velo_y * e3_y
                                      + flow_velo_z * e3_z;

                  cs_real_t u_times_e3[3] = {flow_e3 * e3_x,
                                             flow_e3 * e3_y,
                                             flow_e3 * e3_z};

                  /* Update of the flow seen velocity */

                  particle.velocity_seen[0] = v_n_e1[0] + u_times_e2[0] + u_times_e3[0];
                  particle.velocity_seen[1] = v_n_e1[1] + u_times_e2[1] + u_times_e3[1];
                  particle.velocity_seen[2] = v_n_e1[2] + u_times_e2[2] + u_times_e3[2];
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
                                &error);

        face_num = -face_num;

        if (error == 1) {

          const char msg[]
            = "Error during the particle displacement (Boundary face)";
          _manage_error(failsafe_mode,
                        particle,
                        &error,
                        &n_failed_particles,
                        &failed_particle_weight,
                        msg);

          move_particle = CS_LAGR_PART_MOVE_OFF;
          particle_state = CS_LAGR_PART_ERR;

        }
        else if (indian == 1) { /* Particle moves to the neighbor cell
                                   through the current face "face_num" */

          /* particle / boundary condition interaction
             1 - modify particle.cur_cell_num : 0 or boundary_cell_num
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
            = _bdy_treatment(&prev_particle,
                             &particle,
                             face_num,
                             boundary_stat,
                             bdy_conditions->b_face_zone_num[face_num-1]-1,
                             failsafe_mode,
                             &move_particle,
                             &n_failed_particles,
                             &failed_particle_weight,
                             iensi3,
                             nvisbr,
                             inbr,
                             inbrbd,
                             iflm,
                             iflmbd,
                             iang,
                             iangbd,
                             ivit,
                             ivitbd,
                             nusbor,
                             iusb);

          if (scheme_order == 2)
            particle.switch_order_1 = CS_LAGR_SWITCH_ON;

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

  particle.last_face_num = new_face_num;

  /* Return pointers */

  *p_prev_particle = prev_particle;
  *p_particle = particle;
  *p_n_failed_particles = n_failed_particles;
  *p_failed_particle_weight = failed_particle_weight;

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Exchange counters on the number of particles to send and to receive
 *
 * parameters:
 *  halo        -->  pointer to a cs_halo_t structure
 *  lag_halo    -->  pointer to a cs_lagr_halo_t structure
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
 *  halo        -->  pointer to a cs_halo_t structure
 *  lag_halo    -->  pointer to a cs_lagr_halo_t structure
 *----------------------------------------------------------------------------*/

static void
_exchange_particles(const cs_halo_t          *halo,
                    cs_lagr_halo_t           *lag_halo)
{
  cs_lnum_t  shift;

  cs_lagr_particle_t  *recv_buf = NULL, *send_buf = NULL;

  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

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
        recv_buf = &(lag_halo->recv_buf->particles[shift]);

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

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(cs_glob_mpi_comm);

    /* Send data to distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank] != local_rank) {

        shift = lag_halo->send_shift[rank];
        if (lag_halo->send_count[rank] == 0)
          send_buf = NULL;
        else
          send_buf = &(lag_halo->send_buf->particles[shift]);

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

      cs_lnum_t  i;

      cs_lnum_t  recv_shift = lag_halo->recv_shift[local_rank_id];
      cs_lnum_t  send_shift = lag_halo->send_shift[local_rank_id];

      assert(   lag_halo->recv_count[local_rank_id]
             == lag_halo->send_count[local_rank_id]);

      for (i = 0; i < lag_halo->send_count[local_rank_id]; i++)
        lag_halo->recv_buf->particles[recv_shift + i] =
          lag_halo->send_buf->particles[send_shift + i];

    }
  }

}

/*----------------------------------------------------------------------------
 * Update a cs_lagr_particle_set structure after a parallel and/or periodic
 * exchange between communicating ranks.
 *
 * parameters:
 *  n_recv_particles  -->  number of particles received during the exchange
 *  lag_halo          -->  pointer to a cs_lagr_halo_t structure
 *  set               -->  set of particles to update
 *----------------------------------------------------------------------------*/

static void
_update_particle_set(cs_lnum_t                n_recv_particles,
                     cs_lagr_halo_t          *lag_halo,
                     cs_lagr_particle_set_t  *set)
{
  cs_lnum_t  i, new_id;

  for (i = 0; i < n_recv_particles; i++) {

    cs_lagr_particle_t  new_part = lag_halo->recv_buf->particles[i];

    new_id = set->first_free_id;

    if (set->particles[set->first_free_id].next_id != -1) {
      set->first_free_id = set->particles[set->first_free_id].next_id;
    }
    else {
      set->first_free_id = set->first_free_id + 1 ;
      set->particles[set->first_free_id].next_id = -1;
    }

    /* Add new_part at the beginning of the "used list"
       Update first_used_id */

    if (set->first_used_id != -1)
      set->particles[set->first_used_id].prev_id = new_id;

    new_part.prev_id = -1;
    new_part.next_id = set->first_used_id;

    set->first_used_id = new_id;
    set->particles[new_id] = new_part;
  }

}

/*----------------------------------------------------------------------------
 * Synchronize particle displacement over the ranks.
 *
 * parameters:
 *  lag_halo  -->  pointer to a cs_lagr_halo_t structure
 *  set       -->  set of particles to update
 *----------------------------------------------------------------------------*/

static void
_sync_particle_sets(cs_lagr_halo_t           *lag_halo,
                    cs_interface_set_t       *face_ifs,
                    cs_lagr_particle_set_t   *prev_set,
                    cs_lagr_particle_set_t   *cur_set)
{
  cs_lnum_t  i, j, k, tr_id, rank, shift, ghost_id;
  cs_real_t vect_in[3];
  cs_real_t matrix[3][4];

  cs_lnum_t  n_recv_particles = 0;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const fvm_periodicity_t *periodicity = mesh->periodicity;

  for (i = 0; i < halo->n_c_domains; i++){
    n_recv_particles += lag_halo->recv_count[i];
    lag_halo->send_count[i] = 0;
  }

  /* Fill send_buf for previous particle set */

  for (i = 0, j = cur_set->first_used_id; i < cur_set->n_particles; i++) {

    cs_lagr_particle_t  cur_part = cur_set->particles[j];
    cs_lagr_particle_t  prev_part = prev_set->particles[j];

    if (cur_part.state == CS_LAGR_PART_TO_SYNC) {

      ghost_id = cur_part.cur_cell_num - halo->n_local_elts - 1;
      rank = lag_halo->rank[ghost_id];
      tr_id = lag_halo->transform_id[ghost_id];
      shift = lag_halo->send_shift[rank] + lag_halo->send_count[rank];

      if (tr_id >= 0) { /* Periodicity treatment */

        fvm_periodicity_type_t  perio_type =
          fvm_periodicity_get_type(periodicity, tr_id);

        int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, tr_id);

        fvm_periodicity_get_matrix(periodicity, rev_id, matrix);

        /* Apply transformation to the coordinates in any case */

        _apply_vector_transfo(matrix,
                              prev_part.coord,
                              prev_part.coord);

        /* Apply rotation to velocity vectors in case of rotation */

        if (perio_type >= FVM_PERIODICITY_ROTATION) {

          /* Rotation of the velocity */

          for (k = 0; k < 3; k++)
            vect_in[k] = prev_part.velocity[k];

          _apply_vector_rotation(matrix,
                                 vect_in[0],
                                 vect_in[1],
                                 vect_in[2],
                                 &prev_part.velocity[0],
                                 &prev_part.velocity[1],
                                 &prev_part.velocity[2]);

          /* Rotation of the velocity seen */

          for (k = 0; k < 3; k++)
            vect_in[k] = prev_part.velocity_seen[k];

          _apply_vector_rotation(matrix,
                                 vect_in[0],
                                 vect_in[1],
                                 vect_in[2],
                                 &prev_part.velocity_seen[0],
                                 &prev_part.velocity_seen[1],
                                 &prev_part.velocity_seen[2]);

        } /* Specific treatment in case of rotation for the velocities */

      } /* End of periodicity treatment */

      lag_halo->send_buf->particles[shift] = prev_part;
      lag_halo->send_count[rank] += 1;
      /* Pick out the particle from the particle_set */

      _remove_particle(prev_set, prev_part, j);

    } /* TO_SYNC */

    j = cur_part.next_id;

  } /* End of loop on particles */

    /* Exchange particles */
  _exchange_particles(halo, lag_halo);

  /* Update the particle set after the exchange of particles between ranks */

  _update_particle_set(n_recv_particles, lag_halo, prev_set);

  /* Exchange current particle set */

  for (i = 0; i < halo->n_c_domains; i++)
    lag_halo->send_count[i] = 0;

  /* Fill send_buf  */

  for (i = 0, j = cur_set->first_used_id; i < cur_set->n_particles; i++) {

    cs_lagr_particle_t  cur_part = cur_set->particles[j];

    if (cur_part.state == CS_LAGR_PART_TO_SYNC) {

      ghost_id = cur_part.cur_cell_num - halo->n_local_elts - 1;
      rank = lag_halo->rank[ghost_id];
      tr_id = lag_halo->transform_id[ghost_id];
      cur_part.cur_cell_num = lag_halo->dist_cell_num[ghost_id];

      shift = lag_halo->send_shift[rank] + lag_halo->send_count[rank];

      /* Update if needed last_face_num */

      if (tr_id >= 0) // Same initialization as in previous algo.
      {
        cur_part.last_face_num = 0;

      } else {

        if (cs_glob_n_ranks > 1) {

          assert(face_ifs != NULL);

          {
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

              id = cs_search_binary(n_entities, cur_part.last_face_num - 1, local_num);

              if (id == -1)
                bft_error(__FILE__, __LINE__, 0,
                          _(" Cannot find the relative distant face num.\n"));

              dist_num = cs_interface_get_match_ids(interface);
              cur_part.last_face_num = dist_num[id] + 1;
            }

          }

        }
      }


      if (tr_id >= 0) { /* Periodicity treatment */

        fvm_periodicity_type_t  perio_type =
          fvm_periodicity_get_type(periodicity, tr_id);

        int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, tr_id);
                fvm_periodicity_get_matrix(periodicity, rev_id, matrix);

        /* Apply transformation to the coordinates in any case */

        _apply_vector_transfo(matrix,
                              cur_part.coord,
                              cur_part.coord);

        /* Apply rotation to velocity vectors in case of rotation */

        if (perio_type >= FVM_PERIODICITY_ROTATION) {

          /* Rotation of the velocity */

          for (k = 0; k < 3; k++)
            vect_in[k] = cur_part.velocity[k];

          _apply_vector_rotation(matrix,
                                 vect_in[0],
                                 vect_in[1],
                                 vect_in[2],
                                 &cur_part.velocity[0],
                                 &cur_part.velocity[1],
                                 &cur_part.velocity[2]);

          /* Rotation of the velocity seen */

          for (k = 0; k < 3; k++)
            vect_in[k] = cur_part.velocity_seen[k];

          _apply_vector_rotation(matrix,
                                 vect_in[0],
                                 vect_in[1],
                                 vect_in[2],
                                 &cur_part.velocity_seen[0],
                                 &cur_part.velocity_seen[1],
                                 &cur_part.velocity_seen[2]);

        } /* Specific treatment in case of rotation for the velocities */

      } /* End of periodicity treatment */

      lag_halo->send_buf->particles[shift] = cur_part;
      lag_halo->send_count[rank] += 1;
      /* Pick out the particle from the particle_set */

      _remove_particle(cur_set, cur_part, j);

    } /* TO_SYNC */

    j = cur_part.next_id;

  } /* End of loop on particles */

    /* Exchange particles */
  _exchange_particles(halo, lag_halo);

  /* Update the particle set after the exchange of particles between ranks */

  _update_particle_set(n_recv_particles, lag_halo, cur_set);
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

  cs_real_t delta_weight;

  cs_lnum_t  n_recv_particles = 0, n_send_particles = 0;
  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_halo_t  *lag_halo = builder->halo;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;

  assert(set->n_particles == prev_set->n_particles);
  assert(set->first_used_id == prev_set->first_used_id);

  /* Check this because this assumption is done. */

  /* Initialization */

  for (i = 0; i < halo->n_c_domains; i++) {
    lag_halo->send_count[i] = 0;
    lag_halo->recv_count[i] = 0;
  }

  /* Loop on particles to count number of particles to send on each rank */

  for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

    cs_lagr_particle_t  cur_part = set->particles[j];

    if (cur_part.state == CS_LAGR_PART_TO_SYNC) {

      ghost_id = cur_part.cur_cell_num - mesh->n_cells - 1;
      assert(ghost_id >= 0);
      lag_halo->send_count[lag_halo->rank[ghost_id]] += 1;

    }

    j = cur_part.next_id;

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

  /* Resize particle set only if needed */

  _resize_particle_set(&(lag_halo->send_buf), n_send_particles);
  _resize_particle_set(&(lag_halo->recv_buf), n_recv_particles);

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

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate cs_lagr_particle_set_t structure and initialize useful buffers.
 *
 * parameters:
 *   n_particles_max <--  local max. number of particles
 *   iphyla          <--  kind of physics used for the lagrangian approach
 *   nvls            <--  number of user-defined variables
 *   nbclst          <--  number of stat. class to study sub-set of particles
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagbeg, LAGBEG)(const cs_int_t  *n_particles_max,
                          const cs_int_t  *iphyla,
                          const cs_int_t  *nvls,
                          const cs_int_t  *nbclst)
{
  cs_lnum_t  i;

  /* Initialize global parameter relative to the lagrangian module */

  cs_glob_lagr_param.physic_mode = *iphyla;
  cs_glob_lagr_param.n_user_variables = *nvls;
  cs_glob_lagr_param.n_stat_classes = *nbclst;

  /* Initialize particle set : prev and current */

  _particle_set = _create_particle_set(*n_particles_max);
  _prev_particle_set = _create_particle_set(*n_particles_max);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER CREATION\n");
  _dump_particle_set(_particle_set);
  bft_printf("\n PREV PARTICLE SET AFTER CREATION\n");
  _dump_particle_set(_prev_particle_set);
#endif

  /* Initialization */

  for (i = 0; i < _particle_set->n_particles_max; i++) {

    _particle_set->particles[i].switch_order_1 = CS_LAGR_SWITCH_OFF;
    _prev_particle_set->particles[i].switch_order_1 = CS_LAGR_SWITCH_OFF;

    _particle_set->particles[i].state = CS_LAGR_PART_TO_SYNC;
    _prev_particle_set->particles[i].state = CS_LAGR_PART_TO_SYNC;

#if defined(HAVE_MPI)
    _particle_set->particles[i].cell_rank = cs_glob_rank_id;
    _prev_particle_set->particles[i].cell_rank = cs_glob_rank_id;
#endif

  }

  /* Initialize builder */

  _particle_track_builder = _init_track_builder(*n_particles_max);

  /* Create all useful MPI_Datatype */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _CS_MPI_PARTICLE = _define_particle_datatype();
    _CS_MPI_COAL_PARTICLE = _define_coal_particle_datatype();
    _CS_MPI_HEAT_PARTICLE = _define_heat_particle_datatype();
    _CS_MPI_AUX_PARTICLE = _define_aux_particle_datatype();
    }
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
CS_PROCF (prtget, PRTGET)(const cs_lnum_t   *nbpmax,  /* n_particles max. */
                          const cs_lnum_t   *nbpart,  /* number of current particles */
                          const cs_real_t   *dnbpar,  /* particle total weight */
                          cs_lnum_t          liste[],
                          cs_lnum_t         *nbvis,
                          const cs_real_t    ettp[],
                          const cs_real_t    ettpa[],
                          const cs_lnum_t    itepa[],
                          const cs_real_t    tepa[],
                          const cs_lnum_t    ibord[],
                          const cs_lnum_t    indep[],
                          const cs_lnum_t   *jisor,
                          const cs_lnum_t   *jgnum,
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
                          cs_lnum_t         *idepst
)
{
  cs_lnum_t  i, id;

  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;

  assert(*nbpmax == set->n_particles_max); /* Up to now, we don't manage
                                              a mofification of nbpmax */
  set->n_particles =*nbpart;
  prev_set->n_particles =*nbpart;

  set->weight = 0.0;
  prev_set->weight =0.0;

  set->n_part_out = 0;
  prev_set->n_part_out = 0;

  set->n_part_dep = 0;
  prev_set->n_part_dep = 0;

  set->weight_out = 0.0;
  prev_set->weight_out = 0.0;

  set->weight_dep = 0.0;
  prev_set->weight_dep = 0.0;

  set->n_failed_part = 0;
  prev_set->n_failed_part = 0;

  set->weight_failed = 0.0;
  prev_set->weight_failed = 0.0;

  /* When we receive particles from the FORTRAN, we keep a compact
     storage of particles */

  if (*nbpart > 0) {
    set->first_used_id = 0;
    prev_set->first_used_id = 0;
  }
  else {
    set->first_used_id = -1;
    prev_set->first_used_id = -1;
  }

  set->first_free_id = *nbpart;
  set->particles[set->first_free_id].next_id = -1;

  prev_set->first_free_id = *nbpart;
  prev_set->particles[prev_set->first_free_id].next_id = -1;

  /* Fill set and prev_set structures */

  for (i = 0; i < *nbpart; i++) {

    cs_lagr_particle_t  cur_part = set->particles[i];
    cs_lagr_particle_t  prev_part = prev_set->particles[i];

    if (i > 0) {
      cur_part.prev_id = i-1;
      prev_part.prev_id = i-1;
    }
    else { /* Not defined */
      cur_part.prev_id = -1;
      prev_part.prev_id = -1;
    }

    if (i < *nbpart - 1) {
      cur_part.next_id = i+1;
      prev_part.next_id = i+1;
    }
    else {  /* Not defined */
      cur_part.next_id = -1;
      prev_part.next_id = -1;
    }

    /* Global number (not true in parallel => MPI_Scan to do) Pb in POST-PROCESSING */
    cur_part.global_num = i + 1;
    prev_part.global_num = i + 1;

    cur_part.cur_cell_num = itepa[i + (*jisor-1) * (*nbpmax)];
    prev_part.cur_cell_num = indep[i];

    cur_part.global_num = itepa[i + (*jgnum-1) * (*nbpmax)];
    prev_part.global_num = itepa[i + (*jgnum-1) * (*nbpmax)];

    if (cur_part.cur_cell_num < 0)
      cur_part.state = CS_LAGR_PART_STICKED;
    else if (cur_part.cur_cell_num == 0)
      cur_part.state = CS_LAGR_PART_TO_DELETE;
    else
      cur_part.state = CS_LAGR_PART_TO_SYNC;

    cur_part.last_face_num = 0;
    prev_part.last_face_num = 0;

#if defined(HAVE_MPI)
    cur_part.cell_rank = cs_glob_rank_id;
    prev_part.cell_rank = cs_glob_rank_id;
#endif

    cur_part.switch_order_1 = ibord[i];
    prev_part.switch_order_1 = ibord[i];

    id = (*jrpoi-1) * (*nbpmax) + i;
    cur_part.stat_weight = tepa[id];
    prev_part.stat_weight = tepa[id];

    id = (*jrtsp-1) * (*nbpmax) + i;
    cur_part.residence_time = tepa[id];
    prev_part.residence_time = tepa[id];

    id = (*jmp-1) * (*nbpmax) + i;
    cur_part.mass = ettp[id];
    prev_part.mass = ettpa[id];

    id = (*jdp-1) * (*nbpmax) + i;
    cur_part.diameter =  ettp[id];
    prev_part.diameter = ettpa[id];

    /* Coordinates of the particle */

    id = (*jxp-1) * (*nbpmax) + i;
    cur_part.coord[0] =  ettp[id];
    prev_part.coord[0] =  ettpa[id];

    id = (*jyp-1) * (*nbpmax) + i;
    cur_part.coord[1] =  ettp[id];
    prev_part.coord[1] =  ettpa[id];

    id = (*jzp-1) * (*nbpmax) + i;
    cur_part.coord[2] =  ettp[id];
    prev_part.coord[2] =  ettpa[id];

    /* Velocity of the particle */

    id = (*jup-1) * (*nbpmax) + i;
    cur_part.velocity[0] =  ettp[id];
    prev_part.velocity[0] =  ettpa[id];

    id = (*jvp-1) * (*nbpmax) + i;
    cur_part.velocity[1] =  ettp[id];
    prev_part.velocity[1] =  ettpa[id];

    id = (*jwp-1) * (*nbpmax) + i;
    cur_part.velocity[2] =  ettp[id];
    prev_part.velocity[2] =  ettpa[id];

    /* Velocity seen by the fluid */

    id = (*juf-1) * (*nbpmax) + i;
    cur_part.velocity_seen[0] =  ettp[id];
    prev_part.velocity_seen[0] =  ettpa[id];

    id = (*jvf-1) * (*nbpmax) + i;
    cur_part.velocity_seen[1] =  ettp[id];
    prev_part.velocity_seen[1] =  ettpa[id];

    id = (*jwf-1) * (*nbpmax) + i;
    cur_part.velocity_seen[2] =  ettp[id];
    prev_part.velocity_seen[2] =  ettpa[id];

    id = (*jtaux-1) * (*nbpmax) + i;
    cur_part.taup_aux =  ettp[id];

    /* Default visualization information set to off */
    cur_part.visualized = -1;

    /* Data needed if the deposition model is activated */
    if (*idepst > 0) {

      id = (*jryplu-1) * (*nbpmax) + i;
      cur_part.yplus = tepa[id];

      id = (*jrinpf-1) * (*nbpmax) + i;
      cur_part.interf = tepa[id];

      id = (*jdfac -1) * (*nbpmax) + i;
      cur_part.close_face_id = itepa[id] - 1;

      id = (*jimark -1) * (*nbpmax) + i;
      cur_part.marko_val = itepa[id];

    }
    else {

      cur_part.yplus = 10000;
      cur_part.close_face_id = -1;
      cur_part.marko_val = -1;

    }

    /* Update structures */

    set->particles[i] = cur_part;
    prev_set->particles[i] = prev_part;

  }


#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER PRTGET\n");
  _dump_particle_set(set);
  bft_printf("\n PREV PARTICLE SET AFTER PRTGET\n");
  _dump_particle_set(prev_set);
#endif
}

/*----------------------------------------------------------------------------
 * Put variables and parameters associated to each particles into FORTRAN
 * arrays.
 *
 * parameters:
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (prtput, PRTPUT)(const cs_int_t   *nbpmax,  /* n_particles max. */
                          cs_int_t         *nbpart,  /* number of current particles */
                          cs_real_t        *dnbpar,  /* particle total weight */
                          cs_int_t         *nbpout,  /* number of outgoing particles */
                          cs_real_t        *dnbpou,  /* outgoing particle total weight */
                          cs_int_t         *nbperr,  /* number of failed particles */
                          cs_real_t        *dnbper,  /* failed particles total weight */
                          cs_int_t         *nbpdep,  /* number of depositing particles */
                          cs_real_t        *dnbdep,  /* depositing particles total weight */
                          cs_int_t          liste[],
                          cs_int_t         *nbvis,
                          cs_real_t         ettp[],
                          cs_real_t         ettpa[],
                          cs_int_t          itepa[],
                          cs_real_t         tepa[],
                          cs_int_t          ibord[],
                          const cs_int_t   *jisor,
                          const cs_int_t   *jgnum,
                          const cs_int_t   *jrpoi,
                          const cs_int_t   *jrtsp,
                          const cs_int_t   *jdp,
                          const cs_int_t   *jmp,
                          const cs_int_t   *jxp,
                          const cs_int_t   *jyp,
                          const cs_int_t   *jzp,
                          const cs_int_t   *jup,
                          const cs_int_t   *jvp,
                          const cs_int_t   *jwp,
                          const cs_int_t   *juf,
                          const cs_int_t   *jvf,
                          const cs_int_t   *jwf,
                          const cs_int_t   *jtaux,
                          const cs_int_t   *jryplu,
                          const cs_int_t   *jrinpf,
                          const cs_int_t   *jdfac,
                          const cs_int_t   *jimark,
                          cs_int_t         *idepst)
{
  cs_lnum_t  i, j, k, id , nbp;

  cs_lagr_particle_set_t  *set = _particle_set;

  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;

  assert(*nbpmax == set->n_particles_max);

  j = set->first_used_id;

  nbp = 0;
  k = 0;

  for ( i = 0; (i < set-> n_particles )&( j != -1); i++) {

    nbp++;
    cs_lagr_particle_t  cur_part = set->particles[j];
    cs_lagr_particle_t  prev_part = prev_set->particles[j];

    if  (   cur_part.state == CS_LAGR_PART_TREATED
         || cur_part.state == CS_LAGR_PART_STICKED) {
      itepa[ (*jisor-1) * (*nbpmax) + i] = cur_part.cur_cell_num;
    }
    else if (cur_part.state == CS_LAGR_PART_OUT) {
      itepa[ (*jisor-1) * (*nbpmax) + i] = 0;
    }
    else {
      assert(cur_part.state = CS_LAGR_PART_ERR);

      itepa[ (*jisor-1) * (*nbpmax) + i] = 0;

    }
    itepa[ (*jgnum-1) * (*nbpmax) + i] = cur_part.global_num;


    // Data needed if the deposition model is activated

    if (*idepst > 0) {

      tepa[ (*jryplu - 1) * (*nbpmax) + i] = cur_part.yplus;

      tepa[ (*jrinpf - 1) * (*nbpmax) + i] = cur_part.interf;

      itepa[(*jdfac - 1) * (*nbpmax) + i]  = cur_part.close_face_id + 1;
      itepa[(*jimark - 1) * (*nbpmax) + i]  = cur_part.marko_val;

    }

    ibord[i] = cur_part.switch_order_1;
    // Not useful for prev_part

    id = (*jrpoi-1) * (*nbpmax) + i;
    tepa[id] = cur_part.stat_weight;
    tepa[id] = prev_part.stat_weight;

    id = (*jrtsp-1) * (*nbpmax) + i;
    tepa[id] = cur_part.residence_time;
    tepa[id] = prev_part.residence_time;

    id = (*jmp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.mass;
    ettpa[id] = prev_part.mass;

    id = (*jdp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.diameter;
    ettpa[id] = prev_part.diameter;

    // Coordinates of the particle

    id = (*jxp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.coord[0];
    ettpa[id] = prev_part.coord[0];

    id = (*jyp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.coord[1];
    ettpa[id] = prev_part.coord[1];

    id = (*jzp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.coord[2];
    ettpa[id] = prev_part.coord[2];

    // Velocity of the particle

    id = (*jup-1) * (*nbpmax) + i;
    ettp[id] = cur_part.velocity[0];
    ettpa[id] = prev_part.velocity[0];

    id = (*jvp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.velocity[1];
    ettpa[id] = prev_part.velocity[1];

    id = (*jwp-1) * (*nbpmax) + i;
    ettp[id] = cur_part.velocity[2];
    ettpa[id] = prev_part.velocity[2];

    // Velocity seen by the fluid

    id = (*juf-1) * (*nbpmax) + i;
    ettp[id] = cur_part.velocity_seen[0];
    ettpa[id] = prev_part.velocity_seen[0];

    id = (*jvf-1) * (*nbpmax) + i;
    ettp[id] = cur_part.velocity_seen[1];
    ettpa[id] = prev_part.velocity_seen[1];

    id = (*jwf-1) * (*nbpmax) + i;
    ettp[id] = cur_part.velocity_seen[2];
    ettpa[id] = prev_part.velocity_seen[2];

    id = (*jtaux-1) * (*nbpmax) + i;
    ettp[id] = cur_part.taup_aux;

    /* Next particle id to treat */


    j = cur_part.next_id;

    assert(cur_part.next_id == prev_part.next_id);

  } /* End of loop on particles */

  /* New number of particles */
  *nbpart= set->n_particles;

  /* New weight */
  *dnbpar= set->weight;

  /* Number of exiting particles */
  *nbpout = set->n_part_out;

  /* weight of exiting particles */
  *dnbpou = set->weight_out;

  /* Number of depositing particles */
  *nbpdep = set->n_part_dep;

  /* weight of depositing particles */
  *dnbdep = set->weight_dep;

  /* Number of failed particles */
  *nbperr = set->n_failed_part;

  /* weight of failed particles */
  *dnbper = set->weight_failed;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER PRTPUT\n");
  _dump_particle_set(set);
  bft_printf("\n PREV PARTICLE SET AFTER PRTPUT\n");
  _dump_particle_set(prev_set);
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
                          const cs_int_t     iusmoy[],
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
 *  p_n_particles     <->  pointer to the number of particles
 *  scheme_order    -->  current order of the scheme used for Lagrangian
 *----------------------------------------------------------------------------*/

void
CS_PROCF (dplprt, DPLPRT)(cs_lnum_t        *p_n_particles,
                          cs_real_t        *p_parts_weight,
                          cs_lnum_t        *p_scheme_order,
                          cs_real_t         boundary_stat[],
                          const cs_lnum_t  *iensi3,
                          const cs_lnum_t  *nvisbr,
                          const cs_lnum_t  *inbr,
                          const cs_lnum_t  *inbrbd,
                          const cs_lnum_t  *iflm,
                          const cs_lnum_t  *iflmbd,
                          const cs_lnum_t  *iang,
                          const cs_lnum_t  *iangbd,
                          const cs_lnum_t  *ivit,
                          const cs_lnum_t  *ivitbd,
                          const cs_lnum_t  *nusbor,
                          cs_lnum_t         iusb[],
                          cs_real_t         visc_length[],
                          cs_real_t         dlgeo[],
                          cs_real_t         rtp[],
                          const cs_lnum_t  *iu,
                          const cs_lnum_t  *iv,
                          const cs_lnum_t  *iw,
                          cs_lnum_t        *idepst)
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
  cs_lagr_particle_set_t  *set = _particle_set;
  cs_lagr_particle_set_t  *prev_set = _prev_particle_set;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  const cs_lnum_t  failsafe_mode = 0; /* If 1 : stop as soon as an error is
                                         detected */

  assert(builder != NULL);
  assert(set != NULL && prev_set != NULL);

  /* Main loop on  particles : global propagation */

  while ( _continue_displacement() ) {

    n_delete_particles = 0;
    n_failed_particles = 0;

    r_weight = 0.0;
    tot_weight = 0.0;
    failed_particle_weight = 0.0;

    assert(set->first_free_id != -1);

    /* Local propagation */

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {
      cs_lagr_particle_t*  cur_part = &set->particles[j];
      cs_lagr_particle_t*  prev_part = &prev_set->particles[j];

      if (cur_part->state == CS_LAGR_PART_TO_SYNC)
        cur_part->state = _local_propagation(prev_part,
                                             cur_part,
                                             n_displacement_steps,
                                             scheme_order,
                                             failsafe_mode,
                                             boundary_stat,
                                             &n_failed_particles,
                                             &failed_particle_weight,
                                             iensi3,
                                             nvisbr,
                                             inbr,
                                             inbrbd,
                                             iflm,
                                             iflmbd,
                                             iang,
                                             iangbd,
                                             ivit,
                                             ivitbd,
                                             nusbor,
                                             iusb,
                                             visc_length,
                                             dlgeo,
                                             rtp, iu, iv ,iw,
                                             idepst);

      prev_part->state = cur_part->state;
      j = cur_part->next_id;

      assert(cur_part->next_id == prev_part->next_id);

    } /* End of loop on particles */

    /* Update of the particle set structure. Delete particles. */

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

      cs_lagr_particle_t  cur_part = set->particles[j];
      cs_lagr_particle_t  prev_part = prev_set->particles[j];

      // FIXME: assert(cur_part.state == prev_part.state);

      if (   cur_part.state == CS_LAGR_PART_TO_DELETE
          || cur_part.state == CS_LAGR_PART_OUT
          || cur_part.state == CS_LAGR_PART_ERR) {

        _remove_particle(set, cur_part, j);
        _remove_particle(prev_set, prev_part, j);
        n_delete_particles++;
        r_weight += cur_part.stat_weight;

      }
      else {

        tot_weight += cur_part.stat_weight;

      }

      /*
        cur_part.next_id modified inside _remove_particle() has no effect
        to the next line. As cur_part is a parameter of _remove_particle(),
        it's only a copy which goes through the function.
      */

      j = cur_part.next_id;

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

    if (mesh->n_init_perio > 0 || cs_glob_n_ranks > 1) {

      /* Synchronisation of a selection of particles for parallelism and
         periodicity. Number of particles on the local rank may change. */

      _lagr_halo_sync();
    }

    n_displacement_steps++;

  } /* End of while (global displacement) */

  /* Deposition sub-model additional loop */

  if (*idepst > 0) {

    for (i = 0, j = set->first_used_id; i < set->n_particles; i++) {

      cs_lagr_particle_t*  cur_part = &set->particles[j];

      _test_wall_cell(cur_part,visc_length,dlgeo);

      if (cur_part->yplus < 100.e0) {

        /* Todo : specific treatment */
      }

      j = cur_part->next_id;

    }
  }

  /* Returns pointers */

  *p_n_particles =  set->n_particles;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Delete cs_lagr_particle_set_t structure and delete other useful buffers.
 *----------------------------------------------------------------------------*/

void
cs_lagr_destroy(void)
{
  /* Destroy particle sets */

  _prev_particle_set = _destroy_particle_set(_prev_particle_set);
  _particle_set = _destroy_particle_set(_particle_set);

  /* Destroy builder */
  _particle_track_builder = _destroy_track_builder(_particle_track_builder);

  /* Destroy boundary condition structure */

  _lagr_bdy_conditions = _destroy_bdy_cond_struct(_lagr_bdy_conditions);

  /* Delete MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)  _delete_particle_datatypes();
#endif
}

/*----------------------------------------------------------------------------*/

/* Delete local macro definitions */

END_C_DECLS
