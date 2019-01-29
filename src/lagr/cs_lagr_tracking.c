/*============================================================================
 * Methods for particle localization
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
#include "cs_boundary_zone.h"
#include "cs_physical_constants.h"
#include "cs_geom.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_random.h"
#include "cs_rotation.h"
#include "cs_search.h"
#include "cs_timer_stats.h"
#include "cs_turbomachinery.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr.h"
#include "cs_lagr_deposition_model.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_post.h"
#include "cs_lagr_clogging.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_geom.h"

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

/* State where a particle can be.
   (order is chosen so as to make tests simpler;
   inside domain first, outside after) */

typedef enum {
  CS_LAGR_PART_TO_SYNC,
  CS_LAGR_PART_TREATED,
  CS_LAGR_PART_STUCK,
  CS_LAGR_PART_MERGED,
  CS_LAGR_PART_OUT,
  CS_LAGR_PART_ERR
} cs_lagr_tracking_state_t;

enum {
  CS_LAGR_PART_MOVE_OFF = 0,
  CS_LAGR_PART_MOVE_ON  = 1
};

/* Tracking error types */

typedef enum {

  CS_LAGR_TRACKING_OK,
  CS_LAGR_TRACKING_ERR_MAX_LOOPS,
  CS_LAGR_TRACKING_ERR_LOST_PIC

} cs_lagr_tracking_error_t;

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

  cs_real_t  start_coords[3];       /* starting coordinates for
                                       next displacement */

  cs_lnum_t  last_face_num;         /* last face number encountered */

  cs_lagr_tracking_state_t  state;  /* current state */

} cs_lagr_tracking_info_t;

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
  size_t      extents;        /* Extents for particle set */

  cs_lnum_t  *send_count;     /* number of particles to send to
                                 each communicating rank */
  cs_lnum_t  *recv_count;     /* number of particles to receive from
                                 each communicating rank */

  cs_lnum_t  *send_shift;
  cs_lnum_t  *recv_shift;

  unsigned char  *send_buf;

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

  cs_lagr_halo_t      *halo;   /* Lagrangian halo structure */

  cs_interface_set_t  *face_ifs;

} cs_lagr_track_builder_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Global variable for the current subroutines */

static  cs_lagr_track_builder_t  *_particle_track_builder = NULL;

static  int            _max_propagation_loops = 100;

/* MPI datatype associated to each particle "structure" */

#if defined(HAVE_MPI)
static  MPI_Datatype  _cs_mpi_particle_type;
#endif

/* Global Lagragian module parameters and associated pointer
   Should move to cs_lagr.c */


/*=============================================================================
 * Private function definitions
 *============================================================================*/

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

  BFT_MALLOC(lagr_halo->send_buf,
             lagr_halo->send_buf_size * extents,
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

    BFT_FREE(*halo);
  }

}

/*----------------------------------------------------------------------------
 * Resize a halo's buffers.
 *
 * parameters:
 *   lag_halo         <->  pointer to a cs_lagr_halo_t structure
 *   n_send_particles <-- number of particles to send
 *----------------------------------------------------------------------------*/

static void
_resize_lagr_halo(cs_lagr_halo_t  *lag_halo,
                  cs_lnum_t        n_send_particles)
{
  cs_lnum_t n_halo = lag_halo->send_buf_size;

  size_t tot_extents = lag_halo->extents;

  /* If increase is required */

  if (n_halo < n_send_particles) {
    if (n_halo < CS_LAGR_MIN_COMM_BUF_SIZE)
      n_halo = CS_LAGR_MIN_COMM_BUF_SIZE;
    while (n_halo < n_send_particles)
      n_halo *= 2;
    lag_halo->send_buf_size = n_halo;
    BFT_REALLOC(lag_halo->send_buf, n_halo*tot_extents, unsigned char);
  }

  /* If decrease is allowed, do it progressively, and with a wide
     margin, so as to avoid re-increasing later if possible */

  else if (n_halo > n_send_particles*16) {
    n_halo /= 8;
    lag_halo->send_buf_size = n_halo;
    BFT_REALLOC(lag_halo->send_buf, n_halo*tot_extents, unsigned char);
  }

  /* Otherwise, keep current size */
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

      cs_lnum_t cell_id = mesh->i_face_cells[i][j];

      if (cell_id < mesh->n_cells) {

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

  builder->face_ifs = NULL;

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
 * Manage detected errors
 *
 * parameters:
 *   failsafe_mode            <-- indicate if failsafe mode is used
 *   particle                 <-- pointer to particle data
 *   attr_map                 <-- pointer to attribute map
 *   error_type               <-- error code
 *   msg                      <-- error message
 *----------------------------------------------------------------------------*/

static void
_manage_error(cs_lnum_t                       failsafe_mode,
              void                           *particle,
              const cs_lagr_attribute_map_t  *attr_map,
              cs_lagr_tracking_error_t        error_type)
{
  cs_real_t *prev_part_coord
    = cs_lagr_particle_attr_n(particle, attr_map, 1, CS_LAGR_COORDS);
  cs_real_t *part_coord
    = cs_lagr_particle_attr(particle, attr_map, CS_LAGR_COORDS);

  const cs_real_t  *prev_location
    = ((const cs_lagr_tracking_info_t *)particle)->start_coords;

  cs_real_t d0 = cs_math_3_distance(part_coord, prev_part_coord);
  cs_real_t d1 = cs_math_3_distance(part_coord, prev_location);

  cs_lagr_particle_set_real(particle, attr_map, CS_LAGR_TR_TRUNCATE, d1/d0);

  if (error_type == CS_LAGR_TRACKING_ERR_LOST_PIC)
    cs_lagr_particle_set_real(particle, attr_map, CS_LAGR_TR_TRUNCATE, 2.0);

  if (failsafe_mode == 1) {
    switch (error_type) {
    case CS_LAGR_TRACKING_ERR_MAX_LOOPS:
      bft_error(__FILE__, __LINE__, 0,
                _("Max number of loops reached in particle displacement."));
      break;
    case CS_LAGR_TRACKING_ERR_LOST_PIC:
      bft_error(__FILE__, __LINE__, 0,
                _("Particle lost in local_propagation: it has been removed"));
      break;
    default:
      break;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute the contribution of a particle to the boundary mass flux.
 *
 * Used for the deposition model.
 *
 * parameters:
 *   particle_set     <-- pointer to particle set structure
 *   particle_id      <-- pointer to particle id
 *   sign             <-- -1 to remove contribution, 1 to add it
 *   b_face_surf      <-- boundary face surface
 *   part_b_mass_flux <-> particle mass flux array
 *----------------------------------------------------------------------------*/

static void
_b_mass_contribution(const cs_lagr_particle_set_t   *particles,
                     cs_lnum_t                       particle_id,
                     const cs_real_t                 sign,
                     const cs_real_t                 b_face_surf[],
                     cs_real_t                       part_b_mass_flux[])
{
  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  const unsigned char *particle
    = particles->p_buffer + p_am->extents * particle_id;

  cs_lnum_t depo_flag
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG);

  if (   depo_flag == CS_LAGR_PART_ROLLING
      || depo_flag == CS_LAGR_PART_DEPOSITED) {

    cs_lnum_t neighbor_face_id
      = cs_lagr_particle_get_lnum(particle, p_am,
                                  CS_LAGR_NEIGHBOR_FACE_ID);

    assert(neighbor_face_id > -1);

    cs_real_t cur_stat_weight
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
    cs_real_t cur_mass
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
    cs_real_t face_area = b_face_surf[neighbor_face_id];

    part_b_mass_flux[neighbor_face_id]
      += sign * cur_stat_weight * cur_mass / face_area;

  }
}

/*----------------------------------------------------------------------------
 * Handle particles moving to internal deposition face
 *
 * parameters:
 *   particles       <-- pointer to particle set
 *   particle        <-> particle data for current particle
 *   face_id         <-- index of the treated face
 *   t_intersect     <-- used to compute the intersection of the trajectory and
 *                       the face
 *   p_move_particle <-- particle moves?
 *
 * returns:
 *   particle state
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_internal_treatment(cs_lagr_particle_set_t    *particles,
                    void                      *particle,
                    cs_lnum_t                  face_id,
                    double                     t_intersect,
                    int                       *p_move_particle)
{
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const double bc_epsilon = 1.e-2;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  cs_real_t  disp[3], face_normal[3], intersect_pt[3];

  int move_particle = *p_move_particle;

  cs_lagr_tracking_state_t  particle_state = CS_LAGR_PART_TO_SYNC;

  cs_lagr_internal_condition_t *internal_conditions
    = cs_glob_lagr_internal_conditions;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  cs_real_t  *particle_coord
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_velocity
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);
  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

  cs_real_t particle_stat_weight
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
  cs_real_t particle_mass
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

  assert(internal_conditions != NULL);

  for (int k = 0; k < 3; k++)
    face_normal[k] = fvq->i_face_normal[3*face_id+k];

  cs_real_t face_area  = fvq->i_face_surf[face_id];

  cs_real_t  face_norm[3] = {face_normal[0]/face_area,
                             face_normal[1]/face_area,
                             face_normal[2]/face_area};

  cs_lnum_t  cell_id
    = cs_lagr_particle_get_cell_id(particle, p_am);
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  for (int k = 0; k < 3; k++)
    disp[k] = particle_coord[k] - p_info->start_coords[k];

  assert(! (fabs(disp[0]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[1]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[2]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15));

  for (int k = 0; k < 3; k++)
    intersect_pt[k] = disp[k]*t_intersect + p_info->start_coords[k];

  if (internal_conditions->i_face_zone_id[face_id] == CS_LAGR_OUTLET
   || internal_conditions->i_face_zone_id[face_id] == CS_LAGR_INLET ) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_OUT;

    for (int k = 0; k < 3; k++)
      particle_coord[k] = intersect_pt[k];
  }
  else if (internal_conditions->i_face_zone_id[face_id] == CS_LAGR_DEPO_DLVO) {

    cs_real_t particle_diameter
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

    cs_real_t uxn = particle_velocity[0] * face_norm[0];
    cs_real_t vyn = particle_velocity[1] * face_norm[1];
    cs_real_t wzn = particle_velocity[2] * face_norm[2];

    cs_real_t energ = 0.5 * particle_mass * (uxn+vyn+wzn) * (uxn+vyn+wzn);

    /* No energy barrier yet for internal deposition */
    cs_real_t  energt = 0.;

     /* Deposition criterion: E_kin > E_barr */
    if (energ >= energt * 0.5 * particle_diameter) {

      cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);
      cs_real_3_t vect_cen;
      for (int k = 0; k < 3; k++)
        vect_cen[k] = (cell_cen[k] - intersect_pt[k]);

      for (int k = 0; k < 3; k++) {
        particle_velocity[k] = 0.0;
      }
      /* Force the particle on the intersection but in the original cell */
      for (int k = 0; k < 3; k++) {
        particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
        particle_velocity_seen[k] = 0.0;
      }
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                cell_id +1);
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID,
                                face_id);
      // The particle is not treated yet: the motion is now imposed
      move_particle = CS_LAGR_PART_MOVE_OFF;
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_IMPOSED_MOTION);

      /* Specific treatment in case of particle resuspension modeling */

      particle_state = CS_LAGR_PART_TREATED;

      particles->n_part_dep += 1;
      particles->weight_dep += particle_stat_weight;

    }
  }
  /* FIXME: JBORD* (user-defined boundary condition) not yet implemented
     nor defined by a macro */
  else if (internal_conditions->i_face_zone_id[face_id] != -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Internal condition %d not recognized.\n"),
              internal_conditions->i_face_zone_id[face_id]);

  /* Return pointer */

  *p_move_particle = move_particle;

  /* TODO internal statts ?... */

  return particle_state;
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
                    cs_real_t                 *face_norm,
                    double                     t_intersect,
                    int                        boundary_zone,
                    int                       *p_move_particle)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const double pi = cs_math_pi;
  const double bc_epsilon = 1.e-2;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_real_t  tmp;
  cs_real_t  disp[3], face_normal[3], intersect_pt[3];

  cs_real_t  compo_vel[3] = {0.0, 0.0, 0.0};
  cs_real_t  norm_vel = 0.0;

  int move_particle = *p_move_particle;

  cs_lnum_t  face_id = face_num - 1;
  cs_lagr_tracking_state_t  particle_state = CS_LAGR_PART_TO_SYNC;

  cs_lagr_zone_data_t  *bdy_conditions = cs_lagr_get_boundary_conditions();

  cs_real_t  energt = 0.;
  cs_lnum_t  contact_number = 0;
  cs_real_t  *surface_coverage = NULL;
  cs_real_t* deposit_height_mean = NULL;
  cs_real_t* deposit_height_var = NULL;
  cs_real_t* deposit_diameter_sum = NULL;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  cs_real_t  *particle_coord
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_velocity
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);
  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

  cs_real_t particle_stat_weight
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
  cs_real_t particle_mass
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const char b_type = cs_glob_lagr_boundary_conditions->elt_type[face_id];

  assert(bdy_conditions != NULL);

  for (int k = 0; k < 3; k++)
    disp[k] = particle_coord[k] - p_info->start_coords[k];

  for (int k = 0; k < 3; k++) {
    face_normal[k] = fvq->b_face_normal[3*face_id+k];
  }

  cs_real_t face_area  = fvq->b_face_surf[face_id];


  cs_lnum_t  cell_id
    = cs_lagr_particle_get_cell_id(particle, p_am);
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  assert(! (fabs(disp[0]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[1]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[2]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15));

  /* Save particle impacting velocity */
  if (   cs_glob_lagr_boundary_interactions->iangbd > 0
      || cs_glob_lagr_boundary_interactions->ivitbd > 0) {
    norm_vel = cs_math_3_norm(particle_velocity);
    for (int k = 0; k < 3; k++)
      compo_vel[k] = particle_velocity[k];
  }

  for (int k = 0; k < 3; k++)
    intersect_pt[k] = disp[k]*t_intersect + p_info->start_coords[k];

  if (   b_type == CS_LAGR_OUTLET
      || b_type == CS_LAGR_INLET
      || b_type == CS_LAGR_DEPO1) {

    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_OUT;

    if (b_type == CS_LAGR_DEPO1) {
      particles->n_part_dep += 1;
      particles->weight_dep += particle_stat_weight;
      if (cs_glob_lagr_model->deposition == 1)
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                  CS_LAGR_PART_DEPOSITED);
    }

    /* FIXME: For post-processing by trajectory purpose */

    for (int k = 0; k < 3; k++)
      particle_coord[k] = intersect_pt[k];
  }

  else if (b_type == CS_LAGR_DEPO2) {

    move_particle = CS_LAGR_PART_MOVE_OFF;

    cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);
    cs_real_3_t vect_cen;
    for (int k = 0; k < 3; k++) {
      vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
      particle_velocity[k] = 0.0;
      particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
    }

    particles->n_part_dep += 1;
    particles->weight_dep += particle_stat_weight;

    /* Specific treatment in case of particle resuspension modeling */

    cs_lnum_t *cell_num = cs_lagr_particle_attr(particle,
                                                p_am,
                                                CS_LAGR_CELL_NUM);

    if (cs_glob_lagr_model->deposition == 1)
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_DEPOSITED);

    if (cs_glob_lagr_model->resuspension == 0) {

      *cell_num = - *cell_num;

      for (int k = 0; k < 3; k++)
        particle_velocity_seen[k] = 0.0;

      particle_state = CS_LAGR_PART_STUCK;

    } else {

      *cell_num = cs_glob_mesh->b_face_cells[face_id] + 1;

      particle_state = CS_LAGR_PART_TREATED;

    }

  }

  else if (b_type == CS_LAGR_DEPO_DLVO) {

    cs_real_t particle_diameter
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id] + 1);

    cs_real_t uxn = particle_velocity[0] * face_norm[0];
    cs_real_t vyn = particle_velocity[1] * face_norm[1];
    cs_real_t wzn = particle_velocity[2] * face_norm[2];

    cs_real_t energ = 0.5 * particle_mass * (uxn+vyn+wzn) * (uxn+vyn+wzn);
    cs_real_t min_porosity;
    cs_real_t limit;

    if (cs_glob_lagr_model->clogging) {

      /* If the clogging modeling is activated,                 */
      /* computation of the number of particles in contact with */
      /* the depositing particle                                */

      surface_coverage = &bound_stat[cs_glob_lagr_boundary_interactions->iscovc * n_b_faces + face_id];
      deposit_height_mean = &bound_stat[cs_glob_lagr_boundary_interactions->ihdepm * n_b_faces + face_id];
      deposit_height_var = &bound_stat[cs_glob_lagr_boundary_interactions->ihdepv * n_b_faces + face_id];

      deposit_diameter_sum = &bound_stat[cs_glob_lagr_boundary_interactions->ihsum * n_b_faces + face_id];

      contact_number = cs_lagr_clogging_barrier(particle,
                                                p_am,
                                                face_id,
                                                &energt,
                                                surface_coverage,
                                                &limit,
                                                &min_porosity);

      if (contact_number == 0 && cs_glob_lagr_model->roughness > 0) {
        cs_lagr_roughness_barrier(particle,
                                  p_am,
                                  face_id,
                                  &energt);
      }
    }
    else {

      if (cs_glob_lagr_model->roughness > 0)
        cs_lagr_roughness_barrier(particle,
                                  p_am,
                                  face_id,
                                  &energt);

      else if (cs_glob_lagr_model->roughness == 0) {
        cs_lagr_barrier(particle,
                        p_am,
                        face_id,
                        &energt);
      }

    }

     /* Deposition criterion: E_kin > E_barr */
    if (energ > energt * 0.5 * particle_diameter) {

      cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);
      cs_real_3_t vect_cen;
      for (int k = 0; k < 3; k++)
        vect_cen[k] = (cell_cen[k] - intersect_pt[k]);

      /* The particle deposits*/
      if (!cs_glob_lagr_model->clogging && !cs_glob_lagr_model->resuspension) {
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                  CS_LAGR_PART_DEPOSITED);

        move_particle = CS_LAGR_PART_MOVE_OFF;

        /* Set negative value for current cell number */
        cs_lagr_particle_set_lnum
          (particle, p_am, CS_LAGR_CELL_NUM,
           - cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM));

        particles->n_part_dep += 1;
        particles->weight_dep += particle_stat_weight;

        particle_state = CS_LAGR_PART_STUCK;
      }

      if (!cs_glob_lagr_model->clogging && cs_glob_lagr_model->resuspension > 0) {

        move_particle = CS_LAGR_PART_MOVE_OFF;
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                  CS_LAGR_PART_DEPOSITED);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                  cs_glob_mesh->b_face_cells[face_id] + 1);

        /* The particle is replaced towards the cell center
         * (and not using the normal vector face_norm)
         * to avoid problems where the replacement is out of the cell.
         * This is valid only for star-shaped cells !!! */
        for (int k = 0; k < 3; k++) {
          particle_velocity[k] = 0.0;
          particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
        }
        particles->n_part_dep += 1;
        particles->weight_dep += particle_stat_weight;
        particle_state = CS_LAGR_PART_TREATED;

      }

      if (cs_glob_lagr_model->clogging) {

        bound_stat[cs_glob_lagr_boundary_interactions->inclgt
                   * n_b_faces + face_id] += particle_stat_weight;
        *deposit_diameter_sum += particle_diameter;

        cs_real_t particle_height
          = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_HEIGHT);

        cs_real_t depositing_radius = particle_diameter * 0.5;

        if (contact_number == 0) {

          /* The surface coverage increases if the particle
             has deposited on a naked surface */

          *surface_coverage += (pi * pow(depositing_radius,2))
                               *  particle_stat_weight / face_area;

          *deposit_height_mean +=   particle_height* pi * pow(depositing_radius,2)
                                 *  particle_stat_weight / face_area;
          *deposit_height_var +=   pow(particle_height * pi
                                 * particle_stat_weight / face_area, 2)
                                 * pow(depositing_radius,4);

          bound_stat[cs_glob_lagr_boundary_interactions->inclg
                     * n_b_faces + face_id] += particle_stat_weight;

          /* The particle is replaced towards the cell center
           * (and not using the normal vector face_norm)
           * to avoid problems where the replacement is out of the cell.
           * This is valid only for star-shaped cells !!! */
          for (int k = 0; k < 3; k++) {
            particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
            particle_velocity[k] = 0.0;
            particle_velocity_seen[k] = 0.0;
          }

          move_particle = CS_LAGR_PART_MOVE_OFF;
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                    CS_LAGR_PART_DEPOSITED);
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                    cs_glob_mesh->b_face_cells[face_id] + 1);
          cs_lagr_particle_set_lnum(particle, p_am,CS_LAGR_NEIGHBOR_FACE_ID ,
                                    face_id);

          particles->n_part_dep += 1;
          particles->weight_dep += particle_stat_weight;
          particle_state = CS_LAGR_PART_TREATED;
        }
        else {

          cs_lnum_t i;
          cs_real_t random = -1;
          cs_real_t scov_cdf;
          void *cur_part = NULL;

          /* We choose randomly a deposited particle to interact
             with the depositing one (according to the relative surface coverage
             of each class of particles) */

          cs_random_uniform(1, &random);
          cs_real_t scov_rand =  (random * (*surface_coverage));

          scov_cdf = 0.;

          for (i = 0; i < particles->n_particles; i++) {

            /* FIXME wrong test (don't know what the intent was,
               but compiler warns result may not be as intended)

               was:

               if (CS_LAGR_PART_TREATED
                   <= _get_tracking_info(particles, i)->state
                   <= CS_LAGR_PART_TO_SYNC)
                 continue;

               evaluates to expression below:
            */

            if (  CS_LAGR_PART_TREATED
                > _get_tracking_info(particles, i)->state)
              continue;

            cur_part = (void *)(particles->p_buffer + p_am->extents * i);

            cs_lnum_t cur_part_depo
              = cs_lagr_particle_get_lnum(cur_part, p_am,
                                          CS_LAGR_DEPOSITION_FLAG);

            cs_lnum_t cur_part_close_face_id
              = cs_lagr_particle_get_lnum(cur_part, p_am,
                                          CS_LAGR_NEIGHBOR_FACE_ID);

            cs_real_t cur_part_stat_weight
              = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_STAT_WEIGHT);

            cs_real_t cur_part_diameter
              = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_DIAMETER);

            if (   (cur_part_depo == CS_LAGR_PART_DEPOSITED)
                && (cur_part_close_face_id == face_id)) {
              scov_cdf +=   (pi * pow(cur_part_diameter,2) / 4.)
                          *  cur_part_stat_weight / face_area;
              if (scov_cdf >= scov_rand)
                break;
            }
          }

          cs_lnum_t cur_part_close_face_id
            = cs_lagr_particle_get_lnum(cur_part, p_am,
                                        CS_LAGR_NEIGHBOR_FACE_ID);

          cs_lnum_t particle_close_face_id
            = cs_lagr_particle_get_lnum(particle, p_am,
                                        CS_LAGR_NEIGHBOR_FACE_ID);

          if (cur_part_close_face_id != face_id) {
            bft_error(__FILE__, __LINE__, 0,
                      _(" Error in %s: in the face number %d \n"
                        "no deposited particle found to form a cluster \n"
                        "using the surface coverage %e (scov_cdf %e) \n"
                        "The particle used thus belongs to another face (%d) \n"),
                      __func__,
                      particle_close_face_id, *surface_coverage,
                      scov_cdf, cur_part_close_face_id);
          }

          /* The depositing particle is merged with the existing one */
          /* Statistical weight obtained conserving weight*mass*/
          cs_real_t cur_part_stat_weight
            = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_STAT_WEIGHT);

          cs_real_t cur_part_mass
            = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_MASS);

          cs_real_t cur_part_diameter
            = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_DIAMETER);

          cs_real_t cur_part_height
            = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_HEIGHT);

          cs_real_t cur_part_cluster_nb_part
            = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_CLUSTER_NB_PART);

          cs_real_t particle_cluster_nb_part
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART);

          *deposit_height_mean -=   cur_part_height*pi*pow(cur_part_diameter, 2)
                                  * cur_part_stat_weight / (4.0*face_area);
          *deposit_height_var -=   pow(cur_part_height*pi
                                  * cur_part_stat_weight / (4.0*face_area),2)
                                  * pow(cur_part_diameter, 4);

          if (*surface_coverage >= limit) {
            cs_lagr_particle_set_real(cur_part, p_am, CS_LAGR_HEIGHT,
                                      cur_part_height
                                      +  (  pow(particle_diameter, 3)
                                          / cs_math_sq(cur_part_diameter)
                                          * particle_stat_weight
                                          / cur_part_stat_weight));
          }
          else {
            *surface_coverage -= (pi * pow(cur_part_diameter,2)/4.)
              * cur_part_stat_weight / face_area;

            cs_lagr_particle_set_real(cur_part, p_am, CS_LAGR_DIAMETER,
                                      pow (  pow(cur_part_diameter,3)
                                           + pow(particle_diameter,3)
                                           * particle_stat_weight
                                             / cur_part_stat_weight , 1./3.));

            cur_part_diameter = cs_lagr_particle_get_real(cur_part, p_am,
                                                          CS_LAGR_DIAMETER);

            *surface_coverage +=   (pi * pow(cur_part_diameter,2)/4.)
                                 * cur_part_stat_weight / face_area;

            cs_lagr_particle_set_real(cur_part, p_am, CS_LAGR_HEIGHT,
                                      cur_part_diameter);
          }

          cs_lagr_particle_set_real(cur_part, p_am, CS_LAGR_MASS,
                                    cur_part_mass + particle_mass
                                    * particle_stat_weight / cur_part_stat_weight);
          cs_lagr_particle_set_real(cur_part, p_am, CS_LAGR_CLUSTER_NB_PART,
                                    cur_part_cluster_nb_part+particle_cluster_nb_part
                                    * particle_stat_weight / cur_part_stat_weight);

          move_particle = CS_LAGR_PART_MOVE_OFF;
          particle_state = CS_LAGR_PART_OUT;
          particles->n_part_dep += 1;
          particles->weight_dep += particle_stat_weight;

          cur_part_height   = cs_lagr_particle_get_real(cur_part, p_am,
                                                        CS_LAGR_HEIGHT);

          *deposit_height_mean +=   cur_part_height*pi*pow(cur_part_diameter,2)
                                  * cur_part_stat_weight / (4.0*face_area);
          *deposit_height_var +=    pow(cur_part_height*pi
                                  * cur_part_stat_weight / (4.0*face_area),2)
                                  * pow(cur_part_diameter,4);
        }

      }

    }
    else {

      /*The particle does not deposit:
        It 'rebounds' on the energy barrier*/

      move_particle = CS_LAGR_PART_MOVE_ON;
      particle_state = CS_LAGR_PART_TO_SYNC;
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                cs_glob_mesh->b_face_cells[face_id] + 1);
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_IN_FLOW);

      cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);
      cs_real_3_t vect_cen;
      for (int k = 0; k < 3; k++) {
        vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
        p_info->start_coords[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
      }

      /* Modify the ending point. */

      for (int k = 0; k < 3; k++)
        disp[k] = particle_coord[k] - intersect_pt[k];

      tmp = 2. * cs_math_3_dot_product(disp, face_norm);

      for (int k = 0; k < 3; k++)
        particle_coord[k] -= tmp * face_norm[k];

      /* Modify particle velocity and velocity seen */

      tmp = 2. * cs_math_3_dot_product(particle_velocity, face_norm);

      for (int k = 0; k < 3; k++)
        particle_velocity[k] -= tmp * face_norm[k];

      tmp = 2. * cs_math_3_dot_product(particle_velocity_seen, face_norm);

      for (int k = 0; k < 3; k++)
        particle_velocity_seen[k] -= tmp * face_norm[k];
    }
  }

  else if (b_type == CS_LAGR_REBOUND || b_type == CS_LAGR_SYM) {

    move_particle = CS_LAGR_PART_MOVE_ON;
    particle_state = CS_LAGR_PART_TO_SYNC;
    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                              cs_glob_mesh->b_face_cells[face_id] + 1);

    cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);
    cs_real_3_t vect_cen;
    for (int k = 0; k < 3; k++) {
      vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
      p_info->start_coords[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
    }

    /* Modify the ending point. */

    for (int k = 0; k < 3; k++)
      disp[k] = particle_coord[k] - intersect_pt[k];

    tmp = 2. * cs_math_3_dot_product(disp, face_norm);

    for (int k = 0; k < 3; k++)
      particle_coord[k] -= tmp * face_norm[k];

    /* Modify particle velocity and velocity seen */

    tmp = 2. * cs_math_3_dot_product(particle_velocity, face_norm);

    for (int k = 0; k < 3; k++)
      particle_velocity[k] -= tmp * face_norm[k];

    tmp = 2. * cs_math_3_dot_product(particle_velocity_seen, face_norm);

    for (int k = 0; k < 3; k++)
      particle_velocity_seen[k] -= tmp * face_norm[k];

  }

  else if (b_type == CS_LAGR_FOULING) {

    /* Fouling of the particle, if its properties make it possible and
       with respect to a probability
       HERE if  Tp     > TPENC
           if viscp <= VISREF ==> Probability of fouling equal to 1
           if viscp  > VISREF ==> Probability equal to TRAP = 1-VISREF/viscp
                              ==> Fouling if VNORL is between TRAP et 1. */

    cs_real_t  random = -1, viscp = -1, trap = -1;

    /* Selection of the fouling coefficient*/

    const cs_lnum_t p_coal_id
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_COAL_ID);
    const cs_lnum_t n_layers = p_am->count[0][CS_LAGR_TEMPERATURE];
    const cs_real_t *particle_temp
      = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_TEMPERATURE);

    cs_real_t  temp_ext_part = particle_temp[n_layers - 1];
    cs_real_t  tprenc_icoal
      = cs_glob_lagr_encrustation->tprenc[p_coal_id];
    cs_real_t  visref_icoal
      = cs_glob_lagr_encrustation->visref[p_coal_id];
    cs_real_t  enc1_icoal
      = cs_glob_lagr_encrustation->enc1[p_coal_id];
    cs_real_t  enc2_icoal
      = cs_glob_lagr_encrustation->enc2[p_coal_id];

    if (temp_ext_part > tprenc_icoal+cs_physical_constants_celsius_to_kelvin) {

      /* Coal viscosity*/
      tmp = (  (1.0e7*enc1_icoal)
             / ((temp_ext_part-150.e0-cs_physical_constants_celsius_to_kelvin)
               *(temp_ext_part-150.e0-cs_physical_constants_celsius_to_kelvin)))
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
        cs_random_uniform(1, &random);
        trap = 1.e0- (visref_icoal / viscp);
      }

      if (   (viscp <= visref_icoal)
          || (viscp >= visref_icoal  &&  random >= trap)) {

        move_particle = CS_LAGR_PART_MOVE_OFF;
        particle_state = CS_LAGR_PART_OUT;

        /* Recording for log/lagrangian.log */
        particles->n_part_fou += 1;
        particles->weight_fou += particle_stat_weight;

        /* Recording for statistics */
        if (cs_glob_lagr_boundary_interactions->iencnbbd > 0) {
          bound_stat[  cs_glob_lagr_boundary_interactions->iencnb
                     * n_b_faces + face_id]
            += particle_stat_weight;
        }
        if (cs_glob_lagr_boundary_interactions->iencmabd > 0) {
          bound_stat[  cs_glob_lagr_boundary_interactions->iencma
                     * n_b_faces + face_id]
            += particle_stat_weight * particle_mass / face_area;
        }
        if (cs_glob_lagr_boundary_interactions->iencdibd > 0) {
          bound_stat[  cs_glob_lagr_boundary_interactions->iencdi
                     * n_b_faces + face_id]
            +=   particle_stat_weight
               * cs_lagr_particle_get_real(particle, p_am,
                                           CS_LAGR_SHRINKING_DIAMETER);
        }
        if (cs_glob_lagr_boundary_interactions->iencckbd > 0) {
          if (particle_mass > 0) {
            const cs_real_t *particle_coal_mass
              = cs_lagr_particle_attr_const(particle, p_am,
                                            CS_LAGR_COAL_MASS);
            const cs_real_t *particle_coke_mass
              = cs_lagr_particle_attr_const(particle, p_am,
                                            CS_LAGR_COKE_MASS);
            for (int k = 0; k < n_layers; k++) {
              bound_stat[  cs_glob_lagr_boundary_interactions->iencck
                         * n_b_faces + face_id]
                +=   particle_stat_weight
                   * (particle_coal_mass[k] + particle_coke_mass[k])
                   / particle_mass;
            }
          }
        }

        /* FIXME: For post-processing by trajectory purpose */

        for (int k = 0; k < 3; k++) {
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

    cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);
    cs_real_3_t vect_cen;
    for (int k = 0; k < 3; k++) {
      vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
      p_info->start_coords[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
    }

      /* Modify the ending point. */

      for (int k = 0; k < 3; k++)
        disp[k] = particle_coord[k] - intersect_pt[k];

      tmp = 2. * cs_math_3_dot_product(disp, face_norm);

      for (int k = 0; k < 3; k++)
        particle_coord[k] -= tmp * face_norm[k];

      /* Modify particle velocity and velocity seen */

      tmp = 2. * cs_math_3_dot_product(particle_velocity, face_norm);

      for (int k = 0; k < 3; k++)
        particle_velocity[k] -= tmp * face_norm[k];

      tmp = 2. * cs_math_3_dot_product(particle_velocity_seen, face_norm);

      for (int k = 0; k < 3; k++) {
        particle_velocity_seen[k] -= tmp * face_norm[k];
        // particle_velocity_seen[k] = 0.0; //FIXME
      }

    }

  }

  /* FIXME: JBORD* (user-defined boundary condition) not yet implemented
     nor defined by a macro */
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Boundary condition %d not recognized.\n"),
              b_type);

  /* Ensure some fields are updated */

  if (p_am->size[CS_LAGR_DEPOSITION_FLAG] > 0) {

    cs_lnum_t depo_flag
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG);

    if (   depo_flag == CS_LAGR_PART_ROLLING
        || depo_flag == CS_LAGR_PART_DEPOSITED) {
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                cs_glob_mesh->b_face_cells[face_id] + 1);
      cs_lagr_particle_set_lnum(particle, p_am,CS_LAGR_NEIGHBOR_FACE_ID ,
                                face_id);
    }

  }

  /* Return pointer */

  *p_move_particle = move_particle;

  /* Update per-zone flow rate measure for exiting particles */

  if (particle_state == CS_LAGR_PART_OUT) {
    int n_stats = cs_glob_lagr_model->n_stat_classes + 1;

    cs_real_t fr =   particle_stat_weight
                   * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

    bdy_conditions->particle_flow_rate[boundary_zone*n_stats] -= fr;

    if (n_stats > 1) {
      int class_id
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_STAT_CLASS);
      if (class_id > 0 && class_id < n_stats)
        bdy_conditions->particle_flow_rate[  boundary_zone*n_stats
                                           + class_id] -= fr;
    }

  }

  /* FIXME: Post-treatment not yet implemented... */

  if  (   b_type == CS_LAGR_DEPO1
       || b_type == CS_LAGR_DEPO2
       || b_type == CS_LAGR_DEPO_DLVO
       || b_type == CS_LAGR_REBOUND
       || b_type == CS_LAGR_FOULING) {

    /* Number of particle-boundary interactions  */
    if (cs_glob_lagr_boundary_interactions->has_part_impact_nbr > 0)
      bound_stat[cs_glob_lagr_boundary_interactions->inbr * n_b_faces + face_id]
        += particle_stat_weight;

    /* Particle impact angle and velocity*/
    if (cs_glob_lagr_boundary_interactions->iangbd > 0) {
      cs_real_t imp_ang = acos(cs_math_3_dot_product(compo_vel, face_normal)
                               / (face_area * norm_vel));
      bound_stat[cs_glob_lagr_boundary_interactions->iang * n_b_faces + face_id]
        += imp_ang * particle_stat_weight;
    }

    if (cs_glob_lagr_boundary_interactions->ivitbd > 0)
      bound_stat[cs_glob_lagr_boundary_interactions->ivit * n_b_faces + face_id]
        += norm_vel * particle_stat_weight;

    /* User statistics management. By defaut, set to zero */
    if (cs_glob_lagr_boundary_interactions->nusbor > 0)
      for (int n1 = 0; n1 < cs_glob_lagr_boundary_interactions->nusbor; n1++)
        bound_stat[cs_glob_lagr_boundary_interactions->iusb[n1] * n_b_faces + face_id] = 0.0;
  }

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Move a particle as far as possible while remaining on a given rank.
 *
 * parameters:
 *   particle                 <-> pointer to particle data
 *   p_am                     <-- particle attribute map
 *   displacement_step_id     <-- id of displacement step
 *   failsafe_mode            <-- with (0) / without (1) failure capability
 *   b_face_zone_id           <-- boundary face zone id
 *   visc_length              <-- viscous layer thickness
 *
 * returns:
 *   a state associated to the status of the particle (treated, to be deleted,
 *   to be synchonised)
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_local_propagation(void                           *particle,
                   const cs_lagr_attribute_map_t  *p_am,
                   int                             displacement_step_id,
                   int                             failsafe_mode,
                   const int                       b_face_zone_id[],
                   const cs_real_t                 visc_length[],
                   const cs_field_t               *u)
{
  cs_lnum_t  i;
  cs_real_t  disp[3];
  cs_real_t  null_yplus;

  cs_lnum_t  *neighbor_face_id;
  cs_real_t  *particle_yplus;

  int move_particle = CS_LAGR_PART_MOVE_ON;
  cs_lagr_tracking_state_t  particle_state = CS_LAGR_PART_TO_SYNC;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  cs_lnum_t  *cell_face_idx = builder->cell_face_idx;
  cs_lnum_t  *cell_face_lst = builder->cell_face_lst;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  cs_real_t  *particle_coord
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);
  cs_real_t  *prev_location = p_info->start_coords;

  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

  for (int k = 0; k < 3; k++)
    disp[k] = particle_coord[k] - prev_location[k];

  cs_lnum_t  cell_id
    = cs_lagr_particle_get_cell_id(particle, p_am);

  const cs_real_3_t *vtx_coord
    = (const cs_real_3_t *)(cs_glob_mesh->vtx_coord);
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  /* Dimension less test: no movement? */
  cs_real_t inv_ref_length = 1./pow(cell_vol[cell_id], 1./3.);
  if (fabs(disp[0] * inv_ref_length) < 1e-15 &&
      fabs(disp[1] * inv_ref_length) < 1e-15 &&
      fabs(disp[2] * inv_ref_length) < 1e-15 ) {
    move_particle = CS_LAGR_PART_MOVE_OFF;
    particle_state = CS_LAGR_PART_TREATED;
  }

  if (lagr_model->deposition > 0) {
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

  for (int n_loops = displacement_step_id;
       move_particle == CS_LAGR_PART_MOVE_ON;
       n_loops++) {

    cell_id
      = cs_lagr_particle_get_cell_id(particle, p_am);

    assert(cell_id < mesh->n_cells);
    assert(cell_id > -1);

    if (n_loops > _max_propagation_loops) { /* Manage error */

      _manage_error(failsafe_mode,
                    particle,
                    p_am,
                    CS_LAGR_TRACKING_ERR_MAX_LOOPS);

      move_particle  = CS_LAGR_PART_MOVE_OFF;
      particle_state = CS_LAGR_PART_TREATED;

      cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);

      for (int k = 0; k < 3; k++) {
        prev_location[k] = cell_cen[k];
        particle_coord[k] = cell_cen[k];
        disp[k] = 0.;
      }

      cs_lnum_t n_rep
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_TR_REPOSITION);
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_TR_REPOSITION, n_rep+1);

      return particle_state;

    }

    /* Treatment for depositing particles */

    if (lagr_model->deposition > 0 && *particle_yplus < 0.) {

      cs_lagr_test_wall_cell(particle, p_am, visc_length,
                             particle_yplus, neighbor_face_id);

      if (*particle_yplus < 100.) {

        cs_real_t flow_velo_x, flow_velo_y, flow_velo_z;

        const cs_real_3_t *rot_m
          = (const cs_real_3_t *)(  cs_glob_lagr_b_face_proj
                                  + *neighbor_face_id);

        flow_velo_x = u->val[cell_id*3];
        flow_velo_y = u->val[cell_id*3 + 1];
        flow_velo_z = u->val[cell_id*3 + 2];

        /* e1 (normal) vector coordinates */
        cs_real_t e1_x = cs_glob_lagr_b_u_normal[*neighbor_face_id][0];
        cs_real_t e1_y = cs_glob_lagr_b_u_normal[*neighbor_face_id][1];
        cs_real_t e1_z = cs_glob_lagr_b_u_normal[*neighbor_face_id][2];

        /* e2 vector coordinates */
        cs_real_t e2_x = rot_m[1][0];
        cs_real_t e2_y = rot_m[1][1];
        cs_real_t e2_z = rot_m[1][2];

        /* e3 vector coordinates */
        cs_real_t e3_x = rot_m[2][0];
        cs_real_t e3_y = rot_m[2][1];
        cs_real_t e3_z = rot_m[2][2];

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

    bool restart = false;

  reloop_cen:;

    cs_lnum_t exit_face = 0; /* > 0 for interior faces,
                                < 0 for boundary faces */

    double adist_min = 2.;

    double t_intersect = -1;
    cs_real_3_t face_norm;

    int n_in = 0;
    int n_out = 0;

    const cs_real_t  *next_location
      = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_COORDS);

    /* Loop on faces to see if the particle trajectory crosses it*/
    for (i = cell_face_idx[cell_id];
         i < cell_face_idx[cell_id+1] && move_particle == CS_LAGR_PART_MOVE_ON;
         i++) {

      cs_lnum_t face_id, vtx_start, vtx_end, n_vertices;
      const cs_lnum_t *face_connect;
      const cs_real_t *face_cog;

      /* Outward normal: always well oriented for external faces, depend on the
       * connectivity for internal faces */
      int reorient_face = 1;

      cs_lnum_t face_num = cell_face_lst[i];

      if (face_num > 0) {

        /* Interior face */

        face_id = face_num - 1;
        if (cell_id == mesh->i_face_cells[face_id][1])
          reorient_face = -1;
        vtx_start = mesh->i_face_vtx_idx[face_id];
        vtx_end = mesh->i_face_vtx_idx[face_id+1];
        n_vertices = vtx_end - vtx_start;

        face_connect = mesh->i_face_vtx_lst + vtx_start;
        face_cog = i_face_cog[face_id];

      }
      else {

        assert(face_num < 0);

        /* Boundary faces */

        face_id = -face_num - 1;
        vtx_start = mesh->b_face_vtx_idx[face_id];
        vtx_end = mesh->b_face_vtx_idx[face_id+1];
        n_vertices = vtx_end - vtx_start;

        face_connect = mesh->b_face_vtx_lst + vtx_start;
        face_cog = b_face_cog[face_id];

      }

      /*
        adimensional distance estimation of face intersection
        (2 if no chance of intersection)
      */

      int n_crossings[2] = {0, 0};

      double t = cs_geom_segment_intersect_face(reorient_face,
                                                n_vertices,
                                                face_connect,
                                                vtx_coord,
                                                face_cog,
                                                prev_location,
                                                next_location,
                                                n_crossings,
                                                face_norm);

      n_in += n_crossings[0];
      n_out += n_crossings[1];

      /* Store the nearest intesection from the O point...*/
      if (t < adist_min) {
        exit_face = face_num;
        t_intersect = t;
        adist_min = t_intersect;
      }

    }

    /* We test here if the particle is truly within the current cell
     * (meaning n_in = n_out > 0 )
     * If there is a problem (pb with particle strictly // and on the face ?),
     * the particle initial position is replaced at the cell center
     * and we continue the trajectory analysis. */

    bool test_in = (n_in == 0 && n_out == 0);

    if ((n_in != n_out || test_in)
        && (move_particle == CS_LAGR_PART_MOVE_ON)) {
      cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);

      for (int k = 0; k < 3; k++)
        prev_location[k] = cell_cen[k];

      cs_lnum_t n_rep
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_TR_REPOSITION);
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_TR_REPOSITION, n_rep+1);

      if (!(restart))
        restart = true;
      else {
        _manage_error(failsafe_mode,
                      particle,
                      p_am,
                      CS_LAGR_TRACKING_ERR_LOST_PIC);

        move_particle  = CS_LAGR_PART_MOVE_OFF;
        particle_state = CS_LAGR_PART_TREATED;
        return particle_state;
      }
      goto reloop_cen;
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

      /* Deposition on internal faces? */

      /* particle / internal condition interaction
         1 - modify particle cell_num : 0 or boundary_cell_num
         2 -

         O -->  *         *  <-- D
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
        = _internal_treatment(cs_glob_lagr_particle_set,
                              particle,
                              face_id,
                              t_intersect,
                              &move_particle);

      if (move_particle != CS_LAGR_PART_MOVE_OFF) {

        /* Now update cell */

        if (cell_id == c_id1) {
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                    c_id2 + 1);
          cell_id = c_id2;
        }

        else {
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_NUM,
                                    c_id1 + 1);
          cell_id = c_id1;
        }

        /* Particle changes rank */

        if (cell_id >= mesh->n_cells) {

          particle_state = CS_LAGR_PART_TO_SYNC;
          move_particle = CS_LAGR_PART_MOVE_OFF;

          /* Specific treatment for the particle deposition model */

          if (lagr_model->deposition > 0 && *particle_yplus < 100.) {

            /* Marking of particles */

            *particle_yplus = - *particle_yplus;
            //FIXME... a projection was made, it
            //is safer to project when we use it !

          }

        } /* end of case where particle changes rank */

        else if (lagr_model->deposition > 0) {

          /* Specific treatment for the particle deposition model */

          cs_real_t save_yplus = *particle_yplus;

          /* Wall cell detection */

          cs_lagr_test_wall_cell(particle, p_am, visc_length,
                                 particle_yplus, neighbor_face_id);

          if (save_yplus < 100.) {

            cs_real_t intersect_pt[3], vect_cen[3];
            cs_real_t bc_epsilon = 1e-8;
            cs_real_t *cell_cen = fvq->cell_cen + (3*cell_id);

            for (int k = 0; k < 3; k++)
              intersect_pt[k]
                =   (particle_coord[k] - prev_location[k])*t_intersect
                  + prev_location[k];

            for (int k = 0; k < 3; k++) {
              vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
              particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
            }

            if (*particle_yplus < 100.0) {

              cs_real_t flow_velo_x = u->val[cell_id*3];
              cs_real_t flow_velo_y = u->val[cell_id*3 + 1];
              cs_real_t flow_velo_z = u->val[cell_id*3 + 2];

              /* The particle is still in the boundary layer */

              const cs_real_3_t *rot_m
                = (const cs_real_3_t *)(  cs_glob_lagr_b_face_proj
                                        + *neighbor_face_id);

              /* e1 (normal) vector coordinates */
              cs_real_t e1_x = cs_glob_lagr_b_u_normal[*neighbor_face_id][0];
              cs_real_t e1_y = cs_glob_lagr_b_u_normal[*neighbor_face_id][1];
              cs_real_t e1_z = cs_glob_lagr_b_u_normal[*neighbor_face_id][2];

              /* e2 vector coordinates */
              cs_real_t e2_x = rot_m[1][0];
              cs_real_t e2_y = rot_m[1][1];
              cs_real_t e2_z = rot_m[1][2];

              /* e3 vector coordinates */
              cs_real_t e3_x = rot_m[2][0];
              cs_real_t e3_y = rot_m[2][1];
              cs_real_t e3_z = rot_m[2][2];

              cs_real_t old_fl_seen_norm
                =   particle_velocity_seen[0] * e1_x
                  + particle_velocity_seen[1] * e1_y
                  + particle_velocity_seen[2] * e1_z;
              //FIXME check that it is the new normal to use...

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
          } /* end of case for y+ <100 */

        } /* end of case for deposition model */

      }

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
        = _boundary_treatment(cs_glob_lagr_particle_set,
                              particle,
                              face_num,
                              face_norm,
                              t_intersect,
                              b_face_zone_id[face_num-1],
                              &move_particle);

      if (cs_glob_lagr_time_scheme->t_order == 2)
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_REBOUND_ID, 0);

      assert(   move_particle == CS_LAGR_PART_MOVE_ON
             || move_particle == CS_LAGR_PART_MOVE_OFF);

      p_info->last_face_num = -face_num;

    } /* end if exit_face < 0 */

  } /* End of while : local displacement */

  assert(   move_particle != CS_LAGR_PART_MOVE_ON
         || particle_state != CS_LAGR_PART_TO_SYNC);

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
 * Exchange particles
 *
 * parameters:
 *  halo      <-- pointer to a cs_halo_t structure
 *  lag_halo  <-> pointer to a cs_lagr_halo_t structure
 *  particles <-- set of particles to update
 *----------------------------------------------------------------------------*/

static void
_exchange_particles(const cs_halo_t         *halo,
                    cs_lagr_halo_t          *lag_halo,
                    cs_lagr_particle_set_t  *particles)
{
  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

  const size_t tot_extents = lag_halo->extents;

  cs_lnum_t  n_recv_particles = 0;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int  rank;
    int  request_count = 0;
    const int  local_rank = cs_glob_rank_id;

    /* Receive data from distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      cs_lnum_t shift =   particles->n_particles
                        + lag_halo->recv_shift[rank];

      if (lag_halo->recv_count[rank] > 0) {

        if (halo->c_domain_rank[rank] != local_rank) {
          void  *recv_buf = particles->p_buffer + tot_extents*shift;
          n_recv_particles += lag_halo->recv_count[rank];
          MPI_Irecv(recv_buf,
                    lag_halo->recv_count[rank],
                    _cs_mpi_particle_type,
                    halo->c_domain_rank[rank],
                    halo->c_domain_rank[rank],
                    cs_glob_mpi_comm,
                    &(lag_halo->request[request_count++]));
        }
        else
          local_rank_id = rank;

      }

    }

    /* We wait for posting all receives
       (often recommended in the past, apparently not anymore) */

#if 0
    MPI_Barrier(cs_glob_mpi_comm);
#endif

    /* Send data to distant ranks */

    for (rank = 0; rank < halo->n_c_domains; rank++) {

      /* If this is not the local rank */

      if (   halo->c_domain_rank[rank] != local_rank
          && lag_halo->send_count[rank] > 0) {
        cs_lnum_t shift = lag_halo->send_shift[rank];
        void  *send_buf = lag_halo->send_buf + tot_extents*shift;
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

      cs_lnum_t  recv_shift =   particles->n_particles
                              + lag_halo->recv_shift[local_rank_id];
      cs_lnum_t  send_shift = lag_halo->send_shift[local_rank_id];

      assert(   lag_halo->recv_count[local_rank_id]
             == lag_halo->send_count[local_rank_id]);

      n_recv_particles += lag_halo->send_count[local_rank_id];

      for (cs_lnum_t i = 0; i < lag_halo->send_count[local_rank_id]; i++)
        memcpy(particles->p_buffer + tot_extents*(recv_shift + i),
               lag_halo->send_buf + tot_extents*(send_shift + i),
               tot_extents);

    }
  }

  /* Update particle count and weight */

  cs_real_t tot_weight = 0.;

  for (cs_lnum_t i = 0; i < n_recv_particles; i++) {

    cs_real_t cur_part_stat_weight
      = cs_lagr_particles_get_real(particles,
                                   particles->n_particles + i,
                                   CS_LAGR_STAT_WEIGHT);

    tot_weight += cur_part_stat_weight;

  }

  particles->n_particles += n_recv_particles;
  particles->weight += tot_weight;
}

/*----------------------------------------------------------------------------
 * Determine particle halo sizes
 *
 * parameters:
 *   mesh      <-- pointer to associated mesh
 *   lag_halo  <-> pointer to particle halo structure to update
 *   particles <-- set of particles to update
 *----------------------------------------------------------------------------*/

static void
_lagr_halo_count(const cs_mesh_t               *mesh,
                 cs_lagr_halo_t                *lag_halo,
                 const cs_lagr_particle_set_t  *particles)
{
  cs_lnum_t  i, ghost_id;

  cs_lnum_t  n_recv_particles = 0, n_send_particles = 0;

  const cs_halo_t  *halo = mesh->halo;

  /* Initialization */

  for (i = 0; i < halo->n_c_domains; i++) {
    lag_halo->send_count[i] = 0;
    lag_halo->recv_count[i] = 0;
  }

  /* Loop on particles to count number of particles to send on each rank */

  for (i = 0; i < particles->n_particles; i++) {

    if (_get_tracking_info(particles, i)->state == CS_LAGR_PART_TO_SYNC) {

      ghost_id =   cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_NUM)
                 - mesh->n_cells - 1;

      assert(ghost_id >= 0);
      lag_halo->send_count[lag_halo->rank[ghost_id]] += 1;

    }

  } /* End of loop on particles */

  /* Exchange counters */

  _exchange_counter(halo, lag_halo);

  for (i = 0; i < halo->n_c_domains; i++) {
    n_recv_particles += lag_halo->recv_count[i];
    n_send_particles += lag_halo->send_count[i];
  }

  lag_halo->send_shift[0] = 0;
  lag_halo->recv_shift[0] = 0;

  for (i = 1; i < halo->n_c_domains; i++) {

    lag_halo->send_shift[i] =  lag_halo->send_shift[i-1]
                             + lag_halo->send_count[i-1];

    lag_halo->recv_shift[i] =  lag_halo->recv_shift[i-1]
                             + lag_halo->recv_count[i-1];

  }

  /* Resize particle set and/or halo only if needed */

  cs_lagr_particle_set_resize(particles->n_particles + n_recv_particles);
  _resize_lagr_halo(lag_halo, n_send_particles);
}

/*----------------------------------------------------------------------------
 * Update particle sets, including halo synchronization.
 *
 * parameters:
 *   particles <-> set of particles to update
 *
 * returns:
 *   1 if displacement needs to continue, 0 if finished
 *----------------------------------------------------------------------------*/

static int
_sync_particle_set(cs_lagr_particle_set_t  *particles)
{
  cs_lnum_t  i, k, tr_id, rank, shift, ghost_id;
  cs_real_t matrix[3][4];

  cs_lnum_t  n_recv_particles = 0;
  cs_lnum_t  particle_count = 0;

  cs_lnum_t  n_merged_particles = 0;

  cs_lnum_t  n_exit_particles = 0;
  cs_lnum_t  n_failed_particles = 0;

  cs_real_t  exit_weight = 0.0;
  cs_real_t  fail_weight = 0.0;
  cs_real_t  tot_weight = 0.0;

  cs_lagr_track_builder_t  *builder = _particle_track_builder;
  cs_lagr_halo_t  *lag_halo = builder->halo;

  const cs_lagr_attribute_map_t *p_am = particles->p_am;
  const size_t extents = particles->p_am->extents;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const fvm_periodicity_t *periodicity = mesh->periodicity;
  const cs_interface_set_t  *face_ifs = builder->face_ifs;

  int continue_displacement = 0;

  if (halo != NULL) {

    _lagr_halo_count(mesh, lag_halo, particles);

    for (i = 0; i < halo->n_c_domains; i++) {
      n_recv_particles += lag_halo->recv_count[i];
      lag_halo->send_count[i] = 0;
    }
  }

  /* Loop on particles, transferring particles to synchronize to send_buf
     for particle set, and removing particles that otherwise exited the domain */

  for (i = 0; i < particles->n_particles; i++) {

    cs_lagr_tracking_state_t cur_part_state
      = _get_tracking_info(particles, i)->state;

    cs_real_t cur_part_stat_weight
      = cs_lagr_particles_get_real(particles, i, CS_LAGR_STAT_WEIGHT);

    /* Particle changes domain */

    if (cur_part_state == CS_LAGR_PART_TO_SYNC) {

      continue_displacement = 1;

      ghost_id =   cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_NUM)
                 - halo->n_local_elts - 1;
      rank = lag_halo->rank[ghost_id];
      tr_id = lag_halo->transform_id[ghost_id];

      cs_lagr_particles_set_lnum(particles, i, CS_LAGR_CELL_NUM,
                                 lag_halo->dist_cell_num[ghost_id]);

      shift = lag_halo->send_shift[rank] + lag_halo->send_count[rank];

      /* Update if needed last_face_num */

      if (tr_id >= 0) { /* Same initialization as in previous algorithm */

        _tracking_info(particles, i)->last_face_num = 0;

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
                    _get_tracking_info(particles, i)->last_face_num  - 1,
                    local_num);

            if (id == -1)
              bft_error(__FILE__, __LINE__, 0,
                        _(" Cannot find the relative distant face num.\n"));

            dist_num = cs_interface_get_match_ids(interface);

            _tracking_info(particles, i)->last_face_num = dist_num[id] + 1;

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
                              cs_lagr_particles_attr(particles, i,
                                                     CS_LAGR_COORDS));

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              _tracking_info(particles, i)->start_coords);

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              cs_lagr_particles_attr_n(particles, i, 1,
                                                       CS_LAGR_COORDS));

        /* Apply rotation to velocity vectors in case of rotation */

        if (perio_type >= FVM_PERIODICITY_ROTATION) {

          /* Rotation of the velocity */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(particles, i,
                                                        CS_LAGR_VELOCITY));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr_n(particles, i, 1,
                                                          CS_LAGR_VELOCITY));

          /* Rotation of the velocity seen */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr(particles, i,
                                                        CS_LAGR_VELOCITY_SEEN));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
                                 cs_lagr_particles_attr_n(particles, i, 1,
                                                          CS_LAGR_VELOCITY_SEEN));

        } /* Specific treatment in case of rotation for the velocities */

      } /* End of periodicity treatment */

      memcpy(lag_halo->send_buf + extents*shift,
             particles->p_buffer + extents*i,
             extents);

      lag_halo->send_count[rank] += 1;

      /* Remove the particle from the local set (do not copy it) */

    } /* TO_SYNC */

    /* Particle remains in domain */

    else if (cur_part_state == CS_LAGR_PART_MERGED) {
      n_merged_particles++;
    }

    else if (cur_part_state < CS_LAGR_PART_OUT) {

      if (particle_count < i)
        memcpy(particles->p_buffer + p_am->extents*particle_count,
               particles->p_buffer + p_am->extents*i,
               p_am->extents);

      particle_count += 1;
      tot_weight += cur_part_stat_weight;

    }

    /* Particle exits domain */

    else if (cur_part_state < CS_LAGR_PART_ERR) {
      n_exit_particles++;
      exit_weight += cur_part_stat_weight;
    }

    else {
      n_failed_particles++;
      fail_weight += cur_part_stat_weight;
    }

  } /* End of loop on particles */

  particles->n_particles = particle_count;
  particles->weight = tot_weight;

  particles->n_part_out += n_exit_particles;
  particles->weight_out += exit_weight;

  particles->n_failed_part += n_failed_particles;
  particles->weight_failed += fail_weight;

  particles->n_part_merged += n_merged_particles;

  /* Exchange particles, then update set */

  if (halo != NULL)
    _exchange_particles(halo, lag_halo, particles);

  cs_parall_max(1, CS_INT_TYPE, &continue_displacement);

  return continue_displacement;
}

/*----------------------------------------------------------------------------
 * Prepare for particle movement phase
 *
 * parameters:
 *   particles        <-> pointer to particle set structure
 *   part_b_mass_flux <-> particle mass flux array, or NULL
 *----------------------------------------------------------------------------*/

static void
_initialize_displacement(cs_lagr_particle_set_t  *particles,
                         cs_real_t                part_b_mass_flux[])
{
  cs_lnum_t  i;

  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;

  const cs_lagr_attribute_map_t  *am = particles->p_am;

  const cs_real_t  *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  /* Initialize builder if needed */

  if (_particle_track_builder == NULL)
    _particle_track_builder
      = _init_track_builder(particles->n_particles_max,
                            particles->p_am->extents);

  assert(am->lb >= sizeof(cs_lagr_tracking_info_t));

  /* Info for rotor-stator cases; the time step should actually
     be based on the global (non-Lagrangian) time step in case
     the two differ. Currently both are in lock-step; if one
     becomes a mutiplier of the other, the rotor-movement update
     would need to be done at a frequency matching that multiplier. */

  const int *cell_rotor_num = NULL;
  cs_real_34_t  *rot_m = NULL;

  if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
    cell_rotor_num = cs_turbomachinery_get_cell_rotor_num();
    cs_real_t dt = cs_glob_lagr_time_step->dtp;
    rot_m = cs_turbomachinery_get_rotation_matrices(-1.0 * dt);
  }

  /* Prepare tracking info */

  for (i = 0; i < particles->n_particles; i++) {

    cs_lnum_t cur_part_cell_num
      = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_NUM);
    cs_real_t r_truncate
      = cs_lagr_particles_get_real(particles, i, CS_LAGR_TR_TRUNCATE);

    if (am->size[CS_LAGR_DEPOSITION_FLAG] > 0) {
      if (    cs_lagr_particles_get_lnum(particles, i, CS_LAGR_DEPOSITION_FLAG)
           == CS_LAGR_PART_TO_DELETE) {
        _tracking_info(particles, i)->state = CS_LAGR_PART_OUT;
      }
    }

    if (cur_part_cell_num < 0) {
      cs_real_t particle_weight
        = cs_lagr_particles_get_real(particles, i, CS_LAGR_STAT_WEIGHT);
      if (particle_weight <= 0.) {
        _tracking_info(particles, i)->state = CS_LAGR_PART_MERGED;
      }
      else {
        _tracking_info(particles, i)->state = CS_LAGR_PART_STUCK;
      }
    }
    else if (r_truncate > 1.9) /* from previous displacement */
      _tracking_info(particles, i)->state = CS_LAGR_PART_ERR;
    else {
      _tracking_info(particles, i)->state = CS_LAGR_PART_TO_SYNC;
      if (am->size[CS_LAGR_DEPOSITION_FLAG] > 0) {
        if (    cs_lagr_particles_get_lnum(particles, i, CS_LAGR_DEPOSITION_FLAG)
            == CS_LAGR_PART_DEPOSITED)
          _tracking_info(particles, i)->state = CS_LAGR_PART_TREATED;
      }
    }

    _tracking_info(particles, i)->last_face_num = 0;

    /* Coordinates of the particle */

    cs_real_t *prv_part_coord
      = cs_lagr_particles_attr_n(particles, i, 1, CS_LAGR_COORDS);

    _tracking_info(particles, i)->start_coords[0] = prv_part_coord[0];
    _tracking_info(particles, i)->start_coords[1] = prv_part_coord[1];
    _tracking_info(particles, i)->start_coords[2] = prv_part_coord[2];

    if (cell_rotor_num != NULL) {
      cs_lnum_t  cell_id = CS_ABS(cur_part_cell_num) - 1;
      int r_num = cell_rotor_num[cell_id];
      if (r_num > 0) {
        _apply_vector_transfo((const cs_real_t (*)[4])rot_m[r_num],
                              _tracking_info(particles, i)->start_coords);
      }
    }

    /* Just after injection, reduce displacment so as to simulate
       continuous injection */

    cs_real_t res_time = cs_lagr_particles_get_real(particles, i,
                                                    CS_LAGR_RESIDENCE_TIME);

    if (res_time < 0) {
      cs_real_t fraction =   (cs_glob_lagr_time_step->dtp + res_time)
                           / cs_glob_lagr_time_step->dtp;
      cs_real_t *part_coord
        = cs_lagr_particles_attr(particles, i, CS_LAGR_COORDS);
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_real_t d = part_coord[j] - prv_part_coord[j];
        part_coord[j] = prv_part_coord[j] + fraction*d;
      }
    }

    cs_lagr_particles_set_real(particles, i, CS_LAGR_TR_TRUNCATE, 0);
    cs_lagr_particles_set_lnum(particles, i, CS_LAGR_TR_REPOSITION, 0);

    /* Data needed if the deposition model is activated */
    if (   lagr_model->deposition <= 0
        && am->size[CS_LAGR_DEPOSITION_FLAG] > 0)
      cs_lagr_particles_set_lnum(particles, i, CS_LAGR_DEPOSITION_FLAG,
                                 CS_LAGR_PART_IN_FLOW);

    /* Remove contribution from deposited or rolling particles
       to boundary mass flux at the beginning of their movement. */

    else if (lagr_model->deposition > 0 && part_b_mass_flux != NULL)
      _b_mass_contribution(particles,
                           i,
                           -1.0,
                           b_face_surf,
                           part_b_mass_flux);

  }

  BFT_FREE(rot_m);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Particle set after %s\n", __func__);
  cs_lagr_particle_set_dump(particles);
#endif
}

/*----------------------------------------------------------------------------
 * Update particle set structures: compact array.
 *
 * parameters:
 *   particles        <-> pointer to particle set structure
 *   part_b_mass_flux <-> particle mass flux array, or NULL
 *----------------------------------------------------------------------------*/

static void
_finalize_displacement(cs_lagr_particle_set_t  *particles,
                       cs_real_t                part_b_mass_flux[])
{
  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;
  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;
  const cs_real_t  *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  const cs_lnum_t n_particles = particles->n_particles;

  cs_lnum_t *cell_idx;
  unsigned char *swap_buffer;
  size_t swap_buffer_size = p_am->extents * ((size_t)n_particles);

  BFT_MALLOC(cell_idx, n_cells+1, cs_lnum_t);
  BFT_MALLOC(swap_buffer, swap_buffer_size, unsigned char);

  /* Cell index (count first) */

  for (cs_lnum_t i = 0; i < n_cells+1; i++)
    cell_idx[i] = 0;

  /* Copy unordered particle data to buffer */

  for (cs_lnum_t i = 0; i < n_particles; i++) {

    cs_lnum_t cur_part_state = _get_tracking_info(particles, i)->state;

    assert(   cur_part_state < CS_LAGR_PART_OUT
           && cur_part_state != CS_LAGR_PART_TO_SYNC);

    cs_lnum_t cell_num = cs_lagr_particles_get_lnum(particles, i,
                                                    CS_LAGR_CELL_NUM);

    assert(cell_num != 0);

    if (cell_num < 0)
      cell_num = - cell_num;

    cs_lnum_t cell_id = cell_num - 1;

    memcpy(swap_buffer + p_am->extents*i,
           particles->p_buffer + p_am->extents*i,
           p_am->extents);

    cell_idx[cell_id+1] += 1;

  }

  /* Convert count to index */

  for (cs_lnum_t i = 1; i < n_cells; i++)
    cell_idx[i+1] += cell_idx[i];

  assert(n_particles == cell_idx[n_cells]);

  /* Now copy particle data and update some statistics */

  const cs_lnum_t p_extents = particles->p_am->extents;
  const cs_lnum_t cell_num_displ = particles->p_am->displ[0][CS_LAGR_CELL_NUM];

  for (cs_lnum_t i = 0; i < n_particles; i++) {

    cs_lnum_t cell_num
      = *((const cs_lnum_t *)(swap_buffer + p_extents*i + cell_num_displ));

    assert(cell_num != 0);

    if (cell_num < 0)
      cell_num = - cell_num;

    cs_lnum_t cell_id = cell_num - 1;

    cs_lnum_t particle_id = cell_idx[cell_id];

    cell_idx[cell_id] += 1;

    memcpy(particles->p_buffer + p_am->extents*particle_id,
           swap_buffer + p_am->extents*i,
           p_am->extents);

    /* Add contribution from deposited or rolling particles
       to boundary mass flux at the end of their movement. */

    if (lagr_model->deposition > 0 && part_b_mass_flux != NULL)
      _b_mass_contribution(particles,
                           particle_id,
                           1.0,
                           b_face_surf,
                           part_b_mass_flux);

  }

  BFT_FREE(swap_buffer);
  BFT_FREE(cell_idx);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Particle set after %s\n", __func__);
  cs_lagr_particle_set_dump(particles);
#endif

  if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT)
    _particle_track_builder = _destroy_track_builder(_particle_track_builder);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize particle tracking subsystem
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_tracking_initialize(void)
{
  /* Initialize particle set */

  cs_lagr_particle_set_create();

  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n PARTICLE SET AFTER CREATION\n");
  cs_lagr_particle_set_dump(p_set);
#endif

  /* Initialization */

  for (cs_lnum_t i = 0; i < p_set->n_particles_max; i++)
    _tracking_info(p_set, i)->state = CS_LAGR_PART_TO_SYNC;

  /* Create all useful MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _cs_mpi_particle_type = _define_particle_datatype(p_set->p_am);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply one particle movement step.
 *
 * \param[in]  visc_length     viscous layer thickness
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_tracking_particle_movement(const cs_real_t  visc_length[])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;

  int  displacement_step_id = 0;
  int  continue_displacement = 1;

  cs_lagr_particle_set_t  *particles = cs_glob_lagr_particle_set;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;

  const cs_lnum_t  failsafe_mode = 0; /* If 1 : stop as soon as an error is
                                         detected */

  const cs_field_t *u = cs_glob_lagr_extra_module->vel;

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_real_t *part_b_mass_flux = NULL;

  int t_stat_id = cs_timer_stats_id_by_name("particle_displacement_stage");

  int t_top_id = cs_timer_stats_switch(t_stat_id);

  if (cs_glob_lagr_boundary_interactions->iflmbd) {
    assert(cs_glob_lagr_boundary_interactions->iflm >= 0);
    part_b_mass_flux = bound_stat
      + (cs_glob_lagr_boundary_interactions->iflm * mesh->n_b_faces);
  }

  assert(particles != NULL);

  const int *b_face_zone_id = cs_boundary_zone_face_class_id();

  /* particles->n_part_new: handled in injection step */

  particles->weight = 0.0;
  particles->n_part_out = 0;
  particles->n_part_dep = 0;
  particles->n_part_fou = 0;
  particles->weight_out = 0.0;
  particles->weight_dep = 0.0;
  particles->weight_fou = 0.0;
  particles->n_failed_part = 0;
  particles->weight_failed = 0.0;

  _initialize_displacement(particles,
                           part_b_mass_flux);

  /* Main loop on particles: global propagation */

  while (continue_displacement) {

    /* Local propagation */

    for (cs_lnum_t i = 0; i < particles->n_particles; i++) {

      unsigned char *particle = particles->p_buffer + p_am->extents * i;

      /* Local copies of the current and previous particles state vectors
         to be used in case of the first pass of _local_propagation fails */

      cs_lagr_tracking_state_t cur_part_state
        = _get_tracking_info(particles, i)->state;

      if (cur_part_state == CS_LAGR_PART_TO_SYNC) {

        /* Main particle displacement stage */

        cur_part_state = _local_propagation(particle,
                                            p_am,
                                            displacement_step_id,
                                            failsafe_mode,
                                            b_face_zone_id,
                                            visc_length,
                                            u);

        _tracking_info(particles, i)->state = cur_part_state;

      }

    } /* End of loop on particles */

    /* Update of the particle set structure. Delete exited particles,
       update for particles which change domain. */

    continue_displacement = _sync_particle_set(particles);

#if 0
    bft_printf("\n Particle set after sync\n");
    cs_lagr_particle_set_dump(particles);
#endif

    /*  assert(j == -1);  After a loop on particles, next_id of the last
        particle must not be defined */

    displacement_step_id++;

  } /* End of while (global displacement) */

  /* Deposition sub-model additional loop */

  if (lagr_model->deposition > 0) {

    for (cs_lnum_t i = 0; i < particles->n_particles; i++) {

      unsigned char *particle = particles->p_buffer + p_am->extents * i;

      cs_lnum_t *cur_neighbor_face_id
        = cs_lagr_particle_attr(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
      cs_real_t *cur_part_yplus
        = cs_lagr_particle_attr(particle, p_am, CS_LAGR_YPLUS);

      cs_lagr_test_wall_cell(particle, p_am, visc_length,
                             cur_part_yplus, cur_neighbor_face_id);

      /* Modification of MARKO pointer */
      if (*cur_part_yplus > 100.0)
        cs_lagr_particles_set_lnum(particles, i, CS_LAGR_MARKO_VALUE, CS_LAGR_COHERENCE_STRUCT_BULK);

      else {

        if (*cur_part_yplus <
            cs_lagr_particles_get_real(particles, i, CS_LAGR_INTERF)) {

          if (cs_lagr_particles_get_lnum(particles, i, CS_LAGR_MARKO_VALUE) < 0)
            cs_lagr_particles_set_lnum(particles,
                                       i,
                                       CS_LAGR_MARKO_VALUE,
                                       CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF);
          else
            cs_lagr_particles_set_lnum(particles,
                                       i,
                                       CS_LAGR_MARKO_VALUE,
                                       CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF);

        }
        else {

          if (cs_lagr_particles_get_lnum(particles, i, CS_LAGR_MARKO_VALUE) < 0)
            cs_lagr_particles_set_lnum(particles,
                                       i,
                                       CS_LAGR_MARKO_VALUE,
                                       CS_LAGR_COHERENCE_STRUCT_DEGEN_SWEEP);

          else if (cs_lagr_particles_get_lnum(particles, i, CS_LAGR_MARKO_VALUE)
              == CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF
              ||
              cs_lagr_particles_get_lnum(particles, i, CS_LAGR_MARKO_VALUE)
              == CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF)
            cs_lagr_particles_set_lnum(particles,
                                       i,
                                       CS_LAGR_MARKO_VALUE,
                                       CS_LAGR_COHERENCE_STRUCT_DEGEN_EJECTION);
        }

      }

    }
  }

  /* Internal deposition: additional loop */
  if (cs_glob_porous_model == 3) {
    cs_real_t *covered_surface = NULL;
    BFT_MALLOC(covered_surface, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);

    /* Initialization */
    for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells_with_ghosts; cell_id++)
      covered_surface[cell_id] = 0.;

    for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_i_faces ; face_id++) {
      /* Internal face flagged as internal deposition */
      if (cs_glob_lagr_internal_conditions->i_face_zone_id[face_id] >= 0) {
        for (cs_lnum_t j = 0; j < 3; j++)
          fvq->i_f_face_normal[3*face_id+j] = fvq->i_face_normal[3*face_id+j];
      }
    }

    if (lagr_model->deposition == 1) {

      for (cs_lnum_t ip = 0; ip < particles->n_particles; ip++) {

        cs_lnum_t cell_num
          = cs_lagr_particles_get_lnum(particles, ip, CS_LAGR_CELL_NUM);

        if (   cell_num >=0
            && (   cs_lagr_particles_get_lnum(particles, ip, CS_LAGR_DEPOSITION_FLAG)
              == CS_LAGR_PART_IMPOSED_MOTION)) {

          cs_lnum_t cell_id = cell_num - 1;

          covered_surface[cell_id] += cs_math_pi * 0.25
            * pow(cs_lagr_particles_get_real(particles, ip, CS_LAGR_DIAMETER),2.)
            * cs_lagr_particles_get_real(particles, ip, CS_LAGR_FOULING_INDEX)
            * cs_lagr_particles_get_real(particles, ip, CS_LAGR_STAT_WEIGHT);

          /* Loop over internal faces of the current particle faces
           * NB: usefull for resuspension, the last face_id is stored.
           * face_id is unique in many cases.
           * */
          for (cs_lnum_t i = _particle_track_builder->cell_face_idx[cell_id];
              i < _particle_track_builder->cell_face_idx[cell_id+1] ;
              i++ ) {

            cs_lnum_t face_num = _particle_track_builder->cell_face_lst[i];

            if (face_num > 0) {

              cs_lnum_t face_id = face_num - 1;

              /* Internal face flagged as internal deposition */
              if (cs_glob_lagr_internal_conditions->i_face_zone_id[face_id] >= 0)
                cs_lagr_particles_set_lnum(particles, ip, CS_LAGR_NEIGHBOR_FACE_ID, face_id);

            }
          }

        }
      }
    }

    /* Synchronization */
    if (mesh->halo != NULL)
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, covered_surface);


    /* Compute fluid section and clip it to 0 if negative */
    for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_i_faces; face_id++) {
      cs_real_t temp = 0.;

      /* Internal face flagged as internal deposition */
      if (cs_glob_lagr_internal_conditions->i_face_zone_id[face_id] >= 0) {
        cs_lnum_t cell_id1 = mesh->i_face_cells[face_id][0];;
        cs_lnum_t cell_id2 = mesh->i_face_cells[face_id][1];
        /* Remove from the particle area from fluid section */
        for (cs_lnum_t id = 0; id < 3; id++)
          fvq->i_f_face_normal[3*face_id + id] -=
            (covered_surface[cell_id1] + covered_surface[cell_id2])
            * fvq->i_face_normal[3*face_id + id]
            / fvq->i_face_surf[face_id];


        /* If S_fluid . S is negative, that means we removed too much surface
         * to fluid surface */
        for (cs_lnum_t j = 0; j < 3; j++)
          temp += fvq->i_f_face_normal[3*face_id+j] * fvq->i_face_normal[3*face_id+j];

        if (temp <= 0.) {
          for (cs_lnum_t j = 0; j < 3; j++)
            fvq->i_f_face_normal[3*face_id+j] = 0.;
        }
        fvq->i_f_face_surf[face_id] = cs_math_3_norm(fvq->i_f_face_normal + 3*face_id);
      }
    }
    /* Free memory */
    BFT_FREE(covered_surface);
  }

  _finalize_displacement(particles, part_b_mass_flux);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize Lagrangian module.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_tracking_finalize(void)
{
  if (cs_glob_lagr_particle_set == NULL)
    return;

  /* Destroy particle set */

  cs_lagr_particle_finalize();

  /* Destroy builder */
  _particle_track_builder = _destroy_track_builder(_particle_track_builder);

  /* Destroy internal condition structure*/

  cs_lagr_finalize_internal_cond();

  /* Destroy the structure dedicated to dlvo modeling */

  if (cs_glob_lagr_model->dlvo)
    cs_lagr_dlvo_finalize();

  /* Destroy the structure dedicated to clogging modeling */

  if (cs_glob_lagr_model->clogging)
    cs_lagr_clogging_finalize();

  /* Destroy the structure dedicated to roughness surface modeling */

  if (cs_glob_lagr_model->roughness)
    cs_lagr_roughness_finalize();

  /* Delete MPI_Datatypes */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)  _delete_particle_datatypes();
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine the number of the closest wall face from the particle
 *        as well as the corresponding wall normal distance (y_p^+)
 *
 * Used for the deposition model.
 *
 * \param[in]   particle     particle attributes for current time step
 * \param[in]   p_am         pointer to attributes map for current time step
 * \param[in]   visc_length  viscous layer thickness
 * \param[out]  yplus        associated yplus value
 * \param[out]  face_id      associated neighbor wall face, or -1
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_test_wall_cell(const void                     *particle,
                       const cs_lagr_attribute_map_t  *p_am,
                       const cs_real_t                 visc_length[],
                       cs_real_t                      *yplus,
                       cs_lnum_t                      *face_id)
{
  cs_lnum_t cell_num
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM);

  if (cell_num < 0) return;

  cs_lnum_t  *cell_b_face_idx = cs_glob_mesh_adjacencies->cell_b_faces_idx;
  cs_lnum_t  *cell_b_faces = cs_glob_mesh_adjacencies->cell_b_faces;
  cs_lnum_t cell_id = cell_num - 1;

  *yplus = 10000;
  *face_id = -1;

  cs_lnum_t  start = cell_b_face_idx[cell_id];
  cs_lnum_t  end =  cell_b_face_idx[cell_id + 1];

  for (cs_lnum_t i = start; i < end; i++) {
    cs_lnum_t f_id = cell_b_faces[i];

    assert(cs_glob_lagr_boundary_conditions != NULL);

    const char b_type = cs_glob_lagr_boundary_conditions->elt_type[f_id];

    if (   (b_type == CS_LAGR_DEPO1)
        || (b_type == CS_LAGR_DEPO2)
        || (b_type == CS_LAGR_DEPO_DLVO)) {

      cs_real_t x_face = cs_glob_lagr_b_u_normal[f_id][0];
      cs_real_t y_face = cs_glob_lagr_b_u_normal[f_id][1];
      cs_real_t z_face = cs_glob_lagr_b_u_normal[f_id][2];

      cs_real_t offset_face = cs_glob_lagr_b_u_normal[f_id][3];
      const cs_real_t  *particle_coord
        = cs_lagr_particle_attr_const(particle, p_am, CS_LAGR_COORDS);

      cs_real_t dist_norm = CS_ABS(  particle_coord[0] * x_face
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointers to cell face connectivity used in particle tracking.
 *
 * \param[out]  cell_face_idx  cell face index
 * \param[out]  cell_face_lst  cell face connectivity (signed 1-to-n based,
 *                             negative for boundary faces, positive for
 *                             interior faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_get_cell_face_connectivity(cs_lnum_t  **cell_face_idx,
                                   cs_lnum_t  **cell_face_lst)
{
  cs_lagr_track_builder_t  *builder = _particle_track_builder;

  if (builder != NULL) {
    *cell_face_idx = builder->cell_face_idx;
    *cell_face_lst = builder->cell_face_lst;
  }
  else {
    cs_lagr_particle_set_t  *particles = cs_glob_lagr_particle_set;
    if (particles != NULL) {
      _particle_track_builder
        =_init_track_builder(particles->n_particles_max,
                             particles->p_am->extents);
      builder = _particle_track_builder;
      *cell_face_idx = builder->cell_face_idx;
      *cell_face_lst = builder->cell_face_lst;
    }
    else {
      *cell_face_idx = NULL;
      *cell_face_lst = NULL;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
