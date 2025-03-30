/*============================================================================
 * Methods for particle localization
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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

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
#include "bft/bft_mem.h"

#include "fvm/fvm_periodicity.h"

#include "base/cs_base.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_physical_constants.h"
#include "mesh/cs_geom.h"
#include "base/cs_halo.h"
#include "base/cs_interface.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_porous_model.h"
#include "base/cs_random.h"
#include "base/cs_rotation.h"
#include "base/cs_search.h"
#include "base/cs_timer_stats.h"
#include "base/cs_turbomachinery.h"
#include "turb/cs_turbulence_model.h"

#include "base/cs_field.h"
#include "base/cs_field_pointer.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_car.h"
#include "lagr/cs_lagr_clogging.h"
#include "lagr/cs_lagr_coupling.h"
#include "lagr/cs_lagr_deposition_model.h"
#include "lagr/cs_lagr_dlvo.h"
#include "lagr/cs_lagr_orientation.h"
#include "lagr/cs_lagr_event.h"
#include "lagr/cs_lagr_particle.h"
#include "lagr/cs_lagr_porosity.h"
#include "lagr/cs_lagr_prototypes.h"
#include "lagr/cs_lagr_post.h"
#include "lagr/cs_lagr_resuspension.h"
#include "lagr/cs_lagr_roughness.h"
#include "lagr/cs_lagr_sde.h"
#include "lagr/cs_lagr_sde_model.h"
#include "lagr/cs_lagr_stat.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_tracking.h"

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

  cs_lnum_t  last_face_id;          /* last face id encountered
                                       (interior first, boundary next) */

  cs_lagr_tracking_state_t  state;  /* current state */

  int tracking_step_id;             /* 1 if cell_wise_integ == 0
                                         if cell_wise_integ == 1
                                       0 track deterministic ghost particle
                                       1 track true stochastic integ particle
                                       2 final location reached */

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
  cs_lnum_t  *dist_cell_id;   /* local cell id on distant ranks */
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

  cs_lagr_halo_t      *halo;   /* Lagrangian halo structure */

  cs_interface_set_t  *face_ifs;

} cs_lagr_track_builder_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Global variable for the current subroutines */

static  cs_lagr_track_builder_t  *_particle_track_builder = nullptr;

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

  CS_MALLOC(cs_type, tot_extents, cs_datatype_t);

  for (i = 0; i < tot_extents; i++)
    cs_type[i] = CS_CHAR;

  /* Map tracking info */

  size_t attr_start, attr_end;

  attr_start = offsetof(cs_lagr_tracking_info_t, start_coords);
  attr_end = attr_start + 3*sizeof(cs_real_t);
  for (i = attr_start; i < attr_end; i++)
    cs_type[i] = CS_REAL_TYPE;

  attr_start = offsetof(cs_lagr_tracking_info_t, last_face_id);
  attr_end = attr_start + sizeof(cs_lnum_t);
  for (i = attr_start; i < attr_end; i++)
    cs_type[i] = CS_LNUM_TYPE;

  attr_start = offsetof(cs_lagr_tracking_info_t, last_face_id);
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

  CS_MALLOC(blocklengths, count, int);
  CS_MALLOC(types, count, MPI_Datatype);
  CS_MALLOC(displacements, count, MPI_Aint);

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

  CS_FREE(displacements);
  CS_FREE(types);
  CS_FREE(blocklengths);
  CS_FREE(cs_type);

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
  cs_lnum_t  *cell_id = nullptr;
  cs_lagr_halo_t  *lagr_halo = nullptr;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;
  const cs_lnum_t  n_halo_cells = halo->n_elts[CS_HALO_EXTENDED];

  CS_MALLOC(lagr_halo, 1, cs_lagr_halo_t);

  assert(n_halo_cells == halo->index[2*halo->n_c_domains]);
  assert(n_halo_cells == mesh->n_ghost_cells);

  lagr_halo->extents = extents;
  lagr_halo->n_cells = n_halo_cells;

  /* Allocate buffers to enable the exchange between communicating ranks */

  CS_MALLOC(lagr_halo->send_shift, halo->n_c_domains, cs_lnum_t);
  CS_MALLOC(lagr_halo->send_count, halo->n_c_domains, cs_lnum_t);
  CS_MALLOC(lagr_halo->recv_shift, halo->n_c_domains, cs_lnum_t);
  CS_MALLOC(lagr_halo->recv_count, halo->n_c_domains, cs_lnum_t);

  lagr_halo->send_buf_size = CS_LAGR_MIN_COMM_BUF_SIZE;

  CS_MALLOC(lagr_halo->send_buf,
            lagr_halo->send_buf_size * extents,
            unsigned char);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    cs_lnum_t  request_size = 2 * halo->n_c_domains;

    CS_MALLOC(lagr_halo->request, request_size, MPI_Request);
    CS_MALLOC(lagr_halo->status,  request_size, MPI_Status);

  }
#endif

  /* Fill rank */

  CS_MALLOC(lagr_halo->rank, n_halo_cells, cs_lnum_t);

  for (rank = 0; rank < halo->n_c_domains; rank++) {

    for (i = halo->index[2*rank]; i < halo->index[2*rank+2]; i++)
      lagr_halo->rank[halo_cell_id++] = rank;

  }

  assert(halo_cell_id == n_halo_cells);

  /* Fill transform_id */

  CS_MALLOC(lagr_halo->transform_id, n_halo_cells, cs_lnum_t);

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

  /* Fill dist_cell_id */

  CS_MALLOC(lagr_halo->dist_cell_id, n_halo_cells, cs_lnum_t);

  CS_MALLOC(cell_id, mesh->n_cells_with_ghosts, cs_lnum_t);

  for (i = 0; i < mesh->n_cells_with_ghosts; i++)
    cell_id[i] = i;

  cs_halo_sync_num(halo, CS_HALO_EXTENDED, cell_id);

  for (i = 0; i < n_halo_cells; i++)
    lagr_halo->dist_cell_id[i] = cell_id[mesh->n_cells + i];

  /* Free memory */

  CS_FREE(cell_id);

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
  if (*halo != nullptr) {

    cs_lagr_halo_t *h = *halo;

    CS_FREE(h->rank);
    CS_FREE(h->transform_id);
    CS_FREE(h->dist_cell_id);

    CS_FREE(h->send_shift);
    CS_FREE(h->send_count);
    CS_FREE(h->recv_shift);
    CS_FREE(h->recv_count);

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {
      CS_FREE(h->request);
      CS_FREE(h->status);
    }
#endif

    CS_FREE(h->send_buf);
    CS_FREE(*halo);
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
    CS_REALLOC(lag_halo->send_buf, n_halo*tot_extents, unsigned char);
  }

  /* If decrease is allowed, do it progressively, and with a wide
     margin, so as to avoid re-increasing later if possible */

  else if (n_halo > n_send_particles*16) {
    n_halo /= 8;
    lag_halo->send_buf_size = n_halo;
    CS_REALLOC(lag_halo->send_buf, n_halo*tot_extents, unsigned char);
  }

  /* Otherwise, keep current size */
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

  if (n_particles_max == 0)
    return nullptr;

  cs_lagr_track_builder_t  *builder = nullptr;
  CS_MALLOC(builder, 1, cs_lagr_track_builder_t);

  /* Ensure a cell->face connectivity is defined */

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  if (ma->cell_i_faces == nullptr)
    cs_mesh_adjacencies_update_cell_i_faces();

  /* Define a cs_lagr_halo_t structure to deal with parallelism and
     periodicity */

  if (cs_glob_mesh->n_init_perio > 0 || cs_glob_n_ranks > 1)
    builder->halo = _create_lagr_halo(extents);
  else
    builder->halo = nullptr;

  /* Define an interface set on interior faces for keeping up-to-date
     the last_face_id value across ranks. Not used in serial mode */

  builder->face_ifs = nullptr;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    builder->face_ifs = cs_interface_set_create(mesh->n_i_faces,
                                                nullptr,
                                                mesh->global_i_face_num,
                                                nullptr,
                                                0,
                                                nullptr,
                                                nullptr,
                                                nullptr);

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
 *   nullptr
 *----------------------------------------------------------------------------*/

static cs_lagr_track_builder_t *
_destroy_track_builder(cs_lagr_track_builder_t  *builder)
{
  if (builder == nullptr)
    return builder;

  /* Destroy the cs_lagr_halo_t structure */

  _delete_lagr_halo(&(builder->halo));
  cs_interface_set_destroy(&(builder->face_ifs));

  /* Destroy the builder structure */

  CS_FREE(builder);

  return nullptr;
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
    = cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, attr_map, 1,
                                                 CS_LAGR_COORDS);
  cs_real_t *part_coord
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, attr_map,
                                               CS_LAGR_COORDS);

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
              _("Max number of local loops reached in particle displacement."));
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
 * Integrate SDEs associated to a particle
 *
 * parameters:
 *   particles       <-> pointer to particle set
 *   p_id            <-- particle index in set
 *   nor             <-- current step id (for 2nd order scheme)
 *   dt_part         <-- remaining time step associated to the particle
 *   taup            <-- dynamic characteristic time
 *   tlag            <-- fluid characteristic time
 *   piil            <-- term in integration of U-P SDEs
 *   force_p         <-- forces per mass unit on particles (m/s^2)
 *   bx              <-- turbulence characteristics
 *   tempct          <-- thermal characteristic time
 *   tsfext          --> info for return coupling source terms
 *   vislen          <-- nu/u* = y/y+
 *   beta            <-- proportional to the gradient of T_lag
 *   vagaus          <-- gaussian random variables
 *   br_gaus         <-- gaussian random variables
 *   cpgd1           <-> devolatilization term 1
 *   cpgd2           <-> devolatilization term 2
 *   cpght           <-> term for heterogeneous combustion (coal with thermal
 *                       return coupling)
 *   nresnew         <--
 *   loc_rebound     <-- save if a specific treatment occured during
 *                       this sub-iteration for cell_wise_integ == 1

 *----------------------------------------------------------------------------*/

static void
_integ_particle_quantities(cs_lagr_particle_set_t          *particles,
                           cs_lnum_t                        p_id,
                           int                              nor,
                           cs_real_t                        dt_part,
                           const cs_real_t                 *taup,
                           const cs_real_3_t               *tlag,
                           const cs_real_3_t               *piil,
                           const cs_real_3_t                force_p,
                           const cs_real_33_t              *bx,
                           cs_real_2_t                      tempct,
                           cs_real_t                       *tsfext,
                           const cs_real_t                  vislen[],
                           const cs_real_3_t                beta,
                           cs_real_3_t                     *vagaus,
                           cs_real_6_t                      br_gaus,
                           cs_real_t                       *cpgd1,
                           cs_real_t                       *cpgd2,
                           cs_real_t                       *cpght,
                           cs_lnum_t                       *nresnew,
                           bool                             loc_rebound)
{
  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;
  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();
  cs_real_t dtp = cs_glob_lagr_time_step->dtp;

  assert(dt_part >= 0.);

  /* test whether a correction step should occur when using second order scheme
   * no correction if a specific treatment as occured at a face during the
   * trajectography step */
  if (   nor == 2
      && (   (cs_glob_lagr_time_scheme->cell_wise_integ == 1 && loc_rebound)
          || (    cs_glob_lagr_time_scheme->cell_wise_integ == 0
               && cs_lagr_particles_get_lnum(particles, p_id,
                                             CS_LAGR_REBOUND_ID) == 0)))
    return;

  if (    !cs_lagr_particles_get_flag(particles, p_id, CS_LAGR_PART_FIXED)
      &&  (  cs_lagr_particles_get_real(particles, p_id, CS_LAGR_RESIDENCE_TIME)
                  > 1.e-12 * dtp
           || cs_glob_lagr_time_scheme->cell_wise_integ == 1)) {

    /* if imposed motion no correction of the velocities */
    if (!(   nor == 2
          && cs_lagr_particles_get_flag(particles, p_id,
                                        CS_LAGR_PART_IMPOSED_MOTION)))
      /* Integration of SDEs: position, fluid and particle velocity */
      cs_lagr_sde(p_id,
                  dt_part,
                  nor,
                  taup,
                  tlag,
                  piil,
                  bx,
                  tsfext,
                  force_p,
                  vislen,
                  beta,
                  vagaus,
                  br_gaus,
                  nresnew);

    /* Integration of SDEs for orientation of spheroids without inertia */
    int iprev;
    if (nor == 1)
      /* Use fields at previous time step    */
      iprev = 1;
    else
      /* Use fields at current time step     */
      iprev = 0;
    if (lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_STOC_MODEL) {

      cs_lagr_orientation_dyn_spheroids(p_id,
                                        iprev,
                                        dt_part,
                                        (const cs_real_33_t *)extra->grad_vel);
    }
    /* Integration of Jeffrey equations for ellipsoids */
    else if (lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_JEFFERY_MODEL) {
      cs_lagr_orientation_dyn_jeffery(p_id,
                                      iprev,
                                      dt_part,
                                      (const cs_real_33_t *)extra->grad_vel);
    }
  }

  /* Integration of SDE related to physical models */
  if (lagr_model->physical_model != CS_LAGR_PHYS_OFF)
    cs_lagr_sde_model(p_id, dt_part, nor, tempct, cpgd1, cpgd2, cpght);

  /* Integration of additional user variables
     -----------------------------------------*/

  if (cs_glob_lagr_model->n_user_variables > 0)
    cs_user_lagr_sde(dt_part, p_id, taup, tlag, tempct, nor);
}

/*----------------------------------------------------------------------------
 * Handle particles moving to internal deposition face
 *
 * parameters:
 *   particles                     <-- pointer to particle set
 *   p_id                          <-- particle id
 *   face_id                       <-- index of the treated face
 *   rel_disp_intersect            <-- used to compute the intersection of the
 *                                     trajectory and the face
 *   next_location                 <--> next position in the trajectography
 *   specific_face_interaction     <--> true if specific interaction at face
 *
 * returns:
 *   particle state
 *----------------------------------------------------------------------------*/

static cs_lagr_tracking_state_t
_internal_treatment(cs_lagr_particle_set_t  *particles,
                    cs_lnum_t                p_id,
                    cs_lnum_t                face_id,
                    double                   rel_disp_intersect,
                    cs_real_t                next_location[],
                    bool                    *specific_face_interaction)
{
  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  int n_phases = extra->n_phases;

  cs_lagr_tracking_state_t  particle_state = CS_LAGR_PART_TO_SYNC;

  cs_lagr_internal_condition_t *internal_conditions
    = cs_glob_lagr_internal_conditions;

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const double bc_epsilon = 1.e-6;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  unsigned char *particle = particles->p_buffer + p_am->extents * p_id;

  cs_real_t  disp[3], intersect_pt[3];
  const cs_mesh_t  *mesh = cs_glob_mesh;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  cs_real_t  *particle_coord
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_velocity
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY);
  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY_SEEN);

  cs_real_t particle_stat_weight
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
  cs_real_t particle_mass
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

  cs_lnum_t cell_id
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);

  /* In most cases, we just have to stop the integration at the face */

  if (   internal_conditions == nullptr
      || internal_conditions->i_face_zone_id[face_id] < 0) {
    if(cs_glob_lagr_time_scheme->cell_wise_integ > 0){

      cs_lnum_t  c_id1 = mesh->i_face_cells[face_id][0];
      cs_lnum_t  c_id2 = mesh->i_face_cells[face_id][1];

      for (int k = 0; k < 3; k++)
        disp[k] = next_location[k] - p_info->start_coords[k];

      for (int k = 0; k < 3; k++)
        intersect_pt[k] = disp[k] * rel_disp_intersect + p_info->start_coords[k];

      cs_real_t *cell_cen;
      if (cell_id == c_id1)
        cell_cen = fvq->cell_cen[c_id2];
      else
        cell_cen = fvq->cell_cen[c_id1];

      cs_real_3_t vect_cen;
      for (int k = 0; k < 3; k++){
        vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
        p_info->start_coords[k] =   intersect_pt[k] + bc_epsilon * vect_cen[k];
      }
    }
    return particle_state;
  }

  /* Now, we know we have some specific interactions to handle */

  *specific_face_interaction = true;
  /* FIXME GB set CS_LAGR_REBOUND_ID, 0 */

  assert(internal_conditions != nullptr);

  const cs_nreal_t *face_norm = fvq->i_face_u_normal[face_id];

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  CS_NO_WARN_IF_UNUSED(cell_vol);

  for (int k = 0; k < 3; k++)
    disp[k] = next_location[k] - p_info->start_coords[k];

  assert(! (fabs(disp[0]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[1]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[2]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15));

  for (int k = 0; k < 3; k++)
    intersect_pt[k] = disp[k] * rel_disp_intersect + p_info->start_coords[k];

  if (   internal_conditions->i_face_zone_id[face_id] == CS_LAGR_OUTLET
           || internal_conditions->i_face_zone_id[face_id] == CS_LAGR_INLET) {

    particle_state = CS_LAGR_PART_OUT;
    p_info->tracking_step_id = 2;

    if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                                -1.);

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

      cs_real_t *cell_cen = fvq->cell_cen[cell_id];
      cs_real_3_t vect_cen;
      for (int k = 0; k < 3; k++)
        vect_cen[k] = (cell_cen[k] - intersect_pt[k]);

      for (int k = 0; k < 3; k++) {
        particle_velocity[k] = 0.0;
      }
      /* Force the particle on the intersection but in the original cell */
      for (int k = 0; k < 3; k++) {
        next_location[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
        particle_coord[k] = next_location[k];
        for (int phase_id = 0; phase_id < n_phases; phase_id++)
          particle_velocity_seen[3 * phase_id + k] = 0.0;
      }

      /* The particle is not treated yet: the motion is now imposed */
      cs_lagr_particles_set_flag(particles, p_id,
                                 CS_LAGR_PART_IMPOSED_MOTION);

      /* Specific treatment in case of particle resuspension modeling */

      particle_state = CS_LAGR_PART_TREATED;
      p_info->tracking_step_id = 2;
      if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                                  -1.);

      particles->n_part_dep += 1;
      particles->weight_dep += particle_stat_weight;

    }
  }
  else if (internal_conditions->i_face_zone_id[face_id] == CS_LAGR_BC_USER){
    cs_lagr_user_internal_interaction(particles,
                                      p_id,
                                      face_id,
                                      face_norm,
                                      intersect_pt,
                                      rel_disp_intersect,
                                      &particle_state);

  /* FIXME: JBORD* (user-defined boundary condition) not yet implemented
   nor defined by a macro */
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Internal condition %d not recognized.\n"),
              internal_conditions->i_face_zone_id[face_id]);

  /* TODO internal face events ?... */

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Add event when particle rolls off an interior face
 *
 * parameters:
 *   particles  <-- pointer to particle set
 *   events     <-> events structure
 *   p_id       <-- particle id
 *   b_face_id  <-- associated boundary face id
 *----------------------------------------------------------------------------*/

static void
_roll_off_event(cs_lagr_particle_set_t    *particles,
                cs_lagr_event_set_t       *events,
                cs_lnum_t                  p_id,
                cs_lnum_t                  b_face_id)
{
  /* Get event id, flushing events if necessary */

  cs_lnum_t event_id = events->n_events;
  if (event_id >= events->n_events_max) {
    cs_lagr_stat_update_event(events,
                              CS_LAGR_STAT_GROUP_TRACKING_EVENT);
    events->n_events = 0;
    event_id = 0;
  }
  events->n_events += 1;

  /* Now set event values */

  cs_lagr_event_init_from_particle(events, particles, event_id, p_id);

  cs_lagr_events_set_lnum(events,
                          event_id,
                          CS_LAGR_E_FACE_ID,
                          b_face_id);

  cs_real_t *p_vel =
    cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, p_id,
                                              CS_LAGR_VELOCITY);
  cs_real_t *e_vel_post =
    cs_lagr_events_attr_get_ptr<cs_real_t>(events, event_id,
                                           CS_LAGR_E_VELOCITY);
  for (int k = 0; k < 3; k++)
    e_vel_post[k] = p_vel[k];

  cs_lnum_t *e_flag = cs_lagr_events_attr_get_ptr<cs_lnum_t>(events,
                                                             event_id,
                                                             CS_LAGR_E_FLAG);
  *e_flag = *e_flag | CS_EVENT_ROLL_OFF;
}

/*----------------------------------------------------------------------------
 * Handle particles moving to boundary
 *
 * parameters:
 *   particles          <-- pointer to particle set
 *   events             <-> events structure
 *   p_id               <-- particle id
 *   face_id            <-- boundary face id
 *   face_norm          <-- unit face (or face subdivision) normal
 *   rel_disp_intersect <-- relative distance (in [0, 1]) of the intersection
 *                          point with the face relative to the initial
 *                           trajectory segment
 *   b_z_id             <-- boundary zone id of the matching face
 *
 * returns:
 *   particle state
 *----------------------------------------------------------------------------*/

static cs_lagr_tracking_state_t
_boundary_treatment(cs_lagr_particle_set_t    *particles,
                    cs_lagr_event_set_t       *events,
                    cs_lnum_t                  p_id,
                    cs_lnum_t                  face_id,
                    cs_real_t                 *face_norm,
                    double                     rel_disp_intersect,
                    int                        b_z_id)
{
  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  int n_phases = extra->n_phases;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const double pi = cs_math_pi;
  const double bc_epsilon = 1.e-6;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  unsigned char *particle = particles->p_buffer + p_am->extents * p_id;

  cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_real_t  tmp;
  cs_real_t  disp[3], intersect_pt[3];

  int  event_flag = 0;

  cs_lagr_tracking_state_t  particle_state = CS_LAGR_PART_TO_SYNC;

  cs_lagr_zone_data_t  *bdy_conditions = cs_lagr_get_boundary_conditions();
  cs_lagr_boundary_interactions_t  *bdy_interactions
    = cs_glob_lagr_boundary_interactions;

  cs_real_t  energt = 0.;
  cs_lnum_t  contact_number = 0;
  cs_real_t  *surface_coverage = nullptr;
  cs_real_t* deposit_height_mean = nullptr;
  cs_real_t* deposit_height_var = nullptr;
  cs_real_t* deposit_diameter_sum = nullptr;

  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  cs_real_t  *particle_coord
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am, CS_LAGR_COORDS);
  cs_real_t  *particle_velocity
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY);
  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY_SEEN);

  cs_real_t particle_stat_weight
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
  cs_real_t particle_mass
    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const char b_type = cs_glob_lagr_boundary_conditions->elt_type[face_id];

  assert(bdy_conditions != nullptr);

  for (int k = 0; k < 3; k++)
    disp[k] = particle_coord[k] - p_info->start_coords[k];

  cs_real_t face_area  = fvq->b_face_surf[face_id];

  cs_lnum_t  cell_id
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  CS_NO_WARN_IF_UNUSED(cell_vol);

  assert(! (fabs(disp[0]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[1]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15 &&
            fabs(disp[2]/pow(cell_vol[cell_id],1.0/3.0)) < 1e-15));

  for (int k = 0; k < 3; k++)
    intersect_pt[k] = disp[k] * rel_disp_intersect + p_info->start_coords[k];

  /* Cancel events for symmetry boundary, as it is not a "true" boundary */

  if (b_type == CS_LAGR_SYM)
    events = nullptr;

  /* Generate event */

  int event_id = -1;

  if (events != nullptr) {

    event_id = events->n_events;
    if (event_id >= events->n_events_max) {
      /* flush events */
      cs_lagr_stat_update_event(events,
                                CS_LAGR_STAT_GROUP_TRACKING_EVENT);
      events->n_events = 0;
      event_id = 0;
    }

    cs_lagr_event_init_from_particle(events, particles, event_id, p_id);

    cs_lagr_events_set_lnum(events,
                            event_id,
                            CS_LAGR_E_FACE_ID,
                            face_id);

    cs_real_t *e_coords
      = cs_lagr_events_attr_get_ptr<cs_real_t>(events,
                                               event_id,
                                    (cs_lagr_event_attribute_t)CS_LAGR_COORDS);
    for (int k = 0; k < 3; k++)
      e_coords[k] = intersect_pt[k];

  }

  if (   b_type == CS_LAGR_OUTLET
      || b_type == CS_LAGR_INLET
      || b_type == CS_LAGR_DEPO1) {

    particle_state = CS_LAGR_PART_OUT;
    p_info->tracking_step_id = 2;
    if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                                -1.);

    if (b_type == CS_LAGR_DEPO1) {
      particles->n_part_dep += 1;
      particles->weight_dep += particle_stat_weight;
      cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_DEPOSITED);
      event_flag = event_flag | CS_EVENT_DEPOSITION;
    }
    else
      event_flag = event_flag | CS_EVENT_OUTFLOW;

    for (int k = 0; k < 3; k++)
      particle_coord[k] = intersect_pt[k];
  }

  else if (b_type == CS_LAGR_DEPO2) {

    p_info->tracking_step_id = 2;
    if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                                  -1.);
    cs_real_t *cell_cen = fvq->cell_cen[cell_id];
    cs_real_3_t vect_cen;
    for (int k = 0; k < 3; k++) {
      vect_cen[k] = (cell_cen[k] - intersect_pt[k]);
      particle_velocity[k] = 0.0;
      particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
    }

    particles->n_part_dep += 1;
    particles->weight_dep += particle_stat_weight;

    /* Specific treatment in case of particle resuspension modeling */

    cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_DEPOSITED);

    if (cs_glob_lagr_model->resuspension == 0) {

      for (int k = 0; k < 3; k++)
        for (int phase_id = 0; phase_id < n_phases; phase_id++)
          particle_velocity_seen[3 * phase_id + k] = 0.0;

      cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_FIXED);
      particle_state = CS_LAGR_PART_STUCK;
      if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                                  -1.);
      p_info->tracking_step_id = 2;

    }
    else {
      particle_state = CS_LAGR_PART_TREATED;
      if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                                  -1.);
      p_info->tracking_step_id = 2;
    }

    event_flag = event_flag | CS_EVENT_DEPOSITION;

  }

  else if (b_type == CS_LAGR_DEPO_DLVO) {

    cs_real_t particle_diameter
      = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

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

      surface_coverage = &bound_stat[  bdy_interactions->iscovc * n_b_faces
                                     + face_id];
      deposit_height_mean = &bound_stat[  bdy_interactions->ihdepm * n_b_faces
                                        + face_id];
      deposit_height_var = &bound_stat[  bdy_interactions->ihdepv * n_b_faces
                                       + face_id];
      deposit_diameter_sum = &bound_stat[bdy_interactions->ihsum * n_b_faces
                                         + face_id];

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

      p_info->tracking_step_id = 2;
      if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                  -1.);
      cs_real_t *cell_cen = fvq->cell_cen[cell_id];
      cs_real_3_t vect_cen;
      for (int k = 0; k < 3; k++)
        vect_cen[k] = (cell_cen[k] - intersect_pt[k]);

      /* The particle deposits*/
      if (!cs_glob_lagr_model->clogging && !cs_glob_lagr_model->resuspension) {
        cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_DEPOSITED);

        particles->n_part_dep += 1;
        particles->weight_dep += particle_stat_weight;

        cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_FIXED);
        particle_state = CS_LAGR_PART_STUCK;
        if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                    -1.);
        p_info->tracking_step_id = 2;
      }

      if (!cs_glob_lagr_model->clogging && cs_glob_lagr_model->resuspension > 0) {

        cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_DEPOSITED);

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
        if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_REMAINING_INTEG_TIME,
                                    -1.);
        p_info->tracking_step_id = 2;
      }

      if (cs_glob_lagr_model->clogging) {

        bound_stat[bdy_interactions->inclgt
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

          bound_stat[bdy_interactions->inclg
                     * n_b_faces + face_id] += particle_stat_weight;

          /* The particle is replaced towards the cell center
           * (and not using the normal vector face_norm)
           * to avoid problems where the replacement is out of the cell.
           * This is valid only for star-shaped cells !!! */
          for (int k = 0; k < 3; k++) {
            particle_coord[k] = intersect_pt[k] + bc_epsilon * vect_cen[k];
            particle_velocity[k] = 0.0;
            for (int phase_id = 0; phase_id < n_phases; phase_id++)
              particle_velocity_seen[3 * phase_id + k] = 0.0;
          }

          cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_DEPOSITED);
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID,
                                    face_id);

          particles->n_part_dep += 1;
          particles->weight_dep += particle_stat_weight;
          particle_state = CS_LAGR_PART_TREATED;
          if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_REMAINING_INTEG_TIME, -1.);
          p_info->tracking_step_id = 2;
        }
        else {

          cs_lnum_t i;
          cs_real_t random = -1;
          cs_real_t scov_cdf;
          void *cur_part = nullptr;

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

            cs_lnum_t cur_part_flag
              = cs_lagr_particle_get_lnum(cur_part, p_am,
                                          CS_LAGR_P_FLAG);

            cs_lnum_t cur_part_close_face_id
              = cs_lagr_particle_get_lnum(cur_part, p_am,
                                          CS_LAGR_NEIGHBOR_FACE_ID);

            cs_real_t cur_part_stat_weight
              = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_STAT_WEIGHT);

            cs_real_t cur_part_diameter
              = cs_lagr_particle_get_real(cur_part, p_am, CS_LAGR_DIAMETER);

            if (   (cur_part_flag & CS_LAGR_PART_DEPOSITED)
                && (cur_part_close_face_id == face_id)) {
              scov_cdf +=   (pi * cs_math_pow2(cur_part_diameter) / 4.)
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
                      _(" Error in %s: in the face number %ld \n"
                        "no deposited particle found to form a cluster \n"
                        "using the surface coverage %e (scov_cdf %e) \n"
                        "The particle used thus belongs to another face (%ld) \n"),
                      __func__,
                      (long)particle_close_face_id, *surface_coverage,
                      scov_cdf, (long)cur_part_close_face_id);
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
                                 * cur_part_stat_weight / (4.0*face_area), 2)
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
                                      pow(  cs_math_pow3(cur_part_diameter)
                                          + cs_math_pow3(particle_diameter)
                                            * particle_stat_weight
                                            / cur_part_stat_weight, 1./3.));

            cur_part_diameter = cs_lagr_particle_get_real(cur_part, p_am,
                                                          CS_LAGR_DIAMETER);

            *surface_coverage +=   (pi * pow(cur_part_diameter,2)/4.)
                                 * cur_part_stat_weight / face_area;

            cs_lagr_particle_set_real(cur_part, p_am, CS_LAGR_HEIGHT,
                                      cur_part_diameter);
          }

          cs_lagr_particle_set_real
            (cur_part, p_am, CS_LAGR_MASS,
             cur_part_mass + (  particle_mass
                              * particle_stat_weight
                              / cur_part_stat_weight));
          cs_lagr_particle_set_real
            (cur_part, p_am, CS_LAGR_CLUSTER_NB_PART,
             cur_part_cluster_nb_part + (  particle_cluster_nb_part
                                         * particle_stat_weight
                                         / cur_part_stat_weight));

          particle_state = CS_LAGR_PART_OUT;
          if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_REMAINING_INTEG_TIME, -1.);
          p_info->tracking_step_id = 2;
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

      event_flag = event_flag | CS_EVENT_DEPOSITION;
    }
    else {

      /*The particle does not deposit:
        It 'rebounds' on the energy barrier*/

      particle_state = CS_LAGR_PART_TO_SYNC;
      cs_lagr_particles_unset_flag(particles, p_id,
                                   CS_LAGR_PART_DEPOSITION_FLAGS);

      cs_real_t *cell_cen = fvq->cell_cen[cell_id];
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

      for (int phase_id = 0; phase_id < n_phases; phase_id++) {
        tmp = 2. * cs_math_3_dot_product(particle_velocity_seen + 3 * phase_id,
                                         face_norm);

        for (int k = 0; k < 3; k++)
          particle_velocity_seen[3 * phase_id + k] -= tmp * face_norm[k];
      }

      event_flag = event_flag | CS_EVENT_REBOUND;
    }
  }

  /* Anelastic rebound */
  else if (b_type == CS_LAGR_REBOUND ) {

    particle_state = CS_LAGR_PART_TO_SYNC;

    cs_real_t *cell_cen = fvq->cell_cen[cell_id];
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

    for (int phase_id = 0; phase_id < n_phases; phase_id++) {
      tmp = 2. * cs_math_3_dot_product(particle_velocity_seen + 3 * phase_id,
                                       face_norm);

      /* Wall function */
      if (extra_i[phase_id].cvar_rij != nullptr) {
        /* Reynolds stress tensor (current value) */
        cs_real_t *r_ij = &(extra_i[phase_id].cvar_rij->vals[0][6*cell_id]);
        /* Component Rnn = ni Rij nj */
        cs_real_t r_nn =
            cs_math_3_sym_33_3_dot_product(face_norm, r_ij, face_norm);

        /* Vector Rij.nj */
        cs_real_t r_in[3];
        cs_math_sym_33_3_product(r_ij, face_norm, r_in);

        /* Modify velocity seen */
        for (int k = 0; k < 3; k++)
          particle_velocity_seen[3 * phase_id + k] -= r_in[k] / r_nn * tmp ;

        /* Modify particle temperature seen */
        if (   cs_glob_lagr_model->physical_model != CS_LAGR_PHYS_OFF
            && extra->temperature != nullptr
            && extra->temperature_turbulent_flux != nullptr) {

          cs_real_t temperature_out =
            cs_lagr_particles_get_real(particles, p_id,
                                       CS_LAGR_TEMPERATURE_SEEN);
          /* Normal thermal turbulent flux */
          cs_real_t *_rit = &extra->temperature_turbulent_flux->val[3*cell_id];
          cs_real_t  r_tn = cs_math_3_dot_product(_rit, face_norm);
          cs_lagr_particles_set_real(particles, p_id, CS_LAGR_TEMPERATURE_SEEN,
                                     temperature_out -  r_tn / r_nn * tmp);
        }
      }

      /* else: for EVM u*^2 / r_nn = 1/sqrt(C0) as algebraic closure for SLM */
      else {
        cs_real_3_t r_in_ov_rnn;
        cs_real_3_t fluid_vel_tau;
        for (int k = 0; k < 3; k++)
          fluid_vel_tau[k] = extra_i[phase_id].vel->val[3*cell_id+k];

        cs_real_t fluid_vel_n =
          cs_math_3_dot_product(fluid_vel_tau, face_norm);
        for (int k = 0; k < 3; k++)
          fluid_vel_tau[k] = fluid_vel_tau[k] - fluid_vel_n * face_norm[k];

        cs_real_3_t dir_tau;
        cs_math_3_normalize(fluid_vel_tau, dir_tau);

        for (int k = 0; k < 3; k++) {
          r_in_ov_rnn[k] = face_norm[k] - dir_tau[k] / sqrt(cs_turb_crij_c0);
          particle_velocity_seen[3 * phase_id + k] -= r_in_ov_rnn[k] * tmp;
        }

        if (   cs_glob_lagr_model->physical_model != CS_LAGR_PHYS_OFF
            && extra->temperature != nullptr
            && extra->temperature_variance != nullptr
            && (extra->tstar != nullptr || extra->grad_tempf != nullptr)
            && extra_i[phase_id].cvar_k != nullptr) {
          cs_real_t temperature_out =
            cs_lagr_particles_get_real(particles, p_id,
                                       CS_LAGR_TEMPERATURE_SEEN);

          /* Modify temperature seen with algebraic solution from SLM--IEM model
           * in the neutral limit (=close to walls, could be poor further away)*/
          cs_real_t rnt_ov_rnn =
            -sqrt(cs_turb_crij_ct * (1.5 * cs_turb_crij_c0 + 1.) /
              (cs_turb_crij_c0 * (cs_turb_crij_ct + 0.75 * cs_turb_crij_c0 + 1.))
                * extra->temperature_variance->val[cell_id]
                / extra_i[phase_id].cvar_k->val[cell_id]);

          /* Estimate the direction of the flux */
          if (extra->grad_tempf != nullptr ) {
            cs_real_t grad_n = cs_math_3_dot_product(extra->grad_tempf[cell_id],
                                                     face_norm);
            if (cs::abs(grad_n) < cs_math_epzero)
              rnt_ov_rnn = 0;
            else if (grad_n < 0.)
              rnt_ov_rnn *= -1.;
          }
          else if (extra->tstar != nullptr) { //define only at walls
            if (cs::abs(extra->tstar->val[face_id]) < cs_math_epzero)
              rnt_ov_rnn = 0;
            else if (extra->tstar->val[face_id] < 0.)
              rnt_ov_rnn *= -1.;
          }
          else
            rnt_ov_rnn = 0;

          cs_lagr_particles_set_real(particles, p_id, CS_LAGR_TEMPERATURE_SEEN,
                                     temperature_out -  rnt_ov_rnn * tmp);
        }
        else if (   cs_glob_lagr_model->physical_model != CS_LAGR_PHYS_OFF
                 && extra->temperature != nullptr
                 && extra->tstar != nullptr
                 && extra->ustar != nullptr) {
          if (extra->ustar->val[face_id] > cs_math_epzero) {
            cs_real_t temperature_out =
              cs_lagr_particles_get_real(particles, p_id,
                                         CS_LAGR_TEMPERATURE_SEEN);
            cs_real_t rnt_ov_rnn = - extra->tstar->val[face_id] /
              (sqrt(cs_turb_crij_c0) * extra->ustar->val[face_id]);
            cs_lagr_particles_set_real(particles, p_id, CS_LAGR_TEMPERATURE_SEEN,
                                       temperature_out -  rnt_ov_rnn * tmp);
          }
        }
      }
    }

    event_flag = event_flag | CS_EVENT_REBOUND;

  }

  /* Elastic rebound */
  else if (b_type == CS_LAGR_SYM){

    particle_state = CS_LAGR_PART_TO_SYNC;

    cs_real_t *cell_cen = fvq->cell_cen[cell_id];
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

    for (int phase_id = 0; phase_id < n_phases; phase_id++) {
      tmp = 2. * cs_math_3_dot_product(particle_velocity_seen + 3 * phase_id,
                                       face_norm);

      for (int k = 0; k < 3; k++)
        particle_velocity_seen[3 * phase_id + k] -= tmp * face_norm[k];
    }

    event_flag = event_flag | CS_EVENT_REBOUND;
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
      = cs_lagr_particle_attr_get_const_ptr<cs_real_t>(particle, p_am,
                                                       CS_LAGR_TEMPERATURE);

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

        event_flag = event_flag | CS_EVENT_FOULING;

        particle_state = CS_LAGR_PART_OUT;
        if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
          cs_lagr_particle_set_real(particle, p_am,
                                    CS_LAGR_REMAINING_INTEG_TIME, -1.);
        p_info->tracking_step_id = 2;

        /* Recording for log/lagrangian.log */
        particles->n_part_fou += 1;
        particles->weight_fou += particle_stat_weight;

        /* Recording for statistics */
        /* FIXME: For post-processing by trajectory purpose */

        for (int k = 0; k < 3; k++) {
          particle_coord[k] = intersect_pt[k];
          particle_velocity[k] = 0.0;
          for (int phase_id = 0; phase_id < n_phases; phase_id++)
            particle_velocity_seen[3 * phase_id + k] = 0.0;
        }
      }
    }

    /*--> if there is no fouling, then it is an elastic rebound */
    if (particle_state == CS_LAGR_PART_TO_SYNC) {

      cs_real_t *cell_cen = fvq->cell_cen[cell_id];
      cs_real_t vect_cen[3];
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

      for (int phase_id = 0; phase_id < n_phases; phase_id++) {
        tmp = 2. * cs_math_3_dot_product(particle_velocity_seen + 3 * phase_id,
                                         face_norm);

        /*TODO set anelastic rebound */
        for (int k = 0; k < 3; k++) {
          particle_velocity_seen[3 * phase_id + k] -= tmp * face_norm[k];
          // particle_velocity_seen[k] = 0.0; //FIXME
        }
      }

      event_flag = event_flag | CS_EVENT_REBOUND;
    }

  }

  else if (b_type == CS_LAGR_BC_USER)
    cs_lagr_user_boundary_interaction(particles,
                                      p_id,
                                      face_id,
                                      face_norm,
                                      intersect_pt,
                                      rel_disp_intersect,
                                      b_z_id,
                                      &event_flag,
                                      &particle_state);

  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Boundary condition %d not recognized.\n"),
              b_type);

  /* Ensure some fields are updated */

  if (p_am->size[CS_LAGR_NEIGHBOR_FACE_ID] > 0) {

    int depo_flag
      = cs_lagr_particles_get_flag(particles, p_id,
                                   CS_LAGR_PART_DEPOSITED | CS_LAGR_PART_ROLLING);
    if (depo_flag) {
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID,
                                face_id);

      if (depo_flag & CS_LAGR_PART_ROLLING)
        event_flag = event_flag | CS_EVENT_ROLL_ON;
    }
  }

  if (events != nullptr) {
    cs_lnum_t *e_flag =
      cs_lagr_events_attr_get_ptr<cs_lnum_t>(events, event_id, CS_LAGR_E_FLAG);

    cs_real_t *e_vel_post =
      cs_lagr_events_attr_get_ptr<cs_real_t>(events, event_id,
                                             CS_LAGR_E_VELOCITY);

    *e_flag = *e_flag | event_flag;

    for (int k = 0; k < 3; k++)
      e_vel_post[k] = particle_velocity[k];

    events->n_events += 1;
  }

  /* Update per-zone flow rate measure for exiting particles */

  if (particle_state == CS_LAGR_PART_OUT) {
    int n_stats = cs_glob_lagr_model->n_stat_classes + 1;

    cs_real_t fr =   particle_stat_weight
                   * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

    bdy_conditions->particle_flow_rate[b_z_id*n_stats] -= fr;

    if (n_stats > 1) {
      int class_id
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_STAT_CLASS);
      if (class_id > 0 && class_id < n_stats)
        bdy_conditions->particle_flow_rate[  b_z_id*n_stats
                                           + class_id] -= fr;
    }
  }

  if  (   b_type == CS_LAGR_DEPO1
       || b_type == CS_LAGR_DEPO2
       || b_type == CS_LAGR_DEPO_DLVO
       || b_type == CS_LAGR_REBOUND
       || b_type == CS_LAGR_FOULING) {

    /* Number of particle-boundary interactions  */
    if (bdy_interactions->has_part_impact_nbr > 0)
      bound_stat[bdy_interactions->inbr * n_b_faces + face_id]
        += particle_stat_weight;

  }

  return particle_state;
}

/*----------------------------------------------------------------------------
 * Move a particle as far as possible while remaining on a given rank.
 *
 * parameters:
 *   particles                    <-> pointer to particle set
 *   events                   <-> events structure
 *   p_id                     <-- particle id
 *   failsafe_mode            <-- with (0) / without (1) failure capability
 *   b_face_zone_id           <-- boundary face zone id
 *   visc_length              <-- viscous layer thickness
 *   resol_sde                <--  true  the sdes will be resolved
 *                                 false only the trajectography is done
 *   nresnew
 * returns:
 *   a state associated to the status of the particle (treated, to be deleted,
 *   to be synchonised)
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_local_propagation(cs_lagr_particle_set_t         *particles,
                   cs_lagr_event_set_t            *events,
                   cs_lnum_t                       p_id,
                   int                             failsafe_mode,
                   const int                       b_face_zone_id[],
                   const cs_real_t                 visc_length[],
                   const bool                      resol_sde,
                   cs_lnum_t                      *nresnew)
{

  cs_real_t fluid_vel[3] = {0, 0, 0};
  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  int n_phases = extra->n_phases;

  cs_lagr_tracking_state_t  particle_state = CS_LAGR_PART_TO_SYNC;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;

  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *)fvq->b_face_cog;

  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  unsigned char *particle = particles->p_buffer + p_am->extents * p_id;
  cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

  /* integrated location */
  cs_real_t  *particle_coord
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am, CS_LAGR_COORDS);

  /* starting point for the tracking*/
  cs_real_t  *prev_location = p_info->start_coords;

  /* location tracked  may be different of the integrated position */
  cs_real_3_t  next_location;

  /* Usefull to avoid oscillation at face */
  cs_lnum_t save_old_cell_id = -1;
  cs_real_t dt_incremented_in_subiter;

  int cell_wise_integ = cs_glob_lagr_time_scheme->cell_wise_integ;

  /* Useful for cell_wise_integ == 1
   * when the deterministic virtual partner bounces back on a face */
  cs_real_t sum_t_intersect = 0.;
  bool save_specific_face_interaction = false;

 /* trajectory step for (2nd order scheme)
  * 1: 1st order scheme or prediction step for second order moment
  * 2: correction step for second order moment */
  int nor = 1;

  cs_real_t  *particle_velocity_seen
    = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY_SEEN);

  const cs_real_3_t *vtx_coord
    = (const cs_real_3_t *)(cs_glob_mesh->vtx_coord);

  /* Neighbor face id (allow test even without attribute */
  cs_lnum_t  null_face_id = -1;
  cs_lnum_t  *neighbor_face_id = &null_face_id;
  if (p_am->size[CS_LAGR_NEIGHBOR_FACE_ID] > 0)
    neighbor_face_id
      = cs_lagr_particle_attr_get_ptr<cs_lnum_t>(particle, p_am,
                                                 CS_LAGR_NEIGHBOR_FACE_ID);

  /* Particle y+  (allow tes even without attribute */
  cs_real_t  null_yplus = 0.;
  cs_real_t  *particle_yplus = &null_yplus;
  if (p_am->size[CS_LAGR_YPLUS] > 0)
    particle_yplus
      = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am, CS_LAGR_YPLUS);

  /* particle properties */
  cs_real_2_t tempct;
  cs_real_t cpgd1, cpgd2, cpght, tsfext;
  cs_real_3_t beta = { 0., 0., 0.};
  cs_real_3_t force_p;
  cs_real_6_t br_gaus;
  cs_real_t *taup = nullptr;
  cs_real_3_t *tlag = nullptr;
  cs_real_3_t *piil = nullptr;
  cs_real_3_t *vagaus = nullptr;
  cs_real_3_t *null_vagaus = nullptr;
  cs_real_33_t *bx = nullptr;
  if (resol_sde && cell_wise_integ == 1) {
    CS_MALLOC(taup, n_phases, cs_real_t);
    CS_MALLOC(tlag, n_phases, cs_real_3_t);
    CS_MALLOC(piil, n_phases, cs_real_3_t);
    CS_MALLOC(vagaus, n_phases + 2, cs_real_3_t);
    CS_MALLOC(null_vagaus, n_phases + 2, cs_real_3_t);
    CS_MALLOC(bx, n_phases, cs_real_33_t);
    for (int phase_id = 0; phase_id < n_phases + 2; phase_id++) {
      for (int i = 0; i < 3; i++)
        null_vagaus[phase_id][i] = 0.;
    }
  }

 /* limit number of loop to avoid infinite loop */
  cs_lnum_t max_propagation_loops
    = cs_glob_lagr_time_scheme->max_track_propagation_loops;
  /* particle displaced at init can be moved from any cell to any cell*/
  if (!resol_sde)
    max_propagation_loops = cs::max(mesh->n_cells, max_propagation_loops);

  /* local loops on the cell within a given proc
   * if either the deterministic virtual partiner or the stochastic particle
   * is blocked on face_treatment the particle is blocked *
   * current tracking step associated to the particle
   *  0 : ghost trajectory associated to deterministic virtual partner
   *  1 : trajectory step associated to the integrated stochastic particle
   *  2 : final location reached
   */
  for (int n_loops_local = 0;
      particle_state == CS_LAGR_PART_TO_SYNC;
      n_loops_local++) {

    cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                   CS_LAGR_CELL_ID);
    cs_lagr_particle_set_lnum_n(particle, p_am, 1, CS_LAGR_CELL_ID, cell_id);

    assert(cell_id < mesh->n_cells);
    assert(cell_id > -1);

    /* Manage error */
    if (n_loops_local > max_propagation_loops) {

      _manage_error(failsafe_mode,
                    particle,
                    p_am,
                    CS_LAGR_TRACKING_ERR_MAX_LOOPS);

      particle_state = CS_LAGR_PART_TREATED;

      cs_real_t *cell_cen = fvq->cell_cen[cell_id];

      for (int k = 0; k < 3; k++) {
        prev_location[k] = cell_cen[k];
        particle_coord[k] = cell_cen[k];
      }

      cs_lnum_t n_rep
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_TR_REPOSITION);
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_TR_REPOSITION, n_rep+1);

      return particle_state;
    }

    cs_real_t dt_part = -1.;
    if (p_info->tracking_step_id  == 0 && sum_t_intersect < cs_math_epzero) {
      /* Integration of deterministic virtual partner for cell_wise_integ
       * No new integration if the last face crossed this iteration is
       * is boundary condition resulting on a rebound*/
      dt_part =   cs_glob_lagr_time_step->dtp
                * cs_lagr_particle_get_real(particle,p_am,
                                            CS_LAGR_REMAINING_INTEG_TIME);
      for (int phase_id = 0; phase_id < n_phases; phase_id ++)
        cs_lagr_car(1, /* iprev = 1: Use fields at previous time step    */
                    phase_id,
                    p_id,
                    nor,
                    dt_part,
                    &taup[phase_id],
                    tlag[phase_id],
                    piil[phase_id],
                    bx[phase_id],
                    tempct,
                    beta,
                    vagaus,
                    br_gaus);

      cs_lagr_get_force_p(dt_part, p_id, taup, tlag, piil, bx,
                          tsfext, vagaus, force_p);

      cs_real_6_t null_brgaus = {0. ,0. ,0. ,0. , 0. , 0.};

      /* save the integrated fields */
      cs_sde_vels_pos_1_st_order_time_integ(p_id,
                                            dt_part,
                                            nor,
                                            taup,
                                            tlag,
                                            piil,
                                            bx,
                                            null_vagaus,
                                            null_brgaus,
                                            force_p,
                                            beta);
      dt_incremented_in_subiter = 0.;

    }
    /* track the integrated position of the stochastic particle */
    for(int k = 0 ; k < 3 ; k++)
      next_location[k] = particle_coord[k];

    /* Dimension less test: no movement ? */
    const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
    cs_real_t inv_ref_length = 1./pow(cell_vol[cell_id], 1./3.);
    cs_real_t  disp[3];
    for (int k = 0; k < 3; k++)
      disp[k] = next_location[k] - prev_location[k];
    if (   fabs(disp[0] * inv_ref_length) < 1e-15
        && fabs(disp[1] * inv_ref_length) < 1e-15
        && fabs(disp[2] * inv_ref_length) < 1e-15) {


      if (p_info->tracking_step_id == 0) {
        /* deterministic virtual partner tracked */
        dt_part =   cs_glob_lagr_time_step->dtp
                  * cs_lagr_particle_get_real(particle,p_am,
                                              CS_LAGR_REMAINING_INTEG_TIME);
        _integ_particle_quantities(particles,
                                   p_id,
                                   nor, // 1
                                   dt_part,
                                   taup,
                                   tlag,
                                   piil,
                                   force_p,
                                   bx,
                                   tempct,
                                   &tsfext,
                                   visc_length,
                                   beta,
                                   vagaus,
                                   br_gaus,
                                   &cpgd1,
                                   &cpgd2,
                                   &cpght,
                                   nresnew,
                                   save_specific_face_interaction);

        cs_lagr_particles_set_real(particles, p_id,
                                   CS_LAGR_REMAINING_INTEG_TIME, -1);
        p_info->tracking_step_id = 1; /* now track stochastic particle */
      }
      else if (p_info->tracking_step_id == 1) {
        /* stochastic particle is currently tracked */
        if (save_specific_face_interaction) {
          /* Restart new deterministic virtual partner trajectory now that
           * local effect of specific treatment are seen */
          for (int i = 0; i < 3; i++)
            prev_location[i] = particle_coord[i];
          save_specific_face_interaction = false;
          p_info->tracking_step_id = 0;
          /* Update previous state */
          cs_lagr_particles_current_to_previous(particles, p_id);
        }

        else { /* stochastic particle final location reached */
          particle_state = CS_LAGR_PART_TREATED;
          p_info->tracking_step_id = 2;
        }
      }

      continue;
    }
    /* Treatment for depositing particles */

    if (lagr_model->deposition > 0 && *particle_yplus < 0.) {

      cs_lagr_test_wall_cell(particle, p_am, visc_length,
                             particle_yplus, neighbor_face_id);

      if (*particle_yplus < 100.) {

        for (cs_lnum_t phase_id = 0; phase_id < n_phases; phase_id++) {
          for (int i = 0; i < 3; i++)
            fluid_vel[i] = extra_i[phase_id].vel->val[3 * cell_id + i];
          cs_real_t fluid_vel_proj[3];

          /* normal vector coordinates */
          const cs_nreal_t *normal = b_face_u_normal[*neighbor_face_id];

          /* (V . n) * n  */
          cs_real_t v_dot_n =
            cs_math_3_dot_product(particle_velocity_seen + 3 * phase_id,
                                  normal);

          /* tangential projection to the wall:
           * (Id -n (x) n) U */
          cs_math_3_orthogonal_projection(normal, fluid_vel, fluid_vel_proj);

          /* Update of the flow seen velocity */

          for (int i = 0; i < 3; i++)
            particle_velocity_seen[3 * phase_id + i] = v_dot_n * normal[i]
                                                     + fluid_vel_proj[i];
        }
      }
    }

    /* Loop on faces connected to the current cell */

    bool restart = false;

    reloop_cen:;

    cs_lnum_t exit_face = -1; /* id in [interior faces + boundary faces] */

    double adist_min = 2.;
    /* exit_face which originated the new computing */
    cs_real_t t_intersect = 1.;
    cs_real_t rel_disp_intersect = -1;
    cs_real_3_t face_norm={0,0,0};

    int n_in = 0;
    int n_out = 0;

    /* Loop on faces to see if the particle trajectory crosses it*/

    const cs_lnum_t n_cell_i_faces =   ma->cell_cells_idx[cell_id+1]
                                     - ma->cell_cells_idx[cell_id];
    const cs_lnum_t n_cell_b_faces =   ma->cell_b_faces_idx[cell_id+1]
                                     - ma->cell_b_faces_idx[cell_id];

    cs_lnum_t n_cell_faces = n_cell_i_faces + n_cell_b_faces;

    if (ma->cell_hb_faces_idx != nullptr) {
      n_cell_faces +=   ma->cell_hb_faces_idx[cell_id+1]
                      - ma->cell_hb_faces_idx[cell_id];
    }

    for (cs_lnum_t i = 0;
         i < n_cell_faces && particle_state == CS_LAGR_PART_TO_SYNC;
         i++) {

      cs_lnum_t face_id, vtx_start, vtx_end, n_vertices;
      const cs_lnum_t *face_connect;
      const cs_real_t *face_cog;

      /* Outward normal: always well oriented for external faces, depend on the
       * connectivity for internal faces */
      int reorient_face = 1;

      if (i < n_cell_i_faces) { /* Interior face */

        face_id = ma->cell_i_faces[ma->cell_cells_idx[cell_id] + i];

        if (cell_id == mesh->i_face_cells[face_id][1])
          reorient_face = -1;
        vtx_start = mesh->i_face_vtx_idx[face_id];
        vtx_end = mesh->i_face_vtx_idx[face_id+1];
        n_vertices = vtx_end - vtx_start;

        face_connect = mesh->i_face_vtx_lst + vtx_start;
        face_cog = i_face_cog[face_id];
      }
      else { /* Boundary faces */

        cs_lnum_t j = i - n_cell_i_faces;
        if (j < n_cell_b_faces)
          face_id = ma->cell_b_faces[ma->cell_b_faces_idx[cell_id] + j];

        else {
          assert(ma->cell_hb_faces_idx != nullptr);
          j -= n_cell_b_faces;
          face_id = ma->cell_hb_faces[ma->cell_hb_faces_idx[cell_id] + j];
        }

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
                                                (cs_real_t*)next_location,
                                                n_crossings,
                                                face_norm);

      n_in += n_crossings[0];
      n_out += n_crossings[1];

      /* Store the nearest intesection from the O point...*/
      if (t < adist_min && t >= 0) {
        exit_face = face_id;
        if (i >= n_cell_i_faces)
          exit_face += mesh->n_i_faces;
        rel_disp_intersect = t;
        adist_min = rel_disp_intersect;
      }
    }

    /* We test here if the particle is truly within the current cell
     * (meaning n_in = n_out > 0 )
     * If there is a problem (pb with particle strictly // and on the face ?),
     * the particle initial position is replaced at the cell center
     * and we continue the trajectory analysis. */

    bool test_in = (n_in == 0 && n_out == 0);

    if ((n_in != n_out || test_in)
        && (particle_state == CS_LAGR_PART_TO_SYNC)) {
      cs_real_t *cell_cen = fvq->cell_cen[cell_id];

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
        particle_state = CS_LAGR_PART_TREATED;
        if (cell_wise_integ)
          cs_lagr_particle_set_real(particle,p_am,CS_LAGR_REMAINING_INTEG_TIME,
                                                  -1.);
        return particle_state;
      }
      goto reloop_cen;
    }

    bool specific_face_interaction = false;
    int prev_tracking_step_id = p_info->tracking_step_id;

    if (lagr_model->deposition && exit_face >= 0) {
      cs_lnum_t b_face_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                      CS_LAGR_NEIGHBOR_FACE_ID);
      if (b_face_id > -1) {
        if (cs_lagr_particles_get_flag(particles, p_id, CS_LAGR_PART_ROLLING))
          _roll_off_event(particles, events, p_id, b_face_id);
      }
    }

    cs_real_t remain_time = 1.;
    t_intersect = cs::min(cs::abs(rel_disp_intersect), 1.);
    /* may be use for cell-wise integration when rebound at boundary occurs */
    sum_t_intersect += t_intersect * (1. - sum_t_intersect);

    if (p_info->tracking_step_id == 0) {
      /* we have tracked the deterministic virtual partner, the remaining time
       * is obtained thank to the linear interp. based on the previous state */
      remain_time = cs_lagr_particles_get_real(particles, p_id,
                                               CS_LAGR_REMAINING_INTEG_TIME);
      dt_part = remain_time * sum_t_intersect * cs_glob_lagr_time_step->dtp;
    }

    if (exit_face >= 0 && exit_face < mesh->n_i_faces) { /* Particle moves to the
                                                           neighbor cell */

      cs_lnum_t face_id = exit_face;

      cs_lnum_t  c_id1 = mesh->i_face_cells[face_id][0];
      cs_lnum_t  c_id2 = mesh->i_face_cells[face_id][1];

      p_info->last_face_id = exit_face;

      /* Neighbor face id needs update */

      if (p_am->size[CS_LAGR_NEIGHBOR_FACE_ID] > 0)
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID, -1);

      /* Deposition on internal faces? */

      /* particle / internal condition interaction
         1 - modify particle cell_num: 0 or boundary_cell_num
         2 -

         O -->  *         *  <-- D
         \       /
         \     /
         \   /
         \ /
         ------------------      boundary condition
         K
      */
      particle_state
        = _internal_treatment(particles,
                              p_id,
                              face_id,
                              rel_disp_intersect,
                              next_location,
                              &specific_face_interaction);

      if (particle_state == CS_LAGR_PART_TO_SYNC) {

        if (cell_id == c_id1)
          cell_id = c_id2;
        else
          cell_id = c_id1;

        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_CELL_ID, cell_id);

        /* Particle changes rank */

        if (cell_id >= mesh->n_cells) {

          particle_state = CS_LAGR_PART_TO_SYNC_NEXT;

          /* Specific treatment for the particle deposition model:
           * mark it so that it will be treated in the next iteration */

          *particle_yplus = - 1.;

        } /* end of case where particle changes rank */

        /* Update values for depositing particles */

        else if (lagr_model->deposition > 0) {

          cs_lagr_test_wall_cell(particle, p_am, visc_length,
                                 particle_yplus, neighbor_face_id);

          if (*particle_yplus < 100.) {

            for (cs_lnum_t phase_id = 0; phase_id < n_phases; phase_id++) {
              for (cs_lnum_t i = 0; i < 3; i++)
                fluid_vel[i] = extra_i[phase_id].vel->val[3 * cell_id + i];

              cs_real_t fluid_vel_proj[3];

              /* normal vector coordinates */
              const cs_nreal_t *normal = b_face_u_normal[*neighbor_face_id];

              /* (V . n) * n  */
              cs_real_t v_dot_n =
                cs_math_3_dot_product(particle_velocity_seen + 3 * phase_id,
                                      normal);

              /* tangential projection to the wall:
               * (Id -n (x) n) U */
              cs_math_3_orthogonal_projection(normal, fluid_vel, fluid_vel_proj);

              /* Update of the flow seen velocity */

              for (cs_lnum_t i = 0; i < 3; i++)
                particle_velocity_seen[3 * phase_id + i] = v_dot_n * normal[i]
                                                         + fluid_vel_proj[i];
            }
          }

        } /* end of case for deposition model */
      }
    }
    else if (exit_face >=mesh->n_i_faces) {
      /* Particle moves to the boundary through the current face "face_num"
         particle / boundary condition interaction
         1 - modify particle cell_num : 0 or boundary_cell_num
         2 -

         P -->  *         *  <-- Q
         \       /
         \     /
         \   /
         \ /
         ------------------      boundary condition
         K
      */

      cs_lnum_t face_id = exit_face - mesh->n_i_faces;

      particle_state
        = _boundary_treatment(particles,
                              events,
                              p_id,
                              face_id,
                              face_norm,
                              rel_disp_intersect,
                              b_face_zone_id[face_id]);

      if (cs_glob_lagr_time_scheme->t_order == 2)
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_REBOUND_ID, 0);

      specific_face_interaction = true;

      /* update  next position*/
      for(int k = 0; k < 3 ; k++){
        next_location[k] = particle_coord[k];
      }
      p_info->last_face_id = face_id;

    } /* end for test on exit_face */
    else if (exit_face < 0) {
      /* no face crossed with trajectography location corresponding
       * to the integrated position */
      if (p_info->tracking_step_id == 1) {
        /* stochastic particle is currently tracked */
        if (save_specific_face_interaction) {
          /* Restart new deterministic virtual partner trajectory now that
           * local effect of specific treatment are seen */
          for (int i = 0; i < 3; i++)
            prev_location[i] = particle_coord[i];
          save_specific_face_interaction = false;
          p_info->tracking_step_id = 0;
          sum_t_intersect = 0.;
          /* Update previous state */
          cs_lagr_particles_current_to_previous(particles, p_id);
        }

        else { /* stochastic particle final location reached */
          particle_state = CS_LAGR_PART_TREATED;
          p_info->tracking_step_id = 2;
        }
      }
      else if (p_info->tracking_step_id == 0) {
        /* the deterministic virtual partner remains in the cell the new
         * location to track is the final stochastic particle position*/
        for (int i = 0; i < 3; i++)
          next_location[i] = particle_coord[i];

        p_info->tracking_step_id = 1;
      }
    }

    if (prev_tracking_step_id == 0) {
      assert(dt_part > 0);

      /* FIXME improve condition */
      /* Made to avoid spurious oscillation around specific points */
      if (sum_t_intersect < 1./max_propagation_loops
          && exit_face >= 0
          && cell_id != save_old_cell_id) {
        specific_face_interaction = true;
      }

      if (   !specific_face_interaction
          && (cell_id != save_old_cell_id || exit_face < 0) ) {
        /* if no rebound nor leaving, nor deposition nor oscillation around a
         * given face.
         * Integrate of the stochastic particle over estimated residence time */

        _integ_particle_quantities(particles,
                                   p_id,
                                   nor, // 1
                                   dt_part,
                                   taup,
                                   tlag,
                                   piil,
                                   force_p,
                                   bx,
                                   tempct,
                                   &tsfext,
                                   visc_length,
                                   beta,
                                   vagaus,
                                   br_gaus,
                                   &cpgd1,
                                   &cpgd2,
                                   &cpght,
                                   nresnew,
                                   save_specific_face_interaction);

        /* set back the new cell_id to the particle */

        if (cs_glob_lagr_time_scheme->t_order == 2) {
          nor = 2;

          for (int phase_id =0; phase_id < n_phases; phase_id++)
            cs_lagr_car(0, /* iprev = 0: Use fields at current time step    */
                        phase_id,
                        p_id,
                        nor,
                        dt_part,
                        &taup[phase_id],
                        tlag[phase_id],
                        piil[phase_id],
                        bx[phase_id],
                        tempct,
                        beta,
                        vagaus,
                        br_gaus);

          cs_lagr_get_force_p(dt_part, p_id, taup, tlag, piil, bx,
                              tsfext, vagaus, force_p);

          _integ_particle_quantities(particles,
                                     p_id,
                                     nor, // 2
                                     dt_part, // rel_integ_time
                                     taup,
                                     tlag,
                                     piil,
                                     force_p,
                                     bx,
                                     tempct,
                                     &tsfext,
                                     visc_length,
                                     beta,
                                     vagaus,
                                     br_gaus,
                                     &cpgd1,
                                     &cpgd2,
                                     &cpght,
                                     nresnew,
                                     save_specific_face_interaction);
        }
      }

      /* If we track the deterministic virtual partner and there is a rebound
       * we track the stochastic particle so it sees the rebound condition */
      if (specific_face_interaction && particle_state == CS_LAGR_PART_TO_SYNC)
        save_specific_face_interaction = true;

      /* Increment quantities associated to the previous cell if required */

      cs_real_t dt_incr = dt_part - dt_incremented_in_subiter;

      assert(dt_incr > 0.);
      /* Increment two-way-coupling source terms */
      if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING)
        cs_lagr_coupling_increment_part_contrib(particles,
                                                p_id,
                                                dt_incr,
                                                specific_face_interaction,
                                                taup[0],
                                                force_p,
                                                tempct);

      /* Increment stat in the previous cell */
      if(    cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt)
        cs_lagr_stat_update_all_incr(particles, p_id,
                                     dt_incr /cs_glob_lagr_time_step->dtp);

      if (  (cell_id == save_old_cell_id && exit_face >= 0)
          || specific_face_interaction) {
        /* save the incremented time but continue tracking virtual partner
         * until next face to avoid spurious oscillation around a face */
        dt_incremented_in_subiter = dt_part;

        save_old_cell_id = cs_lagr_particle_get_lnum_n(particle, p_am, 1,
                                                       CS_LAGR_CELL_ID);
        cs_lagr_particle_set_lnum_n(particle, p_am, 1,CS_LAGR_CELL_ID, cell_id);
      }

      else {
        save_old_cell_id = cs_lagr_particle_get_lnum_n(particle, p_am, 1,
                                                       CS_LAGR_CELL_ID);

        /* Update previous state after integration */
        cs_lagr_particles_current_to_previous(particles, p_id);

        /* update remaining time */
        if (exit_face < 0) {
          /* deterministic local partner has reached its final location no
           * new integration on the position occur */
          cs_lagr_particle_set_real(particle, p_am,
                                    CS_LAGR_REMAINING_INTEG_TIME, -1.);
          /* search directly the final stochastic particle location
           * without intermediate step */
          save_specific_face_interaction = false;
        }
        else {
          /*  continue tracking the deterministic virtual partner*/
          remain_time = remain_time * (1. - sum_t_intersect);

          cs_lagr_particles_set_real(particles,
                                     p_id,
                                     CS_LAGR_REMAINING_INTEG_TIME,
                                     remain_time);
        }
        sum_t_intersect = 0.;

        /* temporary track stochastic particle to apply specific interaction
         * on the stochastic particle */
        if (save_specific_face_interaction)
          p_info->tracking_step_id = 1;

      }
    }
  } /* End of while : local displacement */

  CS_FREE(taup);
  CS_FREE(tlag);
  CS_FREE(piil);
  CS_FREE(bx);
  CS_FREE(vagaus);
  CS_FREE(null_vagaus);

  assert(particle_state != CS_LAGR_PART_TO_SYNC);
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
                  CS_MPI_LNUM,
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
                  CS_MPI_LNUM,
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
 *  halo           <-- pointer to a cs_halo_t structure
 *  lag_halo       <-> pointer to a cs_lagr_halo_t structure
 *  particles          <-- set of particles to update
 *  particle_range <-> start and past-the-end ids of tracked particles
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER < 1800
#pragma optimization_level 1 /* Bug with O2 or above with icc 17.0.0 20160721 */
#endif
#endif

static void
_exchange_particles(const cs_halo_t         *halo,
                    cs_lagr_halo_t          *lag_halo,
                    cs_lagr_particle_set_t  *particles,
                    cs_lnum_t                particle_range[2])
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

      for (cs_lnum_t i = 0; i < lag_halo->send_count[local_rank_id]; i++) {
        memcpy(particles->p_buffer + tot_extents*(recv_shift + i),
               lag_halo->send_buf + tot_extents*(send_shift + i),
               tot_extents);

      }
    }
  }

  /* Update particle count and weight */

  cs_real_t tot_weight = 0.;

  for (cs_lnum_t i = 0; i < n_recv_particles; i++) {

    cs_lnum_t j = particles->n_particles + i;

    cs_real_t cur_part_stat_weight
      = cs_lagr_particles_get_real(particles, j, CS_LAGR_STAT_WEIGHT);

    tot_weight += cur_part_stat_weight;
    _tracking_info(particles, j)->state = CS_LAGR_PART_TO_SYNC;
  }

  particles->n_particles += n_recv_particles;
  particles->weight += tot_weight;

  particle_range[1] += n_recv_particles;
}

/*----------------------------------------------------------------------------
 * Determine particle halo sizes
 *
 * parameters:
 *   mesh           <-- pointer to associated mesh
 *   lag_halo       <-> pointer to particle halo structure to update
 *   particles      <-- set of particles to update
 *   particle_range <-- start and past-the-end ids of tracked particles
 *----------------------------------------------------------------------------*/

static void
_lagr_halo_count(const cs_mesh_t               *mesh,
                 cs_lagr_halo_t                *lag_halo,
                 const cs_lagr_particle_set_t  *particles,
                 const cs_lnum_t                particle_range[2])
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

  for (i = particle_range[0]; i < particle_range[1]; i++) {

    if (_get_tracking_info(particles, i)->state == CS_LAGR_PART_TO_SYNC_NEXT) {

      ghost_id =   cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_ID)
                 - mesh->n_cells;

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
 *   particles      <-> set of particles to update
 *   particle_range <-> start and past-the-end ids of tracked particles
 *
 * returns:
 *   1 if displacement needs to continue, 0 if finished
 *----------------------------------------------------------------------------*/

static int
_sync_particle_set(cs_lagr_particle_set_t  *particles,
                   cs_lnum_t                particle_range[2])
{
  cs_lnum_t  i, k, tr_id, rank, shift, ghost_id;
  cs_real_t matrix[3][4];

  cs_lnum_t  particle_count = 0;

  cs_lnum_t  n_merged_particles = 0;

  cs_lnum_t  n_exit_particles = 0;
  cs_lnum_t  n_failed_particles = 0;

  cs_real_t  exit_weight = 0.0;
  cs_real_t  merged_weight = 0.0;
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

  if (halo != nullptr) {

    _lagr_halo_count(mesh, lag_halo, particles, particle_range);

    for (i = 0; i < halo->n_c_domains; i++) {
      lag_halo->send_count[i] = 0;
    }
  }

  /* Loop on particles, transferring particles to synchronize to send_buf
     for particle set, and removing particles that otherwise exited the domain */

  for (i = particle_range[0]; i < particle_range[1]; i++) {

    cs_lagr_tracking_state_t cur_part_state
      = _get_tracking_info(particles, i)->state;

    cs_real_t cur_part_stat_weight
      = cs_lagr_particles_get_real(particles, i, CS_LAGR_STAT_WEIGHT);

    /* Particle changes domain */

    if (cur_part_state == CS_LAGR_PART_TO_SYNC_NEXT) {

      continue_displacement = 1;

      ghost_id =   cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_ID)
                 - halo->n_local_elts;
      rank = lag_halo->rank[ghost_id];
      tr_id = lag_halo->transform_id[ghost_id];

      cs_lagr_particles_set_lnum(particles, i, CS_LAGR_CELL_ID,
                                 lag_halo->dist_cell_id[ghost_id]);

      shift = lag_halo->send_shift[rank] + lag_halo->send_count[rank];

      /* Update if needed last_face_num */

      if (tr_id >= 0) { /* Same initialization as in previous algorithm */

        _tracking_info(particles, i)->last_face_id = -1;

      }

      else {

        if (cs_glob_n_ranks > 1) {

          assert(face_ifs != nullptr);

          int  distant_rank;

          const int search_rank = halo->c_domain_rank[rank];
          const cs_interface_t  *interface = nullptr;
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

            cs_lnum_t n_entities = cs_interface_size(interface);
            const cs_lnum_t *local_num = cs_interface_get_elt_ids(interface);

            cs_lnum_t id = cs_search_binary
                             (n_entities,
                              _get_tracking_info(particles, i)->last_face_id,
                              local_num);

            if (id == -1)
              bft_error(__FILE__, __LINE__, 0,
                        _(" Cannot find the relative distant face id.\n"));

            const cs_lnum_t *dist_ids = cs_interface_get_match_ids(interface);

            _tracking_info(particles, i)->last_face_id = dist_ids[id];

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
          cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, i, CS_LAGR_COORDS));

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
                              _tracking_info(particles, i)->start_coords);

        _apply_vector_transfo((const cs_real_t (*)[4])matrix,
          cs_lagr_particles_attr_n_get_ptr<cs_real_t>(particles, i, 1,
                                                      CS_LAGR_COORDS));

        /* Apply rotation to velocity vectors in case of rotation */

        if (perio_type >= FVM_PERIODICITY_ROTATION) {

          /* Rotation of the velocity */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
            cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, i,
                                                      CS_LAGR_VELOCITY));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
              cs_lagr_particles_attr_n_get_ptr<cs_real_t>(particles, i, 1,
                                                          CS_LAGR_VELOCITY));

          /* Rotation of the velocity seen */

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
            cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, i,
                                                        CS_LAGR_VELOCITY_SEEN));

          _apply_vector_rotation((const cs_real_t (*)[4])matrix,
            cs_lagr_particles_attr_n_get_ptr<cs_real_t>(particles, i, 1,
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
      merged_weight += cur_part_stat_weight;
    }

    else if (cur_part_state < CS_LAGR_PART_OUT) {

      if (particle_count + particle_range[0] < i) {
        memcpy(particles->p_buffer + p_am->extents*
                                        (particle_count +  particle_range[0]),
               particles->p_buffer + p_am->extents*i,
               p_am->extents);
      }

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

  particles->n_particles += particle_count - particle_range[1]
                                           + particle_range[0];
  particle_range[1] = particle_range[0] + particle_count;
  /* Only the weight of within particle_range if different from the whole set*/
  particles->weight = tot_weight;

  particles->n_part_out += n_exit_particles;
  particles->weight_out += exit_weight;

  particles->n_failed_part += n_failed_particles;
  particles->weight_failed += fail_weight;

  particles->n_part_merged += n_merged_particles;
  particles->weight_merged += merged_weight;

  /* Exchange particles, then update set */

  if (halo != nullptr)
    _exchange_particles(halo, lag_halo, particles, particle_range);

  cs_parall_max(1, CS_INT_TYPE, &continue_displacement);

  return continue_displacement;
}

/*----------------------------------------------------------------------------
 * Prepare for particle movement phase
 *
 * parameters:
 *   particles      <->  pointer to particle set structure
 *   particle_range <--  start and past-the-end ids of tracked particles
 *   resol_sde      <--  true  the sdes will be resolved
 *                       false only the trajectography is done
 *----------------------------------------------------------------------------*/

static void
_initialize_displacement(cs_lagr_particle_set_t  *particles,
                         const cs_lnum_t          particle_range[2],
                         const bool               resol_sde)
{
  /* Initialize builder if needed */

  if (_particle_track_builder == nullptr)
    _particle_track_builder
      = _init_track_builder(particles->n_particles_max,
                            particles->p_am->extents);

  assert(particles->p_am->lb >= sizeof(cs_lagr_tracking_info_t));

  /* Info for rotor-stator cases; the time step should actually
     be based on the global (non-Lagrangian) time step in case
     the two differ. Currently both are in lock-step; if one
     becomes a mutiplier of the other, the rotor-movement update
     would need to be done at a frequency matching that multiplier. */

  const int *cell_rotor_num = nullptr;
  cs_real_34_t  *rot_m = nullptr;

  if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
    cell_rotor_num = cs_turbomachinery_get_cell_rotor_num();
    cs_real_t dt = cs_glob_lagr_time_step->dtp;
    rot_m = cs_turbomachinery_get_rotation_matrices(-1.0 * dt);
  }

  /* Prepare tracking info */
  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    cs_lnum_t cur_part_cell_id
      = cs_lagr_particles_get_lnum(particles, p_id, CS_LAGR_CELL_ID);
    cs_real_t r_truncate
      = cs_lagr_particles_get_real(particles, p_id, CS_LAGR_TR_TRUNCATE);

    int p_flag = cs_lagr_particles_get_lnum(particles, p_id, CS_LAGR_P_FLAG);

    if (p_flag & CS_LAGR_PART_TO_DELETE) {
      _tracking_info(particles, p_id)->state = CS_LAGR_PART_OUT;
    }

    else if (p_flag & CS_LAGR_PART_FIXED) {
      cs_real_t particle_weight
        = cs_lagr_particles_get_real(particles, p_id, CS_LAGR_STAT_WEIGHT);
      if (particle_weight <= 0.) /* == 0, use <= to avoid warning */
        _tracking_info(particles, p_id)->state = CS_LAGR_PART_MERGED;
      else
        _tracking_info(particles, p_id)->state = CS_LAGR_PART_STUCK;
    }

    else if (r_truncate > 1.9){ /* from previous displacement */
      _tracking_info(particles, p_id)->state = CS_LAGR_PART_ERR;
    }
    else {
      _tracking_info(particles, p_id)->state = CS_LAGR_PART_TO_SYNC;
      int depo_flag = p_flag & CS_LAGR_PART_DEPOSITION_FLAGS;
      if (depo_flag == CS_LAGR_PART_DEPOSITED) { /* deposited, no other flag */
        _tracking_info(particles, p_id)->state = CS_LAGR_PART_TREATED;
      }
    }
    if (   cs_glob_lagr_time_scheme->cell_wise_integ == 1
        && _tracking_info(particles, p_id)->state != CS_LAGR_PART_TO_SYNC)
      cs_lagr_particles_set_real(particles, p_id, CS_LAGR_REMAINING_INTEG_TIME,
                                                  -1);

    _tracking_info(particles, p_id)->last_face_id = -1;

    /* Coordinates of the particle */

    cs_real_t *prv_part_coord
      = cs_lagr_particles_attr_n_get_ptr<cs_real_t>(particles, p_id, 1,
                                                    CS_LAGR_COORDS);
    _tracking_info(particles, p_id)->start_coords[0] = prv_part_coord[0];
    _tracking_info(particles, p_id)->start_coords[1] = prv_part_coord[1];
    _tracking_info(particles, p_id)->start_coords[2] = prv_part_coord[2];


    if (cell_rotor_num != nullptr) {
      int r_num = cell_rotor_num[cur_part_cell_id];
      if (r_num > 0) {
        _apply_vector_transfo((const cs_real_t (*)[4])rot_m[r_num],
                              _tracking_info(particles, p_id)->start_coords);
      }
    }

    cs_lagr_particles_set_real(particles, p_id, CS_LAGR_TR_TRUNCATE, 0);
    cs_lagr_particles_set_lnum(particles, p_id, CS_LAGR_TR_REPOSITION, 0);

  }


  if (cs_glob_lagr_time_scheme->cell_wise_integ == 1 && resol_sde)
  {
    /* track ghost trajectory associated to deterministic virtual partner*/
    for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++)
      _tracking_info(particles, p_id)->tracking_step_id = 0;
  }
  else {
    /* track trajectory associated directly to the stochastic particle*/
    for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++)
      _tracking_info(particles, p_id)->tracking_step_id = 1;
  }
  CS_FREE(rot_m);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Particle set after %s\n", __func__);
  cs_lagr_particle_set_dump(particles);
#endif
}

/*----------------------------------------------------------------------------
 * Update particle set structures: compact array.
 *
 * parameters:
 *   particles      <-> pointer to particle set structure
 *----------------------------------------------------------------------------*/

static void
_finalize_displacement(cs_lagr_particle_set_t  *particles)
{
  const cs_lagr_attribute_map_t  *p_am = particles->p_am;
  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;

  const cs_lnum_t n_particles = particles->n_particles;

  cs_lnum_t *cell_idx;
  unsigned char *swap_buffer;
  size_t swap_buffer_size = p_am->extents * ((size_t)n_particles);

  CS_MALLOC(cell_idx, n_cells+1, cs_lnum_t);
  CS_MALLOC(swap_buffer, swap_buffer_size, unsigned char);

  /* Cell index (count first) */

  for (cs_lnum_t i = 0; i < n_cells+1; i++)
    cell_idx[i] = 0;

  /* Copy unordered particle data to buffer */

  for (cs_lnum_t i = 0; i < n_particles; i++) {

    cs_lnum_t cur_part_state = _get_tracking_info(particles, i)->state;

    CS_NO_WARN_IF_UNUSED(cur_part_state);
    assert(   cur_part_state < CS_LAGR_PART_OUT
           && cur_part_state != CS_LAGR_PART_TO_SYNC);

    cs_lnum_t cell_id = cs_lagr_particles_get_lnum(particles, i,
                                                   CS_LAGR_CELL_ID);

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
  const cs_lnum_t cell_num_displ = particles->p_am->displ[0][CS_LAGR_CELL_ID];

  for (cs_lnum_t i = 0; i < n_particles; i++) {

    cs_lnum_t cell_id
      = *((const cs_lnum_t *)(swap_buffer + p_extents*i + cell_num_displ));

    assert(cell_id > -1);

    cs_lnum_t particle_id = cell_idx[cell_id];

    cell_idx[cell_id] += 1;

    memcpy(particles->p_buffer + p_am->extents*particle_id,
           swap_buffer + p_am->extents*i,
           p_am->extents);

  }

  CS_FREE(swap_buffer);
  CS_FREE(cell_idx);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Particle set after %s\n", __func__);
  cs_lagr_particle_set_dump(particles);
#endif

  if (cs_glob_mesh->time_dep == CS_MESH_TRANSIENT_CONNECT)
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
 * \brief Integrate or not SDEs associated to the particle and apply one
 * trajectography step to track the displacement.
 *
 * \param[in]     visc_length     viscous layer thickness
 * \param[in,out] particle_range  start and past-the-end ids of tracked
 * \param[in,out] particle_range  start and past-the-end ids of tracked
 * \param[in]     resol_sde       true  the sdes will be resolved
 *                                false only the trajectography step is done
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_integ_track_particles(const cs_real_t  visc_length[],
                              cs_lnum_t        particle_range[2],
                              const bool       resol_sde)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;

  int  displacement_step_id = 0;
  int  continue_displacement = 1;

  cs_lagr_particle_set_t  *particles = cs_glob_lagr_particle_set;
  cs_lagr_event_set_t     *events = nullptr;

  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  int n_phases = extra->n_phases;

  const cs_lagr_attribute_map_t  *p_am = particles->p_am;

  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;

  const cs_lnum_t  failsafe_mode = 0; /* If 1 : stop as soon as an error is
                                         detected */

  int t_stat_id = cs_timer_stats_id_by_name("particle_displacement_stage");
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  int t_top_id = cs_timer_stats_switch(t_stat_id);

  if (cs_lagr_stat_is_active(CS_LAGR_STAT_GROUP_TRACKING_EVENT)) {
    events = cs_lagr_event_set_boundary_interaction();
    /* Event set "expected" size: n boundary faces*2 */
    cs_lnum_t events_min_size = mesh->n_b_faces * 2;
    if (events->n_events_max < events_min_size)
      cs_lagr_event_set_resize(events, events_min_size);
  }

  int cell_wise_integ = cs_glob_lagr_time_scheme->cell_wise_integ;

  assert(particles != nullptr);

  const int *b_face_zone_id = cs_boundary_zone_face_class_id();

  _initialize_displacement(particles, particle_range, resol_sde);

  int nresnew = 0;

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  /* particle properties */
  cs_real_2_t tempct;
  cs_real_t cpgd1, cpgd2, cpght, tsfext;
  cs_real_3_t beta = { 0., 0., 0.};
  cs_real_3_t force_p;
  cs_real_6_t br_gaus;
  cs_real_t *taup = nullptr;
  cs_real_3_t *tlag = nullptr;
  cs_real_3_t *piil = nullptr;
  cs_real_3_t *vagaus = nullptr;
  cs_real_33_t *bx = nullptr;
  if (resol_sde && cell_wise_integ == 0) {
    CS_MALLOC(taup, n_phases, cs_real_t);
    CS_MALLOC(tlag, n_phases, cs_real_3_t);
    CS_MALLOC(piil, n_phases, cs_real_3_t);
    CS_MALLOC(vagaus, n_phases + 2, cs_real_3_t);
    CS_MALLOC(bx, n_phases, cs_real_33_t);
  }

  int nor = 1;

  cs_real_t   *list_taup = nullptr;
  cs_real_3_t *list_force_p = nullptr;
  cs_real_2_t *list_tempct = nullptr;

  if (resol_sde && cell_wise_integ == 0) {
    cs_real_t dtp = cs_glob_lagr_time_step->dtp;
    /* if two-way coupling used: save taup force_p and temcpt for coupling
     * Reinit of t_st_... after cs_lagr_coupling as lagr_st_vel used in piil*/
    if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
        && cs_glob_lagr_time_scheme->t_order == 1) {
      CS_MALLOC(list_taup, particles->n_particles, cs_real_t);
      CS_MALLOC(list_force_p, particles->n_particles, cs_real_3_t);
      CS_MALLOC(list_tempct, particles->n_particles, cs_real_2_t);
    }
    for (cs_lnum_t p_id = 0; p_id < particles->n_particles; p_id++) {
      unsigned char *particle = particles->p_buffer + p_am->extents * p_id;
      cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;

      if (cs_lagr_particles_get_flag(particles, p_id, CS_LAGR_PART_FIXED)) {
        _tracking_info(particles, p_id)->state = CS_LAGR_PART_STUCK;
        continue;
      }
      for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
        cs_lagr_car(1, /* iprev = 1: Use fields at previous time step    */
                    phase_id,
                    p_id,
                    1, // nor
                    dtp,
                    &(taup[phase_id]),
                    tlag[phase_id],
                    piil[phase_id],
                    bx[phase_id],
                    tempct,
                    beta,
                    vagaus,
                    br_gaus);

        /* Save bx values associated with particles for next pass */
        if (cs_glob_lagr_time_scheme->t_order == 2) {
          cs_real_t *jbx1 =
            cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, p_id,
                                                      CS_LAGR_TURB_STATE_1);

          for (cs_lnum_t i = 0; i < 3; i++)
            jbx1[9 * phase_id + i] = bx[phase_id][i][0];
        }
      }

      cs_lagr_get_force_p(dtp, p_id, taup, tlag, piil, bx,
                          tsfext, vagaus, force_p);

      /* Integrate the particle associated SDEs */
      _integ_particle_quantities(particles,
                                 p_id,
                                 nor, // 1
                                 dtp,
                                 taup,
                                 tlag,
                                 piil,
                                 force_p,
                                 bx,
                                 tempct,
                                 &tsfext,
                                 visc_length,
                                 beta,
                                 vagaus,
                                 br_gaus,
                                 &cpgd1,
                                 &cpgd2,
                                 &cpght,
                                 &nresnew,
                                 true);

      /* save field particle quantities if still needed for coupling */
      if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
          && cs_glob_lagr_time_scheme->t_order == 1) {
        list_taup[p_id] = taup[0];
        for (int i = 0; i < 3; i++)
          list_force_p[p_id][i] = force_p[i];
        list_tempct[p_id][0] = tempct[0];
        list_tempct[p_id][1] = tempct[1];
      }

      cs_real_t  *particle_coord
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_COORDS);
      /* starting point for the tracking*/
      cs_real_t  *prev_location = p_info->start_coords;

      cs_real_t  disp[3];
      /* Assert if there is a displacement */
      for (int k = 0; k < 3; k++)
        disp[k] = particle_coord[k] - prev_location[k];

      cs_lnum_t cell_id = cs_lagr_particles_get_lnum(particles, p_id,
                                                    CS_LAGR_CELL_ID);
      cs_real_t inv_ref_length = 1. / pow(cell_vol[cell_id], 1./3.);

      if (   fabs(disp[0] * inv_ref_length) < 1e-15
          && fabs(disp[1] * inv_ref_length) < 1e-15
          && fabs(disp[2] * inv_ref_length) < 1e-15) {
        p_info->state = CS_LAGR_PART_TREATED;
      }
    }
    if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
        && cs_glob_lagr_time_scheme->t_order == 1) {
      cs_lagr_coupling_initialize();
      for (cs_lnum_t p_id = 0; p_id < particles->n_particles; p_id++) {
        cs_lagr_coupling_increment_part_contrib(particles,
                                                p_id,
                                                dtp,
                                                false,
                                                list_taup[p_id],
                                                list_force_p[p_id],
                                                list_tempct[p_id]);

      }
      CS_FREE(list_taup);
      CS_FREE(list_force_p);
      CS_FREE(list_tempct);
    }
  } /* end resol_sde */

  /* limit number of ranks and periodicity crossed to avoid infinite loop */
  int  max_perio_or_rank_crossed =
    cs_glob_lagr_time_scheme->max_perio_or_rank_crossed;
  if (!resol_sde)
    max_perio_or_rank_crossed = cs::max(cs_glob_n_ranks,
                                        max_perio_or_rank_crossed);
  /* Main loop on particles: global propagation */
  while (continue_displacement) {

    /* Local propagation */
    for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++)
    {
      unsigned char *particle = particles->p_buffer + p_am->extents * p_id;
      cs_lagr_tracking_info_t *p_info = (cs_lagr_tracking_info_t *)particle;
      cs_lagr_tracking_state_t cur_part_state = p_info->state;

      /* assert if there is still particles to track and if so propage them*/
      if (cur_part_state == CS_LAGR_PART_TO_SYNC) {

        /* Main particle displacement stage */
        cur_part_state =
          (cs_lagr_tracking_state_t)_local_propagation(particles,
                                                       events,
                                                       p_id,
                                                       failsafe_mode,
                                                       b_face_zone_id,
                                                       visc_length,
                                                       resol_sde,
                                                       &nresnew);

        if(displacement_step_id == max_perio_or_rank_crossed -1) {
          cs_real_t  *particle_coord
            = cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, p_id,
                                                        CS_LAGR_COORDS);
          cs_lnum_t cell_id = cs_lagr_particles_get_lnum(particles, p_id,
                                                         CS_LAGR_CELL_ID);
          _manage_error(failsafe_mode,
                       particle,
                       p_am,
                       CS_LAGR_TRACKING_ERR_MAX_LOOPS);

          p_info->state = CS_LAGR_PART_TREATED;
          cs_real_t *cell_cen = fvq->cell_cen[cell_id];

          cs_lnum_t n_rep =
            cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_TR_REPOSITION);
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_TR_REPOSITION,
                                    n_rep+1);

          for (int k = 0; k < 3; k++) {
            p_info->start_coords[k] = cell_cen[k];
            particle_coord[k] = cell_cen[k];
          }
        }
        p_info->state = cur_part_state;

      }
    } /* End of loop on particles */
    /* Update of the particle set structure. Delete exited particles,
     update for particles which change domain. */

    continue_displacement = _sync_particle_set(particles, particle_range);

#if 0
    bft_printf("\n Particle set after sync\n");
    cs_lagr_particle_set_dump(particles);
#endif

    /*  assert(j == -1);  After a loop on particles, next_id of the last
        particle must not be defined */

    displacement_step_id++;

    /* Update size of the particle set if it has changed */
    particles->n_particles += nresnew;
  } /* End of while (global displacement) */

  if (   cell_wise_integ == 0
      && cs_glob_lagr_time_scheme->t_order == 2
      && resol_sde) {

    if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {
      CS_MALLOC(list_taup, particles->n_particles, cs_real_t);
      CS_MALLOC(list_force_p, particles->n_particles, cs_real_3_t);
      CS_MALLOC(list_tempct, particles->n_particles, cs_real_2_t);
    }

    cs_real_t dtp = cs_glob_lagr_time_step->dtp;
    nor = 2;

    for (cs_lnum_t p_id = 0; p_id < particles->n_particles; p_id++) {
      bool has_rebound_occured =
        (cs_lagr_particles_get_lnum(particles, p_id, CS_LAGR_REBOUND_ID) == 0 );

      for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
      /*   Retrieve bx values associated with particles from previous pass */
        cs_real_t *jbx1 =
          cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, p_id,
                                                    CS_LAGR_TURB_STATE_1);

        for (cs_lnum_t i = 0; i < 3; i++)
          bx[phase_id][i][0] = jbx1[3 * phase_id + i];

        cs_lagr_car(0, /* iprev = 0: Use fields at current time step    */
                    phase_id,
                    p_id,
                    nor,
                    dtp,
                    &taup[phase_id],
                    tlag[phase_id],
                    piil[phase_id],
                    bx[phase_id],
                    tempct,
                    beta,
                    vagaus,
                    br_gaus);
      }

      cs_lagr_get_force_p(dtp, p_id, taup, tlag, piil, bx,
                          tsfext, vagaus, force_p);

      if (!has_rebound_occured)
        _integ_particle_quantities(particles,
                                   p_id,
                                   nor, // 2
                                   dtp,
                                   taup,
                                   tlag,
                                   piil,
                                   force_p,
                                   bx,
                                   tempct,
                                   &tsfext,
                                   visc_length,
                                   beta,
                                   vagaus,
                                   br_gaus,
                                   &cpgd1,
                                   &cpgd2,
                                   &cpght,
                                   &nresnew,
                                   true);

      /* save field particle quantities if still needed for coupling */
      if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {
        list_taup[p_id] = taup[0];
        for (int i = 0; i < 3; i++)
          list_force_p[p_id][i] = force_p[i];
        list_tempct[p_id][0] = tempct[0];
        list_tempct[p_id][1] = tempct[1];
      }

      /* Increment stat in the previous cell */
      if(    cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt)
        cs_lagr_stat_update_all_incr(particles, p_id, 1);
    }

    if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

      cs_lagr_coupling_initialize();

      for (cs_lnum_t p_id = 0; p_id < particles->n_particles; p_id++) {
        bool has_rebound_occured =
          (cs_lagr_particles_get_lnum(particles, p_id, CS_LAGR_REBOUND_ID) == 0);

        cs_lagr_coupling_increment_part_contrib(particles,
                                                p_id,
                                                dtp,
                                                has_rebound_occured,
                                                list_taup[p_id],
                                                list_force_p[p_id],
                                                list_tempct[p_id]);
      }
    }
  }
  CS_FREE(list_taup);
  CS_FREE(list_force_p);
  CS_FREE(list_tempct);
  CS_FREE(taup);
  CS_FREE(tlag);
  CS_FREE(piil);
  CS_FREE(bx);
  CS_FREE(vagaus);

  /* Deposition sub-model additional loop */

  if (lagr_model->deposition > 0) {

    for (cs_lnum_t i = particle_range[0]; i < particle_range[1]; i++) {

      unsigned char *particle = particles->p_buffer + p_am->extents * i;

      cs_lnum_t *neighbor_face_id
        = cs_lagr_particle_attr_get_ptr<cs_lnum_t>(particle, p_am,
                                                   CS_LAGR_NEIGHBOR_FACE_ID);
      cs_real_t *particle_yplus
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_YPLUS);

      cs_lagr_test_wall_cell(particle, p_am, visc_length,
                             particle_yplus, neighbor_face_id);

      /* Modification of MARKO pointer */
      if (*particle_yplus > 100.0)
        cs_lagr_particles_set_lnum(particles, i,
                                   CS_LAGR_MARKO_VALUE,
                                   CS_LAGR_COHERENCE_STRUCT_BULK);

      else {

        if (  *particle_yplus
            < cs_lagr_particles_get_real(particles, i, CS_LAGR_INTERF)) {

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

  /* If the whole set of particles is tracked rearrange particles */
  if (particle_range[1] - particle_range[0] == particles->n_particles)
    _finalize_displacement(particles);

  if (   cs_glob_porous_model == 3
      && lagr_model->deposition == 1)
  cs_lagr_porosity();

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
  if (cs_glob_lagr_particle_set == nullptr)
    return;

  /* Destroy event structures */

  cs_lagr_event_finalize();

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
  cs_lnum_t cell_id
    = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);

  *yplus = 10000;
  *face_id = -1;

  cs_lnum_t  *cell_b_face_idx = cs_glob_mesh_adjacencies->cell_b_faces_idx;
  cs_lnum_t  *cell_b_faces = cs_glob_mesh_adjacencies->cell_b_faces;
  const cs_nreal_3_t *restrict b_face_u_normal
    = cs_glob_mesh_quantities->b_face_u_normal;
  const cs_real_3_t *restrict b_face_cog = cs_glob_mesh_quantities->b_face_cog;

  const cs_real_t  *particle_coord
    = cs_lagr_particle_attr_get_const_ptr<cs_real_t>(particle, p_am,
                                                     CS_LAGR_COORDS);

  cs_lnum_t  start = cell_b_face_idx[cell_id];
  cs_lnum_t  end =  cell_b_face_idx[cell_id + 1];

  for (cs_lnum_t i = start; i < end; i++) {
    cs_lnum_t f_id = cell_b_faces[i];

    assert(cs_glob_lagr_boundary_conditions != nullptr);

    const char b_type = cs_glob_lagr_boundary_conditions->elt_type[f_id];

    if (   (b_type == CS_LAGR_DEPO1)
        || (b_type == CS_LAGR_DEPO2)
        || (b_type == CS_LAGR_DEPO_DLVO)) {

      /* normal vector coordinates */
      const cs_nreal_t *normal = b_face_u_normal[f_id];

      /* [(x_f - x_p) . n ] / L */
      cs_real_t dist_norm = cs::abs(
          cs_math_3_distance_dot_product(b_face_cog[f_id],
                                         particle_coord,
                                         normal)) / visc_length[f_id];
      if (dist_norm  < *yplus) {
        *yplus = dist_norm;
        *face_id = f_id;
      }
    }

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
