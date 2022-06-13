/*============================================================================
 * Functions dealing with ghost cells
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_order.h"

#include "cs_interface.h"
#include "cs_rank_neighbors.h"

#include "fvm_periodicity.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"

#if defined(HAVE_CUDA)
#include "cs_halo_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/* Remarks:
 *
 * The current available mode for MPI-3 RMA uses "get" semantics.
 * A "put" semantics variant could easily be added, either:
 * - Using MPI_Win_create_dynamic and attaching the halo section of the
 *   current array to that window (which would seem cumbersome as this also
 *   requires exchanging base addresses obtained with MPI_Get_Address for
 *   each array, but could be amortized for iterative algorithms).
 * - Using a fixed receive buffer, and copying data to the tail section of
 *   the array afterwards (as is done for accelerators when the MPI library
 *   used is not no accelerator-aware). This would add an extra copy, but
 *   be much simpler.
 *
 * It may also be useful to allow this setting on a "per halo" basis, as
 * some publications report better performance with RMA for large data,
 * and better performance with P2P for small data, so in uses such
 * as multigrid solvers, either may be preferred for different levels.
*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#if defined(OMPI_MAJOR_VERSION)
  #include <mpi-ext.h>
#endif

#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  #define _CS_MPI_DEVICE_SUPPORT 1
#else
  #if defined(_CS_MPI_DEVICE_SUPPORT)
    #undef _CS_MPI_DEVICE_SUPPORT
  #endif
#endif

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure to maintain halo exchange state */

struct _cs_halo_state_t {

  /* Current synchronization state */

  cs_halo_type_t  sync_mode;      /* Standard or extended */
  cs_datatype_t   data_type;      /* Datatype */
  int             stride;         /* Number of values per location */

  cs_alloc_mode_t var_location;   /* Allocation info for exchanged variable */

  void        *send_buffer_cur;   /* Send buffer used for current progress
                                     (either _send_buffer or passed by caller) */

  int       n_requests;        /* Number of MPI requests */
  int       local_rank_id;     /* Id of halo for own rank, -1 if not present */

  /* Buffers for synchronization;
     receive buffers only needed for some communication modes */

  size_t       send_buffer_size;  /* Size of send buffer, in bytes */
  size_t       recv_buffer_size;  /* Size of receive buffer, in bytes */

  void        *send_buffer;       /* Send buffer (maintained by this object) */
  void        *recv_buffer;       /* Recv. buffer (maintained by this object) */

#if defined(HAVE_MPI)

  int          request_size;      /* Size of requests and status arrays */

  MPI_Request  *request;          /* Array of MPI requests */
  MPI_Status   *status;           /* Array of MPI status */

  MPI_Win       win;              /* MPI-3 RMA window */

#endif

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of defined halos */
static int _n_halos = 0;

/* Allocation mode for arrays which might be used on accelerator device
   Note that an alternative option would be to use shared memory with
   prefetching. We will need to do performance comparisons first, but
   in the case of similar performance, going for the shared approach
   would be preferred for its "more generic" aspect. */
static cs_alloc_mode_t _halo_buffer_alloc_mode = CS_ALLOC_HOST_DEVICE_PINNED;

/* Should we use barriers after posting receives ? */
static int _halo_use_barrier = false;

/* Default halo state handler */
static cs_halo_state_t *_halo_state = NULL;

/* Halo communications mode */
static int _halo_comm_mode = CS_HALO_COMM_P2P;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if an array of global numbers is ordered.
 *
 * \param[in]  list    optional list (1 to n numbering) of selected entities
 *                     (or NULL if all nb_ent are selected). This list may
 *                     contain element numbers in any order
 * \param[in]  nb_ent  number of entities considered
 *
 * \return  1 if ordered, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

static int
_order_int_test(const int  list[],
                size_t     nb_ent)
{
  size_t i = 0;

  /* If numbering is explicit */

  if (list != NULL) {
    for (i = 1 ; i < nb_ent ; i++) {
      if (list[i] < list[i-1])
        break;
    }
  }
  else
    i = nb_ent;

  if (i == nb_ent || nb_ent == 0)
    return 1;
  else
    return 0;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Update state request arrays so as to be usable with a given halo.
 *
 * This function should be called at the end of any halo creation,
 * so that buffer sizes are increased if necessary.
 *
 * parameters:
 *   halo       <-- pointer to cs_halo_t structure.
 *   halo_state <-> pointer to halo state structure.
 *---------------------------------------------------------------------------*/

static void
_update_requests(const cs_halo_t  *halo,
                 cs_halo_state_t  *hs)
{
  if (halo == NULL)
    return;

  int n_requests = halo->n_c_domains*2;

  if (n_requests > hs->request_size) {
    hs->request_size = n_requests;
    BFT_REALLOC(hs->request, hs->request_size, MPI_Request);
    BFT_REALLOC(hs->status, hs->request_size,  MPI_Status);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange send shift in send buffer for one-sided get.
 *
 * \param[in]  halo  halo structure to update
 */
/*----------------------------------------------------------------------------*/

static void
_exchange_send_shift(cs_halo_t  *halo)
{
  MPI_Comm comm = cs_glob_mpi_comm;
  MPI_Request *request = NULL;
  MPI_Status *status = NULL;

  BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
  BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);

  BFT_REALLOC(halo->c_domain_s_shift, halo->n_c_domains, cs_lnum_t);

  /* Exchange local range with neighbor ranks */

  const int local_rank = CS_MAX(cs_glob_rank_id, 0);

  for (int i = 0; i < halo->n_c_domains; i++) {
    int rank_id = halo->c_domain_rank[i];
    MPI_Irecv(halo->c_domain_s_shift + i,
              1,
              CS_MPI_LNUM,
              rank_id,
              local_rank,
              comm,
              &(request[i]));
  }

  for (int i = 0; i < halo->n_c_domains; i++) {
    int rank_id = halo->c_domain_rank[i];
    MPI_Isend(halo->send_index + 2*i,
              1,
              CS_MPI_LNUM,
              rank_id,
              rank_id,
              comm,
              &(request[halo->n_c_domains + i]));
  }

  MPI_Waitall(halo->n_c_domains*2, request, status);

  BFT_FREE(request);
  BFT_FREE(status);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Local copy from halo send buffer to destination array.
 *
 * This allows pariodicity data which may be present on the local rank to
 * be exchanged without any MPI call.
 *
 * Data is untyped; only its size is given, so this function may also
 * be used to synchronize interleaved multidimendsional data, using
 * size = element_size*dim (assuming a homogeneous environment, at least
 * as far as data encoding goes).
 *
 * parameters:
 *   halo          <-- pointer to halo structure
 *   local_rank_id <-- id of local rank
 *   sync_mode     <-- synchronization mode (standard or extended)
 *   size          <-- datatype size
 *   var_location  <-- Allocation info for exchanged variable (host/device)
 *   send_buf      <-> pointer to send data
 *   val           <-> pointer to destination data
 *----------------------------------------------------------------------------*/

static void
_sync_local(const cs_halo_t  *halo,
            int               local_rank_id,
            cs_halo_type_t    sync_mode,
            size_t            size,
            cs_alloc_mode_t   var_location,
            const void       *send_buf,
            void             *val)
{
  cs_lnum_t end_shift = (sync_mode == CS_HALO_EXTENDED) ? 2 : 1;

  unsigned char *_val = val;
  unsigned char *recv
    = _val + (halo->n_local_elts + halo->index[2*local_rank_id]) * size;

  cs_lnum_t start = halo->send_index[2*local_rank_id]*size;
  cs_lnum_t length = (  halo->send_index[2*local_rank_id + end_shift]
                      - halo->send_index[2*local_rank_id]);

  size_t count = length * size;

  if (var_location == CS_ALLOC_HOST) {
    const unsigned char *buffer = (const unsigned char *)send_buf;
    const unsigned char *_buffer = buffer + start;
    memcpy(recv, _buffer, count);
  }

#if defined(HAVE_ACCEL)
  else {
    const unsigned char *buffer
      = (const unsigned char *)cs_get_device_ptr_const(send_buf);
    const unsigned char *_buffer = buffer + start;
    cs_copy_d2d(recv, _buffer, count);
  }
#endif
}

#if defined(HAVE_MPI)
#if (MPI_VERSION >= 3)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Launch update of ghost values in case of parallelism
 *        for one-sided communication.
 *
 * The cs_halo_sync_pack function should have been called before this function,
 * using the same hs argument.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       val         pointer to variable value array
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

static void
_halo_sync_start_one_sided(const cs_halo_t  *halo,
                           void             *val,
                           cs_halo_state_t  *hs)
{
  cs_lnum_t end_shift = (hs->sync_mode == CS_HALO_EXTENDED) ? 2 : 1;
  cs_lnum_t stride = hs->stride;
  size_t elt_size = cs_datatype_size[hs->data_type] * stride;
  size_t n_loc_elts = halo->n_local_elts;

  unsigned char *restrict _val = val;
  unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

  MPI_Datatype mpi_datatype = cs_datatype_to_mpi[hs->data_type];

  const int local_rank = CS_MAX(cs_glob_rank_id, 0);

  /* Get data from distant ranks */

  if (_halo_comm_mode == CS_HALO_COMM_RMA_GET) {

    /* Use active target synchronization */

    if (halo->c_domain_group != MPI_GROUP_NULL) {
      /* Start RMA exposure epoch */
      MPI_Win_post(halo->c_domain_group,
                   MPI_MODE_NOPUT,          /* program assertion */
                   hs->win);

      /* Access Epoch */
      MPI_Win_start(halo->c_domain_group,
                    0,                      /* program assertion */
                    hs->win);
    }
    else {
      MPI_Win_fence(0, hs->win);
    }

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      int  length = (  halo->index[2*rank_id + end_shift]
                     - halo->index[2*rank_id]) * stride;

      if (halo->c_domain_rank[rank_id] != local_rank) {

        if (length > 0) {
          cs_lnum_t start = halo->index[2*rank_id]*elt_size;
          unsigned char *dest = _val_dest + start;
          MPI_Aint displacement = halo->c_domain_s_shift[rank_id]*elt_size;

          MPI_Get(dest,
                  length,                        /* origin count */
                  mpi_datatype,                  /* origin datatype */
                  halo->c_domain_rank[rank_id],  /* target rank */
                  displacement,                  /* target displacement */
                  length,                        /* target count */
                  mpi_datatype,                  /* target datatype */
                  hs->win);
        }

      }
      else
        hs->local_rank_id = rank_id;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize update of ghost values in case of parallelism
 *        for one-sided communication.
 *
 * The cs_halo_sync_pack function should have been called before this function,
 * using the same hs argument.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       val         pointer to variable value array
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

static void
_halo_sync_complete_one_sided(const cs_halo_t  *halo,
                              void             *val,
                              cs_halo_state_t  *hs)
{
  /* Use active target synchronization */

  /* Access Epoch */
  if (halo->c_domain_group != MPI_GROUP_NULL) {
    MPI_Win_complete(hs->win);

    /* Complete RMA exposure epoch */
    MPI_Win_wait(hs->win);
  }
  else {
    MPI_Win_fence(0, hs->win);
  }

  /* Copy local values in case of periodicity */

  if (hs->local_rank_id > -1) {
    size_t elt_size = cs_datatype_size[hs->data_type] * hs->stride;
    _sync_local(halo, hs->local_rank_id, hs->sync_mode, elt_size,
                hs->var_location, hs->send_buffer_cur, val);
  }
}

#endif /* (MPI_VERSION >= 3) */
#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a halo structure given an interface set.
 *
 * \param[in]  ifs  pointer to a cs_interface_set structure
 *
 * \return  pointer to created cs_halo_t structure
 */
/*----------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create(const cs_interface_set_t  *ifs)
{
  cs_lnum_t  i, tmp_id, perio_lst_size;

  cs_lnum_t  loc_id = -1;

  cs_halo_t  *halo = NULL;

  const cs_interface_t  *interface = NULL;

  BFT_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = cs_interface_set_size(ifs);
  halo->n_transforms = 0;

  halo->periodicity = cs_interface_set_periodicity(ifs);
  halo->n_rotations = 0;

  halo->n_local_elts = 0;

  for (i = 0; i < CS_HALO_N_TYPES; i++) {
    halo->n_send_elts[i] = 0;
    halo->n_elts [i] = 0;
  }

  BFT_MALLOC(halo->c_domain_rank, halo->n_c_domains, int);

  /* Check if cs_glob_rank_id belongs to interface set in order to
     order ranks with local rank at first place */

  for (i = 0; i < halo->n_c_domains; i++) {

    interface = cs_interface_set_get(ifs, i);
    halo->c_domain_rank[i] = cs_interface_rank(interface);

    if (cs_glob_rank_id == cs_interface_rank(interface))
      loc_id = i;

  } /* End of loop on ranks */

  if (loc_id > 0) {

    tmp_id = halo->c_domain_rank[loc_id];
    halo->c_domain_rank[loc_id] = halo->c_domain_rank[0];
    halo->c_domain_rank[0] = tmp_id;

  }

  /* Order ranks */

  if (   halo->n_c_domains > 2
      && _order_int_test(&(halo->c_domain_rank[1]),
                         halo->n_c_domains-1) == 0) {

    cs_lnum_t  *order = NULL;
    cs_gnum_t  *buffer = NULL;

    BFT_MALLOC(order, halo->n_c_domains - 1, cs_lnum_t);
    BFT_MALLOC(buffer, halo->n_c_domains - 1, cs_gnum_t);

    for (i = 1; i < halo->n_c_domains; i++)
      buffer[i-1] = (cs_gnum_t)halo->c_domain_rank[i];

    cs_order_gnum_allocated(NULL,
                            buffer,
                            order,
                            halo->n_c_domains - 1);

    for (i = 0; i < halo->n_c_domains - 1; i++)
      halo->c_domain_rank[i+1] = (cs_lnum_t)buffer[order[i]];

    BFT_FREE(buffer);
    BFT_FREE(order);

  } /* End of ordering ranks */

  CS_MALLOC_HD(halo->send_index, 2*halo->n_c_domains + 1, cs_lnum_t,
               _halo_buffer_alloc_mode);
  BFT_MALLOC(halo->index, 2*halo->n_c_domains + 1, cs_lnum_t);

  for (i = 0; i < 2*halo->n_c_domains + 1; i++) {
    halo->send_index[i] = 0;
    halo->index[i] = 0;
  }

  halo->send_perio_lst = NULL;
  halo->perio_lst = NULL;

  if (halo->periodicity != NULL) {

    halo->n_transforms = fvm_periodicity_get_n_transforms(halo->periodicity);

    for (i = 0; i < halo->n_transforms; i++) {
      if (   fvm_periodicity_get_type(halo->periodicity, i)
          >= FVM_PERIODICITY_ROTATION)
        halo->n_rotations += 1;
    }

    /* We need 2 values per transformation and there are n_transforms
       transformations. For each rank, we need a value for standard and
       extended halo. */

    perio_lst_size = 2*halo->n_transforms * 2*halo->n_c_domains;

    BFT_MALLOC(halo->send_perio_lst, perio_lst_size, cs_lnum_t);
    BFT_MALLOC(halo->perio_lst, perio_lst_size, cs_lnum_t);

    for (i = 0; i < perio_lst_size; i++) {
      halo->send_perio_lst[i] = 0;
      halo->perio_lst[i] = 0;
    }

  }

  halo->send_list = NULL;

#if defined(HAVE_MPI)
  halo->c_domain_group = MPI_GROUP_NULL;
  halo->c_domain_s_shift = NULL;
#endif

  _n_halos += 1;

  return halo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ready halo for use.
 *
 * This function should be called after building a halo using the
 * cs_halo_create_function and defined locally.
 * It is called automatically by cs_halo_create_from_ref and
 * cs_halo_create_from_rank_neighbors so does not need to be called again
 * using these functions.
 *
 * \param[in]  halo  pointer to halo structure
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_create_complete(cs_halo_t  *halo)
{
#if defined(HAVE_MPI)

  /* Make buffer available on device if relevant */
  cs_sync_h2d(halo->send_index);
  cs_sync_h2d(halo->send_list);

  /* Create group for one-sided communication */
  if (_halo_comm_mode > CS_HALO_COMM_P2P) {
    const int local_rank = CS_MAX(cs_glob_rank_id, 0);
    int n_group_ranks = 0;
    int *group_ranks = NULL;
    BFT_MALLOC(group_ranks, halo->n_c_domains + 1, int);
    for (int i = 0; i < halo->n_c_domains; i++) {
      if (halo->c_domain_rank[i] < local_rank)
        group_ranks[n_group_ranks++] = halo->c_domain_rank[i];
    }
    group_ranks[n_group_ranks++] = local_rank;
    for (int i = 0; i < halo->n_c_domains; i++) {
      if (halo->c_domain_rank[i] > local_rank)
        group_ranks[n_group_ranks++] = halo->c_domain_rank[i];
    }

    if (_order_int_test(group_ranks, n_group_ranks)) {

      MPI_Group glob_group;
      MPI_Comm_group(cs_glob_mpi_comm, &glob_group);
      MPI_Group_incl(glob_group,
                     n_group_ranks,
                     group_ranks,
                     &(halo->c_domain_group));
      MPI_Group_free(&glob_group);

    }

    BFT_FREE(group_ranks);
  }

  /* Exchange shifts for one-sided communication */
  if (_halo_comm_mode == CS_HALO_COMM_RMA_GET)
    _exchange_send_shift(halo);

  if (_halo_state == NULL)
    _halo_state = cs_halo_state_create();

#endif /* defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a halo structure, given a reference halo.
 *
 * \param[in]  ref  pointer to reference halo
 *
 * \return  pointer to created cs_halo_t structure
 */
/*----------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create_from_ref(const cs_halo_t  *ref)
{
  cs_lnum_t  i;

  cs_halo_t  *halo = NULL;

  BFT_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = ref->n_c_domains;
  halo->n_transforms = ref->n_transforms;

  halo->periodicity = ref->periodicity;
  halo->n_rotations = ref->n_rotations;

  halo->n_local_elts = 0;

  BFT_MALLOC(halo->c_domain_rank, halo->n_c_domains, int);

  for (i = 0; i < halo->n_c_domains; i++)
    halo->c_domain_rank[i] = ref->c_domain_rank[i];

  CS_MALLOC_HD(halo->send_index, 2*halo->n_c_domains + 1, cs_lnum_t,
               _halo_buffer_alloc_mode);
  BFT_MALLOC(halo->index, 2*halo->n_c_domains + 1, cs_lnum_t);

  for (i = 0; i < 2*halo->n_c_domains + 1; i++) {
    halo->send_index[i] = 0;
    halo->index[i] = 0;
  }

  halo->send_perio_lst = NULL;
  halo->perio_lst = NULL;

  if (halo->n_transforms > 0) {

    cs_lnum_t  perio_lst_size = 2*halo->n_transforms * 2*halo->n_c_domains;

    BFT_MALLOC(halo->send_perio_lst, perio_lst_size, cs_lnum_t);
    BFT_MALLOC(halo->perio_lst, perio_lst_size, cs_lnum_t);

    for (i = 0; i < perio_lst_size; i++) {
      halo->send_perio_lst[i] = 0;
      halo->perio_lst[i] = 0;
    }

  }

  halo->send_list = NULL;

#if defined(HAVE_MPI)
  halo->c_domain_group = MPI_GROUP_NULL;
  halo->c_domain_s_shift = NULL;
#endif

  _n_halos += 1;

  cs_halo_create_complete(halo);

  return halo;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a halo structure from distant element distant ranks and ids.
 *
 * \remark  This function does not handle periodicity. For most matrix-vector,
 *          products and similar operations, periodicity of translation an
 *          even rotation could be handled with no specific halo information,
 *          simply by assigning an equivalence between two periodic elements.
 *          For rotation, this would require also applying a rotation through
 *          the matrix coefficients (this would have the advantage of being
 *          compatible with external libraries). An alternative would be
 *          to add rotation information to a given halo as a second stage,
 *          through a specialized operator which can be added in the future.
 *
 * \param[in]  rn              associated rank neighbors info
 * \param[in]  n_local_elts    number of elements for local rank
 * \param[in]  n_distant_elts  number of distant elements for local rank
 * \param[in]  elt_rank_id     distant element rank index in rank neighbors,
 *                             ordered by rank (size: n_distant_elts)
 * \param[in]  elt_id          distant element id (at distant rank),
 *                             ordered by rank (size: n_distant_elts)
 *
 * \return  pointer to created cs_halo_t structure
 */
/*----------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create_from_rank_neighbors(const cs_rank_neighbors_t  *rn,
                                   cs_lnum_t                   n_local_elts,
                                   cs_lnum_t                   n_distant_elts,
                                   const int                   elt_rank_id[],
                                   const cs_lnum_t             elt_id[])
{
  cs_halo_t  *halo = NULL;

  BFT_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = 0;
  halo->n_transforms = 0;

  halo->n_rotations = 0;

  halo->periodicity = NULL;
  halo->send_perio_lst = NULL;
  halo->perio_lst = NULL;

#if defined(HAVE_MPI)
  halo->c_domain_group = MPI_GROUP_NULL;
  halo->c_domain_s_shift = NULL;
#endif

  halo->n_local_elts = n_local_elts;

  for (int i = 0; i < CS_HALO_N_TYPES; i++) {
    halo->n_send_elts[i] = 0;
    halo->n_elts [i] = n_distant_elts;
  }

  /* Count elements for each rank;
     check they are are ordered lexicographically */

  cs_lnum_t *rank_count;
  BFT_MALLOC(rank_count, rn->size*2, cs_lnum_t);
  for (int i = 0; i < rn->size; i++)
    rank_count[i] = 0;

  int rank_prev = -1;
  int elt_prev = -1;
  for (cs_lnum_t i = 0; i < n_distant_elts; i++) {
    int rank_id = elt_rank_id[i];
    if (   rank_id < rank_prev
        || (rank_id == rank_prev && elt_id[i] <= elt_prev))
      bft_error
        (__FILE__, __LINE__, 0,
         "%s:\n"
         "  Rank and distant element ids passed to this function must\n"
         "  be lexicographically ordered; this is not the case here.",
         __func__);
    rank_count[rank_id] += 1;
    rank_prev = rank_id;
    elt_prev = elt_id[i];
  }

  /* Now exchange counts with neighboring elements */

  MPI_Comm comm = cs_glob_mpi_comm;
  MPI_Request *request = NULL;
  MPI_Status *status = NULL;

  BFT_MALLOC(request, rn->size*2, MPI_Request);
  BFT_MALLOC(status, rn->size*2, MPI_Status);

  /* Exchange local range with neighbor ranks */

  int request_count = 0;
  const int local_rank = CS_MAX(cs_glob_rank_id, 0);

  for (int i = 0; i < rn->size; i++) {
    MPI_Irecv(rank_count + rn->size + i,
              1,
              CS_MPI_LNUM,
              rn->rank[i],
              local_rank,
              comm,
              &(request[request_count++]));
  }

  for (int i = 0; i < rn->size; i++) {
    MPI_Isend(rank_count + i,
              1,
              CS_MPI_LNUM,
              rn->rank[i],
              rn->rank[i],
              comm,
              &(request[request_count++]));
  }

  MPI_Waitall(request_count, request, status);

  /* Now build send and receive indexes to exchange data;
     the receive index can be directly assigned to the halo;
     also check if cs_glob_rank_id belongs to interface set in order to
     order ranks with local rank at first place */

  int        loc_r_index = -1;
  cs_lnum_t  r_displ = 0, loc_r_displ = 0;
  cs_lnum_t  recv_count = 0, send_count = 0;

  halo->n_c_domains = 0;
  for (int i = 0; i < rn->size; i++) {
    if (rank_count[i] + rank_count[rn->size + i] > 0) {
      halo->n_c_domains += 1;
      if (rn->rank[i] == local_rank) {
        loc_r_index = i;
        loc_r_displ = r_displ;
        assert(rank_count[i] == rank_count[rn->size + i]);
      }
      r_displ += rank_count[i];
      recv_count += rank_count[rn->size + i];
    }
  }

  BFT_MALLOC(halo->c_domain_rank, halo->n_c_domains, int);

  CS_MALLOC_HD(halo->send_list, recv_count, cs_lnum_t,
               _halo_buffer_alloc_mode);
  CS_MALLOC_HD(halo->send_index, 2*halo->n_c_domains + 1, cs_lnum_t,
               _halo_buffer_alloc_mode);
  BFT_MALLOC(halo->index, halo->n_c_domains*2+1, cs_lnum_t);

  halo->n_c_domains = 0;
  send_count = 0;
  recv_count = 0;

  halo->index[0] = 0;
  halo->send_index[0] = 0;

  if (loc_r_index > -1) {
    halo->c_domain_rank[0] = local_rank;
    cs_lnum_t  l_count = rank_count[loc_r_index];
    for (cs_lnum_t i = 0; i < l_count; i++)
      halo->send_list[i] = elt_id[loc_r_displ + i];
    send_count += l_count;
    recv_count += l_count;
    halo->n_c_domains = 1;
    for (int j = 1; j < 3; j++) {
      halo->index[j] = recv_count;
      halo->send_index[j] = send_count;
    }
  }

  for (int i = 0; i < rn->size; i++) {
    if (   rank_count[i] + rank_count[rn->size + i] > 0
        && rn->rank[i] != local_rank) {
      halo->c_domain_rank[halo->n_c_domains] = rn->rank[i];
      recv_count += rank_count[i];
      send_count += rank_count[rn->size + i];
      for (int j = 1; j < 3; j++) {
        halo->index[halo->n_c_domains*2 + j] = recv_count;
        halo->send_index[halo->n_c_domains*2 + j] = send_count;
      }
      halo->n_c_domains += 1;
    }
  }

  BFT_FREE(rank_count);

  for (int i = 0; i < CS_HALO_N_TYPES; i++)
    halo->n_send_elts[i] = send_count;

  /* Now send lists to matching ranks (reverse send and receive) */

  request_count = 0;

  for (int i = 0; i < halo->n_c_domains; i++) {
    int rank_id = halo->c_domain_rank[i];
    if (rank_id == local_rank) continue;
    cs_lnum_t r_shift = halo->send_index[2*i];
    cs_lnum_t r_size  = halo->send_index[2*i+1] - r_shift;
    if (r_size > 0)
      MPI_Irecv(halo->send_list + r_shift,
                r_size,
                CS_MPI_LNUM,
                rank_id,
                local_rank,
                comm,
                &(request[request_count++]));
  }

  for (int i = 0; i < halo->n_c_domains; i++) {
    int rank_id = halo->c_domain_rank[i];
    if (rank_id == local_rank) continue;
    cs_lnum_t s_shift = halo->index[2*i];
    cs_lnum_t s_size  = halo->index[2*i+1] - s_shift;
    if (s_shift < loc_r_displ) { /* case with local rank first */
      assert(halo->c_domain_rank[0] == local_rank);
      s_shift -= halo->index[2];
    }
    if (s_size > 0)
      MPI_Isend(elt_id + s_shift,
                s_size,
                CS_MPI_LNUM,
                rank_id,
                rank_id,
                comm,
                &(request[request_count++]));
  }

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  _n_halos += 1;

  cs_halo_create_complete(halo);

  return halo;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * brief Destroy a halo structure.
 *
 * \param[in, out]  halo  pointer to pointer to cs_halo structure to destroy.
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_destroy(cs_halo_t  **halo)
{
  if (halo == NULL)
    return;

  if (*halo == NULL)
    return;

  cs_halo_t  *_halo = *halo;

#if defined(HAVE_MPI)
  if (_halo->c_domain_group != MPI_GROUP_NULL)
    MPI_Group_free(&(_halo->c_domain_group));

  BFT_FREE(_halo->c_domain_s_shift);
#endif

  BFT_FREE(_halo->c_domain_rank);

  CS_FREE_HD(_halo->send_list);
  CS_FREE_HD(_halo->send_index);
  BFT_FREE(_halo->index);

  BFT_FREE(_halo->send_perio_lst);
  BFT_FREE(_halo->perio_lst);

  BFT_FREE(*halo);

  _n_halos -= 1;

  /* Delete default state if no halo remains */

  if (_n_halos == 0)
    cs_halo_state_destroy(&_halo_state);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a halo state structure.
 *
 * \return  pointer to created cs_halo_state_t structure.
 */
/*----------------------------------------------------------------------------*/

cs_halo_state_t *
cs_halo_state_create(void)
{
  cs_halo_state_t *hs;
  BFT_MALLOC(hs, 1, cs_halo_state_t);

  cs_halo_state_t hs_ini = {
    .sync_mode = CS_HALO_STANDARD,
    .data_type = CS_DATATYPE_NULL,
    .stride = 0,
    .var_location = CS_ALLOC_HOST,
    .send_buffer_cur = NULL,
    .n_requests = 0,
    .local_rank_id = -1,
    .send_buffer_size = 0,
    .recv_buffer_size = 0,
    .send_buffer = NULL,
    .recv_buffer = NULL
#if defined(HAVE_MPI)
    ,
    .request_size = 0,
    .request = NULL,
    .status = NULL,
    .win = MPI_WIN_NULL

#endif
  };

  *hs = hs_ini;

  return hs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a halo state structure.
 *
 * \param[in, out]  halo_state  pointer to pointer to cs_halo_state
 *                              structure to destroy.
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_state_destroy(cs_halo_state_t  **halo_state)
{
  if (halo_state != NULL) {
    cs_halo_state_t *hs = *halo_state;

#if defined(HAVE_MPI)
#if (MPI_VERSION >= 3)
    if (hs->win != MPI_WIN_NULL) {
      MPI_Win_free(&(hs->win));
      hs->win = MPI_WIN_NULL;
    }
#endif
#endif

    CS_FREE_HD(hs->send_buffer);
    CS_FREE_HD(hs->recv_buffer);

#if defined(HAVE_MPI)
    BFT_FREE(hs->request);
    BFT_FREE(hs->status);
#endif

    BFT_FREE(*halo_state);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to default halo state structure.
 *
 * \return]  halo  pointer to pointer to cs_halo structure to destroy.
 */
/*----------------------------------------------------------------------------*/

cs_halo_state_t *
cs_halo_state_get_default(void)
{
  return _halo_state;
}

/*----------------------------------------------------------------------------
 * Apply local cells renumbering to a halo
 *
 * parameters:
 *   halo        <-- pointer to halo structure
 *   new_cell_id <-- array indicating old -> new cell id (0 to n-1)
 *---------------------------------------------------------------------------*/

void
cs_halo_renumber_cells(cs_halo_t        *halo,
                       const cs_lnum_t   new_cell_id[])
{
  if (halo != NULL) {

    const cs_lnum_t n_elts = halo->n_send_elts[CS_HALO_EXTENDED];

    for (cs_lnum_t j = 0; j < n_elts; j++)
      halo->send_list[j] = new_cell_id[halo->send_list[j]];

    cs_sync_h2d(halo->send_list);

  }
}

/*----------------------------------------------------------------------------
 * Apply ghost cells renumbering to a halo
 *
 * parameters:
 *   halo        <-- pointer to halo structure
 *   old_cell_id <-- array indicating new -> old cell id (0 to n-1)
 *---------------------------------------------------------------------------*/

void
cs_halo_renumber_ghost_cells(cs_halo_t        *halo,
                             const cs_lnum_t   old_cell_id[])
{
  if (halo == NULL)
    return;

  /* Reverse update from distant cells */

  cs_lnum_t *send_buf, *recv_buf;

  BFT_MALLOC(send_buf, halo->n_send_elts[1], cs_lnum_t);
  BFT_MALLOC(recv_buf, halo->n_elts[1], cs_lnum_t);

  for (int i = 0; i < halo->n_c_domains; i++) {
    cs_lnum_t start = halo->index[2*i];
    cs_lnum_t end = halo->index[2*i+2];
    cs_lnum_t shift = halo->n_local_elts + halo->index[2*i];
    for (cs_lnum_t j = start; j < end; j++) {
      recv_buf[j] = old_cell_id[halo->n_local_elts + j] - shift;
      assert(recv_buf[j] >= 0 && recv_buf[j] < (end - start));
    }
  }

  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int rank_id;
    int request_count = 0;
    const int local_rank = cs_glob_rank_id;

    MPI_Request  *request;
    MPI_Status   *status;

    BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);

    /* Receive data from distant ranks */

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      cs_lnum_t start = halo->send_index[2*rank_id];
      cs_lnum_t length = (  halo->send_index[2*rank_id + 2]
                          - halo->send_index[2*rank_id]);

      if (halo->c_domain_rank[rank_id] != local_rank) {
        if (length > 0)
          MPI_Irecv(send_buf + start,
                    length,
                    CS_MPI_LNUM,
                    halo->c_domain_rank[rank_id],
                    local_rank,
                    cs_glob_mpi_comm,
                    &(request[request_count++]));
      }
      else
        local_rank_id = rank_id;

    }

    /* We wait for posting all receives (often recommended) */

    if (_halo_use_barrier)
      MPI_Barrier(cs_glob_mpi_comm);

    /* Send data to distant ranks */

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank_id] != local_rank) {

        cs_lnum_t start = halo->index[2*rank_id];
        cs_lnum_t length = (  halo->index[2*rank_id + 2]
                            - halo->index[2*rank_id]);

        if (length > 0)
          MPI_Isend(recv_buf + start,
                    length,
                    CS_MPI_LNUM,
                    halo->c_domain_rank[rank_id],
                    halo->c_domain_rank[rank_id],
                    cs_glob_mpi_comm,
                    &(request[request_count++]));

      }

    }

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    BFT_FREE(request);
    BFT_FREE(status);

  }

#endif /* defined(HAVE_MPI) */

  /* Copy local values if present */

  if (local_rank_id > -1) {

    cs_lnum_t *recv = recv_buf + halo->index[2*local_rank_id];

    cs_lnum_t start = halo->send_index[2*local_rank_id];
    cs_lnum_t length = (  halo->send_index[2*local_rank_id + 2]
                        - halo->send_index[2*local_rank_id]);

    for (cs_lnum_t j = 0; j < length; j++)
      send_buf[j+start] = recv[j];

  }

  BFT_FREE(recv_buf);

  /* Now apply renumbering to send list */

  for (int i = 0; i < halo->n_c_domains; i++) {
    cs_lnum_t start = halo->send_index[2*i];
    cs_lnum_t end = halo->send_index[2*i+2];
    for (cs_lnum_t j = start; j < end; j++)
      send_buf[j] = halo->send_list[start + send_buf[j]];
    for (cs_lnum_t j = start; j < end; j++)
      halo->send_list[j] = send_buf[j];
  }

  cs_sync_h2d(halo->send_list);

  BFT_FREE(send_buf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize halo state prior to packing halo data to send.
 *
 * A local state handler may be provided, or the default state handler will
 * be used.
 *
 * This function is included in \ref cs_halo_sync_pack, but may be called
 * separately for specific implementations, such as for accelerator devices.
 *
 * A local state and/or buffer may be provided, or the default (global) state
 * and buffer will be used. If provided explicitely,
 * the buffer must be of sufficient size.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       data_type   data type
 * \param[in]       stride      number of (interlaced) values by entity
 * \param[out]      send_buf    pointer to send buffer, NULL for global buffer
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 *
 * \return  pointer to halo send buffer
 */
/*----------------------------------------------------------------------------*/

void *
cs_halo_sync_pack_init_state(const cs_halo_t  *halo,
                             cs_halo_type_t    sync_mode,
                             cs_datatype_t     data_type,
                             int               stride,
                             void             *send_buf,
                             cs_halo_state_t  *hs)
{
  void *_send_buffer = send_buf;

  if (halo == NULL)
    return _send_buffer;

  cs_halo_state_t  *_hs = (hs != NULL) ? hs : _halo_state;

  if (_send_buffer == NULL) {
    size_t send_buffer_size = cs_halo_pack_size(halo, data_type, stride);

    if (send_buffer_size > _hs->send_buffer_size) {
      cs_alloc_mode_t alloc_mode = cs_check_device_ptr(halo->send_list);

      _hs->send_buffer_size = send_buffer_size;

#if defined(HAVE_MPI)
#if (MPI_VERSION >= 3)
      if (_hs->win != MPI_WIN_NULL) {
        MPI_Win_free(&(_hs->win));
        _hs->win = MPI_WIN_NULL;
      }
#endif
#endif

      CS_FREE_HD(_hs->send_buffer);
      CS_MALLOC_HD(_hs->send_buffer,
                   _hs->send_buffer_size,
                   char,
                   alloc_mode);

#if defined(HAVE_MPI)
#if (MPI_VERSION >= 3)
      if (_halo_comm_mode == CS_HALO_COMM_RMA_GET)
        MPI_Win_create(_hs->send_buffer,
                       _hs->send_buffer_size,
                       1,   /* displacement unit */
                       MPI_INFO_NULL,
                       MPI_COMM_WORLD,
                       &(_hs->win));
#endif
#endif
    }

    _send_buffer = _hs->send_buffer;
  }

  _hs->var_location = CS_ALLOC_HOST;
  _hs->send_buffer_cur = _send_buffer;

  _hs->sync_mode = sync_mode;
  _hs->data_type = data_type;
  _hs->stride = stride;

  return _send_buffer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Pack halo data to send into dense buffer.
 *
 * A local state handler may be provided, or the default state handler will
 * be used.
 *
 * A local state and/or buffer may be provided, or the default (global) state
 * and buffer will be used. If provided explicitely,
 * the buffer must be of sufficient size.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       data_type   data type
 * \param[in]       stride      number of (interlaced) values by entity
 * \param[in]       val         pointer to variable value array
 * \param[out]      send_buf    pointer to send buffer, NULL for global buffer
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_pack(const cs_halo_t  *halo,
                  cs_halo_type_t    sync_mode,
                  cs_datatype_t     data_type,
                  int               stride,
                  void             *val,
                  void             *send_buf,
                  cs_halo_state_t  *hs)
{
  if (halo == NULL)
    return;

  void *_send_buffer = cs_halo_sync_pack_init_state(halo,
                                                    sync_mode,
                                                    data_type,
                                                    stride,
                                                    send_buf,
                                                    hs);

  cs_lnum_t end_shift = 0;

  if (sync_mode == CS_HALO_STANDARD)
    end_shift = 1;

  else if (sync_mode == CS_HALO_EXTENDED)
    end_shift = 2;

  /* Assemble buffers for halo exchange; avoid threading for now, as dynamic
     scheduling led to slightly higher cost here in some tests,
     and even static scheduling might lead to false sharing for small
     halos. */

  if (data_type == CS_REAL_TYPE) {

    cs_real_t *buffer = (cs_real_t *)_send_buffer;
    cs_real_t *var = val;

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      cs_lnum_t p_start = halo->send_index[2*rank_id]*stride;
      size_t start = halo->send_index[2*rank_id];
      size_t length = (  halo->send_index[2*rank_id + end_shift]
                       - halo->send_index[2*rank_id]);

      if (stride == 3) { /* Unroll loop for this case */
        for (size_t i = 0; i < length; i++) {
          buffer[p_start + i*3]
            = var[(halo->send_list[start + i])*3];
          buffer[p_start + i*3 + 1]
            = var[(halo->send_list[start + i])*3 + 1];
          buffer[p_start + i*3 + 2]
            = var[(halo->send_list[start + i])*3 + 2];
        }
      }
      else {
        size_t _stride = stride;
        for (size_t i = 0; i < length; i++) {
          size_t r_start = halo->send_list[start + i] * stride;
          for (size_t j = 0; j < _stride; j++)
            buffer[p_start + i*_stride + j] = var[r_start + j];
        }
      }

    }

  }

  else {

    unsigned char *buffer = (unsigned char *)_send_buffer;

    size_t elt_size = cs_datatype_size[data_type] * stride;

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      cs_lnum_t p_start = halo->send_index[2*rank_id]*elt_size;
      size_t start = halo->send_index[2*rank_id];
      size_t length = (  halo->send_index[2*rank_id + end_shift]
                       - halo->send_index[2*rank_id]);

      unsigned char *restrict _val = val;
      unsigned char *_buffer = buffer + p_start;

      for (size_t i = 0; i < length; i++) {
        size_t r_start = halo->send_list[start + i] * elt_size;
        for (size_t j = 0; j < elt_size; j++)
          _buffer[i*elt_size + j] = _val[r_start + j];
      }

    }

  }
}

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Pack halo data to send into dense buffer on accelerator device.
 *
 * A local state handler may be provided, or the default state handler will
 * be used.
 *
 * A local state and/or buffer may be provided, or the default (global) state
 * and buffer will be used. If provided explicitely,
 * the buffer must be of sufficient size.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       data_type   data type
 * \param[in]       stride      number of (interlaced) values by entity
 * \param[in]       val         pointer to variable value array (on device)
 * \param[out]      send_buf    pointer to send buffer (on device),
 *                              NULL for global buffer
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_pack_d(const cs_halo_t  *halo,
                    cs_halo_type_t    sync_mode,
                    cs_datatype_t     data_type,
                    int               stride,
                    void             *val,
                    void             *send_buf,
                    cs_halo_state_t  *hs)
{
  if (halo == NULL)
    return;

  cs_halo_state_t  *_hs = (hs != NULL) ? hs : _halo_state;

  void *_send_buf = cs_halo_sync_pack_init_state(halo,
                                                 sync_mode,
                                                 data_type,
                                                 stride,
                                                 send_buf,
                                                 _hs);

  void *_send_buf_d = cs_get_device_ptr(_send_buf);

#if defined(HAVE_CUDA)

  cs_halo_cuda_pack_send_buffer_real(halo,
                                     sync_mode,
                                     stride,
                                     val,
                                     _send_buf_d);

#else

  cs_halo_sync_pack(halo,
                    sync_mode,
                    data_type,
                    stride,
                    val,
                    send_buf,
                    _hs);

#endif

  /* As device pointer is passed, cs_check_device_ptr will provide
     the correct result only in case of CS_ALLOC_HOST_DEVICE_SHARED
     (where the pointers are identical). */

  _hs->var_location = cs_check_device_ptr(val);
  if (_hs->var_location != CS_ALLOC_HOST_DEVICE_SHARED)
    _hs->var_location = CS_ALLOC_DEVICE;
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Launch update array of values in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * The cs_halo_sync_pack function should have been called before this function,
 * using the same hs argument.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       val         pointer to variable value array
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_start(const cs_halo_t  *halo,
                   void             *val,
                   cs_halo_state_t  *hs)
{
  if (halo == NULL)
    return;

  cs_halo_state_t  *_hs = (hs != NULL) ? hs : _halo_state;

#if (MPI_VERSION >= 3)
  if (_halo_comm_mode > CS_HALO_COMM_P2P) {
    _halo_sync_start_one_sided(halo, val, _hs);
    return;
  }
#endif

  cs_lnum_t end_shift = (_hs->sync_mode == CS_HALO_EXTENDED) ? 2 : 1;
  cs_lnum_t stride = _hs->stride;
  size_t elt_size = cs_datatype_size[_hs->data_type] * stride;
  size_t n_loc_elts = halo->n_local_elts;

  unsigned char *restrict _val = val;
  unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

  unsigned char *buffer = (unsigned char *)(_hs->send_buffer_cur);

  if (_hs->var_location > CS_ALLOC_HOST) {
#   if defined(_CS_MPI_DEVICE_SUPPORT)
    /* For CUDA-aware MPI, directly work with buffer on device */
    buffer = cs_get_device_ptr(buffer);
# else
    /* For host-based MPI, copy or prefetch buffer */
    cs_sync_d2h(buffer);

    /* When array passed is defined on device but is not shared, use separate
       (smaller) CPU buffer for receive (as we cannot know whether a matching
       host array is present without complexifying the API);
       this will be copied back to device at the next step */
    if (_hs->var_location != CS_ALLOC_HOST_DEVICE_SHARED) {
      size_t recv_size = halo->n_elts[_hs->sync_mode] * elt_size;
      if (_hs->recv_buffer_size < recv_size) {
        _hs->recv_buffer_size = recv_size;
        CS_FREE_HD(_hs->recv_buffer);
        CS_MALLOC_HD(_hs->recv_buffer, _hs->recv_buffer_size, unsigned char,
                     CS_ALLOC_HOST_DEVICE_PINNED);
      }
      _val_dest = _hs->recv_buffer;
    }
#endif
  }

#if defined(HAVE_MPI)

  _update_requests(halo, _hs);

  MPI_Datatype mpi_datatype = cs_datatype_to_mpi[_hs->data_type];

  int request_count = 0;
  const int local_rank = CS_MAX(cs_glob_rank_id, 0);

  /* Receive data from distant ranks */

  for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    cs_lnum_t length = (  halo->index[2*rank_id + end_shift]
                        - halo->index[2*rank_id]) * stride;

    if (halo->c_domain_rank[rank_id] != local_rank) {

      if (length > 0) {
        size_t start = (size_t)(halo->index[2*rank_id]);
        unsigned char *dest = _val_dest + start*elt_size;

        MPI_Irecv(dest,
                  length*_hs->stride,
                  mpi_datatype,
                  halo->c_domain_rank[rank_id],
                  halo->c_domain_rank[rank_id],
                  cs_glob_mpi_comm,
                  &(_hs->request[request_count++]));
      }

    }
    else
      _hs->local_rank_id = rank_id;
  }

  /* We may wait for posting all receives (sometimes recommended) */

  if (_halo_use_barrier)
    MPI_Barrier(cs_glob_mpi_comm);

  /* Send data to distant ranks */

  for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    cs_lnum_t start = halo->send_index[2*rank_id]*elt_size;
    cs_lnum_t length = (  halo->send_index[2*rank_id + end_shift]
                        - halo->send_index[2*rank_id]);

    if (halo->c_domain_rank[rank_id] != local_rank && length > 0)
      MPI_Isend(buffer + start,
                length*stride,
                mpi_datatype,
                halo->c_domain_rank[rank_id],
                local_rank,
                cs_glob_mpi_comm,
                &(_hs->request[request_count++]));

  }

  _hs->n_requests = request_count;

#endif /* defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wait for completion of update array of values in case of
 *  parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * The cs_halo_sync_start function should have been called before this function,
 * using the same hs argument.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       val         pointer to variable value array
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_wait(const cs_halo_t  *halo,
                  void             *val,
                  cs_halo_state_t  *hs)
{
  if (halo == NULL)
    return;

  cs_halo_state_t  *_hs = (hs != NULL) ? hs : _halo_state;

#if (MPI_VERSION >= 3)
  if (_halo_comm_mode > CS_HALO_COMM_P2P) {
    _halo_sync_complete_one_sided(halo, val, _hs);
    return;
  }
#endif

#if defined(HAVE_MPI)

  /* Wait for all exchanges */

  if (_hs->n_requests > 0)
    MPI_Waitall(_hs->n_requests, _hs->request, _hs->status);

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_ACCEL)
#if !defined(_CS_MPI_DEVICE_SUPPORT)
  if (_hs->var_location > CS_ALLOC_HOST) {

    size_t n_loc_elts = halo->n_local_elts;
    size_t n_elts = (   _hs->sync_mode
                     == CS_HALO_EXTENDED) ? halo->n_elts[1] : halo->n_elts[0];
    size_t elt_size = cs_datatype_size[_hs->data_type] * _hs->stride;
    size_t n_bytes = n_elts*elt_size;

    unsigned char *restrict _val = val;
    unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

    if (_hs->var_location == CS_ALLOC_HOST_DEVICE_SHARED)
      cs_prefetch_h2d(_val_dest, n_bytes);
    else
      cs_copy_h2d(_val_dest, _hs->recv_buffer, n_bytes);

  }
#endif
#endif /* defined(HAVE_ACCEL) */

  /* Copy local values in case of periodicity */

  if (_hs->local_rank_id > -1) {
    size_t elt_size = cs_datatype_size[_hs->data_type] * _hs->stride;
    _sync_local(halo, _hs->local_rank_id, _hs->sync_mode, elt_size,
                _hs->var_location, _hs->send_buffer_cur, val);
  }

  /* Cleanup */

  _hs->sync_mode = CS_HALO_STANDARD;
  _hs->data_type = CS_DATATYPE_NULL;
  _hs->stride = 0;
  _hs->send_buffer_cur = NULL;
  _hs->n_requests = 0;
  _hs->local_rank_id  = -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update array of values in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \param[in]   halo        pointer to halo structure
 * \param[in]   sync_mode   synchronization mode (standard or extended)
 * \param[in]   data_type   data type
 * \param[in]   stride      number of (interlaced) values by entity
 * \param[in]   val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync(const cs_halo_t  *halo,
             cs_halo_type_t    sync_mode,
             cs_datatype_t     data_type,
             int               stride,
             void             *val)
{
  if (halo == NULL)
    return;

  cs_halo_sync_pack(halo,
                    sync_mode,
                    data_type,
                    stride,
                    val,
                    NULL,
                    NULL);

  cs_halo_sync_start(halo, val, NULL);

  cs_halo_sync_wait(halo, val, NULL);
}

/*----------------------------------------------------------------------------
 * Update array of any type of halo values in case of parallelism or
 * periodicity.
 *
 * Data is untyped; only its size is given, so this function may also
 * be used to synchronize interleaved multidimendsional data, using
 * size = element_size*dim (assuming a homogeneous environment, at least
 * as far as data encoding goes).
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   size      <-- datatype size
 *   num       <-> pointer to local number value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_untyped(const cs_halo_t  *halo,
                     cs_halo_type_t    sync_mode,
                     size_t            size,
                     void             *val)
{
  cs_halo_sync(halo, sync_mode, CS_CHAR, size, val);
}

/*----------------------------------------------------------------------------
 * Update array of integer halo values in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   num       <-> pointer to local number value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_num(const cs_halo_t  *halo,
                 cs_halo_type_t    sync_mode,
                 cs_lnum_t         num[])
{
  cs_halo_sync(halo, sync_mode, CS_LNUM_TYPE, 1, num);
}

/*----------------------------------------------------------------------------
 * Update array of variable (floating-point) halo values in case of
 * parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   var       <-> pointer to variable value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_var(const cs_halo_t  *halo,
                 cs_halo_type_t    sync_mode,
                 cs_real_t         var[])
{
  cs_halo_sync(halo, sync_mode, CS_REAL_TYPE, 1, var);
}

/*----------------------------------------------------------------------------
 * Update array of strided variable (floating-point) values in case
 * of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   var       <-> pointer to variable value array
 *   stride    <-- number of (interlaced) values by entity
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_var_strided(const cs_halo_t  *halo,
                         cs_halo_type_t    sync_mode,
                         cs_real_t         var[],
                         int               stride)
{
  cs_halo_sync(halo, sync_mode, CS_REAL_TYPE, stride, var);
}

/*----------------------------------------------------------------------------
 * Return MPI_Barrier usage flag.
 *
 * returns:
 *   true if MPI barriers are used after posting receives and before posting
 *   sends, false otherwise
 *---------------------------------------------------------------------------*/

bool
cs_halo_get_use_barrier(void)
{
  return _halo_use_barrier;
}

/*----------------------------------------------------------------------------
 * Set MPI_Barrier usage flag.
 *
 * parameters:
 *   use_barrier <-- true if MPI barriers should be used after posting
 *                   receives and before posting sends, false otherwise.
 *---------------------------------------------------------------------------*/

void
cs_halo_set_use_barrier(bool use_barrier)
{
  _halo_use_barrier = use_barrier;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get default communication mode for halo exchange.
 *
 * \return  allocation mode
 */
/*----------------------------------------------------------------------------*/

cs_halo_comm_mode_t
cs_halo_get_comm_mode(void)
{
  return _halo_comm_mode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set default communication mode for halo exchange.
 *
 * \param[in]  mode  allocation mode to set
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_set_comm_mode(cs_halo_comm_mode_t  mode)
{
  if (mode >= CS_HALO_COMM_P2P && mode <= CS_HALO_COMM_RMA_GET)
    _halo_comm_mode = mode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get default host/device allocation mode for message packing arrays.
 *
 * \return  allocation mode
 */
/*----------------------------------------------------------------------------*/

cs_alloc_mode_t
cs_halo_get_buffer_alloc_mode(void)
{
  return _halo_buffer_alloc_mode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set default host/device allocation mode for message packing arrays.
 *
 * \param[in]  mode  allocation mode to set
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_set_buffer_alloc_mode(cs_alloc_mode_t  mode)
{
  _halo_buffer_alloc_mode = mode;
}

/*----------------------------------------------------------------------------
 * Dump a cs_halo_t structure.
 *
 * parameters:
 *   halo           <-- pointer to cs_halo_t struture
 *   print_level    <--  0 only dimensions and indexes are printed, else (1)
 *                       everything is printed
 *---------------------------------------------------------------------------*/

void
cs_halo_dump(const cs_halo_t  *halo,
             int               print_level)
{
  if (halo == NULL) {
    bft_printf("\n\n  halo: nil\n");
    return;
  }

  bft_printf("\n  halo:         %p\n"
             "  n_transforms:   %d\n"
             "  n_c_domains:    %d\n"
             "  periodicity:    %p\n"
             "  n_rotations:    %d\n"
             "  n_local_elts:   %ld\n",
             (const void *)halo,
             halo->n_transforms, halo->n_c_domains,
             (const void *)halo->periodicity,
             halo->n_rotations, (long)halo->n_local_elts);

  bft_printf("\nRanks on halo frontier:\n");
  for (int i = 0; i < halo->n_c_domains; i++)
    bft_printf("%5d", halo->c_domain_rank[i]);

  for (int halo_id = 0; halo_id < 2; halo_id++) {

    cs_lnum_t  n_elts[2];
    cs_lnum_t  *index = NULL, *list = NULL, *perio_lst = NULL;

    bft_printf("\n    ---------\n");

    if (halo_id == 0) {

      bft_printf("    send_list:\n");
      n_elts[0] = halo->n_send_elts[0];
      n_elts[1] = halo->n_send_elts[1];
      index = halo->send_index;
      list = halo->send_list;
      perio_lst = halo->send_perio_lst;

    }
    else if (halo_id == 1) {

      bft_printf("    halo:\n");
      n_elts[0] = halo->n_elts[0];
      n_elts[1] = halo->n_elts[1];
      index = halo->index;
      list = NULL;
      perio_lst = halo->perio_lst;

    }

    bft_printf("    ---------\n\n");
    bft_printf("  n_ghost_cells:        %ld\n"
               "  n_std_ghost_cells:    %ld\n", (long)n_elts[1], (long)n_elts[0]);

    if (index == NULL)
      return;

    if (halo->n_transforms > 0) {

      const cs_lnum_t  stride = 4*halo->n_c_domains;

      for (int i = 0; i < halo->n_transforms; i++) {

        bft_printf("\nTransformation number: %d\n", i+1);

        for (int j = 0; j < halo->n_c_domains; j++) {

          bft_printf("    rank %3d <STD> %5ld %5ld <EXT> %5ld %5ld\n",
                     halo->c_domain_rank[j],
                     (long)perio_lst[i*stride + 4*j],
                     (long)perio_lst[i*stride + 4*j+1],
                     (long)perio_lst[i*stride + 4*j+2],
                     (long)perio_lst[i*stride + 4*j+3]);
        }

      } /* End of loop on perio */

    } /* End if n_perio > 0 */

    for (int i = 0; i < halo->n_c_domains; i++) {

      bft_printf("\n  rank      %d:\n", halo->c_domain_rank[i]);

      if (index[2*i+1] - index[2*i] > 0) {

        bft_printf("\n  Standard halo\n");
        bft_printf("  idx start %ld:          idx end   %ld:\n",
                   (long)index[2*i], (long)index[2*i+1]);

        if (print_level > 0 && list != NULL) {
          bft_printf("\n            idx     elt id\n");
          for (cs_lnum_t j = index[2*i]; j < index[2*i+1]; j++)
            bft_printf("    %10ld %10ld\n", (long)j, (long)list[j]);
        }

      } /* there are elements on standard neighborhood */

      if (index[2*i+2] - index[2*i+1] > 0) {

        bft_printf("\n  Extended halo\n");
        bft_printf("  idx start %ld:          idx end   %ld:\n",
                   (long)index[2*i+1], (long)index[2*i+2]);

        if (print_level > 0 && list != NULL) {
          bft_printf("\n            idx     elt id\n");
          for (long j = index[2*i+1]; j < index[2*i+2]; j++)
            bft_printf("    %10ld %10ld %10ld\n",
                       (long)j, (long)list[j], (long)halo->n_local_elts+j);
        }

      } /* If there are elements on extended neighborhood */

    } /* End of loop on involved ranks */

  } /* End of loop on halos (send_halo/halo) */

  bft_printf("\n\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
