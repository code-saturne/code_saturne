/*============================================================================
 * Functions dealing with ghost cells
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

#include <stdio.h>
#include <string.h>
#include <assert.h>

#if defined(HAVE_NCCL)
#include <nccl.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_profiling.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_base_accel.h"
#if defined(HAVE_CUDA)
#include "base/cs_base_cuda.h"
#endif
#include "base/cs_interface.h"
#include "base/cs_mem.h"
#include "base/cs_order.h"
#include "base/cs_rank_neighbors.h"

#include "fvm/fvm_periodicity.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"

#if defined(HAVE_CUDA)
#include "base/cs_halo_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

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
  cs_alloc_mode_t send_buffer_location;  /* Allocation mode for send buffer */

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

#if defined(HAVE_CUDA)

  cudaStream_t  stream;           /* Associated CUDA stream */

#endif
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of defined halos */
static int _n_halos = 0;

/* Allocation mode for arrays which might be used on accelerator device.
   Best performance is expected either with page-locked (pinned) host memory
   and separate device memory, or shared (managed) memory with prefetching.
   We should run performance comparisons, but in the case of similar
   performance, going for the shared approach would be preferred for its
   other advantages (simplicity most of all, and pageable device memory). */

#if defined(HAVE_ACCEL)
static cs_alloc_mode_t _halo_buffer_alloc_mode = CS_ALLOC_HOST_DEVICE_PINNED;
#else
static cs_alloc_mode_t _halo_buffer_alloc_mode = CS_ALLOC_HOST;
#endif

/* Should we use barriers after posting receives ? */
static int _halo_use_barrier = false;

/* Default halo state handler */
static cs_halo_state_t *_halo_state = nullptr;

/* Halo communications mode */
static cs_halo_comm_mode_t _halo_comm_mode = CS_HALO_COMM_P2P;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if an array of global numbers is ordered.
 *
 * \param[in]  list    optional list (1 to n numbering) of selected entities
 *                     (or null if all nb_ent are selected). This list may
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

  if (list != nullptr) {
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
  if (halo == nullptr)
    return;

  int n_requests = halo->n_c_domains*2;

  if (n_requests > hs->request_size) {
    hs->request_size = n_requests;
    CS_REALLOC(hs->request, hs->request_size, MPI_Request);
    CS_REALLOC(hs->status, hs->request_size,  MPI_Status);
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
  MPI_Request *request = nullptr;
  MPI_Status *status = nullptr;

  CS_MALLOC(request, halo->n_c_domains*2, MPI_Request);
  CS_MALLOC(status, halo->n_c_domains*2, MPI_Status);

  CS_REALLOC(halo->c_domain_s_shift, halo->n_c_domains, cs_lnum_t);

  /* Exchange local range with neighbor ranks */

  const int local_rank = cs::max(cs_glob_rank_id, 0);

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

  CS_FREE(request);
  CS_FREE(status);
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

  unsigned char *_val = (unsigned char *)val;
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
 * \param[in, out]  hs          pointer to halo state, null for global state
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

  unsigned char *restrict _val = (unsigned char *)val;
  unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

  MPI_Datatype mpi_datatype = cs_datatype_to_mpi[hs->data_type];

  const int local_rank = cs::max(cs_glob_rank_id, 0);

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
 * \param[in, out]  hs          pointer to halo state, null for global state
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

#if defined(HAVE_NCCL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Launch update array of values with NCCL in case of parallelism.
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
 * \param[in, out]  hs          pointer to halo state
 */
/*----------------------------------------------------------------------------*/

static void
_sync_start_nccl(const cs_halo_t  *halo,
                 void             *val,
                 cs_halo_state_t  *hs)
{
  if (halo == nullptr)
    return;

  cs_lnum_t end_shift = (hs->sync_mode == CS_HALO_EXTENDED) ? 2 : 1;
  cs_lnum_t stride = hs->stride;
  size_t elt_size = cs_datatype_size[hs->data_type] * stride;
  size_t n_loc_elts = halo->n_local_elts;

  unsigned char *restrict _val = (unsigned char *)val;
  unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

  unsigned char *buffer = (unsigned char *)(hs->send_buffer_cur);
  buffer = (unsigned char *)cs_get_device_ptr(buffer);

  ncclDataType_t nccl_datatype = cs_datatype_to_nccl[hs->data_type];

  const int local_rank = cs::max(cs_glob_rank_id, 0);

  ncclGroupStart();

  /* Receive data from distant ranks */

  for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    cs_lnum_t length = (  halo->index[2*rank_id + end_shift]
                        - halo->index[2*rank_id]);

    if (halo->c_domain_rank[rank_id] != local_rank) {

      if (length > 0) {
        size_t start = (size_t)(halo->index[2*rank_id]);
        unsigned char *dest = _val_dest + start*elt_size;

        ncclRecv(dest,
                 length*stride,
                 nccl_datatype,
                 halo->c_domain_rank[rank_id],
                 cs_glob_nccl_comm,
                 hs->stream);
      }

    }
    else
      hs->local_rank_id = rank_id;
  }

  /* Send data to distant ranks */

  for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    cs_lnum_t start = halo->send_index[2*rank_id]*elt_size;
    cs_lnum_t length = (  halo->send_index[2*rank_id + end_shift]
                        - halo->send_index[2*rank_id]);

    if (halo->c_domain_rank[rank_id] != local_rank && length > 0)
      ncclSend(buffer + start,
               length*stride,
               nccl_datatype,
               halo->c_domain_rank[rank_id],
               cs_glob_nccl_comm,
               hs->stream);

  }

  ncclGroupEnd();
}

#endif // defined(HAVE_NCCL)

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public C++ function definitions
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

  cs_halo_t  *halo = nullptr;

  const cs_interface_t  *interface = nullptr;

  CS_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = cs_interface_set_size(ifs);
  halo->n_transforms = 0;

  halo->periodicity = cs_interface_set_periodicity(ifs);
  halo->n_rotations = 0;

  halo->n_local_elts = 0;

  for (i = 0; i < CS_HALO_N_TYPES; i++) {
    halo->n_send_elts[i] = 0;
    halo->n_elts [i] = 0;
  }

  CS_MALLOC(halo->c_domain_rank, halo->n_c_domains, int);

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

    cs_lnum_t  *order = nullptr;
    cs_gnum_t  *buffer = nullptr;

    CS_MALLOC(order, halo->n_c_domains - 1, cs_lnum_t);
    CS_MALLOC(buffer, halo->n_c_domains - 1, cs_gnum_t);

    for (i = 1; i < halo->n_c_domains; i++)
      buffer[i-1] = (cs_gnum_t)halo->c_domain_rank[i];

    cs_order_gnum_allocated(nullptr,
                            buffer,
                            order,
                            halo->n_c_domains - 1);

    for (i = 0; i < halo->n_c_domains - 1; i++)
      halo->c_domain_rank[i+1] = (cs_lnum_t)buffer[order[i]];

    CS_FREE(buffer);
    CS_FREE(order);

  } /* End of ordering ranks */

  CS_MALLOC_HD(halo->send_index, 2*halo->n_c_domains + 1, cs_lnum_t,
               _halo_buffer_alloc_mode);
  CS_MALLOC(halo->index, 2*halo->n_c_domains + 1, cs_lnum_t);

  for (i = 0; i < 2*halo->n_c_domains + 1; i++) {
    halo->send_index[i] = 0;
    halo->index[i] = 0;
  }

  halo->send_perio_lst = nullptr;
  halo->perio_lst = nullptr;

  if (halo->periodicity != nullptr) {

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

    CS_MALLOC(halo->send_perio_lst, perio_lst_size, cs_lnum_t);
    CS_MALLOC(halo->perio_lst, perio_lst_size, cs_lnum_t);

    for (i = 0; i < perio_lst_size; i++) {
      halo->send_perio_lst[i] = 0;
      halo->perio_lst[i] = 0;
    }

  }

  halo->send_list = nullptr;

  halo->std_send_block_size = 256;  /* Fixed size for now */
  halo->n_std_send_blocks   = 0;
  halo->std_send_blocks = nullptr;

#if defined(HAVE_MPI)
  halo->c_domain_group = MPI_GROUP_NULL;
  halo->c_domain_s_shift = nullptr;
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
  cs_alloc_mode_t halo_alloc_mode = cs_check_device_ptr(halo->send_index);
  if (halo_alloc_mode > CS_ALLOC_HOST) {
    cs_sync_h2d_start(halo->send_index);
    cs_sync_h2d_start(halo->send_list);
    if (halo_alloc_mode == CS_ALLOC_HOST_DEVICE_PINNED)
      cs_mem_hd_async_wait();
  }

  /* Create group for one-sided communication */
  if (_halo_comm_mode > CS_HALO_COMM_P2P) {
    const int local_rank = cs::max(cs_glob_rank_id, 0);
    int n_group_ranks = 0;
    int *group_ranks = nullptr;
    CS_MALLOC(group_ranks, halo->n_c_domains + 1, int);
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

    CS_FREE(group_ranks);
  }

  /* Exchange shifts for one-sided communication */
  if (_halo_comm_mode == CS_HALO_COMM_RMA_GET)
    _exchange_send_shift(halo);

  /* Build block send info for packing of standard exchange
     (not needed if the halo only has a standard component,
     as we can use the whole send list in that case) */

  if (halo->n_send_elts[1] > halo->n_send_elts[0]) {
    const cs_lnum_t block_size = halo->std_send_block_size;

    cs_lnum_t n_blocks_ref
      = (halo->n_send_elts[1] % block_size)
          ? halo->n_send_elts[1]/block_size + 1
          : halo->n_send_elts[1]/block_size;
    cs_lnum_t n_blocks = 0;

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
      size_t n = (  halo->send_index[2*rank_id + 1]
                  - halo->send_index[2*rank_id]);
      n_blocks +=  (n % block_size) ?  n/block_size + 1 : n/block_size;
    }

    if (n_blocks < n_blocks_ref) {
      halo->n_std_send_blocks = n_blocks;

      cs_lnum_t *send_blocks;
      CS_MALLOC_HD(send_blocks, halo->n_std_send_blocks*2, cs_lnum_t,
                   cs_check_device_ptr(halo->send_list));

      n_blocks = 0;

      for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        cs_lnum_t e_id = halo->send_index[2*rank_id+1];

        for (cs_lnum_t s_id = halo->send_index[2*rank_id];
             s_id < e_id;
             s_id += block_size) {
          send_blocks[n_blocks*2] = s_id;
          send_blocks[n_blocks*2+1] = cs::min(s_id + block_size, e_id);
          n_blocks += 1;
        }

      }

      halo->std_send_blocks = send_blocks;
    }
  }

#endif /* defined(HAVE_MPI) */

  if (_halo_state == nullptr)
    _halo_state = cs_halo_state_create();
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
  cs_halo_t  *halo = nullptr;

  cs_alloc_mode_t halo_alloc_mode = cs_check_device_ptr(ref->send_index);

  CS_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = ref->n_c_domains;
  halo->n_transforms = ref->n_transforms;

  halo->periodicity = ref->periodicity;
  halo->n_rotations = ref->n_rotations;

  halo->n_local_elts = 0;
  halo->n_send_elts[0] = 0;
  halo->n_send_elts[1] = 0;

  CS_MALLOC(halo->c_domain_rank, halo->n_c_domains, int);

  for (int i = 0; i < halo->n_c_domains; i++)
    halo->c_domain_rank[i] = ref->c_domain_rank[i];

  CS_MALLOC_HD(halo->send_index, 2*halo->n_c_domains + 1, cs_lnum_t,
               halo_alloc_mode);
  CS_MALLOC(halo->index, 2*halo->n_c_domains + 1, cs_lnum_t);

  for (int i = 0; i < 2*halo->n_c_domains + 1; i++) {
    halo->send_index[i] = 0;
    halo->index[i] = 0;
  }

  halo->send_perio_lst = nullptr;
  halo->perio_lst = nullptr;

  if (halo->n_transforms > 0) {

    cs_lnum_t  perio_lst_size = 2*halo->n_transforms * 2*halo->n_c_domains;

    CS_MALLOC(halo->send_perio_lst, perio_lst_size, cs_lnum_t);
    CS_MALLOC(halo->perio_lst, perio_lst_size, cs_lnum_t);

    for (int i = 0; i < perio_lst_size; i++) {
      halo->send_perio_lst[i] = 0;
      halo->perio_lst[i] = 0;
    }

  }

  halo->send_list = nullptr;

  halo->std_send_block_size = 256;  /* Fixed size for now */
  halo->n_std_send_blocks   = 0;
  halo->std_send_blocks = nullptr;

#if defined(HAVE_MPI)
  halo->c_domain_group = MPI_GROUP_NULL;
  halo->c_domain_s_shift = nullptr;
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
 * \param[in]  rn              associated rank neighbors info
 * \param[in]  n_local_elts    number of elements for local rank
 * \param[in]  n_distant_elts  number of distant elements for local rank
 * \param[in]  elt_rank_idx    distant element rank index in rank neighbors,
 *                             ordered by rank (size: n_distant_elts)
 * \param[in]  elt_id          distant element id (at distant rank),
 *                             ordered by rank (size: n_distant_elts)
 * \param[in]  elt_tr_id       distant element transform id (-1 for
 *                             non-periodic elements), null if non-periodic
 * \param[in]  periodicity     optional periodicity, or null
 *
 * \return  pointer to created cs_halo_t structure
 */
/*----------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create_from_rank_neighbors
(
  const cs_rank_neighbors_t  *rn,
  cs_lnum_t                   n_local_elts,
  cs_lnum_t                   n_distant_elts,
  const int                   elt_rank_idx[],
  const cs_lnum_t             elt_id[],
  const int16_t               elt_tr_id[],
  const fvm_periodicity_t    *periodicity
)
{
  cs_halo_t  *halo = nullptr;

  CS_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = 0;
  halo->n_transforms = 0;

  halo->n_rotations = 0;

  halo->periodicity = periodicity;
  halo->send_perio_lst = nullptr;
  halo->perio_lst = nullptr;

  if (periodicity != nullptr) {
    halo->n_transforms = fvm_periodicity_get_n_transforms(periodicity);
  }

#if defined(HAVE_MPI)
  halo->c_domain_group = MPI_GROUP_NULL;
  halo->c_domain_s_shift = nullptr;
#endif

  halo->n_local_elts = n_local_elts;

  for (int i = 0; i < CS_HALO_N_TYPES; i++) {
    halo->n_send_elts[i] = 0;
    halo->n_elts [i] = n_distant_elts;
  }

  halo->std_send_block_size = 256;  /* Fixed size for now */
  halo->n_std_send_blocks   = 0;
  halo->std_send_blocks = nullptr;

  /* Count elements for each rank;
     check they are are ordered lexicographically */

  int tr_mult = (elt_tr_id != nullptr) ? halo->n_transforms + 1 : 1;

  cs_lnum_t rank_tr_count_size = rn->size * tr_mult;
  cs_lnum_t *rank_tr_count;
  CS_MALLOC(rank_tr_count, rank_tr_count_size*2, cs_lnum_t);
  for (int i = 0; i < rank_tr_count_size; i++)
    rank_tr_count[i] = 0;

  bool ordered = true;

  /* Check that the input data is correctly ordered. */

  if (elt_tr_id == nullptr) {
    int rank_idx_prev = -1;
    cs_lnum_t elt_prev = -1;
    for (cs_lnum_t i = 0; i < n_distant_elts; i++) {
      int rank_idx = elt_rank_idx[i];
      if (   rank_idx < rank_idx_prev
          || (rank_idx == rank_idx_prev && elt_id[i] <= elt_prev))
        ordered = false;
      rank_tr_count[rank_idx] += 1;
      rank_idx_prev = rank_idx;
      elt_prev = elt_id[i];
    }
  }
  else {
    int rank_tr_idx_prev = -1;
    cs_lnum_t elt_prev = -1;
    for (cs_lnum_t i = 0; i < n_distant_elts; i++) {
      int rank_idx = elt_rank_idx[i];
      int rank_tr_idx = rank_idx*tr_mult + elt_tr_id[i] + 1;
      assert(rank_tr_idx > -1 && rank_tr_idx < rank_tr_count_size);
      if (   rank_tr_idx < rank_tr_idx_prev
          || (rank_tr_idx == rank_tr_idx_prev && elt_id[i] <= elt_prev))
        ordered = false;
      rank_tr_count[rank_tr_idx] += 1;
      rank_tr_idx_prev = rank_tr_idx;
      elt_prev = elt_id[i];
    }
  }

  if (ordered == false)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s:\n"
       "  Rank and distant element ids passed to this function must\n"
       "  be lexicographically ordered; this is not the case here.",
       __func__);

  /* Now exchange counts with neighboring elements */

  MPI_Comm comm = cs_glob_mpi_comm;
  MPI_Request *request = nullptr;
  MPI_Status *status = nullptr;

  CS_MALLOC(request, rn->size*2, MPI_Request);
  CS_MALLOC(status, rn->size*2, MPI_Status);

  /* Exchange local range with neighbor ranks */

  int request_count = 0;
  const int local_rank = cs::max(cs_glob_rank_id, 0);

  for (int i = 0; i < rn->size; i++) {
    MPI_Irecv(rank_tr_count + rank_tr_count_size + tr_mult*i,
              tr_mult,
              CS_MPI_LNUM,
              rn->rank[i],
              local_rank,
              comm,
              &(request[request_count++]));
  }

  for (int i = 0; i < rn->size; i++) {
    MPI_Isend(rank_tr_count + tr_mult*i,
              tr_mult,
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
    const cs_lnum_t *tr_count = rank_tr_count + tr_mult*i;
    cs_lnum_t r_count = 0, s_count = 0;
    for (int j = 0; j < tr_mult; j++) {
      r_count += tr_count[j];
      s_count += tr_count[rank_tr_count_size + j];
    }
    if (r_count + s_count > 0) {
      halo->n_c_domains += 1;
      if (rn->rank[i] == local_rank) {
        loc_r_index = i;
        loc_r_displ = r_displ;
        assert(r_count == s_count);
      }
      r_displ += r_count;
      recv_count += s_count;
    }
  }

  CS_MALLOC(halo->c_domain_rank, halo->n_c_domains, int);

  CS_MALLOC_HD(halo->send_list, recv_count, cs_lnum_t,
               _halo_buffer_alloc_mode);
  CS_MALLOC_HD(halo->send_index, 2*halo->n_c_domains + 1, cs_lnum_t,
               _halo_buffer_alloc_mode);
  CS_MALLOC(halo->index, halo->n_c_domains*2+1, cs_lnum_t);

  send_count = 0;
  recv_count = 0;

  halo->index[0] = 0;
  halo->send_index[0] = 0;

  if (tr_mult > 1) {
    const int n_transforms = halo->n_transforms;
    const cs_lnum_t  perio_lst_size = 2*n_transforms * 2*halo->n_c_domains;

    CS_MALLOC(halo->perio_lst, perio_lst_size, cs_lnum_t);
    CS_MALLOC(halo->send_perio_lst, perio_lst_size, cs_lnum_t);

    cs_lnum_t *perio_lst = halo->perio_lst;
    cs_lnum_t *send_perio_lst = halo->send_perio_lst;
    for (cs_lnum_t i = 0; i < perio_lst_size; i++) {
      perio_lst[i] = 0;
      send_perio_lst[i] = 0;
    }
  }

  /* Build indexes. */

  int d_idx = 0;

  for (int rn_idx = 0; rn_idx < rn->size; rn_idx++) {
    int i = rn_idx;
    if (loc_r_index > -1) { // local halo first
      if (rn_idx == 0)
        i = loc_r_index;
      else if (rn_idx <= loc_r_index)
        i = rn_idx - 1;
    }

    const cs_lnum_t *r_tr_count = rank_tr_count + tr_mult*i;
    const cs_lnum_t *s_tr_count = r_tr_count + rank_tr_count_size;
    cs_lnum_t r_count = 0, s_count = 0;
    for (int j = 0; j < tr_mult; j++) {
      r_count += r_tr_count[j];
      s_count += s_tr_count[j];
    }
    if (r_count + s_count == 0)
      continue;

    halo->c_domain_rank[d_idx] = rn->rank[i];
    recv_count += r_count;
    send_count += s_count;
    for (int j = 1; j < 3; j++) {
      halo->index[d_idx*2 + j] = recv_count;
      halo->send_index[d_idx*2 + j] = send_count;
    }

    if (tr_mult > 1) {
      cs_lnum_t *perio_lst = halo->perio_lst;
      cs_lnum_t *send_perio_lst = halo->send_perio_lst;
      cs_lnum_t r_shift = halo->index[d_idx*2] + r_tr_count[0];
      cs_lnum_t s_shift = halo->send_index[d_idx*2] + s_tr_count[0];
      for (int j = 1; j < tr_mult; j++) {
        int tr_id = j-1;
        int tr_shift = 4 * (halo->n_c_domains*tr_id + d_idx);
        perio_lst[tr_shift] = r_shift;
        perio_lst[tr_shift+1] = r_tr_count[j];
        r_shift += r_tr_count[j];
        perio_lst[tr_shift+2] = r_shift;
        send_perio_lst[tr_shift] = s_shift;
        send_perio_lst[tr_shift+1] = s_tr_count[j];
        s_shift += s_tr_count[j];
        send_perio_lst[tr_shift+2] = s_shift;
      }
    }

    d_idx += 1;
  }

  assert(d_idx == halo->n_c_domains);
  CS_FREE(rank_tr_count);

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

  /* Local (self) halo first if present */

  if (loc_r_index > -1) {
    cs_lnum_t s_shift = halo->index[0];
    cs_lnum_t s_size  = halo->index[1] - s_shift;
    memcpy(halo->send_list, elt_id + loc_r_displ, s_size * sizeof(cs_lnum_t));
  }

  MPI_Waitall(request_count, request, status);

  CS_FREE(request);
  CS_FREE(status);

  _n_halos += 1;

  cs_halo_create_complete(halo);

  return halo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a halo structure from distant element distant ranks and ids.
 *
 * \param[in]  rn              associated rank neighbors info
 * \param[in]  n_local_elts    number of elements for local rank
 * \param[in]  n_distant_elts  number of distant elements for local rank
 * \param[in]  elt_rank_idx    distant element rank index in rank neighbors,
 *                             ordered by rank (size: n_distant_elts)
 * \param[in]  elt_id          distant element id (at distant rank),
 *                             ordered by rank (size: n_distant_elts)
 *
 * \return  pointer to created cs_halo_t structure
 */
/*----------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create_from_rank_neighbors
(
  const cs_rank_neighbors_t  *rn,
  cs_lnum_t                   n_local_elts,
  cs_lnum_t                   n_distant_elts,
  const int                   elt_rank_idx[],
  const cs_lnum_t             elt_id[]
)
{
  return cs_halo_create_from_rank_neighbors(rn,
                                            n_local_elts,
                                            n_distant_elts,
                                            elt_rank_idx,
                                            elt_id,
                                            nullptr,
                                            nullptr);
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
  if (halo == nullptr)
    return;

  if (*halo == nullptr)
    return;

  cs_halo_t  *_halo = *halo;

#if defined(HAVE_MPI)
  if (_halo->c_domain_group != MPI_GROUP_NULL)
    MPI_Group_free(&(_halo->c_domain_group));

  CS_FREE(_halo->c_domain_s_shift);
#endif

  CS_FREE(_halo->c_domain_rank);

  CS_FREE_HD(_halo->send_list);
  CS_FREE_HD(_halo->send_index);
  CS_FREE(_halo->index);

  CS_FREE(_halo->send_perio_lst);
  CS_FREE(_halo->perio_lst);

  CS_FREE(_halo->std_send_blocks);

  CS_FREE(*halo);

  _n_halos -= 1;

  /* Delete default state if no halo remains */

  if (_n_halos == 0)
    cs_halo_state_destroy(&_halo_state);
}

BEGIN_C_DECLS

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
  CS_MALLOC(hs, 1, cs_halo_state_t);

  cs_halo_state_t hs_ini = {
    .sync_mode = CS_HALO_STANDARD,
    .data_type = CS_DATATYPE_NULL,
    .stride = 0,
    .var_location = CS_ALLOC_HOST,
    .send_buffer_location = CS_ALLOC_HOST,
    .send_buffer_cur = nullptr,
    .n_requests = 0,
    .local_rank_id = -1,
    .send_buffer_size = 0,
    .recv_buffer_size = 0,
    .send_buffer = nullptr,
    .recv_buffer = nullptr
#if defined(HAVE_MPI)
    ,
    .request_size = 0,
    .request = nullptr,
    .status = nullptr,
    .win = MPI_WIN_NULL

#endif

#if defined(HAVE_CUDA)
    ,
    .stream = cs_cuda_get_stream(0)

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
  if (halo_state != nullptr) {
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
    CS_FREE(hs->request);
    CS_FREE(hs->status);
#endif

    CS_FREE(*halo_state);
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
  if (halo != nullptr) {

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
  if (halo == nullptr)
    return;

  /* Reverse update from distant cells */

  cs_lnum_t *send_buf, *recv_buf;

  CS_MALLOC(send_buf, halo->n_send_elts[1], cs_lnum_t);
  CS_MALLOC(recv_buf, halo->n_elts[1], cs_lnum_t);

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

    CS_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    CS_MALLOC(status, halo->n_c_domains*2, MPI_Status);

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

    CS_FREE(request);
    CS_FREE(status);

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

  CS_FREE(recv_buf);

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

  CS_FREE(send_buf);
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
 * \param[out]      send_buf    pointer to send buffer, null for global buffer
 * \param[in, out]  hs          pointer to halo state, null for global state
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

  if (halo == nullptr)
    return _send_buffer;

  cs_halo_state_t  *_hs = (hs != nullptr) ? hs : _halo_state;

  if (_send_buffer == nullptr) {
    size_t send_buffer_size = cs_halo_pack_size(halo, data_type, stride);

    if (send_buffer_size > _hs->send_buffer_size) {
      cs_alloc_mode_t alloc_mode = _halo_buffer_alloc_mode;

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
      _hs->send_buffer_location = alloc_mode;

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
    else
      _hs->send_buffer_location = cs_check_device_ptr(_hs->send_buffer);

    _send_buffer = _hs->send_buffer;
  }
  else {
    _hs->send_buffer_location = cs_check_device_ptr(_send_buffer);
  }

  _hs->var_location = CS_ALLOC_HOST;
  _hs->send_buffer_cur = _send_buffer;

  _hs->sync_mode = sync_mode;
  _hs->data_type = data_type;
  _hs->stride = stride;


#if defined(HAVE_CUDA)
  _hs->stream = cs_cuda_get_stream(0);
#endif

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
 * \param[out]      send_buf    pointer to send buffer, null for global buffer
 * \param[in, out]  hs          pointer to halo state, null for global state
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
  if (halo == nullptr)
    return;

  void *_send_buffer = cs_halo_sync_pack_init_state(halo,
                                                    sync_mode,
                                                    data_type,
                                                    stride,
                                                    send_buf,
                                                    hs);

  const size_t block_size = halo->std_send_block_size;
  const cs_lnum_t *send_list = halo->send_list;
  const cs_lnum_t *send_blocks = nullptr;

  size_t n_send = 0, n_blocks = 0;

  /* If we do not need to send the full halo, use blocks
     to possibly allow for threading */

  if (sync_mode == CS_HALO_STANDARD && halo->n_std_send_blocks > 0) {
    n_send = halo->std_send_block_size * halo->n_std_send_blocks;
    n_blocks = halo->n_std_send_blocks;
    send_blocks = halo->std_send_blocks;
  }
  else {
    n_send = halo->n_send_elts[1];
    n_blocks = (n_send % block_size) ? n_send/block_size + 1 : n_send/block_size;
  }

  #pragma omp parallel for  if (n_send > CS_THR_MIN)
  for (size_t b_id = 0; b_id < n_blocks; b_id++) {

    size_t s_id, e_id;
    if (send_blocks != nullptr) {
      s_id = halo->std_send_blocks[b_id*2];
      e_id = halo->std_send_blocks[b_id*2 + 1];
    }
    else {
      s_id = b_id*block_size;
      e_id = (b_id+1)*block_size;
      if (e_id > n_send)
        e_id = n_send;
    }

    if (data_type == CS_REAL_TYPE) {

      cs_real_t *buffer = (cs_real_t *)_send_buffer;
      cs_real_t *var = (cs_real_t *)val;

      if (stride == 1) {
        for (size_t i = s_id; i < e_id; i++)
          buffer[i] = var[send_list[i]];
      }
      if (stride == 3) { /* Unroll loop for this case */
        for (size_t i = s_id; i < e_id; i++) {
          const size_t j = send_list[i];
          buffer[i*3]     = var[j*3];
          buffer[i*3 + 1] = var[j*3 + 1];
          buffer[i*3 + 2] = var[j*3 + 2];
        }
      }
      else {
        size_t _stride = stride;
        for (size_t i = s_id; i < e_id; i++) {
          size_t j_s = send_list[i]*_stride;
          for (size_t k = 0; k < _stride; k++)
            buffer[i*_stride + k] = var[j_s + k];
        }
      }

    }

    else {

      unsigned char *buffer = (unsigned char *)_send_buffer;
      unsigned char *var = (unsigned char *)val;

      size_t elt_size = cs_datatype_size[data_type] * stride;
      for (size_t i = s_id; i < e_id; i++) {
        size_t i_s = i*elt_size;
        size_t j_s = send_list[i]*elt_size;
        for (size_t k = 0; k < elt_size; k++)
          buffer[i_s + k] = var[j_s + k];
      }

    }

  } /* End of loop on blocks */
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
 *                              null for global buffer
 * \param[in, out]  hs          pointer to halo state, null for global state
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
  if (halo == nullptr)
    return;

  cs_halo_state_t  *_hs = (hs != nullptr) ? hs : _halo_state;

#if defined(HAVE_CUDA)

  void *_send_buf = cs_halo_sync_pack_init_state(halo,
                                                 sync_mode,
                                                 data_type,
                                                 stride,
                                                 send_buf,
                                                 _hs);

  void *val_host_ptr = cs_cuda_get_host_ptr(val);
  void *_send_buf_d = (send_buf != nullptr) ?
    send_buf : cs_get_device_ptr(_send_buf);

  cs_halo_cuda_pack_send_buffer(halo,
                                _hs->stream,
                                sync_mode,
                                data_type,
                                stride,
                                val,
                                _send_buf_d);

  /* We do not try to optimize for CS_ALLOC_HOST_DEVICE or
     CS_ALLOC_HOST_DEVICE_PINNED here as this would require tracking
     the host pointer in the halo state, would possibly lead to
     lower performance when that memory is not pinned, and is not the
     expected dominant future paradigm (using managed memory with
     CS_ALLOC_HOST_DEVICE_SHARED being that one).
     So when synchronizing halos for device arrays using separate host
     and device memory (and when CUDA-aware MPI is not available),
     we will use a local (pinned) memory buffer to receive halo
     data, then copy it back to the device memory. */

  _hs->var_location = CS_ALLOC_HOST_DEVICE_SHARED;

  if (val_host_ptr != val && val != nullptr)
    _hs->var_location = CS_ALLOC_DEVICE;

#else // defined(HAVE_CUDA)

  cs_halo_sync_pack(halo,
                    sync_mode,
                    data_type,
                    stride,
                    val,
                    send_buf,
                    _hs);

  /* As device pointer is passed, cs_check_device_ptr will provide
     the correct result only in case of CS_ALLOC_HOST_DEVICE_SHARED
     (where the pointers are identical). */

  _hs->var_location = cs_check_device_ptr(val);
  if (_hs->var_location != CS_ALLOC_HOST_DEVICE_SHARED)
    _hs->var_location = CS_ALLOC_DEVICE;

#endif
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
 * \param[in, out]  hs          pointer to halo state, null for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_start(const cs_halo_t  *halo,
                   void             *val,
                   cs_halo_state_t  *hs)
{
  if (halo == nullptr)
    return;

  cs_halo_state_t  *_hs = (hs != nullptr) ? hs : _halo_state;

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

  unsigned char *restrict _val = (unsigned char *)val;
  unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

  unsigned char *buffer = (unsigned char *)(_hs->send_buffer_cur);

#if defined(HAVE_ACCEL)

  if (_hs->var_location > CS_ALLOC_HOST) {

#if defined(HAVE_NCCL)

    if (cs_glob_nccl_comm != nullptr) {
      _sync_start_nccl(halo, val, _hs);
      return;
    }

#endif

    /* For CUDA-aware MPI, directly work with buffer on device */

    if (cs_mpi_device_support)
      buffer = (unsigned char *)cs_get_device_ptr(buffer);

    /* For host-based MPI, copy or prefetch buffer */

    else if (cs_glob_n_ranks > 1) {
      size_t pack_size = halo->n_send_elts[CS_HALO_EXTENDED] * elt_size;
      size_t recv_size = halo->n_elts[CS_HALO_EXTENDED] * elt_size;

      /* When array passed is defined on device but is not shared, use separate
         (smaller) CPU buffer for receive (as we cannot know whether a matching
         host array is present without complexifying the API);
         this will be copied back to device at the next step */
      if (_hs->send_buffer_location != CS_ALLOC_HOST_DEVICE_SHARED) {
        void *d_buffer = cs_get_device_ptr(buffer);
        cs_copy_d2h(buffer, d_buffer, pack_size);
      }
      else {
        cs_prefetch_d2h(buffer, pack_size);
      }

      if (_hs->var_location != CS_ALLOC_HOST_DEVICE_SHARED) {
        if (_hs->recv_buffer_size < recv_size) {
          _hs->recv_buffer_size = recv_size;
          CS_FREE_HD(_hs->recv_buffer);
          CS_MALLOC_HD(_hs->recv_buffer, _hs->recv_buffer_size, unsigned char,
                       CS_ALLOC_HOST_DEVICE_PINNED);
        }
        _val_dest = (unsigned char *)_hs->recv_buffer;
      }
    }
  }

#endif /* defined(HAVE_ACCEL) */

#if defined(HAVE_MPI)

  _update_requests(halo, _hs);

  MPI_Datatype mpi_datatype = cs_datatype_to_mpi[_hs->data_type];

  int request_count = 0;
  const int local_rank = cs::max(cs_glob_rank_id, 0);

  /* Receive data from distant ranks */

  for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    cs_lnum_t length = (  halo->index[2*rank_id + end_shift]
                        - halo->index[2*rank_id]);

    if (halo->c_domain_rank[rank_id] != local_rank) {

      if (length > 0) {
        size_t start = (size_t)(halo->index[2*rank_id]);
        unsigned char *dest = _val_dest + start*elt_size;

        MPI_Irecv(dest,
                  length*stride,
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

#else /* defined(HAVE_MPI) */

  const int local_rank = 0;

  /* Receive data from distant ranks */

  for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
    if (halo->c_domain_rank[rank_id] == local_rank)
      _hs->local_rank_id = rank_id;
  }

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
 * \param[in, out]  hs          pointer to halo state, null for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_wait(const cs_halo_t  *halo,
                  void             *val,
                  cs_halo_state_t  *hs)
{
  if (halo == nullptr)
    return;

  cs_halo_state_t  *_hs = (hs != nullptr) ? hs : _halo_state;

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

  if (_hs->var_location > CS_ALLOC_HOST) {

    int device_support = cs_mpi_device_support;

#if defined(HAVE_NCCL)

    if (cs_glob_nccl_comm != nullptr) {
      CS_CUDA_CHECK(cudaStreamSynchronize(_hs->stream));
      CS_CUDA_CHECK(cudaGetLastError());

      device_support = 1;
    }

#endif

    if (device_support == 0) {

      size_t n_loc_elts = halo->n_local_elts;
      size_t n_elts = (   _hs->sync_mode
                       == CS_HALO_EXTENDED) ? halo->n_elts[1] : halo->n_elts[0];
      size_t elt_size = cs_datatype_size[_hs->data_type] * _hs->stride;
      size_t n_bytes = n_elts*elt_size;

      if (n_elts > 0) {
        unsigned char *restrict _val = (unsigned char *)val;
        unsigned char *restrict _val_dest = _val + n_loc_elts*elt_size;

        if (_hs->var_location == CS_ALLOC_HOST_DEVICE_SHARED)
          cs_prefetch_h2d(_val_dest, n_bytes);
        else
          cs_copy_h2d(_val_dest, _hs->recv_buffer, n_bytes);
      }

    }

  }

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
  _hs->send_buffer_cur = nullptr;
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
  if (halo == nullptr)
    return;

  cs_halo_sync_pack(halo,
                    sync_mode,
                    data_type,
                    stride,
                    val,
                    nullptr,
                    nullptr);

  cs_halo_sync_start(halo, val, nullptr);

  cs_halo_sync_wait(halo, val, nullptr);
}

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update array of values on device in case of parallelism
 *        or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \param[in]   halo        pointer to halo structure
 * \param[in]   sync_mode   synchronization mode (standard or extended)
 * \param[in]   data_type   data type
 * \param[in]   stride      number of (interlaced) values by entity
 * \param[in]   val         pointer to variable value array (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_d(const cs_halo_t  *halo,
               cs_halo_type_t    sync_mode,
               cs_datatype_t     data_type,
               int               stride,
               void             *val)
{
  if (halo == nullptr)
    return;

  cs_halo_state_t  *hs = _halo_state;

  cs_halo_sync_pack_d(halo,
                      sync_mode,
                      data_type,
                      stride,
                      val,
                      nullptr,
                      hs);

  cs_halo_sync_start(halo, val, hs);
  cs_halo_sync_wait(halo, val, hs);
}

#endif /* defined(HAVE_ACCEL) */

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
  if (halo == nullptr) {
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
    cs_lnum_t  *index = nullptr, *list = nullptr, *perio_lst = nullptr;

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
      list = nullptr;
      perio_lst = halo->perio_lst;

    }

    bft_printf("    ---------\n\n");
    bft_printf("  n_ghost_cells:        %ld\n"
               "  n_std_ghost_cells:    %ld\n", (long)n_elts[1], (long)n_elts[0]);

    if (index == nullptr)
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

        if (print_level > 0 && list != nullptr) {
          bft_printf("\n            idx     elt id\n");
          for (cs_lnum_t j = index[2*i]; j < index[2*i+1]; j++)
            bft_printf("    %10ld %10ld\n", (long)j, (long)list[j]);
        }

      } /* there are elements on standard neighborhood */

      if (index[2*i+2] - index[2*i+1] > 0) {

        bft_printf("\n  Extended halo\n");
        bft_printf("  idx start %ld:          idx end   %ld:\n",
                   (long)index[2*i+1], (long)index[2*i+2]);

        if (print_level > 0 && list != nullptr) {
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

/*============================================================================
 * Public C++ function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update array of values in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync(const cs_halo_t        *halo,
             cs_halo_type_t          sync_mode,
             [[maybe_unused]] bool   on_device,
             T                       val[])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  if (halo == nullptr)
    return;

  cs_datatype_t datatype = cs_datatype_from_type<T>();

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_pack_d(halo,
                        sync_mode,
                        datatype,
                        1,
                        val,
                        nullptr,
                        nullptr);
  else
#endif
    cs_halo_sync_pack(halo,
                      sync_mode,
                      datatype,
                      1,
                      val,
                      nullptr,
                      nullptr);

  cs_halo_sync_start(halo, val, nullptr);

  cs_halo_sync_wait(halo, val, nullptr);
}

// Force instanciation

template void
cs_halo_sync(const cs_halo_t  *halo,
             cs_halo_type_t    sync_mode,
             bool              on_device,
             cs_real_t         val[]);

template void
cs_halo_sync(const cs_halo_t  *halo,
             cs_halo_type_t    sync_mode,
             bool              on_device,
             cs_lnum_t         val[]);

/*----------------------------------------------------------------------------*/

template <int Stride, typename T>
void
cs_halo_sync(const cs_halo_t       *halo,
             cs_halo_type_t         sync_mode,
             [[maybe_unused]] bool  on_device,
             T                      val[][Stride])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  if (halo == nullptr)
    return;

  cs_datatype_t datatype = cs_datatype_from_type<T>();

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_pack_d(halo,
                        sync_mode,
                        datatype,
                        Stride,
                        val,
                        nullptr,
                        nullptr);
  else
#endif
    cs_halo_sync_pack(halo,
                      sync_mode,
                      datatype,
                      Stride,
                      val,
                      nullptr,
                      nullptr);

  cs_halo_sync_start(halo, val, nullptr);

  cs_halo_sync_wait(halo, val, nullptr);
}

// Force instanciation

template void
cs_halo_sync(const cs_halo_t  *halo,
             cs_halo_type_t    sync_mode,
             bool              on_device,
             cs_real_t         val[][3]);

template void
cs_halo_sync(const cs_halo_t  *halo,
             cs_halo_type_t    sync_mode,
             bool              on_device,
             cs_real_t         val[][6]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update array of values in case of parallelism or periodicity,
 *        using the standard neighborhood.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      Stride      number of (interlaced) values by entity
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync(const cs_halo_t  *halo,
             bool              on_device,
             T                 val[])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  cs_halo_sync(halo, CS_HALO_STANDARD, on_device, val);
}

// Force instanciation

template void
cs_halo_sync(const cs_halo_t  *halo,
             bool              on_device,
             cs_real_t         val[]);

/*----------------------------------------------------------------------------*/

template <int Stride, typename T>
void
cs_halo_sync(const cs_halo_t  *halo,
             bool              on_device,
             T                 val[][Stride])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  cs_halo_sync(halo, CS_HALO_STANDARD, on_device, val);
}

// Force instanciation

template void
cs_halo_sync(const cs_halo_t  *halo,
             bool              on_device,
             cs_real_t         val[][3]);

template void
cs_halo_sync(const cs_halo_t  *halo,
             bool              on_device,
             cs_real_t         val[][6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update ghost cell values of a spatial vector field,
 *        including rotational periodicity if present.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync_r(const cs_halo_t       *halo,
               cs_halo_type_t         sync_mode,
               [[maybe_unused]]bool   on_device,
               T                      val[][3])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  if (halo == nullptr)
    return;

  cs_datatype_t datatype = cs_datatype_from_type<T>();

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_pack_d(halo, sync_mode, datatype, 3, val,
                        nullptr, nullptr);
  else
#endif
    cs_halo_sync_pack(halo, sync_mode, datatype, 3, val,
                      nullptr, nullptr);

  cs_halo_sync_start(halo, val, nullptr);

  cs_halo_sync_wait(halo, val, nullptr);

  if (halo->n_rotations == 0)
    return;

  /* Rotation if needed */

  // TODO: implement this on GPU instead of syncing.
#if defined(HAVE_ACCEL)
  if (on_device)
    cs_sync_d2h((void  *)val);
#endif

  assert(datatype == CS_REAL_TYPE);  // TODO: use templated type below

  cs_halo_perio_sync_var_vect(halo, sync_mode, (cs_real_t *)val, 3);

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_sync_h2d((void  *)val);
#endif
}

// Force instanciation

template void
cs_halo_sync_r(const cs_halo_t  *halo,
               cs_halo_type_t    sync_mode,
               bool              on_device,
               cs_real_t         val[][3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update ghost cell values of a spatial vector field,
 *        including rotational periodicity if present,
 *        using the standard neighborhood.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync_r(const cs_halo_t  *halo,
               bool              on_device,
               T                 val[][3])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  cs_halo_sync_r(halo, CS_HALO_STANDARD, on_device, val);
}

// Force instanciation

template void
cs_halo_sync_r(const cs_halo_t  *halo,
               bool              on_device,
               cs_real_t         val[][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update ghost cell values of a symmetric tensor field,
 *        including rotational periodicity if present.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync_r(const cs_halo_t       *halo,
               cs_halo_type_t         sync_mode,
               [[maybe_unused]]bool   on_device,
               T                      val[][6])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  if (halo == nullptr)
    return;

  cs_datatype_t datatype = cs_datatype_from_type<T>();

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_pack_d(halo, sync_mode, datatype, 6, val,
                        nullptr, nullptr);
  else
#endif
    cs_halo_sync_pack(halo, sync_mode, datatype, 6, val,
                      nullptr, nullptr);

  cs_halo_sync_start(halo, val, nullptr);

  cs_halo_sync_wait(halo, val, nullptr);

  if (halo->n_rotations == 0)
    return;

  /* Rotation if needed */

  // TODO: implement this on GPU instead of syncing.
#if defined(HAVE_ACCEL)
  if (on_device)
    cs_sync_d2h((void  *)val);
#endif

  assert(datatype == CS_REAL_TYPE);  // TODO: use templated type below

  cs_halo_perio_sync_var_sym_tens(halo, sync_mode, (cs_real_t *)val);

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_sync_h2d((void  *)val);
#endif
}

// Force instanciation

template void
cs_halo_sync_r(const cs_halo_t  *halo,
               cs_halo_type_t    sync_mode,
               bool              on_device,
               cs_real_t         val[][6]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update ghost cell values of a symmetric tensor field,
 *        including rotational periodicity if present,
 *        using the standard neighborhood.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync_r(const cs_halo_t  *halo,
               bool              on_device,
               T                 val[][6])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  cs_halo_sync_r(halo, CS_HALO_STANDARD, on_device, val);
}

// Force instanciation

template void
cs_halo_sync_r(const cs_halo_t  *halo,
               bool              on_device,
               cs_real_t         val[][6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update ghost cell values of a non-symmetric tensor field,
 *        including rotational periodicity if present.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * \tparam[in]      T           value type
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       on_device   run on accelerated device if possible
 * \param[in, out]  val         pointer to variable value array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_sync_r(const cs_halo_t       *halo,
               cs_halo_type_t         sync_mode,
               [[maybe_unused]]bool   on_device,
               T                      val[][3][3])
{
  CS_PROFILE_FUNC_RANGE();
  CS_PROFILE_MARK_LINE();

  if (halo == nullptr)
    return;

  cs_datatype_t datatype = cs_datatype_from_type<T>();

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_pack_d(halo, sync_mode, datatype, 6, val,
                        nullptr, nullptr);
  else
#endif
    cs_halo_sync_pack(halo, sync_mode, datatype, 6, val,
                      nullptr, nullptr);

  cs_halo_sync_start(halo, val, nullptr);

  cs_halo_sync_wait(halo, val, nullptr);

  if (halo->n_rotations == 0)
    return;

  /* Rotation if needed */

  // TODO: implement this on GPU instead of syncing.
#if defined(HAVE_ACCEL)
  if (on_device)
    cs_sync_d2h((void  *)val);
#endif

  assert(datatype == CS_REAL_TYPE);  // TODO: use templated type below

  cs_halo_perio_sync_var_tens(halo, sync_mode, (cs_real_t *)val);

#if defined(HAVE_ACCEL)
  if (on_device)
    cs_sync_h2d((void  *)val);
#endif
}

// Force instanciation

template void
cs_halo_sync_r(const cs_halo_t  *halo,
               cs_halo_type_t    sync_mode,
               bool              on_device,
               cs_real_t         val[][3][3]);

/*----------------------------------------------------------------------------*/
