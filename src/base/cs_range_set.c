/*============================================================================
 * Operations related to handling of an owning rank for distributed entities.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_interface.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_range_set.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_range_set.c

  \brief Operations related to handling of an owning rank for distributed
         entities.

  Global element id ranges are assigned to each rank, and global ids are
  defined by a parallel scan type operation counting elements on parallel
  interfaces only once. Each element will appear inside one rank's range and
  outside the range of all other ranks.

  Ranges across different ranks are contiguous.

  This allows building distribution information such as that used in many
  external libraries, such as PETSc, HYPRE, and may also simplify many
  internal operations, where it is needed that elements have a unique owner
  rank, and are ghosted on others (such as linear solvers operating on
  elements which may be on parallel boundaries, such as vertices, edges,
  and faces).

  Elements and their periodic matches allways have distinct global ids;
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define global ids and a partitioning of data based on local ranges
 *        for elements which may be shared across ranks.
 *
 * Global id ranges are assigned to each rank of the interface set's associated
 * communicator, and global ids are defined by a parallel scan type operation
 * counting elements on parallel interfaces only once. Each element will
 * appear inside one rank's range and outside the range of all other ranks.
 * Ranges across different ranks are contiguous.
 *
 * This allows building distribution information such as that used in many
 * external libraries, such as PETSc, HYPRE, and may also simplify many
 * internal operations, where it is needed that elements have a unique owner
 * rank, and are ghosted on others (such as linear solvers operating on
 * elements which may be on parallel boundaries, such as vertices, edges,
 * and faces.
 *
 * Periodicity is ignored for this operation (i.e. elements and their periodic
 * matches are considered distinct).
 *
 * \param[in]   ifs          pointer to interface set structure
 * \param[in]   n_elts       number of elements
 * \param[in]   balance      try to balance shared elements across ranks ?
 * \param[in]   g_id_base    global id base index (usually 0, but 1
 *                           could be used to generate an IO numbering)
 * \param[out]  l_range      global id range assigned to local rank:
 *                           [start, past-the-end[
 * \param[out]  g_id         global id assigned to elements
 */
/*----------------------------------------------------------------------------*/

static void
_interface_set_partition_ids(const cs_interface_set_t  *ifs,
                             cs_lnum_t                  n_elts,
                             bool                       balance,
                             cs_gnum_t                  g_id_base,
                             cs_gnum_t                  l_range[2],
                             cs_gnum_t                 *g_id)
{
  /* Use OpenMP pragma in case of first touch */

# pragma omp parallel for
  for (cs_lnum_t i = 0; i < n_elts; i++)
    g_id[i] = 0;

  /* Second stage: mark elements which are not only local,
     with the corresponding min or max rank, +1
     (as g_id is used as a work array first and cannot have
     negative values, mark with rank + 1 instead of rank */

  int l_rank = CS_MAX(cs_glob_rank_id, 0);

  int ifs_size = cs_interface_set_size(ifs);

  for (int i = 0; i < ifs_size; i++) {

    const cs_interface_t *itf = cs_interface_set_get(ifs, i);

    cs_lnum_t   _tr_index[2] = {0, 0};
    const cs_lnum_t *tr_index = cs_interface_get_tr_index(itf);

    if (tr_index == NULL) {
      _tr_index[1] = cs_interface_size(itf);
      tr_index = _tr_index;
    }

    int itf_rank = cs_interface_rank(itf);

    cs_lnum_t start_id = tr_index[0];
    cs_lnum_t end_id = tr_index[1];

    cs_gnum_t max_rank = CS_MAX(cs_glob_rank_id, itf_rank) + 1;

    const cs_lnum_t *elt_ids = cs_interface_get_elt_ids(itf);

    /* In case of balancing algorithm, assign 1st half
       of elements to lowest rank, 2nd half to highest rank */

    if (balance) {
      cs_gnum_t min_rank= CS_MIN(l_rank, itf_rank) + 1;
      cs_lnum_t mid_id = (start_id + end_id) / 2;
      for (cs_lnum_t j = start_id; j < mid_id; j++) {
        cs_lnum_t k = elt_ids[j];
        if (g_id[k] == 0)
          g_id[k] = min_rank;
        else if (min_rank < g_id[k])
          g_id[k] = min_rank;
      }
      start_id = mid_id;
    }

    for (cs_lnum_t j = start_id; j < end_id; j++) {
      cs_lnum_t k = elt_ids[j];
      g_id[k] = CS_MAX(g_id[k], max_rank);
    }
  }

  /* For balancing option, elements belonging to 2 ranks
     should have a final value, but those belonging to 3
     might have inconsistant values between ranks, so
     toke highest rank for those (should cause only a
     slight imbalance) */

  if (balance)
    cs_interface_set_max(ifs,
                         n_elts,
                         1,
                         true,
                         CS_GNUM_TYPE,
                         g_id);

  /* Now count and mark fo global elements */

  l_range[0] = 0;
  l_range[1] = 0;

  cs_gnum_t l_rank_mark = cs_glob_rank_id + 1;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    if (g_id[i] == 0 || g_id[i] == l_rank_mark)
      l_range[1] += 1;
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t n_local = l_range[1];
    MPI_Scan(&n_local, l_range + 1, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
    l_range[0] = l_range[1] - n_local;
  }
#endif

  /* Mark with 1-based global id, 0 for non owned */

  cs_gnum_t g_id_next = l_range[0] + 1;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    if (g_id[i] == 0 || g_id[i] == l_rank_mark) {
      g_id[i] = g_id_next;
      g_id_next++;
    }
    else
      g_id[i] = 0;
  }

  cs_interface_set_max(ifs,
                       n_elts,
                       1,
                       true,
                       CS_GNUM_TYPE,
                       g_id);

  /* Now assign to correct base */

  if (g_id_base != 1) {
    int64_t g_id_shift = (int64_t)g_id_base - 1;
    for (cs_lnum_t i = 0; i < n_elts; i++)
      g_id[i] += g_id_shift;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define global ids and a partitioning of data based on local ranges
 *        for elements which may be shared across ranks.
 *
 * Global id ranges are assigned to each rank of the interface set's associated
 * communicator, and global ids are defined by a parallel scan type operation
 * counting elements on parallel interfaces only once. Each element will
 * appear inside one rank's range and outside the range of all other ranks.
 * Ranges across different ranks are contiguous.
 *
 * This allows building distribution information such as that used in many
 * external libraries, such as PETSc, HYPRE, and may also simplify many
 * internal operations, where it is needed that elements have a unique owner
 * rank, and are ghosted on others (such as linear solvers operating on
 * elements which may be on parallel boundaries, such as vertices, edges,
 * and faces.
 *
 * Periodicity is ignored for this operation (i.e. elements and their periodic
 * matches are considered distinct).
 *
 * \param[in]       ifs        pointer to interface set structure
 * \param[in]       stride     number of (interlaced) values by entity
 * \param[in]       datatype   type of data considered
 * \param[in]       l_range    global id range assigned to local rank:
 *                             [start, past-the-end[
 * \param[in]       g_id       global id assigned to elements
 * \param[in, out]  val        values buffer
 */
/*----------------------------------------------------------------------------*/

static void
_interface_set_zero_out_of_range(const cs_interface_set_t  *ifs,
                                 cs_datatype_t              datatype,
                                 cs_lnum_t                  stride,
                                 const cs_gnum_t            l_range[2],
                                 const cs_gnum_t           *g_id,
                                 void                      *val)
{
  int ifs_size = cs_interface_set_size(ifs);

  switch(datatype) {

  case CS_FLOAT:
    {
      if (stride > 1) {
        float *v = val;
        for (int i = 0; i < ifs_size; i++) {
          const cs_interface_t *itf = cs_interface_set_get(ifs, i);
          cs_lnum_t        n_elts = cs_interface_size(itf);
          const cs_lnum_t *elt_ids = cs_interface_get_elt_ids(itf);
          for (cs_lnum_t j = 0; j < n_elts; j++) {
            cs_lnum_t k = elt_ids[j];
            if (g_id[k] < l_range[0] || g_id[k] >= l_range[1]) {
              for (cs_lnum_t l = 0; l < stride; l++)
                v[k*stride + l] = 0;
            }
          }
        }
      }
      else {
        float *v = val;
        for (int i = 0; i < ifs_size; i++) {
          const cs_interface_t *itf = cs_interface_set_get(ifs, i);
          cs_lnum_t        n_elts = cs_interface_size(itf);
          const cs_lnum_t *elt_ids = cs_interface_get_elt_ids(itf);
          for (cs_lnum_t j = 0; j < n_elts; j++) {
            cs_lnum_t k = elt_ids[j];
            if (g_id[k] < l_range[0] || g_id[k] >= l_range[1])
              v[k] = 0;
          }
        }
      }
    }
    break;

  case CS_DOUBLE:
    {
      if (stride > 1) {
        double *v = val;
        for (int i = 0; i < ifs_size; i++) {
          const cs_interface_t *itf = cs_interface_set_get(ifs, i);
          cs_lnum_t        n_elts = cs_interface_size(itf);
          const cs_lnum_t *elt_ids = cs_interface_get_elt_ids(itf);
          for (cs_lnum_t j = 0; j < n_elts; j++) {
            cs_lnum_t k = elt_ids[j];
            if (g_id[k] < l_range[0] || g_id[k] >= l_range[1]) {
              for (cs_lnum_t l = 0; l < stride; l++)
                v[k*stride + l] = 0;
            }
          }
        }
      }
      else {
        double *v = val;
        for (int i = 0; i < ifs_size; i++) {
          const cs_interface_t *itf = cs_interface_set_get(ifs, i);
          cs_lnum_t        n_elts = cs_interface_size(itf);
          const cs_lnum_t *elt_ids = cs_interface_get_elt_ids(itf);
          for (cs_lnum_t j = 0; j < n_elts; j++) {
            cs_lnum_t k = elt_ids[j];
            if (g_id[k] < l_range[0] || g_id[k] >= l_range[1])
              v[k] = 0;
          }
        }
      }
    }
    break;

  default:
    {
      cs_lnum_t stride_size = cs_datatype_size[datatype]*stride;
      unsigned char *v = val;
      for (int i = 0; i < ifs_size; i++) {
        const cs_interface_t *itf = cs_interface_set_get(ifs, i);
        cs_lnum_t        n_elts = cs_interface_size(itf);
        const cs_lnum_t *elt_ids = cs_interface_get_elt_ids(itf);
        for (cs_lnum_t j = 0; j < n_elts; j++) {
          cs_lnum_t k = elt_ids[j];
          if (g_id[k] < l_range[0] || g_id[k] >= l_range[1])
            memset(v + k*stride_size, 0, stride_size);
        }
      }
    }
    break;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define global ids and a partitioning of data based on local ranges
 *        for elements which may be shared across ranks or have halo elements.
 *
 * This is a utility function, allowing a similar call for cases where matching
 * elements or parallel ranks are identified using an interface set (for
 * elements which may be on rank boundaries, such as vertices or faces),
 * elements with an associated a halo (such as for cells), or neither
 * (in the single-rank case).
 *
 * Global id ranges are assigned to each rank, and global ids are defined
 * by a parallel scan type operation counting elements on parallel
 * interfaces only once. Each element will appear inside one rank's range
 * and outside the range of all other ranks.
 * Ranges across different ranks are contiguous.
 *
 * This allows building distribution information such as that used in many
 * external libraries, such as PETSc, HYPRE, and may also simplify many
 * internal operations, where it is needed that elements have a unique owner
 * rank, and are ghosted on others (such as linear solvers operating on
 * elements which may be on parallel boundaries, such as vertices, edges,
 * and faces).
 *
 * Elements and their periodic matches allways have distinct global ids;
 *
 * \param[in]   ifs          pointer to interface set structure, or NULL
 * \param[in]   halo         pointer to halo structure, or NULL
 * \param[in]   n_elts       number of elements
 * \param[in]   balance      try to balance shared elements across ranks ?
 *                           (for elements shared across an interface set)
 * \param[in]   g_id_base    global id base index (usually 0, but 1
 *                           could be used to generate an IO numbering)
 * \param[out]  l_range      global id range assigned to local rank:
 *                           [start, past-the-end[
 * \param[out]  g_id         global id assigned to elements
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_define(const cs_interface_set_t  *ifs,
                    const cs_halo_t           *halo,
                    cs_lnum_t                  n_elts,
                    bool                       balance,
                    cs_gnum_t                  g_id_base,
                    cs_gnum_t                  l_range[2],
                    cs_gnum_t                 *g_id)
{
  assert(halo == NULL || ifs == NULL);

  if (ifs != NULL)
    _interface_set_partition_ids(ifs,
                                 n_elts,
                                 balance,
                                 g_id_base,
                                 l_range,
                                 g_id);

  else {

    l_range[0] = g_id_base;
    l_range[1] = g_id_base + n_elts;

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {
      cs_gnum_t loc_shift = n_elts;
      MPI_Scan(&loc_shift, l_range+1, 1, CS_MPI_GNUM, MPI_SUM,
               cs_glob_mpi_comm);
      l_range[1] += g_id_base;
      l_range[0] = l_range[1] - loc_shift;
    }
#endif

#   pragma omp parallel for
    for (cs_lnum_t i = 0; i < n_elts; i++)
      g_id[i] = (cs_gnum_t)i + l_range[0];

    if (halo != NULL)
      cs_halo_sync_untyped(halo,
                           CS_HALO_EXTENDED,
                           sizeof(cs_gnum_t),
                           g_id);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a range set (with associated range and global ids) for the
 *        partitioning of data based on local ranges for elements which may
 *        be shared across ranks or have halo elements.
 *
 * Global id ranges are assigned to each rank of the interface set's associated
 * communicator, and global ids are defined by a parallel scan type operation
 * counting elements on parallel interfaces only once. Each element will
 * appear inside one rank's range and outside the range of all other ranks.
 * Ranges across different ranks are contiguous.
 *
 * The range set maintains pointers to the optional interface set and halo
 * structures, but does not copy them, so those structures should have a
 * lifetime at least as long as the returned range set.
 *
 * \param[in]   ifs          pointer to interface set structure, or NULL
 * \param[in]   halo         pointer to halo structure, or NULL
 * \param[in]   n_elts       number of elements
 * \param[in]   balance      try to balance shared elements across ranks?
 *                           (for elements shared across an interface set)
 * \param[in]   g_id_base    global id base index (usually 0, but 1
 *                           could be used to generate an IO numbering)
 *
 * \return  pointer to created range set structure
 */
/*----------------------------------------------------------------------------*/

cs_range_set_t *
cs_range_set_create(const cs_interface_set_t  *ifs,
                    const cs_halo_t           *halo,
                    cs_lnum_t                  n_elts,
                    bool                       balance,
                    cs_gnum_t                  g_id_base)
{
  cs_gnum_t  *g_id;
  cs_gnum_t  l_range[2];

  BFT_MALLOC(g_id, n_elts, cs_gnum_t);

  cs_range_set_define(ifs,
                      halo,
                      n_elts,
                      balance,
                      g_id_base,
                      l_range,
                      g_id);

  cs_range_set_t *rs = cs_range_set_create_from_shared(ifs,
                                                       halo,
                                                       n_elts,
                                                       l_range,
                                                       g_id);

  rs->_g_id = g_id;

  return rs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a range set (with associated range and global ids) from
 *        an existing partition of data based on local ranges for elements
 *        which may be shared across ranks or have halo elements.
 *
 * The optional interface set, halo, and global element id array are only
 * shared by the range set, not copied, so they should have a lifetime at
 * least as long as the returned range set.
 *
 * \param[in]  ifs      pointer to interface set structure, or NULL
 * \param[in]  halo     pointer to halo structure, or NULL
 * \param[in]  n_elts   number of elements
 * \param[in]  l_range  global id range assigned to local rank:
 *                      [start, past-the-end[
 * \param[in]  g_id     global id assigned to elements
 *
 * \return  pointer to created range set structure
 */
/*----------------------------------------------------------------------------*/

cs_range_set_t *
cs_range_set_create_from_shared(const cs_interface_set_t  *ifs,
                                const cs_halo_t           *halo,
                                cs_lnum_t                  n_elts,
                                cs_gnum_t                  l_range[2],
                                cs_gnum_t                 *g_id)
{
  cs_range_set_t *rs;
  BFT_MALLOC(rs, 1, cs_range_set_t);

  rs->n_elts[0] = 0;
  if (l_range[1] > l_range[0])
    rs->n_elts[0] = l_range[1] - l_range[0];

  rs->n_elts[1] = n_elts;

  rs->l_range[0] = l_range[0];
  rs->l_range[1] = l_range[1];

  rs->ifs= ifs;
  rs->halo= halo;

  rs->g_id = g_id;
  rs->_g_id = NULL;

  return rs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a range set structure.
 *
 * \param[in, out]  rs  pointer to pointer to structure to destroy
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_destroy(cs_range_set_t  **rs)
{
  if (rs != NULL) {
    cs_range_set_t  *_rs = *rs;
    if (_rs != NULL) {
      BFT_FREE(_rs->_g_id);
      BFT_FREE(*rs);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set values of a given array to zero for indexes of elements
 *        outside the local range.
 *
 * If an interface set used to define the range set is available, it may be
 * used to accelerate this operation, as only elements on that interface need
 * to be checked.
 *
 * \param[in]       rs        pointer to range set structure, or NULL
 * \param[in]       datatype  type of data considered
 * \param[in]       stride    number of values per entity (interlaced)
 * \param[in, out]  val       pointer to array values
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_zero_out_of_range(const cs_range_set_t  *rs,
                               cs_datatype_t          datatype,
                               cs_lnum_t              stride,
                               void                  *val)
{
  if (rs == NULL)
    return;

  else if (rs->ifs != NULL) {
    _interface_set_zero_out_of_range(rs->ifs,
                                     datatype,
                                     stride,
                                     rs->l_range,
                                     rs->g_id,
                                     val);
    return;
  }

  cs_lnum_t start_id = 0;
  cs_lnum_t n_elts = rs->n_elts[1];
  const cs_gnum_t  l_range[2] = {rs->l_range[0], rs->l_range[1]};
  const cs_gnum_t *g_id = rs->g_id;

  if (rs->halo != NULL)
    start_id = rs->halo->n_local_elts;

  switch (datatype) {

  case CS_CHAR:
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        char *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  case CS_FLOAT:
#   pragma omp parallel for if(n_elts - start_id > CS_THR_MIN)
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        double *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  case CS_DOUBLE:
#   pragma omp parallel for if(n_elts - start_id > CS_THR_MIN)
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        double *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  case CS_INT32:
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        int32_t *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  case CS_INT64:
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        int64_t *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  case CS_UINT32:
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        uint32_t *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  case CS_UINT64:
    for (cs_lnum_t i = start_id; i < n_elts; i++) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1]) {
        uint64_t *v = val;
        for (cs_lnum_t j = 0; j < stride; j++)
          v[i*stride + j] = 0;
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Called %s with unhandled datatype (%d)."),
              __func__, (int)datatype);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize values elements associated with a range set, using
 *        either a halo or an interface set.
 *
 * \param[in]       rs        pointer to range set structure, or NULL
 * \param[in]       datatype  type of data considered
 * \param[in]       stride    number of values per entity (interlaced)
 * \param[in, out]  val       values buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_sync(const cs_range_set_t  *rs,
                  cs_datatype_t          datatype,
                  cs_lnum_t              stride,
                  void                  *val)
{
  if (rs == NULL)
    return;

  else if (rs->ifs != NULL) {
    cs_lnum_t n_elts = rs->n_elts[1];
    _interface_set_zero_out_of_range(rs->ifs,
                                     datatype,
                                     stride,
                                     rs->l_range,
                                     rs->g_id,
                                     val);
    cs_interface_set_sum(rs->ifs, n_elts, stride, true, datatype, val);
  }

  else if (rs->halo != NULL) {
    if (datatype == CS_REAL_TYPE) {
      if (stride == 1)
        cs_halo_sync_var(rs->halo, CS_HALO_STANDARD, val);
      else
        cs_halo_sync_var_strided(rs->halo, CS_HALO_STANDARD, val, stride);
    }
    else {
      size_t d_size = cs_datatype_size[datatype]*stride;
      cs_halo_sync_untyped(rs->halo, CS_HALO_STANDARD, d_size, val);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gather element values associated with a range set to a compact set.
 *
 * \param[in]   rs        pointer to range set structure, or NULL
 * \param[in]   datatype  type of data considered
 * \param[in]   stride    number of values per entity (interlaced)
 * \param[in]   src_val   source values buffer
 * \param[out]  dest_val  destination values buffer (may be identical to
 *                        src_val, in which case operation is "in-place")
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_gather(const cs_range_set_t  *rs,
                    cs_datatype_t          datatype,
                    cs_lnum_t              stride,
                    const void            *src_val,
                    void                  *dest_val)
{
  if (rs == NULL)
    return;

  else if (rs->halo != NULL)
    return;

  cs_lnum_t n_elts = rs->n_elts[1];
  cs_lnum_t d_size = cs_datatype_size[datatype]*stride;
  const cs_gnum_t  l_range[2] = {rs->l_range[0], rs->l_range[1]};
  const cs_gnum_t  *g_id = rs->g_id;

  /* Case with overlapping source and destination */

  if (src_val == dest_val) {

    unsigned char *src = dest_val;
    unsigned char *dest = dest_val;

    cs_lnum_t i = 0;

    /* First set of compact values */

    while (i < n_elts) {
      if (g_id[i] < l_range[0] || g_id[i] >= l_range[1])
        break;
      i++;
    }

    src += i*d_size;
    dest += i*d_size;

    while (i < n_elts) {
      if (g_id[i] >= l_range[0] && g_id[i] < l_range[1]) {
        memcpy(dest, src, d_size);
        dest += d_size;
      }
      src += d_size;
      i++;
    }

  }

  /* Case with non-overlapping values */

  else { /* src_val != dest_val */

    const unsigned char *src = src_val;
    unsigned char *dest = dest_val;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      if (g_id[i] >= l_range[0] && g_id[i] < l_range[1]) {
        memcpy(dest, src, d_size);
        dest += d_size;
      }
      src += d_size;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Scatter element values associated with a range set to the full set.
 *
 * This includes parallel synchronization when the range set is associated
 * with a halo or interface set structure.
 *
 * \param[in]   rs        pointer to range set structure, or NULL
 * \param[in]   datatype  type of data considered
 * \param[in]   stride    number of values per entity (interlaced)
 * \param[in]   src_val   source values buffer
 * \param[out]  dest_val  destination values buffer (may be identical to
 *                        src_val, in which case operation is "in-place")
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_scatter(const cs_range_set_t  *rs,
                     cs_datatype_t          datatype,
                     cs_lnum_t              stride,
                     const void            *src_val,
                     void                  *dest_val)
{
  if (rs == NULL)
    return;

  else if (rs->halo != NULL) {
    cs_range_set_sync(rs, datatype, stride, dest_val);
    return;
  }

  cs_lnum_t d_size = cs_datatype_size[datatype]*stride;
  const cs_gnum_t  l_range[2] = {rs->l_range[0], rs->l_range[1]};
  const cs_gnum_t  *g_id = rs->g_id;

  /* Case with overlapping source and destination */

  if (src_val == dest_val) {

    cs_lnum_t i_d = rs->n_elts[1] - 1;
    cs_lnum_t i_s = rs->n_elts[0] - 1;

    /* First set of compact values */

    unsigned char *src = (unsigned char *)dest_val + i_s*d_size;
    unsigned char *dest = (unsigned char *)dest_val + i_d*d_size;

    while (i_d >= 0) {

      if (g_id[i_d] >= l_range[0] && g_id[i_d] < l_range[1]) {
        memmove(dest, src, d_size);
        src -= d_size;
        i_s--;
      }
      dest -= d_size;
      i_d--;

      if (i_d == i_s)
        break;

    }

  }

  /* Case with non-overlapping values */

  else { /* src_val != dest_val */

    cs_lnum_t n_elts = rs->n_elts[1];

    const unsigned char *src = src_val;
    unsigned char *dest = dest_val;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      if (g_id[i] >= l_range[0] && g_id[i] < l_range[1]) {
        memcpy(dest, src, d_size);
        src += d_size;
      }
      dest += d_size;
    }

  }

  /* Now synchronize values */

  cs_range_set_sync(rs, datatype, stride, dest_val);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
