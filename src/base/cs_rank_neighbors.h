#ifndef __CS_RANK_NEIGHBORS_H__
#define __CS_RANK_NEIGHBORS_H__

/*============================================================================
 * Management of parallel rank neighbors.
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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_block_dist.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Exchange algorithm choice */

typedef enum {

  CS_RANK_NEIGHBORS_PEX,
  CS_RANK_NEIGHBORS_NBX,
  CS_RANK_NEIGHBORS_CRYSTAL_ROUTER

} cs_rank_neighbors_exchange_t;

#if defined(HAVE_MPI)

/* Rank neighbors structure */

typedef struct {

  int   size;   /*!< number of neighboring ranks */
  int  *rank;   /*!< neighboring rank ids (always ordered)  */

} cs_rank_neighbors_t;

#endif

/*=============================================================================
 * Global variables
 *============================================================================*/

extern const char  *cs_rank_neighbors_exchange_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a rank neighbors structure base on a list of element ranks
 *
 * \param[in]  n_elts    number of elements
 * \param[in]  elt_rank  element rank
 *
 * \return  pointer to new rank neighborhood.
 */
/*----------------------------------------------------------------------------*/

cs_rank_neighbors_t *
cs_rank_neighbors_create(size_t     n_elts,
                         const int  elt_rank[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a rank neighborhood structure
 *
 * \param[in, out]  n   pointer to associated rank neighborhood
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_destroy(cs_rank_neighbors_t  **n);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Given a list of element ranks, determine the associated element
 *        rank indexes in a rank neighborhood structure.
 *
 * The elt_rank and elt_rank_index may be identical,
 * in which case it is updated.
 *
 * \param[in]   n               pointer to rank neighborhood structure
 * \param[in]   n_elts          number of elements
 * \param[in]   elt_rank        element rank (size: n_elts)
 * \param[out]  elt_rank_index  element rank index in neighborhood
 *                              (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_to_index(const cs_rank_neighbors_t  *n,
                           size_t                      n_elts,
                           const int                   elt_rank[],
                           int                        *elt_rank_index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Symmetrize a rank neighborhood structure.
 *
 * This is a collective operation ,which ensures that if rank i has j among
 * its neighbors, then j will also jave i among its neighbors.
 *
 * \param[in, out]  n     pointer to rank neighborhood structure
 * \param[in]       comm  associated communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_symmetrize(cs_rank_neighbors_t  *n,
                             MPI_Comm              comm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Given a list of element rank indexes, count occurences for
 *        a rank neighborhood structure.
 *
 * \param[in]   n               pointer to rank neighborhood structure
 * \param[in]   n_elts          number of elements
 * \param[in]   elt_rank_index  element rank index in neighborhood
 *                              (size: n_elts)
 * \param[out]  elt_rank_count  element rank count in neighborhood
 *                              (size: n->size)
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_count(const cs_rank_neighbors_t  *n,
                        size_t                      n_elts,
                        const int                  *elt_rank_index,
                        cs_lnum_t                  *elt_rank_count);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange send and receive counts for rank neighborhoods.
 *
 * This allocates the n_recv ranks neighborhood structure and the
 * recv_count counts array, which the caller is responsible for freeing.
 *
 * \param[in]   n_send       pointer to rank neighborhood used for sending
 * \param[out]  n_recv       pointer to rank neighborhood used for receiving
 * \param[in]   send_count   pointer to rank neighborhood used for sending
 * \param[in]   recv_count   pointer to rank neighborhood used for sending
 * \param[in]   comm         associated communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_sync_count(const cs_rank_neighbors_t   *n_send,
                             cs_rank_neighbors_t        **n_recv,
                             const cs_lnum_t             *send_count,
                             cs_lnum_t                  **recv_count,
                             MPI_Comm                     comm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange send and receive counts for rank neighborhoods,
 *        using a given method.
 *
 * This allocates the n_recv ranks neighborhood structure and the
 * recv_count counts array, which the caller is responsible for freeing.
 *
 * \param[in]   n_send         pointer to rank neighborhood used for sending
 * \param[out]  n_recv         pointer to rank neighborhood used for receiving
 * \param[in]   send_count     pointer to rank neighborhood used for sending
 * \param[in]   recv_count     pointer to rank neighborhood used for sending
 * \param[in]   exchange_type  exchange type
 * \param[in]   comm           associated communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_sync_count_m(const cs_rank_neighbors_t      *n_send,
                               cs_rank_neighbors_t           **n_recv,
                               const cs_lnum_t                *send_count,
                               cs_lnum_t                     **recv_count,
                               cs_rank_neighbors_exchange_t    exchange_type,
                               MPI_Comm                        comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current type of rank neighbors collective algorithm choice.
 *
 * \return  current type of rank neighbors collective algorithm choice
 */
/*----------------------------------------------------------------------------*/

cs_rank_neighbors_exchange_t
cs_rank_neighbors_get_exchange_type(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current type of rank neighbors collective algorithm choice.
 *
 * \param  t  type of rank neighbors collective algorithm choice
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_set_exchange_type(cs_rank_neighbors_exchange_t  t);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RANK_NEIGHBORS_H__ */
