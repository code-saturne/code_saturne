#ifndef __CS_CRYSTAL_ROUTER_H__
#define __CS_CRYSTAL_ROUTER_H__

/*============================================================================
 * Crystal Router parallel exchange algorithm based operations.
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

/*!
 * Crystal router ordering and metadata flags.
 */

#define CS_CRYSTAL_ROUTER_USE_DEST_ID         (1 << 0)

#define CS_CRYSTAL_ROUTER_ADD_SRC_ID          (1 << 1)
#define CS_CRYSTAL_ROUTER_ADD_SRC_RANK        (1 << 2)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Opaque crystal router structure */

#if defined(HAVE_MPI)

typedef struct _cs_crystal_router_t  cs_crystal_router_t;

#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Crystal Router for strided data.
 *
 * If the flags constant contains \ref CS_CRYSTAL_ROUTER_USE_DEST_ID,
 * data exchanged will be ordered by the array passed to the
 * \c dest_id argument. For \c n total values received on a rank
 * (as given by \ref cs_crystal_router_n_elts), those destination ids
 * must be in the range [0, \c n[.
 *
 * If the flags bit mask matches \ref CS_CRYSTAL_ROUTER_ADD_SRC_ID,
 * source ids are added. If it matches \ref CS_CRYSTAL_ROUTER_ADD_SRC_RANK,
 * source rank metadata is added.
 *
 * \param[in]  n_elts            number of elements
 * \param[in]  stride            number of values per entity (interlaced)
 * \param[in]  datatype          type of data considered
 * \param[in]  flags             add destination id ?
 * \param[in]  elt               element values
 * \param[in]  dest_id           element destination id, or NULL
 * \param[in]  dest_rank         destination rank for each element
 * \param[in]  comm              associated MPI communicator
 *
 * \return  pointer to new Crystal Router structure.
 */
/*----------------------------------------------------------------------------*/

cs_crystal_router_t *
cs_crystal_router_create_s(size_t            n_elts,
                           int               stride,
                           cs_datatype_t     datatype,
                           int               flags,
                           const void       *elt,
                           const cs_lnum_t  *dest_id,
                           const int         dest_rank[],
                           MPI_Comm          comm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Crystal Router for indexed data.
 *
 * If the flags constant contains \ref CS_CRYSTAL_ROUTER_USE_DEST_ID,
 * data exchanged will be ordered by the array passed to the
 * \c dest_id argument. For \c n total values received on a rank
 * (as given by \ref cs_crystal_router_n_elts), those destination ids
 * must be in the range [0, \c n[.
 *
 * If the flags bit mask matches \ref CS_CRYSTAL_ROUTER_ADD_SRC_ID,
 * source ids are added. If it matches \ref CS_CRYSTAL_ROUTER_ADD_SRC_RANK,
 * source rank metadata is added.
 *
 * \param[in]  n_elts            number of elements
 * \param[in]  stride            number of values per entity (interlaced)
 * \param[in]  datatype          type of data considered
 * \param[in]  flags             add destination id ?
 * \param[in]  elt_idx           element values start and past-the-last index
 * \param[in]  elt               element values
 * \param[in]  dest_id           element destination id, or NULL
 * \param[in]  dest_rank         destination rank for each element
 * \param[in]  comm              associated MPI communicator
 *
 * \return  pointer to new Crystal Router structure.
 */
/*----------------------------------------------------------------------------*/

cs_crystal_router_t *
cs_crystal_router_create_i(size_t            n_elts,
                           cs_datatype_t     datatype,
                           int               flags,
                           const cs_lnum_t  *eld_idx,
                           const void       *elt,
                           const cs_lnum_t  *dest_id,
                           const int         dest_rank[],
                           MPI_Comm          comm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a Crystal Router.
 *
 * \param[in, out]  cr   pointer to associated Crystal Router
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_destroy(cs_crystal_router_t  **cr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange data with a Crystal Router.
 *
 * Order of data from a same source rank is preserved.
 *
 * \param[in, out]  cr  pointer to associated Crystal Router
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_exchange(cs_crystal_router_t  *cr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements associated with Crystal Router.
 *
 * The number of elements is the number of elements received after exchange.
 *
 * \param[in]  cr  pointer to associated Crystal Router
 *
 * \return  number of elements associated with distributor.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_crystal_router_n_elts(const cs_crystal_router_t  *cr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data elements associated with a Crystal Router.
 *
 * Depending on the creation options, some elements may be available or not.
 *
 * For each element type, a pointer to a pointer to a buffer is provided.
 * If the pointer to the buffer is non-null, the buffer must be of sufficient
 * size for the number of elements returned by \ref cs_crystal_router_n_elts.
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * Note also that if the destination id is provided in the Crystal Router,
 * it will be applied automatically to all extracted arrays, except the
 * dest_id array itself, which is always in receive order. If the Crystal
 * Router does not contain destination id info but the \c dest_id
 * argument points to a non-NULL value, the provided id will be used to
 * order extracted data. This allows saving the destination id on the receive
 * side, and not re-sending it (saving bandwidth) for subsequent calls
 * with a similar Crystal Router.
 *
 * \param[in]       cr          pointer to associated Crystal Router
 * \param[out]      src_rank    pointer to (pointer to) source rank array,
 *                              or NULL
 * \param[in, out]  dest_id     pointer to (pointer to) destination id array,
 *                              or NULL
 * \param[out]      src_id      pointer to (pointer to) source id array, or NULL
 * \param[out]      data_index  pointer to (pointer to) destination index,
 *                              or NULL
 * \param[out]      data        pointer to (pointer to) destination data,
 *                              or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_get_data(cs_crystal_router_t   *cr,
                           int                  **src_rank,
                           cs_lnum_t            **dest_id,
                           cs_lnum_t            **src_id,
                           cs_lnum_t            **data_index,
                           void                 **data);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log performance information relative to Crystal Router exchange.
 *
 * Call count is based on Crystal router creation calls, as the number
 * of exchange and data access calls should be identical.
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_log_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CRYSTAL_ROUTER_H__ */
