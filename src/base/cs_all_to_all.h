#ifndef __CS_ALL_TO_ALL_H__
#define __CS_ALL_TO_ALL_H__

/*============================================================================
 * All-to-all parallel data exchange.
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
#include "cs_rank_neighbors.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * All-to-all distributor ordering and metadata flags.
 */

#define CS_ALL_TO_ALL_USE_DEST_ID             (1 << 0)
#define CS_ALL_TO_ALL_ORDER_BY_SRC_RANK       (1 << 1)

#define CS_ALL_TO_ALL_NO_REVERSE              (1 << 2)
#define CS_ALL_TO_ALL_NEED_SRC_RANK           (1 << 3)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* All-to-all algorithm choice */

typedef enum {

  CS_ALL_TO_ALL_MPI_DEFAULT,
  CS_ALL_TO_ALL_HYBRID,
  CS_ALL_TO_ALL_CRYSTAL_ROUTER

} cs_all_to_all_type_t;

/* Opaque all-to-all distribution structure */

#if defined(HAVE_MPI)

typedef struct _cs_all_to_all_t  cs_all_to_all_t;

#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor based on destination rank.
 *
 * This is a collective operation on communicator comm.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_USE_DEST_ID,
 * data exchanged will be ordered by the array passed to the
 * \c dest_id argument. For \c n total values received on a rank
 * (as given by \ref cs_all_to_all_n_elts_dest), those destination ids
 * must be in the range [0, \c n[.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,
 * data exchanged will be ordered by source rank (this is incompatible
 * with \ref CS_ALL_TO_ALL_USE_DEST_ID.
 *
 * \attention
 * The \c dest_rank and \c dest_id arrays are only referenced by
 * the distributor, not copied, and must remain available throughout
 * the distributor's lifetime. They may be fully transferred to
 * the structure if not needed elsewhere using the
 * \ref cs_all_to_all_transfer_dest_rank and
 * \ref cs_all_to_all_transfer_dest_id functions.
 *
 * \param[in]  n_elts       number of elements
 * \param[in]  flags        sum of ordering and metadata flag constants
 * \param[in]  dest_id      element destination id (required if flags
 *                          contain \ref CS_ALL_TO_ALL_USE_DEST_ID),
 *                          or NULL
 * \param[in]  dest_rank    destination rank for each element
 * \param[in]  comm         associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create(size_t            n_elts,
                     int               flags,
                     const cs_lnum_t  *dest_id,
                     const int         dest_rank[],
                     MPI_Comm          comm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor for elements whose destination
 * rank is determined from global numbers and block distribution information.
 *
 * This is a collective operation on communicator comm.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_USE_DEST_ID,
 * data exchanged will be ordered by global element number.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,
 * data exchanged will be ordered by source rank (this is incompatible
 * with \ref CS_ALL_TO_ALL_USE_DEST_ID.
 *
 * \param[in]  n_elts       number of elements
 * \param[in]  flags        sum of ordering and metadata flag constants
 * \param[in]  src_gnum     global source element numbers
 * \param[in]  bi           destination block distribution info
 * \param[in]  comm         associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_from_block(size_t                 n_elts,
                                int                    flags,
                                const cs_gnum_t       *src_gnum,
                                cs_block_dist_info_t   bi,
                                MPI_Comm               comm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy an all-to-all distributor.
 *
 * \param[in, out]  d   pointer to associated all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_destroy(cs_all_to_all_t **d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer ownership of destination rank to an all-to-all distributor.
 *
 * The dest_rank array should be the same as the one used for the creation of
 * the distributor.
 *
 * \param[in, out]  d          pointer to associated all-to-all distributor
 * \param[in, out]  dest_rank  pointer to element destination rank
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_transfer_dest_rank(cs_all_to_all_t   *d,
                                 int              **dest_rank);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer ownership of destination ids to an
 *        all-to-all distributor.
 *
 * The dest_id array should be the same as the one used for the creation of
 * the distributor.
 *
 * \param[in, out]  d        pointer to associated all-to-all distributor
 * \param[in, out]  dest_id  pointer to element destination id
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_transfer_dest_id(cs_all_to_all_t   *d,
                               cs_lnum_t        **dest_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements associated with all-to-all distributor.
 *
 * The number of elements is the number of elements received after exchange.
 *
 * If no exchange has been done yet (depending on the communication protocol),
 * metadata will be exchanged by this call, so it is a collective operation.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 *
 * \return  number of elements associated with distributor.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_all_to_all_n_elts_dest(cs_all_to_all_t  *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Communicate array data using all-to-all distributor.
 *
 * If a destination buffer is provided, it should be of sufficient size for
 * the number of elements returned by \ref cs_all_to_all_n_elts_dest
 * (multiplied by stride and datatype size).
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * If used in reverse mode, data is still communicated from src_data
 * to dest_buffer or an internal buffer, but communication direction
 * (i.e. source and destination ranks) are reversed.
 *
 * This is obviously a collective operation, and all ranks must provide
 * the same datatype, stride, and reverse values.
 *
 * \param[in, out]  d          pointer to associated all-to-all distributor
 * \param[in]       datatype   type of data considered
 * \param[in]       stride     number of values per entity (interlaced),
 * \param[in]       reverse    if true, communicate in reverse direction
 * \param[in]       src_data   source data
 * \param[out]      dest_data  pointer to destination data, or NULL
 *
 * \return pointer to destination data (dest_buffer if non-NULL)
 */
/*----------------------------------------------------------------------------*/

void *
cs_all_to_all_copy_array(cs_all_to_all_t   *d,
                         cs_datatype_t      datatype,
                         int                stride,
                         bool               reverse,
                         const void        *src_data,
                         void              *dest_data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Communicate local index using all-to-all distributor.
 *
 * If a destination buffer is provided, it should be of sufficient size for
 * the number of elements returned by \ref cs_all_to_all_n_elts_dest.
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * If used in reverse mode, data is still communicated from src_index
 * to dest_index or an internal buffer, but communication direction
 * (i.e. source and destination ranks) are reversed.
 *
 * This is obviously a collective operation, and all ranks must provide
 * the same value for the reverse parameter.
 *
 * \param[in, out]  d           pointer to associated all-to-all distributor
 * \param[in]       reverse     if true, communicate in reverse direction
 * \param[in]       src_index   source index
 * \param[out]      dest_index  pointer to destination index, or NULL
 *
 * \return pointer to destination data (dest_buffer if non-NULL)
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_all_to_all_copy_index(cs_all_to_all_t  *d,
                         bool              reverse,
                         const cs_lnum_t  *src_index,
                         cs_lnum_t        *dest_index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Communicate local index using all-to-all distributor.
 *
 * If a destination buffer is provided, it should be of sufficient size for
 * the number of elements indicated by
 * dest_index[\ref cs_all_to_all_n_elts_dest "cs_all_to_all_n_elts_dest(d)"];
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * If used in reverse mode, data is still communicated from src_index
 * to dest_index or an internal buffer, but communication direction
 * (i.e. source and destination ranks) are reversed.
 *
 * This is obviously a collective operation, and all ranks must provide
 * the same value for the reverse parameter.
 *
 * \param[in, out]  d           pointer to associated all-to-all distributor
 * \param[in]       datatype    type of data considered
 * \param[in]       reverse     if true, communicate in reverse direction
 * \param[in]       src_index   source index
 * \param[in]       src_data    source data
 * \param[in]       dest_index  destination index
 * \param[out]      dest_data   pointer to destination data, or NULL
 *
 * \return pointer to destination data (dest_buffer if non-NULL)
 */
/*----------------------------------------------------------------------------*/

void *
cs_all_to_all_copy_indexed(cs_all_to_all_t  *d,
                           cs_datatype_t     datatype,
                           bool              reverse,
                           const cs_lnum_t  *src_index,
                           const void       *src_data,
                           const cs_lnum_t  *dest_index,
                           void             *dest_data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get array of source element ranks associated with an
 *        all-to-all distributor.
 *
 * This function should be called only after \ref cs_all_to_all_copy_array,
 * and allocates and returns an array of source element ranks matching the
 * exchanged data elements.
 *
 * It should also be called only if the distributor creation flags match
 * CS_ALL_TO_ALL_NEED_SRC_RANK or CS_ALL_TO_ALL_ORDER_BY_SRC_RANK.
 *
 * The returned data is owned by the caller, who is responsible for freeing
 * it when no longer needed.
 *
 * If source ranks are not available (depending on the distributor
 * creation function and options), the matching pointer will
 * be set to NULL.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 *
 * \return  array of source ranks (or NULL)
 */
/*----------------------------------------------------------------------------*/

int *
cs_all_to_all_get_src_rank(cs_all_to_all_t  *d);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current type of all-to-all distributor algorithm choice.
 *
 * \return  current type of all-to-all distributor algorithm choice
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_type_t
cs_all_to_all_get_type(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current type of all-to-all distributor algorithm choice.
 *
 * \param  t  type of all-to-all distributor algorithm choice to select
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_set_type(cs_all_to_all_type_t  t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current type of hybrid all-to-all distributor parameters.
 *
 * \param[out]  rne_type  type of metadata exchange algorithm, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_get_hybrid_parameters(cs_rank_neighbors_exchange_t  *rne_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current type of all-to-all distributor algorithm choice.
 *
 * \param[in]  rne_type  type of metadata exchange algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_set_hybrid_parameters(cs_rank_neighbors_exchange_t  rne_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log performance information relative to instrumented all-to-all
 * distribution.
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_log_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ALL_TO_ALL_H__ */
