#ifndef __CS_ALL_TO_ALL_H__
#define __CS_ALL_TO_ALL_H__

/*============================================================================
 * All-to-all parallel data exchange.
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

/* All-to-all algorithm choice */

typedef enum {

  CS_ALL_TO_ALL_MPI_DEFAULT,
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

/*----------------------------------------------------------------------------
 * Create an all-to-all distributor for strided data.
 *
 * parameters:
 *   n_elts      <-- number of elements
 *   stride      <-- number of values per entity (interlaced)
 *   datatype    <-- type of data considered
 *   elt         <-- element values
 *   dest_rank   <-- destination rank for each element
 *   comm        <-- associated MPI communicator
 *
 * returns:
 *   pointer to new all-to-all distributor
 *---------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_s(size_t          n_elts,
                       int             stride,
                       cs_datatype_t   datatype,
                       void           *elt,
                       const int       dest_rank[],
                       MPI_Comm        comm);

/*----------------------------------------------------------------------------
 * Create an all-to-all distributor for strided data with additional metadata.
 *
 * This variant allows optional tracking of destination ids or global
 * numbers associated with elements, as well as their source ids.
 *
 * In cases where those arrays are required and already available, this
 * may avoid the need for a specific element values buffer mixing actual
 * data values and numbering metadata. It also makes extraction of the
 * metadata easier using cs_all_to_all_get_id_pointers().
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL depending on elt_id values)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   elt              <-- element values
 *   dest_id          <-- element destination id, global number, or NULL
 *   dest_rank        <-- destination rank for each element
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   pointer to new all-to-all distributor
 *---------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_with_ids_s(size_t            n_elts,
                                int               stride,
                                cs_datatype_t     datatype,
                                cs_datatype_t     dest_id_datatype,
                                bool              add_src_id,
                                void             *elt,
                                void             *dest_id,
                                const int         dest_rank[],
                                MPI_Comm          comm);

/*----------------------------------------------------------------------------
 * Create an all-to-all distributor for strided data with additional metadata,
 * with destination rank determined from global numbers and block distribution
 * information.
 *
 * This variant allows optional tracking of destination ids or global
 * numbers associated with elements, as well as their source ids.
 *
 * In cases where those arrays are required and already available, this
 * may avoid the need for a specific element values buffer mixing actual
 * data values and numbering metadata. It also makes extraction of the
 * metadata easier using cs_all_to_all_get_id_pointers().
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   elt              <-- element values
 *   elt_gnum         <-- global element numbers
 *   bi               <-- destination block distribution info
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   pointer to new all-to-all distributor
 *---------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_from_block_s(size_t                 n_elts,
                                  int                    stride,
                                  cs_datatype_t          datatype,
                                  cs_datatype_t          dest_id_datatype,
                                  bool                   add_src_id,
                                  void                  *elt,
                                  const cs_gnum_t       *elt_gnum,
                                  cs_block_dist_info_t   bi,
                                  MPI_Comm               comm);

/*----------------------------------------------------------------------------
 * Destroy an all-to-all distributor.
 *
 * parameters:
 *   d <-> pointer to associated all-to-all distributor
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_destroy(cs_all_to_all_t **d);

/*----------------------------------------------------------------------------
 * Exchange data with an all-to-all distributor.
 *
 * Order of data from a same source rank is preserved.
 *
 * parameters:
 *   d <-> pointer to associated all-to-all distributor
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_exchange(cs_all_to_all_t  *d);

/*----------------------------------------------------------------------------
 * Sort stride crystal router data by source rank.
 *
 * parameters:
 *   d <-> pointer to associated all-to-all distributor
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_sort_by_source_rank(cs_all_to_all_t  *d);

/*----------------------------------------------------------------------------
 * Get number of elements associated with all-to-all distributor.
 *
 * The number of elements is the number of elements received after exchange.
 *
 * parameters:
 *   d <-- pointer to associated all-to-all distributor
 *---------------------------------------------------------------------------*/

cs_lnum_t
cs_all_to_all_n_elts(const cs_all_to_all_t  *d);

/*----------------------------------------------------------------------------
 * Swap source and destination ranks of all-to-all distributor.
 *
 * parameters:
 *   d <-> associated all-to-all distributor pointer
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_swap_src_dest(cs_all_to_all_t  *d);

/*----------------------------------------------------------------------------
 * Get pointer to data elements associated with an all-to-all distributor.
 *
 * This allows modification and/or extraction of those elements.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * parameters:
 *   d           <-- pointer to associated all-to-all distributor
 *   data_stride --> stride (in bytes) between data items
 *   data_ptr    --> pointer to data items
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_get_data_pointer(cs_all_to_all_t   *d,
                               size_t            *data_stride,
                               unsigned char    **data);

/*----------------------------------------------------------------------------
 * Get pointer to ranks of elements associated with an all-to-all distributor.
 *
 * This allows modification and/or extraction of those ranks.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * parameters:
 *   d           <-- pointer to associated all-to-all distributor
 *   rank_stride --> stride (in integers) between rank values
 *   src_rank    --> pointer to source rank values (or NULL)
 *   dest_rank   --> pointer to destination rank values (or NULL)
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_get_rank_pointers(cs_all_to_all_t   *d,
                                size_t            *rank_stride,
                                int              **src_rank,
                                int              **dest_rank);

/*----------------------------------------------------------------------------
 * Get pointer to source or destination rank element ids associated with an
 * all-to-all distributor.
 *
 * If a requested type of id is not available (depending on the all-to-all
 * distributor creation function and options), the matching pointer will
 * be set to NULL.
 *
 * This allows modification and/or extraction of those ids, though it is
 * intended primarily for identification.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * parameters:
 *   d         <-- pointer to associated all-to-all distributor
 *   id_stride --> stride (in integers) between id items
 *   dest_id   --> pointer to destination ids (or NULL)
 *   dest_id   --> pointer to source ids (or NULL)
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_get_id_pointers(cs_all_to_all_t   *d,
                              size_t            *id_stride,
                              cs_lnum_t        **dest_id,
                              cs_lnum_t        **src_id);

/*----------------------------------------------------------------------------
 * Get pointer to element global numbers associated with an all-to-all
 * distributor.
 *
 * If this data is not available (depending on the all-to-all distributor
 * creation function and options), the matching pointer will
 * be set to NULL.
 *
 * This allows modification and/or extraction of those numbers.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * parameters:
 *   d           <-- pointer to associated all-to-all distributor
 *   gnum_stride --> stride (in integers) between element global numbers
 *   gnum        --> pointer to global numbers
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_get_gnum_pointer(cs_all_to_all_t   *d,
                               size_t            *gnum_stride,
                               cs_gnum_t        **gnum);

/*----------------------------------------------------------------------------
 * Get current type of all-to-all distributor algorithm choice.
 *
 * returns:
 *   current type of all-to-all distributor algorithm choice
 *---------------------------------------------------------------------------*/

cs_all_to_all_type_t
cs_all_to_all_get_type(void);

/*----------------------------------------------------------------------------
 * Set current type of all-to-all distributor algorithm choice.
 *
 * parameters:
 *   t <-- type of all-to-all distributor algorithm choice to select
 *---------------------------------------------------------------------------*/

void
cs_all_to_all_set_type(cs_all_to_all_type_t  t);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Log performance information relative to instrumented all-to-all
 * distribution.
 *----------------------------------------------------------------------------*/

void
cs_all_to_all_log_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ALL_TO_ALL_H__ */
