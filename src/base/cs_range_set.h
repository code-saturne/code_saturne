#ifndef __CS_RANGE_SET_H__
#define __CS_RANGE_SET_H__

/*============================================================================
 * Operations related to handling of an owning rank for distributed entities.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "fvm_group.h"
#include "fvm_selector.h"
#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_numbering.h"

#include "cs_mesh_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*! Structure for range set management */

typedef struct {

  cs_lnum_t  n_elts[2];             /*!< Number of associated local elements
                                         (0: local range, 1: total) */

  cs_gnum_t  l_range[2];            /*!< global id range assigned to local
                                         rank: [start, past-the-end[ */
  const cs_gnum_t  *g_id;           /*!< global id assigned to elements
                                         (possibly shared) */
  cs_gnum_t        *_g_id;          /*!< global id assigned to elements
                                         (private) */

  const cs_interface_set_t  *ifs;   /*!< Associated interface set, or NULL */
  const cs_halo_t           *halo;  /*!< Associated halo, or NULL */

} cs_range_set_t;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
                    cs_gnum_t                 *g_id);

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
 * \param[in]   balance      try to balance shared elements across ranks ?
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
                    cs_gnum_t                  g_id_base);

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
                                cs_gnum_t                 *g_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a range set structure.
 *
 * \param[in, out]  rs  pointer to pointer to structure to destroy
 */
/*----------------------------------------------------------------------------*/

void
cs_range_set_destroy(cs_range_set_t  **rs);

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
                               void                  *val);

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
                  void                  *val);

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
                    void                  *dest_val);

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
                     void                  *dest_val);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RANGE_SET_H__ */
