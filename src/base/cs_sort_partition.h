#ifndef __FVM_SORT_PARTITION_H__
#define __FVM_SORT_PARTITION_H__

/*============================================================================
 * Data partitioning for parallel sort or order operations.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief function pointer for conversion of a double precision value in
 * range [0, 1] to a given element.
 *
 * \remark for simplicity, calling functions assume the input is shared
 *  between \ref cs_sort_partition_s_to_elt_t and
 * \ref cs_sort_partition_compare_t functions (as they use both).
 *
 * \param[in]   s      coordinate between 0 and 1
 * \param[out]  elt    pointer to element
 * \param[in]   input  pointer to optional (untyped) value or structure.
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_sort_partition_s_to_elt_t) (double        s,
                                void         *elt,
                                const void   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief function pointer for comparison of 2 elements.
 *
 * This function is the same type as that used by qsort_r.
 *
 * \remark for simplicity, calling functions assume the input is shared
 *  between \ref cs_sort_partition_s_to_elt_t and
 * \ref cs_sort_partition_compare_t functions (as they use both).
 *
 * \param[in]  elt1   coordinate between 0 and 1
 * \param[in]  elt2   pointer to optional (untyped) value or structure.
 * \param[in]  input  pointer to optional (untyped) value or structure.
 *
 * \return < 0 if elt1 < elt2, 0 if elt1 == elt2, > 0 if elt1 > elt2
 */
/*----------------------------------------------------------------------------*/

typedef int
(cs_sort_partition_compare_t) (const void  *elt1,
                               const void  *elt2,
                               const void  *input);

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine to which rank data elements should be sent for parallel
 *        sorting or ordering.
 *
 * \param[in]   sampling_factor  number of samples per rank
 * \param[in]   elt_size         size associated with each element
 * \param[in]   n_elts           number of elements to be indexed
 * \param[in]   elts             array of elements
 * \param[in]   weight           optional weight of each element, or NULL
 * \param[in]   order            ordering array
 * \param[out]  dest_rank_id     destination rank id (size: n_elts)
 * \param[in]   s_to_elt         coordinate to element conversion function
 * \param[in]   compare          comparison function
 * \param[in]   f_input          optional input to s_to_elt and compare, or NULL
 * \param[in]   comm             MPI communicator on which we build the
 *                               global index
 */
/*----------------------------------------------------------------------------*/

void
cs_sort_partition_dest_rank_id(cs_lnum_t                      sampling_factor,
                               size_t                         elt_size,
                               cs_lnum_t                      n_elts,
                               const void                    *elts,
                               const cs_lnum_t               *weight,
                               const cs_lnum_t                order[],
                               int                            dest_rank_id[],
                               cs_sort_partition_s_to_elt_t   s_to_elt,
                               cs_sort_partition_compare_t    compare,
                               const void                    *f_input,
                               MPI_Comm                       comm);

#endif /* if HAVE_MPI */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_SORT_PARTITION_H__ */
