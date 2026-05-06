#pragma once

/*============================================================================
 * Various base algorithms.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"

#include "base/cs_dispatch.h"
#include "base/cs_halo.h"
#include "alge/cs_matrix.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

namespace cs {

/*=============================================================================
 * Public inline functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given local id in a given array of
 *        ordered local ids, when the id might not be present.
 *
 * \param[in]  l_id_array size  array_size
 * \param[in]  l_id             local id to search for
 * \param[in]  l_id_array       ordered unique local ids array
 *
 * \return  index of l_id in l_id_array, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_lnum_t
l_binary_search(cs_lnum_t        l_id_array_size,
                cs_lnum_t        l_id,
                const cs_lnum_t  l_id_array[])
{
  if (l_id_array_size < 1)
    return -1;

  cs_lnum_t start_id = 0;
  cs_lnum_t end_id = l_id_array_size - 1;
  cs_lnum_t mid_id = (end_id -start_id) / 2;
  while (start_id < end_id) {
    if (l_id_array[mid_id] < l_id)
      start_id = mid_id + 1;
    else if (l_id_array[mid_id] > l_id)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }
  if (l_id_array[mid_id] != l_id)
    mid_id = -1;

  return mid_id;
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*--------------------------------------------------------------------------*/
/*
 * \brief Transform a count to an index in-place.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * \param[in]       n     number of elements
 * \param[in, out]  a <-> count in, index out (size: n+1)
 */
/*--------------------------------------------------------------------------*/

void
count_to_index(cs_dispatch_context  &ctx,
               cs_lnum_t             n,
               cs_lnum_t             a[]);

/*----------------------------------------------------------------------------*/

} // namespace cs
