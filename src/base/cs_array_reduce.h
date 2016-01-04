#ifndef __CS_ARRAY_REDUCE_H__
#define __CS_ARRAY_REDUCE_H__

/*============================================================================
 * Common array reduction operations.
 *============================================================================*/

/*
  This file is part of the Code_Saturne Kernel, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1998-2016 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Kernel is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Kernel is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Kernel; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute sums of an n-dimensional cs_real_t array's components.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior.
 *
 * Note that for locations with elements shared across ranks
 * (such as interior faces and vertices), sums may be incorrect as
 * contributions from multiple ranks may be counted several times.
 *
 * parameters:
 *   n_elts     <-- number of local elements
 *   v_elt_list <-- optional list of parent elements on which values
 *                  are defined, or NULL
 *   dim        <-- local array dimension (max: 9)
 *   v          <-- pointer to array values
 *   vsum       --> resulting sum array (size: dim, or 4 if dim = 3)
 *----------------------------------------------------------------------------*/

void
cs_array_reduce_sum_l(cs_lnum_t         n_elts,
                      int               dim,
                      const cs_lnum_t  *v_elt_list,
                      const cs_real_t   v[],
                      double            vsum[]);

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of an
 * n-dimensional cs_real_t array's components.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior.
 *
 * Note that for locations with elements shared across ranks
 * (such as interior faces and vertices), sums may be incorrect as
 * contributions from multiple ranks may be counted several times.
 *
 * parameters:
 *   n_elts     <-- number of local elements
 *   dim        <-- local array dimension (max: 9)
 *   v_elt_list <-- optional list of parent elements on which values
 *                  are defined, or NULL
 *   v          <-- pointer to array values
 *   vmin       --> resulting min array (size: dim, or 4 if dim = 3)
 *   vmax       --> resulting max array (size: dim, or 4 if dim = 3)
 *   vsum       --> resulting sum array (size: dim, or 4 if dim = 3)
 *----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_stats_l(cs_lnum_t         n_elts,
                               int               dim,
                               const cs_lnum_t  *v_elt_list,
                               const cs_real_t   v[],
                               double            vmin[],
                               double            vmax[],
                               double            vsum[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local stats (minima, maxima, sum, weighted sum) of
 * an n-dimensional cs_real_t array's components for a given mesh location.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_elts     <-- number of local elements
 *   dim        <-- local array dimension (max: 9)
 *   v_elt_list <-- optional list of parent elements on which values
 *                  are defined, or NULL
 *   w_elt_list <-- optional list of parent elements on which weights
 *                  are defined, or NULL; if v_elt_list is defined
 *                  (ie. non-NULL),then w_elt_list = v_elt_list is assumed,
 *                  so this parameter is ignored
 *   v          <-- pointer to array values
 *   w          <-- pointer to weights
 *   vmin       --> resulting min array (size: dim, or 4 if dim = 3)
 *   vmax       --> resulting max array (size: dim, or 4 if dim = 3)
 *   vsum       --> resulting sum array (size: dim, or 4 if dim = 3)
 *   wsum       --> resulting weighted sum array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_stats_l_w(cs_lnum_t         n_elts,
                                 int               dim,
                                 const cs_lnum_t  *v_elt_list,
                                 const cs_lnum_t  *w_elt_list,
                                 const cs_real_t   v[],
                                 const cs_real_t   w[],
                                 double            vmin[],
                                 double            vmax[],
                                 double            vsum[],
                                 double            wsum[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ARRAY_REDUCE_H__ */
