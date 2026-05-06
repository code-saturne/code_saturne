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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_algorithm.h"

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_algorithm.cpp
        Various base algorithms.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*--------------------------------------------------------------------------*/
/*!
 * \brief Transform a count to an index in place, OpenMP threaded version.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * \param[in]       n     number of elements
 * \param[in, out]  a <-> count in, index out (size: n+1)
 */
/*--------------------------------------------------------------------------*/

#if defined(HAVE_OPENMP)

static void
_count_to_index_inplace_omp(int        n_threads,
                            cs_lnum_t  n,
                            cs_lnum_t  a[])
{
  assert(a != nullptr && n > 0);

  constexpr int cl_size_max = 128;  // Max expected cache line size
  constexpr int n_t_noalloc = 32;

  // Space partial counts by stride to try to avoid false sharing.
  constexpr cs_lnum_t stride = (cl_size_max / sizeof(cs_lnum_t))
                               + (cl_size_max % sizeof(cs_lnum_t) == 0) ? 0 : 1;
  cs_lnum_t partial_sum_[cl_size_max * n_t_noalloc];
  cs_lnum_t *partial_sum = partial_sum_;
  if (n_threads > n_t_noalloc)
    CS_MALLOC(partial_sum, n_threads*stride, cs_lnum_t);

  // Compute partial sums (only reading from a, not writing to it)

  #pragma omp parallel shared(partial_sum) num_threads(n_threads)
  {
    int t_id = omp_get_thread_num();
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                           &t_s_id, &t_e_id);

    cs_lnum_t l_count = 0;
    for (cs_lnum_t i = t_s_id; i < t_e_id; i++)
      l_count += a[i];

    partial_sum[t_id*stride] = l_count;
  }

  // Serial inclusive scan of partial sums.

  {
    cs_lnum_t s = partial_sum[0];
    for (int j = 1; j < n_threads; j++) {
      s += partial_sum[j*stride];
      partial_sum[j*stride] = s;
    }
  }

  // Now finalize computation of index.

  #pragma omp parallel shared(partial_sum) num_threads(n_threads)
  {
    int t_id = omp_get_thread_num();
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                           &t_s_id, &t_e_id);

    cs_lnum_t s = (t_id == 0) ? 0 : partial_sum[t_id-1];

    for (cs_lnum_t i = t_s_id; i < t_e_id; i++) {
      cs_lnum_t c = a[i];
      a[i] = s;
      s += c;
    }

    if (t_id == n_threads-1)  // Last thread handles past-the end value
      a[n] = s;
  }

  if (partial_sum != partial_sum_)
    CS_FREE(partial_sum);
}

#endif // defined(HAVE_OPENMP)

/*--------------------------------------------------------------------------*/
/*!
 * \brief Transform a count to an index in place, Serial version.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * \param[in]       n     number of elements
 * \param[in, out]  a <-> count in, index out (size: n+1)
 */
/*--------------------------------------------------------------------------*/

static void
_count_to_index_inplace_serial(cs_lnum_t  n,
                               cs_lnum_t  a[])
{
  assert(a != nullptr && n > 0);

  cs_lnum_t s = 0;

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_lnum_t c = a[i];
    a[i] = s;
    s += c;
  }

  a[n] = s;
}

namespace cs {

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*--------------------------------------------------------------------------*/
/*!
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
               cs_lnum_t             a[])
{
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  if (ctx.use_gpu()) {
    // TODO: call GPU code for this and return;
  }
#endif

#if defined(HAVE_OPENMP)

  int n_threads = ctx.n_cpu_threads();
  if (n_threads == -1)
    n_threads = cs_parall_n_threads(n, CS_THR_MIN);

  if (n_threads > 1) {
    _count_to_index_inplace_omp(n_threads, n, a);
    return;
  }

#endif

  _count_to_index_inplace_serial(n, a);
}

/*----------------------------------------------------------------------------*/

} // namespace cs
