#ifndef CS_DEBUG_H
#define CS_DEBUG_H

/*============================================================================
 * Debug macros and utilities
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

END_C_DECLS

/*=============================================================================
 * Templated function definitions
 *============================================================================*/

#ifdef __cplusplus

/* Compute the unit in the last place (ULP) */
template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, T>::type
cs_diff_ulp(T x,
            T y)
{
  // Since `epsilon()` is the gap size (ULP, unit in the last place)
  // of floating-point numbers in interval [1, 2), we can scale it to
  // the gap size in interval [2^e, 2^{e+1}), where `e` is the exponent
  // of `x` and `y`.

  // If `x` and `y` have different gap sizes (which means they have
  // different exponents), we take the smaller one. Taking the bigger
  // one is also reasonable, I guess.
  const T m = fmin(std::fabs(x), std::fabs(y));

  // Subnormal numbers have fixed exponent, which is `min_exponent - 1`.
  const int exp = m < std::numeric_limits<T>::min()
    ? std::numeric_limits<T>::min_exponent - 1
    : std::ilogb(m);

  // We divide the absolute difference by
  // the epsilon times the exponent (1 ulp)
  return std::fabs(x - y) / std::ldexp(std::numeric_limits<T>::epsilon(), exp);
}

#endif // __cplusplus

/*----------------------------------------------------------------------------*/

#endif /* CS_DEBUG_H */
