/*============================================================================
 * Math functions for CUDA.
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

#pragma once

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base_cuda.h"

/*=============================================================================
 * Device function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the absolute value of a real value.
 *
 * \param[in]  x  value
 *
 * \return absolute value of the given value
 */
/*----------------------------------------------------------------------------*/

__device__ static double
cs_math_abs_cuda(double  x)
{
  return fabs(x);
}

__device__ static float
cs_math_abs_cuda(float  x)
{
  return fabsf(x);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of two vectors of 3 real values.
 *
 * \param[in]   u   vector of 3 real values
 * \param[in]   v   vector of 3 real values
 *
 * \return the resulting dot product u.v.
 */
/*----------------------------------------------------------------------------*/

template <typename T>
__device__ static cs_real_t
cs_math_3_dot_product_cuda(const T  u[3],
                           const T  v[3])
{
  T prod = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

  return prod;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Normalise a vector of 3 real values.
 *
 * To normalize in-place, \p vin and \p vout may point to the same array.
 *
 * \param[in]     vin           vector
 * \param[out]    vout          normalized vector
 */
/*----------------------------------------------------------------------------*/

template <typename T>
__device__ static void
cs_math_3_normalize_cuda(const T  in[3],
                         T        out[3])
{
  T norm = sqrt(  in[0]*in[0]
                + in[1]*in[1]
                + in[2]*in[2]);

  T inverse_norm =  1. / norm;

  out[0] = inverse_norm * in[0];
  out[1] = inverse_norm * in[1];
  out[2] = inverse_norm * in[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the square norm of a vector of 3 real values.
 *
 * \param[in]     v             vector of 3 real values
 *
 * \return square norm of v.
 */
/*----------------------------------------------------------------------------*/

template <typename T>
__device__ static T
cs_math_3_square_norm_cuda(const T  in[3])
{
  T norm = in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
  return norm;
}

/*----------------------------------------------------------------------------*/
