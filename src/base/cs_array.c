/*============================================================================
 * Array handling utilities.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_array.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_array.c
        Array handling utilities.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy an array (ref) into another array (dest) and apply an
 *        inderection at the same time. Array with stride > 1 are assumed to be
 *        interlaced. The indirectino is applied to ref
 *
 * \param[in]      n_elts   number of elements in the array
 * \param[in]      stride   number of values for each element
 * \param[in]      elt_ids  indirection list
 * \param[in]      ref      reference values to copy
 * \param[in, out] dest     array storing values after applying the indirection
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_copy_with_indirection(cs_lnum_t        n_elts,
                                    int              stride,
                                    const cs_lnum_t  elt_ids[],
                                    const cs_real_t  ref[],
                                    cs_real_t        dest[])
{
  if (n_elts < 1)
    return;

  assert(stride > 0);
  assert(ref != NULL && dest != NULL);

  if (elt_ids == NULL)
    memcpy(dest, ref, sizeof(cs_real_t)*stride*n_elts);

  else {

    if (stride == 1) {

#     pragma omp parallel for if (n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_elts; i++)
        dest[i] = ref[elt_ids[i]];

    }
    else {

#     pragma omp parallel for if (n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_elts; i++) {

        const cs_real_t  *_ref = ref + elt_ids[i]*stride;
        cs_real_t  *_dest = dest + i*stride;

        for (int k = 0; k < stride; k++)
          _dest[k] = _ref[k];

      }

    } /* stride > 1 */

  } /* elt_ids != NULL */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy real values from an array to another of the same dimensions.
 *
 * \param[in]   n_elts  number of associated elements
 * \param[in]   dim     associated dimension
 * \param[in]   src     source array values (size: n_elts*dim)
 * \param[out]  dest    destination array values (size: n_elts*dim)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_copy(cs_lnum_t        n_elts,
                   cs_lnum_t        dim,
                   const cs_real_t  src[],
                   cs_real_t        dest[restrict])
{
  const cs_lnum_t _n_elts = dim * n_elts;

# pragma omp parallel for if (_n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < _n_elts; ii++)
    dest[ii] = src[ii];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant scalar value to an array.
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      ref_val  value to assign
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_scalar(cs_lnum_t  n_elts,
                         cs_real_t  ref_val,
                         cs_real_t  a[])
{
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++)
    a[ii] = ref_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant vector to an array of stride 3 which is interlaced
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      ref_val  vector to assign
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_vector(cs_lnum_t         n_elts,
                         const cs_real_t   ref_val[3],
                         cs_real_t        *a)
{
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    cs_real_t  *_a = a + 3*ii;
    _a[0] = ref_val[0], _a[1] = ref_val[1], _a[2] = ref_val[2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant vector of size 6 (optimzed way to define a
 *        symmetric tensor) to an array (of stride 6) which is interlaced
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      ref_val  vector to assign
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_symm_tensor(cs_lnum_t         n_elts,
                              const cs_real_t   ref_val[6],
                              cs_real_t        *a)
{
  if (n_elts < 1)
    return;

  assert(ref_val != NULL);

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    cs_real_t  *_a = a + 6*ii;
    _a[0] = ref_val[0], _a[1] = ref_val[1], _a[2] = ref_val[2];
    _a[3] = ref_val[3], _a[4] = ref_val[4], _a[5] = ref_val[5];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant 3x3 tensor to an array (of stride 9) which is
 *        interlaced
 *
 * \param[in]      n_elts    number of elements
 * \param[in]      ref_tens  tensor to assign
 * \param[in, out] a         array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_tensor(cs_lnum_t         n_elts,
                         const cs_real_t   ref_tens[3][3],
                         cs_real_t        *a)
{
  if (n_elts < 1)
    return;

  assert(ref_tens != NULL);

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

    cs_real_t  *_a = a + 9*ii;

    _a[0] = ref_tens[0][0], _a[1] = ref_tens[0][1], _a[2] = ref_tens[0][2];
    _a[3] = ref_tens[1][0], _a[4] = ref_tens[1][1], _a[5] = ref_tens[1][2];
    _a[6] = ref_tens[2][0], _a[7] = ref_tens[2][1], _a[8] = ref_tens[2][2];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant value to an array (deprecated function).
 *
 * \param[in]   n_elts  number of associated elements
 * \param[in]   dim     associated dimension
 * \param[in]   v       value to assign
 * \param[out]  a       array values (size: n_elts*dim]
 */
/*----------------------------------------------------------------------------*/

void
cs_array_set_value_real(cs_lnum_t  n_elts,
                        cs_lnum_t  dim,
                        cs_real_t  v,
                        cs_real_t  a[])
{
  cs_array_real_set_scalar(dim*n_elts, v, a);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
