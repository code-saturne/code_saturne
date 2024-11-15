/*============================================================================
 * Array handling utilities.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
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
  \file cs_array.cpp
        Array handling utilities.
*/

static const size_t  size_of_real = sizeof(cs_real_t);
static const size_t  size_of_lnum = sizeof(cs_lnum_t);
static const size_t  size_of_flag = sizeof(cs_flag_t);

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
 * \brief Assign true to all elements of an array. Case of an array of booleans
 *
 * \param[in]      size    total number of elements to set
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_bool_fill_true(cs_lnum_t  size,
                        bool       a[])
{
# pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < size; i++)
    a[i] = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign false to all elements of an array. Case of an array of booleans
 *
 * \param[in]      size    total number of elements to set
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_bool_fill_false(cs_lnum_t  size,
                         bool       a[])
{
# pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < size; i++)
    a[i] = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign zero to all elements of an array. Case of a cs_flag_t array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array of flags to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_flag_fill_zero(cs_lnum_t  size,
                        cs_flag_t  a[])
{
  if (cs_glob_n_threads > 1) {
#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      a[i] = 0;
  }
  else
    memset(a, 0, size*size_of_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign zero to all elements of an array. Case of a cs_lnum_t array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_lnum_fill_zero(cs_lnum_t  size,
                        cs_lnum_t  a[])
{
  if (cs_glob_n_threads > 1) {
#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      a[i] = 0;
  }
  else
    memset(a, 0, size*size_of_lnum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy values from an array of cs_lnum_t type to another of the same
 *        dimensions.
 *
 * \param[in]  size  number of elements * dimension
 * \param[in]  src   source array values (size: n_elts*dim)
 * \param[out] dest  destination array values (size: n_elts*dim)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_lnum_copy(cs_lnum_t       size,
                   const cs_lnum_t src[],
                   cs_lnum_t       dest[])
{
  if (cs_glob_n_threads > 1) {

    cs_lnum_t *restrict _dest = dest;

#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < size; ii++)
      _dest[ii] = src[ii];

  }
  else
    memcpy(dest, src, size_of_lnum*size);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign the value "num" to all elements of an array. Case of a
 *        cs_lnum_t array.
 *
 * \param[in]      size    total number of elements to set
 * \param[in]      num     value to set
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_lnum_set_value(cs_lnum_t  size,
                        cs_lnum_t  num,
                        cs_lnum_t  a[])
{
# pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < size; i++)
    a[i] = num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign the value "num" to an array on a selected subset of elements.
 *        if elt_ids is null, then one recovers the function
 *        \ref cs_array_lnum_set_value
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      num      value to set
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_lnum_set_value_on_subset(cs_lnum_t        n_elts,
                                  const cs_lnum_t  elt_ids[],
                                  cs_lnum_t        num,
                                  cs_lnum_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_lnum_set_value(n_elts, num, a);

  else {
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++)
      a[elt_ids[ii]] = num;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign zero to all elements of an array. Case of a int array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_int_fill_zero(cs_lnum_t  size,
                       int        a[])
{
  if (cs_glob_n_threads > 1) {
#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      a[i] = 0;
  }
  else
    memset(a, 0, size*size_of_lnum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign the value "num" to all elements of an array. Case of a
 *        int array.
 *
 * \param[in]      size    total number of elements to set
 * \param[in]      num     value to set
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_int_set_value(cs_lnum_t  size,
                       int        num,
                       int        a[])
{
# pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < size; i++)
    a[i] = num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign the value "num" to an array on a selected subset of elements.
 *        if elt_ids is null, then one recovers the function
 *        \ref cs_array_int_set_value
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      num      value to set
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_int_set_value_on_subset(cs_lnum_t        n_elts,
                                 const cs_lnum_t  elt_ids[],
                                 int              num,
                                 int              a[])
{
  if (elt_ids == nullptr)
    cs_array_int_set_value(n_elts, num, a);

  else {
    int *restrict _a = a;
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++)
      _a[elt_ids[ii]] = num;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy an array ("ref") into another array ("dest") on possibly only a
 *        part of the array(s). Array with stride > 1 are assumed to be
 *        interlaced. The subset of elements on which working is defined by
 *        "elt_ids". The way to apply the subset is defined with the parameter
 *        "mode" as follows:
 *        - Only the "ref" array if mode = 0 (CS_ARRAY_SUBSET_IN)
 *        - Only the "dest" array if mode = 1 (CS_ARRAY_SUBSET_OUT)
 *        - Both "ref" and "dest" arrays if mode = 2 (CS_ARRAY_SUBSET_INOUT)
 *
 *        It elt_ids is null or mode < 0 (CS_ARRAY_SUBSET_NULL), then the
 *        behavior is as \ref cs_array_real_copy
 *
 *        One assumes that all arrays are allocated with a correct size.
 *
 * \param[in]      n_elts   number of elements in the array
 * \param[in]      stride   number of values for each element
 * \param[in]      mode     apply the subset ids to which array(s)
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref      reference values to copy
 * \param[in, out] dest     array storing values after applying the indirection
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_copy_subset(cs_lnum_t         n_elts,
                          int               stride,
                          const cs_lnum_t   elt_ids[],
                          int               mode,
                          const cs_real_t   ref[],
                          cs_real_t         dest[])
{
  if (n_elts < 1)
    return;

  assert(stride > 0);
  assert(ref != nullptr && dest != nullptr);

  if (elt_ids == nullptr)
    cs_array_real_copy(n_elts*stride, ref, dest);

  else {

    switch (mode) {

    case CS_ARRAY_SUBSET_IN: /* Indirection is applied to ref */
      if (stride == 1) {

        cs_real_t *restrict _dest = dest;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++)
          _dest[i] = ref[elt_ids[i]];

      }
      else {

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {

          const cs_real_t  *_ref = ref + elt_ids[i]*stride;
          cs_real_t  *restrict _dest = dest + i*stride;

          for (int k = 0; k < stride; k++)
            _dest[k] = _ref[k];

        }

      } /* stride > 1 */
      break;

    case CS_ARRAY_SUBSET_OUT: /* Indirection is applied to dest */
      if (stride == 1) {

        cs_real_t *restrict _dest = dest;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++)
          _dest[elt_ids[i]] = ref[i];

      }
      else {

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {

          const cs_real_t  *_ref = ref + i*stride;
          cs_real_t  *restrict _dest = dest + elt_ids[i]*stride;

          for (int k = 0; k < stride; k++)
            _dest[k] = _ref[k];

        }

      } /* stride > 1 */
      break;

    case CS_ARRAY_SUBSET_INOUT: /* Indirection is applied to ref/dest */
      if (stride == 1) {

        cs_real_t *restrict _dest = dest;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _dest[elt_id] = ref[elt_id];
        }

      }
      else {

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {

          const cs_lnum_t  shift = elt_ids[i]*stride;
          const cs_real_t  *_ref = ref + shift;
          cs_real_t  *restrict _dest = dest + shift;

          for (int k = 0; k < stride; k++)
            _dest[k] = _ref[k];

        }

      } /* stride > 1 */
      break;

    default: /* No indirection */
      memcpy(dest, ref, size_of_real*stride*n_elts);
      break;

    } /* Switch on the indirection mode */

  } /* elt_ids != nullptr */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy real values from an array to another of the same dimensions.
 *
 * \param[in]   size    number of elements * dimension
 * \param[in]   src     source array values (size: n_elts*dim)
 * \param[out]  dest    destination array values (size: n_elts*dim)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_copy(cs_lnum_t         size,
                   const cs_real_t   src[],
                   cs_real_t         dest[])
{
  if (cs_glob_n_threads > 1) {

    cs_real_t *restrict _dest = dest;

#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < size; ii++)
      _dest[ii] = src[ii];

  }
  else
    memcpy(dest, src, size_of_real*size);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Multiply each value by a scaling factor s.t. dest *= scaling_factor
 *        If elt_ids is non-null, one applies an indirection.
 *        A stride can also be applied. One assumes an interlaced array.
 *
 * \param[in]  n_elts          number of elements
 * \param[in]  stride          number of values for each element
 * \param[in]  elt_ids         list of ids in the subset or null (size: n_elts)
 * \param[in]  scaling_factor  value of the scaling factor
 * \param[out] dest            destination array values
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_scale(cs_lnum_t         n_elts,
                    int               stride,
                    const cs_lnum_t  *elt_ids,
                    cs_real_t         scaling_factor,
                    cs_real_t         dest[])
{
  if (elt_ids == nullptr) {

    cs_real_t *restrict _dest = dest;

#   pragma omp parallel for if (n_elts*stride > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts*stride; ii++)
      _dest[ii] *= scaling_factor;

  }
  else {

    if (stride == 1) {

      cs_real_t *restrict _dest = dest;

#     pragma omp parallel for if (n_elts > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_elts; ii++)
        _dest[elt_ids[ii]] *= scaling_factor;

    }
    else {

      assert(stride > 0);
#     pragma omp parallel for if (n_elts > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
        cs_real_t  *restrict _dest = dest + stride*elt_ids[ii];
        for (int k = 0; k < stride; k++)
          _dest[k] *= scaling_factor;
      }

    }

  } /* elt_ids */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add in place an array s.t. r += l_add
 *
 * \param[in]  n_elts   number of elements
 * \param[in]  l_add    array to add
 * \param[out] r        destination array values
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_padd(cs_lnum_t       n_elts,
                   const cs_real_t l_add[],
                   cs_real_t       r[])
{
  if (n_elts < 1)
    return;

  cs_real_t *restrict _r = r;
#pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++)
    _r[ii] += l_add[ii];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant value of dim "stride" to an interlaced array
 *        sharing the same stride
 *
 * \param[in]      n_elts    number of elements
 * \param[in]      stride    number of values for each element
 * \param[in]      ref_val   list of values to assign (size: stride)
 * \param[in, out] a         array to set (size: n_elts*stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_value(cs_lnum_t        n_elts,
                        int              stride,
                        const cs_real_t  ref_val[],
                        cs_real_t        a[])
{
  if (n_elts < 1)
    return;

  assert(ref_val != nullptr);

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

    cs_real_t  *restrict _a = a + stride*ii;
    for (int k = 0; k < stride; k++)
      _a[k] = ref_val[k];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a weighted constant value of dim "stride" to an interlaced
 *        array sharing the same stride. Apply a weight for each element. This
 *        weight is constant for each component of an element.
 *
 * \param[in]      n_elts    number of elements
 * \param[in]      stride    number of values for each element
 * \param[in]      ref_val   list of values to assign (size: stride)
 * \param[in]      weight    values of the weight to apply (size: n_elts)
 * \param[in, out] a         array to set (size: n_elts*stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_wvalue(cs_lnum_t        n_elts,
                         int              stride,
                         const cs_real_t  ref_val[],
                         const cs_real_t  weight[],
                         cs_real_t        a[])
{
  if (n_elts < 1)
    return;

  assert(ref_val != nullptr && weight != nullptr);

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

    const cs_real_t  w = weight[ii];
    cs_real_t  *restrict _a = a + stride*ii;

    for (int k = 0; k < stride; k++)
      _a[k] = w*ref_val[k];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant value of dim "stride" to an interlaced array
 *        sharing the same stride. Only a subset of elements are considered.
 *        If elt_ids is null, then one recovers the function
 *        \ref cs_array_real_set_value
 *
 * \param[in]      n_elts    number of elements
 * \param[in]      stride    number of values for each element
 * \param[in]      elt_ids   list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val   list of values to assign (size: stride)
 * \param[in, out] a         array to set (size >= n_elts * stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_value_on_subset(cs_lnum_t        n_elts,
                                  int              stride,
                                  const cs_lnum_t  elt_ids[],
                                  const cs_real_t  ref_val[],
                                  cs_real_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_value(n_elts, stride, ref_val, a);

  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

      cs_real_t  *restrict _a = a + stride*elt_ids[ii];
      for (int k = 0; k < stride; k++)
        _a[k] = ref_val[k];

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a weighted constant value of dim "stride" to an interlaced
 *        array sharing the same stride. Only a subset of elements are
 *        considered.  If elt_ids is null, then one recovers the function \ref
 *        cs_array_real_set_wvalue Apply a weight for each element. This
 *        weight is constant for each component of an element.
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      stride   number of values for each element
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val  list of values to assign (size: stride)
 * \param[in]      weight   values of the weight to apply (size >= n_elts)
 * \param[in, out] a        array to set (size >= n_elts*stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_wvalue_on_subset(cs_lnum_t        n_elts,
                                   int              stride,
                                   const cs_lnum_t  elt_ids[],
                                   const cs_real_t  ref_val[],
                                   const cs_real_t  weight[],
                                   cs_real_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_wvalue(n_elts, stride, ref_val, weight, a);

  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

      const cs_lnum_t  elt_id = elt_ids[ii];
      const cs_real_t  w = weight[elt_id];
      cs_real_t  *restrict _a = a + stride*elt_id;

      for (int k = 0; k < stride; k++)
        _a[k] = w*ref_val[k];

    }

  }
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
 * \brief Assign a weighted constant scalar value to an array.
 *        The weight array has the same size as the array "a".
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      ref_val  value to assign
 * \param[in]      weight   values of the weight to apply
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_wscalar(cs_lnum_t        n_elts,
                          cs_real_t        ref_val,
                          const cs_real_t  weight[],
                          cs_real_t        a[])
{
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++)
    a[ii] = ref_val * weight[ii];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant scalar value to an array on a selected subset of
 *        elements. If elt_ids is null, then one recovers the function
 *        cs_array_real_set_scalar
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val  value to assign
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_scalar_on_subset(cs_lnum_t        n_elts,
                                   const cs_lnum_t  elt_ids[],
                                   cs_real_t        ref_val,
                                   cs_real_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_scalar(n_elts, ref_val, a);

  else {

    cs_real_t *restrict _a = a;
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++)
      _a[elt_ids[ii]] = ref_val;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a weighted constant scalar value to an array on a selected
 *        subset of elements. If elt_ids is null, then one recovers the function
 *        cs_array_real_set_wscalar
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val  value to assign
 * \param[in]      weight   values of weights to apply
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_wscalar_on_subset(cs_lnum_t        n_elts,
                                    const cs_lnum_t  elt_ids[],
                                    cs_real_t        ref_val,
                                    const cs_real_t  weight[],
                                    cs_real_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_wscalar(n_elts, ref_val, weight, a);

  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
      const cs_lnum_t  elt_id = elt_ids[ii];
      a[elt_id] = ref_val * weight[elt_id];
    }

  }
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
                         cs_real_t         a[])
{
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    cs_real_t  *restrict _a = a + 3*ii;
    _a[0] = ref_val[0], _a[1] = ref_val[1], _a[2] = ref_val[2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a weighted constant vector value to an interlaced array (of
 *        stride 3). The array of weights has the same size as the array "a".
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      ref_val  vector to assign
 * \param[in]      weight   values of the weight to apply
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_wvector(cs_lnum_t        n_elts,
                          const cs_real_t  ref_val[3],
                          const cs_real_t  weight[],
                          cs_real_t        a[])
{
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    const cs_real_t  w = weight[ii];
    cs_real_t  *_a = a + 3*ii;
    _a[0] = w*ref_val[0], _a[1] = w*ref_val[1], _a[2] = w*ref_val[2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant vector to an interlaced array (of stride 3) on a
 *        selected subset of elements. If elt_ids is null, then one recovers the
 *        function cs_array_real_set_vector
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val  vector to assign
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_vector_on_subset(cs_lnum_t        n_elts,
                                   const cs_lnum_t  elt_ids[],
                                   const cs_real_t  ref_val[3],
                                   cs_real_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_vector(n_elts, ref_val, a);

  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
      cs_real_t  *_a = a + 3*elt_ids[ii];
      _a[0] = ref_val[0], _a[1] = ref_val[1], _a[2] = ref_val[2];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a weighted constant vector value to an interlaced array (of
 *        stride 3). The subset selection is given by elt_ids. If null, then
 *        one recovers the function \ref cs_array_real_set_wvector
 *        The array of weights has the same size as the array "a".
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val  vector to assign
 * \param[in]      weight   values of the weight to apply
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_wvector_on_subset(cs_lnum_t        n_elts,
                                    const cs_lnum_t  elt_ids[],
                                    const cs_real_t  ref_val[3],
                                    const cs_real_t  weight[],
                                    cs_real_t        a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_wvector(n_elts, ref_val, weight, a);

  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
      const cs_lnum_t  id = elt_ids[ii];
      const cs_real_t  w = weight[id];
      cs_real_t  *restrict _a = a + 3*id;
      _a[0] = w*ref_val[0], _a[1] = w*ref_val[1], _a[2] = w*ref_val[2];
    }

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
cs_array_real_set_tensor(cs_lnum_t        n_elts,
                         const cs_real_t  ref_tens[3][3],
                         cs_real_t        a[])
{
  if (n_elts < 1)
    return;

  assert(ref_tens != nullptr);

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

    cs_real_t  *restrict _a = a + 9*ii;

    _a[0] = ref_tens[0][0], _a[1] = ref_tens[0][1], _a[2] = ref_tens[0][2];
    _a[3] = ref_tens[1][0], _a[4] = ref_tens[1][1], _a[5] = ref_tens[1][2];
    _a[6] = ref_tens[2][0], _a[7] = ref_tens[2][1], _a[8] = ref_tens[2][2];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant 3x3 tensor to an interlaced array (of stride 9) on
 *        a subset of elements. If elt_ids is null, then one recovers the
 *        function cs_array_real_set_tensor
 *
 * \param[in]      n_elts    number of elements
 * \param[in]      elt_ids   list of ids defining the subset or nullptr
 * \param[in]      ref_tens  tensor to assign
 * \param[in, out] a         array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_tensor_on_subset(cs_lnum_t         n_elts,
                                   const cs_lnum_t   elt_ids[],
                                   const cs_real_t   ref_tens[3][3],
                                   cs_real_t         a[])
{
  if (elt_ids == nullptr)
    cs_array_real_set_tensor(n_elts, ref_tens, a);

  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {

      cs_real_t  *restrict _a = a + 9*elt_ids[ii];

      _a[0] = ref_tens[0][0], _a[1] = ref_tens[0][1], _a[2] = ref_tens[0][2];
      _a[3] = ref_tens[1][0], _a[4] = ref_tens[1][1], _a[5] = ref_tens[1][2];
      _a[6] = ref_tens[2][0], _a[7] = ref_tens[2][1], _a[8] = ref_tens[2][2];

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign zero to all elements of an array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_fill_zero(cs_lnum_t  size,
                        cs_real_t  a[])
{
  if (cs_glob_n_threads > 1) {

#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      a[i] = 0;

  }
  else
    memset(a, 0, size*size_of_real);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant value to an array (deprecated function).
 *
 * \param[in]   n_elts  number of associated elements
 * \param[in]   dim     associated dimension
 * \param[in]   v       value to assign
 * \param[out]  a       array values (size: n_elts*dim)
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
