#ifndef __CS_ARRAY_H__
#define __CS_ARRAY_H__

/*============================================================================
 * Array handling utilities.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_base_accel.h"
#if defined (__NVCC__)
#include "base/cs_array_cuda.h"
#include "base/cs_base_cuda.h"
#endif

#include "base/cs_dispatch.h"
#include "base/cs_mdspan.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Define the way to apply a subset. In case of a copy, the subset is related
   to the reference (input array) and/or the destination array (output
   array) */

#define CS_ARRAY_SUBSET_NULL  -1
#define CS_ARRAY_SUBSET_IN     0
#define CS_ARRAY_SUBSET_OUT    1
#define CS_ARRAY_SUBSET_INOUT  2

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public inline function prototypes
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(__cplusplus)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values to all elements of multiple arrays. ref_val is input
 *        as a pointer or an array
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      n_elts  total number of elements to set
 * \param[in]      ref_val value to assign
 * \param[out]     arrays  arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_value(const cs_lnum_t  n_elts,
                    const T         *ref_val,
                    Arrays&&...      arrays)
{
  /* Expand the parameter pack */
  T* array_ptrs[] = {arrays ... };

#if defined (__NVCC__)
  bool is_available_on_device = cs_check_device_ptr(ref_val);
  for (T* array : array_ptrs)
    is_available_on_device =  is_available_on_device
                           && (cs_check_device_ptr(array)
                               == CS_ALLOC_HOST_DEVICE_SHARED);

  if (is_available_on_device) {
    cudaStream_t stream_ = cs_cuda_get_stream(0);
    cs_arrays_set_value<T, stride>(stream_,
                                   n_elts,
                                   ref_val,
                                   arrays...);
    return;
  }
#endif

  auto set_value = [=](cs_lnum_t i_elt) {
    for (T* array : array_ptrs)
      memcpy(array + i_elt*stride, ref_val, stride*sizeof(T));
  };

  /* Loop through each index and assign values */
# pragma omp parallel for  if (n_elts >= CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++)
    set_value(i);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values to all elements of multiple arrays. ref_val is input
 *        as a scalar
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      n_elts  total number of elements to set
 * \param[in]      ref_val value to assign
 * \param[out]     arrays  arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_value(const cs_lnum_t  n_elts,
                    const T          ref_val,
                    Arrays&&...      arrays)
{

  /* Expand the parameter pack */
  T* array_ptrs[] = {arrays ... };

#if defined (__NVCC__)
  bool is_available_on_device = true;
  for (T* array : array_ptrs)
    is_available_on_device =  is_available_on_device
                           && (cs_check_device_ptr(array)
                               == CS_ALLOC_HOST_DEVICE_SHARED);

  if (is_available_on_device) {
    cudaStream_t stream_ = cs_cuda_get_stream(0);
    cs_arrays_set_value<T, stride>(stream_,
                                   n_elts,
                                   ref_val,
                                   arrays...);
    return;
  }
#endif

  auto set_value = [=](cs_lnum_t i_elt) {
    for (T* array : array_ptrs)
      for (size_t k = 0; k < stride; k++) {
        array[i_elt*stride + k] = ref_val;
      }
  };

  /* Loop through each index and assign values */
# pragma omp parallel for  if (n_elts >= CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++)
    set_value(i);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values on a selected subset of elements to multiple arrays.
 *        ref_val is input as a pointer or an array
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      n_elts   total number of elements to set
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val  value to assign
 * \param[out]     arrays   arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_value_on_subset(const cs_lnum_t  n_elts,
                              const cs_lnum_t *elt_ids,
                              const T         *ref_val,
                              Arrays&&...      arrays)
{
  if (n_elts < 1)
    return;

  if (elt_ids == NULL)
    cs_arrays_set_value<T, stride>(n_elts, ref_val, arrays...);
  else {
    /* Expand the parameter pack */
    T* array_ptrs[] = {arrays ... };

    auto set_value = [=](cs_lnum_t i_elt) {
      for (T* array : array_ptrs)
        memcpy(array + i_elt*stride, ref_val, stride*sizeof(T));
    };

    /* Loop through each index on the subset and assign values */
#   pragma omp parallel for  if (n_elts >= CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      set_value(elt_ids[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values on a selected subset of elements to multiple arrays.
 *        ref_val is input as a scalar
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      n_elts  total number of elements to set
 * \param[in]      elt_ids  list of ids in the subset or null (size: n_elts)
 * \param[in]      ref_val value to assign
 * \param[out]     arrays  arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_value_on_subset(const cs_lnum_t  n_elts,
                              const cs_lnum_t *elt_ids,
                              const T          ref_val,
                              Arrays&&...      arrays)
{
  if (n_elts < 1)
    return;

  if (elt_ids == NULL)
    cs_arrays_set_value<T, stride>(n_elts, ref_val, arrays...);
  else {

    /* Expand the parameter pack */
    T* array_ptrs[] = {arrays ... };

    auto set_value = [=](cs_lnum_t i_elt) {
      for (T* array : array_ptrs)
        for (size_t k = 0; k < stride; k++) {
          array[i_elt*stride + k] = ref_val;
        }
    };

   /* Loop through each index on the subset and assign values */
#   pragma omp parallel for  if (n_elts >= CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      set_value(elt_ids[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy values from an array to another of the same dimensions.
 *
 * Template parmeters.
 *                 T       type name
 *
 * \param[in]   size    number of elements * dimension
 * \param[in]   src     source array values
 * \param[out]  dest    destination array values
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_array_copy(const cs_lnum_t  size,
              const T*         src,
              T*               dest)
{
#if defined (__NVCC__)
  bool is_available_on_device =  cs_check_device_ptr(src)
                              && cs_check_device_ptr(dest);

  if (is_available_on_device) {
    cudaStream_t stream_ = cs_cuda_get_stream(0);
    cs_array_copy<T>(stream_, size, src, dest);
    return;
  }
#endif

# pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < size; ii++)
    dest[ii] = src[ii];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the difference diff = x - y. All arrays have the same
 * dimension.
 *
 *
 * Template parmeters.
 *                 T       type name
 *
 * \param[in]   size    number of elements * dimension
 * \param[in]   x       x array values
 * \param[in]   y       y array values
 * \param[out]  diff    difference array values
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_array_difference(const cs_lnum_t size, const T *x, const T *y, T *diff)
{
  cs_array_copy(size, x, diff);

#pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < size; ii++)
    diff[ii] -= y[ii];
}

#endif // __cplusplus

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign true to all elements of an array. Case of an array of booleans
 *
 * \param[in]      size    total number of elements to set
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_bool_fill_true(cs_lnum_t  size,
                        bool       a[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign false to all elements of an array. Case of an array of booleans
 *
 * \param[in]      size    total number of elements to set
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_bool_fill_false(cs_lnum_t  size,
                         bool       a[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign zero to all elements of an array. Case of a cs_flag_t array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array of flags to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_flag_fill_zero(cs_lnum_t  size,
                        cs_flag_t  a[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign zero to all elements of an array. Case of a cs_lnum_t array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_lnum_fill_zero(cs_lnum_t  size,
                        cs_lnum_t  a[]);

/*----------------------------------------------------------------------------*/
/*
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
                   cs_lnum_t       dest[]);

/*----------------------------------------------------------------------------*/
/*
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
                        cs_lnum_t  a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                  cs_lnum_t        a[]);
/*----------------------------------------------------------------------------*/
/*
 * \brief Assign zero to all elements of an array. Case of a int array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_int_fill_zero(cs_lnum_t  size,
                       int        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                       int        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                 int              a[]);

/*----------------------------------------------------------------------------*/
/*
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
                          cs_real_t         dest[]);

/*----------------------------------------------------------------------------*/
/*
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
                   cs_real_t         dest[]);

/*----------------------------------------------------------------------------*/
/*
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
                    cs_real_t         dest[]);

/*----------------------------------------------------------------------------*/
/*
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
                   cs_real_t       r[]);

/*----------------------------------------------------------------------------*/
/*
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
                        cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                         cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                  cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                   cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                         cs_real_t  a[]);

/*----------------------------------------------------------------------------*/
/*
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
                          cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                   cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                    cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                         cs_real_t         a[]);

/*----------------------------------------------------------------------------*/
/*
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
                          cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                   cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                    cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                         cs_real_t        a[]);

/*----------------------------------------------------------------------------*/
/*
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
                                   cs_real_t         a[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign zero to all elements of an array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_fill_zero(cs_lnum_t  size,
                        cs_real_t  a[]);

/*----------------------------------------------------------------------------*/
/*
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
                        cs_real_t  a[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*----------------------------------------------------------------------------*/

#if defined(__cplusplus)

namespace cs {

/*----------------------------------------------------------------------------*/
/*!
 * \class Templated data array class.
 *
 * @tparam T : data type
 * @tparam N : number of dimensions (int)
 * @tparam L : memory layout (cs::layout)
 */
/*----------------------------------------------------------------------------*/

template<class T, int N = 1, layout L = layout::right>
class array {

public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method leading to "empty container".
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array():
    _extent{0},
    _offset{0},
    _size(0),
    _owner(true),
    _data(nullptr),
    _mode(cs_alloc_mode)
  {
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using only size.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    cs_lnum_t     size, /*!<[in] size of array */
    cs_alloc_mode_t alloc_mode = cs_alloc_mode, /*!<[in] Memory allocation mode */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
 __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  :
    _extent{0},
    _offset{0},
    _size(size),
    _owner(true),
    _data(nullptr),
    _mode(alloc_mode)
  {
    /* Only usable for array */
    static_assert(N == 1,
                  "Instantiation using only total size only possible for "
                  " cs::array<T,1> or cs::array<T>");
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using dimensions.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    const cs_lnum_t(&dims)[N], /*!<[in] Array of dimensions sizes */
    cs_alloc_mode_t alloc_mode = cs_alloc_mode, /*!<[in] Memory allocation mode */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
 __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  :
    _extent{0},
    _offset{0},
    _owner(true),
    _data(nullptr),
    _mode(alloc_mode)
  {
    set_size_(dims);
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for 2D array based on sizes
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    cs_lnum_t     size1, /*!<[in] First size of array */
    cs_lnum_t     size2, /*!<[in] Second size of array */
    cs_alloc_mode_t alloc_mode = cs_alloc_mode, /*!<[in] Memory allocation mode */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
 __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  :
    _extent{0},
    _offset{0},
    _owner(true),
    _data(nullptr),
    _mode(alloc_mode)
  {
    /* Only usable for array */
    static_assert(N == 2,
                  "Instantiation using (size1, size2) only possible for "
                  " cs::array<T,2>");
    cs_lnum_t tmp_size[N] = {size1, size2};
    set_size_(tmp_size);
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for 3D array based on sizes
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    cs_lnum_t     size1, /*!<[in] First size of array */
    cs_lnum_t     size2, /*!<[in] Second size of array */
    cs_lnum_t     size3, /*!<[in] Third size of array */
    cs_alloc_mode_t alloc_mode = cs_alloc_mode, /*!<[in] Memory allocation mode */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
 __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  :
    _extent{0},
    _offset{0},
    _owner(true),
    _data(nullptr),
    _mode(alloc_mode)
  {
    /* Only usable for array */
    static_assert(N == 3,
                  "Instantiation using (size1, size2, size3) only possible "
                  " for cs::array<T,3>");
    cs_lnum_t tmp_size[N] = {size1, size2, size3};
    set_size_(tmp_size);
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for non owner version
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    cs_lnum_t  size,  /*!<[in] Size of array */
    T         *data,  /*!<[in] Pointer to data array */
    cs_alloc_mode_t alloc_mode = cs_alloc_mode, /*!<[in] Memory allocation mode */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
 __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  :
    _extent{0},
    _offset{0},
    _size(size),
    _owner(false),
    _mode(alloc_mode),
    _data(data)
  {
    static_assert(N == 1,
                  "Instantiation using only total size only possible for "
                  "cs::array<T,1> or cs::array<T>");
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using copy. May be a shallow copy.
   *
   * Important !
   * The default mode is shallow copy, since we want for lambda captures
   * to have a lightweight copy. Otherwise, with the deep copy as default,
   * this may lead to segmentation faults and unnecessary allocations.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    array&      other, /*!<[in] Reference of data array to copy */
    bool        deep_copy=false, /*!<[in] Make a deep copy (owner) or not. */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {

    set_size_(other._extent);
    _mode = other._mode;

    /* If shallow copy new instance is not owner. Otherwise same ownership
     * as original instance since we copy it.
     */
    _owner = (deep_copy) ? other._owner : false;

    if (_owner) {
      allocate_(file_name, line_number);
      copy_data(other._data);
    }
    else {
      CS_UNUSED(file_name);
      CS_UNUSED(line_number);
      _data = other._data;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Move constructor.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array
  (
    array&& other /*!<[in] Original reference to move */
  )
  : array()
  {
    swap(*this, other);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  ~array()
  {
    clear();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Class swap operator used for assignment or move.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  friend void
  swap
  (
    array& first, /*!<[in] First instance to swap */
    array& second /*!<[in] Second instance to swap */
  )
  {
    using std::swap;

    /* Swap the different members between the two references. */
    swap(first._extent, second._extent);
    swap(first._offset, second._offset);
    swap(first._size, second._size);
    swap(first._owner, second._owner);
    swap(first._data, second._data);
    swap(first._mode, second._mode);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Assignment operator.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array& operator=(array other)
  {
    swap(*this, other);

    return *this;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Clear data (empty container).
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  clear()
  {
    if (_owner) {
      CS_FREE(_data);
    }
    else {
      _data = nullptr;
    }
    set_size_(0);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get span view of array, same dimensions as array.
   *
   * \return mdspan view with same number of dimensions as array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan<T, N, L>
  view()
  {
    if (N == 1)
      return mdspan<T,N,L>(_size, _data);
    else
      return mdspan<T,N,L>(_extent, _data);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get span view of array (default), same dimensions as array.
   *
   * \return mdspan view with same number of dimensions as array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan<T, N, L>
  get_span()
  {
    if (N == 1)
      return mdspan<T,N,L>(_size, _data);
    else
      return mdspan<T,N,L>(_extent, _data);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get span view of array with a given array of dimensions (extent).
   *        If total size of span is different than that of the array an error
   *        is raised.
   *
   * \return mdspan view with given dimensions.
   */
  /*--------------------------------------------------------------------------*/

  template<int _N_, layout _L_ = L>
  CS_F_HOST_DEVICE
  mdspan<T, _N_, _L_>
  get_mdspan
  (
    const cs_lnum_t(&dims)[_N_] /*!<[in] Array of dimensions sizes */
  )
  {
    /* sanity check for total size */
    cs_lnum_t s = 1;
    for (int i = 0; i < _N_; i++)
      s *= dims[i];

    if (s != _size)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: requested span has total size of %ld instead of %ld "
                  "for this array.\n"),
                __func__, s, _size);

    return mdspan<T,_N_,_L_>(dims, _data);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set all values of the data array to a constant value.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_to_val
  (
    T               val,         /*!<[in] Value to set to entire data array. */
    const cs_lnum_t n_vals = -1 /*!<[in] Number of values to copy.
                                          If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    cs_dispatch_context ctx;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = val;
    });

    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set all values of the data array to a constant value while providing
   *        a dispatch context. It is up to the call to synchronize the context
   *        after this call.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_to_val
  (
    cs_dispatch_context &ctx,        /*!< Reference to dispatch context */
    T                    val,        /*!<[in] Value to set to entire data array. */
    const cs_lnum_t      n_vals = -1 /*!<[in] Number of values to copy.
                                              If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    /* No wait here since context is passed as argument */
    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = val;
    });
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set all values of the data array to 0.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  zero()
  {
    T _zero = static_cast<T>(0);

    cs_dispatch_context ctx;

    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = _zero;
    });

    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set all values of the data array to 0.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  zero
  (
    cs_dispatch_context &ctx /*!< Reference to dispatch context */
  )
  {
    T _zero = static_cast<T>(0);

    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = _zero;
    });
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Initializer method for empty containers.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  set_empty()
  {
    set_size_(0);
    _owner = false;
    _data = nullptr;
    _mode = cs_alloc_mode;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void reshape
  (
    cs_lnum_t    new_size,      /*!<[in] New size */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    static_assert(N == 1,
                  "Method reshape(size, ...) only possible for "
                  " cs::array<T> or cs::array<T,1>");
    assert(new_size >= 0);

    /* If same dimensions, nothing to do ... */
    if (new_size == _size)
      return;

    if (_owner) {
      set_size_(new_size);
      reallocate_(file_name, line_number);
    }

  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array and copy old data
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  reshape_and_copy
  (
    cs_lnum_t       new_size,     /*!<[in] New size */
    cs_lnum_t       size_to_keep, /*!<[in] Size of data to keep (-1 for all) */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    static_assert(N==1,
                  "Method reshape_and_copy(size, ...) only possible for "
                  " cs::array<T> or cs::array<T,1>");
    assert(new_size >= 0);
    assert(size_to_keep <= new_size && size_to_keep <= _size);

    /* If same dimensions, nothing to do ... */
    if (new_size == _size)
      return;

    if (!(size_to_keep <= new_size && size_to_keep <= _size))
      bft_error(__FILE__, __LINE__, 0,
                "%s: Data cannot be saved when new size is smaller than size to keep.\n",
                __func__);

    if (_owner) {
      /* If new size is larger, realloc is sufficient. */
      set_size_(new_size);
      reallocate_(file_name, line_number);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array and do not keep old data.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void reshape
  (
    const cs_lnum_t(&dims)[N], /*!<[in] Array of dimensions sizes */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    /* If same dimensions, nothing to do ... */
    bool same_dim = true;
    for (int i = 0; i < N; i++)
      if (dims[i] != _extent[i])
        same_dim = false;

    if (same_dim)
      return;

    if (_owner) {
      set_size_(dims);
      reallocate_(file_name, line_number);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array and copy old data
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void reshape_and_copy
  (
    const cs_lnum_t(&dims)[N],   /*!<[in] Array of dimensions sizes */
    cs_lnum_t size_to_keep = -1, /*!<[in] Size to keep when copying data */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    assert(size_to_keep <= dims[N-1] && size_to_keep <= _extent[N-1]);

    /* If same dimensions, nothing to do ... */
    bool same_dim = true;
    for (int i = 0; i < N; i++)
      if (dims[i] != _extent[i])
        same_dim = false;

    if (same_dim)
      return;

    if (_owner) {
      /* check that size_to_keep is either >0 or == -1 (keep all) */
      if (size_to_keep != 0) {
        /* Temporary copy */
        array<T,N,L> tmp(*this, false);

        /* Update instance size */
        set_size_(dims);
        reallocate_(file_name, line_number);

        /* Local copies for the context loop */
        cs_lnum_t new_o[N], new_e[N], old_o[N], old_e[N], max_e[N];

        for (int i = 0; i < N; i++) {
          new_o[i] = _offset[i];
          new_e[i] = _extent[i];
          old_o[i] = tmp._offset[i];
          old_e[i] = tmp._extent[i];

          max_e[i] = (new_e[i] > old_e[i]) ? old_e[i] : new_e[i];
        }
        if (max_e[N-1] > size_to_keep && size_to_keep > 0)
          max_e[N-1] = size_to_keep;

        T* old_data = tmp._data;
        T* new_data = _data;

        cs_lnum_t loop_size = tmp._size;

        /* Loop using dispatch */
        cs_dispatch_context ctx;

        ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
          cs_lnum_t idx[N];
          cs_lnum_t dummy = e_id;

          for (int i = N-1; i >= 1; i--) {
            idx[i] = dummy % old_e[i];
            dummy = (dummy - idx[i]) / old_e[i];
          }
          idx[0] = dummy;

          bool all_ok = true;
          for (int i = 0; i < N; i++) {
            if (idx[i] >= max_e[i]) {
              all_ok = false;
              break;
            }
          }

          if (all_ok) {
            cs_lnum_t o_id = 0;
            cs_lnum_t n_id = 0;

            for (int i = 0; i < N; i++) {
              o_id += idx[i] * old_o[i];
              n_id += idx[i] * new_o[i];
            }

            new_data[n_id] = old_data[o_id];
          }
        });

        ctx.wait();
      }
      else {
        reshape(dims, file_name, line_number);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array 2D (no explicit copy, only realloc)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  reshape
  (
    const cs_lnum_t size1,       /*!<[in] Size along first dimension */
    const cs_lnum_t size2,       /*!<[in] Size along second dimension */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    static_assert(N == 2,
                  "Method reshape(size1, size2, ...) only possible for "
                  " cs::array<T,2>");
    cs_lnum_t tmp_size[N] = {size1, size2};

    reshape(tmp_size, file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array and copy old data (2D)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  reshape_and_copy
  (
    const cs_lnum_t size1,  /*!<[in] Size along first dimension */
    const cs_lnum_t size2,  /*!<[in] Size along second dimension */
    cs_lnum_t size_to_keep, /*!<[in] Size to keep when copying data  (-1 for all) */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    static_assert(N == 2,
                  "Method reshape_and_copy(size1, size2, ...) only possible for "
                  " cs::array<T,2>");
    cs_lnum_t tmp_size[N] = {size1, size2};

    reshape_and_copy(tmp_size, size_to_keep, file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array 3D (no explicit copy, only realloc)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  reshape
  (
    const cs_lnum_t size1,       /*!<[in] Size along first dimension */
    const cs_lnum_t size2,       /*!<[in] Size along second dimension */
    const cs_lnum_t size3,       /*!<[in] Size along third dimension */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    static_assert(N == 3,
                  "Method reshape(size1, size2, size3, ...) only possible for "
                  " cs::array<T,3>");
    cs_lnum_t tmp_size[N] = {size1, size2, size3};

    reshape(tmp_size, file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array and copy old data (3D)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  reshape_and_copy
  (
    const cs_lnum_t size1,  /*!<[in] Size along first dimension */
    const cs_lnum_t size2,  /*!<[in] Size along second dimension */
    const cs_lnum_t size3,  /*!<[in] Size along third dimension */
    cs_lnum_t size_to_keep, /*!<[in] Size to keep when copying data  (-1 for all) */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  )
  {
    static_assert(N == 3,
                  "Method reshape_and_copy(size1, size2, size3, ...) only "
                  "possible for cs::array<T,3>");
    cs_lnum_t tmp_size[N] = {size1, size2, size3};

    reshape_and_copy(tmp_size, size_to_keep, file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set memory allocation mode.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_alloc_mode
  (
    cs_alloc_mode_t mode /*!<[in] Memory allocation mode. */
  )
  {
    if (_mode != mode) {
      /* Since allocation mode is changed, deallocate data */
      if (_size != 0)
        clear();

      _mode = mode;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on index.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T*
  sub_array
  (
    cs_lnum_t i /*!<[in] index of subarray */
  )
  {
    if (L == layout::right)
      return _data + i*_offset[0];
    else
      return _data + i*_offset[N-1];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on two indexes.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T*
  sub_array
  (
    cs_lnum_t i, /*!<[in] index of subarray */
    cs_lnum_t j  /*!<[in] index of subarray */
  )
  {
    static_assert(N > 2,
                  "sub_array(i,j) can only be called for N>2");
    if (L == layout::right)
      return _data + i*_offset[0] + j*_offset[1];
    else
      return _data + i*_offset[N-1] + j*_offset[N-2];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on three indexes.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T*
  sub_array
  (
    cs_lnum_t i, /*!<[in] index of subarray */
    cs_lnum_t j, /*!<[in] index of subarray */
    cs_lnum_t k  /*!<[in] index of subarray */
  )
  {
    static_assert(N > 3,
                  "sub_array(i,j,k) can only be called for N>3");
    if (L == layout::right)
      return _data + i*_offset[0] + j*_offset[1] + k*_offset[2];
    else
      return _data + i*_offset[N-1] + j*_offset[N-2] + k*_offset[N-3];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on index. (const)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T*
  sub_array
  (
    cs_lnum_t i /*!<[in] index of subarray */
  ) const
  {
    if (L == layout::right)
      return _data + i*_offset[0];
    else
      return _data + i*_offset[N-1];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on index.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T*
  sub_array
  (
    cs_lnum_t i, /*!<[in] index of subarray */
    cs_lnum_t j  /*!<[in] index of subarray */
  ) const
  {
    static_assert(N > 2,
                  "sub_array(i,j) can only be called for N>2");
    if (L == layout::right)
      return _data + i*_offset[0] + j*_offset[1];
    else
      return _data + i*_offset[N-1] + j*_offset[N-2];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on three indexes.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T*
  sub_array
  (
    cs_lnum_t i, /*!<[in] index of subarray */
    cs_lnum_t j, /*!<[in] index of subarray */
    cs_lnum_t k  /*!<[in] index of subarray */
  ) const
  {
    static_assert(N > 3,
                  "sub_array(i,j,k) can only be called for N>3");
    if (L == layout::right)
      return _data + i*_offset[0] + j*_offset[1] + k*_offset[2];
    else
      return _data + i*_offset[N-1] + j*_offset[N-2] + k*_offset[N-3];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter to data array
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T*
  data()
  {
    return _data;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Const getter to data array
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T*
  data() const
  {
    return _data;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator[]
  (
    cs_lnum_t i /*!<[in] Index of value to get */
  )
  {
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator[]
  (
    cs_lnum_t i /*!<[in] Index of value to get */
  ) const
  {
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator()
  (
    cs_lnum_t i /*!<[in] Index of value to get */
  )
  {
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator()
  (
    cs_lnum_t i /*!<[in] Index of value to get */
  ) const
  {
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the (i,j)-th value couple.
   *
   * \returns reference to the (i,j)-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator()
  (
    cs_lnum_t i, /*!<[in] Index along first dimension */
    cs_lnum_t j  /*!<[in] Index along second dimension */
  )
  {
    static_assert(N == 2,
                  "Operator (i,j) can only be called for cs::array<T,2>");

    if (L == layout::right)
      return _data[i*_offset[0] + j];
    else if (L == layout::left)
      return _data[i + j*_offset[1]];
    else
      return _data[i*_offset[0] + j*_offset[1]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the (i,j)-th value couple.
   *
   * \returns const reference to the (i,j)-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator()
  (
    cs_lnum_t i, /*!<[in] Index along first dimension */
    cs_lnum_t j  /*!<[in] Index along second dimension */
  ) const
  {
    static_assert(N == 2,
                  "Operator (i,j) can only be called for cs::array<T,2>");

    if (L == layout::right)
      return _data[i*_offset[0] + j];
    else if (L == layout::left)
      return _data[i + j*_offset[1]];
    else
      return _data[i*_offset[0] + j*_offset[1]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the (i,j,k)-th value tuple.
   *
   * \returns const reference to the (i,j,k)-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator()
  (
    cs_lnum_t i, /*!<[in] Index along first dimension */
    cs_lnum_t j, /*!<[in] Index along second dimension */
    cs_lnum_t k  /*!<[in] Index along third dimension */
  )
  {
    static_assert(N == 3,
                  "Operator (i,j,k) can only be called for cs::array<T,3>");

    if (L == layout::right)
      return _data[i*_offset[0] + j*_offset[1] + k];
    else if (L == layout::left)
      return _data[i + j*_offset[1] + k*_offset[2]];
    else
      return _data[i*_offset[0] + j*_offset[1] + k*_offset[2]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the (i,j,k)-th value tuple.
   *
   * \returns const reference to the (i,j,k)-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  T& operator()
  (
    cs_lnum_t i, /*!<[in] Index along first dimension */
    cs_lnum_t j, /*!<[in] Index along second dimension */
    cs_lnum_t k  /*!<[in] Index along third dimension */
  ) const
  {
    static_assert(N == 3,
                  "Operator (i,j,k) can only be called for cs::array<T,3>");

    if (L == layout::right)
      return _data[i*_offset[0] + j*_offset[1] + k];
    else if (L == layout::left)
      return _data[i + j*_offset[1] + k*_offset[2]];
    else
      return _data[i*_offset[0] + j*_offset[1] + k*_offset[2]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for size per dimension (extent).
   *
   * \returns value of size for given dimension of data array (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  cs_lnum_t
  extent
  (
    int i /*!<[in] dimension id */
  )
  {
    return _extent[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for offset per dimension.
   *
   * \returns value of offset- given dimension (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  cs_lnum_t
  offset
  (
    int i /*!<[in] dimension id */
  )
  {
    return _offset[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for total size.
   *
   * \returns value for total size of data array (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  cs_lnum_t
  size()
  {
    return _size;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for owner status.
   *
   * \returns True if owner, false otherwise (bool)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  bool
  owner()
  {
    return _owner;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for allocation mode.
   *
   * \returns memory allocation mode (cs_alloc_mode_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  cs_alloc_mode_t
  mode()
  {
    return _mode;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Copy data from raw pointer, we suppose that data size is same
   *        as that of the array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    T               *data,       /*!<[in] Pointer to copy */
    const cs_lnum_t  n_vals = -1 /*!<[in] Number of values to copy.
                                          If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    cs_dispatch_context ctx;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = data[e_id];
    });

    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Copy data from another array, we suppose that data size is same
   *        as that of the array. An assert test the sizes in debug.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    array&          other,      /*!<[in] Reference to another array object */
    const cs_lnum_t n_vals = -1 /*!<[in] Number of values to copy.
                                         If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);
    assert(other._size == _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    cs_dispatch_context ctx;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = other._data[e_id];
    });

    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Copy data from a span, we suppose that data size is same
   *        as that of the array. An assert test the sizes in debug.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    mdspan<T, N, L>& span,       /*!<[in] Reference to a span object */
    const cs_lnum_t  n_vals = -1 /*!<[in] Number of values to copy.
                                          If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);
    assert(span._size == _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    cs_dispatch_context ctx;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = span._data[e_id];
    });

    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Copy data from raw pointer, we suppose that data size is same
   *        as that of the array. A dispatch_context is provided, hence
   *        no implicit synchronization which should be done by the caller.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    cs_dispatch_context &ctx,        /*!<[in] Reference to dispatch context  */
    T                   *data,       /*!<[in] Pointer to copy */
    const cs_lnum_t      n_vals = -1 /*!<[in] Number of values to copy.
                                         If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);
    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = data[e_id];
    });
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Copy data from raw pointer, we suppose that data size is same
   *        as that of the array. A dispatch_context is provided, hence
   *        no implicit synchronization which should be done by the caller.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    cs_dispatch_context &ctx,        /*!<[in] Reference to dispatch context  */
    array               &other,      /*!<[in] Reference to another array object */
    const cs_lnum_t      n_vals = -1 /*!<[in] Number of values to copy.
                                              If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);
    assert(other.size() == _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = other._data[e_id];
    });
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Copy data from a span, we suppose that data size is same
   *        as that of the array. A dispatch_context is provided, hence
   *        no implicit synchronization which should be done by the caller.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    cs_dispatch_context &ctx,        /*!<[in] Reference to dispatch context  */
    mdspan<T,N,L>       &span,       /*!<[in] Reference to a span object */
    const cs_lnum_t      n_vals = -1 /*!<[in] Number of values to copy.
                                         If -1, default, we use array size */
  )
  {
    assert(n_vals <= _size);
    assert(span._size == _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = span._data[e_id];
    });
  }

private:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set size of array based on global size.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  set_size_
  (
    const cs_lnum_t size /*!<[in] Global size to set for the array */
  )
  {
    _size = size;
    for (int i = 0; i < N; i++) {
      _extent[i] = 0;
      _offset[i] = 0;
    }
    if (N > 0) {
      _extent[0] = _size;
      _offset[0] = 1;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set size of array based on dimensions.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  set_size_
  (
    const cs_lnum_t dims[] /*!<[in] Array of dimensions' sizes */
  )
  {
    _size = (N > 0) ? 1 : 0;
    for (int i = 0; i < N; i++) {
      _extent[i] = dims[i];
      _size *= dims[i];
    }

    /* Compute offset values for getters */

    if (L == layout::right) {
      /* Version for Layout right */
      for (int i = 0; i < N-1; i++) {
        _offset[i] = 1;
        for (int j = i + 1; j < N; j++)
          _offset[i] *= dims[j];
      }
      _offset[N-1] = 1;
    }
    else if (L == layout::left) {
      for (int i = N-1; i >= 1; i--) {
        _offset[i] = 1;
        for (int j = i - 1; j >= 0; j--)
          _offset[i] *= dims[j];
      }
      _offset[0] = 1;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Private allocator
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  allocate_
  (
    const char *file_name,  /*!<[in] Caller file (for log) */
    const int   line_number /*!<[in] Caller line (for log) */
  )
  {
    const char *_ptr_name = "cs::array._data";
    _data = static_cast<T *>(cs_mem_malloc_hd(_mode,
                                              _size,
                                              sizeof(T),
                                              _ptr_name,
                                              file_name,
                                              line_number));
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Private re-allocator
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  reallocate_
  (
    const char *file_name,  /*!<[in] Caller file (for log) */
    const int   line_number /*!<[in] Caller line (for log) */
  )
  {
    /* If not owner no-op */
    if (_owner) {
      const char *_ptr_name = "cs::array._data";
      _data = static_cast<T *>(cs_mem_realloc_hd(_data,
                                                 _mode,
                                                 _size,
                                                 sizeof(T),
                                                 _ptr_name,
                                                 file_name,
                                                 line_number));
    }
  };

  /*===========================================================================
   * Private members
   *==========================================================================*/

  cs_lnum_t       _extent[N];
  cs_lnum_t       _offset[N];
  cs_lnum_t       _size;
  bool            _owner;
  T*              _data;
  cs_alloc_mode_t _mode;

};

} /* namespace cs */

/*--------------------------------------------------------------------------*/
/* Usefull aliases without namespace */
/*--------------------------------------------------------------------------*/

template<class T>
using cs_array = cs::array<T,1,cs::layout::right>;

template<class T, cs::layout L = cs::layout::right>
using cs_array_2d = cs::array<T,2,L>;

template<class T, cs::layout L = cs::layout::right>
using cs_array_3d = cs::array<T,3,L>;

template<class T, int N, cs::layout L = cs::layout::right>
using cs_mdarray = cs::array<T,N,L>;

template<class T, int N>
using cs_mdarray_r = cs::array<T,N,cs::layout::right>;

template<class T, int N>
using cs_mdarray_l = cs::array<T,N,cs::layout::left>;

#endif /* __cplusplus */

#endif /* __CS_ARRAY_H__ */
