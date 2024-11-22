#ifndef __CS_ARRAY_H__
#define __CS_ARRAY_H__

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

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_base_accel.h"
#if defined (__NVCC__)
#include "cs_array_cuda.h"
#include "cs_base_cuda.h"
#endif

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

#endif /* __CS_ARRAY_H__ */
