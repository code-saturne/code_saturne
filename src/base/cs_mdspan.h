#ifndef __CS_MDSPAN_H__
#define __CS_MDSPAN_H__

/*============================================================================
 * mdspan class and handling utilities.
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

#include "base/cs_dispatch.h"

/*----------------------------------------------------------------------------*/

#if defined(__cplusplus)

namespace cs {

/*! Enum with layout options */
enum class
layout {
  right,   /*!< Layout right */
  left,    /*!< Layout left */
  unknown  /*!< Layout unknown */
};

/*----------------------------------------------------------------------------*/
/*!
 * \class Templated mdspan class. (non owner)
 *
 * @tparam T : data type
 * @tparam N : number of dimensions (int)
 * @tparam L : memory layout (cs::layout)
 */
/*----------------------------------------------------------------------------*/

template<class T, int N, layout L = layout::right>
class mdspan {

/*=============================================================================
 * Public methods
 *============================================================================*/

public:

  /* Make array class friend */
  template<class _T_, int _N_, layout _L_> friend class array;

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method leading to "empty container".
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan():
    _extent{0},
    _offset{0},
    _size(0),
    _data(nullptr)
  {}

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using only global size (works for N=1).
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  mdspan
  (
    T       *data,   /*!<[in] data pointer (raw) */
    Args...  indices /*!<[in] total size of data (1D) */
  )
  {
    check_operator_args_(indices...);
    set_size_(indices...);
    _data = data;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using dimensions.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan
  (
    T                     *data, /*!<[in] data pointer (raw) */
    const cs_lnum_t(&dims)[N]    /*!<[in] array of sizes along dimensions */
  )
  {
    set_size_(dims);
    _data = data;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using copy. May be a shallow copy.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan
  (
    mdspan& other
  )
  {
    set_size_(other._extent);
    _data = other._data;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Move constructor.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan
  (
    mdspan&& other /*!<[in] reference to other instance */
  )
  : mdspan()
  {
    swap(*this, other);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Desstructor method using dimensions.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  ~mdspan()
  {
    _data = nullptr;
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
    mdspan& first, /*!<[in] reference to first instance to swap */
    mdspan& second /*!<[in] reference to second instance to swap */
  )
  {
    using std::swap;

    swap(first._extent, second._extent);
    swap(first._offset, second._offset);
    swap(first._size, second._size);
    swap(first._data, second._data);
  }

  /*===========================================================================
   * Operators
   *==========================================================================*/

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Assignment operator.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan& operator=(mdspan other)
  {
    swap(*this, other);

    return *this;
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

  template<typename Id1>
  CS_F_HOST_DEVICE
  inline
  std::enable_if_t<cs::always_true<Id1>::value && N==1, T&>
  operator()
  (
    Id1 i /*!<[in] Index of value to get */
  )
  {
    check_operator_args_(i);
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  template<typename Id1>
  CS_F_HOST_DEVICE
  inline
  std::enable_if_t<cs::always_true<Id1>::value && N==1, T&>
  operator()
  (
    Id1 i /*!<[in] Index of value to get */
  ) const
  {
    check_operator_args_(i);
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the (i,j)-th value couple.
   *
   * \returns reference to the (i,j)-th value
   */
  /*--------------------------------------------------------------------------*/

  template<typename Id1, typename Id2>
  CS_F_HOST_DEVICE
  inline
  std::enable_if_t<cs::always_true<Id1, Id2>::value && N==2, T&>
  operator()
  (
    cs_lnum_t i, /*!<[in] Index along first dimension */
    cs_lnum_t j  /*!<[in] Index along second dimension */
  )
  {
    check_operator_args_(i,j);

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

  template<typename Id1, typename Id2>
  CS_F_HOST_DEVICE
  inline
  std::enable_if_t<cs::always_true<Id1, Id2>::value && N==2, T&>
  operator()
  (
    Id1 i, /*!<[in] Index along first dimension */
    Id2 j  /*!<[in] Index along second dimension */
  ) const
  {
    check_operator_args_(i,j);

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

  template<typename... Args>
  CS_F_HOST_DEVICE
  inline
  std::enable_if_t<cs::always_true<Args...>::value && (N!=2) && (N!=1), T&>
  operator()
  (
    Args... indices /*!<[in] Input arguments (parameter pack) */
  )
  {
    check_operator_args_(indices...);

    return _data[data_offset_(indices...)];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the (i,j,k)-th value tuple.
   *
   * \returns const reference to the (i,j,k)-th value
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  inline
  std::enable_if_t<cs::always_true<Args...>::value && (N!=2) && (N!=1), T&>
  operator()
  (
    Args... indices /*!<[in] Input arguments (parameter pack) */
  ) const
  {
    check_operator_args_(indices...);

    return _data[data_offset_(indices...)];
  }

  /*===========================================================================
   * Getters
   *==========================================================================*/

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter for total size.
   *
   * \return total size.
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
   * \brief Getter for extent along a given dimension.
   *
   * \return extent.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline
  cs_lnum_t
  extent
  (
    int i
  )
  {
    return _extent[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter for data raw pointer
   *
   * \return raw pointer
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
   * \brief Getter for data raw pointer with recast (casts T* to U*)
   *
   * @tparam U : data type used to cast T* into U*
   *
   * \return raw pointer (U*)
   */
  /*--------------------------------------------------------------------------*/

  template<typename U>
  CS_F_HOST_DEVICE
  U*
  data()
  {
    return reinterpret_cast<U*>(_data);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter for a subspan based on first dimension
   *
   * \return pointer
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  mdspan<T, N-sizeof...(Args), L>
  sub_view
  (
    Args... indices /*!<[in] Input arguments (parameter pack) */
  )
  {
    check_sub_function_args_(indices...);

    constexpr int n_idx = sizeof...(Args);

    cs_lnum_t dims[N - n_idx];

    if (L == layout::right) {
      for (int i = 0; i < N - n_idx; i++)
        dims[i] = _extent[i+n_idx];
    }
    else if (L == layout::left) {
      for (int i = 0; i < N - n_idx; i++)
        dims[i] = _extent[i];
    }

    return mdspan<T,N-n_idx,L>(_data + contiguous_data_offset_(indices...), dims);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get sub array based on index.
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  T*
  sub_array
  (
    Args... indices /*!<[in] Input arguments (parameter pack) */
  )
  {
    check_sub_function_args_(indices...);

    return _data + contiguous_data_offset_(indices...);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set all values of the data array to a constant value.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_to_val
  (
    T               val,        /*!<[in] Value to set to entire data array. */
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
   * \brief Set subset of values of the data array to a constant value.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_to_val_on_subset
  (
    T               val,      /*!<[in] Value to set to entire data array. */
    const cs_lnum_t n_elts,   /*!<[in] Number of values to set. */
    const cs_lnum_t elt_ids[] /*!<[in] list of ids in the subset or null (size:n_elts) */
  )
  {
    assert(n_elts <= _size && n_elts >= 0);

    if (n_elts < 1)
      return;

    cs_dispatch_context ctx;

    if (elt_ids == nullptr)
      set_to_val(ctx, val, n_elts);
    else {
      ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
          _data[elt_ids[e_id]] = val;
      });
    }

    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set a subset of values of the data array to a constant value while providing
   *        a dispatch context. It is up to the call to synchronize the context
   *        after this call.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_to_val_on_subset
  (
    cs_dispatch_context &ctx,      /*!< Reference to dispatch context */
    T                    val,      /*!<[in] Value to set to entire data array. */
    const cs_lnum_t      n_elts,   /*!<[in] Number of values to set. */
    const cs_lnum_t      elt_ids[] /*!<[in] list of ids in the subset or null (size:n_elts) */
  )
  {
    assert(n_elts <= _size && n_elts >= 0);

    if (n_elts < 1)
      return;

    /* No wait here since context is passed as argument */
    if (elt_ids == nullptr)
      set_to_val(ctx, val, n_elts);
    else {
      ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
          _data[elt_ids[e_id]] = val;
      });
    }
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
   * \brief Copy data from another mdspan, we suppose that data size is same
   *        as that of the array. An assert test the sizes in debug.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    mdspan&         other,      /*!<[in] Reference to another mdspan object */
    const cs_lnum_t n_vals = -1 /*!<[in] Number of values to copy.
                                         If -1, default, we use array size */
  )
  {
    assert(other._size == _size);

    assert(n_vals <= _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    cs_dispatch_context ctx;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = other._data[e_id];
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
   * \brief Copy data from another mdspan, we suppose that data size is same
   *        as that of the array. An assert test the sizes in debug.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  copy_data
  (
    cs_dispatch_context &ctx,        /*!<[in] Reference to dispatch context  */
    mdspan              &other,      /*!<[in] Reference to another mdspan object */
    const cs_lnum_t      n_vals = -1 /*!<[in] Number of values to copy.
                                              If -1, default, we use array size */
  )
  {
    assert(other.size() == _size);

    assert(n_vals <= _size);

    const cs_lnum_t loop_size = (n_vals == -1) ? _size : n_vals;

    ctx.parallel_for(loop_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = other._data[e_id];
    });
  }

private:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Helper function to static check operator input arguments.
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  static inline
  void
  check_operator_args_
  (
    Args... /*!<[in] Input arguments (parameter pack) */
  )
  {
    static_assert(sizeof...(Args) == N, "Wrong number of arguments");
    static_assert(sizeof...(Args) != 0, "No input arguments provided...");
    static_assert(cs::are_integral<Args...>::value, "Non integral input arguments.");
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Helper function to static check sub-function input arguments.
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  static inline
  void
  check_sub_function_args_
  (
    Args... /*!<[in] Input arguments (parameter pack) */
  )
  {
    static_assert(sizeof...(Args) < N);
    static_assert(sizeof...(Args) != 0, "No input arguments provided...");
    static_assert(cs::are_integral<Args...>::value);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Helper function to compute value offset.
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  inline
  cs_lnum_t
  data_offset_
  (
    Args... indices /*!<[in] Input arguments (parameter pack) */
  )
  {
    static_assert(sizeof...(Args) <= N && sizeof...(Args) > 0);

    constexpr int n_idx = sizeof...(Args);

    cs_lnum_t _indices[n_idx] = {indices...};

    cs_lnum_t retval = 0;
    for (int i = 0; i < n_idx; i++)
      retval +=_indices[i] * _offset[i];

    return retval;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Helper function to compute value offset.
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  inline
  cs_lnum_t
  contiguous_data_offset_
  (
    Args... indices /*!<[in] Input arguments (parameter pack) */
  )
  {
    static_assert(sizeof...(Args) <= N && sizeof...(Args) > 0);

    constexpr int n_idx = sizeof...(Args);

    cs_lnum_t _indices[n_idx] = {indices...};

    cs_lnum_t retval = 0;
    if (L == layout::right) {
      for (int i = 0; i < n_idx; i++)
        retval +=_indices[i] * _offset[i];
    }
    else if (L == layout::left) {
      for (int i = 0; i < n_idx; i++)
        retval +=_indices[i] * _offset[N-1-i];
    }

    return retval;
  }

  /*===========================================================================
   * Private methods
   *==========================================================================*/

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set size of array based on dimensions.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  set_size_
  (
    const cs_lnum_t(&dims)[N] /*!<[in] array of sizes along dimensions */
  )
  {
    _size = (N > 0) ? 1 : 0;
    for (int i = 0; i < N; i++) {
      _extent[i] = dims[i];
      _size *= dims[i];
      _offset[i] = 1;
    }

    /* Compute offset values for getters */

    if (L == layout::right) {
      /* Version for Layout right */
      for (int i = 0; i < N-1; i++) {
        for (int j = i + 1; j < N; j++)
          _offset[i] *= dims[j];
      }
    }
    else if (L == layout::left) {
      for (int i = N-1; i >= 1; i--) {
        for (int j = i - 1; j >= 0; j--)
          _offset[i] *= dims[j];
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set size of N indices (variadic template)
   */
  /*--------------------------------------------------------------------------*/

  template<typename... Args>
  CS_F_HOST_DEVICE
  void
  set_size_
  (
    Args... dims /*!<[in] Array of dimensions' sizes */
  )
  {
    check_operator_args_(dims...);

    cs_lnum_t loc_dims[N] = {dims...};

    set_size_(loc_dims);
  }

  /*===========================================================================
   * Private members
   *==========================================================================*/

  cs_lnum_t  _extent[N];
  cs_lnum_t  _offset[N];
  cs_lnum_t  _size;
  T*         _data;
};

} /* namespace cs */

template<class T, cs::layout L = cs::layout::right>
using cs_span = cs::mdspan<T, 1, L>;

template<class T, cs::layout L = cs::layout::right>
using cs_span_2d = cs::mdspan<T, 2, L>;

template<class T, cs::layout L = cs::layout::right>
using cs_span_3d = cs::mdspan<T, 3, L>;

template<class T, int N>
using cs_mdspan_r = cs::mdspan<T, N, cs::layout::right>;

template<class T, int N>
using cs_mdspan_l = cs::mdspan<T, N, cs::layout::left>;

#endif /* __cplusplus */

/*--------------------------------------------------------------------------*/

#endif /* __CS_MDSPAN_H__ */
