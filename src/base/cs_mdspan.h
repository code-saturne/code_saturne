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

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method leading to "empty container".
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan():
    _extent({0}),
    _offset({0}),
    _size(0),
    _data(nullptr)
  {}

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using only global size (works for N=1).
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan
  (
    cs_lnum_t  size, /*!<[in] total size of data (1D) */
    T         *data  /*!<[in] data pointer (raw) */
  )
  {
    set_size_(size);
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
    const cs_lnum_t(&dims)[N],  /*!<[in] array of sizes along dimensions */
    T                     *data /*!<[in] data pointer (raw) */
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
                  "Operator (i,j) can only be called for cs::mdspan<T,2>");

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
                  "Operator (i,j) can only be called for cs::mdspan<T,2>");

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
                  "Operator (i,j) can only be called for cs::mdspan<T,3>");

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
                  "Operator (i,j) can only be called for cs::mdspan<T,3>");

    if (L == layout::right)
      return _data[i*_offset[0] + j*_offset[1] + k];
    else if (L == layout::left)
      return _data[i + j*_offset[1] + k*_offset[2]];
    else
      return _data[i*_offset[0] + j*_offset[1] + k*_offset[2]];
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
   * \brief Getter for a subspan based on first dimension
   *
   * \return pointer
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  mdspan<T, N-1, L>
  sub_span
  (
    const cs_lnum_t idx /*!<[in] index of sub span */
  )
  {
    static_assert(N > 1, "Subspan method can only be called for dimension > 1");
    cs_lnum_t dims[N-1];

    cs_lnum_t offset = 0;

    if (L == layout::right) {
      for (int i = 0; i < N-1; i++)
        dims[i] = _extent[i+1];

      offset = _offset[0];
    }
    else {
      for (int i = 0; i < N-1; i++)
        dims[i] = _extent[i];

      offset = _offset[N-1];
    }

    return mdspan<T,N-1,L>(dims, _data + idx*offset);
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
    const cs_lnum_t i /*!<[in] index of sub span */
  )
  {
    if (L == layout::right)
      return _data + i*_offset[0];
    else
      return _data + i*_offset[N-1];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set all values of the data array to a constant value.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void set_to_val
  (
    T val /*!<[in] Value to set to entire data array. */
  )
  {
    cs_dispatch_context ctx;

    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = val;
    });

    ctx.wait();
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
    T *data
  )
  {
    cs_dispatch_context ctx;

    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
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
    mdspan& other
  )
  {
    assert(other._size == _size);

    cs_dispatch_context ctx;

    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
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
    cs_dispatch_context &ctx,
    T                   *data
  )
  {
    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
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
    cs_dispatch_context &ctx,
    mdspan              &other
  )
  {
    assert(other.size() == _size);

    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = other._data[e_id];
    });
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
    cs_dispatch_context &ctx, /*!< Reference to dispatch context */
    T                    val  /*!<[in] Value to set to entire data array. */
  )
  {
    /* No wait here since context is passed as argument */
    ctx.parallel_for(_size, [=] CS_F_HOST_DEVICE (cs_lnum_t e_id) {
        _data[e_id] = val;
    });
  }

private:

  /*===========================================================================
   * Private methods
   *==========================================================================*/

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set size of array based on global size.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  set_size_
  (
    const cs_lnum_t size
  )
  {
    _size = size;
    for (int i = 0; i < N; i++) {
      _extent[i] = 0;
      _offset[i] = 0;
    }

    _extent[0] = _size;
    _offset[0] = 1;
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

  /*===========================================================================
   * Private members
   *==========================================================================*/

  cs_lnum_t  _extent[N];
  cs_lnum_t  _offset[N];
  cs_lnum_t  _size;
  T*         _data;
};

template<class T, layout L = layout::right>
using span = mdspan<T, 1, L>;

template<class T, int N>
using mdspan_r = mdspan<T, N, layout::right>;

template<class T, int N>
using mdspan_l = mdspan<T, N, layout::left>;

} /* namespace cs */

#endif /* __cplusplus */

/*--------------------------------------------------------------------------*/

#endif /* __CS_MDSPAN_H__ */
