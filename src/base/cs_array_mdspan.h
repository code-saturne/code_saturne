#ifndef __CS_ARRAY_SPAN_H__
#define __CS_ARRAY_SPAN_H__

/*============================================================================
 * Templated 2D array object
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_defs.h"

#include "base/cs_array.h"
#include "base/cs_mem.h"

#if defined(__cplusplus)

#include <array>

/*=============================================================================
 * BUILTIN TEST
 *============================================================================*/

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

/*=============================================================================
 * Public C++ class template
 *============================================================================*/

namespace cs {

namespace array {


/*----------------------------------------------------------------------------*/
/*!
 * \class Templated data array class.
 */
/*----------------------------------------------------------------------------*/

template<class T>
class data_array {

public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method leading to "empty container".
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  data_array():
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
  data_array
  (
    cs_lnum_t     size,
    cs_alloc_mode_t alloc_mode = cs_alloc_mode,
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
    _size(size),
    _owner(true),
    _data(nullptr),
    _mode(alloc_mode)
  {
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for non owner version
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  data_array
  (
    cs_lnum_t  size,
    T         *data,
    cs_alloc_mode_t alloc_mode = cs_alloc_mode,
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
    _size(size),
    _owner(false),
    _mode(alloc_mode),
    _data(data)
  {
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using copy. May be a shallow copy.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  data_array
  (
    data_array& other,
    bool        shallow_copy=false,
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
    _size = other._size;
    _mode = other._mode;

    /* If shallow copy new instance is not owner. Otherwise same ownership
     * as original instance since we copy it.
     */
    _owner = (shallow_copy) ? false : other._owner;

    if (_owner) {
      allocate_(file_name, line_number);
      cs_array_copy<T>(_size, other._data, _data);
    }
    else {
      _data = other._data;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Move constructor.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  data_array
  (
    data_array&& other /*!<[in] Original reference to move */
  )
  : data_array()
  {
    swap(*this, other);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  ~data_array()
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
    data_array& first,
    data_array& second
  )
  {
    using std::swap;

    /* Swap the different members between the two references. */
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
  data_array& operator=(data_array other)
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

    _size = 0;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Change pointers/size of an existing container.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  point_to
  (
    data_array& other /*!<[in] Other instance to which we want to point
                               to (shallow copy) */
  )
  {
    clear();
    *this = other.view();
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
    cs_arrays_set_value<T,1>(_size, val, _data);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Initializer method for empty containers.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  empty()
  {
    _size = 0;
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
  void resize
  (
    cs_lnum_t       new_size,         /*!<[in] New size */
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
    assert(new_size >= 0);

    /* If same dimensions, nothing to do ... */
    if (new_size == _size)
      return;

    if (_owner) {
      clear();
      _size = new_size;
      allocate_(file_name, line_number);
    }

  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void resize
  (
    cs_lnum_t       new_size,     /*!<[in] New size */
    bool            copy_data,    /*!<[in] Copy data from old pointer to new
                                           array. Default is false. */
    cs_lnum_t       size_to_keep, /*!<[in] Size of data to keep */
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
    assert(new_size >= 0);
    assert(size_to_keep <= new_size && size_to_keep <= _size);

    /* If same dimensions, nothing to do ... */
    if (new_size == _size)
      return;

    if (copy_data && !(size_to_keep <= new_size && size_to_keep <= _size))
      bft_error(__FILE__, __LINE__, 0,
                "%s: Data cannot be saved when new size is smaller than size to keep.\n",
                __func__);

    if (_owner) {
      if (copy_data) {
        /* If new size is larger, realloc is sufficient. */
        _size = new_size;
        reallocate_(file_name, line_number);
      }
      else {
        resize(new_size, file_name, line_number);
      }
    }

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
      if (_size != 0)
        clear();

      _mode = mode;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter to data array
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T*
  vals()
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
  vals() const
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
  T& operator[]
  (
    int i
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
  const T& operator[]
  (
    int i
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
  T& operator()
  (
    int i
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
  const T& operator()
  (
    int i
  ) const
  {
    return _data[i];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for total size.
   *
   * \returns value for total size of data array (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  cs_lnum_t size()
  {
    return _size;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for allocation mode.
   *
   * \returns memory allocation mode (cs_alloc_mode_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  cs_alloc_mode_t mode()
  {
    return _mode;
  }

private:

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
    const char *_ptr_name = "[data_array]._data";
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
      const char *_ptr_name = "[data_array]._data";
      _data = static_cast<T *>(cs_mem_realloc_hd(_data,
                                                 _mode,
                                                 _size,
                                                 sizeof(T),
                                                 _ptr_name,
                                                 file_name,
                                                 line_number));
    }
  };

  /*--------------------------------------------------------------------------*/
  /* Private members */
  /*--------------------------------------------------------------------------*/

  cs_lnum_t       _size;
  bool            _owner;
  T*              _data;
  cs_alloc_mode_t _mode;

};

/*----------------------------------------------------------------------------*/
/*!
 * \class Templated data array class.
 */
/*----------------------------------------------------------------------------*/

template <class T, int N>
class span {
public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method leading to "empty container".
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span() :
    _dim({0}),
    _offset({0}),
    _size(0),
    _owner(true),
    _data(nullptr),
    _mode(cs_alloc_mode)
  {}

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using only dimension.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span
  (
    const cs_lnum_t(&dims)[N],
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
    _owner(true),
    _data(nullptr),
    _mode(cs_alloc_mode)
  {
    set_size_(dims);
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method with specified allocation method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span
  (
    const cs_lnum_t(&dims)[N],
    cs_alloc_mode_t alloc_mode, /*!<[in] Memory allocation mode */
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
    _owner(true),
    _data(nullptr),
    _mode(alloc_mode)
  {
    set_size_(dims);
    allocate_(file_name, line_number);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for non owner version based on raw pointer
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span
  (
    const cs_lnum_t(&dims)[N],
    T*                      data,              /*!<[in] Pointer to data array */
    cs_alloc_mode_t alloc_mode = cs_alloc_mode /*!<[in] Memory allocation mode,
                                                        default is HOST. */
  )
  :
    _owner(false),
    _data(data),
    _mode(alloc_mode)
  {
    set_size_(dims);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for non owner from data_array class
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span
  (
    const cs_lnum_t(&dims)[N],
    data_array<T>& array,
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
    _owner(false)
  {
    // Check that the span is correct
    set_size_(dims);

    if (_size != array.size())
      bft_error(__FILE__, __LINE__, 0,
                _("%s: The dimensions provided for the span lead to size %d "
                  "which is different than size %d of the data array.\n"
                  "Called from file \"%s\" L%d\n"),
                __func__, _size, array.size(), file_name, line_number);

    _mode = array.mode();
    _data = array.vals();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using copy. May be a shallow copy.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span
  (
    span& other,              /*!<[in] Instance to copy */
    bool  shallow_copy=false, /*!<[in] Do a shallow copy or not */
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
    set_size_(other._dim);
    _mode = other._mode;

    /* If shallow copy new instance is not owner. Otherwise same ownership
     * as original instance since we copy it.
     */
    _owner = (shallow_copy) ? false : other._owner;

    if (_owner) {
      allocate_(file_name, line_number);
      cs_array_copy<T>(_size, other._data, _data);
    }
    else {
      _data = other._data;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Move constructor.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span
  (
    span&& other /*!<[in] Original reference to move */
  )
  : span()
  {
    swap(*this, other);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  ~span()
  {
    clear();
  }

  CS_F_HOST_DEVICE
  friend void
  swap
  (
    span& first, /*!<[in,out] First class instance */
    span& second /*!<[in,out] Second class instance */
  )
  {
    using std::swap;
    /* Swap the different members between the two references. */
    swap(first._dim, second._dim);
    swpa(first._offset, second._offset);
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
  span& operator=(span other)
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
    cs_lnum_t dummy_size[N] = {0};
    set_size_(dummy_size);
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
    cs_lnum_t dummy_size[N] = {0};
    set_size_(dummy_size);

    _owner = false;
    _data = nullptr;
    _mode = cs_alloc_mode;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Change values so as to point
   */
  /*--------------------------------------------------------------------------*/


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
    cs_arrays_set_value<T,1>(_size, val, _data);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void resize
  (
    const cs_lnum_t(&dims)[N],
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
    bool same_dim = false;
    for (int i = 0; i < N; i++)
      if (dims[i] == _dim[i])
        same_dim = true;

    if (same_dim)
      return;

    if (_owner) {
      clear();
      set_size_(dims);
      allocate_(file_name, line_number);
    }

  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array while keeping some of the old data.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void resize
  (
    cs_lnum_t d1,
    bool      copy_data,    /*!<[in] Copy data from old pointer to new
                                array. Default is false. */
    cs_lnum_t       size_to_keep, /*!<[in] Size of data to keep */
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
    static_assert(N == 1);

    assert(size_to_keep <= d1 && size_to_keep <= _dim[0]);

    /* If same dimensions, nothing to do ... */
    if (d1 == _dim[0])
      return;

    if (_owner) {
      if (copy_data) {
        set_size_({d1});
        reallocate_(file_name, line_number);
      }
      else {
        resize({d1}, file_name, line_number);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array while keeping some of the old data.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void resize
  (
    cs_lnum_t d1,
    cs_lnum_t d2,
    bool      copy_data,    /*!<[in] Copy data from old pointer to new
                                array. Default is false. */
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
    static_assert(N == 2);

    /* If same dimensions, nothing to do ... */
    if (d1 == _dim[0])
      return;

    if (_owner) {
      if (copy_data) {
        /* Temporary copy */
        span<T,N> tmp(*this, false);

        /* Update this instance sizes */
        set_size_({d1, d2});
        allocate_(file_name, line_number);

        /* Copy what can be copied */
        for (cs_lnum_t i = 0; i < d1 && i < tmp.dim(0); i++) {
          for (cs_lnum_t j = 0; j < d2 && j < tmp.dim(1); j++) {
            _data[i*_offset[0] + j*_offset[1]] = tmp(i,j);
          }
        }
      }
      else {
        resize({d1,d2}, file_name, line_number);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array while keeping some of the old data.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void resize
  (
    cs_lnum_t d1,
    cs_lnum_t d2,
    cs_lnum_t d3,
    bool      copy_data,    /*!<[in] Copy data from old pointer to new
                                array. Default is false. */
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
    static_assert(N == 3);

    /* If same dimensions, nothing to do ... */
    if (d1 == _dim[0] && d2 == _dim[1] && d3 == _dim[2])
      return;

    if (_owner) {
      if (copy_data) {
        /* Temporary copy */
        span<T,N> tmp(*this, false);

        /* Update this instance sizes */
        set_size_({d1, d2, d3});
        allocate_(file_name, line_number);

        /* Copy what can be copied */
        for (cs_lnum_t i = 0; i < d1 && i < tmp.dim(0); i++) {
          for (cs_lnum_t j = 0; j < d2 && j < tmp.dim(1); j++) {
            for (cs_lnum_t k = 0; k < d3 && k < tmp.dim(2); k++) {
              this(i,j,k) = tmp(i,j,k);
            }
          }
        }
      }
      else {
        resize({d1,d2,d3}, file_name, line_number);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get a view (non owner) of the array.
   *
   * \return temporary instance which is a view (non owner) of the data
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  span view()
  {
    /* Use the shallow copy constructor to return a temporary instance */
    return span(*this, true);
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
      /* Changing allocation mode requires freeing the memory ... */
      if (_size != 0)
        clear();

      _mode = mode;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  cs_lnum_t
  size()
  {
    return _size;
  }

  CS_F_HOST_DEVICE
  cs_alloc_mode_t
  mode()
  {
    return _mode;
  }

  CS_F_HOST_DEVICE
  cs_lnum_t
  dim
  (
    int i
  )
  {
    assert(i >= 0 && i < N);
    return _dim[i];
  }

  CS_F_HOST_DEVICE
  cs_lnum_t
  offset
  (
    int i
  )
  {
    assert(i >= 0 && i < N);
    return _offset(i);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T& operator()
  (
    cs_lnum_t i
  )
  {
    static_assert(N == 1);
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
  const T& operator()
  (
    cs_lnum_t i
  ) const
  {
    static_assert(N == 1);
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
  T& operator()
  (
    cs_lnum_t i,
    cs_lnum_t j
  )
  {
    static_assert(N == 2);
    return _data[i*_offset[0] + j*_offset[1]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T& operator()
  (
    cs_lnum_t i,
    cs_lnum_t j
  ) const
  {
    static_assert(N == 2);
    return _data[i*_offset[0] + j*_offset[1]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded () operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T& operator()
  (
    cs_lnum_t i,
    cs_lnum_t j,
    cs_lnum_t k
  )
  {
    static_assert(N == 3);
    return _data[i*_offset[0] + j*_offset[1] + k*_offset[2]];
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith value (val[i]).
   *
   * \returns raw pointer to the i-th value
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T& operator()
  (
    cs_lnum_t i,
    cs_lnum_t j,
    cs_lnum_t k
  ) const
  {
    static_assert(N == 3);
    return _data[i*_offset[0] + j*_offset[1] + k*_offset[2]];
  }

private:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Private set dimensions method
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  set_size_
  (
    const cs_lnum_t dims[]
  )
  {
    _size = (N > 0) ? 1 : 0;
    for (int i = 0; i < N; i++) {
      _dim[i] = dims[i];
      _size *= dims[i];
    }

    /* Compute offset values for getters */

    /* Version for Layout right */
    for (int i = 0; i < N-1; i++) {
      _offset[i] = 1;
      for (int j = i + 1; j < N; j++)
        _offset[i] *= dims[j];
    }
    _offset[N-1] = 1;
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
    const char *ptr_name = "cs::array::span._data";
    _data = static_cast<T *>(cs_mem_malloc_hd(_mode,
                                              _size,
                                              sizeof(T),
                                              ptr_name,
                                              file_name,
                                              line_number));
  };

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
      const char *_ptr_name = "cs::array::span.data";
      _data = static_cast<T *>(cs_mem_realloc_hd(_data,
                                                 _mode,
                                                 _size,
                                                 sizeof(T),
                                                 _ptr_name,
                                                 file_name,
                                                 line_number));
    }
  };

  /*--------------------------------------------------------------------------*/
  /* Private members */
  /*--------------------------------------------------------------------------*/

  cs_lnum_t       _dim[N];
  cs_lnum_t       _offset[N];
  cs_lnum_t       _size;
  bool            _owner;
  T*              _data;
  cs_alloc_mode_t _mode;
};

} /* namespace array */
} /* namespace cs */

/*----------------------------------------------------------------------------*/

#endif // __cplusplus

/*----------------------------------------------------------------------------*/

#endif // __CS_ARRAY_SPAN_H__
