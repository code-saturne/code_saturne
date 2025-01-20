#ifndef __CS_ARRAY_2DSPAN_H__
#define __CS_ARRAY_2DSPAN_H__

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
#include "base/cs_mem.h"

#if defined(__cplusplus)

/*=============================================================================
 * Public C++ class template
 *============================================================================*/

namespace cs {

/*----------------------------------------------------------------------------*/
/*!
 * \class Templated 2DSpan array class.
 */
/*----------------------------------------------------------------------------*/

template <class T>
class array_2dspan {
public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method leading to "empty container".
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan() :
    _dim1(0),
    _dim2(0),
    _size(0),
    _is_owner(true),
    _full_array(nullptr),
    _mode(CS_ALLOC_HOST)
  {
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using only dimension.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan
  (
    cs_lnum_t dim1, /*!<[in] First dimension size */
    cs_lnum_t dim2  /*!<[in] Second dimension size */
  )
  :
    _dim1(dim1),
    _dim2(dim2),
    _is_owner(true),
    _mode(CS_ALLOC_HOST)
  {
    /* Sanity check for both dimensions */
    assert(_dim1 > 0 && _dim2 > 0);
    _size = dim1*dim2;

    allocate_();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method with allocation method provided.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan
  (
    cs_lnum_t       dim1,      /*!<[in] First dimension size */
    cs_lnum_t       dim2,      /*!<[in] Second dimension size */
    cs_alloc_mode_t alloc_mode /*!<[in] Memory allocation mode */
  )
  :
    _dim1(dim1),
    _dim2(dim2),
    _is_owner(true),
    _mode(alloc_mode)
  {
    /* Sanity check for both dimensions */
    assert(_dim1 > 0 && _dim2 > 0);
    _size = dim1*dim2;

    allocate_();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for non owner version
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan
  (
    cs_lnum_t       dim1,                      /*!<[in] First dimension size */
    cs_lnum_t       dim2,                      /*!<[in] Second dimension size */
    T*              data_array,                /*!<[in] Pointer to data array */
    cs_alloc_mode_t alloc_mode = CS_ALLOC_HOST /*!<[in] Memory allocation mode,
                                                        default is HOST. */
  )
  :
    _dim1(dim1),
    _dim2(dim2),
    _is_owner(false),
    _mode(alloc_mode)
  {
    /* Sanity check for both dimensions */
    assert(_dim1 > 0 && _dim2 > 0);

    _size = dim1*dim2;

    _full_array = data_array;

    allocate_();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method using copy. May be a shallow copy.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan
  (
   array_2dspan& other,             /*!<[in] Instance to copy */
   bool          shallow_copy=false /*!<[in] Do a shallow copy or not */
  )
  {
    _dim1 = other._dim1;
    _dim2 = other._dim2;
    _size = other._size;
    _mode = other._mode;

    /* If shallow copy new instance is not owner. Otherwise same ownership
     * as original instance since we copy it.
     */
    _is_owner = (shallow_copy) ? false : other._is_owner;

    if (_is_owner) {
      allocate_();
      for (cs_lnum_t e_id = 0; e_id < _size; e_id++)
        _full_array[e_id] = other._full_array[e_id];
    }
    else {
      _full_array = other._full_array;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Move constructor.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan
  (
    array_2dspan&& other /*!<[in] Original reference to move */
  )
    : array_2dspan()
  {
    swap(*this, other);
  }
  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  ~array_2dspan()
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
    array_2dspan& first, /*!<[in,out] First class instance */
    array_2dspan& second /*!<[in,out] Second class instance */
  )
  {
    using std::swap;
    /* Swap the different members between the two references. */
    swap(first._dim1, second._dim1);
    swap(first._dim2, second._dim2);
    swap(first._size, second._size);
    swap(first._is_owner, second._is_owner);
    swap(first._mode, second._mode);
    swap(first._full_array, second._full_array);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Assignment operator.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  array_2dspan& operator=(array_2dspan other)
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
    if (_is_owner) {
      CS_FREE(_full_array);
    }
    else {
      _full_array = nullptr;
    }

    _dim1 = 0;
    _dim2 = 0;
    _size = 0;
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
    _dim1 = _dim2 = _size = 0;
    _is_owner = false;
    _full_array = nullptr;
    _mode = CS_ALLOC_HOST;
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
    array_2dspan& other /*!<[in] Other instance to which we want to point
                                 to (shallow copy) */
  )
  {
    clear();
    *this = array_2dspan(other, true);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Resize data array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void resize
  (
    cs_lnum_t       dim1, /*!<[in] First dimension size */
    cs_lnum_t       dim2  /*!<[in] Second dimension size */
  )
  {
    assert(dim1 > 0 && dim2 > 0);

    /* If same dimensions, nothing to do ... */
    if (dim1 == _dim1 && dim2 == _dim2)
      return;

    if (_is_owner) {
      clear();

      _dim1 = dim1;
      _dim2 = dim2;
      _size = dim1 * dim2;

      allocate_();
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
    cs_lnum_t       dim1,         /*!<[in] First dimension size */
    cs_lnum_t       dim2,         /*!<[in] Second dimension size */
    bool            copy_data,    /*!<[in] Copy data from old pointer to new
                                           array. Default is false. */
    cs_lnum_t       size_to_keep  /*!<[in] Size of data to keep */
  )
  {
    assert(dim1 > 0 && dim2 > 0);
    assert(size_to_keep <= dim2 && size_to_keep <= _dim2);

    /* If same dimensions, nothing to do ... */
    if (dim1 == _dim1 && dim2 == _dim2)
      return;

    if (copy_data && dim1 != _dim1)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Data cannot be saved with new dim1 is different from previous.\n",
                __func__);

    if (_is_owner) {
      if (copy_data) {
        T *tmp = nullptr;
        CS_MALLOC_HD(tmp, _size, T, _mode);
        for (cs_lnum_t e_id = 0; e_id < _size; e_id++)
          tmp[e_id] = _full_array[e_id];

        cs_lnum_t old_dim2 = _dim2;

        resize(dim1, dim2);

        /* We loop on "_dim1" since dim1 = _dim1 */
        for (cs_lnum_t i = 0; i < _dim1; i++) {
          for (cs_lnum_t j = 0; j < size_to_keep; j++)
            _full_array[i*_dim2 + j] = tmp[i*old_dim2 + j];
        }

        CS_FREE(tmp);
      }
      else {
        resize(dim1, dim2);
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
   * \brief Getter function for full array.
   *
   * \returns raw pointer to the full data array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T *vals()
  {
    return _full_array;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Const getter function for full array.
   *
   * \returns const raw pointer to the full data array.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T *vals() const
  {
    return _full_array;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded [] operator to access the ith subarray (val[i][...]).
   *
   * \returns raw pointer to the i-th sub-array
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  T* operator[]
  (
    int i /*!<[in] sub-array index to access */
  )
  {
    return _full_array + i * _dim2;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Overloaded const [] operator to access the ith subarray (val[i][...]).
   *
   * \returns const raw pointer to the i-th sub-array
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  const T* operator[]
  (
    int i /*!<[in] sub-array index to access */
  ) const
  {
    return _full_array + i * _dim2;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for first dimension size.
   *
   * \returns value for first dimension size (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  cs_lnum_t dim1()
  {
    return _dim1;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for second dimension size.
   *
   * \returns value for second dimension size (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  cs_lnum_t dim2()
  {
    return _dim2;
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

private:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Private allocator
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  void
  allocate_()
  {
    /* Initialize total size of data array and allocate it if owner */
    if (_is_owner) {
      CS_MALLOC_HD(_full_array, _size, T, _mode);
    }
  };

  /*--------------------------------------------------------------------------*/
  /* Private members */
  /*--------------------------------------------------------------------------*/

  cs_lnum_t       _dim1;       /*!< First dimension size */
  cs_lnum_t       _dim2;       /*!< Second dimension size */
  cs_lnum_t       _size;       /*!< Total size of data array */
  bool            _is_owner;   /*!< Is this array owner of data */
  T*              _full_array; /*!< Full data array */
  cs_alloc_mode_t _mode;       /*!< Data allocation mode */
};

}

/*----------------------------------------------------------------------------*/

#endif // __cplusplus

/*----------------------------------------------------------------------------*/

#endif // __CS_ARRAY_2DSPAN_H__
