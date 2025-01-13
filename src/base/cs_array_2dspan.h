#ifndef __CS_2D_ARRAY_H__
#define __CS_2D_ARRAY_H__

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

#ifdef __cplusplus

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

    if (shallow_copy) {
      _is_owner = false;
      _full_array = other._full_array;
    }
    else {
      _is_owner = true;
      allocate_();
      for (cs_lnum_t e_id = 0; e_id < _size; e_id++)
        _full_array[e_id] = other._full_array[e_id];
    }
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  ~array_2dspan()
  {
    clear();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Clear data (empty container).
   */
  /*--------------------------------------------------------------------------*/

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
   * \brief Resize data array.
   */
  /*--------------------------------------------------------------------------*/

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
                "%s: Data cannot be saved with new dim1 is different from previous.\n");

    cs_lnum_t new_size = dim1 * dim2;

    if (_is_owner) {
      if (copy_data) {
        T *tmp = nullptr;
        CS_MALLOC_HD(tmp, _size, T, _mode);
        for (cs_lnum_t e_id = 0; e_id < _size; e_id++)
          tmp[e_id] = _full_array[e_id];

        CS_FREE(_full_array);
        CS_MALLOC_HD(_full_array, new_size, T, _mode);

        for (cs_lnum_t i = 0; i < _dim1; i++) {
          for (cs_lnum_t j = 0; j < size_to_keep; j++)
            _full_array[i*dim2 + j] = tmp[i*_dim2 + j];
        }
        CS_FREE(tmp);
      }
      else {
        CS_FREE(_full_array);
        CS_MALLOC_HD(_full_array, new_size, T, _mode);
      }
    }

    _dim1 = dim1;
    _dim2 = dim2;
    _size = dim1 * dim2;
  }
  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set memory allocation mode.
   */
  /*--------------------------------------------------------------------------*/

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

  T*              _full_array; /*!< Full data array */
  cs_lnum_t       _dim1;       /*!< First dimension size */
  cs_lnum_t       _dim2;       /*!< Second dimension size */
  cs_lnum_t       _size;       /*!< Total size of data array */
  bool            _is_owner;   /*!< Is this array owner of data */
  cs_alloc_mode_t _mode;       /*!< Data allocation mode */
};

}

/*----------------------------------------------------------------------------*/

#endif // __cplusplus

/*----------------------------------------------------------------------------*/

#endif // __CS_2D_ARRAY_H__
