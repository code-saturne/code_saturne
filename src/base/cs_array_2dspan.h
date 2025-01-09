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
class cs_array_2dspan {
public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method.
   */
  /*--------------------------------------------------------------------------*/

  cs_array_2dspan
  (
    cs_lnum_t dim1, /*!<[in] First dimension size */
    cs_lnum_t dim2  /*!<[in] Second dimension size */
  )
  :
    _dim1(dim1),
    _dim2(dim2),
    _size(dim1*dim2),
    _is_owner(true),
    _mode(CS_ALLOC_HOST)
  {
    /* Sanity check for both dimensions */
    assert(_dim1 > 0 && _dim2 > 0);

    allocate_arrays();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method with allocation method provided.
   */
  /*--------------------------------------------------------------------------*/

  cs_array_2dspan
  (
    cs_lnum_t       dim1,      /*!<[in] First dimension size */
    cs_lnum_t       dim2,      /*!<[in] Second dimension size */
    cs_alloc_mode_t alloc_mode /*!<[in] Memory allocation mode */
  )
  :
    _dim1(dim1),
    _dim2(dim2),
    _size(dim1*dim2),
    _is_owner(true),
    _mode(alloc_mode)
  {
    /* Sanity check for both dimensions */
    assert(_dim1 > 0 && _dim2 > 0);

    allocate_arrays();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method for non owner version
   */
  /*--------------------------------------------------------------------------*/

  cs_array_2dspan
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
    _size(dim1*dim2),
    _is_owner(false),
    _mode(alloc_mode)
  {
    /* Sanity check for both dimensions */
    assert(_dim1 > 0 && _dim2 > 0);

    _full_array = data_array;

    allocate_arrays();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  ~cs_array_2dspan()
  {
    if (_is_owner) {
      CS_FREE(_sub_array);
      CS_FREE(_full_array);
    }
    else {
      _sub_array = nullptr;
      _full_array = nullptr;
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
    return _sub_array[i];
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
    return _sub_array[i];
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
  allocate_arrays()
  {
    /* Initialize total size of data array and allocate it if owner */
    if (_is_owner) {
      CS_MALLOC_HD(_full_array, _size, T, _mode);
    }

    /* Allocate and initialize sub array pointer */
    CS_MALLOC_HD(_sub_array, _dim1, T*, _mode);
    for (cs_lnum_t i = 0; i < _dim1; i++)
      _sub_array[i] = _full_array + i * _dim2;
  };

  /*--------------------------------------------------------------------------*/
  /* Private members */
  /*--------------------------------------------------------------------------*/

  T*              _full_array; /*!< Full data array */
  T**             _sub_array;  /*!< Array of sub-arrays pointers */
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
