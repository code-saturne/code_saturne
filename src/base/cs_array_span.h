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
  template <class T, size_t N>
  class span {
    public:

      CS_F_HOST_DEVICE
      span() :
        _dim({0}),
        _size(0),
        _is_owner(true),
        _data(nullptr),
        _mode(cs_alloc_mode)

      CS_F_HOST_DEVICE
      template<cs_lnum_t ...Args>
      span(Args ...args):
        _is_owner(true),
        _full_array(nullptr),
        _mode(cs_alloc_mode)
      {
        static_assert(sizeof...(Args) == N);
        cs_lnum_t dims[] = {args ...};
        for (int i = 0; i < N; i++)
          _dim[i] = dims[i];

        _size = _dim[0];
        for (int i = 1; i < N; i++)
          _size *= _dim[i];

        _data = static_cast<T *>(cs_mem_malloc_hd(_mode,
                                                  _size,
                                                  sizeof(T),
                                                  "_full_array",
                                                  __FILE__,
                                                  __LINE__));
      }


    private:
      cs_lnum_t       _dim[N];
      cs_lnum_t       _size;
      bool            _is_owner;
      T*              _data;
      cs_alloc_mode_t _mode;
  }
}
}

/*----------------------------------------------------------------------------*/

#endif // __cplusplus

/*----------------------------------------------------------------------------*/

#endif // __CS_ARRAY_SPAN_H__
