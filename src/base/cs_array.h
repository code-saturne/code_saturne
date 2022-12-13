#ifndef __CS_ARRAY_H__
#define __CS_ARRAY_H__

/*============================================================================
 * Array handling utilities.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign zero to all elements of an array.
 *
 * \param[in]      size    total number of elements to set to zero
 * \param[in, out] a       array to set
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_array_real_fill_zero(cs_lnum_t  size,
                        cs_real_t  a[])
{
  memset(a, 0, size*sizeof(cs_real_t));
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy an array (ref) into another array (dest) and apply an
 *        inderection at the same time. Array with stride > 1 are assumed to be
 *        interlaced. The indirectino is applied to ref
 *
 * \param[in]      n_elts   number of elements in the array
 * \param[in]      stride   number of values for each element
 * \param[in]      elt_ids  indirection list
 * \param[in]      ref      reference values to copy
 * \param[in, out] dest     array storing values after applying the indirection
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_copy_with_indirection(cs_lnum_t        n_elts,
                                    int              stride,
                                    const cs_lnum_t  elt_ids[],
                                    const cs_real_t  ref[],
                                    cs_real_t        dest[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy real values from an array to another of the same dimensions.
 *
 * \param[in]   n_elts  number of associated elements
 * \param[in]   dim     associated dimension
 * \param[in]   src     source array values (size: n_elts*dim]
 * \param[out]  dest    destination array values (size: n_elts*dim]
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_copy(cs_lnum_t        n_elts,
                   cs_lnum_t        dim,
                   const cs_real_t  src[],
                   cs_real_t        dest[restrict]);

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
                         cs_real_t  a[]);

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
                         cs_real_t        *a);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant vector of size 6 (optimzed way to define a
 *        symmetric tensor) to an array (of stride 6) which is interlaced
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      ref_val  vector to assign
 * \param[in, out] a        array to set
 */
/*----------------------------------------------------------------------------*/

void
cs_array_real_set_symm_tensor(cs_lnum_t         n_elts,
                              const cs_real_t   ref_val[6],
                              cs_real_t        *a);

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
cs_array_real_set_tensor(cs_lnum_t         n_elts,
                         const cs_real_t   ref_tens[3][3],
                         cs_real_t        *a);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant value to an array (deprecated function).
 *
 * \param[in]   n_elts  number of associated elements
 * \param[in]   dim     associated dimension
 * \param[in]   v       value to assign
 * \param[out]  a       array values (size: n_elts*dim]
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
