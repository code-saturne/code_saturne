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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_array.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_array.c
        Array handling utilities.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

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
cs_array_copy_real(cs_lnum_t        n_elts,
                   cs_lnum_t        dim,
                   const cs_real_t  src[],
                   cs_real_t        dest[restrict])
{
  const cs_lnum_t _n_elts = dim * n_elts;

# pragma omp parallel for if (_n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < _n_elts; ii++)
    dest[ii] = src[ii];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant value to an array.
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
                        cs_real_t  a[])
{
  const cs_lnum_t _n_elts = dim * n_elts;

# pragma omp parallel for if (_n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < _n_elts; ii++)
    a[ii] = v;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
