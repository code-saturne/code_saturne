/*============================================================================
 * Mathematical base functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_math.c
           Mathematical base functions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_inv_cramer
 *----------------------------------------------------------------------------*/

void CS_PROCF (symmetric_matrix_inverse, SYMMETRIC_MATRIX_INVERSE)
(
  const cs_real_6_t s,
  cs_real_6_t       sout
)
{
  cs_math_sym_33_inv_cramer(s,
                            sout);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_product
 *----------------------------------------------------------------------------*/

void CS_PROCF (symmetric_matrix_product, SYMMETRIC_MATRIX_PRODUCT)
(
  const cs_real_6_t s1,
  const cs_real_6_t s2,
  cs_real_6_t       sout
)
{
  cs_math_sym_33_product(s1,
                         s2,
                         sout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
