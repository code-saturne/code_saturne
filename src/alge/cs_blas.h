#ifndef __CS_BLAS_H__
#define __CS_BLAS_H__

/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * External library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_ESSL)
#include <essl.h>

#elif defined(HAVE_MKL)
#include <mkl_cblas.h>

#elif defined(HAVE_ACML)
#include <acml.h>

#elif defined(HAVE_ATLAS) || defined(HAVE_CBLAS)
#include <cblas.h>

#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/* Return the dot product of 2 vectors: x.y */

double CS_PROCF(csdot, CSDOT)(const cs_int_t  *n,
                              const double    *x,
                              const double    *y);

/*============================================================================
 *  Public function prototypes or wrapper macros
 *============================================================================*/

/* For the IBM ESSL library, function prototypes are defined in essl.h,
   with legacy blas names <name> mapped to esv<name>;
   for the AMD ACML library, prototypes are defined in acml.h.
   In both cases, purely input array arguments are not defined as const,
   so a cast must be used.
   For the Intel MKL library, function prototypes are defined in mkl_cblas.h,
   with standard legacy C BLAS names */

/*----------------------------------------------------------------------------
 * Constant times a vector plus a vector: y <-- ax + y
 *
 * parameters:
 *   n <-- size of arrays x and y
 *   a <-- multiplier for x
 *   x <-- array of floating-point values
 *   y <-- array of floating-point values
 *----------------------------------------------------------------------------*/

#if defined(HAVE_ACML) || defined(HAVE_ESSL)

#define cs_axpy(_n, _a, _x, _y) \
  daxpy(_n, _a, (double *)_x, 1, (double *)_y, 1)

#elif defined(HAVE_ATLAS) || defined(HAVE_CBLAS) || defined(HAVE_MKL)

#define cs_axpy(_n, _a, _x, _y) \
  cblas_daxpy(_n, _a, _x, 1, _y, 1)

#else

void
cs_axpy(cs_lnum_t      n,
        double         a,
        const double  *x,
        double        *restrict y);

#endif /* BLAS defined */

/*----------------------------------------------------------------------------
 * Return the dot product of 2 vectors: x.y
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n <-- size of arrays x and y
 *   x <-- array of floating-point values
 *   y<-- array of floating-point values
 *
 * returns:
 *   dot product
 *----------------------------------------------------------------------------*/

double
cs_dot(cs_lnum_t      n,
       const double  *x,
       const double  *y);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLAS_H__ */
