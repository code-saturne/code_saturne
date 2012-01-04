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

#elif defined(HAVE_CBLAS)
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

#if    defined(HAVE_CBLAS) || defined(HAVE_ESSL) \
    || defined (HAVE_MKL)  || defined (HAVE_ACML)
#define HAVE_BLAS 1
#endif

/*----------------------------------------------------------------------------
 * Definition of some C BLAS functions depending on external libraries
 *----------------------------------------------------------------------------*/

/*
  - Constant times a vector plus a vector: y <-- ax + y

    void daxpy(int            n,
               double         a,
               const double  *x,
               double        *y);

  - Return the dot product of 2 vectors: x.y

    double ddot(int            n,
                const double  *x,
                const double  *y);

*/

/* For the IBM ESSL library, function prototypes are defined in essl.h,
   with legacy blas names <name> mapped to esv<name>;
   for the AMD ACML library, prototypes are defined in acml.h.
   In both cases, purely input array arguments are not defined as const,
   so a cast must be used. */

#if defined(HAVE_ESSL) || defined(HAVE_ACML)

#define cs_axpy(_n, _a, _x, _y) \
  daxpy(_n, _a, (double *)_x, 1, (double *)_y, 1)

#define cs_dot(_n, _x, _y) \
  ddot(_n, (double *)_x, 1, (double *)_y, 1)

/* For the Intel MKL library, function prototypes are defined in mkl_cblas.h,
   with standard legacy C BLAS names */

#elif defined(HAVE_MKL) || defined(HAVE_CBLAS)

#define cs_axpy(_n, _a, _x, _y) \
  cblas_daxpy(_n, _a, _x, 1, _y, 1)

#define cs_dot(_n, _x, _y) \
  cblas_ddot(_n, _x, 1, _y, 1)

#define HAVE_CBLAS 1

/* Otherwise, if the legacy C BLAS names are not defined, we define
   some double precision legacy BLAS 1 prototypes */

#else

/* Constant times a vector plus a vector: y <-- ax + y */

void cs_axpy(cs_lnum_t      n,
             double         a,
             const double  *x,
             double        *restrict y);

/* Return the dot product of 2 vectors: x.y */

double cs_dot(cs_lnum_t      n,
              const double  *x,
              const double  *y);

#endif /* BLAS type */

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/* Return the dot product of 2 vectors: x.y */

double CS_PROCF(csdot, CSDOT)(const cs_int_t  *n,
                              const double    *x,
                              const double    *y);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLAS_H__ */
