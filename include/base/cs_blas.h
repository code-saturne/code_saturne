/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_BLAS_H__
#define __CS_BLAS_H__

/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * External library headers
 *----------------------------------------------------------------------------*/

#if defined(_CS_HAVE_ESSL)
#include <essl.h>

#elif defined(_CS_HAVE_MKL)
#include <mkl_cblas.h>

#elif defined(_CS_HAVE_CBLAS)
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

#if    defined(_CS_HAVE_CBLAS) || defined(_CS_HAVE_F77BLAS) \
    || defined(_CS_HAVE_ESSL) || defined (_CS_HAVE_MKL)
#define _CS_HAVE_BLAS 1
#endif

/*----------------------------------------------------------------------------
 * Definition of some C BLAS functions depending on external libraries
 *----------------------------------------------------------------------------*/

/*
  - Sum of the absolute values of a vector

    double dasum(int            n,
                 const double  *x,
                 int            incx);

  - Constant times a vector plus a vector: y <-- ax + y

    void daxpy(int            n,
               double         a,
               const double  *x,
               int            incx,
               double        *y,
               int            incy);

  - Copy a vector x to a vector y: y <-- x

    void dcopy(int            n,
               const double  *x,
               int            incx,
               double        *y,
               int            incy);

  - Return the dot product of 2 vectors: x.y

    double ddot(int            n,
                const double  *x,
                int            incx,
                const double  *y,
                int            incy);

  - Return the euclidean norm of a vector

    double dnrm2(int            n,
                 const double  *x,
                 int            incx);

  - Scales a vector by a constant: x <-- ax

    void dscal(int     n,
               double  a,
               double *x,
               int     incx);

  - Interchange vectors

    void dswap(int      n,
               double  *x,
               int      incx,
               double  *y,
               int      incy);

  - Finds the index of element having max absolute value

    int idamax(int            n,
               const double  *x,
               int            incx);
*/

/* For the IBM ESSL library, function prototypes are defined in essl.h,
   with legacy blas names <name> mapped to esv<name> */

#if defined(_CS_HAVE_ESSL)

#define cblas_dasum  dasum
#define cblas_daxpy  daxpy
#define cblas_dcopy  dcopy
#define cblas_ddot   ddot
#define cblas_dnrm2  dnrm2
#define cblas_dscal  dscal
#define cblas_dswap  dswap

/* For the Intel MKL library, function prototypes are defined in mkl_cblas.h,
   with standard legacy C BLAS names */

#elif defined(_CS_HAVE_MKL)

/* Otherwise, if the legacy C BLAS names are not defined, we define double
   precision legacy BLAS 1 prototypes */

#define _CS_HAVE_CBLAS 1

#elif !defined(_CS_HAVE_CBLAS)

/* Sum of the absolute values of a vector */

double cblas_dasum(cs_int_t       n,
                   const double  *x,
                   cs_int_t       incx);

/* Constant times a vector plus a vector: y <-- ax + y */

void cblas_daxpy(cs_int_t       n,
                 double         a,
                 const double  *x,
                 cs_int_t       incx,
                 double        *restrict y,
                 cs_int_t       incy);

/* Copy a vector x to a vector y: y <-- x */

void cblas_dcopy(cs_int_t       n,
                 const double  *x,
                 cs_int_t       incx,
                 double        *restrict y,
                 cs_int_t       incy);

/* Return the dot product of 2 vectors: x.y */

double cblas_ddot(cs_int_t       n,
                  const double  *x,
                  cs_int_t       incx,
                  const double  *y,
                  cs_int_t       incy);

/* Return the euclidean norm of a vector */

double cblas_dnrm2(cs_int_t       n,
                   const double  *x,
                   cs_int_t       incx);

/* Scales a vector by a constant: x <-- ax */

void cblas_dscal(cs_int_t   n,
                 double     a,
                 double    *x,
                 cs_int_t   incx);

/* Interchange vectors */

void cblas_dswap(cs_int_t   n,
                 double    *restrict x,
                 cs_int_t   incx,
                 double    *restrict y,
                 cs_int_t   incy);

/* Finds the index of element having max absolute value */

cs_int_t cblas_idamax(cs_int_t       n,
                      const double  *x,
                      cs_int_t       incx);

#endif /* !defined(_CS_HAVE_CBLAS) */

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

double CS_PROCF(dasum, DASUM)(const cs_int_t  *n,
                              const double    *x,
                              const cs_int_t  *incx);

/* Constant times a vector plus a vector: y <-- ax + y */

void CS_PROCF(daxpy, DAXPY)(const cs_int_t  *n,
                            const double    *a,
                            const double    *x,
                            const cs_int_t  *incx,
                            double          *restrict y,
                            const cs_int_t  *incy);

/* Copy a vector x to a vector y: y <-- x */

void CS_PROCF(dcopy, DCOPY)(const cs_int_t  *n,
                            const double    *x,
                            const cs_int_t  *incx,
                            double          *restrict y,
                            const cs_int_t  *incy);

/* Return the dot product of 2 vectors: x.y */

double CS_PROCF(ddot, DDOT)(const cs_int_t  *n,
                            const double    *x,
                            const cs_int_t  *incx,
                            const double    *y,
                            const cs_int_t  *incy);

/* Return the euclidean norm of a vector */

double CS_PROCF(dnrm2, DNRM2)(const cs_int_t  *n,
                              const double    *x,
                              const cs_int_t  *incx);

/* Scales a vector by a constant: x <-- ax */

void CS_PROCF(dscal, DSCAL)(const cs_int_t  *n,
                            const double    *a,
                            double          *x,
                            const cs_int_t  *incx);

/* Interchange vectors */

void CS_PROCF(dswap, DSWAP)(const cs_int_t  *n,
                            double          *restrict x,
                            const cs_int_t  *incx,
                            double          *restrict y,
                            const cs_int_t  *incy);

/* Finds the index of element having max absolute value */

cs_int_t CS_PROCF(idamax, IDAMAX)(const cs_int_t  *n,
                                  const double    *x,
                                  const cs_int_t  *incx);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLAS_H__ */
