/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
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

/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * External library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

#if defined(_CS_HAVE_CBLAS)
#include <cblas.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function definitions (see BLAS reference)
 *============================================================================*/

#if     defined(_CS_HAVE_F77BLAS) && !defined(_CS_HAVE_CBLAS) \
    && !defined(_CS_HAVE_ESSL) && !defined (_CS_HAVE_MKL)

/* If we have F77 BLAS but no C BLAS functions, define C BLAS as F77 wrappers */
/*----------------------------------------------------------------------------*/

/* Sum of the absolute values of a vector */

double cblas_dasum(cs_int_t       n,
                   const double  *x,
                   cs_int_t       incx)
{
  return CS_PROCF(dasum, DASUM)(&n, x, &incx);
}

/* Constant times a vector plus a vector: y <-- ax + y */

void cblas_daxpy(cs_int_t       n,
                 double         a,
                 const double  *x,
                 cs_int_t       incx,
                 double        *restrict y,
                 cs_int_t       incy)
{
  CS_PROCF(daxpy, DAXPY)(&n, &a, x, &incx, y, &incy);
}

/* Copy a vector x to a vector y: y <-- x */

void cblas_dcopy(cs_int_t       n,
                 const double  *x,
                 cs_int_t       incx,
                 double        *restrict y,
                 cs_int_t       incy)
{
  CS_PROCF(dcopy, DCOPY)(&n, x, &incx, y, &incy);
}

/* Return the dot product of 2 vectors: x.y */

double cblas_ddot(cs_int_t       n,
                  const double  *x,
                  cs_int_t       incx,
                  const double  *y,
                  cs_int_t       incy)
{
  return CS_PROCF(ddot, DDOT)(&n, x, &incx, y, &incy);
}

/* Return the euclidean norm of a vector */

double cblas_dnrm2(cs_int_t       n,
                   const double  *x,
                   cs_int_t       incx)
{
  return CS_PROCF(dnrm2, DNRM2)(&n, x, &incx);
}

/* Scales a vector by a constant: x <-- ax */

void cblas_dscal(cs_int_t   n,
                 double     a,
                 double    *x,
                 cs_int_t   incx)
{
  CS_PROCF(dscal, DSCAL)(&n, &a, x, &incx);
}

/* Interchange vectors */

void cblas_dswap(cs_int_t   n,
                 double    *restrict x,
                 cs_int_t   incx,
                 double    *restrict y,
                 cs_int_t   incy)
{
  CS_PROCF(dswap, DSWAP)(&n, x, &incx, y, &incy);
}

/* Finds the index of element having max absolute value */

cs_int_t cblas_idamax(cs_int_t       n,
                      const double  *x,
                      cs_int_t       incx)
{
  return CS_PROCF(idamax, IDAMAX)(&n, x, &incx);
}

#endif  /*     defined(_CS_HAVE_F77BLAS) && !defined(_CS_HAVE_CBLAS) \
           && !defined(_CS_HAVE_ESSL) && !defined (_CS_HAVE_MKL) */

#if    !defined(_CS_HAVE_F77BLAS) \
    && !defined(_CS_HAVE_ESSL) && !defined (_CS_HAVE_MKL)

/* If we have no F77 BLAS functions, define F77 BLAS as C wrappers */
/*-----------------------------------------------------------------*/

/* Sum of the absolute values of a vector */

double CS_PROCF(dasum, DASUM)(const cs_int_t  *n,
                              const double    *x,
                              const cs_int_t  *incx)
{
  return cblas_dasum(*n, x, *incx);
}

/* Constant times a vector plus a vector: y <-- ax + y */

void CS_PROCF(daxpy, DAXPY)(const cs_int_t  *n,
                            const double    *a,
                            const double    *x,
                            const cs_int_t  *incx,
                            double          *restrict y,
                            const cs_int_t  *incy)
{
  cblas_daxpy(*n, *a, x, *incx, y, *incy);
}

/* Copy a vector x to a vector y: y <-- x */

void CS_PROCF(dcopy, DCOPY)(const cs_int_t  *n,
                            const double    *x,
                            const cs_int_t  *incx,
                            double          *restrict y,
                            const cs_int_t  *incy)
{
  cblas_dcopy(*n, x, *incx, y, *incy);
}

/* Return the dot product of 2 vectors: x.y */

double CS_PROCF(ddot, DDOT)(const cs_int_t  *n,
                            const double    *x,
                            const cs_int_t  *incx,
                            const double    *y,
                            const cs_int_t  *incy)
{
  return cblas_ddot(*n, x, *incx, y, *incy);
}

/* Return the euclidean norm of a vector */

double CS_PROCF(dnrm2, DNRM2)(const cs_int_t  *n,
                              const double    *x,
                              const cs_int_t  *incx)
{
  return cblas_dnrm2(*n, x, *incx);
}

/* Scales a vector by a constant: x <-- ax */

void CS_PROCF(dscal, DSCAL)(const cs_int_t  *n,
                            const double    *a,
                            double          *x,
                            const cs_int_t  *incx)
{
  cblas_dscal(*n, *a, x, *incx);
}

/* Interchange vectors */

void CS_PROCF(dswap, DSWAP)(const cs_int_t  *n,
                            double          *restrict x,
                            const cs_int_t  *incx,
                            double          *restrict y,
                            const cs_int_t  *incy)
{
  cblas_dswap(*n, x, *incx, y, *incy);
}

/* Finds the index of element having max absolute value */

cs_int_t CS_PROCF(idamax, IDAMAX)(const cs_int_t  *n,
                                  const double    *x,
                                  const cs_int_t  *incx)
{
  return cblas_idamax(*n, x, *incx);
}

#endif /*    !defined(_CS_HAVE_F77BLAS) \
          && !defined(_CS_HAVE_ESSL) && !defined (_CS_HAVE_MKL) */

#if !defined(_CS_HAVE_BLAS)

/* If we have no external BLAS, define fallback legacy C BLAS */
/*------------------------------------------------------------*/

/* Return the sum of absolute values */

double cblas_dasum(cs_int_t       n,
                   const double  *x,
                   cs_int_t       incx)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);

  double     sum = 0;

  if (n < 0) return sum;

  if (inc_x == 1) {
    for (i = 0; i < n; i++) {
      sum += CS_ABS(x[i]);
    }
  }
  else {
    cs_int_t j;
    for (i = 0, j = 0; i < n; i++, j += inc_x) {
      sum += CS_ABS(x[j]);
    }
  }

  return sum;
}

/* Constant times a vector plus a vector: y <-- ax + y */

void cblas_daxpy(cs_int_t       n,
                 double         a,
                 const double  *x,
                 cs_int_t       incx,
                 double        *restrict y,
                 cs_int_t       incy)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);
  cs_int_t   inc_y = CS_ABS(incy);

  if (n < 0) return;

  if (inc_x == 1 && inc_y == 1) {
    for (i = 0; i < n; i++) {
      y[i] += (a * x[i]);
    }
  }
  else {
    cs_int_t j, k;
    for (i = 0, j = 0, k = 0; i < n; i++, j += inc_x, k += inc_y) {
      y[k] += (a * x[j]);
    }
  }
}

/* Copy a vector x to a vector y: y <-- x */

void cblas_dcopy(cs_int_t       n,
                 const double  *x,
                 cs_int_t       incx,
                 double        *restrict y,
                 cs_int_t       incy)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);
  cs_int_t   inc_y = CS_ABS(incy);

  if (n < 0) return;

  if (inc_x == 1 && inc_y == 1) {
    for (i = 0; i < n; i++) {
      y[i] = x[i];
    }
  }
  else {
    cs_int_t j, k;
    for (i = 0, j = 0, k = 0; i < n; i++, j += inc_x, k += inc_y) {
      y[k] = x[j];
    }
  }
}

/* Return the dot product of 2 vectors: x.y */

double cblas_ddot(cs_int_t       n,
                  const double  *x,
                  cs_int_t       incx,
                  const double  *y,
                  cs_int_t       incy)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);
  cs_int_t   inc_y = CS_ABS(incy);

  double     sum = 0;

  if (n < 0) return sum;

  if (inc_x == 1 && inc_y == 1) {
    for (i = 0; i < n; i++) {
      sum += (x[i] * y[i]);
    }
  }
  else {
    cs_int_t j, k;
    for (i = 0, j = 0, k = 0; i < n; i++, j += inc_x, k += inc_y) {
      sum += (x[j] * y[k]);
    }
  }

  return sum;
}

/* Return the euclidean norm of a vector */

double cblas_dnrm2(cs_int_t         n,
                   const double    *x,
                   cs_int_t         incx)
{
  cs_int_t   i, j;
  cs_int_t   inc_x = CS_ABS(incx);
  double     absx, temp;
  double     scale = 1e-18; /* Close to zero */
  double     ssq = 0.0;

  if (n < 0) return ssq;

  for (i = 0, j = 0; i < n; i++, j += inc_x) {
    absx = CS_ABS(x[j]);
    if (scale < absx) {
      temp = scale / absx;
      scale = absx;
      ssq = 1.0 + (ssq * (temp*temp));
    }
    else {
      temp = absx / scale;
      ssq += (temp * temp);
    }
  }

  return (scale * sqrt(ssq));
}

/* Scales a vector by a constant: x <-- ax */

void cblas_dscal(cs_int_t         n,
                 double           a,
                 double          *x,
                 cs_int_t         incx)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);

  if (n < 0) return;

  if (inc_x == 1) {
    for (i = 0; i < n; i++) {
      x[i] *= a;
    }
  }
  else {
    cs_int_t j;
    for (i = 0, j = 0; i < n; i++, j += inc_x) {
      x[j] *= a;
    }
  }
}

/* Interchange vectors */

void cblas_dswap(cs_int_t   n,
                 double    *restrict x,
                 cs_int_t   incx,
                 double    *restrict y,
                 cs_int_t   incy)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);
  cs_int_t   inc_y = CS_ABS(incy);

  double     temp;

  if (n < 0) return;

  if (inc_x == 1 && inc_y == 1) {
    for (i = 0; i < n; i++) {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }
  }
  else {
    cs_int_t j, k;
    for (i = 0, j = 0, k = 0; i < n; i++, j += inc_x, k += inc_y) {
      temp = x[j];
      x[j] = y[k];
      y[k] = temp;
    }
  }
}

/* Finds the index of element having max absolute value */

cs_int_t cblas_idamax(cs_int_t       n,
                      const double  *x,
                      cs_int_t       incx)
{
  cs_int_t   i;
  cs_int_t   inc_x = CS_ABS(incx);
  double     dmax;

  cs_int_t   index_max = 0;

  if (n < 1 || inc_x < 0) return index_max;

  index_max = 1;
  dmax = CS_ABS(x[0]);

  if (inc_x == 1) {
    for (i = 0; i < n; i++) {
      if (CS_ABS(x[i]) > dmax) {
        index_max = i+1;
        dmax= CS_ABS(x[i]);
      }
    }
  }
  else {
    cs_int_t j;
    for (i = 0, j = 0; i < n; i++, j += inc_x) {
      if (CS_ABS(x[j]) > dmax) {
        index_max = i+1;
        dmax= CS_ABS(x[j]);
      }
    }
  }

  return index_max;
}

#endif /* !defined(_CS_HAVE_CBLAS) && !defined(_CS_HAVE_F77BLAS) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
