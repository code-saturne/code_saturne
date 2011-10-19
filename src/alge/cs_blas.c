/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
 * External library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

#if defined(HAVE_CBLAS)
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

BEGIN_C_DECLS

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

/* Return the dot product of 2 vectors: x.y */

double CS_PROCF(csdot, CSDOT)(const cs_int_t  *n,
                              const double    *x,
                              const double    *y)
{
  return cs_dot(*n, x, y);
}

#if !defined(HAVE_BLAS)

/* If we have no external BLAS, use local implementation */
/*-------------------------------------------------------*/

/* Constant times a vector plus a vector: y <-- ax + y */

void cs_axpy(cs_lnum_t      n,
             double         a,
             const double  *x,
             double        *restrict y)
{
  cs_lnum_t  i;

  if (n < 1)
    return;

# pragma omp parallel for
  for (i = 0; i < n; i++)
    y[i] += (a * x[i]);
}

/* Return the dot product of 2 vectors: x.y */

double cs_dot(cs_lnum_t      n,
              const double  *x,
              const double  *y)
{
  cs_lnum_t  i;
  double     s = 0;

  if (n < 1)
    return sum;

# pragma omp parallel for reduction(+: s)
  for (i = 0; i < n; i++) {
    s += (x[i] * y[i]);
  }

  return s;
}

#endif /* !defined(HAVE_BLAS) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
