#ifndef __CS_SLES_PC_PRIV_H__
#define __CS_SLES_PC_PRIV_H__

/*============================================================================
 * Sparse Linear Equation Solver preconditioner, private elements
 *
 * These elements are shared between host and device implementations,
 * but are not accessible to calling code.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_log.h"
#include "cs_halo_perio.h"
#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure for Jacobi or polynomial preconditioner */
/*---------------------------------------------------*/

typedef struct {

#if defined(HAVE_ACCEL)
  bool                 accelerated;       /* Use accelerated version ? */
#endif

  int                  poly_degree;       /* 0: Jacobi, > 0: polynomial */
  cs_lnum_t            n_rows;            /* Number of associated rows */
  cs_lnum_t            n_cols;            /* Number of associated columns */

  cs_lnum_t            n_aux;             /* Size of auxiliary data */

  const cs_matrix_t   *a;                 /* Pointer to associated matrix */
  const cs_real_t     *ad_inv;            /* pointer to diagonal inverse */
  cs_real_t           *_ad_inv;           /* private pointer to
                                             diagonal inverse */

  cs_real_t           *aux;               /* Auxiliary data */

} cs_sles_pc_poly_t;


/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_PC_PRIV_H__ */
