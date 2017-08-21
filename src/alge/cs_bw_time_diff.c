/*============================================================================
 * Time operators.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

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

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_bw_time_diff.c

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_backward_differentiation_in_time(const int     field_id,
                                    cs_real_t    *exp_part,
                                    cs_real_t    *imp_part);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

void
cs_backward_differentiation_in_time(const int     field_id,
                                    cs_real_t    *exp_part,
                                    cs_real_t    *imp_part)
{
  cs_lnum_t iel;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *rho = CS_F_(rho)->val;

  const cs_field_t *f = cs_field_by_id(field_id);
  const int dim = f->dim;

  if (dim == 3) {
    cs_real_3_t *exp_part_3 = (cs_real_3_t *)(exp_part);
    cs_real_33_t *imp_part_33 = (cs_real_33_t *)(imp_part);
    for (iel = 0; iel < n_cells; iel++) {
      for (int jj = 0; jj < 3; jj++) {
        exp_part_3[iel][jj] += rho[iel]*cell_vol[iel]/dt[iel]
                             *(f->vals[1][dim*iel + jj]
                             - 0.5*f->vals[2][dim*iel + jj]);
        imp_part_33[iel][jj][jj] += -0.5*rho[iel]
                                  *cell_vol[iel]/dt[iel];
      }
    }
  }
  else {
    for (iel = 0; iel < n_cells; iel++) {
      exp_part[iel] += rho[iel]*cell_vol[iel]/dt[iel]
                               *(f->vals[1][dim*iel]
                               - 0.5*f->vals[2][dim*iel]);
      imp_part[iel] += -0.5*rho[iel]*cell_vol[iel]/dt[iel];
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
