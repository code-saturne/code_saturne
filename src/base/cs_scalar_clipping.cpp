/*============================================================================
 * Clip scalar (i.e. convected/diffused) fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_field.h"
#include "base/cs_log_iteration.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_scalar_clipping.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_scalar_clipping.cpp

  \brief Clipping scalar field.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran function prototypes for subroutines from field.f90.
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_scalar_clipping(int   id);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Clip a scalar variable field.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id       <-- field id
 *----------------------------------------------------------------------------*/

void
cs_f_scalar_clipping(int   id)
{
  cs_field_t *f = cs_field_by_id(id);

  cs_scalar_clipping(f);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip a scalar variable field.
 *
 * \param[in]  f  pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_scalar_clipping(cs_field_t  *f)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const int kscavr = cs_field_key_id("first_moment_id");
  const int kclvfl = cs_field_key_id_try("variance_clipping");
  const int kscmin = cs_field_key_id_try("min_scalar_clipping");
  const int kscmax = cs_field_key_id_try("max_scalar_clipping");

  int kclipp = cs_field_key_id("clipping_id");

  /* Post-process clippings ? */
  int clip_f_id = cs_field_get_key_int(f, kclipp);
  cs_real_t *cpro_f_clipped = NULL;
  if (clip_f_id > -1) {
    cpro_f_clipped = cs_field_by_id(clip_f_id)->val;
    cs_array_real_fill_zero(n_cells, cpro_f_clipped);
  }

  cs_real_t *cvar_scal = f->val;

  int variance_id = cs_field_get_key_int(f, kscavr);

  /* Logging and clippings
   * --------------------- */

  /* Variances are always clipped at positive values */

  /* Compute min and max values */
  cs_real_t *vmin;
  cs_real_t *vmax;

  CS_MALLOC(vmin, f->dim, cs_real_t);
  CS_MALLOC(vmax, f->dim, cs_real_t);

  for (cs_lnum_t i = 0; i < f->dim; i++) {

    vmin[i] = HUGE_VAL;
    vmax[i] = -HUGE_VAL;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      vmin[i] = cs::min(cvar_scal[f->dim * c_id + i], vmin[i]);
      vmax[i] = cs::max(cvar_scal[f->dim * c_id + i], vmax[i]);
    }
  }

 /* Clipping of non-variance scalars
  * And first clippings for the variances
  * (usually 0 for min and +grand for max) */

  cs_lnum_t *iclmax;
  cs_lnum_t *iclmin;
  CS_MALLOC(iclmin, f->dim, cs_lnum_t);
  CS_MALLOC(iclmax, f->dim, cs_lnum_t);

  for (cs_lnum_t i = 0; i < f->dim; i++) {
    iclmin[i] = 0;
    iclmax[i] = 0;
  }

  const cs_real_t scminp = cs_field_get_key_double(f, kscmin);
  const cs_real_t scmaxp = cs_field_get_key_double(f, kscmax);

  if (scmaxp > scminp && f->dim == 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (cvar_scal[c_id] > scmaxp) {
        iclmax[0] += 1;
        if (cpro_f_clipped != NULL)
          cpro_f_clipped[c_id] = cvar_scal[c_id] - scmaxp;

        cvar_scal[c_id] = scmaxp;
      }
      if (cvar_scal[c_id] < scminp) {
         iclmin[0] += 1;
         if (cpro_f_clipped != NULL)
           cpro_f_clipped[c_id] = cvar_scal[c_id] - scminp;
         cvar_scal[c_id] = scminp;
      }
    }
  }

  /* Clipping of max of variances
   * Based on associated scalar values (or 0 at min.) */
  /* Clipping of variances */

  if (variance_id > -1) {

    const int iclvfl = cs_field_get_key_int(f, kclvfl);

    if (iclvfl == 1) {

      iclmax[0] = 0;

      /* Get the corresponding scalar */
      cs_field_t *fl = cs_field_by_id(variance_id);
      cs_real_t *cvar_scav = fl->val;

      /* Get the min clipping of the corresponding scalar */
      const cs_real_t scmin = cs_field_get_key_double(fl, kscmin);
      const cs_real_t scmax = cs_field_get_key_double(fl, kscmax);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t vfmax = (cvar_scav[c_id]-scmin)*(scmax-cvar_scav[c_id]);
        if (cvar_scal[c_id] > vfmax) {
          iclmax[0] += 1;
          cvar_scal[c_id] = vfmax;
        }
      }
    }

  }

  /* Save counts and extrema for delayed (batched) global reduction */

  cs_log_iteration_clipping_field(f->id,
                                  iclmin[0],
                                  iclmax[0],
                                  vmin,
                                  vmax,
                                  iclmin,
                                  iclmax);


  CS_FREE(vmin);
  CS_FREE(vmax);
  CS_FREE(iclmin);
  CS_FREE(iclmax);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
