/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User operations called at the end of each time step.
 *
 * This function has a very general purpose, although it is recommended to
 * handle mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations([[maybe_unused]] cs_domain_t  *domain)
{
  /*! [extra_clipping] */

  /* Get total number of fields */
  const int n_fields = cs_field_n_fields();

  /* Loop over all fields */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);

    /* Filter fields of type variable (i.e. solved variables) */
    if (f->type & CS_FIELD_VARIABLE) {

      /* Retrieve solving information (which contains globally reduced
         clipping statistics for this iteration) */
      const cs_solving_info_t *sinfo = cs_field_get_solving_info_const(f);

      if (sinfo == nullptr)
        continue;

      /* Check if clipping occurred on this field */
      if (sinfo->n_clip_min > 0 || sinfo->n_clip_max > 0) {
        bft_printf
          ("cs_user_extra_operations: field '%s' clipped:\n"
           "  clips to min: %llu (pre-clip min: %14.5e)\n"
           "  clips to max: %llu (pre-clip max: %14.5e)\n",
           f->name,
           (unsigned long long)sinfo->n_clip_min, sinfo->min_pre_clip[0],
           (unsigned long long)sinfo->n_clip_max, sinfo->max_pre_clip[0]);

        /* For multi-component fields (e.g. vectors, tensors),
           per-component clipping counts are also available */
        if (f->dim > 1) {
          int n_comp = cs::min(f->dim, CS_SOLVING_INFO_MAX_DIM);
          for (int c_id = 0; c_id < n_comp; c_id++) {
            if (   sinfo->n_clip_min_comp[c_id] > 0
                || sinfo->n_clip_max_comp[c_id] > 0) {
              bft_printf
                ("    component %d: %llu min, %llu max "
                 "(range: [%14.5e, %14.5e])\n",
                 c_id,
                 (unsigned long long)sinfo->n_clip_min_comp[c_id],
                 (unsigned long long)sinfo->n_clip_max_comp[c_id],
                 sinfo->min_pre_clip[c_id],
                 sinfo->max_pre_clip[c_id]);
            }
          }
        }
      }
    }
  }

  /*! [extra_clipping] */
}

/*----------------------------------------------------------------------------*/
