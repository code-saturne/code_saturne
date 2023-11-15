/*============================================================================
 * Initialize fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_array.h"
#include "cs_assert.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_scalar_clipping.h"
#include "cs_turbulence_init.h"
#include "cs_vof.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_initialize_fields.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_initialize_fields.c

  \brief Various field initializations.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computed variable and property initialization.
 *
 * The time step, the indicator of wall distance computation are also
 * initialized just before reading a restart file or use the user
 * initializations.
 */
/*----------------------------------------------------------------------------*/

void
cs_initialize_fields_stage_0(void)
{
  const cs_mesh_t *m = cs_glob_mesh;

  const int n_fields = cs_field_n_fields();

  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kdflim = cs_field_key_id("diffusion_limiter_id");

  /* Reference value for P*
     ---------------------- */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           cs_glob_fluid_properties->pred0,
                           CS_F_(p)->val);

  /* VoF
     --- */

  /* Note that the initialization here should not be strictly
     necessary, as cs_scalar_clipping should handle it, but
     doing so here will avoid activating the clipping
     counters; maybe a better manner of doing things would
     be for the counters to be deactivated at this step,
     which comes before user settings */

  if (cs_glob_vof_parameters->vof_model > 0) {
    const int kscmin = cs_field_key_id_try("min_scalar_clipping");

    cs_field_t *f_vf = cs_field_by_name("void_fraction");
    cs_real_t clvfmn = cs_field_get_key_double(f_vf, kscmin);

    cs_array_real_set_scalar(m->n_cells_with_ghosts, clvfmn, f_vf->val);
  }

  /* Turbulence
     ---------- */

  cs_turbulence_init_by_ref_quantities();

  /* Clip scalar variables (except k-epsilon, see above)
     --------------------------------------------------- */

  /* Clipping of non-variance scalars */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (   (! (f->type & CS_FIELD_VARIABLE))
        || (f->type & CS_FIELD_CDO)) /* CDO variables handled elsewhere */
      continue;
    if (cs_field_get_key_int(f, keysca) < 0)
      continue;
    if (cs_field_get_key_int(f, kscavr) >= 0)  /* is a variance */
      continue;

    cs_scalar_clipping(f);
  }

  /* Clipping of variances  */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (   (! (f->type & CS_FIELD_VARIABLE))
        || (f->type & CS_FIELD_CDO)) /* CDO variables handled elsewhere */
      continue;
    if (cs_field_get_key_int(f, keysca) < 0)
      continue;

    if (cs_field_get_key_int(f, kscavr) >= 0)
      cs_scalar_clipping(f);
  }

  /* ALE
     --- */

  if (cs_glob_ale > 0) {
    cs_field_t *f_coord0 = cs_field_by_name("vtx_coord0");
    cs_array_real_copy(m->n_vertices*3, m->vtx_coord, f_coord0->val);
  }

  /* Specific field initializations
     ------------------------------ */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (! (f->type & CS_FIELD_VARIABLE))
      continue;

    if (f->type & CS_FIELD_CDO) /* CDO variables handled elsewhere */
      continue;

    /* Diffusion limiters;
       Note that allowing a default value initialization for 1D-fields
       would allow handling this upon field definition */

    const int dl_id = cs_field_get_key_int(f, kdflim);
    if (dl_id > -1) {
      cs_field_t *f_dl = cs_field_by_id(dl_id);
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 1.0, f_dl->val);
    }

  }

  /* Previous value initializations
     ------------------------------

     This is not restricted to variables
     (for example Cp, Cv, and mass fluxes may be handled here also,
     depending on the time stepping options). */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (f->type & CS_FIELD_CDO) /* CDO variables handled elsewhere */
      continue;

    for (int n_p = f->n_time_vals; n_p > 1; n_p--)
      cs_field_current_to_previous(f);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
