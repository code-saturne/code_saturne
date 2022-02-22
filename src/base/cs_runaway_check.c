/*============================================================================
 * Runaway (diverging) computation detection.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_runaway_check.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_runaway_check

  \brief  Runaway (diverging) computation detection.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/* Error threshold data */
/*----------------------*/

typedef struct {

  int     f_id;         /* Associated field id, or -1 */
  int     nt_test;      /* Last tested time step */
  double  max_allowed;  /* Maximum allowed */
  double  max_reached;  /* Maximum reached */

} cs_log_check_bounds_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_log_check_bounds_t  _check_bounds = {.f_id = -1,
                                               .nt_test = -1,
                                               .max_allowed = HUGE_VAL,
                                               .max_reached = -HUGE_VAL};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update bounds check.
 *----------------------------------------------------------------------------*/

static void
_update_check_bounds(cs_log_check_bounds_t *cb)
{
  if (cb->nt_test >= cs_glob_time_step->nt_cur)
    return;

  const cs_field_t *f = cs_field_by_id(cb->f_id);
  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
  const cs_lnum_t _n_elts = n_elts[0];
  int _dim = (f->dim == 3) ? 4 : f->dim;
  if (_dim > 6)
    _dim = 6;

  cs_real_t vmax = {-HUGE_VAL};

  if (f->dim == 3) {
    const cs_real_3_t *val = (const cs_real_3_t *)(f->val);
    for (cs_lnum_t i = 0; i < _n_elts; i++) {
      cs_real_t v = cs_math_3_square_norm(val[i]);
      if (v > vmax)
        vmax = v;
    }
    if (vmax > 0)
      vmax = sqrt(vmax);
  }
  else if (f->dim == 6) {
    const cs_real_6_t *val = (const cs_real_6_t *)(f->val);
    for (cs_lnum_t i = 0; i < _n_elts; i++) {
      cs_real_t v = val[i][0] + val[i][1] + val[i][2];
      if (v > vmax)
        vmax = v;
    }
  }
  else {
    const cs_real_t *val = f->val;
    const cs_lnum_t n = _n_elts*((cs_lnum_t)(f->dim));
    for (cs_lnum_t i = 0; i < n; i++) {
      if (val[i] > vmax)
        vmax = val[i];
    }
  }

  cs_parall_max(1, CS_REAL_TYPE, &vmax);

  cb->nt_test = cs_glob_time_step->nt_cur;
  cb->max_reached = vmax;

  if (cb->max_reached > cb->max_allowed) {
    bft_printf(_("\n"
                 "Error (runaway computation):\n"
                 "-----\n"
                 "  Maximum allowed value of %g exceeded for field %s.\n"),
               cb->max_allowed, f->name);
    int nt_cur = cs_glob_time_step->nt_cur;
    cs_time_step_define_nt_max(nt_cur);
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that defined field bounds are not exceeded.
 *
 * \return 0 if no bounds are exceeded, 1 if bounds are exceeded.
 */
/*----------------------------------------------------------------------------*/

int
cs_runaway_check(void)
{
  int retval = 0;

  if (_check_bounds.f_id > -1) {

    _update_check_bounds(&_check_bounds);

    if (_check_bounds.max_reached > _check_bounds.max_allowed) {
      int nt_cur = cs_glob_time_step->nt_cur;
      cs_time_step_define_nt_max(nt_cur);
      retval += 1;
    }

  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define maximum value for a field, beyon which computation is aborted.
 *
 * Currently, only one field is handled, so calling this multiple times
 * replaced the previous setting. Using a negative field id removes
 * this check.
 */
/*----------------------------------------------------------------------------*/

void
cs_runaway_check_define_field_max(int        f_id,
                                  cs_real_t  max_allowed)
{
  _check_bounds.f_id = f_id;
  _check_bounds.max_allowed = max_allowed;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that defined field bounds are not exceeded.
 */
/*----------------------------------------------------------------------------*/

void
cs_runaway_check_finalize(void)
{
  int n_checks = 1, n_errors = 0;

  /* Note: we may turn this into a loop if we allow for more than one
     bound checks */

  for (int i = 0; i < n_checks; i++) {
    cs_log_check_bounds_t  *cb = &_check_bounds + i;
    if (cb->max_reached > cb->max_allowed)
      n_errors += 1;
  }

  if (n_errors > 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Runaway computation:\n"
                "  At least one field exceeded allowed bounds (see log)."));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
