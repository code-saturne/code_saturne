/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_parameters.c
        General parameters and options management.
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_var_cal_opt_t _var_cal_opt =
{
  0,     /* iwarni */
  1,     /* iconv  */
  1,     /* istat  */
  1,     /* idiff  */
  1,     /* idifft */
  1,     /* idften */
  0,     /* iswdyn */
  1,     /* ischcv */
  1,     /* isstpc */
  100,   /* nswrgr */
  1,     /* nswrsm */
  0,     /* imrgra */
  -1,    /* imligr */
  1,     /* ircflu */
  1.,    /* thetav */
  1.,    /* blencv */
  1.e-8, /* epsilo */
  1.e-7, /* epsrsm */
  1.e-5, /* epsrgr */
  1.5,   /* climgr */
  0.,    /* extrag */
  1.     /* relaxv */
};

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log default values of the structure */

static void
_log_func_var_opt_cal(const void *t)
{
  const cs_var_cal_opt_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "iwarni", _t->iwarni);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "iconv ", _t->iconv );
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "istat ", _t->istat );
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "idiff ", _t->idiff );
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "idifft", _t->idifft);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "idften", _t->idften);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "iswdyn", _t->iswdyn);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "ischcv", _t->ischcv);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "isstpc", _t->isstpc);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "nswrgr", _t->nswrgr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "nswrsm", _t->nswrsm);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "imrgra", _t->imrgra);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "imligr", _t->imligr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "ircflu", _t->ircflu);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "thetav", _t->thetav);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "blencv", _t->blencv);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "epsilo", _t->epsilo);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "epsrsm", _t->epsrsm);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "epsrgr", _t->epsrgr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "climgr", _t->climgr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "extrag", _t->extrag);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "relaxv", _t->relaxv);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define general field keys.
 *
 * A recommened practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_keys(void)
{
  cs_field_define_key_int("inner_mass_flux_id", -1, 0);
  cs_field_define_key_int("boundary_mass_flux_id", -1, 0);

  cs_field_define_key_int("variable_id", -1, 0); /* inverse of ivarfl(ivar) */
  cs_field_define_key_int("post_id", -1, 0);     /* inverse of the ipp array */

  cs_field_define_key_int("scalar_diffusivity_id", -1, 0);
  cs_field_define_key_int("diffusivity_tensor", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("scalar_id", -1, 0); /* inverse of isca(iscal) */
  cs_field_define_key_int("drift_scalar_model", 0, 0);
  cs_field_define_key_int("scalar_class", -1, 0);
  cs_field_define_key_int("first_moment_id", -1, 0); // old iscavr(iscal)
  cs_field_define_key_double("min_scalar_clipping", -1.e12, 0);
  cs_field_define_key_double("max_scalar_clipping", 1.e12, 0);

  cs_field_define_key_int("property_id", -1, 0); /* inverse of iprpfl(iprop) */

  cs_field_define_key_struct("var_cal_opt",
                             &_var_cal_opt,
                             _log_func_var_opt_cal,
                             sizeof(cs_var_cal_opt_t),
                             CS_FIELD_VARIABLE);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read general restart info.
 *
 * This updates the previous time step info.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_read_restart_info(void)
{
  if (cs_restart_present()) {
    cs_restart_t *r
      = cs_restart_create("main", "restart", CS_RESTART_MODE_READ);
    cs_restart_read_time_step_info(r);
    r = cs_restart_destroy(r);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
