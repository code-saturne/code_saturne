/*============================================================================
 * General parameters management for groundwater flow module.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <stdarg.h>
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

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_gwf_parameters.c
        General parameters and options management for ground water flow module.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* soilwater partition model structure */

static cs_gwf_soilwater_partition_t
_gwf_soilwater_partition =
{
  .kinetic =  0, /* 0: Kd model, 1: EK model */
  .ikd     = -1,
  .idel    = -1,
  .ikp     = -1,
  .ikm     = -1,
  .isorb   = -1
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log default values of the structure */

static void
_log_func_gwf_soilwater_partition(const void *t)
{
  const char fmt_i[] = N_("      %-19s  %-4d\n");
  const cs_gwf_soilwater_partition_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "kinetic ", _t->kinetic);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "ikd     ", _t->ikd);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "idel    ", _t->idel);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "ikp     ", _t->ikp);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "ikm     ", _t->ikm);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "isorb   ", _t->isorb);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define field key for soilwater partition model.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_parameters_define_field_key_soilwater_partition(void)
{
  cs_field_define_key_struct("gwf_soilwater_partition",
                             &_gwf_soilwater_partition,
                             _log_func_gwf_soilwater_partition,
                             sizeof(cs_gwf_soilwater_partition_t),
                             0);
 }

/*----------------------------------------------------------------------------*/

END_C_DECLS
