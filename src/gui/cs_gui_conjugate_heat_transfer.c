/*============================================================================
 * Management of the GUI parameters file: conjugate heat transfer
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_mesh.h"
#include "cs_prototypes.h"
#include "cs_syr_coupling.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_conjugate_heat_transfer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the value to the asociated markup.
 *
 * parameters:
 *   keyword  <-- label of the markup
 *   number   <-- number of the syrthes coupling
 *----------------------------------------------------------------------------*/

static char*
_get_syrthes_coupling(const char*  keyword,
                      int          number)
{
  char* value = NULL;
  char *path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "conjugate_heat_transfer",
                        "external_coupling");
  cs_xpath_add_element_num(&path, "syrthes", number);
  cs_xpath_add_element(&path, keyword);
  cs_xpath_add_function_text(&path);
  value = cs_gui_get_text_value(path);
  BFT_FREE(path);
  return value;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * syrthes_name argument is ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument.
 *
 *----------------------------------------------------------------------------*/

void
cs_gui_syrthes_coupling(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  int zone_id;
  int verbosity = 0;
  int visualization = 1;
  bool allow_nonmatching = false;
  double tolerance = 0.1;
  char* syrthes_name = NULL;
  char* syrthes_verbosity = NULL;
  char* syrthes_visu = NULL;
  char* syrthes_tolerance = NULL;
  char* syrthes_nonmatching = NULL;
  char* projection_axis = NULL;
  char* boundary_criteria = NULL;
  char* volume_criteria = NULL;

  int n_couplings
    = cs_gui_get_tag_count("/conjugate_heat_transfer/external_coupling/syrthes",
                           1);

  for (zone_id = 0; zone_id < n_couplings; zone_id++) {

    syrthes_name        = _get_syrthes_coupling("syrthes_name",       zone_id+1);
    syrthes_tolerance   = _get_syrthes_coupling("tolerance",          zone_id+1);
    syrthes_verbosity   = _get_syrthes_coupling("verbosity",          zone_id+1);
    syrthes_visu        = _get_syrthes_coupling("visualization",      zone_id+1);
    projection_axis     = _get_syrthes_coupling("projection_axis",    zone_id+1);
    syrthes_nonmatching = _get_syrthes_coupling("allow_nonmatching",  zone_id+1);
    boundary_criteria   = _get_syrthes_coupling("selection_criteria", zone_id+1);
    volume_criteria     = _get_syrthes_coupling("volume_criteria",    zone_id+1);

    if (syrthes_verbosity != NULL)
      verbosity = atoi(syrthes_verbosity);

    if (syrthes_visu != NULL)
      visualization = atoi(syrthes_visu);

    if (syrthes_tolerance != NULL)
      tolerance = atof(syrthes_tolerance);

    if (syrthes_nonmatching != NULL) {
      if (atoi(syrthes_nonmatching))
        allow_nonmatching = true;
    }

    cs_syr_coupling_define(syrthes_name,
                           boundary_criteria,
                           volume_criteria,
                           *projection_axis,
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           visualization);

#if _XML_DEBUG_
    bft_printf("==>uisyrc\n");
    bft_printf("--syrthes_name          = %s\n", syrthes_name);
    bft_printf("--boundary_criteria     = %s\n", boundary_criteria);
    bft_printf("--volume_criteria       = %s\n", volume_criteria);
    bft_printf("--syrthes_verbosity     = %s\n", syrthes_verbosity);
    bft_printf("--syrthes_visualization = %s\n", syrthes_visu);
    bft_printf("--syrthes_nonmatching   = %s\n", syrthes_nonmatching);
    bft_printf("--syrthes_tolerance     = %s\n", syrthes_tolerance);
    bft_printf("--projection_axis       = %s\n", projection_axis);
#endif
    BFT_FREE(syrthes_name);
    BFT_FREE(syrthes_verbosity);
    BFT_FREE(syrthes_visu);
    BFT_FREE(syrthes_tolerance);
    BFT_FREE(syrthes_nonmatching);
    BFT_FREE(projection_axis);
    BFT_FREE(boundary_criteria);
    BFT_FREE(volume_criteria);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
