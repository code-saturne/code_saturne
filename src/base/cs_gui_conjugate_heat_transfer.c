/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Management of the GUI parameters file: conjugate heat transfer
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvm_selector.h"

/*----------------------------------------------------------------------------
 * MEI library headers
 *----------------------------------------------------------------------------*/

#ifdef HAVE_MEI
#include "mei_evaluate.h"
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

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
 * parameter:
 *   keyword    --> label of the markup
 *   number     --> number of the syrthes coupling
 *----------------------------------------------------------------------------*/

static char*
_get_syrthes_coupling(const char* keyword, const int number)
{
    char* value = NULL;
    char *path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "thermophysical_models",
                                    "conjugate_heat_transfer",
                                    "external_coupling");
    cs_xpath_add_element_num(&path, "syrthes", number);
    cs_xpath_add_element(&path, keyword);
    cs_xpath_add_function_text(&path);
    value = cs_gui_get_attribute_value(path);
    BFT_FREE(path);
    return value;
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * syrthes_app_num and syrthes_name arguments are ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument, then
 * on the syrthes_app_num argument. If syrthes_name is empty, matching will
 * be based on syrthes_app_num only.
 *
 * subroutine uisyrc
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uisyrc, UISYRC) (void)
{
    int izone;
    int verbosity = 0;
    int app_num = 0;
    char* syrthes_name = NULL;
    char* syrthes_app_num = NULL;
    char* projection_axis = NULL;
    char* boundary_criteria = NULL;
    char* volume_criteria = NULL;

    int coupling =
      cs_gui_get_tag_number("/conjugate_heat_transfer/external_coupling/syrthes",
                            1);

    for (izone=0 ; izone < coupling ; izone++)
    {
      syrthes_name      = _get_syrthes_coupling("syrthes_name",       izone+1);
      syrthes_app_num   = _get_syrthes_coupling("syrthes_app_num",    izone+1);
      projection_axis   = _get_syrthes_coupling("projection_axis",    izone+1);
      boundary_criteria = _get_syrthes_coupling("selection_criteria", izone+1);

      app_num = atoi(syrthes_app_num);

      cs_syr_coupling_define(app_num,
                             syrthes_name,
                             boundary_criteria,
                             volume_criteria,
                             *projection_axis,
                             verbosity);

#if _XML_DEBUG_
      bft_printf("==>uisyrc\n");
      bft_printf("--syrthes_name      = %s\n", syrthes_name);
      bft_printf("--syrthes_app_num   = %s\n", syrthes_app_num);
      bft_printf("--boundary_criteria = %s\n", boundary_criteria);
      bft_printf("--projection_axis   = %s\n", projection_axis);
#endif
      BFT_FREE(syrthes_name);
      BFT_FREE(syrthes_app_num);
      BFT_FREE(projection_axis);
      BFT_FREE(boundary_criteria);

    }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
