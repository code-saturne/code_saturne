/*============================================================================
 * Define conjuguate heat transfer couplings with the SYRTHES code
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_sat_coupling.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define couplings with other instances of Code_Saturne.
 *
 * This is done by calling the cs_sat_coupling_define() function for each
 * coupling to add.
 *
 * The arguments to cs_sat_coupling_define are:
 *   saturne_name          <-- matching Code_Saturne application name
 *   volume_sup_criteria   <-- cell selection criteria for support
 *   boundary_sup_criteria <-- boundary face selection criteria for support
 *                             (not functional)
 *   volume_cpl_criteria   <-- cell selection criteria for coupled cells
 *   boundary_cpl_criteria <-- boundary face selection criteria for coupled
 *                             faces
 *   verbosity             <-- verbosity level
 *
 * In the case of only 2 Code_Saturne instances, the 'saturne_name' argument
 * is ignored, as there is only one matching possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * Code_Saturne instances based on the 'saturne_name' argument.
 *----------------------------------------------------------------------------*/

void
cs_user_saturne_coupling(void)
{
  int  verbosity = 1;
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /*-------------------------------------------------------------------------
   * Example 1: coupling with instance "SATURNE_01".
   *
   * - coupled faces of groups "3" or "4"
   * - all cells available as location support
   *-------------------------------------------------------------------------*/

  cs_sat_coupling_define("SATURNE_01",
                         "all[]",
                         NULL,
                         NULL,
                         "3 or 4",
                         verbosity);

  /*-------------------------------------------------------------------------
   * Example 2: coupling with instance "SATURNE_03".
   *
   * - coupled faces of groups "coupled_faces"
   * - coupled cells (every cell overlapping the distant mesh)
   * - all cells available as location support
   *-------------------------------------------------------------------------*/

  cs_sat_coupling_define("SATURNE_03",
                         "all[]",
                         NULL,
                         "all[]",
                         "coupled_faces",
                         verbosity);
}
