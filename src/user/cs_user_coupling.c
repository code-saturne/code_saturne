/*============================================================================
 * Code couplings definition with SYRTHES and Code_Saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of Code_Saturne
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
#include "cs_syr_coupling.h"

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
 * Define couplings with SYRTHES code.
 *
 * This is done by calling the cs_syr_coupling_define() function for each
 * coupling to add.
 *
 * The arguments to cs_syr_coupling_define are:
 *   syrthes_name      <-- matching SYRTHES application name
 *   boundary_criteria <-- surface selection criteria
 *   volume  _criteria <-- volume selection criteria
 *   projection_axis   <-- ' ' : standard 3D coupling
 *                         'x', 'y', or 'z': projection axis for coupling
 *                                           with 2D SYRTHES.
 *   verbosity         <-- verbosity level
 *   plot              <-- visualization level
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * 'syrthes_name' argument is ignored, as there is only one matching
 * possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances based on the 'syrthes_name' argument.
 *----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling(void)
{
  int  verbosity = 1, plot = 1;
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /*-------------------------------------------------------------------------
   * Example 1:
   *
   * Boundary faces of group '3' coupled with instance named 'SYRTHES_01'.
   *-------------------------------------------------------------------------*/

  cs_syr_coupling_define("SYRTHES_01",
                         "3",             /* boundary criteria */
                         NULL,            /* volume_criteria */
                         ' ',             /* projection_axis */
                         verbosity,
                         plot);

  /*-------------------------------------------------------------------------
   * Example 2:
   *
   * Boundary faces of group 'Wall' coupled with 2D SYRTHES instance
   * named 'SYRTHES_02'.
   *-------------------------------------------------------------------------*/

  cs_syr_coupling_define("SYRTHES_02",
                         "Wall",          /* boundary criteria */
                         NULL,            /* volume_criteria */
                         'z',             /* projection_axis */
                         verbosity,
                         plot);

  /*-------------------------------------------------------------------------
   * Example 3:
   *
   * Cells in box with corners (0, 0, 0) and (1, 1, 1) coupled with
   * SYRTHES instance named 'Solid' (volume coupling).
   *-------------------------------------------------------------------------*/

  cs_syr_coupling_define("Solid",
                         NULL,                          /* boundary */
                         "box[0., 0., 0., 1., 1., 1.]", /* volume */
                         ' ',                           /* projection */
                         verbosity,
                         plot);

  /* By default, conservativity forcing flag is switched off (value 0)
     If one wants to switch on the conservativity forcing flag:

     cs_syr_coupling_set_conservativity(1);
  */

  /* Only for a volume coupling:
      By default, implicit treatment is done. You can switch to
      an explicit treatment by using the following function:

     cs_syr_coupling_set_explicit_treatment();
  */

}

/*----------------------------------------------------------------------------
 * Define couplings with other instances of Code_Saturne.
 *
 * This is done by calling the cs_sat_coupling_define() function for each
 * coupling to add.
 *
 * The arguments to cs_sat_coupling_define are:
 *   saturne_name          <-- matching Code_Saturne application name
 *   boundary_cpl_criteria <-- boundary face selection criteria for coupled
 *                             faces
 *   volume_cpl_criteria   <-- cell selection criteria for coupled cells
 *   boundary_sup_criteria <-- boundary face selection criteria for support
 *                             (not functional)
 *   volume_sup_criteria   <-- cell selection criteria for support
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
                         "3 or 4",
                         NULL,
                         NULL,
                         "all[]",
                         verbosity);

  /*-------------------------------------------------------------------------
   * Example 2: coupling with instance "SATURNE_03".
   *
   * - coupled faces of groups "coupled_faces"
   * - coupled cells (every cell overlapping the distant mesh)
   * - all cells available as location support
   *-------------------------------------------------------------------------*/

  cs_sat_coupling_define("SATURNE_03",
                         "coupled_faces",
                         "all[]",
                         NULL,
                         "all[]",
                         verbosity);
}

END_C_DECLS
