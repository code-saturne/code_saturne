/*============================================================================
 * Code couplings definition with SYRTHES and Code_Saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of Code_Saturne
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_coupling.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_coupling-saturne.c
 *
 * \brief Code couplings definition with SYRTHES and Code_Saturne.
 *
 * See \subpage user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with other instances of Code_Saturne.
 *
 * This is done by calling the \ref cs_sat_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_saturne_coupling(void)
{
  /*! [coupling_saturne_1] */
  {
    int  verbosity = 1;

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

  }
  /*! [coupling_saturne_1] */

  /*! [coupling_saturne_2] */
  {

    int  verbosity = 1;

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
  /*! [coupling_saturne_2] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
