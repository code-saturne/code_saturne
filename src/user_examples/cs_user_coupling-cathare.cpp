/*============================================================================
 * Code couplings definition with CATHARE and code_saturne.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_coupling.c
 *
 * \brief Code couplings definition with SYRTHES and code_saturne.
 *
 * See \ref user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with CATHARE code.
 *
 * This is done by calling the \ref cs_sys_coupling_add function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cathare_coupling(void)
{

  /*! [coupling_cathare_1] */
  {
    /* Add a coupling with cathare using only one phase */
    cs_sys_coupling_add("CAT-PIPE", // Name of cathare instance
                        1);         // Number of coupled phases

    /* Couple the inlet zone only */
    cs_sys_cpl_t *catcpl = cs_sys_coupling_by_name("CAT-PIPE");
    cs_zone_t *zone = cs_boundary_zone_by_name("inlet");

    cs_sys_coupling_add_cplbc(catcpl,               // pointer to coupling structure
                              CS_SYS_CPL_BC_INLET,  // Type of BC
                              zone,                 // coupled zone
                              "z < 0.1",            // homogenized 1D volume selection criteria
                              "TUBE1",              // name of cathare element
                              3,                    // Scalar cell just outside boundary
                              4,                    // Scalar cell just after the boundary
                              1);                   // Number of elements, 1 for inlet/outlet
  }
  /*! [coupling_cathare_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
