/*============================================================================
 * Code couplings definition with SYRTHES and Code_Saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of Code_Saturne
 *============================================================================*/

/* VERS */

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
 * \file cs_user_coupling.c
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
 * \brief Define global options for couplings.
 *
 * These options allow defining the time step synchronization policy,
 * as well as a time step multiplier.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_coupling(void)
{
  /*! [coupling_1] */
  {
    /*-------------------------------------------------------------------------
     * Example for time step multiplier for external couplings.
     *
     * The apparent time step for the current instance times (as viewed by
     * coupled codes) is equal to the true time step times this multiplier.
     *
     * When coupling with SYRTHES, it is recommended to use the same multiplier
     * here as for the thermal variable time step (this is not automated,
     * to allow for more advanced combinations if necessary, so the user
     * should ensure this when using a time step multiplier).
     *-------------------------------------------------------------------------*/

    cs_coupling_set_ts_multiplier(10.);

  }
  /*! [coupling_1] */
}

END_C_DECLS
