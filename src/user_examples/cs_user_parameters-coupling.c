/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

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
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

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
 * \file cs_user_parameters-coupling.c
 *
 * \brief Code coupling parameters examples.
 *
 * See \ref user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  /*-------------------------------------------------------------------------
   * Ensure coupled codes use the smallest of their defined time steps.
   *-------------------------------------------------------------------------*/

  /*! [coupling_ts] */
  cs_coupling_set_sync_flag(PLE_COUPLING_TS_MIN);
  /*! [coupling_ts] */

  /*-------------------------------------------------------------------------
   * Example time step multiplier for external couplings.
   *
   * The apparent time step for the current instance times (as viewed by
   * coupled codes) is equal to the true time step times this multiplier.
   *
   * When coupling with SYRTHES, it is recommended to use the same multiplier
   * here as for the thermal variable time step (this is not automated,
   * to allow for more advanced combinations if necessary, so the user
   * should ensure this when using a time step multiplier).
   *-------------------------------------------------------------------------*/
  {
    /*! [coupling_1] */
    cs_coupling_set_ts_multiplier(10.);
    /*! [coupling_1] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
