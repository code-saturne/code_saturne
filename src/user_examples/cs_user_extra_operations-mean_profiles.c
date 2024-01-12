/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
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
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  User Header
 *----------------------------------------------------------------------------*/

#include "cs_user_profile.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.

 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_initialize(cs_domain_t *domain)
{
  CS_UNUSED(domain);

  /*![Initialize]*/

  /* Initialize a  mean temperature profile over z */
  cs_real_t v_dir[3] = { 0.0, 0.0, 1.0 };

#if defined(HAVE_MEDCOUPLING)
  user_profile_t *profile
    = user_create_profile("T_vertical_profile", /* name */
                          "temperature",  /* field*/
                          "all[]",        /* cell selection */
                          v_dir,          /* profile direction */
                          10,             /* number of layers */
                          "PARABOLIC",    /* progression law */
                          1.5,            /* geometric progression */
                          "MASS",         /* mass, volume or no: weight */
                          "MEDCOUPLING"); /* method used to intersect volume */

#else
  user_profile_t *profile
    = user_create_profile("T_vertical_profile", /* name */
                          "temperature",   /* field*/
                          "all[]",         /* cell selection */
                          v_dir,           /* profile direction */
                          10,              /* number of layers */
                          "CONSTANT",      /* progression law */
                          1.0,             /* geometric progression */
                          "MASS",          /* mass, volume or no: weight */
                          "STL");          /* method used to intersect volume */

#endif

  /* Calculate once for each cell percent lying in each layers */

  user_compute_cell_volume_per_layer(profile);
  /*![Initialize]*/
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t *domain)
{
  CS_UNUSED(domain);

  /* Mean profile calculation and results output */
   /*![generate]*/
  user_profile_t *profile = user_profile_get_by_name("T_vertical_profile");
  user_profile_compute(profile);

  user_profile_output(profile, /* pointer to profile structure */
                      15);     /* interval (time step) */
  /*![generate]*/
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of the calculation.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.

 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_finalize(cs_domain_t *domain)
{
  CS_UNUSED(domain);

  /*![finalize]*/
  /* Mean profile calculation and output at last time step */

  int interval = 1;  /* Time step interval between outputs */

  user_profiles_compute_all();
  user_profiles_output_all(interval);
  user_profiles_histogram_ot_output_all(interval);

  /* Free profiles memory */

  user_free_profiles();

  /*![finalize]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
