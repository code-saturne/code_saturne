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
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_thermal_model.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"
#include "cs_atmo_aerosol_ssh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function initializes the external aerosol code
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_initialize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_initialize();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_finalize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_finalize();

  BFT_FREE(cs_glob_atmo_chemistry->aero_conc_file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with aerosol concentrations
 *        and numbers from the external aerosol code..
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_aero(cs_real_t* array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_get_aero(array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_gas(cs_real_t* array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_get_gas(array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_time_advance(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_time_advance();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
