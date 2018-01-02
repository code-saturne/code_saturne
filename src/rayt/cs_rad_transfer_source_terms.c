/*============================================================================
 * Radiation solver source term integration for thermal scalar
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_parameters.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_rad_transfer_source_terms.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for fortran API
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
 * \brief Implicit and explicit radiative source terms for thermal scalar.
 *
 * \param[inout]  smbrs   work array for right hand side
 * \param[inout]  rovsdt  work array for unsteady term
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_source_terms(cs_real_t  smbrs[],
                             cs_real_t  rovsdt[])
{
  if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
      || cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    /* Implicit part   */
    cs_field_t *f_tsri = cs_field_by_name("rad_st_implicit");
    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {
      f_tsri->val[iel] = CS_MAX( -f_tsri->val[iel], 0.0);
      rovsdt[iel] += f_tsri->val[iel] * cs_glob_mesh_quantities->cell_vol[iel];
    }

    /* Explicit part   */
    cs_field_t *f_tsre = cs_field_by_name("rad_st");
    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {
      smbrs[iel] += f_tsre->val[iel] * cs_glob_mesh_quantities->cell_vol[iel];
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
