/*============================================================================
 * Compressible models data
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cf_model.c
        Compressible models data.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* main compressible model structure */

static cs_cf_model_t  _cf_model =
{
  .ithvar          =  10000,
  .hgn_relax_eq_st = -1
};

const cs_cf_model_t  *cs_glob_cf_model = &_cf_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_cf_model_get_pointers(int  **ithvar,
                           int  **hgn_relax_eq_st);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compressible model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ithvar           --> pointer to cs_glob_cf_model->ithvar
 *   hgn_relax_eq_st  --> pointer to cs_glob_cf_model->hgn_relax_eq_st
 *----------------------------------------------------------------------------*/

void
cs_f_cf_model_get_pointers(int  **ithvar,
                           int  **hgn_relax_eq_st)
{
  *ithvar           = &(_cf_model.ithvar);
  *hgn_relax_eq_st  = &(_cf_model.hgn_relax_eq_st);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to compressible model global structure cs_glob_cf_model
 */
/*----------------------------------------------------------------------------*/

cs_cf_model_t *
cs_get_glob_cf_model(void)
{
  return &_cf_model;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
