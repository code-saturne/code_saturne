/*============================================================================
 * Compressible models data
 *============================================================================*/

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
/*----------------------------------------------------------------------------*/

/*!
  \defgroup compressible Compressible models options

  @addtogroup compressible
  @{

  \struct cs_cf_model_t

  \brief Compressible model general options descriptor

  Members of these fluid properties are publicly accessible, to allow for
  concise syntax, as they are expected to be used in many places.

  \var  cs_cf_model_t::ieos
        indicator of equation of state
        -  CS_EOS_IDEAL_GAS: ideal gas with a constant adiabatic coefficient
        -  CS_EOS_STIFFENED_GAS: stiffened gas
        -  CS_EOS_GAS_MIX: mix of ideal gas
        -  CS_EOS_HOMOGENEOUS_TWO_PHASE: two-phase homogeneous model only,
           each phase follows a stiffened gas law.

  \var  cs_cf_model_t::ithvar
        indicator for thermodynamic variables initialization

  \var  cs_cf_model_t::psginf
        stiffened gas limit pressure (zero in perfect gas) for single phase
        model in Pa

  \var  cs_cf_model_t::gammasg
        stiffened gas polytropic coefficient (dimensionless) for single phase
        model

  \defgroup comp_homogeneous Homogeneous two-phase compressible model options

  @addtogroup comp_homogeneous
  @{

  \var  cs_cf_model_t::hgn_relax_eq_st
        source term step indicator for two-phase homogeneous model:
        - -1 disabled
        -  0 enabled

  @}

  @}

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* main compressible model structure */

static cs_cf_model_t  _cf_model =
{
  .ieos            = -1,
  .ithvar          =  10000,
  .psginf          = 0.,
  .gammasg         = 1.4,
  .hgn_relax_eq_st = -1
};

const cs_cf_model_t  *cs_glob_cf_model = &_cf_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_cf_model_get_pointers(int    **ieos,
                           int    **ithvar,
                           double **psginf,
                           double **gammasg,
                           int    **hgn_relax_eq_st);

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
 *   ieos             --> pointer to cs_glob_cf_model->ieos
 *   ithvar           --> pointer to cs_glob_cf_model->ithvar
 *   psginf           --> pointer to cs_glob_cf_model->psginf
 *   gammasg          --> pointer to cs_glob_cf_model->gammasg
 *   hgn_relax_eq_st  --> pointer to cs_glob_cf_model->hgn_relax_eq_st
 *----------------------------------------------------------------------------*/

void
cs_f_cf_model_get_pointers(int    **ieos,
                           int    **ithvar,
                           double **psginf,
                           double **gammasg,
                           int    **hgn_relax_eq_st)
{
  *ieos             = &(_cf_model.ieos);
  *ithvar           = &(_cf_model.ithvar);
  *psginf           = &(_cf_model.psginf);
  *gammasg          = &(_cf_model.gammasg);
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
