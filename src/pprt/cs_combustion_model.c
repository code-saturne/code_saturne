/*============================================================================
 * Combustion model parameters.
 *============================================================================*/

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"
#include "cs_base.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_combustion_gas.h"
#include "cs_coal.h"

#include "cs_combustion_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_model.c
        Combustion  model selection parameters.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_coal_model_map(void);

void
cs_f_ppthch_get_pointers(int     **ngaze,
                         int     **ngazg,
                         int     **nato,
                         int     **nrgaz,
                         int     **iic,
                         int     **npo,
                         double  **wmole,
                         double  **wmolg,
                         double  **xco2,
                         double  **xh2o,
                         double  **ckabs1,
                         double  **fs,
                         double  **th,
                         double  **cpgazg);

void
cs_f_coincl_get_pointers(int     **isoot,
                         bool    **use_janaf,
                         double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot,
                         double  **hinfue,
                         double  **hinoxy,
                         double  **pcigas,
                         double  **tinfue,
                         double  **tinoxy);

void
cs_f_ppcpfu_get_pointers(double  **oxyo2,
                         double  **oxyn2,
                         double  **oxyh2o,
                         double  **oxyco2);

void
cs_f_combustion_model_get_pointers(double  **srrom);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global physical model flags.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ngaze  --> pointer to number of elementary species
 *   ngazg  --> pointer to number of global species
 *   nato   --> pointer to number of atomic species
 *   iic    --> pointer to rank of C in gas composition
 *   wmole  --> pointer to molar mass of elementary gas components
 *   wmolg  --> pointer to molar mass of global species
 *   xco2   --> pointer to molar coefficient of co2
 *   xh2o   --> pointer to molar coefficient of h2o
 *   ckabs1 --> pointer to absorption coefficient of gas mixture
 *   fs     --> pointer to mixing rate at the stoichiometry
 *----------------------------------------------------------------------------*/

void
cs_f_ppthch_get_pointers(int     **ngaze,
                         int     **ngazg,
                         int     **nato,
                         int     **nrgaz,
                         int     **iic,
                         int     **npo,
                         double  **wmole,
                         double  **wmolg,
                         double  **xco2,
                         double  **xh2o,
                         double  **ckabs1,
                         double  **fs,
                         double  **th,
                         double  **cpgazg)
{
  *ckabs1 = NULL;

  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *ngaze  = &(cm->n_gas_el_comp);
    *ngazg  = &(cm->n_gas_species);
    *nato   = &(cm->n_atomic_species);
    *nrgaz  = &(cm->n_reactions);
    *npo    = &(cm->n_tab_points);
    *iic    = &(cm->iic);

    *wmole  = cm->wmole;
    *wmolg  = cm->wmolg;
    *xco2   = &(cm->xco2);
    *xh2o   = &(cm->xh2o);
    *fs     = cm->fs;
    *th     = cm->th;
    *cpgazg = (double *)cm->cpgazg;

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *ngaze  = &(cm->n_gas_el_comp);
    *ngazg  = &(cm->n_gas_species);
    *nato   = &(cm->n_atomic_species);
    *nrgaz  = &(cm->n_reactions);
    *npo    = &(cm->n_tab_points);

    *wmole  = cm->wmole;
    *xco2   = &(cm->xco2);
    *xh2o   = &(cm->xh2o);
    *ckabs1 = &(cm->ckabs0);
    *th     = cm->th;

  }
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global physical model (coincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   coefeg --> pointer to conversion coefficients
 *   compog --> pointer to conversion coefficients
 *----------------------------------------------------------------------------*/

void
cs_f_coincl_get_pointers(int     **isoot,
                         bool    **use_janaf,
                         double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot,
                         double  **hinfue,
                         double  **hinoxy,
                         double  **pcigas,
                         double  **tinfue,
                         double  **tinoxy)
{
  *isoot  = NULL;
  *use_janaf = NULL;
  *coefeg = NULL;
  *compog = NULL;
  *xsoot  = NULL;
  *rosoot = NULL;
  *hinfue = NULL;
  *tinfue = NULL;
  *hinoxy = NULL;
  *tinoxy = NULL;
  *pcigas = NULL;

  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *isoot  = &(cm->isoot);
    *use_janaf  = &(cm->use_janaf);
    *coefeg = &(cm->coefeg[0][0]);
    *compog = &(cm->compog[0][0]);
    *xsoot  = &(cm->xsoot);
    *rosoot = &(cm->rosoot);
    *hinfue = &(cm->hinfue);
    *tinfue = &(cm->tinfue);
    *hinoxy = &(cm->hinoxy);
    *tinoxy = &(cm->tinoxy);
    *pcigas = &(cm->pcigas);

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *pcigas = &(cm->pcigas);

  }
}

/*----------------------------------------------------------------------------
 * Get pointers to members of combustion model (ppcpfu).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ieqco2 --> pointer to cm->ieqco2
 *   oxyo2  --> pointer to cm->oxyo2
 *   oxyn2  --> pointer to cm->oxyn2
 *   oxyh2o --> pointer to cm->oxyh2o
 *   oxyco2 --> pointer to cm->oxyco2
 *----------------------------------------------------------------------------*/

void
cs_f_ppcpfu_get_pointers(double  **oxyo2,
                         double  **oxyn2,
                         double  **oxyh2o,
                         double  **oxyco2)
{
  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *oxyo2 = NULL;
    *oxyn2 =  cm->oxyn2;
    *oxyh2o = cm->oxyh2o;
    *oxyco2 = cm->oxyco2;

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *oxyo2 =  cm->oxyo2;
    *oxyn2 =  cm->oxyn2;
    *oxyh2o = cm->oxyh2o;
    *oxyco2 = cm->oxyco2;

  }
}

/*----------------------------------------------------------------------------
 * Get pointers to generic physical model pointers (pincl)
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_combustion_model_get_pointers(double  **srrom)
{
  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *srrom = &(cm->srrom);

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *srrom = &(cm->srrom);

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
