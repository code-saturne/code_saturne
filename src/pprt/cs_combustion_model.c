/*============================================================================
 * Combustion model parameters.
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

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_combustion_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_physical_model.c
        Specific physical models selection.
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

/*! Combustion model parameters structure */

cs_combustion_model_t
  _combustion_model = {.gas = {.iic = 0,
                               .xsoot = 0.,
                               .rosoot = 0.},
                       .coal = {.nclacp = 0},
                       .fuel = {.nclafu = 0},
                       .isoot = -1,
                       .ckabs0 = 0,
                       .xco2 = -1,
                       .xh2o = -1};

cs_combustion_model_t
  *cs_glob_combustion_model = &_combustion_model;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_combustion_model_get_pointers(int  **isoot);

void
cs_f_ppthch_get_pointers(int     **iic,
                         double  **wmole,
                         double  **wmolg,
                         double  **xco2,
                         double  **xh2o,
                         double  **ckabs1);

void
cs_f_coincl_get_pointers(double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot);

void
cs_f_cpincl_get_pointers(int     **ico2,
                         int     **ih2o,
                         int     **nclacp,
                         int     **ichcor,
                         double  **xashch,
                         double  **diam20,
                         double  **dia2mn,
                         double  **rho20,
                         double  **rho2mn,
                         double  **xmp0,
                         double  **xmasch);

void
cs_f_fuel_get_pointers(int     **nclafu);

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
 *   isoot  --> pointer to issot in combustion model
 *----------------------------------------------------------------------------*/

void
cs_f_combustion_model_get_pointers(int  **isoot)
{
  *isoot = &(cs_glob_combustion_model->isoot);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global physical model flags.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iic    --> pointer to rank of C in gas composition
 *   wmole  --> pointer to molar mass of elementary gas components
 *   wmolg  --> pointer to molar mass of global species
 *   xco2   --> pointer to molar coefficient of co2
 *   xh2o   --> pointer to molar coefficient of h2o
 *   ckabs1 --> pointer to absorption coefficient of gas mixture
 *----------------------------------------------------------------------------*/

void
cs_f_ppthch_get_pointers(int     **iic,
                         double  **wmole,
                         double  **wmolg,
                         double  **xco2,
                         double  **xh2o,
                         double  **ckabs1)
{
  *iic    = &(cs_glob_combustion_model->gas.iic);
  *wmole  = cs_glob_combustion_model->wmole;
  *wmolg  = cs_glob_combustion_model->gas.wmolg;
  *xco2   = &(cs_glob_combustion_model->xco2);
  *xh2o   = &(cs_glob_combustion_model->xh2o);
  *ckabs1 = &(cs_glob_combustion_model->ckabs0);
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
cs_f_coincl_get_pointers(double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot)
{
  *coefeg = &(cs_glob_combustion_model->gas.coefeg[0][0]);
  *compog = &(cs_glob_combustion_model->gas.compog[0][0]);
  *xsoot  = &(cs_glob_combustion_model->gas.xsoot);
  *rosoot = &(cs_glob_combustion_model->gas.rosoot);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ico2   --> pointer to cs_glob_combustion_model->ico2
 *   ih2o   --> pointer to cs_glob_combustion_model->ih2o
 *   nclacp --> pointer to cs_glob_combustion_model->coal.nclacp
 *   ichcor --> pointer to cs_glob_combustion_model->coal.ichcor
 *   xashch --> pointer to cs_glob_combustion_model->coal.xashch
 *   diam20 --> pointer to cs_glob_combustion_model->coal.diam20
 *   dia2mn --> pointer to cs_glob_combustion_model->coal.dia2mn
 *   rho20  --> pointer to cs_glob_combustion_model->coal.rho20
 *   rho2mn --> pointer to cs_glob_combustion_model->coal.rho2mn
 *   xmp0   --> pointer to cs_glob_combustion_model->coal.xmp0
 *   xmasch --> pointer to cs_glob_combustion_model->coal.xmash
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_get_pointers(int     **ico2,
                         int     **ih2o,
                         int     **nclacp,
                         int     **ichcor,
                         double  **xashch,
                         double  **diam20,
                         double  **dia2mn,
                         double  **rho20,
                         double  **rho2mn,
                         double  **xmp0,
                         double  **xmasch)
{
  *ico2   = &(cs_glob_combustion_model->ico2);
  *ih2o   = &(cs_glob_combustion_model->ih2o);
  *nclacp = &(cs_glob_combustion_model->coal.nclacp);
  *ichcor = cs_glob_combustion_model->coal.ichcor;
  *xashch = cs_glob_combustion_model->coal.xashch;
  *diam20 = cs_glob_combustion_model->coal.diam20;
  *dia2mn = cs_glob_combustion_model->coal.dia2mn;
  *rho20  = cs_glob_combustion_model->coal.rho20;
  *rho2mn = cs_glob_combustion_model->coal.rho2mn;
  *xmp0   = cs_glob_combustion_model->coal.xmp0;
  *xmasch = cs_glob_combustion_model->coal.xmasch;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the fuel combustion model (cs_fuel_incl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   nclafu --> pointer to cs_glob_combustion_model->fuel.nclafu
 *----------------------------------------------------------------------------*/

void
cs_f_fuel_get_pointers(int     **nclafu)
{
  *nclafu = &(cs_glob_combustion_model->fuel.nclafu);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
