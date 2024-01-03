/*============================================================================
 * Combustion model parameters.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

cs_combustion_gas_model_t _glob_gas_model = {.iic = 0,
                                             .hinfue = 0,
                                             .xsoot = 0.,
                                             .rosoot = 0.};

/*! Combustion model parameters structure */

cs_combustion_model_t
  _combustion_model = {.gas = &_glob_gas_model,
                       .coal = NULL,
                       .n_gas_el_comp = 0,
                       .n_gas_species = 0,
                       .n_atomic_species = 0,
                       .n_reactions = 0,
                       .idrift = 0,
                       .ieqco2 = 0,
                       .ieqnox = 0,
                       .isoot = -1,
                       .ckabs0 = 0,
                       .diftl0 = -1,
                       .xco2 = -1,
                       .xh2o = -1,
                       .hinoxy = 0};

cs_combustion_model_t
  *cs_glob_combustion_model = &_combustion_model;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_coal_model_map(void);

void
cs_f_combustion_model_get_pointers(int  **isoot);

void
cs_f_ppthch_get_pointers(int     **ngaze,
                         int     **ngazg,
                         int     **nato,
                         int     **nrgaz,
                         int     **iic,
                         int     **npo,
                         double  **wmole,
                         double  **wmolg,
                         double  **diftl0,
                         double  **xco2,
                         double  **xh2o,
                         double  **ckabs1,
                         double  **fs,
                         double  **th,
                         double  **cpgazg);

void
cs_f_coincl_get_pointers(double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot,
                         double  **hinfue,
                         double  **hinoxy,
                         double  **pcigas,
                         double  **tinfue,
                         double  **tinoxy);

void
cs_f_cpincl_comb_get_pointers(int     **ico2,
                              int     **ih2o,
                              int     **io2,
                              int     **in2);

void
cs_f_cpincl_coal_get_pointers(int     **ncharb,
                              int     **nclacp,
                              int     **nclpch,
                              int     **ich,
                              int     **ick,
                              int     **iash,
                              int     **iwat,
                              double  **ehsoli,
                              double  **wmols,
                              double  **eh0sol,
                              int     **ichcor,
                              double  **cp2ch,
                              double  **xashch,
                              double  **xwatch,
                              double  **diam20,
                              double  **dia2mn,
                              double  **rho20,
                              double  **rho2mn,
                              double  **xmp0,
                              double  **xmasch);

void
cs_f_ppcpfu_get_pointers(int     **idrift,
                         int     **ieqco2,
                         int     **ieqnox,
                         double  **oxyo2,
                         double  **oxyn2,
                         double  **oxyh2o,
                         double  **oxyco2);

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
 *   isoot  --> pointer to isoot in combustion model
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
 *   ngaze  --> pointer to number of elementary species
 *   ngazg  --> pointer to number of global species
 *   nato   --> pointer to number of atomic species
 *   iic    --> pointer to rank of C in gas composition
 *   wmole  --> pointer to molar mass of elementary gas components
 *   wmolg  --> pointer to molar mass of global species
 *   diftl0 --> pointer to molecular diffusivity for the enthalpy.
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
                         double  **diftl0,
                         double  **xco2,
                         double  **xh2o,
                         double  **ckabs1,
                         double  **fs,
                         double  **th,
                         double  **cpgazg)
{
  *ngaze  = &(cs_glob_combustion_model->n_gas_el_comp);
  *ngazg  = &(cs_glob_combustion_model->n_gas_species);
  *nato   = &(cs_glob_combustion_model->n_atomic_species);
  *nrgaz  = &(cs_glob_combustion_model->n_reactions);
  *npo    = &(cs_glob_combustion_model->n_tab_points);
  *iic    = &(cs_glob_combustion_model->gas->iic);

  *wmole  = cs_glob_combustion_model->wmole;
  *wmolg  = cs_glob_combustion_model->gas->wmolg;
  *diftl0 = &(cs_glob_combustion_model->diftl0);
  *xco2   = &(cs_glob_combustion_model->xco2);
  *xh2o   = &(cs_glob_combustion_model->xh2o);
  *ckabs1 = &(cs_glob_combustion_model->ckabs0);
  *fs     = cs_glob_combustion_model->gas->fs;
  *th     = cs_glob_combustion_model->th;
  *cpgazg = (double *)cs_glob_combustion_model->gas->cpgazg;
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
                         double  **rosoot,
                         double  **hinfue,
                         double  **hinoxy,
                         double  **pcigas,
                         double  **tinfue,
                         double  **tinoxy)
{
  *coefeg = &(cs_glob_combustion_model->gas->coefeg[0][0]);
  *compog = &(cs_glob_combustion_model->gas->compog[0][0]);
  *xsoot  = &(cs_glob_combustion_model->gas->xsoot);
  *rosoot = &(cs_glob_combustion_model->gas->rosoot);
  *hinfue = &(cs_glob_combustion_model->gas->hinfue);
  *tinfue = &(cs_glob_combustion_model->gas->tinfue);
  *hinoxy = &(cs_glob_combustion_model->hinoxy);
  *tinoxy = &(cs_glob_combustion_model->tinoxy);
  *pcigas = &(cs_glob_combustion_model->pcigas);
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
 *   io2    --> pointer to cs_glob_combustion_model->io2
 *   in2    --> pointer to cs_glob_combustion_model->in2
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_comb_get_pointers(int     **ico2,
                              int     **ih2o,
                              int     **io2,
                              int     **in2)
{
  *io2   = &(cs_glob_combustion_model->io2);
  *in2   = &(cs_glob_combustion_model->in2);
  *ico2   = &(cs_glob_combustion_model->ico2);
  *ih2o   = &(cs_glob_combustion_model->ih2o);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ncharb --> pointer to cs_glob_combustion_model->coal->n_coals
 *   nclacp --> pointer to cs_glob_combustion_model->coal->nclacp
 *   nclpch --> pointer to cs_glob_combustion_model->coal->n_classes_per_coal
 *   ich    --> pointer to cs_glob_combustion_model->coal->ich
 *   ick    --> pointer to cs_glob_combustion_model->coal->ick
 *   iash   --> pointer to cs_glob_combustion_model->coal->iash
 *   iwat   --> pointer to cs_glob_combustion_model->coal->iwat
 *   ehsoli --> pointer to cs_glob_combustion_model->coal->ehsoli
 *   wmols  --> pointer to cs_glob_combustion_model->coal->wmols
 *   eh0sol --> pointer to cs_glob_combustion_model->coal->eh0sol
 *   ichcor --> pointer to cs_glob_combustion_model->coal->ichcor
 *   cp2ch  --> pointer to cs_glob_combustion_model->coal->cp2ch
 *   xashch --> pointer to cs_glob_combustion_model->coal->xashch
 *   diam20 --> pointer to cs_glob_combustion_model->coal->diam20
 *   dia2mn --> pointer to cs_glob_combustion_model->coal->dia2mn
 *   rho20  --> pointer to cs_glob_combustion_model->coal->rho20
 *   rho2mn --> pointer to cs_glob_combustion_model->coal->rho2mn
 *   xmp0   --> pointer to cs_glob_combustion_model->coal->xmp0
 *   xmasch --> pointer to cs_glob_combustion_model->coal->xmash
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_coal_get_pointers(int     **ncharb,
                              int     **nclacp,
                              int     **nclpch,
                              int     **ich,
                              int     **ick,
                              int     **iash,
                              int     **iwat,
                              double  **ehsoli,
                              double  **wmols,
                              double  **eh0sol,
                              int     **ichcor,
                              double  **cp2ch,
                              double  **xashch,
                              double  **xwatch,
                              double  **diam20,
                              double  **dia2mn,
                              double  **rho20,
                              double  **rho2mn,
                              double  **xmp0,
                              double  **xmasch)
{
  cs_coal_model_t *coal = cs_glob_combustion_model->coal;

  *ncharb = &(coal->n_coals);
  *nclacp = &(coal->nclacp);
  *nclpch = coal->n_classes_per_coal;

  *ich    = coal->ich;
  *ick    = coal->ick;
  *iash   = coal->iash;
  *iwat   = coal->iwat;
  *ehsoli = (double *)coal->ehsoli;
  *wmols  = coal->wmols;
  *eh0sol = coal->eh0sol;

  *ichcor = coal->ichcor;
  *cp2ch  = coal->cp2ch;
  *xashch = coal->xashch;
  *xwatch = coal->xwatch;
  *diam20 = coal->diam20;
  *dia2mn = coal->dia2mn;
  *rho20  = coal->rho20;
  *rho2mn = coal->rho2mn;
  *xmp0   = coal->xmp0;
  *xmasch = coal->xmasch;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of combustion model (ppcpfu).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   idrift --> pointer to cs_glob_combustion_model->idrift
 *   ieqco2 --> pointer to cs_glob_combustion_model->ieqco2
 *   ieqnox --> pointer to cs_glob_combustion_model->ieqnox
 *   oxyo2  --> pointer to cs_glob_combustion_model->oxyo2
 *   oxyn2  --> pointer to cs_glob_combustion_model->oxyn2
 *   oxyh2o --> pointer to cs_glob_combustion_model->oxyh2o
 *   oxyco2 --> pointer to cs_glob_combustion_model->oxyco2
 *----------------------------------------------------------------------------*/

void
cs_f_ppcpfu_get_pointers(int     **idrift,
                         int     **ieqco2,
                         int     **ieqnox,
                         double  **oxyo2,
                         double  **oxyn2,
                         double  **oxyh2o,
                         double  **oxyco2)
{
  *idrift = &(cs_glob_combustion_model->idrift);
  *ieqco2 = &(cs_glob_combustion_model->ieqco2);
  *ieqnox = &(cs_glob_combustion_model->ieqnox);
  *oxyo2 =  cs_glob_combustion_model->oxyo2;
  *oxyn2 =  cs_glob_combustion_model->oxyn2;
  *oxyh2o = cs_glob_combustion_model->oxyh2o;
  *oxyco2 = cs_glob_combustion_model->oxyco2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize combustion model based on active physical models
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_initialize(void)
{
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0) {
    if (_combustion_model.coal == NULL) {
      _combustion_model.coal = cs_coal_model_create();
      cs_f_coal_model_map();
    }

    cs_base_at_finalize(cs_combustion_finalize);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize combustion model based on active physical models
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_finalize(void)
{
  if (_combustion_model.coal != NULL)
    cs_coal_model_destroy(&_combustion_model.coal);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the combustion module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_log_setup(void)
{
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_LW] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0) {

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "Combustion module options\n"
                    "-------------------------\n\n"));

    /* Gas combistion only */

    if (cs_glob_combustion_model->coal == NULL) {

      switch(cs_glob_combustion_model->isoot) {
      case -1:
        /* No Soot model */
        cs_log_printf(CS_LOG_SETUP,
                      _("    isoot:    -1 (No Soot model)\n\n"));
        break;
      case 0:
        /* constant fraction of product Xsoot */
        cs_log_printf(CS_LOG_SETUP,
                      _("    isoot:     0 (Constant soot yield)\n\n"));
        cs_log_printf(CS_LOG_SETUP,
                      _("  Parameters for the soot model:\n"
                        "    xsoot:  %14.5e (Fraction of product - Used only\n"
                        "            if the soot yield is not defined in the\n"
                        "            thermochemistry data file)\n"
                        "    rosoot: %14.5e (Soot density)\n\n"),
                      cs_glob_combustion_model->gas->xsoot,
                      cs_glob_combustion_model->gas->rosoot);
        break;
      case 1:
        /* 2 equations model of Moss et al. */
        cs_log_printf
          (CS_LOG_SETUP,
           _("    isoot:     1 (2 equations model of Moss et al.)\n\n"));
        cs_log_printf(CS_LOG_SETUP,
                      _("  Parameter for the soot model:\n"
                        "    rosoot: %14.5e (Soot density)\n\n"),
                      cs_glob_combustion_model->gas->rosoot);
        break;
      default:
        break;
      }

      const char *ipthrm_value_str[] = {N_("0 (no mean pressure computation)"),
                                        N_("1 (mean pressure computation)")};
      cs_log_printf(CS_LOG_SETUP,
                    _("    ipthrm:    %s\n\n"),
                    _(ipthrm_value_str[cs_glob_fluid_properties->ipthrm]));

    }

    cs_log_printf(CS_LOG_SETUP,
                  _("  Molecular diffusivity for the Enthalpy:\n"
                    "    diftl0: %14.5e\n"),
                    cs_glob_combustion_model->diftl0);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
