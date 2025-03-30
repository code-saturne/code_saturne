/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_log.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_parameters_check.h"
#include "base/cs_restart.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_it.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_control.h"
#include "base/cs_timer.h"

#include "comb/cs_coal.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"

#include "rayt/cs_rad_transfer.h"
#include "rayt/cs_rad_transfer_dir.h"
#include "rayt/cs_rad_transfer_fields.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "rayt/cs_rad_transfer_options.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_rad_transfer_options.cpp */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Radiation solver initialization
 *
 *  1) Default initialization
 *  2) User parameter input reading
 *  3) Coherency controls against particular physics models
 *  4) User parameters controls
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_options(void)
{
  /* Pointer with shorter syntax ... */
  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  /* nrphas: for pulverized coal and fuel combustion:
   *         nrphas = 1 (gas) + number of classes (particles or droplets) */

  /* For pulverized coal combustion:   */

  if (cs_glob_coal_model != nullptr) {
    rt_params->nrphas += cs_glob_coal_model->nclacp;
  }

  /* Default initializations
   * ----------------------- */

  rt_params->restart = (cs_restart_present()) ? 1 : 0;

  /* Explicit radiative source term computation mode
   * - 0 => Semi-analytic (mandatory if transparent)
   * - 1 => Conservative
   * - 2 => Corrected semi-analytic (to be conservative)
   *
   * REMARK: if transparent, idiver is automatically reset to -1
   *         in _rad_transfer_solve. */

  if (   rt_params->atmo_model == CS_RAD_ATMO_3D_NONE
      && (rt_params->idiver == 0 || rt_params->idiver == 1)) {
    cs_parameters_error(CS_WARNING,
                        _("in Radiative module"),
                        "Using the semi-analytic (idiver = 0)\n"
                        "or conservative (idiver = 1) explicit source term\n"
                        "computation mode is not recommended.\n");
  }

  if (rt_params->imoadf == 1)
    rt_params->nwsgg = 8;

  else if (rt_params->imoadf == 2)
    rt_params->nwsgg = 50;

  if (rt_params->imfsck == 1)
    rt_params->nwsgg = 7;

  if (rt_params->imrcfsk == 1)
    rt_params->nwsgg = 10;

  /* Add bands for Direct Solar, diFUse solar and InfraRed
   * and for solar: make the distinction between UV-visible (absorbed by O3)
   * and Solar IR (SIR) absobed by H2O
   * if activated */
  {
    if (rt_params->atmo_model != CS_RAD_ATMO_3D_NONE)
      rt_params->nwsgg = 0;

    /* Fields for atmospheric Direct Solar (DR) model
     * (SIR only if SUV is not activated) */
    if (rt_params->atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR) {
      rt_params->atmo_dr_id = rt_params->nwsgg;
      rt_params->nwsgg++;
    }

    /* Fields for atmospheric Direct Solar (DR) model
     * (SUV band) */
    if (rt_params->atmo_model & CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND) {
      rt_params->atmo_dr_o3_id = rt_params->nwsgg;
      rt_params->nwsgg++;
    }

    /* Fields for atmospheric diFfuse Solar (DF) model
     * (SIR only if SUV is activated) */
    if (rt_params->atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR) {
      rt_params->atmo_df_id = rt_params->nwsgg;
      rt_params->nwsgg++;
    }

    /* Fields for atmospheric diFfuse Solar (DF) model
     * (SUV band) */
    if (rt_params->atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND) {
      rt_params->atmo_df_o3_id = rt_params->nwsgg;
      rt_params->nwsgg++;
    }

    /* Fields for atmospheric infrared absorption model */
    if (rt_params->atmo_model & CS_RAD_ATMO_3D_INFRARED) {
      rt_params->atmo_ir_id = rt_params->nwsgg;
      rt_params->nwsgg++;
    }
  }

  /* In case of imfsck == 2, spectral radiative properties depend on
   * user provided radiation library coupled with flamelet library.
   *
   * Spectral number must be prescribed in cs_user_model
   * if (rt_params->imfsck == 2)
       rt_params->nwsgg = user prescribed in cs_user_model;
  */

  /* Coherency check with thermal model */

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Radiative module"),
                                "cs_glob_rad_transfer_params->type",
                                cs_glob_rad_transfer_params->type,
                                0, 3);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Radiative module"),
                                "cs_glob_rad_transfer_params->imgrey",
                                cs_glob_rad_transfer_params->imgrey,
                                0, 3);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Radiative module"),
                                "cs_glob_rad_transfer_params->imfsck",
                                cs_glob_rad_transfer_params->imfsck,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Radiative module"),
                                "cs_glob_rad_transfer_params->imrcfsk",
                                cs_glob_rad_transfer_params->imrcfsk,
                                0, 2);

  /* Verifications */

  if (rt_params->type > CS_RAD_TRANSFER_NONE) {

    /* --> NFREQR */

    if (rt_params->time_control.interval_nt <= 0)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("in Radiative module"),
         _("Thermal model resolution frequency"
           " (cs_glob_rad_transfer_params->time_control.interval_nt)\n"
           "must be > 0, and not %d.\n"),
         rt_params->time_control.interval_nt);

    /* --> i_quadrature     */

    if (rt_params->type == CS_RAD_TRANSFER_DOM) {
      cs_parameters_is_in_range_int
        (CS_ABORT_DELAYED,
         _("in Radiative module"),
         _("The quadrature type number"
           " (cs_glob_rad_transfer_params->i_quadrature)"),
         rt_params->i_quadrature,
         1, 11);
    }

    /* --> NDIREC */

    if (   rt_params->type == CS_RAD_TRANSFER_DOM
        && rt_params->i_quadrature == 6) {
      if (rt_params->ndirec < 2)
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("in Radiative module"),
           _("Tn quadrature parameter n must be > 1, and not %d.\n"),
           rt_params->ndirec);
    }

    /* --> IDIVER */

    cs_parameters_is_in_range_int
      (CS_ABORT_DELAYED,
       _("in Radiative module"),
       _("Computation mode parameter"
         " (cs_glob_rad_transfer_params->idiver"),
       rt_params->idiver,
       -1, 3);

    cs_parameters_error_barrier();
  }

  else
    return;

  /* Quadrature initialization */

  cs_rad_transfer_dir();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log radiative module options in setup file.
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_log_setup(void)
{
  if (cs_glob_rad_transfer_params->type <= CS_RAD_TRANSFER_NONE)
    return;

  /* Now add Lagrangian setup info */

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Radiative thermal transfer options\n"
                  "----------------------------------\n\n"));

  cs_log_printf
    (CS_LOG_SETUP,
     _("  Continuous phase:\n"
       "    type:          %s\n"),
     cs_rad_transfer_model_name[cs_glob_rad_transfer_params->type]);

  const char *restart_value_str[] = {N_("0 (no restart)"),
                                     N_("1 (restart)")};

  cs_log_printf(CS_LOG_SETUP,
                _("    restart:       %s\n"),
                _(restart_value_str[cs_glob_rad_transfer_params->restart]));

  char buf[128];
  cs_time_control_get_description(&(cs_glob_rad_transfer_params->time_control),
                                  buf,
                                  128);

  cs_log_printf(CS_LOG_SETUP,
                _("    time control:      %s\n"), buf);

  if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_DOM) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("    i_quadrature:  %s\n"),
       _(cs_rad_transfer_quadrature_name
           [cs_glob_rad_transfer_params->i_quadrature]));
    if (cs_glob_rad_transfer_params->i_quadrature == 6)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    ndirec:       %d\n"),
         cs_glob_rad_transfer_params->ndirec);
  }

  const char *imgrey_value_str[]
    = {N_("0 (do not use the grey body model)"),
       N_("1 (Classical grey body model - Modak)"),
       N_("2 (New grey body model - Magnussen, Kent and Honnery models")};

  const char *imoadf_value_str[]
    = {N_("0 (no AFD model)"),
       N_("1 (ADF model with 8 wavelength intervals)"),
       N_("2 (ADF model with 50 wavelength intervals)")};

  const char *imfsck_value_str[]
    = {N_("0 (no FSCK model)"),
       N_("1 (FSCK model activated)"),
       N_("2 (FSCK model activated with tabulated properties)")};

  const char *imrcfsk_value_str[]
    = {N_("0 (no RCFSK model)"),
       N_("1 (RCFSK model activated)")};

  const char *idiver_value_str[]
    = {N_("-1 (no renormalization)"),
       N_("0 (semi-analytic radiative S.T. calculation;\n"
          "                   "
          "   compulsory with transparent media)"),
       N_("1 (conservative radiative S.T. calculation)"),
       N_("2 (semi-analytic radiative S.T. calculation,\n"
          "                   "
          "   corrected for global conservation)")};

  cs_log_printf(CS_LOG_SETUP,
                _("    idiver:        %s\n"),
                _(idiver_value_str[cs_glob_rad_transfer_params->idiver+1]));

  cs_log_printf(CS_LOG_SETUP,
                _("    imgrey:        %s\n"),
                _(imgrey_value_str[cs_glob_rad_transfer_params->imgrey]));

  const char *iimpar_value_str[]
    = {N_("0 (do not log wall temperature)"),
       N_("1 (standard wall temperature log)"),
       N_("2 (detailed wall temperature compute log)")};

  cs_log_printf(CS_LOG_SETUP,
                _("    iimpar:        %s\n"),
                _(iimpar_value_str[cs_glob_rad_transfer_params->iimpar]));

  const char *iimlum_value_str[]
    = {N_("0 (no solver logging)"),
       N_("1 (standard solver log)"),
       N_("2 (detailed solver logging)")};

  int iimlum = cs::min(2, cs::max(cs_glob_rad_transfer_params->verbosity, 0));

  cs_log_printf(CS_LOG_SETUP,
                  _("    verbosity:        %s\n"),
                _(iimlum_value_str[iimlum]));

  cs_log_printf(CS_LOG_SETUP,
                _("    imoadf:        %s\n"),
                _(imoadf_value_str[cs_glob_rad_transfer_params->imoadf]));

  cs_log_printf(CS_LOG_SETUP,
                _("    imfsck:        %s\n"),
                _(imfsck_value_str[cs_glob_rad_transfer_params->imfsck]));

  cs_log_printf(CS_LOG_SETUP,
                _("    imrcfsk:       %s\n"),
                _(imrcfsk_value_str[cs_glob_rad_transfer_params->imrcfsk]));

  if (cs_glob_rad_transfer_params->atmo_model) {

    if (  cs_glob_rad_transfer_params->atmo_model
        & CS_RAD_ATMO_3D_DIRECT_SOLAR)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Direct solar atmospheric 3D model on\n"
           "      band id = %d\n"),
         cs_glob_rad_transfer_params->atmo_dr_id);

    if (  cs_glob_rad_transfer_params->atmo_model
        & CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND)
      cs_log_printf(CS_LOG_SETUP,
                    _("    Direct solar O3 atmospheric 3D model on\n"
                      "      band id = %d\n"),
                    cs_glob_rad_transfer_params->atmo_dr_o3_id);

    if (cs_glob_rad_transfer_params->atmo_model & CS_RAD_ATMO_3D_DIFFUSE_SOLAR)
      cs_log_printf(CS_LOG_SETUP,
                    _("    Diffuse solar atmospheric 3D model on\n"
                      "      band id = %d\n"),
                    cs_glob_rad_transfer_params->atmo_df_id);

    if (cs_glob_rad_transfer_params->atmo_model
        & CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND)
      cs_log_printf(CS_LOG_SETUP,
                    _("    Diffuse solar O3 atmospheric 3D model on\n"
                      "      band id = %d\n"),
                    cs_glob_rad_transfer_params->atmo_df_o3_id);

    if (cs_glob_rad_transfer_params->atmo_model & CS_RAD_ATMO_3D_INFRARED)
      cs_log_printf(CS_LOG_SETUP,
                    _("    Infra-red atmospheric 3D model on\n"
                      "      band id = %d\n"),
                    cs_glob_rad_transfer_params->atmo_ir_id);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
