/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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

#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_restart.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_time_control.h"
#include "cs_timer.h"
#include "cs_thermal_model.h"

#include "cs_combustion_model.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"

#include "cs_rad_transfer.h"
#include "cs_rad_transfer_fields.h"
#include "cs_rad_transfer_dir.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_options.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_rad_transfer_options.c */

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

  /* -> nrphas: for pulverized coal and fuel combustion:
   *            nrphas = 1 (gas) + number of classes (particles or droplets) */

  /* -> For pulverized coal and fuel combustion:   */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0)
    rt_params->nrphas += cs_glob_combustion_model->coal->nclacp;

  /* Default initializations
   * ----------------------- */

  /* ->  Restart computation (read restart) */

  rt_params->restart
    = (cs_restart_present()) ? 1 : 0;

  /* ->  Radiation solver call frequency */

  cs_time_control_init_by_time_step(&(rt_params->time_control),
                                    - 1,     /* nt_start */
                                    -1,      /* nt_end */
                                    1,       /* interval */
                                    true,    /* at start */
                                    false);  /* at end */

  /* ->  Explicit radiative source term computation mode
   *            = 0 => Semi-analytic (mandatory if transparent)
   *            = 1 => Conservative
   *            = 2 => Corrected semi-analytic (to be conservative)
   *     REMARK: if transparent, idiver = -1 automatically in raydom  */

  if (rt_params->atmo_model == CS_RAD_ATMO_3D_NONE)
    rt_params->idiver = 2;

  if (rt_params->imoadf == 1)
    rt_params->nwsgg = 8;

  else if (rt_params->imoadf == 2)
    rt_params->nwsgg = 50;

  if (rt_params->imfsck == 1)
    rt_params->nwsgg = 7;

  /* Add bands for Direct Solar, diFUse solar and InfraRed
   * and for solar: make the distinction between UV-visible (absorbed by O3)
   * and Solar IR (SIR) absobed by H2O
   * if activated */
  {
    if (rt_params->atmo_model
        != CS_RAD_ATMO_3D_NONE)
      rt_params->nwsgg = 0;

    /* Fields for atmospheric Direct Solar (DR) model
     * (SIR only if SUV is not activated) */
    if (rt_params->atmo_model
        & CS_RAD_ATMO_3D_DIRECT_SOLAR) {
      rt_params->atmo_dr_id = rt_params->nwsgg;
      rt_params->nwsgg++;
    }

    /* Fields for atmospheric Direct Solar (DR) model
     * (SUV band) */
    if (rt_params->atmo_model
        & CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND) {
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
   * Spectral number must be prescribed in usppmod
   * if (rt_params->imfsck == 2)
       rt_params->nwsgg = user prescribed in usppmod;
  */

  /* Coherency check with thermal model */

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Radiative module"),
                                "cs_glob_rad_transfer_params->type",
                                cs_glob_rad_transfer_params->type,
                                0, 3);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Radiative module"),
                                "cs_glob_rad_transfer_params->imodak",
                                cs_glob_rad_transfer_params->imodak,
                                0, 3);

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

  const char *imodak_value_str[]
    = {N_("0 (do not use Modak)"),
       N_("1 (Modak absorption coefficient)"),
       N_("2 (Magnussen, Kent and Honnery models")};

  const char *imoadf_value_str[]
    = {N_("0 (no AFD model)"),
       N_("1 (ADF model with 8 wavelength intervals)"),
       N_("2 (ADF model with 50 wavelength intervals)")};

  const char *imfsck_value_str[]
    = {N_("0 (no FSCK model)"),
       N_("1 (FSCK model activated)"),
       N_("2 (FSCK model activated with tabulated properties)")};

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
                  _("    imodak:        %s\n"),
                _(imodak_value_str[cs_glob_rad_transfer_params->imodak]));

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

  int iimlum = CS_MIN(2, CS_MAX(cs_glob_rad_transfer_params->verbosity, 0));

  cs_log_printf(CS_LOG_SETUP,
                  _("    verbosity:        %s\n"),
                _(iimlum_value_str[iimlum]));

  cs_log_printf(CS_LOG_SETUP,
                  _("    imoadf:        %s\n"),
                _(imoadf_value_str[cs_glob_rad_transfer_params->imoadf]));

  cs_log_printf(CS_LOG_SETUP,
                  _("    imfsck:        %s\n"),
                _(imfsck_value_str[cs_glob_rad_transfer_params->imfsck]));

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
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Direct solar O3 atmospheric 3D model on\n"
           "      band id = %d\n"),
         cs_glob_rad_transfer_params->atmo_dr_o3_id);

    if (  cs_glob_rad_transfer_params->atmo_model
        & CS_RAD_ATMO_3D_DIFFUSE_SOLAR)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Diffuse solar atmospheric 3D model on\n"
           "      band id = %d\n"),
         cs_glob_rad_transfer_params->atmo_df_id);

   if (  cs_glob_rad_transfer_params->atmo_model
        & CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Diffuse solar O3 atmospheric 3D model on\n"
           "      band id = %d\n"),
         cs_glob_rad_transfer_params->atmo_df_o3_id);

    if (  cs_glob_rad_transfer_params->atmo_model
        & CS_RAD_ATMO_3D_INFRARED)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Infra-red atmospheric 3D model on\n"
           "      band id = %d\n"),
         cs_glob_rad_transfer_params->atmo_ir_id);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
