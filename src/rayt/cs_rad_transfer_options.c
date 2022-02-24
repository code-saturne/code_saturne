/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
#include "cs_gui_output.h"

#include "cs_combustion_model.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"

#include "cs_gui_radiative_transfer.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_property_fields.h"
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
    rt_params->nrphas = 1 + cs_glob_combustion_model->coal.nclacp;
  else if (cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
    rt_params->nrphas = 1 + cs_glob_combustion_model->fuel.nclafu;
  else
    rt_params->nrphas = 1;

  int nmodpp = 0;
  for (int ipp = 2; ipp < CS_N_PHYSICAL_MODEL_TYPES; ipp++) {
    if (cs_glob_physical_model_flag[ipp] != -1)
      nmodpp++;
  }

  /* ====================================================================
   * 1. Default initializations
   *    ^^^^^^^
   * ==================================================================== */

  /* ->  Absorption coefficent computation
   *      IMODAK = 0 : without modak
   *               1 : with modak    */

  rt_params->imodak       = 0;

  /* ->  IMOADF = 0 ADF models are not used
   *            = 1 ADF Model with wavelength interval of a 8
   *            = 2 ADF Model with wavelength interval of a 50*/

  rt_params->imoadf       = 0;

  /* ->  IMFSCK = 0 switch off FSCK model
   *            = 1 switch on FSCK model */

  rt_params->imfsck       = 0;

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

  /* ->  Quadrature number and Tn parameter */

  rt_params->i_quadrature = 1;
  rt_params->ndirec       = 3;

  /* ->  Cell percentage for which it is admitted that optical length
   *     is greater than one for P-1 model*/

  rt_params->xnp1mx       = 10.0;

  /* ->  Explicit radiative source term computation mode
   *            = 0 => Semi-analytic (mandatory if transparent)
   *            = 1 => Conservative
   *            = 2 => Corrected semi-analytic (to be conservative)
   *     REMARK: if transparent, idiver = -1 automatically in raydom  */

  rt_params->idiver       = 2;

  /* -> Wall temperature verbosity */

  rt_params->iimpar       = 1;

  /* Radiance resolution verbosity */

  rt_params->verbosity = 0;

  /* -> Number of iterations used to solve the ETR.
   *    Must at least be one to make the standard models work. */

  /* User parameters  */

  cs_gui_radiative_transfer_parameters();
  cs_user_radiative_transfer_parameters();

  if (rt_params->imoadf == 1)
    rt_params->nwsgg = 8;

  else if (rt_params->imoadf == 2)
    rt_params->nwsgg = 50;

  if (rt_params->imfsck == 1)
    rt_params->nwsgg = 7;

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
                                0, 2);

  if (   rt_params->type == CS_RAD_TRANSFER_DOM
      || rt_params->type == CS_RAD_TRANSFER_P1)
    cs_parameters_is_in_range_int
      (CS_ABORT_DELAYED,
       _("in Radiative module"),
       _("Thermal model option (cs_glob_thermal model->itherm)"),
       cs_glob_thermal_model->itherm,
       CS_THERMAL_MODEL_TEMPERATURE, CS_THERMAL_MODEL_TOTAL_ENERGY);

  cs_parameters_error_barrier();

  /* Verifications */

  if (rt_params->type > CS_RAD_TRANSFER_NONE) {

    /* Property fields */
    cs_rad_transfer_prp();

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
       0, 3);

    cs_parameters_error_barrier();
  }

  else
    return;

  /* Quadrature initialization */

  cs_rad_transfer_dir();

  /* Postprocessing */

  cs_gui_radiative_transfer_postprocess();
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
       N_("1 (Modak absorption coefficient)")};

  const char *imoadf_value_str[]
    = {N_("0 (no AFD model)"),
       N_("1 (ADF model with 8 wavelength intervals)"),
       N_("2 (ADF model with 50 wavelength intervals)")};

  const char *imfsck_value_str[]
    = {N_("0 (no FSCK model)"),
       N_("1 (FSCK model activated)"),
       N_("2 (FSCK model activated with tabulated properties)")};

  const char *idiver_value_str[]
    = {N_("0 (semi-analytic radiative S.T. calculation;\n"
          "                   "
          "   compulsory with transparent media)"),
       N_("1 (conservative radiative S.T. calculation)"),
       N_("2 (semi-analytic radiative S.T. calculation,\n"
          "                   "
          "   corrected for global conservation)")};

  cs_log_printf(CS_LOG_SETUP,
                  _("    idiver:        %s\n"),
                _(idiver_value_str[cs_glob_rad_transfer_params->idiver]));

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

  cs_log_printf(CS_LOG_SETUP,
                  _("    iimlum:        %s\n"),
                _(iimlum_value_str[cs_glob_rad_transfer_params->verbosity]));

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
        & CS_RAD_ATMO_3D_DIFFUSE_SOLAR)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    Diffuse solar atmospheric 3D model on\n"
           "      band id = %d\n"),
         cs_glob_rad_transfer_params->atmo_df_id);

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
