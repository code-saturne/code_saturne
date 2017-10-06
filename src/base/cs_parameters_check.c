/*============================================================================
 * Check computation parameters after user modification.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_1d_wall_thermal.h"
#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_lagr.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_rad_transfer.h"
#include "cs_restart_default.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_stokes_model.h"
#include "cs_syr_coupling.h"
#include "cs_wall_functions.h"
#include "cs_convection_diffusion.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_parameters_check.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/


/*!
  \file cs_parameters_check.c
        Parameters and options management check.
*/

/*----------------------------------------------------------------------------*/

/*!
  \enum cs_parameter_error_behavior_t

  \brief File acces modes

  \var CS_WARNING
       Warn only
  \var CS_ABORT_DELAYED
       Abort when \ref cs_parameters_error_barrier is called.
  \var CS_FILE_MODE_APPEND
       Abort immediately

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Counter for parameter checking errors */

static int  _param_check_errors = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * Build section description based on field label.
 *----------------------------------------------------------------------------*/

inline static char *
_field_section_desc(cs_field_t *f, const char *section_desc_b)
{
  const char *f_label = cs_field_get_label(f);

  /* 2 stands for the terminal character and one blank */
  int s_size =  cs_log_strlen(section_desc_b)
              + cs_log_strlen(f_label) + 2;

  char *section_desc = NULL;
  BFT_MALLOC(section_desc, s_size, char);
  snprintf(section_desc, s_size, "%s %s", section_desc_b, f_label);

  return section_desc;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print general parameters error or warning info.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param [in] format        format string, as printf() and family.
 * \param [in] ...           variable arguments based on format string.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error(cs_parameter_error_behavior_t   err_behavior,
                    const char                     *section_desc,
                    const char                     *format,
                    ...)
{
  cs_parameters_error_header(err_behavior, section_desc);

  int log_id = CS_LOG_DEFAULT;

  va_list  arg_ptr;
  va_start(arg_ptr, format);

  cs_log_vprintf(log_id, format, arg_ptr);

  va_end(arg_ptr);

  cs_parameters_error_footer(err_behavior);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print header for a given parameters error message type.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_header(cs_parameter_error_behavior_t   err_behavior,
                           const char                     *section_desc)
{
  const int err_type_id = (err_behavior <= CS_WARNING) ? 0 : 1;
  const char *error_type[] = {N_("Warning"),
                              N_("Error")};

  int log_id = CS_LOG_DEFAULT;

  if (section_desc != NULL)
    cs_log_printf(log_id,
                  "\n%s %s\n",
                  _(error_type[err_type_id]),
                  section_desc);
  else
    cs_log_printf(log_id, "\n%s\n", _(error_type[err_type_id]));
  size_t l = cs_log_strlen(_(error_type[err_type_id]));
  char underline[81];

  for (size_t i = 0; i < 80 && i < l; i++)
    underline[i] = '-';
  underline[CS_MIN(l,80)] = '\0';
  cs_log_printf(log_id, "%s\n", underline);

  if (err_behavior > CS_WARNING)
    _param_check_errors++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print footer for a given parameters error message type.
 *
 * \param[in]  err_behavior  warn or abort ?
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_footer(cs_parameter_error_behavior_t   err_behavior)
{
  if (err_behavior == CS_ABORT_IMMEDIATE)
    bft_error
      (__FILE__, __LINE__, 0,
       _("\nCheck your data and parameters (GUI and user subroutines)."));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword has values in a specified range.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  range_l       range lower bound (included)
 * \param[in]  range_u       range upper bound (excluded)
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_in_range_int(cs_parameter_error_behavior_t   err_behavior,
                              const char                     *section_desc,
                              const char                     *param_name,
                              int                             param_value,
                              int                             range_l,
                              int                             range_u)
{
  if (param_value < range_l || param_value >= range_u) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be in range [%d, %d].\n"),
                  param_name, param_value, range_l, range_u-1);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword has not values in a specified
 * range.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  range_l       range lower bound (included)
 * \param[in]  range_u       range upper bound (excluded)
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_not_in_range_int(cs_parameter_error_behavior_t   err_behavior,
                                  const char                     *section_desc,
                                  const char                     *param_name,
                                  int                             param_value,
                                  int                             range_l,
                                  int                             range_u)
{
  if (param_value >= range_l || param_value < range_u) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be out of range [%d, %d].\n"),
                  param_name, param_value, range_l, range_u-1);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword has values in a specified list.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  enum_size     size of possible enumeration
 * \param[in]  enum_values   optional list of enumerated values, or NULL
 *                           (in which case {0, ... enum_sizes-1} assumed
 * \param[in]  enum_names    optional list of value names, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_in_list_int(cs_parameter_error_behavior_t   err_behavior,
                             const char                     *section_desc,
                             const char                     *param_name,
                             int                             param_value,
                             int                             enum_size,
                             const int                      *enum_values,
                             const char                     *enum_names[])
{
  /* Check if we are in the defined list */

  if (enum_values != NULL) {
    for (int i = 0; i < enum_size; i++) {
      if (param_value == enum_values[i])
        return;
    }
  }
  else if (param_value >= 0 && param_value < enum_size)
    return;

  /* If we are not, report error */

  cs_parameters_error_header(err_behavior, section_desc);

  int log_id = CS_LOG_DEFAULT;

  if (enum_names != NULL) {
    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be one of:\n"),
                  param_name, param_value);
    for (int i = 0; i < enum_size; i++)
      cs_log_printf(log_id, "  %s\n", enum_names[i]);
  }
  else if (enum_values != NULL) {
    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be one of:\n"),
                  param_name, param_value);
    for (int i = 0; i < enum_size; i++)
      cs_log_printf(log_id, "  %d\n", enum_values[i]);
  }
  else {
    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be in range [%d, %d].\n"),
                  param_name, param_value, 0, enum_size-1);
  }

  cs_parameters_error_footer(err_behavior);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword does not have values in a
 * specified list.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  enum_size     size of possible enumeration
 * \param[in]  enum_values   optional list of enumerated values, or NULL
 *                           (in which case {0, ... enum_sizes-1} assumed
 * \param[in]  enum_names    optional list of value names, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_not_in_list_int(cs_parameter_error_behavior_t   err_behavior,
                                 const char                     *section_desc,
                                 const char                     *param_name,
                                 int                             param_value,
                                 int                             enum_size,
                                 const int                      *enum_values,
                                 const char                     *enum_names[])
{
  /* Check if we are in the defined list */

  int in = 0;
  if (enum_values != NULL) {
    for (int i = 0; i < enum_size; i++) {
      if (param_value == enum_values[i]) {
        in = 1;
        break;
      }
    }
  }
  else if (param_value >= 0 && param_value < enum_size)
    in = 1;

  /* If we are not, report error */

  if (in == 1) {
    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    if (enum_names != NULL) {
      cs_log_printf(log_id,
                    _("Parameter: %s = %d\n"
                      "while its value must not be one of:\n"),
                    param_name, param_value);
      for (int i = 0; i < enum_size; i++)
        cs_log_printf(log_id, "  %s\n", enum_names[i]);
    }
    else if (enum_values != NULL) {
      cs_log_printf(log_id,
                    _("Parameter: %s = %d\n"
                      "while its value must not be one of:\n"),
                    param_name, param_value);
      for (int i = 0; i < enum_size; i++)
        cs_log_printf(log_id, "  %d\n", enum_values[i]);
    }
    else {
      cs_log_printf(log_id,
                    _("Parameter: %s = %d\n"
                      "while its value must be out of range [%d, %d].\n"),
                    param_name, param_value, 0, enum_size-1);
    }

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword is equal to a specified value
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  std_value     compulsory or recommended parameter's value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_equal_int(cs_parameter_error_behavior_t   err_behavior,
                           const char                     *section_desc,
                           const char                     *param_name,
                           int                             param_value,
                           int                             std_value)
{
  if (param_value != std_value) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    if (err_behavior > CS_WARNING) {
      cs_log_printf(log_id,
                    _("Parameter: %s = %d\n"
                      "while its value must be equal to %d.\n"),
                    param_name, param_value, std_value);
    } else {
      cs_log_printf(log_id,
                    _("Parameter: %s = %d\n"
                      "while its recommended value is equal to %d.\n"),
                    param_name, param_value, std_value);
    }

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword is not equal to a specified value
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  fbd_value     forbidden value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_not_equal_int(cs_parameter_error_behavior_t   err_behavior,
                               const char                     *section_desc,
                               const char                     *param_name,
                               int                             param_value,
                               int                             fbd_value)
{
  if (param_value == fbd_value) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "which is a forbidden value.\n"),
                  param_name, param_value);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword is strictly positive
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_positive_int(cs_parameter_error_behavior_t   err_behavior,
                              const char                     *section_desc,
                              const char                     *param_name,
                              int                             param_value)
{
  if (param_value <= 0) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be strictly positive.\n"),
                  param_name, param_value);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given int keyword is greater than a specified value.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  ib_value     inferior bound value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_greater_int(cs_parameter_error_behavior_t   err_behavior,
                             const char                     *section_desc,
                             const char                     *param_name,
                             int                             param_value,
                             int                             ib_value)
{
  if (param_value < ib_value) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be greater than %d.\n"),
                  param_name, param_value, ib_value);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given double keyword has values in a specified range.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  range_l       range lower bound (included)
 * \param[in]  range_u       range upper bound (included)
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_in_range_double(cs_parameter_error_behavior_t   err_behavior,
                                 const char                     *section_desc,
                                 const char                     *param_name,
                                 double                          param_value,
                                 double                          range_l,
                                 double                          range_u)
{
  if (param_value < range_l || param_value > range_u) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %-5.3g\n"
                    "while its value must be in range [%-5.3g, %-5.3g].\n"),
                  param_name, param_value, range_l, range_u);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given double keyword has values in a specified list.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  enum_size     size of possible enumeration
 * \param[in]  enum_values   list of enumerated values
 * \param[in]  enum_names    optional list of value names, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_in_list_double(cs_parameter_error_behavior_t   err_behavior,
                                const char                     *section_desc,
                                const char                     *param_name,
                                double                          param_value,
                                int                             enum_size,
                                const double                   *enum_values,
                                const char                     *enum_names[])
{
  /* Check if we are in the defined list */

  if (enum_values != NULL) {
    for (int i = 0; i < enum_size; i++) {
      if (CS_ABS(param_value - enum_values[i]) > cs_math_epzero)
        return;
    }
  }

  /* If we are not, report error */

  cs_parameters_error_header(err_behavior, section_desc);

  int log_id = CS_LOG_DEFAULT;

  if (enum_names != NULL) {
    cs_log_printf(log_id,
                  _("Parameter: %s = %-5.3g\n"
                    "while its value must be one of:\n"),
                  param_name, param_value);
    for (int i = 0; i < enum_size; i++)
      cs_log_printf(log_id, "  %s\n", enum_names[i]);
  }
  else if (enum_values != NULL) {
    cs_log_printf(log_id,
                  _("Parameter: %s = %-5.3g\n"
                    "while its value must be one of:\n"),
                  param_name, param_value);
    for (int i = 0; i < enum_size; i++)
      cs_log_printf(log_id, "  %-5.3g\n", enum_values[i]);
  }

  cs_parameters_error_footer(err_behavior);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given double keyword is equal to a specified value.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  std_value     compulsory or recommended parameter's value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_equal_double(cs_parameter_error_behavior_t   err_behavior,
                              const char                     *section_desc,
                              const char                     *param_name,
                              double                          param_value,
                              double                          std_value)
{
  if (CS_ABS(param_value-std_value) > cs_math_epzero) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    if (err_behavior > CS_WARNING) {
      cs_log_printf(log_id,
                    _("Parameter: %s = %-5.3g\n"
                      "while its value must be equal to %-5.3g.\n"),
                    param_name, param_value, std_value);
    } else {
      cs_log_printf(log_id,
                    _("Parameter: %s = %-5.3g\n"
                      "while its recommended value is equal to %-5.3g.\n"),
                    param_name, param_value, std_value);
    }

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given double keyword is greater than a specified value.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  ib_value     inferior bound value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_greater_double(cs_parameter_error_behavior_t   err_behavior,
                                const char                     *section_desc,
                                const char                     *param_name,
                                double                          param_value,
                                double                          ib_value)
{
  if (param_value < ib_value) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %-5.3g\n"
                    "while its value must be greater than %-5.3g.\n"),
                  param_name, param_value, ib_value);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Abort if the the parameter errors count is nonzero.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_barrier(void)
{
  cs_lnum_t n_errors = _param_check_errors;
  cs_parall_counter_max(&n_errors, 1);

  if (n_errors > 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%d parameter error(s) reported.\n"
         "\n"
         "Read error messages above for details, then\n"
         "check your data and parameters (GUI and user subroutines)."),
       n_errors);

  _param_check_errors = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check computation parameters after user modification.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_check(void)
{
  int n_fields = cs_field_n_fields();
  cs_var_cal_opt_t var_cal_opt;
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int keyvar = cs_field_key_id("variable_id");
  const int kcpsyr = cs_field_key_id("syrthes_coupling");
  const int kivisl = cs_field_key_id("scalar_diffusivity_id");
  const int kvisls0 = cs_field_key_id("scalar_diffusivity_ref");
  const int restart_file_key_id = cs_field_key_id("restart_file");
  const int key_limiter = cs_field_key_id("limiter_choice");

  cs_field_t *f_pot = NULL;
  if (cs_glob_physical_model_flag[CS_GROUNDWATER] > 0)
    f_pot = CS_F_(head);
  else
    f_pot = CS_F_(p);

  cs_field_t *f_th = cs_thermal_model_field();

  char *f_desc = NULL;

  /*--------------------------------------------------------------------------
   * listing log frequency
   *--------------------------------------------------------------------------*/

  cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_log_frequency",
                               cs_glob_log_frequency,
                               -1);

  /*--------------------------------------------------------------------------
   * Computation parameters
   *--------------------------------------------------------------------------*/

  /* Thermal model */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_thermal_model->itherm",
                                cs_glob_thermal_model->itherm,
                                CS_THERMAL_MODEL_NONE, CS_THERMAL_MODEL_N_TYPES);

  /* Constant or variable rho and viscosity */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_fluid_properties->irovar",
                                cs_glob_fluid_properties->irovar,
                                0, 2);
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_fluid_properties->ivivar",
                                cs_glob_fluid_properties->irovar,
                                0, 2);

  /* Equations definition, time scheme, convective scheme */
  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      f_desc = _field_section_desc(f, "while reading numerical "
                                      "parameters for variable");

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.iconv (convection flag)",
                                    var_cal_opt.iconv,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.istat (unsteadiness flag)",
                                    var_cal_opt.istat,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.idiff (diffusion flag)",
                                    var_cal_opt.idiff,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.idifft (turbulent diffusion "
                                                                       "flag)",
                                    var_cal_opt.idifft,
                                    0, 2);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "var_cal_opt.thetav (theta-scheme)",
                                       var_cal_opt.thetav,
                                       0., 1.);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "var_cal_opt.blencv (2nd order scheme "
                                       "share for convection)",
                                       var_cal_opt.blencv,
                                       0., 1.);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "var_cal_opt.blend_st (2nd order scheme "
                                       "share for convection)",
                                       var_cal_opt.blend_st,
                                       0., 1.);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.ischcv (2nd order scheme "
                                    "choice)",
                                    var_cal_opt.ischcv,
                                    0, 4);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.isstpc (slope test)",
                                    var_cal_opt.isstpc,
                                    -1, 4);

      BFT_FREE(f_desc);
    }
  }

  /* check if NVD scheme for thermal scalar is not one of the VOF schemes */
  if (f_th != NULL) {
    cs_field_get_key_struct(f_th, key_cal_opt_id, &var_cal_opt);
    if (var_cal_opt.isstpc >= 3) { /* NVD scheme on thermal scalar? */
      cs_nvd_type_t limiter_choice = cs_field_get_key_int(f_th, key_limiter);

      f_desc = _field_section_desc(f_th, "while reading numerical "
                                         "parameters for variable");

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "NVD scheme",
                                    limiter_choice,
                                    CS_NVD_GAMMA, CS_NVD_VOF_HRIC);
    }
  }

  /* thetav for pressure must be equal to 1 */
  cs_field_get_key_struct(f_pot, key_cal_opt_id, &var_cal_opt);
  f_desc = _field_section_desc(f_pot, "while reading numerical "
                                      "parameters for variable");

  cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                _(f_desc),
                                "var_cal_opt.thetav (theta-scheme)",
                                var_cal_opt.thetav,
                                1.);

  BFT_FREE(f_desc);

  /* In LES, additional consistency checkings are needed *
   * Only a warning for non standard parameters, but stop if
   * there is more than 5% of upwind.
   * Centered scheme with/without slope test, nwsrsm */
  if (cs_glob_turb_model->itytur == 4) {
    cs_field_get_key_struct(CS_F_(u), key_cal_opt_id, &var_cal_opt);
    f_desc = _field_section_desc(CS_F_(u), "in LES, while reading time "
                                           "scheme parameters for variable");

    cs_parameters_is_equal_double(CS_WARNING,
                                  _(f_desc),
                                  "var_cal_opt.thetav (theta-scheme)",
                                  var_cal_opt.thetav,
                                  0.5);

    BFT_FREE(f_desc);

    f_desc = _field_section_desc(CS_F_(u), "in LES, while reading "
                                           "convection scheme "
                                           "parameters for variable");

    cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                     _(f_desc),
                                     "var_cal_opt.blencv (2nd order scheme "
                                                     "share for convection)",
                                     var_cal_opt.blencv,
                                     0.95, 1.);

    cs_parameters_is_equal_double(CS_WARNING,
                                  _(f_desc),
                                  "var_cal_opt.blencv (2nd order scheme "
                                                  "share for convection)",
                                  var_cal_opt.blencv,
                                  1.);

    cs_parameters_is_equal_int(CS_WARNING,
                               _(f_desc),
                               "var_cal_opt.isstpc (slope test)",
                               var_cal_opt.isstpc,
                               1);

    BFT_FREE(f_desc);

    for (int f_id = 0 ; f_id < n_fields ; f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int isca = cs_field_get_key_int(f, keysca);
      if (isca > 0) {
        cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
        f_desc = _field_section_desc(f, "in LES, while reading time "
                                        "scheme parameters for variable");

        cs_parameters_is_equal_double(CS_WARNING,
                                      _(f_desc),
                                      "var_cal_opt.thetav (theta-scheme)",
                                      var_cal_opt.thetav,
                                      0.5);

        BFT_FREE(f_desc);

        f_desc = _field_section_desc(f, "in LES, while reading "
                                        "convection scheme "
                                        "parameters for variable");

        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _(f_desc),
                                         "var_cal_opt.blencv (2nd order "
                                           "scheme share for convection)",
                                         var_cal_opt.blencv,
                                         0.95, 1.);

        cs_parameters_is_equal_double(CS_WARNING,
                                      _(f_desc),
                                      "var_cal_opt.blencv (2nd order scheme "
                                                      "share for convection)",
                                      var_cal_opt.blencv,
                                      1.);

        cs_parameters_is_equal_int(CS_WARNING,
                                   _(f_desc),
                                   "var_cal_opt.isstpc (slope test)",
                                   var_cal_opt.isstpc,
                                   0);

        BFT_FREE(f_desc);
      }
    }
  }

  /* navsto sub-iterations
   * It must be a integer greater or equal to 1
   * For the moment, we forbid nterup > 1 with
   * estimators, weight matrix (reinforced U-P coupling), hydrostatic pressure
   * and steady algorithm. */
  cs_parameters_is_positive_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_piso->nterup (Navier-Stokes sub-iterations)",
                                cs_glob_piso->nterup);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_time_step_options->idtvar",
                                cs_glob_time_step_options->idtvar,
                                -1, 3);

  /* Steady Algorithm */
  if (cs_glob_time_step_options->idtvar < 0) {
    for (int f_id = 0 ; f_id < n_fields ; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int ivar = cs_field_get_key_int(f, keyvar);
      if (ivar > 0) {
        cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
        f_desc = _field_section_desc(f, "With the steady algorithm "
                                        "(SIMPLE), "
                                        "while reading the relaxation "
                                        "coefficient for variable");

        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _(f_desc),
                                         "var_cal_opt.relaxv (relax. coef.)",
                                         var_cal_opt.relaxv,
                                         0., 1.);

        BFT_FREE(f_desc);
      }
    }

    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("The steady algorithm is not compatible "
                                 "with the lagrangian module which is "
                                 "time-dependant by nature"),
                               "cs_glob_lagr_time_scheme->iilagr",
                               cs_glob_lagr_time_scheme->iilagr,
                               0);

    int les_iturb[3] = {40, 41, 42};
    cs_parameters_is_not_in_list_int(CS_ABORT_DELAYED,
                                     _("The steady algorithm is not compatible "
                                       "with L.E.S. turbulence modelling "
                                       "which is time-dependant by nature"),
                                     "cs_glob_turb_model->iturb",
                                     cs_glob_turb_model->iturb,
                                     3,
                                     les_iturb,
                                     NULL);
  }

  /* Gradient reconstruction */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_space_disc->imrgra",
                                cs_glob_space_disc->imrgra,
                                -10, 11);

  int imrgrl = CS_ABS(cs_glob_space_disc->imrgra);
  imrgrl = imrgrl%10;

  /* We check the non-orthogonality angle of the extended neighborhood
   * in case of a least squares gradient method */
  if (imrgrl == 3 || imrgrl == 6 || imrgrl == 9) {
    cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                     _("while reading gradient "
                                       "reconstruction parameters for "
                                       "least squares method on reduced "
                                       "extended neighborhood"),
                                     "cs_glob_space_disc->anomax "
                                     "(max. non-orthogonality angle in radians)",
                                     cs_glob_space_disc->anomax,
                                     0., cs_math_pi*0.5);
  }

  /* Numbers of sweeps don't have to be checked : they are simply
   * integers (negative if one wants to be sure to never enter the loops */
  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int ivar = cs_field_get_key_int(f, keyvar);
    if (ivar > 0) {
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      f_desc = _field_section_desc(f, "Wile reading numerical parameters "
                                      " for variable");

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.imligr "
                                    "(gradient limitation method)",
                                    var_cal_opt.imligr,
                                    -1, 2);

      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _(f_desc),
                                      "var_cal_opt.climgr "
                                      "(gradient limitation coeffcient)",
                                      var_cal_opt.climgr,
                                      1.);

      /* Extrag can be different from zero only for the pressure
         and in that case equal to 1 */
      if (f_id == f_pot->id) {
        const cs_real_t extrag_vals[2] = {0., 1.};
        cs_parameters_is_in_list_double(CS_ABORT_DELAYED,
                                        _(f_desc),
                                        "var_cal_opt.extrag",
                                        var_cal_opt.extrag,
                                        2,
                                        extrag_vals,
                                        NULL);
      } else {
        cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                      _(f_desc),
                                      "var_cal_opt.extrag",
                                      var_cal_opt.extrag,
                                      0.);
      }

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "var_cal_opt.ircflu (fluxes "
                                    "reconstruction)",
                                    var_cal_opt.ircflu,
                                    0, 2);

      BFT_FREE(f_desc);

      if (   var_cal_opt.ischcv == 0
          && CS_ABS(var_cal_opt.blencv) > cs_math_epzero) {
        f_desc = _field_section_desc(f, "Second order linear upwind "
                                        "enabled for variable ");

        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "var_cal_opt.ircflu (fluxes reconstruction)",
                                   var_cal_opt.ircflu,
                                   1);

        BFT_FREE(f_desc);
      }
    }
  }

  /* Time stepping */
  /* The number of time step could be negative : no test */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_time_step_options->inpdt0",
                                cs_glob_time_step_options->inpdt0,
                                0, 2);

  if (cs_glob_time_step_options->idtvar >= 0) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_time_step_options->dtref",
                                    cs_glob_time_step_options->dtref,
                                    0.);
  }

  if (cs_glob_time_step_options->idtvar > 0) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_time_step_options->varrdt",
                                    cs_glob_time_step_options->varrdt,
                                    0.);

    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_time_step_options->dtmax",
                                    cs_glob_time_step_options->dtmax,
                                    0.);

    cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                     _("while reading input data"),
                                     "cs_glob_time_step_options->dtmin",
                                     cs_glob_time_step_options->dtmin,
                                     0., cs_glob_time_step_options->dtmax);
  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_time_step_options->iptlro "
                                "(time-step clipping related to buyoancy)",
                                cs_glob_time_step_options->iptlro,
                                0, 2);

  /* If constant time step, the time step is left unmodified
     but we highlight the local criteria exceedances */
  if (cs_glob_time_step_options->idtvar == 0) {
    cs_parameters_is_equal_int(CS_WARNING,
                               _("while reading input data,\n"
                                 "a time-step clipping is enabled but a "
                                 "constant time-step is set.\nPotential local "
                                 "criteria violations will be highlighted\n"
                                 "but the time-step won't be clipped."),
                               "cs_glob_time_step_options->iptlro",
                               cs_glob_time_step_options->iptlro,
                               0);
  } else if (cs_glob_time_step_options->idtvar == -1) {
    cs_parameters_is_equal_int(CS_WARNING,
                               _("while reading input data,\n"
                                 "a time-step clipping is enabled but a "
                                 "steady-state algorithm (SIMPLE) is set.\n "
                                 "This setting will be ignored."),
                               "cs_glob_time_step_options->iptlro",
                               cs_glob_time_step_options->iptlro,
                               0);
  }

  /* Turbulence */

  /* Model */
  const int iturb_vals[14] = {0,          /* laminar */
                              10,         /* mixing length */
                              20, 21,     /* k-epsilon */
                              30, 31, 32, /* RSM */
                              40, 41, 42, /* LES */
                              50, 51,     /* v2f */
                              60,         /* k-omega */
                              70};        /* Spallart-Allmaras */

  cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_turb_model->iturb",
                               cs_glob_turb_model->iturb,
                               14,
                               iturb_vals,
                               NULL);

  /* Rotation curvature correction for eddy-viscosity models */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_turb_rans_model->irccor",
                                cs_glob_turb_rans_model->irccor,
                                0, 2);

  /* Rotation curvature correction compatible only with RANS eddy-viscosity
     models */
  if (cs_glob_turb_rans_model->irccor == 1) {
    const int iturb_evm_vals[6] = {20, 21,     /* k-epsilon */
                                   50, 51,     /* v2f */
                                   60,         /* k-omega */
                                   70};        /* Spallart-Allmaras */

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "rotation curvature correction is only "
                                   "compatible with eddy viscosity turbulence "
                                   "models"),
                                 "cs_glob_turb_model->iturb",
                                 cs_glob_turb_model->iturb,
                                 6,
                                 iturb_evm_vals,
                                 NULL);
  }

  /* Delayed Detached Eddy Simulation model compatible only with RANS
     k-omega SST model */
  if (cs_glob_turb_rans_model->iddes == 1) {
    const int iturb_ddes_vals[1] = {60} ;     /* k-omega */

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "Delayed Detached Eddy Simulation is only "
                                   "compatible with k-omega SST model "
                                   "(iturb=60)"),
                                 "cs_glob_turb_model->iturb",
                                 cs_glob_turb_model->iturb,
                                 1,
                                 iturb_ddes_vals,
                                 NULL);

  }

  /* In lagrangian with two-way coupling, k-omega SST is forbidden (not
     properly implemented) */
  if (cs_glob_lagr_time_scheme->iilagr == 2) {
    cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "two way coupling in lagrangian modelling "
                                     "is not compatible with k-omega SST "
                                     "turbulence model"),
                                   "cs_glob_turb_model->iturb",
                                   cs_glob_turb_model->iturb,
                                   60);
  }

  /* Vortex method for LES */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_turb_les_model->ivrtex",
                                cs_glob_turb_les_model->ivrtex,
                                0, 2);

  if (   cs_glob_turb_model->itytur != 4
      && cs_glob_turb_les_model->ivrtex == 1) {
    cs_parameters_is_equal_int(CS_WARNING,
                               _("while reading input data,\n"
                                 "synthetic vortex method only for LES , "
                                 "this setting will be ignored"),
                               "cs_glob_turb_les_model->ivrtex",
                               cs_glob_turb_les_model->ivrtex,
                               0);

    cs_turb_les_model_t *les_model = cs_get_glob_turb_les_model();
    les_model->ivrtex = 0;
  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_wall_functions->iwallf",
                                cs_glob_wall_functions->iwallf,
                                0, 7);

  if (cs_glob_wall_functions->iwallf > 2) {
    const int itytur_vals[5] = {2, 3, 5, 6, 7};

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "two scales wall function not compatible "
                                   "with a laminar, mixing length or LES "
                                   "computation"),
                                 "cs_glob_turb_model->itytur",
                                 cs_glob_turb_model->itytur,
                                 5,
                                 itytur_vals,
                                 NULL);
  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_wall_functions->iwallt",
                                cs_glob_wall_functions->iwallt,
                                0, 2);

  /* Specific k-epsilon, v2f and k-omega */
  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 5
      || cs_glob_turb_model->itytur == 6) {
    /* iclkep option not available in k-omega */
    if (cs_glob_turb_model->iturb != 60) {
      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_turb_rans_model->iclkep",
                                    cs_glob_turb_rans_model->iclkep,
                                    0, 2);
    }

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->ikecou",
                                  cs_glob_turb_rans_model->ikecou,
                                  0, 2);

    /* En k-eps a prod lin et en v2f on force IKECOU a 0 */
    if (   cs_glob_turb_model->iturb == 21
        || cs_glob_turb_model->itytur == 5) {
      cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "with k-epsilon LP (iturb=21) or v2f model "
                                   "(iturb=50/51)"),
                                 "cs_glob_turb_rans_model->ikecou",
                                 cs_glob_turb_rans_model->ikecou,
                                 0);
    }

    /* En stationnaire on force IEKCOU a 0 */
    if (cs_glob_time_step_options->idtvar < 0) {
      cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "with steady-state algorithm"),
                                 "cs_glob_turb_rans_model->ikecou",
                                 cs_glob_turb_rans_model->ikecou,
                                 0);
    }

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->igrhok",
                                  cs_glob_turb_rans_model->igrhok,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->igrake",
                                  cs_glob_turb_rans_model->igrake,
                                  0, 2);

    /* if relaxv(ik) was modified by the user, and if ikecou=1, we warn the
       user that relaxv settings won't have any effect,
       otherwise check that relaxv is in range [O,1] (already done in steady) */

    if (cs_glob_turb_model->itytur == 6) { /* k-omega */
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(omg), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxvo = var_cal_opt.relaxv;

      if (   CS_ABS(relaxvk + 999.) > cs_math_epzero
          || CS_ABS(relaxvo + 999.) > cs_math_epzero) {
        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "modifications of relaxation coefficients "
                                     "for turbulence variables k and omega "
                                     "will be ignored"),
                                   "cs_glob_turb_rans_model->ikecou",
                                   cs_glob_turb_rans_model->ikecou,
                                   0);
      }

      if (   cs_glob_turb_rans_model->ikecou == 0
          && cs_glob_time_step_options->idtvar >= 0) {
        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _("while reading input data "
                                           "for turbulent variable k"),
                                         "var_cal_opt.relaxv",
                                         relaxvk,
                                         0, 1);
        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _("while reading input data "
                                           "for turbulent variable omega"),
                                         "var_cal_opt.relaxv",
                                         relaxvo,
                                         0, 1);
      }
    } else { /* k-epsilon */
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxve = var_cal_opt.relaxv;

      if (   CS_ABS(relaxvk + 999.) > cs_math_epzero
          || CS_ABS(relaxve + 999.) > cs_math_epzero) {
        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "modifications of relaxation coefficients "
                                     "for turbulence variables k and epsilon "
                                     "will be ignored"),
                                   "cs_glob_turb_rans_model->ikecou",
                                   cs_glob_turb_rans_model->ikecou,
                                   0);
      }

      if (   cs_glob_turb_rans_model->ikecou == 0
          && cs_glob_time_step_options->idtvar >= 0) {
        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _("while reading input data "
                                           "for turbulent variable k"),
                                         "var_cal_opt.relaxv",
                                         relaxvk,
                                         0, 1);
        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _("while reading input data "
                                           "for turbulent variable epsilon"),
                                         "var_cal_opt.relaxv",
                                         relaxve,
                                         0, 1);
      }
    }
  }

  /* Specifique Rij-epsilon */
  if (cs_glob_turb_model->itytur == 3) {
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->irijnu",
                                  cs_glob_turb_rans_model->irijnu,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->irijrb",
                                  cs_glob_turb_rans_model->irijrb,
                                  0, 2);

    /* wall echo and specific implicitation of the diffusion of epsilon only
       in Rij LRR */
    if (cs_glob_turb_model->iturb == 30) {
      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_turb_rans_model->irijec",
                                    cs_glob_turb_rans_model->irijec,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_turb_rans_model->idifre",
                                    cs_glob_turb_rans_model->idifre,
                                    0, 2);
    }

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->igrari",
                                  cs_glob_turb_rans_model->igrari,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->iclsyr",
                                  cs_glob_turb_rans_model->iclsyr,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->iclptr",
                                  cs_glob_turb_rans_model->iclptr,
                                  0, 2);
  }

  /* Specifique LES */
  if (cs_glob_turb_model->itytur == 4) {
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_les_model->idries",
                                  cs_glob_turb_les_model->idries,
                                  0, 2);


    if (   cs_glob_turb_model->iturb == 41
        || cs_glob_turb_model->iturb == 42) {
      cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "Van Driest near wall damping not "
                                   "compatible with dynamic Smagorinsky or "
                                   "Wale models"),
                                 "cs_glob_turb_les_model->idries",
                                 cs_glob_turb_les_model->idries,
                                 0);
    }

    /* The reduction of the extended neighborhood can degrade the results of the
       LES dynamic model */
    if (cs_glob_turb_model->iturb == 41) {
      int imrgra_cmp = CS_ABS(cs_glob_turb_model->iturb);
      switch(imrgra_cmp) {
      case 3:
      case 6:
      case 9:
        cs_parameters_error
          (CS_WARNING,
           _("while reading input data"),
           _("A reduction of the extended neighborhood was selected for the\n"
             "calculation of the gradients by least squares.\n"
             "This will also be applied to the averaging in the "
             "selected LES Dynamic model.\n"
             "The computation will run, but the averaging of the Smagorinsky\n"
             "constant can be degraded, "
             "as it uses the same reduced neighborhood.\n"
             "Use fully extended neighborhood or directly "
             "define the averaging of the\n"
             "dynamic Smagorinsky constant via the ussmag subroutine."));
        break;
      default:
        break;
      }
    }
  }

  /* Stokes */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_stokes_model->iprco",
                                cs_glob_stokes_model->iprco,
                                0, 2);

  if (cs_glob_stokes_model->iprco == 1) {
    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_stokes_model->irevmc",
                               cs_glob_stokes_model->irevmc,
                               0);

    cs_real_t arakfr = cs_glob_stokes_model->arak;
    if (cs_glob_time_step_options->idtvar < 0) {
      cs_field_get_key_struct(CS_F_(u), key_cal_opt_id, &var_cal_opt);
      arakfr *= var_cal_opt.relaxv;
    }

    cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                     _("while reading input data"),
                                     "cs_glob_stokes_model->arak",
                                     arakfr,
                                     0., 1.);
  }

  /* U-P reinforced coupling */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_stokes_model->ipucou",
                                cs_glob_stokes_model->ipucou,
                                0, 2);

  cs_field_get_key_struct(CS_F_(u), key_cal_opt_id, &var_cal_opt);
  /* steady or variable time step time algorithm not compatible with theta
     scheme with theta different from 1 for the velocity */
  if (cs_glob_time_step_options->idtvar != 0) {
    cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                  _("while reading time scheme parameters,\n"
                                    "theta-scheme with theta different from 1 "
                                    "for the velocity\n"
                                    "only compatible with constant time step "
                                    "unsteady algorithm"),
                                  "var_cal_opt.thetav",
                                  var_cal_opt.thetav,
                                  1);
  }
  /* U-P reinforced coupling not compatible with theta scheme with theta
     different from 1 for the velocity */
  if (cs_glob_stokes_model->ipucou == 1) {
    cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                  _("while reading time scheme parameters,\n"
                                    "theta-scheme with theta different from 1 "
                                    "for the velocity\n"
                                    "not compatible with reinforced "
                                    "velocity-pressure coupling"),
                                  "var_cal_opt.thetav",
                                  var_cal_opt.thetav,
                                  1);
  }

  /* frozen velocity field */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_stokes_model->iccvfg",
                                cs_glob_stokes_model->iccvfg,
                                0, 2);

  /* face viscosity interpolation */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_space_disc->imvisf",
                                cs_glob_space_disc->imvisf,
                                0, 2);

  /* thermal SYRTHES coupling */

  /* check if values of key syrthes_coupling are realistic */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    f_desc = _field_section_desc(f, "while reading parameters for "
                                    "field ");
    int icpsyr = cs_field_get_key_int(f, kcpsyr);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _(f_desc),
                                  "key 'syrthes coupling'",
                                  icpsyr,
                                  0, 2);

    BFT_FREE(f_desc);
  }

  /* the number of coupled scalars is counted */
  int nbsccp = 0, n_coupl;
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    nbsccp += cs_field_get_key_int(f, kcpsyr);
  }

  /* On regarde s'il y a du couplage */
  n_coupl = cs_syr_coupling_n_couplings();

  /* if no coupling with SYRTHES */
  if (n_coupl == 0) {
    /* and scalars are defined as coupled with SYRTHES */
    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("Inconsistency in SYRTHES coupling settings,\n"
                                 "no coupling with SYRTHES defined but some "
                                 "scalars are defined as coupled"),
                               "number of coupled scalars",
                               nbsccp,
                               0);
  } else { /* if coupling with SYRTHES */
    /* and more than one scalar is defined as coupled */
    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("Inconsistency in SYRTHES coupling settings,\n"
                                 "more than one "
                                 "scalars are defined as coupled"),
                               "number of coupled scalars",
                               nbsccp,
                               1);

    /* and the coupled scalar is not the thermal scalar
       or if the thermal scalar is not the temperature
       if compressible is not enabled */
    for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
      cs_field_t  *f = cs_field_by_id(f_id);
      int icpsyr = cs_field_get_key_int(f, kcpsyr);
      if (icpsyr == 1) {
        /* compressible */
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > 0) {
          cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                     _("Inconsistency in SYRTHES coupling "
                                       "settings,\n"
                                       "with the compressible module, the "
                                       "coupled scalar must be the energy "
                                       "(use thermal scalar number "
                                       "cs_glob_thermal_model->iscalt)"),
                                     "number of coupled scalar",
                                     cs_field_get_key_int(f, keysca),
                                     cs_glob_thermal_model->iscalt);
        } else { /* not compressible */
          cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                     _("Inconsistency in SYRTHES coupling "
                                       "settings,\n"
                                       "the coupled scalar must be the "
                                       "the thermal scalar \n"
                                       "(use thermal scalar number "
                                       "cs_glob_thermal_model->iscalt)"),
                                     "number of coupled scalar",
                                     cs_field_get_key_int(f, keysca),
                                     cs_glob_thermal_model->iscalt);

          cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                     _("Inconsistency in SYRTHES coupling "
                                       "settings,\n"
                                       "the thermal scalar must be "
                                       "the temperature"),
                                     "cs_glob_thermal_model->itherm",
                                     cs_glob_thermal_model->itherm,
                                     CS_THERMAL_MODEL_TEMPERATURE);
        }
      }
    }
  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_stokes_model->itbrrb",
                                cs_glob_stokes_model->itbrrb,
                                0, 2);

  /* Dynamic relaxation option */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      if (var_cal_opt.iswdyn >= 1) {
        f_desc = _field_section_desc(f, "Dynamic relaxation enabled for "
                                        "variable ");
        /* The number of reconstruction sweeps is set to 20 at least */
        cs_parameters_is_greater_int(CS_WARNING,
                                     _(f_desc),
                                     "var_cal_opt.nswrsm",
                                     var_cal_opt.nswrsm,
                                     19);

        BFT_FREE(f_desc);

        if (var_cal_opt.nswrsm < 20) {
          var_cal_opt.nswrsm = 20;
          int log_id = CS_LOG_DEFAULT;
          cs_log_printf(log_id,
                        _("The calculation continues with nswrsm = 20.\n"));
          cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
        }

        if (f->id == CS_F_(u)->id) {
          cs_parameters_is_equal_int(CS_WARNING,
                                     _("Dynamic relaxation enabled for "
                                       "variable velocity.\n"
                                       "The slope test must be disabled."),
                                     "var_cal_opt.isstpc",
                                     var_cal_opt.isstpc,
                                     1);

          if (var_cal_opt.isstpc == 0) {
            var_cal_opt.isstpc = 1;
            int log_id = CS_LOG_DEFAULT;
            cs_log_printf(log_id,
                          _("The calculation continues with isstpc = 1 for "
                            "the variable velocity.\n"));
            cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
          }
        }
      }
    }
  }

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      f_desc = _field_section_desc(f, "Number of RHS reconstruction "
                                      "sweeps for variable ");

      cs_parameters_is_positive_int(CS_WARNING,
                                    _(f_desc),
                                    "var_cal_opt.nswrsm",
                                    var_cal_opt.nswrsm);

      BFT_FREE(f_desc);

      if (var_cal_opt.nswrsm <= 0) {
        var_cal_opt.nswrsm = 1;
        int log_id = CS_LOG_DEFAULT;
        cs_log_printf(log_id,
                      _("The calculation continues with nswrsm = 1 "
                        "for this variable.\n"));
        cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
      }
    }
  }

  /* physical properties reference values */
  cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                  _("while reading reference density value"),
                                  "cs_glob_fluid_properties->ro0",
                                  cs_glob_fluid_properties->ro0,
                                  0.);

  cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                  _("while reading reference dynamic "
                                    "viscosity value"),
                                  "cs_glob_fluid_properties->viscl0",
                                  cs_glob_fluid_properties->viscl0,
                                  0.);

  /* check variances */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      f_desc = _field_section_desc(f, "while reading parameters for "
                                      "variance ");

      cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                     _(f_desc),
                                     "associated fluctuating field scalar id",
                                     cs_parameters_iscavr(f),
                                     -1);

      int fscavr_id = cs_field_get_key_int(f, kscavr);
      if (fscavr_id > -1) {
        cs_field_t *fluct_f = cs_field_by_id(fscavr_id);

        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "associated fluctuating field is already a "
                                   "variance.\n Associated fluctuating field "
                                   "scalar id for this variance",
                                   cs_parameters_iscavr(fluct_f),
                                   0);
      }

      BFT_FREE(f_desc);
    }
  }

  /* check positivity of reference diffusivities */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_real_t visls0 = cs_field_get_key_double(f, kvisls0);
      f_desc = _field_section_desc(f, "reference diffusivity value for "
                                      "variable ");

      int diff_id = cs_field_get_key_int(f, kivisl);
      int isca = cs_field_get_key_int(f, keysca);

      if (isca > -1 && diff_id == -1) {
        cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                        _(f_desc),
                                        "key scalar_diffusivity_ref",
                                        visls0,
                                        0.);
      }

      BFT_FREE(f_desc);
    }
  }

  /* Turbulence */

  /* uref is needed for automatic initialisation and boundary conditions
     of turbulence variables check that it is positive
     and warn the user if not */
  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 3
      || cs_glob_turb_model->itytur == 5
      || cs_glob_turb_model->iturb == 60
      || cs_glob_turb_model->iturb == 70) {
    cs_parameters_is_greater_double(CS_WARNING,
                                    _("Reference velocity is used for the "
                                      "automatic initialisation and boundary "
                                      "conditions of the turbulence "
                                      "variables\n(if not overridden by "
                                      "user or via a restart file)"),
                                    "uref (reference velocity)",
                                    cs_glob_turb_ref_values->uref,
                                    0.);
  }

  if (cs_glob_turb_model->iturb == 10) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_turb_rans_model->xlomlg "
                                    "(mixing length)",
                                    cs_glob_turb_rans_model->xlomlg,
                                    0.);
  }

  /* LES (check also in Fortran for now because constants are duplicated) */
  if (cs_glob_turb_model->itytur == 4) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("LES dynamic model constant"),
                                    "cs_turb_xlesfl",
                                    cs_turb_xlesfl,
                                    0.);

    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("LES dynamic model constant"),
                                    "cs_turb_ales",
                                    cs_turb_ales,
                                    0.);

    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("LES dynamic model constant"),
                                    "cs_turb_bles",
                                    cs_turb_bles,
                                    0.);

    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("LES dynamic model constant"),
                                    "cs_turb_csmago",
                                    cs_turb_csmago,
                                    0.);

    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("LES dynamic model constant"),
                                    "cs_turb_cwale",
                                    cs_turb_cwale,
                                    0.);

    if (cs_glob_turb_les_model->idries == 1) {
      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _("LES dynamic model constant"),
                                      "cs_turb_cdries",
                                      cs_turb_cdries,
                                      0.);
    }

    if (cs_glob_turb_model->iturb == 41) {
      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _("LES dynamic model constant"),
                                      "cs_turb_xlesfd",
                                      cs_turb_xlesfd,
                                      0.);

      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _("LES dynamic model constant"),
                                      "cs_turb_smagmx",
                                      cs_turb_smagmx,
                                      0.);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _("LES dynamic model constant"),
                                       "cs_turb_smagmn (must be < "
                                       "cs_turb_smagmx)",
                                       cs_turb_smagmn,
                                       0., cs_turb_smagmx);
    }
  }

  /*--------------------------------------------------------------
   * Rotating frame and unsteady rotor/stator coupling.
   *--------------------------------------------------------------*/

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("flag for consideration of rotation"),
                                "cs_glob_physical_constants->icorio",
                                cs_glob_physical_constants->icorio,
                                0, 2);

  /* Unsteady rotor/stator coupling is not compatible with the Lagrangian
     module */
  int model = cs_turbomachinery_get_model();

  if (model == CS_TURBOMACHINERY_TRANSIENT) {
    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("transient turbomachinery module not "
                                 "compatible "
                                 "with the lagrangian module"),
                               "cs_glob_lagr_time_scheme->iilagr",
                               cs_glob_lagr_time_scheme->iilagr,
                               0);
  }

  if (model > CS_TURBOMACHINERY_NONE) {
    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("turbomachinery module not compatible "
                                 "with computation in relative reference "
                                 "frame.\n"
                                 "If necessary, use the turbomachinery module "
                                 "with several rotors."),
                               "cs_glob_physical_constants->icorio",
                               cs_glob_physical_constants->icorio,
                               0);
  }

  /* Check the consistency of the restart_file key */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    cs_restart_file_t r_id = cs_field_get_key_int(f, restart_file_key_id);
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                  "restart_file",
                                  r_id,
                                  CS_RESTART_DISABLED,
                                  CS_RESTART_N_RESTART_FILES);

    if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_NONE) {
      cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "the lagrangian checkpoint cannot "
                                     "be used\n"
                                     "without having the lagrangian module "
                                     "enabled."),
                                     "restart_file",
                                     r_id,
                                     CS_RESTART_LAGR_STAT);

      cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "the lagrangian checkpoint cannot "
                                     "be used\n"
                                     "without having the lagrangian module "
                                     "enabled."),
                                     "restart_file",
                                     r_id,
                                     CS_RESTART_LAGR);
    }

    if (cs_glob_1d_wall_thermal->nfpt1d == 0)
      cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "the 1D wall thermal checkpoint cannot "
                                     "be used\n"
                                     "without having the 1D wall thermal module "
                                     "enabled."),
                                     "restart_file",
                                     r_id,
                                     CS_RESTART_1D_WALL_THERMAL);

  }

  /* stop the calculation if needed once all checks have been done */
  cs_parameters_error_barrier();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
