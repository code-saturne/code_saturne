/*============================================================================
 * Check computation parameters after user modification.
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
#include "cs_ale.h"
#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_ibm.h"
#include "cs_lagr.h"
#include "cs_les_balance.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_porosity_from_scan.h"
#include "cs_rad_transfer.h"
#include "cs_restart_default.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_syr_coupling.h"
#include "cs_wall_functions.h"
#include "cs_convection_diffusion.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"

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
  \var CS_ABORT_IMMEDIATE
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
 *
 * The caller is responsible for freeing the returned array.
 *----------------------------------------------------------------------------*/

inline static char *
_field_section_desc(cs_field_t  *f,
                    const char  *section_desc_b)
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
  int n_errors = _param_check_errors;
  cs_parall_sum(1, CS_INT_TYPE, &n_errors);

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
  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kclvfl = cs_field_key_id("variance_clipping");
  const int keyvar = cs_field_key_id("variable_id");
  const int kcpsyr = cs_field_key_id("syrthes_coupling");
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisls0 = cs_field_key_id("diffusivity_ref");
  const int restart_file_key_id = cs_field_key_id("restart_file");
  const int key_limiter = cs_field_key_id("limiter_choice");

  if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_ONLY)
    return; /* Avoid the detection of false setting errors when using
               CDO schemes */

  cs_field_t *f_pot = NULL;
  if (cs_glob_physical_model_flag[CS_GROUNDWATER] > 0) {
    f_pot = CS_F_(head);
    if (cs_glob_velocity_pressure_param->iphydr != 0) {
      cs_velocity_pressure_param_t *_vp_param
        = cs_get_glob_velocity_pressure_param();
      _vp_param->iphydr = 0;
    }
  }
  else
    f_pot = CS_F_(p);

  const cs_velocity_pressure_model_t *vp_model
    = cs_glob_velocity_pressure_model;
  cs_velocity_pressure_param_t *vp_param
    = cs_get_glob_velocity_pressure_param();

  cs_field_t *f_th = cs_thermal_model_field();

  char *f_desc = NULL;

  int list_01[2] = {0, 1};

  /*--------------------------------------------------------------------------
   * log frequency and checkpoint options
   *--------------------------------------------------------------------------*/

  cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_log_frequency",
                               cs_glob_log_frequency,
                               -1);

  cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_restart_auxiliary->read_auxiliary",
                               cs_glob_restart_auxiliary->read_auxiliary,
                               2,
                               list_01,
                               NULL);

  cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_restart_auxiliary->write_auxiliary",
                               cs_glob_restart_auxiliary->write_auxiliary,
                               2,
                               list_01,
                               NULL);

  /*--------------------------------------------------------------
   * Error estimators
   *--------------------------------------------------------------*/

  int n_ns_error_estimators = 0;

  {
    const char *name[] = {"est_error_pre_2",
                          "est_error_der_2",
                          "est_error_cor_2",
                          "est_error_tot_2"};

    for (int i = 0; i < 4; i++) {
      const cs_field_t *f = cs_field_by_name_try(name[i]);
      if (f != NULL)
        n_ns_error_estimators += 1;
    }

    if (n_ns_error_estimators > 0) {
      const char *ee_active
        = N_("One or several error estimates are activated for Navier-Stokes");

      if (cs_glob_time_scheme->iccvfg == 1)
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("while reading input data"),
           _("%s\n"
             "with frozen velocity field."), _(ee_active));

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > 0)
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("while reading input data"),
           _("%s\n"
             "this is not compatible with the compressible flow model."),
           _(ee_active));

      if (vp_param->nterup > 1)
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("%s\n"
             "this is not compatible with sub-iterations\n"
             "(here cs_glob_velocity_pressure_param->nterup = %d."),
           _(ee_active), vp_param->nterup);
    }
  }

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

  const int icp = cs_field_id_by_name("specific_heat");
  if (   cs_glob_physical_model_flag[CS_COOLING_TOWERS] > 0
        && icp == -1)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("while reading input data"),
         _("Cooling towers model requires variable specific_heat field.\n"));

  /* Equations definition, time scheme, convective scheme */
  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      f_desc = _field_section_desc(f, "while reading numerical "
                                      "parameters for variable");

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param iconv (convection flag)",
                                    eqp->iconv,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param istat (unsteadiness flag)",
                                    eqp->istat,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param idircl (reinforce matrix diag)",
                                    eqp->idircl,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param idiff (diffusion flag)",
                                    eqp->idiff,
                                    0, 2);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param idifft (turbulent diffusion "
                                                                       "flag)",
                                    eqp->idifft,
                                    0, 2);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "equation param thetav (theta-scheme)",
                                       eqp->thetav,
                                       0., 1.);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "equation param blencv (2nd order scheme "
                                       "share for convection)",
                                       eqp->blencv,
                                       0., 1.);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "equation param blend_st (2nd order scheme "
                                       "share for convection)",
                                       eqp->blend_st,
                                       0., 1.);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param ischcv (2nd order scheme "
                                    "choice)",
                                    eqp->ischcv,
                                    0, 5);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param isstpc (limiter type)",
                                    eqp->isstpc,
                                    0, 3);

      BFT_FREE(f_desc);
    }
  }

  /* check if NVD scheme for thermal scalar is not one of the VOF schemes */
  if (f_th != NULL) {
    cs_equation_param_t *eqp = cs_field_get_equation_param(f_th);
    if (eqp->ischcv >= 4) { /* NVD scheme on thermal scalar? */
      cs_nvd_type_t limiter_choice = cs_field_get_key_int(f_th, key_limiter);

      f_desc = _field_section_desc(f_th, "while reading numerical "
                                         "parameters for variable");

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "NVD scheme",
                                    limiter_choice,
                                    CS_NVD_GAMMA, CS_NVD_VOF_HRIC);

      BFT_FREE(f_desc);
    }
  }

  /* thetav for pressure must be equal to 1 */
  {
    cs_equation_param_t *eqp = cs_field_get_equation_param(f_pot);
    f_desc = _field_section_desc(f_pot, "while reading numerical "
                                        "parameters for variable");

    cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                  _(f_desc),
                                  "equation param thetav (theta-scheme)",
                                  eqp->thetav,
                                  1.);

    BFT_FREE(f_desc);
  }

  /* 2nd order in time (rho, visc, N.S source terms, velocity theta)
     is assumed to be incompatible with:
     - error estimateurs
     - ipucou
     - iphydr = 2
     - local or variable time step */

  if (CS_F_(vel) != NULL) {
    int key_t_ext_id = cs_field_key_id("time_extrapolated");
    cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(vel));

    const char *tds_err_str
      = N_("Some options are incompatible with the time discretization scheme\n"
           "\n"
           " A second order time-scheme was requested:\n"
           "\n"
           " Velocity:                    theta = %g\n"
           " Navier-Stokes source terms:  isno2t = %d\n"
           "                              thetsn = %g\n"
           " Density:                     iroext = %d\n"
           " Viscosity:                   iviext = %d\n"
           "                              thetvi = %d\n"
           "\n"
           "This is not compatible with:\n"
           "- error estimators (%d)\n"
           "- reinforced U-P coupling (ipucou): %d\n"
           "- non-constant time step (idtvar): %d.");

    int iroext = 0, iviext = 0;
    if (CS_F_(rho) != NULL)
      iroext = cs_field_get_key_int(CS_F_(rho), key_t_ext_id);
    if (CS_F_(mu) != NULL)
      iviext = cs_field_get_key_int(CS_F_(mu), key_t_ext_id);

    const cs_time_scheme_t *t_sch = cs_glob_time_scheme;

    if (   fabs(eqp->thetav-1.0) > 1e-3
        || t_sch->thetvi > 0
        || t_sch->thetsn > 0
        || t_sch->isno2t > 0
        || iroext > 0
        || iviext > 0) {

      if (   n_ns_error_estimators > 0
          || vp_param->ipucou == 1
          || vp_param->iphydr == 2
          || cs_glob_time_step_options->idtvar != 0)
        cs_parameters_error(CS_ABORT_DELAYED,
                            _("while reading input data"),
                            _(tds_err_str),
                            eqp->theta,
                            t_sch->isno2t,
                            t_sch->thetsn,
                            iroext,
                            iviext,
                            t_sch->thetvi,
                            n_ns_error_estimators,
                            vp_param->ipucou,
                            cs_glob_time_step_options->idtvar);
    }

    if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] > 0
        && vp_param->nterup > 1)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("while reading input data"),
         _("Pressure-Velocity coupling with sub-iterations\n"
           "(cs_glob_velocity_pressure_param->nterup = %d.\n"
           "is not compatible with the compressible flow model."),
         vp_param->nterup);
  }

  /* In LES, additional consistency checkings are needed *
   * Only a warning for non standard parameters, but stop if
   * there is more than 5% of upwind.
   * Centered scheme with/without slope test, nwsrsm */
  const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();
  assert(turb_model != NULL);

  if (turb_model->type == CS_TURB_LES) {
    cs_equation_param_t *eqp_v = cs_field_get_equation_param(CS_F_(vel));
    f_desc = _field_section_desc(CS_F_(vel), "in LES, while reading time "
                                 "scheme parameters for variable");

    cs_parameters_is_equal_double(CS_WARNING,
                                  _(f_desc),
                                  "equation param thetav (theta-scheme)",
                                  eqp_v->thetav,
                                  0.5);

    BFT_FREE(f_desc);

    f_desc = _field_section_desc(CS_F_(vel), "in LES, while reading "
                                 "convection scheme "
                                 "parameters for variable");

    cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                     _(f_desc),
                                     "equation param blencv (2nd order scheme "
                                     "share for convection)",
                                     eqp_v->blencv,
                                     0.95, 1.);

    cs_parameters_is_equal_double(CS_WARNING,
                                  _(f_desc),
                                  "equation param blencv (2nd order scheme "
                                  "share for convection)",
                                  eqp_v->blencv,
                                  1.);

    cs_parameters_is_equal_int(CS_WARNING,
                               _(f_desc),
                               "equation param isstpc (limiter type)",
                               eqp_v->isstpc,
                               1);

    BFT_FREE(f_desc);

    for (int f_id = 0 ; f_id < n_fields ; f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int isca = cs_field_get_key_int(f, keysca);
      if (isca > 0) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        f_desc = _field_section_desc(f, "in LES, while reading time "
                                        "scheme parameters for variable");

        cs_parameters_is_equal_double(CS_WARNING,
                                      _(f_desc),
                                      "equation param thetav (theta-scheme)",
                                      eqp->thetav,
                                      0.5);

        BFT_FREE(f_desc);

        f_desc = _field_section_desc(f, "in LES, while reading "
                                        "convection scheme "
                                        "parameters for variable");

        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _(f_desc),
                                         "equation param blencv (2nd order "
                                           "scheme share for convection)",
                                         eqp->blencv,
                                         0.95, 1.);

        cs_parameters_is_equal_double(CS_WARNING,
                                      _(f_desc),
                                      "equation param blencv (2nd order scheme "
                                                      "share for convection)",
                                      eqp->blencv,
                                      1.);

        cs_parameters_is_equal_int(CS_WARNING,
                                   _(f_desc),
                                   "equation param isstpc (limiter type)",
                                   eqp->isstpc,
                                   0);

        BFT_FREE(f_desc);
      }
    }

    cs_parameters_is_equal_int(CS_WARNING,
                               _("Checkpoint settings with LES computation."),
                               "cs_glob_restart_auxiliary->read_auxiliary",
                               cs_glob_restart_auxiliary->read_auxiliary,
                               1);

    cs_parameters_is_equal_int(CS_WARNING,
                               _("Checkpoint settings with LES computation."),
                               "cs_glob_restart_auxiliary->write_auxiliary",
                               cs_glob_restart_auxiliary->write_auxiliary,
                               1);
  }

  /* navsto sub-iterations
   * It must be a integer greater or equal to 1
   * For the moment, we forbid nterup > 1 with
   * estimators, weight matrix (reinforced U-P coupling), hydrostatic pressure
   * and steady algorithm. */
  cs_parameters_is_positive_int
    (CS_ABORT_DELAYED,
     _("while reading input data"),
     "cs_glob_velocity_pressure_param->nterup (Navier-Stokes sub-iterations)",
     vp_param->nterup);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_time_step_options->idtvar",
                                cs_glob_time_step_options->idtvar,
                                -1, 3);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_time_step_options->itpcol",
                                vp_param->itpcol,
                                -1, 2);

  if (cs_glob_time_scheme->iccvfg == 1) {
    cs_parameters_is_equal_int(CS_WARNING,
                              _("while reading input data,\n"
                                "cs_glob_velocity_pressure_param->nterup "
                                "(Navier-Stokes sub-iterations)\n"
                                "is incompatible with fozen velocity field\n"
                                "it will be set to 1."),
                              "nterup",
                              vp_param->nterup,
                              1);
    vp_param->nterup = 1;
  }

  if (vp_param->ipucou == 1) {
    cs_parameters_is_equal_int(CS_WARNING,
                              _("while reading input data,\n"
                                "cs_glob_velocity_pressure_param->nterup "
                                "(Navier-Stokes sub-iterations)\n"
                                "is not compatible with reinforced "
                                "velocity-pressure coupling (ipucou=1)\n"
                                "it will be set to 1."),
                              "nterup",
                              vp_param->nterup,
                              1);
    vp_param->nterup = 1;
  }

  if (cs_glob_time_step_options->idtvar == -1) {
    cs_parameters_is_equal_int(CS_WARNING,
                              _("while reading input data,\n"
                                "cs_glob_velocity_pressure_param->nterup "
                                "(Navier-Stokes sub-iterations)\n"
                                "is not compatible with steady algorithm "
                                "(idtvar=-1)\n"
                                "it will be set to 1."),
                              "nterup",
                              vp_param->nterup,
                              1);
    vp_param->nterup = 1;
  }

  /* Steady Algorithm */
  if (cs_glob_time_step_options->idtvar < 0) {
    for (int f_id = 0 ; f_id < n_fields ; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int ivar = cs_field_get_key_int(f, keyvar);
      if (ivar > 0) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        f_desc = _field_section_desc(f, "With steady algorithm (SIMPLE), while "
                                        "reading numerical parameters for "
                                        "variable");

        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _(f_desc),
                                         "equation param relaxv (relax. coef.)",
                                         eqp->relaxv,
                                         0., 1.);

        cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                      _(f_desc),
                                      "equation param isstpc (limiter type)",
                                      eqp->isstpc,
                                      0, 2);

        cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                      _(f_desc),
                                      "equation param ischcv",
                                      eqp->ischcv,
                                      0, 3);

        const int kiflux = cs_field_key_id("inner_flux_id");
        const int kbflux = cs_field_key_id("boundary_flux_id");
        int i_flux_id = cs_field_get_key_int(f, kiflux);
        int b_flux_id = cs_field_get_key_int(f, kbflux);

        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "key inner_flux_id (inner flux field id)",
                                   i_flux_id,
                                   -1);

        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "key boundary_flux_id (boundary flux field id)",
                                   b_flux_id,
                                   -1);

        BFT_FREE(f_desc);
      }
    }

    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("The steady algorithm is not compatible "
                                 "with the Lagrangian module which is "
                                 "time-dependant by nature"),
                               "cs_glob_lagr_time_scheme->iilagr",
                               cs_glob_lagr_time_scheme->iilagr,
                               CS_LAGR_OFF);

    int les_iturb[3]
      = {CS_TURB_LES_SMAGO_CONST, CS_TURB_LES_SMAGO_DYN, CS_TURB_LES_WALE};
    cs_parameters_is_not_in_list_int(CS_ABORT_DELAYED,
                                     _("The steady algorithm is not compatible "
                                       "with L.E.S. turbulence modelling "
                                       "which is time-dependant by nature"),
                                     "cs_glob_turb_model->iturb",
                                     turb_model->iturb,
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

  /* Numbers of sweeps don't need to be checked: they are simply
   * integers (negative if one wants to be sure to never enter the loops */
  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int ivar = cs_field_get_key_int(f, keyvar);
    if (ivar > 0) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      f_desc = _field_section_desc(f, "Wile reading numerical parameters "
                                      " for variable");

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param imligr "
                                    "(gradient limitation method)",
                                    eqp->imligr,
                                    -1, 2);

      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _(f_desc),
                                      "equation param climgr "
                                      "(gradient limitation coeffcient)",
                                      eqp->climgr,
                                      1.);

      cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                    _(f_desc),
                                    "equation param ircflu (fluxes "
                                    "reconstruction)",
                                    eqp->ircflu,
                                    0, 2);

      BFT_FREE(f_desc);

      if (   eqp->ischcv == 0
          && CS_ABS(eqp->blencv) > cs_math_epzero) {
        f_desc = _field_section_desc(f, "Second order linear upwind "
                                        "enabled for variable ");

        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "equation param ircflu (fluxes reconstruction)",
                                   eqp->ircflu,
                                   1);

        BFT_FREE(f_desc);
      }
    }
  }

  /* Time stepping */
  /* The number of time step could be negative : no test */

  if (cs_glob_time_step_options->idtvar >= 0) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_time_step->dt_ref",
                                    cs_glob_time_step->dt_ref,
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
                                0, 3);

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
  }
  else if (cs_glob_time_step_options->idtvar == -1) {
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
  const int iturb_vals[16] = {CS_TURB_NONE,   /* laminar */
                              CS_TURB_MIXING_LENGTH,
                              CS_TURB_K_EPSILON,
                              CS_TURB_K_EPSILON_LIN_PROD,
                              CS_TURB_K_EPSILON_LS,
                              CS_TURB_K_EPSILON_QUAD,
                              CS_TURB_RIJ_EPSILON_LRR,
                              CS_TURB_RIJ_EPSILON_SSG,
                              CS_TURB_RIJ_EPSILON_EBRSM,
                              CS_TURB_LES_SMAGO_CONST,
                              CS_TURB_LES_SMAGO_DYN,
                              CS_TURB_LES_WALE,
                              CS_TURB_V2F_PHI,
                              CS_TURB_V2F_BL_V2K,
                              CS_TURB_K_OMEGA,
                              CS_TURB_SPALART_ALLMARAS};

  cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                               _("while reading input data"),
                               "cs_glob_turb_model->iturb",
                               turb_model->iturb,
                               16,
                               iturb_vals,
                               NULL);

  /* Rotation curvature correction for eddy-viscosity models */
  assert(cs_glob_turb_rans_model != NULL);
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_turb_rans_model->irccor",
                                cs_glob_turb_rans_model->irccor,
                                0, 2);

  /* Rotation curvature correction compatible only with RANS eddy-viscosity
     models */
  if (cs_glob_turb_rans_model->irccor == 1) {
    const int iturb_evm_vals[6] = {CS_TURB_K_EPSILON,
                                   CS_TURB_K_EPSILON_LIN_PROD,
                                   CS_TURB_V2F_PHI,
                                   CS_TURB_V2F_BL_V2K,
                                   CS_TURB_K_OMEGA,
                                   CS_TURB_SPALART_ALLMARAS};

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "rotation curvature correction is only "
                                   "compatible with eddy viscosity turbulence "
                                   "models"),
                                 "cs_glob_turb_model->iturb",
                                 turb_model->iturb,
                                 6,
                                 iturb_evm_vals,
                                 NULL);
  }

  /* Hybrid RANS/LES models (DES, DDES, SAS) only with k-omega SST model */
  if (   turb_model->hybrid_turb > 0
      && turb_model->hybrid_turb < 4) {
    const int iturb_ddes_vals[1] = {CS_TURB_K_OMEGA};

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "Hybrid RANS/LES model is only "
                                   "compatible with k-omega SST model "
                                   "(iturb=60)"),
                                 "cs_glob_turb_model->iturb",
                                 turb_model->iturb,
                                 1,
                                 iturb_ddes_vals,
                                 NULL);

  }

  /* HTLES model only with k-omega SST and BL-v2/k model */
  if (turb_model->hybrid_turb == 4) {
    const int iturb_htles_vals[2] = {CS_TURB_K_OMEGA, CS_TURB_V2F_BL_V2K};

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "HTLES model is only compatible"
                                   "with k-omega SST model (iturb=60)"
                                   "and BL-v2/k model (iturb=51)"),
                                 "cs_glob_turb_model->iturb",
                                 turb_model->iturb,
                                 2,
                                 iturb_htles_vals,
                                 NULL);

  }

  /* In Lagrangian with two-way coupling, k-omega SST is forbidden (not
     properly implemented) */
  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {
    cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "two way coupling in Lagrangian modelling "
                                     "is not compatible with k-omega SST "
                                     "turbulence model"),
                                   "cs_glob_turb_model->iturb",
                                   turb_model->iturb,
                                   CS_TURB_K_OMEGA);
  }

  if (cs_glob_lagr_time_scheme->extended_t_scheme == 1) {

    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("while reading input data,\n"
                                 "extended scheme in Lagrangian modelling "
                                 "is only compatible with "
                                 "spherical particles"),
                               "cs_glob_lagr_model->shape",
                               cs_glob_lagr_model->shape,
                               CS_LAGR_SHAPE_SPHERE_MODEL);

    cs_parameters_is_equal_int(CS_WARNING,
                               _("while reading input data,\n"
                                 "extended scheme in Lagrangian modelling "
                                 "is usefull only if turbulent dispersion "
                                 "is considered"),
                               "cs_glob_lagr_model->idistu",
                               cs_glob_lagr_model->idistu,
                               1);
  }

  /* LES balance */
  if (   cs_glob_turb_model->type != CS_TURB_LES
      && cs_glob_les_balance->i_les_balance != 0) {
    cs_parameters_is_equal_int(CS_WARNING,
                               _("while reading input data,\n"
                                 "LES balance only for LES , "
                                 "this setting will be ignored"),
                               "cs_glob_les_balance->i_les_balance",
                               cs_glob_les_balance->i_les_balance,
                               0);

    cs_les_balance_t *les_balance = cs_get_glob_les_balance();
    les_balance->i_les_balance = 0;
  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_wall_functions->iwallf",
                                cs_glob_wall_functions->iwallf,
                                0, 8);

  if (cs_glob_wall_functions->iwallf > 2) {
    const int itytur_vals[5] = {2, 3, 5, 6, 7};

    cs_parameters_is_in_list_int(CS_ABORT_DELAYED,
                                 _("while reading input data,\n"
                                   "two scales wall function not compatible "
                                   "with a laminar, mixing length or LES "
                                   "computation"),
                                 "cs_glob_turb_model->itytur",
                                 turb_model->itytur,
                                 5,
                                 itytur_vals,
                                 NULL);
  }

  /* Specific k-epsilon, v2f and k-omega */
  if (   turb_model->itytur == 2
      || turb_model->itytur == 5
      || turb_model->itytur == 6) {
    /* iclkep option not available in k-omega */
    if (turb_model->iturb != CS_TURB_K_OMEGA) {
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

    /* In k-eps with lin prod/LS/Quad and in v2f we force ikecou to 0 */
    if (   turb_model->iturb == CS_TURB_K_EPSILON_LIN_PROD
        || turb_model->iturb == CS_TURB_K_EPSILON_LS
        || turb_model->iturb == CS_TURB_K_EPSILON_QUAD
        || turb_model->itytur == 5) {
      cs_parameters_is_equal_int
        (CS_ABORT_DELAYED,
         _("while reading input data,\n"
           "with k-epsilon LP (iturb=CS_TURB_K_EPSILON_LIN_PROD),\n"
           "k-epsilon LS (iturb=CS_TURB_K_EPSILON_LS),\n"
           "k-epislon quadratic (iturb=CS_TURB_K_EPSILON_QUAD),\n"
           "or v2f model (iturb=CS_TURB_V2F_PHI, CS_TURB_V2F_BL_V2K)."),
         "cs_glob_turb_rans_model->ikecou",
         cs_glob_turb_rans_model->ikecou,
         0);
    }

    /* In steady mode force IEKCOU to 0 */
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
  }

  /* if relaxv of k, epsilon or omega was modified by the user, and if ikecou=1,
     we warn the user that relaxv settings won't have any effect,
     otherwise check that relaxv is in range [O,1] (already done in steady) */

  /* relaxv takes the value 1. in modini if not modified by the user
     if idtvar >= 0 */

  if (   (   turb_model->itytur == 2
          || turb_model->itytur == 6)
      && cs_glob_time_step_options->idtvar >= 0) {
    cs_field_t *f_eo = (turb_model->itytur == 2) ? CS_F_(eps) : CS_F_(omg);
    int f_ids[2] = {CS_F_(k)->id, f_eo->id};

    for (int ii = 0; ii < 2; ii++) {
      cs_field_t *f = cs_field_by_id(f_ids[ii]);
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      if (cs_glob_turb_rans_model->ikecou == 1) {
        f_desc = _field_section_desc(f,
                                     "while reading input data,\n"
                                     "modifications of relaxation "
                                     "coefficients will be ignored with ikecou=1 "
                                     " for variable");

        cs_parameters_is_equal_double(CS_WARNING,
                                      _(f_desc),
                                      "equation param relaxv",
                                      eqp->relaxv,
                                      1.);
        BFT_FREE(f_desc);
      } else { /* ikecou = 0 */
        f_desc = _field_section_desc(f, "while reading numerical "
                                        "parameters for variable");

        cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                         _(f_desc),
                                         "equation param relaxv",
                                         eqp->relaxv,
                                         0, 1);
        BFT_FREE(f_desc);
      }
    }
  }

  /* Check that relaxv is in [0,1] for Spallart-Allmaras nu variable
     (already done for steady) */
  if (   turb_model->iturb == CS_TURB_SPALART_ALLMARAS
      && cs_glob_time_step_options->idtvar >= 0) {
      cs_field_t *f = CS_F_(nusa);
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      f_desc = _field_section_desc(f, "while reading numerical "
                                      "parameters for variable");

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _(f_desc),
                                       "equation param relaxv",
                                       eqp->relaxv,
                                       0, 1);
      BFT_FREE(f_desc);
  }

  /* checks for RSM models */
  if (turb_model->order == CS_TURB_SECOND_ORDER) {
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->irijnu",
                                  cs_glob_turb_rans_model->irijnu,
                                  0, 3);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_rans_model->irijrb",
                                  cs_glob_turb_rans_model->irijrb,
                                  0, 2);

    /* wall echo and specific implicitation of the diffusion of epsilon only
       in Rij LRR */
    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
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
  if (turb_model->type == CS_TURB_LES) {
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_turb_les_model->idries",
                                  cs_glob_turb_les_model->idries,
                                  0, 2);


    if (   turb_model->iturb == CS_TURB_LES_SMAGO_DYN
        || turb_model->iturb == CS_TURB_LES_WALE) {
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
    if (turb_model->iturb == CS_TURB_LES_SMAGO_DYN) {
      int imrgra_cmp = CS_ABS(turb_model->iturb);
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
             "dynamic Smagorinsky constant via the\n"
             "cs_user_physical_properties_turb_viscosity function."));
        break;
      default:
        break;
      }
    }
  }

  /* Stokes */

  /* transposed velocity gradient term and secondary viscosity */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_velocity_pressure_model->ivisse",
                                vp_model->ivisse,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_velocity_pressure_param->iprco",
                                vp_param->iprco,
                                0, 2);

  if (vp_param->iprco == 1) {
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_velocity_pressure_param->irevmc",
                                  vp_param->irevmc,
                                  0, 2);

    cs_real_t arakfr = vp_param->arak;
    if (cs_glob_time_step_options->idtvar < 0) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(vel));
      arakfr *= eqp->relaxv;
    }

    cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                     _("while reading input data"),
                                     "cs_glob_velocity_pressure_param->arak",
                                     arakfr,
                                     0., 1.);
  }

  /* U-P reinforced coupling */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_velocity_pressure_param->ipucou",
                                vp_param->ipucou,
                                0, 2);

  /* Dilatable algorithm: 0 to 5 */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_velocity_pressure_model->idilat",
                                vp_model->idilat,
                                0, 6);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_velocity_pressure_param->iphydr",
                                vp_param->iphydr,
                                0, 3);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_porous_model",
                                cs_glob_porous_model,
                                0, 4);

  if (cs_glob_porosity_from_scan_opt->compute_porosity_from_scan)
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data,\n"
                                    "Porosity from scan has been activated"
                                    "and cs_glob_porous_model changed\n"),
                                  "cs_glob_porous_model",
                                  cs_glob_porous_model,
                                  3, 4);

  if (cs_glob_porosity_ibm_opt->porosity_mode > 0)
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data,\n"
                                    "Porosity from immersed boundaries has"
                                    "been activated (experimental) and "
                                    "cs_glob_porous_model changed\n"),
                                  "cs_glob_porous_model",
                                  cs_glob_porous_model,
                                  3, 4);

  if (cs_glob_porous_model == 3) { /* integral formulation */
    const int imrgra_vals[5] ={0, 4, 5, 6, 7};

    cs_parameters_is_in_list_int(CS_WARNING,
                               _("while reading porous model,\n"
                                 "integral formulation "
                                 "(cs_glob_porous_model=3) "
                                 "not compatible \n"
                                 "with gradient calculation method: "
                                 "least squares"),
                               "cs_glob_space_disc->imrgra",
                               cs_glob_space_disc->imrgra,
                               5,
                               imrgra_vals,
                               NULL);
  }

  cs_equation_param_t *eqp_v = cs_field_get_equation_param(CS_F_(vel));
  /* steady or variable time step time algorithm not compatible with theta
     scheme with theta different from 1 for the velocity */
  if (cs_glob_time_step_options->idtvar != 0) {
    cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                  _("while reading time scheme parameters,\n"
                                    "theta-scheme with theta different from 1 "
                                    "for the velocity\n"
                                    "only compatible with constant time step "
                                    "unsteady algorithm"),
                                  "equation param thetav",
                                  eqp_v->thetav,
                                  1);
  }
  /* U-P reinforced coupling not compatible with theta scheme with theta
     different from 1 for the velocity */
  if (vp_param->ipucou == 1) {
    cs_parameters_is_equal_double(CS_ABORT_DELAYED,
                                  _("while reading time scheme parameters,\n"
                                    "theta-scheme with theta different from 1 "
                                    "for the velocity\n"
                                    "not compatible with reinforced "
                                    "velocity-pressure coupling"),
                                  "equation param thetav",
                                  eqp_v->thetav,
                                  1);
  }

  /* frozen velocity field */
  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_time_scheme->iccvfg",
                                cs_glob_time_scheme->iccvfg,
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

  /* Check if there is coupling */
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
  }
  else { /* if coupling with SYRTHES */
    /* and more than one scalar is defined as coupled */
    cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                               _("Inconsistency in SYRTHES coupling settings,\n"
                                 "more than one "
                                 "scalars are defined as coupled"),
                               "number of coupled scalars",
                               nbsccp,
                               1);

    const cs_field_t *tf = cs_thermal_model_field();
    const char none_defined[] = "<none_defined>";
    const char *tf_name = none_defined;
    if (tf)
      tf_name = tf->name;

    /* and the coupled scalar is not the thermal scalar
       or if the thermal scalar is not the temperature
       if compressible is not enabled */
    for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
      cs_field_t  *f = cs_field_by_id(f_id);
      int icpsyr = cs_field_get_key_int(f, kcpsyr);
      if (icpsyr == 1 && f != tf) {
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("Inconsistency in SYRTHES coupling "
             "settings,\n"
             "the coupled scalar must be the "
             "the thermal scalar (%s), not %s."),
           tf_name, f->name);
      }
    }
  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("while reading input data"),
                                "cs_glob_space_disc->itbrrb",
                                cs_glob_space_disc->itbrrb,
                                0, 2);

  /* Dynamic relaxation option */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (eqp->iswdyn >= 1) {
        if (f->id == CS_F_(vel)->id) {
          cs_parameters_is_equal_int(CS_WARNING,
                                     _("Dynamic relaxation enabled for "
                                       "variable velocity.\n"
                                       "The slope test must be disabled."),
                                     "equation param isstpc",
                                     eqp->isstpc,
                                     1);

          if (eqp->isstpc == 0) {
            eqp->isstpc = 1;
            int log_id = CS_LOG_DEFAULT;
            cs_log_printf(log_id,
                          _("The calculation continues with isstpc = 1 for "
                            "the variable velocity.\n"));
          }
        }
      }
    }
  }

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      f_desc = _field_section_desc(f, "Number of RHS reconstruction "
                                      "sweeps for variable ");

      cs_parameters_is_positive_int(CS_WARNING,
                                    _(f_desc),
                                    "equation param nswrsm",
                                    eqp->nswrsm);

      BFT_FREE(f_desc);

      if (eqp->nswrsm <= 0) {
        eqp->nswrsm = 1;
        int log_id = CS_LOG_DEFAULT;
        cs_log_printf(log_id,
                      _("The calculation continues with nswrsm = 1 "
                        "for this variable.\n"));
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

  if (cs_glob_vof_parameters->vof_model > 0) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                  _("while reading reference surface "
                                    "tension value"),
                                  "cs_glob_vof_parameters->sigma_s",
                                  cs_glob_vof_parameters->sigma_s,
                                  0.);
  }

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
      int vr_clip = cs_field_get_key_int(f, kclvfl);

      if (fscavr_id > -1) {
        cs_field_t *fluct_f = cs_field_by_id(fscavr_id);

        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "associated fluctuating field is already a "
                                   "variance.\n Associated fluctuating field "
                                   "scalar id for this variance",
                                   cs_parameters_iscavr(fluct_f),
                                   0);

        cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                      _(f_desc),
                                      "variance_clipping",
                                      vr_clip,
                                      0, 3);

      }
      else
        cs_parameters_is_equal_int(CS_ABORT_DELAYED,
                                   _(f_desc),
                                   "variance_clipping",
                                   vr_clip,
                                   -1);

      BFT_FREE(f_desc);
    }
  }

  /* check positivity of reference diffusivities */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_real_t visls_0 = cs_field_get_key_double(f, kvisls0);
      f_desc = _field_section_desc(f, "reference diffusivity value for "
                                      "variable ");

      int diff_id = cs_field_get_key_int(f, kivisl);
      int isca = cs_field_get_key_int(f, keysca);
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      if (isca > 0 && diff_id == -1 && eqp->idiff > 0) {
        cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                        _(f_desc),
                                        "key diffusivity_ref",
                                        visls_0,
                                        0.);
      }

      BFT_FREE(f_desc);
    }
  }

  /* Turbulence */

  /* uref is needed for automatic initialisation and boundary conditions
     of turbulence variables check that it is positive
     and warn the user if not */
  if (   turb_model->type == CS_TURB_RANS
      && turb_model->order == CS_TURB_FIRST_ORDER) {
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

  if (turb_model->iturb == CS_TURB_MIXING_LENGTH) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _("while reading input data"),
                                    "cs_glob_turb_rans_model->xlomlg "
                                    "(mixing length)",
                                    cs_glob_turb_rans_model->xlomlg,
                                    0.);
  }

  /* LES (check also in Fortran for now because constants are duplicated) */
  if (turb_model->type == CS_TURB_LES) {
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

    if (turb_model->iturb == CS_TURB_LES_SMAGO_DYN) {
      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _("LES dynamic model constant"),
                                      "cs_turb_xlesfd",
                                      cs_turb_xlesfd,
                                      0.);

      cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                      _("LES dynamic model constant"),
                                      "cs_turb_csmago_max",
                                      cs_turb_csmago_max,
                                      0.);

      cs_parameters_is_in_range_double(CS_ABORT_DELAYED,
                                       _("LES dynamic model constant"),
                                       "cs_turb_csmago_min (must be < "
                                       "cs_turb_csmago_max)",
                                       cs_turb_csmago_min,
                                       0., cs_turb_csmago_max);
    }
  }

  /*--------------------------------------------------------------
   * ALE
   *--------------------------------------------------------------*/

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("flag for ALE method"),
                                "cs_glob_ale or iale in Fortran",
                                cs_glob_ale,
                                0, 3);

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

  /* Check the consistency of time extrapolation key */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);

    int t_ext = cs_field_get_key_int(f, cs_field_key_id("time_extrapolated"));

    f_desc = _field_section_desc(f, "while reading parameters for "
                                    "field ");
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _(f_desc),
                                  "key 'time_extrapolated'",
                                  t_ext,
                                  -1,
                                  3);

    BFT_FREE(f_desc);
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
                                     "the Lagrangian checkpoint cannot "
                                     "be used\n"
                                     "without having the Lagrangian module "
                                     "enabled."),
                                     "restart_file",
                                     r_id,
                                     CS_RESTART_LAGR_STAT);

      cs_parameters_is_not_equal_int(CS_ABORT_DELAYED,
                                   _("while reading input data,\n"
                                     "the Lagrangian checkpoint cannot "
                                     "be used\n"
                                     "without having the Lagrangian module "
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
