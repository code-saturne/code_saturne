#ifndef __CS_PARAMETERS_CHECK_H__
#define __CS_PARAMETERS_CHECK_H__

/*============================================================================
 * Log field and other array statistics at relevant time steps.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * type definitions
 *============================================================================*/

/* Parameter check behavior when an error is detected */

typedef enum {

  CS_WARNING,
  CS_ABORT_DELAYED,
  CS_ABORT_IMMEDIATE

} cs_parameter_error_behavior_t;

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

#if defined(__GNUC__)

void
cs_parameters_error(cs_parameter_error_behavior_t   err_behavior,
                    const char                     *section_desc,
                    const char                     *format,
                    ...)
  __attribute__((format(printf, 3, 4)));

#else

void
cs_parameters_error(cs_parameter_error_behavior_t   err_behavior,
                    const char                     *section_desc,
                    const char                     *format,
                    ...);

#endif

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
                    ...);

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
                           const char                     *section_desc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print footer for a given parameters error message type.
 *
 * \param[in]  err_behavior  warn or abort ?
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_footer(cs_parameter_error_behavior_t   err_behavior);

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
                              int                             range_u);

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
                                  int                             range_u);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword has values in a specified range.
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
                             const char                     *enum_names[]);

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
                                 const char                     *enum_names[]);

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
                           int                             std_value);

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
                               int                             fbd_value);

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
                              int                             param_value);

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
                             int                             ib_value);

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
                                 double                          range_u);

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
                                const char                     *enum_names[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given double keyword is equal to a specified value.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  std_value     compulsory parameter's value
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_equal_double(cs_parameter_error_behavior_t   err_behavior,
                              const char                     *section_desc,
                              const char                     *param_name,
                              double                          param_value,
                              double                          std_value);

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
                                double                          ib_value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Abort if the the parameter errors count is nonzero.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_barrier(void);

/*----------------------------------------------------------------------------
!
 * \brief Check data settings
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_check(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMETERS_CHECK_H__ */
