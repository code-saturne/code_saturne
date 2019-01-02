#ifndef __CS_FLAG_CHECK_H__
#define __CS_FLAG_CHECK_H__

/*============================================================================
 * Mesh element flag checking and error handling.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for and handle errors with an associated element flag
 *
 * It is assumed that element flags are usually positive integers, and that
 * in case of a detected error, their signs have been set to a negative value.
 *
 * A minimum allowed value may be specified, so for example 0 may be
 * considered a valid or invalid flag depending on that minimum.
 *
 * This function exits silently if no such marked elements are present in the
 * computational domain.
 *
 * Otherwise, it logs information on the first detected error location, and
 * outputs postprocessing visualization information to assist debugging.
 *
 * If the error status (i.e. negative flag) is known locally but not
 * globally, use \ref cs_flag_check.
 *
 * Currently supported locations are CS_MESH_LOCATION_CELLS and
 * CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * \param[in]  err_elt_descr    description fro first element with error
 * \param[in]  flag_descr       flag type description
 * \param[in]  flag_label       field label for flag postprocessing
 * \param[in]  error_mesh_name  postprocessing mesh name for elements with error
 * \param[in]  valid_mesh_name  postprocessing mesh name for valid elements
 * \param[in]  location_id      associated mesh location
 * \param[in]  min_flag         minimum allowed flag
 * \param[in]  elt_flag         current element flag
 *
 * \return 0 in case no error was detected, 1 in case of errors.
 */
/*----------------------------------------------------------------------------*/

int
cs_flag_check(const char   *err_elt_descr,
              const char   *flag_descr,
              const char   *flag_label,
              const char   *error_mesh_name,
              const char   *valid_mesh_name,
              int           location_id,
              int           min_flag,
              const int     elt_flag[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle an error with an associated element flag
 *
 * This function logs information on the first detected error location, and
 * outputs postprocessing visualization information to assist debugging.
 *
 * It is assumed that element flags are usually positive integers, and that
 * in case of a detected error, their signs have been set to a negative value.
 *
 * A minimum allowed value may be specified, so for example 0 may be
 * considered a valid or invalid flag depending on that minimum.
 *
 * This function should be called when the error status has been previously
 * checked, and all ranks know that an error is present.
 *
 * If the error status (i.e. negative flag) is known locally but not
 * globally, use \ref cs_flag_check.
 *
 * Currently supported locations are CS_MESH_LOCATION_CELLS and
 * CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * \param[in]  err_elt_descr    description fro first element with error
 * \param[in]  flag_descr       flag type description
 * \param[in]  flag_label       field label for flag postprocessing
 * \param[in]  error_mesh_name  postprocessing mesh name for elements with error
 * \param[in]  valid_mesh_name  postprocessing mesh name for valid elements
 * \param[in]  location_id      associated mesh location
 * \param[in]  min_flag         minimum allowed flag
 * \param[in]  elt_flag         current element flag
 */
/*----------------------------------------------------------------------------*/

void
cs_flag_check_error_info(const char   *err_elt_descr,
                         const char   *flag_descr,
                         const char   *flag_label,
                         const char   *error_mesh_name,
                         const char   *valid_mesh_name,
                         int           location_id,
                         int           min_flag,
                         const int     elt_flag[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FLAG_CHECK_H__ */
