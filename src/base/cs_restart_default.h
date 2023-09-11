#ifndef __CS_RESTART_DEFAULT_H__
#define __CS_RESTART_DEFAULT_H__

/*============================================================================
 * Checkpoint/restart handling for default application.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"
#include "cs_map.h"
#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Restart file ids */
/*-------------------*/

typedef enum {

  CS_RESTART_DISABLED = -1,       /*!< no values to save */
  CS_RESTART_MAIN = 0,            /*!< save values in main restart file */
  CS_RESTART_AUXILIARY = 1,       /*!< save values in auxiliary restart file */
  CS_RESTART_RAD_TRANSFER = 2,    /*!< save values in radiative transfer
                                       restart file */
  CS_RESTART_LAGR = 3,            /*!< save values in lagrangian restart file */
  CS_RESTART_LAGR_STAT = 4,       /*!< save values in restart file for
                                       lagrangian statistics */
  CS_RESTART_1D_WALL_THERMAL = 5, /*!< save values in 1D wall thermal restart
                                       file */
  CS_RESTART_LES_INFLOW = 6,      /*!< save values in LES inflow restart file */
  CS_RESTART_N_RESTART_FILES = 7  /*!< Number of types of restart file */

} cs_restart_file_t ;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read field metadata from checkpoint.
 *
 * Old ids associated to each field are determined for future use.
 * Note that when reading legacy files (code_saturne version 3.3 and below),
 * the old id will actually be the old scalar id (-1 for others).
 *
 * parameters:
 *   r             <-> associated restart file pointer
 *   old_field_map --> name to id map of fields in restart file
 *----------------------------------------------------------------------------*/

void
cs_restart_read_field_info(cs_restart_t          *r,
                           cs_map_name_to_id_t  **old_field_map);

/*----------------------------------------------------------------------------
 * Write field metadata to checkpoint.
 *
 * parameters:
 *   r <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_write_field_info(cs_restart_t  *r);

/*----------------------------------------------------------------------------
 * Read variables from checkpoint.
 *
 * parameters:
 *   r             <-> associated restart file pointer
 *   old_field_map <-- name to id map of fields in restart file
 *   t_id_flag     <-- -1: all time values; 0: current values;
 *                      > 0: previous values
 *   read_flag     <-- optional flag to track fields read, or NULL;
 *                     set to sum of 2^time_id for fields read (size: n_fields)
 *----------------------------------------------------------------------------*/

void
cs_restart_read_variables(cs_restart_t               *r,
                          const cs_map_name_to_id_t  *old_field_map,
                          int                         t_id_flag,
                          int                         read_flag[]);

/*----------------------------------------------------------------------------
 * Write variables to checkpoint.
 *
 * parameters:
 *   r          <-> associated restart file pointer
 *   t_id_flag  <-- -1: all time values; 0: current values;
 *                  > 0: previous values
 *   write_flag <-- optional flag to track fields written, or NULL;
 *                  set to sum of 2^time_id for fields written (size: n_fields)
*----------------------------------------------------------------------------*/

void
cs_restart_write_variables(cs_restart_t  *r,
                           int            t_id_flag,
                           int            write_flag[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read notebook parameters from main checkpoint.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_notebook_variables(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write notebook parameters to main checkpoint.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_notebook_variables(cs_restart_t  *r);

/*----------------------------------------------------------------------------
 * Read fields depending on others from checkpoint.
 *
 * Old ids associate to each field are determined for future use.
 * Note that when reading legacy files (code_saturne version 3.3 and below),
 * the old id will actually be the old scalar id (-1 for others).
 *
 * parameters:
 *   r             <-> associated restart file pointer
 *   old_field_map <-- name to id map of fields in restart file
 *   key           <-- key for field association
 *   read_flag     <-- optional flag to track fields read, or NULL;
 *                     set to sum of 2^time_id for fields read, -1 for fields
 *                     failed to read (size: n_fields)
 *----------------------------------------------------------------------------*/

void
cs_restart_read_linked_fields(cs_restart_t               *r,
                              const cs_map_name_to_id_t  *old_field_map,
                              const char                 *key,
                              int                         read_flag[]);

/*----------------------------------------------------------------------------
 * Write fields depending on others to checkpoint.
 *
 * Write field metadata to main checkpoint.
 *
 * parameters:
 *   r          <-> associated restart file pointer
 *   key        <-- key for field association
 *   write_flag <-- optional flag to track fields written, or NULL;
 *                  set to sum of 2^time_id for fields written (size: n_fields)
 *
 * returns:
 *   number of fields written
 *----------------------------------------------------------------------------*/

int
cs_restart_write_linked_fields(cs_restart_t  *r,
                               const char    *key,
                               int            write_flag[]);

/*----------------------------------------------------------------------------
 * Read boundary condition coefficients for all fields from checkpoint.
 *
 * parameters:
 *   r <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_read_bc_coeffs(cs_restart_t  *r);

/*----------------------------------------------------------------------------
 * Write boundary condition coefficients for all fields to checkpoint.
 *
 * parameters:
 *   r <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_write_bc_coeffs(cs_restart_t  *r);

/*----------------------------------------------------------------------------
 * Read field values from checkpoint.
 *
 * If the values are not found using the default rules based on the
 * field's name, its name itself, or a "restart_rename" keyed string value,
 * an old name may be used for compatibility with older files.
 * For cell-based fields, the old name base is appended automatically with
 * "_ce_phase01", except for scalars, where the name uses a different scheme,
 * based on "scalaire_ce_%04" % s_num;
 *
 * parameters:
 *   r    <-> associated restart file pointer
 *   f_id <-- field id
 *   t_id <-- time id (0 for current, 1 for previous, ...)
 *
 * returns:
 *   CS_RESTART_SUCCESS in case of success, CS_RESTART_ERR_... otherwise
 *----------------------------------------------------------------------------*/

int
cs_restart_read_field_vals(cs_restart_t  *r,
                           int            f_id,
                           int            t_id);

/*----------------------------------------------------------------------------
 *  Write field values to checkpoint.
 *
 * parameters:
 *   r     <-> associated restart file pointer
 *   f_id  <-- field id
 *   t_id  <-- time id (0 for current, 1 for previous, ...)
 *----------------------------------------------------------------------------*/

void
cs_restart_write_field_vals(cs_restart_t  *r,
                            int            f_id,
                            int            t_id);

/*----------------------------------------------------------------------------
 * Read restart time step info.
 *
 * parameters:
 *   r  <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_read_time_step_info(cs_restart_t  *r);

/*----------------------------------------------------------------------------
 * Loop over all fields and save them in the restart file which id is
 * passed in argument if it matches their "restart_file" key value.
 *
 * parameters:
 *   r      <-> associated restart file pointer
 *   r_id   <-- value of the key "restart_file"
 *----------------------------------------------------------------------------*/

void
cs_restart_write_fields(cs_restart_t        *r,
                        cs_restart_file_t    r_id);

/*----------------------------------------------------------------------------
 * Loop over all fields and read them in the restart file which id is
 * passed in argument if it matches their "restart_file" key value.
 *
 * parameters:
 *   r      <-> associated restart file pointer
 *   r_id   <-- value of the key "restart_file"
 *----------------------------------------------------------------------------*/

void
cs_restart_read_fields(cs_restart_t       *r,
                       cs_restart_file_t   r_id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set restart file values for fields when those values cannot
 *        be determined at field definition time.
 *
 * This is needed when the need for restart data depends on various
 * combinations of settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_set_auxiliary_field_options(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Initialize fields read status array
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_initialize_fields_read_status(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Finalize fields read status array
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_finalize_fields_read_status(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Get checkpoint read status for a field based on its id
 *
 * \param[in] f_id  field id
 *
 * \returns 0 if field read action failed, 1 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_get_field_read_status(const int f_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RESTART_DEFAULT_H__ */
