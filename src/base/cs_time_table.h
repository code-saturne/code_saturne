#ifndef __CS_TIME_TABLE_H__
#define __CS_TIME_TABLE_H__

/*============================================================================
 * Tabulation handling for code_saturne
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_time_table_t cs_time_table_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get column id corresponding to a given header/label.
 *
 * \param[in] table pointer to time table
 * \param[in] name  name of the header/column
 *
 * \return index of the corresponding column. -1 if not found.
 */
/*----------------------------------------------------------------------------*/

int
cs_time_table_column_id_by_name(const cs_time_table_t *table,
                                const char            *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Try to get time table based on name.
 *
 * \param[in] name Name of time table
 *
 * \returns pointer to time table, NULL if not found.
 */
/*----------------------------------------------------------------------------*/

cs_time_table_t *
cs_time_table_by_name_try(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get time table based on name.
 *
 * \param[in] name Name of time table
 *
 * \returns pointer to time table, raises error if not found.
 */
/*----------------------------------------------------------------------------*/

cs_time_table_t *
cs_time_table_by_name(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set time offset value for a time table
 *
 * \param[in] table        Pointer to time table structure
 * \param[in] time_offset  Value of time offset for the table.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_set_offset(cs_time_table_t *table,
                         cs_real_t        time_offset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set headers labels for a time table. Number of headers needs to be
 * equal to the number of columns, otherwise an error is raised.
 *
 * \param[in] table      Pointer to time table structure
 * \param[in] n_headers  Number of headers to define
 * \param[in] headers    Array of headers strings
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_set_headers(cs_time_table_t *table,
                          const int        n_headers,
                          const char     **headers);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time table from a CSV file.
 *
 * \param[in] name                  Name of the table to be created
 * \param[in] file_name             Path to CSV file
 * \param[in] separator             Separator used in the CSV file
 * \param[in] n_headers             Number of header lines to be ignored during parsing
 * \param[in] n_columns             Number of columns to read. -1 for all
 * \param[in] col_idx               Array containing indexes of columns to read if n_columns != -1
 * \param[in] ignore_missing_tokens Should we ignore missing tokens
 *
 * \returns pointer to newly created time table
 */
/*----------------------------------------------------------------------------*/

cs_time_table_t *
cs_time_table_from_csv_file(const char  *name,
                            const char  *file_name,
                            const char  *separator,
                            const int    n_headers,
                            const int    n_columns,
                            const int   *col_idx,
                            const bool   ignore_missing_tokens);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time table from a CSV file. We suppose that all columns are
 * to be read, no headers line to skip and that missing tokens are to be ignored.
 *
 * \param[in] name                  Name of the table to be created
 * \param[in] file_name             Path to CSV file
 * \param[in] separator             Separator used in the CSV file
 *
 * \returns pointer to newly created time table
 */
/*----------------------------------------------------------------------------*/

cs_time_table_t *
cs_time_table_from_csv_file_simple(const char *name,
                                   const char *file_name,
                                   const char *separator);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time table from a CSV file. We suppose that all columns are
 * to be read and that missing tokens are to be ignored.
 *
 * \param[in] name                  Name of the table to be created
 * \param[in] file_name             Path to CSV file
 * \param[in] separator             Separator used in the CSV file
 * \param[in] n_headers             Number of header lines to be ignored during parsing
 *
 * \returns pointer to newly created time table
 */
/*----------------------------------------------------------------------------*/

cs_time_table_t *
cs_time_table_from_csv_file_simple_headers(const char *name,
                                           const char *file_name,
                                           const char *separator,
                                           const int   n_headers);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the column id for time based on an index
 *
 * \param[in] table   Pointer to time table structure
 * \param[in] col_id  Index of column which is to be used as time
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_set_time_col_id(cs_time_table_t *table,
                              const int        col_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the column id for time based on a label
 *
 * \param[in] table      Pointer to time table structure
 * \param[in] time_label Label to identify index of column which is to be used as time
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_set_time_from_label(cs_time_table_t *table,
                                  const char      *time_label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset time table time value (force search from beginning of table).
 *
 * \param[in] table  Pointer to time table structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_reset_position(cs_time_table_t *table);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update time coefficients used for data interpolation.
 *
 * \param[in] table            Pointer to time table structure
 * \param[in] time             Time value
 * \param[in] reset_time_value Reset current time value (bool)
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_update_position(cs_time_table_t *table,
                              cs_real_t        time,
                              bool             reset_time_value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute value using a given abscissa and a specific column
 *
 * \param[in] name            Name of the used time table
 * \param[in] t               Time for which we seek values
 * \param[in] col             Index of column used for computation
 * \param[in] overwrite_prev  Start search of value using first value
 *
 * \returns Interpolated value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_time_table_compute_time_value(const char *name,
                                 cs_real_t   t,
                                 const int   col,
                                 bool        overwrite_prev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute value using a given abscissa and a label
 *
 * \param[in] name            Name of the used time table
 * \param[in] t               Time for which we seek values
 * \param[in] label           Label of column used for computation
 * \param[in] overwrite_prev  Start search of value using first value
 *
 * \returns Interpolated value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_time_table_compute_time_value_by_label(const char *name,
                                          cs_real_t   t,
                                          const char *label,
                                          bool        overwrite_prev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute values for several columns of a time table for a given
 *        abscissa
 *
 * \param[in]  name            Name of the time table to use
 * \param[in]  t               Time for which we seek values
 * \param[in]  n_cols          Number of values to compute
 * \param[in]  cols            Array with the indices of columns for computation
 * \param[in]  overwrite_prev  Start search of value using first value
 * \param[out] retvals         Array of output values
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_compute_n_time_values(const char *name,
                                    cs_real_t   t,
                                    const int   n_cols,
                                    const int   cols[],
                                    bool        overwrite_prev,
                                    cs_real_t  *retvals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute values for several columns of a time table for a given
 *        abscissa
 *
 * \param[in]  name            Name of the time table to use
 * \param[in]  t               Time for which we seek values
 * \param[in]  n_cols          Number of values to compute
 * \param[in]  labels          Array with the labels of columns for computation
 * \param[in]  overwrite_prev  Start search of value using first value
 * \param[out] retvals         Array of output values
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_compute_n_time_values_by_label(const char *name,
                                             cs_real_t   t,
                                             const int   n_cols,
                                             const char *labels[],
                                             bool        overwrite_prev,
                                             cs_real_t  *retvals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute value from time table based on label and current time.
 * Positions are not updated for the table.
 *
 * \param[in] name            Name of the used time table
 * \param[in] label           Label of column used for computation
 *
 * \returns Interpolated value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
CS_TIME_TABLE(const char *name,
              const char *label)
{
  return cs_time_table_compute_time_value_by_label(name,
                                                   cs_glob_time_step->t_cur,
                                                   label,
                                                   false);
}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all data structures related to datasets
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_destroy_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIME_TABLE_H__ */
