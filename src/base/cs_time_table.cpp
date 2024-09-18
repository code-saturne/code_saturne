/*============================================================================
 * Tabulation handling for code_saturne
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_file.h"
#include "cs_file_csv_parser.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_time_table.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Time table structure */

struct _cs_time_table_t {
 char *name;                /* Table name */

 char      **headers;       /* Column headers */
 cs_real_t **columns;       /* Data columns */

 cs_real_t time_offset;     /* Offset for computation */
 int time_col_id;           /* Index of time column */

 int n_cols;                /* Number of columns */
 int n_rows;                /* Number of rows */

 cs_double_int_t coeffs[2]; /* Coefficients used for interpolation */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of time tables defined */
static int _n_time_tables = 0;

/* Array of time tables */
static cs_time_table_t **_time_tables = nullptr;

/*============================================================================
 * Private functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information concerning a time table, for debugging
 *
 * \param[in] t pointer to time table
 */
/*----------------------------------------------------------------------------*/

#if 0

static void
_print_time_table(cs_time_table_t *t)
{
  bft_printf("------------------------------\n");
  bft_printf(" PRINTING TIME TABLE INFORMATION \n");
  bft_printf("------------------------------\n");

  bft_printf("[NAME]  \"%s\"\n", t->name);
  bft_printf("[NROWS] \"%d\"\n", t->n_rows);
  bft_printf("[NCOLS] \"%d\"\n", t->n_cols);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get time table by name (internal function)
 *
 * \param[in] name  name of the time table
 *
 * \return pointer to corresponding time table. nullptr if not found
 */
/*----------------------------------------------------------------------------*/

static cs_time_table_t *
_time_table_by_name_try(const char *name)
{
  if (name == nullptr || strcmp(name, "") == 0)
    bft_error(__FILE__, __LINE__, 0,
              "Error: Empty time table name.\n");

  cs_time_table_t *retval = nullptr;

  for (int i = 0; i < _n_time_tables; i++) {
    if (strcmp(_time_tables[i]->name, name) == 0) {
      retval = _time_tables[i];
      break;
    }
  }

  return retval;
}

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

static int
_time_table_column_id_by_name(const cs_time_table_t *table,
                              const char            *name)
{
  assert(table != nullptr);
  assert(name != nullptr);
  assert(table->headers != nullptr);

  if (table->headers == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: table \"%s\" has no defined headers.\n"), table->name);

  int retval = -1;
  for (int i = 0; i < table->n_cols; i++) {
    if (strcmp(name, table->headers[i]) == 0) {
      retval = i;
      break;
    }
  }

  return retval;
}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate a new time table pointer and returns it
 *
 * \param[in] name  name of the new time table
 *
 * \return pointer to newly created time table
 */
/*----------------------------------------------------------------------------*/

static cs_time_table_t *
_time_table_create(const char *name)
{

  cs_time_table_t *retval = _time_table_by_name_try(name);

  if (retval != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Error: time table \"%s\" allready exists.\n",
              name);

  int new_id = _n_time_tables;

  BFT_REALLOC(_time_tables, _n_time_tables+1, cs_time_table_t *);

  BFT_MALLOC(retval, 1, cs_time_table_t);

  BFT_MALLOC(retval->name, strlen(name)+1, char);
  strcpy(retval->name, name);

  retval->n_rows      = 0;
  retval->n_cols      = 0;
  retval->time_col_id = 0;
  retval->time_offset = 0.;

  for (int i = 0; i < 2; i++) {
    retval->coeffs[i].id  = 0;
    retval->coeffs[i].val = 0.;
  }

  retval->headers = nullptr;
  retval->columns = nullptr;

  _time_tables[new_id] = retval;
  _n_time_tables += 1;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a given time table
 *
 * \param[in] t pointer to time table which is to be freed.
 */
/*----------------------------------------------------------------------------*/

static void
_free_time_table(cs_time_table_t *t)
{
  assert(t != nullptr);

  for (int i = 0; i < t->n_cols; i++) {
    BFT_FREE(t->columns[i]);
    if (t->headers != nullptr)
      BFT_FREE(t->headers[i]);
  }
  BFT_FREE(t->columns);
  BFT_FREE(t->headers);

  BFT_FREE(t->name);

  BFT_FREE(t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute value for column using the current defined time.
 *
 * \param[in] table   Pointer to time table structure
 * \param[in] col_id  Column id for which to compute time value
 *
 * \return time value (cs_real_t)
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_time_table_compute_value(cs_time_table_t *table,
                          const int        col_id)
{
  assert(table != nullptr);
  cs_real_t *_vals = table->columns[col_id];
  cs_double_int_t *c = table->coeffs;

  cs_real_t retval = _vals[c[0].id] * c[0].val + _vals[c[1].id] * c[1].val;

  return retval;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Try to get time table based on name.
 *
 * \param[in] name Name of time table
 *
 * \returns pointer to time table, nullptr if not found.
 */
/*----------------------------------------------------------------------------*/

cs_time_table_t *
cs_time_table_by_name_try(const char *name)
{
  cs_time_table_t *retval = _time_table_by_name_try(name);

  return retval;
}

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
cs_time_table_by_name(const char *name)
{
  cs_time_table_t *retval = _time_table_by_name_try(name);

  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Error: time table \"%s\" does not exist.\n",
              name);

  return retval;
}

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
                         cs_real_t        time_offset)
{
  assert(table != nullptr);

  table->time_offset = time_offset;
}

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
                          const char     **headers)
{
  assert(n_headers > 0 && headers != nullptr);

  if (n_headers != table->n_cols)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: table \"%s\" has %d columns but you defined only %d"
                " headers"),
              table->name, table->n_cols, n_headers);

  BFT_MALLOC(table->headers, n_headers, char *);
  for (int i = 0; i < n_headers; i++) {
    char *_h = nullptr;
    BFT_MALLOC(_h, strlen(headers[i]) + 1, char);

    strcpy(_h, headers[i]);
    _h[strlen(headers[i])] = '\0';

    table->headers[i] = _h;
  }
}

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
                            const bool   ignore_missing_tokens)
{
  assert(name != nullptr);
  assert(file_name != nullptr);
  assert(separator != nullptr);

  cs_time_table_t *t = cs_time_table_by_name_try(name);
  if (t != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: time table \"%s\" allready exists.\n"),
              name);

  int _n_rows = 0;
  int _n_cols = 0;

  char ***_data = cs_file_csv_parse(file_name,
                                    separator,
                                    n_headers,
                                    n_columns,
                                    col_idx,
                                    ignore_missing_tokens,
                                    &_n_rows,
                                    &_n_cols);

  t = _time_table_create(name);

  t->n_rows = _n_rows;
  t->n_cols = _n_cols;

  BFT_MALLOC(t->columns, _n_cols, cs_real_t *);
  for (int i = 0; i < _n_cols; i++)
    BFT_MALLOC(t->columns[i], _n_rows, cs_real_t);

  for (int ir = 0; ir < _n_rows; ir++) {
    char **_row = _data[ir];
    for (int ic = 0; ic < _n_cols; ic++)
      t->columns[ic][ir] = atof(_row[ic]);
  }

  // Free data which is no longer needed.
  for (int i = 0; i < _n_rows; i++) {
    for (int j = 0; j < _n_cols; j++)
      BFT_FREE(_data[i][j]);
    BFT_FREE(_data[i]);
  }
  BFT_FREE(_data);

  return t;
}

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
                                   const char *separator)
{
  cs_time_table_t *retval = cs_time_table_from_csv_file(name,
                                                        file_name,
                                                        separator,
                                                        0,
                                                        -1,
                                                        nullptr,
                                                        true);

  return retval;
}

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
                                           const int   n_headers)
{
  cs_time_table_t *retval = cs_time_table_from_csv_file(name,
                                                        file_name,
                                                        separator,
                                                        n_headers,
                                                        -1,
                                                        nullptr,
                                                        true);

  return retval;
}

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
                              const int        col_id)
{
  assert(table != nullptr);
  assert(col_id > -1 && col_id < table->n_cols);
  table->time_col_id = col_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the column id for time based on a label
 *
 * \param[in] table       Pointer to time table structure
 * \param[in] time_label  Label to identify index of column used as time
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_set_time_from_label(cs_time_table_t *table,
                                  const char      *time_label)
{
  assert(table != nullptr);

  const int t_id = _time_table_column_id_by_name(table, time_label);

  if (t_id > -1)
    table->time_col_id = t_id;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Error: table \"%s\" has no column with header \"%s\"\n"),
              table->name, time_label);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset time table time value (force search from beginning of table).
 *
 * \param[in] table  Pointer to time table structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_reset_position(cs_time_table_t *table)
{
  assert(table != nullptr);

  for (int i = 0; i < 2; i++) {
    table->coeffs[i].id  = 0;
    table->coeffs[i].val = 0.;
  }
}

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
                              bool             reset_time_value)
{
  // Reset time
  if (reset_time_value)
    cs_time_table_reset_position(table);

  const cs_real_t *time_vals = table->columns[table->time_col_id];
  const int n_rows = table->n_rows;
  cs_double_int_t *coeffs = table->coeffs;

  const int t0_id = coeffs[0].id;

  /* Compute time for interpolation using the table defined offset value */
  cs_real_t _time = time + table->time_offset;

  if (_time < time_vals[0]) {
    coeffs[0].id  = 0;
    coeffs[1].id  = 0;
    coeffs[0].val = 1.;
    coeffs[1].val = 0.;
  }
  else if (_time > time_vals[n_rows - 1]) {
    coeffs[0].id = n_rows - 1;
    coeffs[1].id = n_rows - 1;
    coeffs[0].val = 1.;
    coeffs[1].val = 0.;
  }
  else {
    for (int i = t0_id; i < n_rows - 1; i++) {
      if (_time >= time_vals[i] && _time < time_vals[i+1]) {
        coeffs[1].id = i + 1;
        coeffs[1].val = (_time - time_vals[i]) / (time_vals[i+1] - time_vals[i]);

        coeffs[0].id = i;
        coeffs[0].val = 1. - coeffs[1].val;

        break;
      }
    }
  }
}

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
                                 bool        overwrite_prev)
{
  cs_time_table_t *table = cs_time_table_by_name(name);

  cs_time_table_update_position(table, t, overwrite_prev);

  return _time_table_compute_value(table, col);
}

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
                                          bool        overwrite_prev)
{
  assert(name != nullptr);
  assert(label != nullptr);

  cs_time_table_t *table = cs_time_table_by_name(name);

  int _id = _time_table_column_id_by_name(table, label);

  cs_time_table_update_position(table, t, overwrite_prev);

  return _time_table_compute_value(table, _id);
}

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
                                    cs_real_t  *retvals)
{
  cs_time_table_t *table = cs_time_table_by_name(name);

  cs_time_table_update_position(table, t, overwrite_prev);

  for (int i = 0; i < n_cols; i++)
    retvals[i] = _time_table_compute_value(table, cols[i]);
}

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
                                             cs_real_t  *retvals)
{
  cs_time_table_t *table = cs_time_table_by_name(name);
  cs_time_table_update_position(table, t, overwrite_prev);

  for (int i = 0; i < n_cols; i++) {
    int _id = _time_table_column_id_by_name(table, labels[i]);
    retvals[i] = _time_table_compute_value(table, _id);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all data structures related to datasets
 */
/*----------------------------------------------------------------------------*/

void
cs_time_table_destroy_all(void)
{
  for (int i = 0; i < _n_time_tables; i++)
    _free_time_table(_time_tables[i]);

  BFT_FREE(_time_tables);
  _n_time_tables = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
