#ifndef __CS_FILE_CSV_PARSER_H__
#define __CS_FILE_CSV_PARSER_H__

/*============================================================================
 * Read data from CSV files
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Parse a csv file and export to a dataset (char ***).
 *
 * The caller is responsible for freeing the dataset when not needed anymore.
 *
 * \param[in] file_name              Name of the file to read
 * \param[in] separator              Separator (int)
 * \param[in] n_headers              Number of headers (to ignore during import)
 * \param[in] n_columns              Number of columns to read.
 *                                   -1 if all columns are to be read
 * \param[in] col_idx                Array of indices of columns to read
 *                                   (if n_columns != -1)
 * \param[in] ignore_missing_tokens  Ignore missing tokens (NULL)
 * \param[in] n_rows                 Pointer to number of rows in file
 * \param[in] n_cols                 Pointer to number of columns in file
 *
 * \returns Pointer to newly created dataset.
 */
/*----------------------------------------------------------------------------*/

char ***
cs_file_csv_parse(const char  *file_name,
                  const char  *separator,
                  const int    n_headers,
                  const int    n_columns,
                  const int   *col_idx,
                  const bool   ignore_missing_tokens,
                  int         *n_rows,
                  int         *n_cols);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FILE_CSV_PARSER_H__ */
