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
#include "cs_file.h"

#include "cs_file_csv_parser.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

// Define max char array size based on C99 specifications
#define CS_MAX_STR_SIZE 65535 / sizeof(char)

/*============================================================================
 * Private functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a token based on the difference between two strings (beginning)
 *
 * \param[in] str1  Longest string which should contain the token at its
 *                  beginning
 * \param[in] str2  Shortest string. If nullptr treated as 0 length string.
 *
 * \returns token
 */
/*----------------------------------------------------------------------------*/

static char *
_get_token(const char *str1,
           const char *str2)
{
  char *t = nullptr;

  size_t l_1 = strlen(str1);
  size_t l_2 = (str2 != nullptr) ? strlen(str2) : 0;

  if (l_1 > l_2) {
    size_t l_t = l_1 - l_2;
    BFT_MALLOC(t, l_t + 1, char);
    memcpy(t, str1, l_t);
    t[l_t] = '\0';
  }

  return t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Parse a line using a given separator and returns the tokens.
 *
 * \param[in]      line                  line to parse (char *)
 * \param[in]      separator             separator (int)
 * \param[in, out] n_tokens              number of tokens
 * \param[in]      keep_missing_tokens   true or false
 *
 * \return array of strings, each string corresponding to a token
 */
/*----------------------------------------------------------------------------*/

static char **
_parse_line(char          *line,
            const char    *separator,
            int           *n_tokens,
            bool           keep_missing_tokens)
{
  char **tokens = nullptr;

  int _nt = 0;

  char *n = nullptr;
  char *c = line;

  while (true) {

    int isep = (int)separator[0];
    n = strchr(c, isep);

    if (n != nullptr) {

      char *_t = _get_token(c, n);
      if (keep_missing_tokens || (_t != nullptr && strlen(_t) > 0) ) {
        BFT_REALLOC(tokens, _nt + 1, char *);
        tokens[_nt] = _t;
        _nt += 1;
      }
      else
        BFT_FREE(_t);
      c = n + 1;

    }
    else if (c != nullptr && strlen(c) > 0 && strcmp(c, "\n")) {

      /* Check if there is still a token left (line not ending with the
         separator) */

      char *_t = _get_token(c, n);
      if (keep_missing_tokens || (_t != nullptr && strlen(_t) > 0) ) {
        BFT_REALLOC(tokens, _nt + 1, char *);
        BFT_FREE(_t);
        tokens[_nt] = _get_token(c, n);
        _nt += 1;
      }
      else
        BFT_FREE(_t);
      break;

    }
    else
      break;

  }

  *n_tokens = _nt;

  return tokens;
}

#if 0

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the i-th token of a given line.
 *
 * \param[in] line       Line to parse
 * \param[in] separator  Separator
 * \param[in] position   Position (int) of the searched token
 *
 * \return token
 */
/*----------------------------------------------------------------------------*/

static char *
_get_token_from_line(char    *line,
                     int      separator,
                     int      position)
{
  int _nt = 0;

  char *n = nullptr;
  char *c = line;

  char *retval = nullptr;

  while (true) {
    n = strchr(c, separator);

    if (n != nullptr) {
      if (_nt == position) {
        size_t l_c = strlen(c);
        size_t l_n = strlen(n);

        if (l_c > l_n) {
          size_t l_t = l_c - l_n;
          BFT_MALLOC(retval, l_t + 1, char);
          memcpy(retval, c, l_t);
          retval[l_t] = '\0';
        }
        break;
      }
      _nt += 1;
    }
    else
      break;
  }

  return retval;
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a line is empty
 *
 * \param[in] line  text line read from file to check
 *
 * \return 1 if line is empty, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

static int
_check_if_line_is_empty(char *line)
{
  int retval = 0;

  if (line == nullptr)
    retval = 1;
  else if (strcmp(line, "\n") == 0 ||
           strcmp(line, "\r\n") == 0 ||
           strcmp(line, "\t") == 0)
    retval = 1;
  else if (strlen(line) <= 1)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count the number of lines in the file
 *
 * \param[in] file_name Name of the file to open
 *
 * \returns number of lines
 */
/*----------------------------------------------------------------------------*/

static int
_count_lines_in_file(const char *file_name)
{
  FILE *f = fopen(file_name, "r");

  char line[CS_MAX_STR_SIZE] = "";

  int n_lines = 0;

  while (fgets(line, CS_MAX_STR_SIZE, f)) {
    if (_check_if_line_is_empty(line) == 0)
      n_lines += 1;
  }

  fclose(f);

  return n_lines;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
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
 * \param[in] ignore_missing_tokens  Ignore missing tokens (nullptr)
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
                  int         *n_cols)
{
  if (cs_file_isreg(file_name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              "Error: file \"%s\" does not exist.\n", file_name);

  if (n_columns > -1 && col_idx == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Error: a subset of the columns is requested without providing "
              "the list of columns");

  if (separator == nullptr || strcmp(separator, "") == 0)
    bft_error(__FILE__, __LINE__, 0,
              "Error: empty separator provided.\n");

  bool keep_missing_tokens = !(ignore_missing_tokens);
  int _n_pts_max = _count_lines_in_file(file_name);

  const int _n_rows = _n_pts_max - n_headers;

  *n_rows = _n_rows;

  char ***retval = nullptr;
  BFT_MALLOC(retval, _n_rows, char **);
  for (int i = 0; i < _n_rows; i++)
    retval[i] = nullptr;

  FILE *f = fopen(file_name, "r");
  char line[CS_MAX_STR_SIZE] = "";

  int row = 0;

  int read_lines = 0;
  while (fgets(line, CS_MAX_STR_SIZE, f)) {

    if (read_lines >= n_headers) {
      if (_check_if_line_is_empty(line)) {
        bft_printf(_("Warning: line #%d of file \"%s\" is empty and "
                     "will be ignored. Please consider removing it from the "
                     "file.\n"),
                   read_lines + 1, file_name);
      }
      else {
        char **tokens = nullptr;
        int    n_tokens = 0;

        tokens = _parse_line(line, separator, &n_tokens, keep_missing_tokens);

        assert(n_columns <= n_tokens);

        int _n_col = (n_columns == -1) ? n_tokens : n_columns;

        if (retval[0] == nullptr) {
          *n_cols = _n_col;
          for (int i = 0; i < _n_rows; i++) {
            BFT_MALLOC(retval[i], _n_col, char *);
            for (int j = 0; j < _n_col; j++)
              BFT_MALLOC(retval[i][j], 120, char);
          }
        }

        for (int col = 0; col < _n_col; col++) {
          int tk_id = (n_columns == -1) ? col : col_idx[col];
          strncpy(retval[row][col], tokens[tk_id], 120);
        }

        row += 1;

        // Free tokens
        for (int i = 0; i < n_tokens; i++)
          BFT_FREE(tokens[i]);
        BFT_FREE(tokens);
      }
    }
    read_lines += 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
