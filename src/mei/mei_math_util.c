/*!
 * \file mei_math_util.c
 *
 * \brief Provides mathemathical functions facilities
 */

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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <assert.h>

#if defined(HAVE_UNISTD_H)
# include <unistd.h>
#endif

/*----------------------------------------------------------------------------
 * Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mei_math_util.h"

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local macro definitions
 *-----------------------------------------------------------------------------*/

#undef SIZE_MAX
#define SIZE_MAX 1000

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static mei_user_data_t **data = NULL;  /* array of pointers on structures
                                          that contains user data sets */

static int data_length = 0;            /* number of user data sets */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump the content of a mei_user_data structure.
 *
 * \param [in] d structure that contains the user data for 1D interpolation
 */
/*----------------------------------------------------------------------------*/

#if 0
static void
_user_data_dump(const mei_user_data_t *d)
{
  bft_printf("\n\nDUMP OF THE MEI_USER_DATA STRUCTURE: %p\n\n",
             (const void *)d);

  bft_printf("  name:              %s\n",   d->name);
  bft_printf("  description:\n%s\n",        d->description);
  bft_printf("  number of columns: %i\n",   d->ncols);
  bft_printf("  number of lines:   %i\n\n", d->nrows);

  for (int i = 0; i < d->nrows; i++) {
    bft_printf("\nline #%i: ", i);
    for (int j = 0; j < d->ncols; j++) {
      bft_printf("%f ", d->values[i][j]);
    }
  }

  bft_printf("\n\nEND OF DUMP OF MEI_USER_DATA STRUCTURE\n\n");
  bft_printf_flush();
  return;
}
#endif

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Compare two strings.
 *
 * \param [in] s1 first string
 * \param [in] s2 second string
 *
 * \return  1 if the strings are equal, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

static int
_user_data_strcmp(const char *s1,
                  const char *s2)
{
  if (s1 == NULL || s2 == NULL) return 0;
  if ( strlen(s1) != strlen(s2)) return 0;
  if (!strncmp(s1, s2, strlen(s1))) return 1;
  return 0;
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Read a single user data set.
 *
 * \param [in] filename name of the file that contents the data set
 */
/*-----------------------------------------------------------------------------*/

static void
_user_data_reader(const char *filename)
{
  FILE *fr;
  int nb_col = 0;
  int nb_line_tot = 0;
  int nb_line = 0;
  int row;
  char line[SIZE_MAX];
  char *string_tok = NULL;
  char *saveptr = NULL;
  char *buff = NULL;
  char *ext = NULL;
  char *separator = NULL;
  char *end = NULL;

  assert(filename);

  /* Verification of the existence of the file while opening it */

  int file_descriptor = open(filename, O_RDONLY);

  if (file_descriptor ==  -1)
    bft_error(__FILE__, __LINE__, 0,
              _("The user data file: %s needed by interp1d is not found.\n"
                "Add this data file to the list of users file to copy \n"
                "in the computational directory.\n"), filename);
  else
    /* If file exists, close it. It will be reopened. */
    close(file_descriptor);

  /* Searching the file extension */

  BFT_MALLOC(buff, strlen(filename) + 1, char);
  saveptr = buff;  /* keep the pointer for free memory */
  strcpy(buff, filename);

  /* first call of strtok for initialization */

  ext = strtok(buff, ".");

  do {
    ext = strtok(NULL, ".");

    if (_user_data_strcmp(ext, "dat")) {
      BFT_MALLOC(separator, 4, char);
      strcpy(separator, " \t");
      break;

    }
    else if (_user_data_strcmp(ext, "csv")) {
      BFT_MALLOC(separator, 3, char);
      strcpy(separator, ",");
      break;
    }

  } while (ext != NULL);

  if (!(_user_data_strcmp(ext, "dat") || _user_data_strcmp(ext, "csv")))
    bft_error(__FILE__, __LINE__, 0,
              _("interp1d: the file extention expected "
                "is dat or csv.  Extension found: %s\n"), ext);

  BFT_FREE(saveptr);

  /* opening of the file */

  fr = fopen(filename, "rt");

  /* Position at the beginning of the file */

  fseek(fr, 0, SEEK_SET);

  /* 1- searching the number of lines */

  nb_line_tot = 0;
  nb_line = 0;

  while (fgets(line, SIZE_MAX, fr) != NULL) {
    nb_line_tot++;
    if (strncmp(line, "#", 1) && strncmp(line, "\n", 1))
      nb_line++;
  }

  if (nb_line < 2)
    bft_error(__FILE__, __LINE__, 0,
              _("At least two lines are expected in %s.\n"), filename);

  data[data_length-1]->nrows = nb_line;

  BFT_MALLOC(data[data_length-1]->values,
             data[data_length-1]->nrows,
             double*);

  /* 2- searching the number of columns */

  fseek(fr, 0, SEEK_SET);

  nb_col = 0;
  fgets(line, SIZE_MAX ,fr);

  while(!strncmp(line, "\n", 1) || !strncmp(line, "#", 1))
    fgets(line, SIZE_MAX ,fr);

  BFT_MALLOC(buff, SIZE_MAX, char);
  saveptr = buff;

  string_tok = strtok(line, separator);

  while (string_tok != NULL) {
    if (strncmp(string_tok, "\n", 1))
      nb_col += 1;
    string_tok = strtok(NULL, separator);
  }
  BFT_FREE(saveptr);

  if (nb_col < 2)
    bft_error(__FILE__, __LINE__, 0,
              _("At least two columns are expected in %s.\n"), filename);

  data[data_length-1]->ncols = nb_col;

  for (int i = 0; i < data[data_length-1]->nrows; i++)
    BFT_MALLOC(data[data_length-1]->values[i],
               data[data_length-1]->ncols,
               double);

  /* 3- data storing */

  fseek(fr, 0, SEEK_SET);

  row = 0;

  for (int i = 0; i < nb_line_tot; i++) {

    fgets(line, SIZE_MAX, fr);

    if (!_user_data_strcmp(line, "\n")) {

      if (!strncmp(line, "#", 1)) {
        if (data[data_length-1]->description == NULL) {
          BFT_MALLOC(data[data_length-1]->description,
                     strlen(line)+1,
                     char);
          strcpy(data[data_length-1]->description, line);
        }
        else {
          BFT_REALLOC(data[data_length-1]->description,
                      strlen(data[data_length-1]->description)+strlen(line)+1,
                      char);
          strcat(data[data_length-1]->description, line);
        }
      }
      else {

        BFT_MALLOC(buff, SIZE_MAX, char);
        saveptr = buff;
        string_tok = strtok(line, separator);

        if (string_tok != NULL) {
          row += 1;
          data[data_length-1]->values[row-1][0] = strtod(string_tok, &end);

          if (!_user_data_strcmp(end, ""))
            bft_error(__FILE__, __LINE__, 0,
                      _("In the file %s, line %i, there is a wrong value :%s\n"),
                      filename, i+1, string_tok);

          for (int j = 2; j < nb_col + 1; j++) {
            string_tok = strtok(NULL, separator);

            if (string_tok == NULL || _user_data_strcmp(string_tok, "\n")) {
              bft_error(__FILE__, __LINE__, 0,
                        _("The file %s does not have a correct "
                          "structure at line %i.\n"),
                        filename, i);
            }
            else {
              data[data_length - 1]->values[row - 1][j - 1]
                = strtod(string_tok, &end);

              if (   !_user_data_strcmp(end, "")
                  && !_user_data_strcmp(end, "\n"))
                bft_error(__FILE__, __LINE__, 0,
                          _("In the file %s, line %i, "
                            "there is a wrong value :%s\n"),
                          filename, i, string_tok);
            }
          }
        }
        BFT_FREE(saveptr);
      }
    }
  }

  fclose(fr);

  assert(row == nb_line);

#if 0
  _user_data_dump(data[data_length-1]);
#endif

  BFT_FREE(separator);

  return;
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Create a single user data set in the array data.
 *
 * \param [in] filename name of the file that contents the data set
 */
/*-----------------------------------------------------------------------------*/

static void
_user_data_create(const char *filename)
{
  data_length += 1;

  if (data_length == 1)
    BFT_MALLOC(data, data_length, mei_user_data_t*);
  else
    BFT_REALLOC(data, data_length, mei_user_data_t*);

  BFT_MALLOC(data[data_length - 1], 1, mei_user_data_t);

  BFT_MALLOC(data[data_length - 1]->name, strlen(filename) + 1, char);
  strcpy(data[data_length - 1]->name, filename);

  BFT_MALLOC(data[data_length-1]->description, strlen(" ") + 1, char);
  strcpy(data[data_length-1]->description, "");

  data[data_length-1]->ncols = -1;
  data[data_length-1]->nrows = -1;

  _user_data_reader(filename);
  return;
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Return the 1D interpolation if a value.
 *
 * \param [in] d data set
 * \param [in] values data set
 * \param [in] c1 column number of the file for abscisse
 * \param [in] c2 column number of the file for ordinate
 * \param [in] x variable to interpolate
 *
 * \return interpolated value
 */
/*-----------------------------------------------------------------------------*/

static double
_user_data_interp(const mei_user_data_t *d,
                  const int c1,
                  const int c2,
                  const double x)
{
  int k = 0;
  int position = -1;
  double y;
  int row = d->nrows;
  double **values = d->values;

  /* verification of the abscissa colomn */

  for (int i = 0; i < row - 1; i++) {
    if (values[i + 1][c1 - 1] < values[i][c1 - 1])
      bft_error(__FILE__, __LINE__, 0,
                _("Abscissa colomn is not in the rigth order.\n"));
  }

  /* interpolation */

  if (x > values[row - 1][c1 - 1]) {
    /* if the x value is after the table */
    y = values[row - 2][c2 - 1] +
        (x - values[row - 2][c1 - 1]) *
        (values[row - 1][c2 - 1] - values[row - 2][c2 - 1]) /
        (values[row - 1][c1 - 1] - values[row - 2][c1 - 1]);
  }
  else if (x < values[0][c1 - 1]) {
    /* if the x value is before the table */
    y = values[0][c2-1] + (x - values[0][c1 - 1]) *
        (values[1][c2 - 1] - values[0][c2 - 1]) /
        (values[1][c1 - 1] - values[0][c1 - 1]);
  }
  else {

    while (position < 0 && k <= row - 1) {
        if (x > values[k][c1 - 1])
          k = k + 1;
        else
          position = k - 1;
    }

    y = values[position][c2 - 1] +
        (x - values[position][c1 - 1]) *
        (values[position + 1][c2 - 1] - values[position][c2 - 1]) /
        (values[position + 1][c1 - 1] - values[position][c1 - 1]);
  }

  return y;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Return the 1D interpolation if a value.
 *
 * \param [in] filename name of file of data
 * \param [in] c1 column number of the file for abscisse
 * \param [in] c2 column number of the file for ordinate
 * \param [in] x variable to interpolate
 *
 * \return interpolated value
 */
/*-----------------------------------------------------------------------------*/

double
mei_interp1d(const char  *filename,
             int          c1,
             int          c2,
             double       x)
{
  int data_index = -1;

  if (data_length > 0) {
    /* Data exists: search if data are alerady read */
    for (int i = 0; i < data_length; i++) {
      if (_user_data_strcmp(data[i]->name, filename))
        data_index = i;
    }
    /* Data exists but required data are not found */
    if (data_index == -1) {
      _user_data_create(filename);
      data_index = data_length - 1;
    }
  }
  else {
    /* No data yet */
    _user_data_create(filename);
    data_index = 0;
  }

  return _user_data_interp(data[data_index], c1, c2, x);
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all user data set for 1D interpolation.
 */
/*-----------------------------------------------------------------------------*/

void
mei_data_free(void)
{
  for(int i = 0; i < data_length - 1; i++) {
    BFT_FREE(data[i]->name);
    BFT_FREE(data[i]->description);
    for (int j = 0; j < data[i]->nrows; j++)
      BFT_FREE(data[i]->values[i]);
    BFT_FREE(data[i]->values);
    BFT_FREE(data[i]);
  }

  BFT_FREE(data);

  data_length = 0;
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


