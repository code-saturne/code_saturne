/*============================================================================
 * Program logging information
 *============================================================================*/

/*
  This file is part of the Code_Saturne Kernel, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1998-2011 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Kernel is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Kernel is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Kernel; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*-----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_log.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static bool  _cs_log_atexit_set = false;

static FILE* _cs_log[] = {NULL};
static const char* _cs_log_name[] = {"performance.log"};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Close all log files.
 *----------------------------------------------------------------------------*/

static void
_close_log_files(void)
{
  int i;

  for (i = 0; i < CS_LOG_N_TYPES; i++) {
    if (_cs_log[i] != NULL)
      fclose(_cs_log[i]);
  }
}

/*----------------------------------------------------------------------------
 * Open log file.
 *
 * parameters:
 *   log <-- log file type
 *----------------------------------------------------------------------------*/

static void
_open_log(cs_log_t log)
{
  if (cs_glob_rank_id < 1 && _cs_log[log] == NULL) {
    _cs_log[log] = fopen(_cs_log_name[log], "w");
    if (_cs_log[log] == NULL)
      bft_error(__FILE__, __LINE__, errno,
                _("Error opening log file: %s"),
                _cs_log_name[log]);
    if (_cs_log_atexit_set == false) {
      atexit(_close_log_files);
      _cs_log_atexit_set = true;
    }
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count printable length of a character string.
 *
 * This should also include UTF-8 strings.
 *
 * \param[in]  str  pointer to printable string
 *
 * \return  printable length of character string
 *----------------------------------------------------------------------------*/

size_t
cs_log_strlen(const char  *str)
{
  static int mode_utf8 = -1;

  int l = 0;
  int retval = 0;

  if (mode_utf8 == -1) {
    char *lang = getenv("LANG");
    mode_utf8 = 0;
    if (lang != NULL) {
      if (   strcmp(lang + strlen(lang) - 5, "UTF-8") == 0
          || strcmp(lang + strlen(lang) - 4, "utf8") == 0)
        mode_utf8 = 1;
    }
  }

  if (str != NULL) {

    l = strlen(str);

    if (mode_utf8 == 0)
      retval = l;

    else if (mode_utf8 == 1) {

      int i;
      bool multibyte = false;

      for (i = 0; i < l; i++) {

        char c = str[i];

        if (multibyte == false || (c < 0x80 || c > 0xBF)) {

          multibyte = false;

          if (c <= 0x7F) {
            retval++;
          }
          else if (c >= 0xC0 || c <= 0xFD) {
            multibyte = true;
            retval++;
          }
        }
      }
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print log info to a given log type.
 *
 * The format and variable arguments are similar to those of the printf()
 * type functions.
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  log     log file type
 * \param[in]  format  format string, as printf() and family.
 * \param[in]  ...     variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\0' used
 *         to end output strings
 */
/*----------------------------------------------------------------------------*/

int
cs_log_printf(cs_log_t     log,
              const char  *format,
              ...)
{
  int  retval;
  va_list  arg_ptr;

  if (cs_glob_rank_id > 0)
    return 0;
  else if (_cs_log[log] == NULL)
    _open_log(log);

  va_start(arg_ptr, format);

  retval = vfprintf(_cs_log[log], format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print a string to a given width info to a given type.
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  log    log file type
 * \param[in]  str    string to print.
 * \param[in]  width  minimum width to which the string must be printed.
 *
 * \return number of characters printed, not counting the trailing '\0' used
 *         to end output strings
 */
/*----------------------------------------------------------------------------*/

int
cs_log_print_padded_str(cs_log_t     log,
                        const char  *str,
                        int          width)
{
  int l = cs_log_strlen(str);
  int  retval = 0;

  if (cs_glob_rank_id > 0)
    return 0;
  else if (_cs_log[log] == NULL)
    _open_log(log);

  if (str != NULL)
    retval += fprintf(_cs_log[log], "%s", str);

  if (width > l)
    fprintf(_cs_log[log], "%-*s", width-l, "");

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flush output of a log file.
 *
 * In parallel, output is only handled by rank 0.
 *
 * If the argument is set to CS_LOG_N_TYPES, all log files are flushed.
 *
 * \param[in]  log  log file type
 *
 * \return 0 upon successful completion 0 is returned. Otherwise,
 *           EOF is returned and  errno  is  set  to indicate the error.
 */
/*----------------------------------------------------------------------------*/

int
cs_log_printf_flush(cs_log_t log)
{
  int i;
  int retval = 0;

  if (log < CS_LOG_N_TYPES)
    retval = fflush(_cs_log[log]);

  else {
    for (i = 0; i < CS_LOG_N_TYPES; i++) {
      retval = fflush(_cs_log[log]);
      if (retval != 0)
        break;
    }
  }

  return retval;
}


/*-----------------------------------------------------------------------------*/

END_C_DECLS
