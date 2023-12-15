/*============================================================================
 * Program logging information
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_log.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local variable definitions
 *----------------------------------------------------------------------------*/

static bool  _cs_log_atexit_set = false;

static FILE* _cs_log[] = {NULL, NULL, NULL, NULL};
static const char* _cs_log_name[] = {"",
                                     "setup.log",
                                     "performance.log",
                                     "warnings.log"};

static bool  _cs_log_default_active = true;

int cs_glob_log_frequency = 1;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_log_frequency_get_pointer(int  **ntlist);

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

/*----------------------------------------------------------------------------*
 * Pad a string so that its printable length is the required length.
 *
 * This allows pretty-printing with UTF-8 strings, whose actual length may be
 * larger than their printable length in the presence of multibyte characters.
 *
 * If either the printable length of the string is longer than the target
 * width or the actual length is long than the destination buffer's size,
 * it is truncated.
 *
 * parameters:
 *   dest    -->  pointer to destination buffer
 *   str      <-- pointer to printable string
 *   width    <-- desired printed length
 *   destsize <-- destination buffer size
 *   align    <-- 1: left, 0: right
 *----------------------------------------------------------------------------*/

static void
_log_strpad(char        *dest,
            const char  *src,
            size_t       width,
            size_t       destsize,
            int          align)
{
  size_t i, j;
  size_t pad_l = 0, pad_r = 0, p_len = 0, c_len = 0;
  size_t _destsize = destsize - 1;

  static int mode_utf8 = -1;

  assert(dest != NULL && destsize > 0);

  if (mode_utf8 == -1) {
    char *lang = getenv("LANG");
    mode_utf8 = 0;
    if (lang != NULL) {
      if (   strcmp(lang + strlen(lang) - 5, "UTF-8") == 0
          || strcmp(lang + strlen(lang) - 4, "utf8") == 0)
        mode_utf8 = 1;
    }
  }

  if (src != NULL) {
    if (mode_utf8 == 0) {
      p_len = strlen(src);
      if (p_len > _destsize)
        p_len = _destsize;
      c_len = p_len;
    }
    else { /* UTF-8 case */
      for (i = 0;
           i < _destsize && p_len < width;
           i++) {
        unsigned char c = src[i];
        if (c == '\0') {
          c_len = i;
          break;
        }
        else if (c < 0x80 || c > 0xBF) { /* Single byte or first byte in UTF-8 */
          p_len++;
          c_len = i+1;
        }
      }
    }
  }

  if (p_len < width && c_len < _destsize) {
    size_t pad = width - p_len;
    if (c_len + pad > _destsize)
      pad = _destsize - c_len;
    if (align == 0)
      pad_r = pad;
    else
      pad_l = pad;
  }

  j = 0;
  for (i = 0; i < pad_l; i++)
    dest[j++] = ' ';
  for (i = 0; i < c_len; i++)
    dest[j++] = src[i];
  for (i = 0; i < pad_r; i++)
    dest[j++] = ' ';

  dest[j] = '\0';
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointer to log frequency (ntlist in Fortran).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ntlist  --> pointer to ntlist
 *----------------------------------------------------------------------------*/

void
cs_f_log_frequency_get_pointer(int     **ntlist)
{
  *ntlist = &cs_glob_log_frequency;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update "active" or "inactive" flag of default log.
 *
 * This does not prevent output to the log file, but the flag can be queried
 * using \ref cs_log_default_is_active, so in most cases, this status
 * should be checked before logging output
 *
 * \param[in]  activate  true to activate, false to deactivate.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_default_activate(bool  activate)
{
  _cs_log_default_active = activate;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update "active" or "inactive" flag of default log.
 *
 * This does not prevent output to the log file, but the flag can be queried
 * using \
 *
 * \return  true if active, false otherwise.
 */
/*----------------------------------------------------------------------------*/

bool
cs_log_default_is_active(void)
{
  return _cs_log_default_active;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count printable length of a character string.
 *
 * This should also include UTF-8 strings.
 *
 * \param[in]  str  pointer to printable string
 *
 * \return  printable length of character string.
 */
/*----------------------------------------------------------------------------*/

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
          else {
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
 * \brief Pad a string so that its printable length is the required length.
 *
 * This allows pretty-printing with UTF-8 strings, whose actual length may be
 * larger than their printable length in the presence of multibyte characters.
 *
 * If either the printable length of the string is longer than the target
 * width or the actual length is long than the destination buffer's size,
 * it is truncated.
 *
 * \param[out] dest      pointer to destination buffer
 * \param[in]  src       pointer to printable string
 * \param[in]  width     desired printed length
 * \param[in]  destsize  destination buffer size
 */
/*----------------------------------------------------------------------------*/

void
cs_log_strpad(char        *dest,
              const char  *src,
              size_t       width,
              size_t       destsize)
{
  _log_strpad(dest, src, width, destsize, 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Pad a string on the left so that its printable length is
 * the required length.
 *
 * This allows pretty-printing with UTF-8 strings, whose actual length may be
 * larger than their printable length in the presence of multibyte characters.
 *
 * If either the printable length of the string is longer than the target
 * width or the actual length is long than the destination buffer's size,
 * it is truncated.
 *
 * \param[out] dest      pointer to destination buffer
 * \param[in]  src       pointer to printable string
 * \param[in]  width     desired printed length
 * \param[in]  destsize  destination buffer size
 */
/*----------------------------------------------------------------------------*/

void
cs_log_strpadl(char        *dest,
               const char  *src,
               size_t       width,
               size_t       destsize)
{
  _log_strpad(dest, src, width, destsize, 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Pretty-print int-32 based bit field to string
 *
 * \param[in]   code  value to print
 * \param[out]  buf   output buffer (must be at least 33 bytes).
 */
/*----------------------------------------------------------------------------*/

void
cs_log_binary_pp_int32(int32_t     code,
                       char     buf[33])
{
  int i;
  int32_t n = code;

  for (i = 0; i < 33; i++)
    buf[i] = ' ';
  buf[32] = '\0';
  buf[31] = '0';

  i = 31;
  while (n && i > -1) {
    if (n & 1)
      buf[i] = '1';
    else
      buf[i] = '0';
    n = n >> 1;
    i--;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print log info to a given log type.
 *
 * The format and variable arguments are similar to those of the vprintf()
 * type functions.
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  log      log file type
 * \param[in]  format   format string, as printf() and family.
 * \param[in]  arg_ptr  variable arguments list pointer
 *
 * \return number of characters printed, not counting the trailing '\0' used
 *         to end output strings
 */
/*----------------------------------------------------------------------------*/

int
cs_log_vprintf(cs_log_t     log,
               const char  *format,
               va_list      arg_ptr)
{
  int  retval;

  if (cs_glob_rank_id > 0)
    return 0;

  if (log != CS_LOG_DEFAULT) {

    if (_cs_log[log] == NULL && log != CS_LOG_DEFAULT)
      _open_log(log);

    retval = vfprintf(_cs_log[log], format, arg_ptr);

  }

  else {

    bft_printf_proxy_t *_printf_proxy = bft_printf_proxy_get();

    retval = _printf_proxy(format, arg_ptr);

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

  if (log != CS_LOG_DEFAULT) {

    if (_cs_log[log] == NULL && log != CS_LOG_DEFAULT)
      _open_log(log);

    va_start(arg_ptr, format);

    retval = vfprintf(_cs_log[log], format, arg_ptr);

    va_end(arg_ptr);

  }

  else {

    bft_printf_proxy_t *_printf_proxy = bft_printf_proxy_get();

    va_start(arg_ptr, format);

    retval = _printf_proxy(format, arg_ptr);

    va_end(arg_ptr);

  }

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

  if (log < CS_LOG_N_TYPES) {
    if (log == CS_LOG_DEFAULT)
      retval = bft_printf_flush();
    else if (_cs_log[log] != NULL)
      retval = fflush(_cs_log[log]);
  }

  else {
    for (i = 0; i < CS_LOG_N_TYPES; i++) {
      if (_cs_log[i] != NULL)
        retval = fflush(_cs_log[i]);
      if (retval != 0)
        break;
    }
    retval = bft_printf_flush();
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print a separator line in a log file
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  log  log file type
 */
/*----------------------------------------------------------------------------*/

void
cs_log_separator(cs_log_t log)
{
  int i;
  char separator[81];

  for (i = 0; i < 80; i++)
    separator[i] = '-';
  separator[80] = '\0';

  cs_log_printf(log, "%s\n", separator);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output timing data array header to a given log.
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  log           log file type
 * \param[in]  indent        indentation before first column
 * \param[in]  header_title  title for optional header line
 * \param[in]  calls         true if calls column is to be used
 */
/*----------------------------------------------------------------------------*/

void
cs_log_timer_array_header(cs_log_t     log,
                          int          indent,
                          const char  *header_title,
                          bool         calls)
{
  int title_width = 80 - 16 - indent;
  char tmp_s[4][64] =  {"", "", "", ""};

  /* Available width for title */

  if (calls)
    title_width -= 10; /* 1 field, 1 space + 9 digits */

  /* Header line if requested */

  assert(header_title != NULL);

  if (strlen(header_title) > 0)
    cs_log_strpad(tmp_s[0], _(header_title), title_width, 64);
  else
    cs_log_strpad(tmp_s[0], "", title_width, 64);

  cs_log_strpadl(tmp_s[2], _("time"), 12, 64);

  if (calls) {
    cs_log_strpadl(tmp_s[1], _("calls"), 9, 64);
    cs_log_printf(log,
                  "%*s%s %s %s\n",
                  indent, " ",
                  tmp_s[0], tmp_s[1], tmp_s[2]);
  }
  else
    cs_log_printf(log,
                  "%*s%s %s\n",
                  indent, " ", tmp_s[0], tmp_s[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output timing data block to a given log.
 *
 * If the optional array of call counters is used, only lines
 * with a number of calls greater than 0 are logged.
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  log           log file type
 * \param[in]  indent        indentation before first column
 * \param[in]  n_lines       number of lines in array, excluding header
 * \param[in]  line_titles   array of titles for data lines
 * \param[in]  calls         optional array of call counters, or NULL
 * \param[in]  time_count    array of time counters
 */
/*----------------------------------------------------------------------------*/

void
cs_log_timer_array(cs_log_t                   log,
                   int                        indent,
                   int                        n_lines,
                   const char                *line_titles[],
                   const unsigned             calls[],
                   const cs_timer_counter_t   time_count[])
{
  int i;
  int title_width = 80 - 16 - indent;
  char tmp_s[4][64] =  {"", "", "", ""};

  /* Available width for title */

  if (calls != NULL)
    title_width -= 10; /* 1 field, 1 space + 9 digits */

  /* Data lines */

  for (i = 0; i < n_lines; i++) {
    double wtime = (time_count[i]).nsec * 1.e-9;
    if (line_titles != NULL)
      cs_log_strpad(tmp_s[0], _(line_titles[i]), title_width, 64);
    else
      cs_log_strpad(tmp_s[0], "", title_width, 64);
    if (calls != NULL) {
      if (calls[i] > 0)
        cs_log_printf(log,
                      "%*s%s %9u %12.3f\n",
                      indent, " ", tmp_s[0], calls[i], wtime);
    }
    else
      cs_log_printf(log,
                    "%*s%s %12.3f\n",
                    indent, " ", tmp_s[0], wtime);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print a warning message to the warnings.log file and a copy to the
 * default log file.
 *
 * The format and variable arguments are similar to those of the printf()
 * type functions.
 *
 * In parallel, output is only handled by rank 0.
 *
 * \param[in]  format  format string, as printf() and family.
 * \param[in]  ...     variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\0' used
 *         to end output strings
 */
/*----------------------------------------------------------------------------*/

int
cs_log_warning(const char *format,
               ...)
{
  int retval = 0;
  va_list arg_ptr;

  if (cs_glob_rank_id > 0)
    return 0;

  va_start(arg_ptr, format);

  /* Print to warning log */
  cs_log_separator(CS_LOG_WARNINGS);
  retval = cs_log_vprintf(CS_LOG_WARNINGS, format, arg_ptr);
  cs_log_separator(CS_LOG_WARNINGS);
  cs_log_printf(CS_LOG_WARNINGS, "\n");

  /* Print copy to main log */
  cs_log_printf(CS_LOG_DEFAULT, "\nWarning message :\n");
  cs_log_printf(CS_LOG_DEFAULT, "-----------------\n\n");
  cs_log_vprintf(CS_LOG_DEFAULT, format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

/*-----------------------------------------------------------------------------*/



END_C_DECLS
