#ifndef __CS_LOG_H__
#define __CS_LOG_H__

/*============================================================================
 * Program timing information
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
#include "cs_timer.h"
#include "stdarg.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/* code_saturne log file types */

typedef enum {

  CS_LOG_DEFAULT,      /* Default (main) log */
  CS_LOG_SETUP,        /* Calculation setup and options log */
  CS_LOG_PERFORMANCE,  /* Performance log */
  CS_LOG_WARNINGS,     /* Warnings log */
  CS_LOG_N_TYPES       /* Number of log file types */

} cs_log_t;

extern int cs_glob_log_frequency;

/*============================================================================
 * Public macros
 *============================================================================*/

/*============================================================================
 * Public function prototypes
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
cs_log_default_activate(bool  activate);

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
cs_log_default_is_active(void);

/*----------------------------------------------------------------------------
 * Count printable length of a character string.
 *
 * This should also include UTF-8 strings.
 *
 * parameters:
 *   str <-- pointer to printable string
 *
 * returns:
 *   printable length of character string
 *----------------------------------------------------------------------------*/

size_t
cs_log_strlen(const char  *s);

/*----------------------------------------------------------------------------
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
 *   dest     --> pointer to destination buffer
 *   str      <-- pointer to printable string
 *   width    <-- desired printed length
 *   destsize <-- destination buffer size
 *----------------------------------------------------------------------------*/

void
cs_log_strpad(char        *dest,
              const char  *src,
              size_t       width,
              size_t       destsize);

/*----------------------------------------------------------------------------
 * Pad a string on the left so that its printable length is
 * the required length.
 *
 * This allows pretty-printing with UTF-8 strings, whose actual length may be
 * larger than their printable length in the presence of multibyte characters.
 *
 * If either the printable length of the string is longer than the target
 * width or the actual length is long than the destination buffer's size,
 * it is truncated.
 *
 * parameters:
 *   dest     --> pointer to destination buffer
 *   str      <-- pointer to printable string
 *   width    <--  desired printed length
 *   destsize <-- destination buffer size
 *----------------------------------------------------------------------------*/

void
cs_log_strpadl(char        *dest,
               const char  *src,
               size_t       width,
               size_t       destsize);

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
                       char     buf[33]);

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
               va_list      arg_ptr);

/*----------------------------------------------------------------------------
 * Print log info to a given log type.
 *
 * The format and variable arguments are similar to those of the printf()
 * type functions.
 *
 * In parallel, output is only handled by rank 0.
 *
 * parameters:
 *   format <-- format string, as printf() and family.
 *   ...    <-- variable arguments based on format string.
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0' used
 *   to end output strings
 *----------------------------------------------------------------------------*/

#if defined(__GNUC__)

int
cs_log_printf(cs_log_t     log,
              const char  *format,
              ...)
  __attribute__((format(printf, 2, 3)));

#else

int
cs_log_printf(cs_log_t     log,
              const char  *format,
              ...);

#endif


/*----------------------------------------------------------------------------
 * Flush for output of cs_log_printf() with modifiable behavior.
 *
 * If the argument is set to CS_LOG_N_TYPES, all log files are flushed.
 *
 * returns:
 *   0 upon successful completion 0 is returned. Otherwise, EOF is returned
 *   and  errno  is  set  to indicate the error.
 *----------------------------------------------------------------------------*/

int
cs_log_printf_flush(cs_log_t log);

/*----------------------------------------------------------------------------
 * Print a separator line in a log file
 *
 * In parallel, output is only handled by rank 0.
 *
 * parameters:
 *   log <-- log file type
 *----------------------------------------------------------------------------*/

void
cs_log_separator(cs_log_t log);

/*----------------------------------------------------------------------------
 * Output timing data block to a given log.
 *
 * If the optional array of call counters is used, only lines
 * with a number of calls greater than 0 are logged.
 *
 * In parallel, output is only handled by rank 0.
 *
 * parameters:
 *   log          <-- log file type
 *   indent       <-- indentation before first column
 *   header_title <-- title for optional header line
 *   calls        <-- true if calls column is to be used
 *----------------------------------------------------------------------------*/

void
cs_log_timer_array_header(cs_log_t     log,
                          int          indent,
                          const char  *header_title,
                          bool         calls);

/*----------------------------------------------------------------------------
 * Output timing data block to a given log.
 *
 * If the optional array of call counters is used, only lines
 * with a number of calls greater than 0 are logged.
 *
 * In parallel, output is only handled by rank 0.
 *
 * parameters:
 *   log          <-- log file type
 *   indent       <-- indentation before first column
 *   n_lines      <-- number of lines in array, excluding header
 *   line_titles  <-- array of titles for data lines
 *   calls        <-- optional array of call counters, or NULL
 *   time_count   <-- array of time counters
 *----------------------------------------------------------------------------*/

void
cs_log_timer_array(cs_log_t                   log,
                   int                        indent,
                   int                        n_lines,
                   const char                *line_titles[],
                   const unsigned             calls[],
                   const cs_timer_counter_t   time_count[]);

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
cs_log_warning(const char *format,
               ...);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LOG_H__ */
