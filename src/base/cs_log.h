#ifndef __CS_LOG_H__
#define __CS_LOG_H__

/*============================================================================
 * Program timing information
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/* Code_Saturne log file types */

typedef enum {

  CS_LOG_PERFORMANCE,  /* Performance log */
  CS_LOG_N_TYPES       /* Number of log file types */

} cs_log_t;


/*============================================================================
 * Public macros
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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

int
cs_log_printf(cs_log_t     log,
              const char  *format,
              ...);

/*----------------------------------------------------------------------------
 * Print a string to a given width info to a given type.
 *
 * In parallel, output is only handled by rank 0.
 *
 * parameters:
 *   log   <-- log file type
 *   str   <-- string to print
 *   width <-- minimum width to which the string must be printed
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0' used
 *   to end output strings
 *----------------------------------------------------------------------------*/

int
cs_log_print_padded_str(cs_log_t     log,
                        const char  *str,
                        int          width);

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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LOG_H__ */
