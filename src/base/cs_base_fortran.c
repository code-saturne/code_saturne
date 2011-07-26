/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Initialization and handling of Fortran-related mechanisms
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_file.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_base_fortran.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print a message to standard output
 *----------------------------------------------------------------------------*/

static int
_bft_printf(const char     *const format,
            va_list         arg_ptr)
{
 cs_int_t  line;
 cs_int_t  msgsize;

 /* Buffer for printing from C code: print to a character string, which will
    be printed to file by Fortran code.
    If Fortran output is completely replaced by C output in the future,
    we will be able to remove this step, but it is currently necessary
    so as to ensure that the same output files may be used and output
    remains ordered. */

#undef CS_BUF_PRINT_F_SIZE
#define CS_BUF_PRINT_F_SIZE 16384

 static char cs_buf_print_f[CS_BUF_PRINT_F_SIZE];

 /* Write to buffer */

#if (__STDC_VERSION__ < 199901L)
  msgsize = vsprintf (cs_buf_print_f, format, arg_ptr);
#else
  msgsize = vsnprintf (cs_buf_print_f, CS_BUF_PRINT_F_SIZE, format, arg_ptr);
#endif

  line = __LINE__ - 1;

  if (msgsize == -1 || msgsize > CS_BUF_PRINT_F_SIZE - 1) {
    fprintf(stderr,
            _("Fatal error: bft_printf() called on a message of size %d\n"
              "whereas the print buffer is of size %d."),
            msgsize, CS_BUF_PRINT_F_SIZE);

    /* Try to force segmentation fault (to call signal handlers);
       as stack has most likely been corrupted, this is the most
       "similar" error that allows for portable handling. */
    {
      int *_force_err = NULL;
      *_force_err = 0;
    }
    cs_exit(EXIT_FAILURE);
  }

  /* Effective output by Fortran code */

  CS_PROCF (csprnt, CSPRNT) (cs_buf_print_f, &msgsize);

  return msgsize;
}

/*----------------------------------------------------------------------------
 * Flush log output buffer
 *----------------------------------------------------------------------------*/

static int
_bft_printf_flush(void)
{
  CS_PROCF (csflsh, CSFLSH) ();

  return 0;
}

/*----------------------------------------------------------------------------
 * Close Fortran log files.
 *----------------------------------------------------------------------------*/

static void
_close_log_files(void)
{
  CS_PROCF(csclli, CSCLLI)();
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a directory, or check it exists.
 *
 * Fortran interface
 *
 * SUBROUTINE CSMKDR (DIRNAM, DIRLEN)
 * *****************
 *
 * CHARACTER*       DIRNAM      : --> : Directory name
 * INTEGER          DIRLEN      : --> : Directory name length
 *----------------------------------------------------------------------------*/

void CS_PROCF (csmkdr, CSMKDR)
(
 const char       *dirnam,
 const cs_int_t   *dirlen
)
{
  char    *bufname;

  /* Handle name for C API */

  bufname = cs_base_string_f_to_c_create(dirnam, *dirlen);

  if (cs_file_mkdir_default(bufname) == 1)
    bft_error(__FILE__, __LINE__, 0,
              _("The directory %s cannot be created"), bufname);

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&bufname);
}

/*----------------------------------------------------------------------------
 * Copy a Fortan string buffer to a C string buffer
 *
 * The aim of this function is to aviod issues with Fortran array bounds
 * checking when compilers such as icc 11 consider a character array from C
 * as an array of 1-character length strings.
 *
 * Fortran interface
 *
 * SUBROUTINE CSSF2C (LEN, CSTR, FSTR)
 * *****************
 *
 * INTEGER          LEN         : --> : String length
 * CHARACTER*       FSTR        : --> : Fortran string
 * CHARACTER*       CSTR        : <-- : C string
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssf2c, CSSF2C)
(
 const cs_int_t   *len,
 const char       *fstr,
 char             *cstr
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  memcpy(cstr, fstr, *len);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This is necessary for good consistency of messages output from C or
 * from Fortran, and to handle parallel and serial logging options.
 *
 * parameters:
 *   r0_log_flag <-- redirection for rank 0 log;
 *                   0: not redirected; 1: redirected to "listing" file
 *   rn_log_flag <-- redirection for ranks > 0 log:
 *                   0: not redirected; 1: redirected to "listing_n*" file;
 *                   2: redirected to "/dev/null" (suppressed)
 *----------------------------------------------------------------------------*/

void
cs_base_fortran_bft_printf_set(int r0_log_flag,
                               int rn_log_flag)
{
  cs_int_t _rank_id = cs_glob_rank_id;
  cs_int_t _n_ranks = cs_glob_n_ranks;
  cs_int_t _ilisr0 = r0_log_flag;
  cs_int_t _ilisrn = rn_log_flag;

  bft_printf_proxy_set(_bft_printf);
  bft_printf_flush_proxy_set(_bft_printf_flush);
  ple_printf_function_set(_bft_printf);

  /* Open Fortran log files */

  CS_PROCF(csopli, CSOPLI)(&_rank_id,
                           &_n_ranks,
                           &_ilisr0,
                           &_ilisrn);

  /* Close Fortran log files upon exit */

  atexit(_close_log_files);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
