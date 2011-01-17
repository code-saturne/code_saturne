//============================================================================
//
//     This file is part of the Code_Saturne CFD tool.
//
//     Copyright (C) 2006-2011 EDF S.A., France
//
//     contact: saturne-support@edf.fr
//
//     The Code_Saturne CFD tool is free software; you can redistribute it
//     and/or modify it under the terms of the GNU General Public License
//     as published by the Free Software Foundation; either version 2 of
//     the License, or (at your option) any later version.
//
//     The Code_Saturne CFD tool is distributed in the hope that it will be
//     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with the Code_Saturne Kernel; if not, write to the
//     Free Software Foundation, Inc.,
//     51 Franklin St, Fifth Floor,
//     Boston, MA  02110-1301  USA
//
//============================================================================

//============================================================================
// Definitions, Global variables, and basic functions
//============================================================================

// System headers

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Library headers

#include "cfd_proxy_defs.h"

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//=============================================================================
// Global variables
//=============================================================================

int cfd_proxy_glob_base_rank = - 1;           // Parallel rank; -1 if serial

char cfd_proxy_glob_build_date[] = __DATE__;  // Build date

#if defined(HAVE_MPI)
int cfd_proxy_glob_have_mpi = 1;              // Indicator for MPI support
#else
int cfd_proxy_glob_have_mpi = 0;
#endif

#if defined(HAVE_SOCKET)
int  cfd_proxy_glob_have_socket = 1;          // Indicator for socket support
#else
int  cfd_proxy_glob_have_socket = 0;
#endif

static void  *_cfd_proxy_glob_def_component[1] = {NULL};
int           cfd_proxy_glob_n_components = 1;
void        **cfd_proxy_glob_component = _cfd_proxy_glob_def_component;

// Sizes associated with each type

size_t  cfd_proxy_glob_type_size[] = {1,               // CFD_PROXY_TYPE_char
                                      sizeof(int),     // CFD_PROXY_TYPE_bool
                                      sizeof(int),     // CFD_PROXY_TYPE_int
                                      sizeof(float),   // CFD_PROXY_TYPE_float
                                      sizeof(double),  // CFD_PROXY_TYPE_float
                                      0};

//=============================================================================
// Private functions
//=============================================================================

//=============================================================================
// Public functions
//=============================================================================

//----------------------------------------------------------------------------
// Error handler
//----------------------------------------------------------------------------

void
cfd_proxy_error(const char  *filename,
                int          linenum,
                int          sys_err_code,
                const char  *format,
                ...)
{
  va_list  arg_ptr;

  fflush(stdout); fflush(stderr);

  fprintf(stderr, _("\n"
                    "Error executing CFD_Proxy\n"
                    "=========================\n")) ;

  if (sys_err_code != 0)
    fprintf(stderr, _("\nSystem error: %s\n"), strerror(sys_err_code));

  fprintf(stderr, _("\n%s:%d: Fatal error.\n\n"), filename, linenum);

  va_start(arg_ptr, format);

  vfprintf(stderr, format, arg_ptr);

  va_end(arg_ptr);

  fprintf(stderr, "\n\n");

  fflush(stderr);
}

//-----------------------------------------------------------------------------
// Print output (wrapper or replacement for printf)
//-----------------------------------------------------------------------------

int
cfd_proxy_printf(const char *format,
                 ...)
{
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  return vprintf(format, arg_ptr);

  va_end(arg_ptr);
}

//-----------------------------------------------------------------------------
// Flush output
//-----------------------------------------------------------------------------

int
cfd_proxy_printf_flush(void)
{
  return fflush(stdout);
}

//-----------------------------------------------------------------------------
// Print warning
//-----------------------------------------------------------------------------

void
cfd_proxy_warn(void)
{
  cfd_proxy_printf(_("\n"
                     "Warning (CFD_Proxy)\n"
                     "===================\n"));
  cfd_proxy_printf_flush();
}

// Allocate memory for ni elements of size bytes.
//
// parameters:
//   ni        <-- number of items.
//   size      <-- element size.
//   var_name  <-- allocated variable name string.
//   file_name <-- name of calling source file.
//   line_num  <-- line number in calling source file.
//
// returns:
//   pointer to allocated memory.

void *
cfd_proxy_malloc(size_t       ni,
                 size_t       size,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num)
{
  void       *p_loc;
  size_t      alloc_size = ni * size;

  if (ni == 0)
    return NULL;

  // Allocate memory and check return

  p_loc = malloc(alloc_size);

  if (p_loc == NULL)
    cfd_proxy_error(file_name, line_num, errno,
                    _("Failure to allocate \"%s\" (%lu bytes)"),
                    var_name, (unsigned long)alloc_size);

  return p_loc;
}

// Reallocate memory for ni elements of size bytes.
//
// parameters:
//   ptr       <-- pointer to previous memory location
//   ni        <-- number of items.
//   size      <-- element size.
//   var_name  <-- allocated variable name string.
//   file_name <-- name of calling source file.
//   line_num   -> line number in calling source file
//
// returns:
//   pointer to allocated memory.

void *
cfd_proxy_realloc(void        *ptr,
                  size_t       ni,
                  size_t       size,
                  const char  *var_name,
                  const char  *file_name,
                  int          line_num)
{
  void      *p_loc;

  size_t new_size = ni * size;

  p_loc = realloc(ptr, new_size);

  if (p_loc == NULL && new_size != 0)
    cfd_proxy_error(file_name, line_num, errno,
                    _("Failure to reallocate \"%s\" (%lu bytes)"),
                    var_name, (unsigned long)new_size);

  return p_loc;
}

// Free allocated memory.
//
// parameters:
//   ptr       <-- pointer to previous memory location
//
// returns:
//   NULL pointer.

void *
cfd_proxy_free(void  *ptr)
{
  if (ptr != NULL)
    free(ptr);

  return NULL;
}

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

