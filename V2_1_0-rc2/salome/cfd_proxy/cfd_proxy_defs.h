#ifndef _CFD_PROXY_DEF_H_
#define _CFD_PROXY_DEF_H_

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

#include "cs_config.h"

/*============================================================================
 * Definitions that may not always be provided directly by the system
 *============================================================================*/

#include <stdarg.h>

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if HAVE_STDDEF_H
# include <stddef.h>
#else
# include <stdlib.h>
#endif

/*
 * Usually stdint.h is included by inttypes.h, but only inttypes.h exists
 * on certain systems, such as Tru64Unix
 */

#if HAVE_STDINT_H
# include <stdint.h>
#elif HAVE_INTTYPES_H
# include <inttypes.h>
#endif

/* C99 _Bool type */

#if HAVE_STDBOOL_H
# include <stdbool.h>
#else
# if !HAVE__BOOL
#  ifdef __cplusplus
typedef bool _Bool;
#  else
typedef unsigned char _Bool;
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
# define __bool_tru_false_are_defined 1
#endif

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//============================================================================
// Type definitions
//============================================================================

// Enum for datatype description

typedef enum {
  CFD_PROXY_TYPE_char,
  CFD_PROXY_TYPE_bool,
  CFD_PROXY_TYPE_int,
  CFD_PROXY_TYPE_float,
  CFD_PROXY_TYPE_double,
  CFD_PROXY_TYPE_void
} cfd_proxy_type_t;

//=============================================================================
// Macro definitions
//=============================================================================

// System name

#if defined(__sgi__) || defined(__sgi) || defined(sgi)
#define CFD_PROXY_ARCH_IRIX_64

#elif defined(__hpux__) || defined(__hpux) || defined(hpux)
#define CFD_PROXY_ARCH_HP_UX

#elif defined(__linux__) || defined(__linux) || defined(linux)
#define CFD_PROXY_ARCH_Linux

#elif defined(__sun__) || defined(__sun) || defined(sun)
#define CFD_PROXY_ARCH_SunOS

#elif defined(__uxpv__) || defined(__uxpv) || defined(uxpv)
#define CFD_PROXY_ARCH_UNIX_System_V

#endif

// "Classical" macros

#define CFD_PROXY_ABS(a)     ((a) <  0  ? -(a) : (a))
#define CFD_PROXY_MIN(a,b)   ((a) > (b) ?  (b) : (a))
#define CFD_PROXY_MAX(a,b)   ((a) < (b) ?  (b) : (a))

// Internationalization macros (using gettext() or some similar function).

#define _(String) String
#define N_(String) String
#define textdomain(String) String
#define gettext(String) String
#define dgettext(Domain,String) String
#define dcgettext(Domain,String,Type) String
#define bindtextdomain(Domain,Directory) Domain

// Memory allocation macros

// Allocate memory for _ni items of type _type.
//
// parameters:
//   _ptr  --> pointer to allocated memory.
//   _ni   <-- number of items.
//   _type <-- element type.

#define CFDP_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) cfd_proxy_malloc(_ni, sizeof(_type), \
                                  #_ptr, __FILE__, __LINE__)

// Reallocate memory for _ni items of type _type.
//
// parameters:
//   _ptr  <->  pointer to allocated memory.
//   _ni   <-- number of items.
//   _type <-- element type.

#define CFDP_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) cfd_proxy_realloc(_ptr, _ni, sizeof(_type), \
                                   #_ptr, __FILE__, __LINE__)

// Free allocated memory.
//
// The freed pointer is set to NULL to avoid accidental reuse.
//
// parameters:
//   _ptr  <->  pointer to allocated memory.

#define CFDP_FREE(_ptr) \
cfd_proxy_free(_ptr), _ptr = NULL

//=============================================================================
// Global variables
//=============================================================================

extern int cfd_proxy_glob_base_rank;      // Parallel rank; -1 if serial

extern char cfd_proxy_glob_build_date[];  // Build date

extern int  cfd_proxy_glob_have_mpi;      // Indicator for MPI support
extern int  cfd_proxy_glob_have_socket;   // Indicator for socket support

extern int     cfd_proxy_glob_n_components;
extern void  **cfd_proxy_glob_component;

extern size_t  cfd_proxy_glob_type_size[];  // Size associated with each type

//=============================================================================
// Function prototypes
//=============================================================================

//----------------------------------------------------------------------------
// Error handler
//----------------------------------------------------------------------------

void
cfd_proxy_error(const char  *filename,
                int          linenum,
                int          sys_err_code,
                const char  *format,
                ...);

//-----------------------------------------------------------------------------
// Print output (wrapper or replacement for printf)
//-----------------------------------------------------------------------------

int
cfd_proxy_printf(const char *format,
                 ...);

//-----------------------------------------------------------------------------
// Flush output
//-----------------------------------------------------------------------------

int
cfd_proxy_printf_flush(void);

//-----------------------------------------------------------------------------
// Print warning
//-----------------------------------------------------------------------------

void
cfd_proxy_warn(void);

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
                 int          line_num);

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
                  int          line_num);

// Free allocated memory.
//
// parameters:
//   ptr       <-- pointer to previous memory location
//
// returns:
//   NULL pointer.

void *
cfd_proxy_free(void  *ptr);

//-----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _CFD_PROXY_DEF_H_ */
