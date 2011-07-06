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
// Load a shared library and configure its settings.
//============================================================================

#include "cs_config.h"

// System headers

#include <errno.h>
#include <stdio.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

// SALOME headers

#include <calcium.h>

// Local headers

#include "cfd_proxy_defs.h"

#include "cfd_proxy_loader.h"

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


//============================================================================
//                          Function prototypes
//============================================================================

// CALCIUM Function pointer types

typedef int
(cs_calcium_connect_t) (void  *component,
                        char  *s);

typedef int
(cs_calcium_disconnect_t) (void  *component,
                           int    cont);

typedef int
(cs_calcium_read_int_t)(void    *component,
                        int      time_dep,
                        float   *min_time,
                        float   *max_time,
                        int     *iteration,
                        char    *var_name,
                        int      n_val_max,
                        int     *n_val_read,
                        int      val[]);

typedef int
(cs_calcium_read_float_t)(void    *component,
                          int      time_dep,
                          float   *min_time,
                          float   *max_time,
                          int     *iteration,
                          char    *var_name,
                          int      n_val_max,
                          int     *n_val_read,
                          float    val[]);

typedef int
(cs_calcium_read_double_t)(void    *component,
                           int      time_dep,
                           double  *min_time,
                           double  *max_time,
                           int     *iteration,
                           char    *var_name,
                           int      n_val_max,
                           int     *n_val_read,
                           double   val[]);

typedef int
(cs_calcium_nb_read_int_t)(void    *component,
                           int      time_dep,
                           float   *min_time,
                           float   *max_time,
                           int     *iteration,
                           char    *var_name,
                           int      n_val_max,
                           int     *n_val_read,
                           int      val[]);

typedef int
(cs_calcium_nb_read_float_t)(void    *component,
                             int      time_dep,
                             float   *min_time,
                             float   *max_time,
                             int     *iteration,
                             char    *var_name,
                             int      n_val_max,
                             int     *n_val_read,
                             float    val[]);

typedef int
(cs_calcium_nb_read_double_t)(void    *component,
                              int      time_dep,
                              double  *min_time,
                              double  *max_time,
                              int     *iteration,
                              char    *var_name,
                              int      n_val_max,
                              int     *n_val_read,
                              double   val[]);

typedef int
(cs_calcium_write_int_t)(void    *component,
                         int      time_dep,
                         float    cur_time,
                         int      iteration,
                         char    *var_name,
                         int      n_val,
                         int      val[]);

typedef int
(cs_calcium_write_float_t)(void    *component,
                           int      time_dep,
                           float    cur_time,
                           int      iteration,
                           char    *var_name,
                           int      n_val,
                           float    val[]);

typedef int
(cs_calcium_write_double_t)(void    *component,
                            int      time_dep,
                            double   cur_time,
                            int      iteration,
                            char    *var_name,
                            int      n_val,
                            double   val[]);

// Exit function prototype

typedef void (cs_exit_t)(int status);
typedef void (cs_base_exit_set_t)(cs_exit_t *exit_func);

// Prototypes of functions used to set CALCIUM function pointers

typedef void
(cs_calcium_set_connection_funcs_t)(cs_calcium_connect_t     *cp_cd_func,
                                    cs_calcium_disconnect_t  *cp_fin_func);

typedef void
(cs_calcium_set_int_rw_funcs_t)(cs_calcium_read_int_t     *cp_len_func,
                                cs_calcium_write_int_t    *cp_een_func);

typedef void
(cs_calcium_set_float_rw_funcs_t)(cs_calcium_read_float_t     *cp_lre_func,
                                  cs_calcium_write_float_t    *cp_ere_func);

typedef void
(cs_calcium_set_double_rw_funcs_t)(cs_calcium_read_double_t     *cp_ldb_func,
                                   cs_calcium_write_double_t    *cp_edb_func);


// Main function prototype

typedef int (cs_main_t)(int argc, char **argv);

//============================================================================
//                      Local structure definitions
//============================================================================

typedef struct {

  void       *handle;                      // Handle to shared library

  cs_main_t  *cs_main;                     // Main function

  bool        cs_base_exit;                // Is cs_base_exit set ?
  bool        cs_calcium_connection_funcs; // Are connect funcs set ?
  bool        cs_calcium_set_int_rw;       // Are int r/w funcs set ?
  bool        cs_calcium_set_float_rw;     // Are float r/w funcs set ?
  bool        cs_calcium_set_double_rw;    // Are double r/w funcs set ?

} _cfd_proxy_library_t ;

//============================================================================
// Global variables
//============================================================================

static bool                   _cfd_proxy_exit_env_is_set = false;
static jmp_buf                _cfd_proxy_exit_env;

static _cfd_proxy_library_t  *_cfd_proxy_library = NULL;
static bool                   _cfd_proxy_loader_has_run = false;
static int                    _cfd_proxy_loader_exit_status = 0;

//============================================================================
//                      Private Function Definitions
//============================================================================

//----------------------------------------------------------------------------
// Exit function
//----------------------------------------------------------------------------

static void
_cfd_proxy_cs_base_exit(int status)
{
  cfd_proxy_printf(_("Finished with exit status %d\n"), status);

  _cfd_proxy_loader_has_run = true;
  _cfd_proxy_loader_exit_status = status;

  if (_cfd_proxy_exit_env_is_set == true)
    longjmp(_cfd_proxy_exit_env, status);
}

//----------------------------------------------------------------------------
// Set a shared library function pointer
//----------------------------------------------------------------------------

static void *
_set_dl_function_pointer(void        *handle,
                         const char  *name,
                         bool         errors_are_verbose)
{
  void  *retval = NULL;
  char  *error = NULL;

  dlerror();    // Clear any existing error

  retval = dlsym(handle, name);
  error = dlerror();

  if (error != NULL && errors_are_verbose)
    cfd_proxy_printf(_("Error calling dlsym: %s\n"), dlerror());

  return retval;
}

//============================================================================
//                      Public Function Definitions
//============================================================================

//----------------------------------------------------------------------------
// Load a shared library
//----------------------------------------------------------------------------

int
cfd_proxy_loader_init(const char  *lib_name)
{
#if defined(HAVE_DLOPEN)

  _cfd_proxy_library_t *l = NULL;
  int retval = 0;

  CFDP_MALLOC(l, 1, _cfd_proxy_library_t);

  // Initialize structure

  l->handle = NULL;
  l->cs_main = NULL;

  l->cs_base_exit = false;
  l->cs_calcium_connection_funcs = false;
  l->cs_calcium_set_int_rw = false;
  l->cs_calcium_set_float_rw = false;
  l->cs_calcium_set_double_rw = false;

  // Load library

  l->handle = dlopen(lib_name, RTLD_LAZY);

  if (l->handle == NULL) {
    cfd_proxy_printf(_("Error loading %s: %s\n"), lib_name, dlerror());
    retval = 1;
  }

  else {

    // Functions used to set function pointers to functions
    // provided by proxy

    cs_base_exit_set_t  *cs_base_exit_set = NULL; // Set exit function
    cs_calcium_set_connection_funcs_t *cs_calcium_set_connection_funcs = NULL;
    cs_calcium_set_int_rw_funcs_t *cs_calcium_set_int_rw_funcs = NULL;
    cs_calcium_set_float_rw_funcs_t *cs_calcium_set_float_rw_funcs = NULL;
    cs_calcium_set_double_rw_funcs_t *cs_calcium_set_double_rw_funcs = NULL;

    cs_base_exit_set = _set_dl_function_pointer(l->handle,
                                                "cs_base_exit_set",
                                                true);

    cs_calcium_set_connection_funcs
      = _set_dl_function_pointer(l->handle,
                                 "cs_calcium_set_connection_funcs",
                                 true);

    cs_calcium_set_int_rw_funcs
      = _set_dl_function_pointer(l->handle,
                                 "cs_calcium_set_int_rw_funcs",
                                 true);

    cs_calcium_set_float_rw_funcs
      = _set_dl_function_pointer(l->handle,
                                 "cs_calcium_set_float_rw_funcs",
                                 true);

    cs_calcium_set_double_rw_funcs
      = _set_dl_function_pointer(l->handle,
                                 "cs_calcium_set_double_rw_funcs",
                                 true);

    if (cs_base_exit_set != NULL) {
      cs_base_exit_set(_cfd_proxy_cs_base_exit);
      l->cs_base_exit = true;
    }
    else
      retval = 2;

    if (cs_calcium_set_connection_funcs != NULL) {
      cs_calcium_set_connection_funcs(cp_cd, cp_fin);
      l->cs_calcium_connection_funcs = true;
    }
    else
      retval = 2;

    if (cs_calcium_set_int_rw_funcs != NULL) {
      cs_calcium_set_int_rw_funcs(cp_len, cp_een);
      l->cs_calcium_set_int_rw = true;
    }
    else
      retval = 2;

    if (cs_calcium_set_float_rw_funcs != NULL) {
      cs_calcium_set_float_rw_funcs(cp_lre, cp_ere);
      l->cs_calcium_set_float_rw = true;
    }
    else
      retval = 2;

    if (cs_calcium_set_double_rw_funcs != NULL) {
      cs_calcium_set_double_rw_funcs(cp_ldb, cp_edb);
      l->cs_calcium_set_double_rw = true;
    }
    else
      retval = 2;

    // Set main function

    l->cs_main =_set_dl_function_pointer(l->handle,
                                         "cs_main",
                                         true);

    if (l->cs_main == NULL)
      retval = 2;
  }

  // Handle errors

  if (retval == 2)
    dlclose(l->handle);

  if (retval != 0)
    CFDP_FREE(l);

  _cfd_proxy_library = l;

  return retval;

#else /* if !defined(HAVE_DLOPEN) */

  cfd_proxy_printf(_("Shared library support not available.\n"
                     "Unable to load: %s\n"), lib_name);

  return 1;

#endif
}

//----------------------------------------------------------------------------
// Unload a shared library and free associated structure
//----------------------------------------------------------------------------

void
cfd_proxy_loader_finalize(void)
{
  _cfd_proxy_library_t *l = _cfd_proxy_library;

  if (l == NULL)
    return;

#if defined(HAVE_DLOPEN)
  dlclose(l->handle);
#endif

  CFDP_FREE(l);
  _cfd_proxy_library = l;

  _cfd_proxy_loader_has_run = false;
  _cfd_proxy_loader_exit_status = 0;
}

//----------------------------------------------------------------------------
// Run a shared library's main call sequence
//----------------------------------------------------------------------------

int
cfd_proxy_loader_run(int argc, char **argv)
{
  _cfd_proxy_library_t *l = _cfd_proxy_library;

  if (l->cs_main != NULL) {

    _cfd_proxy_exit_env_is_set = true;

    if (setjmp (_cfd_proxy_exit_env) == 0) {

      _cfd_proxy_exit_env_is_set = true;

      l->cs_main(argc, argv);

    }
    else
      _cfd_proxy_exit_env_is_set = false;
  }

  return _cfd_proxy_loader_exit_status;
}

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

