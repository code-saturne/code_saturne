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
// Main API functions for CFD_Proxy component
//============================================================================

#include "cs_config.h"

//----------------------------------------------------------------------------
// System headers
//----------------------------------------------------------------------------

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "cfd_proxy_defs.h"
#include "cfd_proxy_comm.h"
#include "cfd_proxy_child.h"
#include "cfd_proxy_loader.h"

#include "cfd_proxy_api.h"

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//============================================================================
// Local type and structure definitions
//============================================================================

typedef enum {

  CFD_PROXY_UNDEF,    // Mode not defined yet
  CFD_PROXY_LOAD,     // Load a shared library
  CFD_PROXY_SPAWN     // Execute a sub-process

} cfd_proxy_mode_t;

// Main CFD_proxy object

typedef struct {

  cfd_proxy_mode_t  mode;        // CFD code execution mode

  int               n_args;
  int               n_launcher_args;

  char             *dirname;
  char             *filename;
  char            **args;
  char            **launcher_args;

} _cfd_proxy_t ;

//============================================================================
// Static global variables
//============================================================================

static _cfd_proxy_t _proxy = {CFD_PROXY_UNDEF, 0, 0, NULL, NULL, NULL, NULL};

//============================================================================
// Private function definitions
//============================================================================

//----------------------------------------------------------------------------
// Set arguments in a structure.
//
// parameters:
//   n_args   <-- number of arguments defined
//   args     <-- array of argument strings.
//   s_n_args <-> pointer to structure's n_args
//   s_args   <-> pointer to structure's args
//----------------------------------------------------------------------------

static void
_set_args_by_list(int           n_args,
                  const char   *args[],
                  int          *s_n_args,
                  char       ***s_args)
{
  int i;
  size_t args_size = 0;

  int _n_args = *s_n_args;
  char **_args = *s_args;

  if (_args != NULL) {
    CFDP_FREE(_args[0]);
    CFDP_FREE(_args);
  }

  for (i = 0; i < n_args; i++) {
    if (args[i] == NULL)
      break;
    args_size += (strlen(args[i]) + 1);
    _n_args = i+1;
  }

  CFDP_MALLOC(_args, _n_args, char*);
  CFDP_MALLOC(_args[0], args_size, char);

  if (_n_args == 0)
    return;

  for (i = 1, args_size = strlen(args[0]) + 1;
       i < _n_args;
       i++) {
    _args[i] = _args[0] + args_size;
    args_size += (strlen(args[i]) + 1);
  }

  for (i = 0, args_size; i < _n_args; i++)
    strcpy(_args[i], args[i]);

  // Set return values

  *s_n_args = _n_args;
  *s_args = _args;
}

//----------------------------------------------------------------------------
// Set arguments in a structure, using a complete arguments string.
//
// parameters:
//   args     <-- argument string
//   s_n_args <-> pointer to structure's n_args
//   s_args   <-> pointer to structure's args
//----------------------------------------------------------------------------

static void
_set_args_by_string(const char   *args,
                    int          *s_n_args,
                    char       ***s_args)
{
  size_t i, tok_len, start_quote_id;
  char *s;
  int protected; // 0 for unprotected, 1 for protected char, 2 for substring,
                 // 3 for protected char in substring
  size_t args_size = 0, s_size = 0;

  int _n_args = *s_n_args;
  char **_args = *s_args;

  if (_args != NULL) {
    CFDP_FREE(_args[0]);
    CFDP_FREE(_args);
  }

  if (args == NULL) {
    *s_n_args = 0;
    *s_args = _args;
    return;
  }

  args_size = strlen(args);

  // Estimate number of arguments (may overestimate, not underestimate)

  for (i = 0, _n_args = 1; i < args_size; i++) {
    char c = args[i];
    if (c == ' ' || c == '\t' || c == '\n'|| c == '\r')
      _n_args++;
  }

  // Prepare tokenization of the arguments string

  CFDP_MALLOC(_args, _n_args, char*);
  CFDP_MALLOC(_args[0], args_size + 1, char);

  s = _args[0];

  // Tokenize

  _n_args = 0;

  i = 0;                      // String position marker
  start_quote_id = args_size; // Start position marker for quoted strings
                              // (unset if j == j, set if j < args_size)

  protected = 0;
  tok_len = 0;

  while (i < args_size) {

    char c = args[i];

    // Regular case, where previous character was not protected

    if (protected == 0) {

      // Protection character

      if (c == '\\')
        protected = 1;

      // Fully protected string

      else if (c == '"' || c == '\'') {
        protected = 2;
        start_quote_id = i;
      }

      // Whitespace

      else if (c == ' ' || c == '\t' || c == '\n'|| c == '\r') {
        if (tok_len > 0) { // Finish previous token
          _args[_n_args] = s + s_size;
          s[s_size + tok_len] = '\0';
          _n_args += 1;
          s_size += tok_len+1;
          tok_len = 0;
        }
      }

      // Regular characters (token)

      else {
        s[s_size + tok_len] = args[i];
        if (tok_len == 0)
          _args[_n_args] = s + s_size;
        tok_len++;
      }

    }

    // Cases where previous character was protected

    else if (protected == 1) {

      protected = 0;
      s[s_size + tok_len] = args[i];
      if (tok_len == 0)
        _args[_n_args] = s + s_size;
      tok_len++;

    }

    else if (protected == 2) {

      // Single protection character

      if (c == '\\')
        protected = 3;

      // End of string protection

      else if (c == args[start_quote_id]) {
        protected = 0;
        start_quote_id = args_size;
      }

      else {
        s[s_size + tok_len] = args[i];
        if (tok_len == 0)
          _args[_n_args] = s + s_size;
        tok_len++;
      }

    }

    else { // if (protected == 3)

      s[s_size + tok_len] = args[i];
      if (tok_len == 0)
        _args[_n_args] = s + s_size;
      tok_len++;
      protected = 2;

    }

    i+= 1;

  } /* End of loop in infix string characters */

  if (tok_len > 0) { /* Finish previous token */
    _args[_n_args] = s + s_size;
    s[s_size + tok_len] = '\0';
    _n_args += 1;
    s_size += tok_len+1;
    tok_len = 0;
  }

  if (protected == 1)
    cfd_proxy_error(__FILE__, __LINE__, 0,
                    _("Error tokenizing expression:\n"
                      "%s\n"
                      "Missing character after \\\n"),
                    args);
  else if (protected >= 2)
    cfd_proxy_error(__FILE__, __LINE__, 0,
                    _("Error tokenizing expression:\n"
                      "%s\n"
                      "Missing closing quote for subexpression:\n"
                      "%s\n"),
                    args, args + start_quote_id);

  // Resize to adjusted size (do not resize _args[0], as _args[i]
  // would then need to be updated, and extra size is limited.

  CFDP_REALLOC(_args, _n_args, char *);

  for (i = s_size; i < args_size; i++)
    s[i] = '\0';

  // Set return values

  *s_n_args = _n_args;
  *s_args = _args;
}

//----------------------------------------------------------------------------
// Set working directory.
//
// parameters:
//   path <-- new working directory
//
// returns:
//   0 on success, -1 on error
//----------------------------------------------------------------------------

static int
_change_dir(const char  *path)
{
  int retval = 0;

  if (path != NULL) {

    retval = chdir(path);

    if (retval != 0)
      cfd_proxy_error(__FILE__, __LINE__, errno,
                      _("Error setting the working directory to:\n"
                        "%s"),
                      path);
  }

  return retval;
}

//============================================================================
// Public functions
//============================================================================

//----------------------------------------------------------------------------
// Set the proxy's associated supervisable SALOME component.
//
// Multiple components may be asociated with the proxy if necessary,
// using a different id for each (0 for the first). This may be useful
// in threaded cases used for multiple couplings.
//
// parameters:
//   component    <-- pointer of type Superv_Component_i* to the
//                    supervisable SALOME component
//   component_id <-- id of component (0 by default, >= 0 if multiple
//                    components are managed in the same process space,
//                    which may be possible with multiple threads)
//----------------------------------------------------------------------------

void
cfd_proxy_set_component(void  *component,
                        int    component_id)
{
  if (component_id < 0)
    return;

  if (component_id > cfd_proxy_glob_n_components) {

    if (cfd_proxy_glob_n_components == 1)
      // pointed to stack value, see cfd_proxy_defs.c
      cfd_proxy_glob_component = NULL;

    CFDP_REALLOC(cfd_proxy_glob_component, component_id+1, void *);

  }

  cfd_proxy_glob_component[component_id] = component;
}

//----------------------------------------------------------------------------
// Set command-line arguments for execution, using an argument array.
//
// Only user arguments need to be defined (i.e. the executable file
// at argv[0] is set automatically, so only argv[1] to argv[argc-1]
// need to be defined here.
//
// parameters:
//   n_args <-- number of arguments defined
//   args   <-- array of argument strings.
//----------------------------------------------------------------------------

void
cfd_proxy_set_args_by_list(int          n_args,
                           const char  *args[])
{
  _set_args_by_list(n_args, args, &(_proxy.n_args), &(_proxy.args));
}

//----------------------------------------------------------------------------
// Set command-line arguments for execution, using a single string.
//
// Arguments are separated by whitespace. Whitespaces may be protected
// using couples of " or ' quotes, or single \ escape characters.
//
// Only user arguments need to be defined (i.e. the executable file
// at argv[0] is set automatically, so only argv[1] to argv[argc-1]
// need to be defined here.
//
// parameters:
//   args <-- arguments string.
//----------------------------------------------------------------------------

void
cfd_proxy_set_args(const char  *args)
{
  _set_args_by_string(args, &(_proxy.n_args), &(_proxy.args));
}

//----------------------------------------------------------------------------
// Set working directory.
//
// parameters:
//   path <-- new working directory
//----------------------------------------------------------------------------

void
cfd_proxy_set_dir(const char  *path)
{
  if (path != NULL) {
    CFDP_REALLOC(_proxy.dirname, strlen(path) + 1, char);
    strcpy(_proxy.dirname, path);
  }
}

//----------------------------------------------------------------------------
// Set shared library (also setting the shared library proxy mode)
//
// parameters:
//   filename <-- name of dynamic library file
//----------------------------------------------------------------------------

void
cfd_proxy_set_lib(const char  *filename)
{
#if defined(HAVE_DLOPEN)

  size_t alloc_size = strlen(filename) + 1;

  // Avoid valgrind warning on Linux which seems to be due to dlopen
  // reading in multiples of size_t by allocating a slightly larger buffer.

  if (alloc_size % sizeof(size_t))
    alloc_size += sizeof(size_t) - (alloc_size % sizeof(size_t));

  _proxy.mode = CFD_PROXY_LOAD;

  CFDP_REALLOC(_proxy.filename, alloc_size, char);
  memset(_proxy.filename, '\0', alloc_size);

  strcpy(_proxy.filename, filename);

#endif
}

//----------------------------------------------------------------------------
// Set executable (also setting the child/IPC proxy mode)
//
// parameters:
//   filename <-- name of executable file
//----------------------------------------------------------------------------

void
cfd_proxy_set_exe(const char  *filename)
{
  _proxy.mode = CFD_PROXY_SPAWN;

  CFDP_REALLOC(_proxy.filename, strlen(filename) + 1, char);
  strcpy(_proxy.filename, filename);
}

//----------------------------------------------------------------------------
// Define intermediate launcher and associated arguments (only of use for
// the child/IPC proxy mode), using an argument array.
//
// This allows running the child executable through another process, such
// as mpiexec (for parallel runs), or a debugger.
//
// parameters:
//   n_launcher_args  <-- number of arguments defined
//   launcher_args    <-- array of string arguments (the first of which
//                        should be the launcher executable file name)
//----------------------------------------------------------------------------

void
cfd_proxy_set_launcher_by_list(int          n_launcher_args,
                               const char  *launcher_args[])
{
  _set_args_by_list(n_launcher_args,
                    launcher_args,
                    &(_proxy.n_launcher_args),
                    &(_proxy.launcher_args));
}

//----------------------------------------------------------------------------
// Define intermediate launcher and associated arguments (only of use for
// the child/IPC proxy mode), using a single string.
//
// This allows running the child executable through another process, such
// as mpiexec (for parallel runs), or a debugger.
//
// Arguments are separated by whitespace. Whitespaces may be protected
// using couples of " or ' quotes, or single \ escape characters.
//
// parameters:
//   launcher_args <-- string of launcher arguments (the first of which
//                     should be the launcher executable file name)
//----------------------------------------------------------------------------

void
cfd_proxy_set_launcher(const char  *launcher_args)
{
  _set_args_by_string(launcher_args,
                      &(_proxy.n_launcher_args),
                      &(_proxy.launcher_args));
}

//----------------------------------------------------------------------------
// Run full execution of the CFD code
//
// returns:
//   execution's return value (or prior error code if execution impossible);
//----------------------------------------------------------------------------

int
cfd_proxy_run_all()
{
  int n_args = 0;
  int arg_id = 0;
  char *old_wd = NULL;

  int retval = 0;

  extern char **environ;

  if (_proxy.dirname != NULL) {

    size_t old_wd_size = 128;
    char *wd = NULL;

    // Save old working directory

    CFDP_MALLOC(old_wd, old_wd_size, char);

    wd = getcwd(old_wd, old_wd_size);
    while (wd == NULL && errno == ERANGE) {
      old_wd_size *= 2;
      CFDP_REALLOC(old_wd, old_wd_size, char);
      wd = getcwd(old_wd, old_wd_size);
    }
    if (wd != NULL)
      CFDP_REALLOC(old_wd, strlen(old_wd) + 1, char);
    else {
      cfd_proxy_warn();
      cfd_proxy_printf
        (_("Could not obtain and save current working directory.\n"));
    }

    // Switch to working directory

    retval = _change_dir(_proxy.dirname);
    if (retval != 0)
      return retval;
  }

  // Handle shared library mode
  //---------------------------

#if defined(HAVE_DLOPEN)

  if (_proxy.mode == CFD_PROXY_LOAD) {

    char **_argv = NULL;

    CFDP_MALLOC(_argv, 1 + _proxy.n_args + 1, char *);

    _argv[n_args++] = _proxy.filename;

    for (arg_id = 0; arg_id < _proxy.n_args; arg_id++)
      _argv[n_args++] = _proxy.args[arg_id];

    _argv[n_args] = NULL;

    retval = cfd_proxy_loader_init(_proxy.filename);

    if (retval == 0)
      retval = cfd_proxy_loader_run(1 + _proxy.n_args, _argv);

    cfd_proxy_loader_finalize();

    CFDP_FREE(_argv);
  }

#endif // defined(HAVE_DLOPEN)

  // Handle spawn mode
  //------------------

#if defined(HAVE_POSIX_SPAWN) || defined(HAVE_FORK_EXECVE)

  if (_proxy.mode == CFD_PROXY_SPAWN) {

    cfd_proxy_child_t *child = NULL;
    char **child_argv = NULL;

    CFDP_MALLOC(child_argv,
                _proxy.n_launcher_args + 1 + _proxy.n_args + 1,
                char *);

    for (arg_id = 0; arg_id < _proxy.n_launcher_args; arg_id++)
      child_argv[n_args++] = _proxy.launcher_args[arg_id];

    child_argv[n_args++] = _proxy.filename;

    for (arg_id = 0; arg_id < _proxy.n_args; arg_id++)
      child_argv[n_args++] = _proxy.args[arg_id];

    child_argv[n_args] = NULL;

    child = cfd_proxy_child_start(child_argv[0],
                                  child_argv, environ,
                                  CFD_PROXY_COMM_TYPE_SOCKET, 10);

    if (child != NULL) {

      cfd_proxy_child_forward_all(child);

      retval = cfd_proxy_child_stop(&child);

    }

    else
      retval = 1;

    CFDP_FREE(child_argv);
  }

#endif // if defined(HAVE_POSIX_SPAWN) || defined(HAVE_FORK_EXECVE)

  if (old_wd != NULL) {
    _change_dir(old_wd);
    CFDP_FREE(old_wd);
  }

  return retval;
}

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

