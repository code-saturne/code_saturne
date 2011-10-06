#ifndef _CFD_PROXY_COUPLING_H_
#define _CFD_PROXY_COUPLING_H_

//============================================================================
// Main API functions for CFD_Proxy component
//============================================================================

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

//----------------------------------------------------------------------------
// System headers
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Local headers
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//----------------------------------------------------------------------------
// Structure definitions
//----------------------------------------------------------------------------

typedef struct _cfd_proxy_coupling_t cfd_proxy_coupling_t;

//============================================================================
// Public function prototypes
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
                        int    component_id);

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
                           const char  *args[]);

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
//   args <-- arguments string
//----------------------------------------------------------------------------

void
cfd_proxy_set_args(const char  *args);

//----------------------------------------------------------------------------
// Set working directory.
//
// parameters:
//   path <-- new working directory
//----------------------------------------------------------------------------

void
cfd_proxy_set_dir(const char  *path);

//----------------------------------------------------------------------------
// Set shared library (also setting the shared library proxy mode)
//
// parameters:
//   filename <-- name of dynamic library file
//----------------------------------------------------------------------------

void
cfd_proxy_set_lib(const char  *filename);

//----------------------------------------------------------------------------
// Set executable (also setting the child/IPC proxy mode)
//
// parameters:
//   filename <-- name of executable file
//----------------------------------------------------------------------------

void
cfd_proxy_set_exe(const char  *filename);

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
                               const char  *launcher_args[]);

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
cfd_proxy_set_launcher(const char  *launcher_args);

//----------------------------------------------------------------------------
// Run full execution of the CFD code
//
// returns:
//   execution's return value (or prior error code if execution impossible);
//----------------------------------------------------------------------------

int
cfd_proxy_run_all(void);

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CFD_PROXY_COUPLING_H_ */
