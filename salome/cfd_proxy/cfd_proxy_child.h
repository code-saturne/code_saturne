#ifndef _CFD_PROXY_CHILD_H_
#define _CFD_PROXY_CHILD_H_

//============================================================================
// Spawn a child process and establish a connection.
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

#include "cfd_proxy_defs.h"

#include "cfd_proxy_comm.h"

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

typedef struct _cfd_proxy_child_t cfd_proxy_child_t;

//============================================================================
// Global variables
//============================================================================

//============================================================================
// Public function prototypes
//============================================================================

//----------------------------------------------------------------------------
// Spawn a child process and establish a connection
//
// returns:
//   child process handle in case of success, NULL in case of error.
//----------------------------------------------------------------------------

cfd_proxy_child_t *
cfd_proxy_child_start(const char                   *path,
                      char                   *const argv[restrict],
                      char                   *const envp[restrict],
                      cfd_proxy_comm_type_t         comm_type,
                      int                           comm_verbosity);

//----------------------------------------------------------------------------
// End connection with a child process and free associated structure
//----------------------------------------------------------------------------

int
cfd_proxy_child_stop(cfd_proxy_child_t **child);

//----------------------------------------------------------------------------
// Forward all calls from the client and their responses.
//----------------------------------------------------------------------------

void
cfd_proxy_child_forward_all(cfd_proxy_child_t  *child);

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CFD_PROXY_CHILD_H_ */
