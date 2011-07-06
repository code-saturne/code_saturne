#ifndef _CFD_PROXY_LOADER_H_
#define _CFD_PROXY_LOADER_H_

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

//----------------------------------------------------------------------------
// System headers
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Local headers
//----------------------------------------------------------------------------

#include "cfd_proxy_defs.h"

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

//============================================================================
// Global variables
//============================================================================

//============================================================================
// Public function prototypes
//============================================================================

//----------------------------------------------------------------------------
// Load a shared library
//----------------------------------------------------------------------------

int
cfd_proxy_loader_init(const char  *lib_name);

//----------------------------------------------------------------------------
// Unload a shared library and free associated structure
//----------------------------------------------------------------------------

void
cfd_proxy_loader_finalize(void);

//----------------------------------------------------------------------------
// Run a shared library's main call sequence
//----------------------------------------------------------------------------

int
cfd_proxy_loader_run(int argc, char **argv);

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CFD_PROXY_LOADER_H_ */
