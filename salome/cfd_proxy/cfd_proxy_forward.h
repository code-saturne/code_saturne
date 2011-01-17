#ifndef _CFD_PROXY_FORWARD_H_
#define _CFD_PROXY_FORWARD_H_
 
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
// Forwarding of function calls
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
// Message types
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Macro definitions
//----------------------------------------------------------------------------

#define CFD_PROXY_FORWARD_CMD_ABORT                       "cmd:abort"
#define CFD_PROXY_FORWARD_CMD_STOP                         "cmd:stop"
#define CFD_PROXY_FORWARD_FORWARD                           "forward"

#define CFD_PROXY_FORWARD_PROC_ALL                                -1

#define CFD_PROXY_FORWARD_L_SEC_NAME                              32

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
// Forward a call from the client and its response.
//
// parameters:
//   comm   <-> communicator
//   r      <-> container for temporary communication request data
//
// returns:
//   0 if request was forwarded correctly, -1 for end-of-file, -2 for error
//----------------------------------------------------------------------------

int
cfd_proxy_forward(cfd_proxy_comm_t          *comm,
                  cfd_proxy_comm_request_t  *r);

//----------------------------------------------------------------------------
// Forward all calls from the client and their responses.
//----------------------------------------------------------------------------

void
cfd_proxy_forward_all(cfd_proxy_comm_t  *comm);

//-----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CFD_PROXY_FORWARD_H_ */
