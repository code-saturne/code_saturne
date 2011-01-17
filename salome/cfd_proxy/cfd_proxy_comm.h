#ifndef _CFD_PROXY_COMM_H_
#define _CFD_PROXY_COMM_H_

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
// Definitions of base communication functions
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
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//----------------------------------------------------------------------------
// Message types
//----------------------------------------------------------------------------

typedef enum {

  CFD_PROXY_COMM_MSG_NONE,       // Not a message
  CFD_PROXY_COMM_MSG_ABORT,      // Emergency stop
  CFD_PROXY_COMM_MSG_STOP,       // End of communication
  CFD_PROXY_COMM_MSG_FORWARD,    // Message to forward
  CFD_PROXY_COMM_MSG_OTHER       // Undefined message type

} cfd_proxy_comm_msg_t;

typedef enum {

  CFD_PROXY_COMM_TYPE_SOCKET,    // Communicate through sockets
  CFD_PROXY_COMM_TYPE_NULL       // Null communicator

} cfd_proxy_comm_type_t;

//----------------------------------------------------------------------------
// Macro definitions
//----------------------------------------------------------------------------

#define CFD_PROXY_COMM_CMD_ABORT                       "cmd:abort"
#define CFD_PROXY_COMM_CMD_STOP                         "cmd:stop"
#define CFD_PROXY_COMM_FORWARD                           "forward"

#define CFD_PROXY_COMM_PROC_ALL                                -1

#define CFD_PROXY_COMM_L_SEC_NAME                              32

//----------------------------------------------------------------------------
// Structure definitions
//----------------------------------------------------------------------------

typedef struct _cfd_proxy_comm_t cfd_proxy_comm_t;

// Public structure used to save data from a message header, simplifying
// the transfer of this data to processing functions

typedef struct {

  char                   func_name[32]; // Function name

  int                    comp_id;       // Component id

  int                    n_ints;        // Number of integer values
  int                    n_doubles;     // Number of double-precision values
  int                    n_strings;     // Number of string values

  int                    ints_size;     // Size of array containing ints
  int                    doubles_size;  // Size of array containing doubles
  int                    strings_size;  // Size of array containing strings

  int                   *int_vals;      // Array of integer values
  double                *double_vals;   // Array of double-precision values
  char                  *string_vals;   // Array of string values

} cfd_proxy_comm_request_t;

//============================================================================
// Global variables
//============================================================================

//============================================================================
// Public function prototypes
//============================================================================

//----------------------------------------------------------------------------
// Initialize a communicator
//
// parameters:
//   type  <-- communication type
//   echo  <-- echo on main output
//----------------------------------------------------------------------------

cfd_proxy_comm_t *
cfd_proxy_comm_initialize(cfd_proxy_comm_type_t  type,
                          int                    echo);

//----------------------------------------------------------------------------
// Finalize a communicator
//----------------------------------------------------------------------------

cfd_proxy_comm_t *
cfd_proxy_comm_finalize(cfd_proxy_comm_t *comm);

//----------------------------------------------------------------------------
// Establish a communicator connection
//
// parameters:
//   comm          <-> communicator
//   magic_string  <-- magic string for verification
//
// returns:
//   0 in case of success, -1 otherwise;
//----------------------------------------------------------------------------

int
cfd_proxy_comm_connect(cfd_proxy_comm_t  *comm,
                       const char        *magic_string);

//----------------------------------------------------------------------------
// Get a communicator's name
//
// This function returns a pointer to an existing name, so the string
// returned should not be freed by the user.
//----------------------------------------------------------------------------

const char *
cfd_proxy_comm_get_name(const cfd_proxy_comm_t *comm);

//----------------------------------------------------------------------------
// Get a communicator's associated connection key
//----------------------------------------------------------------------------

int
cfd_proxy_comm_get_key(const cfd_proxy_comm_t *comm);

//----------------------------------------------------------------------------
// Initialize a function request structure.
//
// parameters:
//   r <-> pointer to function request structure
//----------------------------------------------------------------------------

void
cfd_proxy_comm_init_request(cfd_proxy_comm_request_t *r);

//----------------------------------------------------------------------------
// Finalize a function request structure.
//
// parameters:
//   r <-> pointer to function request structure
//----------------------------------------------------------------------------

void
cfd_proxy_comm_finalize_request(cfd_proxy_comm_request_t *r);

//----------------------------------------------------------------------------
// Write a record to a client.
//
// parameters:
//   comm    <-- communicator
//   rec     <-- pointer to data to write
//   size    <-- size of each data element, in bytes
//   count   <-- number of data elements
//
// returns:
//   0 if request was written correctly, -2 for error
//----------------------------------------------------------------------------

int
cfd_proxy_comm_write(cfd_proxy_comm_t  *comm,
                     const void        *rec,
                     size_t             size,
                     size_t             count);

//----------------------------------------------------------------------------
// Read a record from a proxy.
//
// parameters:
//   comm    <-- communicator
//   rec     --> pointer to data to write
//   size    <-- size of each data element, in bytes
//   count   <-- number of data elements
//
// returns:
//   0 if request was read correctly, -1 for end-of-file, -2 for error
//----------------------------------------------------------------------------

int
cfd_proxy_comm_read(cfd_proxy_comm_t  *comm,
                    void              *rec,
                    size_t             size,
                    size_t             count);

//----------------------------------------------------------------------------
// Read a function-relay response from a proxy.
//
// parameters:
//   comm <-> communicator
//   r    <-> pointer to function request structure
//
// returns:
//   0 if request was read correctly, -1 for end-of-file, -2 for error
//----------------------------------------------------------------------------*/

int
cfd_proxy_comm_read_request(cfd_proxy_comm_t          *comm,
                            cfd_proxy_comm_request_t  *r);

//----------------------------------------------------------------------------
// Send a function-relay response to a proxy.
//
// parameters:
//   comm        <-- communicator
//   retcode     <-- function return code
//   comp_id     <-- associated component id
//   n_ints      <-- number of integer arguments
//   n_doubles   <-- number of floating-point arguments
//   n_strings   <-- number of string arguments
//   int_vals    <-- integer argument values
//   double_vals <-- floating-point argument values
//   string_vals <-- string argument values
//
// returns:
//   0 if response was written correctly, -2 for error
//----------------------------------------------------------------------------

int
cfd_proxy_comm_write_response(cfd_proxy_comm_t  *comm,
                              int                retcode,
                              int                n_ints,
                              int                n_doubles,
                              int                n_strings,
                              const int          int_vals[],
                              const double       double_vals[],
                              const char        *string_vals[]);

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CFD_PROXY_COMM_H_ */
