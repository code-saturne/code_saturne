/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_PROXY_COMM_H__
#define __CS_PROXY_COMM_H__

/*============================================================================
 * Base communication functions for use with proxy
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_PROXY_COMM_CMD_ABORT                       "cmd:abort"
#define CS_PROXY_COMM_CMD_STOP                         "cmd:stop"
#define CS_PROXY_COMM_FORWARD                           "forward"

#define CS_PROXY_COMM_L_SEC_NAME                              32

/*============================================================================
 * Type and structure definitions
 *============================================================================*/

typedef enum {
  CS_PROXY_TYPE_string,
  CS_PROXY_TYPE_int,
  CS_PROXY_TYPE_float,
  CS_PROXY_TYPE_double,
} cs_proxy_type_t;

typedef enum {

  CS_PROXY_COMM_MSG_NONE,       /* Not a message */
  CS_PROXY_COMM_MSG_ABORT,      /* Emergency stop */
  CS_PROXY_COMM_MSG_STOP,       /* End of communication */
  CS_PROXY_COMM_MSG_FORWARD,    /* Message to forward */
  CS_PROXY_COMM_MSG_OTHER       /* Undefined message type */

} cs_proxy_comm_msg_t;

typedef enum {

  CS_PROXY_COMM_TYPE_SOCKET,    /* Communicate through sockets */
  CS_PROXY_COMM_TYPE_NULL       /* Null communicator */

} cs_proxy_comm_type_t;

typedef struct _cs_proxy_comm_t cs_proxy_comm_t;

/* Public structure used to save data from a message header, simplifying
   the transfer of this data to processing functions */

typedef struct {

  cs_proxy_comm_msg_t msg_type;
  int  sec_num;
  char sec_name[CS_PROXY_COMM_L_SEC_NAME + 1];
  int  n_elts;
  cs_type_t elt_type;

} cs_proxy_comm_msg_header_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Establish a connection to a proxy.
 *
 * parameters:
 *   port_name     <-- name of server port (host:port for IP sockets)
 *   key           <-- key for authentification
 *   type          <-- communication type
 *----------------------------------------------------------------------------*/

void
cs_proxy_comm_initialize(const char           *port_name,
                         int                   key,
                         cs_proxy_comm_type_t  type);

/*----------------------------------------------------------------------------
 * Finalize a connection to a proxy.
 *----------------------------------------------------------------------------*/

void
cs_proxy_comm_finalize(void);

/*----------------------------------------------------------------------------
 * Write a record to a proxy.
 *
 * parameters:
 *   rec     <-- pointer to data to write
 *   size    <-- size of each data element, in bytes
 *   count   <-- number of data elements
 *----------------------------------------------------------------------------*/

void
cs_proxy_comm_write(const void  *rec,
                    size_t       size,
                    size_t       count);

/*----------------------------------------------------------------------------
 * Read a record from a proxy.
 *
 * parameters:
 *   rec     --> pointer to data to write
 *   size    <-- size of each data element, in bytes
 *   count   <-- number of data elements
 *----------------------------------------------------------------------------*/

void
cs_proxy_comm_read(void    *rec,
                   size_t   size,
                   size_t   count);

/*----------------------------------------------------------------------------
 * Send a function-relay request to a proxy.
 *
 * parameters:
 *   func_name   <-- name of function associated with request
 *   comp_id     <-- associated component id
 *   n_ints      <-- number of integer arguments
 *   n_doubles   <-- number of floating-point arguments
 *   n_strings   <-- number of string arguments
 *   int_vals    <-- integer argument values
 *   double_vals <-- floating-point argument values
 *   string_vals <-- string argument values
 *----------------------------------------------------------------------------*/

void
cs_proxy_comm_write_request(const char      *func_name,
                            int              comp_id,
                            int              n_ints,
                            int              n_doubles,
                            int              n_strings,
                            const int        int_vals[],
                            const double     double_vals[],
                            const char      *string_vals[]);

/*----------------------------------------------------------------------------
 * Read a function-relay response from a proxy.
 *
 * Return value arrays must be large enough to receive all values.
 *
 * parameters:
 *   n_ints      <-- number of integer arguments
 *   n_doubles   <-- number of floating-point arguments
 *   n_strings   <-- number of string arguments
 *   int_vals    --> integer argument values
 *   double_vals --> floating-point argument values
 *   string_vals --> string argument values
 *
 * returns:
 *   the relayed function's return value
 *----------------------------------------------------------------------------*/

int
cs_proxy_comm_read_response(int        n_ints,
                            int        n_doubles,
                            int        n_strings,
                            int        int_vals[],
                            double     double_vals[],
                            char      *string_vals[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROXY_COMM_H__ */
