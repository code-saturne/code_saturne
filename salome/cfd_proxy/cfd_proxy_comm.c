//============================================================================
// Definitions of base communication functions
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

#include "cs_config.h"

#if defined(__linux__)
# define _GNU_SOURCE 1
#endif

// System headers

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(HAVE_SOCKET)
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

#if defined(HAVE_POLL)
#include <poll.h>
#endif

// Local headers

#include "cfd_proxy_defs.h"
#include "cfd_proxy_comm.h"

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//============================================================================
//  Constants and Macros
//============================================================================

#define CFD_PROXY_COMM_FILE_NUM_LEN        4

#define NE_COMM_SOCKET_HEADER       "NE_comm_socket"

#define NE_COMM_L_HOSTNAME          256
#define NE_COMM_L_NOM_MAX           256

#define CFD_PROXY_COMM_L_TYPE_NAME         2
#define CFD_PROXY_COMM_L_SEC_NUM           4

// If SSIZE_MAX is not defined by the sytem headers, we take the minimum value
// required by POSIX (for low level reads/writes with sockets).

#if !defined(SSIZE_MAX)
# define SSIZE_MAX  32767
#endif

//============================================================================
// Structure definition
//============================================================================

struct _cfd_proxy_comm_t {

  char                  *name;          // Communicator name

#if defined(HAVE_SOCKET)
  int                    socket;        // Socket number
  int                    server_socket; // Server socket number
  struct sockaddr_in     sock_addr;     // Server socket structure
#endif

  int                    key;           // Key associated with connection
  bool                   swap_endian;   // Force big-endian communications

  cfd_proxy_comm_type_t  type;          // Communicator type
  int                    echo;          // Log (printout) level of communications

  int                    n_sec_elts;    // Number of elements in a section
                                        // (for read mode)
  int                    status;        // Error status (0 if no error,
                                        // -1 for end-of-file, -2 for error)
};

//============================================================================
// Global variables
//============================================================================

//============================================================================
// Private function definitions
//============================================================================

// Convert data from "little-endian" to "big-endian" or the reverse.
//
// The memory areas pointed to by src and dest should overlap either
// exactly or not at all.
//
// parameters:
//   dest <-> pointer to converted data location.
//   src  <-- pointer to source data location.
//   size <-- size of each item of data in bytes.
//   ni   <-- number of data items.

static void
_swap_endian(void *const          dest,
             const void    *const src,
             const size_t         size,
             const size_t         ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)dest;
  const unsigned char  *psrc = (const unsigned char *)src;

  for (i = 0 ; i < ni ; i++) {

    shift = i * size;

    for (ib = 0 ; ib < (size / 2) ; ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }

  }

  if (dest != src && size == 1)
    memcpy(dest, src, ni);
}

//----------------------------------------------------------------------------
// Destroy a communicator
//----------------------------------------------------------------------------

static void
_comm_destroy(cfd_proxy_comm_t  **comm)
{
  cfd_proxy_comm_t  *_comm = *comm;

  CFDP_FREE(_comm->name);
  CFDP_FREE(*comm);
}

#if defined(HAVE_SOCKET)

//----------------------------------------------------------------------------
// Read a record from an interface socket
//----------------------------------------------------------------------------

static void
_comm_read_sock(cfd_proxy_comm_t  *comm,
                void              *rec,
                size_t             size,
                size_t             count)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  size_t   n_bytes;
  ssize_t  ret;
  char    *_rec = rec;

  assert(rec  != NULL);
  assert(comm != NULL);
  assert(comm->socket > -1);

  n_bytes = size * count;

  // Read record from socket

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CFD_PROXY_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    ret = read(comm->socket, _rec + start_id, n_loc);

    if (ret <= 4) {
      if (strncmp(_rec + start_id, "EOF", ret) == 0)
        comm->status = -1;
    }
    else if (ret < 1 && n_loc > 0) {
      cfd_proxy_error(__FILE__, __LINE__, errno,
                      _("Communication %s:\n"
                        "Error receiving data through socket"),
                comm->name);
      comm->status = -2;
    }

    if (comm->status != 0)
      break;

    start_id += ret;

  }

  if (comm->swap_endian == true && size > 1)
    _swap_endian(rec, rec, size, count);
}

//----------------------------------------------------------------------------
// Write a record to an interface socket
//----------------------------------------------------------------------------

static void
_comm_write_sock(cfd_proxy_comm_t  *comm,
                 const void        *rec,
                 size_t             size,
                 size_t             count)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  size_t   n_bytes;
  ssize_t  ret;

  char        *_rec_swap = NULL;
  const char  *_rec = rec;

  assert(rec  != NULL);
  assert(comm != NULL);
  assert(comm->socket > -1);

  // Determine associated size

  n_bytes = size * count;

  // Convert if "little-endian"

  if (comm->swap_endian == true && size != 1) {
    CFDP_MALLOC(_rec_swap, n_bytes, char);
    _swap_endian(_rec_swap, rec, size, count);
  }

  // Write record to socket

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CFD_PROXY_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    if (_rec_swap == NULL)
      ret = write(comm->socket, _rec + start_id, n_loc);
    else
      ret = write(comm->socket, _rec_swap + start_id, n_loc);

    if (ret < 1 && n_loc > 0) {
      cfd_proxy_error(__FILE__, __LINE__, errno,
                      _("Communication %s:\n"
                        "Error sending data through socket"),
                comm->name);
      comm->status = -2;
      break;
    }

    start_id += ret;
  }

  if (_rec_swap != NULL)
    CFDP_FREE(_rec_swap);
}

//----------------------------------------------------------------------------
//  Initialize socket communication
//----------------------------------------------------------------------------

static void _comm_sock_init(cfd_proxy_comm_t  *comm)
{
  int  port_num;

#define _MAX_HOSTNAME_LEN 128
  char   namestr[_MAX_HOSTNAME_LEN + 1];

#if defined(__linux__)
  socklen_t sock_len;
#else
  int       sock_len;  // size_t by SUS-v2 standard, but according to the
                       // Linux gethostbyname man page, the standard is
                       // bad, we should have an int (or socklen_t).
#endif

  struct sockaddr_in  *sock_addr = &(comm->sock_addr);
  struct hostent      *host_ent;

  // Create server socket

  comm->server_socket = socket(AF_INET, SOCK_STREAM, 0);

  if (comm->server_socket == -1) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
                    _("Error initializing server socket.\n"));
    comm->status = -1;
  }

  // Prepare for use

  sock_len = sizeof(comm->sock_addr);

  memset((char *) sock_addr, 0, sock_len);

  sock_addr->sin_family = AF_INET;
  sock_addr->sin_addr.s_addr = INADDR_ANY;
  sock_addr->sin_port = 0;

  if (comm->swap_endian == true) {
    _swap_endian(&(sock_addr->sin_addr.s_addr),
                 &(sock_addr->sin_addr.s_addr),
                 sizeof(sock_addr->sin_addr.s_addr),
                 1);
    _swap_endian(&(sock_addr->sin_port),
                 &(sock_addr->sin_port),
                 sizeof(sock_addr->sin_port),
                 1);
  }

  if (gethostname(namestr, _MAX_HOSTNAME_LEN) < 0)
    cfd_proxy_error(__FILE__, __LINE__, errno,
                    _("Error obtaining machine name."));
  namestr[_MAX_HOSTNAME_LEN] = '\0';

  host_ent = gethostbyname(namestr);

  if (host_ent == NULL)
    host_ent = gethostbyname("localhost");

  if (host_ent == NULL) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
                    _("Error initializing server socket.\n"));
    comm->status = -1;
  }

  memcpy(host_ent->h_addr, &(sock_addr->sin_addr), host_ent->h_length);

  if (bind(comm->server_socket, (struct sockaddr *)sock_addr, sock_len) != 0) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("Error initializing server socket.\n"));
    comm->status = -1;
  }

  if (listen(comm->server_socket, 1) < 0) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("Error initializing server socket.\n"));
    comm->status = -1;
  }

  // Get service number

  if (getsockname(comm->server_socket,
                  (struct sockaddr *)sock_addr,
                  &sock_len) != 0) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("Error initializing server socket.\n"));
    comm->status = -1;
  }

  port_num = sock_addr->sin_port;
  if (comm->swap_endian == true) {
    _swap_endian(&(sock_addr->sin_port),
                 &(sock_addr->sin_port),
                 sizeof(sock_addr->sin_port), 1);
    port_num = sock_addr->sin_port;
    _swap_endian(&(sock_addr->sin_port),
                 &(sock_addr->sin_port),
                 sizeof(sock_addr->sin_port), 1);
  }

  // Build name for IP socket

  CFDP_MALLOC(comm->name, strlen(namestr) + sizeof(port_num)*8 + 2, char);
  sprintf(comm->name, "%s:%d", namestr, port_num);
  CFDP_REALLOC(comm->name, strlen(comm->name) + 1, char);

#undef _MAX_HOSTNAME_LEN
}

//----------------------------------------------------------------------------
// Close an interface socket
//----------------------------------------------------------------------------

static void
_comm_sock_disconnect(cfd_proxy_comm_t  *comm)
{
  if (close(comm->socket) != 0) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("Communication %s:\n"
                "Error closing socket.\n"),
              comm->name);
    comm->status = -1;
  }

  comm->socket = -1;
}

//----------------------------------------------------------------------------
// Connection for socket initialization
//----------------------------------------------------------------------------

static void
_comm_sock_connect(cfd_proxy_comm_t  *comm,
                   const char        *magic_string)
{
  int32_t  read_key;
  int      poll_ret = 1;
  int      poll_timeout = 10000; // in milliseconds

#if defined(CFD_PROXY_ARCH_Linux)
  socklen_t  sock_len;
#else
  size_t     sock_len;
#endif

  char *cmp_string = NULL;
  size_t string_len = strlen(magic_string);

  // Connection to server socket

  sock_len = sizeof(comm->sock_addr);

#if defined(HAVE_POLL)
  {
    struct pollfd poll_event;

    poll_event.fd = comm->server_socket;
    poll_event.events = POLLIN | POLLPRI;

    poll_ret = poll(&poll_event, 1, poll_timeout);
  }
#endif

  if (poll_ret <= 0) {
    cfd_proxy_printf
      (_("[error]\n"
         "Socket communication: timeout (%f s) connecting client.\n"),
       (double)(poll_timeout/1000.));
    comm->status = -1;
    return;
  }

  comm->socket = accept(comm->server_socket,
                        (struct sockaddr *)&(comm->sock_addr),
                        &sock_len);

  if (comm->socket < 0) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("[error]\n"
                "Socket communication: error connecting client."));
    comm->status = -1;
    return;
  }

  // Verify connection

  if (read(comm->socket, &read_key, 4) < 4) {
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("[error]\n"
                "Error in socket communication\n"));
    comm->status = -1;
    return;
  }

  if (comm->swap_endian == true)
    _swap_endian(&read_key, &read_key, 4, 1);

  if (read_key != comm->key) {
    _comm_sock_disconnect(comm);
    cfd_proxy_error(__FILE__, __LINE__, errno,
              _("[error]\n"
                "Socket connection attempt with incorrect key (%d)\n"),
              read_key);
    comm->status = -1;
    return;
  }

  // Read then write "magic string"

  CFDP_MALLOC(cmp_string, string_len + 1, char);

  _comm_read_sock(comm, cmp_string, 1, string_len);
  cmp_string[string_len] = '\0';

  _comm_write_sock(comm, magic_string, 1, string_len);

  // If the magic string does not correspond, we have an error

  if (strcmp(cmp_string, magic_string) != 0) {
    cfd_proxy_error(__FILE__, __LINE__, 0,
              _("Error reading: \"%s\".\n"
                "The interface format version is incompatible.\n"
                "The magic string indicates the interface format version:\n"
                "magic string read:     \"%s\"\n"
                "expected magic string: \"%s\""),
              comm->name, cmp_string, magic_string);
    _comm_sock_disconnect(comm);
    comm->status = -1;
    return;
  }

  CFDP_FREE(cmp_string);

  return;
}

#endif /* defined(HAVE_SOCKET) */

//============================================================================
// Public function definitions
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
                          int                    echo)
{
  unsigned  int_endian;

  cfd_proxy_comm_t *comm;

  CFDP_MALLOC(comm, 1, cfd_proxy_comm_t);

  // Initialize fields

  srandom(time(NULL));

  comm->name = NULL;
  comm->type = type;
  comm->key = random();

#if defined(HAVE_SOCKET)

  comm->socket = -1;
  comm->server_socket = -1;
  memset((char *) &(comm->sock_addr), 0, sizeof(comm->sock_addr));

#endif

  comm->echo = echo;
  comm->n_sec_elts = 0;

  comm->status = 0;

  // Test if system is big-endian

  comm->swap_endian = false; // Use "big-endian" mode to communicate

  int_endian = 0;
  *((char *) (&int_endian)) = '\1';

  if (int_endian == 1)
    comm->swap_endian = true;

#if defined(DEBUG) && !defined(NDEBUG)
  else {
    int_endian = 0;
    *((char *) (&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert (int_endian == 1);
  }
#endif

  // Info on interface creation

  if (comm->name != NULL)
    cfd_proxy_printf(_("\n"
                       "Initializing communication:  %s ..."),
                     comm->name);
  else
    cfd_proxy_printf(_("Initializing communication ..."));
  cfd_proxy_printf_flush();

  // Initialize interface

  if (type == CFD_PROXY_COMM_TYPE_SOCKET) {

#if defined(HAVE_SOCKET)

    _comm_sock_init(comm);
    if (comm->status != 0) {
      _comm_destroy(&comm);
      return NULL;
    }

#else

    cfd_proxy_error
      (__FILE__, __LINE__, 0,
       _("Library compiled without sockets support, so the communicator\n"
         "type argument to cfd_proxy_comm_initialize() must different\n"
         "from CFD_PROXY_COMM_TYPE_SOCKET (%d)."),
       (int)CFD_PROXY_COMM_TYPE_SOCKET);

#endif

  }

  cfd_proxy_printf("[ok]\n");
  cfd_proxy_printf_flush();

  // End

  return comm;
}

//----------------------------------------------------------------------------
// Finalize a communicator
//----------------------------------------------------------------------------

cfd_proxy_comm_t *
cfd_proxy_comm_finalize(cfd_proxy_comm_t *comm)
{
  cfd_proxy_printf("\n");

  // Info on closing interface files

  cfd_proxy_printf(_("Closing communication: %s\n"),
             comm->name);

#if defined(HAVE_SOCKET)
  if (comm->socket > -1)
    _comm_sock_disconnect(comm);
#endif

  CFDP_FREE(comm->name);

  CFDP_FREE(comm);

  return NULL;
}

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
                       const char        *magic_string)
{
  int retval = 0;

  if (comm->name != NULL)
    cfd_proxy_printf(_("Connecting client on:  %s ..."), comm->name);
  else
    cfd_proxy_printf(_("Connecting client ..."));

  cfd_proxy_printf_flush();

#if defined(HAVE_SOCKET)

  _comm_sock_connect(comm, magic_string);
  retval = comm->status;

#endif

  if (retval == 0)
    cfd_proxy_printf("[ok]\n");
  cfd_proxy_printf_flush();

  return retval;
}

//----------------------------------------------------------------------------
// Get a communicator's name
//
// This function returns a pointer to an existing name, so the string
// returned should not be freed by the user.
//----------------------------------------------------------------------------

const char *
cfd_proxy_comm_get_name(const cfd_proxy_comm_t *comm)
{
  return comm->name;
}

//----------------------------------------------------------------------------
// Get a communicator's associated connection key
//----------------------------------------------------------------------------

int
cfd_proxy_comm_get_key(const cfd_proxy_comm_t *comm)
{
  return comm->key;
}

//----------------------------------------------------------------------------
// Initialize a function request structure.
//
// parameters:
//   r <-> pointer to function request structure
//----------------------------------------------------------------------------

void
cfd_proxy_comm_init_request(cfd_proxy_comm_request_t *r)
{
  int i;

  for (i = 0; i < 32; i++)
    r->func_name[i] = '\0';

  r->comp_id = -1;

  r->n_ints = 0;
  r->n_doubles = 0;
  r->n_strings = 0;

  r->ints_size = 8;
  r->doubles_size = 4;
  r->strings_size = 512;

  CFDP_MALLOC(r->int_vals, r->ints_size, int);
  CFDP_MALLOC(r->double_vals, r->doubles_size, double);
  CFDP_MALLOC(r->string_vals, r->strings_size, char);
}

//----------------------------------------------------------------------------
// Finalize a function request structure.
//
// parameters:
//   r <-> pointer to function request structure
//----------------------------------------------------------------------------

void
cfd_proxy_comm_finalize_request(cfd_proxy_comm_request_t *r)
{
  int i;

  for (i = 0; i < 32; i++)
    r->func_name[i] = '\0';

  r->comp_id = -1;

  r->n_ints = 0;
  r->n_doubles = 0;
  r->n_strings = 0;

  r->ints_size = 0;
  r->doubles_size = 0;
  r->strings_size = 0;

  CFDP_FREE(r->int_vals);
  CFDP_FREE(r->double_vals);
  CFDP_FREE(r->string_vals);
}

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
                     size_t             count)
{
  assert(comm != NULL);

#if defined(HAVE_SOCKET)
  if (comm->socket > -1)
    _comm_write_sock(comm, rec, size, count);
#endif

  return comm->status;
}

//----------------------------------------------------------------------------
// Read a record from a client.
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
                    size_t             count)
{
  assert(comm != NULL);

#if defined(HAVE_SOCKET)
  if (comm->socket > -1)
    _comm_read_sock(comm, rec, size, count);
#endif

  return comm->status;
}

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
                            cfd_proxy_comm_request_t  *r)
{
  int i;
  char _base_header[256];
  char *_header = _base_header;

  size_t s_len = 0;
  size_t s_pos = 0;
  size_t header_size = 0;
  size_t header_pos = 0;
  size_t block_size = 256;

  // Read initial part of message

  cfd_proxy_comm_read(comm, _header, 1, block_size);

  if (comm->status < 0)
    return (comm->status);

  strncpy(r->func_name, _header, 32);

  // Unpack integer values

  if (comm->swap_endian == true)
    _swap_endian(_header + 32, _header + 32, 4, 5);

  {
    int32_t *_val_p;

    _val_p = (int32_t *)(_header + 32);
    r->comp_id = *_val_p;
    _val_p = (int32_t *)(_header + 36);
    header_size = *_val_p;
    _val_p = (int32_t *)(_header + 40);
    r->n_ints = *_val_p;
    _val_p = (int32_t *)(_header + 44);
    r->n_doubles = *_val_p;
    _val_p = (int32_t *)(_header + 48);
    r->n_strings = *_val_p;
  }

  header_pos = 52;

  if (header_size > 256) {
    block_size = 256*((header_size/256) + 1);
    CFDP_MALLOC(_header, block_size, char);
    memcpy(_header, _base_header, 256);
    cfd_proxy_comm_read(comm, _header+256, 1, block_size-256);
    if (comm->status < 0) {
      CFDP_FREE(_header);
      return (comm->status);
    }
  }

  // Unpack integer values

  if (r->n_ints > r->ints_size) {
    r->ints_size = r->n_ints;
    CFDP_REALLOC(r->int_vals, r->ints_size, int);
  }

  if (comm->swap_endian == true)
    _swap_endian(_header + header_pos, _header + header_pos, 4, r->n_ints);

  for (i = 0; i < r->n_ints; i++) {
    int32_t *_val_p = (int32_t *)(_header + header_pos);
    r->int_vals[i] = *_val_p;
    header_pos += 4;
  }

  // Unpack double values

  assert(sizeof(double) == 8);

  if (r->n_doubles > r->doubles_size) {
    r->doubles_size = r->n_doubles;
    CFDP_REALLOC(r->double_vals, r->doubles_size, double);
  }

  if (comm->swap_endian == true)
    _swap_endian(_header + header_pos, _header + header_pos, 8, r->n_doubles);

  for (i = 0; i < r->n_doubles; i++) {
    double *_val_p = (double *)(_header + header_pos);
    r->double_vals[i] = *_val_p;
    header_pos += 8;
  }

  // Unpack strings

  for (i = 0; i < r->n_strings; i++) {
    s_len += strlen(_header + header_pos) + 1;
    if ((size_t)(r->strings_size) < s_len) {
      r->strings_size = s_len;
      CFDP_REALLOC(r->string_vals, r->strings_size, char);
    }
    strcpy(r->string_vals + s_pos, _header + header_pos);
    s_pos = s_len;
    header_pos += strlen(_header + header_pos) + 1;
  }

  // Free allocated memory

  if (_header != _base_header)
    CFDP_FREE(_header);

  return (comm->status);
}

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
                              const char        *string_vals[])
{
  int i;
  char _base_header[256];
  char *_header = _base_header;

  size_t header_pos = 0;
  size_t header_size = 0;
  size_t block_size = 256;

  /* Compute request length, allocate if necessary */

  header_size = 4 + 4 + 4*n_ints + 8*n_doubles;

  for (i = 0; i < n_strings; i++)
    header_size += strlen(string_vals[i]) + 1;

  if (header_size > 256) {
    block_size = 256*((header_size/256) + 1);
    CFDP_MALLOC(_header, block_size, char);
  }

  // Form request header

  memset(_header, block_size, '\0');

  // Pack integer values

  {
    int32_t *_val_p;

    _val_p = (int32_t *)(_header);
    *_val_p = retcode;
    _val_p = (int32_t *)(_header + 4);
    *_val_p = header_size;
    _val_p = (int32_t *)(_header + 8);
    *_val_p = n_ints;
    _val_p = (int32_t *)(_header + 12);
    *_val_p = n_doubles;
    _val_p = (int32_t *)(_header + 16);
    *_val_p = n_strings;

    header_pos = 20; // 5*4

    for (i = 0; i < n_ints; i++) {
      _val_p = (int32_t *)(_header + header_pos);
      *_val_p = int_vals[i];
      header_pos += 4;
    }
  }

  if (comm->swap_endian == true)
    _swap_endian(_header, _header, 4, header_pos/4);

  // Pack double values

  assert(sizeof(double) == 8);

  {
    int header_pos_start = header_pos;
    double *_val_p;

    for (i = 0; i < n_doubles; i++) {
      _val_p = (double *)(_header + header_pos);
      *_val_p = double_vals[i];
      header_pos += 8;
    }

    if (comm->swap_endian == true)
     _swap_endian(_header + header_pos_start,
                  _header + header_pos_start,
                  8,
                  (header_pos - header_pos_start)/8);

  }

  /* Pack strings */

  for (i = 0; i < n_strings; i++) {

    strcpy(_header + header_pos, string_vals[i]);
    header_pos += strlen(string_vals[i]);
    _header[header_pos] = '\0';

  }

  /* Write message */

  cfd_proxy_comm_write(comm, _header, 1, block_size);

  /* Free allocated memory */

  if (_header != _base_header)
    CFDP_FREE(_header);

  return (comm->status);
}

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

