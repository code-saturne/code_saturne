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

/*============================================================================
 * Base communication functions for use with proxy
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/* System and BFT headers */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_SOCKET)
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

#include <bft_file.h>
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/* Local headers */

#include "cs_base.h"

#include "cs_proxy_comm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_PROXY_COMM_MAGIC_STRING       "CFD_Proxy_comm_socket"

#define CS_PROXY_COMM_FILE_NUM_LEN        4

#define CS_PROXY_COMM_HOSTNAME_L        256
#define CS_PROXY_COMM_NAME_L            256

#define CS_PROXY_COMM_L_TYPE_NAME         2
#define CS_PROXY_COMM_L_SEC_NUM           4

/* If SSIZE_MAX is not defined by the sytem headers, we take the minimum value
   required by POSIX (for low level reads/writes with sockets). */

#if !defined(SSIZE_MAX)
# define SSIZE_MAX  32767
#endif

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

struct _cs_proxy_comm_t {

  char                  *port_name;     /* Port name (hostname:socket
                                           for IP sockets) */

#if defined(HAVE_SOCKET)
  int                    socket;        /* Socket number */
#endif

  cs_bool_t              swap_endian;   /* Force big-endian communications */

  cs_proxy_comm_type_t   type;          /* Communicator type */

  int                    n_sec_elts;    /* Number of elements in a section
                                           (for read mode) */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_proxy_comm_t *_cs_glob_proxy_comm = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Read a record from an interface socket
 *----------------------------------------------------------------------------*/

static void _comm_read_sock(const cs_proxy_comm_t  *comm,
                            void                   *rec,
                            size_t                  size,
                            size_t                  count)
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

  /* Read record from socket */

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    ret = read(comm->socket, _rec + start_id, n_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error receiving data through socket."),
                comm->port_name);

    start_id += ret;

  }

  if (comm->swap_endian == true && size > 1)
    bft_file_swap_endian(rec, rec, size, count);
}

/*----------------------------------------------------------------------------
 * Write a record to an interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_write_sock(const cs_proxy_comm_t  *comm,
                 const void             *rec,
                 size_t                  size,
                 size_t                  count)
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

  /* Determine associated size */

  n_bytes = size * count;

  /* Convert if "little-endian" */

  if (comm->swap_endian == true && size != 1) {
    BFT_MALLOC(_rec_swap, n_bytes, char);
    bft_file_swap_endian(_rec_swap, rec, size, count);
  }

  /* Write record to socket */

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    if (_rec_swap == NULL)
      ret = write(comm->socket, _rec + start_id, n_loc);
    else
      ret = write(comm->socket, _rec_swap + start_id, n_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error sending data through socket."),
                comm->port_name);

    start_id += ret;
  }

  if (_rec_swap != NULL)
    BFT_FREE(_rec_swap);
}

/*----------------------------------------------------------------------------
 * Close an interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_sock_disconnect(cs_proxy_comm_t  *comm)
{
  if (close(comm->socket) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Communication %s:\n"
                "Error closing socket."),
              comm->port_name);

  comm->socket = -1;
}

/*----------------------------------------------------------------------------
 * Connection for socket initialization
 *----------------------------------------------------------------------------*/

static void
_comm_sock_connect(cs_proxy_comm_t  *comm)
{
  int  port_num, id;

  char  *host_name = NULL;

#if defined(_CS_ARCH_Linux)
  socklen_t  sock_len;
#else
  size_t     sock_len;
#endif

  struct sockaddr_in   sock_addr;
  struct hostent      *host_ent;

  /* Decode comm->port_name string */

  for (id = strlen(comm->port_name) - 1;
       id > 0 && comm->port_name[id] != ':'; id--);

  port_num = atoi(comm->port_name + id + 1);

  BFT_MALLOC(host_name, id + 1, char);
  strncpy(host_name, comm->port_name, id);
  host_name[id] = '\0';

  /* Establish communication with CFD_Proxy */
  /*----------------------------------------*/

  /* Create socket interface descriptor */

  comm->socket = socket(AF_INET, SOCK_STREAM, 0);

  if (comm->socket == -1)
    bft_error(__FILE__, __LINE__, errno,
              _("Error initializing socket communication."));

  /* Prepare connection */

  sock_len = sizeof(sock_addr);

  memset((char *) &sock_addr, 0, sock_len);

  sock_addr.sin_family = AF_INET;
  sock_addr.sin_addr.s_addr = inet_addr(host_name);

  if (sock_addr.sin_addr.s_addr == INADDR_NONE) {
    host_ent = gethostbyname(host_name);

    if (host_ent == NULL)
      host_ent = gethostbyname("localhost");

    if (host_ent == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Socket communication: host \"%s\" unknown."),
                host_name);

    memcpy(&sock_addr.sin_addr, host_ent->h_addr_list[0], host_ent->h_length);
  }

  sock_addr.sin_port = port_num;

  if (comm->swap_endian == true)
    bft_file_swap_endian((char *)&(sock_addr.sin_port),
                         (char *)&(sock_addr.sin_port),
                         sizeof(sock_addr.sin_port),
                         1);

  if (connect(comm->socket,
              (struct sockaddr *)&sock_addr, sock_len) < 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Socket communication: error connecting to\n"
                "%s (port %d)."), host_name, port_num);

  /* Free temporary string */

  BFT_FREE(host_name);
}

/*----------------------------------------------------------------------------
 * Initialize a socket for communication
 *----------------------------------------------------------------------------*/

static void
_comm_sock_handshake(cs_proxy_comm_t  *comm,
                     const char       *magic_string,
                     int               key)
{
  int   len = strlen(magic_string);
  char  keybuf[4] = {'\0', '\0', '\0', '\0'};
  char *str_cmp = NULL;

  /* Send key */

  assert(sizeof(int) == 4 || sizeof(short) == 4);

  if (sizeof(int) == 4)
    *((int *)&keybuf) = key;
  else if (sizeof(short) == 4)
    *((short *)&keybuf) = key;

  _comm_write_sock(comm, keybuf, 4, 1);

  /* Write "magic string" */

  _comm_write_sock(comm, magic_string, 1, len);

  /* Read same magic string */

  BFT_MALLOC(str_cmp, len + 1, char);

  _comm_read_sock(comm, str_cmp, 1, len);
  bft_printf("compare : %s\n", strcmp);

  if (strncmp(str_cmp, magic_string, len))
    bft_error(__FILE__, __LINE__, 0, _("Handshake with proxy failed."));

  BFT_FREE(str_cmp);
}

#endif /* defined(HAVE_SOCKET) */

/*----------------------------------------------------------------------------
 * Establish a communicator connection
 *
 * parameters:
 *   port_name     <-- name of server port (host:port for IP sockets)
 *   key           <-- key for authentification
 *   type          <-- communication type
 *
 * returns:
 *   pointer to initialized communicator;
 *----------------------------------------------------------------------------*/

static cs_proxy_comm_t *
_comm_initialize(const char           *port_name,
                 int                   key,
                 cs_proxy_comm_type_t  type)
{
  unsigned  int_endian;

  int retval = 0;

  cs_proxy_comm_t *comm = NULL;

  BFT_MALLOC(comm, 1, cs_proxy_comm_t);

  /* Initialize fields */

  BFT_MALLOC(comm->port_name, strlen(port_name) + 1, char);
  strcpy(comm->port_name, port_name);

  comm->type = type;

  comm->n_sec_elts = 0;

  /* Test if system is big-endian */

  comm->swap_endian = false; /* Use "big-endian" mode to communicate */

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

  /* Info on interface creation */

  if (comm->port_name != NULL)
    bft_printf(_("Connecting to proxy:  %s ..."), comm->port_name);
  else
    bft_printf(_("Connecting to proxy ..."));
  bft_printf_flush();

  /* Initialize interface */

  if (type == CS_PROXY_COMM_TYPE_SOCKET) {

#if defined(HAVE_SOCKET)

    _comm_sock_connect(comm);
    _comm_sock_handshake(comm, CS_PROXY_COMM_MAGIC_STRING, key);

#else

    bft_printf("\n");
    bft_error
      (__FILE__, __LINE__, 0,
       _("Library compiled without sockets support, so the communicator\n"
         "type argument to cs_proxy_comm_initialize() must be different\n"
         "from CS_PROXY_COMM_TYPE_SOCKET (%d)."),
       (int)CS_PROXY_COMM_TYPE_SOCKET);

#endif

  }

  if (retval == 0)
    bft_printf("[ok]\n");
  else {
    BFT_FREE(comm);
  }

  bft_printf_flush();

  /* End */

  return comm;
}

/*----------------------------------------------------------------------------
 * Finalize a communicator
 *----------------------------------------------------------------------------*/

static cs_proxy_comm_t *
_comm_finalize(cs_proxy_comm_t *comm)
{
  if (comm != NULL) {

    bft_printf("\n");

    /* Info on closing interface files */

    bft_printf(_("Closing communication: %s\n"),
               comm->port_name);

#if defined(HAVE_SOCKET)
    if (comm->socket > -1)
      _comm_sock_disconnect(comm);
#endif

    BFT_FREE(comm->port_name);

    BFT_FREE(comm);
  }

  return NULL;
}

/*============================================================================
 * Public function definitions
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
                         cs_proxy_comm_type_t  type)
{
  if (cs_glob_rank_id <= 0)
    _cs_glob_proxy_comm = _comm_initialize(port_name,
                                           key,
                                           type);
}

/*----------------------------------------------------------------------------
 * Finalize a connection to a proxy.
 *----------------------------------------------------------------------------*/

void
cs_proxy_comm_finalize(void)
{
  if (cs_glob_rank_id <= 0)
    _cs_glob_proxy_comm = _comm_finalize(_cs_glob_proxy_comm);
}

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
                    size_t       count)
{
  cs_proxy_comm_t *comm = _cs_glob_proxy_comm;

  assert(comm != NULL);

#if defined(HAVE_SOCKET)
  if (comm->socket > -1)
    _comm_write_sock(comm, rec, size, count);
#endif
}

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
                   size_t   count)
{
  cs_proxy_comm_t *comm = _cs_glob_proxy_comm;

  assert(comm != NULL);

#if defined(HAVE_SOCKET)
  if (comm->socket > -1)
    _comm_read_sock(comm, rec, size, count);
#endif
}

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
                            const char      *string_vals[])
{
  int i;
  char _base_header[256];
  char *_header = _base_header;

  size_t header_size = 0;
  size_t block_size = 256;

  cs_proxy_comm_t *comm = _cs_glob_proxy_comm;

  /* Compute request length, allocate if necessary */

  header_size = 32 + 4 + 4*n_ints + 8*n_doubles;

  for (i = 0; i < n_strings; i++)
    header_size += strlen(string_vals[i]) + 1;

  if (header_size > 256) {
    block_size = 256*((header_size/256) + 1);
    BFT_MALLOC(_header, block_size, char);
  }

  /* Form request header */

  memset(_header, 0, block_size);

  strncpy(_header, func_name, 32);

  /* Pack integer values */

#if (_CS_STDC_VERSION < 199901L)

  assert(sizeof(int) == 4);

  {
    int *_val_p;

    _val_p = (int *)(_header + 32);
    *_val_p = comp_id;
    _val_p = (int *)(_header + 36);
    *_val_p = header_size;
    _val_p = (int *)(_header + 40);
    *_val_p = n_ints;
    _val_p = (int *)(_header + 44);
    *_val_p = n_doubles;
    _val_p = (int *)(_header + 48);
    *_val_p = n_strings;

    header_size = 52; /* 32 + 5*4 */

    for (i = 0; i < n_ints; i++) {
      _val_p = (int *)(_header + header_size);
      *_val_p = int_vals[i];
      header_size += 4;
    }
  }

#else

  {
    int32_t *_val_p;

    _val_p = (int32_t *)(_header + 32);
    *_val_p = comp_id;
    _val_p = (int32_t *)(_header + 36);
    *_val_p = header_size;
    _val_p = (int32_t *)(_header + 40);
    *_val_p = n_ints;
    _val_p = (int32_t *)(_header + 44);
    *_val_p = n_doubles;
    _val_p = (int32_t *)(_header + 48);
    *_val_p = n_strings;

    header_size = 52; /* 32 + 5*4 */

    for (i = 0; i < n_ints; i++) {
      _val_p = (int32_t *)(_header + header_size);
      *_val_p = int_vals[i];
      header_size += 4;
    }
  }

#endif

  if (comm->swap_endian == true)
    bft_file_swap_endian(_header+32, _header+32, 4, (header_size - 32)/4);

  /* Pack double values (in C99, int32_t could avoid the sizeof() condition) */

  assert(sizeof(double) == 8);

  {
    int header_size_start = header_size;
    double *_val_p;

    for (i = 0; i < n_doubles; i++) {
      _val_p = (double *)(_header + header_size);
      *_val_p = double_vals[i];
      header_size += 8;
    }

    if (comm->swap_endian == true)
      bft_file_swap_endian(_header + header_size_start,
                           _header + header_size_start,
                           8,
                           (header_size - header_size_start)/8);

  }

  /* Pack strings */

  for (i = 0; i < n_strings; i++) {

    strcpy(_header + header_size, string_vals[i]);
    header_size += strlen(string_vals[i]);
    _header[header_size] = '\0';

  }

  /* Write message */

  cs_proxy_comm_write(_header, 1, block_size);

  /* Free allocated memory */

  if (_header != _base_header)
    BFT_FREE(_header);
}

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
                            char      *string_vals[])
{
  int i;
  char _base_header[256];
  char *_header = _base_header;

  size_t header_size = 0;
  size_t header_pos = 0;
  size_t block_size = 256;

  int retval = 0;

  cs_proxy_comm_t *comm = _cs_glob_proxy_comm;

  /* Read initial part of message */

  cs_proxy_comm_read(_header, 1, block_size);

  /* Unpack base info (return code, header_size, n_args/type) */

  if (comm->swap_endian == true)
    bft_file_swap_endian(_header, _header, 4, 5);

#if (_CS_STDC_VERSION < 199901L)

  assert(sizeof(int) == 4);

  {
    int *_val_p;

    _val_p = (int *)(_header);
    retval = *_val_p;
    _val_p = (int *)(_header + 4);
    header_size = *_val_p;
    _val_p = (int *)(_header + 8);
    n_ints = *_val_p;
    _val_p = (int *)(_header + 12);
    n_doubles = *_val_p;
    _val_p = (int *)(_header + 16);
    n_strings = *_val_p;
  }

#else

  {
    int32_t *_val_p;

    _val_p = (int32_t *)(_header);
    retval = *_val_p;
    _val_p = (int32_t *)(_header + 4);
    header_size = *_val_p;
    _val_p = (int32_t *)(_header + 8);
    n_ints = *_val_p;
    _val_p = (int32_t *)(_header + 12);
    n_doubles = *_val_p;
    _val_p = (int32_t *)(_header + 16);
    n_strings = *_val_p;
  }

#endif

  if (header_size > 256) {
    block_size = 256*((header_size/256) + 1);
    BFT_MALLOC(_header, block_size, char);
    memcpy(_header, _base_header, 256);
    cs_proxy_comm_read(_header+256, 1, block_size-256);
  }

  if (retval != 0) {
    if (_header != _base_header)
      BFT_FREE(_header);
    return retval;
  }

  header_pos = 20; /* 5*4 */

  /* Unpack integer values */

  if (comm->swap_endian == true)
    bft_file_swap_endian(_header + header_pos,
                         _header + header_pos,
                         4,
                         n_ints);

#if (_CS_STDC_VERSION < 199901L)

  assert(sizeof(int) == 4);

  {
    int *_val_p;

    for (i = 0; i < n_ints; i++) {
      _val_p = (int *)(_header + header_pos);
      int_vals[i] = *_val_p;
      header_pos += 4;
    }
  }

#else

  {
    int32_t *_val_p;

    for (i = 0; i < n_ints; i++) {
      _val_p = (int32_t *)(_header + header_pos);
      int_vals[i] = *_val_p;
      header_pos += 4;
    }
  }

#endif

  /* Unpack double values (in C99, int32_t could avoid the sizeof() condition) */

  assert(sizeof(double) == 8);

  if (comm->swap_endian == true)
    bft_file_swap_endian(_header + header_pos,
                         _header + header_pos,
                         8,
                         n_doubles);

  {
    double *_val_p;

    for (i = 0; i < n_doubles; i++) {
      _val_p = (double *)(_header + header_pos);
      double_vals[i] = *_val_p;
      header_pos += 8;
    }
  }

  /* Unpack strings */

  for (i = 0; i < n_strings; i++) {
    strcpy(string_vals[i], _header + header_pos);
    header_pos += strlen(string_vals[i]) + 1;
  }

  /* Free allocated memory */

  if (_header != _base_header)
    BFT_FREE(_header);

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
