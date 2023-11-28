/*============================================================================
 * Interactive control management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_UNISTD_H) && defined(HAVE_ACCESS)
#include <unistd.h>
#endif

#if defined(HAVE_SOCKET)
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_file.h"
#include "cs_log.h"
#include "cs_notebook.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_resource.h"
#include "cs_restart.h"
#include "cs_time_plot.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_control.h"

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_CONTROL_COMM_MAGIC_STRING       "CFD_control_comm_socket"

#define CS_CONTROL_COMM_FILE_NUM_LEN        4

#define CS_CONTROL_COMM_HOSTNAME_L        256
#define CS_CONTROL_COMM_NAME_L            256

#define CS_CONTROL_COMM_L_TYPE_NAME         2
#define CS_CONTROL_COMM_L_SEC_NUM           4

/* If SSIZE_MAX is not defined by the sytem headers, we take the minimum value
   required by POSIX (for low level reads/writes with sockets). */

#if !defined(SSIZE_MAX)
# define SSIZE_MAX  32767
#endif

#define CS_CONTROL_TIME_COMM  1

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_CONTROL_COMM_TYPE_SOCKET,    /* Communicate through sockets */
  CS_CONTROL_COMM_TYPE_NULL       /* Null communicator */

} cs_control_comm_type_t;

/* Communication handler structure */

typedef struct {

  char                   *port_name;         /* Port name (hostname:socket
                                                for IP sockets) */

#if defined(HAVE_SOCKET)
  FILE                   *recv_log;          /* If active, log of date read */
  FILE                   *trace;             /* If active, trace */

  int                     socket;            /* Socket number; */
#endif

  bool                    connected;         /* Currently connected ? */
  bool                    swap_endian;       /* Force big-endian communication */

  cs_control_comm_type_t  type;              /* Communicator type */

  bool                    errors_are_fatal;  /* If true, abort in case of error;
                                                otherwise, disconnect */

} cs_control_comm_t;

/* Saving of queued commands */

typedef struct {

  size_t                  buf_idx[4];    /* Buffer index (0: next command,
                                            1: partial read start; 2: end,
                                            3: size) */
  char                   *buf;           /* Buffer */

} cs_control_queue_t;

/* TODO: use control queue after advance for remaining commands
   * commands may return response in case of controller sockets
   * socket errors should give warnings, not errors to allow dirty disconnect */

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_control.c
 *
 *  \brief Handle control file usable for interactive change of stop,
 *         post-processing or checkpoint behavior.
 */

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_control_queue_t *_cs_glob_control_queue = NULL;
cs_control_comm_t *_cs_glob_control_comm = NULL;

static double  _control_file_wt_interval = 0.;
static double  _control_file_wt_last = -1.;

static int     _control_advance_steps = -1;
static int     _flush_nt = -1;

static int     _n_input_notebook_vars = 0;
static int     _n_output_notebook_vars = 0;

static int     *_input_notebook_vars = NULL;
static int     *_output_notebook_vars = NULL;

/* Timer statistics */

static cs_timer_counter_t   _control_t_tot = {.nsec = 0};
static cs_timer_counter_t   _control_comm_t_tot = {.nsec = 0};
static cs_timer_counter_t   _control_recv_t_tot = {.nsec = 0};
static cs_timer_counter_t   _control_send_t_tot = {.nsec = 0};

static long _n_reads = 0;
static cs_timer_counter_t   _control_read_t_min = {.nsec = 0};
static cs_timer_counter_t   _control_read_t_max = {.nsec = 0};
static cs_timer_counter_t   _control_read_t_sum = {.nsec = 0};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a control queue
 *
 * returns:
 *   pointer to initialized control queue
 *----------------------------------------------------------------------------*/

static cs_control_queue_t *
_queue_initialize(void)
{
  cs_control_queue_t *queue = NULL;

  BFT_MALLOC(queue, 1, cs_control_queue_t);

  queue->buf = NULL;

  queue->buf_idx[0] = 0;
  queue->buf_idx[1] = 0;
  queue->buf_idx[2] = 0;
  queue->buf_idx[3] = 0;

  return queue;
}

/*----------------------------------------------------------------------------
 * Finalize a queue
 *----------------------------------------------------------------------------*/

static void
_queue_finalize(cs_control_queue_t  **queue)
{
  if (queue != NULL) {
    if (*queue == NULL)
      return;
    cs_control_queue_t  *_queue = *queue;
    BFT_FREE(_queue->buf);
    BFT_FREE(*queue);
  }
}

#if defined(HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Write string for buffer characters, including non-printable ones
 *
 * parameters:
 *   f   <-- FILE
 *   rec <-- record to trace
 *   n   <-- number of elements
 *----------------------------------------------------------------------------*/

static void
_trace_buf(FILE        *f,
           const char   rec[],
           size_t       n)
{
  for (size_t j = 0; j < n; j++) {
    char c = rec[j];
    if (isprint(c))
      fprintf(f, "%c", c);
    else {
      switch(c) {
      case '\r':
        fprintf(f, "\\r");
        break;
      case '\n':
        fprintf(f, "\\n");
        break;
      case '\t':
        fprintf(f, "\\t");
        break;
      default:
        fprintf(f, "'\\%u'", (int)c);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *
 * parameters:
 *   dest <-- pointer to converted data location.
 *   src  --> pointer to source data location.
 *   size <-- size of each item of data in bytes.
 *   ni   <-- number of data items.
 *----------------------------------------------------------------------------*/

static void
_swap_endian(void        *dest,
             const void  *src,
             size_t       size,
             size_t       ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)dest;
  const unsigned char  *psrc = (const unsigned char *)src;

  for (i = 0; i < ni; i++) {

    shift = i * size;

    for (ib = 0; ib < (size / 2); ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }

  }

  if (dest != src && size == 1)
    memcpy(dest, src, ni);
}

/*----------------------------------------------------------------------------
 * Close an interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_sock_disconnect(cs_control_comm_t  *comm)
{
  if (close(comm->socket) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Communication %s:\n"
                "Error closing socket."),
              comm->port_name);

  comm->socket = -1;
  comm->connected = false;
}

/*----------------------------------------------------------------------------
 * Read a record from an interface socket into a given buffer
 *
 * parameters:
 *   comm    <-- pointer to communicator
 *   rec     <-- buffer for record, or NULL for default
 *   n_bytes <-- number of elements to read
 *----------------------------------------------------------------------------*/

static void
_comm_read_sock_r0(const cs_control_comm_t  *comm,
                   void                     *rec,
                   size_t                    n_bytes)
{
  char    *_rec = rec;

  assert(rec  != NULL);
  assert(comm != NULL);
  assert(comm->socket > -1);

  size_t start_id = 0;

  /* Read record from socket */

  while (start_id < n_bytes) {

    size_t end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);
    size_t n_loc = end_id - start_id;

    if (comm->trace != NULL) {
      fprintf(comm->trace, "-- reading up to %d values...\n",
              (int)n_bytes);
    }

    ssize_t ret = read(comm->socket, _rec + start_id, n_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error receiving data through socket."),
                comm->port_name);
    else if (comm->trace != NULL) {
      fprintf(comm->trace, "   read %d bytes: [", (int)ret);
      _trace_buf(comm->trace, _rec + start_id, ret);
      fprintf(comm->trace, "]\n");
      fflush(comm->trace);
    }

    start_id += ret;

  }
}

/*----------------------------------------------------------------------------
 * Read a record from an interface socket into a given buffer
 *
 * parameters:
 *   comm  <-- pointer to communicator
 *   queue <-- pointer to queue, or NULL
 *   rec   <-- buffer for record, or NULL for default
 *   size  <-- element size
 *   count <-- number of elements to read
 *----------------------------------------------------------------------------*/

static void
_comm_read_sock(const cs_control_comm_t  *comm,
                cs_control_queue_t       *queue,
                void                     *rec,
                size_t                    size,
                size_t                    count)
{
  char    *_rec = rec;

  assert(rec  != NULL);
  assert(comm != NULL);

  size_t n_bytes = size * count;

  size_t start_id = 0;

  /* Some data might already be buffered in queue */

  if (queue != NULL) {
    if (queue->buf_idx[2] > queue->buf_idx[1]) {

      start_id = queue->buf_idx[2] - queue->buf_idx[1];
      size_t n_remain = 0;
      if (start_id > n_bytes) {
        n_remain = start_id - n_bytes;
        start_id = n_bytes;
      }

      size_t s_id = queue->buf_idx[1];
      memcpy(rec, queue->buf + s_id, start_id);

      if (n_remain > 0)
        memmove(queue->buf + s_id,
                queue->buf + s_id + n_bytes,
                n_remain);
      queue->buf_idx[2] = queue->buf_idx[1] + n_remain;

    }
  }

  /* Read record from socket */

  if (cs_glob_rank_id <= 0) {

    assert(comm->socket > -1);

    while (start_id < n_bytes) {

      size_t end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);
      size_t n_loc = end_id - start_id;

      if (comm->trace != NULL) {
        fprintf(comm->trace, "-- reading up to %d bytes, "
                " %d of %d bytes already buffered...\n",
                (int)n_loc, (int)start_id, (int)n_bytes);
      }

      ssize_t ret = read(comm->socket, _rec + start_id, n_loc);

      if (ret < 1)
        bft_error(__FILE__, __LINE__, errno,
                  _("Communication %s:\n"
                    "Error receiving data through socket."),
                  comm->port_name);
      else if (comm->trace != NULL) {
        fprintf(comm->trace, "   read %d bytes: [", (int)ret);
        _trace_buf(comm->trace, _rec + start_id, ret);
        fprintf(comm->trace, "]\n");
        fflush(comm->trace);
      }

      start_id += ret;

    }

    if (comm->trace != NULL && size > 1) {
      for (size_t i = 0; i < count; i++) {
        fprintf(comm->trace, "    ");
        for (size_t j = 0; j < size; j++)
          fprintf(comm->trace, " %x", (unsigned)_rec[i*size + j]);
        fprintf(comm->trace, "\n");
      }
    }

    if (comm->swap_endian == true && size > 1)
      _swap_endian(rec, rec, size, count);
  }

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0) {
    MPI_Bcast(rec, n_bytes, MPI_BYTE, 0, cs_glob_mpi_comm);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Write a record to an interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_write_sock(cs_control_comm_t  *comm,
                 const void         *rec,
                 size_t              size,
                 size_t              count)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  ssize_t  ret;

  char        *_rec_swap = NULL;
  const char  *_rec = rec;

  assert(rec  != NULL);
  assert(comm != NULL);

  if (comm->socket < 0)
    return;

  cs_timer_t t0 = cs_timer_time();

  /* Determine associated size */

  size_t n_bytes = size * count;

  /* Convert if "little-endian" */

  if (comm->swap_endian == true && size != 1) {
    BFT_MALLOC(_rec_swap, n_bytes, char);
    _swap_endian(_rec_swap, rec, size, count);
    _rec = _rec_swap;
  }

  if (comm->trace != NULL) {
    if (size == 1) {
      fprintf(comm->trace, "-- write %d bytes: [",
              (int)count);
      _trace_buf(comm->trace, _rec, n_bytes);
      fprintf(comm->trace, "]...\n");
    }
    else {
      fprintf(comm->trace, "-- write %d values of size %d:\n",
              (int)count, (int)size);
      for (size_t i = 0; i < count; i++) {
        fprintf(comm->trace, "    ");
        for (size_t j = 0; j < size; j++)
          fprintf(comm->trace, " %x", (unsigned)(_rec[i*size + j]));
        fprintf(comm->trace, "\n");
      }
    }
    fflush(comm->trace);
  }

  /* Write record to socket */

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    ret = write(comm->socket, _rec + start_id, n_loc);

    if (ret < 1) {
      if (comm->errors_are_fatal)
        bft_error(__FILE__, __LINE__, errno,
                  _("Communication %s:\n"
                    "Error sending data through socket."),
                  comm->port_name);
      else {
        bft_printf(_("Communication %s:\n"
                     "Error sending data through socket."),
                   comm->port_name);
        _comm_sock_disconnect(comm);
      }
    }
    else if (comm->trace != NULL) {
      fprintf(comm->trace, "   wrote %d bytes\n", (int)ret);
      fflush(comm->trace);
    }

    start_id += ret;
  }

  if (_rec_swap != NULL)
    BFT_FREE(_rec_swap);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_control_send_t_tot, &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Prepare queue for reading of record into a given buffer
 *
 * parameters:
 *   queue <-- pointer to queue
 *   comm  <-> pointer to communicator
 *
 * return: size of data initially present
 *----------------------------------------------------------------------------*/

static size_t
_comm_read_sock_to_queue_prepare(cs_control_queue_t  *queue,
                                 cs_control_comm_t   *comm)
{
  if (queue->buf == NULL) {
    queue->buf_idx[0] = 0,
    queue->buf_idx[1] = 0,
    queue->buf_idx[2] = 0,
    queue->buf_idx[3] = 32767;
    BFT_MALLOC(queue->buf, queue->buf_idx[3]+1, char);
  }

  assert(comm != NULL);

  if (queue->buf_idx[0] < queue->buf_idx[1]) {
    bft_error(__FILE__, __LINE__, errno,
              "%s:\n"
              "  queue must be empty before reading additional data "
              "through socket.", __func__);
    return 0;
  }

  /* Move previously read but not consumed data to beginning of queue */

  ssize_t n_prv = queue->buf_idx[2] - queue->buf_idx[1];
  if (n_prv > 0) {
    memmove(queue->buf, queue->buf + queue->buf_idx[1], n_prv);
    queue->buf_idx[2] -= queue->buf_idx[1];
  }
  else
    queue->buf_idx[2] = 0;

  queue->buf_idx[1] = 0;
  queue->buf_idx[0] = 0;

  return (size_t)n_prv;
}

/*----------------------------------------------------------------------------
 * Read a record from an interface socket into a given buffer
 *
 * parameters:
 *   queue <-- pointer to queue
 *   comm  <-> pointer to communicator
 *----------------------------------------------------------------------------*/

static void
_comm_read_sock_to_queue(cs_control_queue_t  *queue,
                         cs_control_comm_t   *comm)
{
  /* Read record from socket */

  ssize_t  n_loc_max = queue->buf_idx[3];

  size_t start_id = queue->buf_idx[2] - queue->buf_idx[1];

  while (queue->buf_idx[1] < 1) {

    n_loc_max = queue->buf_idx[3] - start_id;

    if (comm->trace != NULL) {
      fprintf(comm->trace, "-- reading up to %d bytes...\n",
              (int)n_loc_max);
    }

#if CS_CONTROL_TIME_COMM != 0
    cs_timer_t t0 = cs_timer_time();
#endif

    ssize_t ret = read(comm->socket, queue->buf + start_id, n_loc_max);

#if CS_CONTROL_TIME_COMM != 0
    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_t td = cs_timer_diff(&t0, &t1);
#if CS_CONTROL_TIME_COMM > 1
    bft_printf("%s: read %d bytes in %g seconds\n",
               __func__, (int)ret, (double)(td.nsec)*1e-9);
#endif
    _control_read_t_sum.nsec += td.nsec;
    if (   _control_read_t_min.nsec == 0
        || _control_read_t_min.nsec > td.nsec)
      _control_read_t_min.nsec = td.nsec;
    if (_control_read_t_max.nsec < td.nsec)
      _control_read_t_max.nsec = td.nsec;
    _n_reads += 1;
#endif

    if (comm->trace != NULL) {
      fprintf(comm->trace, "   read %d bytes: [", (int)ret);
      _trace_buf(comm->trace, queue->buf + start_id, ret);
      fprintf(comm->trace, "]\n");
      fflush(comm->trace);
    }

    if (ret < 1 && start_id == 0) {
      if (comm->errors_are_fatal)
        bft_error(__FILE__, __LINE__, errno,
                  _("Communication %s:\n"
                    "Error receiving data through socket."),
                  comm->port_name);
      else {
        bft_printf(_("Communication %s:\n"
                     "Error receiving data through socket."),
                   comm->port_name);
        _comm_sock_disconnect(comm);
      }
    }

    queue->buf_idx[2] = start_id + ret;

    /* Check for last end of string */
    for (size_t i = queue->buf_idx[2]; i > start_id; i--) {
      if (queue->buf[i-1] == '\0') {
        queue->buf_idx[1] = i;
        break;
      }
    }

    start_id += ret;

    /* If we do not have a complete line, continue reading if possible */
    if (ret == n_loc_max) {
      if (queue->buf_idx[1] < 1) {
        queue->buf_idx[3] *= 2;
        BFT_REALLOC(queue->buf, queue->buf_idx[3], char);
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Synchronize queue across MPI ranks, and update connect status.
 *
 * parameters:
 *   start_id <-- previous start id in queue
 *   queue    <-- pointer to queue
 *   comm     <-> pointer to communicator
 *
 * return:
 *   number of elements usable prior to next read (i.e. buffered elements)
 *----------------------------------------------------------------------------*/

static size_t
_comm_read_sock_to_queue_sync(size_t               start_id,
                              cs_control_queue_t  *queue,
                              cs_control_comm_t   *comm)
{
#if defined(HAVE_MPI)

  if (cs_glob_rank_id >= 0) {

    long buf[4] = {queue->buf_idx[1],
                   queue->buf_idx[2],
                   queue->buf_idx[3],
                   0};

#if defined(HAVE_SOCKET)
    buf[3] = comm->socket;
#endif

    MPI_Bcast(buf, 4, MPI_LONG,  0, cs_glob_mpi_comm);

    if (buf[3] < 0)
      comm->connected = false;

    if (buf[2] > (long)queue->buf_idx[3]) {
      queue->buf_idx[3] = buf[2];
      BFT_REALLOC(queue->buf, queue->buf_idx[3], char);
    }

    queue->buf_idx[2] = buf[1];
    queue->buf_idx[1] = buf[0];

    if (queue->buf_idx[2] > start_id) {
      size_t count = queue->buf_idx[2] - start_id;
      MPI_Bcast(queue->buf + start_id, count, MPI_CHAR, 0,
                cs_glob_mpi_comm);
    }

  }

#endif /* defined(HAVE_MPI) */

  return queue->buf_idx[1];
}

/*----------------------------------------------------------------------------
 * Connection for socket initialization
 *----------------------------------------------------------------------------*/

static void
_comm_sock_connect(cs_control_comm_t  *comm)
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

  /* Establish communication with client */
  /*-------------------------------------*/

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
    _swap_endian((char *)&(sock_addr.sin_port),
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
 *
 * return:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

static int
_comm_sock_handshake(cs_control_comm_t  *comm,
                     const char         *magic_string,
                     const char         *key)
{
  int retval = 0;

  int   len = strlen(magic_string);
  char *str_cmp = NULL;

  /* Send key */

  _comm_write_sock(comm, key, 1, strlen(key));

  /* Write "magic string" */

  _comm_write_sock(comm, magic_string, 1, len);

  /* Read same magic string */

  BFT_MALLOC(str_cmp, len + 1, char);

  _comm_read_sock_r0(comm, str_cmp, len);
  str_cmp[len] = '\0';

  if (strncmp(str_cmp, magic_string, len)) {
    bft_error(__FILE__, __LINE__, 0, _("Handshake with client failed."));
    retval = 1;
  }

  BFT_FREE(str_cmp);

  return retval;
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

static cs_control_comm_t *
_comm_initialize(const char             *port_name,
                 const char             *key,
                 cs_control_comm_type_t  type)
{
  int retval = 0;

  cs_control_comm_t *comm = NULL;

  BFT_MALLOC(comm, 1, cs_control_comm_t);

  /* Initialize fields */

  if (cs_glob_rank_id <= 0) {
    BFT_MALLOC(comm->port_name, strlen(port_name) + 1, char);
    strcpy(comm->port_name, port_name);
  }
  else
    comm->port_name = NULL;

  comm->type = type;
  comm->errors_are_fatal = true;
  comm->connected = false;

  /* Test if system is big-endian */

  comm->swap_endian = false; /* Use "big-endian" mode to communicate */

  unsigned int_endian = 0;
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

#if defined(HAVE_SOCKET)
  comm->recv_log = NULL;
  comm->trace = NULL;

  comm->socket = -1;
#endif

  /* Info on interface creation */

  if (cs_glob_rank_id <= 0) {

    if (comm->port_name != NULL)
      bft_printf(_("Connecting to client:  %s ..."), comm->port_name);
    else
      bft_printf(_("Connecting to client ..."));
    bft_printf_flush();

    /* Initialize interface */

    if (type == CS_CONTROL_COMM_TYPE_SOCKET) {

#if defined(HAVE_SOCKET)

      const char *p;
      p = getenv("CS_CONTROL_RECV_LOG");
      if (p != NULL)
        comm->recv_log = fopen(p, "w");
      p = getenv("CS_CONTROL_COMM_TRACE");
      if (p != NULL)
        comm->trace = fopen(p, "w");

      _comm_sock_connect(comm);
      retval = _comm_sock_handshake(comm, CS_CONTROL_COMM_MAGIC_STRING, key);

#else

      bft_printf("\n");
      bft_error
        (__FILE__, __LINE__, 0,
         _("Library compiled without sockets support, so the communicator\n"
           "type argument to cs_control_comm_initialize() must be different\n"
           "from CS_CONTROL_COMM_TYPE_SOCKET (%d)."),
         (int)CS_CONTROL_COMM_TYPE_SOCKET);

#endif

    }

    if (retval == 0)
      bft_printf("[ok]\n");

    bft_printf_flush();

  }

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0) {
    MPI_Bcast(&retval, 1, MPI_INT, 0, cs_glob_mpi_comm);
  }
#endif

  if (retval != 0)
    BFT_FREE(comm);
  else
    comm->connected = true;

  /* End */

  return comm;
}

/*----------------------------------------------------------------------------
 * Finalize a communicator
 *----------------------------------------------------------------------------*/

static void
_comm_finalize(cs_control_comm_t  **comm)
{
  if (comm != NULL) {

    if (*comm == NULL)
      return;

    cs_control_comm_t  *_comm = *comm;

    if (cs_glob_rank_id <= 0) {
      bft_printf(_("\nClosing communication: %s\n"),
                 _comm->port_name);
    }

#if defined(HAVE_SOCKET)
    if (_comm->socket > -1) {
      _comm_sock_disconnect(_comm);
      if (_comm->recv_log != NULL)
        fclose(_comm->recv_log);
      if (_comm->trace != NULL)
        fclose(_comm->trace);
    }
#endif

    BFT_FREE(_comm->port_name);

    BFT_FREE(*comm);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Establish a connection to a client.
 *
 * \param[in]  port_name  name of server port (host:port for IP sockets)
 * \param[in]  key        key for authentification
 * \param[in]  type       communication type
 */
/*----------------------------------------------------------------------------*/

static void
_control_comm_initialize(const char              *port_name,
                         const char              *key,
                         cs_control_comm_type_t   type)
{
  _cs_glob_control_comm = _comm_initialize(port_name, key, type);

  _control_advance_steps = 1;

  if (_cs_glob_control_queue == NULL)
    _cs_glob_control_queue = _queue_initialize();

  /* Call controller again right away */
  cs_control_check_file();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read data from a client into a command queue
 *
 * The function updates a pointer (view) to the data.
 *
 * \return number of useable elements read
 *         (i.e. text elements before next read is required)
 */
/*----------------------------------------------------------------------------*/

static size_t
_control_comm_read_to_queue(void)
{
  size_t retval = 0;
  cs_control_queue_t *queue = _cs_glob_control_queue;
  cs_control_comm_t *comm = _cs_glob_control_comm;

  cs_timer_t t0 = cs_timer_time();

  /* If no communicator, simply update queue for possible
     remaining operations */

  if (comm == NULL) {
    return retval;
  }
  else if (queue != NULL) {
    size_t n_prv = _comm_read_sock_to_queue_prepare(queue, comm);

#if defined(HAVE_SOCKET)
    if (comm->socket > -1)
      _comm_read_sock_to_queue(queue, comm);
#endif

    retval = _comm_read_sock_to_queue_sync(n_prv, queue, comm);
    if (comm->connected == false)
      _comm_finalize(&comm);
    _cs_glob_control_comm = comm;
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_control_recv_t_tot, &t0, &t1);

  return retval;
}

/*----------------------------------------------------------------------------
 * Read next value, expecting integer
 *
 * parameters:
 *   cur_line <-- pointer to line buffer
 *   s        <-> current string position
 *   val      --> integer read
 *
 * returns:
 *   number of integers read (1 for success, 0 otherwise)
 *----------------------------------------------------------------------------*/

static int
_read_next_int(const char   *cur_line,
               const char  **s,
               int          *val)
{
  int n_val = 0;

  const char *p = *s;
  while (*p != '\0' && *p != ' ' && *p != '\t')
    p++;
  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  n_val = sscanf(*s, "%i", val);

  if (n_val == 0)
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

  return n_val;
}

/*----------------------------------------------------------------------------
 * Read next optional value, expecting integer
 *
 * parameters:
 *   s        <-> current string position
 *   val     --> integer read
 *
 * returns:
 *   number of integers read (1 for success, 0 otherwise)
 *----------------------------------------------------------------------------*/

static int
_read_next_opt_int(const char  **s,
                   int          *val)
{
  int n_val = 0;

  const char *p = *s;
  while (*p != '\0' && *p != ' ' && *p != '\t')
    p++;
  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  n_val = sscanf(*s, "%i", val);

  return n_val;
}

/*----------------------------------------------------------------------------
 * Read next value, expecting double precision floating point value
 *
 * parameters:
 *   skip_prev <-- skip previous string ?
 *   cur_line  <-- pointer to line buffer
 *   s         <-> current string position
 *   val       --> value read
 *
 * returns:
 *   number of values read (1 for success, 0 otherwise)
 *----------------------------------------------------------------------------*/

static int
_read_next_double(bool          skip_prev,
                  const char   *cur_line,
                  const char  **s,
                  double       *val)
{
  int n_val = 0;

  const char *p = *s;

  if (skip_prev) {
    while (*p != '\0' && *p != ' ' && *p != '\t')
      p++;
  }

  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  n_val = sscanf(*s, "%lg", val);

  if (n_val == 0)
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

  return n_val;
}

/*----------------------------------------------------------------------------
 * Read next value, expecting string
 *
 * parameters:
 *   skip_prev <-- skip previous string ?
 *   s         <-> current string position
 *   s_val     --> pointer to string read
 *----------------------------------------------------------------------------*/

static void
_read_next_string(bool    skip_prev,
                  char  **s,
                  char  **s_val)
{
  *s_val = NULL;

  char *p = *s;
  if (skip_prev) {
    while (*p != '\0' && *p != ' ' && *p != '\t')
      p++;
  }
  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  *s_val = *s;

  while (*p != '\0' && *p != ' ' && *p != '\t')
    p++;
  if (*p != '\0') {
    *p = '\0';
    p++;
  }
  *s = p;
}

/*----------------------------------------------------------------------------
 * Handle command file line relative to checkpointing
 *
 * parameters:
 *   cur_line <-- pointer to line buffer
 *   s        <-> pointer to current position in line
 *----------------------------------------------------------------------------*/

static void
_control_checkpoint(const char   *cur_line,
                    const char  **s)
{
  *s += 11; /* shift in string by length of "checkpoint_" part */

  if (strncmp(*s, "time_step ", 10) == 0) {
    int nt;
    if (_read_next_int(cur_line, s, &nt) > 0) {
      cs_restart_checkpoint_set_next_ts(nt);
      bft_printf("  %-32s %12d\n",
                 "checkpoint_time_step", nt);
    }
  }
  else if (strncmp(*s, "time_value ", 11) == 0) {
    double t;
    if (_read_next_double(true, cur_line, s, &t) > 0) {
      cs_restart_checkpoint_set_next_tv(t);
      bft_printf("  %-32s %12.5g\n",
                 "checkpoint_time_value", t);
    }
  }
  else if (strncmp(*s, "wall_time ", 10) == 0) {
    double wt;
    if (_read_next_double(true, cur_line, s, &wt) > 0) {
      cs_restart_checkpoint_set_next_wt(wt);
      bft_printf("  %-32s %12.5g\n",
                 "checkpoint_wall_time", wt);
    }
  }
  else if (strncmp(*s, "time_step_interval ", 19) == 0) {
    int nt;
    if (_read_next_int(cur_line, s, &nt) > 0) {
      cs_restart_checkpoint_set_interval(nt, -1., -1.);
      bft_printf("  %-32s %12d\n",
                 "checkpoint_time_step_interval", nt);
    }
  }
  else if (strncmp(*s, "time_value_interval ", 20) == 0) {
    double t;
    if (_read_next_double(true, cur_line, s, &t) > 0) {
      if (t > 0) {
        cs_restart_checkpoint_set_interval(-1, t, -1.);
        bft_printf("  %-32s %12.5g\n",
                   "checkpoint_time_value_interval", t);
      }
      else
        bft_printf("  %-32s %12.5g %s\n",
                   "checkpoint_time_value_interval", t, _("ignored"));
    }
  }
  else if (strncmp(*s, "wall_time_interval ", 19) == 0) {
    double wt;
    if (_read_next_double(true, cur_line, s, &wt) > 0) {
      if (wt > 0) {
        cs_restart_checkpoint_set_interval(-1, -1., wt);
        bft_printf("  %-32s %12.5g\n",
                   "checkpoint_wall_time_interval", wt);
      }
      else
         bft_printf("  %-32s %12.5g %s\n",
                   "checkpoint_wall_time_interval", wt, _("ignored"));
    }
  }
  else
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);
}

/*----------------------------------------------------------------------------
 * Handle command file line relative to notebook
 *
 * parameters:
 *   ts       <-- pointer to time step status
 *   cur_line <-> pointer to the current line
 *   s        <-> pointer to current position in line
 *----------------------------------------------------------------------------*/

static void
_control_notebook(const cs_time_step_t   *ts,
                  char                   *cur_line,
                  char                  **s)
{
  CS_UNUSED(ts);

  *s += 9; /* shift in string by length of "notebook_" part */

  bool ignored = true;

  /* Set specifically at this current time */

  if (strncmp(*s, "set ", 4) == 0) {
    char *name;
    double val = 0.;
    _read_next_string(true, s, &name);
    if (_read_next_double(false, cur_line, (const char **)s, &val) == 1) {
      int editable;
      cs_notebook_parameter_is_present(name, &editable);
      if (editable) {
        cs_notebook_parameter_set_value(name, val);
        bft_printf("  %-32s \"%s\" set to %12.5g\n",
                   "notebook", name, val);
      }
      ignored = false;
    }
  }

  /* Get specifically at this current time */

  else if (strncmp(*s, "get ", 4) == 0) {
    char *name;
    double val = 0.;
    _read_next_string(true, s, &name);
    int id = cs_notebook_parameter_get_id(name);
    if (id > -1) {
      val = cs_notebook_parameter_value_by_name(name);
      bft_printf("  %-32s \"%s\" get : %12.5g\n",
                 "notebook", name, val);
      ignored = false;

#if defined(HAVE_SOCKET)
      if (_cs_glob_control_comm != NULL) {
        char reply[20] = "\0";
        sprintf(reply, "get: %.3f", val);
        _comm_write_sock(_cs_glob_control_comm, reply, 1, strlen(reply) + 1);
      }
#endif
    }
  }

  /* Define for systematic exchange */

  else if (strncmp(*s, "add_input ", 10) == 0) {
    char *name;
    double val = 0.;
    _read_next_string(true, s, &name);

    int id = cs_notebook_parameter_get_id(name);
    int editable = 0;
    cs_notebook_parameter_is_present(name, &editable);
    if (! editable)
      id = -1;

    _n_input_notebook_vars += 1;
    BFT_REALLOC(_input_notebook_vars, _n_input_notebook_vars, int);
    _input_notebook_vars[_n_input_notebook_vars - 1] = id;

    if (id < 0)
      bft_printf("  %-32s \"%s\" input does not match "
                 "an editable notebook variable\n",
                 "notebook", name);

    if (_read_next_double(false, cur_line, (const char **)s, &val) == 1) {
      if (id > -1) {
        cs_notebook_parameter_set_value(name, val);
        bft_printf("  %-32s \"%s\" set to %12.5g\n",
                   "notebook", name, val);
      }
    }
  }

  else if (strncmp(*s, "add_output ", 10) == 0) {
    char *name;
    double val = 0.;
    _read_next_string(true, s, &name);

    int id = cs_notebook_parameter_get_id(name);

    _n_output_notebook_vars += 1;
    BFT_REALLOC(_output_notebook_vars, _n_output_notebook_vars, int);
    _output_notebook_vars[_n_output_notebook_vars - 1] = id;

    if (id > -1) {
      val = cs_notebook_parameter_value_by_name(name);
      bft_printf("  %-32s \"%s\" get : %12.5g\n",
                 "notebook", name, val);
      ignored = false;

#if defined(HAVE_SOCKET)
      if (_cs_glob_control_comm != NULL) {
        char reply[20] = "\0";
        sprintf(reply, "get: %.3f", val);
        _comm_write_sock(_cs_glob_control_comm, reply, 1, strlen(reply) + 1);
      }
#endif
    }

  }

  if (ignored)
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);
}

/*----------------------------------------------------------------------------
 * Handle command file line relative to postprocessing
 *
 * parameters:
 *   ts       <-- pointer to time step status
 *   s        <-> pointer to current position in line
 *----------------------------------------------------------------------------*/

static void
_control_postprocess(const cs_time_step_t   *ts,
                     char                   *cur_line,
                     const char            **s)
{
  *s += 12; /* shift in string by length of "postprocess_" part */

  if (strncmp(*s, "time_step ", 10) == 0) {
    int nt = 0, writer_id = 0;
    if (_read_next_int(cur_line, s, &nt) > 0) {
      if (_read_next_opt_int(s, &writer_id) == 0)
        writer_id = 0;
      if (nt >= 0)
        nt = CS_MAX(nt, ts->nt_cur);
      else
        nt = CS_MAX(nt, -ts->nt_cur);
      cs_post_add_writer_t_step(writer_id, nt);
      bft_printf("  %-32s %12d %12d\n",
                 "postprocess_time_step", nt, writer_id);
    }
  }
  else if (strncmp(*s, "time_value ", 11) == 0) {
    int writer_id = 0;
    double t = 0.;
    if (_read_next_double(true, cur_line, s, &t) > 0) {
      if (_read_next_opt_int(s, &writer_id) == 0)
        writer_id = 0;
      if (t >= 0)
        t = CS_MAX(t, ts->t_cur);
      else
        t = CS_MAX(t, -ts->t_cur);
      bft_printf("  %-32s %12.5g %12d\n",
                 "postprocess_time_value", t, writer_id);
      cs_post_add_writer_t_value(writer_id, t);
    }
  }
  else
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

}

/*----------------------------------------------------------------------------
 * Parse control file or queue
 *
 * For some commands (such as advance), this will return before reading
 * the full buffer, so the read size is returned so as to allow resuming
 * later.
 *
 * parameters:
 *   name         <-- name of control system
 *   buffer       <-- pointer to file contents buffer
 *   f_size       <-- size of buffer for file and connector
 *   control_comm <-- control communicator, or NULL
 *   queue        <-- pointer to queue, or NULL
 *
 * returns:
 *   0 if buffer completely handled, index of unhandled part otherwise
 *----------------------------------------------------------------------------*/

static size_t
_parse_control_buffer(const char          *name,
                      char                *buffer,
                      long                 f_size,
                      cs_control_comm_t   *control_comm,
                      cs_control_queue_t  *queue)
{
  size_t retval = 0;

  int nt_max;

  char *s;
  char *cur_line = buffer, *next_line = NULL;
  const cs_time_step_t  *ts = cs_glob_time_step;

  if (name != NULL)
    bft_printf
      (_("\n"
         " Options set or changed by \"%s\":\n"
         " -------------------------\n\n"), name);

  /* Loop on buffer's lines */

  /* Note that when parsing lines, we do not use strtok() type functions, to
     avoid modifying a line (so that log/error/warning messages may simply use
     that line); hence also tests using strncp on strings appended with a
     whitespace (always present in case of arguments) rather than strcmp. */

  while (cur_line != NULL) {

    int retcode = 0;

    /* Prepare current and next line for parsing */

    next_line = cur_line;

    while (   *next_line != '\0'
           && *next_line != '\n' && *next_line != '\r' && *next_line != '\f')
      next_line++;

    *next_line = '\0'; /* marks end of cur_line */
    next_line += 1;

    if (next_line >= (buffer + f_size))
      next_line = NULL;
    else {
      while (    *next_line != '\0'
             && (*next_line == '\n' || *next_line == '\r' || *next_line == '\f'))
        next_line++;
    }

    /* Check for keywords given in control_file and store the related values */

    size_t l = strlen(cur_line);

    if (control_comm != NULL) {
      if (control_comm->recv_log != NULL)
        fprintf(control_comm->recv_log, "%s\n", cur_line);
    }

    for (size_t i = 0; i < l; i++) {
      if (cur_line[i] == '#') {
        cur_line[i] = '\0';
        break;
      }
    }

    for (s = cur_line; *s == ' ' || *s == '\t'; s++)
      *s = ' ';

    if (*s == '\0') {
      cur_line = next_line;
      continue;
    }

    /* Calculation end options
       default with no keyword; max_time_step */

    nt_max = -1;
    if (sscanf(s, "%i", &nt_max) > 0)
      nt_max = CS_MAX(nt_max, 0);
    else if (strncmp(s, "max_time_step ", 14) == 0) {
      if (_read_next_int(cur_line, (const char **)&s, &nt_max) > 0)
        nt_max = CS_MAX(nt_max, 0);
    }
    else if (strncmp(s, "time_step_limit ", 16) == 0) {
      if (_read_next_int(cur_line, (const char **)&s, &nt_max) > 0) {
        nt_max = CS_MAX(nt_max, 0);
        if (ts->nt_max > -1)
          nt_max = CS_MIN(nt_max + ts->nt_prev, ts->nt_max);
        else
          nt_max = nt_max + ts->nt_prev;
      }
    }

    if (nt_max > -1) {
      nt_max = CS_MAX(nt_max, ts->nt_cur);
      cs_time_step_define_nt_max(nt_max);
      bft_printf("  %-32s %12d (%s %d)\n",
                 "max_time_step", ts->nt_max, _("current:"), ts->nt_cur);
    }
    else if (strncmp(s, "max_time_value ", 15) == 0) {
      double t_max;
      if (_read_next_double(true, cur_line, (const char **)&s, &t_max) > 0)
        t_max = CS_MAX(t_max, ts->t_cur);
      cs_time_step_define_t_max(t_max);
      bft_printf("  %-32s %12.5g (%s %12.5g)\n",
                 "max_time_value", ts->t_max, _("current:"), ts->t_cur);
    }
    else if (strncmp(s, "max_wall_time ", 14) == 0) {
      double wt_max;
      if (_read_next_double(true, cur_line, (const char **)&s, &wt_max) > 0)
        wt_max = CS_MAX(wt_max, cs_timer_wtime());
      bft_printf("  %-32s %12.5g (%s %12.5g)\n",
                 "max_wall_time", wt_max, _("current:"),
                 cs_resource_get_wt_limit());
      cs_resource_set_wt_limit(wt_max);
    }

    /* Control file check interval */

    else if (strncmp(s, "control_file_wtime_interval ", 28) == 0) {
      double wt;
      if (_read_next_double(true, cur_line, (const char **)&s, &wt) > 0)
        _control_file_wt_interval = wt;
    }

    /* Checkpointing options */

    else if (strncmp(s, "checkpoint_", 11) == 0)
      _control_checkpoint(cur_line, (const char **)&s);

    /* Postprocessing options */

    else if (strncmp(s, "postprocess_", 12) == 0)
      _control_postprocess(ts, cur_line, (const char **)&s);

    /* Notebook options */

    else if (strncmp(s, "notebook_", 9) == 0)
      _control_notebook(ts, cur_line, &s);

    /* Force flush of logs */

    else if (strncmp(s, "flush", 5) == 0) {
      int nt = -1;
      if (_read_next_int(cur_line, (const char **)&s, &nt) > 0) {
        if (nt > -1)
          _flush_nt = CS_MAX(nt, ts->nt_cur);
        else
          _flush_nt = -1;
      }
      else
        _flush_nt = ts->nt_cur;
      bft_printf(_("  flush logs and time plots at time step %12d\n"),
                 _flush_nt);
    }

    /* Connect/disconnect request */

    else if (strncmp(s, "connect ", 8) == 0) {
      char *port_name, *key;
      _read_next_string(true, &s, &port_name);
      _read_next_string(false, &s, &key);
      _control_comm_initialize(port_name, key,
                               CS_CONTROL_COMM_TYPE_SOCKET);
    }

    else if (strncmp(s, "disconnect ", 11) == 0) {
      _comm_finalize(&_cs_glob_control_comm);
      break;
    }

    /* Advance time */

    else if (strncmp(s, "advance ", 8) == 0) {

      int n = -1;
      if (_read_next_opt_int((const char **)&s, &n) == 0)
        n = 1;
      if (_control_advance_steps <= 0)
        _control_advance_steps = n;
      else
        _control_advance_steps += n;
      cur_line = next_line;

#if defined(HAVE_SOCKET)
      if (control_comm != NULL) {
        char reply[5] = "0\0";
        _comm_write_sock(control_comm, reply, 1, 2);
      }
#endif

      /* Use binary exchanges for multiple variables here,
         to minimize number of messages */

      int n_notebook_vals = CS_MAX(_n_input_notebook_vars,
                                   _n_output_notebook_vars);

      if (n_notebook_vals > 0) {
        double *notebook_vals = NULL;

        BFT_MALLOC(notebook_vals, n_notebook_vals, double);

        /* Send output variables from previous time step,
           get input variables for new time step. */

        cs_notebook_get_values(_n_output_notebook_vars,
                               _output_notebook_vars,
                               notebook_vals);

#if defined(HAVE_SOCKET)
        if (control_comm != NULL) {
          if (_n_output_notebook_vars > 0)
            _comm_write_sock(control_comm, notebook_vals,
                             sizeof(double), _n_output_notebook_vars);
          if (_n_input_notebook_vars > 0)
            _comm_read_sock(control_comm, queue, notebook_vals,
                            sizeof(double), _n_input_notebook_vars);
        }
#endif

        /* Get input variables */

        cs_notebook_set_values(_n_input_notebook_vars,
                               _input_notebook_vars,
                               notebook_vals);

        BFT_FREE(notebook_vals);
      }

      break;

    }

    /* Unhandled lines */

    else {
      bft_printf(_("   ignored control command: \"%s\"\n"), cur_line);
      retcode = -1;
    }

    /* Prepare for next line */

    cur_line = next_line;

#if defined(HAVE_SOCKET)
    if (control_comm != NULL) {
      char reply[5];
      snprintf(reply, 5, "%d", retcode);
      reply[4] = '\0';
      _comm_write_sock(control_comm, reply, 1, strlen(reply) + 1);
    }
#endif

  } /* End of loop on lines */

  if (name != NULL) {

    /* Empty control_file equivalent to flush request */

    if (f_size < 1) {
      _flush_nt = ts->nt_cur;
      bft_printf(_("  flush logs and time plots at time step %12d\n"),
                 _flush_nt);
    }

    bft_printf
      (_("\n"
         " Finished reading \"%s\".\n\n"), name);
  }

  if (next_line != NULL) {
    ptrdiff_t i = next_line - buffer;
    retval = i;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize controller structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_control_finalize(void)
{
  _n_input_notebook_vars = 0;
  _n_output_notebook_vars = 0;

  BFT_FREE(_input_notebook_vars);
  BFT_FREE(_output_notebook_vars);

  cs_timer_t t0 = cs_timer_time();

  _comm_finalize(&_cs_glob_control_comm);
  _queue_finalize(&_cs_glob_control_queue);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_control_t_tot, &t0, &t1);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Total elapsed time for controller:  %.3f s\n"),
                _control_t_tot.nsec*1e-9);

  if (_control_comm_t_tot.nsec > 0)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Useful elapsed time for controller:  %.3f s\n"),
                  _control_comm_t_tot.nsec*1e-9);

  if (_control_recv_t_tot.nsec > 0) {
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Communication times for controller:\n"
                    "                            send:    %.3f s\n"
                    "                            receive: %.3f s\n"),
                  _control_send_t_tot.nsec*1e-9,
                  _control_recv_t_tot.nsec*1e-9);
    if (_n_reads > 0)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("\n"
                      "Communication times detail:\n"
                      "                            reads:   %lu\n"
                      "                            read time min: %.3e s\n"
                      "                                      max: %.3e s\n"
                      "                                      mean: %.3e s\n"),
                    (unsigned long)_n_reads,
                    _control_read_t_min.nsec*1e-9,
                    _control_read_t_max.nsec*1e-9,
                    _control_read_t_sum.nsec/_n_reads*1e-9);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the presence of a control file and deal with the interactive
 *        control.
 */
/*----------------------------------------------------------------------------*/

void
cs_control_check_file(void)
{
  bool time_comm = (_cs_glob_control_comm != NULL) ? true : false;

  cs_timer_t t0 = cs_timer_time();

  long f_size = -1;
  char *buffer = NULL;
  const cs_time_step_t  *ts = cs_glob_time_step;

  const char path[] = "control_file";

  /* Test existence and size of file */

  if (cs_glob_rank_id <= 0) {

    if (   _control_file_wt_interval <= 0.
        ||(    cs_timer_wtime() - _control_file_wt_last
           >= _control_file_wt_interval)) {

#if defined(HAVE_UNISTD_H) && defined(HAVE_ACCESS)

      /* Test existence of file using access() before stat(),
         as this may be less costly on some filesytems
         (such as on LUSTRE, due to metadata handling aspects). */

      if (access(path, F_OK) == 0)
        f_size = cs_file_size(path);

#else

      f_size = cs_file_size(path);

#endif

    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0)
    MPI_Bcast(&f_size,  1, MPI_LONG,  0, cs_glob_mpi_comm);
#endif

  /* If file exists, handle it */

  if (f_size >= 0) {

    BFT_MALLOC(buffer, f_size + 1, char);

    if (cs_glob_rank_id <= 0) {

      size_t r_size = 0;
      FILE *control_file = fopen("control_file", "r");

      if (control_file != NULL) {
        r_size = fread(buffer, 1, f_size, control_file);
        buffer[r_size] = '\0'; /* precaution for partial read */
        fclose(control_file);
        remove("control_file");
      }
      else
        bft_printf
          (_("\n"
             " Warning: error opening %s (ignored):\n"
             " --------\n"
             "   \"%s\"\n\n"), path, strerror(errno));

      _control_file_wt_last = cs_timer_wtime();

    }

#if defined(HAVE_MPI)
    if (cs_glob_rank_id >= 0)
      MPI_Bcast(buffer, f_size + 1, MPI_CHAR, 0, cs_glob_mpi_comm);
#endif

    /* Now all ranks have required buffer */

    _parse_control_buffer("control_file", buffer, f_size, NULL, NULL);

    BFT_FREE(buffer);
  }

  /* Test control queue and connection second */

  if (_control_advance_steps > 0) {

    _control_advance_steps -= 1;

#if defined(HAVE_SOCKET)
    char reply[13] = "Iteration OK";
    if (_cs_glob_control_comm != NULL)
      _comm_write_sock(_cs_glob_control_comm, reply, 1, strlen(reply)+1);
#endif
  }

  if (_cs_glob_control_queue != NULL) {

    while (_control_advance_steps < 1) {

      /* If the queue is empty, receive new commands */

      if (   _cs_glob_control_queue->buf_idx[0]
          >= _cs_glob_control_queue->buf_idx[1]) {

        size_t n = _control_comm_read_to_queue();

        if (n == 0 && _cs_glob_control_comm == NULL) {
          _queue_finalize(&_cs_glob_control_queue);
          break;
        }

      }

      /* Handle queued elements */

      while (  _cs_glob_control_queue->buf_idx[0]
             < _cs_glob_control_queue->buf_idx[1]) {

        size_t idx = _parse_control_buffer(NULL,
                                           _cs_glob_control_queue->buf,
                                           _cs_glob_control_queue->buf_idx[1],
                                           _cs_glob_control_comm,
                                           _cs_glob_control_queue);

        if (idx == 0)
          _cs_glob_control_queue->buf_idx[0]
            = _cs_glob_control_queue->buf_idx[1];
        else
          _cs_glob_control_queue->buf_idx[0] = idx;

      }

    }

  }

  if (_flush_nt == ts->nt_cur) {
    _flush_nt = -1;
    cs_log_printf_flush(CS_LOG_N_TYPES);
    bft_printf_flush();
    cs_time_plot_flush_all();
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_control_t_tot, &t0, &t1);

  if (time_comm)
    cs_timer_counter_add_diff(&_control_comm_t_tot, &t0, &t1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
