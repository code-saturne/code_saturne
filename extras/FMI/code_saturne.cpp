/*============================================================================
 * Functional Mock-up Unit main methods implementation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.
  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C++/C library headers
 *----------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string.h>
#include <array>
#include <limits.h>
#include <sys/time.h>
#include <sys/time.h>
#include <pthread.h>

#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "code_saturne.h"

/*----------------------------------------------------------------------------*/

using namespace std;

#define CS_DRY_RUN 0   // Set to 1 to for simulation mode with no actual
                       // connection to code_saturne.

#define CS_TIMING 0    // Log some timings.

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct sockaddr_in sockaddr_in;
typedef struct sockaddr sockaddr;

struct thread_data {
  sockaddr_in address;
  int         addrlen;
};

/* Log handling types */
/*--------------------*/

/* As the generated code does not provide for detecting whether logging is on,
   and the generator does not provide for defining allowed categories,
   we use an additional layer here, with mapping to default categories,
   and a lower-level log mechanism for debugging. */

typedef enum {

  CS_LOG_EVENTS,   /* Events */
  CS_LOG_WARNING,  /* Warnings */
  CS_LOG_ERROR,    /* Error */
  CS_LOG_FATAL,    /* Fatal error */
  CS_LOG_ALL,      /* All messages */
  CS_LOG_COMM,     /* Log communication (low-level) */
  CS_LOG_LAUNCH,   /* Log system calls */
  CS_LOG_TRACE     /* Other trace logs (low-level) */

} cs_log_types_t;

/* Saving of queued reads */

typedef struct {

  size_t                  buf_idx[3];    /* Buffer index
                                            (0: partial read start;
                                            1: end,
                                            2: size) */
  char                   *buf;           /* Buffer */

} cs_control_queue_t;

/* Save value references of read and written variables */

typedef struct {

  int  n_input;
  int  n_output;

  int  input_max;
  int  output_max;

  int  input_max_size;
  int  output_max_size;

  int     *input_ids;
  int     *output_ids;

  double  *input_vals;
  double  *output_vals;

} cs_variables_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static int _n_iter = 0;
static int _master_socket = -1;
static int _cs_socket = -1;
static int _cs_swap_endian = 0;

static fmi2ComponentEnvironment  _component_environment = nullptr;
static fmi2String                _instance_name = "[unnamed]";

static FILE *tracefile = NULL;

/* Mapping to default log categories
   categories up to "logAll" are standardized, the rest are local. */

static const char *_log_categories[] = {"logEvents",
                                        "logStatusWarning",
                                        "logStatusError",
                                        "logStatusFatal",
                                        "logAll",
                                        "logComm",
                                        "logLaunch",
                                        "logTrace"};

static const char *_log_prefix[] = {"[event]   ",
                                    "[warning] ",
                                    "[error]   ",
                                    "[fatal]   ",
                                    "[]        ",
                                    "[comm]    ",
                                    "[launch]  ",
                                    "[trace]   "};

/* Active logging: 0: inactive, 1, log using FMI, -1: low-level trace */

static int _log_active[] = {-1,  /* events */
                            -1,  /* warning */
                            -1,  /* error */
                            -1,  /* fatal */
                            -1,  /* all */
                             0,  /* comm */
                            -1,  /* launch */
                            -1}; /* trace */

/* Read_queue */

cs_control_queue_t *_cs_glob_control_queue = NULL;

/* Structure for grouped notebook value settings */

cs_variables_t *_cs_variables = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log messages. The generated code does not provide for querying of
   logging activation, so we add a layer at least allowing easier mapping. */

void _cs_log(fmi2Status       status,
             cs_log_types_t   log_type,
             const char      *text)
{
  if (_log_active[log_type] > 0)
    log(_component_environment, _instance_name,
        status, _log_categories[log_type], text, nullptr);
  else if (_log_active[log_type] < 0)
    cout << _log_prefix[log_type] << _instance_name << ": " << text << endl;
}

void _cs_log(fmi2Status       status,
             cs_log_types_t   log_type,
             string           text)
{
  if (_log_active[log_type] > 0)
    log(_component_environment, _instance_name,
        status, _log_categories[log_type], text.c_str(), nullptr);
  else if (_log_active[log_type] < 0)
    cout << _log_prefix[log_type] << _instance_name << ": " << text << endl;
}

void _cs_log(fmi2ComponentEnvironment  environment,
             fmi2String                id,
             fmi2Status                status,
             cs_log_types_t            log_type,
             fmi2String                text,
             void*                     pointers)
{
  if (_log_active[log_type] > 0)
    log(environment, id, status, _log_categories[log_type], text, pointers);
  else if (_log_active[log_type] < 0)
    cout << _log_prefix[log_type] << id << ": " << text << endl;
}

void _cs_log(fmi2ComponentEnvironment  environment,
             fmi2String                id,
             fmi2Status                status,
             cs_log_types_t            log_type,
             string                    text,
             void*                     pointers)
{
  if (_log_active[log_type] > 0)
    log(environment, id, status, _log_categories[log_type], text.c_str(),
        pointers);
  else if (_log_active[log_type] < 0)
    cout << _log_prefix[log_type] << id << ": " << text << endl;
}

/*----------------------------------------------------------------------------*/

/* Generate a random integer between min and max */

static int
_randrange(int   min,
           int   max)
{
  return min + rand() / (RAND_MAX / (max - min + 1) + 1);
}

/*----------------------------------------------------------------------------*/

/* Generate a random key which will be later exchanged with code_saturne */

static int
_generate_key(void)
{
  time_t t;
  srand(static_cast<unsigned int>(time(&t)));

  int key = _randrange(0, static_cast<int>(2e31));

  return key;
}

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
 * Initialize a control queue
 *
 * returns:
 *   pointer to initialized control queue
 *----------------------------------------------------------------------------*/

static cs_control_queue_t *
_queue_initialize(void)
{
  cs_control_queue_t *queue = NULL;

  queue = (cs_control_queue_t *)malloc(sizeof(cs_control_queue_t));

  queue->buf = NULL;

  queue->buf_idx[0] = 0;
  queue->buf_idx[1] = 0;
  queue->buf_idx[2] = 32767;
  queue->buf = (char *)malloc(queue->buf_idx[2]+1);

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
    free(_queue->buf);
    free(*queue);
  }
}

/*----------------------------------------------------------------------------
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * parameters:
 *   buf  <-> pointer to data location.
 *   size <-- size of each item of data in bytes.
 *   ni   <-- number of data items.
 *----------------------------------------------------------------------------*/

static void
_swap_endian(void        *buf,
             size_t       size,
             size_t       ni)
{
  unsigned char  *pbuf = (unsigned char *)buf;

  for (size_t i = 0; i < ni; i++) {

    size_t shift = i * size;

    for (size_t ib = 0; ib < (size / 2); ib++) {

      unsigned char tmpswap = *(pbuf + shift + ib);
      *(pbuf + shift + ib) = *(pbuf + shift + (size - 1) - ib);
      *(pbuf + shift + (size - 1) - ib) = tmpswap;

    }

  }
}

/*----------------------------------------------------------------------------*/

/* Send a data on a socket */

static void
_send_sock(int       sock,
           char     *buffer,
           size_t    size,
           size_t    ni)
{
  size_t start_id = 0;
  size_t n_bytes = size*ni;

#if CS_DRY_RUN == 1
  return;
#endif

  if (_cs_swap_endian && size > 1)
    _swap_endian(buffer, size, ni);

  if (tracefile != NULL) {
    if (size == 1) {
      fprintf(tracefile, "== send %d bytes: [", (int)n_bytes);
      _trace_buf(tracefile, buffer, n_bytes);
      fprintf(tracefile, "]...\n");
    }
    else {
      fprintf(tracefile, "== send %d values of size %d:\n", (int)ni, (int)size);
      for (size_t i = 0; i < ni; i++) {
        fprintf(tracefile, "    ");
        for (size_t j = 0; j < size; j++)
          fprintf(tracefile, " %x", (unsigned)buffer[i*size + j]);
        fprintf(tracefile, "\n");
      }
    }
    fflush(tracefile);
  }

  while (start_id < n_bytes) {

    size_t end_id = start_id + SSIZE_MAX;
    if (n_bytes < end_id)
      end_id = n_bytes;
    size_t n_loc = end_id - start_id;

    ssize_t ret = send(sock, buffer+start_id, n_loc, 0);

    if (ret < 1) {
      string s0 = "Error sending buffer ";
      string s2 = strerror(errno);
      _cs_log(fmi2Fatal, CS_LOG_FATAL, s0 + s2);
      exit(EXIT_FAILURE);
    }
    else if (tracefile != NULL) {
      fprintf(tracefile, "   sent %d bytes\n", (int)ret);
      fflush(tracefile);
    }

    start_id += ret;
  }

  if (_cs_swap_endian && size > 1)  /* restore endiannes */
    _swap_endian(buffer, size, ni);
}

/*----------------------------------------------------------------------------*/

/* Send a message on a socket */

static void
_send_sock_str(int          sock,
               const char  *str)
{
#if CS_DRY_RUN == 1
  return;
#endif

  size_t n = strlen(str)+1;

  if (tracefile != NULL) {
    fprintf(tracefile, "== send %d bytes: [", (int)n);
    _trace_buf(tracefile, str, n);
    fprintf(tracefile, "]...\n");
    fflush(tracefile);
  }

  ssize_t ret = send(sock, str, n, 0);

  if (ret < 1) {
    string s0 = "Error sending ";
    string s1 = str;
    string s2 = strerror(errno);
    _cs_log(fmi2Fatal, CS_LOG_FATAL, s0 + s1 + s2);
    exit(EXIT_FAILURE);
  }
  else if (tracefile != NULL) {
    fprintf(tracefile, "   sent %d bytes\n", (int)ret);
    fflush(tracefile);
  }
}

/*----------------------------------------------------------------------------*/

/* Receive data message (fixed size) on a socket */

static void
_recv_sock(int                  socket,
           char                *buffer,
           cs_control_queue_t  *queue,
           size_t               size,
           size_t               ni)
{
  size_t start_id = 0;
  size_t n_bytes = size*ni;

  /* Cleaning buffer */
  memset(buffer, 0, n_bytes);

#if CS_DRY_RUN == 1
  buffer[0] = '\0';
  return;
#endif

  /* Part of the message may already have been read to queue */

  if (queue != NULL) {
    ssize_t n_prv = queue->buf_idx[1] - queue->buf_idx[0];
    if (n_prv > 0) {
      if ((size_t)n_prv > n_bytes)
        n_prv = n_bytes;
      memcpy(buffer, queue->buf + queue->buf_idx[0], n_prv);
      queue->buf_idx[0] += n_prv;
      start_id = n_prv;
    }
  }

  while (start_id < n_bytes) {

    size_t end_id = start_id + SSIZE_MAX;
    if (n_bytes < end_id)
      end_id = n_bytes;
    size_t n_loc = end_id - start_id;

    if (tracefile != NULL) {
      fprintf(tracefile, "== receiving up to %d bytes, %d of %d bytes already buffered...\n",
              (int)n_loc, (int)start_id, (int)n_bytes);
    }

    ssize_t ret = recv(socket, buffer + start_id, n_loc, 0);

    if (ret < 1) {
      string s0 = "Error receiving data: ";
      string s1;
      if (errno != 0)
        s1 = strerror(errno);
      else
        s1 = "code_saturne may have disconnected.";
      _cs_log(fmi2Fatal, CS_LOG_FATAL, s0 + s1);
      exit(EXIT_FAILURE);
    }
    else if (tracefile != NULL) {
      fprintf(tracefile, "   received %d bytes: [", (int)ret);
      if (size == 1)
        _trace_buf(tracefile, buffer, ret);
      fprintf(tracefile, "]\n");
      fflush(tracefile);
    }

    start_id += ret;

  }

  if (tracefile != NULL && size > 1) {
    for (size_t i = 0; i < ni; i++) {
      fprintf(tracefile, "    ");
      for (size_t j = 0; j < size; j++)
        fprintf(tracefile, " %x", (unsigned)buffer[i*size + j]);
      fprintf(tracefile, "\n");
    }
  }

  if (_cs_swap_endian && size > 1)
    _swap_endian(buffer, size, ni);
}

/*----------------------------------------------------------------------------*/

/* Receive message (text expected) on a socket */

static char *
_recv_sock_with_queue(int                  socket,
                      cs_control_queue_t  *queue,
                      size_t               min_size)
{
  size_t start_id = queue->buf_idx[1] - queue->buf_idx[0];
  ssize_t cut_id = -1;

  if (tracefile != NULL) {
    fprintf(tracefile, "_recv_sock_with_queue: %d %d\n", (int)start_id, (int)min_size);
  }

  /* Move previously read but not consumed data to beginning of queue */

  if (start_id > 0) {
    memmove(queue->buf, queue->buf + queue->buf_idx[0], start_id);

    queue->buf_idx[1] -= queue->buf_idx[0];
    queue->buf_idx[0] = 0;

    for (size_t i = 0; i < start_id; i++) {
      if (queue->buf[i] == '\0') {
        cut_id = i;
        break;
      }
    }
  }

  /* Read records from socket until complete string is read */

  ssize_t n_loc_max = queue->buf_idx[2];

  while (cut_id < 0) {

    n_loc_max = queue->buf_idx[2] - start_id;

    if (n_loc_max < 1) {
      /* If we do not have a complete line, continue reading if possible */
      queue->buf_idx[2] *= 2;
      queue->buf = (char *)realloc(queue->buf, queue->buf_idx[2]);
      n_loc_max = queue->buf_idx[2] - start_id;
    }

    if (tracefile != NULL) {
      fprintf(tracefile, "== receiving up to %d bytes...\n",
              (int)n_loc_max);
    }

#if CS_DRY_RUN == 1
    queue->buf[start_id] =$ '\0';
    return;
#endif

    ssize_t ret = (socket != NULL) ?
      read(socket, queue->buf + start_id, n_loc_max) : 0;

    if (ret < 1) {
      string s0 = "Error receiving data: ";
      string s1;
      if (errno != 0)
        s1 = strerror(errno);
      else
        s1 = "code_saturne may have disconnected.";
      _cs_log(fmi2Fatal, CS_LOG_FATAL, s0 + s1);
      exit(EXIT_FAILURE);
    }

    if (tracefile != NULL) {
      fprintf(tracefile, "   received %d bytes: [", (int)ret);
      _trace_buf(tracefile, queue->buf + start_id, ret);
      fprintf(tracefile, "]\n");
      fflush(tracefile);
    }

    queue->buf_idx[1] = start_id + ret;

    /* Check for end of string */

    for (ssize_t i = 0; i < ret; i++) {
      if (queue->buf[start_id + i] == '\0') {
        cut_id = start_id + i;
        break;
      }
    }

  }

  /* A full string has now been read, and data from following
     reads may already be present also */

  queue->buf_idx[0] = cut_id + 1;

  /* Clean end of line */

  for (ssize_t i = cut_id; i > 0; i--) {
    char c = queue->buf[i];
    if (c == ' '  || c == '\t' || c == '\n' || c == '\r' || c == '\f') {
      queue->buf[i] = '\0';
      cut_id = i;
    }
    else
      break;
  }

  /* Set pointers to remaining part for next call */

  if (queue->buf_idx[0] > queue->buf_idx[1]) {
    queue->buf_idx[0] = 0;
    queue->buf_idx[1] = 0;
  }

  /* Return pointer to usable string in buffer
     (always at the beginning of the buffer). */

  return queue->buf;
}

/*----------------------------------------------------------------------------*/

/* Communication with code_saturne handshake at the beginning
 * (first connection) and reception of messages from CS after */

static void
_comm_with_saturne(int   key)
{
  char *key_buffer = nullptr;
  char *magic_buffer = nullptr;
  string magic_string = "CFD_control_comm_socket";
  string key_s = to_string(key);
  size_t len_key = key_s.size();

  /* Allocation of buffers */
  key_buffer = new char[len_key+1]();
  magic_buffer = new char[magic_string.size()+1]();

  /* code_saturne socket */
  _recv_sock(_cs_socket, key_buffer, NULL, 1, len_key);

  if (strncmp(key_s.c_str(), key_buffer, len_key) != 0) {
    _cs_log(fmi2Fatal, CS_LOG_ERROR,
            "wrong key received (socket handshake)");
    exit(EXIT_FAILURE);
  }

  _recv_sock(_cs_socket, magic_buffer, NULL, 1, magic_string.size());

  _send_sock_str(_cs_socket, magic_buffer);

  _cs_log(fmi2OK, CS_LOG_ALL,
          "Connection between FMU client and code_saturne established");

  /* Iteration OK */
  char buf_rcv[13];
  _recv_sock(_cs_socket, buf_rcv, NULL,1, 13);

  delete[] key_buffer;
  delete[] magic_buffer;

  _cs_glob_control_queue = _queue_initialize();
}

/*----------------------------------------------------------------------------*/

/* Configuration of the master socket (server) */

static int
_configure_server(sockaddr_in  *address)
{
  int master_socket;
  int opt = 1;
  socklen_t optlen = sizeof(opt);

  /* Create a master socket */
  if ((master_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
    _cs_log(fmi2Fatal, CS_LOG_FATAL, "socket() failed");
    exit(EXIT_FAILURE);
  }

  /* Set master socket to allow multiple connections */
  if (setsockopt(master_socket, SOL_SOCKET, SO_REUSEADDR, &opt, optlen) < 0) {
    _cs_log(fmi2Fatal, CS_LOG_FATAL, "setsockopt() failed");
    exit(EXIT_FAILURE);
  }

  /* Master socket configuration */
  address->sin_family = AF_INET;
  address->sin_addr.s_addr = htonl(INADDR_ANY);
  address->sin_port = 0;

  /* Bind the socket to any available localhost port*/
  if (bind(master_socket,
           (sockaddr *)address,
           sizeof(*address)) < 0) {
    _cs_log(fmi2Fatal, CS_LOG_FATAL, "bind() failed");
    exit(EXIT_FAILURE);
  }

  return master_socket;
}

/*----------------------------------------------------------------------------*/

/* Writing of control_files for code_saturne which is deleted after reading */

static void
_write_control_file(string path,
                    string hostname,
                    int    port,
                    int    key)
{
  /* Write control_file for CS */
  ofstream control_file(path + "/control_file");

  if (control_file.is_open()) {
    control_file << "connect " << hostname << ":" << port << " " << key;
    control_file.close();
  }
  else {
    _cs_log(fmi2Fatal, CS_LOG_FATAL, "Error writing control_file.");
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*/

static int
_init_server(string        path,
             sockaddr_in  *address)
{
  /* Master socket */
  int port;
  sockaddr_in foo;
  socklen_t foosize = sizeof(foo);
  string hostname = "127.0.0.1";

  /* Generation of the key */
  int key = _generate_key();

  /* Server configuration */
  _master_socket = _configure_server(address);

  /* Get the actual port */
  getsockname(_master_socket, (sockaddr *)&foo, &foosize);
  port = ntohs(foo.sin_port);

  string s0 = "Listener on port ";
  string s1 = to_string(port);
  _cs_log(fmi2OK, CS_LOG_ALL, s0 + s1);

  _write_control_file(path, hostname, port, key);

  return key;
}

/*----------------------------------------------------------------------------*/

static void *
_start_server(void *data)
{
  struct thread_data *t_data;
  t_data = (struct thread_data *)data;
  sockaddr_in address = t_data->address;
  int addrlen = t_data->addrlen;

  /* Try to specify maximum of 1 pending connection *
     for the master socket                          */
  if (listen(_master_socket, 1) < 0) {
    _cs_log(fmi2Fatal, CS_LOG_FATAL, "listen() failed.");
    exit(EXIT_FAILURE);
  }

  /* Test if system is big-endian */
  _cs_swap_endian = 0;
  unsigned int_endian = 0;
  *((char *) (&int_endian)) = '\1';
  if (int_endian == 1)
    _cs_swap_endian = 1;

  /* Create socket */

  if ((_cs_socket = accept(_master_socket,
                           (sockaddr *)&address,
                           (socklen_t*)&addrlen)) < 0) {
    _cs_log(fmi2Fatal, CS_LOG_FATAL, "accept() failed.");
    exit(EXIT_FAILURE);
  }

  /* Inform user of socket number - used in send and receive commands */

  string s = "New connection, socket fd is " + to_string(_cs_socket);
  _cs_log(fmi2OK, CS_LOG_COMM, s);

  s = " ip is: " + string(inet_ntoa(address.sin_addr));
  _cs_log(fmi2OK, CS_LOG_COMM, s);

  s = " port: " + to_string(ntohs(address.sin_port));;
  _cs_log(fmi2OK, CS_LOG_COMM, s);

  pthread_exit(nullptr);
}

/*----------------------------------------------------------------------------*/

static void
*_start_cs(void *resu)
{
  string &res = *(static_cast<string*>(resu));
  string cmd;

  /* Wait for the server to be listening */
  sleep(1);

  /* Run Code_saturne client in background */
  cmd = "cd " + res + " && ./run_solver &";
  system(cmd.c_str());

  pthread_exit(nullptr);
}

/*----------------------------------------------------------------------------
 * Initialize a variables grouping structure
 *
 * returns:
 *   pointer to initialized control queue
 *----------------------------------------------------------------------------*/

static cs_variables_t *
_variables_initialize(void)
{
  cs_variables_t *v = NULL;

  v = (cs_variables_t *)malloc(sizeof(cs_variables_t));

  v->n_input = 0;
  v->n_output = 0;

  v->input_max = 0;
  v->output_max = 0;

  v->input_max_size = 0;
  v->output_max_size = 0;

  v->input_ids = NULL;
  v->output_ids = NULL;

  v->input_vals = NULL;
  v->output_vals = NULL;

  return v;
}

/*----------------------------------------------------------------------------
 * Finalize a variables grouping structure
 *----------------------------------------------------------------------------*/

static void
_variables_finalize(cs_variables_t  **v)
{
  if (v != NULL) {
    if (*v == NULL)
      return;
    cs_variables_t  *_v = *v;

    if (_v->input_max_size > 0) {
      free(_v->input_ids);
      free(_v->input_vals);
    }
    if (_v->output_max_size > 0) {
      free(_v->output_ids);
      free(_v->output_vals);
    }

    free(*v);
  }
}

/*----------------------------------------------------------------------------
 * Add or set a variables input;
 *
 * return 0 if variable was already present, 1 if inserted
 *----------------------------------------------------------------------------*/

static int
_variables_add_input(cs_variables_t  *v,
                     int              id,
                     double           val)
{
  int retval = 0;

  assert(id >= 0);

  if (id >= v->input_max) {

    if (id >= v->input_max_size) {

      if (v->input_max_size == 0)
        v->input_max_size = 2;
      while (id >= v->input_max_size)
        v->input_max_size *= 2;

      v->input_ids = (int *)realloc(v->input_ids,
                                    v->input_max_size*sizeof(int));
      v->input_vals = (double *)realloc(v->input_vals,
                                        v->input_max_size*sizeof(double));

      for (int i = v->input_max; i < v->input_max_size; i++) {
        v->input_ids[i] = -1;
        v->input_vals[i] = 0;
      }
    }

    v->input_max = id + 1;
  }

  if (v->input_ids[id] < 0) {
    retval = 1;
    v->input_ids[id] = v->n_input;
    v->n_input += 1;
  }

  v->input_vals[id] = val;

  return retval;
}

/*----------------------------------------------------------------------------
 * Add or get a variables output;
 *
 * return 0 if variable was already present, 1 if inserted
 *----------------------------------------------------------------------------*/

static int
_variables_add_output(cs_variables_t  *v,
                      int              id)
{
  int retval = 0;

  assert(id >= 0);

  if (id >= v->output_max) {

    if (id >= v->output_max_size) {

      if (v->output_max_size == 0)
        v->output_max_size = 2;
      while (id >= v->output_max_size)
        v->output_max_size *= 2;

      v->output_ids = (int *)realloc(v->output_ids,
                                     v->output_max_size*sizeof(int));
      v->output_vals = (double *)realloc(v->output_vals,
                                         v->output_max_size*sizeof(double));

      for (int i = v->output_max; i < v->output_max_size; i++) {
        v->output_ids[i] = -1;
        v->output_vals[i] = 0;
      }
    }

    v->output_max = id + 1;
  }

  if (v->output_ids[id] < 0) {
    retval = 1;
    v->output_ids[id] = v->n_output;
    v->n_output += 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/

static void
_advance(int   n)
{
#if CS_TIMING == 1
  struct timeval  tv_time_0, tv_time_1;
  (void)gettimeofday(&tv_time_0, NULL);
#endif

  cs_variables_t *v = _cs_variables;

  string buffer = "advance " + to_string(n);

  string s = "Advancing " + to_string(n) + " iterations";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

#if CS_DRY_RUN == 1
  {
    for (int i = 0; i < v->output_max; i++) {
      static long _counter = 0;
      _counter += 1;
      int id = v->output_ids[i];
      if (id < 0)
        continue;
      v->output_vals[id] = (double)(_counter % 50) * 0.02;
    }
  }

  return;
#endif

  _send_sock_str(_cs_socket, buffer.c_str());

  /* receive 0 */
  char *buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);
  if (strcmp(buf_rcv, "0") != 0) {
    s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
    _cs_log(fmi2Warning, CS_LOG_WARNING, s);
  }

  /* Get output notebook values */

  if (v->n_output > 0) {
    string s = "Get " + to_string(v->n_output) + " output variables";
    _cs_log(fmi2OK, CS_LOG_TRACE, s);
    _recv_sock(_cs_socket, (char *)_cs_variables->output_vals,
               _cs_glob_control_queue,
               sizeof(double), v->n_output);
  }

  /* Send input notebook values */

  if (v->n_input > 0) {
    string s = "Send " + to_string(v->n_input) + " input variables";
    _send_sock(_cs_socket, (char *)_cs_variables->input_vals,
               sizeof(double), v->n_input);
  }

  /* receive iteration OK */
  if (n > 0) {

    buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);
    if (strcmp(buf_rcv, "Iteration OK") != 0) {
      s =   string(__func__) + ": expected \"Iteration OK\", not: "
          + string(buf_rcv);
      _cs_log(fmi2Warning, CS_LOG_WARNING, s);
    }

#if CS_TIMING == 1
    (void)gettimeofday(&tv_time_1, NULL);

    long usec =   (tv_time_1.tv_sec - tv_time_0.tv_sec) * (long)1000000
                + (tv_time_1.tv_usec - tv_time_0.tv_usec);
    double sec = (double)usec / 1000000.;

    cout << "Time to advance " << _instance_name << ": " << sec << endl;
#endif

  }
}

/*----------------------------------------------------------------------------*/

/* Send 'notebook_add_input value' to the server */

static void
_set_notebook_variable(int          sock,
                       const char  *variable,
                       double       value)
{
  string var_s(variable);

  string buffer = "notebook_add_input " + var_s + " " + to_string(value);

  string s = string(__func__) + ": sending " + string(variable);
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  _send_sock_str(sock, buffer.c_str());

#if CS_DRY_RUN == 1
  return;
#endif

  s = string(__func__) + ": waiting for reply";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  // Receive 0
  char *buf_rcv = _recv_sock_with_queue(sock, _cs_glob_control_queue, 0);

  if (strcmp(buf_rcv, "0") != 0) {
    s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
    _cs_log(fmi2Warning, CS_LOG_WARNING, s);
  }
  else {
    s = string(__func__) + ": variable set.";
    _cs_log(fmi2OK, CS_LOG_TRACE, s);
  }
}

/*----------------------------------------------------------------------------*/

/* Send 'notebook_get variable' to the server */

static double
_get_notebook_variable(int         sock,
                       const char *variable)
{
  char *eptr = nullptr;
  double val = 0.;
  string var_s(variable);

  string buffer = "notebook_add_output " + var_s;

  string s_log = string(__func__) + ": send query for " + string(variable);
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

#if CS_DRY_RUN == 1

  static long _counter = 0;
  _counter += 1;
  val = (double)(_counter % 50) * 0.02;

#else

  _send_sock_str(sock, buffer.c_str());

  s_log = string(__func__) + ": waiting for " + string(variable);
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  /* Received get: val */
  char *buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);

  if (strncmp(buf_rcv, "get: ", 5) != 0) {
    s_log =   string(__func__) + ": unexpected reply; " + string(buf_rcv);
    _cs_log(fmi2Error, CS_LOG_ERROR, s_log);
  }

  val = strtod(buf_rcv + 5, &eptr);

  s_log =   string(__func__) + ": retrieved "
          + string(variable) + " (" + to_string(val) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  /* Received 0 : OK */
  buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);
  if (strcmp(buf_rcv, "0") != 0) {
    s_log =   string(__func__) + ": unexpected return code " + string(buf_rcv);
    _cs_log(fmi2Error, CS_LOG_ERROR, s_log);
  }

#endif

  return val;
}

/*----------------------------------------------------------------------------*/

/* Send 'disconnect ' to the server */

static void
_disconnect(void)
{
  char buffer[13] = "disconnect ";

  _send_sock_str(_cs_socket, buffer);

  _cs_log(fmi2OK, CS_LOG_TRACE, "Disconnecting the controller...");
}

/*----------------------------------------------------------------------------*/

/* Equivalent of "system" but is able to get stdout */

string
exec_popen(const char* cmd)
{
  array<char, 128> buffer;
  string result;
  unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);

  if (!pipe) {
    throw runtime_error("popen() failed!");
  }

  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

/*----------------------------------------------------------------------------*/

void code_saturne::setResourceLocation(const char* fmiResourceLocation)
{
  string s0(__func__);
  string s_rl(fmiResourceLocation);

  /* As the generated code provides no other entry under
     fmi2Instantiate than this function, also set some other global
     pointers here (not the cleanest solution we would dream of, but
     a workaround the current limited scope of user-generated and
     auto-generated code boundaries). */

  _instance_name = (fmi2String)malloc(sizeof(id) + 1);
  strcpy(_instance_name, id);

  /* Also set id to saved string, as it seems the string to which
     id was assigned may have a temporary lifetime
     work around FMU generator bug */

  id = _instance_name;

  _cs_log(fmi2OK, CS_LOG_TRACE, s0 + s_rl);
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::init()
{
  /* As the code generator used provides no other earlier entry point,
     global pointer to component environment (for logging) here. */

  _component_environment = component;

  string s_casename(casename);
  string s_run_id(run_id);
  string s_code_saturne(code_saturne);
  string resu = s_casename + "/RESU/" + s_run_id;
  string cmd;
  int addrlen;
  sockaddr_in address;
  pthread_t thread_server, thread_cs;
  pthread_attr_t attr;
  void *status;

  _cs_log(fmi2OK, CS_LOG_TRACE, __func__);

  /* Just in case, to avoid comma separators */
  setlocale(LC_NUMERIC, "C");

  /* Set tracefile if requested here */

  const char *p = getenv("CS_FMU_COMM_TRACE");
  if (p != NULL) {
    if (strcmp(p, "stdout") != 0 && strlen(p) > 0)
      tracefile = fopen(p, "w");
    else
      tracefile = stdout;
  }

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* Initialize code_saturne calculation */
  cmd = "cd " + s_casename + " && " + s_code_saturne +
        " run --initialize --force --id='" + s_run_id + "'";

#if CS_DRY_RUN == 1

  _cs_log(fmi2OK, CS_LOG_LAUNCH, cmd);
  _cs_socket = 444444;

#else

  _cs_log(fmi2OK, CS_LOG_LAUNCH, exec_popen(cmd.c_str()).c_str());

  int key = _init_server(resu, &address);

  /* Accept the incoming connection */
  addrlen = sizeof(address);

  _cs_log(component, id, fmi2OK, CS_LOG_ALL, "Waiting for connection...",
          nullptr);

  struct thread_data data;
  data.address = address;
  data.addrlen = addrlen;

  pthread_create(&thread_server, nullptr, _start_server, static_cast<void*>(&data));
  pthread_create(&thread_cs, nullptr, _start_cs, static_cast<void*>(&resu));

  pthread_attr_destroy(&attr);
  pthread_join(thread_server, &status);
  pthread_join(thread_cs, &status);

  _comm_with_saturne(key);

#endif

  /* Sleep 3 seconds just to make sure the connection is well established */
  sleep(3);

  /* Now exchange input and output values which are already declared */

  if (_cs_variables == NULL)
    _cs_variables = _variables_initialize();

  for (int i = 0; i < _cs_variables->input_max; i++) {
    if (_cs_variables->input_ids[i] > -1)
      _set_notebook_variable(_cs_socket,
                             realVariables[i].name,
                             realVariables[i].getter(this));
  }
  for (int i = 0; i < _cs_variables->output_max; i++) {
    if (_cs_variables->output_ids[i] > -1) {
      double value
        = _get_notebook_variable(_cs_socket,
                                 realVariables[i].name);
      realVariables[i].setter(this, value);
    }
  }

  /* Increase the number of maximum time steps, as FMI is controlling
     the time stepping, not the code. */

  {
    _send_sock_str(_cs_socket, "max_time_step 9999999");

    char *buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);

    string s = string(__func__);
    if (strcmp(buf_rcv, "0") != 0) {
      s += ": unexpected return code: " + string(buf_rcv);
      _cs_log(fmi2Warning, CS_LOG_WARNING, s);
    }
    else {
      s += ": set max time step to 9999999 for FMI control.";
      _cs_log(fmi2OK, CS_LOG_TRACE, s);
    }
  }

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::doStep(double step)
{
  string s = "----- " + string(__func__) + "(" + to_string(step) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  cs_variables_t *v = _cs_variables;

  /* Update input notebook values */

  for (int i = 0; i < v->input_max; i++) {
    int id = v->input_ids[i];
    if (id < 0)
      continue;
    double val = realVariables[i].getter(this);
    v->input_vals[id] = val;
  }

  /* Advance 1 iteration */
  _advance(1);

  /* Update output notebook values */

  for (int i = 0; i < v->output_max; i++) {
    int id = v->output_ids[i];
    if (id < 0)
      continue;
    double val = v->output_vals[id];
    realVariables[i].setter(this, val);
  }

  _n_iter++;

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::terminate()
{
  string s_casename(casename);
  string s_run_id(run_id);
  string resu = s_casename + "/RESU/" + s_run_id;

  _send_sock_str(_cs_socket, "max_time_step 0");
  _advance(0);

  _variables_finalize(&_cs_variables);

  _queue_finalize(&_cs_glob_control_queue);

  /* Send disconnect to the server */
  _disconnect();

#if CS_DRY_RUN == 1

  _cs_socket = -1;

#else

  shutdown(_cs_socket, SHUT_RDWR);
  close(_cs_socket);

  /* Stop the computation manually */
  ofstream control_file(resu + "/control_file");
  if (control_file.is_open()) {
    control_file << 1;
    control_file.close();
  }
  else {
    _cs_log(fmi2Error, CS_LOG_ERROR, "Error writing control_file.");
  }

#endif

  if (tracefile != NULL) {
    if (tracefile != stdout && tracefile != stderr) {
      fclose(tracefile);
      tracefile = NULL;
    }
  }

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::reset()
{
  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

double code_saturne::getReal(int fmi2ValueReference)
{
  if (_cs_variables == NULL)
    _cs_variables = _variables_initialize();

  string s = string(__func__) + "(" + to_string(fmi2ValueReference) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  if (realVariables[fmi2ValueReference].causality() != input) {
    int first = _variables_add_output(_cs_variables, fmi2ValueReference);
    if (first && _cs_socket >= 0) {
      double value
        = _get_notebook_variable(_cs_socket,
                                 realVariables[fmi2ValueReference].name);
      realVariables[fmi2ValueReference].setter(this, value);
    }
  }

  return realVariables[fmi2ValueReference].getter(this);
}

/*----------------------------------------------------------------------------*/

void code_saturne::setReal(int fmi2ValueReference, double value)
{
  if (_cs_variables == NULL)
    _cs_variables = _variables_initialize();

  string s = string(__func__) + "(" + to_string(fmi2ValueReference)
    + ", " + to_string(value) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  if (realVariables[fmi2ValueReference].causality() == input) {
    realVariables[fmi2ValueReference].setter(this, value);

    int first = _variables_add_input(_cs_variables,
                                     fmi2ValueReference,
                                     value);
    if (first && _cs_socket >= 0) {
      _set_notebook_variable(_cs_socket,
                             realVariables[fmi2ValueReference].name,
                             value);
    }
  }
}

/*----------------------------------------------------------------------------*/

int code_saturne::getInteger(int fmi2ValueReference)
{
  return integerVariables[fmi2ValueReference].getter(this);
}

/*----------------------------------------------------------------------------*/

void code_saturne::setInteger(int fmi2ValueReference, int value)
{
  integerVariables[fmi2ValueReference].setter(this, value);
}

/*----------------------------------------------------------------------------*/

bool code_saturne::getBoolean(int fmi2ValueReference)
{
  return booleanVariables[fmi2ValueReference].getter(this);
}

/*----------------------------------------------------------------------------*/

void code_saturne::setBoolean(int fmi2ValueReference, bool value)
{
  booleanVariables[fmi2ValueReference].setter(this, value);
}

/*----------------------------------------------------------------------------*/

const char* code_saturne::getString(int fmi2ValueReference)
{
  return stringVariables[fmi2ValueReference].getter(this);
}

/*----------------------------------------------------------------------------*/

void code_saturne::setString(int fmi2ValueReference, const char* value)
{
  string s = string(__func__) + "(" + to_string(fmi2ValueReference)
    + ", " + string(value) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  /* We need to make a deep-copy of a provided string (to avoid a bug
     in the default generated-code). As the initial strings probably
     point to read-only strings, we cannot free them. A map could be
     used to determine which variables have been set or not, but it would
     need to represent a global variable and not be part of the
     code_saturne class, since the code_saturne.h file is generated.

     Since we only expect to use string variables for a few fixed
     environment  values, we accept a small (one-time per string
     in the model) memory leak when setting them, to avoid needlessly
     complex code. */

  char *p_var = NULL;

#if 0
  p_var = stringVariables[fmi2ValueReference].getter(this);
  if (p_var != NULL)
    free(p_var);
#endif

  size_t l = strlen(value);
  p_var = (char *)(malloc(l+1));
  strncpy(p_var, value, l);
  p_var[l] = '\0';

  stringVariables[fmi2ValueReference].setter(this, p_var);
}

/*----------------------------------------------------------------------------*/
