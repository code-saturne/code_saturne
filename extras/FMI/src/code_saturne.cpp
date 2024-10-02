/*============================================================================
 * Functional Mock-up Unit main methods implementation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.
  Copyright (C) 1998-2024 EDF S.A.

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

#include <map>
#include <functional>
#include <iostream>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include <src/callbacks.h>

typedef enum {
    exact,
    approx,
    calculated
} fmi_initial_t;

typedef enum {
    constant,
    fixed,
    tunable,
    discrete,
    continuous
} fmi_variability_t;

typedef enum {
    parameter,
    calculatedParamater,
    local,
    independent,
    input,
    output
} fmi_causality_t;

typedef enum {
    Real,
    Integer,
    Boolean,
    String
} fmi_type_t;

typedef struct {
  int                variable_index;  // Index in global variables list
  int                value_reference;
  fmi_causality_t    causality;
  fmi_variability_t  variability;
  fmi_initial_t      initial;
  double             start;
} cs_real_var_t;

typedef struct {
  int                variable_index;  // Index in global variables list
  int                value_reference;
  fmi_causality_t    causality;
  fmi_variability_t  variability;
  fmi_initial_t      initial;
  int                start;
} cs_integer_var_t;

typedef struct {
  int                variable_index;  // Index in global variables list
  int                value_reference;
  fmi_causality_t    causality;
  fmi_variability_t  variability;
  fmi_initial_t      initial;
  bool               start;
} cs_bool_var_t;


typedef struct {
  int                variable_index;  // Index in global variables list
  int                value_reference;
  fmi_causality_t    causality;
  fmi_variability_t  variability;
  fmi_initial_t      initial;
  const char *       start;
} cs_string_var_t;

#include "code_saturne_fmu_variables.h"

/*----------------------------------------------------------------------------*/

using namespace std;

#define CS_DRY_RUN 0   // Set to 1 to for simulation mode with no actual
                       // connection to code_saturne.

#define CS_TIMING 0    // Log some timings.

/*----------------------------------------------------------------------------
 * Macro used to silence "unused argument" warnings.
 *
 * This is useful when a function must match a given function pointer
 * type, but does not use all possible arguments.
 *----------------------------------------------------------------------------*/

#define CS_UNUSED(x) (void)(x)

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

  CS_LOG_EVENT,    /* Event */
  CS_LOG_WARNING,  /* Warning */
  CS_LOG_ERROR,    /* Error */
  CS_LOG_FATAL,    /* Fatal error */
  CS_LOG_PENDING,  /* Pending (asynchronous cases) */
  CS_LOG_LAUNCH,   /* Log system calls */
  CS_LOG_COMM,     /* Log communication (low-level) */
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

/* Save state of read and written variables */

typedef struct {

  double  *vals;

} cs_state_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static std::map<int, int> _real_variable_reference_map;
static std::map<int, int> _integer_variable_reference_map;
static std::map<int, int> _bool_variable_reference_map;
static std::map<int, int> _string_variable_reference_map;

int *_state_p = nullptr;
static std::map<fmi2FMUstate, fmi2FMUstate> _states;

static double *_real_variable_values = nullptr;
static int *_integer_variable_values = nullptr;
static bool *_bool_variable_values = nullptr;
static char **_string_variable_values = nullptr;

static int _n_iter = 0;
static int _master_socket = -1;
static int _cs_socket = -1;
static int _cs_swap_endian = 0;

static fmi2ComponentEnvironment  _component_environment = nullptr;
static fmi2Char                  _instance_name_init[] = "[unnamed]";
static fmi2Char *                _instance_name = _instance_name_init;

static FILE *tracefile = nullptr;

/* Mapping to default log categories
   categories up to "logStatusPending" are standardized, the rest are local.
   The "logAll" is handled separately */

static const size_t _n_log_categories = 8;
static const char *_log_categories[] = {"logEvents",
                                        "logStatusWarning",
                                        "logStatusError",
                                        "logStatusFatal",
                                        "logStatusPending",
                                        "logLaunch",
                                        "logComm",
                                        "logTrace"};

static const char *_log_prefix[] = {"[event]   ",
                                    "[warning] ",
                                    "[error]   ",
                                    "[fatal]   ",
                                    "[pending] ",
                                    "[launch]  ",
                                    "[comm]    ",
                                    "[trace]   "};

/* Active logging: 0: inactive, 1, log using FMI, -1: low-level trace,
   - 2: low level trace if FMI logging not active, changed to 1
   if FMI debug logging activated.

   Actual initialization to default values in fmi2Instantiate
   (preferred to static initialization in case FMI is restarted) */

static int _log_active[] = {-1, -1, -1, -1, -1, -1, -1, -1};

/* Read_queue */

cs_control_queue_t *_cs_glob_control_queue = nullptr;

/* Structure for grouped notebook value settings */

cs_variables_t *_cs_variables = nullptr;

/* Serialization (state) data */

size_t _serialized_size = 0;
char  *_serialized_data = nullptr;

/* Variables from generated template */

fmi2Component *_component = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log messages. We can also send some messages directly to stdout if
   not logged by the FMI environment. */

void
_cs_log(fmi2Status       status,
        cs_log_types_t   log_type,
        const char      *text)
{
  if (_log_active[log_type] > 0)
    log(_component_environment, _instance_name,
        status, _log_categories[log_type], text, nullptr);
  else if (_log_active[log_type] < 0)
    cout << _log_prefix[log_type] << _instance_name << ": " << text << endl;
}

void
_cs_log(fmi2Status       status,
        cs_log_types_t   log_type,
        string           text)
{
  if (_log_active[log_type] > 0)
    log(_component_environment, _instance_name,
        status, _log_categories[log_type], text.c_str(), nullptr);
  else if (_log_active[log_type] < 0)
    cout << _log_prefix[log_type] << _instance_name << ": " << text << endl;
}

void
_cs_log(fmi2ComponentEnvironment  environment,
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

void
_cs_log(fmi2ComponentEnvironment  environment,
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
  cs_control_queue_t *queue = nullptr;

  queue = (cs_control_queue_t *)malloc(sizeof(cs_control_queue_t));

  queue->buf = nullptr;

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
  if (queue != nullptr) {
    if (*queue == nullptr)
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

  if (tracefile != nullptr) {
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
    else if (tracefile != nullptr) {
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

  if (tracefile != nullptr) {
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
  else if (tracefile != nullptr) {
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

  if (queue != nullptr) {
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

    if (tracefile != nullptr) {
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
    else if (tracefile != nullptr) {
      fprintf(tracefile, "   received %d bytes: [", (int)ret);
      if (size == 1)
        _trace_buf(tracefile, buffer, ret);
      fprintf(tracefile, "]\n");
      fflush(tracefile);
    }

    start_id += ret;

  }

  if (tracefile != nullptr && size > 1) {
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

  if (tracefile != nullptr) {
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

    if (tracefile != nullptr) {
      fprintf(tracefile, "== receiving up to %d bytes...\n",
              (int)n_loc_max);
    }

#if CS_DRY_RUN == 1
    queue->buf[start_id] = '\0';
    return queue->buf;
#endif

    ssize_t ret = (socket != 0) ?
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

    if (tracefile != nullptr) {
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
  _recv_sock(_cs_socket, key_buffer, nullptr, 1, len_key);

  if (strncmp(key_s.c_str(), key_buffer, len_key) != 0) {
    _cs_log(fmi2Fatal, CS_LOG_FATAL,
            "wrong key received (socket handshake)");
    exit(EXIT_FAILURE);
  }

  _recv_sock(_cs_socket, magic_buffer, nullptr, 1, magic_string.size());

  _send_sock_str(_cs_socket, magic_buffer);

  _cs_log(fmi2OK, CS_LOG_COMM,
          "Connection between FMU client and code_saturne established");

  /* Iteration OK */
  char buf_rcv[13];
  _recv_sock(_cs_socket, buf_rcv, nullptr,1, 13);

  delete[] key_buffer;
  delete[] magic_buffer;

  _cs_glob_control_queue = _queue_initialize();
}

/*----------------------------------------------------------------------------*/

/* Configuration of the master socket (client) */

static int
_configure_client(sockaddr_in  *address)
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
_init_client(string        path,
             sockaddr_in  *address)
{
  /* Master socket */
  int port;
  sockaddr_in foo;
  socklen_t foosize = sizeof(foo);
  string hostname = "127.0.0.1";

  /* Generation of the key */
  int key = _generate_key();

  /* Client configuration */
  _master_socket = _configure_client(address);

  /* Get the actual port */
  getsockname(_master_socket, (sockaddr *)&foo, &foosize);
  port = ntohs(foo.sin_port);

  string s0 = "Listener on port ";
  string s1 = to_string(port);
  _cs_log(fmi2OK, CS_LOG_COMM, s0 + s1);

  _write_control_file(path, hostname, port, key);

  return key;
}

/*----------------------------------------------------------------------------*/

static void *
_start_client(void *data)
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

  string s = "New connection, socket fd: " + to_string(_cs_socket);
  s += "; ip: " + string(inet_ntoa(address.sin_addr));
  s += "; port: " + to_string(ntohs(address.sin_port));;

  _cs_log(fmi2OK, CS_LOG_COMM, s);

  pthread_exit(nullptr);
}

/*----------------------------------------------------------------------------*/

static void
*_start_cs(void *resu)
{
  string &res = *(static_cast<string*>(resu));
  string cmd;

  /* Wait for the client to be listening */
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

static void
_map_initialize(void)
{
  _real_variable_values = (double *)malloc(_n_real_vars*sizeof(double));
  for (int i = 0; i < _n_real_vars; i++) {
    int j = _real_vars[i].value_reference;
    _real_variable_reference_map[j] = i;
    _real_variable_values[i] = _real_vars[i].start;
  }
  _integer_variable_values = (int *)malloc(_n_integer_vars*sizeof(int));
  for (int i = 0; i < _n_integer_vars; i++) {
    int j = _integer_vars[i].value_reference;
    _integer_variable_reference_map[j] = i;
    _integer_variable_values[i] = _integer_vars[i].start;
  }
  _bool_variable_values = (bool *)malloc(_n_bool_vars*sizeof(bool));
  for (int i = 0; i < _n_bool_vars; i++) {
    int j = _bool_vars[i].value_reference;
    _bool_variable_reference_map[j] = i;
    _bool_variable_values[i] = _bool_vars[i].start;
  }
  _string_variable_values = (char **)malloc(_n_string_vars*sizeof(char *));
  for (int i = 0; i < _n_string_vars; i++) {
    int j = _string_vars[i].value_reference;
    _string_variable_reference_map[j] = i;
    size_t l = strlen(_string_vars[i].start);
    _string_variable_values[i] = (char *)malloc(l + 1);
    strncpy(_string_variable_values[i], _string_vars[i].start, l);
  }
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
  cs_variables_t *v = nullptr;

  v = (cs_variables_t *)malloc(sizeof(cs_variables_t));

  v->n_input = 0;
  v->n_output = 0;

  v->input_max = 0;
  v->output_max = 0;

  v->input_max_size = 0;
  v->output_max_size = 0;

  v->input_ids = nullptr;
  v->output_ids = nullptr;

  v->input_vals = nullptr;
  v->output_vals = nullptr;

  return v;
}

/*----------------------------------------------------------------------------
 * Finalize a variables grouping structure
 *----------------------------------------------------------------------------*/

static void
_variables_finalize(cs_variables_t  **v)
{
  if (v != nullptr) {
    if (*v == nullptr)
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
  (void)gettimeofday(&tv_time_0, nullptr);
#endif

  cs_variables_t *v = _cs_variables;

  string buffer = "advance " + to_string(n);

  string s = "Advancing " + to_string(n) + " iteration(s)";
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
    (void)gettimeofday(&tv_time_1, nullptr);

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

/* Query serialized snapshot from the server */

static void
_get_serialized_snapshot(int     sock)
{
  string buffer = "snapshot_get_serialized";

  string s_log = string(__func__) + ": send query for " + buffer;
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  size_t restart_size[1] = {0};

  /* Add input and output variables to end of snapshot */

  cs_variables_t *v = _cs_variables;
  size_t n_add = (v->n_input + v->n_output) * sizeof(double);

  if (_serialized_data != nullptr)
    free(_serialized_data);

#if CS_DRY_RUN == 1

  _serialized_size = n_add;
  _serialized_data = malloc(n_add);

#else

  _send_sock_str(sock, buffer.c_str());

  s_log = string(__func__) + ": waiting for " + buffer;
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  /* Received get: val */
  char *buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);

  if (strncmp(buf_rcv, "serialized_snapshot", 19) != 0) {
    s_log =   string(__func__) + ": unexpected reply; " + string(buf_rcv);
    _cs_log(fmi2Error, CS_LOG_ERROR, s_log);
  }

  /* Data size */
  _recv_sock(_cs_socket, (char *)restart_size, _cs_glob_control_queue,
             sizeof(size_t), 1);

  _serialized_size = restart_size[0] + n_add;
  _serialized_data = (char *)malloc(_serialized_size);

  if (restart_size[0] > 0)
    _recv_sock(_cs_socket, _serialized_data, _cs_glob_control_queue,
               1, restart_size[0]);

  s_log = string(__func__) + ": received snapshot";
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  /* receive 0 */
  buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);
  if (strcmp(buf_rcv, "0") != 0) {
    string s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
    _cs_log(fmi2Warning, CS_LOG_WARNING, s);
  }

#endif

  char *p = _serialized_data + restart_size[0];
  memcpy(p, v->input_vals, v->n_input*sizeof(double));
  p += v->n_input*sizeof(double);
  memcpy(p, v->output_vals, v->n_output*sizeof(double));
}

/*----------------------------------------------------------------------------*/

/* Send 'disconnect ' to the server */

static void
_disconnect(void)
{
  char buffer[13] = "disconnect ";

  _send_sock_str(_cs_socket, buffer);

  _cs_log(fmi2OK, CS_LOG_COMM, "Disconnecting the controller...");
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

/*============================================================================
 * code_saturne FMU method definitions
 *============================================================================*/

void
setResourceLocation(const char* fmiResourceLocation)
{
  string s0(__func__);
  string s_rl(fmiResourceLocation);

  _cs_log(fmi2OK, CS_LOG_TRACE, s0 + s_rl);
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2EnterInitializationMode(fmi2Component  component)
{
  /* As the code generator used provides no other earlier entry point,
     global pointer to component environment (for logging) here. */

  _component_environment = component;

  string s_casename, s_run_id, s_code_saturne;

  for (int i = 0; i < _n_string_vars; i++) {
    if (strcmp(_string_names[i], "code_saturne") == 0)
      s_code_saturne = _string_vars[i].start;
    else if (strcmp(_string_names[i], "casename") == 0)
      s_casename = _string_vars[i].start;
    else if (strcmp(_string_names[i], "run_id") == 0)
      s_run_id = _string_vars[i].start;
  }

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
  if (p != nullptr) {
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
  _cs_glob_control_queue = _queue_initialize();

#else

  _cs_log(fmi2OK, CS_LOG_LAUNCH, string("launch: \"") + cmd + string("\""));
  _cs_log(fmi2OK, CS_LOG_LAUNCH,
          string("output: \"") + exec_popen(cmd.c_str()) + string("\""));

  int key = _init_client(resu, &address);

  /* Accept the incoming connection */
  addrlen = sizeof(address);

  _cs_log(component, _instance_name, fmi2OK, CS_LOG_COMM,
          "Waiting for connection...",
          nullptr);

  struct thread_data data;
  data.address = address;
  data.addrlen = addrlen;

  pthread_create(&thread_server, nullptr, _start_client,
                 static_cast<void*>(&data));
  pthread_create(&thread_cs, nullptr, _start_cs,
                 static_cast<void*>(&resu));

  pthread_attr_destroy(&attr);
  pthread_join(thread_server, &status);
  pthread_join(thread_cs, &status);

  _comm_with_saturne(key);

#endif

  /* Sleep 3 seconds just to make sure the connection is well established */
  sleep(3);

  /* Now exchange input and output values which are already declared */

  if (_cs_variables == nullptr)
    _cs_variables = _variables_initialize();

  for (int i = 0; i < _cs_variables->input_max; i++) {
    if (_cs_variables->input_ids[i] > -1) {
      int var_index = _real_variable_reference_map[i];
      _set_notebook_variable(_cs_socket,
                             _real_names[var_index],
                             _real_variable_values[var_index]);
    }
  }
  for (int i = 0; i < _cs_variables->output_max; i++) {
    if (_cs_variables->output_ids[i] > -1) {
      int var_index = _real_variable_reference_map[i];
      double value
        = _get_notebook_variable(_cs_socket,
                                 _real_names[var_index]);
      _real_variable_values[var_index] = value;
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

fmi2Status
fmi2DoStep(fmi2Component  c,
           fmi2Real       currentComunicationTime,
           fmi2Real       stepSize,
           fmi2Boolean    noSetFmuFmuStatePriorToCurrentPoint)
{
  CS_UNUSED(c);
  CS_UNUSED(currentComunicationTime);
  CS_UNUSED(noSetFmuFmuStatePriorToCurrentPoint);

  string s = "----- " + string(__func__) + "(" + to_string(stepSize) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  cs_variables_t *v = _cs_variables;

  /* Update input notebook values */

  for (int i = 0; i < v->input_max; i++) {
    int id = v->input_ids[i];
    if (id < 0)
      continue;
    int var_index = _real_variable_reference_map[i];
    v->input_vals[id] = _real_variable_values[var_index];
  }

  /* Advance 1 iteration */
  _advance(1);

  /* Update output notebook values */

  for (int i = 0; i < v->output_max; i++) {
    int id = v->output_ids[i];
    if (id < 0)
      continue;
    int var_index = _real_variable_reference_map[i];
    _real_variable_values[var_index] = v->output_vals[id];
  }

  _n_iter++;

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status fmi2Terminate(fmi2Component  c)
{
  CS_UNUSED(c);

  string s_casename, s_run_id;

  for (int i = 0; i < _n_string_vars; i++) {
    if (strcmp(_string_names[i], "casename") == 0)
      s_casename = _string_vars[i].start;
    else if (strcmp(_string_names[i], "run_id") == 0)
      s_run_id = _string_vars[i].start;
  }

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

  if (tracefile != nullptr) {
    if (tracefile != stdout && tracefile != stderr) {
      fclose(tracefile);
      tracefile = nullptr;
    }
  }

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2Reset(fmi2Component  c)
{
  CS_UNUSED(c);

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Real
fmi2GetReal(fmi2ValueReference  reference)
{
  if (_cs_variables == nullptr)
    _cs_variables = _variables_initialize();

  string s = string(__func__) + "(" + to_string(reference) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  int var_index = _real_variable_reference_map[reference];

  if (_real_vars[var_index].causality != input) {
    int first = _variables_add_output(_cs_variables, reference);
    if (first && _cs_socket >= 0) {
      double value = _get_notebook_variable(_cs_socket,
                                            _real_names[var_index]);
      _real_variable_values[var_index] = value;
    }
  }

  return _real_variable_values[var_index];
}

/*----------------------------------------------------------------------------*/

void
fmi2SetReal(fmi2ValueReference  reference,
            fmi2Real            value)
{
  if (_cs_variables == nullptr)
    _cs_variables = _variables_initialize();

  string s = string(__func__) + "(" + to_string(reference)
    + ", " + to_string(value) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  int var_index = _real_variable_reference_map[reference];

  if (_real_vars[var_index].causality == input) {
    _real_variable_values[var_index] = value;

    int first = _variables_add_input(_cs_variables,
                                     reference,
                                     value);
    if (first && _cs_socket >= 0) {
      _set_notebook_variable(_cs_socket,
                             _real_names[var_index],
                             value);
    }
  }
}

/*----------------------------------------------------------------------------*/

fmi2Integer
fmi2GetInteger(fmi2ValueReference  reference)
{
  int var_index = _integer_variable_reference_map[reference];

  return _integer_variable_values[var_index];
}

/*----------------------------------------------------------------------------*/

void
fmi2SetInteger(fmi2ValueReference  reference,
               fmi2Integer         value)
{
  int var_index = _integer_variable_reference_map[reference];

  _integer_variable_values[var_index] = value;
}

/*----------------------------------------------------------------------------*/

fmi2Integer
fmi2GetBoolean(fmi2ValueReference  reference)
{
  int var_index = _bool_variable_reference_map[reference];

  return _bool_variable_values[var_index];
}

/*----------------------------------------------------------------------------*/

void
fmi2SetBoolean(fmi2ValueReference  reference,
               fmi2Boolean         value)
{
  int var_index = _bool_variable_reference_map[reference];

  _bool_variable_values[var_index] = value;
}

/*----------------------------------------------------------------------------*/

fmi2String
fmi2GetString(fmi2ValueReference  reference)
{
  int var_index = _string_variable_reference_map[reference];

  return _string_variable_values[var_index];
}

/*----------------------------------------------------------------------------*/

void fmi2SetString(fmi2ValueReference  reference,
                   const char*         value)
{
  string s = string(__func__) + "(" + to_string(reference)
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

  int var_index = _string_variable_reference_map[reference];

  char *p_var = _string_variable_values[var_index];
  if (p_var != nullptr)
    free(p_var);

  size_t l = strlen(value);
  p_var = (char *)(malloc(l+1));
  strncpy(p_var, value, l);
  p_var[l] = '\0';

  _string_variable_values[var_index] = p_var;
}

/*----------------------------------------------------------------------------*/

fmi2Status fmi2GetFMUstate(fmi2Component  c,
                           fmi2FMUstate*  state)
{
  CS_UNUSED(c);

  char s[128];
  snprintf(s, 127, "%s: called for %p; only client state is partially saved",
           __func__, (void *)*state);
  _cs_log(fmi2Warning, CS_LOG_WARNING, s);

  /* Not a valid pointer to actual data here,
     just a distinct pointer per state... */

  cs_variables_t *v = _cs_variables;
  cs_state_t *csst = nullptr;

  if (*state == nullptr) {
    csst = (cs_state_t *)malloc(sizeof(cs_state_t));
    csst->vals = (double *)malloc(sizeof(double) * (v->n_input+v->n_output));
    fmi2FMUstate fs = (fmi2FMUstate)csst;
    _states[fs] = fs;

    *state = fs;
  }
  else {
    if (_states.count(*state) > 0)
      csst = (cs_state_t *)_states[*state];
    else {
      snprintf(s, 127, "%s: called for %p, not previously defined",
               __func__, (void *)*state);
      _cs_log(fmi2Error, CS_LOG_ERROR, s);
      return fmi2Error;
    }
  }

  {
    int j = 0;
    for (int i = 0; i < v->n_input; i++)
      csst->vals[j++] = v->input_vals[i];
    for (int i = 0; i < v->n_output; i++)
      csst->vals[j++] = v->output_vals[i];
  }

  return fmi2Warning;
}

/*----------------------------------------------------------------------------*/

fmi2Status fmi2SetFMUstate(fmi2Component  c,
                           fmi2FMUstate   state)
{
  CS_UNUSED(c);

  char s[128];
  snprintf(s, 127, "%s: called for %p; only client state is partially set",
           __func__, (void *)state);
  _cs_log(fmi2Warning, CS_LOG_WARNING, s);

  cs_variables_t *v = _cs_variables;
  cs_state_t *csst = nullptr;

  if (_states.count(state) > 0)
    csst = (cs_state_t *)_states[state];
  else {
    snprintf(s, 127, "%s: called for %p, not previously defined",
             __func__, (void *)state);
    _cs_log(fmi2Error, CS_LOG_ERROR, s);
    return fmi2Error;
  }

  {
    int j = 0;
    for (int i = 0; i < v->n_input; i++)
      v->input_vals[i] = csst->vals[j++];
    for (int i = 0; i < v->n_output; i++)
      v->output_vals[i] = csst->vals[j++];
  }

  return fmi2Warning;
}

/*----------------------------------------------------------------------------*/

fmi2Status fmi2FreeFMUstate(fmi2Component  c,
                            fmi2FMUstate*  state)
{
  CS_UNUSED(c);

  if (*state == nullptr)
    return fmi2OK;

  char s[128];

  if (_states.count(*state) > 0) {
    snprintf(s, 127, "%s: called for %p\n",
             __func__, (void *)*state);
    _cs_log(fmi2OK, CS_LOG_TRACE, s);
    cs_state_t *csst = (cs_state_t *)_states[*state];
    free(csst->vals);
    free(csst);
    _states.erase(*state);
    *state = nullptr;
  }
  else {
    snprintf(s, 127, "%s: called for %p, which is not present\n",
           __func__, (void *)*state);
    _cs_log(fmi2Error, CS_LOG_WARNING, s);
  }

  return fmi2Warning;
}

/*----------------------------------------------------------------------------*/

fmi2Status fmi2SerializedFMUstateSize(fmi2Component  c,
                                      fmi2FMUstate   state,
                                      size_t*        stateSize)
{
  CS_UNUSED(c);

  _get_serialized_snapshot(_cs_socket);

  *stateSize = _serialized_size;
  char s[128];
  snprintf(s, 127, "%s: called for %p (%ld bytes)",
           __func__, (void *)state, (long)_serialized_size);
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2SerializeFMUstate(fmi2Component  c,
                      fmi2FMUstate   state,
                      fmi2Byte       serializedState[],
                      size_t         serializedStateSize)
{
  CS_UNUSED(c);

  char s[128];

  if (serializedStateSize != _serialized_size) {
    snprintf(s, 127,
             "%s: called for %p\n"
             "with size %d (%d expected)",
             __func__, (void *)state,
             (int)serializedStateSize, (int)_serialized_size);
    _cs_log(fmi2Error, CS_LOG_TRACE, s);
    return fmi2Error;
  }

  memcpy(serializedState, _serialized_data, _serialized_size);

  _serialized_size = 0;
  free(_serialized_data);

  snprintf(s, 127, "%s: called for %p", __func__, (void *)state);
  _cs_log(fmi2Error, CS_LOG_TRACE, s);
  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2DeSerializeFMUstate(fmi2Component   c,
                        const fmi2Byte  serializedState[],
                        size_t          size,
                        fmi2FMUstate*   state)
{
  CS_UNUSED(c);
  CS_UNUSED(state);  // Assumed to be handled by caller only.

  char s[128];
  snprintf(s, 127, "%s: called for %p\n", __func__, (void *)state);
  _cs_log(fmi2OK, CS_LOG_ERROR, s);

  /* Add input and output variables to end of snapshot */

  cs_variables_t *v = _cs_variables;

  char *p = (char *)serializedState;
  p += size - v->n_output*sizeof(double);
  memcpy(v->output_vals, p, v->n_output*sizeof(double));
  p -= v->n_input*sizeof(double);
  memcpy(v->input_vals, p, v->n_input*sizeof(double));

  size_t restart_size = size - (v->n_output + v->n_input)*sizeof(double);

  /* Now send rest of snapshot to server */

  string buffer = "snapshot_load_serialized " + to_string(restart_size);

  string s_log = string(__func__) + ": send command : " + buffer;
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  _send_sock_str(_cs_socket, buffer.c_str());

  _send_sock(_cs_socket, (char *)serializedState, 1, restart_size);

  return fmi2OK;
}

/*============================================================================
 * General FMI function definitions
 *============================================================================*/

const char*
fmi2GetTypesPlatform()
{
  return "default";
}

const char*
fmi2GetVersion()
{
  return "2.0";
}

fmi2Status
fmi2SetDebugLogging(fmi2Component     c,
                    fmi2Boolean       loggingOn,
                    size_t            nCategories,
                    const fmi2String  categories[])
{
  CS_UNUSED(c);

  if (loggingOn) {
    if (nCategories == 0) {
      for (size_t i = 0; i < _n_log_categories; i++) {
        if (_log_active[i] == -2)
          _log_active[i] = 1;
      }
    }

    // Categories not provided/read in XML yet
    else {
      for (size_t i = 0; i < nCategories; i++) {
        const fmi2String s = categories[i];
        if (s == nullptr)
          continue;
        else if (strcmp(s, "logAll") == 0) {
          for (size_t k = 0; k <= _n_log_categories; k++) {
            if (_log_active[k] == -2)
              _log_active[k] = 1;
          }
          continue;
        }
        else {
          size_t j = 0;
          while (j < _n_log_categories) {
            if (strcmp(s, _log_categories[j]) != 0) {
              _log_active[j] = 1;
              break;
            }
            j++;
          }
          if (j >= _n_log_categories)
            cout << string(__func__) + ": category '" + string(s)  \
              + "' ignored" << endl;
        }
      }
    }
  }

  /* If FMI logging off, keep low level log for now */
  else {
    for (size_t i = 0; i < _n_log_categories; i++) {
      if (_log_active[i] != 0) {
        _log_active[i] = -1;
      }
    }
  }

  return fmi2OK;
}

fmi2Component
fmi2Instantiate(fmi2String   instanceName,
                fmi2Type     fmuType,
                fmi2String   fmuGUID,
                fmi2String   fmuResourceLocation,
                const fmi2CallbackFunctions*  callbacks,
                fmi2Boolean  visible,
                fmi2Boolean  loggingOn)
{
  CS_UNUSED(fmuType);
  CS_UNUSED(fmuGUID);
  CS_UNUSED(visible);

  _log_active[CS_LOG_EVENT] = -2;
  _log_active[CS_LOG_WARNING] = -2;
  _log_active[CS_LOG_ERROR] = -2;
  _log_active[CS_LOG_FATAL] = -2;
  _log_active[CS_LOG_PENDING] = -2;
  _log_active[CS_LOG_LAUNCH] = -2;
  _log_active[CS_LOG_COMM] = -2;
  _log_active[CS_LOG_TRACE] = -1;

  setFunctions(callbacks);

  _map_initialize();

  _instance_name = (fmi2Char *)malloc(sizeof(instanceName) + 1);
  strcpy(_instance_name, instanceName);

  setResourceLocation(fmuResourceLocation);
  _component = (fmi2Component*) fmiAlloc(sizeof(fmi2Component));

  if (loggingOn) {
    for (size_t i = 0; i < _n_log_categories; i++) {
      if (_log_active[i] == -2) {
        _log_active[i] = 1;
      }
    }
  }

  return _component;
}

fmi2Status
fmi2SetupExperiment(fmi2Component  c,
                    fmi2Boolean    toleranceDefined,
                    fmi2Real       tolerance,
                    fmi2Real       startTime,
                    fmi2Boolean    stopTimeDefined,
                    fmi2Real       stopTime)
{
  CS_UNUSED(c);
  CS_UNUSED(toleranceDefined);
  CS_UNUSED(tolerance);
  CS_UNUSED(startTime);
  CS_UNUSED(stopTimeDefined);
  CS_UNUSED(stopTime);

  return fmi2OK;
}

fmi2Status
fmi2ExitInitializationMode(fmi2Component  c)
{
  CS_UNUSED(c);

  return fmi2OK;
}

void
fmi2FreeInstance(fmi2Component  c)
{
  CS_UNUSED(c);

  fmiFree(c);
}

fmi2Status
fmi2GetReal(fmi2Component             c,
            const fmi2ValueReference  valueReference[],
            size_t                    numberOfValues,
            fmi2Real                  values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = fmi2GetReal(valueReference[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetInteger(fmi2Component             c,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               fmi2Integer               values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = fmi2GetInteger(valueReference[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetBoolean(fmi2Component             c,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               fmi2Boolean               values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = fmi2GetBoolean(valueReference[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetString(fmi2Component             c,
              const fmi2ValueReference  valueReference[],
              size_t                    numberOfValues,
              fmi2String                values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = fmi2GetString(valueReference[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetReal(fmi2Component             c,
            const fmi2ValueReference  valueReference[],
            size_t                    numberOfValues,
            const fmi2Real            values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      fmi2SetReal(valueReference[i], values[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetInteger(fmi2Component             c,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               const fmi2Integer         values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      fmi2SetInteger(valueReference[i], values[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetBoolean(fmi2Component             c,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               const fmi2Boolean         values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      fmi2SetBoolean(valueReference[i], values[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetString(fmi2Component             c,
              const fmi2ValueReference  valueReference[],
              size_t                    numberOfValues,
              const fmi2String          values[])
{
  CS_UNUSED(c);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      fmi2SetString(valueReference[i], values[i]);
    return fmi2OK;
  } catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetDirectionalDerivative(fmi2Component             c,
                             const fmi2ValueReference  unknownValueReferences[],
                             size_t                    numberOfUnknowns,
                             const fmi2ValueReference  knownValueReferences[],
                             fmi2Integer               numberOfKnowns,
                             fmi2Real                  knownDifferential[],
                             fmi2Real                  unknownDifferential[])
{
  CS_UNUSED(c);
  CS_UNUSED(unknownValueReferences);
  CS_UNUSED(numberOfUnknowns);
  CS_UNUSED(knownValueReferences);
  CS_UNUSED(numberOfKnowns);
  CS_UNUSED(knownDifferential);
  CS_UNUSED(unknownDifferential);

  return fmi2Error;
}

fmi2Status
fmi2SetRealInputDerivatives(fmi2Component             c,
                            const fmi2ValueReference  valueReferences[],
                            size_t                    numberOfValueReferences,
                            fmi2Integer               orders[],
                            const fmi2Real            values[])
{
  CS_UNUSED(c);
  CS_UNUSED(valueReferences);
  CS_UNUSED(numberOfValueReferences);
  CS_UNUSED(orders);
  CS_UNUSED(values);

  return fmi2Error;
}

fmi2Status
fmi2GetRealOutputDerivatives(fmi2Component             c,
                             const fmi2ValueReference  valueReference[],
                             size_t                    numberOfValues,
                             const fmi2Integer         order[],
                             fmi2Real                  values[])
{
  CS_UNUSED(c);
  CS_UNUSED(valueReference);
  CS_UNUSED(numberOfValues);
  CS_UNUSED(order);
  CS_UNUSED(values);

  return fmi2Error;
}

fmi2Status
fmi2CancelStep(fmi2Component  c)
{
  CS_UNUSED(c);

  return fmi2Warning;
}

fmi2Status
fmi2GetStatus(fmi2Component         c,
              const fmi2StatusKind  kind,
              fmi2Status*           status)
{
  CS_UNUSED(c);
  CS_UNUSED(kind);
  CS_UNUSED(status);

  return fmi2Error;
}

fmi2Status
fmi2GetRealStatus(fmi2Component         c,
                  const fmi2StatusKind  kind,
                  fmi2Real*             value)
{
  CS_UNUSED(c);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

fmi2Status
fmi2GetIntegerStatus(fmi2Component         c,
                     const fmi2StatusKind  kind,
                     fmi2Integer*          value)
{
  CS_UNUSED(c);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

fmi2Status
fmi2GetBooleanStatus(fmi2Component         c,
                     const fmi2StatusKind  kind,
                     fmi2Boolean*          value)
{
  CS_UNUSED(c);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

fmi2Status
fmi2GetStringStatus(fmi2Component         c,
                    const fmi2StatusKind  kind,
                    fmi2String*           value)
{
  CS_UNUSED(c);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

/*----------------------------------------------------------------------------*/
