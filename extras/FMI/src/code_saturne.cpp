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
#include <thread>

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

  bool    *input_init;

} cs_variables_t;

/* Save state of read and written variables */

typedef struct {

  double  *vals;

} cs_state_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
    *queue = nullptr;
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

/* Equivalent of "system" but is able to get stdout */

static string
_exec_popen(const char* cmd)
{
  array<char, 128> buffer;
  string result;
  unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);

  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }

  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

/*============================================================================
 * Base class for code_saturne FMU component.
 *============================================================================*/

class cs_fmu {

public:

  bool _initialization_mode = false;

  std::map<int, int> _real_variable_reference_map;
  std::map<int, int> _integer_variable_reference_map;
  std::map<int, int> _bool_variable_reference_map;
  std::map<int, int> _string_variable_reference_map;

  int *_state_p = nullptr;
  std::map<fmi2FMUstate, fmi2FMUstate> _states;

  double *_real_variable_values = nullptr;
  int *_integer_variable_values = nullptr;
  bool *_bool_variable_values = nullptr;
  char **_string_variable_values = nullptr;

  int _n_iter = 0;
  int _master_socket = -1;
  int _cs_socket = -1;
  int _cs_swap_endian = 0;

  sockaddr_in  _client_address;

  fmi2ComponentEnvironment  _component_environment = nullptr;
  fmi2Char                 *_instance_name;
  fmi2Char                 *_resource_location;

  FILE *_tracefile = nullptr;

  /* Active logging: 0: inactive, 1, log using FMI, -1: low-level trace,
     - 2: low level trace if FMI logging not active, changed to 1
     if FMI debug logging activated.

   Actual initialization to default values in fmi2Instantiate
   (preferred to static initialization in case FMI is restarted) */

  int _log_active[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

  /* Read_queue */

  cs_control_queue_t *_control_queue = nullptr;

  /* Structure for grouped notebook value settings */

  cs_variables_t *_cs_variables = nullptr;

  /* Serialization (state) data */

  size_t _serialized_size = 0;
  char  *_serialized_data = nullptr;

  //==========================================================================
  // Default contructor and destructor
  //==========================================================================

  cs_fmu(fmi2String   instanceName,
         fmi2String   fmuResourceLocation) {

    _instance_name = new fmi2Char[strlen(instanceName) + 1];
    strcpy(_instance_name, instanceName);

    _resource_location = new fmi2Char[strlen(fmuResourceLocation) + 1];
    strcpy(_resource_location, fmuResourceLocation);

    _log_active[CS_LOG_EVENT] = -2;
    _log_active[CS_LOG_WARNING] = -2;
    _log_active[CS_LOG_ERROR] = -2;
    _log_active[CS_LOG_FATAL] = -2;
    _log_active[CS_LOG_PENDING] = -2;
    _log_active[CS_LOG_LAUNCH] = -2;
    _log_active[CS_LOG_COMM] = -2;
    _log_active[CS_LOG_TRACE] = -1;

    _cs_variables = _variables_initialize();

 }

  ~cs_fmu() {

    delete[] _instance_name;
    delete[] _resource_location;

    _variables_finalize();

  }

  //==========================================================================
  // Member functions
  //==========================================================================

  /* Log messages. We can also send some messages directly to stdout if
     not logged by the FMI environment. */

  void
  _log(fmi2Status                status,
       cs_log_types_t            log_type,
       const char               *text)
  {
    if (_log_active[log_type] > 0)
      log(_component_environment, _instance_name,
          status, _log_categories[log_type], text, nullptr);
    else if (_log_active[log_type] < 0)
      cout << _log_prefix[log_type] << _instance_name << ": " << text << endl;
  }

  void
  _log(fmi2Status                status,
       cs_log_types_t            log_type,
       string                    text)
  {
    if (_log_active[log_type] > 0)
      log(_component_environment, _instance_name,
          status, _log_categories[log_type], text.c_str(), nullptr);
    else if (_log_active[log_type] < 0)
      cout << _log_prefix[log_type] << _instance_name << ": " << text << endl;
  }

  void
  _log(fmi2ComponentEnvironment  environment,
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
  _log(fmi2ComponentEnvironment  environment,
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

  /* Send data on a socket */

  void
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

    if (_tracefile != nullptr) {
      if (size == 1) {
        fprintf(_tracefile, "== send %d bytes: [", (int)n_bytes);
        _trace_buf(_tracefile, buffer, n_bytes);
        fprintf(_tracefile, "]...\n");
      }
      else {
        fprintf(_tracefile, "== send %d values of size %d:\n", (int)ni, (int)size);
        for (size_t i = 0; i < ni; i++) {
          fprintf(_tracefile, "    ");
          for (size_t j = 0; j < size; j++)
            fprintf(_tracefile, " %x", (unsigned)buffer[i*size + j]);
          fprintf(_tracefile, "\n");
        }
      }
      fflush(_tracefile);
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
        _log(fmi2Fatal, CS_LOG_FATAL, s0 + s2);
        throw std::runtime_error(s0 + s2);
      }
      else if (_tracefile != nullptr) {
        fprintf(_tracefile, "   sent %d bytes\n", (int)ret);
        fflush(_tracefile);
      }

      start_id += ret;
    }

    if (_cs_swap_endian && size > 1)  /* restore endiannes */
      _swap_endian(buffer, size, ni);
  }

  /*----------------------------------------------------------------------------*/

  /* Send a message on a socket */

  void
  _send_sock_str(int          sock,
                 const char  *str)
  {
#if CS_DRY_RUN == 1
    return;
#endif

    size_t n = strlen(str)+1;

    if (_tracefile != nullptr) {
      fprintf(_tracefile, "== send %d bytes: [", (int)n);
      _trace_buf(_tracefile, str, n);
      fprintf(_tracefile, "]...\n");
      fflush(_tracefile);
    }

    ssize_t ret = send(sock, str, n, 0);

    if (ret < 1) {
      string s0 = "Error sending ";
      string s1 = str;
      string s2 = strerror(errno);
      string err_str = s0 + s1 + s2;
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }
    else if (_tracefile != nullptr) {
      fprintf(_tracefile, "   sent %d bytes\n", (int)ret);
      fflush(_tracefile);
    }
  }

  /*--------------------------------------------------------------------------*/

  /* Receive data message (fixed size) on a socket */

  void
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

      if (_tracefile != nullptr) {
        fprintf(_tracefile,
                "== receiving up to %d bytes, %d of %d bytes already buffered...\n",
                (int)n_loc,
                (int)start_id,
                (int)n_bytes);
      }

      ssize_t ret = recv(socket, buffer + start_id, n_loc, 0);

      if (ret < 1) {
        string s0 = "Error receiving data: ";
        string s1;
        if (errno != 0)
          s1 = strerror(errno);
        else
          s1 = "code_saturne may have disconnected.";
        _log(fmi2Fatal, CS_LOG_FATAL, s0 + s1);
        throw std::runtime_error(s0 + s1);
      }
      else if (_tracefile != nullptr) {
        fprintf(_tracefile, "   received %d bytes: [", (int)ret);
        if (size == 1)
          _trace_buf(_tracefile, buffer, ret);
        fprintf(_tracefile, "]\n");
        fflush(_tracefile);
      }

      start_id += ret;

    }

    if (_tracefile != nullptr && size > 1) {
      for (size_t i = 0; i < ni; i++) {
        fprintf(_tracefile, "    ");
        for (size_t j = 0; j < size; j++)
          fprintf(_tracefile, " %x", (unsigned)buffer[i*size + j]);
        fprintf(_tracefile, "\n");
      }
    }

    if (_cs_swap_endian && size > 1)
      _swap_endian(buffer, size, ni);
  }

  /*--------------------------------------------------------------------------*/

  /* Receive message (text expected) on a socket */

  char *
  _recv_sock_with_queue(int                  socket,
                        cs_control_queue_t  *queue,
                        size_t               min_size)
  {
    size_t start_id = queue->buf_idx[1] - queue->buf_idx[0];
    ssize_t cut_id = -1;

    if (_tracefile != nullptr) {
      fprintf(_tracefile,
              "_recv_sock_with_queue: %d %d\n",
              (int)start_id,
              (int)min_size);
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

      if (_tracefile != nullptr) {
        fprintf(_tracefile, "== receiving up to %d bytes...\n",
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
        _log(fmi2Fatal, CS_LOG_FATAL, s0 + s1);
        throw std::runtime_error(s0 + s1);
      }

      if (_tracefile != nullptr) {
        fprintf(_tracefile, "   received %d bytes: [", (int)ret);
        _trace_buf(_tracefile, queue->buf + start_id, ret);
        fprintf(_tracefile, "]\n");
        fflush(_tracefile);
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

  /*--------------------------------------------------------------------------*/

  /* Communication with code_saturne handshake at the beginning
   * (first connection) and reception of messages from CS after */

  void
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
      const char err_str[] = "wrong key received (socket handshake)";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }

    _recv_sock(_cs_socket, magic_buffer, nullptr, 1, magic_string.size());

    _send_sock_str(_cs_socket, magic_buffer);

    _log(fmi2OK, CS_LOG_COMM,
         "Connection between FMU client and code_saturne established");

    /* Iteration OK */
    char buf_rcv[13];
    _recv_sock(_cs_socket, buf_rcv, nullptr,1, 13);

    delete[] key_buffer;
    delete[] magic_buffer;

    _control_queue = _queue_initialize();
  }

  /*--------------------------------------------------------------------------*/

  /* Configuration of the master socket (client) */

  int
  _configure_client(void)
  {
    int master_socket;
    int opt = 1;
    socklen_t optlen = sizeof(opt);

    /* Create a master socket */
    if ((master_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
      const char err_str[] = "socket() failed";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }

    /* Set master socket to allow multiple connections */
    if (setsockopt(master_socket, SOL_SOCKET, SO_REUSEADDR, &opt, optlen) < 0) {
      const char err_str[] = "setsockopt() failed";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }

    /* Master socket configuration */
    _client_address.sin_family = AF_INET;
    _client_address.sin_addr.s_addr = htonl(INADDR_ANY);
    _client_address.sin_port = 0;

    /* Bind the socket to any available localhost port*/
    if (bind(master_socket,
             (sockaddr *)&_client_address,
             sizeof(_client_address)) < 0) {
      const char err_str[] = "bind() failed";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }

    return master_socket;
  }

  /*--------------------------------------------------------------------------*/

  /* Writing of control_files for code_saturne which is deleted after reading */

  void
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
      const char err_str[] = "Error writing control_file.";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }
  }

  /*--------------------------------------------------------------------------*/

  int
  _init_client(string   path)
  {
    /* Master socket */
    int port;
    sockaddr_in foo;
    socklen_t foosize = sizeof(foo);
    string hostname = "127.0.0.1";

    /* Generation of the key */
    int key = _generate_key();

    /* Client configuration */
    _master_socket = _configure_client();

    /* Get the actual port */
    getsockname(_master_socket, (sockaddr *)&foo, &foosize);
    port = ntohs(foo.sin_port);

    string s0 = "Listener on port ";
    string s1 = to_string(port);
    _log(fmi2OK, CS_LOG_COMM, s0 + s1);

    _write_control_file(path, hostname, port, key);

    return key;
  }

  /*--------------------------------------------------------------------------*/

  void *
  _start_client(void)
  {
    sockaddr_in address = _client_address;
    int addrlen = sizeof(address);

    /* Try to specify maximum of 1 pending connection *
       for the master socket                          */
    if (listen(_master_socket, 1) < 0) {
      const char err_str[] = "listen() failed.";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
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
      const char err_str[] = "accept() failed.";
      _log(fmi2Fatal, CS_LOG_FATAL, err_str);
      throw std::runtime_error(err_str);
    }

    /* Inform user of socket number - used in send and receive commands */

    string s = "New connection, socket fd: " + to_string(_cs_socket);
    s += "; ip: " + string(inet_ntoa(address.sin_addr));
    s += "; port: " + to_string(ntohs(address.sin_port));;

    _log(fmi2OK, CS_LOG_COMM, s);

    pthread_exit(nullptr);
  }

  /*--------------------------------------------------------------------------*/

  void
  _start_cs(void)
  {
    string resu = _get_resu_dir();
    string cmd;

    /* Wait for the client to be listening */
    sleep(1);

    /* Run Code_saturne client in background */
    cmd = "cd " + resu + " && ./run_solver &";
    system(cmd.c_str());
  }

  /*--------------------------------------------------------------------------
   * Initialize a variables grouping structure
   *
   * returns:
   *   pointer to initialized control queue
   *--------------------------------------------------------------------------*/

  cs_variables_t *
  _variables_initialize(void)
  {
    cs_variables_t *v = nullptr;

    v = new cs_variables_t;

    v->n_input = 0;
    v->n_output = 0;

    v->input_max = 0;
    v->output_max = 0;

    v->input_max_size = 0;
    v->output_max_size = 0;

    v->input_ids = nullptr;
    v->output_ids = nullptr;

    v->input_init = nullptr;

    v->input_vals = nullptr;
    v->output_vals = nullptr;

    return v;
  }

  /*--------------------------------------------------------------------------
   * Finalize a variables grouping structure
   *--------------------------------------------------------------------------*/

  void
  _variables_finalize(void)
  {
    if (_cs_variables != nullptr) {

      if (_cs_variables->input_max_size > 0) {
        free(_cs_variables->input_ids);
        free(_cs_variables->input_vals);
        free(_cs_variables->input_init);
      }
      if (_cs_variables->output_max_size > 0) {
        free(_cs_variables->output_ids);
        free(_cs_variables->output_vals);
      }

      delete _cs_variables; _cs_variables = nullptr;
    }
  }

  /*--------------------------------------------------------------------------
   * Add or set a variables input;
   *--------------------------------------------------------------------------*/

  void
  _variables_add_input(cs_variables_t  *v,
                       int              id,
                       double           val)
  {
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

        v->input_init = (bool *)realloc(v->input_init,
                                        v->input_max_size*sizeof(bool));

        for (int i = v->input_max; i < v->input_max_size; i++) {
          v->input_ids[i] = -1;
          v->input_vals[i] = 0;
          v->input_init[i] = false;
        }
      }

      v->input_max = id + 1;
    }

    if (v->input_ids[id] < 0) {
      v->input_ids[id] = v->n_input;
      v->n_input += 1;
    }

    v->input_vals[id] = val;
  }

  /*--------------------------------------------------------------------------
   * Set a variables input;
   *--------------------------------------------------------------------------*/

  int
  _variables_set_input(cs_variables_t  *v,
                       int              id,
                       double           val)
  {
    int retval = 0;

    assert(id >= 0);

    bool unknown = false;
    if (id >= v->input_max)
      unknown = true;
    else if (v->input_ids[id] < 0)
      unknown = true;
    if (unknown) {
      string s = "Input variable id " + to_string(id) + " unknown in model";
      _log(fmi2Fatal, CS_LOG_FATAL, s);
      throw std::runtime_error(s);
    }

    if (v->input_init[id] == false) {
      v->input_init[id] = true;
      retval = 1;
    }

    v->input_vals[id] = val;

    return retval;
  }

  /*--------------------------------------------------------------------------
   * Add an output variable;
   *--------------------------------------------------------------------------*/

  void
  _variables_add_output(cs_variables_t  *v,
                        int              id)
  {
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
      v->output_ids[id] = v->n_output;
      v->n_output += 1;
    }
  }

  /*--------------------------------------------------------------------------
   * Initialize a variables grouping structure
   *
   * returns:
   *   pointer to initialized control queue
   *--------------------------------------------------------------------------*/

  void
  _map_initialize(void)
  {
    _real_variable_values = (double *)malloc(_n_real_vars*sizeof(double));
    for (int i = 0; i < _n_real_vars; i++) {
      int j = _real_vars[i].value_reference;
      _real_variable_reference_map[j] = i;
      _real_variable_values[i] = _real_vars[i].start;
      if (_real_vars[i].causality == input)
        _variables_add_input(_cs_variables, j, _real_vars[i].start);
      else
        _variables_add_output(_cs_variables, j);
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

  /*--------------------------------------------------------------------------*/

  void
  _advance(int   n)
  {
#if CS_TIMING == 1
    struct timeval  tv_time_0, tv_time_1;
    (void)gettimeofday(&tv_time_0, nullptr);
#endif

    cs_variables_t *v = _cs_variables;

    string buffer = "advance " + to_string(n);

    string s = "Advancing " + to_string(n) + " iteration(s)";
    _log(fmi2OK, CS_LOG_TRACE, s);

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
    char *buf_rcv = _recv_sock_with_queue(_cs_socket, _control_queue, 0);
    if (strcmp(buf_rcv, "0") != 0) {
      s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
      _log(fmi2Warning, CS_LOG_WARNING, s);
    }

    /* Get output notebook values */

    if (v->n_output > 0) {
      string s = "Get " + to_string(v->n_output) + " output variables";
      _log(fmi2OK, CS_LOG_TRACE, s);
      _recv_sock(_cs_socket, (char *)_cs_variables->output_vals,
                 _control_queue,
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

      buf_rcv = _recv_sock_with_queue(_cs_socket, _control_queue, 0);
      if (strcmp(buf_rcv, "Iteration OK") != 0) {
        s =   string(__func__) + ": expected \"Iteration OK\", not: "
            + string(buf_rcv);
        _log(fmi2Warning, CS_LOG_WARNING, s);
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

  /*--------------------------------------------------------------------------*/

  /* Send 'notebook_add_input value' to the server */

  void
  _set_notebook_variable(int          sock,
                         const char  *variable,
                         double       value)
  {
    string var_s(variable);

    string buffer = "notebook_add_input " + var_s + " " + to_string(value);

    string s = string(__func__) + ": sending " + string(variable);
    _log(fmi2OK, CS_LOG_TRACE, s);

    _send_sock_str(sock, buffer.c_str());

#if CS_DRY_RUN == 1
    return;
#endif

    s = string(__func__) + ": waiting for reply";
    _log(fmi2OK, CS_LOG_TRACE, s);

    // Receive 0
    char *buf_rcv = _recv_sock_with_queue(sock, _control_queue, 0);

    if (strcmp(buf_rcv, "0") != 0) {
      s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
      _log(fmi2Warning, CS_LOG_WARNING, s);
    }
    else {
      s = string(__func__) + ": variable set.";
      _log(fmi2OK, CS_LOG_TRACE, s);
    }
  }

  /*--------------------------------------------------------------------------*/

  /* Send 'notebook_get variable' to the server */

  double
  _get_notebook_variable(int         sock,
                         const char *variable)
  {
    char *eptr = nullptr;
    double val = 0.;
    string var_s(variable);

    string buffer = "notebook_add_output " + var_s;

    string s_log = string(__func__) + ": send query for " + string(variable);
    _log(fmi2OK, CS_LOG_TRACE, s_log);

#if CS_DRY_RUN == 1

    static long _counter = 0;
    _counter += 1;
    val = (double)(_counter % 50) * 0.02;

#else

    _send_sock_str(sock, buffer.c_str());

    s_log = string(__func__) + ": waiting for " + string(variable);
    _log(fmi2OK, CS_LOG_TRACE, s_log);

    /* Received get: val */
    char *buf_rcv = _recv_sock_with_queue(_cs_socket, _control_queue, 0);

    if (strncmp(buf_rcv, "get: ", 5) != 0) {
      s_log =   string(__func__) + ": unexpected reply; " + string(buf_rcv);
      _log(fmi2Error, CS_LOG_ERROR, s_log);
      throw std::runtime_error(s_log);
    }

    val = strtod(buf_rcv + 5, &eptr);

    s_log =   string(__func__) + ": retrieved "
            + string(variable) + " (" + to_string(val) + ")";
    _log(fmi2OK, CS_LOG_TRACE, s_log);

    /* Received 0 : OK */
    buf_rcv = _recv_sock_with_queue(_cs_socket, _control_queue, 0);
    if (strcmp(buf_rcv, "0") != 0) {
      s_log =   string(__func__) + ": unexpected return code " + string(buf_rcv);
      _log(fmi2Error, CS_LOG_ERROR, s_log);
      throw std::runtime_error(s_log);
    }

#endif

    return val;
  }

  /*--------------------------------------------------------------------------*/

  /* Query serialized snapshot from the server */

  void
  _get_serialized_snapshot(int     sock)
  {
    string buffer = "snapshot_get_serialized";

    string s_log = string(__func__) + ": send query for " + buffer;
    _log(fmi2OK, CS_LOG_TRACE, s_log);

    size_t restart_size[1] = {0};

    /* Add input and output variables to end of snapshot */

    cs_variables_t *v = _cs_variables;
    size_t n_add = (v->n_input + v->n_output) * sizeof(double);

    if (_serialized_data != nullptr) {
      free(_serialized_data);
      _serialized_data = nullptr;
    }

#if CS_DRY_RUN == 1

    _serialized_size = n_add;
    _serialized_data = malloc(n_add);

#else

    _send_sock_str(sock, buffer.c_str());

    s_log = string(__func__) + ": waiting for " + buffer;
    _log(fmi2OK, CS_LOG_TRACE, s_log);

    /* Received get: val */
    char *buf_rcv = _recv_sock_with_queue(_cs_socket, _control_queue, 0);

    if (strncmp(buf_rcv, "serialized_snapshot", 19) != 0) {
      s_log =   string(__func__) + ": unexpected reply; " + string(buf_rcv);
      _log(fmi2Error, CS_LOG_ERROR, s_log);
      throw std::runtime_error(s_log);
    }

    /* Data size */
    _recv_sock(_cs_socket, (char *)restart_size, _control_queue,
               sizeof(size_t), 1);

    _serialized_size = restart_size[0] + n_add;
    _serialized_data = (char *)malloc(_serialized_size);
    assert(_serialized_data != nullptr || _serialized_size == 0);

    if (restart_size[0] > 0)
      _recv_sock(_cs_socket, _serialized_data, _control_queue,
                 1, restart_size[0]);

    s_log = string(__func__) + ": received snapshot";
    _log(fmi2OK, CS_LOG_TRACE, s_log);

    /* receive 0 */
    buf_rcv = _recv_sock_with_queue(_cs_socket, _control_queue, 0);
    if (strcmp(buf_rcv, "0") != 0) {
      string s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
      _log(fmi2Warning, CS_LOG_WARNING, s);
    }

#endif

    char *p = _serialized_data + restart_size[0];
    memcpy(p, v->input_vals, v->n_input*sizeof(double));
    p += v->n_input*sizeof(double);
    memcpy(p, v->output_vals, v->n_output*sizeof(double));
  }

  /*--------------------------------------------------------------------------*/

  /* Send 'disconnect ' to the server */

  void
  _disconnect(void)
  {
    char buffer[13] = "disconnect ";

    _send_sock_str(_cs_socket, buffer);

    _log(fmi2OK, CS_LOG_COMM, "Disconnecting the controller...");
  }

  /*--------------------------------------------------------------------------*/

  string
  _get_resu_dir() {

    string s_casename, s_run_id;
    for (int i = 0; i < _n_string_vars; i++) {
      if (strcmp(_string_names[i], "casename") == 0)
        s_casename = _string_vars[i].start;
      else if (strcmp(_string_names[i], "run_id") == 0)
        s_run_id = _string_vars[i].start;
    }

    string resu = s_casename + "/RESU/" + s_run_id;
    return resu;

  }

  /*--------------------------------------------------------------------------*/

  void
  _enter_initialization_mode(fmi2Component  component) {

    _initialization_mode = true;
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

    _log(fmi2OK, CS_LOG_TRACE, __func__);

    /* Just in case, to avoid comma separators */
    setlocale(LC_NUMERIC, "C");

    /* Set tracefile if requested here */

    const char *p = getenv("CS_FMU_COMM_TRACE");
    if (p != nullptr) {
      if (strcmp(p, "stdout") != 0 && strlen(p) > 0)
        _tracefile = fopen(p, "w");
      else
        _tracefile = stdout;
    }

    /* Initialize code_saturne calculation */
    cmd = "cd " + s_casename + " && " + s_code_saturne +
          " run --initialize --force --id='" + s_run_id + "'";

#if CS_DRY_RUN == 1

    _log(fmi2OK, CS_LOG_LAUNCH, cmd);
    _cs_socket = 444444;
    _control_queue = _queue_initialize();

#else

    _log(fmi2OK, CS_LOG_LAUNCH, string("launch: \"") + cmd + string("\""));
    _log(fmi2OK, CS_LOG_LAUNCH,
         string("output: \"") + _exec_popen(cmd.c_str()) + string("\""));

    int key = _init_client(resu);

    _log(component, _instance_name, fmi2OK, CS_LOG_COMM,
         "Waiting for connection...",
         nullptr);

    std::thread t_server(&cs_fmu::_start_client, this);
    std::thread t_cs(&cs_fmu::_start_cs, this);

    t_server.join();
    t_cs.join();

    _comm_with_saturne(key);

#endif

    /* Sleep 3 seconds just to make sure the connection is well established */
    sleep(3);

    /* Now exchange input and output values which are already declared */

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

      char *buf_rcv = _recv_sock_with_queue(_cs_socket,
                                            _control_queue,
                                            0);

      string s = string(__func__);
      if (strcmp(buf_rcv, "0") != 0) {
        s += ": unexpected return code: " + string(buf_rcv);
        _log(fmi2Warning, CS_LOG_WARNING, s);
      }
      else {
        s += ": set max time step to 9999999 for FMI control.";
        _log(fmi2OK, CS_LOG_TRACE, s);
      }
    }
  }

  /*--------------------------------------------------------------------------*/

  fmi2Real
  _get_real(fmi2ValueReference  reference)
  {
    if (_cs_variables == nullptr)
      _cs_variables = _variables_initialize();

    string s = string(__func__) + "(" + to_string(reference) + ")";
    _log(fmi2OK, CS_LOG_TRACE, s);

    int var_index = _real_variable_reference_map[reference];

    return _real_variable_values[var_index];
  }

  /*--------------------------------------------------------------------------*/

  void
  _set_real(fmi2ValueReference  reference,
            fmi2Real            value)
  {
    if (_cs_variables == nullptr)
      _cs_variables = _variables_initialize();

    string s = string(__func__) + "(" + to_string(reference)
      + ", " + to_string(value) + ")";
    _log(fmi2OK, CS_LOG_TRACE, s);

    int var_index = _real_variable_reference_map[reference];

    /* TODO: it should be possible to call _set_notebook_variable
     *       for a all relevant input variables upon initialization,
     *       just after connection. This would remove the need for the
     *       _variables_set_input function and asociated
     *       input_init arrays. We do need to make sure though that values
     *       are set before a time step, and not after.
     *
     *       For this, we should check whether "first" occurs during
     *       the initializatrion mode, in which case we could replace
     *       the check on "first" and underlying initialization tracking
     *       with a simpler check on the initialization mode. */

    if (_real_vars[var_index].causality == input) {
      _real_variable_values[var_index] = value;

      int first = _variables_set_input(_cs_variables,
                                       reference,
                                       value);
      if (first && _cs_socket >= 0) {
        _set_notebook_variable(_cs_socket,
                               _real_names[var_index],
                               value);
      }
    }
  }

  /*--------------------------------------------------------------------------*/

  fmi2Integer
  _get_integer(fmi2ValueReference  reference)
  {
    int var_index = _integer_variable_reference_map[reference];

    return _integer_variable_values[var_index];
  }

  /*--------------------------------------------------------------------------*/

  void
  _set_integer(fmi2ValueReference  reference,
               fmi2Integer         value)
  {
    int var_index = _integer_variable_reference_map[reference];

    _integer_variable_values[var_index] = value;
  }

  /*--------------------------------------------------------------------------*/

  fmi2Integer
  _get_boolean(fmi2ValueReference  reference)
  {
    int var_index = _bool_variable_reference_map[reference];

    return _bool_variable_values[var_index];
  }

  /*--------------------------------------------------------------------------*/

  void
  _set_boolean(fmi2ValueReference  reference,
               fmi2Boolean         value)
  {
    int var_index = _bool_variable_reference_map[reference];

    _bool_variable_values[var_index] = value;
  }

  /*--------------------------------------------------------------------------*/

  fmi2String
  _get_string(fmi2ValueReference  reference)
  {
    int var_index = _string_variable_reference_map[reference];

    return _string_variable_values[var_index];
  }

  /*--------------------------------------------------------------------------*/

  void
  _set_string(fmi2ValueReference  reference,
              const char*         value)
  {
    string s = string(__func__) + "(" + to_string(reference)
      + ", " + string(value) + ")";
    _log(fmi2OK, CS_LOG_TRACE, s);

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

};

/*============================================================================
 * code_saturne FMU method definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2EnterInitializationMode(fmi2Component  component)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);
  fmi2Status  status = fmi2OK;

  try {
    /* As the code generator used provides no other earlier entry point,
       global pointer to component environment (for logging) here. */

    c->_enter_initialization_mode(component);
  }
  catch (const std::runtime_error& e) {
    cout << "An exception occured: " << e.what() << "\n";
    status = fmi2Fatal;
  }

  return status;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2ExitInitializationMode(fmi2Component  component)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  c->_initialization_mode = false;

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2DoStep(fmi2Component  component,
           fmi2Real       currentComunicationTime,
           fmi2Real       stepSize,
           fmi2Boolean    noSetFmuFmuStatePriorToCurrentPoint)
{
  CS_UNUSED(currentComunicationTime);
  CS_UNUSED(noSetFmuFmuStatePriorToCurrentPoint);

  cs_fmu *c = static_cast<cs_fmu *>(component);

  fmi2Status  status = fmi2OK;

  try {
    string s = "----- " + string(__func__) + "(" + to_string(stepSize) + ")";
    c->_log(fmi2OK, CS_LOG_TRACE, s);

    cs_variables_t *v = c->_cs_variables;

    /* Update input notebook values */

    for (int i = 0; i < v->input_max; i++) {
      int id = v->input_ids[i];
      if (id < 0)
        continue;
      int var_index = c->_real_variable_reference_map[i];
      v->input_vals[id] = c->_real_variable_values[var_index];
    }

    /* Advance 1 iteration */
    c->_advance(1);

    /* Update output notebook values */

    for (int i = 0; i < v->output_max; i++) {
      int id = v->output_ids[i];
      if (id < 0)
        continue;
      int var_index = c->_real_variable_reference_map[i];
      c->_real_variable_values[var_index] = v->output_vals[id];
    }

    c->_n_iter++;
  }

  catch (const std::runtime_error& e) {
    cout << "An exception occured: " << e.what() << "\n";
    status = fmi2Fatal;
  }

  return status;
}

/*----------------------------------------------------------------------------*/

fmi2Status fmi2Terminate(fmi2Component  component)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  fmi2Status  status = fmi2OK;

  try {
    string s_casename, s_run_id;

    for (int i = 0; i < _n_string_vars; i++) {
      if (strcmp(_string_names[i], "casename") == 0)
        s_casename = _string_vars[i].start;
      else if (strcmp(_string_names[i], "run_id") == 0)
        s_run_id = _string_vars[i].start;
    }

    string resu = s_casename + "/RESU/" + s_run_id;

    c->_send_sock_str(c->_cs_socket, "max_time_step 0");
    c->_advance(0);

    c->_variables_finalize();

    _queue_finalize(&(c->_control_queue));

    /* Send disconnect to the server */
    c->_disconnect();

#if CS_DRY_RUN == 1

    _cs_socket = -1;

#else

    shutdown(c->_cs_socket, SHUT_RDWR);
    close(c->_cs_socket);

    /* Stop the computation manually */
    ofstream control_file(resu + "/control_file");
    if (control_file.is_open()) {
      control_file << 1;
      control_file.close();
    }
    else {
      string s_log = string(__func__) + ": error writing control_file";
      c->_log(fmi2Error, CS_LOG_ERROR, s_log);
      throw std::runtime_error(s_log);
    }

#endif

    if (c->_tracefile != nullptr) {
      if (c->_tracefile != stdout && c->_tracefile != stderr) {
        fclose(c->_tracefile);
        c->_tracefile = nullptr;
      }
    }
  }

  catch (const std::runtime_error& e) {
    cout << "An exception occured: " << e.what() << "\n";
    status = fmi2Fatal;
  }

  return status;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2Reset(fmi2Component  component)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  c->_log(fmi2Error, CS_LOG_ERROR, "fmi2Reset called but not handled.");

  return fmi2Error;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2GetReal(fmi2Component             component,
            const fmi2ValueReference  valueReference[],
            size_t                    numberOfValues,
            fmi2Real                  values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = c->_get_real(valueReference[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetInteger(fmi2Component             component,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               fmi2Integer               values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = c->_get_integer(valueReference[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetBoolean(fmi2Component             component,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               fmi2Boolean               values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = c->_get_boolean(valueReference[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2GetString(fmi2Component             component,
              const fmi2ValueReference  valueReference[],
              size_t                    numberOfValues,
              fmi2String                values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      values[i] = c->_get_string(valueReference[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2SetReal(fmi2Component             component,
            const fmi2ValueReference  valueReference[],
            size_t                    numberOfValues,
            const fmi2Real            values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      c->_set_real(valueReference[i], values[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetInteger(fmi2Component             component,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               const fmi2Integer         values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      c->_set_integer(valueReference[i], values[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetBoolean(fmi2Component             component,
               const fmi2ValueReference  valueReference[],
               size_t                    numberOfValues,
               const fmi2Boolean         values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      c->_set_boolean(valueReference[i], values[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

fmi2Status
fmi2SetString(fmi2Component             component,
              const fmi2ValueReference  valueReference[],
              size_t                    numberOfValues,
              const fmi2String          values[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  try {
    for (size_t i = 0; i < numberOfValues; i++)
      c->_set_string(valueReference[i], values[i]);
    return fmi2OK;
  }
  catch (...) {
    return fmi2Fatal;
  }
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2GetFMUstate(fmi2Component  component,
                fmi2FMUstate*  state)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  char s[128];
  snprintf(s, 127, "%s: called for %p; only client state is partially saved",
           __func__, (void *)*state);
  c->_log(fmi2Warning, CS_LOG_WARNING, s);

  /* Not a valid pointer to actual data here,
     just a distinct pointer per state... */

  cs_variables_t *v = c->_cs_variables;
  cs_state_t *csst = nullptr;

  if (*state == nullptr) {
    csst = (cs_state_t *)malloc(sizeof(cs_state_t));
    csst->vals = (double *)malloc(sizeof(double) * (v->n_input+v->n_output));
    fmi2FMUstate fs = (fmi2FMUstate)csst;
    c->_states[fs] = fs;

    *state = fs;
  }
  else {
    if (c->_states.count(*state) > 0)
      csst = (cs_state_t *)c->_states[*state];
    else {
      snprintf(s, 127, "%s: called for %p, not previously defined",
               __func__, (void *)*state);
      c->_log(fmi2Error, CS_LOG_ERROR, s);
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

fmi2Status
fmi2SetFMUstate(fmi2Component  component,
                fmi2FMUstate   state)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  char s[128];
  snprintf(s, 127, "%s: called for %p; only client state is partially set",
           __func__, (void *)state);
  c->_log(fmi2Warning, CS_LOG_WARNING, s);

  cs_variables_t *v = c->_cs_variables;
  cs_state_t *csst = nullptr;

  if (c->_states.count(state) > 0)
    csst = (cs_state_t *)c->_states[state];
  else {
    snprintf(s, 127, "%s: called for %p, not previously defined",
             __func__, (void *)state);
    c->_log(fmi2Error, CS_LOG_ERROR, s);
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

fmi2Status
fmi2FreeFMUstate(fmi2Component  component,
                 fmi2FMUstate*  state)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  if (*state == nullptr)
    return fmi2OK;

  char s[128];

  if (c->_states.count(*state) > 0) {
    snprintf(s, 127, "%s: called for %p\n",
             __func__, (void *)*state);
    c->_log(fmi2OK, CS_LOG_TRACE, s);
    cs_state_t *csst = (cs_state_t *)c->_states[*state];
    free(csst->vals);
    free(csst);
    c->_states.erase(*state);
    *state = nullptr;
  }
  else {
    snprintf(s, 127, "%s: called for %p, which is not present\n",
           __func__, (void *)*state);
    c->_log(fmi2Error, CS_LOG_WARNING, s);
  }

  return fmi2Warning;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2SerializedFMUstateSize(fmi2Component  component,
                           fmi2FMUstate   state,
                           size_t*        stateSize)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  c->_get_serialized_snapshot(c->_cs_socket);

  *stateSize = c->_serialized_size;
  char s[128];
  snprintf(s, 127, "%s: called for %p (%ld bytes)",
           __func__, (void *)state, (long)c->_serialized_size);
  c->_log(fmi2OK, CS_LOG_TRACE, s);

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2SerializeFMUstate(fmi2Component  component,
                      fmi2FMUstate   state,
                      fmi2Byte       serializedState[],
                      size_t         serializedStateSize)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  char s[128];

  if (serializedStateSize != c->_serialized_size) {
    snprintf(s, 127,
             "%s: called for %p\n"
             "with size %d (%d expected)",
             __func__, (void *)state,
             (int)serializedStateSize, (int)c->_serialized_size);
    c->_log(fmi2Error, CS_LOG_TRACE, s);
    return fmi2Error;
  }

  memcpy(serializedState, c->_serialized_data, c->_serialized_size);

  c->_serialized_size = 0;
  free(c->_serialized_data);
  c->_serialized_data = nullptr;

  snprintf(s, 127, "%s: called for %p", __func__, (void *)state);
  c->_log(fmi2Error, CS_LOG_TRACE, s);
  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status
fmi2DeSerializeFMUstate(fmi2Component   component,
                        const fmi2Byte  serializedState[],
                        size_t          size,
                        fmi2FMUstate*   state)
{
  CS_UNUSED(state);  // Assumed to be handled by caller only.

  cs_fmu *c = static_cast<cs_fmu *>(component);

  char s[128];
  snprintf(s, 127, "%s: called for %p\n", __func__, (void *)state);
  c->_log(fmi2OK, CS_LOG_ERROR, s);

  /* Add input and output variables to end of snapshot */

  cs_variables_t *v = c->_cs_variables;

  char *p = (char *)serializedState;
  p += size - v->n_output*sizeof(double);
  memcpy(v->output_vals, p, v->n_output*sizeof(double));
  p -= v->n_input*sizeof(double);
  memcpy(v->input_vals, p, v->n_input*sizeof(double));

  size_t restart_size = size - (v->n_output + v->n_input)*sizeof(double);

  /* Now send rest of snapshot to server */

  string buffer = "snapshot_load_serialized " + to_string(restart_size);

  string s_log = string(__func__) + ": send command : " + buffer;
  c->_log(fmi2OK, CS_LOG_TRACE, s_log);

  c->_send_sock_str(c->_cs_socket, buffer.c_str());

  c->_send_sock(c->_cs_socket, (char *)serializedState, 1, restart_size);

  return fmi2OK;
}

fmi2Status
fmi2SetDebugLogging(fmi2Component     component,
                    fmi2Boolean       loggingOn,
                    size_t            nCategories,
                    const fmi2String  categories[])
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  if (loggingOn) {
    if (nCategories == 0) {
      for (size_t i = 0; i < _n_log_categories; i++) {
        if (c->_log_active[i] == -2)
          c->_log_active[i] = 1;
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
            if (c->_log_active[k] == -2)
              c->_log_active[k] = 1;
          }
          continue;
        }
        else {
          size_t j = 0;
          while (j < _n_log_categories) {
            if (strcmp(s, _log_categories[j]) != 0) {
              c->_log_active[j] = 1;
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
      if (c->_log_active[i] != 0) {
        c->_log_active[i] = -1;
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
  cs_fmu *c = new cs_fmu(instanceName,
                         fmuResourceLocation);

  CS_UNUSED(fmuType);
  CS_UNUSED(fmuGUID);
  CS_UNUSED(visible);

  setFunctions(callbacks);

  c->_map_initialize();

  if (loggingOn) {
    for (size_t i = 0; i < _n_log_categories; i++) {
      if (c->_log_active[i] == -2) {
        c->_log_active[i] = 1;
      }
    }
  }

  /* Only use low level trace if requested. */

  bool trace_to_stdout = false;
  const char s[] = "CS_FMU_TRACE";
  if (getenv(s) != nullptr) {
    int i = atoi(getenv(s));
    if (i > 0)
      trace_to_stdout = true;
  }
  if (trace_to_stdout == false) {
    for (size_t i = 0; i < _n_log_categories; i++) {
      if (c->_log_active[i] < 0) {
        c->_log_active[i] = 0;
      }
    }
  }

  return (fmi2Component*)c;
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

void
fmi2FreeInstance(fmi2Component  component)
{
  cs_fmu *c = static_cast<cs_fmu *>(component);

  delete c;
}

fmi2Status
fmi2GetDirectionalDerivative(fmi2Component             component,
                             const fmi2ValueReference  unknownValueReferences[],
                             size_t                    numberOfUnknowns,
                             const fmi2ValueReference  knownValueReferences[],
                             fmi2Integer               numberOfKnowns,
                             fmi2Real                  knownDifferential[],
                             fmi2Real                  unknownDifferential[])
{
  CS_UNUSED(component);
  CS_UNUSED(unknownValueReferences);
  CS_UNUSED(numberOfUnknowns);
  CS_UNUSED(knownValueReferences);
  CS_UNUSED(numberOfKnowns);
  CS_UNUSED(knownDifferential);
  CS_UNUSED(unknownDifferential);

  return fmi2Error;
}

fmi2Status
fmi2SetRealInputDerivatives(fmi2Component             component,
                            const fmi2ValueReference  valueReferences[],
                            size_t                    numberOfValueReferences,
                            fmi2Integer               orders[],
                            const fmi2Real            values[])
{
  CS_UNUSED(component);
  CS_UNUSED(valueReferences);
  CS_UNUSED(numberOfValueReferences);
  CS_UNUSED(orders);
  CS_UNUSED(values);

  return fmi2Error;
}

fmi2Status
fmi2GetRealOutputDerivatives(fmi2Component             component,
                             const fmi2ValueReference  valueReference[],
                             size_t                    numberOfValues,
                             const fmi2Integer         order[],
                             fmi2Real                  values[])
{
  CS_UNUSED(component);
  CS_UNUSED(valueReference);
  CS_UNUSED(numberOfValues);
  CS_UNUSED(order);
  CS_UNUSED(values);

  return fmi2Error;
}

fmi2Status
fmi2CancelStep(fmi2Component  component)
{
  CS_UNUSED(component);

  return fmi2Error;
}

fmi2Status
fmi2GetStatus(fmi2Component         component,
              const fmi2StatusKind  kind,
              fmi2Status*           status)
{
  CS_UNUSED(component);
  CS_UNUSED(kind);
  CS_UNUSED(status);

  return fmi2Error;
}

fmi2Status
fmi2GetRealStatus(fmi2Component         component,
                  const fmi2StatusKind  kind,
                  fmi2Real*             value)
{
  CS_UNUSED(component);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

fmi2Status
fmi2GetIntegerStatus(fmi2Component         component,
                     const fmi2StatusKind  kind,
                     fmi2Integer*          value)
{
  CS_UNUSED(component);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

fmi2Status
fmi2GetBooleanStatus(fmi2Component         component,
                     const fmi2StatusKind  kind,
                     fmi2Boolean*          value)
{
  CS_UNUSED(component);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

fmi2Status
fmi2GetStringStatus(fmi2Component         component,
                    const fmi2StatusKind  kind,
                    fmi2String*           value)
{
  CS_UNUSED(component);
  CS_UNUSED(kind);
  CS_UNUSED(value);

  return fmi2Error;
}

/*----------------------------------------------------------------------------*/
