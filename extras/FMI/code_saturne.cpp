/*============================================================================
 * Functional Mock-up Unit main methods implementation
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
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string.h>
#include <array>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <sys/time.h>
#include <pthread.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "code_saturne.h"

/*----------------------------------------------------------------------------*/

using namespace std;

#define CS_DRY_RUN 0   // set to 1 to for simulation mode with no actual
                       // connection to code_saturne.

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

  size_t                  buf_idx[4];    /* Buffer index (0: next command,
                                            1: partial read start; 2: end,
                                            3: size) */
  char                   *buf;           /* Buffer */

} cs_control_queue_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static int _n_iter = 0;
static int _master_socket = -1;
static int _cs_socket = -1;

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
                            0,   /* comm */
                            -1,  /* launch */
                            -1}; /* trace */

/* Read_queue */

cs_control_queue_t *_cs_glob_control_queue = NULL;

/* Map for delayed notebook value settings */

std::map<int, double> _queued_notebook_send;

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
  queue->buf_idx[2] = 0;
  queue->buf_idx[3] = 32767;
  queue->buf = (char *)malloc(queue->buf_idx[3]+1);

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

/*----------------------------------------------------------------------------*/

/* Send a message on a socket */

static void
_send_sock(int          sock,
           const char  *str)
{
#if CS_DRY_RUN == 1
  return;
#endif

  size_t n = strlen(str)+1;

  if (tracefile != NULL) {
    fprintf(tracefile, "== send %d values: [", (int)n);
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

/* Reception of a message on a socket */

static void
_recv_sock(int       socket,
           char     *buffer,
           size_t    len_buffer)
{
  /* Cleaning buffer */
  memset(buffer, 0, len_buffer);

  if (tracefile != NULL) {
    fprintf(tracefile, "== receiving up to %d bytes...\n",
            (int)len_buffer);
  }

#if CS_DRY_RUN == 1
  buffer[0] = '\0';
  return;
#endif

  ssize_t ret = recv(socket, buffer, len_buffer, 0);

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
    _trace_buf(tracefile, buffer, ret);
    fprintf(tracefile, "]\n");
    fflush(tracefile);
  }
}

/*----------------------------------------------------------------------------*/

/* Reception of a message on a socket */

static char *
_recv_sock_with_queue(int                  socket,
                      cs_control_queue_t  *queue,
                      size_t               min_size)
{
  /* Read record from socket (skip if data still in queue) */

  ssize_t  n_loc_max = queue->buf_idx[3];

  size_t start_id = queue->buf_idx[2] - queue->buf_idx[1];

  if (queue->buf_idx[0] == 0) {

    /* Move previously read but not consumed data to beginning of queue */

    ssize_t n_prv = queue->buf_idx[2] - queue->buf_idx[1];
    if (n_prv > 0) {
      memmove(queue->buf, queue->buf + queue->buf_idx[1], n_prv);
    }

    while (true) {

      n_loc_max = queue->buf_idx[3] - start_id;

      if (tracefile != NULL) {
        fprintf(tracefile, "== receiving up to %d bytes...\n",
                (int)n_loc_max);
    }

#if CS_DRY_RUN == 1
      queue->buf[start_id] = '\0';
      return;
#endif

      ssize_t  ret = read(socket, queue->buf + start_id, n_loc_max);

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

      /* Check for end of command (end of line if not escaped
         by continuation character such as "/" or ","
         min_size should usually be 0, but could be set to a known size
         to handle binary data in the future. */
      size_t cut_id = start_id + ret - 1;
      bool escape = false;
      queue->buf_idx[2] = cut_id;
      while (cut_id > min_size && queue->buf[cut_id] != '\0') {
        if (queue->buf[cut_id] == ',' || queue->buf[cut_id] == '\\')
          escape = true;
        else if (queue->buf[cut_id] == '\n') {
          if (escape == false)
            break;
          else
            escape = false;
        }
        cut_id -= 1;
      }
      queue->buf_idx[1] = cut_id;
      queue->buf[cut_id] = '\0';

      start_id += ret;

      /* If we do not have a complete line, continue reading if possible */
      if (ret == n_loc_max && cut_id < 1) {
        queue->buf_idx[3] *= 2;
        queue->buf = (char *)realloc(queue->buf, queue->buf_idx[3]);
      }
      else
        break;

    }

  }  /* End of actual receive */

  /* Return pointer to next requested data to buffer,
     after updating next pointer in buffer */

  char *next_p = queue->buf + queue->buf_idx[0];
  size_t n_cur = strlen(next_p);

  queue->buf_idx[0] += n_cur;

  if (queue->buf_idx[0] >= queue->buf_idx[1])
    queue->buf_idx[0] = 0;

  return next_p;
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
  _recv_sock(_cs_socket, key_buffer, len_key);

  if (strncmp(key_s.c_str(), key_buffer, len_key) != 0) {
    _cs_log(fmi2Fatal, CS_LOG_ERROR,
            "wrong key received (socket handshake)");
    exit(EXIT_FAILURE);
  }

  _recv_sock(_cs_socket, magic_buffer, magic_string.size());

  _send_sock(_cs_socket, magic_buffer);

  _cs_log(fmi2OK, CS_LOG_ALL,
          "Connection between FMU client and code_saturne established");

  /* Iteration OK */
  char buf_rcv[13];
  _recv_sock(_cs_socket, buf_rcv, 13);

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

/*----------------------------------------------------------------------------*/

static void
_advance(int   n)
{
  string buffer = "advance " + to_string(n);

  string s = "Advancing " + to_string(n) + " iterations";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

#if CS_DRY_RUN == 1
  return;
#endif

  _send_sock(_cs_socket, buffer.c_str());

  /* receive 0 */
  char *buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);
  if (strcmp(buf_rcv, "0") != 0) {
    s = string(__func__) + ": unexpected return code: " + string(buf_rcv);
    _cs_log(fmi2Warning, CS_LOG_WARNING, s);
  }

  /* receive iteration OK */
  buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);
  if (strcmp(buf_rcv, "Iteration OK") != 0) {
    s =   string(__func__) + ": expected \"Iteration OK\", not: "
        + string(buf_rcv);
    _cs_log(fmi2Warning, CS_LOG_WARNING, s);
  }
}

/*----------------------------------------------------------------------------*/

/* Queue notebook variable to send to the server */

static void
_queue_notebook_variable(int       fmi2ValueReference,
                         double    value)
{
  _queued_notebook_send[fmi2ValueReference] = value;
}

/*----------------------------------------------------------------------------*/

/* Send 'notebook_set variable value' to the server */

static void
_set_notebook_variable(int          sock,
                       const char  *variable,
                       double       value)
{
  string var_s(variable);

  string buffer = "notebook_set " + var_s + " " + to_string(value);

  string s = string(__func__) + ": sending " + string(variable);
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  _send_sock(sock, buffer.c_str());

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
  string buffer = "notebook_get " + var_s;

  string s_log = string(__func__) + ": send query for " + string(variable);
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

#if CS_DRY_RUN == 1

  static long _counter = 0;
  _counter += 1;
  val = (double)(_counter % 50) * 0.02;

#else

  _send_sock(sock, buffer.c_str());

  s_log = string(__func__) + ": waiting for " + string(variable);
  _cs_log(fmi2OK, CS_LOG_TRACE, s_log);

  /* Received get: val */
  char *buf_rcv = _recv_sock_with_queue(_cs_socket, _cs_glob_control_queue, 0);

  char *s = buf_rcv;
  s += 4;
  val = strtod(s, &eptr);

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

  _send_sock(_cs_socket, buffer);

  _cs_log(fmi2Error, CS_LOG_TRACE, "Disconnecting the controller...");
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

  _instance_name = malloc(sizeof(id) + 1);
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
        " run --param setup.xml --initialize --id='" + s_run_id + "'";

#if CS_DRY_RUN == 1

  _cs_log(fmi2OK, CS_LOG_LAUNCH, cmd);

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

  /* Now send queued values */

  if (_queued_notebook_send.size() > 0) {
    for (const auto& n : _queued_notebook_send) {
      _set_notebook_variable(_cs_socket,
                             realVariables[n.first].name,
                             n.second);
    }
    _queued_notebook_send.clear();
  }

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::doStep(double step)
{
  string s = "----- " + string(__func__) + "(" + to_string(step) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  /* Advance 1 iterations */
  _advance(1);

  _n_iter++;

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::terminate()
{
  string s_casename(casename);
  string s_run_id(run_id);
  string resu = s_casename + "/RESU/" + s_run_id;

  _queue_finalize(&_cs_glob_control_queue);

  /* Send disconnect to the server */
  _disconnect();

  _cs_log(fmi2OK, CS_LOG_ALL, "Disconnecting the controller...");

#if CS_DRY_RUN != 1

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
  string s = string(__func__) + "(" + to_string(fmi2ValueReference) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  if (realVariables[fmi2ValueReference].causality() != input) {
    return _get_notebook_variable(_cs_socket, realVariables[fmi2ValueReference].name);
  }
  else {
    return realVariables[fmi2ValueReference].getter(this);
  }
}

/*----------------------------------------------------------------------------*/

void code_saturne::setReal(int fmi2ValueReference, double value)
{
  string s = string(__func__) + "(" + to_string(fmi2ValueReference)
    + ", " + to_string(value) + ")";
  _cs_log(fmi2OK, CS_LOG_TRACE, s);

  if (realVariables[fmi2ValueReference].causality() == input) {
    if (_cs_socket >= 0)
      _set_notebook_variable(_cs_socket,
                             realVariables[fmi2ValueReference].name,
                             value);
    else
      _queue_notebook_variable(fmi2ValueReference, value);
    realVariables[fmi2ValueReference].setter(this, value);
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
