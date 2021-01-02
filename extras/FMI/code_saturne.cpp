/*============================================================================
 * Functional Mock-up Unit main methods implementation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
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

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct sockaddr_in sockaddr_in;
typedef struct sockaddr sockaddr;

struct thread_data {
  sockaddr_in address;
  int         addrlen;
};

/*============================================================================
 * Static global variables
 *============================================================================*/

static int _n_iter = 0;
static int _master_socket = -1;
static int _cs_socket = -1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

/*----------------------------------------------------------------------------*/

/* Send a message on a socket */

static void
_send_sock(int         sock,
           const char *string)
{
  if (send(sock, string, strlen(string)+1, 0) > 0) {
    cout << string << " sent" << endl;
  } else {
    cout << "Error sending " << string << endl;
    exit(1);
  }
}

/*----------------------------------------------------------------------------*/

/* Reception of a message on a socket */

static void
_recv_sock(int    socket,
           char  *buffer,
           size_t len_buffer)
{
  /* Cleaning buffer */
  memset(buffer, 0, len_buffer);

  if (recv(socket, buffer, len_buffer, 0) != -1) {
    cout << "Received: " << buffer << endl;
  } else {
    cout << "Error receiving data" << endl;
    exit(1);
  }
}

/*----------------------------------------------------------------------------*/

/* Communication with code_saturne handshake at the beginning
 * (first connection) and reception of messages from CS after */

static void
_comm_with_saturne(int key)
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
    cout << "wrong key received" << endl;
    exit(1);
  }

  _recv_sock(_cs_socket, magic_buffer, magic_string.size());

  _send_sock(_cs_socket, magic_buffer);

  cout << "Connection between server and code_saturne established !" << endl;

  /* Iteration OK */
  char buf_rcv[13];
  _recv_sock(_cs_socket, buf_rcv, 13);

  delete[] key_buffer;
  delete[] magic_buffer;
}

/*----------------------------------------------------------------------------*/

/* Configuration of the master socket (server) */

static int
_configure_server(sockaddr_in *address)
{
  int master_socket;
  int opt = 1;
  socklen_t optlen = sizeof(opt);

  /* Create a master socket */
  if ((master_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
    cout << "socket failed" << endl;
    exit(EXIT_FAILURE);
  }

  /* Set master socket to allow multiple connections */
  if (setsockopt(master_socket, SOL_SOCKET, SO_REUSEADDR, &opt, optlen) < 0) {
    cout << "setsockopt" << endl;
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
    cout << "bind failed" << endl;
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
  } else {
    cout << "Error writing control_file." << endl;
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*/

static int
_init_server(string       path,
             sockaddr_in *address)
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
  cout << "Listener on port " <<  port << endl;

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
    cout << "listen" << endl;
    exit(EXIT_FAILURE);
  }

  if ((_cs_socket = accept(_master_socket,
                           (sockaddr *)&address,
                           (socklen_t*)&addrlen)) < 0) {
    cout << "accept";
    exit(EXIT_FAILURE);
  }

  /* Inform user of socket number - used in send and receive commands */
  cout << "New connection, socket fd is " << _cs_socket;
  cout << " ip is: " << inet_ntoa(address.sin_addr);
  cout << " port: " << ntohs(address.sin_port) << endl;

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
  char buf_rcv[13];
  string buffer = "advance " + to_string(n);

  cout << "Advancing " << n << " iterations" << endl;

  _send_sock(_cs_socket, buffer.c_str());

  /* receive 0 */
  _recv_sock(_cs_socket, buf_rcv, 2);

  /* receive iteration OK */
  _recv_sock(_cs_socket, buf_rcv, 13);
}

/*----------------------------------------------------------------------------*/

/* Send 'notebook_set variable value' to the server */

static void
_set_notebook_variable(int         sock,
                       const char *variable,
                       double      value)
{
  char buf_rcv[2];
  string var_s(variable);

  string buffer = "notebook_set " + var_s + " " + to_string(value);

  cout << "Sending " << buffer << endl;

  _send_sock(sock, buffer.c_str());

  cout << "Setting " << variable << " to " << value << endl;

  /* Receive 0 */
  _recv_sock(sock, buf_rcv, 2);

  cout << "Variable " << variable << " set." << endl;
}

/*----------------------------------------------------------------------------*/

/* Send 'notebook_get variable' to the server */

static double
_get_notebook_variable(int         sock,
                       const char *variable)
{
  char *eptr = nullptr;
  char buf_rcv[20];
  double val = 0.;
  string var_s(variable);
  string buffer = "notebook_get " + var_s;

  cout << "Trying to get " << variable << endl;

  _send_sock(sock, buffer.c_str());

  cout << "Getting " << variable << endl;

  /* Received get: val */
  _recv_sock(sock, buf_rcv, 20);

  char *s = buf_rcv;
  s += 4;
  val = strtod(s, &eptr);

  /* Received 0 : OK */
  _recv_sock(sock, buf_rcv, 2);

  cout << "Variable " << variable << " retrieved." << endl;

  return val;
}

/*----------------------------------------------------------------------------*/

/* Send 'disconnect ' to the server */

static void
_disconnect(void)
{
  char buffer[13] = "disconnect ";

  _send_sock(_cs_socket, buffer);

  cout << "Disconnecting the controller..." << endl;
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
  logInfo(fmiResourceLocation);
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::init()
{
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

  /* Just in case, to avoid comma separators */
  setlocale(LC_NUMERIC, "C");

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* Initialize code_saturne calculation */
  cmd = "cd " + s_casename + " && " + s_code_saturne +
        " run --param setup.xml --initialize --id='" + s_run_id + "'";
  logInfo(exec_popen(cmd.c_str()).c_str());

  int key = _init_server(resu, &address);

  /* Accept the incoming connection */
  addrlen = sizeof(address);
  logInfo("Waiting for connections ...");

  struct thread_data data;
  data.address = address;
  data.addrlen = addrlen;

  pthread_create(&thread_server, nullptr, _start_server, static_cast<void*>(&data));
  pthread_create(&thread_cs, nullptr, _start_cs, static_cast<void*>(&resu));

  pthread_attr_destroy(&attr);
  pthread_join(thread_server, &status);
  pthread_join(thread_cs, &status);

  _comm_with_saturne(key);

  /* Sleep 3 seconds just to make sure the connection is well established */
  sleep(3);

  return fmi2OK;
}

/*----------------------------------------------------------------------------*/

fmi2Status code_saturne::doStep(double step)
{
  cout << "-----------------------" << endl;

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

  /* Send disconnect to the server */
  _disconnect();

  shutdown(_cs_socket, SHUT_RDWR);
  close(_cs_socket);

  /* Stop the computation manually */
  ofstream control_file(resu + "/control_file");
  if (control_file.is_open()) {
    control_file << 1;
    control_file.close();
  } else {
    cout << "Error writing control_file" << endl;
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
  if (realVariables[fmi2ValueReference].causality() != input) {
    return _get_notebook_variable(_cs_socket, realVariables[fmi2ValueReference].name);
  } else {
    return realVariables[fmi2ValueReference].getter(this);
  }
}

/*----------------------------------------------------------------------------*/

void code_saturne::setReal(int fmi2ValueReference, double value)
{
  if (realVariables[fmi2ValueReference].causality() == input) {
    _set_notebook_variable(_cs_socket, realVariables[fmi2ValueReference].name, value);
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
  stringVariables[fmi2ValueReference].setter(this, value);
}

/*----------------------------------------------------------------------------*/
