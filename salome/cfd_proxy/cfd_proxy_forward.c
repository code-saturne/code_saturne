//============================================================================
// Forwarding of function calls
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

// System and SALOME headers

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <calcium.h>

// Local headers

#include "cfd_proxy_defs.h"
#include "cfd_proxy_comm.h"
#include "cfd_proxy_forward.h"

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//============================================================================
//  Constants and Macros
//============================================================================

//============================================================================
// Structure definition
//============================================================================

//============================================================================
// Global variables
//============================================================================

//============================================================================
// Private function definitions
//============================================================================

//----------------------------------------------------------------------------
// Forward a CALCIUM write call
//
// parameters:
//   comm   <-> communicator
//   r      <-> container for temporary communication request data
//   type   <-- data type
//
// returns:
//   0 if request was forwarded correctly, -1 for end-of-file, -2 for error
//----------------------------------------------------------------------------

static int
_forward_calcium_write(cfd_proxy_comm_t          *comm,
                       cfd_proxy_comm_request_t  *r,
                       cfd_proxy_type_t           type)
{
  void *component = cfd_proxy_glob_component[r->comp_id];

  int time_dep = r->int_vals[0];
  int iteration = r->int_vals[1];
  int n_vals = r->int_vals[2];
  double cur_time = r->double_vals[0];
  float _cur_time = cur_time;
  char *var_name = r->string_vals;

  size_t type_size = cfd_proxy_glob_type_size[type];

  int comm_retval = 0, retval = 0;
  char *buffer = NULL;

  if (n_vals > 0) {
    CFDP_MALLOC(buffer, n_vals*type_size, char);
    comm_retval = cfd_proxy_comm_read(comm, buffer, type_size, n_vals);
    if (comm_retval != 0) {
      CFDP_FREE(buffer);
      return comm_retval;
    }
  }

  switch(type) {
  case CFD_PROXY_TYPE_int:
    retval = cp_een(component, time_dep, _cur_time, iteration,
                    var_name, n_vals, (int *)buffer);
    break;
  case CFD_PROXY_TYPE_float:
    retval = cp_ere(component, time_dep, _cur_time, iteration,
                    var_name, n_vals, (float *)buffer);
    break;
  case CFD_PROXY_TYPE_double:
    retval = cp_edb(component, time_dep, cur_time, iteration,
                    var_name, n_vals, (double *)buffer);
    break;
  default:
    assert(0);
  }

  comm_retval = cfd_proxy_comm_write_response(comm, retval, 0, 0, 0,
                                              NULL, NULL, NULL);

  if (n_vals > 0)
    CFDP_FREE(buffer);

  return comm_retval;
}

//----------------------------------------------------------------------------
// Forward a CALCIUM read call
//
// parameters:
//   comm     <-> communicator
//   r        <-> container for temporary communication request data
//   type     <-- data type
//
// returns:
//   0 if request was forwarded correctly, -2 for error
//----------------------------------------------------------------------------

static int
_forward_calcium_read(cfd_proxy_comm_t          *comm,
                      cfd_proxy_comm_request_t  *r,
                      cfd_proxy_type_t           type)
{
  void *component = cfd_proxy_glob_component[r->comp_id];

  int time_dep = r->int_vals[0];
  int iteration = r->int_vals[1];
  int n_vals_max = r->int_vals[2];
  int n_vals = 0;
  double min_time = r->double_vals[0];
  double max_time = r->double_vals[1];
  float _min_time = min_time;
  float _max_time = max_time;
  char *var_name = r->string_vals;

  int int_out[2];
  double double_out[1];

  size_t type_size = cfd_proxy_glob_type_size[type];

  int comm_retval = 0;
  int retval = 0;
  char *buffer = NULL;

  if (n_vals_max > 0)
    CFDP_MALLOC(buffer, n_vals_max*type_size, char);

  switch(type) {
  case CFD_PROXY_TYPE_int:
    retval = cp_len(component, time_dep, &_min_time, &_max_time, &iteration,
                    var_name, n_vals_max, &n_vals, (int *)buffer);
    min_time = _min_time;
    break;
  case CFD_PROXY_TYPE_float:
    retval = cp_lre(component, time_dep, &_min_time, &_max_time, &iteration,
                    var_name, n_vals_max, &n_vals, (float *)buffer);
    min_time = _min_time;
    break;
  case CFD_PROXY_TYPE_double:
    retval = cp_ldb(component, time_dep, &min_time, &max_time, &iteration,
                    var_name, n_vals_max, &n_vals, (double *)buffer);
    break;
  default:
    assert(0);
  }

  int_out[0] = iteration;
  int_out[1] = n_vals;
  double_out[0] = min_time;

  comm_retval = cfd_proxy_comm_write_response(comm, retval, 2, 1, 0,
                                              int_out, double_out, NULL);

  if (n_vals > 0 && retval == 0 && comm_retval == 0)
    comm_retval = cfd_proxy_comm_write(comm, buffer, type_size, n_vals);

  if (n_vals_max > 0)
    CFDP_FREE(buffer);

  return comm_retval;
}

//============================================================================
// Public function definitions
//============================================================================

//----------------------------------------------------------------------------
// Forward a call from the client and its response.
//
// parameters:
//   comm   <-> communicator
//   r      <-> container for temporary communication request data
//
// returns:
//   0 if request was forwarded correctly, -1 for end-of-file, -2 for error
//----------------------------------------------------------------------------

int
cfd_proxy_forward(cfd_proxy_comm_t          *comm,
                  cfd_proxy_comm_request_t  *r)
{
  int comm_retval = 0;

  comm_retval = cfd_proxy_comm_read_request(comm, r);

  if (comm_retval != 0)
    return comm_retval;

  assert(   r->comp_id >= 0
         && r->comp_id < cfd_proxy_glob_n_components);

  // Handle function depending on function name

  if (strncmp(r->func_name, "cp_e", 4) == 0) {

    if (strncmp(r->func_name, "cp_een", 32) == 0)
      comm_retval = _forward_calcium_write(comm, r, CFD_PROXY_TYPE_int);
    else if (strncmp(r->func_name, "cp_ere", 32) == 0)
      comm_retval = _forward_calcium_write(comm, r, CFD_PROXY_TYPE_float);
    else if (strncmp(r->func_name, "cp_edb", 32) == 0)
      comm_retval = _forward_calcium_write(comm, r, CFD_PROXY_TYPE_double);

  }
  else if (strncmp(r->func_name, "cp_l", 4) == 0) {

    if (strncmp(r->func_name, "cp_len", 32) == 0)
      comm_retval = _forward_calcium_read(comm, r, CFD_PROXY_TYPE_int);
    else if (strncmp(r->func_name, "cp_lre", 32) == 0)
      comm_retval = _forward_calcium_read(comm, r, CFD_PROXY_TYPE_float);
    else if (strncmp(r->func_name, "cp_ldb", 32) == 0)
      comm_retval = _forward_calcium_read(comm, r, CFD_PROXY_TYPE_double);

  }

  else if (strncmp(r->func_name, "cp_cd", 32) == 0) {

    void *component = cfd_proxy_glob_component[r->comp_id];
    char instance_name[256]; // Normally, INSTANCE_LEN = 72
    const char *string_vals[] = {instance_name};
    int retval = 0;

    memset(instance_name, '\0', sizeof(instance_name));

    retval = cp_cd(component, instance_name);

    comm_retval = cfd_proxy_comm_write_response(comm, retval, 0, 0, 1,
                                                NULL, NULL, string_vals);

  }

  else if (strncmp(r->func_name, "cp_fin", 32) == 0) {

    void *component = cfd_proxy_glob_component[r->comp_id];
    int   cont = r->int_vals[0];
    int retval = 0;

    retval = cp_fin(component, cont);

    comm_retval = cfd_proxy_comm_write_response(comm, retval, 0, 0, 0,
                                                NULL, NULL, NULL);

  }

  return comm_retval;
}

//----------------------------------------------------------------------------
// Forward all calls from the client and their responses.
//----------------------------------------------------------------------------

void
cfd_proxy_forward_all(cfd_proxy_comm_t  *comm)
{
  cfd_proxy_comm_request_t r;

  int retval = 0;

  if (comm == NULL)
    return;

  cfd_proxy_comm_init_request(&r);

  while (retval == 0) {

    retval = cfd_proxy_forward(comm, &r);

  }

  cfd_proxy_comm_finalize_request(&r);
}

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

