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
 * Basic CALCIUM-mappable functions for code coupling using SALOME's YACS
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_proxy_comm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_calcium.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Maximum string Lengths (based on CALCIUM's limits) */

#define CS_CALCIUM_INSTANCE_LEN 72
#define CS_CALCIUM_VARIABLE_LEN 144

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* CALCIUM datatypes */

typedef enum {

  CALCIUM_integer,      /* Integer values */
  CALCIUM_real,         /* Floating-point values */
  CALCIUM_double,       /* Double-precision floating-point values */
  CALCIUM_complex,      /* Complex values (not used by Code_Saturne) */
  CALCIUM_string,       /* character string */
  CALCIUM_logical       /* Logical values */

} cs_calcium_datatype_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Use communication with proxy ? */

static int _cs_glob_calcium_comm_proxy = 0;

/* Verbosity (none if -1, headers if 0,
   headers + n first and last elements if > 0 */

static int _cs_glob_calcium_n_echo = -1;

/* Pointer of type Superv_Component_i* to the supervisable SALOME component */

static void *_cs_glob_calcium_component[8] = {NULL, NULL, NULL, NULL,
                                              NULL, NULL, NULL, NULL};

/* Map from enumerated values to SALOME's Calcium API defined values */

static int _cs_calcium_timedep_map[3] = {40,   /* CP_TEMPS      = 40 */
                                         41,   /* CP_ITERATION  = 41 */
                                         42};  /* CP_SEQUENTIAL = 42 */

static int _cs_calcium_continuation_map[2] = {20,   /* CP_CONT   = 20 */
                                              21};  /* CP_ARRET  = 21 */

/* Calcium datatype names */

static const char *cs_calcium_datatype_name[] = {"integer", "real", "double",
                                                 "complex", "string",
                                                 "logical"};
static const char *cs_calcium_timedep_name[] = {"T", "I", "S"};

/* Function pointers */

static cs_calcium_connect_t         *_cs_glob_calcium_connect = NULL;
static cs_calcium_disconnect_t      *_cs_glob_calcium_disconnect = NULL;
static cs_calcium_read_int_t        *_cs_glob_calcium_read_int = NULL;
static cs_calcium_read_float_t      *_cs_glob_calcium_read_float = NULL;
static cs_calcium_read_double_t     *_cs_glob_calcium_read_double = NULL;
static cs_calcium_write_int_t       *_cs_glob_calcium_write_int = NULL;
static cs_calcium_write_float_t     *_cs_glob_calcium_write_float = NULL;
static cs_calcium_write_double_t    *_cs_glob_calcium_write_double = NULL;

/*=============================================================================
 * Local function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print (part of) an array
 *
 * parameters:
 *   datatype <-- section data type
 *   n_echo   <-- number of beginning and end values to print
 *   n_val    <-- number of values in array
 *   val      <-- array values
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_body(cs_calcium_datatype_t   datatype,
                   int                     n_echo,
                   int                     n_val,
                   const void             *val)
{
  int start_id, end_id, id;

  if (n_val == 0) return;

  assert(val != NULL);

  start_id = 0;

  if (n_echo * 2 < n_val) {

    end_id = n_echo;
    bft_printf(_("    %d first and last elements:\n"), n_echo);

  }
  else {

    end_id = n_val;
    bft_printf(_("    elements:\n"));

  }

  do {

    switch(datatype) {

    case CALCIUM_integer:
      {
        const int *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : %12d\n", id + 1, *(_val + id));
      }
      break;

    case CALCIUM_real:
      {
        const float *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : %12.5e\n", id + 1,
                     (double)(*(_val + id)));
      }
      break;

    case CALCIUM_double:
      {
        const double *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : %14.7e\n", id + 1, *(_val + id));
      }
      break;

    case CALCIUM_complex:
      {
        const float *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : (%12.5e, %12.5e)\n", id + 1,
                     (double)(*(_val + 2*id)), (double)(*(_val + 2*id + 1)));

      }
      break;

    case CALCIUM_string:
      {
        const char *const *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : '%s\n", id + 1, _val[id]);
      }
      break;

    default:

      assert(0);

    } /* End of switch on element type */

    if (end_id < n_val) {

      bft_printf(_("    ..........   ............\n"));

      start_id = n_val - n_echo;
      end_id = n_val;

    }
    else {

      assert(end_id == n_val);
      end_id = n_val + 1;

    }

  } while (end_id <= n_val);

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print message indicating that we are ready to read data
 *
 * parameters:
 *   comp_id    <-- coupled component id
 *   var_name   <-- variable name
 *   time_dep   <-- time dependency
 *   min_time   <-- time interval low
 *   max_time   <-- time interval high
 *   iteration  <-- iteration step
 *   datatype   <-- section data type
 *   n_max_vals <-- maximum number of values to read
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_pre_read(int                     comp_id,
                       const char             *var_name,
                       cs_calcium_timedep_t    time_dep,
                       double                  min_time,
                       double                  max_time,
                       int                     iteration,
                       cs_calcium_datatype_t   datatype,
                       int                     n_max_vals)
{
  if (_cs_glob_calcium_n_echo < 0)
    return;

  assert(var_name != NULL);

  if (_cs_glob_calcium_component[comp_id] != NULL)
    bft_printf(_("\nComponent %d [%p]:\n"),
               comp_id, _cs_glob_calcium_component[comp_id]);
  else
    bft_printf(_("\nComponent %d:\n"), comp_id);

  bft_printf(_("Reading up to %d values of type %s, time_dependency %s\n"
               "              (min/max time %f/%f, iteration %d) ..."),
             n_max_vals, cs_calcium_datatype_name[datatype],
             cs_calcium_timedep_name[time_dep],
             min_time, max_time, iteration);
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print message indicating that we are finished reading data, and optionnaly
 * print a part of the corresponding data.
 *
 * parameters:
 *   min_time  <-- time interval low
 *   iteration <-- iteration step
 *   datatype  <-- section data type
 *   n_val     <-- number of values in array
 *   val       <-- array values
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_post_read(double                 min_time,
                        int                    iteration,
                        cs_calcium_datatype_t  datatype,
                        int                    n_val,
                        const void            *val)
{
  if (_cs_glob_calcium_n_echo < 0)
    return;

  bft_printf(_("[ok]\n"
               "Read          %d values (min time %f, iteration %d).\n"),
             n_val, min_time, iteration);

  _calcium_echo_body(datatype, _cs_glob_calcium_n_echo, n_val, val);
}

/*----------------------------------------------------------------------------
 * Print message indicating that we are ready to write data
 *
 * parameters:
 *   comp_id   <-- coupled component id
 *   var_name  <-- variable name
 *   time_dep  <-- time dependency
 *   cur_time  <-- current time
 *   iteration <-- iteration step
 *   datatype  <-- section data type
 *   n_vals    <-- number of values to read
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_pre_write(int                     comp_id,
                        const char             *var_name,
                        cs_calcium_timedep_t    time_dep,
                        double                  cur_time,
                        int                     iteration,
                        cs_calcium_datatype_t   datatype,
                        int                     n_vals)
{
  if (_cs_glob_calcium_n_echo < 0)
    return;

  assert(var_name != NULL);

  if (_cs_glob_calcium_component[comp_id] != NULL)
    bft_printf(_("\nComponent %d [%p]:\n"),
               comp_id, _cs_glob_calcium_component[comp_id]);
  else
    bft_printf(_("\nComponent %d:\n"), comp_id);

  bft_printf(_("Writing %d values of type %s, time_dependency %s\n"
               "        (time %f, iteration %d) ..."),
             n_vals, cs_calcium_datatype_name[datatype],
             cs_calcium_timedep_name[time_dep], cur_time, iteration);
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print message indicating that we are finished writing data, and optionnaly
 * print a part of the corresponding data.
 *
 * parameters:
 *   datatype <-- section data type
 *   n_val    <-- number of values in array
 *   val      <-- array values
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_post_write(cs_calcium_datatype_t  datatype,
                         int                    n_val,
                         const void            *val)
{
  if (_cs_glob_calcium_n_echo < 0)
    return;

  bft_printf(_("[ok]\n"));

  _calcium_echo_body(datatype, _cs_glob_calcium_n_echo, n_val, val);
}

/*----------------------------------------------------------------------------
 * Call CALCIUM-type connection through inter-process communication
 * with a proxy.
 *
 * parameters:
 *   comp_id <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   s       --> name of calling instance
 *               (CS_CALCIUM_INSTANCE_LEN chars max)
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

static int
_proxy_comm_connect(int   comp_id,
                    char *s)
{
  int  retval = 0;

  char *string_out[] = {s};

  cs_proxy_comm_write_request("cp_cd", comp_id,
                              0, 0, 0, NULL, NULL, NULL);

  retval = cs_proxy_comm_read_response(0, 0, 1, NULL, NULL, string_out);

  return retval;
}

/*----------------------------------------------------------------------------
 * Call CALCIUM-type disconnection through inter-process communication
 * with a proxy.
 *
 * parameters:
 *   comp_id <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   cont    <-- continuation directive
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

static int
_proxy_comm_disconnect(int  comp_id,
                       int  cont)
{
  int  retval = 0;

  const int int_in[] = {cont};

  cs_proxy_comm_write_request("cp_fin", comp_id,
                              1, 0, 0, int_in, NULL, NULL);

  retval = cs_proxy_comm_read_response(0, 0, 0,
                                       NULL, NULL, NULL);

  return retval;
}

/*----------------------------------------------------------------------------
 * Base function for reading values, relaying the CALCIUM-type call through
 * inter-process communication with a proxy.
 *
 * parameters:
 *   func_name  <-- associated function name
 *   blocking   <-- use blocking variant (1) or non-blocking (0)
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   min_time   <-> lower bound of read interval
 *   max_time   <-- upper bound of read interval
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val_max  <-- maximum number of values to read
 *   n_val_read <-- maximum number of values to read
 *   val        --> values read
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

static int
_proxy_comm_read_any(const char  *func_name,
                     int          blocking,
                     int          comp_id,
                     int          time_dep,
                     double      *min_time,
                     double      *max_time,
                     int         *iteration,
                     const char  *var_name,
                     int          n_val_max,
                     int         *n_val_read,
                     size_t       type_size,
                     void        *val)
{
  int  retval = 0;

  int int_out[2];
  double double_out[1];

  const int int_in[] = {time_dep, *iteration, n_val_max};
  const double double_in[] = {*min_time, *max_time};
  const char *string_in[] = {var_name};

  cs_proxy_comm_write_request(func_name, comp_id,
                              3, 2, 1,
                              int_in, double_in, string_in);

  retval = cs_proxy_comm_read_response(2, 1, 0,
                                       int_out, double_out, NULL);

  if (blocking == 0 && retval == 13) /* 13: CPATTENTE */
    return retval;

  *min_time = double_out[0];
  *iteration = int_out[0];
  *n_val_read = int_out[1];

  if (*n_val_read > 0)
    cs_proxy_comm_read(val, type_size, *n_val_read);

  return retval;
}

/*----------------------------------------------------------------------------
 * Base function for writing values, relaying the CALCIUM-type call through
 * inter-process communication with a proxy.
 *
 * parameters:
 *   func_name  <-- associated function name
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   cur_time   <-- current time
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val      <-- maximum number of values to write
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

static int
_proxy_comm_write_any(const char  *func_name,
                      int          comp_id,
                      int          time_dep,
                      double       cur_time,
                      int          iteration,
                      const char  *var_name,
                      int          n_val,
                      size_t       type_size,
                      void        *val)
{
  int  retval = 0;

  const int int_in[] = {time_dep, iteration, n_val};
  const double double_in[] = {cur_time};
  const char *string_in[] = {var_name};

  cs_proxy_comm_write_request(func_name, comp_id,
                              3, 1, 1,
                              int_in, double_in, string_in);

  if (n_val > 0)
    cs_proxy_comm_write(val, type_size, n_val);

  retval = cs_proxy_comm_read_response(0, 0, 0, NULL, NULL, NULL);

  return retval;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Connection
 *
 * parameters:
 *   comp_id <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   s       --> name of calling instance
 *               (CS_CALCIUM_INSTANCE_LEN chars max)
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_connect(int   comp_id,
                   char *s)
{
  int retval = 0;

  if (_cs_glob_calcium_connect != NULL)
    retval = _cs_glob_calcium_connect(_cs_glob_calcium_component[comp_id],
                                      s);

  else if (_cs_glob_calcium_comm_proxy)
    retval = _proxy_comm_connect(comp_id, s);

  return retval;
}

/*----------------------------------------------------------------------------
 * Disconnection
 *
 * parameters:
 *   comp_id <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   cont    --> continuation directive (continue with last values or stop)
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_disconnect(int                       comp_id,
                      cs_calcium_continuation_t cont)
{
  int _cont = _cs_calcium_continuation_map[cont];
  int retval = 0;

  if (_cs_glob_calcium_connect != NULL)
    retval = _cs_glob_calcium_disconnect(_cs_glob_calcium_component[comp_id],
                                         _cont);

  else if (_cs_glob_calcium_comm_proxy)
    retval = _proxy_comm_disconnect(comp_id, _cont);

  return retval;
}

/*----------------------------------------------------------------------------
 * Read integer values, blocking until they are available.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   min_time   <-> lower bound of read interval
 *   max_time   <-- upper bound of read interval
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val_max  <-- maximum number of values to read
 *   n_val_read <-- maximum number of values to read
 *   val        --> values read
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_read_int(int                    comp_id,
                    cs_calcium_timedep_t   time_dep,
                    double                *min_time,
                    double                *max_time,
                    int                   *iteration,
                    const char            *var_name,
                    int                    n_val_max,
                    int                   *n_val_read,
                    int                    val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  _time_dep = _cs_calcium_timedep_map[time_dep];
  float _min_time = *min_time;
  float _max_time = *max_time;
  int  _retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_read(comp_id,
                         _var_name, time_dep, *min_time, *max_time, *iteration,
                         CALCIUM_integer, n_val_max);

  if (_cs_glob_calcium_read_int != NULL) {
    _retval = _cs_glob_calcium_read_int(_cs_glob_calcium_component[comp_id],
                                        _time_dep,
                                        &_min_time,
                                        &_max_time,
                                        iteration,
                                        _var_name,
                                        n_val_max,
                                        n_val_read,
                                        val);
    *min_time = _min_time;
    *max_time = _max_time;
  }

  else if (_cs_glob_calcium_comm_proxy)
    _retval = _proxy_comm_read_any("cp_len", 1, comp_id,
                                   _time_dep, min_time, max_time, iteration,
                                   _var_name, n_val_max, n_val_read,
                                   sizeof(int), val);

  _calcium_echo_post_read(*min_time, *iteration,
                          CALCIUM_integer, *n_val_read, val);

  return _retval;
}

/*----------------------------------------------------------------------------
 * Read single-precision floating-point values,
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   min_time   <-> lower bound of read interval
 *   max_time   <-- upper bound of read interval
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val_max  <-- maximum number of values to read
 *   n_val_read <-- maximum number of values to read
 *   val        --> values read
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_read_float(int                    comp_id,
                      cs_calcium_timedep_t   time_dep,
                      double                *min_time,
                      double                *max_time,
                      int                   *iteration,
                      const char            *var_name,
                      int                    n_val_max,
                      int                   *n_val_read,
                      float                  val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  _time_dep = _cs_calcium_timedep_map[time_dep];
  float _min_time = *min_time;
  float _max_time = *max_time;
  int  _retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_read(comp_id,
                         _var_name, time_dep, *min_time, *max_time, *iteration,
                         CALCIUM_real, n_val_max);

  if (_cs_glob_calcium_read_float != NULL) {
    _retval = _cs_glob_calcium_read_float(_cs_glob_calcium_component[comp_id],
                                          _time_dep,
                                          &_min_time,
                                          &_max_time,
                                          iteration,
                                          _var_name,
                                          n_val_max,
                                          n_val_read,
                                          val);
    *min_time = _min_time;
    *max_time = _max_time;
  }

  else if (_cs_glob_calcium_comm_proxy)
    _retval = _proxy_comm_read_any("cp_lre", 1, comp_id,
                                   _time_dep, min_time, max_time, iteration,
                                   _var_name, n_val_max, n_val_read,
                                   sizeof(float), val);

  _calcium_echo_post_read(*min_time, *iteration,
                          CALCIUM_real, *n_val_read, val);

  return _retval;
}

/*----------------------------------------------------------------------------
 * Read double-precision floating-point values,
 * blocking until they are available.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   min_time   <-> lower bound of read interval
 *   max_time   <-- upper bound of read interval
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val_max  <-- maximum number of values to read
 *   n_val_read <-- maximum number of values to read
 *   val        --> values read
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_read_double(int                    comp_id,
                       cs_calcium_timedep_t   time_dep,
                       double                *min_time,
                       double                *max_time,
                       int                   *iteration,
                       const char            *var_name,
                       int                    n_val_max,
                       int                   *n_val_read,
                       double                 val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  _time_dep = _cs_calcium_timedep_map[time_dep];
  int  _retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_read(comp_id,
                         _var_name, time_dep, *min_time, *max_time, *iteration,
                         CALCIUM_double, n_val_max);

  if (_cs_glob_calcium_read_double != NULL)
    _retval = _cs_glob_calcium_read_double(_cs_glob_calcium_component[comp_id],
                                           _time_dep,
                                           min_time,
                                           max_time,
                                           iteration,
                                           _var_name,
                                           n_val_max,
                                           n_val_read,
                                           val);

  else if (_cs_glob_calcium_comm_proxy)
    _retval = _proxy_comm_read_any("cp_ldb", 1, comp_id,
                                   _time_dep, min_time, max_time, iteration,
                                   _var_name, n_val_max, n_val_read,
                                   sizeof(double), val);

  _calcium_echo_post_read(*min_time, *iteration,
                          CALCIUM_double, *n_val_read, val);

  return _retval;
}

/*----------------------------------------------------------------------------
 * Write integer values.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   cur_time   <-- current time
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_int(int                    comp_id,
                     cs_calcium_timedep_t   time_dep,
                     double                 cur_time,
                     int                    iteration,
                     const char            *var_name,
                     int                    n_val,
                     const int              val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  _time_dep = _cs_calcium_timedep_map[time_dep];
  int  *_val = NULL;

  int retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_write(comp_id,
                          _var_name, time_dep, cur_time, iteration,
                          CALCIUM_integer, n_val);

  BFT_MALLOC(_val, n_val, int);
  memcpy(_val, val, n_val * sizeof(int));

  if (_cs_glob_calcium_write_int != NULL)
    retval = _cs_glob_calcium_write_int(_cs_glob_calcium_component[comp_id],
                                        _time_dep,
                                        cur_time,
                                        iteration,
                                        _var_name,
                                        n_val,
                                        _val);

  else if (_cs_glob_calcium_comm_proxy)
    _proxy_comm_write_any("cp_een", comp_id, _time_dep,
                          cur_time, iteration, _var_name,
                          n_val, sizeof(int), _val);

  BFT_FREE(_val);

  _calcium_echo_post_write(CALCIUM_integer, n_val, val);

  return retval;
}

/*----------------------------------------------------------------------------
 * Write float values.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   cur_time   <-- current time
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_float(int                    comp_id,
                       cs_calcium_timedep_t   time_dep,
                       double                 cur_time,
                       int                    iteration,
                       const char            *var_name,
                       int                    n_val,
                       const float            val[])
{
  char  _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int   _time_dep = _cs_calcium_timedep_map[time_dep];
  float *_val = NULL;

  int retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_write(comp_id,
                          _var_name, time_dep,  cur_time, iteration,
                          CALCIUM_real, n_val);

  BFT_MALLOC(_val, n_val, float);
  memcpy(_val, val, n_val * sizeof(float));

  if (_cs_glob_calcium_write_float != NULL)
    retval = _cs_glob_calcium_write_float(_cs_glob_calcium_component[comp_id],
                                          _time_dep,
                                          cur_time,
                                          iteration,
                                          _var_name,
                                          n_val,
                                          _val);

  else if (_cs_glob_calcium_comm_proxy)
    _proxy_comm_write_any("cp_ere", comp_id, _time_dep,
                          cur_time, iteration, _var_name,
                          n_val, sizeof(float), _val);

  BFT_FREE(_val);

  _calcium_echo_post_write(CALCIUM_real, n_val, val);

  return retval;
}

/*----------------------------------------------------------------------------
 * Write double-precision float values.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   cur_time   <-- current time
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_double(int                    comp_id,
                        cs_calcium_timedep_t   time_dep,
                        double                 cur_time,
                        int                    iteration,
                        const char            *var_name,
                        int                    n_val,
                        const double           val[])
{
  char   _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int    _time_dep = _cs_calcium_timedep_map[time_dep];
  double *_val = NULL;

  int retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_write(comp_id,
                          _var_name, time_dep,  cur_time, iteration,
                          CALCIUM_double, n_val);

  BFT_MALLOC(_val, n_val, double);
  memcpy(_val, val, n_val * sizeof(double));

  if (_cs_glob_calcium_write_double != NULL)
    retval = _cs_glob_calcium_write_double(_cs_glob_calcium_component[comp_id],
                                           _time_dep,
                                           cur_time,
                                           iteration,
                                           _var_name,
                                           n_val,
                                           _val);

  else if (_cs_glob_calcium_comm_proxy)
    _proxy_comm_write_any("cp_edb", comp_id, _time_dep,
                          cur_time, iteration, _var_name,
                          n_val, sizeof(double), _val);

  BFT_FREE(_val);

  _calcium_echo_post_write(CALCIUM_double, n_val, val);

  return retval;
}

/*----------------------------------------------------------------------------
 * Assign a component and its id
 *
 * parameters:
 *   comp_id <-- id of component (0 to n-1, Code_Saturne local)
 *   comp    <-- pointer to component
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_component(int    comp_id,
                         void  *comp)
{
  assert(comp_id > -1 && comp_id < 8); /* Current limit, easily made dynamic */

  _cs_glob_calcium_component[comp_id] = comp;
}

/*----------------------------------------------------------------------------
 * Set proxy IPC communication mode
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_comm_proxy(void)
{
  _cs_glob_calcium_comm_proxy = 1;
}

/*----------------------------------------------------------------------------
 * Set the CALCIUM-mappable function's verbosity
 *
 * parameters:
 *   n_echo <-- verbosity (none if -1, headers if 0,
 *              headers + n first and last elements if > 0.
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_verbosity(int  n_echo)
{
  _cs_glob_calcium_n_echo = n_echo;
}

/*----------------------------------------------------------------------------
 * Set connect and disconnect functions
 *
 * parameters:
 *   cp_cd_func  <-- pointer to cp_cd function or equivalent
 *   cp_fin_func <-- pointer to cp_fin function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_connection_funcs(cs_calcium_connect_t     *cp_cd_func,
                                cs_calcium_disconnect_t  *cp_fin_func)
{
  _cs_glob_calcium_connect = cp_cd_func;
  _cs_glob_calcium_disconnect = cp_fin_func;
}

/*----------------------------------------------------------------------------
 * Set read, non-blocking read, and write functions for integers
 *
 * parameters:
 *   cp_len_func  <-- pointer to cp_len function or equivalent
 *   cp_een_func  <-- pointer to cp_een function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_int_rw_funcs(cs_calcium_read_int_t     *cp_len_func,
                            cs_calcium_write_int_t    *cp_een_func)
{
  _cs_glob_calcium_read_int = cp_len_func;
  _cs_glob_calcium_write_int = cp_een_func;
}

/*----------------------------------------------------------------------------
 * Set read, non-blocking read, and write functions for floats
 *
 * parameters:
 *   cp_lre_func  <-- pointer to cp_lre function or equivalent
 *   cp_ere_func  <-- pointer to cp_ere function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_float_rw_funcs(cs_calcium_read_float_t     *cp_lre_func,
                              cs_calcium_write_float_t    *cp_ere_func)
{
  _cs_glob_calcium_read_float = cp_lre_func;
  _cs_glob_calcium_write_float = cp_ere_func;
}

/*----------------------------------------------------------------------------
 * Set read, non-blocking read, and write functions for doubles
 *
 * parameters:
 *   cp_ldb_func  <-- pointer to cp_ldb function or equivalent
 *   cp_edb_func  <-- pointer to cp_edb function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_double_rw_funcs(cs_calcium_read_double_t     *cp_ldb_func,
                               cs_calcium_write_double_t    *cp_edb_func)
{
  _cs_glob_calcium_read_double = cp_ldb_func;
  _cs_glob_calcium_write_double = cp_edb_func;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
