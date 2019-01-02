/*============================================================================
 * Basic CALCIUM-mappable functions for code coupling using SALOME's YACS
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_calcium.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef int
(cs_calcium_yacsinit_t)(void);

typedef int
(cs_calcium_connect_t)(void  *component,
                       char  *s);

typedef int
(cs_calcium_disconnect_t)(void  *component,
                          int    cont);

typedef int
(cs_calcium_read_int_t)(void    *component,
                        int      time_dep,
                        float   *min_time,
                        float   *max_time,
                        int     *iteration,
                        char    *var_name,
                        int      n_val_max,
                        int     *n_val_read,
                        int      val[]);

typedef int
(cs_calcium_read_float_t)(void    *component,
                          int      time_dep,
                          float   *min_time,
                          float   *max_time,
                          int     *iteration,
                          char    *var_name,
                          int      n_val_max,
                          int     *n_val_read,
                          float    val[]);

typedef int
(cs_calcium_read_double_t)(void    *component,
                           int      time_dep,
                           double  *min_time,
                           double  *max_time,
                           int     *iteration,
                           char    *var_name,
                           int      n_val_max,
                           int     *n_val_read,
                           double   val[]);

typedef int
(cs_calcium_write_int_t)(void    *component,
                         int      time_dep,
                         float    cur_time,
                         int      iteration,
                         char    *var_name,
                         int      n_val,
                         int      val[]);

typedef int
(cs_calcium_write_float_t)(void    *component,
                           int      time_dep,
                           float    cur_time,
                           int      iteration,
                           char    *var_name,
                           int      n_val,
                           float    val[]);

typedef int
(cs_calcium_write_double_t)(void    *component,
                            int      time_dep,
                            double   cur_time,
                            int      iteration,
                            char    *var_name,
                            int      n_val,
                            double   val[]);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Verbosity (none if -1, headers if 0,
   headers + n first and last elements if > 0 */

static int _cs_calcium_n_echo = -1;

/* Pointer of type Superv_Component_i* to the supervisable SALOME component */

static void *_cs_calcium_component[8] = {NULL, NULL, NULL, NULL,
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

/* YACS dynamic library, initialization, and specific error handling */

static void  *_cs_calcium_yacslib = NULL;

static cs_calcium_yacsinit_t  *_cs_calcium_yacsinit = NULL;

/* Calcium function pointers */

static cs_calcium_connect_t         *_cs_calcium_connect = NULL;
static cs_calcium_disconnect_t      *_cs_calcium_disconnect = NULL;
static cs_calcium_read_int_t        *_cs_calcium_read_int = NULL;
static cs_calcium_read_float_t      *_cs_calcium_read_float = NULL;
static cs_calcium_read_double_t     *_cs_calcium_read_double = NULL;
static cs_calcium_write_int_t       *_cs_calcium_write_int = NULL;
static cs_calcium_write_float_t     *_cs_calcium_write_float = NULL;
static cs_calcium_write_double_t    *_cs_calcium_write_double = NULL;

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
  if (_cs_calcium_n_echo < 0)
    return;

  assert(var_name != NULL);

  if (_cs_calcium_component[comp_id] != NULL)
    bft_printf(_("\nComponent %d [%p], port %s:\n"),
               comp_id, _cs_calcium_component[comp_id], var_name);
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
  if (_cs_calcium_n_echo < 0)
    return;

  bft_printf(_("[ok]\n"
               "Read          %d values (min time %f, iteration %d).\n"),
             n_val, min_time, iteration);

  _calcium_echo_body(datatype, _cs_calcium_n_echo, n_val, val);
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
  if (_cs_calcium_n_echo < 0)
    return;

  assert(var_name != NULL);

  if (_cs_calcium_component[comp_id] != NULL)
    bft_printf(_("\nComponent %d [%p], port %s:\n"),
               comp_id, _cs_calcium_component[comp_id], var_name);
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
  if (_cs_calcium_n_echo < 0)
    return;

  bft_printf(_("[ok]\n"));

  _calcium_echo_body(datatype, _cs_calcium_n_echo, n_val, val);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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

  if (_cs_calcium_connect != NULL)
    retval = _cs_calcium_connect(_cs_calcium_component[comp_id],
                                 s);

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

  if (_cs_calcium_disconnect != NULL)
    retval = _cs_calcium_disconnect(_cs_calcium_component[comp_id],
                                    _cont);

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

  if (_cs_calcium_read_int != NULL) {
    _retval = _cs_calcium_read_int(_cs_calcium_component[comp_id],
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

  if (_cs_calcium_read_float != NULL) {
    _retval = _cs_calcium_read_float(_cs_calcium_component[comp_id],
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

  if (_cs_calcium_read_double != NULL)
    _retval = _cs_calcium_read_double(_cs_calcium_component[comp_id],
                                      _time_dep,
                                      min_time,
                                      max_time,
                                      iteration,
                                      _var_name,
                                      n_val_max,
                                      n_val_read,
                                      val);

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

  if (_cs_calcium_write_int != NULL)
    retval = _cs_calcium_write_int(_cs_calcium_component[comp_id],
                                   _time_dep,
                                   cur_time,
                                   iteration,
                                   _var_name,
                                   n_val,
                                   _val);

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

  if (_cs_calcium_write_float != NULL)
    retval = _cs_calcium_write_float(_cs_calcium_component[comp_id],
                                     _time_dep,
                                     cur_time,
                                     iteration,
                                     _var_name,
                                     n_val,
                                     _val);

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

  if (_cs_calcium_write_double != NULL)
    retval = _cs_calcium_write_double(_cs_calcium_component[comp_id],
                                      _time_dep,
                                      cur_time,
                                      iteration,
                                      _var_name,
                                      n_val,
                                      _val);

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

  _cs_calcium_component[comp_id] = comp;
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
  _cs_calcium_n_echo = n_echo;
}

/*----------------------------------------------------------------------------
 * Load YACS and corresponding Calcium functions.
 *
 * parameters:
 *   lib_path <-- path to shared library containing the yacsinit() function.
 *----------------------------------------------------------------------------*/

void
cs_calcium_load_yacs(const char *lib_path)
{
#if defined(HAVE_DLOPEN)

  /* Load symbols from shared library */

  _cs_calcium_yacslib = cs_base_dlopen(lib_path);

  /* Function pointers need to be double-casted so as to first convert
     a (void *) type to a memory address and then convert it back to the
     original type. Otherwise, the compiler may issue a warning.
     This is a valid ISO C construction. */

  _cs_calcium_yacsinit = (cs_calcium_yacsinit_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "yacsinit", true);

  _cs_calcium_read_int = (cs_calcium_read_int_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "cp_len", true);

  _cs_calcium_write_int = (cs_calcium_write_int_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "cp_een", true);

  _cs_calcium_read_float = (cs_calcium_read_float_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "cp_lre", true);

  _cs_calcium_write_float = (cs_calcium_write_float_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "cp_ere", true);

  _cs_calcium_read_double = (cs_calcium_read_double_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "cp_ldb", true);

  _cs_calcium_write_double = (cs_calcium_write_double_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_cs_calcium_yacslib, "cp_edb", true);

  if (   _cs_calcium_yacsinit == NULL
      || _cs_calcium_read_int == NULL
      || _cs_calcium_write_int == NULL
      || _cs_calcium_read_float == NULL
      || _cs_calcium_write_float == NULL
      || _cs_calcium_read_double== NULL
      || _cs_calcium_write_double == NULL) {
    cs_base_dlclose(lib_path, _cs_calcium_yacslib);
    _cs_calcium_yacslib = NULL;
  }

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to load: %s\n"), lib_path);

#endif
}

/*----------------------------------------------------------------------------
 * Unload YACS and corresponding Calcium functions
 *----------------------------------------------------------------------------*/

void
cs_calcium_unload_yacs(void)
{
#if defined(HAVE_DLOPEN)

  if (_cs_calcium_yacslib != NULL)
    cs_base_dlclose(NULL, _cs_calcium_yacslib);

  /* Reset function pointers to NULL */

  _cs_calcium_yacslib = NULL;

  _cs_calcium_yacsinit = NULL;

  _cs_calcium_read_int = NULL;
  _cs_calcium_write_int = NULL;
  _cs_calcium_read_float = NULL;
  _cs_calcium_write_float = NULL;
  _cs_calcium_read_double = NULL;
  _cs_calcium_write_double = NULL;

#endif
}

/*----------------------------------------------------------------------------
 * Initialize YACS component and enter event loop.
 *
 * This must be called after cs_calcium_load_yacs().
 *
 * Note that the YACS event loop does not return, sot the YACS component
 * description should ensure that the code's main run() method (or similar)
 * is called in the component body.
 *----------------------------------------------------------------------------*/

void
cs_calcium_start_yacs(void)
{
  if (_cs_calcium_yacslib != NULL)
    _cs_calcium_yacsinit();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
