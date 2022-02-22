/*============================================================================
 * Basic CALCIUM-mappable functions for code coupling
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

/* Maximum string Lengths */

#define CS_CALCIUM_VARIABLE_LEN 127

#if defined(HAVE_MPI)

static MPI_Comm _comm = MPI_COMM_WORLD;

/* MPI tag for file operations */

#define CS_CALCIUM_MPI_TAG  0

#endif

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Verbosity (none if -1, headers if 0,
   headers + n first and last elements if > 0 */

static int _cs_calcium_n_echo = 1;

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
_calcium_echo_body(cs_datatype_t           datatype,
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

    case CS_INT_TYPE:
      {
        const int *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : %12d\n", id + 1, *(_val + id));
      }
      break;

    case CS_DOUBLE:
      {
        const double *_val = val;
        for (id = start_id; id < end_id; id++)
          bft_printf("    %10d : %14.7e\n", id + 1, *(_val + id));
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
 *   rank_id    <-- communicating MPI rank id
 *   var_name   <-- variable name
 *   iteration  <-- iteration step
 *   datatype   <-- section data type
 *   n_max_vals <-- maximum number of values to read
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_pre_read(int                     rank_id,
                       const char             *var_name,
                       int                     iteration,
                       cs_datatype_t           datatype,
                       int                     n_max_vals)
{
  if (_cs_calcium_n_echo < 0)
    return;

  assert(var_name != NULL);

  bft_printf(_("\nRank %d, %s:\n"), rank_id, var_name);

  bft_printf(_("Reading up to %d values of type %s (iteration %d) ..."),
             n_max_vals, cs_datatype_name[datatype],
             iteration);
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print message indicating that we are finished reading data, and optionnaly
 * print a part of the corresponding data.
 *
 * parameters:
 *   iteration <-- iteration step
 *   datatype  <-- section data type
 *   n_val     <-- number of values in array
 *   val       <-- array values
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_post_read(int                    iteration,
                        cs_datatype_t          datatype,
                        int                    n_val,
                        const void            *val)
{
  if (_cs_calcium_n_echo < 0)
    return;

  bft_printf(_("[ok]\n"
               "Read          %d values (iteration %d).\n"),
             n_val, iteration);

  _calcium_echo_body(datatype, _cs_calcium_n_echo, n_val, val);
}

/*----------------------------------------------------------------------------
 * Print message indicating that we are ready to write data
 *
 * parameters:
 *   rank_id   <-- communicating MPI rank id
 *   var_name  <-- variable name
 *   iteration <-- iteration step
 *   datatype  <-- section data type
 *   n_vals    <-- number of values to read
 *----------------------------------------------------------------------------*/

static void
_calcium_echo_pre_write(int                     rank_id,
                        const char             *var_name,
                        int                     iteration,
                        cs_datatype_t           datatype,
                        int                     n_vals)
{
  if (_cs_calcium_n_echo < 0)
    return;

  assert(var_name != NULL);

  bft_printf(_("\nRank %d, %s:\n"), rank_id, var_name);

  bft_printf(_("Writing %d values of type %s (iteration %d) ..."),
             n_vals, cs_datatype_name[datatype],
             iteration);
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
_calcium_echo_post_write(cs_datatype_t          datatype,
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
 * Read integer values, blocking until they are available.
 *
 * parameters:
 *   rank_id    <-- communicating MPI rank id
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
cs_calcium_read_int(int                    rank_id,
                    int                   *iteration,
                    const char            *var_name,
                    int                    n_val_max,
                    int                   *n_val_read,
                    int                    val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  _retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_read(rank_id,
                         _var_name, *iteration,
                         CS_INT_TYPE, n_val_max);

#if defined(HAVE_MPI)

  char var_cmp[CS_CALCIUM_VARIABLE_LEN + 1];
  int meta[3] = {0, 0, 0};

  MPI_Status status;
  MPI_Recv(var_cmp, CS_CALCIUM_VARIABLE_LEN + 1, MPI_CHAR, rank_id,
           CS_CALCIUM_MPI_TAG, _comm, &status);

  if (strncmp(var_cmp, _var_name, CS_CALCIUM_VARIABLE_LEN + 1) != 0) {
    bft_printf("\n"
               "Warning: received %s\n"
               "         expected %s\n",
               _var_name, var_cmp);
    bft_printf_flush();
  }

  MPI_Recv(meta, 3, MPI_INT, rank_id,
           CS_CALCIUM_MPI_TAG, _comm, &status);

  if (meta[0] != *iteration || meta[1] != n_val_max || meta[2] != 4) {
    bft_printf("\n"
               "Warning: received [%d, %d, %d] for %s\n"
               "         expected [%d, %d, %d]\n",
               meta[0], meta[1], meta[2], _var_name,
               *iteration, n_val_max, 4);
    bft_printf_flush();
  }

  MPI_Recv(meta, n_val_max, MPI_INT, rank_id,
           CS_CALCIUM_MPI_TAG, _comm, &status);

  MPI_Get_count(&status, MPI_INT, n_val_read);

#endif

  _calcium_echo_post_read(*iteration, CS_INT_TYPE, *n_val_read, val);

  return _retval;
}

/*----------------------------------------------------------------------------
 * Read double-precision floating-point values,
 * blocking until they are available.
 *
 * parameters:
 *   rank_id    <-- communicating MPI rank id
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
cs_calcium_read_double(int                    rank_id,
                       int                   *iteration,
                       const char            *var_name,
                       int                    n_val_max,
                       int                   *n_val_read,
                       double                 val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  _retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_read(rank_id, _var_name, *iteration,
                         CS_DOUBLE, n_val_max);

#if defined(HAVE_MPI)

  char var_cmp[CS_CALCIUM_VARIABLE_LEN + 1];
  int meta[3] = {0, 0, 0};

  MPI_Status status;
  MPI_Recv(var_cmp, CS_CALCIUM_VARIABLE_LEN + 1, MPI_CHAR, rank_id,
           CS_CALCIUM_MPI_TAG, _comm, &status);

  if (strncmp(var_cmp, _var_name, CS_CALCIUM_VARIABLE_LEN + 1) != 0) {
    bft_printf("\n"
               "Warning: received %s\n"
               "         expected %s\n",
               _var_name, var_cmp);
    bft_printf_flush();
  }

  MPI_Recv(meta, 3, MPI_INT, rank_id,
           CS_CALCIUM_MPI_TAG, _comm, &status);

  if (meta[0] != *iteration || meta[1] != n_val_max || meta[2] != 8) {
    bft_printf("\n"
               "Warning: received [%d, %d, %d] for %s\n"
               "         expected [%d, %d, %d]\n",
               meta[0], meta[1], meta[2], _var_name,
               *iteration, n_val_max, 8);
    bft_printf_flush();
  }

  MPI_Recv(meta, n_val_max, MPI_DOUBLE, rank_id,
           CS_CALCIUM_MPI_TAG, _comm, &status);

  MPI_Get_count(&status, MPI_DOUBLE, n_val_read);

#endif

  _calcium_echo_post_read(*iteration, CS_DOUBLE, *n_val_read, val);

  return _retval;
}

/*----------------------------------------------------------------------------
 * Write integer values.
 *
 * parameters:
 *   rank_id    <-- communicating MPI rank id
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_int(int                    rank_id,
                     int                    iteration,
                     const char            *var_name,
                     int                    n_val,
                     const int              val[])
{
  char _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  int  *_val = NULL;

  int retval = 0;

  memset(_var_name, 0, CS_CALCIUM_VARIABLE_LEN + 1);
  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_write(rank_id,
                          _var_name, iteration, CS_INT_TYPE, n_val);

  BFT_MALLOC(_val, n_val, int);
  memcpy(_val, val, n_val * sizeof(int));

#if defined(HAVE_MPI)

  int meta[3] = {iteration, n_val, 4};

  MPI_Send(_var_name, CS_CALCIUM_VARIABLE_LEN + 1, MPI_CHAR, rank_id,
           CS_CALCIUM_MPI_TAG, _comm);
  MPI_Send(meta, 3, MPI_INT, rank_id, CS_CALCIUM_MPI_TAG, _comm);

  MPI_Send(_val, n_val, MPI_INT, rank_id, CS_CALCIUM_MPI_TAG, _comm);

#endif

  BFT_FREE(_val);

  _calcium_echo_post_write(CS_INT_TYPE, n_val, val);

  return retval;
}

/*----------------------------------------------------------------------------
 * Write double-precision float values.
 *
 * parameters:
 *   rank_id    <-- communicating MPI rank id
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_double(int                    rank_id,
                        int                    iteration,
                        const char            *var_name,
                        int                    n_val,
                        const double           val[])
{
  char   _var_name[CS_CALCIUM_VARIABLE_LEN + 1];
  double *_val = NULL;

  int retval = 0;

  strncpy(_var_name, var_name, CS_CALCIUM_VARIABLE_LEN);

  _calcium_echo_pre_write(rank_id,
                          _var_name, iteration, CS_DOUBLE, n_val);

  BFT_MALLOC(_val, n_val, double);
  memcpy(_val, val, n_val * sizeof(double));

#if defined(HAVE_MPI)

  int meta[3] = {iteration, n_val, 8};

  MPI_Send(_var_name, CS_CALCIUM_VARIABLE_LEN + 1, MPI_CHAR, rank_id,
           CS_CALCIUM_MPI_TAG, _comm);
  MPI_Send(meta, 3, MPI_INT, rank_id, CS_CALCIUM_MPI_TAG, _comm);

  MPI_Send(_val, n_val, MPI_DOUBLE, rank_id, CS_CALCIUM_MPI_TAG, _comm);

#endif

  BFT_FREE(_val);

  _calcium_echo_post_write(CS_DOUBLE, n_val, val);

  return retval;
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
