/*============================================================================
 * Functions dealing with parallelism
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_parall.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

#define CS_PARALL_ARRAY_SIZE  500

typedef struct
{
  double  val;
  int     rank;
} _mpi_double_int_t;

/*============================================================================
 * Global static variables
 *============================================================================*/

#if defined(HAVE_MPI)

/* Minimum recommended scatter/gather buffer size */

static size_t _cs_parall_min_coll_buf_size = 1024*1024*8;

#endif

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Call MPI_Allreduce for a given Code_Saturne datatype and MPI
 * operation on all default communicator processes.
 *
 * parameters:
 *   n         <-- number of values
 *   datatype  <-- matching Code_Saturne datatype
 *   operation <-- MPI operation
 *   val       <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

static void
_cs_parall_allreduce(int             n,
                     cs_datatype_t   datatype,
                     MPI_Op          operation,
                     void           *val)
{
  int             data_size = n*cs_datatype_size[datatype];
  unsigned char  *locval;
  unsigned char  _locval[256];

  if (data_size > 256)
    BFT_MALLOC(locval, data_size, unsigned char);
  else
    locval = _locval;

  memcpy(locval, val, data_size);

  MPI_Allreduce(locval, val, n, cs_datatype_to_mpi[datatype], operation,
                cs_glob_mpi_comm);

  if (locval != _locval)
    BFT_FREE(locval);
}

#endif

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the maximum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * subroutine parcmx (counter)
 * *****************
 *
 * integer          counter       <-> : input = local counter
 *                                      output = global max counter
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcmx, PARCMX)(cs_int_t  *counter)
{
#if defined(HAVE_MPI)

  cs_int_t  global_max;

  assert(sizeof(int) == sizeof(cs_int_t));

  MPI_Allreduce(counter, &global_max, 1, CS_MPI_INT, MPI_MAX,
                cs_glob_mpi_comm);

  *counter = global_max;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the minimum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * subroutine parcmn (counter)
 * *****************
 *
 * integer          counter       <-> : input = local counter
 *                                      output = global min counter
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcmn, PARCMN)(cs_int_t  *counter)
{
#if defined(HAVE_MPI)

  cs_int_t  global_max;

  assert(sizeof(int) == sizeof(cs_int_t));

  MPI_Allreduce(counter, &global_max, 1, CS_MPI_INT, MPI_MIN,
                cs_glob_mpi_comm);

  *counter = global_max;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global sum of a counter (int) for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parcpt (counter)
 * *****************
 *
 * integer          counter     : <-> : input = counter to sum
 *                                      output = global sum
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcpt, PARCPT)(cs_int_t  *counter)
{
#if defined(HAVE_MPI)

  cs_int_t  global_sum;

  assert(sizeof(int) == sizeof(cs_int_t));

  MPI_Allreduce(counter, &global_sum, 1, CS_MPI_INT, MPI_SUM,
                cs_glob_mpi_comm);

  *counter = global_sum;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global sum of a real for the entire domain in case of parellism
 *
 * Fortran Interface :
 *
 * subroutine parsom (var)
 * *****************
 *
 * double precision var         : <-> : input = value to sum
 *                                      output = global sum
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parsom, PARSOM)(cs_real_t  *var)
{
#if defined(HAVE_MPI)

  cs_real_t  global_sum;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (var, &global_sum, 1, CS_MPI_REAL, MPI_SUM,
                 cs_glob_mpi_comm);

  *var = global_sum;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the maximum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parmax (var)
 * *****************
 *
 * double precision var         : <-> : input = local maximum value
 *                                      output = global maximum value
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmax, PARMAX)(cs_real_t  *var)
{
#if defined(HAVE_MPI)

  cs_real_t global_max;

  assert(sizeof(double) == sizeof(cs_real_t));

  MPI_Allreduce(var, &global_max, 1, CS_MPI_REAL, MPI_MAX,
                cs_glob_mpi_comm);

  *var = global_max;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the minimum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parmin (var)
 * *****************
 *
 * double precision var         : <-> : input = local minimum value
 *                                      output = global minimum value
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmin, PARMIN)(cs_real_t  *var)
{
#if defined(HAVE_MPI)

  cs_real_t global_min;

  assert(sizeof(double) == sizeof(cs_real_t));

  MPI_Allreduce(var, &global_min, 1, CS_MPI_REAL, MPI_MIN,
                cs_glob_mpi_comm);

  *var = global_min;

#endif
}

/*----------------------------------------------------------------------------
 * Maximum value of a real and the value of the related variable for the entire
 * domain in case of parallelism.
 *
 * Fortran Interface
 *
 * subroutine parmxl (nbr, var, xyzvar)
 * *****************
 *
 * integer          nbr         : <-- : size of the related variable
 * double precision var         : <-> : input = local max. value
 *                                      output = global max. value
 * double precision xyzvar(nbr) : <-> : input = value related to local max.
 *                                      output = value related to global max.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmxl, PARMXL)(cs_int_t   *nbr,
                          cs_real_t  *var,
                          cs_real_t   xyzvar[])
{
#if defined(HAVE_MPI)

  _mpi_double_int_t  val_in, val_max;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *var;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                cs_glob_mpi_comm);

  *var = val_max.val;

  MPI_Bcast(xyzvar, *nbr, CS_MPI_REAL, val_max.rank, cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Minimum value of a real and the value of the related variable for the entire
 * domain in case of parallelism.
 *
 * Fortran Interface
 *
 * subroutine parmnl (nbr, var, xyzvar)
 * *****************
 *
 * integer          nbr         : <-- : size of the related variable
 * double precision var         : <-> : input = local max. value
 *                                      output = global max. value
 * double precision xyzvar(nbr) : <-> : input = value related to local max.
 *                                      output = value related to global max.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmnl, PARMNL)(cs_int_t   *nbr,
                          cs_real_t  *var,
                          cs_real_t   xyzvar[])
{
#if defined(HAVE_MPI)

  _mpi_double_int_t  val_in, val_min;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *var;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                cs_glob_mpi_comm);

  *var = val_min.val;

  MPI_Bcast(xyzvar, *nbr, CS_MPI_REAL, val_min.rank, cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global sum for each element of an array of int in case of
 * parallelism.
 *
 * Fortran Interface
 *
 * subroutine parism (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * integer          array(*)     : <-> : input = local array
 *                                       output = array of global sum values.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parism, PARISM)(cs_int_t  *n_elts,
                          cs_int_t   array[])
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_int_t  set_sum_array[CS_PARALL_ARRAY_SIZE];

  cs_int_t  *sum_array = NULL;

  assert(sizeof(int) == sizeof(cs_int_t));

  if (CS_PARALL_ARRAY_SIZE < *n_elts) {

    BFT_MALLOC(sum_array, *n_elts, cs_int_t);

    MPI_Allreduce(array, sum_array, *n_elts, CS_MPI_INT, MPI_SUM,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = sum_array[i];

    BFT_FREE(sum_array);

  }
  else {

    MPI_Allreduce(array, set_sum_array, *n_elts, CS_MPI_INT, MPI_SUM,
                  cs_glob_mpi_comm);

    for (i = 0 ; i < *n_elts ; i++)
      array[i] = set_sum_array[i];

  }

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global maximum value for each element of an array of int in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parimx (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * integer          array(*)     : <-> : input = local array
 *                                       output = array of global max. values.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parimx, PARIMX)(cs_int_t  *n_elts,
                          cs_int_t   array[])
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_int_t  set_globmax_array[CS_PARALL_ARRAY_SIZE];

  cs_int_t  *globmax_array = NULL;

  assert(sizeof(int) == sizeof(cs_int_t));

  if (CS_PARALL_ARRAY_SIZE < *n_elts) {

    BFT_MALLOC(globmax_array, *n_elts, cs_int_t);

    MPI_Allreduce(array, globmax_array, *n_elts, CS_MPI_INT, MPI_MAX,
                  cs_glob_mpi_comm);

    for (i = 0 ; i < *n_elts ; i++)
        array[i] = globmax_array[i];

    BFT_FREE(globmax_array);

  }
  else {

    MPI_Allreduce(array, set_globmax_array, *n_elts, CS_MPI_INT, MPI_MAX,
                  cs_glob_mpi_comm);

    for (i = 0 ; i < *n_elts ; i++)
      array[i] = set_globmax_array[i];

  }

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global minimum value for each element of an array of int in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parimn (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * integer          array(*)     : <-> : input = local array
 *                                       output = array of global min. values.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parimn, PARIMN)(cs_int_t  *n_elts,
                          cs_int_t   array[])
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_int_t  set_globmin_array[CS_PARALL_ARRAY_SIZE];

  cs_int_t  *globmin_array = NULL;

  assert(sizeof(int) == sizeof(cs_int_t));

  if (CS_PARALL_ARRAY_SIZE < *n_elts) {

    BFT_MALLOC(globmin_array, *n_elts, cs_int_t);

    MPI_Allreduce (array, globmin_array, *n_elts, CS_MPI_INT, MPI_MIN,
                   cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
        array[i] = globmin_array[i];

    BFT_FREE(globmin_array);

  }
  else {

    MPI_Allreduce(array, set_globmin_array, *n_elts, CS_MPI_INT, MPI_MIN,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = set_globmin_array[i];

  }

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global sum for each element of an array of real in case of
 * parallelism.
 *
 * Fortran Interface
 *
 * subroutine parrsm (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * double precision array(*)     : <-> : input = local array
 *                                       output = array of global sum values.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parrsm, PARRSM)(cs_int_t   *n_elts,
                          cs_real_t   array[])
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_real_t  set_sum_array[CS_PARALL_ARRAY_SIZE];

  cs_real_t  *sum_array = NULL;

  assert(sizeof(double) == sizeof(cs_real_t));

  if (CS_PARALL_ARRAY_SIZE < *n_elts) {

    BFT_MALLOC(sum_array, *n_elts, cs_real_t);

    MPI_Allreduce(array, sum_array, *n_elts, CS_MPI_REAL, MPI_SUM,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
        array[i] = sum_array[i];

    BFT_FREE(sum_array);

  }
  else {

    MPI_Allreduce(array, set_sum_array, *n_elts, CS_MPI_REAL, MPI_SUM,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = set_sum_array[i];

  }

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global maximum value for each element of an array of real in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parrmx (n_elts, array)
 * *****************
 *
 * integer          n_elts        : <-- : size of the array
 * double precision array(*)      : <-> : input = local array
 *                                        output = array of global max. values.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parrmx, PARRMX)(cs_int_t   *n_elts,
                          cs_real_t   array[])
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_real_t  set_globmax_array[CS_PARALL_ARRAY_SIZE];

  cs_real_t  *globmax_array = NULL;

  assert(sizeof(double) == sizeof(cs_real_t));

  if (CS_PARALL_ARRAY_SIZE < *n_elts) {

    BFT_MALLOC(globmax_array, *n_elts, cs_real_t);

    MPI_Allreduce(array, globmax_array, *n_elts, CS_MPI_REAL, MPI_MAX,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = globmax_array[i];

    BFT_FREE(globmax_array);

  }
  else {

    MPI_Allreduce(array, set_globmax_array, *n_elts, CS_MPI_REAL, MPI_MAX,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = set_globmax_array[i];

  }

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global minimum value for each element of an array of real in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parrmn (n_elts, array)
 * *****************
 *
 * integer          n_elts        : <-- : size of the array
 * double precision array(*)      : <-> : input = local array
 *                                        output = array of global min. values.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parrmn, PARRMN)(cs_int_t   *n_elts,
                          cs_real_t   array[])
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_real_t  set_globmin_array[CS_PARALL_ARRAY_SIZE];

  cs_real_t  *globmin_array = NULL;

  assert(sizeof(double) == sizeof(cs_real_t));

  if (CS_PARALL_ARRAY_SIZE < *n_elts) {

    BFT_MALLOC(globmin_array, *n_elts, cs_real_t);

    MPI_Allreduce(array, globmin_array, *n_elts, CS_MPI_REAL, MPI_MIN,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = globmin_array[i];

    BFT_FREE(globmin_array);

  }
  else {

    MPI_Allreduce(array, set_globmin_array, *n_elts, CS_MPI_REAL, MPI_MIN,
                  cs_glob_mpi_comm);

    for (i = 0; i < *n_elts; i++)
      array[i] = set_globmin_array[i];

  }

#endif
}

/*----------------------------------------------------------------------------
 * Broadcast to all the ranks the value of each element of an array of int.
 * (encapsulation of MPI_Bcast())
 *
 * Fortran Interface :
 *
 * subroutine parbci (irank, n_elts, array)
 * *****************
 *
 * integer          irank       : <-- : rank related to the sending process
 * integer          n_elts      : <-- : size of the array
 * integer          array(*)    : <-> : array t broadcast
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbci, PARBCI)(cs_int_t   *irank,
                          cs_int_t   *n_elts,
                          cs_int_t    array[])
{
#if defined(HAVE_MPI)

  MPI_Bcast(array, *n_elts, CS_MPI_INT, *irank, cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Broadcast to all the ranks the value of each element of an array of real.
 * (encapsulation of MPI_Bcast())
 *
 * Fortran Interface :
 *
 * subroutine parbcr (irank, n_elts, array)
 * *****************
 *
 * integer            irank     : <-- : rank related to the sending process
 * integer            n_elts    : <-- : size of the array
 * double precision   array(*)  : <-> : array to broadcast
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbcr, PARBCR)(cs_int_t   *irank,
                          cs_int_t   *n_elts,
                          cs_real_t   array[])
{
#if defined(HAVE_MPI)

  MPI_Bcast(array, *n_elts, CS_MPI_REAL, *irank, cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Build a global array from each local array in each domain. The size of
 * each local array can be different.
 *
 * Fortran Interface :
 *
 * subroutine paragv (nvar, nvargb, var, vargb)
 * *****************
 *
 * integer           n_elts      : <-- : size of the local array
 * integer           n_g_elts    : <-- : size of the global array
 * double precision  array(*)    : <-- : local array
 * double precision  g_array(*)  : --> : global array
 *----------------------------------------------------------------------------*/

void
CS_PROCF (paragv, PARAGV)(cs_int_t   *n_elts,
                          cs_int_t   *n_g_elts,
                          cs_real_t   array[],
                          cs_real_t  *g_array)
{
#if defined(HAVE_MPI)

  int  i;
  int  *count = NULL;
  int  *shift = NULL;

  const int  n_domains = cs_glob_n_ranks;

  assert(sizeof(double) == sizeof(cs_real_t));

  BFT_MALLOC(count, n_domains, cs_int_t);
  BFT_MALLOC(shift, n_domains, cs_int_t);

  MPI_Allgather(n_elts, 1, CS_MPI_INT, count, 1, CS_MPI_INT,
                cs_glob_mpi_comm);

  shift[0] = 0;
  for (i = 1; i < n_domains; i++)
    shift[i] = shift[i-1] + count[i-1];

  assert(*n_g_elts == (shift[n_domains - 1] + count[n_domains - 1]));

  MPI_Allgatherv(array, *n_elts, CS_MPI_REAL,
                 g_array, count, shift, CS_MPI_REAL, cs_glob_mpi_comm);

  BFT_FREE(count);
  BFT_FREE(shift);

#endif
}

/*----------------------------------------------------------------------------
 * Find a node which minimizes a given distance and its related rank.
 * May be used to locate a node among several domains.
 *
 * Fortran Interface :
 *
 * subroutine parfpt (node, ndrang, dis2mn)
 * *****************
 *
 * integer          node        : <-> : local number of the closest node
 * integer          ndrang      : --> : rank id for which the distance is the
 *                                      smallest
 * double precision dis2mn      : <-- : square distance between the closest node
 *                                      and the wanted node.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfpt, PARFPT)(cs_int_t   *node,
                          cs_int_t   *ndrang,
                          cs_real_t  *dis2mn)
{
#if defined(HAVE_MPI)

  cs_int_t buf[2];

  _mpi_double_int_t  val_in, val_min;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *dis2mn;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                cs_glob_mpi_comm);

  *ndrang = cs_glob_rank_id;

  buf[0] = *node;
  buf[1] = *ndrang;

  MPI_Bcast(buf, 2, CS_MPI_INT, val_min.rank, cs_glob_mpi_comm);

  *node = buf[0];
  *ndrang = buf[1];

#endif
}

/*----------------------------------------------------------------------------
 * Return the value associated to a probe.
 *
 * Fortran Interface :
 *
 * subroutine parhis (node, ndrang, var, varcap)
 * *****************
 *
 * integer          node        : <-- : local number of the element related to
 *                                      a measure node
 * integer          ndrang      : <-- : rank of the process owning the closest
 *                                      node from the measure node
 * double precision var(*)      : <-- : values of the variable on local elements
 * double precision varcap      : --> : value of the variable for the element
 *                                      related to the measure node
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parhis, PARHIS)(cs_int_t   *node,
                          cs_int_t   *ndrang,
                          cs_real_t   var[],
                          cs_real_t  *varcap)
{
#if defined(HAVE_MPI)

  assert(sizeof(double) == sizeof(cs_real_t));

  if (*ndrang == cs_glob_rank_id)
    *varcap = var[*node - 1];
  else
    *varcap = 0.0;

  MPI_Bcast(varcap, 1, CS_MPI_REAL, *ndrang, cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Call a barrier in case of parallelism
 *
 * This function should not be necessary in production code,
 * but it may be useful for debugging purposes.
 *
 * Fortran interface :
 *
 * subroutine parbar
 * *****************
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbar, PARBAR)(void)
{
#if defined(HAVE_MPI)

  MPI_Barrier(cs_glob_mpi_comm);

#endif
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_counter(cs_gnum_t   cpt[],
                  const int   n)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, CS_GNUM_TYPE, MPI_SUM, cpt);
}
#endif

/*----------------------------------------------------------------------------
 * Maximum values of a counter on all default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_counter_max(cs_lnum_t   cpt[],
                      const int   n)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, CS_LNUM_TYPE, MPI_MAX, cpt);
}

#endif

/*----------------------------------------------------------------------------
 * Sum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching Code_Saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_sum(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, datatype, MPI_SUM, val);
}

#endif

/*----------------------------------------------------------------------------
 * Maximum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching Code_Saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_max(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, datatype, MPI_MAX, val);
}

#endif

/*----------------------------------------------------------------------------
 * Minimum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching Code_Saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_min(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, datatype, MPI_MIN, val);
}

#endif

/*----------------------------------------------------------------------------
 * Return minimum recommended scatter or gather buffer size.
 *
 * This is used by some internal part to block or scatter/gather algorithms,
 * so as to allow I/O buffer size tuning.
 *
 * returns:
 *   minimum recommended part to block or gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

size_t
cs_parall_get_min_coll_buf_size(void)
{
#if defined(HAVE_MPI)
  return _cs_parall_min_coll_buf_size;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------
 * Define minimum recommended gather buffer size.
 *
 * This is used by some internal part to block or scatter/gather algorithms,
 * so as to allow I/O buffer size tuning.
 *
 * parameters:
 *   minimum recommended part to block or gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

void
cs_parall_set_min_coll_buf_size(size_t buffer_size)
{
#if defined(HAVE_MPI)
  _cs_parall_min_coll_buf_size = buffer_size;
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
