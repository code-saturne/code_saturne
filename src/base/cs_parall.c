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
 * Functions dealing with parallelism
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_perio.h" /* Only needed for PARCOM */

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

#define CS_PARALL_DEBUG_COUNT  0
#define CS_PARALL_ARRAY_SIZE  500


/*============================================================================
 *  Global static variables
 *============================================================================*/

/* Counter of the number of call for each function */

#if CS_PARALL_DEBUG_COUNT
static  cs_int_t n_total_par_calls = 0;
static  cs_int_t n_pargve_calls = 0;
static  cs_int_t n_parcom_calls = 0;
static  cs_int_t n_parcve_calls = 0;
static  cs_int_t n_parcmx_calls = 0;
static  cs_int_t n_parcmn_calls = 0;
static  cs_int_t n_parcpt_calls = 0;
static  cs_int_t n_parsom_calls = 0;
static  cs_int_t n_parmax_calls = 0;
static  cs_int_t n_parmin_calls = 0;
static  cs_int_t n_parmxl_calls = 0;
static  cs_int_t n_parmnl_calls = 0;
static  cs_int_t n_parism_calls = 0;
static  cs_int_t n_parimx_calls = 0;
static  cs_int_t n_parimn_calls = 0;
static  cs_int_t n_parrsm_calls = 0;
static  cs_int_t n_parrmx_calls = 0;
static  cs_int_t n_parrmn_calls = 0;
static  cs_int_t n_parbci_calls = 0;
static  cs_int_t n_parbcr_calls = 0;
static  cs_int_t n_paragv_calls = 0;
static  cs_int_t n_parfpt_calls = 0;
static  cs_int_t n_parhis_calls = 0;
static  cs_int_t n_parcel_calls = 0;
static  cs_int_t n_interface_sr_calls = 0;
#endif

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Public functions definition for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update a buffer on cells in case of parallelism
 *
 * This function copies values of the cells in the standard send_halo
 * (local cells) to ghost cells on distant ranks.
 *
 * Fortran interface :
 *
 * SUBROUTINE PARCOM (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR(NCELET) : <-> : variable on cells, output is an update
 *                                      of VAR(NCEL+1..NCELET)
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcom, PARCOM)(cs_real_t  var[])
{

  if (cs_glob_mesh->have_rotation_perio != 0)
    cs_perio_save_rotation_halo(cs_glob_mesh->halo, CS_HALO_STANDARD, var);

  cs_halo_sync_var(cs_glob_mesh->halo, CS_HALO_STANDARD, var);

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parcom\n",
         cs_glob_rank_id, n_parcom_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Update a buffer on cells in case of parallelism
 *
 * This function copies values of the cells in the entire (i.e. std + ext)
 * send_halo (local cells) to ghost cells on distant ranks.
 *
 * PVAR has to be well allocated => n_cells + n_cells_with_ghosts where
 * n_cells_with_ghosts = n_std_ghost_cells + n_ext_ghost_cells.
 *
 * Fortran interface :
 *
 * SUBROUTINE PARCVE
 * *****************
 *
 * DOUBLE PRECISION  PVAR      : <-> : variable buffer to sync
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcve, PARCVE)(cs_real_t  pvar[])
{

  cs_halo_sync_var(cs_glob_mesh->halo, CS_HALO_EXTENDED, pvar);

#if CS_PARALL_DEBUG_COUNT
  printf ("irang = %d, iappel = %d, tot = %d, parcve\n",
          cs_glob_rank_id, n_parcve_calls++, n_total_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the maximum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARCMX (IND)
 * *****************
 *
 * INTEGER          COUNTER       <-> : input = local counter
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parcmx\n",
         cs_glob_rank_id, n_parcmx_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the minimum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARCMN (IND)
 * *****************
 *
 * INTEGER          COUNTER       <-> : input = local counter
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parcmn\n",
         cs_glob_rank_id, n_parcmn_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global sum of a counter (int) for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARCPT (CPT)
 * *****************
 *
 * INTEGER          COUNTER     : <-> : input = counter to sum
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parcpt\n",
         cs_glob_rank_id, n_parcpt_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global sum of a real for the entire domain in case of parellism
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARSOM (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR         : <-> : input = value to sum
 *                                      output = global sum
 *
 * Returns:
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parsom\n",
         cs_glob_rank_id, n_parsom_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the maximum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARMAX (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR         : <-> : input = local maximum value
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parmax\n",
         cs_glob_rank_id, n_parmax_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the minimum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARMIN (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR         : <-> : input = local minimum value
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parmin\n",
         cs_glob_rank_id, n_parmin_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Maximum value of a real and the value of the related variable for the entire
 * domain in case of parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARMXL (NBR, VAR, XYZVAR)
 * *****************
 *
 * INTEGER          NBR         : --> : size of the related variable
 * DOUBLE PRECISION VAR         : <-> : input = local max. value
 *                                      output = global max. value
 * DOUBLE PRECISION XYZVAR(NBR) : <-> : input = value related to local max.
 *                                      output = value related to global max.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmxl, PARMXL)(cs_int_t   *nbr,
                          cs_real_t  *var,
                          cs_real_t   xyzvar[])
{
#if defined(HAVE_MPI)

  cs_mpi_real_int_t  val_in, val_max;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *var;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_max, 1, CS_MPI_REAL_INT, MPI_MAXLOC,
                cs_glob_mpi_comm);

  *var = val_max.val;

  MPI_Bcast(xyzvar, *nbr, CS_MPI_REAL, val_max.rank, cs_glob_mpi_comm);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parmxl\n",
         cs_glob_rank_id, n_parmxl_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Minimum value of a real and the value of the related variable for the entire
 * domain in case of parallelism.
 *
 * Fortran Interface
 *
 * Interface Fortran :
 *
 * SUBROUTINE PARMNL (NBR, VAR, XYZVAR)
 * *****************
 *
 * INTEGER          NBR         : --> : size of the related variable
 * DOUBLE PRECISION VAR         : <-> : input = local max. value
 *                                      output = global max. value
 * DOUBLE PRECISION XYZVAR(NBR) : <-> : input = value related to local max.
 *                                      output = value related to global max.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmnl, PARMNL)(cs_int_t   *nbr,
                          cs_real_t  *var,
                          cs_real_t   xyzvar[])
{
#if defined(HAVE_MPI)

  cs_mpi_real_int_t  val_in, val_min;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *var;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_min, 1, CS_MPI_REAL_INT, MPI_MINLOC,
                cs_glob_mpi_comm);

  *var = val_min.val;

  MPI_Bcast(xyzvar, *nbr, CS_MPI_REAL, val_min.rank, cs_glob_mpi_comm);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parmnl\n",
         cs_glob_rank_id, n_parmnl_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global sum for each element of an array of int in case of
 * parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARISM (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * INTEGER          ARRAY(*)     : <-> : input = local array
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parism\n",
         cs_glob_rank_id, n_parism_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global maximum value for each element of an array of int in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARIMX (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * INTEGER          ARRAY(*)     : <-> : input = local array
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

#if CS_PARALL_DEBUG_COUNT
  printf ("irang = %d, iappel = %d, tot = %d, parimx\n",
          cs_glob_rank_id, n_parimx_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global minimum value for each element of an array of int in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARIMN (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * INTEGER          ARRAY(*)     : <-> : input = local array
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parimn\n",
         cs_glob_rank_id, n_parimn_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global sum for each element of an array of real in case of
 * parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARRSM (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * DOUBLE PRECISION ARRAY(*)     : <-> : input = local array
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parrsm\n",
         cs_glob_rank_id, n_parrsm_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global maximum value for each element of an array of real in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARRMX (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS        : --> : size of the array
 * DOUBLE PRECISION ARRAY(*)      : <-> : input = local array
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parrmx\n",
         cs_glob_rank_id, n_parrmx_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Compute the global minimum value for each element of an array of real in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARRMN (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS        : --> : size of the array
 * DOUBLE PRECISION ARRAY(*)      : <-> : input = local array
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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parrmn\n",
         cs_glob_rank_id, n_parrmn_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Broadcast to all the ranks the value of each element of an array of int.
 * (encapsulation of MPI_Bcast())
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARBCI (IRANK, N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          IRANK       : --> : rank related to the sending process
 * INTEGER          N_ELTS      : --> : size of the array
 * INTEGER          ARRAY(*)    : <-> : array of int
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbci, PARBCI)(cs_int_t   *irank,
                          cs_int_t   *n_elts,
                          cs_int_t    array[])
{
#if defined(HAVE_MPI)

  MPI_Bcast(array, *n_elts, CS_MPI_INT, *irank, cs_glob_mpi_comm);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parbci\n",
         cs_glob_rank_id, n_parbci_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Broadcast to all the ranks the value of each element of an array of real.
 * (encapsulation of MPI_Bcast())
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARBCR (IRANK, N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER            IRANK     : --> : rank related to the sending process
 * INTEGER            N_ELTS    : --> : size of the array
 * DOUBLE PRECISION   ARRAY(*)  : <-> : array of real
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbcr, PARBCR)(cs_int_t   *irank,
                          cs_int_t   *n_elts,
                          cs_real_t   array[])
{
#if defined(HAVE_MPI)

  MPI_Bcast(array, *n_elts, CS_MPI_REAL, *irank, cs_glob_mpi_comm);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parbcr\n",
         cs_glob_rank_id, n_parbcr_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Build a global array from each local array in each domain. The size of
 * each local array can be different.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARAGV (NVAR, NVARGB, VAR, VARGB)
 * *****************
 *
 * INTEGER           N_ELTS      : --> : size of the local array
 * INTEGER           N_G_ELTS    : --> : size of the global array
 * DOUBLE PRECISION  ARRAY(*)    : --> : local array
 * DOUBLE PRECISION  G_ARRAY(*)  : <-- : global array
 *----------------------------------------------------------------------------*/

void
CS_PROCF (paragv, PARAGV)(cs_int_t   *n_elts,
                          cs_int_t   *n_g_elts,
                          cs_real_t   array[],
                          cs_real_t  *g_array)
{
#if defined(HAVE_MPI)

  cs_int_t  i;
  cs_int_t  *count = NULL;
  cs_int_t  *shift = NULL;

  const cs_int_t  n_domains = cs_glob_mesh->n_domains;

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

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, paragv\n",
         cs_glob_rank_id, n_paragv_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Find a node which minimizes a given distance and its related rank.
 * May be used to locate a node among several domains.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARFPT (NODE, NDRANG, DIS2MN)
 * *****************
 *
 * INTEGER          NODE        : <-> : local number of the closest node
 * INTEGER          NDRANG      : <-- : rank id for which the distance is the
 *                                      smallest
 * DOUBLE PRECISION DIS2MN      : --> : square distance between the closest node
 *                                      and the wanted node.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfpt, PARFPT)(cs_int_t   *node,
                          cs_int_t   *ndrang,
                          cs_real_t  *dis2mn)
{
#if defined(HAVE_MPI)

  cs_mpi_real_int_t  val_in, val_min;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *dis2mn;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_min, 1, CS_MPI_REAL_INT, MPI_MINLOC,
                cs_glob_mpi_comm);

  *ndrang = cs_glob_rank_id;

  MPI_Bcast(node,   1, CS_MPI_INT, val_min.rank, cs_glob_mpi_comm);
  MPI_Bcast(ndrang, 1, CS_MPI_INT, val_min.rank, cs_glob_mpi_comm);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parfpt\n",
         cs_glob_rank_id, n_parfpt_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Return the value associated to a probe.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARHIS (NODE, NDRANG, VAR, VARCAP)
 * *****************
 *
 * INTEGER          NODE        : --> : local number of the element related to
 *                                      a measure node
 * INTEGER          NDRANG      : --> : rank of the process owning the closest
 *                                      node from the measure node
 * DOUBLE PRECISION VAR(*)      : --> : values of the variable on local elements
 * DOUBLE PRECISION VARCAP      : <-- : value of the variable for the element
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

#if CS_PARALL_DEBUG_COUNT
  printf ("irang = %d, iappel = %d, tot = %d, parhis\n",
          cs_glob_rank_id, n_parhis_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Return the global cell number of a local cell.
 * (send to all the processes)
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARCEL (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM        : --> : local cell number
 * INTEGER          RANKID      : --> : rank of the domain (0 to N-1)
 * INTEGER          GNUM        : <-- : global cell number
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcel, PARCEL)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum)
{
#if defined(HAVE_MPI)

  assert(sizeof(double) == sizeof(cs_real_t));

  if (*rankid == cs_glob_rank_id)
    *gnum = cs_glob_mesh->global_cell_num[*lnum - 1];
  else
    *gnum = 0;

  MPI_Bcast(gnum, 1, CS_MPI_INT, *rankid, cs_glob_mpi_comm);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, parcel\n",
         cs_glob_rank_id, n_parcel_calls++, n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------
 * Return the global cell number knowing its related local cell number. No
 * communication is useful.
 * Return the local cell number in serial mode.
 * Return 0 if the local cell number > mesh->n_cells
 * Return 0 if the current rank domain != RANKID
 *
 * Fortran interface :
 *
 * SUBROUTINE PARCLG (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM      : --> : local cell number
 * INTEGER          RANKID    : --> : rank of the current domain (0 to N-1)
 * INTEGER          GNUM      : <-- : global cell number
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parclg, PARCLG)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum)
{
  if (*rankid < 0)
    *gnum = *lnum;

  else if ( (*rankid == cs_glob_rank_id) && (*lnum <= cs_glob_mesh->n_cells) )
    *gnum = cs_glob_mesh->global_cell_num[*lnum - 1];

  else
    *gnum = 0;

}

/*----------------------------------------------------------------------------
 * Return the global internal face number knowing its related local internal
 * face number. No communication is useful.
 * Return the local internal face number in serial mode.
 * Return 0 if the local internal face number > mesh->n_i_faces
 * Return 0 if the current rank domain != RANKID
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARFIG (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM      : --> : local internal face number
 * INTEGER          RANKID    : --> : rank of the current domain (0 to N-1)
 * INTEGER          GNUM      : <-- : global internal face number
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfig, PARFIG)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum)
{
  if (*rankid < 0)
    *gnum = *lnum;

  else if (   (*rankid == cs_glob_rank_id)
           && (*lnum <= cs_glob_mesh->n_i_faces) )
    *gnum = cs_glob_mesh->global_i_face_num[*lnum - 1];

  else
    *gnum = 0;

}

/*----------------------------------------------------------------------------
 * Return the global border face number knowing its related local border face
 * number. No communication is useful.
 * Return the local border face number in serial mode.
 * Return 0 if the local border face number > mesh->n_b_faces
 * Return 0 if the current rank domain != RANKID
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARFBG (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM      : --> : local border face number
 * INTEGER          RANKID    : --> : rank of the current domain (0 to N-1)
 * INTEGER          GNUM      : <-- : global border face number
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfbg, PARFBG)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum)
{
  if (*rankid < 0)
    *gnum = *lnum;

  else if (   (*rankid == cs_glob_rank_id)
           && (*lnum <= cs_glob_mesh->n_b_faces) )
    *gnum = cs_glob_mesh->global_b_face_num[*lnum - 1];

  else
    *gnum = 0;

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the sum of real values for entities belonging to a fvm_interface_t
 * structure.
 *
 * Only the values of entities belonging to the interface are summed.
 *
 * parameters:
 *   interfaces --> pointer to a fvm_interface_set_t structure
 *   var_size   --> number of elements in var buffer
 *   stride     --> number of values (no interlaced) by entity
 *   var        <-> variable buffer
 *----------------------------------------------------------------------------*/

void
cs_parall_interface_sr(fvm_interface_set_t  *interfaces,
                       cs_int_t              var_size,
                       cs_int_t              stride,
                       cs_real_t            *var)
{
#if defined(HAVE_MPI)

  int  request_count;
  int  distant_rank, n_interfaces;
  cs_int_t  id, ii, jj;
  cs_int_t  total_size;

  cs_int_t  count_size = 0;
  fvm_lnum_t  n_entities = 0;
  cs_real_t  *buf = NULL, *send_buf = NULL, *recv_buf = NULL;

  MPI_Request  *request = NULL;
  MPI_Status  *status  = NULL;

  const fvm_lnum_t  *local_num = NULL;
  const fvm_interface_t  *interface = NULL;

  /* Initialize and allocate */

  n_interfaces = fvm_interface_set_size(interfaces);

  for (id = 0; id < n_interfaces; id++) {
    count_size
      += fvm_interface_size(fvm_interface_set_get(interfaces, id));
  }

  total_size = count_size;

  BFT_MALLOC(buf, total_size * stride * 2, cs_real_t);

  BFT_MALLOC(request, n_interfaces * 2, MPI_Request);
  BFT_MALLOC(status,  n_interfaces * 2, MPI_Status);

  /* Send and Receive data from distant ranks with
     non-blocking communications */

  request_count = 0;
  count_size  = 0;

  /* Receive */

  for (id = 0; id < n_interfaces; id++) {

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);

    recv_buf = buf + (count_size * stride);

    MPI_Irecv(recv_buf,
              n_entities * stride,
              CS_MPI_REAL,
              distant_rank,
              distant_rank,
              cs_glob_mpi_comm,
              &(request[request_count++]));

    count_size += n_entities;

  }

  assert(count_size == total_size);

  /* Send */

  for (id = 0 ; id < n_interfaces ; id++) {

    /* Preparation of data to send */

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank   = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);
    local_num = fvm_interface_get_local_num(interface);

    send_buf = buf + (count_size * stride);

    for (ii = 0 ; ii < n_entities ; ii++) {
      for (jj = 0 ; jj < stride ; jj++)
        send_buf[ii*stride + jj] = var[jj*var_size + (local_num[ii] - 1)];
    }

    MPI_Isend(send_buf,
              n_entities * stride,
              CS_MPI_REAL,
              distant_rank,
              (int)cs_glob_rank_id,
              cs_glob_mpi_comm,
              &(request[request_count++]));

    count_size += n_entities;

  }

  assert(count_size == 2*total_size);

  /* Sync after each rank had received all the messages */

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  /* Now we had each part to var */

  count_size = 0;

  for (id = 0 ; id < n_interfaces ; id++) {

    /* Retrieve data */

    interface = fvm_interface_set_get(interfaces, id);
    n_entities = fvm_interface_size(interface);
    local_num = fvm_interface_get_local_num(interface);

    recv_buf = buf + (count_size * stride);

    for (ii = 0 ; ii < n_entities ; ii++) {
      for (jj = 0 ; jj < stride ; jj++) {
        var[jj*var_size + (local_num[ii] - 1)] += recv_buf[ii*stride + jj];
      }
    }

    count_size += n_entities;

  }

  BFT_FREE(buf);

#endif

#if CS_PARALL_DEBUG_COUNT
  printf("irang = %d, iappel = %d, tot = %d, cs_parallel_interface_sr\n",
         cs_glob_rank_id,
         n_interface_sr_calls++,
         n_total_par_calls++);
#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
