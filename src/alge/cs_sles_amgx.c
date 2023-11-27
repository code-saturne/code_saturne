/*============================================================================
 * Sparse Linear Equation Solvers using AmgX
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * AmgX headers
 *----------------------------------------------------------------------------*/

#include <amgx_c.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_log.h"
#include "cs_fp_exception.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_sles_amgx.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_amgx.c

  \brief handling of AmgX-based linear solvers

  \page sles_amgx AmgX-based linear solvers.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_sles_amgx_setup_t {

  AMGX_solver_handle   solver;           /* Linear solver context */
  AMGX_matrix_handle   matrix;           /* Linear system matrix */

  double               r_norm;           /* residual normalization */
  void                 *cctx;            /* convergence context */

} cs_sles_amgx_setup_t;

struct _cs_sles_amgx_t {

  /* Performance data */

  int                  n_setups;           /* Number of times system setup */
  int                  n_solves;           /* Number of times system solved */

  int                  n_iterations_last;  /* Number of iterations for last
                                              system resolution */
  int                  n_iterations_min;   /* Minimum number of iterations
                                              in system resolution history */
  int                  n_iterations_max;   /* Maximum number of iterations
                                              in system resolution history */
  int long long        n_iterations_tot;   /* Total accumulated number of
                                              iterations */

  cs_timer_counter_t   t_setup;            /* Total setup */
  cs_timer_counter_t   t_solve;            /* Total time used */

  /* Setup data */

  char                   *amgx_config_file;
  char                   *amgx_config_string;

  AMGX_Mode               amgx_mode;

  int                     flags;              /* additional option flags */

  AMGX_config_handle      amgx_config;        /* Solver configuration */
  AMGX_resources_handle   amgx_resources;     /* Associated resources */

  cs_sles_amgx_setup_t   *setup_data;

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static int  _n_amgx_systems = 0;

#if defined(HAVE_MPI)
static MPI_Comm _amgx_comm = MPI_COMM_NULL;
#endif

/* TODO: for multi-device configurations, this will need to be adapted */

static int _n_devices = 1;
static int _devices[] = {0};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print function for AmgX.
 *
 * \param[in]  msg     message to print
 * \param[in]  length  length of message to print
 */
/*----------------------------------------------------------------------------*/

static void
_print_callback(const char  *msg,
                int          length)
{
  CS_NO_WARN_IF_UNUSED(length);

  bft_printf("%s", msg);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize AmgX.
 */
/*----------------------------------------------------------------------------*/

static void
_amgx_initialize(void)
{
  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");
  const char warning_fmt[] = N_("\nwarning: %s returned %d.\n"
                                "%s\n");

  AMGX_RC retval = AMGX_RC_OK;

  retval = AMGX_register_print_callback(_print_callback);
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_printf(_(warning_fmt), "AMGX_initialize", (int)retval, err_str);
  }

  retval = AMGX_initialize();
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_initialize", retval, err_str);
  }

  retval = AMGX_initialize_plugins();
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_initialize_plugins", retval, err_str);
  }

  int major, minor;
  AMGX_get_api_version(&major, &minor);
  bft_printf(_("\nAMGX API version %d.%d\n"), major, minor);

  /* Note: if MPI supports GPUDirect, MPI_DIRECT is also allowed */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    MPI_Comm_dup(cs_glob_mpi_comm, &_amgx_comm);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize AmgX.
 */
/*----------------------------------------------------------------------------*/

static void
_amgx_finalize(void)
{
  char err_str[4096];
  const char warning_fmt[] = N_("\nwarning: %s returned %d.\n"
                                "%s\n");

  AMGX_RC retval = AMGX_finalize_plugins();
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_printf(_(warning_fmt), "AMGX_finalize_plugins", (int)retval, err_str);
  }

  retval = AMGX_finalize();
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_printf(_(warning_fmt), "AMGX_finallize", (int)retval, err_str);
  }

#if defined(HAVE_MPI)
  if (_amgx_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&_amgx_comm);
    _amgx_comm = MPI_COMM_NULL;
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Load AmgX solver configuration.
 *
 * \param[in, out]  c   pointer to AmgX solver info and context
 */
/*----------------------------------------------------------------------------*/

static void
_load_amgx_config(cs_sles_amgx_t  *c)
{
  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");
  AMGX_RC retval = AMGX_RC_OK;

  if (c->amgx_config_file == NULL) {
    retval = AMGX_config_create(&(c->amgx_config),
                                cs_sles_amgx_get_config(c));
    if (retval != AMGX_RC_OK) {
      AMGX_get_error_string(retval, err_str, 4096);
      bft_error(__FILE__, __LINE__, 0, _(error_fmt),
                "AMGX_config_create", retval, err_str);
    }
  }
  else {
    retval = AMGX_config_create_from_file(&(c->amgx_config),
                                          c->amgx_config_file);
    if (retval != AMGX_RC_OK) {
      AMGX_get_error_string(retval, err_str, 4096);
      bft_error(__FILE__, __LINE__, 0, _(error_fmt),
                "AMGX_config_create_from_file", retval, err_str);
    }
  }

  retval = AMGX_config_add_parameters(&(c->amgx_config),
                                      "config_version=2, "
                                      "main:monitor_residual=1, "
                                      "main:store_res_history=1");
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_config_add_parameters", retval, err_str);
  }

#if 0
  /* Exception handling can be ensured by AMGX, but by default,
     we prefer to check return codes and used the regular code_saturne
     error handling. */
  AMGX_config_add_parameters(&(c->amgx_config), "exception_handling=1");
#endif

  void *comm_ptr = NULL;
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    comm_ptr = &_amgx_comm;
#endif

  retval = AMGX_resources_create(&c->amgx_resources, c->amgx_config,
                                 comm_ptr,
                                 _n_devices, _devices);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_resources_create", retval, err_str);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup AmgX matrix with at most 1 ring.
 *
 * \param[in, out]  c   pointer to AmgX solver info and context
 * \param[in]       a   associated matrix
 */
/*----------------------------------------------------------------------------*/

static void
_setup_matrix_1_ring(cs_sles_amgx_t     *c,
                     const cs_matrix_t  *a)
{
  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");

  cs_sles_amgx_setup_t *sd = c->setup_data;

  if (sizeof(cs_lnum_t) != sizeof(int))
    bft_error
      (__FILE__, __LINE__, 0,
       _("AmgX bindings are not currently handled for code_saturne builds\n"
         "using long cs_lnum_t types (i.e. --enable-long-lnum)."));

  const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
  const int n_rows = cs_matrix_get_n_rows(a);
  const int db_size = cs_matrix_get_diag_block_size(a);
  const cs_halo_t *halo = cs_matrix_get_halo(a);

  const cs_lnum_t *a_row_index, *a_col_id;

  const cs_real_t *a_val = NULL, *a_d_val = NULL;

  if (cs_mat_type == CS_MATRIX_CSR)
    cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);
  else if (cs_mat_type == CS_MATRIX_MSR)
    cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_val);

  const cs_lnum_t  *row_index = a_row_index;
  int              *_row_index = NULL;
  const cs_lnum_t  *col_id = a_col_id;
  cs_lnum_t        *_col_id = NULL;

  if (sizeof(int) != sizeof(cs_lnum_t)) {
    BFT_MALLOC(_row_index, n_rows, int);
    for (cs_lnum_t i = 0; i < n_rows; i++)
      _row_index[i] = a_row_index[i];
    row_index = _row_index;
    int nnz = row_index[n_rows];
    BFT_MALLOC(_col_id, nnz, int);
    for (cs_lnum_t i = 0; i < nnz; i++)
      _col_id[i] = a_col_id[i];
    col_id = _col_id;
  }

  /* Matrix */

  AMGX_RC retval;

  retval = AMGX_matrix_create(&(sd->matrix),
                              c->amgx_resources,
                              c->amgx_mode);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_matrix_create", retval, err_str);
  }

  if (cs_glob_n_ranks > 1) {

    int *send_sizes, *recv_sizes;
    BFT_MALLOC(send_sizes, halo->n_c_domains, int);
    BFT_MALLOC(recv_sizes, halo->n_c_domains, int);
    for (int i = 0; i < halo->n_c_domains; i++) {
      send_sizes[i] =   halo->send_index[2*i + 1]
                      - halo->send_index[2*i];
      recv_sizes[i] =   halo->index[2*i + 1]
                      - halo->index[2*i];
    }
    int **send_maps, **recv_maps;
    BFT_MALLOC(send_maps, halo->n_c_domains, int *);
    BFT_MALLOC(recv_maps, halo->n_c_domains, int *);

    assert(sizeof(cs_lnum_t) == sizeof(int));

    for (int i = 0; i < halo->n_c_domains; i++) {
      BFT_MALLOC(send_maps[i], send_sizes[i], int);
      int *_send_map = send_maps[i];
      for (int j = 0; j < send_sizes[i]; j++)
        _send_map[j] = halo->send_list[halo->send_index[2*i] + j];
      BFT_MALLOC(recv_maps[i], recv_sizes[i], int);
      int *_recv_map = recv_maps[i];
      for (int j = 0; j < recv_sizes[i]; j++)
        _recv_map[j] = halo->n_local_elts + halo->index[2*i] + j;
    }

    retval = AMGX_matrix_comm_from_maps_one_ring(sd->matrix,
                                                 1, /* allocated_halo_depth */
                                                 halo->n_c_domains,
                                                 halo->c_domain_rank,
                                                 send_sizes,
                                                 (const int **)send_maps,
                                                 recv_sizes,
                                                 (const int **)recv_maps);

    if (retval != AMGX_RC_OK) {
      AMGX_get_error_string(retval, err_str, 4096);
      bft_error(__FILE__, __LINE__, 0, _(error_fmt),
                "AMGX_matrix_comm_from_maps_one_ring", retval, err_str);
    }

    for (int i = 0; i < halo->n_c_domains; i++) {
      BFT_FREE(recv_maps[i]);
      BFT_FREE(send_maps[i]);
    }
    BFT_FREE(recv_sizes);
    BFT_FREE(send_sizes);

  }

  const int b_size = cs_matrix_get_diag_block_size(a);
  const int b_mem_size = b_size*b_size*sizeof(cs_real_t);

  cs_alloc_mode_t amode_row_index = cs_check_device_ptr(row_index);
  cs_alloc_mode_t amode_col_id = cs_check_device_ptr(col_id);
  cs_alloc_mode_t amode_a_val = cs_check_device_ptr(a_val);
  cs_alloc_mode_t amode_a_d_val = cs_check_device_ptr(a_d_val);

  if (amode_row_index < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)row_index, (n_rows+1)*sizeof(int));
  if (amode_col_id < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)col_id, a_row_index[n_rows]*sizeof(int));
  if (amode_a_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)a_val, a_row_index[n_rows]*b_mem_size);
  if (a_d_val != NULL && amode_a_d_val < CS_ALLOC_HOST_DEVICE_PINNED)
      AMGX_pin_memory((void *)a_d_val, n_rows*b_mem_size);

  retval = AMGX_matrix_upload_all(sd->matrix,
                                  n_rows,
                                  cs_matrix_get_n_entries(a),
                                  db_size,
                                  db_size,
                                  row_index,
                                  col_id,
                                  a_val,
                                  a_d_val);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_matrix_upload_all", retval, err_str);
  }

  if (a_d_val != NULL && amode_a_d_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)a_d_val);
  if (amode_a_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)a_val);
  if (amode_col_id < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)col_id);
  if (amode_row_index < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)row_index);

  BFT_FREE(_row_index);
  BFT_FREE(_col_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup AmgX matrix based on distribution info
 *
 * \param[in, out]  c   pointer to AmgX solver info and context
 * \param[in]       a   associated matrix
 */
/*----------------------------------------------------------------------------*/

static void
_setup_matrix_dist(cs_sles_amgx_t     *c,
                   const cs_matrix_t  *a)
{
  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");

  cs_sles_amgx_setup_t *sd = c->setup_data;

  if (sizeof(cs_lnum_t) != sizeof(int))
    bft_error
      (__FILE__, __LINE__, 0,
       _("AmgX bindings are not currently handled for code_saturne builds\n"
         "using long cs_lnumt_ types (i.e. --enable-long-lnum)."));

  const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
  const int n_rows = cs_matrix_get_n_rows(a);
  const int db_size = cs_matrix_get_diag_block_size(a);
  const cs_halo_t *halo = cs_matrix_get_halo(a);

  const cs_lnum_t *a_row_index, *a_col_id;

  const cs_real_t *a_val = NULL, *a_d_val = NULL;

  if (cs_mat_type == CS_MATRIX_CSR)
    cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);
  else if (cs_mat_type == CS_MATRIX_MSR)
    cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_val);

  const cs_lnum_t  *row_index = a_row_index;
  int              *_row_index = NULL;
  cs_gnum_t        *col_gid = NULL;

  cs_alloc_mode_t amode = CS_ALLOC_HOST_DEVICE;

  if (sizeof(int) != sizeof(cs_lnum_t)) {
    CS_MALLOC_HD(_row_index, n_rows, int, amode);
    for (cs_lnum_t i = 0; i < n_rows; i++)
      _row_index[i] = a_row_index[i];
    row_index = _row_index;
  }

  cs_lnum_t nnz = row_index[n_rows];
  CS_MALLOC_HD(col_gid, nnz, cs_gnum_t, amode);

  const cs_gnum_t *grow_id = cs_matrix_get_block_row_g_id(a);

  for (cs_lnum_t j = 0; j < n_rows; j++) {
    for (cs_lnum_t i = a_row_index[j]; i < a_row_index[j+1]; ++i)
      col_gid[i] = grow_id[a_col_id[i]];
  }

  /* Matrix */

  AMGX_RC retval;

  retval = AMGX_matrix_create(&(sd->matrix),
                              c->amgx_resources,
                              c->amgx_mode);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_matrix_create", retval, err_str);
  }

  /* Assume partitioning with continuous blocks */

  cs_gnum_t *partition_offsets;
  CS_MALLOC_HD(partition_offsets, cs_glob_n_ranks + 1, cs_gnum_t, amode);

  /* Gather the number of rows on each rank, use exclusive scan to get the offsets */

  int n_ranks = cs_glob_n_ranks;

  cs_gnum_t n_g_rows = n_rows;
  partition_offsets[0] = 0;
  MPI_Allgather(&n_g_rows, 1, CS_MPI_GNUM, &partition_offsets[1], 1, CS_MPI_GNUM,
                _amgx_comm);
  for (cs_lnum_t i = 2; i < n_ranks + 1; i++) {
    partition_offsets[i] += partition_offsets[i-1];
  }
  n_g_rows = partition_offsets[n_ranks];

  const int b_size = cs_matrix_get_diag_block_size(a);
  const int b_mem_size = b_size*b_size*sizeof(cs_real_t);

  cs_alloc_mode_t amode_row_index = cs_check_device_ptr(row_index);
  cs_alloc_mode_t amode_a_val = cs_check_device_ptr(a_val);
  cs_alloc_mode_t amode_a_d_val = cs_check_device_ptr(a_d_val);

  if (amode_row_index < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)row_index, (n_rows+1)*sizeof(int));
  if (amode < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory(col_gid, a_row_index[n_rows]*sizeof(int));
  if (amode_a_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)a_val, a_row_index[n_rows]*b_mem_size);
  if (a_d_val != NULL && amode_a_d_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)a_d_val, n_rows*b_mem_size);
  if (amode < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory(partition_offsets, (n_ranks+1)*sizeof(cs_gnum_t));

  AMGX_distribution_handle dist;
  AMGX_distribution_create(&dist, c->amgx_config);

  if (sizeof(cs_gnum_t) == sizeof(int32_t))
    AMGX_distribution_set_32bit_colindices(dist, 1);

  retval = AMGX_distribution_set_partition_data(dist,
                                                AMGX_DIST_PARTITION_OFFSETS,
                                                partition_offsets);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_distribution_set_partition_data", retval, err_str);
  }

  retval = AMGX_matrix_upload_distributed(sd->matrix,
                                          n_g_rows, n_rows, nnz,
                                          db_size, db_size,
                                          row_index, col_gid,
                                          a_val, a_d_val,
                                          dist);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_matrix_upload_distributed", retval, err_str);
  }

  AMGX_distribution_destroy(dist);

  if (amode < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory(partition_offsets);
  if (a_d_val != NULL && amode_a_d_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)a_d_val);
  if (amode_a_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)a_val);
  if (amode < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory(col_gid);
  if (amode_row_index < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)row_index);

  CS_FREE_HD(partition_offsets);
  CS_FREE_HD(col_gid);
  CS_FREE_HD(_row_index);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update AmgX setup with new matrix coefficients for same structure.
 *
 * \param[in, out]  c   pointer to AmgX solver info and context
 * \param[in]       a   associated matrix
 */
/*----------------------------------------------------------------------------*/

static void
_setup_update_coeffs(cs_sles_amgx_t     *c,
                     const cs_matrix_t  *a)
{
  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");

  cs_sles_amgx_setup_t *sd = c->setup_data;

  const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
  const int n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t *a_row_index, *a_col_id;

  const cs_real_t *a_val = NULL, *a_d_val = NULL;

  if (cs_mat_type == CS_MATRIX_CSR)
    cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);
  else if (cs_mat_type == CS_MATRIX_MSR)
    cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_val);

  /* Matrix */

  AMGX_RC retval;

  const int b_size = cs_matrix_get_diag_block_size(a);
  const int b_mem_size = b_size*b_size*sizeof(cs_real_t);

  cs_alloc_mode_t amode_a_val = cs_check_device_ptr(a_val);
  cs_alloc_mode_t amode_a_d_val = cs_check_device_ptr(a_d_val);
  if (amode_a_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)a_val, a_row_index[n_rows]*b_mem_size);
  if (a_d_val != NULL && amode_a_d_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)a_d_val, n_rows*b_mem_size);

  retval = AMGX_matrix_replace_coefficients(sd->matrix,
                                            n_rows,
                                            cs_matrix_get_n_entries(a),
                                            a_val,
                                            a_d_val);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_matrix_replace_coefficients", retval, err_str);
  }

  if (a_d_val != NULL && amode_a_d_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)a_d_val);
  if (amode_a_val < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)a_val);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate an AmgX linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_amgx_create.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * Note that this function returns a pointer directly to the solver
 * management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or NULL
 *
 * \return  pointer to newly created AmgX solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_amgx_t *
cs_sles_amgx_define(int           f_id,
                     const char  *name)
{
  cs_sles_amgx_t * c = cs_sles_amgx_create();

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 c,
                                 "cs_sles_amgx_t",
                                 cs_sles_amgx_setup,
                                 cs_sles_amgx_solve,
                                 cs_sles_amgx_free,
                                 cs_sles_amgx_log,
                                 cs_sles_amgx_copy,
                                 cs_sles_amgx_destroy);

  CS_NO_WARN_IF_UNUSED(sc);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create AmgX linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * \return  pointer to associated linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_amgx_t *
cs_sles_amgx_create(void)

{
  cs_sles_amgx_t *c;

  if (_n_amgx_systems < 1)
    _amgx_initialize();

  _n_amgx_systems += 1;

  BFT_MALLOC(c, 1, cs_sles_amgx_t);
  c->n_setups = 0;
  c->n_solves = 0;
  c->n_iterations_last = 0;
  c->n_iterations_min = 0;
  c->n_iterations_max = 0;
  c->n_iterations_tot = 0;

  CS_TIMER_COUNTER_INIT(c->t_setup);
  CS_TIMER_COUNTER_INIT(c->t_solve);

  /* Setup data */

  c->setup_data = NULL;

  c->amgx_config_file = NULL;
  c->amgx_config_string = NULL;

  cs_sles_amgx_set_use_device(c, true);

  c->flags = 0;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create AmgX linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_amgx_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_amgx_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_amgx_copy(const void  *context)
{
  cs_sles_amgx_t *d = NULL;

  if (context != NULL) {
    const cs_sles_amgx_t *c = context;
    d = cs_sles_amgx_create();

    if (c->amgx_config_file != NULL) {
      size_t l = strlen(c->amgx_config_file);
      BFT_MALLOC(d->amgx_config_file, l+1, char);
      strncpy(d->amgx_config_file, c->amgx_config_file, l);
      d->amgx_config_file[l] = '\0';
    }
    if (c->amgx_config_string != NULL) {
      size_t l = strlen(c->amgx_config_string);
      BFT_MALLOC(d->amgx_config_string, l+1, char);
      strncpy(d->amgx_config_string, c->amgx_config_string, l);
      d->amgx_config_string[l] = '\0';
    }

    d->amgx_mode = c->amgx_mode;
    d->flags = c->flags;
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy AmgX linear system solver info and context.
 *
 * \param[in, out]  context  pointer to AmgX solver info and context
 *                           (actual type: cs_sles_amgx_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_destroy(void **context)
{
  cs_sles_amgx_t *c = (cs_sles_amgx_t *)(*context);
  if (c != NULL) {

    /* Free local strings */

    BFT_FREE(c->amgx_config_file);
    BFT_FREE(c->amgx_config_string);

    if (c->n_setups >= 1) {
      char err_str[4096];
      const char warning_fmt[] = N_("\nwarning: %s returned %d.\n"
                                    "%s\n");

      AMGX_RC retval = AMGX_resources_destroy(c->amgx_resources);
      if (retval != AMGX_RC_OK) {
        AMGX_get_error_string(retval, err_str, 4096);
        bft_printf(_(warning_fmt), "AMGX_resources_destroy",
                   (int)retval, err_str);
      }

      retval = AMGX_config_destroy(c->amgx_config);
      if (retval != AMGX_RC_OK) {
        AMGX_get_error_string(retval, err_str, 4096);
        bft_printf(_(warning_fmt), "AMGX_resources_destroy",
                   (int)retval, err_str);
      }
    }

    /* Free structure */

    cs_sles_amgx_free(c);
    BFT_FREE(c);
    *context = c;

    _n_amgx_systems -= 1;
    if (_n_amgx_systems == 0)
      _amgx_finalize();

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief return the configuration for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration strings syntax.
 *
 * \param[in, out]  context  pointer to AmgX solver info and context
 *
 * \return  configuration string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_amgx_get_config(void  *context)
{
  cs_sles_amgx_t  *c = context;

  if (   c->amgx_config_file == NULL
      && c->amgx_config_string == NULL) {

    const char config[] =
      "{"
      "  \"config_version\": 2, "
      "  \"solver\": {"
      "    \"preconditioner\": {"
      "      \"print_grid_stats\": 1, "
      "      \"print_vis_data\": 0, "
      "      \"solver\": \"AMG\", "
      "      \"print_solve_stats\": 0, "
      "      \"interpolator\": \"D2\", "
      "      \"presweeps\": 1, "
      "      \"max_iters\": 1, "
      "      \"monitor_residual\": 0, "
      "      \"store_res_history\": 0, "
      "      \"scope\": \"amg\", "
      "      \"cycle\": \"V\", "
      "      \"postsweeps\": 1 "
      "    }, "
      "    \"solver\": \"FGMRES\", "
      "    \"print_solve_stats\": 0, "
      "    \"solver_verbose\": 0, "
      "    \"obtain_timings\": 0, "
      "    \"max_iters\": 100, "
      "    \"gmres_n_restart\": 20, "
      "    \"convergence\": \"ABSOLUTE\", "
      "    \"scope\": \"main\", "
      "    \"tolerance\": 0.0001, "
      "    \"norm\": \"L2\""
      "  }"
      "}";

    cs_sles_amgx_set_config(context, config);

  }

  return c->amgx_config_string;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the solver configuration for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration strings syntax.
 *
 * If this function is not called, a default configuration will be used.
 *
 * \param[in, out]  context  pointer to AmgX solver info and context
 * \param[in]       config   string defining configuration to use
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_set_config(void        *context,
                        const char  *config)
{
  cs_sles_amgx_t  *c = context;

  size_t l = strlen(config);

  BFT_REALLOC(c->amgx_config_string, l+1, char);
  strncpy(c->amgx_config_string, config, l);
  c->amgx_config_string[l] = '\0';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief return the name of the solver configuration file for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration file syntax.
 *
 * \param[in, out]  context  pointer to AmgX solver info and context
 *
 * \return  configuration file name, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_amgx_get_config_file(void  *context)
{
  cs_sles_amgx_t  *c = context;

  return c->amgx_config_file;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the solver configuration file for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration file syntax.
 *
 * If this function is not called, a default configuration will be used.
 *
 * \param[in, out]  context  pointer to AmgX solver info and context
 * \param[in]       path     path to configuration file
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_set_config_file(void        *context,
                             const char  *path)
{
  cs_sles_amgx_t  *c = context;

  size_t l = strlen(path);

  BFT_REALLOC(c->amgx_config_file, l+1, char);
  strncpy(c->amgx_config_file, path, l);
  c->amgx_config_file[l] = '\0';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query whether an AmgX solver should use the device or host.
 *
 * \param[in]  context  pointer to AmgX solver info and context
 *
 * \return  true for device, false for host only
 */
/*----------------------------------------------------------------------------*/

bool
cs_sles_amgx_get_use_device(void  *context)
{
  cs_sles_amgx_t  *c = context;
  bool use_device = true;

  if (   c->amgx_mode == AMGX_mode_hDDI
      || c->amgx_mode == AMGX_mode_hFFI) {
    use_device = false;
  }

  return use_device;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define whether an AmgX solver should use the device or host.
 *
 * By default, the device will be used, but by calling this function
 * with "use_device = false", only the host will be used.
 *
 * \param[in, out]  context      pointer to AmgX solver info and context
 * \param[in]       use_device   true for devince, false for host only
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_set_use_device(void  *context,
                            bool   use_device)
{
  cs_sles_amgx_t  *c = context;

  if (use_device) {
    if (sizeof(cs_real_t) == sizeof(double))
      c->amgx_mode = AMGX_mode_dDDI;
    else if (sizeof(cs_real_t) == sizeof(float))
      c->amgx_mode = AMGX_mode_dFFI;
  }

  else {  /* To run on host instead of device */
    if (sizeof(cs_real_t) == sizeof(double))
      c->amgx_mode = AMGX_mode_hDDI;
    else if (sizeof(cs_real_t) == sizeof(float))
      c->amgx_mode = AMGX_mode_hFFI;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query additional AmgX solver usage flags.
 *
 * \param[in]  context  pointer to AmgX solver info and context
 *
 * \return  associated flags
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_amgx_get_flags(void  *context)
{
  cs_sles_amgx_t  *c = context;
  return c->flags;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define additional AmgX solver usage flags
 *
 * By default, the device will be used, but by calling this function
 * with "use_device = false", only the host will be used.
 *
 * \param[in, out]  context   pointer to AmgX solver info and context
 * \param[in]       flags     flags (sum/bitwise of) for AmgX usage options.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_set_flags(void  *context,
                       int    flags)
{
  cs_sles_amgx_t  *c = context;
  c->flags = flags;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup AmgX linear equation solver.
 *
 * \param[in, out]  context    pointer to AmgX solver info and context
 *                             (actual type: cs_sles_amgx_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_setup(void               *context,
                   const char         *name,
                   const cs_matrix_t  *a,
                   int                 verbosity)
{
  CS_NO_WARN_IF_UNUSED(verbosity);

  cs_timer_t t0;
  t0 = cs_timer_time();

  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");

  cs_sles_amgx_t  *c = context;
  cs_sles_amgx_setup_t *sd = c->setup_data;

  AMGX_RC retval = AMGX_RC_OK;

  /* Case where setup data is not currently present */

  if (sd == NULL) {

    BFT_MALLOC(c->setup_data, 1, cs_sles_amgx_setup_t);
    sd = c->setup_data;

    /* Load configuration at first call */

    if (c->n_setups < 1)
      _load_amgx_config(c);

    const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
    const int db_size = cs_matrix_get_diag_block_size(a);
    const cs_halo_t *halo = cs_matrix_get_halo(a);

    /* Periodicity is not handled (at least not in serial mode), as the matrix
       is not square due to ghost cells */

    if (halo != NULL) {
      bool have_perio = false;
      if (halo->n_transforms > 0)
        have_perio = true;
      assert(have_perio == false);
    }

    /* TODO: handle periodicity, by renumbering local periodic cells
       so as to use the main (and not ghost) cell id */

    if (   db_size > 1
        || (cs_mat_type != CS_MATRIX_CSR && cs_mat_type != CS_MATRIX_MSR))
      bft_error
        (__FILE__, __LINE__, 0,
         _("Matrix type %s with block size %d for system \"%s\" "
           "is not usable by AmgX.\n"
           "Only block size 1 with CSR or MSR type "
           "is currently supported by AmgX."),
         cs_matrix_get_type_name(a), db_size,
         name);

    /* Number of ghost cell layers depends on solver configuration
       (classical seems to require 2 rings, aggregation 1 ring,
       but it seems AMGX_config_get_default_number_of_rings cannot
       always be trusted, so we need a user-defined flag to allow
       forcing a setting in case of incorrect compatibility detection;
       we could also use the distribution systematically, but prefer
       to have both options for safety) */

    bool use_dist = false;
    if (cs_glob_n_ranks > 1) {
      if (c->flags & CS_SLES_AMGX_PREFER_COMM_FROM_MAPS) {
        int n_rings = 1;
        AMGX_config_get_default_number_of_rings(c->amgx_config, &n_rings);
        if (n_rings > 1)
          use_dist = true;
      }
      else
        use_dist = true;
    }

    if (use_dist == false)
      _setup_matrix_1_ring(c, a);

    else
      _setup_matrix_dist(c, a);

    /* Solver */

    retval = AMGX_solver_create(&(sd->solver),
                                c->amgx_resources,
                                c->amgx_mode,
                                c->amgx_config);

    if (retval != AMGX_RC_OK) {
      AMGX_get_error_string(retval, err_str, 4096);
      bft_error(__FILE__, __LINE__, 0, _(error_fmt),
                "AMGX_solver_create", retval, err_str);
    }

  }

  /* Case where setup data is already present (simple update) */

  else {
    _setup_update_coeffs(c, a);
  }

  retval = AMGX_solver_setup(sd->solver, sd->matrix);

  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_solver_setup", retval, err_str);
  }

  sd->r_norm = -1;

  /* Update return values */
  c->n_setups += 1;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);

#if 0
  /* Check matrix coefficients are correctly set by comparing
     spmv with reference (debug) */
  {
    const int n_rows = cs_matrix_get_n_rows(a);
    const int n_cols_ext = cs_matrix_get_n_columns(a);
    const int  b_size = cs_matrix_get_diag_block_size(a);
    const cs_lnum_t n_r = n_rows*b_size;
    const cs_lnum_t n_c = n_cols_ext*b_size;

    cs_real_t *x, *y, *z;
    CS_MALLOC_HD(x, n_c, cs_real_t, CS_ALLOC_HOST_DEVICE_PINNED);
    CS_MALLOC_HD(y, n_c, cs_real_t, CS_ALLOC_HOST_DEVICE_PINNED);
    CS_MALLOC_HD(z, n_c, cs_real_t, CS_ALLOC_HOST);

    for (cs_lnum_t ii = 0; ii < n_cols_ext; ii++) {
      cs_gnum_t jj = ii*b_size;
      for (cs_lnum_t kk = 0; kk < b_size; kk++) {
        x[ii*b_size+kk] = sin(jj*b_size+kk);
        y[ii*b_size+kk] = 0.0;
      }
    }

    /* Vector */

    AMGX_vector_handle  hx, hy;

    AMGX_vector_create(&hx, c->amgx_resources, c->amgx_mode);
    AMGX_vector_create(&hy, c->amgx_resources, c->amgx_mode);

    if (cs_glob_n_ranks > 1) {
      AMGX_vector_bind(hx, c->setup_data->matrix);
      AMGX_vector_bind(hy, c->setup_data->matrix);
    }

    AMGX_vector_upload(hx, n_rows, b_size, x);
    AMGX_vector_upload(hy, n_rows, b_size, y);
    AMGX_matrix_vector_multiply(sd->matrix, hx, hy);
    AMGX_vector_download(hy, y);

    cs_matrix_vector_multiply(a, x, z);

    CS_FREE_HD(x);
    AMGX_vector_destroy(hx);

    double dmax = 0.0;
    for (cs_lnum_t ii = 0; ii < n_r; ii++) {
      double d = CS_ABS(y[ii] - z[ii]);
      if (d > dmax)
        dmax = d;
    }

    CS_FREE_HD(y);
    AMGX_vector_destroy(hy);

    CS_FREE_HD(z);

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {
      double dmaxg;
      MPI_Allreduce(&dmax, &dmaxg, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);
      dmax = dmaxg;
    }
#endif

    bft_printf("AmgX matrix check: max diff with ref: %g\n", dmax);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call AmgX linear equation solver.
 *
 * \warning The precision, r_norm, and n_iter parameters are ignored here.
 *          the matching configuration options should be set earlier, using
 *          the \ref cs_sles_amgx_set_config function
 *
 * \param[in, out]  context        pointer to AmgX solver info and context
 *                                 (actual type: cs_sles_amgx_t  *)
 * \param[in]       name           pointer to system name
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if NULL)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_amgx_solve(void                *context,
                   const char          *name,
                   const cs_matrix_t   *a,
                   int                  verbosity,
                   double               precision,
                   double               r_norm,
                   int                 *n_iter,
                   double              *residual,
                   const cs_real_t     *rhs,
                   cs_real_t           *vx,
                   size_t               aux_size,
                   void                *aux_vectors)
{
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  AMGX_RC retval = AMGX_RC_OK;

  char err_str[4096];
  const char error_fmt[] = N_("%s returned %d.\n"
                              "%s");

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_amgx_t  *c = context;
  cs_sles_amgx_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    cs_sles_amgx_setup(c, name, a, verbosity);
    sd = c->setup_data;
  }

  AMGX_vector_handle  x, b;

  sd->r_norm = r_norm;

  int       its = -1;
  double    _residual = -1;
  const int n_rows = cs_matrix_get_n_rows(a);
  const int db_size = cs_matrix_get_diag_block_size(a);

  /* Try to set tolerance to normalized value. */
  {
#if 0
    const char config_fmt[] =
      "{"
      "  \"config_version\": 2, "
      "  \"main\": {"
      "    \"convergence\": \"ABSOLUTE\", "
      "    \"tolerance\": %e}"
      "}";

#else
    const char config_fmt[] =
      "config_version=2, main:convergence=ABSOLUTE, main:tolerance=%e";
#endif

    char options[256];
    snprintf(options, 255, config_fmt, precision*r_norm);
    options[255] = '\0';
    retval = AMGX_config_add_parameters(&(c->amgx_config), options);
    if (retval != AMGX_RC_OK) {
      AMGX_get_error_string(retval, err_str, 4096);
      bft_error(__FILE__, __LINE__, 0, _(error_fmt),
                "AMGX_config_add_parameters", retval, err_str);
    }
  }

  /* Vector */

  AMGX_vector_create(&x, c->amgx_resources, c->amgx_mode);
  AMGX_vector_create(&b, c->amgx_resources, c->amgx_mode);

  if (cs_glob_n_ranks > 1) {
    AMGX_vector_bind(x, c->setup_data->matrix);
    AMGX_vector_bind(b, c->setup_data->matrix);
  }

  unsigned int n_bytes = n_rows*db_size*sizeof(cs_real_t);

  cs_alloc_mode_t amode_vx = cs_check_device_ptr(vx);
  cs_alloc_mode_t amode_rhs = cs_check_device_ptr(rhs);
  if (amode_vx < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)vx, n_bytes);
  if (amode_rhs < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_pin_memory((void *)rhs, n_bytes);

  retval = AMGX_vector_upload(x, n_rows, db_size, vx);
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_vector_upload", retval, err_str);
  }

  retval = AMGX_vector_upload(b, n_rows, db_size, rhs);
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_vector_upload", retval, err_str);
  }

  /* Resolution */

  cs_fp_exception_disable_trap();

  retval = AMGX_solver_solve(sd->solver, b, x);
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_solver_solve", retval, err_str);
  }

  cs_fp_exception_restore_trap();

  retval = AMGX_vector_download(x, vx);
  if (retval != AMGX_RC_OK) {
    AMGX_get_error_string(retval, err_str, 4096);
    bft_error(__FILE__, __LINE__, 0, _(error_fmt),
              "AMGX_vector_download", retval, err_str);
  }

  AMGX_vector_destroy(x);
  AMGX_vector_destroy(b);

  AMGX_solver_get_iterations_number(sd->solver, &its);
  AMGX_solver_get_iteration_residual(sd->solver, its-1, 0, &_residual);

  if (amode_vx < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)vx);
  if (amode_rhs < CS_ALLOC_HOST_DEVICE_PINNED)
    AMGX_unpin_memory((void *)rhs);

  AMGX_SOLVE_STATUS  solve_status;
  AMGX_solver_get_status(sd->solver, &solve_status);

  switch(solve_status) {
  case AMGX_SOLVE_SUCCESS:
    cvg = CS_SLES_CONVERGED;
    break;
  case AMGX_SOLVE_FAILED:
    cvg = CS_SLES_DIVERGED;
    break;
  case AMGX_SOLVE_DIVERGED:
    if (its >= c->n_iterations_max)
      cvg = CS_SLES_MAX_ITERATION;
    else
      cvg = CS_SLES_DIVERGED;
  }

  cs_fp_exception_restore_trap();

  *residual = _residual;
  *n_iter = its;

  /* Update return values */

  if (c->n_solves == 0)
    c->n_iterations_min = its;

  c->n_iterations_last = its;
  c->n_iterations_tot += its;
  if (c->n_iterations_min > its)
    c->n_iterations_min = its;
  if (c->n_iterations_max < its)
    c->n_iterations_max = its;
  c->n_solves += 1;
  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free AmgX linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to AmgX solver info and context
 *                           (actual type: cs_sles_amgx_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_free(void  *context)
{
  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_amgx_t  *c  = context;
  cs_sles_amgx_setup_t *sd = c->setup_data;

  if (sd != NULL) {

    AMGX_solver_destroy(sd->solver);
    AMGX_matrix_destroy(sd->matrix);

  }
  if (c->setup_data != NULL)
    BFT_FREE(c->setup_data);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to AmgX solver info and context
 *                       (actual type: cs_sles_amgx_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_log(const void  *context,
                 cs_log_t     log_type)
{
  const cs_sles_amgx_t  *c = context;

  const char m_type[] = "CSR";

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  _("  Solver type:                       AmgX\n"
                    "    Matrix format:                     %s\n"),
                  m_type);

  }
  else if (log_type == CS_LOG_PERFORMANCE) {

    int n_calls = c->n_solves;
    int n_it_min = c->n_iterations_min;
    int n_it_max = c->n_iterations_max;
    int n_it_mean = 0;

    if (n_calls > 0)
      n_it_mean = (int)(  c->n_iterations_tot
                        / ((unsigned long long)n_calls));

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   AmgX\n"
                    "    Matrix format:               %s\n"
                    "  Number of setups:              %12d\n"
                    "  Number of calls:               %12d\n"
                    "  Minimum number of iterations:  %12d\n"
                    "  Maximum number of iterations:  %12d\n"
                    "  Mean number of iterations:     %12d\n"
                    "  Total setup time:              %12.3f\n"
                    "  Total solution time:           %12.3f\n"),
                  m_type,
                  c->n_setups, n_calls, n_it_min, n_it_max, n_it_mean,
                  c->t_setup.nsec*1e-9,
                  c->t_solve.nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on AmgX library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_amgx_library_info(cs_log_t  log_type)
{
  int major = 0, minor = 0;
  AMGX_get_api_version(&major, &minor);

  cs_log_printf(log_type,
                "    AmgX %d.%d\n", major, minor);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
