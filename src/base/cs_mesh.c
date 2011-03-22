/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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
 * Main structure associated to a mesh
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_order.h>
#include <fvm_parall.h>
#include <fvm_interface.h>
#include <fvm_periodicity.h>
#include <fvm_selector.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh_halo.h"
#include "cs_numbering.h"
#include "cs_perio.h"
#include "cs_mesh_quantities.h"
#include "cs_ext_neighborhood.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_MESH_N_SUBS  5

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main mesh */

cs_mesh_t  *cs_glob_mesh = NULL;

/* Pointer on the temporary mesh used for building main mesh */

cs_mesh_builder_t  *cs_glob_mesh_builder = NULL;

/*============================================================================
 * Public FORTRAN function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check necessity of extended mesh from FORTRAN options.
 *
 * Interface Fortran :
 *
 * SUBROUTINE HALTYP (IVOSET)
 * *****************
 *
 * INTEGER          IVOSET      : <-- : Indicator of necessity of extended mesh
 *----------------------------------------------------------------------------*/

extern void
CS_PROCF (haltyp, HALTYP)(const cs_int_t   *ivoset);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update cell family array in case of parallelism.
 *
 * This function aims at copying main values from cells in halo (id between 1
 * and n_cells) to ghost cells on distant ranks (id between n_cells + 1 to
 * n_cells_with_halo).
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_sync_cell_fam(cs_mesh_t   *const mesh)
{
  cs_halo_t  *halo = mesh->halo;

  if (halo == NULL)
    return;

  bft_printf(_("Synchronizing cell families\n"));

  cs_halo_sync_num(halo, CS_HALO_EXTENDED, mesh->cell_family);
}

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *   min       --> minimum
 *   max       --> maximum
 *----------------------------------------------------------------------------*/

static void
_compute_local_minmax(cs_int_t        n_vals,
                      const cs_int_t  var[],
                      cs_int_t       *min,
                      cs_int_t       *max)
{
  cs_int_t  i;
  cs_int_t  _min = var[0], _max = var[0];

  for (i = 1; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

  if (min != NULL)  *min = _min;
  if (max != NULL)  *max = _max;
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a integer vector.
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *----------------------------------------------------------------------------*/

static void
_display_histograms(cs_int_t        n_vals,
                    const cs_int_t  var[])
{
  cs_int_t  i, j, k, val_max, val_min;
  double step;

  cs_int_t count[CS_MESH_N_SUBS];
  cs_int_t n_steps = CS_MESH_N_SUBS;

  /* Compute local min and max */

  if (n_vals == 0) {
    bft_printf(_("    no value\n"));
    return;
  }

  val_max = var[0];
  val_min = var[0];
  _compute_local_minmax(n_vals, var, &val_min, &val_max);

  bft_printf(_("    minimum value =         %10d\n"), (int)val_min);
  bft_printf(_("    maximum value =         %10d\n\n"), (int)val_max);

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  if (val_max - val_min > 0) {

    if (val_max-val_min < n_steps)
      n_steps = CS_MAX(1, floor(val_max-val_min));

    step = (double)(val_max - val_min) / n_steps;

    /* Loop on values */

    for (i = 0; i < n_vals; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < val_min + k*step)
          break;
      }
      count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3d : [ %10d ; %10d [ = %10d\n",
                 i+1,
                 (int)(val_min + i*step),
                 (int)(val_min + j*step),
                 (int)(count[i]));

    bft_printf("    %3d : [ %10d ; %10d ] = %10d\n",
               n_steps,
               (int)(val_min + (n_steps - 1)*step),
               (int)val_max,
               (int)(count[n_steps - 1]));

  }

  else { /* if (val_max == val_min) */
    bft_printf("    %3d : [ %10d ; %10d ] = %10d\n",
               1, (int)(val_min), (int)val_max, (int)n_vals);
  }

}

/*----------------------------------------------------------------------------
 * Write a summary about halo features in listing
 *
 * parameters:
 *   mesh                   <-- pointer to cs_mesh_t structure
 *   interface_time         <-- time elapsed in interface build
 *   halo_time              <-- time elapsed in halo build
 *   ext_neighborhood_time  <-- time elapsed in extended neighborhood
 *                              connectivity building
 *----------------------------------------------------------------------------*/

static void
_print_halo_info(const cs_mesh_t  *mesh,
                 double            interface_time,
                 double            halo_time,
                 double            ext_neighborhood_time)
{
  /* Summary of the computional times */

  bft_printf(_("\n Halo creation times summary\n\n"));

  if (mesh->n_domains > 1 || mesh->n_init_perio > 0)
    bft_printf(_("     Interface creation:                       %.3g s\n"),
               interface_time);

  if (mesh->halo_type == CS_HALO_EXTENDED)
    bft_printf(_("     Extended connectivity creation:            %.3g s\n"),
               ext_neighborhood_time);

  bft_printf(_("     Halo creation:                             %.3g s\n\n"),
             halo_time);


  bft_printf(_("     Total time for halo creation:              %.3g s\n\n"),
             halo_time + interface_time + ext_neighborhood_time);

  bft_printf(" ----------------------------------------------------------\n\n");

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Write a summary about cell neighbor features in listing
 *
 * parameters:
 *   mesh                   <-- pointer to cs_mesh_t structure
 *----------------------------------------------------------------------------*/

static void
_print_cell_neighbor_info(const cs_mesh_t  *mesh)
{
  cs_int_t i, j, k;
  float step;

  int n_steps = CS_MESH_N_SUBS;
  int n_min_neighbors = 0;
  int n_max_neighbors = 0;

  fvm_gnum_t  count[CS_MESH_N_SUBS];

  int  *n_cell_neighbors = NULL;

  /* Summary of the number of cell neighbors */

  BFT_MALLOC(n_cell_neighbors, mesh->n_cells_with_ghosts, int);

  for (i = 0; i < mesh->n_cells_with_ghosts; i++)
    n_cell_neighbors[i] = 0;

  for (i = 0; i < mesh->n_i_faces; i++) {
    cs_int_t c_id_0 = mesh->i_face_cells[i*2] - 1;
    cs_int_t c_id_1 = mesh->i_face_cells[i*2 + 1] - 1;
    n_cell_neighbors[c_id_0] += 1;
    n_cell_neighbors[c_id_1] += 1;
  }

  if (mesh->n_cells > 0)
    n_min_neighbors = n_cell_neighbors[0];

  for (i = 0; i < mesh->n_cells; i++) {
    if (n_cell_neighbors[i] < n_min_neighbors)
      n_min_neighbors = n_cell_neighbors[i];
    else if (n_cell_neighbors[i] > n_max_neighbors)
      n_max_neighbors = n_cell_neighbors[i];
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int n_g_min, n_g_max;

    MPI_Allreduce(&n_min_neighbors, &n_g_min, 1, MPI_INT, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_max_neighbors, &n_g_max, 1, MPI_INT, MPI_MAX,
                  cs_glob_mpi_comm);

    n_min_neighbors = n_g_min;
    n_max_neighbors = n_g_max;
  }

#endif /* defined(HAVE_MPI) */

  bft_printf(_("\n Histogram of the number of interior faces per cell:\n\n"));

  bft_printf(_("    minimum value =         %10d\n"), n_min_neighbors);
  bft_printf(_("    maximum value =         %10d\n\n"), n_max_neighbors);

  /* Define axis subdivisions */

  for (i = 0; i < CS_MESH_N_SUBS; i++)
    count[i] = 0;

  if (n_max_neighbors - n_min_neighbors > 0) {

    if (n_max_neighbors - n_min_neighbors < n_steps)
      n_steps = n_max_neighbors - n_min_neighbors;

    step = (float)(n_max_neighbors - n_min_neighbors) / n_steps;

    /* Loop on values */

    for (i = 0; i < mesh->n_cells; i++) {

      /* Associated subdivision */
      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (n_cell_neighbors[i] < n_min_neighbors + k*step)
          break;
      }
      count[j] += 1;
    }

#if defined(HAVE_MPI)

    if (cs_glob_n_ranks > 1) {

      fvm_gnum_t g_count[CS_MESH_N_SUBS];

      MPI_Allreduce(count, g_count, n_steps, FVM_MPI_GNUM, MPI_SUM,
                    cs_glob_mpi_comm);

      for (i = 0; i < n_steps; i++)
        count[i] = g_count[i];
    }

#endif

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3d : [ %10d ; %10d [ = %10llu\n",
                 i+1,
                 (int)(n_min_neighbors + i*step),
                 (int)(n_min_neighbors + j*step),
                 (unsigned long long)(count[i]));

    bft_printf("    %3d : [ %10d ; %10d ] = %10llu\n",
               n_steps,
               (int)(n_min_neighbors + (n_steps - 1)*step),
               n_max_neighbors,
               (unsigned long long)(count[n_steps - 1]));
  }

  else { /* if (n_min_neighbors == n_max_neighbors) */
    bft_printf("    %3d : [ %10d ; %10d ] = %10llu\n",
               1, n_min_neighbors, n_max_neighbors,
               (unsigned long long)mesh->n_g_cells);

  }

  bft_printf("\n ----------------------------------------------------------\n");

  /* Cleanup */

  BFT_FREE(n_cell_neighbors);
}

/*----------------------------------------------------------------------------
 * Compare periodic couples in global numbering form (qsort function).
 *
 * parameters:
 *   x <-> pointer to first couple
 *   y <-> pointer to second couple
 *
 * returns:
 *   lexicographical
 *----------------------------------------------------------------------------*/

static int _compare_couples(const void *x, const void *y)
{
  int retval = 1;

  const fvm_gnum_t *c0 = x;
  const fvm_gnum_t *c1 = y;

  if (c0[0] < c1[0])
    retval = -1;

  else if (c0[0] == c1[0]) {
    if (c0[1] < c1[1])
      retval = -1;
    else if (c0[1] == c1[1])
      retval = 0;
  }

  return retval;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Define parameters for building an vertices interface set structure on
 * a given mesh.
 *
 * This version is used for distributed parallel mode.
 *
 * parameters:
 *   mesh              <-- pointer to mesh structure
 *   face_ifs          <-- pointer to face interface structure describing
 *                         periodic face couples
 *   p_n_perio_couples --> pointer to the number of periodic couples
 *   p_perio_couples   --> pointer to the periodic couple list
 *----------------------------------------------------------------------------*/

static void
_define_perio_vtx_couples_g(const cs_mesh_t             *mesh,
                            const fvm_interface_set_t   *face_ifs,
                            fvm_lnum_t                  *p_n_perio_couples[],
                            fvm_gnum_t                 **p_perio_couples[])
{
  fvm_lnum_t  *n_perio_couples = NULL;
  fvm_gnum_t  **perio_couples = NULL;

  int i, j;
  fvm_lnum_t k, l;

  int perio_count = 0;
  fvm_lnum_t  *send_index = NULL, *recv_index = NULL;
  fvm_gnum_t  *send_num = NULL, *recv_num = NULL;
  int  *tr_id = NULL;

  const int n_perio = mesh->n_init_perio;
  const int n_interfaces = fvm_interface_set_size(face_ifs);
  const fvm_gnum_t *vtx_gnum = mesh->global_vtx_num;

  /* List direct and reverse transforms */

  BFT_MALLOC(n_perio_couples, n_perio, fvm_lnum_t);
  BFT_MALLOC(perio_couples, n_perio, fvm_gnum_t *);

  for (i = 0; i < n_perio; i++) {
    n_perio_couples[i] = 0;
    perio_couples[i] = NULL;
  }

  BFT_MALLOC(tr_id, n_perio*2, int);

  for (i = 0; i < n_perio*2; i++) {
    int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, i);
    if (i < rev_id) {
      int parent_ids[2];
      fvm_periodicity_get_parent_ids(mesh->periodicity, i, parent_ids);
      if (parent_ids[0] < 0 && parent_ids[1] < 0) {
        tr_id[perio_count*2] = i + 1;
        tr_id[perio_count*2 + 1] = rev_id + 1;
        perio_count++;
      }
    }
  }
  assert(perio_count == n_perio);

  /* Count maximum vertex couples and build exchange index;
     count only couples in direct periodicity direction, as the number
     of faces in positive and reverse direction should match */

  BFT_MALLOC(send_index, n_interfaces + 1, fvm_lnum_t);
  BFT_MALLOC(recv_index, n_interfaces + 1, fvm_lnum_t);
  send_index[0] = 0;
  recv_index[0] = 0;

  for (i = 0; i < n_interfaces; i++) {

    fvm_lnum_t send_size = 0;
    fvm_lnum_t recv_size = 0;

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);

    for (j = 0; j < n_perio; j++) {

      const fvm_lnum_t tr_start = tr_index[tr_id[j*2]];
      const fvm_lnum_t tr_end = tr_index[tr_id[j*2] + 1];
      const fvm_lnum_t rev_start = tr_index[tr_id[j*2+1]];
      const fvm_lnum_t rev_end = tr_index[tr_id[j*2+1] + 1];

      for (k = rev_start; k < rev_end; k++) {
        fvm_lnum_t l_id = loc_num[k] - 1;
        send_size += (  mesh->i_face_vtx_idx[l_id+1]
                      - mesh->i_face_vtx_idx[l_id]);
      }
      for (k = tr_start; k < tr_end; k++) {
        fvm_lnum_t l_id = loc_num[k] - 1;
        fvm_lnum_t n_vtx = (  mesh->i_face_vtx_idx[l_id+1]
                            - mesh->i_face_vtx_idx[l_id]);
        recv_size += n_vtx;
        n_perio_couples[j] += n_vtx;
      }
    }

    send_index[i+1] = send_index[i] + send_size;
    recv_index[i+1] = recv_index[i] + recv_size;
  }

  BFT_MALLOC(send_num, send_index[n_interfaces], fvm_gnum_t);
  BFT_MALLOC(recv_num, recv_index[n_interfaces], fvm_gnum_t);

  /* Prepare send buffer (send reverse transformation values) */

  for (i = 0, j = 0; i < n_interfaces; i++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);

    for (j = 0, l = send_index[i]; j < n_perio; j++) {

      const fvm_lnum_t start_id = tr_index[tr_id[j*2+1]];
      const fvm_lnum_t end_id = tr_index[tr_id[j*2+1] + 1];

      for (k = start_id; k < end_id; k++) {
        fvm_lnum_t v_id;
        fvm_lnum_t f_id = loc_num[k] - 1;
        for (v_id = mesh->i_face_vtx_idx[f_id];
             v_id < mesh->i_face_vtx_idx[f_id + 1];
             v_id++)
          send_num[l++] = vtx_gnum[mesh->i_face_vtx_lst[v_id - 1] - 1];
      }
    }
  }

  /* Exchange global vertex numbers */

  {
    MPI_Request  *request = NULL;
    MPI_Status  *status  = NULL;

    int request_count = 0;

    BFT_MALLOC(request, n_interfaces*2, MPI_Request);
    BFT_MALLOC(status, n_interfaces*2, MPI_Status);

    for (i = 0; i < n_interfaces; i++) {
      const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
      int distant_rank = fvm_interface_rank(face_if);
      MPI_Irecv(recv_num + recv_index[i],
                recv_index[i+1] - recv_index[i],
                FVM_MPI_GNUM,
                distant_rank,
                distant_rank,
                cs_glob_mpi_comm,
                &(request[request_count++]));
    }

    for (i = 0; i < n_interfaces; i++) {
      const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
      int distant_rank = fvm_interface_rank(face_if);
      MPI_Isend(send_num + send_index[i],
                send_index[i+1] - send_index[i],
                FVM_MPI_GNUM,
                distant_rank,
                (int)cs_glob_rank_id,
                cs_glob_mpi_comm,
                &(request[request_count++]));
    }

    MPI_Waitall(request_count, request, status);

    BFT_FREE(request);
    BFT_FREE(status);
  }

  BFT_FREE(send_num);
  BFT_FREE(send_index);

  for (i = 0; i < n_perio; i++)
    BFT_MALLOC(perio_couples[i], n_perio_couples[i]*2, fvm_gnum_t);

  /* Reset couples count */

  for (i = 0; i < n_perio; i++)
    n_perio_couples[i] = 0;

  /* Now fill couples */

  for (i = 0, j = 0; i < n_interfaces; i++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);

    for (j = 0, l = recv_index[i]; j < n_perio; j++) {

      fvm_lnum_t nc = n_perio_couples[j]*2;
      const fvm_lnum_t start_id = tr_index[tr_id[j*2]];
      const fvm_lnum_t end_id = tr_index[tr_id[j*2] + 1];

      for (k = start_id; k < end_id; k++) {
        fvm_lnum_t v_id;
        fvm_lnum_t f_id = loc_num[k] - 1;
        for (v_id = mesh->i_face_vtx_idx[f_id];
             v_id < mesh->i_face_vtx_idx[f_id + 1];
             v_id++) {
          perio_couples[j][nc++] = vtx_gnum[mesh->i_face_vtx_lst[v_id - 1] - 1];
          perio_couples[j][nc++] = recv_num[l++];
        }
      }
      n_perio_couples[j] = nc/2;
    }
  }

  BFT_FREE(recv_num);
  BFT_FREE(recv_index);
  BFT_FREE(tr_id);

  /* Set return values */

  *p_n_perio_couples = n_perio_couples;
  *p_perio_couples = perio_couples;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Define parameters for building an vertices interface set structure on
 * a given mesh.
 *
 * This version is used for single-process mode.
 *
 * parameters:
 *   mesh              <-- pointer to mesh structure
 *   face_ifs          <-- pointer to face interface structure describing
 *                         periodic face couples
 *   p_n_perio_couples --> pointer to the number of periodic couples
 *   p_perio_couples   --> pointer to the periodic couple list
 *----------------------------------------------------------------------------*/

static void
_define_perio_vtx_couples_l(const cs_mesh_t             *mesh,
                            const fvm_interface_set_t   *face_ifs,
                            fvm_lnum_t                  *p_n_perio_couples[],
                            fvm_gnum_t                 **p_perio_couples[])
{
  fvm_lnum_t  *n_perio_couples = NULL;
  fvm_gnum_t  **perio_couples = NULL;

  int i, j;
  fvm_lnum_t k;

  int perio_count = 0;
  int  *tr_id = NULL;

  const int n_perio = mesh->n_init_perio;
  const int n_interfaces = fvm_interface_set_size(face_ifs);
  const fvm_gnum_t *vtx_gnum = mesh->global_vtx_num;

  /* List direct and reverse transforms */

  BFT_MALLOC(n_perio_couples, n_perio, fvm_lnum_t);
  BFT_MALLOC(perio_couples, n_perio, fvm_gnum_t *);

  for (i = 0; i < n_perio; i++) {
    n_perio_couples[i] = 0;
    perio_couples[i] = NULL;
  }

  BFT_MALLOC(tr_id, n_perio*2, int);

  for (i = 0; i < n_perio*2; i++) {
    int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, i);
    if (i < rev_id) {
      int parent_ids[2];
      fvm_periodicity_get_parent_ids(mesh->periodicity, i, parent_ids);
      if (parent_ids[0] < 0 && parent_ids[1] < 0) {
        tr_id[perio_count*2] = i + 1;
        tr_id[perio_count*2 + 1] = rev_id + 1;
        perio_count++;
      }
    }
  }
  assert(perio_count == n_perio);

  /* Count maximum vertex couples and build exchange index;
     count only couples in direct periodicity direction, as the number
     of faces in positive and reverse direction should match */

  for (i = 0; i < n_interfaces; i++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);

    for (j = 0; j < n_perio; j++) {

      const fvm_lnum_t tr_start = tr_index[tr_id[j*2]];
      const fvm_lnum_t tr_end = tr_index[tr_id[j*2] + 1];

      for (k = tr_start; k < tr_end; k++) {
        fvm_lnum_t l_id = loc_num[k] - 1;
        n_perio_couples[j] += (  mesh->i_face_vtx_idx[l_id+1]
                               - mesh->i_face_vtx_idx[l_id]);
      }
    }
  }

  for (i = 0; i < n_perio; i++)
    BFT_MALLOC(perio_couples[i], n_perio_couples[i]*2, fvm_gnum_t);

  /* Reset couples count */

  for (i = 0; i < n_perio; i++)
    n_perio_couples[i] = 0;

  /* Now fill couples */

  for (i = 0, j = 0; i < n_interfaces; i++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, i);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);
    const fvm_lnum_t *dist_num = fvm_interface_get_distant_num(face_if);

    for (j = 0; j < n_perio; j++) {

      fvm_lnum_t nc = n_perio_couples[j]*2;

      const fvm_lnum_t tr_start = tr_index[tr_id[j*2]];
      const fvm_lnum_t tr_end = tr_index[tr_id[j*2] + 1];

      for (k = tr_start; k < tr_end; k++) {
        fvm_lnum_t v_id;
        fvm_lnum_t f1_id = loc_num[k] - 1;
        fvm_lnum_t f2_id = dist_num[k] - 1;
        fvm_lnum_t n_vtx = (  mesh->i_face_vtx_idx[f1_id+1]
                            - mesh->i_face_vtx_idx[f1_id]);
        fvm_lnum_t v1_start = mesh->i_face_vtx_idx[f1_id] - 1;
        fvm_lnum_t v2_start = mesh->i_face_vtx_idx[f2_id] - 1;

        if (vtx_gnum != NULL) {
          for (v_id = 0; v_id < n_vtx; v_id++) {
            perio_couples[j][nc++]
              = vtx_gnum[mesh->i_face_vtx_lst[v1_start + v_id] - 1];
            perio_couples[j][nc++]
              = vtx_gnum[mesh->i_face_vtx_lst[v2_start + v_id] - 1];
          }
        }
        else {
          for (v_id = 0; v_id < n_vtx; v_id++) {
            perio_couples[j][nc++] = mesh->i_face_vtx_lst[v1_start + v_id];
            perio_couples[j][nc++] = mesh->i_face_vtx_lst[v2_start + v_id];
          }
        }
      }
      n_perio_couples[j] = nc/2;
    }
  }

  BFT_FREE(tr_id);

  /* Set return values */

  *p_n_perio_couples = n_perio_couples;
  *p_perio_couples = perio_couples;
}

/*----------------------------------------------------------------------------
 * Define parameters for building an vertices interface set structure on
 * a given mesh.
 *
 * parameters:
 *   mesh              <-- pointer to mesh structure
 *   face_ifs          <-- pointer to face interface structure describing
 *                         periodic face couples
 *   p_n_perio_couples --> pointer to the number of periodic couples
 *   p_perio_couples   --> pointer to the periodic couple list
 *----------------------------------------------------------------------------*/

static void
_define_perio_vtx_couples(const cs_mesh_t             *mesh,
                          const fvm_interface_set_t   *face_ifs,
                          fvm_lnum_t                  *p_n_perio_couples[],
                          fvm_gnum_t                 **p_perio_couples[])
{
  fvm_lnum_t  *n_perio_couples = NULL;
  fvm_gnum_t  **perio_couples = NULL;

  int i;
  fvm_lnum_t k, l;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _define_perio_vtx_couples_g(mesh,
                                face_ifs,
                                p_n_perio_couples,
                                p_perio_couples);

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1)
    _define_perio_vtx_couples_l(mesh,
                                face_ifs,
                                p_n_perio_couples,
                                p_perio_couples);

  /* Finally, remove duplicate couples */

  n_perio_couples = *p_n_perio_couples;
  perio_couples = *p_perio_couples;

  for (i = 0; i < mesh->n_init_perio; i++) {

    if (n_perio_couples[i] > 1) {

      fvm_lnum_t nc = n_perio_couples[i];
      fvm_gnum_t *pc = perio_couples[i];

      qsort(pc, nc, sizeof(fvm_gnum_t) * 2, &_compare_couples);

      for (k = 1, l = 0; k < nc; k++) {
        if ((pc[k*2] != pc[l*2]) && (pc[k*2+1] != pc[l*2+1])) {
          l += 1;
          pc[l*2] = pc[k*2];
          pc[l*2+1] = pc[k*2+1];
        }
      }
      n_perio_couples[i] = l+1;

      BFT_REALLOC(perio_couples[i], n_perio_couples[i]*2, fvm_gnum_t);
    }
  }

  *p_n_perio_couples = n_perio_couples;
  *p_perio_couples = perio_couples;
}

/*----------------------------------------------------------------------------
 * Assign a periodicity number and optionnaly a rank id to halo elements.
 *
 * parameters:
 *   mesh           <-- pointer to mesh structure
 *   halo_perio_num --> periodicity number associated with each halo cell,
 *                      signed for direct/reverse transform.
 *                      (size: mesh->n_ghost_cells)
 *   halo_rank_id   --> halo local rank id associated with each halo cell,
 *                      or NULL (size: mesh->n_ghost_cells)
 *
 * returns:
 *   the number of base periodicities defined
 *----------------------------------------------------------------------------*/

static int
_get_halo_perio_num(const cs_mesh_t  *mesh,
                    int               halo_perio_num[],
                    int               halo_rank_id[])
{
  int         i, rank_id;
  fvm_lnum_t  j;

  int n_perio = 0;

  int *tr_id = NULL;

  const cs_halo_t  *halo = mesh->halo;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  /* Initialization */

  assert(halo != NULL);
  assert(periodicity != NULL);

  /* List direct and non-combined transforms */

  BFT_MALLOC(tr_id, halo->n_transforms, int);
  for (i = 0; i < halo->n_transforms; i++)
    tr_id[i] = -1;

  n_perio = 0;

  for (i = 0; i < halo->n_transforms; i++) {
    int rev_id = fvm_periodicity_get_reverse_id(periodicity, i);
    if (i < rev_id) {
      int parent_ids[2];
      fvm_periodicity_get_parent_ids(periodicity, i, parent_ids);
      if (parent_ids[0] < 0 && parent_ids[1] < 0) {
        tr_id[n_perio*2] = i;
        tr_id[n_perio*2+1] = rev_id;
        n_perio++;
      }
    }
  }
  assert(n_perio == mesh->n_init_perio);

  BFT_REALLOC(tr_id, n_perio*2, int);
  assert(n_perio == mesh->n_init_perio);

  /* Initialize halo values */

  for (j = 0; j < mesh->n_ghost_cells; j++)
    halo_perio_num[j] = 0;

  if (halo_rank_id != NULL) {
    for (j = 0; j < mesh->n_ghost_cells; j++)
      halo_rank_id[j] = -1;
  }

  /* For each periodic face couple, we have 2 global ghost cells: one for
     the direct periodicity, one for the reverse. To match couples only
     once, we ignore the reverse periodicity, leaving exactly 1 ghost cell
     per couple (which may belong to multiple periodic faces). */

  /* Mark ghost cells with their periodicity number */

  for (i = 0; i < n_perio; i++) {

    fvm_lnum_t shift, start, end;

    /* Direct periodicity */

    shift = 4 * halo->n_c_domains * tr_id[i*2];

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start = halo->perio_lst[shift + 4*rank_id];
      end = start + halo->perio_lst[shift + 4*rank_id + 1];

      for (j = start; j < end; j++)
        halo_perio_num[j] = i + 1;

      if (halo_rank_id != NULL) {
        for (j = start; j < end; j++)
          halo_rank_id[j] = rank_id;
      }
    }

    /* Reverse periodicity */

    shift = 4 * halo->n_c_domains * tr_id[i*2 + 1];

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start = halo->perio_lst[shift + 4*rank_id];
      end = start + halo->perio_lst[shift + 4*rank_id + 1];

      for (j = start; j < end; j++)
        halo_perio_num[j] = -(i + 1);

      if (halo_rank_id != NULL) {
        for (j = start; j < end; j++)
          halo_rank_id[j] = rank_id;
      }
    }
  }

  BFT_FREE(tr_id);

  return n_perio;
}

/*----------------------------------------------------------------------------
 * Build a cell -> face search structure for local cells with a periodic
 * face in the direct peridicity direction.
 *
 * The caller is responsible for freeing the arrays allocated and returned
 * by this function once they are no longer needed.
 *
 * parameters:
 *   mesh           <-- pointer to mesh structure
 *   halo_perio_num <-- periodicity number of ghost cells
 *   cell_face_idx  --> index (0 to n-1) of face ids for a given cell
 *   cell_face      --> face ids (0 to n-1) for a given cell
 *----------------------------------------------------------------------------*/

static void
_build_cell_face_perio_match(const cs_mesh_t    *mesh,
                             const int          *halo_perio_num,
                             fvm_lnum_t        **cell_face_idx,
                             fvm_lnum_t        **cell_face)
{
  fvm_lnum_t  i;

  fvm_lnum_t *_cell_face_count, *_cell_face_idx, *_cell_face;

  /* Build cell -> face indexed search array containing only faces
     adjacent to the halo. */

  BFT_MALLOC(_cell_face_count, mesh->n_cells_with_ghosts, fvm_lnum_t);
  BFT_MALLOC(_cell_face_idx, mesh->n_cells_with_ghosts + 1, fvm_lnum_t);

  for (i = 0; i < mesh->n_cells_with_ghosts; i++)
    _cell_face_count[i] = 0;

  for (i = 0; i < mesh->n_i_faces; i++) {
    const fvm_lnum_t c_id0 = mesh->i_face_cells[i*2] - 1;
    const fvm_lnum_t c_id1 = mesh->i_face_cells[i*2+1] - 1;
    if (c_id0 < mesh->n_cells && c_id1 >= mesh->n_cells) {
      if (halo_perio_num[c_id1 - mesh->n_cells] < 0)
        _cell_face_count[c_id0] += 1;
    }
    else if (c_id1 < mesh->n_cells && c_id0 >= mesh->n_cells) {
      if (halo_perio_num[c_id0 - mesh->n_cells] < 0)
      _cell_face_count[c_id1] += 1;
    }
  }

  _cell_face_idx[0] = 0;
  for (i = 0; i < mesh->n_cells_with_ghosts; i++) {
    _cell_face_idx[i+1] = _cell_face_idx[i] + _cell_face_count[i];
    _cell_face_count[i] = 0;
  }

  BFT_MALLOC(_cell_face, _cell_face_idx[mesh->n_cells_with_ghosts], fvm_lnum_t);

  for (i = 0; i < mesh->n_i_faces; i++) {
    const fvm_lnum_t c_id0 = mesh->i_face_cells[i*2] - 1;
    const fvm_lnum_t c_id1 = mesh->i_face_cells[i*2+1] - 1;
    if (c_id0 < mesh->n_cells && c_id1 >= mesh->n_cells) {
      if (halo_perio_num[c_id1 - mesh->n_cells] < 0) {
        _cell_face[_cell_face_idx[c_id0] + _cell_face_count[c_id0]] = i;
        _cell_face_count[c_id0] += 1;
      }
    }
    else if (c_id1 < mesh->n_cells && c_id0 >= mesh->n_cells) {
      if (halo_perio_num[c_id0 - mesh->n_cells] < 0) {
        _cell_face[_cell_face_idx[c_id1] + _cell_face_count[c_id1]] = i;
        _cell_face_count[c_id1] += 1;
      }
    }
  }

  BFT_FREE(_cell_face_count);

  *cell_face_idx = _cell_face_idx;
  *cell_face = _cell_face;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get global lists of periodic faces.
 *
 * The caller is responsible for freeing the arrays allocated and returned
 * by this function once they are no longer needed.
 *
 * parameters:
 *   mesh                 <-- pointer to mesh structure
 *   n_perio_face_couples --> local number of periodic couples per
 *                          periodicity (size: mesh->n_init_perio)
 *   perio_face_couples   --> arrays of global periodic couple face numbers,
 *                            for each periodicity
 *----------------------------------------------------------------------------*/

static void
_get_perio_faces_g(const cs_mesh_t    *mesh,
                   fvm_lnum_t        **n_perio_face_couples,
                   fvm_gnum_t       ***perio_face_couples)
{
  int         i, rank_id, perio_num;
  fvm_lnum_t  j, k, l;
  fvm_lnum_t  dest_c_id, src_c_id;

  int         n_perio = 0;
  fvm_lnum_t  recv_size = 0;
  fvm_lnum_t *loc_cell_num = NULL;
  fvm_lnum_t *recv_count = NULL, *recv_idx = NULL;
  fvm_lnum_t *send_count = NULL, *send_idx = NULL;
  fvm_lnum_t *cell_face_idx, *cell_face;
  fvm_gnum_t *recv_vals = NULL, *send_vals = NULL;

  int *halo_perio_num = NULL;
  int *halo_rank_id = NULL;

  fvm_lnum_t    *_n_perio_face_couples = NULL;
  fvm_gnum_t    **_perio_face_couples = NULL;

  int           request_count = 0;
  MPI_Request  *halo_request = NULL;
  MPI_Status   *halo_status = NULL;

  const cs_halo_t  *halo = mesh->halo;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  /* Initialization */

  BFT_MALLOC(_n_perio_face_couples, mesh->n_init_perio, fvm_lnum_t);
  BFT_MALLOC(_perio_face_couples, mesh->n_init_perio, fvm_gnum_t *);

  for (i = 0; i < mesh->n_init_perio; i++) {
    _n_perio_face_couples[i] = 0;
    _perio_face_couples[i] = NULL;
  }

  assert(halo != NULL);
  assert(periodicity != NULL);

  /* Mark ghost cells with their periodicity number and halo rank id */

  BFT_MALLOC(halo_perio_num, mesh->n_ghost_cells, int);
  BFT_MALLOC(halo_rank_id, mesh->n_ghost_cells, int);

  n_perio = _get_halo_perio_num(mesh, halo_perio_num, halo_rank_id);

  /* Exchange local cell numbers, so that cell-face reverse matching
     is made easier later */

  BFT_MALLOC(loc_cell_num, mesh->n_cells_with_ghosts, fvm_lnum_t);
  for (j = 0; j < mesh->n_cells; j++)
    loc_cell_num[j] = j+1;
  for (j = mesh->n_cells; j < mesh->n_cells_with_ghosts; j++)
    loc_cell_num[j] = 0;

  cs_halo_sync_num(halo, CS_HALO_STANDARD, loc_cell_num);

  /* Now compute send and receive counts and build indexes */

  BFT_MALLOC(send_count, halo->n_c_domains, fvm_lnum_t);
  BFT_MALLOC(recv_count, halo->n_c_domains, fvm_lnum_t);
  BFT_MALLOC(send_idx, halo->n_c_domains + 1, fvm_lnum_t);
  BFT_MALLOC(recv_idx, halo->n_c_domains + 1, fvm_lnum_t);

  for (i = 0; i < halo->n_c_domains; i++) {
    send_count[i] = 0;
    recv_count[i] = 0;
  }

  for (j = 0; j < mesh->n_i_faces; j++) {
    const fvm_lnum_t h_id0 = mesh->i_face_cells[j*2] - mesh->n_cells - 1;
    const fvm_lnum_t h_id1 = mesh->i_face_cells[j*2+1] - mesh->n_cells - 1;
    if (h_id0 >= 0) {
      perio_num = halo_perio_num[h_id0];
      if (perio_num > 0)
        send_count[halo_rank_id[h_id0]] += 1;
      else if (perio_num < 0)
        recv_count[halo_rank_id[h_id0]] += 1;
    }
    else if (h_id1 >= 0) {
      perio_num = halo_perio_num[h_id1];
      if (perio_num > 0)
        send_count[halo_rank_id[h_id1]] += 1;
      else if (perio_num < 0)
        recv_count[halo_rank_id[h_id1] ]+= 1;
    }
  }

  send_idx[0] = 0;
  recv_idx[0] = 0;
  for (i = 0; i < halo->n_c_domains; i++) {
    send_idx[i+1] = send_idx[i] + send_count[i];
    recv_idx[i+1] = recv_idx[i] + recv_count[i];
    send_count[i] = 0;
  }

  /* Assemble the tuples to send; a tuple contains:
     destination cell id (cell matched by direct periodicity),
     source cell (cell matched by reverse periodicity),
     periodicity number,
     global face number (2nd face in a couple). */

  BFT_MALLOC(send_vals, send_idx[halo->n_c_domains]*4, fvm_gnum_t);

  for (j = 0; j < mesh->n_i_faces; j++) {

    const fvm_lnum_t c_id0 = mesh->i_face_cells[j*2] - 1;
    const fvm_lnum_t c_id1 = mesh->i_face_cells[j*2+1] - 1;

    dest_c_id = -1;
    if (c_id0 >= mesh->n_cells) {
      perio_num = halo_perio_num[c_id0  - mesh->n_cells];
      if (perio_num > 0) {
        dest_c_id = c_id1;
        src_c_id = c_id0;
        rank_id = halo_rank_id[c_id0  - mesh->n_cells];
      }
    }
    else if (c_id1 >= mesh->n_cells) {
      perio_num = halo_perio_num[c_id1 - mesh->n_cells];
      if (perio_num > 0) {
        dest_c_id = c_id0;
        src_c_id = c_id1;
        rank_id = halo_rank_id[c_id1  - mesh->n_cells];
      }
    }
    if (dest_c_id != -1) {
      k = send_idx[rank_id] + send_count[rank_id];
      send_vals[k*4] = loc_cell_num[dest_c_id];
      send_vals[k*4+1] = loc_cell_num[src_c_id];
      send_vals[k*4+2] = perio_num;
      send_vals[k*4+3] = mesh->global_i_face_num[j];
      send_count[rank_id] += 1;
    }
  }

  BFT_FREE(halo_rank_id);

  /* Now exchange face matching data */

  BFT_MALLOC(halo_request, halo->n_c_domains*2, MPI_Request);
  BFT_MALLOC(halo_status, halo->n_c_domains*2, MPI_Status);

  BFT_MALLOC(recv_vals, recv_idx[halo->n_c_domains]*4, fvm_gnum_t);

  /* Receive data from distant ranks */

  for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
    if (recv_count[rank_id] > 0)
      MPI_Irecv(recv_vals + (recv_idx[rank_id]*4),
                recv_count[rank_id]*4,
                FVM_MPI_GNUM,
                halo->c_domain_rank[rank_id],
                halo->c_domain_rank[rank_id],
                cs_glob_mpi_comm,
                &(halo_request[request_count++]));
  }

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
    if (send_count[rank_id] > 0)
      MPI_Isend(send_vals + (send_idx[rank_id]*4),
                send_count[rank_id]*4,
                FVM_MPI_GNUM,
                halo->c_domain_rank[rank_id],
                cs_glob_rank_id,
                cs_glob_mpi_comm,
                &(halo_request[request_count++]));
  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, halo_request, halo_status);

  BFT_FREE(halo_request);
  BFT_FREE(halo_status);

  BFT_FREE(send_vals);
  BFT_FREE(send_count);
  BFT_FREE(send_idx);

  recv_size = recv_idx[halo->n_c_domains];

  BFT_FREE(recv_count);
  BFT_FREE(recv_idx);

  /* Count periodic face couples */

  for (j = 0; j < recv_size; j++) {
    perio_num = recv_vals[j*4+2];
    assert(perio_num > 0);
    _n_perio_face_couples[perio_num - 1] += 1;
  }

  /* Allocate periodic couple lists */

  for (i = 0; i < mesh->n_init_perio; i++) {
    BFT_MALLOC(_perio_face_couples[i],
               _n_perio_face_couples[i]*2,
               fvm_gnum_t);
    _n_perio_face_couples[i] = 0;
  }

  /* Build cell -> face indexed search array containing only faces
     adjacent to the halo. */

  _build_cell_face_perio_match(mesh,
                               halo_perio_num,
                               &cell_face_idx,
                               &cell_face);

  BFT_FREE(halo_perio_num);

  /* Now match faces; faces already matched are marked so as not to be
     re-matched, in case faces have been split and multiple faces thus share
     the same cells (in which case we assume that subfaces are locally ordered
     identically). */

  for (j = 0; j <  recv_size; j++) {

    fvm_lnum_t dest_c_num = recv_vals[j*4];
    fvm_lnum_t src_c_num = recv_vals[j*4+1];
    int perio_id = recv_vals[j*4+2] - 1;
    fvm_gnum_t perio_face_gnum = recv_vals[j*4+3];

    for (k = cell_face_idx[src_c_num-1]; k < cell_face_idx[src_c_num]; k++) {

      l = cell_face[k]; /* Face _id */

      if (l > -1) {
        fvm_lnum_t loc_c_num, dist_c_id;
        if (mesh->i_face_cells[l*2] < mesh->i_face_cells[l*2+1]) {
          loc_c_num = mesh->i_face_cells[l*2];
          dist_c_id = mesh->i_face_cells[l*2+1] - 1;
        }
        else {
          loc_c_num = mesh->i_face_cells[l*2+1];
          dist_c_id = mesh->i_face_cells[l*2] - 1;
        }
        assert(loc_c_num <= mesh->n_cells && dist_c_id >= mesh->n_cells);

        /* If we have a match, update periodic couples */

        if (src_c_num == loc_c_num && dest_c_num == loc_cell_num[dist_c_id]) {
          fvm_lnum_t couple_id = _n_perio_face_couples[perio_id];
          cell_face[k] = -1; /* Mark as used */
          _perio_face_couples[perio_id][couple_id*2]
            = mesh->global_i_face_num[l];
          _perio_face_couples[perio_id][couple_id*2 + 1] = perio_face_gnum;
          _n_perio_face_couples[perio_id] += 1;
        }
      }
    }
  }

  BFT_FREE(cell_face_idx);
  BFT_FREE(cell_face);
  BFT_FREE(recv_vals);
  BFT_FREE(loc_cell_num);

  /* Set return values */

  *n_perio_face_couples = _n_perio_face_couples;
  *perio_face_couples = _perio_face_couples;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Get global lists of periodic faces.
 *
 * This version is used for serial mode.
 *
 * The caller is responsible for freeing the arrays allocated and returned
 * by this function once they are no longer needed.
 *
 * parameters:
 *   mesh                 <-- pointer to mesh structure
 *   n_perio_face_couples --> global number of periodic couples per
 *                            periodicity (size: mesh->n_init_perio)
 *   perio_face_couples   --> arrays of global periodic couple face numbers,
 *                            for each periodicity
 *----------------------------------------------------------------------------*/

static void
_get_perio_faces_l(const cs_mesh_t    *mesh,
                   fvm_lnum_t        **n_perio_face_couples,
                   fvm_gnum_t       ***perio_face_couples)
{
  int         i, perio_num;
  fvm_lnum_t  j, k, l;
  fvm_lnum_t  dest_c_id, src_c_id;

  int         n_perio = 0;
  fvm_lnum_t *loc_cell_num = NULL;
  fvm_lnum_t *cell_face_idx, *cell_face;

  int *halo_perio_num = NULL;

  fvm_lnum_t   *_n_perio_face_couples = NULL;
  fvm_gnum_t  **_perio_face_couples = NULL;

  const cs_halo_t  *halo = mesh->halo;

  /* Initialization */

  BFT_MALLOC(_n_perio_face_couples, mesh->n_init_perio, fvm_lnum_t);
  BFT_MALLOC(_perio_face_couples, mesh->n_init_perio, fvm_gnum_t *);

  for (i = 0; i < mesh->n_init_perio; i++) {
    _n_perio_face_couples[i] = 0;
    _perio_face_couples[i] = NULL;
  }

  assert(halo != NULL);
  assert(halo->n_c_domains == 1);

  /* Mark ghost cells with their periodicity number */

  BFT_MALLOC(halo_perio_num, mesh->n_ghost_cells, int);

  n_perio = _get_halo_perio_num(mesh, halo_perio_num, NULL);

  /* Exchange local cell numbers, so that cell-face reverse matching
     is made easier later */

  BFT_MALLOC(loc_cell_num, mesh->n_cells_with_ghosts, fvm_lnum_t);
  for (j = 0; j < mesh->n_cells; j++)
    loc_cell_num[j] = j+1;
  for (j = mesh->n_cells; j < mesh->n_cells_with_ghosts; j++)
    loc_cell_num[j] = 0;

  cs_halo_sync_num(halo, CS_HALO_STANDARD, loc_cell_num);

  /* Now compute number of periodic couples */

  for (j = 0; j < mesh->n_i_faces; j++) {
    const fvm_lnum_t h_id0 = mesh->i_face_cells[j*2] - mesh->n_cells - 1;
    const fvm_lnum_t h_id1 = mesh->i_face_cells[j*2+1] - mesh->n_cells - 1;
    if (h_id0 >= 0) {
      if (halo_perio_num[h_id0] > 0)
        _n_perio_face_couples[halo_perio_num[h_id0] - 1] += 1;
    }
    else if (h_id1 >= 0) {
      if (halo_perio_num[h_id1] > 0)
        _n_perio_face_couples[halo_perio_num[h_id1] - 1] += 1;
    }
  }

  /* Allocate periodic couple lists */

  for (i = 0; i < mesh->n_init_perio; i++) {
    BFT_MALLOC(_perio_face_couples[i],
               _n_perio_face_couples[i]*2,
               fvm_gnum_t);
    _n_perio_face_couples[i] = 0;
  }

  /* Build cell -> face indexed search array containing only faces
     adjacent to the halo. */

  _build_cell_face_perio_match(mesh,
                               halo_perio_num,
                               &cell_face_idx,
                               &cell_face);

  /* Now match faces; faces already matched are marked so as not to be
     re-matched, in case faces have been split and multiple faces thus share
     the same cells (in which case we assume that subfaces are locally ordered
     identically). */

  for (j = 0; j < mesh->n_i_faces; j++) {

    const fvm_lnum_t c_id0 = mesh->i_face_cells[j*2] - 1;
    const fvm_lnum_t c_id1 = mesh->i_face_cells[j*2+1] - 1;

    dest_c_id = -1;
    if (c_id0 >= mesh->n_cells) {
      perio_num = halo_perio_num[c_id0  - mesh->n_cells];
      if (perio_num > 0) {
        dest_c_id = c_id1;
        src_c_id = c_id0;
      }
    }
    else if (c_id1 >= mesh->n_cells) {
      perio_num = halo_perio_num[c_id1 - mesh->n_cells];
      if (perio_num > 0) {
        dest_c_id = c_id0;
        src_c_id = c_id1;
      }
    }
    if (dest_c_id != -1) {

      fvm_lnum_t dest_c_num = loc_cell_num[dest_c_id];
      fvm_lnum_t src_c_num = loc_cell_num[src_c_id];

      for (k = cell_face_idx[src_c_num-1]; k < cell_face_idx[src_c_num]; k++) {

        l = cell_face[k]; /* Face _id */

        if (l > -1) {
          fvm_lnum_t loc_c_num, dist_c_id;
          if (mesh->i_face_cells[l*2] < mesh->i_face_cells[l*2+1]) {
            loc_c_num = mesh->i_face_cells[l*2];
            dist_c_id = mesh->i_face_cells[l*2+1] - 1;
          }
          else {
            loc_c_num = mesh->i_face_cells[l*2+1];
            dist_c_id = mesh->i_face_cells[l*2] - 1;
          }
          assert(loc_c_num <= mesh->n_cells && dist_c_id >= mesh->n_cells);

          /* If we have a match, update periodic couples */

          if (src_c_num == loc_c_num && dest_c_num == loc_cell_num[dist_c_id]) {
            int perio_id = perio_num - 1;
            fvm_lnum_t couple_id = _n_perio_face_couples[perio_id];
            cell_face[k] = -1; /* Mark as used */
            _perio_face_couples[perio_id][couple_id*2] = l+1;
            _perio_face_couples[perio_id][couple_id*2 + 1] = j+1;
            _n_perio_face_couples[perio_id] += 1;
          }
        }
      }
    }
  }

  if (mesh->global_i_face_num != NULL) {
    for (i = 0; i < mesh->n_init_perio; i++) {
      for (j = 0; j < _n_perio_face_couples[i]*2; j++)
        _perio_face_couples[i][j]
          = mesh->global_i_face_num[_perio_face_couples[i][j] - 1];
    }
  }

  BFT_FREE(halo_perio_num);

  BFT_FREE(cell_face_idx);
  BFT_FREE(cell_face);
  BFT_FREE(loc_cell_num);

  /* Set return values */

  *n_perio_face_couples = _n_perio_face_couples;
  *perio_face_couples = _perio_face_couples;
}

/*----------------------------------------------------------------------------
 * Prepare local processor count for mesh stats.
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   n_group_elts  --> number of cells, interior faces, and boundary faces
 *                     belonging to each group
 *                     (size: mesh->n_groups * 3, interleaved)
 *   n_perio_faces --> number of faces belonging to each group
 *                     (size: mesh->n_init_perio)
 *----------------------------------------------------------------------------*/

static void
_prepare_mesh_group_stats(const cs_mesh_t  *mesh,
                          fvm_gnum_t        n_group_elts[],
                          fvm_gnum_t        n_perio_faces[])
{
  fvm_lnum_t i, j;

  int *i_face_flag = NULL;
  fvm_lnum_t *f_count = NULL;

  /* In case of parallel and/or periodic faces, build a flag so as to count
     purely parallel faces only once, but periodic faces twice.
     The flag is defined to the periodicity number for periodic faces,
     0 for regular or purely parallel faces with an "outside" ghost cell,
     -1 for purely parallel faces with an "inside" ghost cell. */

  if (mesh->halo != NULL) {

    BFT_MALLOC(i_face_flag, mesh->n_i_faces, int);

    if (mesh->n_init_perio > 0) {
      cs_mesh_get_face_perio_num(mesh, i_face_flag);
      for (i = 0; i < mesh->n_i_faces; i++)
        if (i_face_flag[i] < 0)
          i_face_flag[i] = - i_face_flag[i];
    }
    else {
      for (i = 0; i < mesh->n_i_faces; i++)
        i_face_flag[i] = 0;
    }

    for (i = 0; i < mesh->n_i_faces; i++) {
      fvm_lnum_t c_num_0 = mesh->i_face_cells[i*2];
      if (i_face_flag[i] == 0 && c_num_0 > mesh->n_cells)
        i_face_flag[i] = -1;
    }

  }

  /* Set counts to zero */

  BFT_MALLOC(f_count, mesh->n_families*3, fvm_lnum_t);
  for (i = 0; i < mesh->n_families*3; i++)
    f_count[i] = 0;

  for (i = 0; i < mesh->n_init_perio; i++)
    n_perio_faces[i] = 0;

  /* Cell counters */

  for (i = 0; i < mesh->n_cells; i++)
    f_count[(mesh->cell_family[i] - 1) * 3] += 1;

  /* Interior face counters */

  if (i_face_flag != NULL) {
    for (i = 0; i < mesh->n_i_faces; i++) {
      const int flag = i_face_flag[i];
      if (flag > 0)
        n_perio_faces[flag - 1] += 1;
      if (flag > -1)
        f_count[(mesh->i_face_family[i] - 1) * 3 + 1] += 1;
    }
    BFT_FREE(i_face_flag);
  }
  else {
    for (i = 0; i < mesh->n_i_faces; i++)
      f_count[(mesh->i_face_family[i] - 1) * 3 + 1] += 1;
  }

  /* Boundary face counters */

  if (mesh->n_b_faces > 0) {
    for (i = 0; i < mesh->n_b_faces; i++)
      f_count[(mesh->b_face_family[i] - 1) * 3 + 2] += 1;
  }

  /* Now transfer count from group classes to groups */

  for (i = 0; i < mesh->n_groups*3; i++)
    n_group_elts[i] = 0;

  for (i = 0; i < mesh->n_families; i++) {
    for (j = 0; j < mesh->n_max_family_items; j++) {
      int group_id = -1 - mesh->family_item[mesh->n_families*j + i];
      if (group_id >= 0) {
        int k;
        for (k = 0; k < 3; k++)
          n_group_elts[group_id*3 + k] += f_count[i*3 + k];
      }
    }
  }

  BFT_FREE(f_count);
}

/*----------------------------------------------------------------------------
 * Print mesh group statistics.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_print_mesh_group_stats(const cs_mesh_t  *mesh)
{
  int i;
  size_t count_size = mesh->n_groups*3 + mesh->n_init_perio;
  fvm_gnum_t *count = NULL, *n_elt_groups = NULL, *n_perio_faces = NULL;

  if (count_size == 0)
    return;

  BFT_MALLOC(count, count_size, fvm_gnum_t);

  n_elt_groups = count;
  n_perio_faces = count + mesh->n_groups*3;

  _prepare_mesh_group_stats(mesh, n_elt_groups, n_perio_faces);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    fvm_gnum_t *_count = NULL;
    BFT_MALLOC(_count, count_size, fvm_gnum_t);
    MPI_Allreduce(count, _count, count_size, FVM_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    memcpy(count, _count, sizeof(fvm_gnum_t)*count_size);
    BFT_FREE(_count);
  }
#endif

  /* Print perodic faces information */

  if (mesh->n_init_perio > 0) {


    bft_printf(_("\n Periodic faces (which are also interior faces):\n"));

    for (i = 0; i < mesh->n_init_perio; i++)
      bft_printf(_("     Periodicity %2d:          %llu face couples\n"),
                 i+1, (unsigned long long)(n_perio_faces[i]/2));

  }

  /* Print group information */

  if (mesh->n_groups > 0) {

    bft_printf(_("\n Groups:\n"));

    for (i = 0; i < mesh->n_groups; i++) {
      bft_printf("    \"%s\"\n", mesh->group_lst + (mesh->group_idx[i] - 1));
      if (n_elt_groups[i*3] > 0)
        bft_printf(_("       cells:          %12llu\n"),
                   (unsigned long long)(n_elt_groups[i*3]));
      if (n_elt_groups[i*3 + 1] > 0)
        bft_printf(_("       interior faces: %12llu\n"),
                   (unsigned long long)(n_elt_groups[i*3 + 1]));
      if (n_elt_groups[i*3 + 2] > 0)
        bft_printf(_("       boundary faces: %12llu\n"),
                   (unsigned long long)(n_elt_groups[i*3 + 2]));
    }

  }

  BFT_FREE(count);
}

/*----------------------------------------------------------------------------
 * Write a summary about mesh features in log file
 *
 * parameters:
 *   mesh                   <-- pointer to cs_mesh_t structure
 *----------------------------------------------------------------------------*/

static void
_print_mesh_info(cs_mesh_t  *mesh)
{
  cs_halo_t  *halo = mesh->halo;

  cs_int_t  *rank_buffer = NULL;

  /* Summary of cell and ghost cell distribution */

  if (mesh->n_domains > 1) {

    BFT_MALLOC(rank_buffer, mesh->n_domains, cs_int_t);

#if defined(HAVE_MPI)
    MPI_Allgather(&(mesh->n_cells), 1, CS_MPI_INT,
                  rank_buffer     , 1, CS_MPI_INT, cs_glob_mpi_comm);
#endif

    bft_printf(_("\n Histogram of the number of cells per rank:\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

  } /* End if n_domains > 1 */

  else
    bft_printf
      (_(" Number of cells:                                      %d\n"),
       mesh->n_cells);


  bft_printf("\n ----------------------------------------------------------\n");

  if (mesh->n_domains > 1) {

#if defined(HAVE_MPI)
    MPI_Allgather(&(mesh->n_cells_with_ghosts), 1, CS_MPI_INT,
                  rank_buffer, 1, CS_MPI_INT, cs_glob_mpi_comm);
#endif

    bft_printf(_("\n Histogram of the number of standard + halo cells "
                 "per rank:\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

  } /* End if n_domains > 1 */

  else
    bft_printf
      (_(" Number of cells + halo cells:                         %d\n\n"),
       mesh->n_cells_with_ghosts);


  bft_printf("\n ----------------------------------------------------------\n");

  if (halo != NULL) {

    cs_int_t  n_std_ghost_cells = halo->n_elts[CS_HALO_STANDARD];

    if (mesh->n_domains > 1) {

      cs_int_t  n_gcells = mesh->n_ghost_cells;

#if defined(HAVE_MPI)
      MPI_Allgather(&n_gcells  , 1, CS_MPI_INT,
                    rank_buffer, 1, CS_MPI_INT, cs_glob_mpi_comm);
#endif

      bft_printf
        (_("\n Histogram of the number of ghost cells per rank:\n\n"));

      _display_histograms(mesh->n_domains, rank_buffer);

      if (mesh->halo_type == CS_HALO_EXTENDED) {

#if defined(HAVE_MPI)
        MPI_Allgather(&n_std_ghost_cells, 1, CS_MPI_INT,
                      rank_buffer, 1, CS_MPI_INT, cs_glob_mpi_comm);
#endif

        bft_printf
          (_("\n"
             " Histogram of the number of ghost cells\n"
             " in the standard neighborhood per rank:\n\n"));

        _display_histograms(mesh->n_domains, rank_buffer);
      }

    } /* If n_ranks > 1 */

    else {
      bft_printf(_("\n Number of ghost cells:                          %10d\n"),
                 mesh->n_ghost_cells);
      if (mesh->halo_type == CS_HALO_EXTENDED)
        bft_printf(_("   in the standard neighborhood:              %10d\n"),
                   n_std_ghost_cells);
    }

    bft_printf("\n"
               " ----------------------------------------------------------\n");

  } /* End if halo != NULL */

  /* Summary of faces distribution */

  if (mesh->n_domains > 1) {

#if defined(HAVE_MPI)
    MPI_Allgather(&(mesh->n_i_faces), 1, CS_MPI_INT,
                  rank_buffer, 1, CS_MPI_INT, cs_glob_mpi_comm);
#endif

    bft_printf
      (_("\n Histogram of the number of interior faces per rank:\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

  } /* End if n_domains > 1 */

  else
    bft_printf
      (_(" Number of interior faces:                             %d\n"),
       mesh->n_i_faces);

  bft_printf("\n ----------------------------------------------------------\n");

  if (mesh->n_domains > 1) {

#if defined(HAVE_MPI)
    MPI_Allgather(&(mesh->n_b_faces), 1, CS_MPI_INT,
                  rank_buffer, 1, CS_MPI_INT, cs_glob_mpi_comm);
#endif

    bft_printf
      (_("\n Histogram of the number of boundary faces per rank:\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

  } /* End if n_domains > 1 */

  else
    bft_printf
      (_(" Number of boundary faces:                             %d\n"),
       mesh->n_b_faces);

  bft_printf("\n ----------------------------------------------------------\n");

  /* Add cell neighbor info */

  _print_cell_neighbor_info(mesh);

#if defined(HAVE_MPI)

  /* Summary of the number of neighbors */

  if (mesh->n_domains > 1) {

    cs_int_t  n_c_domains = halo->n_c_domains;

    MPI_Allgather(&n_c_domains, 1, CS_MPI_INT,
                  rank_buffer , 1, CS_MPI_INT, cs_glob_mpi_comm);

    bft_printf(_("\n Histogram of the number of neighboring domains "
                 "per rank:\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

    bft_printf("\n"
               " ----------------------------------------------------------\n");

  } /* If n_domains > 1 */

#endif

  /* Cleanup */

  if (mesh->n_domains > 1)
    BFT_FREE(rank_buffer);

  bft_printf_flush();
}

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the group number corresponding to a given name. If the group exists,
 * the number corresponds to the group rank (starting from 1) in the list of
 * the meshe's groups, multiplied by -1. This numbering is that used in
 * family (group class) description array IPRFML(NFML, NPRFML).
 *
 * If the group of the given name is not found, 9999 is returned.
 *
 * Fortran interface:
 *
 * FUNCTION NUMGRP (NAME, LEN)
 * ***************
 *
 * CHARACTER*       NAME        : <-- : Name of the group
 * INTEGER          LEN         : <-- : Group name length
 *----------------------------------------------------------------------------*/

cs_int_t CS_PROCF (numgrp, NUMGRP)
(
 const char       *name,    /* <-- Group name */
 const cs_int_t   *len      /* <-- Name length */
 CS_ARGF_SUPP_CHAINE        /*     (possible 'length' arguments added
                                   by many Fortran compilers) */
)
{
  int  i, lngcmp;

  char  *nomcmp = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  for (i = 0; i < mesh->n_groups; i++) {

    nomcmp = mesh->group_lst + (mesh->group_idx[i] - 1);
    lngcmp = strlen(nomcmp);

    if (lngcmp == *len && strncmp(nomcmp, name, lngcmp) == 0)
      return (-(i + 1));

  }

  return -9999;
}

/*----------------------------------------------------------------------------
 * Update a scalar array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synsca(var)
 * *****************
 *
 * var   : <-> : scalar array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synsca, SYNSCA)
(
 cs_real_t  var[]
)
{
  cs_mesh_sync_var_scal(var);
}

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synvec(var)
 * *****************
 *
 * var1   : <-> : vector component 1 array
 * var2   : <-> : vector component 2 array
 * var3   : <-> : vector component 3 array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synvec, SYNVEC)
(
 cs_real_t  var1[],
 cs_real_t  var2[],
 cs_real_t  var3[]
)
{
  cs_mesh_sync_var_vect(var1, var2, var3);
}

/*----------------------------------------------------------------------------
 * Update a diagonal tensor array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine syndia(var)
 * *****************
 *
 * var11   : <-> : diagonal tensor component 11 array
 * var22   : <-> : diagonal tensor component 22 array
 * var33   : <-> : diagonal tensor component 33 array
 *----------------------------------------------------------------------------*/

void CS_PROCF(syndia, SYNDIA)
(
 cs_real_t  var11[],
 cs_real_t  var22[],
 cs_real_t  var33[]
)
{
  cs_mesh_sync_var_diag(var11, var22, var33);
}

/*----------------------------------------------------------------------------
 * Update a tensor array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synten(var)
 * *****************
 *
 * var11   : <-> : tensor component 11 array
 * var12   : <-> : tensor component 12 array
 * var13   : <-> : tensor component 13 array
 * var21   : <-> : tensor component 21 array
 * var22   : <-> : tensor component 22 array
 * var23   : <-> : tensor component 23 array
 * var31   : <-> : tensor component 31 array
 * var32   : <-> : tensor component 32 array
 * var33   : <-> : tensor component 33 array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synten, SYNTEN)
(
 cs_real_t  var11[],
 cs_real_t  var12[],
 cs_real_t  var13[],
 cs_real_t  var21[],
 cs_real_t  var22[],
 cs_real_t  var23[],
 cs_real_t  var31[],
 cs_real_t  var32[],
 cs_real_t  var33[]
)
{
  cs_mesh_sync_var_tens(var11, var12, var13,
                        var21, var22, var23,
                        var31, var32, var33);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create an empty mesh structure
 *
 * returns:
 *   pointer to created mesh structure
 *----------------------------------------------------------------------------*/

cs_mesh_t *
cs_mesh_create(void)
{
  cs_mesh_t  *mesh = NULL;

  BFT_MALLOC (mesh, 1, cs_mesh_t);

  /* General features */

  mesh->dim = 3;
  mesh->domain_num = cs_glob_rank_id + 1;
  mesh->n_domains = 0;

  /* Local dimensions */

  mesh->n_cells = 0;
  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->n_vertices = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  /* Global dimensions */

  mesh->n_g_cells = 0;
  mesh->n_g_i_faces = 0;
  mesh->n_g_b_faces = 0;
  mesh->n_g_vertices = 0;

  /* Local structures */

  mesh->vtx_coord = NULL;
  mesh->i_face_cells = NULL;
  mesh->b_face_cells = NULL;
  mesh->i_face_vtx_idx = NULL;
  mesh->b_face_vtx_idx = NULL;
  mesh->i_face_vtx_lst = NULL;
  mesh->b_face_vtx_lst = NULL;

  /* Global numbering */

  mesh->global_cell_num = NULL;
  mesh->global_i_face_num = NULL;
  mesh->global_b_face_num = NULL;
  mesh->global_vtx_num = NULL;

  /* Periodic features */

  mesh->n_init_perio = 0;
  mesh->n_transforms = 0;
  mesh->have_rotation_perio = 0;
  mesh->periodicity = NULL;

  /* Halo features */

  mesh->halo_type = CS_HALO_N_TYPES;
  mesh->n_ghost_cells = 0;
  mesh->n_cells_with_ghosts = 0;
  mesh->halo = NULL;

  /* Extended connectivity and neighborhood */

  mesh->n_b_cells = 0;
  mesh->b_cells = NULL;

  mesh->cell_cells_idx = NULL;
  mesh->cell_cells_lst = NULL;
  mesh->gcell_vtx_idx = NULL;
  mesh->gcell_vtx_lst = NULL;

  /* Numbering info */

  mesh->i_face_numbering = NULL;
  mesh->b_face_numbering = NULL;

  /* Group and family features */

  mesh->n_groups = 0;
  mesh->group_idx = NULL;
  mesh->group_lst = NULL;

  mesh->n_max_family_items = 0;
  mesh->n_families = 0;

  mesh->family_item = NULL;
  mesh->cell_family = NULL;
  mesh->i_face_family = NULL;
  mesh->b_face_family = NULL;

  /* Selector features */

  mesh->class_defs = NULL;

  mesh->select_cells = NULL;
  mesh->select_i_faces = NULL;
  mesh->select_b_faces = NULL;

  /* Modification flag */

  mesh->modified = 0;

  return (mesh);
}

/*----------------------------------------------------------------------------
 * Create an empty mesh builder structure.
 *
 * returns:
 *   A pointer to a mesh builder structure
 *----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_create(void)
{
  cs_mesh_builder_t  *mb = NULL;

  BFT_MALLOC(mb, 1, cs_mesh_builder_t);

  mb->n_perio = 0;
  mb->n_perio_couples = NULL;
  mb->perio_couples = NULL;

  return mb;
}

/*----------------------------------------------------------------------------
 * Destroy a mesh structure.
 *
 * mesh       <->  pointer to a mesh structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_mesh_t *
cs_mesh_destroy(cs_mesh_t  *mesh)
{
  BFT_FREE(mesh->vtx_coord);
  BFT_FREE(mesh->i_face_cells);
  BFT_FREE(mesh->b_face_cells);
  BFT_FREE(mesh->i_face_vtx_idx);
  BFT_FREE(mesh->b_face_vtx_idx);
  BFT_FREE(mesh->i_face_vtx_lst);
  BFT_FREE(mesh->b_face_vtx_lst);

  BFT_FREE(mesh->global_cell_num);
  BFT_FREE(mesh->global_i_face_num);
  BFT_FREE(mesh->global_b_face_num);
  BFT_FREE(mesh->global_vtx_num);

  BFT_FREE(mesh->group_idx);
  BFT_FREE(mesh->group_lst);

  BFT_FREE(mesh->family_item);
  BFT_FREE(mesh->cell_family);
  BFT_FREE(mesh->i_face_family);
  BFT_FREE(mesh->b_face_family);

  /* Free periodic structures */

  if (mesh->n_init_perio > 0)
    mesh->periodicity = fvm_periodicity_destroy(mesh->periodicity);

  if (mesh->b_cells != NULL)
    BFT_FREE(mesh->b_cells);

  if (mesh->cell_cells_idx != NULL) {
    BFT_FREE(mesh->cell_cells_idx);
    BFT_FREE(mesh->cell_cells_lst);
  }

  if (mesh->gcell_vtx_idx != NULL) {
    BFT_FREE(mesh->gcell_vtx_idx);
    BFT_FREE(mesh->gcell_vtx_lst);
  }

  /* Free halo structure */

  if (mesh == cs_glob_mesh) {
    cs_perio_free_buffer();
    cs_halo_free_buffer();
  }

  mesh->halo = cs_halo_destroy(mesh->halo);

  /* Free numbering info */

  if (mesh->i_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->i_face_numbering));
  if (mesh->b_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->b_face_numbering));

  /* Free selection structures */

  if (mesh->n_groups > 0) {
    BFT_FREE(mesh->group_idx);
    BFT_FREE(mesh->group_lst);
  }

  if (mesh->select_cells != NULL)
    mesh->select_cells = fvm_selector_destroy(mesh->select_cells);
  if (mesh->select_i_faces != NULL)
    mesh->select_i_faces = fvm_selector_destroy(mesh->select_i_faces);
  if (mesh->select_b_faces != NULL)
    mesh->select_b_faces = fvm_selector_destroy(mesh->select_b_faces);

  /* Destroy group class set after selectors, who reference it */

  if (cs_glob_mesh->class_defs != NULL)
    cs_glob_mesh->class_defs
      = fvm_group_class_set_destroy(cs_glob_mesh->class_defs);

  BFT_FREE(mesh);

  return mesh;
}

/*----------------------------------------------------------------------------
 * Destroy a mesh builder structure
 *
 * mesh_builder <->  pointer to a mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_builder_destroy(cs_mesh_builder_t  **mesh_builder)
{
  if (mesh_builder != NULL) {

    int i;
    cs_mesh_builder_t *_mb = *mesh_builder;

    if (_mb->perio_couples != NULL) {
      for (i = 0; i < _mb->n_perio; i++)
        BFT_FREE(_mb->perio_couples[i]);
    }
    BFT_FREE(_mb->perio_couples);
    BFT_FREE(_mb->n_perio_couples);

    BFT_FREE(*mesh_builder);
  }
}

/*----------------------------------------------------------------------------
 * Renumber vertices.
 *
 * We ensure:
 * If i < j then mesh->global_vtx_num[i] < mesh->global_vtx_num[j]
 * which is not insured by the initial numbering from the pre-processor.
 *
 * parameters:
 *   mesh      <->  pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_order_vertices(cs_mesh_t  *const mesh)
{
  cs_int_t  i, j, size, dim, n_vertices;

  cs_int_t  *tmp_num = NULL;
  cs_real_t  *tmp_coord = NULL;
  fvm_lnum_t  *vertex_order = NULL;
  fvm_lnum_t  *vertex_renum = NULL;
  fvm_gnum_t  *g_vertex_num = NULL;

  assert(mesh != NULL);

  /* No treatment in serial mode */

  if (mesh->global_vtx_num == NULL)
    return;

  dim = mesh->dim;
  n_vertices = mesh->n_vertices;

  /* Compute the new vertex numbering */

  BFT_MALLOC(g_vertex_num, n_vertices, fvm_gnum_t);

  for (i = 0; i < n_vertices; i++)
    g_vertex_num[i] = mesh->global_vtx_num[i];

  vertex_order = fvm_order_local(NULL, g_vertex_num, (size_t)n_vertices);
  BFT_FREE(g_vertex_num);

  vertex_renum = fvm_order_local_renumbering(vertex_order, (size_t)n_vertices);
  BFT_FREE(vertex_order);

  /* Re-define face -> vertices connectivity with the new vertex numbering */

  if (mesh->n_i_faces > 0)
    for (i = 0, size = mesh->i_face_vtx_idx[mesh->n_i_faces] - 1; i < size; i++)
      mesh->i_face_vtx_lst[i] = vertex_renum[mesh->i_face_vtx_lst[i] - 1] + 1;

  if (mesh->n_b_faces > 0)
    for (i = 0, size = mesh->b_face_vtx_idx[mesh->n_b_faces] - 1; i < size; i++)
      mesh->b_face_vtx_lst[i] = vertex_renum[mesh->b_face_vtx_lst[i] - 1] + 1;

  /* Update coordinates */

  BFT_MALLOC(tmp_coord, n_vertices * dim, cs_real_t);

  for (i = 0; i < n_vertices; i++)
    for (j = 0; j < dim; j++)
      tmp_coord[vertex_renum[i]*dim + j] = mesh->vtx_coord[i*dim + j];

  memcpy(mesh->vtx_coord, tmp_coord, (n_vertices * dim * sizeof(cs_real_t)));
  BFT_FREE(tmp_coord);

  /* Update global numbering */

  BFT_MALLOC(tmp_num, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    tmp_num[vertex_renum[i]] = mesh->global_vtx_num[i];

  memcpy(mesh->global_vtx_num, tmp_num, (n_vertices * sizeof(cs_int_t)));

  /* Free memory */

  BFT_FREE(tmp_num);
  BFT_FREE(vertex_renum);
}

/*----------------------------------------------------------------------------
 * Compute global face connectivity size.
 *
 * Faces on simple parallel boundaries are counted only once, but periodic
 * faces are counted twice.
 *
 * parameters:
 *   mesh                   <-- pointer to a cs_mesh_t structure
 *   g_i_face_vertices_size --> global interior face connectivity size, or NULL
 *   g_b_face_vertices_size --> global boundary face connectivity size, or NULL
 *----------------------------------------------------------------------------*/

void
cs_mesh_g_face_vertices_sizes(const cs_mesh_t  *mesh,
                              fvm_gnum_t       *g_i_face_vertices_size,
                              fvm_gnum_t       *g_b_face_vertices_size)
{
  fvm_gnum_t  _g_face_vertices_size[2] = {0, 0};

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    fvm_lnum_t i;
    fvm_gnum_t _l_face_vertices_size[2] = {0, 0};

    /* Without periodicity, count faces bordering halo once */

    if (mesh->n_init_perio == 0) {
      for (i = 0; i < mesh->n_i_faces; i++) {
        if (mesh->i_face_cells[i*2] <= mesh->n_cells)
          _l_face_vertices_size[0] += (  mesh->i_face_vtx_idx[i+1]
                                       - mesh->i_face_vtx_idx[i]);
      }
    }

    /* With periodicity, count faces bordering halo once for purely parallel
       halo */

    else {

      fvm_lnum_t  rank_id, t_id, shift;
      fvm_lnum_t  start = 0, end = 0;
      int *perio_flag = NULL;

      const cs_halo_t *halo = mesh->halo;
      const cs_int_t  n_transforms = halo->n_transforms;

      BFT_MALLOC(perio_flag, mesh->n_ghost_cells, int);
      for (i = 0; i < mesh->n_ghost_cells; i++)
        perio_flag[i] = 0;

      for (t_id = 0; t_id < n_transforms; t_id++) {
        shift = 4 * halo->n_c_domains * t_id;
        for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
          start = halo->perio_lst[shift + 4*rank_id];
          end = start + halo->perio_lst[shift + 4*rank_id + 1];
          for (i = start; i < end; i++)
            perio_flag[i] = 1;
        }
      }

      for (i = 0; i < mesh->n_i_faces; i++) {
        int count_face = 1;
        if (mesh->i_face_cells[i*2] > mesh->n_cells) {
          if (perio_flag[mesh->i_face_cells[i*2] - mesh->n_cells - 1] == 0)
            count_face = 0;
        }
        if (count_face)
          _l_face_vertices_size[0] += (  mesh->i_face_vtx_idx[i+1]
                                       - mesh->i_face_vtx_idx[i]);
      }

      BFT_FREE(perio_flag);
    }

    _l_face_vertices_size[1] = mesh->b_face_vtx_connect_size;

    MPI_Allreduce(&_l_face_vertices_size, &_g_face_vertices_size, 2,
                  FVM_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
  }

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {
    _g_face_vertices_size[0] = mesh->i_face_vtx_connect_size;
    _g_face_vertices_size[1] = mesh->b_face_vtx_connect_size;
  }

  if (g_i_face_vertices_size != NULL)
    *g_i_face_vertices_size = _g_face_vertices_size[0];
  if (g_b_face_vertices_size != NULL)
    *g_b_face_vertices_size = _g_face_vertices_size[1];
}

/*----------------------------------------------------------------------------
 * Compute or update mesh structure members the depend on other members,
 * but whose results may be reused, such as global number of elements
 * (cells, vertices, internal and border faces) and sync cell family.
 *
 * parameters:
 *   mesh  <->  pointer to a cs_mesh_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_update_auxiliary(cs_mesh_t  *mesh)
{
  fvm_lnum_t  i;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    fvm_gnum_t  n_g_elts[4], max_elt_num[4];

    bft_printf(_("\n Global definition of the number of elements "
                 "(cells, vertices, faces...)\n"));

    /* Global dimensions of the mesh */

    max_elt_num[0] = mesh->n_cells;
    MPI_Allreduce(max_elt_num, n_g_elts, 1, FVM_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);

    max_elt_num[1] = 0;
    for (i = 0; i < mesh->n_i_faces; i++) {
      if (mesh->global_i_face_num[i] > max_elt_num[1])
        max_elt_num[1] = mesh->global_i_face_num[i];
    }

    max_elt_num[2] = 0;
    for (i = 0; i < mesh->n_b_faces; i++) {
      if (mesh->global_b_face_num[i] > max_elt_num[2])
        max_elt_num[2] = mesh->global_b_face_num[i];
    }

    max_elt_num[3] = 0;
    for (i = 0; i < mesh->n_vertices; i++) {
      if (mesh->global_vtx_num[i] > max_elt_num[3])
        max_elt_num[3] = mesh->global_vtx_num[i];
    }

    MPI_Allreduce(max_elt_num + 1, n_g_elts + 1, 3, FVM_MPI_GNUM, MPI_MAX,
                  cs_glob_mpi_comm);

    mesh->n_g_cells = n_g_elts[0];
    mesh->n_g_i_faces = n_g_elts[1];
    mesh->n_g_b_faces = n_g_elts[2];
    mesh->n_g_vertices = n_g_elts[3];

  }

#endif

  /* Sync cell family array (also in case of periodicity) */

  _sync_cell_fam(mesh);

  /* Number of boundary cells */

  {
    fvm_lnum_t n_b_cells = 0;
    _Bool *flag = NULL;

    BFT_MALLOC(flag, mesh->n_cells, _Bool);

    for (i = 0; i < mesh->n_cells; i++)
      flag[i] = false;

    for (i = 0; i < mesh->n_b_faces; i++)
      flag[mesh->b_face_cells[i] - 1] = true;

    for (i = 0, n_b_cells = 0; i < mesh->n_cells; i++) {
      if (flag[i] == true)
        n_b_cells++;
    }

    mesh->n_b_cells = n_b_cells;
    BFT_MALLOC(mesh->b_cells, mesh->n_b_cells, cs_int_t);

    for (i = 0, n_b_cells = 0; i < mesh->n_cells; i++) {
      if (flag[i] == true)
        mesh->b_cells[n_b_cells++] = i + 1;
    }

    BFT_FREE(flag);
  }
}

/*----------------------------------------------------------------------------
 * Creation and initialization of halo structures.
 *
 * Treatment of parallel and/or periodic halos for standard and extended
 * ghost cells according to halo type requested by global options.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure
 *   mb    <->  pointer to mesh builder (in case of periodicity)
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_halo(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb)
{
  int  ivoset;
  fvm_lnum_t  i;
  double  t1, t2;
  double  halo_time = 0, interface_time = 0, ext_neighborhood_time = 0;

  int  *perio_num = NULL;
  fvm_lnum_t  *gcell_vtx_idx = NULL, *gcell_vtx_lst = NULL;
  fvm_lnum_t  *n_periodic_couples = NULL;
  fvm_gnum_t  *g_i_face_num = NULL;
  fvm_gnum_t  *g_vertex_num = NULL;
  fvm_gnum_t  **periodic_couples = NULL;

  fvm_interface_set_t  *face_interfaces = NULL;
  fvm_interface_set_t  *vertex_interfaces = NULL;

  const fvm_lnum_t  n_i_faces = mesh->n_i_faces;
  const fvm_lnum_t  n_vertices = mesh->n_vertices;

  /* Choose the type of halo to build according to Fortran options.
     IMRGRA == 3 or 2 OR ITURB == 41 => CS_MESH_HALO_EXTENDED */

  CS_PROCF (haltyp, HALTYP) (&ivoset);

  /* Build halo */

  if (mesh->n_domains > 1 || mesh->n_init_perio > 0) {

    t1 = bft_timer_wtime();

    bft_printf("\n"
               " ----------------------------------------------------------\n");

    if (ivoset == 1) {

      bft_printf(_("\n Halo construction with extended neighborhood\n"
                   " ============================================\n\n"));

      mesh->halo_type = CS_HALO_EXTENDED;

      if (mesh->n_init_perio > 1) {
        bft_printf(_(" Composing periodicities\n"));
        fvm_periodicity_combine(mesh->periodicity, 0);
      }

    }
    else {

      bft_printf(_("\n Halo construction with standard neighborhood\n"
                   " ============================================\n\n"));

      mesh->halo_type = CS_HALO_STANDARD;
    }

    /* Build periodic numbering */

    if (mesh->n_init_perio > 0) {
      BFT_MALLOC(perio_num, mesh->n_init_perio, int);
      for (i = 0; i < mesh->n_init_perio; i++)
        perio_num[i] = i+1;
    }

    bft_printf(_(" Face interfaces creation\n"));

    /* Build purely parallel fvm_interface_set structure */

    g_i_face_num = mesh->global_i_face_num;
    g_vertex_num = mesh->global_vtx_num;

    if (g_i_face_num == NULL) {
      BFT_MALLOC(g_i_face_num, n_i_faces, fvm_gnum_t);
      for (i = 0; i < n_i_faces; i++)
        g_i_face_num[i] = (fvm_gnum_t)i+1;
    }

    if (mb != NULL) {
      face_interfaces = fvm_interface_set_create(n_i_faces,
                                                 NULL,
                                                 g_i_face_num,
                                                 mesh->periodicity,
                                                 mesh->n_init_perio,
                                                 perio_num,
                                                 mb->n_perio_couples,
                      (const fvm_gnum_t *const *)(mb->perio_couples));
      for (i = 0; i < mb->n_perio; i++)
        BFT_FREE(mb->perio_couples[i]);
      BFT_FREE(mb->perio_couples);
      BFT_FREE(mb->n_perio_couples);
      mb->n_perio = 0;
    }
    else
      face_interfaces = fvm_interface_set_create(n_i_faces,
                                                 NULL,
                                                 g_i_face_num,
                                                 mesh->periodicity,
                                                 0,
                                                 NULL,
                                                 NULL,
                                                 NULL);

    if (mesh->global_i_face_num != g_i_face_num)
      BFT_FREE(g_i_face_num);

    if (g_vertex_num == NULL) {
      BFT_MALLOC(g_vertex_num, n_vertices, fvm_gnum_t);
      for (i = 0; i < n_vertices; i++)
        g_vertex_num[i] = (fvm_gnum_t)i+1;
    }

    if (mesh->n_init_perio > 0) {

      mesh->n_transforms = fvm_periodicity_get_n_transforms(mesh->periodicity);

      bft_printf(_(" Definition of periodic vertices\n"));

      _define_perio_vtx_couples(mesh,
                                face_interfaces,
                                &n_periodic_couples,
                                &periodic_couples);

#if 0 /* For debugging purposes */
      for (i = 0; i < mesh->n_init_perio; i++) {
        cs_int_t  j;
        bft_printf("\n\n  Periodicity number: %d\n", perio_num[i]);
        bft_printf("  Number of couples : %d\n", n_periodic_couples[i]);
        for (j = 0; j < n_periodic_couples[i]; j++)
          bft_printf("%12d --> %12d\n",
                     periodic_couples[i][2*j], periodic_couples[i][2*j + 1]);
      }
      fvm_periodicity_dump(mesh->periodicity);
#endif

    }

    bft_printf(_(" Vertex interfaces creation\n"));

    vertex_interfaces = fvm_interface_set_create(n_vertices,
                                                 NULL,
                                                 g_vertex_num,
                                                 mesh->periodicity,
                                                 mesh->n_init_perio,
                                                 perio_num,
                                                 n_periodic_couples,
                      (const fvm_gnum_t *const *)periodic_couples);

    if (mesh->global_vtx_num != g_vertex_num)
      BFT_FREE(g_vertex_num);

    BFT_FREE(perio_num);
    BFT_FREE(n_periodic_couples);

    for (i = 0; i < mesh->n_init_perio; i++)
      BFT_FREE(periodic_couples[i]);
    BFT_FREE(periodic_couples);

#if 0 /* For debugging purposes */
    bft_printf("Dump final vertices interface:\n");
    fvm_interface_set_dump(vertex_interfaces);
#endif

    t2 = bft_timer_wtime();
    interface_time = t2-t1;

    t1 = bft_timer_wtime();

    /* Creation of the cs_halo_t structure. */

    bft_printf(_(" Halo creation\n"));
    bft_printf_flush();

    mesh->halo = cs_halo_create(vertex_interfaces);

    bft_printf(_(" Halo definition\n"));
    bft_printf_flush();

    cs_mesh_halo_define(mesh,
                        face_interfaces,
                        vertex_interfaces,
                        &gcell_vtx_idx,
                        &gcell_vtx_lst);

    if (mesh->n_init_perio > 0)
      cs_perio_update_buffer(mesh->halo);

    fvm_interface_set_destroy(vertex_interfaces);
    fvm_interface_set_destroy(face_interfaces);

    t2 = bft_timer_wtime();
    halo_time = t2-t1;

  } /* end if (mesh->n_domains > 1 || mesh->n_init_perio > 0) */

  /* Define a cell -> cells connectivity for the extended neighborhood
     if necessary */

  if (ivoset == 1) {

    t1 = bft_timer_wtime();
    bft_printf(_(" Extended neighborhood structures definition\n"));
    bft_printf_flush();

    mesh->gcell_vtx_idx = gcell_vtx_idx;
    mesh->gcell_vtx_lst = gcell_vtx_lst;

    cs_ext_neighborhood_define(mesh);

    bft_printf_flush();
    t2 = bft_timer_wtime();
    ext_neighborhood_time = t2-t1;

  }
  else {
    BFT_FREE(gcell_vtx_idx);
    BFT_FREE(gcell_vtx_lst);
  }

  /* Output for listing */

  if (mesh->halo_type != CS_HALO_N_TYPES)
    _print_halo_info(mesh,
                     interface_time,
                     halo_time,
                     ext_neighborhood_time);

  else if (ivoset == 1)
    bft_printf(_("\n Extended connectivity creation (%.3g s)\n"),
               ext_neighborhood_time);

  _print_mesh_info(mesh);
}

/*----------------------------------------------------------------------------
 * Get the global number of ghost cells.
 *
 * parameters:
 *   mesh  <--  pointer to a mesh structure
 *
 * returns:
 *   global number of ghost cells
 *---------------------------------------------------------------------------*/

cs_int_t
cs_mesh_n_g_ghost_cells(cs_mesh_t  *mesh)
{
  cs_int_t  n_g_ghost_cells = 0;

  if (cs_glob_n_ranks == 1)
    n_g_ghost_cells = mesh->n_ghost_cells;

  else {

    assert(cs_glob_n_ranks > 1);

#if defined(HAVE_MPI)
    MPI_Allreduce(&(mesh->n_ghost_cells), &n_g_ghost_cells, 1, MPI_INT,
                  MPI_SUM, cs_glob_mpi_comm);
#endif

  }

  return n_g_ghost_cells;
}

/*----------------------------------------------------------------------------
 * Order family numbers and remove duplicates
 *
 * parameters
 *   mesh <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_clean_families(cs_mesh_t  *mesh)
{
  size_t i, j, gc_id, gc_id_prev;
  int max_val = 0;
  size_t gc_count = 0;
  size_t n_gc = mesh->n_families;
  size_t n_gc_vals = mesh->n_max_family_items;
  size_t size_tot = n_gc * n_gc_vals;
  fvm_gnum_t *interlaced = NULL;
  int *order = NULL, *renum = NULL;

  if (mesh->n_families < 2)
    return;

  /* Build and order interlaced copy with only positive values */

  BFT_MALLOC(interlaced, size_tot, fvm_gnum_t);

  for (i = 0; i < size_tot; i++) {
    int val = mesh->family_item[i];
    if (val > max_val)
      max_val = val;
  }

  for (i = 0; i < n_gc; i++) {
    for (j = 0; j < n_gc_vals; j++) {
      int val = mesh->family_item[j*n_gc + i];
      if (val < 0)
        val = -val + max_val;
      interlaced[i*n_gc_vals + j] = val;
    }
  }

  order = fvm_order_local_s(NULL, interlaced, n_gc_vals, n_gc);

  /* Prepare removal of duplicates and renumbering */

  BFT_MALLOC(renum, n_gc, int);

  gc_id = order[0];
  gc_id_prev = gc_id;
  gc_count = 1;
  renum[gc_id] = 0;

  for (i = 1; i < n_gc; i++) {
    char is_same = '\1';
    gc_id = order[i];
    for (j = 0; j < n_gc_vals; j++) {
      if (   interlaced[gc_id_prev*n_gc_vals + j]
          != interlaced[gc_id*n_gc_vals + j])
        is_same = '\0';
    }
    if (is_same != '\1') {
      gc_id_prev = gc_id;
      gc_count += 1;
    }
    renum[gc_id] = gc_count - 1;
  }

  /* Update definitions */

  mesh->n_families = gc_count;
  BFT_REALLOC(mesh->family_item, gc_count*n_gc_vals, cs_int_t);

  for (i = 0; i < n_gc; i++) {
    gc_id = renum[i];
    for (j = 0; j < n_gc_vals; j++)
      mesh->family_item[j*gc_count + gc_id] = interlaced[i*n_gc_vals + j];
  }

  size_tot = gc_count * n_gc_vals;
  for (i = 0; i < size_tot; i++) {
    int val = mesh->family_item[i];
    if (val > max_val)
      val = -(val - max_val);
    mesh->family_item[i] = val;
  }

  BFT_FREE(interlaced);
  BFT_FREE(order);

  /* Update references */

  if (mesh->cell_family != NULL) {
    for (i = 0; i < (size_t)(mesh->n_cells); i++) {
      int val = mesh->cell_family[i];
      if (val != 0)
        mesh->cell_family[i] = renum[val -1] + 1;
    }
  }

  if (mesh->i_face_family != NULL) {
    for (i = 0; i < (size_t)(mesh->n_i_faces); i++) {
      int val = mesh->i_face_family[i];
      if (val != 0)
        mesh->i_face_family[i] = renum[val - 1] + 1;
    }
  }

  if (mesh->b_face_family != NULL) {
    for (i = 0; i < (size_t)(mesh->n_b_faces); i++) {
      int val = mesh->b_face_family[i];
      if (val != 0)
        mesh->b_face_family[i] = renum[val - 1] + 1;
    }
  }

  BFT_FREE(renum);
}

/*----------------------------------------------------------------------------
 * Add colors to group class descriptions for group names which are
 * convertible to integers.
 *
 * parameters
 *   mesh <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_build_colors(cs_mesh_t  *mesh)
{
  cs_int_t  i, j;
  int  n_max_family_items = mesh->n_max_family_items;
  int  n_colors = 0;

  if (mesh->n_groups == 0)
    return;

  /* Counting pass */

  for (i = 0; i < mesh->n_families; i++) {
    int n_add = 0;
    int n_family_colors = 0;
    int n_null = 0;
    for (j = 0; j < mesh->n_max_family_items; j++) {
      int gc_prop = mesh->family_item[mesh->n_families*j + i];
      if (gc_prop < 0) {
        int color_id;
        int group_id = - mesh->family_item[mesh->n_families*j + i] - 1;
        const char *group_name = mesh->group_lst + mesh->group_idx[group_id] - 1;
        if (sscanf(group_name, "%d", &color_id) == 1) {
          n_colors += 1;
          n_family_colors += 1;
        }
      }
      else if (gc_prop == 0)
        n_null += 1;
    }
    n_add = n_family_colors - n_null;
    if (n_add > 0)
      n_max_family_items = CS_MAX(mesh->n_max_family_items + n_add,
                                  n_max_family_items);
  }

  /* Update structure */

  if (n_colors > 0) {

    int *family_items = NULL;

    /* Increase maximum number of definitions and pad it necessary */

    if (n_max_family_items > mesh->n_max_family_items) {
      BFT_REALLOC(mesh->family_item,
                  mesh->n_families*n_max_family_items,
                  cs_int_t);
      for (i = mesh->n_max_family_items;
           i < n_max_family_items;
           i++) {
        for (j = 0; j < mesh->n_families; j++)
          mesh->family_item[mesh->n_families*i + j] = 0;
      }
      mesh->n_max_family_items = n_max_family_items;
    }

    /* Now insert color numbers before group numbers */

    BFT_MALLOC(family_items, mesh->n_max_family_items, int);

    for (i = 0; i < mesh->n_families; i++) {

      int n_family_items = 0;
      int n_family_colors = 0;

      for (j = 0; j < mesh->n_max_family_items; j++) {
        int gc_prop = mesh->family_item[mesh->n_families*j + i];
        if (gc_prop != 0)
          family_items[n_family_items++] = gc_prop;
        if (gc_prop < 0) {
          int color_id;
          int group_id = - mesh->family_item[mesh->n_families*j + i] - 1;
          const char *group_name
            = mesh->group_lst + mesh->group_idx[group_id] - 1;
          if (sscanf(group_name, "%d", &color_id) == 1) {
            mesh->family_item[mesh->n_families*n_family_colors + i] = color_id;
            n_family_colors += 1;
          }
        }
      }

      if (n_family_colors > 0) {
        for (j = 0; j < n_family_items; j++)
          mesh->family_item[mesh->n_families*(n_family_colors+j) + i]
            = family_items[j];
      }

    }

    BFT_FREE(family_items);

  }
}

/*----------------------------------------------------------------------------
 * Define group classes for a mesh based on its family definitions.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_group_classes(cs_mesh_t  *mesh)
{
  int  i, j;
  int  grp_nbr, grp_num, grp_idx;

  char **group = NULL;

  if (mesh->class_defs != NULL)
    mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

  mesh->class_defs = fvm_group_class_set_create();

  /* Construction of the fvm_group_class structure */

  BFT_MALLOC(group, mesh->n_max_family_items, char*);

  for (i = 0; i < mesh->n_families; i++) {

    grp_nbr  = 0;

    for (j = 0; j <  mesh->n_max_family_items; j++) {

      if (mesh->family_item[j * mesh->n_families + i] < 0) {
        /* Fortran formulation */
        grp_num = -mesh->family_item[j*mesh->n_families + i] -1;
        grp_idx = mesh->group_idx[grp_num];
        group[grp_nbr++] = mesh->group_lst + grp_idx -1;
      }

    }

    fvm_group_class_set_add(mesh->class_defs,
                            grp_nbr,
                            (const char **)group);

  } /* End of loop on families */

  BFT_FREE(group);
}

/*----------------------------------------------------------------------------
 * Assign selectors to global mesh.
 *
 * Should be called once the mesh is fully built.
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_selectors(void)
{
  if (cs_glob_mesh->class_defs == NULL)
    cs_mesh_init_group_classes(cs_glob_mesh);

  /* Construction of the selectors */

  cs_glob_mesh->select_cells
    = fvm_selector_create(cs_glob_mesh->dim,
                          cs_glob_mesh->n_cells,
                          cs_glob_mesh->class_defs,
                          cs_glob_mesh->cell_family,
                          1,
                          cs_glob_mesh_quantities->cell_cen,
                          NULL);

  cs_glob_mesh->select_b_faces
    = fvm_selector_create(cs_glob_mesh->dim,
                          cs_glob_mesh->n_b_faces,
                          cs_glob_mesh->class_defs,
                          cs_glob_mesh->b_face_family,
                          1,
                          cs_glob_mesh_quantities->b_face_cog,
                          cs_glob_mesh_quantities->b_face_normal);

  cs_glob_mesh->select_i_faces
    = fvm_selector_create(cs_glob_mesh->dim,
                          cs_glob_mesh->n_i_faces,
                          cs_glob_mesh->class_defs,
                          cs_glob_mesh->i_face_family,
                          1,
                          cs_glob_mesh_quantities->i_face_cog,
                          cs_glob_mesh_quantities->i_face_normal);

}

/*----------------------------------------------------------------------------
 * Update a scalar array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var  <->  scalar array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_scal(cs_real_t  *var)
{
  const cs_halo_t  *halo = cs_glob_mesh->halo;

  if (halo == NULL) return;

  cs_halo_sync_var(halo, CS_HALO_STANDARD, var);

  if (cs_glob_mesh->n_init_perio > 0)
    cs_perio_sync_var_scal(halo,
                           CS_HALO_STANDARD,
                           CS_PERIO_ROTA_COPY,
                           var);
}

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var1  <->  vector component 1 array
 *   var2  <->  vector component 2 array
 *   var3  <->  vector component 3 array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_vect(cs_real_t  *var1,
                      cs_real_t  *var2,
                      cs_real_t  *var3)
{
  const cs_halo_t  *halo = cs_glob_mesh->halo;

  if (halo == NULL) return;

  cs_halo_sync_var(halo, CS_HALO_STANDARD, var1);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var2);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var3);

  if (cs_glob_mesh->n_init_perio > 0)
    cs_perio_sync_var_vect(halo,
                           CS_HALO_STANDARD,
                           CS_PERIO_ROTA_COPY,
                           var1, var2, var3);
}

/*----------------------------------------------------------------------------
 * Update a diagonal tensor array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var11  <->  diagonal tensor component 11 array
 *   var22  <->  diagonal tensor component 22 array
 *   var33  <->  diagonal tensor component 33 array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_diag(cs_real_t  *var11,
                      cs_real_t  *var22,
                      cs_real_t  *var33)
{
  const cs_halo_t  *halo = cs_glob_mesh->halo;

  if (halo == NULL) return;

  cs_halo_sync_var(halo, CS_HALO_STANDARD, var11);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var22);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var33);

  if (cs_glob_mesh->n_init_perio > 0)
    cs_perio_sync_var_diag(halo,
                           CS_HALO_STANDARD,
                           var11, var22, var33);
}

/*----------------------------------------------------------------------------
 * Update a tensor array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var11  <->  tensor component 11 array
 *   var12  <->  tensor component 12 array
 *   var13  <->  tensor component 13 array
 *   var21  <->  tensor component 21 array
 *   var22  <->  tensor component 22 array
 *   var23  <->  tensor component 23 array
 *   var31  <->  tensor component 31 array
 *   var32  <->  tensor component 32 array
 *   var33  <->  tensor component 33 array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_tens(cs_real_t  *var11,
                      cs_real_t  *var12,
                      cs_real_t  *var13,
                      cs_real_t  *var21,
                      cs_real_t  *var22,
                      cs_real_t  *var23,
                      cs_real_t  *var31,
                      cs_real_t  *var32,
                      cs_real_t  *var33)
{
  const cs_halo_t  *halo = cs_glob_mesh->halo;

  if (halo == NULL) return;

  cs_halo_sync_var(halo, CS_HALO_STANDARD, var11);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var12);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var13);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var21);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var22);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var23);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var31);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var32);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, var33);

  if (cs_glob_mesh->n_init_perio > 0)
    cs_perio_sync_var_tens(halo,
                           CS_HALO_STANDARD,
                           var11, var12, var13,
                           var21, var22, var33,
                           var31, var32, var33);
}

/*----------------------------------------------------------------------------
 * Get global lists of periodic face couples.
 *
 * In parallel, each face couple may appear on only one rank.
 *
 * The caller is responsible for freeing the arrays allocated and returned
 * by this function once they are no onger needed.
 *
 * parameters:
 *   mesh                 <-- pointer to mesh structure
 *   n_perio_face_couples --> global number of periodic couples per
 *                            periodicity (size: mesh->n_init_perio)
 *   perio_face_couples   --> arrays of global periodic couple face numbers,
 *                            for each periodicity
 *----------------------------------------------------------------------------*/

void
cs_mesh_get_perio_faces(const cs_mesh_t    *mesh,
                        fvm_lnum_t        **n_perio_face_couples,
                        fvm_gnum_t       ***perio_face_couples)
{
  assert (mesh != NULL);

  if (mesh->n_init_perio == 0)
    return;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _get_perio_faces_g(mesh,
                       n_perio_face_couples,
                       perio_face_couples);

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1)
    _get_perio_faces_l(mesh,
                       n_perio_face_couples,
                       perio_face_couples);

  /* Ensure couples are sorted (not normally necessary, but safe) */

  if (n_perio_face_couples != NULL) {
    int i;
    for (i = 0; i < mesh->n_init_perio; i++) {
      if ((*n_perio_face_couples)[i] > 1)
        qsort((*perio_face_couples)[i],
              (*n_perio_face_couples)[i],
              sizeof(fvm_gnum_t) * 2,
              &_compare_couples);
    }
  }
}

/*----------------------------------------------------------------------------
 * Build global cell numbering array extended to ghost cell values.
 *
 * If the blank_perio flag is nonzero, periodic ghost cell numbers
 * are set to zero instead of the value of the matching cell.
 *
 * The caller is responsible for freeing the returned array when it
 * is no longer useful.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   blank_perio <-- flag to zeroe periodic cell values
 *----------------------------------------------------------------------------*/

fvm_gnum_t *
cs_mesh_get_cell_gnum(const cs_mesh_t  *mesh,
                      int               blank_perio)
{
  fvm_lnum_t i;
  fvm_gnum_t *cell_gnum = NULL;

  assert (mesh != NULL);

  /* Allocate array */

  BFT_MALLOC(cell_gnum, mesh->n_cells_with_ghosts, fvm_gnum_t);

  /* Build global cell numbering including parallel halos */

  for (i = 0; i < mesh->n_cells; i++)
    cell_gnum[i] = mesh->global_cell_num[i];
  for (i = mesh->n_cells; i < mesh->n_cells_with_ghosts; i++)
    cell_gnum[i] = 0;

  if (mesh->halo != NULL) {

    cs_halo_sync_untyped(mesh->halo,
                         CS_HALO_EXTENDED,
                         sizeof(fvm_gnum_t),
                         cell_gnum);

    if (blank_perio) {

      const cs_halo_t *halo = mesh->halo;

      fvm_lnum_t  rank_id, t_id, shift;
      fvm_lnum_t  start = 0, end = 0;

      const cs_int_t  n_transforms = halo->n_transforms;
      const cs_int_t  n_elts = halo->n_local_elts;

      for (t_id = 0; t_id < n_transforms; t_id++) {

        shift = 4 * halo->n_c_domains * t_id;

        for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

          start = halo->perio_lst[shift + 4*rank_id];
          end = start + halo->perio_lst[shift + 4*rank_id + 1];
          for (i = start; i < end; i++)
            cell_gnum[n_elts+i] = 0;

          start = halo->perio_lst[shift + 4*rank_id + 2];
          end = start + halo->perio_lst[shift + 4*rank_id + 3];
          for (i = start; i < end; i++)
            cell_gnum[n_elts+i] = 0;

        } /* End of loop on ranks */

      } /* End of loop on transformation */

    } /* End for blank_perio */

  }

  /* Return global cell number */

  return cell_gnum;
}

/*----------------------------------------------------------------------------
 * Mark interior faces with the number of their associated periodic
 * transform id.
 *
 * parameters:
 *   mesh      <-- pointer to mesh structure
 *   perio_num --> periodicity number associated with each face, signed for
 *                 direct/reverse transform, 0 for non-periodic faces
 *                 (size: mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_mesh_get_face_perio_num(const cs_mesh_t  *mesh,
                           int               perio_num[])
{
  fvm_lnum_t i;

  for (i = 0; i < mesh->n_i_faces; i++)
    perio_num[i] = 0;

  if (mesh->n_init_perio > 0) {

    int *halo_perio_num = NULL;

    /* Mark ghost cells with their periodicity number */

    BFT_MALLOC(halo_perio_num, mesh->n_ghost_cells, int);

    _get_halo_perio_num(mesh, halo_perio_num, NULL);

    for (i = 0; i < mesh->n_i_faces; i++) {
      const fvm_lnum_t h_id0 = mesh->i_face_cells[i*2] - mesh->n_cells - 1;
      const fvm_lnum_t h_id1 = mesh->i_face_cells[i*2+1] - mesh->n_cells - 1;
      if (h_id0 >= 0) {
        if (halo_perio_num[h_id0] != 0)
          perio_num[i] = halo_perio_num[h_id0];
      }
      else if (h_id1 >= 0) {
        if (halo_perio_num[h_id1] != 0)
          perio_num[i] = halo_perio_num[h_id1];
      }
    }

    BFT_FREE(halo_perio_num);
  }
}

/*----------------------------------------------------------------------------
 * Print information on a mesh structure.
 *
 * parameters:
 *   mesh  <--  pointer to mesh structure.
 *   name  <--  associated name.
 *----------------------------------------------------------------------------*/

void
cs_mesh_print_info(const cs_mesh_t  *mesh,
                   const char       *name)
{
  if (mesh->n_g_vertices > 0) {

    cs_int_t  i, vtx_id;
    cs_int_t  dim = mesh->dim;
    cs_real_t  min_xyz[3] = { 1.e127,  1.e127,  1.e127};
    cs_real_t  max_xyz[3] = {-1.e127, -1.e127, -1.e127};

    for (vtx_id = 0 ; vtx_id < mesh->n_vertices ; vtx_id++) {

      for (i = 0 ; i < dim ; i++) {

        if (mesh->vtx_coord[vtx_id*dim + i] < min_xyz[i])
          min_xyz[i] = mesh->vtx_coord[vtx_id*dim + i];

        if (mesh->vtx_coord[vtx_id*dim + i] > max_xyz[i])
          max_xyz[i] = mesh->vtx_coord[vtx_id*dim + i];

      }

    } /* End of loop on vertices */

#if defined(HAVE_MPI)

    if (cs_glob_n_ranks > 1) {
      cs_real_t  g_min_xyz[3];
      cs_real_t  g_max_xyz[3];
      MPI_Allreduce(min_xyz, g_min_xyz, dim, CS_MPI_REAL, MPI_MIN,
                    cs_glob_mpi_comm);
      MPI_Allreduce(max_xyz, g_max_xyz, dim, CS_MPI_REAL, MPI_MAX,
                    cs_glob_mpi_comm);
      for (i = 0 ; i < dim ; i++) {
        min_xyz[i] = g_min_xyz[i];
        max_xyz[i] = g_max_xyz[i];
      }
    }

#endif

    bft_printf(_("\n"
                 " Mesh coordinates:               minimum    and maximum\n"
                 "                       X : %14.7e %14.7e\n"
                 "                       Y : %14.7e %14.7e\n"
                 "                       Z : %14.7e %14.7e\n"),
               min_xyz[0], max_xyz[0], min_xyz[1], max_xyz[1],
               min_xyz[2], max_xyz[2]);

  }

  bft_printf(_(" %s\n"
               "     Number of cells:          %llu\n"
               "     Number of interior faces: %llu\n"
               "     Number of boundary faces: %llu\n"
               "     Number of vertices:       %llu\n"),
             name,
             (unsigned long long)(mesh->n_g_cells),
             (unsigned long long)(mesh->n_g_i_faces),
             (unsigned long long)(mesh->n_g_b_faces),
             (unsigned long long)(mesh->n_g_vertices));

  _print_mesh_group_stats(mesh);
}

/*----------------------------------------------------------------------------
 * Dump of a mesh structure.
 *
 * parameters:
 *   mesh  <--  pointer to mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_dump(const cs_mesh_t  *mesh)
{
  cs_int_t  i, j;

  bft_printf("\n\nDUMP OF THE MESH STRUCTURE: %p\n\n",mesh);

  bft_printf("space dim :        %d\n"
             "n_domains :        %d\n",
             "domain_num:        %d\n",
             mesh->dim, mesh->n_domains, mesh->domain_num);

  bft_printf("\nLocal dimensions:\n"
             "n_cells:                  %d\n"
             "n_cells_with_ghosts:      %d\n"
             "n_vertices:               %d\n"
             "n_i_faces:                %d\n"
             "n_b_faces:                %d\n",
             mesh->n_cells, mesh->n_cells_with_ghosts,
             mesh->n_vertices,
             mesh->n_i_faces, mesh->n_b_faces);

  bft_printf("\nGlobal dimensions:\n"
             "n_g_cells:                %llu\n"
             "n_g_vertices:             %llu\n"
             "n_g_i_faces:              %llu\n"
             "n_g_b_faces:              %llu\n",
             (unsigned long long)mesh->n_g_cells,
             (unsigned long long)mesh->n_g_vertices,
             (unsigned long long)mesh->n_g_i_faces,
             (unsigned long long)mesh->n_g_b_faces);

  bft_printf("\n\n        --------"
             "        Vertices"
             "        --------\n\n");

  bft_printf("\nVertex coordinates:\n");
  for (i = 0; i < mesh->n_vertices; i++)
    bft_printf("   <%3d >  %10.3f        %10.3f        %10.3f\n",
               i+1, mesh->vtx_coord[3*i], mesh->vtx_coord[3*i+1],
               mesh->vtx_coord[3*i+2]);

  if (mesh->n_domains > 1) {
    bft_printf("\nGlobal vertex numbering:\n");
    for (i = 0; i < mesh->n_vertices; i++)
      bft_printf("   <%7d >  %10u\n",
                 i+1, mesh->global_vtx_num[i]);
  }

  bft_printf("\n\n        ---------------------------"
             "        Internal faces connectivity"
             "        ---------------------------\n\n");

  bft_printf("\nInternal faces -> Cells connectivity:\n");
  for (i = 0; i < mesh->n_i_faces; i++)
    bft_printf("   < %7d >  %7d  <---->  %7d\n", i+1,
               mesh->i_face_cells[2*i], mesh->i_face_cells[2*i+1]);

  bft_printf("\nInternal faces -> vertices connectivity:\n");
  for (i = 0; i < mesh->n_i_faces; i++) {
    bft_printf("    < %7d >", i+1);
    for (j = mesh->i_face_vtx_idx[i]-1; j < mesh->i_face_vtx_idx[i+1]-1; j++)
      bft_printf("  %7d ",mesh->i_face_vtx_lst[j]);
    bft_printf("\n");
  }

  if (mesh->global_i_face_num != NULL) {

    bft_printf("\nInternal faces global numbering:\n");
    for (i = 0; i < mesh->n_i_faces; i++)
      bft_printf("   < %7d >  %12d",
                 i+1, mesh->global_i_face_num[i]);
    bft_printf("\n");

  }

  bft_printf("\n\n        -------------------------"
             "        Border faces connectivity"
             "        -------------------------\n\n");

  bft_printf("\nBorder faces -> Cells connectivity:\n");
  for (i = 0; i < mesh->n_b_faces; i++)
    bft_printf("   < %7d >  %7d\n", i+1, mesh->b_face_cells[i]);

  bft_printf("\nBorder faces -> vertices connectivity:\n");
  for (i = 0; i < mesh->n_b_faces; i++) {
    bft_printf("   < %7d >", i+1);
    for (j = mesh->b_face_vtx_idx[i]-1; j < mesh->b_face_vtx_idx[i+1]-1; j++)
      bft_printf("  %7d ",mesh->b_face_vtx_lst[j]);
    bft_printf("\n");
  }

  bft_printf("\n\n        -------------------------"
             "        Cells"
             "        -------------------------\n\n");

  if (mesh->global_cell_num != NULL) {

    bft_printf("\nCell global numbering:\n");
    for (i = 0; i < mesh->n_cells; i++)
      bft_printf("   < %7d >  %12d", i+1, mesh->global_cell_num[i]);
    bft_printf("\n");

  }

  bft_printf("\nNumber of families: %3d\n",mesh->n_families);
  bft_printf("Family of each cell:\n");
  for (i = 0; i < mesh->n_cells_with_ghosts; i++)
    bft_printf("   < %3d >  %5d\n", i+1, mesh->cell_family[i]);

  if (mesh->halo != NULL) {

    cs_halo_t  *halo = mesh->halo;

    bft_printf("\nHalo information: %p\n", halo);

    bft_printf("n_c_domains:              %d\n", halo->n_c_domains);
    bft_printf("n_ghost_cells:            %d\n", mesh->n_ghost_cells);
    bft_printf("n_std_ghost_cells:        %d\n",
               halo->n_elts[CS_HALO_STANDARD]);
    bft_printf("n_ext_ghost_cells:        %d\n",
               halo->n_elts[CS_HALO_EXTENDED] - halo->n_elts[CS_HALO_STANDARD]);

    for (i = 0; i < halo->n_c_domains; i++) {

      bft_printf("\n\nRank id:        %d\n"
                 "Halo index start:        %d        end:        %d\n"
                 "Send index start:        %d        end:        %d\n"
                 "Send cell numbers:\n",
                 halo->c_domain_rank[i],
                 halo->index[2*i], halo->index[2*i+2],
                 halo->send_index[2*i], halo->send_index[2*i+2]);
      for (j = halo->send_index[2*i]; j < halo->send_index[2*i+2]; j++)
        bft_printf("  %10d : %10d\n", j+1, halo->send_list[j]+1);

    } /* End of loop on the frontiers of halo */

    if (mesh->n_init_perio > 0 && halo->perio_lst != NULL) {

      const cs_int_t  n_c_domains = halo->n_c_domains;
      const cs_int_t  n_transforms = mesh->n_transforms;

      bft_printf("\n\nHalo's data in case of periodicity:\n");
      bft_printf("n_transforms:                %d\n",mesh->n_transforms);

      bft_printf("\nData in the standard halo\n");
      for (i = 0; i < n_transforms; i++)
        for (j = 0; j < n_c_domains; j++)
          bft_printf("< rank:%3d >< transform:%2d > start_idx: %5d"
                     "        n_elts: %5d\n",
                     halo->c_domain_rank[j], i,
                     halo->perio_lst[4*n_c_domains*i + 4*j],
                     halo->perio_lst[4*n_c_domains*i + 4*j+1]);

      bft_printf("\nData in the extended halo\n");
      for (i = 0; i < n_transforms; i++)
        for (j = 0; j < n_c_domains; j++)
          bft_printf("< rank:%3d >< transform:%2d >        "
                     "start_idx:  %5d, n_elts:  %5d\n",
                     halo->c_domain_rank[j], i,
                     halo->perio_lst[4*n_c_domains*i + 4*j+2],
                     halo->perio_lst[4*n_c_domains*i + 4*j+3]);

    } /* End if n_perio > 0 */

  } /* End if mesh->halo != NULL */

  if (mesh->cell_cells_idx != NULL) {

    bft_printf("\n\nCell -> cells connectivity for extended neighborhood\n\n");
    for (i = 0; i < mesh->n_cells; i++) {
      bft_printf("< cell num:%3d>        ", i+1);
      for (j = mesh->cell_cells_idx[i]-1; j < mesh->cell_cells_idx[i+1]-1; j++)
        bft_printf("%d        ", mesh->cell_cells_lst[j]);
      bft_printf("\n");
    }

  }

  /* Dump numbering info */

  cs_numbering_dump(mesh->i_face_numbering);
  cs_numbering_dump(mesh->b_face_numbering);

  /* Modification flag */

  bft_printf("\nModification flag:\n");
  bft_printf("modified:         %d\n",mesh->modified);

  bft_printf("\n\nEND OF DUMP OF MESH STRUCTURE\n\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
