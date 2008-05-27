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
 * Functions handling with periodicity.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_periodicity.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_halo.h"
#include "cs_mesh_quantities.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_perio.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update array of element variable values for periodicity in serial mode,
 * prior to geometric transformations.
 *
 * This function exits immediately in parallel mode, as values should in
 * this case have been updated through halo synchronization.
 *
 * parameters:
 *   halo      --> pointer to halo structure
 *   halo_type --> halo synchronization on standard or extended cells
 *   var       <-> pointer to cell variable
 *----------------------------------------------------------------------------*/

static void
_sync_loc_var(const cs_halo_t  *halo,
              cs_halo_type_t    halo_type,
              cs_real_t         var[])
{
  cs_int_t  i, start, length;

  if (cs_glob_base_nbr == 1) {

    cs_real_t *recv_var
      = var + halo->n_local_elts + halo->index[0];

    start = halo->send_index[0];

    if (halo_type == CS_HALO_EXTENDED)
      length =  halo->send_index[2] - halo->send_index[0];
    else
      length =  halo->send_index[1] - halo->send_index[0];

    for (i = 0; i < length; i++)
      recv_var[i] = var[halo->send_list[start + i]];

  }
}

/*----------------------------------------------------------------------------
 * Update array of strided element variable values for same-rank periodicity,
 * prior to geometric transformations.
 *
 * Periodic values obtained from different ranks should already have been
 * updated through parallel synchronization.
 *
 * parameters:
 *   halo      --> pointer to halo structure
 *   halo_type --> halo synchronization on standard or extended cells
 *   var       <-> pointer to cell variable
 *   stride    --> variable stride
 *----------------------------------------------------------------------------*/

static void
_sync_loc_var_strided(const cs_halo_t  *halo,
                      cs_halo_type_t    halo_type,
                      cs_real_t         var[],
                      int               stride)
{
  cs_int_t  i, j, start, length;

  int rank_id = 0;
  const int  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;

  for (rank_id = 0;
       (   rank_id < halo->n_c_domains
        && halo->c_domain_rank[rank_id] != local_rank);
       rank_id++);

  if (halo->c_domain_rank[rank_id] == local_rank) {

    cs_real_t *recv_var
      = var + halo->n_local_elts + halo->index[2*rank_id];

    start = halo->send_index[2*rank_id];

    if (halo_type == CS_HALO_EXTENDED)
      length =  halo->send_index[2*rank_id + 2] - halo->send_index[2*rank_id];
    else
      length =  halo->send_index[2*rank_id + 1] - halo->send_index[2*rank_id];

    for (i = 0; i < length; i++) {
      for (j = 0; j < stride; j++)
        recv_var[i*stride + j] = var[(halo->send_list[start + i])*stride + j];
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a given
 * vector.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   in_cell_id   --> cell_id of the parent cell.
 *   out_cell_id  --> cell_id of the generated cell.
 *   xyz          <-> array of coordinates
 *----------------------------------------------------------------------------*/

static void
_apply_vector_transfo(cs_real_t    matrix[3][4],
                      cs_int_t     in_cell_id,
                      cs_int_t     out_cell_id,
                      cs_real_t   *xyz)
{
  cs_int_t  i, j;

  cs_real_t  xyz_a[3 + 1];
  cs_real_t  xyz_b[3];

  /* Define the cell center in homogeneous coordinates before
     transformation */

  for (j = 0; j < 3; j++)
    xyz_a[j] = xyz[in_cell_id*3 + j];
  xyz_a[3] = 1;

  /* Initialize output */

  for (i = 0; i < 3; i++)
    xyz_b[i] = 0.;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      xyz_b[i] += matrix[i][j]*xyz_a[j];

  /* Store updated cell center */

  for (j = 0; j < 3; j++)
    xyz[out_cell_id*3 + j] = xyz_b[j];

}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a rotation to a given vector.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   x_in         --> X coord. of the incoming vector
 *   y_in         --> Y coord. of the incoming vector
 *   z_in         --> Z coord. of the incoming vector
 *   x_out        <-- pointer to the X coord. of the output
 *   y_out        <-- pointer to the Y coord. of the output
 *   z_out        <-- pointer to the Z coord. of the output
 *----------------------------------------------------------------------------*/

static void
_apply_vector_rotation(cs_real_t   matrix[3][4],
                       cs_real_t   x_in,
                       cs_real_t   y_in,
                       cs_real_t   z_in,
                       cs_real_t   *x_out,
                       cs_real_t   *y_out,
                       cs_real_t   *z_out)
{
  *x_out = matrix[0][0] * x_in + matrix[0][1] * y_in + matrix[0][2] * z_in;
  *y_out = matrix[1][0] * x_in + matrix[1][1] * y_in + matrix[1][2] * z_in;
  *z_out = matrix[2][0] * x_in + matrix[2][1] * y_in + matrix[2][2] * z_in;
}

/*----------------------------------------------------------------------------
 * Compute a matrix * tensor * Tmatrix product to apply a rotation to a
 * given tensor
 *
 * parameters:
 *   matrix[3][4]        --> transformation matric in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   in11, in12, in13    --> incoming first line of the tensor
 *   in21, in22, in23    --> incoming second line of the tensor
 *   in31, in32, in33    --> incoming third line of the tensor
 *   out11, out12, out13 <-- pointer to the first line of the output
 *   out21, out22, out23 <-- pointer to the second line of the output
 *   out31, out32, out33 <-- pointer to the third line of the output
 *----------------------------------------------------------------------------*/

static void
_apply_tensor_rotation(cs_real_t   matrix[3][4],
                       cs_real_t   in11,
                       cs_real_t   in12,
                       cs_real_t   in13,
                       cs_real_t   in21,
                       cs_real_t   in22,
                       cs_real_t   in23,
                       cs_real_t   in31,
                       cs_real_t   in32,
                       cs_real_t   in33,
                       cs_real_t   *out11,
                       cs_real_t   *out12,
                       cs_real_t   *out13,
                       cs_real_t   *out21,
                       cs_real_t   *out22,
                       cs_real_t   *out23,
                       cs_real_t   *out31,
                       cs_real_t   *out32,
                       cs_real_t   *out33)
{
  cs_int_t  i, j, k;
  cs_real_t  tensorA[3][3], tensorB[3][3];

  tensorA[0][0] = matrix[0][0]*in11 + matrix[0][1]*in12 + matrix[0][2]*in13;
  tensorA[0][1] = matrix[1][0]*in11 + matrix[1][1]*in12 + matrix[1][2]*in13;
  tensorA[0][2] = matrix[2][0]*in11 + matrix[2][1]*in12 + matrix[2][2]*in13;

  tensorA[1][0] = matrix[0][0]*in21 + matrix[0][1]*in22 + matrix[0][2]*in23;
  tensorA[1][1] = matrix[1][0]*in21 + matrix[1][1]*in22 + matrix[1][2]*in23;
  tensorA[1][2] = matrix[2][0]*in21 + matrix[2][1]*in22 + matrix[2][2]*in23;

  tensorA[2][0] = matrix[0][0]*in31 + matrix[0][1]*in32 + matrix[0][2]*in33;
  tensorA[2][1] = matrix[1][0]*in31 + matrix[1][1]*in32 + matrix[1][2]*in33;
  tensorA[2][2] = matrix[2][0]*in31 + matrix[2][1]*in32 + matrix[2][2]*in33;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tensorB[i][j] = 0.;
      for (k = 0; k < 3; k++)
        tensorB[i][j] += matrix[i][k] * tensorA[k][j];
    }
  }

  *out11 = tensorB[0][0];
  *out22 = tensorB[1][1];
  *out33 = tensorB[2][2];

  if (out12 != NULL) {
    *out12 = tensorB[0][1];
    *out13 = tensorB[0][2];
    *out21 = tensorB[1][0];
    *out23 = tensorB[1][2];
    *out31 = tensorB[2][0];
    *out32 = tensorB[2][1];
  }

}

/*----------------------------------------------------------------------------
 * Test if a halo seems compatible with the main mesh's periodic
 * transformations.
 *
 * If a halo is not compatible, abort with an error message.
 *
 * parameters:
 *   halo --> pointer to halo structure
 *----------------------------------------------------------------------------*/

static void
_test_halo_compatibility(const cs_halo_t  *halo)
{
  assert(halo != NULL);

  if (cs_glob_mesh->n_transforms != halo->n_transforms)
    bft_error(__FILE__, __LINE__, 0,
              _("Les %d transformations périodiques du halo ne correspondent\n"
                "pas à celles du maillage principal (au nombre de %d)\n"),
              halo->n_transforms, (int)(cs_glob_mesh->n_transforms));
}

/*----------------------------------------------------------------------------
 * Update dudxyz and wdudxy for periodic ghost cells.
 *
 * Called by PERMAS.
 *
 * parameters:
 *   h_cell_id --> cell id in halo
 *   cell_id   --> cell id
 *   rom       --> density array
 *   call_id   --> first or second call
 *   phase_id  --> phase id
 *   dudxyz    <-> gradient on the components of the velocity.
 *   wdudxy    <-> associated working array.
 *----------------------------------------------------------------------------*/

static void
_update_dudxyz(cs_int_t          h_cell_id,
               cs_int_t          cell_id,
               const cs_real_t  *rom,
               cs_int_t          call_id,
               cs_int_t          phase_id,
               cs_real_t        *dudxyz,
               cs_real_t        *wdudxy)
{
  cs_int_t  i, j, id;

  cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  stride = n_ghost_cells * 3;

  if (call_id == 1) { /* First call */

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        id = h_cell_id + n_ghost_cells*i + stride*j + 3*stride*phase_id;
        wdudxy[id] = dudxyz[id];
        dudxyz[id] *= rom[cell_id];
      }
    }

  }
  else if (call_id == 2) { /* Second call */

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        id = h_cell_id + n_ghost_cells*i + stride*j + 3*stride*phase_id;
        dudxyz[id] = wdudxy[id];
      }
    }

  } /* End if second call */

}

/*----------------------------------------------------------------------------
 * Update drdxyz and wdrdxy for periodic ghost cells.
 *
 * Called by PERMAS.
 *
 * parameters:
 *   h_cell_id --> cell id in halo
 *   cell_id   --> cell id
 *   rom       --> density array
 *   call_id   --> first or second call
 *   phase_id  --> phase id
 *   drdxyz    <-> Gradient on components of Rij (Reynolds stress tensor)
 *   wdrdxy    <-> associated working array.
 *----------------------------------------------------------------------------*/

static void
_update_drdxyz(cs_int_t          h_cell_id,
               cs_int_t          cell_id,
               const cs_real_t  *rom,
               cs_int_t          call_id,
               cs_int_t          phase_id,
               cs_real_t        *drdxyz,
               cs_real_t        *wdrdxy)
{
  cs_int_t  i, j, id;

  cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  stride = 2*3*n_ghost_cells;

  if (call_id == 1) { /* First call */

    for (i = 0; i < 2*3; i++) {
      for (j = 0; j < 3; j++) {
        id = h_cell_id + n_ghost_cells*i + stride*j +3*stride*phase_id;
        wdrdxy[id] = drdxyz[id];
        drdxyz[id] *= rom[cell_id];
      }
    }

  }
  else if (call_id == 2) { /* Second call */

    for (i = 0; i < 2*3; i++) {
      for (j = 0; j < 3; j++) {
        id = h_cell_id + n_ghost_cells*i + stride*j +3*stride*phase_id;
        drdxyz[id] = wdrdxy[id];
      }
    }

  } /* End if second call */

}

/*----------------------------------------------------------------------------
 * Test if there is/are rotation(s).
 *
 * returns:
 *   CS_TRUE or CS_FALSE;
 *----------------------------------------------------------------------------*/

static cs_bool_t
_test_rota_presence(void)
{
  cs_int_t  t_id;

  cs_bool_t  result = CS_FALSE;
  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  for (t_id = 0; t_id < n_transforms; t_id++) {

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {
      result = CS_TRUE;
      break;
    }

  }

  return result;
}

/*----------------------------------------------------------------------------
 * Exchange buffers for PERINR or PERINU
 *
 * parameters:
 *   strid_c    --> stride on the component
 *   strid_v    --> stride on the variable
 *   strid_p    --> stride on the phase
 *   dxyz       <-> gradient on the variable (dudxy or drdxy)
 *   w1, w2, w3 <-> working buffers
 *----------------------------------------------------------------------------*/

static void
_peinur1(cs_int_t      strid_c,
         cs_int_t      strid_v,
         cs_int_t      strid_p,
         cs_real_t    *dxyz,
         cs_real_t    *w1,
         cs_real_t    *w2,
         cs_real_t    *w3)
{
  cs_int_t  i, t_id, rank_id, shift;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_transforms = mesh->n_transforms;

  _sync_loc_var(mesh->halo, mesh->halo_type, w1);
  _sync_loc_var(mesh->halo, mesh->halo_type, w2);
  _sync_loc_var(mesh->halo, mesh->halo_type, w3);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start_std = halo->perio_lst[shift + 4*rank_id];
      length = halo->perio_lst[shift + 4*rank_id + 1];
      end_std = start_std + length;

      if (mesh->halo_type == CS_HALO_EXTENDED) {

        start_ext = halo->perio_lst[shift + 4*rank_id + 2];
        length = halo->perio_lst[shift + 4*rank_id + 3];
        end_ext = start_ext + length;

      }

      for (i = start_std; i < end_std; i++) {
        dxyz[i + strid_c + strid_v*0 + strid_p] = w1[n_cells + i];
        dxyz[i + strid_c + strid_v*1 + strid_p] = w2[n_cells + i];
        dxyz[i + strid_c + strid_v*2 + strid_p] = w3[n_cells + i];
      }

      if (mesh->halo_type == CS_HALO_EXTENDED) {

        for (i = start_ext; i < end_ext; i++) {
          dxyz[i + strid_c + strid_v*0 + strid_p] = w1[n_cells + i];
          dxyz[i + strid_c + strid_v*1 + strid_p] = w2[n_cells + i];
          dxyz[i + strid_c + strid_v*2 + strid_p] = w3[n_cells + i];
        }

      } /* End if extended halo */

    } /* End of loop on ranks */

  } /* End of loop on transformations */

}

/*----------------------------------------------------------------------------
 * Build indexes on periodicities and ranks.
 *
 * parameters:
 *   p_perio_index       <->  pointer to the periodic index
 *   p_send_buffer       <->  pointer to
 *   p_recv_buffer       <->  pointer to
 *   p_vtx_rank_counter  <->  pointer to
 *----------------------------------------------------------------------------*/

static void
_define_periodic_buffers(const cs_mesh_t  *const mesh,
                         cs_int_t         *p_perio_index[],
                         fvm_gnum_t       *p_send_buffer[],
                         fvm_gnum_t       *p_recv_buffer[],
                         cs_int_t         *p_vtx_rank_counter[])
{
  cs_int_t  i, j, k, per_id, rank_id, fac_id, fac_id1, fac_id2, fac_num;
  cs_int_t  distant_rank, shift, pos_count, neg_count, pos_shift, neg_shift;
  cs_int_t  n_f1_vertices, n_f2_vertices;

  cs_int_t  n_vertices_to_send = 0;
  cs_int_t  n_max_face_vertices = 0, n_face_vertices = 0;
  cs_int_t  *vtx_rank_counter = NULL, *counter = NULL;
  cs_int_t  *perio_index = NULL;
  fvm_gnum_t  *send_buffer = NULL, *recv_buffer = NULL;
  fvm_gnum_t  *f1_vertices = NULL, *f2_vertices = NULL;

  const cs_int_t  n_ranks = mesh->n_domains;
  const cs_int_t  n_perio = mesh->n_init_perio;
  const cs_int_t  perio_size = n_ranks * n_perio;
  const cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;
  const cs_int_t  *per_face_idx = cs_glob_mesh_builder->per_face_idx;
  const cs_int_t  *per_face_lst = cs_glob_mesh_builder->per_face_lst;
  const cs_int_t  *per_rank_lst = cs_glob_mesh_builder->per_rank_lst;
  const cs_int_t  *face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_int_t  *face_vtx_lst = mesh->i_face_vtx_lst;
  const fvm_gnum_t  *g_vtx_num = mesh->global_vtx_num;

  /* Allocate and initialize buffers */

  BFT_MALLOC(vtx_rank_counter, n_ranks, cs_int_t);
  BFT_MALLOC(counter, 2*n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    vtx_rank_counter[i] = 0;

  BFT_MALLOC(perio_index, 2*perio_size+1, cs_int_t);

  for (i = 0; i < 2*perio_size+1; i++)
    perio_index[i] = 0;

  /* First step: count and create index */

  for (per_id = 0; per_id < n_perio; per_id++) {

    /*
      The buffer per_face_lst is composed by pairs of face.
      The first element is the local number of the local face. The sign of the
      face indicates the sign of the periodicty.
      The second element is the local number of the distant face.
      The buffer per_rank_lst owns for each distant face the rank id of the
      distant rank associated to the face.
    */

    for (i = per_face_idx[per_id]; i < per_face_idx[per_id + 1]; i++) {

      fac_id = CS_ABS(per_face_lst[2*i]) - 1;
      n_face_vertices = face_vtx_idx[fac_id+1] - face_vtx_idx[fac_id];
      n_max_face_vertices = CS_MAX(n_max_face_vertices, n_face_vertices);

      if (n_ranks > 1)
        distant_rank = per_rank_lst[i] - 1;
      else
        distant_rank = 0;

      shift = 2*n_perio*distant_rank + 2*per_id;

      if (local_rank == distant_rank) {

        vtx_rank_counter[distant_rank] += 2*n_face_vertices;
        perio_index[shift + 1] += n_face_vertices;
        perio_index[shift + 2] += n_face_vertices;

      }
      else {

        vtx_rank_counter[distant_rank] += n_face_vertices;

        if (per_face_lst[2*i] > 0)
          perio_index[shift + 1] += n_face_vertices;
        else
          perio_index[shift + 2] += n_face_vertices;

      }

    } /* End of loop on couples */

  } /* End of loop on periodicities */

  /* Define perio_index */

  for (i = 0; i < 2*perio_size; i++)
    perio_index[i+1] += perio_index[i];

  /* Count the sum of vertices to send */

  for (i = 0; i < n_ranks; i++)
    n_vertices_to_send += vtx_rank_counter[i];

  /* Send/Receive buffer allocation */

  BFT_MALLOC(send_buffer, 2*n_vertices_to_send, fvm_gnum_t);
  recv_buffer = &(send_buffer[n_vertices_to_send]);

  BFT_MALLOC(f1_vertices, 2*n_max_face_vertices, fvm_gnum_t);
  f2_vertices = &(f1_vertices[n_max_face_vertices]);

  /* Second step:
     Fill send_buffer (for all ranks) and recv_buffer (only for local rank) */

  for (per_id = 0; per_id < n_perio; per_id++) {

    for (rank_id = 0; rank_id < 2*n_ranks; rank_id++)
      counter[rank_id] = 0;

    for (i = per_face_idx[per_id]; i < per_face_idx[per_id+1]; i++) {

      fac_num = per_face_lst[2*i];
      fac_id1 = CS_ABS(fac_num) - 1;

      if (n_ranks > 1)
        distant_rank = per_rank_lst[i] - 1;
      else
        distant_rank = 0;

      pos_shift = perio_index[2*n_perio*distant_rank + 2*per_id];
      pos_count = counter[2*distant_rank];
      neg_shift = perio_index[2*n_perio*distant_rank + 2*per_id+1];
      neg_count = counter[2*distant_rank+1];

      for (k= 0, j = face_vtx_idx[fac_id1]-1;
           j < face_vtx_idx[fac_id1+1]-1; j++, k++) {

        if (g_vtx_num != NULL)
          f1_vertices[k] = (fvm_gnum_t)g_vtx_num[face_vtx_lst[j]-1];
        else
          f1_vertices[k] = (fvm_gnum_t)face_vtx_lst[j];

      }

      n_f1_vertices = k;

      if (distant_rank == local_rank) {

        /* => Fill buffer to send and receive for the first and
              second face of the couple */

        fac_id2 = CS_ABS(per_face_lst[2*i+1]) - 1;

        for (k= 0, j = face_vtx_idx[fac_id2]-1;
             j < face_vtx_idx[fac_id2+1]-1; j++, k++) {

          if (g_vtx_num != NULL)
            f2_vertices[k] = (fvm_gnum_t)g_vtx_num[face_vtx_lst[j]-1];
          else
            f2_vertices[k] = (fvm_gnum_t)face_vtx_lst[j];

        }

        n_f2_vertices = k;
        assert(n_f1_vertices == n_f2_vertices);

        /* Send buffer */

        if (fac_num > 0) {

          for (k = 0; k < n_f1_vertices; k++)
            send_buffer[pos_shift + pos_count++] = f1_vertices[k];

          for (k = 0; k < n_f2_vertices; k++)
            send_buffer[neg_shift + neg_count++] = f2_vertices[k];

        }
        else { /* fac_num < 0 */

          for (k = 0; k < n_f1_vertices; k++)
            send_buffer[neg_shift + neg_count++] = f1_vertices[k];

          for (k = 0; k < n_f2_vertices; k++)
            send_buffer[pos_shift + pos_count++] = f2_vertices[k];

        }

        /* Receive buffer */

        pos_count = counter[2*distant_rank];
        neg_count = counter[2*distant_rank+1];

        if (fac_num > 0) {

          for (k = 0; k < n_f1_vertices; k++)
            recv_buffer[neg_shift + neg_count++] = f1_vertices[k];

          for (k =0; k < n_f2_vertices; k++)
            recv_buffer[pos_shift + pos_count++] = f2_vertices[k];

        }
        else {

          for (k = 0; k < n_f1_vertices; k++)
            recv_buffer[pos_shift + pos_count++] = f1_vertices[k];

          for (k = 0; k < n_f2_vertices; k++)
            recv_buffer[neg_shift + neg_count++] = f2_vertices[k];

        }

      } /* End of treatment of self-periodic domain */

      else { /* local_rank != distant_rank => only send_buffer is filled */

        if (fac_num > 0)
          for (k = 0; k < n_f1_vertices; k++)
            send_buffer[pos_shift + pos_count++] = f1_vertices[k];
        else
          for (k = 0; k < n_f1_vertices; k++)
            send_buffer[neg_shift + neg_count++] = f1_vertices[k];

      } /* Non self-periodic */

      counter[2*distant_rank] = pos_count;
      counter[2*distant_rank + 1] = neg_count;

    } /* End of loop on couples */

  } /* End of loop on periodicity */

  BFT_FREE(counter);
  BFT_FREE(f1_vertices);

  /* Return pointers */

  *p_send_buffer = send_buffer;
  *p_recv_buffer = recv_buffer;
  *p_perio_index = perio_index;
  *p_vtx_rank_counter = vtx_rank_counter;

}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update values of periodic cells.
 *
 * VARIJ stands for the periodic variable to deal with.
 *
 * Several cases are possible:
 *
 *   IDIMTE = 0  : VAR11 is a scalar.
 *   IDIMTE = 1  : VAR11, VAR22, VAR33 is a vector.
 *   IDIMTE = 2  : VARIJ is a 3*3 matrix.
 *   IDIMTE = 21 : VARIJ is a diagonal 3*3 matrix (VAR11, VAR22, VAR33).
 *
 * Translation is always treated. Several treatment can be done for rotation:
 *
 *   ITENSO = 0  : only copy values of elements generated by rotation
 *   ITENSO = 1  : ignore rotation.
 *   ITENSO = 11 : reset values of elements generated by rotation
 *
 * - Periodicity for a scalar (IDIMTE = 0, ITENSO = 0). We update VAR11
 *   for translation or rotation periodicity.
 * - Periodicity for a scalar (IDIMTE = 0, ITENSO = 1). We update VAR11 only
 *   for translation periodicity.
 * - Periodicity for a scalar (IDIMTE = 0, ITENSO = 11). We update VAR11 only
 *   for translation periodicity. VAR11 is reseted for rotation periodicicty.
 *
 *   We use this option to cancel the halo for rotational periodicities
 *   in iterative solvers when solving for vectors and tensors by
 *   increment. This is an approximative solution, which does not seem
 *   worse than another.
 *
 * - with a vector (IDIMTE = 0, ITENSO = 2), we update
 *   VAR11, VAR22, VAR33, for translation only.
 * - with a vector (IDIMTE = 1, ITENSO = *), we update
 *   VAR11, VAR22, VAR33, for translation and rotation.
 * - with a tensor of rank 2 (IDIMTE = 2, ITENSO = *), we update
 *   VAR11, V12, VAR13, VAR21, VAR22, VAR23, VAR31, VAR32, VAR33,
 *   for translation and rotation.
 * - with a tensor or rank 2 (IDIMTE = 21, ITENSO = *) , we update
 *   VAR11, VAR22, VAR33, for translation and rotation (the tensor
 *   is considered diagonal).
 *
 * Fortran API:
 *
 * SUBROUTINE PERCOM
 * *****************
 *
 * INTEGER          IDIMTE        :  -> : variable dimension (maximum 3)
 *                                        0 : scalar (VAR11), or considered
 *                                            scalar
 *                                        1 : vector (VAR11,VAR22,VAR33)
 *                                        2 : tensor of rank 2 (VARIJ)
 *                                       21 : tensor of rank 2 supposed
 *                                            diagonal (VAR11, VAR22, VAR33)
 * INTEGER          ITENSO        :  -> : to define rotation behavior
 *                                        0 : scalar (VAR11)
 *                                        1 : tensor or vector component
 *                                            (VAR11), implicit in
 *                                            translation case
 *                                       11 : same as ITENSO=1 with vector
 *                                            or tensor component cancelled
 *                                            for rotation
 *                                        2 : vector (VAR11, VAR22, VAR33)
 *                                            implicit for rotation
 * DOUBLE PRECISION VAR11(NCELET) :  -  : component 11 of rank 2 tensor
 * DOUBLE PRECISION VAR12(NCELET) :  -  : component 12 of rank 2 tensor
 * DOUBLE PRECISION VAR13(NCELET) :  -  : component 13 of rank 2 tensor
 * DOUBLE PRECISION VAR21(NCELET) :  -  : component 21 of rank 2 tensor
 * DOUBLE PRECISION VAR22(NCELET) :  -  : component 22 of rank 2 tensor
 * DOUBLE PRECISION VAR23(NCELET) :  -  : component 23 of rank 2 tensor
 * DOUBLE PRECISION VAR31(NCELET) :  -  : component 31 of rank 2 tensor
 * DOUBLE PRECISION VAR32(NCELET) :  -  : component 32 of rank 2 tensor
 * DOUBLE PRECISION VAR33(NCELET) :  -  : component 33 of rank 2 tensor
 *----------------------------------------------------------------------------*/

void
CS_PROCF (percom, PERCOM) (const cs_int_t  *idimte,
                           const cs_int_t  *itenso,
                           cs_real_t        var11[],
                           cs_real_t        var12[],
                           cs_real_t        var13[],
                           cs_real_t        var21[],
                           cs_real_t        var22[],
                           cs_real_t        var23[],
                           cs_real_t        var31[],
                           cs_real_t        var32[],
                           cs_real_t        var33[])
{
  cs_bool_t  bool_err = CS_FALSE;

  const cs_halo_t *halo = cs_glob_mesh->halo;

  /* 1. Checking    */
  /*----------------*/

  if (*idimte != 0 && *idimte != 1 && *idimte != 2 && *idimte != 21)
    bool_err = CS_TRUE;

  if (*itenso != 0 && *itenso != 1 && *itenso != 11 && *itenso != 2)
    bool_err = CS_TRUE;


  if (bool_err == CS_TRUE)
    bft_error(__FILE__, __LINE__, 0,
              _("IDIMTE et/ou ITENSO ont des valeurs incohérentes"));

  /* 2. Synchronization  */
  /*---------------------*/

  if (*idimte == 0) {

    /* Input parameter is a scalar. Sync values on periodic cells for
       translation and rotation */

    if (*itenso == 0)
      cs_perio_sync_var_scal(halo,
                             CS_HALO_STANDARD,
                             CS_PERIO_ROTA_COPY,
                             var11);

    /* Input parameter is a scalar. Sync values on periodic cells for
       translation. We ignore the rotation tranformations. */

    else if (*itenso == 1)
      cs_perio_sync_var_scal(halo,
                             CS_HALO_STANDARD,
                             CS_PERIO_ROTA_IGNORE,
                             var11);

    /* Reset elements generated by a rotation. (used in Jacobi for Reynolds
       stresses or to solve velocity )*/

    else if (*itenso == 11)
      cs_perio_sync_var_scal(halo,
                             CS_HALO_STANDARD,
                             CS_PERIO_ROTA_RESET,
                             var11);

    /* Variable is part of a tensor, so exchange is possible only in
       translation; 3 components are exchanged (we may guess that for
       rotations, something has been done before; see PERINR for example). */

    else if (*itenso == 2)
      cs_perio_sync_var_vect(halo,
                             CS_HALO_STANDARD,
                             CS_PERIO_ROTA_IGNORE,
                             var11, var22, var33);

  } /* End of idimte == 0 case */

  /* --> If we want to handle the variable as a vector, we suppose that
     it is (at least) a vector. Translation and rotation are exchanged. */

  else if (*idimte == 1)
    cs_perio_sync_var_vect(halo,
                           CS_HALO_STANDARD,
                           CS_PERIO_ROTA_COPY,
                           var11, var22, var33);

  /* --> If we want to handle the variable as a tensor, we suppose that
     it is a tensor. Translation and rotation are exchanged. */

  else if (*idimte == 2)
    cs_perio_sync_var_tens(halo,
                           CS_HALO_STANDARD,
                           var11, var12, var13,
                           var21, var22, var23,
                           var31, var32, var33);

  /* --> If we want to handle the variable as a tensor, but that
     it is a tensor's diagonal, we suppose that it is a tensor.
     Translation and rotation are exchanged. */

  else if (*idimte == 21)
    cs_perio_sync_var_diag(halo,
                           CS_HALO_STANDARD,
                           var11, var22, var33);

}

/*----------------------------------------------------------------------------
 * Save or restore rotation terms in halo for a cell variable.
 *
 * Only the terms from one variable may be saved at a time, and they
 * must be restored before other terms may be saved.
 *
 * This function may be used to cancel the effect of halo synchronization
 * for for rotational periodicities, by saving halo rotation terms
 * before halo synchronization, and restoring them after synchronization.
 *
 * Fortran API:
 *
 * SUBROUTINE PERSVR
 * *****************
 *
 * INTEGER          IMODE         :  -> : mode
 *                                        0 : save halo rotation terms
 *                                        1 : restore halo rotation terms
 * DOUBLE PRECISION VAR(NCELET)   :  -  : cell variable
 *----------------------------------------------------------------------------*/

void
CS_PROCF (persvr, PERSVR) (const cs_int_t  *mode,
                           cs_real_t        var[])
{
  size_t      sr_count = 0;

  static size_t save_buffer_size = 0;
  static cs_real_t *save_buffer = NULL;

  const cs_halo_t  *halo = cs_glob_mesh->halo;

  if (*mode == 0) {

    save_buffer_size = halo->n_elts[0];
    BFT_MALLOC(save_buffer, save_buffer_size, cs_real_t);

    sr_count = cs_perio_save_rotation_halo(halo,
                                           CS_HALO_STANDARD,
                                           var,
                                           save_buffer);
  }
  else if (*mode == 1) {

    if (save_buffer_size != (size_t)(halo->n_elts[0]))
      bft_error(__FILE__, __LINE__, 0,
                _("La taille du tableau de sauvegarde du halo de rotation (%d)\n"
                  "ne correspond pas à la taille enregistréépas (%d)."),
                (int)(halo->n_elts[0]), (int)(save_buffer_size));

    sr_count = cs_perio_restore_rotation_halo(halo,
                                              CS_HALO_STANDARD,
                                              var,
                                              save_buffer);

    BFT_FREE(save_buffer);
    save_buffer_size = 0;
  }

  assert(sr_count <= (size_t)(halo->n_elts[0]));
}

/*----------------------------------------------------------------------------
 * Periodicity management for INIMAS
 *
 * If INIMAS is called by NAVSTO :
 *    We assume that gradient on ghost cells given by a rotation is known
 *    and is equal to the velocity one for the previous time step.
 * If INIMAS is called by DIVRIJ
 *    We assume that (more justifiable than in previous case) gradient on
 *    ghost cells given by rotation is equal to Rij gradient for the previous
 *    time step.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PERMAS
 * *****************
 *
 * INTEGER          IMASPE      :  -> : suivant l'appel de INIMAS
 *                                          = 1 si appel de RESOLP ou NAVSTO
 *                                          = 2 si appel de DIVRIJ
 * INTEGER          IPHAS       :  -> : numero de phase courante
 * INTEGER          IMASPE      :  -> : indicateur d'appel dans INIMAS
 *                                          = 1 si appel au debut
 *                                          = 2 si appel a la fin
 * DOUBLE PRECISION ROM(NCELET) :  -> : masse volumique aux cellules
 * DOUBLE PRECISION DUDXYZ      :  -> : gradient de U aux cellules halo pour
 *                                      l'approche explicite en periodicite
 * DOUBLE PRECISION DRDXYZ      :  -> : gradient de R aux cellules halo pour
 *                                      l'approche explicite en periodicite
 * DOUBLE PRECISION WDUDXY      :  -  : tableau de travail pour DUDXYZ
 * DOUBLE PRECISION WDRDXY      :  -  : tableau de travail pour DRDXYZ
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3*NPHAS
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3*NPHAS
 *----------------------------------------------------------------------------*/

void
CS_PROCF (permas, PERMAS)(const cs_int_t    *imaspe,
                          const cs_int_t    *iphas,
                          const cs_int_t    *iappel,
                          cs_real_t          rom[],
                          cs_real_t         *dudxyz,
                          cs_real_t         *drdxyz,
                          cs_real_t         *wdudxy,
                          cs_real_t         *wdrdxy)
{
  cs_int_t  i, cell_id, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;
  cs_halo_type_t  halo_type = mesh->halo_type;

  const cs_int_t  phase_id = *iphas - 1;

  if (halo_type == CS_HALO_N_TYPES)
    return;

  if (*iappel == 1)
    _sync_loc_var(mesh->halo, mesh->halo_type, rom);

  for (t_id = 0; t_id < mesh->n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start_std = halo->perio_lst[shift + 4*rank_id];
      length = halo->perio_lst[shift + 4*rank_id + 1];
      end_std = start_std + length;

      if (halo_type == CS_HALO_EXTENDED) {

        start_ext = halo->perio_lst[shift + 4*rank_id + 2];
        length = halo->perio_lst[shift + 4*rank_id + 3];
        end_ext = start_ext + length;

      }

      for (i = start_std; i < end_std; i++) {

        cell_id = mesh->n_cells + i;

        if (*imaspe == 1)
          _update_dudxyz(i, cell_id, rom, *iappel, phase_id, dudxyz, wdudxy);

        if (*imaspe == 2)
          _update_drdxyz(i, cell_id, rom, *iappel, phase_id, drdxyz, wdrdxy);

      } /* End of loop on halo elements */

      if (halo_type == CS_HALO_EXTENDED) {

        for (i = start_ext; i < end_ext; i++) {

          cell_id = mesh->n_cells + i;

          if (*imaspe == 1)
            _update_dudxyz(i, cell_id, rom, *iappel, phase_id, dudxyz, wdudxy);

          if (*imaspe == 2)
            _update_drdxyz(i, cell_id, rom, *iappel, phase_id, drdxyz, wdrdxy);

        } /* End of loop on halo elements */

      } /* End if extended halo */

    } /* End of loop on ranks */

  } /* End of loop on transformation */

}

/*----------------------------------------------------------------------------
 * Process DPDX, DPDY, DPDZ buffers in case of rotation on velocity vector and
 * Reynolds stress tensor.
 *
 * We retrieve the gradient given by PERINU and PERINR (PHYVAR) for the
 * velocity and the Reynolds stress tensor in a buffer on ghost cells. Then
 * we define DPDX, DPDY and DPDZ gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a non-scalar
 * variable because we have to know the all three components in GRADRC.
 *
 * Otherwise, we can implicitly treat values given by translation. There will
 * be replace further in GRADRC.
 *
 * We define IDIMTE to 0 and ITENSO to 2 for the velocity vector and the
 * Reynolds stress tensor. We will still have to apply translation to these
 * variables so we define a tag to not forget to do it.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables in GRADRC. We define IDIMTE to 1 and ITENSO to 0.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PERING
 * *****************
 *
 * INTEGER          NPHAS        :  -> : numero de phase courante
 * INTEGER          IVAR         :  -> : numero de la variable
 * INTEGER          IDIMTE       : <-  : dimension de la variable (maximum 3)
 *                                        0 : scalaire (VAR11), ou assimile
 *                                            scalaire
 *                                        1 : vecteur (VAR11,VAR22,VAR33)
 *                                        2 : tenseur d'ordre 2 (VARIJ)
 *                                       21 : tenseur d'ordre 2 suppose
 *                                            diagonal (VAR11, VAR22, VAR33)
 * INTEGER          ITENSO       : <-  : pour l'explicitation de la rotation
 *                                        0 : scalaire (VAR11)
 *                                        1 : composante de vecteur ou de
 *                                            tenseur (VAR11) implicite pour
 *                                            la translation
 *                                       11 : reprend le traitement ITENSO=1
 *                                            et composante de vecteur ou de
 *                                            tenseur annulee pour la rotation
 *                                        2 : vecteur (VAR11 et VAR22 et VAR33)
 *                                            implicite pour la rotation
 * INTEGER          IPEROT       :  -> : indicateur du nombre de periodicte de
 *                                       rotation
 * INTEGER          IGUPER       :  -> : 0/1 indique qu'on a /n'a pas calcule
 *                                       les gradients dans DUDXYZ
 * INTEGER          IGRPER       :  -> : 0/1 indique qu'on a /n'a pas calcule
 *                                       les gradients dans DRDXYZ
 * INTEGER          IU           :  -> : position de la Vitesse(x,y,z)
 * INTEGER          IV           :  -> : dans RTP, RTPA
 * INTEGER          IW           :  -> :     "                   "
 * INTEGER          ITYTUR       :  -> : turbulence (Rij-epsilon ITYTUR = 3)
 * INTEGER          IR11         :  -> : position des Tensions de Reynolds
 * INTEGER          IR22         :  -> : en Rij dans RTP, RTPA
 * INTEGER          IR33         :  -> :     "                   "
 * INTEGER          IR12         :  -> :     "                   "
 * INTEGER          IR13         :  -> :     "                   "
 * INTEGER          IR23         :  -> :     "                   "
 * DOUBLE PRECISION DPDX(NCELET) : <-> : gradient de IVAR
 * DOUBLE PRECISION DPDY(NCELET) : <-> :    "        "
 * DOUBLE PRECISION DPDZ(NCELET) : <-> :    "        "
 * DOUBLE PRECISION DUDXYZ       :  -> : gradient de U aux cellules halo pour
 *                                       l'approche explicite en periodicite
 * DOUBLE PRECISION DRDXYZ       :  -> : gradient de R aux cellules halo pour
 *                                       l'approche explicite en periodicite
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3*NPHAS
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3*NPHAS
 *----------------------------------------------------------------------------*/

void
CS_PROCF (pering, PERING)(const cs_int_t    *nphas,
                          const cs_int_t    *ivar,
                          cs_int_t          *idimte,
                          cs_int_t          *itenso,
                          const cs_int_t    *iperot,
                          const cs_int_t    *iguper,
                          const cs_int_t    *igrper,
                          const cs_int_t     iu[CS_NPHSMX],
                          const cs_int_t     iv[CS_NPHSMX],
                          const cs_int_t     iw[CS_NPHSMX],
                          const cs_int_t     itytur[CS_NPHSMX],
                          const cs_int_t     ir11[CS_NPHSMX],
                          const cs_int_t     ir22[CS_NPHSMX],
                          const cs_int_t     ir33[CS_NPHSMX],
                          const cs_int_t     ir12[CS_NPHSMX],
                          const cs_int_t     ir13[CS_NPHSMX],
                          const cs_int_t     ir23[CS_NPHSMX],
                          cs_real_t          dpdx[],
                          cs_real_t          dpdy[],
                          cs_real_t          dpdz[],
                          const cs_real_t   *dudxyz,
                          const cs_real_t   *drdxyz)
{
  cs_int_t  i, rank_id, phase_id, t_id, shift;
  cs_int_t  start_std, end_std, length, start_ext, end_ext, tag;
  cs_int_t  d_ph, d_var;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_cells   = mesh->n_cells;
  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  stride1 = n_ghost_cells * 3;
  const cs_int_t  stride2 = stride1 * 3;

  /* We treat the gradient like a vector by default ...
     (i.e. we assume that this is the gradient of a scalar. */

  *idimte = 1;
  *itenso = 0;

  /*
    When there is periodicity of rotation :
      - Test if variable is a vector or a tensor,
      - Retrieve its gradient values for ghost cells and for the previous
        time step without reconstruction
      - Finally, we define IDIMTE and ITENSO for the next call to PERCOM
        Ghost cells without rotation are reset and other ghost cells keep
        their value.
  */

  if (_test_rota_presence()) {

    assert(*iperot > 0);

    tag = 0;

    for (phase_id = 0; phase_id < *nphas; phase_id++) {

      if (   *ivar == iu[phase_id]
          || *ivar == iv[phase_id]
          || *ivar == iw[phase_id]) {

        d_ph = stride2*phase_id;
        tag = 1;

        if (*ivar == iu[phase_id]) d_var = 0;
        if (*ivar == iv[phase_id]) d_var = n_ghost_cells;
        if (*ivar == iw[phase_id]) d_var = 2*n_ghost_cells;

        if (*iguper == 1) { /* dudxyz not computed */

          for (t_id = 0; t_id < n_transforms; t_id++) {

            shift = 4 * halo->n_c_domains * t_id;

            for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

              start_std = halo->perio_lst[shift + 4*rank_id];
              length = halo->perio_lst[shift + 4*rank_id + 1];
              end_std = start_std + length;

              if (mesh->halo_type == CS_HALO_EXTENDED) {

                start_ext = halo->perio_lst[shift + 4*rank_id + 2];
                length = halo->perio_lst[shift + 4*rank_id + 3];
                end_ext = start_ext + length;

              }

              for (i = start_std; i < end_std; i++) {
                dpdx[n_cells + i] = dudxyz[i + d_var + stride1*0 + d_ph];
                dpdy[n_cells + i] = dudxyz[i + d_var + stride1*1 + d_ph];
                dpdz[n_cells + i] = dudxyz[i + d_var + stride1*2 + d_ph];
              }

              if (mesh->halo_type == CS_HALO_EXTENDED) {

                for (i = start_ext; i < end_ext; i++) {
                  dpdx[n_cells + i] = dudxyz[i + d_var + stride1*0 + d_ph];
                  dpdy[n_cells + i] = dudxyz[i + d_var + stride1*1 + d_ph];
                  dpdz[n_cells + i] = dudxyz[i + d_var + stride1*2 + d_ph];
                }

              } /* End if extended halo */

            } /* End of loop on ranks */

          } /* End of loop on transformations */

        } /* End if *iguper == 1 */

      } /* If *ivar == iu or iv or iw */

      else if ((itytur[phase_id] == 3) &&
               (*ivar == ir11[phase_id] || *ivar == ir22[phase_id] ||
                *ivar == ir33[phase_id] || *ivar == ir12[phase_id] ||
                *ivar == ir13[phase_id] || *ivar == ir23[phase_id])) {

        d_ph = 2*stride2*phase_id;
        tag = 1;

        if (*ivar == ir11[phase_id]) d_var = 0;
        if (*ivar == ir22[phase_id]) d_var = n_ghost_cells;
        if (*ivar == ir33[phase_id]) d_var = 2*n_ghost_cells;
        if (*ivar == ir12[phase_id]) d_var = 3*n_ghost_cells;
        if (*ivar == ir13[phase_id]) d_var = 4*n_ghost_cells;
        if (*ivar == ir23[phase_id]) d_var = 5*n_ghost_cells;

        if (*igrper == 1) {

          for (t_id = 0; t_id < n_transforms; t_id++) {

            shift = 4 * halo->n_c_domains * t_id;

            for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

              start_std = halo->perio_lst[shift + 4*rank_id];
              length = halo->perio_lst[shift + 4*rank_id + 1];
              end_std = start_std + length;

              if (mesh->halo_type == CS_HALO_EXTENDED) {

                start_ext = halo->perio_lst[shift + 4*rank_id + 2];
                length = halo->perio_lst[shift + 4*rank_id + 3];
                end_ext = start_ext + length;

              }

              for (i = start_std; i < end_std; i++) {
                dpdx[n_cells + i] = drdxyz[i + d_var + 2*stride1*0 + d_ph];
                dpdy[n_cells + i] = drdxyz[i + d_var + 2*stride1*1 + d_ph];
                dpdz[n_cells + i] = drdxyz[i + d_var + 2*stride1*2 + d_ph];
              }

              if (mesh->halo_type == CS_HALO_EXTENDED) {

                for (i = start_ext; i < end_ext; i++) {
                  dpdx[n_cells + i] = drdxyz[i + d_var + 2*stride1*0 + d_ph];
                  dpdy[n_cells + i] = drdxyz[i + d_var + 2*stride1*1 + d_ph];
                  dpdz[n_cells + i] = drdxyz[i + d_var + 2*stride1*2 + d_ph];
                }

              } /* End if extended halo */

            } /* End of loop on ranks */

          } /* End of loop on transformations */

        } /* End if *igrper == 1 */

      } /* If itytur[phase_id] == 3 and *ivar == irij */

    } /* End of loop on phases */

    if (tag == 1) {
      *idimte = 0;
      *itenso = 2;
    }

  } /* If there is/are rotation(s) */

}

/*----------------------------------------------------------------------------
 * Exchange buffers for PERINU
 *
 * Fortran Interface:
 *
 * SUBROUTINE PEINU1
 * *****************
 *
 * INTEGER          ISOU          :  -> : component of the velocity vector
 * INTEGER          IPHAS         :  -> : current phase number
 * DOUBLE PRECISION DUDXYZ        : <-> : gradient of the velocity vector
 *                                        for ghost cells and for an explicit
 *                                        treatment of the periodicity.
 * DOUBLE PRECISION W1..3(NCELET) :  -  : working buffers
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3*NPHAS
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinu1, PEINU1)(const cs_int_t    *isou,
                          const cs_int_t    *iphas,
                          cs_real_t         *dudxyz,
                          cs_real_t          w1[],
                          cs_real_t          w2[],
                          cs_real_t          w3[])
{
  cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  comp_id = *isou - 1;
  const cs_int_t  phas_id = *iphas - 1;
  const cs_int_t  strid_v = n_ghost_cells * 3;
  const cs_int_t  strid_p = strid_v * 3 * phas_id;
  const cs_int_t  strid_c = n_ghost_cells * comp_id;

  _peinur1(strid_c, strid_v, strid_p, dudxyz, w1, w2, w3);

}

/*----------------------------------------------------------------------------
 * Apply rotation on DUDXYZ tensor.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PEINU2 (VAR)
 * *****************
 *
 * INTEGER          IPHAS         :  -> : current phase number
 * DOUBLE PRECISION DUDXYZ        : <-> : gradient of the velocity vector
 *                                        for ghost cells and for an explicit
 *                                        treatment of the periodicity.
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3*NPHAS
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinu2, PEINU2)(const cs_int_t    *iphas,
                          cs_real_t         *dudxyz)
{
  cs_int_t  i, t_id, rank_id, shift;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;
  fvm_periodicity_type_t  perio_type;
  cs_real_t  matrix[3][4];

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  phase_id = *iphas - 1;
  const cs_int_t  stride = 3 * n_ghost_cells;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  /* Macro pour la position au sein d'un tableau */

#define GET_ID1(i, j, k) ( \
   i + n_ghost_cells*j + stride*k + 3*stride*phase_id)

  if (mesh->halo_type == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  /* Compute the new cell centers through periodicity */

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        for (i = start_std; i < end_std; i++)
          _apply_tensor_rotation(matrix,
                                 dudxyz[GET_ID1(i,0,0)],
                                 dudxyz[GET_ID1(i,0,1)],
                                 dudxyz[GET_ID1(i,0,2)],
                                 dudxyz[GET_ID1(i,1,0)],
                                 dudxyz[GET_ID1(i,1,1)],
                                 dudxyz[GET_ID1(i,1,2)],
                                 dudxyz[GET_ID1(i,2,0)],
                                 dudxyz[GET_ID1(i,2,1)],
                                 dudxyz[GET_ID1(i,2,2)],
                                 &dudxyz[GET_ID1(i,0,0)],
                                 &dudxyz[GET_ID1(i,0,1)],
                                 &dudxyz[GET_ID1(i,0,2)],
                                 &dudxyz[GET_ID1(i,1,0)],
                                 &dudxyz[GET_ID1(i,1,1)],
                                 &dudxyz[GET_ID1(i,1,2)],
                                 &dudxyz[GET_ID1(i,2,0)],
                                 &dudxyz[GET_ID1(i,2,1)],
                                 &dudxyz[GET_ID1(i,2,2)]);

        if (mesh->halo_type == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

          for (i = start_ext; i < end_ext; i++)
            _apply_tensor_rotation(matrix,
                                   dudxyz[GET_ID1(i,0,0)],
                                   dudxyz[GET_ID1(i,0,1)],
                                   dudxyz[GET_ID1(i,0,2)],
                                   dudxyz[GET_ID1(i,1,0)],
                                   dudxyz[GET_ID1(i,1,1)],
                                   dudxyz[GET_ID1(i,1,2)],
                                   dudxyz[GET_ID1(i,2,0)],
                                   dudxyz[GET_ID1(i,2,1)],
                                   dudxyz[GET_ID1(i,2,2)],
                                   &dudxyz[GET_ID1(i,0,0)],
                                   &dudxyz[GET_ID1(i,0,1)],
                                   &dudxyz[GET_ID1(i,0,2)],
                                   &dudxyz[GET_ID1(i,1,0)],
                                   &dudxyz[GET_ID1(i,1,1)],
                                   &dudxyz[GET_ID1(i,1,2)],
                                   &dudxyz[GET_ID1(i,2,0)],
                                   &dudxyz[GET_ID1(i,2,1)],
                                   &dudxyz[GET_ID1(i,2,2)]);

        } /* End if extended halo exists */

      } /* End of loop on ranks */

    } /* End if periodicity is a rotation */

  } /* End of loop on transformation */

}

/*----------------------------------------------------------------------------
 * Exchange buffers for PERINR
 *
 * Fortran Interface
 *
 * SUBROUTINE PEINR1 (VAR)
 * *****************
 *
 * INTEGER          ISOU          : -> : component of the Reynolds stress tensor
 * INTEGER          IPHAS         : -> : current phase number
 * DOUBLE PRECISION DRDXYZ        : -> : gradient of the Reynolds stress tensor
 *                                       for ghost cells and for an explicit
 *                                       treatment of the periodicity.
 * DOUBLE PRECISION W1..3(NCELET) : -  : working buffers
 *
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3*NPHAS
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinr1, PEINR1)(const cs_int_t    *isou,
                          const cs_int_t    *iphas,
                          cs_real_t         *drdxyz,
                          cs_real_t          w1[],
                          cs_real_t          w2[],
                          cs_real_t          w3[])
{
  cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  comp_id = *isou - 1;
  const cs_int_t  phase_id = *iphas - 1;
  const cs_int_t  strid_v = 2 * n_ghost_cells * 3;
  const cs_int_t  strid_p = strid_v * 3 * phase_id;
  const cs_int_t  strid_c = n_ghost_cells * comp_id;

  _peinur1(strid_c, strid_v, strid_p, drdxyz, w1, w2, w3);

}

/*----------------------------------------------------------------------------
 * Apply rotation on the gradient of Reynolds stress tensor
 *
 * Fortran Interface:
 *
 * SUBROUTINE PEINR2 (VAR)
 * *****************
 *
 * INTEGER          IPHAS         :  -> : current phase number
 * DOUBLE PRECISION DRDXYZ        :  -> : gradient of the Reynolds stress tensor
 *                                       for ghost cells and for an explicit
 *                                       treatment of the periodicity.
 *
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3*NPHAS
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinr2, PEINR2)(const cs_int_t    *iphas,
                          cs_real_t         *drdxyz)
{
  cs_int_t  i, i1, i2, j, k, l, m, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;
  fvm_periodicity_type_t  perio_type;

  cs_real_t  matrix[3][4];
  cs_real_t  tensa[3][3][3];
  cs_real_t  tensb[3][3][3];
  cs_real_t  tensc[3][3][3];

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  n_ghost_cells = mesh->n_ghost_cells;
  const cs_int_t  stride = 2 * 3 * n_ghost_cells;
  const cs_int_t  phase_id = *iphas - 1;
  const cs_halo_type_t  halo_type = mesh->halo_type;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

#define GET_ID2(elt_id, i, j)   ( \
 elt_id + n_ghost_cells*i + stride*j + 3*stride*phase_id)

  if (halo_type == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        for (i = start_std; i < end_std; i++) {

          tensa[0][0][0] = drdxyz[GET_ID2(i,0,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,0,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,0,2)] * matrix[0][2];
          tensa[0][1][0] = drdxyz[GET_ID2(i,3,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,3,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,3,2)] * matrix[0][2];
          tensa[0][2][0] = drdxyz[GET_ID2(i,4,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,4,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,4,2)] * matrix[0][2];

          tensa[0][0][1] = drdxyz[GET_ID2(i,0,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,0,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,0,2)] * matrix[1][2];
          tensa[0][1][1] = drdxyz[GET_ID2(i,3,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,3,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,3,2)] * matrix[1][2];
          tensa[0][2][1] = drdxyz[GET_ID2(i,4,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,4,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,4,2)] * matrix[1][2];

          tensa[0][0][2] = drdxyz[GET_ID2(i,0,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,0,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,0,2)] * matrix[2][2];
          tensa[0][1][2] = drdxyz[GET_ID2(i,3,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,3,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,3,2)] * matrix[2][2];
          tensa[0][2][2] = drdxyz[GET_ID2(i,4,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,4,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,4,2)] * matrix[2][2];

          tensa[1][0][0] = drdxyz[GET_ID2(i,3,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,3,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,3,2)] * matrix[0][2];
          tensa[1][1][0] = drdxyz[GET_ID2(i,1,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,1,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,1,2)] * matrix[0][2];
          tensa[1][2][0] = drdxyz[GET_ID2(i,5,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,5,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,5,2)] * matrix[0][2];

          tensa[1][0][1] = drdxyz[GET_ID2(i,3,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,3,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,3,2)] * matrix[1][2];
          tensa[1][1][1] = drdxyz[GET_ID2(i,1,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,1,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,1,2)] * matrix[1][2];
          tensa[1][2][1] = drdxyz[GET_ID2(i,5,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,5,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,5,2)] * matrix[1][2];

          tensa[1][0][2] = drdxyz[GET_ID2(i,3,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,3,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,3,2)] * matrix[2][2];
          tensa[1][1][2] = drdxyz[GET_ID2(i,1,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,1,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,1,2)] * matrix[2][2];
          tensa[1][2][2] = drdxyz[GET_ID2(i,5,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,5,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,5,2)] * matrix[2][2];

          tensa[2][0][0] = drdxyz[GET_ID2(i,4,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,4,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,4,2)] * matrix[0][2];
          tensa[2][1][0] = drdxyz[GET_ID2(i,5,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,5,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,5,2)] * matrix[0][2];
          tensa[2][2][0] = drdxyz[GET_ID2(i,2,0)] * matrix[0][0] +
                           drdxyz[GET_ID2(i,2,1)] * matrix[0][1] +
                           drdxyz[GET_ID2(i,2,2)] * matrix[0][2];

          tensa[2][0][1] = drdxyz[GET_ID2(i,4,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,4,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,4,2)] * matrix[1][2];
          tensa[2][1][1] = drdxyz[GET_ID2(i,5,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,5,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,5,2)] * matrix[1][2];
          tensa[2][2][1] = drdxyz[GET_ID2(i,2,0)] * matrix[1][0] +
                           drdxyz[GET_ID2(i,2,1)] * matrix[1][1] +
                           drdxyz[GET_ID2(i,2,2)] * matrix[1][2];

          tensa[2][0][2] = drdxyz[GET_ID2(i,4,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,4,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,4,2)] * matrix[2][2];
          tensa[2][1][2] = drdxyz[GET_ID2(i,5,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,5,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,5,2)] * matrix[2][2];
          tensa[2][2][2] = drdxyz[GET_ID2(i,2,0)] * matrix[2][0] +
                           drdxyz[GET_ID2(i,2,1)] * matrix[2][1] +
                           drdxyz[GET_ID2(i,2,2)] * matrix[2][2];

          for (j = 0; j < 3; j++)
            for (l = 0; l < 3; l++)
              for (k = 0; k < 3; k++) {
                tensb[j][l][k] = 0.;
                for (m = 0; m < 3; m++)
                  tensb[j][l][k] += matrix[j][m] * tensa[l][m][k];
              }

          for (k = 0; k < 3; k++)
            for (i1 = 0; i1 < 3; i1++)
              for (i2 = 0; i2 < 3; i2++) {
                tensc[k][i1][i2] = 0.;
                for (l = 0; l < 3; l++)
                  tensc[k][i1][i2] += matrix[k][l] * tensb[i1][l][i2];
              }


          drdxyz[GET_ID2(i,0,0)] = tensc[0][0][0];
          drdxyz[GET_ID2(i,0,1)] = tensc[0][0][1];
          drdxyz[GET_ID2(i,0,2)] = tensc[0][0][2];

          drdxyz[GET_ID2(i,3,0)] = tensc[0][1][0];
          drdxyz[GET_ID2(i,3,1)] = tensc[0][1][1];
          drdxyz[GET_ID2(i,3,2)] = tensc[0][1][2];

          drdxyz[GET_ID2(i,4,0)] = tensc[0][2][0];
          drdxyz[GET_ID2(i,4,1)] = tensc[0][2][1];
          drdxyz[GET_ID2(i,4,2)] = tensc[0][2][2];

          drdxyz[GET_ID2(i,1,0)] = tensc[1][1][0];
          drdxyz[GET_ID2(i,1,1)] = tensc[1][1][1];
          drdxyz[GET_ID2(i,1,2)] = tensc[1][1][2];

          drdxyz[GET_ID2(i,5,0)] = tensc[1][2][0];
          drdxyz[GET_ID2(i,5,1)] = tensc[1][2][1];
          drdxyz[GET_ID2(i,5,2)] = tensc[1][2][2];

          drdxyz[GET_ID2(i,2,0)] = tensc[2][2][0];
          drdxyz[GET_ID2(i,2,1)] = tensc[2][2][1];
          drdxyz[GET_ID2(i,2,2)] = tensc[2][2][2];

        } /* End of loop on standard ghost cells */

        if (halo_type == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

          for (i = start_ext; i < end_ext; i++) {

            tensa[0][0][0] = drdxyz[GET_ID2(i,0,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,0,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,0,2)] * matrix[0][2];
            tensa[0][1][0] = drdxyz[GET_ID2(i,3,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,3,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,3,2)] * matrix[0][2];
            tensa[0][2][0] = drdxyz[GET_ID2(i,4,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,4,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,4,2)] * matrix[0][2];

            tensa[0][0][1] = drdxyz[GET_ID2(i,0,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,0,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,0,2)] * matrix[1][2];
            tensa[0][1][1] = drdxyz[GET_ID2(i,3,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,3,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,3,2)] * matrix[1][2];
            tensa[0][2][1] = drdxyz[GET_ID2(i,4,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,4,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,4,2)] * matrix[1][2];

            tensa[0][0][2] = drdxyz[GET_ID2(i,0,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,0,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,0,2)] * matrix[2][2];
            tensa[0][1][2] = drdxyz[GET_ID2(i,3,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,3,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,3,2)] * matrix[2][2];
            tensa[0][2][2] = drdxyz[GET_ID2(i,4,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,4,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,4,2)] * matrix[2][2];

            tensa[1][0][0] = drdxyz[GET_ID2(i,3,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,3,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,3,2)] * matrix[0][2];
            tensa[1][1][0] = drdxyz[GET_ID2(i,1,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,1,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,1,2)] * matrix[0][2];
            tensa[1][2][0] = drdxyz[GET_ID2(i,5,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,5,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,5,2)] * matrix[0][2];

            tensa[1][0][1] = drdxyz[GET_ID2(i,3,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,3,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,3,2)] * matrix[1][2];
            tensa[1][1][1] = drdxyz[GET_ID2(i,1,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,1,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,1,2)] * matrix[1][2];
            tensa[1][2][1] = drdxyz[GET_ID2(i,5,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,5,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,5,2)] * matrix[1][2];

            tensa[1][0][2] = drdxyz[GET_ID2(i,3,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,3,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,3,2)] * matrix[2][2];
            tensa[1][1][2] = drdxyz[GET_ID2(i,1,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,1,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,1,2)] * matrix[2][2];
            tensa[1][2][2] = drdxyz[GET_ID2(i,5,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,5,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,5,2)] * matrix[2][2];

            tensa[2][0][0] = drdxyz[GET_ID2(i,4,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,4,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,4,2)] * matrix[0][2];
            tensa[2][1][0] = drdxyz[GET_ID2(i,5,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,5,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,5,2)] * matrix[0][2];
            tensa[2][2][0] = drdxyz[GET_ID2(i,2,0)] * matrix[0][0] +
                             drdxyz[GET_ID2(i,2,1)] * matrix[0][1] +
                             drdxyz[GET_ID2(i,2,2)] * matrix[0][2];

            tensa[2][0][1] = drdxyz[GET_ID2(i,4,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,4,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,4,2)] * matrix[1][2];
            tensa[2][1][1] = drdxyz[GET_ID2(i,5,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,5,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,5,2)] * matrix[1][2];
            tensa[2][2][1] = drdxyz[GET_ID2(i,2,0)] * matrix[1][0] +
                             drdxyz[GET_ID2(i,2,1)] * matrix[1][1] +
                             drdxyz[GET_ID2(i,2,2)] * matrix[1][2];

            tensa[2][0][2] = drdxyz[GET_ID2(i,4,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,4,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,4,2)] * matrix[2][2];
            tensa[2][1][2] = drdxyz[GET_ID2(i,5,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,5,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,5,2)] * matrix[2][2];
            tensa[2][2][2] = drdxyz[GET_ID2(i,2,0)] * matrix[2][0] +
                             drdxyz[GET_ID2(i,2,1)] * matrix[2][1] +
                             drdxyz[GET_ID2(i,2,2)] * matrix[2][2];

            for (j = 0; j < 3; j++)
              for (l = 0; l < 3; l++)
                for (k = 0; k < 3; k++) {
                  tensb[j][l][k] = 0.;
                  for (m = 0; m < 3; m++)
                    tensb[j][l][k] += matrix[j][m] * tensa[l][m][k];
                }

            for (k = 0; k < 3; k++)
              for (i1 = 0; i1 < 3; i1++)
                for (i2 = 0; i2 < 3; i2++) {
                  tensc[k][i1][i2] = 0.;
                  for (l = 0; l < 3; l++)
                    tensc[k][i1][i2] += matrix[k][l] * tensb[i1][l][i2];
                }

            drdxyz[GET_ID2(i,0,0)] = tensc[0][0][0];
            drdxyz[GET_ID2(i,0,1)] = tensc[0][0][1];
            drdxyz[GET_ID2(i,0,2)] = tensc[0][0][2];

            drdxyz[GET_ID2(i,3,0)] = tensc[0][1][0];
            drdxyz[GET_ID2(i,3,1)] = tensc[0][1][1];
            drdxyz[GET_ID2(i,3,2)] = tensc[0][1][2];

            drdxyz[GET_ID2(i,4,0)] = tensc[0][2][0];
            drdxyz[GET_ID2(i,4,1)] = tensc[0][2][1];
            drdxyz[GET_ID2(i,4,2)] = tensc[0][2][2];

            drdxyz[GET_ID2(i,1,0)] = tensc[1][1][0];
            drdxyz[GET_ID2(i,1,1)] = tensc[1][1][1];
            drdxyz[GET_ID2(i,1,2)] = tensc[1][1][2];

            drdxyz[GET_ID2(i,5,0)] = tensc[1][2][0];
            drdxyz[GET_ID2(i,5,1)] = tensc[1][2][1];
            drdxyz[GET_ID2(i,5,2)] = tensc[1][2][2];

            drdxyz[GET_ID2(i,2,0)] = tensc[2][2][0];
            drdxyz[GET_ID2(i,2,1)] = tensc[2][2][1];
            drdxyz[GET_ID2(i,2,2)] = tensc[2][2][2];

          } /* End of loop on extended ghost cells */

        } /* End if an extended halo exists */

      } /* End of loop on ranks */

    } /* If the transformation is a rotation */

  } /* End of loop on transformations */

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Apply transformation on coordinates.
 *
 * parameters:
 *   halo      <-> halo associated with coordinates to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   coords    --> coordinates on which transformation have to be done.
 *----------------------------------------------------------------------------*/

void
cs_perio_sync_coords(const cs_halo_t *halo,
                     cs_halo_type_t   sync_mode,
                     cs_real_t       *coords)
{
  cs_int_t  i, rank_id, t_id, shift;
  cs_int_t  start_std, start_ext, end_std, end_ext, length;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const fvm_periodicity_t  *periodicity = cs_glob_mesh->periodicity;
  const cs_int_t  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts = halo->n_local_elts;

  if (sync_mode == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  _sync_loc_var_strided(halo, sync_mode, coords, 3);

  /* Compute the new cell centers through periodicity */

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);
    fvm_periodicity_get_matrix(periodicity, t_id, matrix);

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start_std = halo->perio_lst[shift + 4*rank_id];
      length = halo->perio_lst[shift + 4*rank_id + 1];
      end_std = start_std + length;

      if (sync_mode == CS_HALO_EXTENDED) {

        start_ext = halo->perio_lst[shift + 4*rank_id + 2];
        length = halo->perio_lst[shift + 4*rank_id + 3];
        end_ext = start_ext + length;

      }

      /* apply transformation for standard halo */

      for (i = start_std; i < end_std; i++)
        _apply_vector_transfo(matrix, n_elts+i, n_elts+i, coords);

      /* apply transformation for extended halo */

      if (sync_mode == CS_HALO_EXTENDED) {

        for (i = start_ext; i < end_ext; i++)
          _apply_vector_transfo(matrix, n_elts+i, n_elts+i, coords);

      } /* End if extended halo */

    } /* End of loop on ranks */

  } /* End of loop on transformation */

}

/*----------------------------------------------------------------------------
 * Initialize cs_mesh_quantities_t elements for periodicity.
 *
 * Updates cell center and cell family for halo cells.
 *----------------------------------------------------------------------------*/

void
cs_perio_sync_geo(void)
{
  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_real_t  *xyzcen = cs_glob_mesh_quantities->cell_cen;
  cs_real_t  *cell_volume = cs_glob_mesh_quantities->cell_vol;

  const cs_halo_type_t  sync_mode = mesh->halo_type;

  /* Compute the new cell centers through periodicity */

  cs_perio_sync_coords(mesh->halo, sync_mode, xyzcen);

  /* Update cell volume for cells in halo */

  cs_perio_sync_var_scal(mesh->halo,
                         sync_mode,
                         CS_PERIO_ROTA_COPY,
                         cell_volume);
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real scalar between periodic elements.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   rota_mode --> Kind of treatment to do on periodic cells of the halo:
 *                 COPY, IGNORE or RESET
 *   var       <-> scalar to synchronize
 *----------------------------------------------------------------------------*/

void
cs_perio_sync_var_scal(const cs_halo_t *halo,
                       cs_halo_type_t   sync_mode,
                       cs_perio_rota_t  rota_mode,
                       cs_real_t        var[])
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const cs_int_t  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;

  if (sync_mode == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  _sync_loc_var(halo, sync_mode, var);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    /* Reset values if mode is CS_PERIO_ROTA_RESET */

    if (   rota_mode == CS_PERIO_ROTA_RESET
        && perio_type >= FVM_PERIODICITY_ROTATION) {

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

        }

        for (i = start_std; i < end_std; i++)
          var[n_elts + i] = 0.;

        if (sync_mode == CS_HALO_EXTENDED) {
          for (i = start_ext; i < end_ext; i++)
            var[n_elts + i] = 0.;
        }

      } /* End of loop on ranks */

    }

  } /* End of loop on transformations */

}

/*----------------------------------------------------------------------------
 * Synchronize values for a real vector between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   rota_mode --> Kind of treatment to do on periodic cells of the halo:
 *                 COPY, IGNORE or RESET
 *   var_x     <-> component of the vector to update
 *   var_y     <-> component of the vector to update
 *   var_z     <-> component of the vector to update
 *----------------------------------------------------------------------------*/

void
cs_perio_sync_var_vect(const cs_halo_t *halo,
                       cs_halo_type_t   sync_mode,
                       cs_perio_rota_t  rota_mode,
                       cs_real_t        var_x[],
                       cs_real_t        var_y[],
                       cs_real_t        var_z[])
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;
  cs_real_t  x_in, y_in, z_in;

  cs_real_t matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const cs_int_t  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;

  if (sync_mode == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  _sync_loc_var(halo, sync_mode, var_x);
  _sync_loc_var(halo, sync_mode, var_y);
  _sync_loc_var(halo, sync_mode, var_z);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

        }

        if (rota_mode == CS_PERIO_ROTA_COPY) {

          for (i = start_std; i < end_std; i++) {

            x_in = var_x[n_elts + i];
            y_in = var_y[n_elts + i];
            z_in = var_z[n_elts + i];

            _apply_vector_rotation(matrix,
                                   x_in,
                                   y_in,
                                   z_in,
                                   &var_x[n_elts+i],
                                   &var_y[n_elts+i],
                                   &var_z[n_elts+i]);
          }

          if (sync_mode == CS_HALO_EXTENDED) {

            for (i = start_ext; i < end_ext; i++) {

              x_in = var_x[n_elts + i];
              y_in = var_y[n_elts + i];
              z_in = var_z[n_elts + i];

              _apply_vector_rotation(matrix,
                                     x_in,
                                     y_in,
                                     z_in,
                                     &var_x[n_elts+i],
                                     &var_y[n_elts+i],
                                     &var_z[n_elts+i]);

            }

          }

        } /* End if rota_mode == CS_PERIO_ROTA_COPY */

        if (rota_mode == CS_PERIO_ROTA_RESET) {

          for (i = start_std; i < end_std; i++) {
            var_x[n_elts + i] = 0;
            var_y[n_elts + i] = 0;
            var_z[n_elts + i] = 0;
          }

          if (sync_mode == CS_HALO_EXTENDED) {

            for (i = start_ext; i < end_ext; i++) {
              var_x[n_elts + i] = 0;
              var_y[n_elts + i] = 0;
              var_z[n_elts + i] = 0;
            }

          }

        } /* End if rota_mode == CS_PERIO_ROTA_RESET */

      } /* End of loop on ranks */

    } /* End of the treatment of rotation */

  } /* End of loop on transformations */

}

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var11     <-> component of the tensor to update
 *   var12     <-> component of the tensor to update
 *   var13     <-> component of the tensor to update
 *   var21     <-> component of the tensor to update
 *   var22     <-> component of the tensor to update
 *   var23     <-> component of the tensor to update
 *   var31     <-> component of the tensor to update
 *   var32     <-> component of the tensor to update
 *   var33     <-> component of the tensor to update
 *----------------------------------------------------------------------------*/

void
cs_perio_sync_var_tens(const cs_halo_t *halo,
                       cs_halo_type_t   sync_mode,
                       cs_real_t        var11[],
                       cs_real_t        var12[],
                       cs_real_t        var13[],
                       cs_real_t        var21[],
                       cs_real_t        var22[],
                       cs_real_t        var23[],
                       cs_real_t        var31[],
                       cs_real_t        var32[],
                       cs_real_t        var33[])
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;
  cs_real_t  v11, v12, v13, v21, v22, v23, v31, v32, v33;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const cs_int_t  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;

  if (sync_mode == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  _sync_loc_var(halo, sync_mode, var11);
  _sync_loc_var(halo, sync_mode, var12);
  _sync_loc_var(halo, sync_mode, var13);
  _sync_loc_var(halo, sync_mode, var21);
  _sync_loc_var(halo, sync_mode, var22);
  _sync_loc_var(halo, sync_mode, var23);
  _sync_loc_var(halo, sync_mode, var31);
  _sync_loc_var(halo, sync_mode, var32);
  _sync_loc_var(halo, sync_mode, var33);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

        }

        for (i = start_std; i < end_std; i++) {

          v11 = var11[n_elts + i];
          v12 = var12[n_elts + i];
          v13 = var13[n_elts + i];
          v21 = var21[n_elts + i];
          v22 = var22[n_elts + i];
          v23 = var23[n_elts + i];
          v31 = var31[n_elts + i];
          v32 = var32[n_elts + i];
          v33 = var33[n_elts + i];

          _apply_tensor_rotation(matrix,
                                 v11, v12, v13, v21, v22, v23,
                                 v31, v32, v33,
                                 &var11[n_elts + i], &var12[n_elts + i],
                                 &var13[n_elts + i], &var21[n_elts + i],
                                 &var22[n_elts + i], &var23[n_elts + i],
                                 &var31[n_elts + i], &var32[n_elts + i],
                                 &var33[n_elts + i]);

        }

        if (sync_mode == CS_HALO_EXTENDED) {

          for (i = start_ext; i < end_ext; i++) {

            v11 = var11[n_elts + i];
            v12 = var12[n_elts + i];
            v13 = var13[n_elts + i];
            v21 = var21[n_elts + i];
            v22 = var22[n_elts + i];
            v23 = var23[n_elts + i];
            v31 = var31[n_elts + i];
            v32 = var32[n_elts + i];
            v33 = var33[n_elts + i];

            _apply_tensor_rotation(matrix,
                                   v11, v12, v13, v21, v22, v23,
                                   v31, v32, v33,
                                   &var11[n_elts + i], &var12[n_elts + i],
                                   &var13[n_elts + i], &var21[n_elts + i],
                                   &var22[n_elts + i], &var23[n_elts + i],
                                   &var31[n_elts + i], &var32[n_elts + i],
                                   &var33[n_elts + i]);

          }

        } /* End of the treatment of rotation */

      } /* End if halo is extended */

    } /* End of loop on ranks */

  } /* End of loop on transformations for the local rank */

}

/*----------------------------------------------------------------------------
 * Synchronize values for a real diagonal tensor between periodic cells.
 *
 * We only know the diagonal of the tensor.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var11     <-> component of the tensor to update
 *   var22     <-> component of the tensor to update
 *   var33     <-> component of the tensor to update
 *----------------------------------------------------------------------------*/

void
cs_perio_sync_var_diag(const cs_halo_t *halo,
                       cs_halo_type_t   sync_mode,
                       cs_real_t        var11[],
                       cs_real_t        var22[],
                       cs_real_t        var33[])
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;
  cs_real_t  v11, v22, v33;
  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const cs_int_t  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;

  if (sync_mode == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  _sync_loc_var(halo, sync_mode, var11);
  _sync_loc_var(halo, sync_mode, var22);
  _sync_loc_var(halo, sync_mode, var33);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

        }

        for (i = start_std; i < end_std; i++) {

          v11 = var11[n_elts + i];
          v22 = var22[n_elts + i];
          v33 = var33[n_elts + i];

          _apply_tensor_rotation(matrix, v11, 0, 0,
                                 0, v22, 0, 0, 0, v33,
                                 &var11[n_elts + i], NULL, NULL,
                                 NULL, &var22[n_elts + i], NULL,
                                 NULL, NULL, &var33[n_elts + i]);

        }

        if (sync_mode == CS_HALO_EXTENDED) {

          for (i = start_ext; i < end_ext; i++) {

            v11 = var11[n_elts + i];
            v22 = var22[n_elts + i];
            v33 = var33[n_elts + i];

            _apply_tensor_rotation(matrix, v11, 0, 0,
                                   0, v22, 0, 0, 0, v33,
                                   &var11[n_elts + i], NULL, NULL,
                                   NULL, &var22[n_elts + i], NULL,
                                   NULL, NULL, &var33[n_elts + i]);

          }

        } /* End if halo is extended */

      } /* End of loop on ranks */

    } /* End of the treatment of rotation */

  } /* End of loop on transformations */

}

/*----------------------------------------------------------------------------
 * Save rotation terms of a halo to a buffer.
 *
 * parameters:
 *   halo        --> pointer to halo structure
 *   op_type     --> kind of halo treatment (standard or extended)
 *   var         --> variable whose halo rotation terms are to be saved
 *                   (size: halo->n_local_elts + halo->n_elts[opt_type])
 *   save_buffer <-- buffer in which the rotation halo terms are to be saved
 *                   (size: halo->n_elts[op_type])
 *
 * returns:
 *   local number of values saved or restored.
 *----------------------------------------------------------------------------*/

size_t
cs_perio_save_rotation_halo(const cs_halo_t   *halo,
                            cs_halo_type_t     op_type,
                            const cs_real_t    var[],
                            cs_real_t          save_buffer[])
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  size_t save_count = 0;

  const int  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts  = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = mesh->periodicity;

  if (op_type == CS_HALO_N_TYPES)
    return 0;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  /* Loop on transforms */

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = n_elts + halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        for (i = start_std; i < end_std; i++)
          save_buffer[save_count++] = var[i];

        if (op_type == CS_HALO_EXTENDED) {

          start_ext = n_elts + halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

          for (i = start_ext; i < end_ext; i++)
            save_buffer[save_count++] = var[i];

        }

      }

    } /* End if perio_type >= FVM_PERIODICITY_ROTATION) */

  }

  return save_count;
}

/*----------------------------------------------------------------------------
 * Restore rotation terms of a halo from a buffer.
 *
 * parameters:
 *   halo        --> pointer to halo structure
 *   op_type     --> kind of halo treatment (standard or extended)
 *   var         <-> variable whose halo rotation terms are to be restored
 *   save_buffer --> buffer in which the rotation halo terms were saved
 *                   (size: halo->n_elts[op_type])
 *
 * returns:
 *   local number of values saved or restored.
 *----------------------------------------------------------------------------*/

size_t
cs_perio_restore_rotation_halo(const cs_halo_t   *halo,
                               cs_halo_type_t     op_type,
                               cs_real_t          var[],
                               const cs_real_t    save_buffer[])
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  size_t restore_count = 0;

  const int  n_transforms = halo->n_transforms;
  const cs_int_t  n_elts  = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = mesh->periodicity;

  if (op_type == CS_HALO_N_TYPES)
    return 0;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  /* Loop on transforms */

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = n_elts + halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        for (i = start_std; i < end_std; i++)
          var[i] = save_buffer[restore_count++];

        if (op_type == CS_HALO_EXTENDED) {

          start_ext = n_elts + halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

          for (i = start_ext; i < end_ext; i++)
            var[i] = save_buffer[restore_count++];

        }

      }

    } /* End if perio_type >= FVM_PERIODICITY_ROTATION) */

  }

  return restore_count;
}

/*----------------------------------------------------------------------------
 * Define parameters for building an interface set structure on the main mesh.
 *
 * parameters:
 *   p_n_periodic_lists   <-> pointer to the number of periodic lists (may
 *                            be local)
 *   p_periodic_num       <-> pointer to the periodicity number (1 to n)
 *                            associated with each periodic list (primary
 *                            periodicities only)
 *   p_n_periodic_couples <-> pointer to the number of periodic couples
 *   p_periodic_couples   <-> pointer to the periodic couple list
 *----------------------------------------------------------------------------*/

void
cs_perio_define_couples(int         *p_n_periodic_lists,
                        int         *p_periodic_num[],
                        fvm_lnum_t  *p_n_periodic_couples[],
                        fvm_gnum_t  **p_periodic_couples[])
{
  cs_int_t  i, per_id, rank_id, shift;

#if defined(_CS_HAVE_MPI)
  cs_int_t  ierror;
  MPI_Status  status;
#endif

  cs_int_t  n_elts = 0, n_per_elts = 0;

  int  *periodic_num = NULL;
  cs_int_t  *vtx_rank_counter = NULL, *perio_index = NULL;
  fvm_lnum_t  *n_periodic_couples = NULL;
  fvm_gnum_t  *send_buffer = NULL, *recv_buffer = NULL;
  fvm_gnum_t  **periodic_couples = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_int_t  n_perio = mesh->n_init_perio;
  const cs_int_t  n_ranks = mesh->n_domains;
  cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0 : cs_glob_base_rang;

  /* Allocate and initialize buffers */

  BFT_MALLOC(periodic_num, n_perio, int);
  BFT_MALLOC(n_periodic_couples, n_perio, fvm_lnum_t);

  for (i = 0; i < n_perio; i++)
    periodic_num[i] = i+1;

  _define_periodic_buffers(mesh,
                           &perio_index,
                           &send_buffer,
                           &recv_buffer,
                           &vtx_rank_counter);

  /* Send/Receive communication */

#if defined(_CS_HAVE_MPI)
  if (n_ranks > 1) {

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {

      if (rank_id != local_rank) {

        /* If there are elements to send and receive */

        if (vtx_rank_counter[rank_id] > 0) {

          shift = perio_index[2*n_perio*rank_id];
          ierror = MPI_Sendrecv((void *)&send_buffer[shift],
                                vtx_rank_counter[rank_id],
                                FVM_MPI_GNUM, rank_id, 0,
                                (void *)&recv_buffer[shift],
                                vtx_rank_counter[rank_id],
                                FVM_MPI_GNUM, rank_id, 0,
                                cs_glob_base_mpi_comm, &status);

        }

      } /* If local_rank != rank_id */

    } /* End of loop on ranks */

  } /* If n_ranks > 1 */
#endif

#if 0
  {
    cs_int_t  j;

    bft_printf("\n\n\t -------- SEND_BUFFER\t ---------\n");
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      bft_printf("\n\trank_id: %d\n\n", rank_id);
      for (i = 0; i < 2*n_perio; i++) {
        shift = 2*n_perio*rank_id + i;
        bft_printf("Perio_id: %d\t idx_start: %d, idx_end  : %d\n",
                   i, perio_index[shift], perio_index[shift+1]);
        for (j = perio_index[shift]; j < perio_index[shift+1]; j++)
        bft_printf("%5u\n", send_buffer[j]);
      }
    }

    bft_printf("\n\n\t -------- RECV_BUFFER\t ---------\n");
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      bft_printf("\n\trank_id: %d\n\n", rank_id);
      for (i = 0; i < 2*n_perio; i++) {
        shift = 2*n_perio*rank_id + i;
        bft_printf("Perio_id: %d\tidx_start: %d, idx_end : %d\n",
                   i, perio_index[shift], perio_index[shift+1]);
        for (j = perio_index[shift]; j < perio_index[shift+1]; j++)
          bft_printf("%5u\n", recv_buffer[j]);
      }
    }
  }
#endif

  BFT_MALLOC(periodic_couples, n_perio, fvm_gnum_t *);

  for (per_id = 0; per_id < n_perio; per_id++) {

    n_per_elts = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {

      shift = 2*n_perio*rank_id + 2*per_id;
      n_elts = perio_index[shift+1] - perio_index[shift];
      n_per_elts += n_elts;

    } /* End of loop on ranks */

    n_periodic_couples[per_id] = n_per_elts;

  } /* End of loop on periodicities */

  /* Define periodic couples in the global numbering */

  for (per_id = 0; per_id < n_perio; per_id++) {

    BFT_MALLOC(periodic_couples[per_id],
               2*n_periodic_couples[per_id],
               fvm_gnum_t);

    n_per_elts = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {

      shift = 2*n_perio*rank_id + 2*per_id;

      for (i = perio_index[shift]; i < perio_index[shift+1]; i++) {

        periodic_couples[per_id][2*n_per_elts] = send_buffer[i];
        periodic_couples[per_id][2*n_per_elts+1] = recv_buffer[i];
        n_per_elts++;

      }

    } /* End of loop on ranks */

  } /* End of loop on periodicities */

  /* Free memory */

  BFT_FREE(vtx_rank_counter);
  BFT_FREE(perio_index);
  BFT_FREE(send_buffer);

  /* Values to return. Affect pointers */

  *p_n_periodic_lists = n_perio;
  *p_periodic_num = periodic_num;
  *p_n_periodic_couples = n_periodic_couples;
  *p_periodic_couples = periodic_couples;

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
