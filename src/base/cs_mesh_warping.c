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
 * Cut warped faces in serial or parallel with/without periodicity.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
#include <fvm_io_num.h>
#include <fvm_order.h>
#include <fvm_triangulate.h>
#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_connect.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_warping.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create the list of faces to cut in order to respect warping criterion.
 *
 * parameters:
 *   n_faces                 --> number of faces
 *   max_warp_angle          --> criterion above which face is cut
 *   face_warping            --> face warping angle
 *   p_n_warp_faces          <-> pointer to the number of warped faces
 *   p_warp_face_lst         <-> pointer to the warped face list
 *----------------------------------------------------------------------------*/

static void
_select_warped_faces(cs_int_t        n_faces,
                     double          max_warp_angle,
                     double          face_warping[],
                     cs_int_t        *p_n_warp_faces,
                     cs_int_t        *p_warp_face_lst[])
{
  cs_int_t  face_id;

  cs_int_t  n_warp_faces = 0;
  cs_int_t  *warp_face_lst = NULL;

  if (n_faces > 0) {

    for (face_id = 0; face_id < n_faces; face_id++)
      if (face_warping[face_id] >= max_warp_angle)
        n_warp_faces++;

    BFT_MALLOC(warp_face_lst, n_warp_faces, cs_int_t);

    n_warp_faces = 0;

    for (face_id = 0; face_id < n_faces; face_id++)
      if (face_warping[face_id] >= max_warp_angle)
        warp_face_lst[n_warp_faces++] = face_id + 1;

  }

  *p_n_warp_faces = n_warp_faces;
  *p_warp_face_lst = warp_face_lst;

}

/*----------------------------------------------------------------------------
 * Define periodic index for send/receive operations.
 *
 * parameters:
 *   mesh                 --> pointer to a mesh structure
 *   face_warping         --> face warping angle
 *   max_warp_angle       --> criterion above which face is cut
 *   p_domain_to_c_ranks  <-> pointer to a rank indirection array
 *   p_s_rank_index       <-> pointer to the sending index on ranks
 *   p_r_rank_index       <-> pointer to the receiving index on ranks
 *----------------------------------------------------------------------------*/

static void
_define_periodic_index(const cs_mesh_t   *const mesh,
                       double            face_warping[],
                       double            max_warp_angle,
                       cs_int_t          *p_domain_to_c_rank[],
                       cs_int_t          *p_s_rank_index[],
                       cs_int_t          *p_r_rank_index[])
{
  cs_int_t  i, face_id, perio_id, rank_id, shift;
  cs_int_t  n_triangles = 0;

  int  cpt_request = 0;
  cs_int_t  *domain_to_c_rank = NULL;
  cs_int_t  *s_rank_index = NULL, *r_rank_index = NULL;

#if defined(_CS_HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  const cs_int_t  n_domains = mesh->n_domains;
  const cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;
  const cs_halo_t  *halo = mesh->halo;
  const cs_int_t  n_c_ranks = halo->n_c_domains;
  const cs_int_t  *per_face_idx = cs_glob_mesh_builder->per_face_idx;
  const cs_int_t  *per_face_lst = cs_glob_mesh_builder->per_face_lst;
  const cs_int_t  *per_rank_lst = cs_glob_mesh_builder->per_rank_lst;

#if defined(_CS_HAVE_MPI)
  if (halo->n_c_domains*2 > 128) {
    BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);
  }
#endif

  BFT_MALLOC(s_rank_index, n_c_ranks + 1, cs_int_t);
  BFT_MALLOC(r_rank_index, n_c_ranks + 1, cs_int_t);

  for (i = 0; i < n_c_ranks + 1; i++) {
    s_rank_index[i] = 0;
    r_rank_index[i] = 0;
  }

  BFT_MALLOC(domain_to_c_rank, n_domains, cs_int_t);

  for (i = 0, shift = 0; i < n_domains; i++) {

    if (halo->c_domain_rank[shift] == i) {
      domain_to_c_rank[i] = shift;
      shift++;
    }
    else
      domain_to_c_rank[i] = -2;

  } /* End of loop on ranks */

  for (perio_id = 0; perio_id < mesh->n_init_perio; perio_id++) {

    for (i = per_face_idx[perio_id]; i < per_face_idx[perio_id+1]; i++) {

      if (per_face_lst[2*i] > 0) {

        face_id = per_face_lst[2*i] - 1;

        if (face_warping[face_id] > max_warp_angle) {

          if (n_domains > 1)
            rank_id = domain_to_c_rank[per_rank_lst[i]];
          else
            rank_id = 0;

          s_rank_index[rank_id + 1] += n_triangles*3 + 1;

        }

      } /* If this is a direct transformation */

    } /* End of loop on periodic index */

  } /* End of loop on periodicities */

  /* Send/Recv number of elements to exchange */

  for (rank_id = 0; rank_id < n_c_ranks; rank_id++) {

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Irecv(&(r_rank_index[rank_id+1]), 1, MPI_INT,
                halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                cs_glob_base_mpi_comm,
                &(request[cpt_request++]));
#endif

    }
    else
      r_rank_index[rank_id+1] = s_rank_index[rank_id+1];

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Barrier(cs_glob_base_mpi_comm);
#endif

  for (rank_id = 0; rank_id < n_c_ranks; rank_id++) {

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
        MPI_Isend(&(s_rank_index[rank_id+1]), 1, MPI_INT,
                  halo->c_domain_rank[rank_id], local_rank,
                  cs_glob_base_mpi_comm,
                  &(request[cpt_request++]));
#endif

    }

  } /* End of loop on ranks */

  /* Sync after each communicating rank had received all the messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Waitall(cpt_request, request, status);
#endif

  /* Define send/receive index */

  for (i = 0; i < n_c_ranks; i++) {
    s_rank_index[i+1] += s_rank_index[i];
    r_rank_index[i+1] += r_rank_index[i];
  }

  *p_domain_to_c_rank = domain_to_c_rank;
  *p_s_rank_index = s_rank_index;
  *p_r_rank_index = r_rank_index;

#if defined(_CS_HAVE_MPI)
  if (request != _request) {
    BFT_FREE(request);
    BFT_FREE(status);
  }
#endif

}

/*----------------------------------------------------------------------------
 * Fill periodic buffer to send
 *
 * parameters:
 *   mesh                 --> pointer to a mesh structure
 *   face_warping         --> face warping angle
 *   max_warp_angle       --> criterion above which face is cut
 *   p_domain_to_c_ranks  <-> pointer to a rank indirection array
 *   s_rank_index         <-> sending index on ranks
 *   r_rank_index         <-> receiving index on ranks
 *   perio_s_buffer       <-> periodic data to send
 *   perio_r_buffer       <-> periodic data to receive
 *----------------------------------------------------------------------------*/

static void
_fill_perio_buffers(const cs_mesh_t   *mesh,
                    double             face_warping[],
                    double             max_warp_angle,
                    cs_int_t           domain_to_c_rank[],
                    cs_int_t           old_to_new[],
                    cs_int_t           old_face_vtx_idx[],
                    cs_int_t           new_face_vtx_idx[],
                    cs_int_t           new_face_vtx_lst[],
                    cs_int_t           s_rank_index[],
                    cs_int_t           r_rank_index[])
{
  cs_int_t  i, j, k, perio_id, local_face_id,  rank_id, new_face_id;
  cs_int_t  perio_shift, shift;
  cs_int_t  n_face_vertices, n_triangles, n_elts;

  int  cpt_request = 0;
  cs_int_t  *counter = NULL;
  cs_int_t  *perio_s_buffer = NULL, *perio_r_buffer = NULL;

#if defined(_CS_HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  const cs_int_t  n_domains = mesh->n_domains;
  const cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;
  const cs_halo_t  *halo = mesh->halo;
  const cs_int_t  *per_face_idx = cs_glob_mesh_builder->per_face_idx;
  const cs_int_t  *per_face_lst = cs_glob_mesh_builder->per_face_lst;
  const cs_int_t  *per_rank_lst = cs_glob_mesh_builder->per_rank_lst;

#if defined(_CS_HAVE_MPI)
  if (halo->n_c_domains*2 > 128) {
    BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);
  }
#endif

  /* Allocation and intialization */

  BFT_MALLOC(perio_s_buffer, s_rank_index[halo->n_c_domains], cs_int_t);
  BFT_MALLOC(perio_r_buffer, r_rank_index[halo->n_c_domains], cs_int_t);
  BFT_MALLOC(counter, halo->n_c_domains, cs_int_t);

  for (i = 0; i < halo->n_c_domains; i++)
    counter[i] = 0;

  /* Define perio_s_buffer */

  for (perio_id = 0; perio_id < mesh->n_init_perio; perio_id++) {

    for (i = per_face_idx[perio_id]; i < per_face_idx[perio_id+1]; i++) {

      if (per_face_lst[2*i] > 0) {

        local_face_id = per_face_lst[2*i] - 1;

        if (face_warping[local_face_id] > max_warp_angle) {

          if (n_domains > 1)
            rank_id = domain_to_c_rank[per_rank_lst[i]];
          else
            rank_id = 0;

          perio_shift = s_rank_index[rank_id] + counter[rank_id];

          perio_s_buffer[perio_shift++] = per_face_lst[2*i+1];

          n_face_vertices =  old_face_vtx_idx[local_face_id+1]
                           - old_face_vtx_idx[local_face_id];
          n_triangles = n_face_vertices - 2;

          new_face_id = old_to_new[local_face_id];

          for (j = 0; j < n_triangles; j++) {

            for (k = new_face_vtx_idx[new_face_id];
                 k < new_face_vtx_idx[new_face_id+1]; k++)
              perio_s_buffer[perio_shift++] = new_face_vtx_lst[k];

            new_face_id++;

          } /* End of loop on triangles */

          counter[rank_id] += 1 + 3*n_triangles;

        } /* If this face should have been cut */

      } /* If this is a direct transformation */

    } /* End of loop on periodic index */

  } /* End of loop on periodicities */

  /* Exchange buffers */

  for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    n_elts = r_rank_index[rank_id+1] - r_rank_index[rank_id];
    shift = r_rank_index[rank_id];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Irecv(&(perio_r_buffer[shift]), n_elts, CS_MPI_INT,
                halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                cs_glob_base_mpi_comm, &(request[cpt_request++]));
#endif

    }
    else {

      assert(n_elts == (s_rank_index[rank_id+1] - s_rank_index[rank_id]));

      for (i = 0; i < n_elts; i++)
        perio_r_buffer[shift + i] = perio_s_buffer[s_rank_index[rank_id] + i];

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Barrier(cs_glob_base_mpi_comm);
#endif

  for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

    if (halo->c_domain_rank[rank_id] != local_rank) {

      n_elts = s_rank_index[rank_id+1] - s_rank_index[rank_id];
      shift = s_rank_index[rank_id];

#if defined(_CS_HAVE_MPI)
      MPI_Isend(&(perio_s_buffer[shift]), n_elts, MPI_INT,
                halo->c_domain_rank[rank_id], local_rank,
                cs_glob_base_mpi_comm, &(request[cpt_request++]));
#endif

    }

  } /* End of loop on ranks */

  /* Sync after each communicating rank had received all the messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Waitall(cpt_request, request, status);
#endif

  /* Apply perio_r_buffer to the new face->vertices connectivity */

  perio_shift = 0;
  while (perio_shift < r_rank_index[halo->n_c_domains]) {

    local_face_id = perio_r_buffer[perio_shift++];

    n_face_vertices =  old_face_vtx_idx[local_face_id+1]
                     - old_face_vtx_idx[local_face_id];
    n_triangles = n_face_vertices - 2;

    new_face_id = old_to_new[local_face_id];

    for (i = 0; i < n_triangles; i++) {

      for (j = new_face_vtx_idx[new_face_id];
           j < new_face_vtx_idx[new_face_id+1]; j++)
        new_face_vtx_lst[j] = perio_r_buffer[perio_shift++];

      new_face_id++;

    } /* End of loop on triangles */

  } /* End of while */

  /* Free memory */

  BFT_FREE(perio_r_buffer);
  BFT_FREE(perio_s_buffer);
  BFT_FREE(counter);

#if defined(_CS_HAVE_MPI)
  if (request != _request) {
    BFT_FREE(request);
    BFT_FREE(status);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Cut faces with periodic treatment and update connectivity.
 * Only useful for internal faces.
 *
 * parameters:
 *   mesh                    <-> pointer to a mesh structure
 *   face_type               --> internal or border faces
 *   max_warp_angle          --> criterion above which face is cut
 *   face_warping            --> face warping angle
 *   p_n_cut_faces           <-> pointer to the number of cut faces
 *   p_cut_face_lst          <-> pointer to the cut face list
 *   p_n_sub_elt_lst         <-> pointer to the sub-elt count list
 *   p_n_faces               <-> pointer to the number of faces
 *   p_face_vtx_connect_size <-> size of the "face -> vertex" connectivity
 *   p_face_cells            <-> "face -> cells" connectivity
 *   p_face_vtx_idx          <-> pointer on "face -> vertices" connect. index
 *   p_face_vtx_lst          <-> pointer on "face -> vertices" connect. list
 *----------------------------------------------------------------------------*/

static void
_cut_warped_faces_perio(cs_mesh_t       *mesh,
                        double           max_warp_angle,
                        double           face_warping[],
                        cs_int_t        *p_n_cut_faces,
                        cs_int_t        *p_cut_face_lst[],
                        fvm_lnum_t      *p_n_sub_elt_lst[],
                        cs_int_t        *p_n_faces,
                        cs_int_t        *p_face_vtx_connect_size,
                        cs_int_t        *p_face_cells[],
                        cs_int_t        *p_face_family[],
                        cs_int_t        *p_face_vtx_idx[],
                        cs_int_t        *p_face_vtx_lst[])
{
  cs_int_t  i, j, face_id, idx_start, idx_end, old_face_idx;
  cs_int_t  n_triangles, num, shift, face_shift;

  cs_int_t  n_face_vertices = 0, n_max_face_vertices = 0;
  cs_int_t  n_new_faces = 0, n_cut_faces = 0, connect_size = 0;

  fvm_triangulate_state_t  *triangle_state = NULL;
  cs_int_t  *face_connectivity = NULL, *new_face_connectivity = NULL;
  cs_int_t  *new_face_vtx_idx = NULL, *new_face_vtx_lst = NULL;
  cs_int_t  *new_face_cells = NULL, *new_face_family = NULL;
  cs_int_t  *cut_face_lst = NULL;
  cs_int_t  *domain_to_c_rank = NULL, *new_face_shift = NULL;
  cs_int_t  *s_rank_index = NULL, *r_rank_index = NULL;
  fvm_lnum_t  *n_sub_elt_lst = NULL;

  const cs_int_t  dim = mesh->dim;
  const cs_int_t  n_init_faces = *p_n_faces;

  assert(dim == CS_DIM_3);

  BFT_MALLOC(n_sub_elt_lst, n_init_faces, fvm_lnum_t);
  BFT_MALLOC(new_face_shift, n_init_faces, cs_int_t);

  /* First loop: compute sizes */

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    idx_start = (*p_face_vtx_idx)[face_id] - 1;
    idx_end = (*p_face_vtx_idx)[face_id + 1] - 1;

    n_face_vertices = idx_end - idx_start;
    n_max_face_vertices = CS_MAX(n_max_face_vertices, n_face_vertices);

    new_face_shift[face_id] = n_new_faces;

    if (face_warping[face_id] >= max_warp_angle) {

      n_triangles = n_face_vertices - 2;
      connect_size += n_triangles*3;
      n_new_faces += n_triangles;
      n_cut_faces += n_triangles;
      n_sub_elt_lst[face_id] = n_triangles;

    }
    else {

      connect_size += n_face_vertices;
      n_new_faces += 1;
      n_sub_elt_lst[face_id] = 1;

    }

  } /* End of loop on faces */

  *p_n_sub_elt_lst = n_sub_elt_lst;

  if (n_cut_faces == 0) {

    BFT_FREE(new_face_shift);
    return;

  }

  BFT_MALLOC(new_face_vtx_idx, n_new_faces + 1, cs_int_t);
  BFT_MALLOC(new_face_vtx_lst, connect_size, cs_int_t);
  BFT_MALLOC(new_face_cells, 2*n_new_faces, cs_int_t);
  BFT_MALLOC(new_face_family, n_new_faces, cs_int_t);

  BFT_MALLOC(cut_face_lst, n_cut_faces, cs_int_t);

  triangle_state = fvm_triangulate_state_create(n_max_face_vertices);

  BFT_MALLOC(face_connectivity, n_max_face_vertices, cs_int_t);
  BFT_MALLOC(new_face_connectivity, (n_max_face_vertices-2)*3, cs_int_t);

  /* Periodic treatment */

  _define_periodic_index(mesh,
                         face_warping,
                         max_warp_angle,
                         &domain_to_c_rank,
                         &s_rank_index,
                         &r_rank_index);

  /* Second loop : define the new connectivity after triangulation */

  new_face_vtx_idx[0] = 1;
  connect_size = 0;
  n_new_faces = 0;
  n_cut_faces = 0;

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    idx_start = (*p_face_vtx_idx)[face_id] - 1;
    idx_end = (*p_face_vtx_idx)[face_id + 1] - 1;
    n_face_vertices = idx_end - idx_start;

    if (face_warping[face_id] >= max_warp_angle) {

      for (j = 0, i = idx_start; i < idx_end; i++, j++)
        face_connectivity[j] = (*p_face_vtx_lst)[i];

      n_triangles = fvm_triangulate_polygon(dim,
                                            n_face_vertices,
                                            mesh->vtx_coord,
                                            NULL,
                                            face_connectivity,
                                            FVM_TRIANGULATE_ELT_DEF,
                                            new_face_connectivity,
                                            triangle_state);

      assert(n_triangles == n_face_vertices - 2);

      /* Update face -> vertex connectivity */

      shift = 0;

      for (i = 0; i < n_triangles; i++) {

        cut_face_lst[n_cut_faces++] = n_new_faces + 1;

        /* Update "face -> cells" connectivity */

        for (j = 0; j < 2; j++)
          new_face_cells[2*n_new_faces + j] = (*p_face_cells)[2*face_id + j];

        /* Update family for each face */

        new_face_family[n_new_faces] = (*p_face_family)[face_id];

        /* Update "face -> vertices" connectivity */

        n_new_faces++;
        connect_size += 3;
        new_face_vtx_idx[n_new_faces]
          = new_face_vtx_idx[n_new_faces-1] + 3;

      } /* End of loop on triangles */

    }
    else {

      /* Update "face -> cells" connectivity */

      for (j = 0; j < 2; j++)
        new_face_cells[2*n_new_faces + j] = (*p_face_cells)[2*face_id + j];

      /* Update family for each face */

      new_face_family[n_new_faces] = (*p_face_family)[face_id];

      /* Update "face -> vertices" connectivity */

      for (j = 0, i = idx_start; i < idx_end; i++, j++)
        new_face_vtx_lst[connect_size + j] = (*p_face_vtx_lst)[i];

      n_new_faces++;
      connect_size += n_face_vertices;
      new_face_vtx_idx[n_new_faces]
        = new_face_vtx_idx[n_new_faces-1] + n_face_vertices;

    }

  } /* End of loop on internal faces */

  triangle_state = fvm_triangulate_state_destroy(triangle_state);

  BFT_FREE(face_connectivity);
  BFT_FREE(new_face_connectivity);
  BFT_FREE(*p_face_cells);

  /* Define perio_s_buffer and exchange it with communicating rank to
     fill perio_r_buffer */

  _fill_perio_buffers(mesh,
                      face_warping,
                      max_warp_angle,
                      domain_to_c_rank,
                      new_face_shift,
                      (*p_face_vtx_idx),
                      new_face_vtx_idx,
                      new_face_vtx_lst,
                      s_rank_index,
                      r_rank_index);

  /* Get mesh numbering from element numbering */

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    if (face_warping[face_id] >= max_warp_angle) {

      idx_start = (*p_face_vtx_idx)[face_id] - 1;
      idx_end = (*p_face_vtx_idx)[face_id + 1] - 1;

      n_face_vertices = idx_end - idx_start;
      n_triangles = n_face_vertices - 2;

      face_shift = new_face_shift[face_id];

      old_face_idx = (*p_face_vtx_idx)[face_id] - 1;

      for (i = 0; i < n_triangles; i++) {

        for (j = new_face_vtx_idx[face_shift] - 1;
             j < new_face_vtx_idx[face_shift+1] - 1; j++) {

          shift = new_face_vtx_lst[j] - 1;
          num = (*p_face_vtx_lst)[old_face_idx + shift];
          new_face_vtx_lst[j] = num;

        }

        face_shift++;

      } /* End of loop on triangles */

    } /* If the face is warped */

  } /* End of loop on faces */

  BFT_FREE(new_face_shift);
  BFT_FREE(domain_to_c_rank);
  BFT_FREE(s_rank_index);
  BFT_FREE(r_rank_index);
  BFT_FREE(*p_face_vtx_idx);
  BFT_FREE(*p_face_vtx_lst);
  BFT_FREE(*p_face_family);

  /* Define returned pointers */

  *p_face_vtx_idx = new_face_vtx_idx;
  *p_face_vtx_lst = new_face_vtx_lst;
  *p_face_cells = new_face_cells;
  *p_face_family = new_face_family;
  *p_face_vtx_connect_size = connect_size;
  *p_n_faces = n_new_faces;
  *p_n_cut_faces = n_cut_faces;

  if (p_cut_face_lst != NULL)
    BFT_FREE(*p_cut_face_lst);
  *p_cut_face_lst = cut_face_lst;

}

/*----------------------------------------------------------------------------
 * Cut faces if necessary and update connectivity without periodicity
 *
 * parameters:
 *   mesh                    <-> pointer to a mesh structure
 *   face_type               --> internal or border faces
 *   max_warp_angle          --> criterion above which face is cut
 *   face_warping            --> face warping angle
 *   p_n_cut_faces           <-> pointer to the number of cut faces
 *   p_cut_face_lst          <-> pointer to the cut face list
 *   p_n_sub_elt_lst         <-> pointer to the sub-elt count list
 *   p_n_faces               <-> pointer to the number of faces
 *   p_face_num              <-> pointer to the global face numbers
 *   p_face_vtx_connect_size <-> size of the "face -> vertex" connectivity
 *   p_face_cells            <-> "face -> cells" connectivity
 *   p_face_vtx_idx          <-> pointer on "face -> vertices" connect. index
 *   p_face_vtx_lst          <-> pointer on "face -> vertices" connect. list
 *----------------------------------------------------------------------------*/

static void
_cut_warped_faces(cs_mesh_t       *mesh,
                  int              stride,
                  double           max_warp_angle,
                  double           face_warping[],
                  cs_int_t        *p_n_cut_faces,
                  cs_int_t        *p_cut_face_lst[],
                  fvm_lnum_t      *p_n_sub_elt_lst[],
                  cs_int_t        *p_n_faces,
                  cs_int_t        *p_face_vtx_connect_size,
                  cs_int_t        *p_face_cells[],
                  cs_int_t        *p_face_family[],
                  cs_int_t        *p_face_vtx_idx[],
                  cs_int_t        *p_face_vtx_lst[])
{
  cs_int_t  i, j, face_id, idx_start, idx_end, shift;
  cs_int_t  n_triangles;

  cs_int_t  n_face_vertices = 0, n_max_face_vertices = 0;
  cs_int_t  n_new_faces = 0, n_cut_faces = 0, connect_size = 0;

  fvm_triangulate_state_t  *triangle_state = NULL;
  cs_int_t  *face_connectivity = NULL, *new_face_connectivity = NULL;
  cs_int_t  *new_face_vtx_idx = NULL, *new_face_vtx_lst = NULL;
  cs_int_t  *new_face_cells = NULL, *new_face_family = NULL;
  cs_int_t  *cut_face_lst = NULL;
  fvm_lnum_t  *n_sub_elt_lst = NULL;

  const cs_int_t  dim = mesh->dim;
  const cs_int_t  n_init_faces = *p_n_faces;

  assert(stride == 1 || stride ==2);
  assert(dim == CS_DIM_3);

  BFT_MALLOC(n_sub_elt_lst, n_init_faces, fvm_lnum_t);

  /* First loop: count */

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    idx_start = (*p_face_vtx_idx)[face_id] - 1;
    idx_end = (*p_face_vtx_idx)[face_id + 1] - 1;

    n_face_vertices = idx_end - idx_start;
    n_max_face_vertices = CS_MAX(n_max_face_vertices, n_face_vertices);

    if (face_warping[face_id] >= max_warp_angle) {

      n_triangles = n_face_vertices - 2;
      connect_size += n_triangles*3;
      n_new_faces += n_triangles;
      n_cut_faces += n_triangles;
      n_sub_elt_lst[face_id] = n_triangles;

    }
    else {

      connect_size += n_face_vertices;
      n_new_faces += 1;
      n_sub_elt_lst[face_id] = 1;

    }

  } /* End of loop on faces */

  *p_n_sub_elt_lst = n_sub_elt_lst;

  if (n_cut_faces == 0)
    return;

  BFT_MALLOC(new_face_vtx_idx, n_new_faces + 1, cs_int_t);
  BFT_MALLOC(new_face_vtx_lst, connect_size, cs_int_t);
  BFT_MALLOC(new_face_cells, n_new_faces*stride, cs_int_t);
  BFT_MALLOC(new_face_family, n_new_faces, cs_int_t);

  BFT_MALLOC(cut_face_lst, n_cut_faces, cs_int_t);

  triangle_state = fvm_triangulate_state_create(n_max_face_vertices);

  BFT_MALLOC(face_connectivity, n_max_face_vertices, cs_int_t);
  BFT_MALLOC(new_face_connectivity, (n_max_face_vertices-2)*3, cs_int_t);

  /* Second loop : define the new connectivity after triangulation */

  new_face_vtx_idx[0] = 1;
  connect_size = 0;
  n_new_faces = 0;
  n_cut_faces = 0;

  for (face_id = 0; face_id < n_init_faces; face_id++) {

    idx_start = (*p_face_vtx_idx)[face_id] - 1;
    idx_end = (*p_face_vtx_idx)[face_id + 1] - 1;
    n_face_vertices = idx_end - idx_start;

    if (face_warping[face_id] >= max_warp_angle) {

      for (j = 0, i = idx_start; i < idx_end; i++, j++)
        face_connectivity[j] = (*p_face_vtx_lst)[i];

      n_triangles = fvm_triangulate_polygon(dim,
                                            n_face_vertices,
                                            mesh->vtx_coord,
                                            NULL,
                                            face_connectivity,
                                            FVM_TRIANGULATE_MESH_DEF,
                                            new_face_connectivity,
                                            triangle_state);

      assert(n_triangles == n_face_vertices - 2);

      /* Update face -> vertex connectivity */

      shift = 0;

      for (i = 0; i < n_triangles; i++) {

        cut_face_lst[n_cut_faces++] = n_new_faces + 1;

        /* Update "face -> cells" connectivity */

        for (j = 0; j < stride; j++)
          new_face_cells[stride*n_new_faces + j] =
            (*p_face_cells)[stride*face_id + j];

        /* Update family for each face */

        new_face_family[n_new_faces] = (*p_face_family)[face_id];

        /* Update "face -> vertices" connectivity */

        for (j = 0; j < 3; j++)
          new_face_vtx_lst[connect_size + j] = new_face_connectivity[shift++];

        n_new_faces++;
        connect_size += 3;
        new_face_vtx_idx[n_new_faces] = new_face_vtx_idx[n_new_faces-1] + 3;

      } /* End of loop on triangles */

    }
    else {

      /* Update "face -> cells" connectivity */

      for (j = 0; j < stride; j++)
        new_face_cells[stride*n_new_faces + j] =
          (*p_face_cells)[stride*face_id + j];

      /* Update family for each faces */

      new_face_family[n_new_faces] = (*p_face_family)[face_id];

      /* Update "face -> vertices" connectivity */

      for (j = 0, i = idx_start; i < idx_end; i++, j++)
        new_face_vtx_lst[connect_size + j] = (*p_face_vtx_lst)[i];

      n_new_faces++;
      connect_size += n_face_vertices;
      new_face_vtx_idx[n_new_faces] =
        new_face_vtx_idx[n_new_faces-1] + n_face_vertices;

    }

  } /* End of loop on internal faces */

  triangle_state = fvm_triangulate_state_destroy(triangle_state);

  BFT_FREE(face_connectivity);
  BFT_FREE(new_face_connectivity);

  BFT_FREE(*p_face_vtx_idx);
  BFT_FREE(*p_face_vtx_lst);
  BFT_FREE(*p_face_cells);
  BFT_FREE(*p_face_family);

  /* Define returned pointers */

  *p_face_vtx_idx = new_face_vtx_idx;
  *p_face_vtx_lst = new_face_vtx_lst;
  *p_face_cells = new_face_cells;
  *p_face_family = new_face_family;
  *p_face_vtx_connect_size = connect_size;
  *p_n_faces = n_new_faces;
  *p_n_cut_faces = n_cut_faces;

  if (p_cut_face_lst != NULL)
    BFT_FREE(*p_cut_face_lst);
  *p_cut_face_lst = cut_face_lst;

}

/*----------------------------------------------------------------------------
 * Update warped faces global numbers after cutting
 *
 * parameters:
 *   mesh                <-> pointer to a mesh structure
 *   n_faces             --> number of faces
 *   n_init_faces        --> initial number of faces
 *   n_cut_faces         --> number of cut faces
 *   cut_face_lst        --> pointer to the cut face list
 *   n_sub_elt_lst       --> sub-elt count list
 *   n_g_faces           <-> global number of faces
 *   p_global_face_num   <-> pointer to the global face numbers
 *----------------------------------------------------------------------------*/

static void
_update_cut_faces_num(cs_mesh_t       *mesh,
                      cs_int_t         n_faces,
                      cs_int_t         n_init_faces,
                      fvm_lnum_t       n_sub_elt_lst[],
                      fvm_gnum_t      *n_g_faces,
                      fvm_gnum_t     **p_global_face_num)
{
  size_t  size;

  fvm_io_num_t  *new_io_num = NULL, *previous_io_num = NULL;
  const fvm_gnum_t  *global_num = NULL;

  /* Simply update global number of faces in trivial case */

  *n_g_faces = n_faces;

  if (*p_global_face_num == NULL)
    return;

  /* Faces should not have been reordered */

  if (fvm_order_local_test(NULL, *p_global_face_num, n_init_faces) == false)
    bft_error(__FILE__, __LINE__, 0,
              _("On a découpé les faces préalablement renumérotées.\n"
                "Ce cas ne devrait pas se produire, car on devrait découper\n"
                "les entités de maillage avant de les renuméroter."));

  /* Update global number of internal faces and its global numbering */

  if (mesh->n_domains > 1) {

    bft_printf(_("\t%12d faces globalement avant découpage\n"),
               *n_g_faces);

    previous_io_num = fvm_io_num_create(NULL,
                                        *p_global_face_num,
                                        n_init_faces,
                                        0);
    new_io_num = fvm_io_num_create_from_sub(previous_io_num,
                                            n_sub_elt_lst);

    previous_io_num = fvm_io_num_destroy(previous_io_num);

    *n_g_faces = fvm_io_num_get_global_count(new_io_num);

    global_num = fvm_io_num_get_global_num(new_io_num);

    BFT_REALLOC(*p_global_face_num, n_faces, fvm_gnum_t);
    size = sizeof(fvm_gnum_t) * n_faces;
    memcpy(*p_global_face_num, global_num, size);

    bft_printf(_("\t%12d faces globalement après découpage\n\n"),
               *n_g_faces);

    new_io_num = fvm_io_num_destroy(new_io_num);

  }
}

/*----------------------------------------------------------------------------
 * Post-process the warped faces before cutting.
 *
 * parameters:
 *   n_i_warp_faces       --> number of internal warped faces
 *   n_b_warp_faces       --> number of border warped faces
 *   i_warp_face_lst      --> internal warped face list
 *   b_warp_face_lst      --> border warped face list
 *   i_face_warping       --> face warping angle for internal faces
 *   b_face_warping       --> face warping angle for internal faces
 *----------------------------------------------------------------------------*/

static void
_post_before_cutting(cs_int_t        n_i_warp_faces,
                     cs_int_t        n_b_warp_faces,
                     cs_int_t        i_warp_face_lst[],
                     cs_int_t        b_warp_face_lst[],
                     double          i_face_warping[],
                     double          b_face_warping[])
{
  int  n_parent_lists = 2;
  fvm_lnum_t  parent_num_shift[2]  = {0, cs_glob_mesh->n_b_faces};
  fvm_nodal_t  *fvm_mesh = NULL;
  fvm_writer_t  *writer = NULL;

  const cs_int_t  writer_id = -1; /* default writer */

  const void  *var_ptr[2] = {NULL, NULL};

  if (cs_post_existe_writer(writer_id) == false)
    return;

  assert(sizeof(double) == sizeof(cs_real_t));

  fvm_mesh = cs_maillage_extrait_fac_nodal(cs_glob_mesh,
                                           _("Faces gauche a decouper"),
                                           n_i_warp_faces,
                                           n_b_warp_faces,
                                           i_warp_face_lst,
                                           b_warp_face_lst);

  writer = cs_post_get_writer(writer_id);

  fvm_writer_set_mesh_time(writer, -1, 0.0);

  /* Write a mesh from the selected faces */

  fvm_writer_export_nodal(writer, fvm_mesh);

  /* Write the warping field */

  var_ptr[0] = ((const char *)b_face_warping);
  var_ptr[1] = ((const char *)i_face_warping);

  fvm_writer_export_field(writer,
                          fvm_mesh,
                          _("Gauchissement des faces"),
                          1,
                          FVM_WRITER_PER_ELEMENT,
                          FVM_INTERLACE,
                          n_parent_lists,
                          parent_num_shift,
                          FVM_DOUBLE,
                          (int)-1,
                          (double)0.0,
                          (const void **)var_ptr);

  fvm_mesh = fvm_nodal_destroy(fvm_mesh);

}

/*----------------------------------------------------------------------------
 * Post-process the warped faces after cutting.
 *
 * parameters:
 *   n_i_cut_faces       --> number of internal faces generated by cutting
 *   n_b_cut_faces       --> number of border faces generated by cutting
 *   i_cut_face_lst      --> face warping angle for internal faces
 *   b_cut_face_lst      --> face warping angle for internal faces
 *----------------------------------------------------------------------------*/

static void
_post_after_cutting(cs_int_t       n_i_cut_faces,
                    cs_int_t       n_b_cut_faces,
                    cs_int_t       i_cut_face_lst[],
                    cs_int_t       b_cut_face_lst[])
{
  fvm_nodal_t  *fvm_mesh = NULL;
  fvm_writer_t  *writer = NULL;

  const cs_int_t  writer_id = -1; /* default writer */

  if (cs_post_existe_writer(writer_id) == false)
    return;

  fvm_mesh = cs_maillage_extrait_fac_nodal(cs_glob_mesh,
                                           _("Faces gauche apres decoupage"),
                                           n_i_cut_faces,
                                           n_b_cut_faces,
                                           i_cut_face_lst,
                                           b_cut_face_lst);

  writer = cs_post_get_writer(writer_id);

  fvm_writer_set_mesh_time(writer, -1, 0.0);

  /* Write a mesh from the selected faces */

  fvm_writer_export_nodal(writer, fvm_mesh);

  fvm_mesh = fvm_nodal_destroy(fvm_mesh);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Cut warped faces.
 *
 * Update border face connectivity and associated mesh quantities.
 *
 * parameters:
 *   mesh             <-> pointer to mesh structure.
 *   max_warp_angle   --> criterion to know which face to cut
 *   post_tag         --> tag to know if we have to post-treat cut faces.
 *----------------------------------------------------------------------------*/

void
cs_mesh_warping_cut_faces(cs_mesh_t    *mesh,
                          double        max_warp_angle,
                          cs_bool_t     post_tag)
{
  cs_int_t  i;

  cs_int_t  n_i_warp_faces = 0, n_b_warp_faces = 0;
  cs_int_t  n_i_cut_faces = 0, n_b_cut_faces = 0;
  cs_int_t  *i_face_lst = NULL, *b_face_lst = NULL;
  cs_real_t  *i_face_normal = NULL, *b_face_normal = NULL;
  double  *working_array = NULL, *i_face_warping = NULL, *b_face_warping = NULL;
  fvm_lnum_t  *n_i_sub_elt_lst = NULL, *n_b_sub_elt_lst = NULL;
  fvm_gnum_t  n_g_i_warp_faces = 0, n_g_b_warp_faces = 0;

  const cs_int_t  n_init_i_faces = mesh->n_i_faces;
  const cs_int_t  n_init_b_faces = mesh->n_b_faces;

#if 0   /* JB DEBUG */
  cs_mesh_dump(mesh);
#endif

  bft_printf(_("\n\n Découpage des faces gauche demandé\n"
               " ----------------------------------\n\n"
               " Angle maximal autorisé (rad) :\t%7.4f\n\n"), max_warp_angle);

  /* Compute face warping */

  BFT_MALLOC(working_array, n_init_i_faces + n_init_b_faces, double);

  for (i = 0; i < n_init_i_faces + n_init_b_faces; i++)
    working_array[i] = 0.;

  i_face_warping = working_array;
  b_face_warping = working_array + n_init_i_faces;

  cs_mesh_quantities_face_normal(mesh,
                                 &i_face_normal,
                                 &b_face_normal);

  cs_mesh_quality_compute_warping(mesh,
                                  i_face_normal,
                                  b_face_normal,
                                  i_face_warping,
                                  b_face_warping);

  BFT_FREE(i_face_normal);
  BFT_FREE(b_face_normal);

  _select_warped_faces(n_init_i_faces,
                       max_warp_angle,
                       i_face_warping,
                       &n_i_warp_faces,
                       &i_face_lst);

  _select_warped_faces(n_init_b_faces,
                       max_warp_angle,
                       b_face_warping,
                       &n_b_warp_faces,
                       &b_face_lst);

  /* Define the global number of faces which need to be cut */

  if (mesh->n_domains > 1) {
#if defined(_CS_HAVE_MPI)
    MPI_Allreduce(&n_i_warp_faces, &n_g_i_warp_faces, 1, CS_MPI_INT,
                  MPI_SUM, cs_glob_base_mpi_comm);

    MPI_Allreduce(&n_b_warp_faces, &n_g_b_warp_faces, 1, CS_MPI_INT,
                  MPI_SUM, cs_glob_base_mpi_comm);
#endif
  }
  else {
    n_g_i_warp_faces = n_i_warp_faces;
    n_g_b_warp_faces = n_b_warp_faces;
  }

  /* Test if there are faces to cut to continue */

  if (n_g_i_warp_faces == 0 && n_g_b_warp_faces == 0) {

    BFT_FREE(i_face_lst);
    BFT_FREE(b_face_lst);
    BFT_FREE(working_array);

    bft_printf(_("\n Aucune face à découper. "
                 " Revoir le critère si nécessaire.\n"));
    return;

  }

  /* Post-processing management */

  if (post_tag == true)
    _post_before_cutting(n_i_warp_faces,
                         n_b_warp_faces,
                         i_face_lst,
                         b_face_lst,
                         i_face_warping,
                         b_face_warping);

  /* Internal face treatment */
  /* ----------------------- */

  if (mesh->n_init_perio == 0)
    _cut_warped_faces(mesh,
                      2,
                      max_warp_angle,
                      i_face_warping,
                      &n_i_cut_faces,
                      &i_face_lst,
                      &n_i_sub_elt_lst,
                      &mesh->n_i_faces,
                      &mesh->i_face_vtx_connect_size,
                      &mesh->i_face_cells,
                      &mesh->i_face_family,
                      &mesh->i_face_vtx_idx,
                      &mesh->i_face_vtx_lst);

  else
    _cut_warped_faces_perio(mesh,
                            max_warp_angle,
                            i_face_warping,
                            &n_i_cut_faces,
                            &i_face_lst,
                            &n_i_sub_elt_lst,
                            &mesh->n_i_faces,
                            &mesh->i_face_vtx_connect_size,
                            &mesh->i_face_cells,
                            &mesh->i_face_family,
                            &mesh->i_face_vtx_idx,
                            &mesh->i_face_vtx_lst);

  bft_printf(_(" Faces internes:\n"
               "\t%12d faces découpées pour respecter le critère\n\n"
               "\t%12d faces localement avant découpage\n"
               "\t%12d faces localement après découpage\n\n"),
             n_i_warp_faces, n_init_i_faces, mesh->n_i_faces);

  bft_printf(_(" Taille de la nouvelle connectivité face -> sommets :%12d\n\n"),
             mesh->i_face_vtx_connect_size);

  /* Update global number of internal faces and its global numbering */

  _update_cut_faces_num(mesh,
                        mesh->n_i_faces,
                        n_init_i_faces,
                        n_i_sub_elt_lst,
                        &(mesh->n_g_i_faces),
                        &(mesh->global_i_face_num));

  /* Partial memory free */

  BFT_FREE(n_i_sub_elt_lst);

  /* Border face treatment */
  /* --------------------- */

  _cut_warped_faces(mesh,
                    1,
                    max_warp_angle,
                    b_face_warping,
                    &n_b_cut_faces,
                    &b_face_lst,
                    &n_b_sub_elt_lst,
                    &mesh->n_b_faces,
                    &mesh->b_face_vtx_connect_size,
                    &mesh->b_face_cells,
                    &mesh->b_face_family,
                    &mesh->b_face_vtx_idx,
                    &mesh->b_face_vtx_lst);

  bft_printf(_(" Faces de bord:\n"
               "\t%12d faces découpées pour respecter le critère\n\n"
               "\t%12d faces localement avant découpage\n"
               "\t%12d faces localement après découpage\n\n"),
             n_b_warp_faces, n_init_b_faces, mesh->n_b_faces);

  bft_printf(_(" Taille de la nouvelle connectivité face -> sommets :%12d\n\n"),
             mesh->b_face_vtx_connect_size);

  /* Update global number of border faces and its global numbering */

  _update_cut_faces_num(mesh,
                        mesh->n_b_faces,
                        n_init_b_faces,
                        n_b_sub_elt_lst,
                        &(mesh->n_g_b_faces),
                        &(mesh->global_b_face_num));

  /* Partial memory free */

  BFT_FREE(n_b_sub_elt_lst);

  /* Post-treatment of the selected faces */

  if (post_tag == true)
    _post_after_cutting(n_i_cut_faces,
                        n_b_cut_faces,
                        i_face_lst,
                        b_face_lst);

  /* Free memory */

  BFT_FREE(working_array);
  BFT_FREE(i_face_lst);
  BFT_FREE(b_face_lst);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
