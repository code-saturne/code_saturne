/*============================================================================
 * Functions dealing with ghost cells
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_periodicity.h"
#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_ctwr.h"
#include "cs_interface.h"
#include "cs_order.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

typedef struct _vtx_lookup_table {

  cs_int_t   n_vertices;    /* Number of local vertices in the problem domain */
  cs_int_t   n_interfaces;  /* Number of interfaces */

  cs_int_t   *if_ranks;     /* List of ranks */
  cs_int_t   *rank_ids;     /* list of rank ids */

  cs_int_t   *index;        /* index on table (size = n_vertices + 1) */

  cs_int_t   *rank_list;    /* list of ranks on which vertices are linked */

} vtx_lookup_table_t;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*---------------------------------------------------------------------------
 * Fill a look up table structure without periodicity
 *
 * The look up list takes the distant rank value with which a local vertex
 * is linked.
 *
 * parameters:
 *   vtx_lookup  -->  pointer to a vtx_lookup_table_t structure
 *   ifs         -->  pointer to a cs_interface_set_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_vtx_lookup(vtx_lookup_table_t  *vtx_lookup,
                 cs_interface_set_t  *ifs)
{
  cs_int_t  i, vtx_id, rank_id, shift;

  cs_int_t  *counter = NULL;

  const cs_int_t  n_interfaces = cs_interface_set_size(ifs);

  BFT_MALLOC(counter, vtx_lookup->n_vertices, cs_int_t);

  for (i = 0; i < vtx_lookup->n_vertices; i++)
    counter[i] = 0;

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    const cs_interface_t  *interface = cs_interface_set_get(ifs, rank_id);
    const cs_lnum_t  interface_size = cs_interface_size(interface);
    const cs_lnum_t  *local_id = cs_interface_get_elt_ids(interface);

    for (i = 0; i < interface_size; i++) { /* Only parallel vertices */

      vtx_id = local_id[i];
      shift = vtx_lookup->index[vtx_id] + counter[vtx_id];

      vtx_lookup->rank_list[shift] = vtx_lookup->rank_ids[rank_id];
      counter[vtx_id] += 1;

    }

  } /* End of loop on ranks */

  BFT_FREE(counter);

}

/*---------------------------------------------------------------------------
 * Create a vtx_look_up_table_t structure
 *
 * parameters:
 *   n_vertices  -->  number of vertices of the table.
 *   ifs         -->  pointer to a cs_interface_set_t structure
 *
 * returns:
 *   A pointer to the created vtx_lookup_table_t structure
 *---------------------------------------------------------------------------*/

static vtx_lookup_table_t *
_vtx_lookup_create(cs_int_t             n_vertices,
                   cs_interface_set_t  *ifs)
{
  cs_int_t  i, rank_id, tmp_id, interface_size;

  cs_int_t  loc_rank_id = -1;
  vtx_lookup_table_t  *vtx_lookup = NULL;

  const cs_interface_t  *interface = NULL;
  const cs_lnum_t  *local_id = NULL;
  const cs_int_t  n_interfaces = cs_interface_set_size(ifs);

  BFT_MALLOC(vtx_lookup, 1, vtx_lookup_table_t);

  vtx_lookup->n_vertices = n_vertices;
  vtx_lookup->n_interfaces = n_interfaces;

  BFT_MALLOC(vtx_lookup->index, n_vertices + 1, cs_int_t);
  BFT_MALLOC(vtx_lookup->if_ranks, n_interfaces, cs_int_t);
  BFT_MALLOC(vtx_lookup->rank_ids, n_interfaces, cs_int_t);

  for (i = 0; i < n_vertices + 1; i++)
    vtx_lookup->index[i] = 0;

  /* Check if cs_glob_rank_id belongs to the interface set in order to
     arrange if_ranks with local rank at first place */

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    interface = cs_interface_set_get(ifs, rank_id);
    vtx_lookup->if_ranks[rank_id] = cs_interface_rank(interface);
    vtx_lookup->rank_ids[rank_id] = rank_id;

    if (cs_glob_rank_id == cs_interface_rank(interface))
      loc_rank_id = rank_id;

  } /* End of loop on if_ranks */

  /* Define vtx_lookup->if_ranks in right order */

  if (loc_rank_id > 0) {

    tmp_id = vtx_lookup->if_ranks[loc_rank_id];
    vtx_lookup->if_ranks[loc_rank_id] = vtx_lookup->if_ranks[0];
    vtx_lookup->if_ranks[0] = tmp_id;

    vtx_lookup->rank_ids[0] = loc_rank_id;
    vtx_lookup->rank_ids[loc_rank_id] = 0;

    /* Order by increasing numbers */

    if (n_interfaces > 2) {

      cs_lnum_t  *order = NULL;
      cs_gnum_t  *buffer = NULL;
      cs_int_t  *_rank_ids = NULL;

      assert(sizeof(cs_lnum_t) == sizeof(cs_int_t));

      BFT_MALLOC(order, n_interfaces - 1, cs_lnum_t);
      BFT_MALLOC(buffer, n_interfaces - 1, cs_gnum_t);
      BFT_MALLOC(_rank_ids, n_interfaces, cs_int_t);

      _rank_ids[0] = vtx_lookup->rank_ids[0];
      for (i = 1; i < n_interfaces; i++) {
        buffer[i-1] = (cs_gnum_t)vtx_lookup->if_ranks[i];
        _rank_ids[i] = vtx_lookup->rank_ids[i];
      }

      cs_order_gnum_allocated(NULL,
                              buffer,
                              order,
                              n_interfaces-1);

      for (i = 0; i < n_interfaces - 1; i++) {
        vtx_lookup->if_ranks[i+1] = (cs_int_t)buffer[order[i]];
        vtx_lookup->rank_ids[i+1] = _rank_ids[order[i] + 1];
      }

      BFT_FREE(buffer);
      BFT_FREE(order);
      BFT_FREE(_rank_ids);

    } /* End of ordering ranks */

  } /* If rank order has to be changed */

  /* First loop to create index */

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    interface = cs_interface_set_get(ifs, rank_id);
    interface_size = cs_interface_size(interface);
    local_id = cs_interface_get_elt_ids(interface);

    for (i = 0; i < interface_size; i++)
      vtx_lookup->index[local_id[i] + 1] += 1;

  } /* End of loop on if_ranks */

  /* Create index and allocate buffers */

  for (i = 0; i < n_vertices; i++)
    vtx_lookup->index[i+1] += vtx_lookup->index[i];

  BFT_MALLOC(vtx_lookup->rank_list, vtx_lookup->index[n_vertices], cs_int_t);

  /* Second loop to fill table(s) */

  _fill_vtx_lookup(vtx_lookup, ifs);

  return vtx_lookup;
}

/*---------------------------------------------------------------------------
 * Destroy a vtx_lookup structure.
 *
 * parameters:
 *   vtx_lookup -->  pointer to a vtx_lookup_table_t structure
 *
 * returns:
 *   A NULL pointer
 *---------------------------------------------------------------------------*/

static void
_vtx_lookup_destroy(vtx_lookup_table_t  *vtx_lookup)
{

  BFT_FREE(vtx_lookup->if_ranks);
  BFT_FREE(vtx_lookup->rank_ids);
  BFT_FREE(vtx_lookup->index);
  BFT_FREE(vtx_lookup->rank_list);

  BFT_FREE(vtx_lookup);

}

/*---------------------------------------------------------------------------
 * Set checker for this vertex_id according to vtx_lookup features.
 *
 * parameters:
 *   vtx_id       -->  vertex id to deal with
 *   vtx_checker  <->  put a tag in the implied categories
 *   vtx_lookup   -->  pointer to a vtx_lookup_table_t structure
 *---------------------------------------------------------------------------*/

static void
_update_vtx_checker(cs_int_t             vtx_id,
                    cs_int_t            *vtx_checker,
                    vtx_lookup_table_t  *vtx_lookup)
{
  cs_int_t  i, rank_id;

  for (i = vtx_lookup->index[vtx_id];
       i < vtx_lookup->index[vtx_id + 1];
       i++) {
    rank_id = vtx_lookup->rank_list[i];
    vtx_checker[rank_id] += 1;
  } /* End of loop on vtx_lookup */
}

/*---------------------------------------------------------------------------
 *
 * TODO  Build in halo with extrusion
 *
 * parameters:
 *   ct  <-> pointer to a ct structure
 *---------------------------------------------------------------------------*/

static void
_build_halo_with_extrusion(cs_ctwr_zone_t  *ct)
{
  cs_int_t        i,  idx, counter, j;
  cs_int_t       *list_tmp    = NULL;
  cs_halo_t   *halo     = ct->water_halo;
  const cs_int_t  n_c_domains = halo->n_c_domains;

  BFT_MALLOC(list_tmp, halo->n_send_elts[0], cs_int_t);

  for (i = 0; i < halo->n_send_elts[0]; i++)
    list_tmp[i] = halo->send_list[i];

  BFT_REALLOC(halo->send_list, halo->n_send_elts[0]*ct->nelect, cs_int_t);

  counter = 0;
  for (i = 0; i < n_c_domains; i++) {

    for (idx = halo->send_index[i];
         idx < halo->send_index[i + 1];
         idx++) {
      for (j = 0; j < ct->nelect; j++) {
        halo->send_list[counter]  = list_tmp[idx]*ct->nelect + j;
        counter++;
      }
    }
    halo->send_index[i] *= ct->nelect;
  }

  halo->send_index[n_c_domains] *= ct->nelect;

  halo->n_send_elts[0] *= ct->nelect;
  halo->n_send_elts[1] = halo->n_send_elts[0];

  BFT_FREE(list_tmp);
}

/*---------------------------------------------------------------------------
 * Define the elements of send_halo structure.
 *
 * Two main loops. First one for counting number of elements and create index.
 * Second one for filling the ghost cells list.
 *
 * parameters:
 *   mesh           --> pointer to cs_mesh_t structure
 *   interface_set  --> pointer to cs_interface_set_t structure
 *  vtx_faces_idx   --> "vtx -> faces" connectivity index
 *  vtx_faces_lst   --> "vtx  -> faces" connectivity list
 *---------------------------------------------------------------------------*/

static void
_fill_send_halo(cs_ctwr_zone_t      *ct,
                cs_interface_set_t  *interface_set,
                cs_int_t            *vtx_faces_idx,
                cs_int_t            *vtx_faces_lst)
{
  cs_int_t  i, fac_idx, shift, rank_id;
  cs_int_t  vtx_id, n_interfaces;

  cs_halo_t  *halo = ct->water_halo;
  vtx_lookup_table_t *vtx_lookup = NULL;

  cs_int_t  *vtx_checker  = NULL;
  cs_int_t  *counter      = NULL;
  cs_int_t  *faces_tag    = NULL;

  const cs_int_t  n_vertices = fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);
  const cs_int_t  n_faces    = fvm_nodal_get_n_entities(ct->face_sup_mesh, 2);


  /* We should have the faces -> vertices connectivity to continue */

  if (vtx_faces_lst == NULL)
    return;

  /* Create a lookup table to accelerate search in
     cs_interface_set structure */

  vtx_lookup   = _vtx_lookup_create(n_vertices, interface_set);
  n_interfaces =  vtx_lookup->n_interfaces;

  BFT_MALLOC(vtx_checker, n_interfaces, cs_int_t);
  BFT_MALLOC(faces_tag, n_faces * n_interfaces, cs_int_t);

  for (fac_idx = 0; fac_idx < n_faces * n_interfaces; fac_idx++)
      faces_tag[fac_idx] = 0;

  /* First loop to create index and allocate  send_list */

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++) {

    for (i = 0; i < n_interfaces; i++)
      vtx_checker[i] = 0;

    _update_vtx_checker(vtx_id, vtx_checker, vtx_lookup);

    for (fac_idx = vtx_faces_idx[vtx_id];
         fac_idx < vtx_faces_idx[vtx_id + 1];
         fac_idx++) {

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
        if (vtx_checker[rank_id] == 1) {
          if (faces_tag[n_faces*rank_id + vtx_faces_lst[fac_idx]] != 1) {
            halo->send_index[rank_id + 1] += 1;
            faces_tag[n_faces*rank_id + vtx_faces_lst[fac_idx]] = 1;
          }
        }
      }
    }
  } /* End of loop on vertices */

  for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++)
    halo->send_index[rank_id + 1] += halo->send_index[rank_id];

  BFT_MALLOC(halo->send_list, halo->send_index[halo->n_c_domains], cs_int_t);

  /* Initialize counter */

  BFT_MALLOC(counter, n_interfaces, cs_int_t);

  for (i = 0; i < n_interfaces; i++)
    counter[i] = 0;

  /* Second loop to build halo->ghost_cells */

  for (fac_idx = 0; fac_idx < n_faces * n_interfaces; fac_idx++)
    faces_tag[fac_idx] = 0;

  for (vtx_id = 0; vtx_id <n_vertices; vtx_id++) {

    for (i = 0; i < n_interfaces; i++)
      vtx_checker[i] = 0;

    _update_vtx_checker(vtx_id, vtx_checker, vtx_lookup);

    for (fac_idx = vtx_faces_idx[vtx_id];
         fac_idx < vtx_faces_idx[vtx_id + 1];
         fac_idx++) {

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
        if (vtx_checker[rank_id] == 1) {
          if (faces_tag[n_faces*rank_id + vtx_faces_lst[fac_idx]] != 1) {

            shift = halo->send_index[rank_id] + counter[rank_id];
            halo->send_list[shift] = vtx_faces_lst[fac_idx];
            counter[rank_id] += 1;

            faces_tag[n_faces*rank_id + vtx_faces_lst[fac_idx]] = 1;

          }
        }
      }
    }
  }

  /* Free memory */

  BFT_FREE(vtx_checker);
  BFT_FREE(counter);

  /* Destroy the lookup table strcuture */

  _vtx_lookup_destroy(vtx_lookup);

  /* Complete halo definition */

  halo->n_send_elts[CS_HALO_STANDARD] = 0;
  halo->n_send_elts[CS_HALO_EXTENDED] = 0;

  for (i = 0; i < halo->n_c_domains; i++)
    halo->n_send_elts[CS_HALO_STANDARD] += halo->send_index[i+1]
                                         - halo->send_index[i];

  halo->n_send_elts[CS_HALO_EXTENDED] += halo->n_send_elts[CS_HALO_STANDARD];
}

/*---------------------------------------------------------------------------
 * Define a buffer on vertices where vertex belonging to the interface_set
 * are tagged with 1 else 0.
 *
 * parameters:
 *   n_vertices    --> size of the buffer
 *   interface_set --> pointer to a cs_interface_set_t structure
 *   p_vertex_tag  <-> pointer to the tagged buffer
 *---------------------------------------------------------------------------*/

static void
_get_vertex_tag(cs_int_t                   n_vertices,
                const cs_interface_set_t  *interface_set,
                cs_int_t                  *p_vertex_tag[])
{
  cs_int_t  i, j, rank_id;

  cs_int_t  *vertex_tag = NULL;

  const int  ifs_size = cs_interface_set_size(interface_set);

  BFT_MALLOC(vertex_tag, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    vertex_tag[i] = 0;

  for (rank_id = 0; rank_id < ifs_size; rank_id++) {

    const cs_interface_t  *interface = cs_interface_set_get(interface_set,
                                                            rank_id);
    const cs_lnum_t  *local_id = cs_interface_get_elt_ids(interface);
    const cs_lnum_t  if_size = cs_interface_size(interface);

    for (j = 0; j < if_size; j++)
      vertex_tag[local_id[j]] = 1;

  } /* End of loop on ranks */

  *p_vertex_tag = vertex_tag;
}

/*---------------------------------------------------------------------------
 * Exchange number and list of cells constituting send_halo structure for each
 * frontier ranks. Fill the halo structure from these data.
 *
 * parameters:
 *   ct  --> pointer to a ct structure
 *---------------------------------------------------------------------------*/

static void
_fill_halo(cs_ctwr_zone_t  *ct)
{
  cs_int_t  rank_id, i;
  cs_int_t  shift;

#if defined(HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  int request_count = 0;

  cs_int_t  *count = NULL;

  cs_halo_t  *halo = ct->water_halo;

  const  cs_int_t  n_c_domains = halo->n_c_domains;
  const  cs_int_t  local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;

#if defined(HAVE_MPI)
  if (halo->n_c_domains*2 > 128) {
    BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);
  }
#endif

  /* Build index */
  /* ----------- */
  /* Receive data from distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(HAVE_MPI)
      MPI_Irecv(&(halo->index[rank_id+1]), 1, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                halo->c_domain_rank[rank_id],
                cs_glob_mpi_comm,
                &(request[request_count++]));
#endif

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Barrier(cs_glob_mpi_comm);
#endif

  BFT_MALLOC(count, n_c_domains, cs_int_t);

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = rank_id;
    count[shift] =   halo->send_index[rank_id+1]
                   - halo->send_index[rank_id];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(HAVE_MPI)
      MPI_Isend(&(count[shift]), 1, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                local_rank,
                cs_glob_mpi_comm,
                &(request[request_count++]));
#endif

    }
    else {

      halo->index[shift+1] = count[shift];

    }

  } /* End of loop on ranks */

  /* Wait for all exchanges being done */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Waitall(request_count, request, status);
#endif
  request_count = 0;

  BFT_FREE(count);

  /* Build index */
  /*-------------*/

  for (i = 0; i < n_c_domains; i++)
    halo->index[i+1] += halo->index[i];

#if defined(HAVE_MPI)
  if (request != _request) {
    BFT_FREE(request);
    BFT_FREE(status);
  }
#endif

  halo->n_elts[CS_HALO_STANDARD] = 0;
  halo->n_elts[CS_HALO_EXTENDED] = 0;

  for (i = 0; i < n_c_domains; i++) {

    halo->n_elts[CS_HALO_STANDARD] += halo->index[i+1] - halo->index[i];

  }

  halo->n_elts[CS_HALO_EXTENDED] += halo->n_elts[CS_HALO_STANDARD];
}

/*---------------------------------------------------------------------------
 *  TODO
 *---------------------------------------------------------------------------*/

static void
_fill_index_out_halo(cs_ctwr_zone_t  *ct)
{
  cs_int_t  rank_id, i;
  cs_int_t  shift;

#if defined(HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  int request_count = 0;

  cs_int_t  *count = NULL;

  cs_halo_t  *halo = ct->water_halo;

  const cs_int_t n_c_domains = halo->n_c_domains;
  const cs_int_t local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;

#if defined(HAVE_MPI)
  if (halo->n_c_domains*2 > 128) {
    BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);
  }
#endif

  /* Build index */
  /* ----------- */
  /* Receive data from distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(HAVE_MPI)
      MPI_Irecv(&(halo->index[rank_id+1]), 1, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                halo->c_domain_rank[rank_id],
                cs_glob_mpi_comm,
                &(request[request_count++]));
#endif

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Barrier(cs_glob_mpi_comm);
#endif

  BFT_MALLOC(count, n_c_domains, cs_int_t);

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = rank_id;
    count[shift] =   halo->send_index[rank_id+1]
                   - halo->send_index[rank_id];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(HAVE_MPI)
      MPI_Isend(&(count[shift]), 1, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                local_rank,
                cs_glob_mpi_comm,
                &(request[request_count++]));
#endif

    }
    else {

      halo->index[shift+1] = count[shift];

    }

  } /* End of loop on ranks */

  /* Wait for all exchanges being done */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Waitall(request_count, request, status);
#endif
  request_count = 0;

  BFT_FREE(count);

  /* Build index */
  /*-------------*/

  for (i = 0; i < n_c_domains; i++)
    halo->index[i+1] += halo->index[i];

  /* Wait for all exchanges being done */

#if defined(HAVE_MPI)
  if (request_count > 0)
    MPI_Waitall(request_count, request, status);
#endif
  request_count = 0;

  /* Exchange number of elements for each periodicity and for each rank.
     Then build out_halo->perio_lst */

  halo->n_elts[0] = 0;

  for (i = 0; i < n_c_domains; i++)
    halo->n_elts[0] += halo->index[i+1] - halo->index[i];

  halo->n_elts[1] = halo->n_elts[0];
}

/*---------------------------------------------------------------------------
 * Send "ghost cells to distant_num vertices" connectivity on communicating
 * ranks and receive the same kind of connectivity from distant ranks.
 *
 * parameters:
 *   ct                      --> pointer to ct structure
 *   send_gface_dist_vtx_idx <-- "ghost face -> distant vertices" index
 *   send_gface_dist_vtx_lst <-- "ghost face -> distant vertices" list
 *   p_gface_dist_vtx_idx    --> "ghost face -> distant vertices" index
 *   p_gface_dist_vtx_lst    --> "ghost face -> distant vertices" list
 *---------------------------------------------------------------------------*/

static void
_exchange_gface_vtx_connect(cs_ctwr_zone_t *ct,
                            cs_int_t   *send_gface_dist_vtx_idx,
                            cs_int_t   *send_gface_dist_vtx_lst,
                            cs_int_t   *p_gface_dist_vtx_idx[],
                            cs_int_t   *p_gface_dist_vtx_lst[])
{
  cs_int_t  i, j, rank_id;
  cs_int_t  send_start_idx, send_end_idx, start_idx, end_idx;
  cs_int_t  n_send_elts, n_recv_elts;

  cs_int_t  send_buffer_size = 0;

  cs_int_t  *send_idx_buffer   = NULL;
  cs_int_t  *gface_dist_vtx_idx = NULL, *gface_dist_vtx_lst = NULL;
  cs_int_t  *send_buffer = NULL, *recv_buffer = NULL;

  cs_halo_t  *halo = ct->water_halo;

  const cs_int_t  local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;
  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  n_ghost_faces = halo->n_elts[CS_HALO_EXTENDED];

#if defined(HAVE_MPI)
  MPI_Status  status;
#endif

  /* Allocate buffers */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {
    if (halo->c_domain_rank[rank_id] != local_rank) {
      n_send_elts = halo->send_index[rank_id + 1]- halo->send_index[rank_id];
      send_buffer_size = CS_MAX(send_buffer_size, n_send_elts);
    }
  }

  BFT_MALLOC(send_idx_buffer, send_buffer_size, cs_int_t);
  BFT_MALLOC(gface_dist_vtx_idx, n_ghost_faces + 1, cs_int_t);

  for (i = 0; i < n_ghost_faces + 1; i++)
    gface_dist_vtx_idx[i] = 0;

  /* Exchange sizes to define gface_dist_vtx_idx */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    recv_buffer = &(gface_dist_vtx_idx[1 + halo->index[rank_id]]);

    if (halo->c_domain_rank[rank_id] != local_rank) {

      /* Fill send buffer */

      for (i = halo->send_index[rank_id], j = 0;
           i < halo->send_index[rank_id + 1]; i++, j++)
        send_idx_buffer[j] =   send_gface_dist_vtx_idx[i+1]
                             - send_gface_dist_vtx_idx[i];

      n_send_elts =  halo->send_index[rank_id + 1] - halo->send_index[rank_id];
      n_recv_elts =  halo->index[rank_id + 1] - halo->index[rank_id];

#if defined(HAVE_MPI)
      MPI_Sendrecv(&(send_idx_buffer[0]), n_send_elts, CS_MPI_INT,
                     halo->c_domain_rank[rank_id], local_rank,
                     &(recv_buffer[0]), n_recv_elts, CS_MPI_INT,
                     halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                     cs_glob_mpi_comm, &status);
#endif

    } /* If rank != local_rank */

    else {

      for (i = halo->send_index[rank_id], j = 0;
           i < halo->send_index[rank_id + 1]; i++, j++)
        recv_buffer[j] =   send_gface_dist_vtx_idx[i+1]
                         - send_gface_dist_vtx_idx[i];

    } /* rank == local_rank */

  } /* End of loop on if_ranks */

  BFT_FREE(send_idx_buffer);

  /* Define index */

  for (i = 0; i < n_ghost_faces; i++)
    gface_dist_vtx_idx[i+1] += gface_dist_vtx_idx[i];

  BFT_MALLOC(gface_dist_vtx_lst, gface_dist_vtx_idx[n_ghost_faces], cs_int_t);

  /* Exchange lists to define gface_dist_vtx_lst */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    /* Exchange conectivity list */

    send_start_idx = send_gface_dist_vtx_idx[halo->send_index[rank_id]];
    send_end_idx   = send_gface_dist_vtx_idx[halo->send_index[rank_id + 1]];
    n_send_elts  = send_end_idx - send_start_idx;
    send_buffer  = &(send_gface_dist_vtx_lst[send_start_idx]);

    start_idx = gface_dist_vtx_idx[halo->index[rank_id]];
    end_idx   = gface_dist_vtx_idx[halo->index[rank_id + 1]];
    n_recv_elts   = end_idx - start_idx;
    recv_buffer   = &(gface_dist_vtx_lst[start_idx]);

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(HAVE_MPI)
      MPI_Sendrecv(send_buffer, n_send_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], local_rank,
                   recv_buffer, n_recv_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                   cs_glob_mpi_comm, &status);
#endif

    }
    else {

      assert(n_recv_elts == n_send_elts);

      for (i = 0; i < n_send_elts; i++)
        recv_buffer[i] = send_buffer[i];

    }

  } /* End of loop on ranks */

  *p_gface_dist_vtx_idx = gface_dist_vtx_idx;
  *p_gface_dist_vtx_lst = gface_dist_vtx_lst;

}

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

static void
_create_gvtx_faces_connect(cs_ctwr_zone_t      *ct,
                           cs_interface_set_t  *interface_set,
                           cs_int_t            *vtx_faces_idx[],
                           cs_int_t            *vtx_faces_lst[])
{
  cs_int_t   id, nb, index, index_g;

  cs_int_t    *vtx_tag         = NULL;
  cs_lnum_t   *_g_vtx_faces_idx = NULL;
  cs_lnum_t   *_g_vtx_faces_lst = NULL;
  cs_int_t    *_vtx_faces_idx  = NULL;
  cs_int_t    *_vtx_faces_lst  = NULL;

  const cs_int_t  n_vtx
    = (cs_int_t) fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);

  fvm_nodal_get_vertex_elements(ct->face_sup_mesh,
                                2,
                                &(_g_vtx_faces_idx),
                                &(_g_vtx_faces_lst));
  *vtx_faces_idx = NULL;
  *vtx_faces_lst = NULL;

  if (interface_set == NULL)
    return;

  _get_vertex_tag(n_vtx, interface_set, &vtx_tag);

  BFT_MALLOC(_vtx_faces_idx,  n_vtx + 1, cs_int_t);

  _vtx_faces_idx[0] = 0;

  for (id = 0;  id < n_vtx; id++) {

    if (vtx_tag[id] == 1) {
      nb = (cs_int_t) _g_vtx_faces_idx[id + 1] - _g_vtx_faces_idx[id];
      _vtx_faces_idx[id + 1] = _vtx_faces_idx[id] + nb;

    }
    else
      _vtx_faces_idx[id + 1] = _vtx_faces_idx[id];

  }

  BFT_MALLOC( _vtx_faces_lst,  _vtx_faces_idx[n_vtx], cs_int_t);

  if (_vtx_faces_idx[n_vtx] == 0)
    return;

  index = 0;

  for (id = 0; id < n_vtx; id++) {
    if (vtx_tag[id] == 1) {
      for (index_g = _g_vtx_faces_idx[id];
           index_g < _g_vtx_faces_idx[id + 1];
           index_g++) {
        _vtx_faces_lst[index] = (cs_int_t) _g_vtx_faces_lst[index_g];
        index++;
      }
    }
  }

  *vtx_faces_idx = _vtx_faces_idx;
  *vtx_faces_lst = _vtx_faces_lst;

  BFT_FREE(_g_vtx_faces_idx);
  BFT_FREE(_g_vtx_faces_lst);
  BFT_FREE(vtx_tag);
}

/*---------------------------------------------------------------------------
 * Create a local "ghost  faces -> distant vertices" connectivity for
 * in_halo cells.
 *
 * parameters:
 *   mesh                 --> pointer to cs_mesh_t structure
 *   interface_set        --> pointer to cs_interface_set_t structure
 *   g_faces_vtx_idx     --> "faces -> vertices" connectivity index
 *   g_faces_vtx_lst     --> "faces -> vertices" connectivity list
 *   p_g_in_faces_vtx_idx<-- "ghost faces -> distant vertices" connect. index
 *   p_g_in_faces_vtx_lst<-- "ghost faces -> distant vertices" connect. list
 *---------------------------------------------------------------------------*/

static void
_create_in_faces_vtx_connect(cs_ctwr_zone_t      *ct,
                             cs_interface_set_t  *interface_set,
                             cs_int_t            *g_faces_vtx_idx,
                             cs_int_t            *g_faces_vtx_lst,
                             cs_int_t            *p_g_in_faces_vtx_idx[],
                             cs_int_t            *p_g_in_faces_vtx_lst[])
{
  cs_int_t  i, j, rank_id, nb_face_in, idfac,ivtx, shift;

  cs_halo_t  *halo = ct->water_halo;

  cs_int_t   *  g_in_faces_vtx_idx;
  cs_int_t   *  g_in_faces_vtx_lst;

  cs_int_t   *vertex_tag = NULL;
  cs_int_t   *counter    = NULL;

  if (interface_set == NULL)
    return;

  const cs_int_t  local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;
  const cs_int_t  n_vertices = fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);

  const int  ifs_size   = cs_interface_set_size(interface_set);

  cs_interface_set_add_match_ids(interface_set);

  nb_face_in  = halo->send_index[halo->n_c_domains];

  BFT_MALLOC(vertex_tag, n_vertices, cs_int_t);

  BFT_MALLOC(g_in_faces_vtx_idx, nb_face_in + 1, cs_int_t);

  /* first loop to count*/
  for (idfac = 0; idfac <  nb_face_in + 1; idfac++)
    g_in_faces_vtx_idx[idfac] = 0;

  for (rank_id = 0; rank_id < ifs_size; rank_id++) {

    for (i = 0; i < n_vertices; i++)
      vertex_tag[i] = -1;

    if (halo->c_domain_rank[rank_id] != local_rank) {
      const cs_interface_t *interface = cs_interface_set_get(interface_set,
                                                             rank_id);
      const cs_lnum_t  *local_id = cs_interface_get_elt_ids(interface);
      const cs_lnum_t  if_size    = cs_interface_size(interface);

      for (j = 0; j < if_size; j++)
        vertex_tag[local_id[j]] = 0;

      for (idfac = halo->send_index[rank_id];
           idfac < halo->send_index[rank_id + 1];
           idfac++) {
        for (ivtx = g_faces_vtx_idx[halo->send_list[idfac]];
             ivtx < g_faces_vtx_idx[halo->send_list[idfac] + 1];
             ivtx++) {
          if (vertex_tag[g_faces_vtx_lst[ivtx]] != -1)
            g_in_faces_vtx_idx[idfac + 1]++;
        }
      }

    }
  } /* End of loop on ranks */

  BFT_MALLOC(counter,  nb_face_in, cs_int_t);

  for(idfac = 0; idfac < nb_face_in; idfac++) {
    g_in_faces_vtx_idx[idfac + 1] += g_in_faces_vtx_idx[idfac];
    counter[idfac] = 0;
  }

  /* second loop to */
  BFT_MALLOC(g_in_faces_vtx_lst, g_in_faces_vtx_idx[nb_face_in], cs_int_t);

  for (rank_id = 0; rank_id < ifs_size; rank_id++) {

    for (i = 0; i < n_vertices; i++)
      vertex_tag[i] = -1;

    if (halo->c_domain_rank[rank_id] != local_rank) {
      const cs_interface_t *interface = cs_interface_set_get(interface_set,
                                                             rank_id);
      const cs_lnum_t  *local_id   = cs_interface_get_elt_ids(interface);
      const cs_lnum_t  *distant_id = cs_interface_get_match_ids(interface);
      const cs_lnum_t  if_size     = cs_interface_size(interface);

      for (j = 0; j < if_size; j++)
        vertex_tag[local_id[j]] = distant_id[j];


      for(idfac = halo->send_index[rank_id];
          idfac < halo->send_index[rank_id  + 1];
          idfac++) {
        for(ivtx = g_faces_vtx_idx[halo->send_list[idfac]];
            ivtx < g_faces_vtx_idx[halo->send_list[idfac] + 1];
            ivtx++) {
          if (vertex_tag[g_faces_vtx_lst[ivtx]] != -1) {
            shift = g_in_faces_vtx_idx[idfac] + counter[idfac];
            g_in_faces_vtx_lst[shift] = vertex_tag[g_faces_vtx_lst[ivtx]];
            counter[idfac]++;
          }
        }
      }
    }
  } /* End of loop on ranks */

  cs_interface_set_free_match_ids(interface_set);

  *p_g_in_faces_vtx_idx = g_in_faces_vtx_idx;
  *p_g_in_faces_vtx_lst = g_in_faces_vtx_lst;
}

/*---------------------------------------------------------------------------
 *
 *---------------------------------------------------------------------------*/

static void
_update_gfaces_connect(cs_ctwr_zone_t  *ct,
                       cs_int_t        *faces_vtx_idx,
                       cs_int_t        *faces_vtx_lst,
                       cs_int_t        *g_vtx_faces_idx[],
                       cs_int_t        *g_vtx_faces_lst[])
{
  cs_int_t  i, idx, shift;

  cs_halo_t  *halo = ct->water_halo;

  cs_int_t   *new_vtx_faces_idx = NULL;
  cs_int_t   *new_vtx_faces_lst = NULL;
  cs_int_t   *_g_vtx_faces_idx  = NULL;
  cs_int_t   *_g_vtx_faces_lst  = NULL;

  cs_int_t   *counter           = NULL;


  const cs_int_t  n_vertices   = fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);
  const cs_int_t  n_faces      = fvm_nodal_get_n_entities(ct->face_sup_mesh, 2);

  BFT_MALLOC(new_vtx_faces_idx, n_vertices + 1, cs_int_t);

  fvm_nodal_get_vertex_elements(ct->face_sup_mesh,
                                2,
                                &(_g_vtx_faces_idx),
                                &(_g_vtx_faces_lst));

  new_vtx_faces_idx[0]= 0;
  for (i = 0; i < n_vertices; i++)
    new_vtx_faces_idx[i + 1] = _g_vtx_faces_idx[i + 1] - _g_vtx_faces_idx[i];

  for (i = 0; i < halo->n_elts[0]; i++)
     for (idx =  faces_vtx_idx[i];
          idx <  faces_vtx_idx[i + 1]; idx++)
     {
        new_vtx_faces_idx[faces_vtx_lst[idx] + 1] +=1;
     }

  for (i = 0; i < n_vertices; i++) {
    new_vtx_faces_idx[i + 1] += new_vtx_faces_idx[i];
  }
  BFT_MALLOC(new_vtx_faces_lst, new_vtx_faces_idx[n_vertices], cs_int_t);

  BFT_MALLOC(counter, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    counter[i] = 0;

  for (i = 0; i < n_vertices; i++)
    for (idx =  _g_vtx_faces_idx[i];
         idx <  _g_vtx_faces_idx[i + 1]; idx++)
    {
      shift = new_vtx_faces_idx[i] + counter[i];
      new_vtx_faces_lst[shift] =  _g_vtx_faces_lst[idx];
      counter[i]++;
    }

  for (i = 0; i < halo->n_elts[0]; i++)
    for (idx =  faces_vtx_idx[i];
         idx <  faces_vtx_idx[i + 1]; idx++)
    {
        shift = new_vtx_faces_idx[faces_vtx_lst[idx]]
                + counter[faces_vtx_lst[idx]];

        new_vtx_faces_lst[shift] =   i + n_faces;

        counter[faces_vtx_lst[idx]]++;
    }


  *g_vtx_faces_idx = new_vtx_faces_idx;
  *g_vtx_faces_lst = new_vtx_faces_lst;

  BFT_FREE(_g_vtx_faces_idx);
  BFT_FREE(_g_vtx_faces_lst);

  BFT_FREE(counter);

}

/*---------------------------------------------------------------------------
 * st
 *---------------------------------------------------------------------------*/

static void
_create_faces_faces_connect(cs_ctwr_zone_t  *ct,
                            cs_int_t        *faces_vtx_idx,
                            cs_int_t        *faces_vtx_lst,
                            cs_int_t        *vtx_faces_idx,
                            cs_int_t        *vtx_faces_lst)
{
  cs_int_t  i, shift, iface, ivtx, ifext;

  cs_halo_t  *halo = ct->water_halo;

  cs_int_t   *counter = NULL;
  cs_int_t   *face_tag = NULL;

  const cs_int_t  n_faces = ct->nnpsct;
  const cs_int_t  n_faces_with_ghosts = n_faces +  halo->n_elts[0];

  BFT_MALLOC(ct->fac_sup_connect_idx, n_faces + 1, cs_int_t);
  BFT_MALLOC(face_tag, n_faces_with_ghosts, cs_int_t);

  if (vtx_faces_idx == NULL)
    fvm_nodal_get_vertex_elements(ct->face_sup_mesh,
                                 2,
                                 &(vtx_faces_idx),
                                 &(vtx_faces_lst));

  for (iface = 0; iface < n_faces + 1; iface++)
    ct->fac_sup_connect_idx[iface]= 0;

  for (iface = 0; iface < n_faces; iface++) {
    for (i  = 0; i  < n_faces_with_ghosts; i ++)
      face_tag[i] = 0;

    face_tag[iface] = -1;

    for (ivtx = faces_vtx_idx[iface];
         ivtx < faces_vtx_idx[iface + 1];
         ivtx++) {

      for (ifext = vtx_faces_idx[faces_vtx_lst[ivtx]];
           ifext < vtx_faces_idx[faces_vtx_lst[ivtx] + 1];
           ifext++) {
        if (face_tag[vtx_faces_lst[ifext]] != -1) {
          ct->fac_sup_connect_idx[iface + 1] ++;
          face_tag[vtx_faces_lst[ifext]] = -1;
        }
      }
    }

  }

  for (iface = 0; iface < n_faces; iface++)
    ct->fac_sup_connect_idx[iface + 1] += ct->fac_sup_connect_idx[iface];

  BFT_MALLOC(ct->fac_sup_connect_lst,
             ct->fac_sup_connect_idx[n_faces],
             cs_int_t);

  BFT_MALLOC(counter, n_faces, cs_int_t);

  for (i = 0; i < n_faces; i++)
    counter[i] = 0;

  for (iface = 0; iface < n_faces; iface++) {
    for (i  = 0; i  <  n_faces_with_ghosts; i ++)
      face_tag[i] = 0;

    face_tag[iface] = -1;

    for (ivtx = faces_vtx_idx[iface];
         ivtx < faces_vtx_idx[iface + 1];
         ivtx++) {
      for (ifext = vtx_faces_idx[faces_vtx_lst[ivtx]];
           ifext < vtx_faces_idx[faces_vtx_lst[ivtx] + 1];
           ifext++) {
        if (face_tag[vtx_faces_lst[ifext]] != -1) {
          shift = ct->fac_sup_connect_idx[iface] + counter[iface];

          ct->fac_sup_connect_lst[shift] = vtx_faces_lst[ifext];

          face_tag[vtx_faces_lst[ifext]] = -1;
          counter[iface]++;
        }
      }
    }

  }

  BFT_FREE(counter);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*---------------------------------------------------------------------------
 * Reverse "ghost cells -> vertex" connectivity into "vertex -> ghost cells"
 * connectivity for out_halo elements.
 * Build the connectivity list.
 *
 * parameters:
 *   halo            --> pointer to a cs_mesh_halo_t structure
 *   n_vertices      --> number of vertices
 *   rank_id         --> rank number to work with
 *   counter         <-> temporary array to count vertices
 *   checker         <-> temporary array to check vertices
 *   gcell_vtx_idx   --> "ghost cell -> vertices" connectivity index
 *   gcell_vtx_lst   --> "ghost cell -> vertices" connectivity list
 *   vtx_gcells_idx  --> "vertex -> ghost cells" connectivity index
 *   vtx_gcells_lst  <-> "vertex -> ghost cells" connectivity list
 *---------------------------------------------------------------------------*/

void
cs_reverse_vtx_faces_connect(const fvm_nodal_t  *this_nodal,
                             cs_int_t           *faces_vtx_idx[],
                             cs_int_t           *faces_vtx_lst[])
{
  cs_int_t iv, ifac, shift;

  cs_int_t   *_faces_vtx_idx = NULL;
  cs_int_t   *_faces_vtx_lst = NULL;
  cs_int_t   *_vtx_faces_idx = NULL;
  cs_int_t   *_vtx_faces_lst = NULL;

  cs_int_t   *counter          = NULL;
  const cs_int_t  n_vtx = (cs_int_t) fvm_nodal_get_n_entities(this_nodal, 0);
  const cs_int_t  n_fac = (cs_int_t) fvm_nodal_get_n_entities(this_nodal, 2);

  BFT_MALLOC(_faces_vtx_idx,  n_fac + 1, cs_int_t);
  BFT_MALLOC(counter, n_fac, cs_int_t);

  fvm_nodal_get_vertex_elements(this_nodal,
                                2,
                                &(_vtx_faces_idx),
                                &(_vtx_faces_lst));

  for (ifac = 0; ifac < n_fac + 1; ifac++)
    _faces_vtx_idx[ifac] = 0;

  for (ifac = 0; ifac < n_fac; ifac++)
    counter[ifac] = 0;

  for (iv = 0; iv <  n_vtx; iv++) {
    for (ifac = _vtx_faces_idx[iv];
          ifac < _vtx_faces_idx[iv + 1]; ifac++) {

           _faces_vtx_idx[_vtx_faces_lst[ifac] + 1] ++;

    }

  }

  for (ifac = 0; ifac < n_fac; ifac++) {
    _faces_vtx_idx[ifac + 1] +=  _faces_vtx_idx[ifac];
  }

  BFT_MALLOC(_faces_vtx_lst,  _faces_vtx_idx[n_fac], cs_int_t);

  for (iv = 0; iv <  n_vtx; iv++) {
    for (ifac = _vtx_faces_idx[iv];
          ifac < _vtx_faces_idx[iv + 1]; ifac++) {

            shift =  _faces_vtx_idx[_vtx_faces_lst[ifac]]
                    + counter[_vtx_faces_lst[ifac]];

           _faces_vtx_lst[shift] = iv;

           counter[_vtx_faces_lst[ifac]]++;

    }

  }

  BFT_FREE(counter);

  *faces_vtx_idx = _faces_vtx_idx;
  *faces_vtx_lst = _faces_vtx_lst;
}


/*----------------------------------------------------------------------------
 * Define halo structures for internal and distant ghost cells.
 *
 * parameters:
 *   mesh             -->  pointer to cs_mesh_t structure
 *   interface_set    -->  pointer to cs_interface_set_t structure.
 *   p_gcell_vtx_idx  <--  pointer to the connectivity index
 *   p_gcell_vtx_lst  <--  pointer to the connectivity list
 *---------------------------------------------------------------------------*/

void
cs_ctwr_halo_define(cs_ctwr_zone_t      *ct,
                    cs_interface_set_t  *interface_set)
{
  cs_int_t  *vtx_faces_idx = NULL;
  cs_int_t  *vtx_faces_lst = NULL;
  cs_int_t  *g_vtx_faces_idx = NULL;
  cs_int_t  *g_vtx_faces_lst = NULL;
  cs_int_t  *g_faces_vtx_idx = NULL;
  cs_int_t  *g_faces_vtx_lst = NULL;

  cs_int_t  *g_out_faces_vtx_idx = NULL;
  cs_int_t  *g_out_faces_vtx_lst = NULL;
  cs_int_t  *g_in_faces_vtx_idx = NULL;
  cs_int_t  *g_in_faces_vtx_lst = NULL;

  cs_halo_t  *halo = ct->water_halo;

  ct->halo_type = CS_HALO_EXTENDED;

  /*  Define vtx -> faces connectivity for ghost faces */

  _create_gvtx_faces_connect(ct,
                             interface_set,
                             &vtx_faces_idx,
                             &vtx_faces_lst);

  /* Fill cs_halo_t structure for send_halo  */

  bft_printf(_("    Local halo definition\n"));
  bft_printf_flush();

  _fill_send_halo(ct,
                  interface_set,
                  vtx_faces_idx,
                  vtx_faces_lst);


  cs_reverse_vtx_faces_connect(ct->face_sup_mesh,
                               &g_faces_vtx_idx,
                               &g_faces_vtx_lst);

  _fill_index_out_halo(ct);


  if (halo->n_elts[0] > 0) {

    _create_in_faces_vtx_connect(ct,
                                 interface_set,
                                 g_faces_vtx_idx,
                                 g_faces_vtx_lst,
                                 &g_in_faces_vtx_idx,
                                 &g_in_faces_vtx_lst);

    _exchange_gface_vtx_connect(ct,
                                g_in_faces_vtx_idx,
                                g_in_faces_vtx_lst,
                                &g_out_faces_vtx_idx,
                                &g_out_faces_vtx_lst);

    bft_printf(_("    Updating the vertex -> faces connectivity\n"));
    bft_printf_flush();

    _update_gfaces_connect(ct,
                           g_out_faces_vtx_idx,
                           g_out_faces_vtx_lst,
                           &g_vtx_faces_idx,
                           &g_vtx_faces_lst);

    /* Free memory */

    BFT_FREE(g_out_faces_vtx_idx);
    BFT_FREE(g_out_faces_vtx_lst);
    BFT_FREE(g_in_faces_vtx_idx);
    BFT_FREE(g_in_faces_vtx_lst);

  }

  _create_faces_faces_connect(ct,
                              g_faces_vtx_idx,
                              g_faces_vtx_lst,
                              g_vtx_faces_idx,
                              g_vtx_faces_lst);

  _build_halo_with_extrusion(ct);

  /* Fill cs_halo_t structure for out_halo.
     We use the data from in_halo structure */

  bft_printf(_("    Distant halo definition\n"));
  bft_printf_flush();

  _fill_halo(ct);

  BFT_FREE(vtx_faces_idx);
  BFT_FREE(vtx_faces_lst);

  /* Update mesh structure elements bound to halo management */

  ct->nnpsct_with_ghosts = ct->nnpsct + halo->n_elts[0] / ct->nelect;

  cs_halo_update_buffers(halo);

#if 0 /* for debugging purposes */
  cs_halo_dump(halo, 1);
#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
