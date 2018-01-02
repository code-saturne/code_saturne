/*============================================================================
 * Functions dealing with ghost cells
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_base.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_order.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

typedef struct _vtx_lookup_table {

  cs_lnum_t  n_vertices;    /* Number of local vertices in the problem domain */
  cs_lnum_t  n_transforms;  /* Number of transformations */
  cs_lnum_t  n_interfaces;  /* Number of interfaces */
  cs_lnum_t  n_categories;  /* Number of possible categories
                               = n_interfaces * (n_transforms + 1)
                               (1 category for purely parallel elements) */

  cs_lnum_t  *if_ranks;     /* List of ranks */
  cs_lnum_t  *rank_ids;     /* list of rank ids */

  cs_lnum_t  *index;        /* index on table (size = n_vertices + 1) */

  cs_lnum_t  *rank_list;    /* list of ranks on which vertices are linked */
  cs_lnum_t  *type_list;    /* list of type (purelly parallel (=0) or number
                               of the periodicity) featuring a vertex. This
                               list is only allocated when n_perio > 0 */

} vtx_lookup_table_t;


/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 *  Global static variables
 *============================================================================*/

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
 *   vtx_lookup  <--  pointer to a vtx_lookup_table_t structure
 *   ifs         <--  pointer to a fvm_interface_set_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_vtx_lookup(vtx_lookup_table_t        *vtx_lookup,
                 const cs_interface_set_t  *ifs)
{
  cs_lnum_t i, vtx_id, rank_id, shift;

  cs_lnum_t *counter = NULL;

  const cs_lnum_t n_interfaces = cs_interface_set_size(ifs);

  BFT_MALLOC(counter, vtx_lookup->n_vertices, cs_lnum_t);

  for (i = 0; i < vtx_lookup->n_vertices; i++)
    counter[i] = 0;

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    const cs_interface_t  *interface = cs_interface_set_get(ifs, rank_id);
    const cs_lnum_t  interface_size = cs_interface_size(interface);
    const cs_lnum_t  *elt_id = cs_interface_get_elt_ids(interface);

    for (i = 0; i < interface_size; i++) { /* Only parallel vertices */

      vtx_id = elt_id[i];
      shift = vtx_lookup->index[vtx_id] + counter[vtx_id];

      vtx_lookup->rank_list[shift] = vtx_lookup->rank_ids[rank_id];
      counter[vtx_id] += 1;

    }

  } /* End of loop on ranks */

  BFT_FREE(counter);

}

/*---------------------------------------------------------------------------
 * Fill a look up table structure with periodicity.
 *
 * The rank_list takes the distant rank value with which a local vertex
 * is linked.
 * The type_list takes the id number associated to the kind of link between
 * local and distant element. If the link is purely parallel, list gets value
 * 0 otherwise the list gets the number of the transformation (+ or -)
 *
 * parameters:
 *   vtx_look_up  <--  pointer to a vtx_lookup_table_t structure
 *   ifs          <--  pointer to a fvm_interface_set_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_vtx_lookup_with_perio(vtx_lookup_table_t        *vtx_lookup,
                            const cs_interface_set_t  *ifs)
{
  cs_lnum_t i, tr_id, vtx_id, rank_id, shift;

  cs_lnum_t *counter = NULL;

  const fvm_periodicity_t  *periodicity = cs_interface_set_periodicity(ifs);
  const cs_lnum_t n_interfaces = cs_interface_set_size(ifs);
  const cs_lnum_t n_transforms = fvm_periodicity_get_n_transforms(periodicity);

  assert(n_transforms > 0);

  BFT_MALLOC(counter, vtx_lookup->n_vertices, cs_lnum_t);

  for (i = 0; i < vtx_lookup->n_vertices; i++)
    counter[i] = 0;

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    const cs_interface_t  *interface
      = cs_interface_set_get(ifs, vtx_lookup->rank_ids[rank_id]);
    const cs_lnum_t  tr_index_size = cs_interface_get_tr_index_size(interface);
    const cs_lnum_t  *tr_index = cs_interface_get_tr_index(interface);
    const cs_lnum_t  *vtx_ids = cs_interface_get_elt_ids(interface);

    assert(n_transforms + 2 == tr_index_size);
    assert(tr_index != NULL);

    for (i = tr_index[0]; i < tr_index[1]; i++) { /* Only parallel vertices */

      vtx_id = vtx_ids[i];
      shift = vtx_lookup->index[vtx_id] + counter[vtx_id];

      vtx_lookup->rank_list[shift] = rank_id;
      vtx_lookup->type_list[shift] = 0;
      counter[vtx_id] += 1;

    }

    for (tr_id = 0; tr_id < n_transforms; tr_id++) {

      for (i = tr_index[tr_id + 1]; i < tr_index[tr_id + 2]; i++) {

        vtx_id = vtx_ids[i];
        shift = vtx_lookup->index[vtx_id] + counter[vtx_id];
        vtx_lookup->rank_list[shift] = rank_id;
        vtx_lookup->type_list[shift] = tr_id + 1;
        counter[vtx_id] += 1;

      }

    } /* End of loop on transformations */

  } /* End of loop on ranks */

  BFT_FREE(counter);

}

/*---------------------------------------------------------------------------
 * Create a vtx_look_up_table_t structure
 *
 * parameters:
 *   n_vertices  <--  number of vertices of the table.
 *   ifs         <--  pointer to a fvm_interface_set_t structure
 *
 * returns:
 *   A pointer to the created vtx_lookup_table_t structure
 *---------------------------------------------------------------------------*/

static vtx_lookup_table_t *
_vtx_lookup_create(cs_lnum_t                  n_vertices,
                   const cs_interface_set_t  *ifs)
{
  cs_lnum_t i, rank_id, tmp_id, interface_size;

  cs_lnum_t loc_rank_id = -1;
  vtx_lookup_table_t  *vtx_lookup = NULL;

  const cs_interface_t  *interface = NULL;
  const cs_lnum_t  *vtx_ids = NULL;
  const fvm_periodicity_t  *periodicity = cs_interface_set_periodicity(ifs);
  const cs_lnum_t n_transforms = fvm_periodicity_get_n_transforms(periodicity);
  const cs_lnum_t n_interfaces = cs_interface_set_size(ifs);

  BFT_MALLOC(vtx_lookup, 1, vtx_lookup_table_t);

  vtx_lookup->n_vertices = n_vertices;
  vtx_lookup->n_interfaces = n_interfaces;
  vtx_lookup->n_transforms = n_transforms;
  vtx_lookup->n_categories = (n_transforms + 1)*n_interfaces;

  BFT_MALLOC(vtx_lookup->index, n_vertices + 1, cs_lnum_t);
  BFT_MALLOC(vtx_lookup->if_ranks, n_interfaces, cs_lnum_t);
  BFT_MALLOC(vtx_lookup->rank_ids, n_interfaces, cs_lnum_t);

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
      cs_lnum_t  *buffer = NULL;
      cs_lnum_t *_rank_ids = NULL;

      assert(sizeof(cs_lnum_t) == sizeof(cs_lnum_t));

      BFT_MALLOC(order, n_interfaces - 1, cs_lnum_t);
      BFT_MALLOC(buffer, n_interfaces - 1, cs_lnum_t);
      BFT_MALLOC(_rank_ids, n_interfaces , cs_lnum_t);

      _rank_ids[0] = vtx_lookup->rank_ids[0];
      for (i = 1; i < n_interfaces; i++) {
        buffer[i-1] = vtx_lookup->if_ranks[i];
        _rank_ids[i] = vtx_lookup->rank_ids[i];
      }

      cs_order_lnum_allocated(NULL,
                              buffer,
                              order,
                              n_interfaces-1);

      for (i = 0; i < n_interfaces - 1; i++) {
        vtx_lookup->if_ranks[i+1] = buffer[order[i]];
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
    vtx_ids = cs_interface_get_elt_ids(interface);

    for (i = 0; i < interface_size; i++)
      vtx_lookup->index[vtx_ids[i] + 1] += 1;

  } /* End of loop on if_ranks */

  /* Create index and allocate buffers */

  for (i = 0; i < n_vertices; i++)
    vtx_lookup->index[i+1] += vtx_lookup->index[i];

  BFT_MALLOC(vtx_lookup->rank_list, vtx_lookup->index[n_vertices], cs_lnum_t);

  /* Second loop to fill table(s) */

  if (n_transforms > 0) {

    BFT_MALLOC(vtx_lookup->type_list, vtx_lookup->index[n_vertices], cs_lnum_t);
    _fill_vtx_lookup_with_perio(vtx_lookup, ifs);

  }
  else {

    vtx_lookup->type_list = NULL;
    _fill_vtx_lookup(vtx_lookup, ifs);

  }

  return vtx_lookup;
}

/*---------------------------------------------------------------------------
 * Destroy a vtx_lookup structure.
 *
 * parameters:
 *   vtx_lookup <--  pointer to a vtx_lookup_table_t structure
 *---------------------------------------------------------------------------*/

static void
_vtx_lookup_destroy(vtx_lookup_table_t  **vtx_lookup)
{
  vtx_lookup_table_t  *vl = *vtx_lookup;

  BFT_FREE(vl->if_ranks);
  BFT_FREE(vl->rank_ids);
  BFT_FREE(vl->index);
  BFT_FREE(vl->rank_list);

  if (vl->type_list != NULL)
    BFT_FREE(vl->type_list);

  BFT_FREE(*vtx_lookup);
}

/*---------------------------------------------------------------------------
 * Set checker for this vertex_id according to vtx_lookup features.
 *
 * parameters:
 *   vtx_id       <--  vertex id to deal with
 *   vtx_checker  <->  put a tag in the implied categories
 *   vtx_lookup   <--  pointer to a vtx_lookup_table_t structure
 *---------------------------------------------------------------------------*/

static void
_update_vtx_checker(cs_lnum_t            vtx_id,
                    cs_lnum_t           *vtx_checker,
                    vtx_lookup_table_t  *vtx_lookup)
{
  cs_lnum_t i, rank_id, type;

  const cs_lnum_t n_interfaces = vtx_lookup->n_interfaces;
  const cs_lnum_t n_transforms = vtx_lookup->n_transforms;

  for (i = vtx_lookup->index[vtx_id];
       i < vtx_lookup->index[vtx_id + 1]; i++) {

    rank_id = vtx_lookup->rank_list[i];

    if (n_transforms == 0)  /* purely parallel */

      vtx_checker[rank_id] += 1;

    else { /* n_perio > 0 */

      type = vtx_lookup->type_list[i];
      vtx_checker[type*n_interfaces + rank_id] += 1;

    }

  } /* End of loop on vtx_lookup */

}

/*---------------------------------------------------------------------------
 * Update the halo type of each face checker element according to vtx_checker
 * values for this face.
 *
 * parameters:
 *   n_face_vertices <-- number of vertices defining a face
 *   n_categories    <-- size of vtx_checker and face_checker
 *   vtx_checker     <-- number of vertices for each categories for the face
 *   face_checker    --> halo_type value in each categories
 *
 * returns:
 *   maximum halo_type value in face_checker
 *---------------------------------------------------------------------------*/

static cs_halo_type_t
_update_face_checker(cs_lnum_t         n_face_vertices,
                     cs_lnum_t         n_categories,
                     cs_lnum_t        *vtx_checker,
                     cs_halo_type_t   *face_checker)
{
  cs_lnum_t i;

  cs_halo_type_t  ret_type = CS_HALO_N_TYPES;

  for (i = 0; i < n_categories; i++) {

    if (vtx_checker[i] == n_face_vertices) { /* => STANDARD HALO */

      face_checker[i] = CS_HALO_STANDARD;
      ret_type = CS_HALO_STANDARD;

    }
    else {

      if (vtx_checker[i] > 0) {  /* => EXTENDED HALO */

        if (ret_type == CS_HALO_N_TYPES)
          ret_type = CS_HALO_EXTENDED;

        if (face_checker[i] == CS_HALO_N_TYPES)
          face_checker[i] = CS_HALO_EXTENDED;

      }

    }

  } /* End of loop on categories */

  return ret_type;
}

/*---------------------------------------------------------------------------
 * Update counter on number of elements in each category of halos
 * according to face_checker values.
 *
 * parameters:
 *   mesh         <-- pointer to a mesh structure
 *   face_checker <-- halo type in each categories for cell's faces
 *---------------------------------------------------------------------------*/

static void
_count_halo_elements(cs_mesh_t        *mesh,
                     cs_halo_type_t   *face_checker)
{
  cs_lnum_t type_id, rank_id;

  const cs_lnum_t n_transforms = mesh->n_transforms;
  const cs_halo_t  *halo = mesh->halo;
  const cs_lnum_t n_c_domains = halo->n_c_domains;
  const cs_lnum_t stride = 4*n_c_domains;

  for (type_id = 0; type_id < n_transforms + 1; type_id++) {

    for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

      if (face_checker[type_id*n_c_domains + rank_id]
          == CS_HALO_STANDARD) {

        halo->send_index[2*rank_id + 1] += 1;

        if (type_id > 0) /* periodic elements */
          halo->send_perio_lst[stride*(type_id-1) + 4*rank_id + 1] += 1;

      } /* STANDARD HALO */

      else if (face_checker[type_id*n_c_domains + rank_id]
               == CS_HALO_EXTENDED) {

        if (mesh->halo_type == CS_HALO_EXTENDED) {

          halo->send_index[2*rank_id + 2] += 1;

          if (type_id > 0) /* periodic elements */
            halo->send_perio_lst[stride*(type_id-1) + 4*rank_id + 3] += 1;

        }

      } /* EXTENDED HALO */

    } /* End of loop on ranks */

  } /* End of loop on categories */

}

/*---------------------------------------------------------------------------
 * Build halo indexes
 *
 * parameters:
 *   mesh  <-> pointer to a mesh structure
 *---------------------------------------------------------------------------*/

static void
_build_halo_index(cs_mesh_t  *mesh)
{
  cs_lnum_t i, rank_id;
  cs_lnum_t buffer_size, n_halo_elts, n_per_halo_elts;

  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t n_c_domains = halo->n_c_domains;
  const cs_lnum_t n_init_perio = mesh->n_init_perio;
  const cs_lnum_t stride = 4*n_c_domains;
  const cs_lnum_t n_transforms = mesh->n_transforms;
  const cs_lnum_t local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;

  buffer_size = 0;
  for (i = 0; i < 2*n_c_domains; i++)
    buffer_size += halo->send_index[i+1];

  BFT_MALLOC(halo->send_list, buffer_size, cs_lnum_t);

  /* Define parallel and periodic index */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    /* Standard halo */
    /* ------------- */

    /* Count number of periodic elements in standard halo for this rank */

    n_per_halo_elts = 0;
    for (i = 0; i < n_transforms; i++)
      n_per_halo_elts += halo->send_perio_lst[stride*i + 4*rank_id + 1];

    /* Define index */

    n_halo_elts = halo->send_index[2*rank_id+1];
    halo->send_index[2*rank_id+1] += halo->send_index[2*rank_id];

    assert(n_halo_elts >= n_per_halo_elts);

    /* Fill perio_lst buffer */

    if (n_init_perio > 0) {

      if (halo->c_domain_rank[rank_id] == local_rank)
        halo->send_perio_lst[4*rank_id] = halo->send_index[2*rank_id];

      else
        halo->send_perio_lst[4*rank_id] =
          halo->send_index[2*rank_id] + (n_halo_elts - n_per_halo_elts);

      for (i = 0; i < n_transforms - 1; i++)
        halo->send_perio_lst[stride*(i+1) + 4*rank_id]
          =  halo->send_perio_lst[stride*i + 4*rank_id]
           + halo->send_perio_lst[stride*i + 4*rank_id + 1];

    } /* Test if n_perio > 0 */

    /* Extended halo */
    /* ------------- */

    n_per_halo_elts = 0;
    for (i = 0; i < n_transforms; i++)
      n_per_halo_elts += halo->send_perio_lst[stride*i + 4*rank_id+3];

    n_halo_elts = halo->send_index[2*rank_id+2];
    halo->send_index[2*rank_id+2] += halo->send_index[2*rank_id+1];

    assert(n_halo_elts >= n_per_halo_elts);

    if (n_init_perio > 0) {

      if (halo->c_domain_rank[rank_id] == local_rank)
        halo->send_perio_lst[4*rank_id+2] = halo->send_index[2*rank_id+1];

      else
        halo->send_perio_lst[4*rank_id+2] =
          halo->send_index[2*rank_id+1] + (n_halo_elts - n_per_halo_elts);

      for (i = 0; i < n_transforms - 1; i++)
        halo->send_perio_lst[stride*(i+1) + 4*rank_id + 2]
          =  halo->send_perio_lst[stride*i + 4*rank_id + 2]
           + halo->send_perio_lst[stride*i + 4*rank_id + 3];

    } /* Test if n_perio > 0 */

  } /* End of loop on c_domain_rank */

}

/*---------------------------------------------------------------------------
 * Fill ghost cells list (member of cs_halo_t structure)
 *
 * parameters:
 *   mesh          <--  pointer to a mesh structure.
 *   face_checker  <--  halo type of each face of the cell.
 *   cell_id       <--  numbering of the treated cell.
 *   counter       <->  counter on each categories.
 *---------------------------------------------------------------------------*/

static void
_add_halo_elements(cs_mesh_t        *mesh,
                   cs_halo_type_t   *face_checker,
                   cs_lnum_t         cell_id,
                   cs_lnum_t        *counter)
{
  cs_lnum_t i, type_id, shift, c_shift;

  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t n_transforms = mesh->n_transforms;
  const cs_lnum_t n_c_domains = halo->n_c_domains;

    /* How is organized counter:

         -------------------------------------------------
 Paral:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          std ext std ext std ext std ext std ext std ext
          _______ _______ _______ _______ _______ _______
           rank0   rank1   rank2   rank3   rank4   rank5

         -------------------------------------------------
    P1:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          std ext std ext std ext std ext std ext std ext
          _______ _______ _______ _______ _______ _______
           rank0   rank1   rank2   rank3   rank4   rank5

         -------------------------------------------------
   P-1:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          std ext std ext std ext std ext std ext std ext
          _______ _______ _______ _______ _______ _______
           rank0   rank1   rank2   rank3   rank4   rank5

    etc...

    */

  for (type_id = 0; type_id < n_transforms + 1; type_id++) {

    for (i = 0; i < n_c_domains; i++) {

      if (face_checker[type_id*n_c_domains + i] == CS_HALO_STANDARD) {

        c_shift = 2*n_c_domains*type_id + 2*i;

        if (type_id == 0)
          shift = halo->send_index[2*i] + counter[c_shift];
        else
          shift = halo->send_perio_lst[4*n_c_domains*(type_id-1) + 4*i]
                + counter[c_shift];

        halo->send_list[shift] = cell_id;
        counter[c_shift] += 1;

      }
      else if (face_checker[type_id*n_c_domains + i] == CS_HALO_EXTENDED) {

        if (mesh->halo_type == CS_HALO_EXTENDED) {

          c_shift = 2*n_c_domains*type_id + 2*i + 1;

          if (type_id == 0)
            shift = halo->send_index[2*i+1] + counter[c_shift];
          else
            shift = halo->send_perio_lst[4*n_c_domains*(type_id-1) + 4*i + 2]
                  + counter[c_shift];

          halo->send_list[shift] = cell_id;
          counter[c_shift] += 1;

        } /* End of extended halo treatment */

      }

    } /* End of loop on interfaces */

  } /* End of loop on categories */

}

/*---------------------------------------------------------------------------
 * Test if loop has to be continued according halo type and face number.
 *
 * parameters:
 *   mesh    <-- pointer to mesh structure
 *   face_id <-- interior face id
 *---------------------------------------------------------------------------*/

static bool
_test_loop_continues(cs_mesh_t  *mesh,
                     cs_lnum_t   face_id)
{
  bool  choice = false;

  /* Face has to be an internal face */

  if (mesh->halo_type == CS_HALO_STANDARD) {

    if (   mesh->i_face_cells[face_id][0] < 0
        || mesh->i_face_cells[face_id][1] < 0)
      choice = true;
    else
      choice = false;

  }
  else {

    assert(mesh->halo_type == CS_HALO_EXTENDED);
    choice = true;

  }

  return choice;
}

/*---------------------------------------------------------------------------
 * Define the elements of send_halo structure.
 *
 * Two main loops. First one for counting number of elements and create index.
 * Second one for filling the ghost cells list.
 *
 * parameters:
 *   mesh            <-- pointer to cs_mesh_t structure
 *   vertex_ifs      <-- pointer to fvm_vertex_ifs_t structure
 *   gcell_faces_idx <-- "ghost cell -> faces" connectivity index
 *   gcell_faces_lst <-- "ghost cell -> faces" connectivity list
 *---------------------------------------------------------------------------*/

static void
_fill_send_halo(cs_mesh_t                 *mesh,
                const cs_interface_set_t  *vertex_ifs,
                cs_lnum_t                 *gcell_faces_idx,
                cs_lnum_t                 *gcell_faces_lst)
{
  cs_lnum_t i, cell_id, i_fac, i_vtx;
  cs_lnum_t fac_id, vtx_id;
  cs_lnum_t n_face_vertices;

  cs_halo_type_t  type_tag = CS_HALO_N_TYPES;
  cs_halo_type_t  face_type = CS_HALO_N_TYPES;
  cs_halo_type_t  cell_type = CS_HALO_N_TYPES;
  cs_lnum_t n_categories = 0;
  cs_halo_t  *halo = mesh->halo;
  vtx_lookup_table_t  *vtx_lookup = NULL;
  cs_halo_type_t  *cell_tag = NULL;
  cs_halo_type_t  *face_checker = NULL;
  cs_lnum_t *vtx_checker = NULL;
  cs_lnum_t *counter = NULL;

  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_vertices = mesh->n_vertices;
  const cs_lnum_t *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *fac_vtx_lst = mesh->i_face_vtx_lst;

  /* We should have the faces -> vertices connectivity to continue */

  if (mesh->i_face_vtx_lst == NULL)
    return;

  /* Create a lookup table to accelerate search in
     cs_interface_set structure */

  vtx_lookup = _vtx_lookup_create(n_vertices, vertex_ifs);

  n_categories = vtx_lookup->n_categories;

  /* How is organized vtx_checker and face_checker:

            -------------------------
    Paral:  |     |     |     |     |
            -------------------------
             rank0 rank1 rank2 rank3

            -------------------------
       P1:  |     |     |     |     |
            -------------------------
             rank0 rank1 rank2 rank3

            -------------------------
      P-1:  |     |     |     |     |
            -------------------------
             rank0 rank1 rank2 rank3

  */

  if (gcell_faces_idx != NULL && gcell_faces_lst != NULL) {

    BFT_MALLOC(vtx_checker, n_categories, cs_lnum_t);
    BFT_MALLOC(face_checker, n_categories, cs_halo_type_t);
    BFT_MALLOC(cell_tag, n_cells, cs_halo_type_t);

    /* First loop to create index and allocate cell_list */

    for (cell_id = 0; cell_id < n_cells; cell_id++) {

      /* Default initialization */

      cell_type = CS_HALO_N_TYPES;

      for (i = 0; i < n_categories; i++)
        face_checker[i] = CS_HALO_N_TYPES;

      /* Loop on faces of the cell */

      for (i_fac = gcell_faces_idx[cell_id];
           i_fac < gcell_faces_idx[cell_id + 1];
           i_fac++) {

        fac_id = gcell_faces_lst[i_fac];

        if (_test_loop_continues(mesh, fac_id) == true) {

          n_face_vertices = fac_vtx_idx[fac_id + 1] - fac_vtx_idx[fac_id];

          /* Initialize checker */

          for (i = 0; i < n_categories; i++)
            vtx_checker[i] = 0;

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id + 1];
               i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx];

            _update_vtx_checker(vtx_id, vtx_checker, vtx_lookup);

          } /* End of loop on vertices */

          face_type = _update_face_checker(n_face_vertices,
                                           n_categories,
                                           vtx_checker,
                                           face_checker);

          cell_type = CS_MIN(face_type, cell_type);

        } /* If the face is an internal face */

      } /* End of loop on faces */

      _count_halo_elements(mesh, face_checker);

      cell_tag[cell_id] = cell_type;

    } /* End of loop on cells */

    /* Build halo index */

    _build_halo_index(mesh);

    /* Initialize counter */

    BFT_MALLOC(counter, 2*n_categories, cs_lnum_t);

    for (i = 0; i < 2*n_categories; i++)
      counter[i] = 0;

    /* Second loop to build halo->ghost_cells */

    for (cell_id = 0; cell_id < n_cells; cell_id++) {

      if (mesh->halo_type == CS_HALO_STANDARD)
        type_tag = CS_HALO_EXTENDED;
      else if (mesh->halo_type == CS_HALO_EXTENDED)
        type_tag = CS_HALO_N_TYPES;

      if (cell_tag[cell_id] < type_tag) {

        for (i = 0; i < n_categories; i++)
          face_checker[i] = CS_HALO_N_TYPES;

        /* Loop on faces of the cell */

        for (i_fac = gcell_faces_idx[cell_id];
             i_fac < gcell_faces_idx[cell_id+1];
             i_fac++) {

          /* Initialize checker */

          for (i = 0; i < n_categories; i++)
            vtx_checker[i] = 0;

          fac_id = gcell_faces_lst[i_fac];

          if (_test_loop_continues(mesh, fac_id) ==  true) {

            /* Loop on vertices of the face */

            n_face_vertices = fac_vtx_idx[fac_id + 1] - fac_vtx_idx[fac_id];

            for (i_vtx = fac_vtx_idx[fac_id];
                 i_vtx < fac_vtx_idx[fac_id + 1];
                 i_vtx++) {

              vtx_id = fac_vtx_lst[i_vtx];

              _update_vtx_checker(vtx_id,
                                  vtx_checker,
                                  vtx_lookup);

            } /* End of loop on vertices */

            face_type = _update_face_checker(n_face_vertices,
                                             n_categories,
                                             vtx_checker,
                                             face_checker);

          } /* If face is an internal face */

        } /* End of loop on faces */

        _add_halo_elements(mesh,
                           face_checker,
                           cell_id,
                           counter);

      } /* End of test on cell_list */

    } /* End of loop on cells */

    /* Free memory */

    BFT_FREE(vtx_checker);
    BFT_FREE(face_checker);
    BFT_FREE(counter);
    BFT_FREE(cell_tag);

  } /* End if gcell_face_idx and gcell_faces_lst != NULL */

  /* Destroy the lookup table strcuture */

  _vtx_lookup_destroy(&vtx_lookup);

  /* Complete halo definition */

  halo->n_send_elts[CS_HALO_STANDARD] = 0;
  halo->n_send_elts[CS_HALO_EXTENDED] = 0;

  for (i = 0; i < halo->n_c_domains; i++) {

    halo->n_send_elts[CS_HALO_STANDARD] += halo->send_index[2*i+1]
                                         - halo->send_index[2*i];
    halo->n_send_elts[CS_HALO_EXTENDED] += halo->send_index[2*i+2]
                                         - halo->send_index[2*i+1];

  }

  halo->n_send_elts[CS_HALO_EXTENDED] += halo->n_send_elts[CS_HALO_STANDARD];
}

/*---------------------------------------------------------------------------
 * Define a buffer on vertices where vertex belonging to the vertex_ifs
 * are tagged with 1 else 0.
 *
 * parameters:
 *   n_vertices    <-- size of the buffer
 *   vertex_ifs    <-- pointer to a cs_interface_set_t structure
 *   p_vertex_tag  <-> pointer to the tagged buffer
 *---------------------------------------------------------------------------*/

static void
_get_vertex_tag(cs_lnum_t                  n_vertices,
                const cs_interface_set_t  *vertex_ifs,
                cs_lnum_t                 *p_vertex_tag[])
{
  cs_lnum_t i, j, rank_id;

  cs_lnum_t *vertex_tag = NULL;

  const int  ifs_size = cs_interface_set_size(vertex_ifs);

  BFT_MALLOC(vertex_tag, n_vertices, cs_lnum_t);

  for (i = 0; i < n_vertices; i++)
    vertex_tag[i] = 0;

  for (rank_id = 0; rank_id < ifs_size; rank_id++) {

    const cs_interface_t  *interface = cs_interface_set_get(vertex_ifs, rank_id);
    const cs_lnum_t  *vtx_ids = cs_interface_get_elt_ids(interface);
    const cs_lnum_t  if_size = cs_interface_size(interface);

    for (j = 0; j < if_size; j++)
      vertex_tag[vtx_ids[j]] = 1;

  } /* End of loop on ranks */

  *p_vertex_tag = vertex_tag;

}

/*---------------------------------------------------------------------------
 * Compute the number of purely parallel ghost cells for a specific rank.
 *
 * parameters:
 *   mesh      <-- pointer to cs_mesh_t structure
 *   rank      <-- rank on which we want to know the number of purely
 *                 parallel elements.
 *   type      <-- standard or extended
 *   index     <-- index on halo's elements
 *   perio_lst <-- periodic details on halo
 *
 * returns:
 *  Number of purely parallel elements in the halo.
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_get_n_par_ghost_cells(const cs_mesh_t  *mesh,
                       cs_lnum_t         rank,
                       cs_halo_type_t    type,
                       const cs_lnum_t   index[],
                       const cs_lnum_t   perio_lst[])
{
  cs_lnum_t i;
  cs_lnum_t n_per_gcells = 0, n_par_gcells = 0;

  const cs_lnum_t n_transforms = mesh->n_transforms;
  const cs_lnum_t n_c_domains = mesh->halo->n_c_domains;

  if (type == CS_HALO_STANDARD) {

    for (i = 0; i < n_transforms; i++)
      n_per_gcells += perio_lst[4*rank+1 + 4*n_c_domains*i];

    n_par_gcells = index[2*rank+1] - index[2*rank];
    n_par_gcells -= n_per_gcells;

  }
  else if (type == CS_HALO_EXTENDED) {

    for (i = 0; i < n_transforms; i++)
      n_per_gcells += perio_lst[4*rank+3 + 4*n_c_domains*i];

    n_par_gcells = index[2*rank+2] - index[2*rank+1];
    n_par_gcells -= n_per_gcells;

  }

  return n_par_gcells;
}

/*---------------------------------------------------------------------------
 * Exchange number and list of cells constituting send_halo structure for each
 * frontier ranks. Fill the halo structure from these data.
 *
 * parameters:
 *   mesh <-- pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_halo(cs_mesh_t  *mesh)
{
  cs_lnum_t rank_id, i, j;
  cs_lnum_t shift;

#if defined(HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  int request_count = 0;

  cs_lnum_t *count = NULL;

  cs_halo_t  *halo = mesh->halo;

  const  cs_lnum_t n_c_domains = halo->n_c_domains;
  const  cs_lnum_t n_transforms = mesh->n_transforms;
  const  cs_lnum_t local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;

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
      MPI_Irecv(&(halo->index[2*rank_id+1]), 2, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                halo->c_domain_rank[rank_id],
                cs_glob_mpi_comm,
                &(request[request_count++]));
#endif

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Barrier(cs_glob_mpi_comm);
#endif

  BFT_MALLOC(count, 2*n_c_domains, cs_lnum_t);

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = 2*rank_id;
    count[shift] =   halo->send_index[2*rank_id+1]
                   - halo->send_index[2*rank_id];
    count[shift+1] =   halo->send_index[2*rank_id+2]
                     - halo->send_index[2*rank_id+1];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(HAVE_MPI)
      MPI_Isend(&(count[shift]), 2, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                local_rank,
                cs_glob_mpi_comm,
                &(request[request_count++]));
#endif

    }
    else {

      halo->index[shift+1] = count[shift];
      halo->index[shift+2] = count[shift+1];

    }

  } /* End of loop on ranks */

  /* Wait for all exchanges being done */

#if defined(HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Waitall(request_count, request, status);
#endif
  request_count = 0;

  BFT_FREE(count);

  /* Build index */
  /*-------------*/

  for (i = 0; i < 2*n_c_domains; i++)
    halo->index[i+1] += halo->index[i];

  /* Exchange number of elements for each periodicity and for each rank.
     Then build halo->perio_lst */

  if (mesh->n_init_perio > 0) {

    cs_lnum_t n_elts;
    cs_lnum_t *exchange_buffer = NULL;

    /* n_transforms periodicities to deal with and for each sub-periodicity
       2 data. One for standard halo and the other one for extended halo */

    const cs_lnum_t n_elts_to_exchange = 2 * n_transforms;

    BFT_MALLOC(exchange_buffer, 4*n_transforms, cs_lnum_t);

    for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

      if (halo->c_domain_rank[rank_id] != local_rank) {

        /* Fill buffer to send */

        for (i = 0; i < n_transforms; i++) {
          shift = 4*n_c_domains*i + 4*rank_id;
          for (j = 0; j < 2; j++)
            exchange_buffer[2*i+j] = halo->send_perio_lst[shift + 2*j + 1];
        }

#if defined(HAVE_MPI)
        MPI_Sendrecv(&(exchange_buffer[0]), n_elts_to_exchange, CS_MPI_INT,
                     halo->c_domain_rank[rank_id], local_rank,
                     &(exchange_buffer[n_elts_to_exchange]), n_elts_to_exchange,
                     CS_MPI_INT,
                     halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                     cs_glob_mpi_comm, status);
#endif

        /* Put received elements in the periodic structure */

        for (i = 0; i < n_transforms; i++) {
          shift = 4*n_c_domains*i + 4*rank_id;
          for (j = 0; j < 2; j++)
            halo->perio_lst[shift + 2*j + 1] =
              exchange_buffer[n_elts_to_exchange + 2*i + j];
        }

      } /* rank != local_rank */

      else {

        for (i = 0; i < n_transforms; i++) {

          shift = 4*n_c_domains*i + 4*rank_id;
          for (j = 0; j < 2; j++)
            halo->perio_lst[shift + 2*j + 1] =
              halo->send_perio_lst[shift + 2*j + 1];

        } /* End of loop on periodicities */

      } /* local_rank == rank */

    } /* End of loop on communicating ranks */

    BFT_FREE(exchange_buffer);

    /* Build index values for perio_lst */

    for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

      /* Build index for standard ghost cells */

      n_elts = _get_n_par_ghost_cells(mesh,
                                      rank_id,
                                      CS_HALO_STANDARD,
                                      halo->index,
                                      halo->perio_lst);

      halo->perio_lst[4*rank_id] = halo->index[2*rank_id] + n_elts;

      for (i = 0; i < n_transforms - 1; i++) {
        shift = 4*n_c_domains*i + 4*rank_id;
        halo->perio_lst[4*n_c_domains + shift] =
          halo->perio_lst[shift] + halo->perio_lst[shift + 1];
      }

      /* Build index for extended ghost cells */

      n_elts = _get_n_par_ghost_cells(mesh,
                                      rank_id,
                                      CS_HALO_EXTENDED,
                                      halo->index,
                                      halo->perio_lst);

      halo->perio_lst[4*rank_id+2] = halo->index[2*rank_id+1] + n_elts;

      for (i = 0; i < n_transforms - 1; i++) {
        shift = 4*n_c_domains*i + 4*rank_id + 2;
        halo->perio_lst[4*n_c_domains + shift] =
          halo->perio_lst[shift] + halo->perio_lst[shift + 1];
      }

    } /* End of loop on communicating ranks */

  } /* End if n_perio > 0 */

#if defined(HAVE_MPI)
  if (request != _request) {
    BFT_FREE(request);
    BFT_FREE(status);
  }
#endif

  halo->n_elts[CS_HALO_STANDARD] = 0;
  halo->n_elts[CS_HALO_EXTENDED] = 0;

  for (i = 0; i < n_c_domains; i++) {

    halo->n_elts[CS_HALO_STANDARD] += halo->index[2*i+1] - halo->index[2*i];
    halo->n_elts[CS_HALO_EXTENDED] += halo->index[2*i+2] - halo->index[2*i+1];

  }

  halo->n_elts[CS_HALO_EXTENDED] += halo->n_elts[CS_HALO_STANDARD];
}

/*---------------------------------------------------------------------------
 * Compute maximum list buffer size.
 *
 * This is done to avoid a reallocation for each rank and transformation.
 *
 * parameters:
 *   ifs <-- pointer to a cs_interface_set_t structure
 *
 * returns:
 *  max buffer size
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_get_list_buffer_size(const cs_interface_set_t  *ifs)
{
  cs_lnum_t i, j, tr_index_size;

  cs_lnum_t max_lst_size = 0;

  const cs_interface_t  *interface = NULL;
  const cs_lnum_t  *tr_index = NULL;
  const cs_lnum_t ifs_size = cs_interface_set_size(ifs);

  if (ifs == NULL)
    return max_lst_size;

  for (i = 0; i < ifs_size; i++) {

    interface = cs_interface_set_get(ifs, i);
    tr_index = cs_interface_get_tr_index(interface);
    tr_index_size = cs_interface_get_tr_index_size(interface) - 1;

    if (tr_index != NULL)
      for (j = 0; j < tr_index_size; j++)
        max_lst_size = CS_MAX(max_lst_size, tr_index[j+1] - tr_index[j]);
    else
      max_lst_size = CS_MAX(max_lst_size,
                            (cs_lnum_t)cs_interface_size(interface));

  } /* End of loop on interfaces */

  return max_lst_size;
}

/*---------------------------------------------------------------------------
 * Define an index on vertices belonging to this interface for this rank
 * and this transformation.
 *
 * parameters:
 *   ifs                <-- pointer to cs_interface_set_t structure
 *   rank_id            <-- rank number to work with
 *   tr_id              <-- transformation id to work with
 *   vtx_interface_idx  <-> index on vertices which match the criterions
 *---------------------------------------------------------------------------*/

static void
_define_vtx_interface_idx(const cs_interface_set_t  *ifs,
                          cs_lnum_t                  rank_id,
                          cs_lnum_t                  tr_id,
                          cs_lnum_t                  n_vertices,
                          cs_lnum_t                 *vtx_interface_idx)
{
  cs_lnum_t i, j, id;

  /* Initialize index */

  for (i = 0; i < n_vertices + 1; i++)
    vtx_interface_idx[i] = 0;

  for (id = 0; id < cs_interface_set_size(ifs); id++) {

    const cs_interface_t  *interface = cs_interface_set_get(ifs, id);

    if (rank_id == cs_interface_rank(interface)) {

      const cs_lnum_t  *tr_index = cs_interface_get_tr_index(interface);
      const cs_lnum_t  *vtx_ids = cs_interface_get_elt_ids(interface);

      if (tr_index == NULL) {  /*  purely parallel elements */

        for (j = 0; j < (cs_lnum_t)cs_interface_size(interface); j++)
          vtx_interface_idx[vtx_ids[j] + 1] += 1;

      }
      else {

        for (j = tr_index[tr_id]; j < tr_index[tr_id+1]; j++)
          vtx_interface_idx[vtx_ids[j] + 1] += 1;

      }

      /* Create index */

      for (j = 0; j < n_vertices; j++)
        vtx_interface_idx[j+1] += vtx_interface_idx[j];

      break;

    }

  } /* End of loop on interfaces */

}

/*---------------------------------------------------------------------------
 * Define an index on vertices belonging the interface for this rank
 * and this perio.
 * Fill the dist_id_lst which is the list of distant ids associated
 * to local vertices (same rank and same transformation).
 *
 * parameters:
 *   ifs               <-- pointer to cs_interface_set_t structure
 *   rank_id           <-- rank number to work with
 *   tr_id             <-- transformation id to work with
 *   n_vertices        <-- number of vertices
 *   dist_id_lst       <-> list of distant vertex numbers matching criteria
 *   vtx_interface_idx <-> index on vertices matching criteria
 *---------------------------------------------------------------------------*/

static void
_define_dist_id_lst(const cs_interface_set_t  *ifs,
                    int                        rank_id,
                    int                        tr_id,
                    cs_lnum_t                  n_vertices,
                    cs_lnum_t                 *dist_id_lst,
                    cs_lnum_t                 *vtx_interface_idx)
{
  cs_lnum_t  i, j, id, shift;

  /* Initialize index */

  for (i = 0; i < n_vertices + 1; i++)
    vtx_interface_idx[i] = 0;

  for (id = 0; id < cs_interface_set_size(ifs); id++) {

    const cs_interface_t  *interface = cs_interface_set_get(ifs, id);

    if (rank_id == cs_interface_rank(interface)) {

      const cs_lnum_t  *tr_index = cs_interface_get_tr_index(interface);
      const cs_lnum_t  *vtx_ids = cs_interface_get_elt_ids(interface);
      const cs_lnum_t  *dist_vtx_ids = cs_interface_get_match_ids(interface);

      if (tr_index == NULL)
        for (j = 0; j < (cs_lnum_t)cs_interface_size(interface); j++)
          vtx_interface_idx[vtx_ids[j] + 1] += 1;

      else
        for (j = tr_index[tr_id]; j < tr_index[tr_id+1]; j++)
          vtx_interface_idx[vtx_ids[j] + 1] += 1;

      /* Create index */

      for (j = 0; j < n_vertices; j++)
        vtx_interface_idx[j+1] += vtx_interface_idx[j];

      /* There must be only one distant vertex id per vtx_id when
         we handle a specific rank and a specific periodicity.
         So, we don't need a counter to fill dist_id_lst */

      if (tr_index == NULL) {

        for (j = 0; j < (cs_lnum_t)cs_interface_size(interface); j++) {
          shift = vtx_interface_idx[vtx_ids[j]];
          dist_id_lst[shift] = dist_vtx_ids[j];
        }


      }
      else {

        for (j = tr_index[tr_id]; j < tr_index[tr_id+1]; j++) {
          shift = vtx_interface_idx[vtx_ids[j]];
          dist_id_lst[shift] = dist_vtx_ids[j];
        }

      }

      break;

    }

  } /* End of loop on interfaces */

}

/*---------------------------------------------------------------------------
 * Compute the start and end index in ghost cells list for the send_halo
 * elements according to its rank, its periodicity and its type.
 *
 * parameters:
 *   mesh        <-- pointer to cs_mesh_t structure
 *   index       <-- index on halo's elements
 *   perio_lst   <-- periodic details on halo
 *   rank_id     <-- rank number to work with
 *   tr_id       <-- transformation id to work with
 *   type_id     <-- standard or extended
 *   p_start_idx --> pointer on start index
 *   p_end_idx   --> pointer on end index
 *---------------------------------------------------------------------------*/

static void
_get_start_end_idx(const cs_mesh_t    *mesh,
                   const cs_lnum_t    *index,
                   const cs_lnum_t    *perio_lst,
                   cs_lnum_t           rank_id,
                   cs_lnum_t           tr_id,
                   cs_lnum_t           type_id,
                   cs_lnum_t          *p_start_idx,
                   cs_lnum_t          *p_end_idx)
{
  cs_lnum_t i, n_par_gcells, n_per_gcells;
  cs_lnum_t start_idx = -1, end_idx = -1;

  const cs_lnum_t n_c_domains = mesh->halo->n_c_domains;

  if (tr_id == 0) { /* Purelly parallel elements */

    if (type_id == 0) { /* STANDARD HALO */

      n_par_gcells = _get_n_par_ghost_cells(mesh,
                                            rank_id,
                                            CS_HALO_STANDARD,
                                            index,
                                            perio_lst);

      start_idx = index[2*rank_id];
      end_idx = start_idx + n_par_gcells;

    }

    if (type_id == 1) { /* EXTENDED HALO */

      n_par_gcells = _get_n_par_ghost_cells(mesh,
                                            rank_id,
                                            CS_HALO_EXTENDED,
                                            index,
                                            perio_lst);

      start_idx = index[2*rank_id+1];
      end_idx = start_idx + n_par_gcells;

    }

  }
  else { /* Periodic elements */

    i = tr_id - 1;

    if (type_id == 0) { /* STANDARD HALO */

      n_per_gcells = perio_lst[4*rank_id + 1 + 4*n_c_domains*i];
      start_idx = perio_lst[4*rank_id + 4*n_c_domains*i];
      end_idx = start_idx + n_per_gcells;

    }

    if (type_id == 1) { /* EXTENDED HALO */

      n_per_gcells = perio_lst[4*rank_id + 3 + 4*n_c_domains*i];
      start_idx = perio_lst[4*rank_id + 2 + 4*n_c_domains*i];
      end_idx = start_idx + n_per_gcells;

    }

  } /* Parallel or periodic elements */

  *p_start_idx = start_idx;
  *p_end_idx = end_idx;

}

/*---------------------------------------------------------------------------
 * Compute the size of the "ghost cell to distant vertices" connectivity for
 * send_halo elements.
 * This will be use to define "ghost cell to distant vertices" index.
 *
 * parameters:
 *   mesh               <-- pointer to cs_mesh_t structure
 *   ifs                <-- pointer to cs_interface_set_t structure
 *   rank_id            <-- rank number to work with
 *   tr_id              <-- transformation id to work with
 *   cell_faces_idx     <-- "cell -> faces" connectivity index
 *   cell_faces_lst     <-- "cell -> faces" connectivity list
 *   vtx_interface_idx  <-> index on vertices which match the criterions
 *   vtx_tag            <-> tag array on vertices
 *   gcell_dist_vtx_idx <-> "ghost cell -> distant vertices" connectivity index
 *---------------------------------------------------------------------------*/

static void
_count_send_gcell_to_dist_vtx_connect(cs_mesh_t            *mesh,
                                      const cs_interface_set_t  *ifs,
                                      int                   rank_id,
                                      int                   tr_id,
                                      cs_lnum_t            *cell_faces_idx,
                                      cs_lnum_t            *cell_faces_lst,
                                      cs_lnum_t            *vtx_interface_idx,
                                      cs_lnum_t            *vtx_tag,
                                      cs_lnum_t            *gcell_dist_vtx_idx)
{
  cs_lnum_t id, cell_id, i_fac, i_vtx, i_loop;
  cs_lnum_t start_idx, end_idx, vtx_id, fac_id;

  cs_lnum_t n_loops = 0;
  cs_lnum_t n_added_vertices = 0;

  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t n_vertices = mesh->n_vertices;
  const cs_lnum_t *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *fac_vtx_lst = mesh->i_face_vtx_lst;

  const char err_corresp[]
    = N_("Repeated inconsistency in the halo construction.\n"
         "Several local points have the same distant correspondant;\n"
         "this is probably a side effect of a poor quality joining\n"
         "or periodicity.\n"
         "Coordinates of the first impacted point: [%12.5e, %12.5e %12.5e].");

  _define_vtx_interface_idx(ifs,
                            halo->c_domain_rank[rank_id],
                            tr_id,
                            n_vertices,
                            vtx_interface_idx);

  /* Count size of the connectivity and define index */

  if (mesh->halo_type == CS_HALO_STANDARD)
    n_loops = 1;
  else if (mesh->halo_type == CS_HALO_EXTENDED)
    n_loops = 2;

  for (i_loop = 0; i_loop < n_loops; i_loop++) {

    /* Define start and end idx */

    _get_start_end_idx(mesh,
                       halo->send_index,
                       halo->send_perio_lst,
                       rank_id,
                       tr_id,
                       i_loop,
                       &start_idx,
                       &end_idx);

    for (id = start_idx; id < end_idx; id++) {

      cell_id = halo->send_list[id];

      for (i_fac = cell_faces_idx[cell_id];
           i_fac < cell_faces_idx[cell_id+1];
           i_fac++) {

        fac_id = cell_faces_lst[i_fac];

        if (_test_loop_continues(mesh, fac_id) == true) {

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id+1];
               i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx];

            /* If vertex is on the interface for this rank and this
               transformation */

            n_added_vertices =  vtx_interface_idx[vtx_id+1]
                              - vtx_interface_idx[vtx_id];

            if (n_added_vertices > 0) {

              if (n_added_vertices > 1) {

                if (cs_glob_n_ranks > 1)
                  bft_printf("fac_num: %d (%llu)\n"
                             "vtx_num: %d (%llu) - n_added: %d\n",
                             fac_id+1,
                             (unsigned long long)(mesh->global_i_face_num[fac_id]),
                             vtx_id+1,
                             (unsigned long long)(mesh->global_vtx_num[vtx_id]),
                             n_added_vertices);
                else
                  bft_printf("fac_num: %d\n"
                             "vtx_num: %d - n_added: %d\n",
                             fac_id+1, vtx_id+1, n_added_vertices);
                bft_printf_flush();

                bft_error(__FILE__, __LINE__, 0, _(err_corresp),
                          mesh->vtx_coord[vtx_id*3],
                          mesh->vtx_coord[vtx_id*3+1],
                          mesh->vtx_coord[vtx_id*3+2]);

              }

              /* Add this vertex if not already checked */

              if (vtx_tag[vtx_id] != id) {
                vtx_tag[vtx_id] = id;
                gcell_dist_vtx_idx[id+1] += n_added_vertices;
              }

            }

          } /* End of loop on vertices */

        } /* Treatment only for implied faces */

      } /* End of loop on faces */

    } /* End of loop on ghost cells */

  } /* Loop on STANDARD and EXTENDED halo */

}

/*---------------------------------------------------------------------------
 * Fill the "ghost cells to distant vertices" connectivity.
 *
 * parameters:
 *   mesh               <-- pointer to cs_mesh_t structure
 *   ifs                <-- pointer to cs_interface_set_t structure
 *   rank_id            <-- rank number to work with
 *   tr_id              <-- transformation id to work with
 *   cell_faces_idx     <-- "cell -> faces" connectivity index
 *   cell_faces_lst     <-- "cell -> faces" connectivity list
 *   vtx_interface_idx  <-> index on vertices matching criteria
 *   dist_id_lst        <-> list of distant vertex numbers matching criteria
 *   counter            <-> counter on vertices
 *   vtx_tag            <-> tag array on vertices
 *   gcell_dist_vtx_idx <-> "ghost cell -> distant vertices" connectivity index
 *   gcell_dist_vtx_lst <-> "ghost cell -> distant vertices" connectivity list
 *---------------------------------------------------------------------------*/

static void
_fill_send_gcell_to_dist_vtx_connect(cs_mesh_t            *mesh,
                                     const cs_interface_set_t  *ifs,
                                     cs_lnum_t             rank_id,
                                     cs_lnum_t             tr_id,
                                     cs_lnum_t            *cell_faces_idx,
                                     cs_lnum_t            *cell_faces_lst,
                                     cs_lnum_t            *vtx_interface_idx,
                                     cs_lnum_t            *dist_id_lst,
                                     cs_lnum_t            *counter,
                                     cs_lnum_t            *vtx_tag,
                                     cs_lnum_t            *gcell_dist_vtx_idx,
                                     cs_lnum_t            *gcell_dist_vtx_lst)
{
  cs_lnum_t i, id, cell_id, i_fac, i_vtx, i_loop;
  cs_lnum_t shift, vtx_id, fac_id, start_idx, end_idx;

  cs_lnum_t n_loops = 0;
  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t n_vertices = mesh->n_vertices;
  const cs_lnum_t *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *fac_vtx_lst = mesh->i_face_vtx_lst;

  _define_dist_id_lst(ifs,
                      halo->c_domain_rank[rank_id],
                      tr_id,
                      n_vertices,
                      dist_id_lst,
                      vtx_interface_idx);

  /* Fill the "ghost cells to distant vertices" connectivity */

  if (mesh->halo_type == CS_HALO_STANDARD)
    n_loops = 1;
  else if (mesh->halo_type == CS_HALO_EXTENDED)
    n_loops = 2;

  for (i_loop = 0; i_loop < n_loops; i_loop++) {

    /* Define start and end idx */

    _get_start_end_idx(mesh,
                       halo->send_index,
                       halo->send_perio_lst,
                       rank_id,
                       tr_id,
                       i_loop,
                       &start_idx,
                       &end_idx);

    for (id = start_idx; id < end_idx; id++) {

      cell_id = halo->send_list[id];

      for (i_fac = cell_faces_idx[cell_id];
           i_fac < cell_faces_idx[cell_id+1];
           i_fac++) {

        fac_id = cell_faces_lst[i_fac];

        if (_test_loop_continues(mesh, fac_id) == true) {

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id + 1];
               i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx];

            /* If vertex is on the interface for this rank and periodicity */

            if (vtx_interface_idx[vtx_id+1] - vtx_interface_idx[vtx_id] > 0) {

              /* Add this vertex if nont already checked */

              if (vtx_tag[vtx_id] != id) { /* Add this vertex */

                vtx_tag[vtx_id] = id;

                for (i = vtx_interface_idx[vtx_id];
                     i < vtx_interface_idx[vtx_id+1];
                     i++) {

                  shift = gcell_dist_vtx_idx[id] + counter[id];
                  gcell_dist_vtx_lst[shift] = dist_id_lst[i];
                  counter[id] += 1;

                }

              }

            } /* If there is something to fill */

          } /* End of loop on vertices */

        } /* Treatment only for implied faces */

      } /* End of loop on faces */

    } /* End of loop on ghost cells */

  } /* Loop on STANDARD or EXTENDED halo */

}

/*---------------------------------------------------------------------------
 * Create a local "ghost cells -> distant vertices" connectivity for
 * send_halo cells.
 *
 * parameters:
 *   mesh                 <-- pointer to cs_mesh_t structure
 *   interface_set        <-- pointer to cs_interface_set_t structure
 *   cell_faces_idx       <-- "cell -> faces" connectivity index
 *   cell_faces_lst       <-- "cell -> faces" connectivity list
 *   p_gcell_dist_vtx_idx --> "ghost cell -> distant vertices" connect. index
 *   p_gcell_dist_vtx_lst --> "ghost cell -> distant vertices" connect. list
 *---------------------------------------------------------------------------*/

static void
_create_send_gcell_vtx_connect(cs_mesh_t           *mesh,
                               cs_interface_set_t  *interface_set,
                               cs_lnum_t           *cell_faces_idx,
                               cs_lnum_t           *cell_faces_lst,
                               cs_lnum_t           *p_gcell_dist_vtx_idx[],
                               cs_lnum_t           *p_gcell_dist_vtx_lst[])
{
  cs_lnum_t i, id, rank_id;

  cs_lnum_t *gcell_dist_vtx_idx = NULL, *gcell_dist_vtx_lst = NULL;
  cs_lnum_t *vtx_interface_idx = NULL;
  cs_lnum_t *dist_id_lst = NULL;
  cs_lnum_t *vtx_tag = NULL;
  cs_lnum_t *counter = NULL;

  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t max_lst_size = _get_list_buffer_size(interface_set);
  const cs_lnum_t n_ghost_cells = halo->n_send_elts[CS_HALO_EXTENDED];
  const cs_lnum_t n_vertices = mesh->n_vertices;
  const cs_lnum_t n_c_domains = halo->n_c_domains;
  const cs_lnum_t tr_index_size = mesh->n_transforms + 1;

  if (n_ghost_cells == 0)
    return;

  /* Allocate and initialize buffers */

  BFT_MALLOC(gcell_dist_vtx_idx, n_ghost_cells + 1, cs_lnum_t);
  BFT_MALLOC(counter, n_ghost_cells, cs_lnum_t);

  gcell_dist_vtx_idx[0] = 0;
  for (i = 0; i < n_ghost_cells; i++) {
    gcell_dist_vtx_idx[i+1] = 0;
    counter[i] = 0;
  }

  BFT_MALLOC(vtx_tag, n_vertices, cs_lnum_t);

  for (i = 0; i < n_vertices; i++)
    vtx_tag[i] = -1;

  BFT_MALLOC(vtx_interface_idx, n_vertices + 1, cs_lnum_t);
  BFT_MALLOC(dist_id_lst, max_lst_size, cs_lnum_t);

  for (i = 0; i < max_lst_size; i++)
    dist_id_lst[i] = -1;

  /* Loop on each rank belonging to send_halo.
     Create a vertex to ghost cells connectivity for each rank */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    /* Define gcell_dist_vtx_idx */

    for (id = 0; id < tr_index_size; id++)
      _count_send_gcell_to_dist_vtx_connect(mesh,
                                            interface_set,
                                            rank_id,
                                            id,
                                            cell_faces_idx,
                                            cell_faces_lst,
                                            vtx_interface_idx,
                                            vtx_tag,
                                            gcell_dist_vtx_idx);

  } /* End of loop on ranks */

  /* Create gcell_dist_vtx_idx */

  for (i = 0; i < n_ghost_cells; i++)
    gcell_dist_vtx_idx[i+1] += gcell_dist_vtx_idx[i];

  BFT_MALLOC(gcell_dist_vtx_lst, gcell_dist_vtx_idx[n_ghost_cells], cs_lnum_t);

  for (i = 0; i < n_vertices; i++)
    vtx_tag[i] = -1;

  cs_interface_set_add_match_ids(interface_set);

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    for (id = 0; id < tr_index_size; id++)
      _fill_send_gcell_to_dist_vtx_connect(mesh,
                                           interface_set,
                                           rank_id,
                                           id,
                                           cell_faces_idx,
                                           cell_faces_lst,
                                           vtx_interface_idx,
                                           dist_id_lst,
                                           counter,
                                           vtx_tag,
                                           gcell_dist_vtx_idx,
                                           gcell_dist_vtx_lst);

  } /* End of loop on ranks */

  cs_interface_set_free_match_ids(interface_set);

  BFT_FREE(counter);
  BFT_FREE(vtx_tag);
  BFT_FREE(vtx_interface_idx);
  BFT_FREE(dist_id_lst);

  *p_gcell_dist_vtx_idx = gcell_dist_vtx_idx;
  *p_gcell_dist_vtx_lst = gcell_dist_vtx_lst;

}

/*---------------------------------------------------------------------------
 * Send "ghost cells to distant_num vertices" connectivity on communicating
 * ranks and receive the same kind of connectivity from distant ranks.
 *
 * parameters:
 *   mesh                     <--  pointer to cs_mesh_t structure
 *   send_gcell_dist_vtx_idx  -->  "ghost cell -> distant vertices" index
 *   send_gcell_dist_vtx_lst  -->  "ghost cell -> distant vertices" list
 *   p_gcell_dist_vtx_idx     <--  "ghost cell -> distant vertices" index
 *   p_gcell_dist_vtx_lst     <--  "ghost cell -> distant vertices" list
 *---------------------------------------------------------------------------*/

static void
_exchange_gcell_vtx_connect(cs_mesh_t  *mesh,
                            cs_lnum_t  *send_gcell_dist_vtx_idx,
                            cs_lnum_t  *send_gcell_dist_vtx_lst,
                            cs_lnum_t  *p_gcell_dist_vtx_idx[],
                            cs_lnum_t  *p_gcell_dist_vtx_lst[])
{
  cs_lnum_t i, j, rank_id;
  cs_lnum_t send_start_idx, send_end_idx, start_idx, end_idx;
  cs_lnum_t n_send_elts, n_recv_elts;

  cs_lnum_t send_buffer_size = 0;

  cs_lnum_t *send_idx_buffer = NULL;
  cs_lnum_t *gcell_dist_vtx_idx = NULL, *gcell_dist_vtx_lst = NULL;
  cs_lnum_t *send_buffer = NULL, *recv_buffer = NULL;

  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;
  const cs_lnum_t n_c_domains = halo->n_c_domains;
  const cs_lnum_t n_ghost_cells = halo->n_elts[CS_HALO_EXTENDED];

#if defined(HAVE_MPI)
  MPI_Status  status;
#endif

  /* Allocate buffers */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {
    if (halo->c_domain_rank[rank_id] != local_rank) {
      n_send_elts = halo->send_index[2*rank_id+2]- halo->send_index[2*rank_id];
      send_buffer_size = CS_MAX(send_buffer_size, n_send_elts);
    }
  }

  BFT_MALLOC(send_idx_buffer, send_buffer_size, cs_lnum_t);
  BFT_MALLOC(gcell_dist_vtx_idx, n_ghost_cells + 1, cs_lnum_t);

  for (i = 0; i < n_ghost_cells + 1; i++)
    gcell_dist_vtx_idx[i] = 0;

  /* Exchange sizes to define gcell_dist_vtx_idx */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    recv_buffer = &(gcell_dist_vtx_idx[1 + halo->index[2*rank_id]]);

    if (halo->c_domain_rank[rank_id] != local_rank) {

      /* Fill send buffer */

      for (i = halo->send_index[2*rank_id], j = 0;
           i < halo->send_index[2*rank_id+2]; i++, j++)
        send_idx_buffer[j] = (  send_gcell_dist_vtx_idx[i+1]
                              - send_gcell_dist_vtx_idx[i]);

      n_send_elts =  halo->send_index[2*rank_id+2] - halo->send_index[2*rank_id];
      n_recv_elts =  halo->index[2*rank_id+2] - halo->index[2*rank_id];

#if defined(HAVE_MPI)
      MPI_Sendrecv(&(send_idx_buffer[0]), n_send_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], local_rank,
                   &(recv_buffer[0]), n_recv_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                   cs_glob_mpi_comm, &status);
#endif

    } /* If rank != local_rank */

    else {

      for (i = halo->send_index[2*rank_id], j = 0;
           i < halo->send_index[2*rank_id+2]; i++, j++)
        recv_buffer[j] = (  send_gcell_dist_vtx_idx[i+1]
                          - send_gcell_dist_vtx_idx[i]);

    } /* rank == local_rank */

  } /* End of loop on if_ranks */

  BFT_FREE(send_idx_buffer);

  /* Define index */

  for (i = 0; i < n_ghost_cells; i++)
    gcell_dist_vtx_idx[i+1] += gcell_dist_vtx_idx[i];

  BFT_MALLOC(gcell_dist_vtx_lst, gcell_dist_vtx_idx[n_ghost_cells], cs_lnum_t);

  /* Exchange lists to define gcell_dist_vtx_lst */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    /* Exchange conectivity list */

    send_start_idx = send_gcell_dist_vtx_idx[halo->send_index[2*rank_id]];
    send_end_idx = send_gcell_dist_vtx_idx[halo->send_index[2*rank_id+2]];
    n_send_elts = send_end_idx - send_start_idx;
    send_buffer = &(send_gcell_dist_vtx_lst[send_start_idx]);

    start_idx = gcell_dist_vtx_idx[halo->index[2*rank_id]];
    end_idx = gcell_dist_vtx_idx[halo->index[2*rank_id+2]];
    n_recv_elts = end_idx - start_idx;
    recv_buffer = &(gcell_dist_vtx_lst[start_idx]);

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

  *p_gcell_dist_vtx_idx = gcell_dist_vtx_idx;
  *p_gcell_dist_vtx_lst = gcell_dist_vtx_lst;

}

/*---------------------------------------------------------------------------
 * Check mesh->i_face_cells array for ghost cells in standard halo.
 *
 * parameters:
 *   mesh   <->  pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_check_i_face_cells(cs_mesh_t  *mesh)
{
  int  i;

  cs_halo_t  *halo = mesh->halo;

  const cs_lnum_t n_c_domains = halo->n_c_domains;

  for (i = 0; i < mesh->n_i_faces; i++) {

    if (mesh->i_face_cells[i][0] < 0) {
      if (n_c_domains > 1)
        bft_error(__FILE__, __LINE__, 0,
                  " Error detected in interior face -> cells connectivity.\n"
                  " Face %d (%llu) has an incomplete connectivity.\n"
                  " Cell1: %d - Cell2: %d (%llu)",
                  i+1, (unsigned long long)(mesh->global_i_face_num[i]),
                  mesh->i_face_cells[i][0],
                  mesh->i_face_cells[i][1],
                  (unsigned long long)(mesh->global_cell_num[mesh->i_face_cells[i][1]]));
      else /* Serial run */
        bft_error(__FILE__, __LINE__, 0,
                  " Error detected in interior face -> cells connectivity.\n"
                  " Face %d has an incomplete connectivity.\n"
                  " Cell1: %d - Cell2: %d",
                  i+1, mesh->i_face_cells[i][0], mesh->i_face_cells[i][1]);
    }

    if (mesh->i_face_cells[i][1] < 0) {
      if (n_c_domains > 1)
        bft_error(__FILE__, __LINE__, 0,
                  " Error detected in interior face -> cells connectivity.\n"
                  " Face %d (%llu) has an incomplete connectivity.\n"
                  " Cell1: %d (%llu) - Cell2: %d",
                  i+1, (unsigned long long)(mesh->global_i_face_num[i]),
                  mesh->i_face_cells[i][0],
                  (unsigned long long)(mesh->global_cell_num[mesh->i_face_cells[i][0]]),
                  mesh->i_face_cells[i][1]);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Error detected in interior face -> cells connectivity.\n"
                  " Face %d has an incomplete connectivity.\n"
                  " Cell1: %d - Cell2: %d",
                  i+1, mesh->i_face_cells[i][0], mesh->i_face_cells[i][1]);
    }

  }

}

/*---------------------------------------------------------------------------
 * Tag local face by its distant face id and reset previous tag
 *
 * parameters:
 *   prev      <--  previous start index
 *   start     <--  start index in lnum and dnum
 *   end       <--  end index in lnum and dnum
 *   l_ids     <--  elt_ids array in cs_interface_t struct.
 *   d_ids     <--  match_ids array in cs_interface_t struct.
 *   l2d_fids  <->  tag on local faces
 *---------------------------------------------------------------------------*/

static void
_tag_interface_faces(int              prev,
                     int              start,
                     int              end,
                     const cs_lnum_t  l_ids[],
                     const cs_lnum_t  d_ids[],
                     cs_lnum_t        l2d_fids[])
{
  cs_lnum_t  i;

  assert(l2d_fids != NULL);

  if (prev > -1) /* Reset previous face num. */
    for (i = prev; i < start; i++)
      l2d_fids[l_ids[i]] = -1;

  for (i = start; i < end; i++)
    l2d_fids[l_ids[i]] = d_ids[i];

}

/*---------------------------------------------------------------------------
 * Update mesh->i_face_cells array for ghost cells in standard halo.
 *
 * parameters:
 *   mesh            <->  pointer to cs_mesh_t structure
 *   face_ifs        <--  faces interface
 *   cell_faces_idx  <--  "cell -> faces" connectivity index
 *   cell_faces_lst  <--  "cell -> faces" connectivity list
 *---------------------------------------------------------------------------*/

static void
_update_i_face_cells(cs_mesh_t           *mesh,
                     cs_interface_set_t  *face_ifs,
                     cs_lnum_t           *cell_faces_idx,
                     cs_lnum_t           *cell_faces_lst)
{
  int  i, j, gcell_id, face_id, if_id, rank, shift, index_end;
  int  tr_id, start, end, length, n_interfaces, mpi_count, distant_rank;

  int  local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;
  cs_halo_t  *halo = mesh->halo;
  cs_lnum_t *gcell_face_count = NULL, *l2d_fids = NULL;
  cs_lnum_t *send_buffer = NULL, *recv_buffer = NULL, *buffer = NULL;
  cs_lnum_t *send_shift = NULL, *recv_shift = NULL, *halo2ifs_rank = NULL;

  const int  n_init_perio = mesh->n_init_perio;
  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = (cs_glob_rank_id == -1) ? 0 : cs_glob_rank_id;
  const int  n_c_domains = halo->n_c_domains;
  const int  *s_index = halo->send_index;
  const int  *s_gcell_ids = halo->send_list;
  const cs_interface_t  *interface = NULL;
  const cs_lnum_t  *loc_ids = NULL;
  const cs_lnum_t  *dist_ids = NULL;
  const cs_lnum_t  *tr_index = NULL;

#if defined(HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  if (n_c_domains*2 > 128) {
    BFT_MALLOC(request, n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, n_c_domains*2, MPI_Status);
  }
#endif

  /* Get an interface set struct. on faces */

  assert(face_ifs != NULL);
  n_interfaces = cs_interface_set_size(face_ifs);

  /* Allocate auxiliary buffers */

  BFT_MALLOC(l2d_fids, mesh->n_i_faces, cs_lnum_t);
  BFT_MALLOC(halo2ifs_rank, n_c_domains, cs_lnum_t);

  /* If there is an extended neighborhood, n_c_domains can be different
     from n_interfaces (c_rank with only vertices on the interface) */

  for (i = 0; i < n_c_domains; i++)
    halo2ifs_rank[i] = -1; /* No face interface between ranks */

  index_end = 1; /* Only purely parralel faces are taking into account */
  index_end += 2*n_init_perio ; /* 2 transform. by initial periodicity */

  /* Identify interface id and communicating rank */

  for (if_id = 0; if_id < n_interfaces; if_id++) {

    interface = cs_interface_set_get(face_ifs, if_id);
    distant_rank = cs_interface_rank(interface);

    for (i = 0; i < n_c_domains; i++)
      if (halo->c_domain_rank[i] == distant_rank)
        break;
    assert(i < n_c_domains);
    halo2ifs_rank[i] = if_id;

  } /* End of loop on interfaces */

  /* First exchange number of faces linked to each ghost cells */

  BFT_MALLOC(send_buffer, halo->n_send_elts[CS_HALO_EXTENDED], cs_lnum_t);
  BFT_MALLOC(send_shift, n_c_domains + 1, cs_lnum_t);

  send_shift[0] = 0;
  for (i = 0; i < halo->n_send_elts[CS_HALO_EXTENDED]; i++)
    send_buffer[i] = 0;

  /* Loop on communicating ranks to build send_buffer */

  cs_interface_set_add_match_ids(face_ifs);

  for (rank = 0, shift = 0; rank < n_c_domains; rank++) {

    if_id = halo2ifs_rank[rank];

    if (if_id > -1) { /* This communicating rank shares elements
                         in the standard neighborhood */

      interface = cs_interface_set_get(face_ifs, if_id);
      loc_ids = cs_interface_get_elt_ids(interface);
      dist_ids = cs_interface_get_match_ids(interface);
      tr_index = cs_interface_get_tr_index(interface);

      assert(cs_interface_rank(interface) == halo->c_domain_rank[rank]);

      /* Initalize l2d_fids */

      for (i = 0; i < mesh->n_i_faces; i++)
        l2d_fids[i] = -1;

      for (tr_id = 0; tr_id < index_end; tr_id++) {

        _get_start_end_idx(mesh, s_index, halo->send_perio_lst,
                           rank, tr_id, 0, /* STANDARD */
                           &start, &end);

        if (tr_id == 0) { /* Purely parallel faces */

          if (tr_index != NULL) {
            _tag_interface_faces(-1, 0, tr_index[1],
                                 loc_ids, dist_ids, l2d_fids);
            assert(n_init_perio > 0);
          }
          else {
            _tag_interface_faces(-1, 0, cs_interface_size(interface),
                                 loc_ids, dist_ids, l2d_fids);
            assert(n_init_perio == 0);
          }

        }
        else /* Periodic faces */
          _tag_interface_faces(tr_index[tr_id-1],
                               tr_index[tr_id],
                               tr_index[tr_id+1],
                               loc_ids, dist_ids, l2d_fids);

        for (i = start; i < end; i++) {

          gcell_id = s_gcell_ids[i];
          for (j = cell_faces_idx[gcell_id];
               j < cell_faces_idx[gcell_id+1];
               j++) {

            face_id = cell_faces_lst[j];
            if (face_id > -1) {
              if (l2d_fids[face_id] > -1) {
                shift++;
                send_buffer[i] += 1;
              }
            }

          } /* End of loop cell -> face connectivity */

        } /* End of loop on standard neighborhood for this tr_id and rank */

      } /* End of loop on tr_ids */

    } /* if_id > -1 */

    send_shift[rank+1] = shift;

  } /* End of loop on ranks */

  /* Exchange data over the ranks  => build face_count */

  BFT_MALLOC(gcell_face_count, halo->n_elts[CS_HALO_EXTENDED], cs_lnum_t);

  for (i = 0; i < halo->n_elts[CS_HALO_EXTENDED]; i++)
    gcell_face_count[i] = 0;

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    /* Receive data from distant ranks */

    mpi_count = 0;

    for (rank = 0; rank < n_c_domains; rank++) {

      start = halo->index[2*rank];
      length = halo->index[2*rank+1] - halo->index[2*rank];

      if (halo->c_domain_rank[rank] != local_rank) {

        buffer = gcell_face_count + start ;

        MPI_Irecv(buffer,
                  length,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank],
                  halo->c_domain_rank[rank],
                  cs_glob_mpi_comm,
                  &(request[mpi_count++]));

      }
      else
        local_rank_id = rank;

    }

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(cs_glob_mpi_comm);

    /* Send data to distant ranks */

    for (rank = 0; rank < n_c_domains; rank++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank] != local_rank) {

        start = s_index[2*rank];
        length = s_index[2*rank + 1] - s_index[2*rank];
        buffer = send_buffer + start;

        MPI_Isend(buffer,
                  length,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank],
                  local_rank,
                  cs_glob_mpi_comm,
                  &(request[mpi_count++]));

      }

    }

    /* Wait for all exchanges */

    MPI_Waitall(mpi_count, request, status);
  }

#endif /* defined(HAVE_MPI) */

  /* Copy local values in case of periodicity */

  if (n_init_perio > 0 && local_rank_id > -1) {

    buffer = gcell_face_count + halo->index[2*local_rank_id];
    start = s_index[2*local_rank_id];
    length = s_index[2*local_rank_id + 1] - s_index[2*local_rank_id];

    for (i = 0; i < length; i++)
      buffer[i] = send_buffer[start + i];

  }

  /* Build recv_shift */

  BFT_MALLOC(recv_shift, n_c_domains + 1, cs_lnum_t);

  recv_shift[0] = 0;
  for (rank = 0; rank < n_c_domains; rank++) {
    recv_shift[rank+1] = 0;
    for (i = halo->index[2*rank]; i < halo->index[2*rank+1]; i++)
      recv_shift[rank+1] += gcell_face_count[i];
  }

  for (rank = 0; rank < n_c_domains; rank++)
    recv_shift[rank+1] += recv_shift[rank];

  BFT_REALLOC(send_buffer, send_shift[n_c_domains], cs_lnum_t);
  BFT_MALLOC(recv_buffer, recv_shift[n_c_domains], cs_lnum_t);

  /* Build send_buffer */

  for (rank = 0; rank < n_c_domains; rank++) {

    if_id = halo2ifs_rank[rank];

    if (if_id > -1) {

      interface = cs_interface_set_get(face_ifs, if_id);
      loc_ids = cs_interface_get_elt_ids(interface);
      dist_ids = cs_interface_get_match_ids(interface);
      tr_index = cs_interface_get_tr_index(interface);
      shift = send_shift[rank];

      /* Initalize l2d_fids */

      for (i = 0; i < mesh->n_i_faces; i++)
        l2d_fids[i] = -1;

      for (tr_id = 0; tr_id < index_end; tr_id++) {

        _get_start_end_idx(mesh, s_index, halo->send_perio_lst,
                           rank, tr_id, 0, /* STANDARD */
                           &start, &end);

        if (tr_id == 0) { /* Purely parallel faces */

          if (tr_index != NULL) {
            _tag_interface_faces(-1, 0, tr_index[1],
                                 loc_ids, dist_ids, l2d_fids);
            assert(n_init_perio > 0);
          }
          else {
            _tag_interface_faces(-1, 0, cs_interface_size(interface),
                                 loc_ids, dist_ids, l2d_fids);
            assert(n_init_perio == 0);
          }

        }
        else /* Periodic faces */
          _tag_interface_faces(tr_index[tr_id-1],
                               tr_index[tr_id],
                               tr_index[tr_id+1],
                               loc_ids, dist_ids, l2d_fids);

        for (i = start; i < end; i++) {

          gcell_id = s_gcell_ids[i];

          for (j = cell_faces_idx[gcell_id];
               j < cell_faces_idx[gcell_id+1];
               j++) {

            face_id = cell_faces_lst[j];
            if (face_id > -1)
              if (l2d_fids[face_id] > -1)
                send_buffer[shift++] = l2d_fids[face_id];

          } /* End of loop cell -> face connectivity */

        }

      } /* End of loop on standard neighborhood for this communicating rank */

    } /* if_id > -1 */

  } /* End of loop on ranks */

  /* Exchange face ids */

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Exchange data over the ranks  => build face_count */

    /* Receive data from distant ranks */

    mpi_count = 0;

    for (rank = 0; rank < n_c_domains; rank++) {

      start = recv_shift[rank];
      length = recv_shift[rank+1] - recv_shift[rank];

      if (halo->c_domain_rank[rank] != local_rank) {

        buffer = recv_buffer + start ;

        MPI_Irecv(buffer,
                  length,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank],
                  halo->c_domain_rank[rank],
                  cs_glob_mpi_comm,
                  &(request[mpi_count++]));

      }

    }

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(cs_glob_mpi_comm);

    /* Send data to distant ranks */

    for (rank = 0; rank < n_c_domains; rank++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank] != local_rank) {

        start = send_shift[rank];
        length = send_shift[rank + 1] - send_shift[rank];
        buffer = send_buffer + start;

        MPI_Isend(buffer,
                  length,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank],
                  local_rank,
                  cs_glob_mpi_comm,
                  &(request[mpi_count++]));

      }

    }

    /* Wait for all exchanges */

    MPI_Waitall(mpi_count, request, status);
  }

#endif /* defined(HAVE_MPI) */

  /* Copy local values in case of periodicity */

  if (n_init_perio > 0 && local_rank_id > -1) {

    buffer = recv_buffer + recv_shift[local_rank_id];
    start = send_shift[local_rank_id];
    length = send_shift[local_rank_id + 1] - send_shift[local_rank_id];

    for (i = 0; i < length; i++)
      buffer[i] = send_buffer[start + i];

  }

  /* Update face-> cells connect */

  for (rank = 0; rank < n_c_domains; rank++) {

    shift = recv_shift[rank];

    for (i = halo->index[2*rank]; i < halo->index[2*rank+1]; i++) {

      for (j = 0; j < gcell_face_count[i]; j++) {

        face_id = recv_buffer[shift++];

        /* In case of periodicity, a ghost cell may appear several times
           in the same halo. So we can have a face->cell connect up-to-date
           when it is not the first pass */

        if (mesh->i_face_cells[face_id][0] == -1)
          mesh->i_face_cells[face_id][0] = mesh->n_cells + i;
        else if (mesh->i_face_cells[face_id][1] == -1)
          mesh->i_face_cells[face_id][1] = mesh->n_cells + i;

      } /* End of loop on related faces */

    } /* End of loop on standard neighborhood for this communicating rank */

  } /* End of loop on ranks */

  /* Free memory */

  cs_interface_set_free_match_ids(face_ifs);

  BFT_FREE(recv_buffer);
  BFT_FREE(send_buffer);
  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);
  BFT_FREE(gcell_face_count);
  BFT_FREE(halo2ifs_rank);
  BFT_FREE(l2d_fids);

#if defined(HAVE_MPI)
  if (request != _request) {
    BFT_FREE(request);
    BFT_FREE(status);
  }
#endif

  /* Sanity check */

  _check_i_face_cells(mesh);
}

#if 0 /* TODO: check algorithm (deadlock on BG/L on one test case) */

/*---------------------------------------------------------------------------
 * Clean a halo i.e. delete rank(s) with no ghost cells.
 *
 * This case happens when an interface has only extended vertices
 * and the halo is standard.
 *
 * parameters:
 *   mesh <--  pointer to a mesh structure
 *---------------------------------------------------------------------------*/

static void
_clean_halo(cs_mesh_t  *mesh)
{
  cs_lnum_t rank_id, i, j;

  cs_halo_t  *halo = mesh->halo;

  cs_lnum_t n_c_domains = halo->n_c_domains;
  cs_lnum_t n_real_c_domains = 0;
  cs_lnum_t counter = 0;
  cs_lnum_t *new_c_domain_rank = NULL;
  cs_lnum_t *new_perio_lst = NULL;
  cs_lnum_t *new_index = NULL;

  const cs_lnum_t n_transforms = mesh->n_transforms;

  /* Is there something to do ? */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++)
    if (halo->send_index[2*rank_id+2] - halo->send_index[2*rank_id] > 0)
      n_real_c_domains++;

  if (n_real_c_domains == n_c_domains)
    return;

  /* halo->index, halo->perio_lst and n_c_domains need an update */

  BFT_MALLOC(new_c_domain_rank, n_real_c_domains, cs_lnum_t);
  BFT_MALLOC(new_index, 2*n_real_c_domains+1, cs_lnum_t);

  if (n_transforms > 0)
    BFT_MALLOC(new_perio_lst, 4*n_transforms*n_real_c_domains, cs_lnum_t);

  /* Define the new buffers */

  new_index[0] = 0;

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    if (halo->send_index[2*rank_id+2] - halo->send_index[2*rank_id] > 0) {

      new_c_domain_rank[counter] = halo->c_domain_rank[rank_id];
      new_index[2*counter+1] = halo->send_index[2*rank_id+1];
      new_index[2*counter+2] = halo->send_index[2*rank_id+2];

      for (i = 0; i < n_transforms; i++)
        for (j = 0; j < 4; j++)
          new_perio_lst[4*counter + j + 4*n_real_c_domains*i]
            = halo->send_perio_lst[4*rank_id + j + 4*n_c_domains*i];

      counter++;

    } /* If there are elements for this rank */

  } /* End of loop on frontier ranks */

  /* Replace halo->send_perio_lst, halo->send_index and
     halo->c_domain_rank by new ones */

  BFT_FREE(halo->c_domain_rank);
  BFT_FREE(halo->send_index);

  if (n_transforms > 0)
    BFT_FREE(halo->send_perio_lst);

  halo->n_c_domains = n_real_c_domains;
  halo->c_domain_rank = new_c_domain_rank;

  halo->send_index = new_index;
  if (n_transforms > 0)
    halo->send_perio_lst = new_perio_lst;

  /* Reallocation of halo's buffers */

  BFT_REALLOC(halo->index, 2*n_real_c_domains+1, cs_lnum_t);
  BFT_REALLOC(halo->perio_lst, 4*n_transforms*n_real_c_domains, cs_lnum_t);

}

#endif /* #if 0 */

/*----------------------------------------------------------------------------
 * Define cell -> internal faces connectivity for ghost cells.
 *
 * Treatment of parallel and/or periodic halos for standard or extended
 * ghost cells according to halo type building option.
 *
 * parameters:
 *   mesh             <-- pointer to cs_mesh_t structure.
 *   interface_set    <-- pointer to cs_interface_set structure.
 *   p_cell_faces_idx --> pointer to the connectivity index
 *   p_cell_faces_lst --> pointer to the connectivity list
 *----------------------------------------------------------------------------*/

static void
_create_gcell_faces_connect(cs_mesh_t                 *mesh,
                            const cs_interface_set_t  *vertex_ifs,
                            cs_lnum_t                 *p_cell_faces_idx[],
                            cs_lnum_t                 *p_cell_faces_lst[])
{
  cs_lnum_t i, fac_id, i_vtx, id1, id2, shift, vtx_id;

  cs_lnum_t *vtx_tag = NULL;
  cs_lnum_t *cell_buffer = NULL, *cell_tag = NULL, *counter = NULL;
  cs_lnum_t *cell_faces_idx = NULL;
  cs_lnum_t *cell_faces_lst = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_2_t *face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const cs_lnum_t *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *fac_vtx_lst = mesh->i_face_vtx_lst;

  *p_cell_faces_idx = cell_faces_idx;
  *p_cell_faces_lst = cell_faces_lst;

  if (vertex_ifs == NULL)
    return;

  BFT_MALLOC(cell_faces_idx, n_cells+1, cs_lnum_t);
  BFT_MALLOC(cell_buffer, 2*n_cells, cs_lnum_t);

  cell_tag = &(cell_buffer[0]);
  counter = &(cell_buffer[n_cells]);

  cell_faces_idx[0] = 0;
  for (i = 0; i < n_cells; i++) {
    cell_faces_idx[i+1] = 0;
    cell_tag[i] = -1;
  }

  assert(sizeof(cs_lnum_t) == sizeof(cs_lnum_t));

  _get_vertex_tag(mesh->n_vertices, vertex_ifs, &vtx_tag);

  for (fac_id = 0; fac_id < n_i_faces; fac_id++) {

    for (i_vtx = fac_vtx_idx[fac_id];
         i_vtx < fac_vtx_idx[fac_id + 1];
         i_vtx++) {

      vtx_id = fac_vtx_lst[i_vtx];

      if (vtx_tag[vtx_id] == 1) {

        id1 = face_cells[fac_id][0];
        id2 = face_cells[fac_id][1];

        if (id1 < 0) {
          if (cell_tag[id2] != fac_id) {
            cell_tag[id2] = fac_id;
            cell_faces_idx[id2 + 1] += 1;
          }
        }

        if (id2 < 0) {
          if (cell_tag[id1] != fac_id) {
            cell_tag[id1] = fac_id;
            cell_faces_idx[id1 + 1] += 1;
          }
        }

        if (mesh->halo_type == CS_HALO_EXTENDED) {
          if (id1 >= 0 && id2 >= 0) {

            if (cell_tag[id1] != fac_id) {
              cell_tag[id1] = fac_id;
              cell_faces_idx[id1 + 1] += 1;
            }

            if (cell_tag[id2] != fac_id) {
              cell_tag[id2] = fac_id;
              cell_faces_idx[id2 + 1] += 1;
            }

          }
        }

      }

    } /* End of loop on vertices */

  } /* End of loop on internal faces */

  /* Build index */

  for (i = 0; i < n_cells; i++) {
    cell_faces_idx[i+1] += cell_faces_idx[i];
    counter[i] = 0;
    cell_tag[i] = -1;
  }

  BFT_MALLOC(cell_faces_lst, cell_faces_idx[n_cells], cs_lnum_t);

  for (fac_id = 0; fac_id < n_i_faces; fac_id++) {

    for (i_vtx = fac_vtx_idx[fac_id];
         i_vtx < fac_vtx_idx[fac_id + 1];
         i_vtx++) {

      vtx_id = fac_vtx_lst[i_vtx];

      if (vtx_tag[vtx_id] == 1) {

        id1 = face_cells[fac_id][0];
        id2 = face_cells[fac_id][1];

        if (id1 < 0) {
          if (cell_tag[id2] != fac_id) {

            cell_tag[id2] = fac_id;
            shift = cell_faces_idx[id2] + counter[id2];
            cell_faces_lst[shift] = fac_id;
            counter[id2] += 1;

          }
        }

        if (id2 < 0) {
          if (cell_tag[id1] != fac_id) {

            cell_tag[id1] = fac_id;
            shift = cell_faces_idx[id1] + counter[id1];
            cell_faces_lst[shift] = fac_id;
            counter[id1] += 1;

          }
        }

        if (mesh->halo_type == CS_HALO_EXTENDED) {
          if (id1 >= 0 && id2 >= 0) {

            if (cell_tag[id1] != fac_id) {
              cell_tag[id1] = fac_id;
              shift = cell_faces_idx[id1] + counter[id1];
              cell_faces_lst[shift] = fac_id;
              counter[id1] += 1;
            }

            if (cell_tag[id2] != fac_id) {
              cell_tag[id2] = fac_id;
              shift = cell_faces_idx[id2] + counter[id2];
              cell_faces_lst[shift] = fac_id;
              counter[id2] += 1;
            }

          }
        }

      }

    } /* End of loop on vertices */

  } /* End of loop on internal faces */

  BFT_FREE(vtx_tag);
  BFT_FREE(cell_buffer);

  *p_cell_faces_idx = cell_faces_idx;
  *p_cell_faces_lst = cell_faces_lst;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define halo structures for internal and distant ghost cells.
 *
 * parameters:
 *   mesh             <--  pointer to cs_mesh_t structure
 *   face_ifs         <--  pointer to faces interfaces
 *   vertex_ifs       <--  pointer to vertex interfaces
 *   p_gcell_vtx_idx  -->  pointer to the connectivity index
 *   p_gcell_vtx_lst  -->  pointer to the connectivity list
 *---------------------------------------------------------------------------*/

void
cs_mesh_halo_define(cs_mesh_t           *mesh,
                    cs_interface_set_t  *face_ifs,
                    cs_interface_set_t  *vertex_ifs,
                    cs_lnum_t           *p_gcell_vtx_idx[],
                    cs_lnum_t           *p_gcell_vtx_lst[])
{
  cs_lnum_t *send_gcell_vtx_idx = NULL, *send_gcell_vtx_lst = NULL;
  cs_lnum_t *gcell_vtx_idx = NULL, *gcell_vtx_lst = NULL;
  cs_lnum_t *gcell_faces_idx = NULL, *gcell_faces_lst = NULL;
  cs_halo_t  *halo = mesh->halo;

  halo->n_local_elts = mesh->n_cells;

  /*  Define cell -> internal faces connectivity for ghost cells */

  _create_gcell_faces_connect(mesh,
                              vertex_ifs,
                              &gcell_faces_idx,
                              &gcell_faces_lst);

  /* Fill cs_halo_t structure for send_halo  */

  if (mesh->verbosity > 0) {
    bft_printf(_("    Local halo definition\n"));
    bft_printf_flush();
  }

  _fill_send_halo(mesh,
                  vertex_ifs,
                  gcell_faces_idx,
                  gcell_faces_lst);

#if 0
  /* Clean cs_halo_t structure.
     Remove communicating ranks with no ghost cells */

  bft_printf(_("    Halo cleaning\n"));
  bft_printf_flush();

  _clean_halo(mesh);
#endif

  /* Fill cs_halo_t structure for halo.
     We use the data from send_halo structure */

  if (mesh->verbosity > 0) {
    bft_printf(_("    Distant halo creation\n"));
    bft_printf_flush();
  }

  _fill_halo(mesh);

  /* Update mesh structure elements bound to halo management */

  mesh->n_ghost_cells = halo->n_elts[CS_HALO_EXTENDED];
  mesh->n_cells_with_ghosts = mesh->n_cells + mesh->n_ghost_cells;

  /* Update connectivity between internal faces and cells */

  if (cs_mesh_n_g_ghost_cells(mesh) > 0) {

    /* Create local ghost connectivity for send_halo cells and send it.
       Receive ghost cells connectivity for halo cells. */

    _create_send_gcell_vtx_connect(mesh,
                                   vertex_ifs,
                                   gcell_faces_idx,
                                   gcell_faces_lst,
                                   &send_gcell_vtx_idx,
                                   &send_gcell_vtx_lst);

    _exchange_gcell_vtx_connect(mesh,
                                send_gcell_vtx_idx,
                                send_gcell_vtx_lst,
                                &gcell_vtx_idx,
                                &gcell_vtx_lst);

    /* Free memory */

    BFT_FREE(send_gcell_vtx_idx);
    BFT_FREE(send_gcell_vtx_lst);

    /* Define mesh->i_face_cells array for ghost cells in standard halo and
       also ghost cells to ghost cells connectivity for standard and extended
       halo if necessary */

    if (mesh->verbosity > 0) {
      bft_printf(_("    Updating the faces -> cells connectivity\n"));
      bft_printf_flush();
    }

    _update_i_face_cells(mesh, face_ifs, gcell_faces_idx, gcell_faces_lst);

  }

  BFT_FREE(gcell_faces_idx);
  BFT_FREE(gcell_faces_lst);

  *p_gcell_vtx_idx = gcell_vtx_idx;
  *p_gcell_vtx_lst = gcell_vtx_lst;

  /* Update mesh structure elements bound to halo management */

  if (mesh->n_ghost_cells > 0)
    BFT_REALLOC(mesh->cell_family, mesh->n_cells_with_ghosts, cs_lnum_t);

  cs_halo_update_buffers(halo);

#if 0 /* for debugging purposes */
  cs_halo_dump(halo, 1);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
