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
 * Functions dealing with ghost cells
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_periodicity.h>
#include <fvm_interface.h>
#include <fvm_order.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

typedef struct _vtx_lookup_table {

  cs_int_t   n_vertices;    /* Number of local vertices in the problem domain */
  cs_int_t   n_transforms;  /* Number of transformations */
  cs_int_t   n_interfaces;  /* Number of interfaces */
  cs_int_t   n_categories;  /* Number of possible categories
                               = n_interfaces * (n_transforms + 1)
                               (1 category for purely parallel elements) */

  cs_int_t   *if_ranks;     /* List of ranks */
  cs_int_t   *rank_ids;     /* list of rank ids */

  cs_int_t   *index;        /* index on table (size = n_vertices + 1) */

  cs_int_t   *rank_list;    /* list of ranks on which vertices are linked */
  cs_int_t   *type_list;    /* list of type (purelly parallel (=0) or number
                               of the periodicity) featuring a vertex. This
                               list is only allocated when n_perio > 0 */

} vtx_lookup_table_t;


typedef struct _table_int {

  cs_int_t  n_elts;
  cs_int_t  size_max;
  cs_int_t  *values;

} table_int_t;


typedef struct {

  cs_int_t  n_elts;   /* Number of elements of the index */
  cs_int_t  *index;   /* size = n_elts + 1: [ 0 ... n_elts ] */

  cs_int_t  size;     /* size of the list */
  cs_int_t  n_lists;  /* number of lists */
  cs_int_t  **lists;  /* list of the corresponding elements */

} lookup_table_t;

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
 * Create a table_int_t structure and initialize it with a given value and
 * a given size.
 *
 * parameters:
 *   init_value  -->  initial values
 *   size        -->  size of the table
 *
 * returns:
 *   pointer to a table_int_t structure
 *---------------------------------------------------------------------------*/

static table_int_t *
_create_table_int(cs_int_t  init_value,
                  cs_int_t  size)
{
  cs_int_t  i;

  table_int_t  *ret_table = NULL;

  BFT_MALLOC(ret_table, 1, table_int_t);
  BFT_MALLOC(ret_table->values, size, cs_int_t);

  ret_table->size_max = size;
  ret_table->n_elts = 0;

  for (i = 0; i < size; i++)
    ret_table->values[i] = init_value;

  return ret_table;
}

/*---------------------------------------------------------------------------
 * Delete table_int_t structure.
 *
 * parameters:
 *   this_table  -->  a pointer to a table_int_t structure
 *
 * returns:
 *   A NULL pointer
 *---------------------------------------------------------------------------*/

static table_int_t *
_delete_table_int(table_int_t  *this_table)
{
  BFT_FREE(this_table->values);
  this_table->size_max = 0;
  this_table->n_elts = 0;
  BFT_FREE(this_table);

  return NULL;
}

/*---------------------------------------------------------------------------
 * Reset a table_int_t structure. Memory is not freed but structure
 * is reinitialized.
 *
 * parameters:
 *   this_table -->  a pointer to a table_int_t structure
 *   init_value -->  new value given to the elements of the table
 *---------------------------------------------------------------------------*/

static void
_reset_table_int(table_int_t  *this_table,
                 cs_int_t      init_value)
{
  cs_int_t  i;

  for (i = 0; i < this_table->n_elts; i++)
    this_table->values[i] = init_value;
  this_table->n_elts = 0;

}

/*---------------------------------------------------------------------------
 * Find value in table_int_t structure.
 *
 * parameters:
 *   this_table   -->  a pointer to a table_int_t structure
 *   value        -->  value to find
 *
 * returns:
 *  The associated index or -1 if value was not found
 *---------------------------------------------------------------------------*/

static cs_int_t
_find_table_int_value(table_int_t  *this_table,
                      cs_int_t     value)
{
  cs_int_t  i;

  cs_int_t  ret_index = -1;

  for (i = 0; i < this_table->n_elts; i++) {
    if (this_table->values[i] == value)
      break;
  }

  if (i < this_table->n_elts)
    ret_index = i;

  return ret_index;
}

/*---------------------------------------------------------------------------
 * Add value in a table_int_t structure if this value does not already exist
 * in the table.
 *
 * parameters:
 *   this_table   -->  a pointer to a table_int_t structure
 *   value        -->  value to add
 *---------------------------------------------------------------------------*/

static void
_add_table_int_value(table_int_t  *this_table,
                     cs_int_t     value)
{
  cs_int_t  index;

  index = _find_table_int_value(this_table, value);

  if (index == -1) { /* Value does not already exist in the table */

    if (this_table->n_elts == this_table->size_max) {
      this_table->size_max = 2*this_table->n_elts;
      BFT_REALLOC(this_table->values, this_table->size_max, cs_int_t);
    }

    this_table->values[this_table->n_elts] = value;
    this_table->n_elts += 1;

  } /* End of value addition */

}

/*---------------------------------------------------------------------------
 * Add a value in a table_int_t structure. Value can already exist in the
 * table.
 *
 * parameters:
 *   this_table   -->  a pointer to a table_int_t structure
 *   value        -->  value to find
 *---------------------------------------------------------------------------*/

static void
_add_table_int_value_dup(table_int_t  *this_table,
                         cs_int_t     value)
{

  if (this_table->n_elts == this_table->size_max) {
    this_table->size_max = 2*this_table->n_elts;
    BFT_REALLOC(this_table->values, this_table->size_max, cs_int_t);
  }

  this_table->values[this_table->n_elts] = value;
  this_table->n_elts += 1;

}

/*---------------------------------------------------------------------------
 * Raise the value of index "idx" by "count".
 *
 * parameters:
 *   this_table   -->  a pointer to a table_int_t structure
 *   idx          -->  value has the index "idx"
 *   count        -->  value is raised by count
 *---------------------------------------------------------------------------*/

static void
_raise_table_int_value(table_int_t  *this_table,
                       cs_int_t      idx,
                       cs_int_t      count)
{

  if (idx >= this_table->n_elts)
    bft_error(__FILE__, __LINE__, 0,
              _("_raise_table_int_value:\n"
                "Invalid index. table_int structure has a lower n_elts.\n"
                "Index = %d, n_elts = %d and table size max = %d\n"),
              idx, this_table->n_elts, this_table->size_max);

  else
    this_table->values[idx] += count;

}

/*---------------------------------------------------------------------------
 * Get value for element of index "idx" in table_int_t.
 *
 * parameters:
 *
 * returns:
 *   The value of the element of index "idx"
 *---------------------------------------------------------------------------*/

static cs_int_t
_get_table_int_value(table_int_t  *this_table,
                     cs_int_t      idx)
{
  cs_int_t  retval = -9999;

  if (idx >= this_table->n_elts)
    bft_error(__FILE__, __LINE__, 0,
              _("get_table_int_value:\n"
                "Invalid index. table_int structure has a lower n_elts.\n"
                "Index = %d, n_elts = %d and table size max = %d\n"),
              idx, this_table->n_elts, this_table->size_max);

  else
    retval = this_table->values[idx];

  return retval;
}

/*---------------------------------------------------------------------------
 * Delete a lookup_table_t structure
 *
 * this_table  -->  pointer to a lookup_table_t structure
 *
 * returns:
 *   NULL
 *---------------------------------------------------------------------------*/

static lookup_table_t *
_delete_lookup_table(lookup_table_t  *this_table)
{
  cs_int_t  i;

  BFT_FREE(this_table->index);

  for (i = 0; i < this_table->n_lists; i++)
    BFT_FREE(this_table->lists[i]);
  BFT_FREE(this_table->lists);

  BFT_FREE(this_table);

  return NULL;
}

/*---------------------------------------------------------------------------
 * Fill a look up table structure without periodicity
 *
 * The look up list takes the distant rank value with which a local vertex
 * is linked.
 *
 * parameters:
 *   vtx_lookup  -->  pointer to a vtx_lookup_table_t structure
 *   ifs         -->  pointer to a fvm_interface_set_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_vtx_lookup(vtx_lookup_table_t   *vtx_lookup,
                 fvm_interface_set_t  *ifs)
{
  cs_int_t  i, vtx_id, rank_id, shift;

  cs_int_t  *counter = NULL;

  const cs_int_t  n_interfaces = fvm_interface_set_size(ifs);

  BFT_MALLOC(counter, vtx_lookup->n_vertices, cs_int_t);

  for (i = 0; i < vtx_lookup->n_vertices; i++)
    counter[i] = 0;

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    const fvm_interface_t  *interface = fvm_interface_set_get(ifs, rank_id);
    const fvm_lnum_t  interface_size = fvm_interface_size(interface);
    const fvm_lnum_t  *local_num = fvm_interface_get_local_num(interface);

    for (i = 0; i < interface_size; i++) { /* Only parallel vertices */

      vtx_id = local_num[i]-1;
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
 *   vtx_look_up  -->  pointer to a vtx_lookup_table_t structure
 *   ifs          -->  pointer to a fvm_interface_set_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_vtx_lookup_with_perio(vtx_lookup_table_t   *vtx_lookup,
                            fvm_interface_set_t  *ifs)
{
  cs_int_t  i, tr_id, vtx_id, rank_id, shift;

  cs_int_t  *counter = NULL;

  const fvm_periodicity_t  *periodicity = fvm_interface_set_periodicity(ifs);
  const cs_int_t  n_interfaces = fvm_interface_set_size(ifs);
  const cs_int_t  n_transforms = fvm_periodicity_get_n_transforms(periodicity);

  assert(n_transforms > 0);

  BFT_MALLOC(counter, vtx_lookup->n_vertices, cs_int_t);

  for (i = 0; i < vtx_lookup->n_vertices; i++)
    counter[i] = 0;

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    const fvm_interface_t  *interface
      = fvm_interface_set_get(ifs, vtx_lookup->rank_ids[rank_id]);
    const fvm_lnum_t  tr_index_size = fvm_interface_get_tr_index_size(interface);
    const fvm_lnum_t  *tr_index = fvm_interface_get_tr_index(interface);
    const fvm_lnum_t  *local_num = fvm_interface_get_local_num(interface);

    assert(n_transforms + 2 == tr_index_size);
    assert(tr_index != NULL);

    for (i = tr_index[0]; i < tr_index[1]; i++) { /* Only parallel vertices */

      vtx_id = local_num[i]-1;
      shift = vtx_lookup->index[vtx_id] + counter[vtx_id];

      vtx_lookup->rank_list[shift] = rank_id;
      vtx_lookup->type_list[shift] = 0;
      counter[vtx_id] += 1;

    }

    for (tr_id = 0; tr_id < n_transforms; tr_id++) {

      for (i = tr_index[tr_id + 1]; i < tr_index[tr_id + 2]; i++) {

        vtx_id = local_num[i]-1;
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
 *   n_vertices  -->  number of vertices of the table.
 *   ifs         -->  pointer to a fvm_interface_set_t structure
 *
 * returns:
 *   A pointer to the created vtx_lookup_table_t structure
 *---------------------------------------------------------------------------*/

static vtx_lookup_table_t *
_vtx_lookup_create(cs_int_t              n_vertices,
                   fvm_interface_set_t  *ifs)
{
  cs_int_t  i, rank_id, tmp_id, interface_size;

  cs_int_t  loc_rank_id = -1;
  vtx_lookup_table_t  *vtx_lookup = NULL;

  const fvm_interface_t  *interface = NULL;
  const fvm_lnum_t  *local_num = NULL;
  const fvm_periodicity_t  *periodicity = fvm_interface_set_periodicity(ifs);
  const cs_int_t  n_transforms = fvm_periodicity_get_n_transforms(periodicity);
  const cs_int_t  n_interfaces = fvm_interface_set_size(ifs);

  BFT_MALLOC(vtx_lookup, 1, vtx_lookup_table_t);

  vtx_lookup->n_vertices = n_vertices;
  vtx_lookup->n_interfaces = n_interfaces;
  vtx_lookup->n_transforms = n_transforms;
  vtx_lookup->n_categories = (n_transforms + 1)*n_interfaces;

  BFT_MALLOC(vtx_lookup->index, n_vertices + 1, cs_int_t);
  BFT_MALLOC(vtx_lookup->if_ranks, n_interfaces, cs_int_t);
  BFT_MALLOC(vtx_lookup->rank_ids, n_interfaces, cs_int_t);

  for (i = 0; i < n_vertices + 1; i++)
    vtx_lookup->index[i] = 0;

  /* Check if cs_glob_base_rang belongs to the interface set in order to
     arrange if_ranks with local rank at first place */

  for (rank_id = 0; rank_id < n_interfaces; rank_id++) {

    interface = fvm_interface_set_get(ifs, rank_id);
    vtx_lookup->if_ranks[rank_id] = fvm_interface_rank(interface);
    vtx_lookup->rank_ids[rank_id] = rank_id;

    if (cs_glob_base_rang == fvm_interface_rank(interface))
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

      fvm_lnum_t  *order = NULL;
      fvm_gnum_t  *buffer = NULL;
      cs_int_t  *_rank_ids = NULL;

      assert(sizeof(fvm_lnum_t) == sizeof(cs_int_t));

      BFT_MALLOC(order, n_interfaces - 1, fvm_lnum_t);
      BFT_MALLOC(buffer, n_interfaces - 1, fvm_gnum_t);
      BFT_MALLOC(_rank_ids, n_interfaces , cs_int_t);

      _rank_ids[0] = vtx_lookup->rank_ids[0];
      for (i = 1; i < n_interfaces; i++) {
        buffer[i-1] = (fvm_gnum_t)vtx_lookup->if_ranks[i];
        _rank_ids[i] = vtx_lookup->rank_ids[i];
      }

      fvm_order_local_allocated(NULL,
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

    interface = fvm_interface_set_get(ifs, rank_id);
    interface_size = fvm_interface_size(interface);
    local_num = fvm_interface_get_local_num(interface);

    for (i = 0; i < interface_size; i++)
      vtx_lookup->index[(cs_int_t)local_num[i]] += 1;

  } /* End of loop on if_ranks */

  /* Create index and allocate buffers */

  for (i = 0; i < n_vertices; i++)
    vtx_lookup->index[i+1] += vtx_lookup->index[i];

  BFT_MALLOC(vtx_lookup->rank_list, vtx_lookup->index[n_vertices], cs_int_t);

  /* Second loop to fill table(s) */

  if (n_transforms > 0) {

    BFT_MALLOC(vtx_lookup->type_list, vtx_lookup->index[n_vertices], cs_int_t);
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

  if (vtx_lookup->type_list != NULL)
    BFT_FREE(vtx_lookup->type_list);

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
  cs_int_t  i, rank_id, type;

  const cs_int_t  n_interfaces = vtx_lookup->n_interfaces;
  const cs_int_t  n_transforms = vtx_lookup->n_transforms;

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
_update_face_checker(cs_int_t          n_face_vertices,
                     cs_int_t          n_categories,
                     cs_int_t         *vtx_checker,
                     cs_halo_type_t   *face_checker)
{
  cs_int_t  i;

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
 *   mesh         --> pointer to a mesh structure
 *   face_checker --> halo type in each categories for cell's faces
 *---------------------------------------------------------------------------*/

static void
_count_halo_elements(cs_mesh_t        *mesh,
                     cs_halo_type_t   *face_checker)
{
  cs_int_t  type_id, rank_id;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_halo_t  *halo = mesh->halo;
  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  stride = 4*n_c_domains;

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
  cs_int_t  i, rank_id;
  cs_int_t  buffer_size, n_halo_elts, n_per_halo_elts;

  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  n_init_perio = mesh->n_init_perio;
  const cs_int_t  stride = 4*n_c_domains;
  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;

  buffer_size = 0;
  for (i = 0; i < 2*n_c_domains; i++)
    buffer_size += halo->send_index[i+1];

  BFT_MALLOC(halo->send_list, buffer_size, cs_int_t);

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
 *   mesh          -->  pointer to a mesh structure.
 *   face_checker  -->  halo type of each face of the cell.
 *   cell_id       -->  numbering of the treated cell.
 *   counter       <->  counter on each categories.
 *---------------------------------------------------------------------------*/

static void
_add_halo_elements(cs_mesh_t        *mesh,
                   cs_halo_type_t   *face_checker,
                   cs_int_t          cell_id,
                   cs_int_t         *counter)
{
  cs_int_t  i, type_id, shift, c_shift;

  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  n_c_domains = halo->n_c_domains;

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
 *
 *---------------------------------------------------------------------------*/

static cs_bool_t
_test_loop_continues(cs_mesh_t   *mesh,
                     cs_int_t     face_num)
{
  cs_int_t  fac_id = face_num - mesh->n_b_faces - 1;
  cs_bool_t  choice = false;

  /* Face has to be an internal face */

  if (face_num > mesh->n_b_faces) {

    if (mesh->halo_type == CS_HALO_STANDARD) {

      if (   mesh->i_face_cells[2*fac_id] < 1
          || mesh->i_face_cells[2*fac_id+1] < 1 )
        choice = true;
      else
        choice = false;

    }
    else {

      assert(mesh->halo_type == CS_HALO_EXTENDED);
      choice = true;

    }

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
 *   mesh           --> pointer to cs_mesh_t structure
 *   interface_set  --> pointer to fvm_interface_set_t structure
 *   cell_faces_idx --> "cell -> faces" connectivity index
 *   cell_faces_lst --> "cell -> faces" connectivity list
 *---------------------------------------------------------------------------*/

static void
_fill_send_halo(cs_mesh_t            *mesh,
                fvm_interface_set_t  *interface_set,
                cs_int_t             *cell_faces_idx,
                cs_int_t             *cell_faces_lst)
{
  cs_int_t  i, cell_id, i_fac, i_vtx;
  cs_int_t  fac_id, vtx_id, fac_num;
  cs_int_t  n_face_vertices;

  cs_halo_type_t  type_tag = CS_HALO_N_TYPES;
  cs_halo_type_t  face_type = CS_HALO_N_TYPES;
  cs_halo_type_t  cell_type = CS_HALO_N_TYPES;
  cs_int_t  n_categories = 0;
  cs_halo_t  *halo = mesh->halo;
  vtx_lookup_table_t  *vtx_lookup = NULL;
  cs_halo_type_t  *cell_tag = NULL;
  cs_halo_type_t  *face_checker = NULL;
  cs_int_t  *vtx_checker = NULL;
  cs_int_t  *counter = NULL;

  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_int_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  /* We should have the faces -> vertices connectivity to continue */

  if (mesh->i_face_vtx_lst == NULL && mesh->b_face_vtx_lst == NULL)
    return;

  /* Create a lookup table to accelerate search in
     fvm_interface_set structure */

  vtx_lookup = _vtx_lookup_create(n_vertices, interface_set);

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

  if (cell_faces_idx != NULL && cell_faces_lst != NULL) {

    BFT_MALLOC(vtx_checker, n_categories, cs_int_t);
    BFT_MALLOC(face_checker, n_categories, cs_halo_type_t);
    BFT_MALLOC(cell_tag, n_cells, cs_halo_type_t);

    /* First loop to create index and allocate cell_list */

    for (cell_id = 0; cell_id < n_cells; cell_id++) {

      /* Default initialization */

      cell_type = CS_HALO_N_TYPES;

      for (i = 0; i < n_categories; i++)
        face_checker[i] = CS_HALO_N_TYPES;

      /* Loop on faces of the cell */

      for (i_fac = cell_faces_idx[cell_id];
           i_fac < cell_faces_idx[cell_id + 1]; i_fac++) {

        fac_num = CS_ABS(cell_faces_lst[i_fac - 1]);

        if (_test_loop_continues(mesh, fac_num) == true) {

          fac_id = fac_num - mesh->n_b_faces - 1;
          n_face_vertices = fac_vtx_idx[fac_id + 1] - fac_vtx_idx[fac_id];

          /* Initialize checker */

          for (i = 0; i < n_categories; i++)
            vtx_checker[i] = 0;

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id + 1]; i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx - 1] - 1;

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

    BFT_MALLOC(counter, 2*n_categories, cs_int_t);

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

        for (i_fac = cell_faces_idx[cell_id];
             i_fac < cell_faces_idx[cell_id+1]; i_fac++) {

          /* Initialize checker */

          for (i = 0; i < n_categories; i++)
            vtx_checker[i] = 0;

          fac_num = CS_ABS(cell_faces_lst[i_fac - 1]);

          if (_test_loop_continues(mesh, fac_num) ==  true) {

            fac_id = fac_num - mesh->n_b_faces - 1;

            /* Loop on vertices of the face */

            n_face_vertices = fac_vtx_idx[fac_id + 1] - fac_vtx_idx[fac_id];

            for (i_vtx = fac_vtx_idx[fac_id];
                 i_vtx < fac_vtx_idx[fac_id + 1]; i_vtx++) {

              vtx_id = fac_vtx_lst[i_vtx - 1] - 1;

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

  } /* End if cell_face_idx and cell_faces_lst != NULL */

  /* Destroy the lookup table strcuture */

  _vtx_lookup_destroy(vtx_lookup);

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
 * Define a buffer on vertices where vertex belonging to the interface_set
 * are tagged with 1 else 0.
 *
 * parameters:
 *   n_vertices    --> size of the buffer
 *   interface_set --> pointer to a fvm_interface_set_t structure
 *   p_vertex_tag  <-> pointer to the tagged buffer
 *---------------------------------------------------------------------------*/

static void
_get_vertex_tag(cs_int_t                    n_vertices,
                const fvm_interface_set_t  *interface_set,
                cs_int_t                   *p_vertex_tag[])
{
  cs_int_t  i, j, rank_id;

  cs_int_t  *vertex_tag = NULL;

  const int  ifs_size = fvm_interface_set_size(interface_set);

  BFT_MALLOC(vertex_tag, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    vertex_tag[i] = 0;

  for (rank_id = 0; rank_id < ifs_size; rank_id++) {

    const fvm_interface_t  *interface = fvm_interface_set_get(interface_set, rank_id);
    const fvm_lnum_t  *local_num = fvm_interface_get_local_num(interface);
    const fvm_lnum_t  if_size = fvm_interface_size(interface);

    for (j = 0; j < if_size; j++)
      vertex_tag[local_num[j]-1] = 1;

  } /* End of loop on ranks */

  *p_vertex_tag = vertex_tag;

}

/*---------------------------------------------------------------------------
 * Compute the number of purely parallel ghost cells for a specific rank.
 *
 * parameters:
 *   mesh      --> pointer to cs_mesh_t structure
 *   rank      --> rank on which we want to know the number of purely
 *                 parallel elements.
 *   type      --> standard or extended
 *   index     --> index on halo's elements
 *   perio_lst --> periodic details on halo
 *
 * returns:
 *  Number of purely parallel elements in the halo.
 *---------------------------------------------------------------------------*/

static cs_int_t
_get_n_par_ghost_cells(cs_mesh_t       *mesh,
                       cs_int_t         rank,
                       cs_halo_type_t   type,
                       cs_int_t        *index,
                       cs_int_t        *perio_lst)
{
  cs_int_t  i;
  cs_int_t  n_per_gcells = 0, n_par_gcells = 0;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  n_c_domains = mesh->halo->n_c_domains;

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
 *   mesh --> pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_fill_halo(cs_mesh_t  *mesh)
{
  cs_int_t  rank_id, i, j;
  cs_int_t  shift;

#if defined(_CS_HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  int request_count = 0;

  cs_int_t  *count = NULL;

  cs_halo_t  *halo = mesh->halo;

  const  cs_int_t  n_c_domains = halo->n_c_domains;
  const  cs_int_t  n_transforms = mesh->n_transforms;
  const  cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;

#if defined(_CS_HAVE_MPI)
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

#if defined(_CS_HAVE_MPI)
      MPI_Irecv(&(halo->index[2*rank_id+1]), 2, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                halo->c_domain_rank[rank_id],
                cs_glob_base_mpi_comm,
                &(request[request_count++]));
#endif

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Barrier(cs_glob_base_mpi_comm);
#endif

  BFT_MALLOC(count, 2*n_c_domains, cs_int_t);

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = 2*rank_id;
    count[shift] =   halo->send_index[2*rank_id+1]
                   - halo->send_index[2*rank_id];
    count[shift+1] =   halo->send_index[2*rank_id+2]
                     - halo->send_index[2*rank_id+1];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Isend(&(count[shift]), 2, CS_MPI_INT,
                halo->c_domain_rank[rank_id],
                local_rank,
                cs_glob_base_mpi_comm,
                &(request[request_count++]));
#endif

    }
    else {

      halo->index[shift+1] = count[shift];
      halo->index[shift+2] = count[shift+1];

    }

  } /* End of loop on ranks */

  /* Wait for all exchanges being done */

#if defined(_CS_HAVE_MPI)
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

    cs_int_t  n_elts;
    cs_int_t  *exchange_buffer = NULL;

    /* n_transforms periodicities to deal with and for each sub-periodicity
       2 data. One for standard halo and the other one for extended halo */

    const cs_int_t  n_elts_to_exchange = 2 * n_transforms;

    BFT_MALLOC(exchange_buffer, 4*n_transforms, cs_int_t);

    for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

      if (halo->c_domain_rank[rank_id] != local_rank) {

        /* Fill buffer to send */

        for (i = 0; i < n_transforms; i++) {
          shift = 4*n_c_domains*i + 4*rank_id;
          for (j = 0; j < 2; j++)
            exchange_buffer[2*i+j] = halo->send_perio_lst[shift + 2*j + 1];
        }

#if defined(_CS_HAVE_MPI)
        MPI_Sendrecv(&(exchange_buffer[0]), n_elts_to_exchange, CS_MPI_INT,
                     halo->c_domain_rank[rank_id], local_rank,
                     &(exchange_buffer[n_elts_to_exchange]), n_elts_to_exchange,
                     CS_MPI_INT,
                     halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                     cs_glob_base_mpi_comm, status);
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

#if defined(_CS_HAVE_MPI)
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

  halo->n_elts[CS_HALO_EXTENDED] += halo->n_elts[CS_HALO_STANDARD] ;

}

/*---------------------------------------------------------------------------
 * Compute maximum list buffer size.
 *
 * This is done to avoid a reallocation for each rank and transformation.
 *
 * parameters:
 *   ifs --> pointer to a fvm_interface_set_t structure
 *
 * returns:
 *  max buffer size
 *---------------------------------------------------------------------------*/

static cs_int_t
_get_list_buffer_size(fvm_interface_set_t  *ifs)
{
  cs_int_t  i, j, tr_index_size;

  cs_int_t  max_lst_size = 0;

  const fvm_interface_t  *interface = NULL;
  const fvm_lnum_t  *tr_index = NULL;
  const cs_int_t  ifs_size = fvm_interface_set_size(ifs);

  if (ifs == NULL)
    return max_lst_size;

  for (i = 0; i < ifs_size; i++) {

    interface = fvm_interface_set_get(ifs, i);
    tr_index = fvm_interface_get_tr_index(interface);
    tr_index_size = fvm_interface_get_tr_index_size(interface) - 1;

    if (tr_index != NULL)
      for (j = 0; j < tr_index_size; j++)
        max_lst_size = CS_MAX(max_lst_size, tr_index[j+1] - tr_index[j]);
    else
      max_lst_size = CS_MAX(max_lst_size,
                            (cs_int_t)fvm_interface_size(interface));

  } /* End of loop on interfaces */

  return max_lst_size;
}

/*---------------------------------------------------------------------------
 * Define an index on vertices belonging to this interface for this rank
 * and this transformation.
 *
 * parameters:
 *   ifs                --> pointer to fvm_interface_set_t structure
 *   rank_id            --> rank number to work with
 *   tr_id              --> transformation id to work with
 *   vtx_interface_idx  <-> index on vertices which match the criterions
 *---------------------------------------------------------------------------*/

static void
_define_vtx_interface_idx(fvm_interface_set_t  *ifs,
                          cs_int_t              rank_id,
                          cs_int_t              tr_id,
                          cs_int_t              n_vertices,
                          cs_int_t             *vtx_interface_idx)
{
  cs_int_t  i, j, id;

  /* Initialize index */

  for (i = 0; i < n_vertices + 1; i++)
    vtx_interface_idx[i] = 0;

  for (id = 0; id < fvm_interface_set_size(ifs); id++) {

    const fvm_interface_t  *interface = fvm_interface_set_get(ifs, id);

    if (rank_id == fvm_interface_rank(interface)) {

      const fvm_lnum_t  *tr_index = fvm_interface_get_tr_index(interface);
      const fvm_lnum_t  *local_num = fvm_interface_get_local_num(interface);

      if (tr_index == NULL) {  /*  purelly parallel elements */

        for (j = 0; j < (cs_int_t)fvm_interface_size(interface); j++)
          vtx_interface_idx[local_num[j]] += 1;

      }
      else {

        for (j = tr_index[tr_id]; j < tr_index[tr_id+1]; j++)
          vtx_interface_idx[local_num[j]] += 1;

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
 * Fill the dist_num_lst which is the list of distant numbering associated
 * to local vertices (same rank and same transformation).
 *
 * parameters:
 *   ifs               --> pointer to fvm_interface_set_t structure
 *   rank_id           --> rank number to work with
 *   tr_id             --> transformation id to work with
 *   n_vertices        --> number of vertices
 *   dist_num_lst      <-> list of distant vertex numbers matching criteria
 *   vtx_interface_idx <-> index on vertices matching criteria
 *---------------------------------------------------------------------------*/

static void
_define_dist_num_lst(fvm_interface_set_t  *ifs,
                     cs_int_t              rank_id,
                     cs_int_t              tr_id,
                     cs_int_t              n_vertices,
                     cs_int_t             *dist_num_lst,
                     cs_int_t             *vtx_interface_idx)
{
  cs_int_t  i, j, id, shift;

  /* Initialize index */

  for (i = 0; i < n_vertices + 1; i++)
    vtx_interface_idx[i] = 0;

  for (id = 0; id < fvm_interface_set_size(ifs); id++) {

    const fvm_interface_t  *interface = fvm_interface_set_get(ifs, id);

    if (rank_id == fvm_interface_rank(interface)) {

      const fvm_lnum_t  *tr_index = fvm_interface_get_tr_index(interface);
      const fvm_lnum_t  *distant_num = fvm_interface_get_distant_num(interface);
      const fvm_lnum_t  *local_num = fvm_interface_get_local_num(interface);

      if (tr_index == NULL)
        for (j = 0; j < (cs_int_t)fvm_interface_size(interface); j++)
          vtx_interface_idx[local_num[j]] += 1;

      else
        for (j = tr_index[tr_id]; j < tr_index[tr_id+1]; j++)
          vtx_interface_idx[local_num[j]] += 1;

      /* Create index */

      for (j = 0; j < n_vertices; j++)
        vtx_interface_idx[j+1] += vtx_interface_idx[j];

      /* There must be only one distant_num per local_num when
         we treat a specific rank and a specific periodicity.
         So, we don't need a counter to fill dist_num_lst */

      if (tr_index == NULL) {

        for (j = 0; j < (cs_int_t)fvm_interface_size(interface); j++) {
          shift = vtx_interface_idx[local_num[j]-1];
          dist_num_lst[shift] = distant_num[j];
        }


      }
      else {

        for (j = tr_index[tr_id]; j < tr_index[tr_id+1]; j++) {
          shift = vtx_interface_idx[local_num[j]-1];
          dist_num_lst[shift] = distant_num[j];
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
 *   mesh        --> pointer to cs_mesh_t structure
 *   index       --> index on halo's elements
 *   perio_lst   --> periodic details on halo
 *   rank_id     --> rank number to work with
 *   tr_id       --> transformation id to work with
 *   type_id     --> standard or extended
 *   p_start_idx <-- pointer on start index
 *   p_end_idx   <-- pointer on end index
 *---------------------------------------------------------------------------*/

static void
_get_start_end_idx(cs_mesh_t    *mesh,
                   cs_int_t     *index,
                   cs_int_t     *perio_lst,
                   cs_int_t      rank_id,
                   cs_int_t      tr_id,
                   cs_int_t      type_id,
                   cs_int_t     *p_start_idx,
                   cs_int_t     *p_end_idx)
{
  cs_int_t  i, n_par_gcells, n_per_gcells;
  cs_int_t  start_idx, end_idx;

  const cs_int_t  n_c_domains = mesh->halo->n_c_domains;

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
 *   mesh               --> pointer to cs_mesh_t structure
 *   ifs                --> pointer to fvm_interface_set_t structure
 *   rank_id            --> rank number to work with
 *   tr_id              --> transformation id to work with
 *   cell_faces_idx     --> "cell -> faces" connectivity index
 *   cell_faces_lst     --> "cell -> faces" connectivity list
 *   vtx_interface_idx  <-> index on vertices which match the criterions
 *   vtx_tag            <-> tag array on vertices
 *   gcell_dist_vtx_idx <-> "ghost cell -> distant vertices" connectivity index
 *---------------------------------------------------------------------------*/

static void
_count_send_gcell_to_dist_vtx_connect(cs_mesh_t            *mesh,
                                      fvm_interface_set_t  *ifs,
                                      cs_int_t              rank_id,
                                      cs_int_t              tr_id,
                                      cs_int_t             *cell_faces_idx,
                                      cs_int_t             *cell_faces_lst,
                                      cs_int_t             *vtx_interface_idx,
                                      cs_int_t             *vtx_tag,
                                      cs_int_t             *gcell_dist_vtx_idx)
{
  cs_int_t  id, cell_id, i_fac, i_vtx, i_loop;
  cs_int_t  start_idx, end_idx, vtx_id, fac_num, fac_id;

  cs_int_t  n_loops = 0;
  cs_int_t  n_added_vertices = 0;

  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_int_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  const char err_corresp[]
    = N_("Repeated inconsistency in the halo construction.\n"
         "Several local points have the same distant correspondant;\n"
         "this is probably due to a multiple-periodicity construction\n"
         "side effect of the Preprocessor.\n"
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
           i_fac < cell_faces_idx[cell_id+1]; i_fac++) {

        fac_num = CS_ABS(cell_faces_lst[i_fac-1]);

        if (_test_loop_continues(mesh, fac_num) == true) {

          fac_id = fac_num - mesh->n_b_faces - 1;

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id+1]; i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx-1] - 1;

            /* If vertex is on the interface for this rank and this
             transformation */

            n_added_vertices =  vtx_interface_idx[vtx_id+1]
                              - vtx_interface_idx[vtx_id];

            if (n_added_vertices > 0) {

              if (n_added_vertices > 1)
                bft_error(__FILE__, __LINE__, 0, _(err_corresp),
                          mesh->vtx_coord[vtx_id*3],
                          mesh->vtx_coord[vtx_id*3+1],
                          mesh->vtx_coord[vtx_id*3+2]);

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
 *   mesh               --> pointer to cs_mesh_t structure
 *   ifs                --> pointer to fvm_interface_set_t structure
 *   rank_id            --> rank number to work with
 *   tr_id              --> transformation id to work with
 *   cell_faces_idx     --> "cell -> faces" connectivity index
 *   cell_faces_lst     --> "cell -> faces" connectivity list
 *   vtx_interface_idx  <-> index on vertices matching criteria
 *   dist_num_lst       <-> list of distant vertex numbers matching criteria
 *   counter            <-> counter on vertices
 *   vtx_tag            <-> tag array on vertices
 *   gcell_dist_vtx_idx <-> "ghost cell -> distant vertices" connectivity index
 *   gcell_dist_vtx_lst <-> "ghost cell -> distant vertices" connectivity list
 *---------------------------------------------------------------------------*/

static void
_fill_send_gcell_to_dist_vtx_connect(cs_mesh_t            *mesh,
                                     fvm_interface_set_t  *ifs,
                                     cs_int_t              rank_id,
                                     cs_int_t              tr_id,
                                     cs_int_t             *cell_faces_idx,
                                     cs_int_t             *cell_faces_lst,
                                     cs_int_t             *vtx_interface_idx,
                                     cs_int_t             *dist_num_lst,
                                     cs_int_t             *counter,
                                     cs_int_t             *vtx_tag,
                                     cs_int_t             *gcell_dist_vtx_idx,
                                     cs_int_t             *gcell_dist_vtx_lst)
{
  cs_int_t  i, id, cell_id, i_fac, i_vtx, i_loop;
  cs_int_t  shift, vtx_id, fac_num, fac_id, start_idx, end_idx;

  cs_int_t  n_loops = 0;
  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_int_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  _define_dist_num_lst(ifs,
                       halo->c_domain_rank[rank_id],
                       tr_id,
                       n_vertices,
                       dist_num_lst,
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
           i_fac < cell_faces_idx[cell_id+1]; i_fac++) {

        fac_num = CS_ABS(cell_faces_lst[i_fac-1]);

        if (_test_loop_continues(mesh, fac_num) == true) {

          fac_id = fac_num - mesh->n_b_faces - 1;

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id+1]; i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx-1] - 1;

            /* If vertex is on the interface for this rank and periodicity */

            if (vtx_interface_idx[vtx_id+1] - vtx_interface_idx[vtx_id] > 0) {

              /* Add this vertex if nont already checked */

              if (vtx_tag[vtx_id] != id) { /* Add this vertex */

                vtx_tag[vtx_id] = id;

                for (i = vtx_interface_idx[vtx_id];
                     i < vtx_interface_idx[vtx_id+1]; i++) {

                  shift = gcell_dist_vtx_idx[id] + counter[id];
                  gcell_dist_vtx_lst[shift] = dist_num_lst[i];
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
 *   mesh                 --> pointer to cs_mesh_t structure
 *   interface_set        --> pointer to fvm_interface_set_t structure
 *   cell_faces_idx       --> "cell -> faces" connectivity index
 *   cell_faces_lst       --> "cell -> faces" connectivity list
 *   p_gcell_dist_vtx_idx <-- "ghost cell -> distant vertices" connect. index
 *   p_gcell_dist_vtx_lst <-- "ghost cell -> distant vertices" connect. list
 *---------------------------------------------------------------------------*/

static void
_create_send_gcell_vtx_connect(cs_mesh_t            *mesh,
                               fvm_interface_set_t  *interface_set,
                               cs_int_t             *cell_faces_idx,
                               cs_int_t             *cell_faces_lst,
                               cs_int_t             *p_gcell_dist_vtx_idx[],
                               cs_int_t             *p_gcell_dist_vtx_lst[])
{
  cs_int_t  i, id, rank_id;

  cs_int_t  *gcell_dist_vtx_idx = NULL, *gcell_dist_vtx_lst = NULL;
  cs_int_t  *vtx_interface_idx = NULL;
  cs_int_t  *dist_num_lst = NULL;
  cs_int_t  *vtx_tag = NULL;
  cs_int_t  *counter = NULL;

  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  max_lst_size = _get_list_buffer_size(interface_set);
  const cs_int_t  n_ghost_cells = halo->n_send_elts[CS_HALO_EXTENDED];
  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  tr_index_size = mesh->n_transforms + 1;

  if (n_ghost_cells == 0)
    return;

  /* Allocate and initialize buffers */

  BFT_MALLOC(gcell_dist_vtx_idx, n_ghost_cells + 1, cs_int_t);
  BFT_MALLOC(counter, n_ghost_cells, cs_int_t);

  gcell_dist_vtx_idx[0] = 0;
  for (i = 0; i < n_ghost_cells; i++) {
    gcell_dist_vtx_idx[i+1] = 0;
    counter[i] = 0;
  }

  BFT_MALLOC(vtx_tag, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    vtx_tag[i] = -1;

  BFT_MALLOC(vtx_interface_idx, n_vertices + 1, cs_int_t);
  BFT_MALLOC(dist_num_lst, max_lst_size, cs_int_t);

  for (i = 0; i < max_lst_size; i++)
    dist_num_lst[i] = -1;

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

  BFT_MALLOC(gcell_dist_vtx_lst, gcell_dist_vtx_idx[n_ghost_cells], cs_int_t);

  for (i = 0; i < n_vertices; i++)
    vtx_tag[i] = -1;

  /* Fill gcell_dist_vtx_lst */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    for (id = 0; id < tr_index_size; id++)
      _fill_send_gcell_to_dist_vtx_connect(mesh,
                                           interface_set,
                                           rank_id,
                                           id,
                                           cell_faces_idx,
                                           cell_faces_lst,
                                           vtx_interface_idx,
                                           dist_num_lst,
                                           counter,
                                           vtx_tag,
                                           gcell_dist_vtx_idx,
                                           gcell_dist_vtx_lst);

  } /* End of loop on ranks */

  BFT_FREE(counter);
  BFT_FREE(vtx_tag);
  BFT_FREE(vtx_interface_idx);
  BFT_FREE(dist_num_lst);

  *p_gcell_dist_vtx_idx = gcell_dist_vtx_idx;
  *p_gcell_dist_vtx_lst = gcell_dist_vtx_lst;

}

/*---------------------------------------------------------------------------
 * Send "ghost cells to distant_num vertices" connectivity on communicating
 * ranks and receive the same kind of connectivity from distant ranks.
 *
 * parameters:
 *   mesh                     -->  pointer to cs_mesh_t structure
 *   send_gcell_dist_vtx_idx  <--  "ghost cell -> distant vertices" index
 *   send_gcell_dist_vtx_lst  <--  "ghost cell -> distant vertices" list
 *   p_gcell_dist_vtx_idx     -->  "ghost cell -> distant vertices" index
 *   p_gcell_dist_vtx_lst     -->  "ghost cell -> distant vertices" list
 *---------------------------------------------------------------------------*/

static void
_exchange_gcell_vtx_connect(cs_mesh_t  *mesh,
                            cs_int_t   *send_gcell_dist_vtx_idx,
                            cs_int_t   *send_gcell_dist_vtx_lst,
                            cs_int_t   *p_gcell_dist_vtx_idx[],
                            cs_int_t   *p_gcell_dist_vtx_lst[])
{
  cs_int_t  i, j, rank_id;
  cs_int_t  send_start_idx, send_end_idx, start_idx, end_idx;
  cs_int_t  n_send_elts, n_recv_elts;

  cs_int_t  send_buffer_size = 0;

  cs_int_t  *send_idx_buffer = NULL;
  cs_int_t  *gcell_dist_vtx_idx = NULL, *gcell_dist_vtx_lst = NULL;
  cs_int_t  *send_buffer = NULL, *recv_buffer = NULL;

  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;
  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  n_ghost_cells = halo->n_elts[CS_HALO_EXTENDED];

#if defined(_CS_HAVE_MPI)
  MPI_Status  status;
#endif

  /* Allocate buffers */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {
    if (halo->c_domain_rank[rank_id] != local_rank) {
      n_send_elts = halo->send_index[2*rank_id+2]- halo->send_index[2*rank_id];
      send_buffer_size = CS_MAX(send_buffer_size, n_send_elts);
    }
  }

  BFT_MALLOC(send_idx_buffer, send_buffer_size, cs_int_t);
  BFT_MALLOC(gcell_dist_vtx_idx, n_ghost_cells + 1, cs_int_t);

  for (i = 0; i < n_ghost_cells + 1; i++)
    gcell_dist_vtx_idx[i] = 0;

  /* Exchange sizes to define gcell_dist_vtx_idx */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    recv_buffer = &(gcell_dist_vtx_idx[1 + halo->index[2*rank_id]]);

    if (halo->c_domain_rank[rank_id] != local_rank) {

      /* Fill send buffer */

      for (i = halo->send_index[2*rank_id], j = 0;
           i < halo->send_index[2*rank_id+2]; i++, j++)
        send_idx_buffer[j] =  send_gcell_dist_vtx_idx[i+1] - send_gcell_dist_vtx_idx[i];

      n_send_elts =  halo->send_index[2*rank_id+2] - halo->send_index[2*rank_id];
      n_recv_elts =  halo->index[2*rank_id+2] - halo->index[2*rank_id];

#if defined (_CS_HAVE_MPI)
      MPI_Sendrecv(&(send_idx_buffer[0]), n_send_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], local_rank,
                   &(recv_buffer[0]), n_recv_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                   cs_glob_base_mpi_comm, &status);
#endif

    } /* If rank != local_rank */

    else {

      for (i = halo->send_index[2*rank_id], j = 0;
           i < halo->send_index[2*rank_id+2]; i++, j++)
        recv_buffer[j] = send_gcell_dist_vtx_idx[i+1] - send_gcell_dist_vtx_idx[i];

    } /* rank == local_rank */

  } /* End of loop on if_ranks */

  BFT_FREE(send_idx_buffer);

  /* Define index */

  for (i = 0; i < n_ghost_cells; i++)
    gcell_dist_vtx_idx[i+1] += gcell_dist_vtx_idx[i];

  BFT_MALLOC(gcell_dist_vtx_lst, gcell_dist_vtx_idx[n_ghost_cells], cs_int_t);

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

#if defined (_CS_HAVE_MPI)
      MPI_Sendrecv(send_buffer, n_send_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], local_rank,
                   recv_buffer, n_recv_elts, CS_MPI_INT,
                   halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                   cs_glob_base_mpi_comm, &status);
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
 * Reverse "ghost cell -> vertex" connectivity into "vertex -> ghost cells"
 * connectivity for halo elements.
 *
 * Build the connectivity index.
 *
 * parameters:
 *   halo            -->  pointer to a cs_halo_t structure
 *   n_vertices      -->  number of vertices
 *   rank_id         -->  rank number to work with
 *   checker         <->  temporary array to check vertices
 *   gcell_vtx_idx   -->  "ghost cell -> vertices" connectivity index
 *   gcell_vtx_lst   -->  "ghost cell -> vertices" connectivity list
 *   vtx_gcells_idx  <->  "vertex -> ghost cells" connectivity index
 *---------------------------------------------------------------------------*/

static void
_reverse_connectivity_idx(cs_halo_t   *halo,
                          cs_int_t     n_vertices,
                          cs_int_t     rank_id,
                          cs_int_t    *checker,
                          cs_int_t    *gcell_vtx_idx,
                          cs_int_t    *gcell_vtx_lst,
                          cs_int_t    *vtx_gcells_idx)
{
  cs_int_t  i, j, id, vtx_id, start_idx, end_idx;

  /* Initialize index */

  vtx_gcells_idx[0] = 0;
  for (i = 0; i < n_vertices; i++) {
    vtx_gcells_idx[i+1] = 0;
    checker[i] = -1;
  }

  if (rank_id == -1) {
    start_idx = 0;
    end_idx = halo->n_elts[CS_HALO_EXTENDED];
  }
  else { /* Call with rank_id > 1 for standard halo */
    start_idx = halo->index[2*rank_id];
    end_idx = halo->index[2*rank_id+1];
  }

  /* Define index */

  for (id = start_idx; id < end_idx; id++) {

    for (j = gcell_vtx_idx[id]; j < gcell_vtx_idx[id+1]; j++) {

      vtx_id = gcell_vtx_lst[j] - 1;

      if (checker[vtx_id] != id) {
        checker[vtx_id] = id;
        vtx_gcells_idx[vtx_id+1] += 1;
      }

    }

  } /* End of loop of ghost cells */

  for (i = 0; i < n_vertices; i++)
    vtx_gcells_idx[i+1] += vtx_gcells_idx[i];

}

/*---------------------------------------------------------------------------
 * Reverse "ghost cells -> vertex" connectivity into "vertex -> ghost cells"
 * connectivity for halo elements.
 *
 * Build the connectivity list.
 *
 * parameters:
 *   halo            -->  pointer to a cs_halo_t structure
 *   n_vertices      -->  number of vertices
 *   rank_id         -->  rank number to work with
 *   counter         <->  temporary array to count vertices
 *   checker         <->  temporary array to check vertices
 *   gcell_vtx_idx   -->  "ghost cell -> vertices" connectivity index
 *   gcell_vtx_lst   -->  "ghost cell -> vertices" connectivity list
 *   vtx_gcells_idx  -->  "vertex -> ghost cells" connectivity index
 *   vtx_gcells_lst  <->  "vertex -> ghost cells" connectivity list
 *---------------------------------------------------------------------------*/

static void
_reverse_connectivity_lst(cs_halo_t   *halo,
                          cs_int_t     n_vertices,
                          cs_int_t     rank_id,
                          cs_int_t    *counter,
                          cs_int_t    *checker,
                          cs_int_t    *gcell_vtx_idx,
                          cs_int_t    *gcell_vtx_lst,
                          cs_int_t    *vtx_gcells_idx,
                          cs_int_t    *vtx_gcells_lst)
{
  cs_int_t  i, j, id, shift, vtx_id, start_idx, end_idx;

  /* Initialize buffers */

  for (i = 0; i < n_vertices; i++) {
    counter[i] = 0;
    checker[i] = -1;
  }

  if (rank_id == -1) {
    start_idx = 0;
    end_idx = halo->n_elts[CS_HALO_EXTENDED];
  }
  else {
    start_idx = halo->index[2*rank_id];
    end_idx = halo->index[2*rank_id+1];
  }

  /* Fill the connectivity list */

  for (id = start_idx; id < end_idx; id++) {

    for (j = gcell_vtx_idx[id]; j < gcell_vtx_idx[id+1]; j++) {

      vtx_id = gcell_vtx_lst[j] - 1;

      if (checker[vtx_id] != id) {

        checker[vtx_id] = id;
        shift = vtx_gcells_idx[vtx_id] + counter[vtx_id];
        vtx_gcells_lst[shift] = id;
        counter[vtx_id] += 1;

      }

    }

  } /* End of loop of ghost cells */

}

/*---------------------------------------------------------------------------
 * Define a lookup table structure on periodic faces.
 *
 * parameters:
 *   mesh                     --> pointer to cs_mesh_t structure
 *   mesh_builder             --> pointer to cs_mesh_builder_t structure
 *
 * returns:
 *   a pointer to a lookup_table_t structure.
 *---------------------------------------------------------------------------*/

static lookup_table_t *
_define_periodic_lookup_table(cs_mesh_t           *mesh,
                              cs_mesh_builder_t   *mesh_builder)
{
  cs_int_t  i, fst_face, snd_face, shift;

  lookup_table_t  *periodic_table =  NULL;
  cs_int_t  *face_list = NULL, *rank_list = NULL;

  const cs_int_t  n_init_perio = mesh->n_init_perio;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  *per_face_idx = mesh_builder->per_face_idx;
  const cs_int_t  *per_face_lst = mesh_builder->per_face_lst;
  const cs_int_t  *per_rank_lst = mesh_builder->per_rank_lst;

  assert(mesh->n_init_perio > 0);
  assert(per_face_idx != NULL);

  /* Create a lookup table structure */

  BFT_MALLOC(periodic_table, 1, lookup_table_t);

  /* Define the lookup table structure */

  periodic_table->n_elts = n_i_faces;

  BFT_MALLOC(periodic_table->index, n_i_faces + 1, cs_int_t);

  for (i = 0; i < n_i_faces + 1; i++)
    periodic_table->index[i] = 0;

  for (i = per_face_idx[0]; i < per_face_idx[n_init_perio]; i++)
    periodic_table->index[CS_ABS(per_face_lst[2*i])] += 1;

  for (i = 0; i < n_i_faces; i++)
    periodic_table->index[i+1] += periodic_table->index[i];

  periodic_table->size = periodic_table->index[n_i_faces];

  if (per_rank_lst != NULL)
    periodic_table->n_lists = 2;
  else
    periodic_table->n_lists = 1;

  BFT_MALLOC(periodic_table->lists, periodic_table->n_lists, cs_int_t *);

  for (i = 0; i < periodic_table->n_lists; i++)
    BFT_MALLOC(periodic_table->lists[i], periodic_table->size, cs_int_t);

  face_list = periodic_table->lists[0];

  if (per_rank_lst != NULL)
    rank_list = periodic_table->lists[1];

  /* Define face_list and rank_list */

  for (i = per_face_idx[0]; i < per_face_idx[n_init_perio]; i++) {

    fst_face = CS_ABS(per_face_lst[2*i]) - 1;
    snd_face = per_face_lst[2*i+1] - 1;

    /* A face is periodic with only one face */

    assert(  periodic_table->index[fst_face + 1]
           - periodic_table->index[fst_face] == 1);

    shift = periodic_table->index[fst_face];
    face_list[shift] = snd_face;

    if (per_rank_lst != NULL)
      rank_list[shift] = per_rank_lst[i];

  }

  return periodic_table;
}

/*---------------------------------------------------------------------------
 * Define the number of vertices owned by each face of halo cells.
 * This will enable to check the right ghost cell in case of multiple
 * choices.
 *
 * parameters:
 *   mesh                   <->  pointer to cs_mesh_t structure
 *   cell_faces_idx         -->  "cell -> faces" connectivity index
 *   cell_faces_lst         -->  "cell -> faces" connectivity list
 *   perio_table            -->  pointer to the periodic face lookup table
 *   gcell_faces_idx        -->  pointer to "ghost cell -> faces" conn. index
 *   gcell_glob_faces_lst   -->  pointer to "ghost cell -> glob face" list
 *---------------------------------------------------------------------------*/

static void
_define_gcells_connect(cs_mesh_t       *mesh,
                       cs_int_t        *cell_faces_idx,
                       cs_int_t        *cell_faces_lst,
                       cs_int_t        *p_gcell_faces_idx[],
                       cs_int_t        *p_gcell_glob_faces_lst[])
{
  cs_int_t  i, j, cell_id, fac_num, fac_id;
  cs_int_t  rank_id, n_faces, shift, _shift, counter;
  cs_int_t  distant_rank, value;

  cs_int_t  n_elts = 0;
  int  request_count = 0;
  cs_int_t  *rank_shift = NULL, *send_buffer = NULL;
  cs_int_t  *gcell_faces_idx = NULL;
  cs_int_t  *gcell_glob_faces_lst = NULL;
  lookup_table_t  *perio_table = NULL;

  cs_halo_t  *halo = mesh->halo;
  cs_mesh_builder_t  *mesh_builder = cs_glob_mesh_builder;

#if defined(_CS_HAVE_MPI)
  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;
#endif

  const cs_int_t  n_init_perio = mesh->n_init_perio;
  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  local_rank = (cs_glob_base_rang == -1) ? 0:cs_glob_base_rang;
  const cs_int_t  n_gcells = halo->n_elts[CS_HALO_EXTENDED];
  const cs_int_t  n_send_gcells = halo->n_send_elts[CS_HALO_EXTENDED];

#if defined(_CS_HAVE_MPI)
  if (halo->n_c_domains*2 > 128) {
    BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
    BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);
  }
#endif

  BFT_MALLOC(gcell_faces_idx, n_gcells + 1, cs_int_t);
  BFT_MALLOC(send_buffer, n_send_gcells, cs_int_t);

  /* Exchange number of face for each ghost cells */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = halo->index[2*rank_id] + 1;
    n_elts = halo->index[2*rank_id+2] - halo->index[2*rank_id];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Irecv(&(gcell_faces_idx[shift]), n_elts, CS_MPI_INT,
                halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                cs_glob_base_mpi_comm, &(request[request_count++]));
#endif

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Barrier(cs_glob_base_mpi_comm);
#endif

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = halo->send_index[2*rank_id];
    n_elts = halo->send_index[2*rank_id+2] - halo->send_index[2*rank_id];

    for (i = halo->send_index[2*rank_id];
         i < halo->send_index[2*rank_id+2];
         i++) {

      cell_id = halo->send_list[i];
      n_faces = cell_faces_idx[cell_id+1] - cell_faces_idx[cell_id];
      send_buffer[i] = n_faces;

    }

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Isend(&(send_buffer[shift]), n_elts, CS_MPI_INT,
                halo->c_domain_rank[rank_id], local_rank,
                cs_glob_base_mpi_comm, &(request[request_count++]));
#endif

    }
    else {

      for (i = halo->send_index[2*rank_id];
           i < halo->send_index[2*rank_id+2];
           i++)
        gcell_faces_idx[i+1] = send_buffer[i];

    }

  } /* End of loop on ranks */

  /* Wait for all exchanges to finish */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Waitall(request_count, request, status);
#endif
  request_count = 0;

  /* Build index */

  gcell_faces_idx[0] = 0;
  for (i = 0; i < n_gcells; i++)
    gcell_faces_idx[i+1] += gcell_faces_idx[i];

  /* Exchange global face numbering for each face of ghost cells */

  BFT_MALLOC(gcell_glob_faces_lst, gcell_faces_idx[n_gcells], cs_int_t);

  /* Define send_buffer size */

  BFT_MALLOC(rank_shift, n_c_domains + 1, cs_int_t);

  for (i = 0; i < n_c_domains + 1; i++)
    rank_shift[i] = 0;

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    for (i = halo->send_index[2*rank_id];
         i < halo->send_index[2*rank_id+2];
         i++) {

      cell_id = halo->send_list[i];

      for (j = cell_faces_idx[cell_id]; j < cell_faces_idx[cell_id+1]; j++) {

        fac_num = CS_ABS(cell_faces_lst[j-1]);

        if (_test_loop_continues(mesh, fac_num) == true)
          rank_shift[rank_id+1] += 1;

      } /* End of loop on faces */

    } /* End of loop on cells belonging to "in halo" */

  } /* End of loop on ranks */

  for (i = 0; i < n_c_domains; i++)
    rank_shift[i+1] += rank_shift[i];

  n_elts = rank_shift[n_c_domains];

  BFT_REALLOC(send_buffer, n_elts, cs_int_t);

  /* Define a look-up table on periodic faces.
     TODO: a global approach will have to be done to avoid mistakes
     on face -> cells connectivity in case of periodicity . This case
     appears only when a periodic face is on a parallel frontier and when
     several faces (at least two) shared the same n_face_vertices vertices.
     This extremely rarely occurs. */

  if (n_init_perio > 0)
    perio_table = _define_periodic_lookup_table(mesh, mesh_builder);

  /* Receive data */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = gcell_faces_idx[halo->index[2*rank_id]];
    n_elts =  gcell_faces_idx[halo->index[2*rank_id+2]]
            - gcell_faces_idx[halo->index[2*rank_id]];

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Irecv(&(gcell_glob_faces_lst[shift]), n_elts, CS_MPI_INT,
                halo->c_domain_rank[rank_id], halo->c_domain_rank[rank_id],
                cs_glob_base_mpi_comm, &(request[request_count++]));
#endif

    }

  } /* End of loop on ranks */

  /* We wait for receiving all messages */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Barrier(cs_glob_base_mpi_comm);
#endif

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    shift = rank_shift[rank_id];
    n_elts = rank_shift[rank_id+1] - rank_shift[rank_id];
    counter = 0;

    /* Fill send_buffer */

    for (i = halo->send_index[2*rank_id];
         i < halo->send_index[2*rank_id+2];
         i++) {

      cell_id = halo->send_list[i];

      for (j = cell_faces_idx[cell_id]; j < cell_faces_idx[cell_id+1]; j++) {

        fac_num = CS_ABS(cell_faces_lst[j-1]);

        if (_test_loop_continues(mesh, fac_num) == true) {

          fac_id = fac_num - mesh->n_b_faces - 1;

          if (mesh->global_i_face_num != NULL) { /* Parallel treatment */

            if (n_init_perio > 0) { /* Periodic treatment */

              if (  perio_table->index[fac_id+1]
                  - perio_table->index[fac_id] > 0) { /* periodic face */

                distant_rank =
                  perio_table->lists[1][perio_table->index[fac_id]] - 1;

                if (distant_rank == halo->c_domain_rank[rank_id])
                  value = - perio_table->lists[0][perio_table->index[fac_id]];
                else
                  value = 0;

              }
              else /* Face is not periodic */
                value = (cs_int_t)mesh->global_i_face_num[fac_id];

            }
            else /* Non periodic treatment */
              value = (cs_int_t)mesh->global_i_face_num[fac_id];

          }
          else { /* Non parallel treatment */

            if (n_init_perio > 0) { /* Periodic treatment */

              if (  perio_table->index[fac_id+1]
                  - perio_table->index[fac_id] > 0) { /* periodic face */

                value = perio_table->lists[0][perio_table->index[fac_id]] + 1;

              }
              else /* face is not periodic */
                value = fac_id + 1;

            }
            else /* Non periodic treatment */
              value = fac_id + 1;

          }

          send_buffer[shift + counter++] = value;

        } /* End if face has to be treat */

      } /* End of loop on cell->faces connectivity */

    } /* End of loop on cells belonging to "in halo" */

    if (halo->c_domain_rank[rank_id] != local_rank) {

#if defined(_CS_HAVE_MPI)
      MPI_Isend(&(send_buffer[shift]), n_elts, CS_MPI_INT,
                halo->c_domain_rank[rank_id], local_rank,
                cs_glob_base_mpi_comm, &(request[request_count++]));
#endif

    }
    else {

      assert(counter == n_elts);
      _shift = gcell_faces_idx[halo->index[2*rank_id]];

      for (i = 0; i < n_elts; i++)
        gcell_glob_faces_lst[_shift + i] = send_buffer[shift + i];

    }

  } /* End of loop on ranks */

  /* Wait for all exchanges being done */

#if defined(_CS_HAVE_MPI)
  if (mesh->n_domains > 1)
    MPI_Waitall(request_count, request, status);
#endif

  *p_gcell_faces_idx = gcell_faces_idx;
  *p_gcell_glob_faces_lst = gcell_glob_faces_lst;

  /* Free memory */

  BFT_FREE(rank_shift);
  BFT_FREE(send_buffer);

#if defined(_CS_HAVE_MPI)
  if (request != _request) {
    BFT_FREE(request);
    BFT_FREE(status);
  }
#endif

  if (n_init_perio > 0)
    perio_table = _delete_lookup_table(perio_table);

}

/*---------------------------------------------------------------------------
 * Update mesh->i_face_cells array for ghost cells in standard halo.
 *
 * Define face to ghost cells connectivity for standard halo and
 * if necessary, for extended halo.
 *
 * parameters:
 *   mesh               <->  pointer to cs_mesh_t structure
 *   cell_faces_idx     -->  "cell -> faces" connectivity index
 *   cell_faces_lst     -->  "cell -> faces" connectivity list
 *   gcell_vtx_idx      -->  "ghost cell -> vertices" connect. index
 *   gcell_vtx_lst      -->  "ghost cell -> vertices" connect. list
 *---------------------------------------------------------------------------*/

static void
_update_gcells_connect(cs_mesh_t       *mesh,
                       cs_int_t        *cell_faces_idx,
                       cs_int_t        *cell_faces_lst,
                       cs_int_t        *gcell_vtx_idx,
                       cs_int_t        *gcell_vtx_lst)
{
  cs_int_t  i, j, k, i_cel, i_fac, i_vtx;
  cs_int_t  rank_id, vtx_id, fac_num, fac_id, cell_id;
  cs_int_t  n_face_vertices, n_shared_vertices, tag, counter;
  cs_int_t  ghost_cell_num, _glob_face_num, glob_face_num;
  cs_bool_t  ok;

  table_int_t  *cell_checker = NULL, *cell_list = NULL;
  cs_int_t  *vtx_gcells_idx = NULL, *vtx_gcells_lst = NULL;
  cs_int_t  *gcell_faces_idx = NULL, *gcell_glob_faces_lst = NULL;
  cs_int_t  *vtx_buffer = NULL, *vtx_counter = NULL, *vtx_checker = NULL;

  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_vertices = mesh->n_vertices;
  const cs_int_t  n_c_domains = halo->n_c_domains;
  const cs_int_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_int_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  /* Define cell -> face  connectivity for halo to know
     exactly the right cell to use in face -> cells connectivity. */

  _define_gcells_connect(mesh,
                         cell_faces_idx,
                         cell_faces_lst,
                         &gcell_faces_idx,
                         &gcell_glob_faces_lst);

  /* Allocation and initialization of buffers */

  BFT_MALLOC(vtx_buffer, 3*n_vertices + 1, cs_int_t);
  vtx_counter = &(vtx_buffer[0]);
  vtx_checker = &(vtx_buffer[n_vertices]);
  vtx_gcells_idx = &(vtx_buffer[2*n_vertices]);

  cell_checker = _create_table_int(0, 20);
  cell_list = _create_table_int(-1, 20);

  /* Loop on each rank to define ghost cell connectivity */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++) {

    _reverse_connectivity_idx(halo,
                              n_vertices,
                              rank_id,
                              vtx_checker,
                              gcell_vtx_idx,
                              gcell_vtx_lst,
                              vtx_gcells_idx);

    BFT_REALLOC(vtx_gcells_lst, vtx_gcells_idx[n_vertices], cs_int_t);

    _reverse_connectivity_lst(halo,
                              n_vertices,
                              rank_id,
                              vtx_counter,
                              vtx_checker,
                              gcell_vtx_idx,
                              gcell_vtx_lst,
                              vtx_gcells_idx,
                              vtx_gcells_lst);

    for (i = halo->send_index[2*rank_id];
         i < halo->send_index[2*rank_id+2]; i++) {

      i_cel = halo->send_list[i];

      _reset_table_int(cell_list, -1);
      _reset_table_int(cell_checker, 0);

      for (i_fac = cell_faces_idx[i_cel];
           i_fac < cell_faces_idx[i_cel+1]; i_fac++) {

        fac_num = CS_ABS(cell_faces_lst[i_fac-1]);

        if (_test_loop_continues(mesh, fac_num) == true) {

          fac_id = fac_num - mesh->n_b_faces - 1;
          n_face_vertices = fac_vtx_idx[fac_id+1] - fac_vtx_idx[fac_id];

          /* Loop on vertices of the face */

          for (i_vtx = fac_vtx_idx[fac_id];
               i_vtx < fac_vtx_idx[fac_id+1]; i_vtx++) {

            vtx_id = fac_vtx_lst[i_vtx-1] - 1;

            /* For each vertex, we store cell in halo */

            for (j = vtx_gcells_idx[vtx_id];
                 j < vtx_gcells_idx[vtx_id+1]; j++) {

              tag = _find_table_int_value(cell_list, vtx_gcells_lst[j]);

              if (tag == -1) /* New cell */
                _add_table_int_value_dup(cell_checker, 1);
              else
                _raise_table_int_value(cell_checker, tag, 1);

              _add_table_int_value(cell_list, vtx_gcells_lst[j]);

            } /* End of loop on out ghost cells connecting to this vertex */

          } /* End of loop on face's vertices */

          counter = 0;
          for (j = 0; j < cell_checker->n_elts; j++) {
            n_shared_vertices = _get_table_int_value(cell_checker, j);
            if (n_shared_vertices == n_face_vertices)
              counter++;
          }

          if (counter > 0) {

            for (j = 0; j < cell_checker->n_elts; j++) {

              n_shared_vertices = _get_table_int_value(cell_checker, j);

              if (n_shared_vertices == n_face_vertices) { /* Standard Cell */

                cell_id = _get_table_int_value(cell_list, j);

                if (counter > 1) {

                  /* Test if this is the right ghost cell */

                  if (mesh->global_i_face_num != NULL)
                    glob_face_num = mesh->global_i_face_num[fac_id];
                  else
                    glob_face_num = fac_id + 1;

                  ok = false;
                  for (k = gcell_faces_idx[cell_id];
                       k < gcell_faces_idx[cell_id+1]; k++) {

                    _glob_face_num = gcell_glob_faces_lst[k];

                    if (_glob_face_num < 0) {

                      assert(mesh->global_i_face_num != NULL);
                      _glob_face_num = mesh->global_i_face_num[-_glob_face_num];

                    }

                    if (_glob_face_num == glob_face_num) {
                      ok = true;
                      break;
                    }

                  } /* End of loop on faces of the out ghost cell */

                } /* End if counter > 1 */

                else if (counter == 1)
                  ok = true;

                if (ok == true) { /* Update face -> cells connectivity */

                  ghost_cell_num = mesh->n_cells + cell_id + 1;

                  if (mesh->i_face_cells[2*fac_id] < 1)
                    mesh->i_face_cells[2*fac_id] = ghost_cell_num;

                  if (mesh->i_face_cells[2*fac_id+1] < 1)
                    mesh->i_face_cells[2*fac_id+1] = ghost_cell_num;

                }

              } /* If n_shared_vertices == n_face_vertices */

            } /* End of loop on cell sharing vertex/vertices */

          } /* End if counter > 0 */

          _reset_table_int(cell_checker, 0);
          _reset_table_int(cell_list, -1);

        } /* If face is concerned */

      } /* End of loop on faces */

    } /* End of loop on ghost cells */

  } /* End of loop on ranks */

  BFT_FREE(vtx_buffer);
  BFT_FREE(vtx_gcells_lst);
  BFT_FREE(gcell_faces_idx);
  BFT_FREE(gcell_glob_faces_lst);

  cell_checker = _delete_table_int(cell_checker);
  cell_list = _delete_table_int(cell_list);

}

#if 0 /* TODO: check algorithm (deadlock on BG/L on one test case) */

/*---------------------------------------------------------------------------
 * Clean a halo i.e. delete rank(s) with no ghost cells.
 *
 * This case happens when an interface has only extended vertices
 * and the halo is standard.
 *
 * parameters:
 *   mesh -->  pointer to a mesh structure
 *---------------------------------------------------------------------------*/

static void
_clean_halo(cs_mesh_t  *mesh)
{
  cs_int_t  rank_id, i, j;

  cs_halo_t  *halo = mesh->halo;

  cs_int_t  n_c_domains = halo->n_c_domains;
  cs_int_t  n_real_c_domains = 0;
  cs_int_t  counter = 0;
  cs_int_t  *new_c_domain_rank = NULL;
  cs_int_t  *new_perio_lst = NULL;
  cs_int_t  *new_index = NULL;

  const cs_int_t  n_transforms = mesh->n_transforms;

  /* Is there something to do ? */

  for (rank_id = 0; rank_id < n_c_domains; rank_id++)
    if (halo->send_index[2*rank_id+2] - halo->send_index[2*rank_id] > 0)
      n_real_c_domains++;

  if (n_real_c_domains == n_c_domains)
    return;


  /* halo->index, halo->perio_lst and n_c_domains need an update */

  BFT_MALLOC(new_c_domain_rank, n_real_c_domains, cs_int_t);
  BFT_MALLOC(new_index, 2*n_real_c_domains+1, cs_int_t);

  if (n_transforms > 0)
    BFT_MALLOC(new_perio_lst, 4*n_transforms*n_real_c_domains, cs_int_t);

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

  BFT_REALLOC(halo->index, 2*n_real_c_domains+1, cs_int_t);
  BFT_REALLOC(halo->perio_lst, 4*n_transforms*n_real_c_domains, cs_int_t);

}

#endif /* #if 0 */

/*----------------------------------------------------------------------------
 * Define cell -> internal faces connectivity for ghost cells.
 *
 * Treatment of parallel and/or periodic halos for standard or extended
 * ghost cells according to halo type building option.
 *
 * parameters:
 *   mesh             --> pointer to cs_mesh_t structure.
 *   interface_set    --> pointer to fvm_interface_set structure.
 *   p_cell_faces_idx <-- pointer to the connectivity index
 *   p_cell_faces_lst <-- pointer to the connectivity list
 *----------------------------------------------------------------------------*/

static void
_create_gcell_faces_connect(cs_mesh_t            *mesh,
                            fvm_interface_set_t  *interface_set,
                            cs_int_t             *p_cell_faces_idx[],
                            cs_int_t             *p_cell_faces_lst[])
{
  cs_int_t  i, fac_id, i_vtx, id1, id2, shift, vtx_id;

  cs_int_t  *vtx_tag = NULL;
  cs_int_t  *cell_buffer = NULL, *cell_tag = NULL, *counter = NULL;
  cs_int_t  *cell_faces_idx = NULL;
  cs_int_t  *cell_faces_lst = NULL;

  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  n_b_faces = mesh->n_b_faces;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  *face_cells = mesh->i_face_cells;
  const cs_int_t  *fac_vtx_idx = mesh->i_face_vtx_idx;
  const cs_int_t  *fac_vtx_lst = mesh->i_face_vtx_lst;

  *p_cell_faces_idx = cell_faces_idx;
  *p_cell_faces_lst = cell_faces_lst;

  if (interface_set == NULL)
    return;

  BFT_MALLOC(cell_faces_idx, n_cells+1, cs_int_t);
  BFT_MALLOC(cell_buffer, 2*n_cells, cs_int_t);

  cell_tag = &(cell_buffer[0]);
  counter = &(cell_buffer[n_cells]);

  cell_faces_idx[0] = 1;
  for (i = 0; i < n_cells; i++) {
    cell_faces_idx[i+1] = 0;
    cell_tag[i] = -1;
  }

  assert(sizeof(cs_int_t) == sizeof(fvm_lnum_t));

  _get_vertex_tag(mesh->n_vertices, interface_set, &vtx_tag);

  for (fac_id = 0; fac_id < n_i_faces; fac_id++) {

    for (i_vtx = fac_vtx_idx[fac_id];
         i_vtx < fac_vtx_idx[fac_id + 1]; i_vtx++) {

      vtx_id = fac_vtx_lst[i_vtx - 1] - 1;

      if (vtx_tag[vtx_id] == 1) {

        id1 = face_cells[2*fac_id] - 1;
        id2 = face_cells[2*fac_id + 1] - 1;

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

  BFT_MALLOC(cell_faces_lst, cell_faces_idx[n_cells]-1, cs_int_t);

  for (fac_id = 0; fac_id < n_i_faces; fac_id++) {

    for (i_vtx = fac_vtx_idx[fac_id];
         i_vtx < fac_vtx_idx[fac_id + 1]; i_vtx++) {

      vtx_id = fac_vtx_lst[i_vtx - 1] - 1;

      if (vtx_tag[vtx_id] == 1) {

        id1 = face_cells[2*fac_id] - 1;
        id2 = face_cells[2*fac_id + 1] - 1;

        if (id1 < 0) {
          if (cell_tag[id2] != fac_id) {

            cell_tag[id2] = fac_id;
            shift = cell_faces_idx[id2] - 1 + counter[id2];
            cell_faces_lst[shift] = fac_id + 1 + n_b_faces;
            counter[id2] += 1;

          }
        }

        if (id2 < 0) {
          if (cell_tag[id1] != fac_id) {

            cell_tag[id1] = fac_id;
            shift = cell_faces_idx[id1] - 1 + counter[id1];
            cell_faces_lst[shift] = fac_id + 1 + n_b_faces;
            counter[id1] += 1;

          }
        }

        if (mesh->halo_type == CS_HALO_EXTENDED) {
          if (id1 >= 0 && id2 >= 0) {

            if (cell_tag[id1] != fac_id) {
              cell_tag[id1] = fac_id;
              shift = cell_faces_idx[id1] - 1 + counter[id1];
              cell_faces_lst[shift] = fac_id + 1 + n_b_faces;
              counter[id1] += 1;
            }

            if (cell_tag[id2] != fac_id) {
              cell_tag[id2] = fac_id;
              shift = cell_faces_idx[id2] - 1 + counter[id2];
              cell_faces_lst[shift] = fac_id + 1 + n_b_faces;
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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define halo structures for internal and distant ghost cells.
 *
 * parameters:
 *   mesh             -->  pointer to cs_mesh_t structure
 *   interface_set    -->  pointer to fvm_interface_set_t structure.
 *   p_gcell_vtx_idx  <--  pointer to the connectivity index
 *   p_gcell_vtx_lst  <--  pointer to the connectivity list
 *---------------------------------------------------------------------------*/

void
cs_mesh_halo_define(cs_mesh_t            *mesh,
                    fvm_interface_set_t  *interface_set,
                    cs_int_t             *p_gcell_vtx_idx[],
                    cs_int_t             *p_gcell_vtx_lst[])
{
  cs_int_t  *send_gcell_vtx_idx = NULL, *send_gcell_vtx_lst = NULL;
  cs_int_t  *gcell_vtx_idx = NULL, *gcell_vtx_lst = NULL;
  cs_int_t  *gcell_faces_idx = NULL, *gcell_faces_lst = NULL;
  cs_halo_t  *halo = mesh->halo;

  halo->n_local_elts = mesh->n_cells;

  /*  Define cell -> internal faces connectivity for ghost cells */

  _create_gcell_faces_connect(mesh,
                              interface_set,
                              &gcell_faces_idx,
                              &gcell_faces_lst);

  /* Fill cs_halo_t structure for send_halo  */

  bft_printf(_("    Local halo definition\n"));
  bft_printf_flush();

  _fill_send_halo(mesh,
                  interface_set,
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

  bft_printf(_("    Distant halo creation\n"));
  bft_printf_flush();

  _fill_halo(mesh);

  /* Update mesh structure elements bound to halo management */

  mesh->n_ghost_cells = halo->n_elts[CS_HALO_EXTENDED];
  mesh->n_cells_with_ghosts = mesh->n_cells + mesh->n_ghost_cells;

  /* Update connectivity between internal faces and cells */

  if (cs_mesh_n_g_ghost_cells(mesh) > 0) {

    /* Create local ghost connectivity for send_halo cells and send it.
       Receive ghost cells connectivity for halo cells. */

    _create_send_gcell_vtx_connect(mesh,
                                   interface_set,
                                   gcell_faces_idx,
                                   gcell_faces_lst,
                                   &send_gcell_vtx_idx,
                                   &send_gcell_vtx_lst);

    _exchange_gcell_vtx_connect(mesh,
                                send_gcell_vtx_idx,
                                send_gcell_vtx_lst,
                                &gcell_vtx_idx,
                                &gcell_vtx_lst);

    /* Define mesh->i_face_cells array for ghost cells in standard halo and
       also ghost cells to ghost cells connectivity for standard and extended
       halo if necessary */

    bft_printf(_("    Updating the faces -> cells connectivity\n"));
    bft_printf_flush();

    _update_gcells_connect(mesh,
                           gcell_faces_idx,
                           gcell_faces_lst,
                           gcell_vtx_idx,
                           gcell_vtx_lst);

    /* Free memory */

    BFT_FREE(send_gcell_vtx_idx);
    BFT_FREE(send_gcell_vtx_lst);

  }

  BFT_FREE(gcell_faces_idx);
  BFT_FREE(gcell_faces_lst);

  *p_gcell_vtx_idx = gcell_vtx_idx;
  *p_gcell_vtx_lst = gcell_vtx_lst;

  /* Update mesh structure elements bound to halo management */

  if (mesh->n_ghost_cells > 0)
    BFT_REALLOC(mesh->cell_family, mesh->n_cells_with_ghosts, cs_int_t);

  cs_halo_update_buffers(halo);

#if 0 /* for debugging purposes */
  cs_halo_dump(halo, 1);
#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
