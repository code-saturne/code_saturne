/*============================================================================
 * Main structure for handling of interfaces associating mesh entities
 * (such as inter-processor or periodic connectivity between cells, faces,
 * or vertices);
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2006-2010  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_order.h"
#include "fvm_parall.h"
#include "fvm_periodicity.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_interface.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an interface
 *----------------------------------------------------------------------------*/

struct _fvm_interface_t {

  int          rank;           /* Associated rank */

  fvm_lnum_t   size;           /* Number of equivalent elements */

  int          tr_index_size;  /* Size of perio_index */
  fvm_lnum_t  *tr_index;       /* Index of sub-sections in local_num and
                                  distant_num for different transformations;
                                  purely parallel equivalences appear at
                                  position 0, equivalences through periodic
                                  transform i appear at position i+1;
                                  NULL in absence of transformations */

  fvm_lnum_t  *local_num;     /* Local element numbers */
  fvm_lnum_t  *distant_num;   /* Distant element numbers */

};

/*----------------------------------------------------------------------------
 * Structure defining a set of interfaces
 *----------------------------------------------------------------------------*/

struct _fvm_interface_set_t {

  int                        size;         /* Number of interfaces */

  fvm_interface_t          **interfaces;   /* Interface structures array */

  const fvm_periodicity_t   *periodicity;  /* Optional periodicity structure */

};

/*----------------------------------------------------------------------------
 * Local structure defining a temporary list of interfaces
 *----------------------------------------------------------------------------*/

typedef struct {

  int          count;    /* Number of equivalences */
  int         *shift;    /* Index of per-equivalence data in rank[] and num[] */
  int         *rank;     /* Rank associated with each element */
  int         *tr_id;    /* Transformation id associated with each element,
                            + 1, with 0 indicating no transformation
                            (NULL in absence of periodicity) */
  fvm_lnum_t  *num;      /* Local number associated with each element */

} _per_slice_equiv_t;

/*----------------------------------------------------------------------------
 * Local structure defining a temporary list of periodic interfaces
 *----------------------------------------------------------------------------*/

typedef struct {

  int          count;     /* Number of periodic couples */
  fvm_lnum_t  *slice_id;  /* local id in slice */
  int         *tr_id;     /* Transform id associated with each couple */
  int         *shift;     /* Index of per-couple data in rank[] and num[] */
  int         *rank;      /* Ranks associated with periodic elements */
  fvm_lnum_t  *num;       /* Local numbers associated with periodic elements */

} _per_slice_period_t;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of an empty interface between entities of a same type.
 *
 * This interface may be used to identify equivalent vertices or faces using
 * domain splitting, as well as periodic entities (on the same or on
 * distant ranks).
 *
 * returns:
 *  pointer to allocated interface structure
 *----------------------------------------------------------------------------*/

static fvm_interface_t *
_fvm_interface_create(void)
{
  fvm_interface_t  *_interface;

  BFT_MALLOC(_interface, 1, fvm_interface_t);

  _interface->rank = -1;
  _interface->size = 0;

  _interface->tr_index_size = 0;
  _interface->tr_index = NULL;

  _interface->local_num = NULL;
  _interface->distant_num = NULL;

  return _interface;
}

/*----------------------------------------------------------------------------
 * Destruction of an interface.
 *
 * parameters:
 *   this_interface <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static fvm_interface_t *
_fvm_interface_destroy(fvm_interface_t  *this_interface)
{
  if (this_interface != NULL) {
    BFT_FREE(this_interface->tr_index);
    BFT_FREE(this_interface->local_num);
    BFT_FREE(this_interface->distant_num);
    BFT_FREE(this_interface);
  }

  return this_interface;
}

/*----------------------------------------------------------------------------
 * Dump printout of an interface.
 *
 * parameters:
 *   this_interface <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

static void
_fvm_interface_dump(const fvm_interface_t  *this_interface)
{
  int i, section_id;
  fvm_lnum_t start_id, end_id;

  int  tr_index_size = 2;
  fvm_lnum_t  _tr_index[2] = {0, 0};
  const fvm_lnum_t  *tr_index = _tr_index;

  if (this_interface == NULL) {
    bft_printf("  interface: nil\n");
    return;
  }

  bft_printf("  interface:             %p\n"
             "  associated rank:       %d\n"
             "  size:                  %llu\n"
             "  transform index size:  %d\n",
             this_interface,
             this_interface->rank,
             (unsigned long long)(this_interface->size),
             this_interface->tr_index_size);

  if (this_interface->tr_index_size > 0) {
    bft_printf("  transform index:\n");
    for (i = 0; i < this_interface->tr_index_size; i++)
      bft_printf("    %5d %lu\n", i, this_interface->tr_index[i]);
  }


  _tr_index[1] = this_interface->size;

  if (this_interface->tr_index_size > 0) {
    tr_index_size = this_interface->tr_index_size;
    tr_index = this_interface->tr_index;
  }

  for (section_id = 0; section_id < tr_index_size - 1; section_id++) {

    if (section_id == 0)
      bft_printf("\n"
                 "            id      local    distant (parallel)\n");
    else
      bft_printf("\n"
                 "            id      local    distant (transform %d)\n",
                 section_id - 1);

    start_id = tr_index[section_id];
    end_id = tr_index[section_id + 1];

    if (this_interface->distant_num != NULL) {
      for (i = start_id; i < end_id; i++)
        bft_printf("    %10d %10d %10d\n", i,
                   this_interface->local_num[i],
                   this_interface->distant_num[i]);
    }
    else {
      for (i = start_id; i < end_id; i++)
        bft_printf("    %10d %10d\n", i,
                   this_interface->local_num[i]);
    }

  }

  bft_printf("\n");
}

/*----------------------------------------------------------------------------
 * Sort and remove duplicates from periodic couple information.
 *
 * parameters:
 *   n_slice_couples  <-> number of couples in current slice
 *   slice_couples    <-> couple information for this rank: for each couple,
 *                        {global number of local element,
 *                         global number of periodic element,
 *                         transform id}
 *----------------------------------------------------------------------------*/

static void
_sort_periodic_couples(fvm_lnum_t   *n_slice_couples,
                       fvm_gnum_t  **slice_couples)
{
  fvm_lnum_t  i, j, k;

  fvm_lnum_t   n_couples = *n_slice_couples;
  fvm_lnum_t  *order = NULL;
  fvm_gnum_t  *couples = *slice_couples;
  fvm_gnum_t  *couples_tmp = NULL;

  if (n_couples < 1)
    return;

  /* Sort periodic couples by local correspondant */

  BFT_MALLOC(order, n_couples, fvm_lnum_t);
  BFT_MALLOC(couples_tmp, n_couples*3, fvm_gnum_t);

  fvm_order_local_allocated_s(NULL,
                              couples,
                              3,
                              order,
                              n_couples);

  /* Copy to temporary array, ignoring duplicates */

  k = order[0]*3;
  couples_tmp[0] = couples[k];
  couples_tmp[1] = couples[k + 1];
  couples_tmp[2] = couples[k + 2];
  j = 3;

  for (i = 1; i < n_couples; i++) {
    k = order[i] * 3;
    if (   (couples[k]   != couples_tmp[j-3])
        || (couples[k+1] != couples_tmp[j-2])
        || (couples[k+2] != couples_tmp[j-1])) {
      couples_tmp[j]     = couples[k];
      couples_tmp[j + 1] = couples[k + 1];
      couples_tmp[j + 2] = couples[k + 2];
      j += 3;
    }
  }
  n_couples = j / 3;

  BFT_FREE(order);

  /* Resize input/outpout array if duplicates were removed */

  if (n_couples <= *n_slice_couples) {
    BFT_REALLOC(couples, n_couples*3, fvm_gnum_t);
    *n_slice_couples = n_couples;
    *slice_couples = couples;
  }

  /* Copy sorted data to input/output array and free temporary storage */

  memcpy(couples, couples_tmp, sizeof(fvm_gnum_t)*n_couples*3);

  BFT_FREE(couples_tmp);
}

/*----------------------------------------------------------------------------
 * Extract periodicity transform data necessary for periodic combinations.
 *
 * Builds a square transformation combination matrix associating a combination
 * transform id with transform ids of levels lower than that given by
 * the level argument; the number of rows and columns is thus equal to the
 * number of transformations of level lower than the argument.
 * Entries corresponding to impossible combinations are set to -1;
 *
 * The caller is responsible for freeing the tr_combine array returned by
 * this function.
 *
 * parameters:
 *   periodicity  <-- periodicity information
 *   level        --> transform combination level (1 or 2)
 *   n_rows       --> number of rows of combination matrix
 *                    (equals transform_level_index[level])
 *   tr_combine   --> combination id matrix (size n_rows * n_rows)
 *----------------------------------------------------------------------------*/

static void
_transform_combine_info(const fvm_periodicity_t   *periodicity,
                        int                        level,
                        int                       *n_rows,
                        int                      **tr_combine)
{
  int  i;
  int  tr_level_idx[4], parent_id[2];

  int n_vals_1 = 0, n_vals_2 = 0;
  int n_rows_1 = 0, n_rows_2 = 0;
  int *tr_combine_1, *tr_combine_2 = NULL;

  assert(periodicity != NULL);
  assert(level == 1 || level == 2);

  /* Extract base periodicity dimension info */

  fvm_periodicity_get_tr_level_idx(periodicity, tr_level_idx);

  /* We always need the level0 x level0 -> level1 array */

  n_rows_1 = tr_level_idx[1];
  n_vals_1 = n_rows_1*n_rows_1;

  BFT_MALLOC(tr_combine_1, n_vals_1, int);

  for (i = 0; i < n_vals_1; i++)
    tr_combine_1[i] = -1;

  for (i = tr_level_idx[1]; i < tr_level_idx[2]; i++) {

    fvm_periodicity_get_parent_ids(periodicity, i, parent_id);

    assert(parent_id[0] > -1 && parent_id[1] > -1);
    assert(parent_id[0] < n_rows_1 && parent_id[1] < n_rows_1);

    tr_combine_1[parent_id[0]*n_rows_1 + parent_id[1]] = i;
    tr_combine_1[parent_id[1]*n_rows_1 + parent_id[0]] = i;

  }

  /* For level 1 transforms, we are done once return values are set */

  if (level == 1) {

    *n_rows = n_rows_1;
    *tr_combine = tr_combine_1;

  }

  /* Handle level 2 transforms */

  else { /* if (level == 2) */

    int  tr_01, tr_02, tr_12;
    int  comp_id[3];

    n_rows_2 = tr_level_idx[2];
    n_vals_2 = n_rows_2*n_rows_2;

    /* Allocate and populate array */

    BFT_MALLOC(tr_combine_2, n_vals_2, int);

    for (i = 0; i < n_vals_2; i++)
      tr_combine_2[i] = -1;

    for (i = tr_level_idx[2]; i < tr_level_idx[3]; i++) {

      fvm_periodicity_get_components(periodicity, i, comp_id);

      assert(comp_id[0] > -1 && comp_id[1] > -1 && comp_id[2] > -1);

      tr_01 = tr_combine_1[comp_id[0]*n_rows_1 + comp_id[1]];
      tr_02 = tr_combine_1[comp_id[0]*n_rows_1 + comp_id[2]];
      tr_12 = tr_combine_1[comp_id[1]*n_rows_1 + comp_id[2]];

      assert(tr_01 > -1 && tr_02 > -1 && tr_12 > -1);

      tr_combine_2[tr_01*n_rows_2 + comp_id[2]] = i;
      tr_combine_2[comp_id[2]*n_rows_2 + tr_01] = i;

      tr_combine_2[tr_02*n_rows_2 + comp_id[1]] = i;
      tr_combine_2[comp_id[1]*n_rows_2 + tr_02] = i;

      tr_combine_2[tr_12*n_rows_2 + comp_id[0]] = i;
      tr_combine_2[comp_id[0]*n_rows_2 + tr_12] = i;

    }

    BFT_FREE(tr_combine_1);

    /* Set return values */

    *n_rows = n_rows_2;
    *tr_combine = tr_combine_2;
  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Maximum global number associated with an I/O numbering structure
 *
 * parameters:
 *   n_ent      <-- local number of entities
 *   global_num <-- global number (id) associated with each entity
 *   comm       <-- associated MPI communicator
 *
 * returns:
 *   maximum global number associated with the I/O numbering
 *----------------------------------------------------------------------------*/

static fvm_gnum_t
_global_num_max(fvm_lnum_t        n_ent,
                const fvm_gnum_t  global_num[],
                MPI_Comm          comm)
{
  fvm_gnum_t  local_max, global_max;

  /* Get maximum global number value */

  if (n_ent > 0)
    local_max = global_num[n_ent - 1];
  else
    local_max = 0;

  MPI_Allreduce(&local_max, &global_max, 1, FVM_MPI_GNUM, MPI_MAX, comm);

  return global_max;
}

/*----------------------------------------------------------------------------
 * Build temporary equivalence structure for data in a given slice,
 * and associate an equivalence id to received elements (-1 for elements
 * with no correponding elements)
 *
 * parameters:
 *   n_ranks         <-- number of associateed ranks
 *   n_ent_recv      <-- number of entities received
 *   recv_shift      <-- shift in received data per rank (size: n_ranks+1)
 *   recv_global_num <-- global numbering received
 *   recv_num        <-- local numbering received
 *   equiv_id        --> equivalence id for each element (-1 if none)
 *
 * returns:
 *   temporary equivalence structure for slice
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 2 /* Crash with O3 on IA64 with icc 9.1 20070320 */
#endif

static _per_slice_equiv_t
_slice_global_num_to_equiv(int                n_ranks,
                           fvm_lnum_t         n_ent_recv,
                           const fvm_lnum_t   recv_shift[],
                           const fvm_gnum_t   recv_global_num[],
                           const fvm_lnum_t   recv_num[],
                           fvm_lnum_t         equiv_id[])
{
  fvm_lnum_t  i, j;
  int         rank;
  fvm_gnum_t  cur_num, prev_num;

  int         *multiple = NULL;
  fvm_lnum_t  *recv_order = NULL;

  _per_slice_equiv_t  e;

  /* Initialize return structure */

  e.count    = 0;
  e.shift    = NULL;
  e.rank     = NULL;
  e.tr_id    = NULL;
  e.num      = NULL;

  if (n_ent_recv == 0)
    return e;

  /* Determine equivalent elements; requires ordering to loop through buffer
     by increasing number (slice blocks associated with each rank are
     already sorted, but the whole "gathered" slice is not). */

  BFT_MALLOC(recv_order, n_ent_recv, fvm_lnum_t);

  fvm_order_local_allocated(NULL,
                            recv_global_num,
                            recv_order,
                            n_ent_recv);

  /* Loop by increasing number: if two elements have the same global
     number, they are equivalent. We do not increment equivalence counts
     as soon as two elements are equivalent, as three equivalent elements
     should have the same equivalence id. Rather, we increment the
     equivalence counter when the previous element was part of an
     equivalence and the current element is not part of this same
     equivalence. */

  e.count = 0;

  equiv_id[recv_order[0]] = -1;
  prev_num = recv_global_num[recv_order[0]];

  for (i = 1; i < n_ent_recv; i++) {
    cur_num = recv_global_num[recv_order[i]];
    if (cur_num == prev_num) {
      equiv_id[recv_order[i-1]] = e.count;
      equiv_id[recv_order[i]]   = e.count;
    }
    else {
      if (equiv_id[recv_order[i-1]] > -1)
        e.count++;
      equiv_id[recv_order[i]] = -1;
    }
    prev_num = cur_num;
  }
  if (equiv_id[recv_order[n_ent_recv-1]] > -1)
    e.count++;

  BFT_FREE(recv_order);

  /* Count number of elements associated with each equivalence */

  BFT_MALLOC(multiple, e.count, int);

  BFT_MALLOC(e.shift, e.count+1, fvm_lnum_t);

  for (i = 0; i < e.count; multiple[i++] = 0);
  for (i = 0; i < n_ent_recv; i++) {
    if (equiv_id[i] > -1)
      multiple[equiv_id[i]] += 1;
  }

  e.shift[0] = 0;
  for (i = 0; i < e.count; i++)
    e.shift[i+1] = e.shift[i] + multiple[i];

  /* Build equivalence data */

  BFT_MALLOC(e.rank, e.shift[e.count], fvm_lnum_t);
  BFT_MALLOC(e.num, e.shift[e.count], fvm_lnum_t);

  for (i = 0; i < e.count; multiple[i++] = 0);

  for (rank = 0; rank < n_ranks; rank++) {
    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {
      if (equiv_id[i] > -1) {
        j = e.shift[equiv_id[i]] + multiple[equiv_id[i]];
        e.rank[j] = rank;
        e.num[j] = recv_num[i];
        multiple[equiv_id[i]] += 1;
      }
    }
  }

  BFT_FREE(multiple);

  return e;
}

/*----------------------------------------------------------------------------
 * Build global interface data from flat equivalence data
 * (usually prepared and received from distant ranks).
 *
 * parameters:
 *   this_interface_set  <-> pointer to structure that should be updated
 *   tr_index_size       <-- size of transform index (number of transforms
 *                           + 1 for identity + 1 for past-the-end);
 *                           0 or 1 if no transforms are present
 *   n_ent_recv          <-- size of received data
 *   equiv_recv          <-- flat (received) equivalence data; for each
 *                           equivalence, we have:
 *                           {local_number, n_equivalents,
 *                            {distant_number, distant_rank}*n_equivalents}
 *----------------------------------------------------------------------------*/

static void
_interfaces_from_flat_equiv(fvm_interface_set_t  *this_interface_set,
                            int                   tr_index_size,
                            fvm_lnum_t            n_ent_recv,
                            const fvm_lnum_t      equiv_recv[])
{
  fvm_lnum_t i, j, k, l;
  fvm_lnum_t local_num, distant_num, n_sub, n_ent_rank_tr_size;
  int rank, tr_id;

  int _tr_index_size = tr_index_size;
  int _tr_stride = tr_index_size > 1 ? tr_index_size - 1 : 1;
  int max_rank = 0, n_ranks = 0, start_id = 0;
  int recv_step = 1;

  fvm_lnum_t  *n_ent_rank = NULL;
  fvm_lnum_t  *n_ent_rank_tr = NULL;
  int *interface_id = NULL;

  fvm_interface_t *_interface = NULL;

  if (_tr_index_size == 0) {
    _tr_index_size = 1;
    _tr_stride = 1;
  }
  else if (_tr_index_size > 1)
    recv_step = 2;

  /* Compute size of subsections for each rank */

  i = 0;
  while (i < n_ent_recv) {
    i++;
    n_sub = equiv_recv[i++];
    for (j = 0; j < n_sub; j++) {
      i += recv_step;
      rank = equiv_recv[i++];
      if (rank > max_rank)
        max_rank = rank;
    }
  }

  BFT_MALLOC(n_ent_rank, max_rank + 1, fvm_lnum_t);

  for (i = 0; i < max_rank + 1; n_ent_rank[i++] = 0);

  i = 0;
  while (i < n_ent_recv) {
    i++;
    n_sub = equiv_recv[i++];
    for (j = 0; j < n_sub; j++) {
      i += recv_step;
      rank = equiv_recv[i++];
      n_ent_rank[rank] += 1;
    }
  }

  /* Build final data structures */

  n_ranks = 0;
  for (i = 0; i < max_rank + 1; i++) {
    if (n_ent_rank[i] > 0)
      n_ranks++;
  }

  /* (Re-)Allocate structures */

  start_id = this_interface_set->size;

  this_interface_set->size += n_ranks;

  BFT_REALLOC(this_interface_set->interfaces,
              this_interface_set->size,
              fvm_interface_t *);

  for (i = start_id; i < this_interface_set->size; i++)
    this_interface_set->interfaces[i] = _fvm_interface_create();

  /* Initialize rank info and interface id */

  n_ranks = 0;
  BFT_MALLOC(interface_id, max_rank + 1, int);
  for (i = 0; i < max_rank + 1; i++) {
    if (n_ent_rank[i] > 0) {
      interface_id[i] = start_id + n_ranks++;
      (this_interface_set->interfaces[interface_id[i]])->rank = i;
      (this_interface_set->interfaces[interface_id[i]])->size = n_ent_rank[i];
    }
    else
      interface_id[i] = -1;
  }

  BFT_FREE(n_ent_rank);

  /* n_ent_rank_tr will be used as a position counter for new interfaces */

  n_ent_rank_tr_size = (this_interface_set->size - start_id)*_tr_stride;
  BFT_MALLOC(n_ent_rank_tr, n_ent_rank_tr_size, fvm_lnum_t);
  for (i = 0; i < n_ent_rank_tr_size; i++)
    n_ent_rank_tr[i] = 0;

  for (i = start_id; i < this_interface_set->size; i++) {

    _interface = this_interface_set->interfaces[i];

    BFT_MALLOC(_interface->local_num, _interface->size, fvm_lnum_t);
    BFT_MALLOC(_interface->distant_num, _interface->size, fvm_lnum_t);

    if (_tr_index_size > 1) {
      _interface->tr_index_size = _tr_index_size;
      BFT_MALLOC(_interface->tr_index, _interface->tr_index_size, fvm_lnum_t);
      for (k = 0; k < _interface->tr_index_size; k++)
        _interface->tr_index[k] = 0;
    }
    else {
      _interface->tr_index_size = 0;
      _interface->tr_index = NULL;
    }

  }

  /* In absence of transforms, we may build the interface in one pass */

  if (_tr_index_size < 2)  {

    i = 0;
    while (i < n_ent_recv) {

      local_num = equiv_recv[i++];
      n_sub = equiv_recv[i++];

      for (j = 0; j < n_sub; j++) {

        distant_num = equiv_recv[i++];
        rank = equiv_recv[i++];

        _interface = this_interface_set->interfaces[interface_id[rank]];
        k = interface_id[rank] - start_id;

        _interface->local_num[n_ent_rank_tr[k]] = local_num;
        _interface->distant_num[n_ent_rank_tr[k]] = distant_num;
        n_ent_rank_tr[k] += 1;

      }
    }

  }

  /* If we have transforms, build the transform index first */

  else { /* if (_tr_index_size > 1) */

    /* Initial count */

    i = 0;
    while (i < n_ent_recv) {

      i++;
      n_sub = equiv_recv[i++];

      for (j = 0; j < n_sub; j++) {

        i++;
        tr_id = equiv_recv[i++];
        rank = equiv_recv[i++];

        _interface = this_interface_set->interfaces[interface_id[rank]];

        _interface->tr_index[tr_id + 1] += 1;

      }
    }

    /* Build index from initial count */

    for (i = start_id; i < this_interface_set->size; i++) {
      _interface = this_interface_set->interfaces[i];
      _interface->tr_index[0] = 0;
      for (j = 1; j < _tr_index_size; j++)
        _interface->tr_index[j] += _interface->tr_index[j-1];
    }

    /* Now populate the arrays */

    i = 0;
    while (i < n_ent_recv) {

      local_num = equiv_recv[i++];
      n_sub = equiv_recv[i++];

      for (j = 0; j < n_sub; j++) {

        distant_num = equiv_recv[i++];
        tr_id = equiv_recv[i++];
        rank = equiv_recv[i++];

        _interface = this_interface_set->interfaces[interface_id[rank]];
        k = (interface_id[rank] - start_id)*_tr_stride + tr_id;

        l = _interface->tr_index[tr_id] + n_ent_rank_tr[k];

        _interface->local_num[l] = local_num;
        _interface->distant_num[l] = distant_num;

        n_ent_rank_tr[k] += 1;

      }
    }

  }

  /* n_ent_rank will now be used as a position counter for new interfaces */

  BFT_FREE(n_ent_rank_tr);
  BFT_FREE(interface_id);
}

/*----------------------------------------------------------------------------
 * Global ordering associated with an I/O numbering structure.
 *
 * The global_num values should be sorted, but need not be contiguous.
 *
 * parameters:
 *   this_interface_set  <-> pointer to structure that should be updated
 *   n_ent               <-- local number of entities
 *   global_num          <-- global number (id) associated with each entity
 *   comm                <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvm_interface_add_global_equiv(fvm_interface_set_t  *this_interface_set,
                                fvm_lnum_t            n_ent,
                                fvm_gnum_t            global_num[],
                                MPI_Comm              comm)
{
  fvm_gnum_t  global_max;
  fvm_lnum_t  i, j, n_ent_recv, n_ent_send;
  size_t      slice_size;
  int         size, rank;

  _per_slice_equiv_t  e;

  int         *send_count = NULL, *recv_count = NULL;
  int         *send_shift = NULL, *recv_shift = NULL;
  fvm_lnum_t  *equiv_id = NULL;
  fvm_lnum_t  *equiv_send, *equiv_recv = NULL;
  fvm_lnum_t  *recv_num = NULL, *send_num = NULL;
  fvm_gnum_t  *recv_global_num = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  /* Get temporary maximum global number value */

  global_max = _global_num_max(n_ent,
                               global_num,
                               comm);

  /* slice_size = ceil(global_max/size) */

  slice_size = global_max / size;
  if (global_max % size > 0)
    slice_size += 1;

  assert(sizeof(fvm_gnum_t) >= sizeof(fvm_lnum_t));

  BFT_MALLOC(send_count, size, int);
  BFT_MALLOC(recv_count, size, int);

  BFT_MALLOC(send_shift, size + 1, int);
  BFT_MALLOC(recv_shift, size + 1, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < n_ent; i++)
    send_count[(global_num[i] - 1) / slice_size] += 1;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the global_num[] array to the correct processors */

  n_ent_recv = recv_shift[size];

  BFT_MALLOC(recv_global_num, n_ent_recv, fvm_gnum_t);
  BFT_MALLOC(recv_num, n_ent_recv, fvm_lnum_t);

  MPI_Alltoallv(global_num, send_count, send_shift, FVM_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, FVM_MPI_GNUM, comm);

  /* We also need to send the corresponding local numbers */

  n_ent_send = send_shift[size];

  BFT_MALLOC(send_num, n_ent_send, fvm_lnum_t);

  for (i = 0; i < n_ent_send; i++)
    send_num[i] = i+1;

  MPI_Alltoallv(send_num, send_count, send_shift, FVM_MPI_LNUM,
                recv_num, recv_count, recv_shift, FVM_MPI_LNUM, comm);

  BFT_FREE(send_num);

  /* Build equivalence data */

  if (n_ent_recv > 0)
    BFT_MALLOC(equiv_id, n_ent_recv, fvm_lnum_t);

  e = _slice_global_num_to_equiv(size,
                                 n_ent_recv,
                                 recv_shift,
                                 recv_global_num,
                                 recv_num,
                                 equiv_id);

  /* Now free Memory */

  BFT_FREE(recv_num);
  BFT_FREE(recv_global_num);

  /* Now that equivalences are marked, count for each rank; for each
     equivalence, we will need to send the local entity number,
     the number of equivalent entities (e.multiple[...] - 1),
     and the corresponding entity numbers and ranks,
     for a total of 1 + 1 + 2*(e.multiple[...] - 1)
     = 2*(e.multiple[...]) values. */

  for (rank = 0; rank < size; rank++) {
    send_count[rank] = 0;
    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {
      if (equiv_id[i] > -1) {
        size_t e_id = equiv_id[i];
        send_count[rank] += 2*(e.shift[e_id+1] - e.shift[e_id]);
      }
    }
  }

  for (rank = 0; rank < size; rank++)
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];

  /* Now prepare new send buffer */

  n_ent_send = send_shift[size];
  BFT_MALLOC(equiv_send, n_ent_send, fvm_lnum_t);

  for (rank = 0; rank < size; rank++) {

    send_count[rank] = 0; /* reset, will be re-incremented */

    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

      if (equiv_id[i] > -1) {

        fvm_lnum_t *equiv_send_p
          = equiv_send + send_shift[rank] + send_count[rank];

        fvm_lnum_t  e_id = equiv_id[i];
        fvm_lnum_t  k = 2;
        const int         multiple = e.shift[e_id+1] - e.shift[e_id];
        const int        *rank_p = e.rank + e.shift[e_id];
        const fvm_lnum_t *num_p  = e.num  + e.shift[e_id];

        send_count[rank] += 2*multiple;

        for (j = 0; j < multiple; j++) {
          if (rank_p[j] == rank) {
            equiv_send_p[0] = num_p[j];
            equiv_send_p[1] = multiple - 1;
          }
          else {
            equiv_send_p[k++] = num_p[j];
            equiv_send_p[k++] = rank_p[j];
          }
        }

      }
    }
  }

  /* Free temporary (slice) equivalence info */

  e.count = 0;
  BFT_FREE(e.shift);
  BFT_FREE(e.rank);
  if (e.tr_id != NULL)
    BFT_FREE(e.tr_id);
  BFT_FREE(e.num);

  BFT_FREE(equiv_id);

  /* Send prepared slice data to destination rank */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  n_ent_recv = recv_shift[size];

  BFT_MALLOC(equiv_recv, n_ent_recv, fvm_lnum_t);

  MPI_Alltoallv(equiv_send, send_count, send_shift, FVM_MPI_LNUM,
                equiv_recv, recv_count, recv_shift, FVM_MPI_LNUM, comm);

  /* At this stage, MPI operations are finished; we may release
     the corresponding counts and indexes */

  BFT_FREE(equiv_send);

  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);

  /* Add interface */

  _interfaces_from_flat_equiv(this_interface_set,
                              1,
                              n_ent_recv,
                              equiv_recv);

  BFT_FREE(equiv_recv);
}

/*----------------------------------------------------------------------------
 * Exchange communication counts.
 *
 * If n_counts > 1, counts are not interlaced, so count arrays are
 * contiguous: {c1[rank], rank = 0, size}{c2[rank], rank = 0, size}.
 *
 * parameters:
 *   n_counts        <-- number of counts per processor (1 or 2)
 *   size            <-- communicator sizentities
 *   send_count      <-- send element count for each rank [size*n_counts]
 *   recv_count      --> receive element count for each rank [size*n_counts]
 *   comm            <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_exchange_counts(int       n_counts,
                 int       size,
                 int       send_count[],
                 int       recv_count[],
                 MPI_Comm  comm)
{
  int i;
  int *count_tmp = NULL;

  /* Interlace parallel and periodic send counts */

  if (n_counts == 2) {

    BFT_MALLOC(count_tmp, size, int);

    memcpy(count_tmp, send_count + size, size*sizeof(int));

    for (i = size -1; i >= 0; i--)
      send_count[i*2] = send_count[i];

    for (i = 0; i < size; i++)
      send_count[i*2+1] = count_tmp[i];

  }

  MPI_Alltoall(send_count, n_counts, MPI_INT,
               recv_count, n_counts, MPI_INT, comm);


  /* De-interlace parallel and periodic send counts */

  if (n_counts == 2) {

    for (i = 0; i < size; i++) {
      count_tmp[i] = send_count[i*2+1];
      send_count[i] = send_count[i*2];
    }

    memcpy(send_count + size, count_tmp, size*sizeof(int));

    for (i = 0; i < size; i++) {
      count_tmp[i] = recv_count[i*2+1];
      recv_count[i] = recv_count[i*2];
    }

    memcpy(recv_count + size, count_tmp, size*sizeof(int));

    BFT_FREE(count_tmp);

  }

}

/*----------------------------------------------------------------------------
 * Exchange periodic couple info between processors providing the data
 * and processors handling the related global numbering interval slices.
 *
 * parameters:
 *   slice_size          <-- size of the slice handled by each processor
 *   periodicity         <-- periodicity information (NULL if none)
 *   n_periodic_lists    <-- number of periodic lists (may be local)
 *   n_periodic_couples  <-- number of periodic couples associated with
 *                           each periodic list
 *   periodic_couples    <-- array indicating periodic couples (using
 *                           global numberings) for each list
 *   recv_couples        <-- couple information received: for each couple,
 *                           {global number of local element,
 *                            global number of periodic element,
 *                            transform id}
 *   send_count          --> local number of values to send to each rank
 *   comm                <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_count_periodic_couple_exchange(size_t                    slice_size,
                                const fvm_periodicity_t  *periodicity,
                                int                       n_periodic_lists,
                                const fvm_lnum_t          n_periodic_couples[],
                                const fvm_gnum_t   *const periodic_couples[],
                                int                       send_count[],
                                MPI_Comm                  comm)
{
  int         list_id;
  int         size, rank;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  if (periodicity == NULL)
    return;

  /* Compute list sizes to send to distant processors */

  for (list_id = 0; list_id < n_periodic_lists; list_id++) {

    fvm_lnum_t couple_id;
    fvm_gnum_t num_1, num_2;
    int rank_1, rank_2;

    const fvm_lnum_t  _n_periodic_couples = n_periodic_couples[list_id];
    const fvm_gnum_t  *_periodic_couples = periodic_couples[list_id];

    for (couple_id = 0; couple_id < _n_periodic_couples; couple_id++) {
      num_1 = _periodic_couples[couple_id*2];
      num_2 = _periodic_couples[couple_id*2 + 1];
      rank_1 = (num_1 - 1) / slice_size;
      rank_2 = (num_2 - 1) / slice_size;
      send_count[rank_1] += 3;
      send_count[rank_2] += 3;
    }

  }

}

/*----------------------------------------------------------------------------
 * Build eventual couples belonging to combined periodicities.
 *
 * parameters:
 *   slice_size        <-- size of the slice handled by each processor
 *   periodicity       <-- periodicity information (NULL if none)
 *   n_slice_couples   <-> number of couples in current slice
 *   slice_couples     <-> couple information for this rank: for each couple,
 *                         {global number of local element,
 *                          global number of periodic element,
 *                          transform id}
 *   comm              <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_combine_periodic_couples(size_t                     slice_size,
                          const fvm_periodicity_t   *periodicity,
                          fvm_lnum_t                *n_slice_couples,
                          fvm_gnum_t               **slice_couples,
                          MPI_Comm                   comm)
{
  int  size;
  int  i, n_tr, level, n_recv;

  fvm_lnum_t  j, k, start_id, end_id;

  int  *send_count = NULL, *recv_count = NULL;
  int  *send_shift = NULL, *recv_shift = NULL;
  int  *tr_reverse_id = NULL;

  fvm_gnum_t  *send_couples = NULL;

  fvm_lnum_t  _n_slice_couples = *n_slice_couples;
  fvm_gnum_t  *_slice_couples = *slice_couples;

  assert (periodicity != NULL);

  /* Initialization */

  MPI_Comm_size(comm, &size);

  /* Build periodicity related arrays for quick access */

  n_tr = fvm_periodicity_get_n_transforms(periodicity);

  BFT_MALLOC(tr_reverse_id, n_tr, int);

  for (i = 0; i < n_tr; i++)
    tr_reverse_id[i] = fvm_periodicity_get_reverse_id(periodicity, i);

  /* Arrays for MPI_AlltoAll */

  BFT_MALLOC(send_count, size, int);
  BFT_MALLOC(recv_count, size, int);
  BFT_MALLOC(send_shift, size+1, int);
  BFT_MALLOC(recv_shift, size+1, int);

  /* Loop on combination levels */

  for (level = 1;
       level < fvm_periodicity_get_n_levels(periodicity);
       level++) {

    int  n_rows = 0;
    int  *tr_combine = NULL;

    _transform_combine_info(periodicity,
                            level,
                            &n_rows,
                            &tr_combine);

    /* Count values to exchange */

    for (i = 0; i < size; i++)
      send_count[i] = 0;

    start_id = 0;
    end_id = 1;

    while (end_id < _n_slice_couples) {

      if (_slice_couples[start_id*3] == _slice_couples[end_id*3]) {

        end_id++;
        while (end_id < _n_slice_couples) {
          if (_slice_couples[end_id*3] != _slice_couples[start_id*3])
            break;
          end_id++;
        }

        for (j = start_id; j < end_id; j++) {
          for (k = j+1; k < end_id; k++) {

            fvm_gnum_t num_1 = _slice_couples[j*3 + 1];
            fvm_gnum_t num_2 = _slice_couples[k*3 + 1];
            int rank_1 = (num_1 - 1) / slice_size;
            int rank_2 = (num_2 - 1) / slice_size;
            int tr_1 = tr_reverse_id[_slice_couples[j*3 + 2]];
            int tr_2 = _slice_couples[k*3 + 2];

            if (tr_combine[tr_1 *n_rows + tr_2] > -1) {
              send_count[rank_1] += 3;
              send_count[rank_2] += 3;
            }

          }
        }

      }

      start_id = end_id;
      end_id += 1;
    }

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    send_shift[0] = 0;
    recv_shift[0] = _n_slice_couples*3;

    for (i = 0; i < size; i++) {
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
      send_count[i] = 0;
    }

    BFT_MALLOC(send_couples, send_shift[size], fvm_gnum_t);

    /* Now exchange combined couples */

    start_id = 0;
    end_id = 1;

    while (end_id < _n_slice_couples) {

      if (_slice_couples[start_id*3] == _slice_couples[end_id*3]) {

        end_id++;
        while (end_id < _n_slice_couples) {
          if (_slice_couples[end_id*3] != _slice_couples[start_id*3])
            break;
          end_id++;
        }

        for (j = start_id; j < end_id; j++) {
          for (k = j+1; k < end_id; k++) {

            fvm_gnum_t num_1 = _slice_couples[j*3 + 1];
            fvm_gnum_t num_2 = _slice_couples[k*3 + 1];
            int rank_1 = (num_1 - 1) / slice_size;
            int rank_2 = (num_2 - 1) / slice_size;
            int tr_1 = tr_reverse_id[_slice_couples[j*3 + 2]];
            int tr_2 = _slice_couples[k*3 + 2];
            int tr_c = tr_combine[tr_1 *n_rows + tr_2];

            if (tr_c > -1) {

              fvm_gnum_t *_send_couples;

              _send_couples
                = send_couples + send_shift[rank_1] + send_count[rank_1];

              _send_couples[0] = num_1;
              _send_couples[1] = num_2;
              _send_couples[2] = tr_c;

              send_count[rank_1] += 3;

              _send_couples
                = send_couples + send_shift[rank_2] + send_count[rank_2];

              _send_couples[0] = num_2;
              _send_couples[1] = num_1;
              _send_couples[2] = tr_reverse_id[tr_c];

              send_count[rank_2] += 3;

            }

          }
        }

      }

      start_id = end_id;
      end_id += 1;
    }

    BFT_FREE(tr_combine);

    n_recv = recv_shift[size] - (_n_slice_couples*3);

    if (n_recv > 0)
      BFT_REALLOC(_slice_couples, recv_shift[size], fvm_gnum_t);

    MPI_Alltoallv(send_couples, send_count, send_shift, FVM_MPI_GNUM,
                  _slice_couples, recv_count, recv_shift, FVM_MPI_GNUM, comm);

    BFT_FREE(send_couples);

    /* Finally, merge additional couples received for slice with
       existing periodicity information */

    if (n_recv > 0) {

      assert(n_recv % 3 == 0);

      _n_slice_couples += n_recv / 3;

      /* Sort and remove dumplicates to update slice periodicity info */

      _sort_periodic_couples(&_n_slice_couples, &_slice_couples);

      *n_slice_couples = _n_slice_couples;
      *slice_couples = _slice_couples;

    }

  }

  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);
  BFT_FREE(send_count);
  BFT_FREE(recv_count);

  BFT_FREE(tr_reverse_id);
}

/*----------------------------------------------------------------------------
 * Exchange periodic couple info between processors providing the data
 * and processors handling the related global numbering interval slices.
 *
 * _count_periodic_couple_exchange() should have been called first.
 *
 * Note that the array pointed to by slice_couples is allocated here,
 * and must be freed by the calling code.
 *
 * parameters:
 *   slice_size          <-- size of the slice handled by each processor
 *   send_count          <-> local number of values to send to each rank
 *   recv_count          <-- local number of values to receive from each rank
 *   periodicity         <-- periodicity information (NULL if none)
 *   n_periodic_lists    <-- number of periodic lists (may be local)
 *   periodicity_num     <-- periodicity number (1 to n) associated with
 *                           each periodic list (primary periodicities only)
 *   n_periodic_couples  <-- number of periodic couples associated with
 *                           each periodic list
 *   periodic_couples    <-- array indicating periodic couples (using
 *                           global numberings) for each list
 *   n_slice_couples     --> number of couples in current slice
 *   slice_couples       --> couple information for slice: for each couple,
 *                           {global number of local element,
 *                            global number of periodic element,
 *                            transform id}
 *   comm                <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_exchange_periodic_couples(size_t                    slice_size,
                           int                       send_count[],
                           int                       recv_count[],
                           const fvm_periodicity_t  *periodicity,
                           int                       n_periodic_lists,
                           const int                 periodicity_num[],
                           const fvm_lnum_t          n_periodic_couples[],
                           const fvm_gnum_t   *const periodic_couples[],
                           fvm_lnum_t               *n_slice_couples,
                           fvm_gnum_t              **slice_couples,
                           MPI_Comm                  comm)
{
  fvm_lnum_t  i, j;
  int         list_id;
  int         size, rank;

  fvm_lnum_t   _n_slice_couples = 0;
  int         *send_shift = NULL, *recv_shift = NULL;
  fvm_gnum_t  *send_couples = NULL, *recv_couples = NULL;

  /* Initialization */

  *n_slice_couples = 0;
  *slice_couples = NULL;

  MPI_Comm_size(comm, &size);

  BFT_MALLOC(send_shift, size + 1, int);
  BFT_MALLOC(recv_shift, size + 1, int);

  /* build shift arrays for MPI_Alltoallv and reset send_count */

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
    send_count[rank] = 0;
  }

  BFT_MALLOC(send_couples, send_shift[size], fvm_gnum_t);
  BFT_MALLOC(recv_couples, recv_shift[size], fvm_gnum_t);

  _n_slice_couples = recv_shift[size]/3;

  /* Prepare lists to send to distant processors */

  for (list_id = 0; list_id < n_periodic_lists; list_id++) {

    fvm_lnum_t couple_id;
    fvm_gnum_t num_1, num_2;
    int rank_1, rank_2;

    const int external_num = periodicity_num[list_id];
    const int direct_id = fvm_periodicity_get_transform_id(periodicity,
                                                           external_num,
                                                           1);
    const int reverse_id = fvm_periodicity_get_transform_id(periodicity,
                                                            external_num,
                                                            -1);

    const fvm_lnum_t  _n_periodic_couples = n_periodic_couples[list_id];
    const fvm_gnum_t  *_periodic_couples = periodic_couples[list_id];

    assert(direct_id >= 0 && reverse_id >= 0);

    for (couple_id = 0; couple_id < _n_periodic_couples; couple_id++) {

      num_1 = _periodic_couples[couple_id*2];
      num_2 = _periodic_couples[couple_id*2 + 1];

      rank_1 = (num_1 - 1) / slice_size;
      rank_2 = (num_2 - 1) / slice_size;

      i = send_shift[rank_1] + send_count[rank_1];

      send_couples[i] = num_1;
      send_couples[i+1] = num_2;
      send_couples[i+2] = direct_id;

      send_count[rank_1] += 3;

      j = send_shift[rank_2] + send_count[rank_2];

      send_couples[j] = num_2;
      send_couples[j+1] = num_1;
      send_couples[j+2] = reverse_id;

      send_count[rank_2] += 3;

    }

  }

#if defined(DEBUG) && !defined(NDEBUG)
  for (rank = 0; rank < size; rank++) {
    assert(send_shift[rank + 1] - send_shift[rank] == send_count[rank]);
  }
#endif

  /* Exchange data */

  MPI_Alltoallv(send_couples, send_count, send_shift, FVM_MPI_GNUM,
                recv_couples, recv_count, recv_shift, FVM_MPI_GNUM, comm);

  BFT_FREE(recv_shift);
  BFT_FREE(send_shift);

  BFT_FREE(send_couples);

  /* Sort periodic couples by local correspondant, remove duplicates */

  _sort_periodic_couples(&_n_slice_couples,
                         &recv_couples);

  /* Set return values */

  *n_slice_couples = _n_slice_couples;
  *slice_couples = recv_couples;
}

/*----------------------------------------------------------------------------
 * Associate slice ids for periodic couples.
 *
 * If a global number appears multiple times in a slice, the lowest
 * occurence id is returned.
 *
 * parameters:
 *   n_slice_elements    <-- number of elements in slice
 *   order               <-- slice ordering by global number received
 *   slice_global_num    <-- global numbering received
 *   n_slice_couples     <-- number of couples in current slice
 *   stride              <-- stride for global number of local element
 *                           in slice_couples[]
 *   slice_couples       <-- couple information for slice: for each couple,
 *                           {global number of local element,
 *                            global number of periodic element,
 *                            transform id} if stride = 3,
 *                           {global number of local element} if stride = 1
 *   couple_slice_id     --> id in slice of local couple element
 *   comm                <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_periodic_couples_slice_id(fvm_lnum_t         n_slice_elements,
                           const fvm_lnum_t   order[],
                           const fvm_gnum_t   slice_global_num[],
                           fvm_lnum_t         n_slice_couples,
                           int                stride,
                           const fvm_gnum_t   slice_couples[],
                           fvm_lnum_t         couple_slice_id[])
{
  fvm_lnum_t   couple_id;

  /* Initialization */

  assert(stride == 3 || stride == 1);

  if (n_slice_couples == 0)
    return;

  /* Use binary search */

  for (couple_id = 0; couple_id < n_slice_couples; couple_id++) {

    fvm_gnum_t num_cmp;
    fvm_lnum_t start_id = 0;
    fvm_lnum_t end_id = n_slice_elements - 1;
    fvm_lnum_t mid_id = (end_id -start_id) / 2;

    const fvm_gnum_t num_1 = slice_couples[couple_id*stride];

    /* use binary search */

    while (start_id <= end_id) {
      num_cmp = slice_global_num[order[mid_id]];
      if (num_cmp < num_1)
        start_id = mid_id + 1;
      else if (num_cmp > num_1)
        end_id = mid_id - 1;
      else
        break;
      mid_id = start_id + ((end_id -start_id) / 2);
    }

    /* In case of multiple occurences, find lowest one */

    while (mid_id > 0 &&  slice_global_num[order[mid_id-1]] == num_1)
      mid_id--;

    assert(slice_global_num[order[mid_id]] == num_1);

    couple_slice_id[couple_id] = order[mid_id];
  }

}

/*----------------------------------------------------------------------------
 * Find rank associated with a given position in a slice.
 *
 * parameters:
 *   n_slices          <-- number of slices (communicator size)
 *   slice_shift       <-- shift in received data per rank (size: n_ranks+1)
 *   slice_id          <-- local numbering received
 *
 * returns:
 *   originating rank associated with id in slice
 *----------------------------------------------------------------------------*/

static int
_rank_by_slice_id(int                   n_slices,
                  const fvm_lnum_t      slice_shift[],
                  fvm_lnum_t            slice_id)
{
  fvm_lnum_t start_id = 0;
  fvm_lnum_t end_id = n_slices - 1;
  fvm_lnum_t mid_id = (end_id -start_id) / 2;

  /* Use binary search */

  while (start_id <= end_id) {
    if (slice_shift[mid_id + 1] <= slice_id)
      start_id = mid_id + 1;
    else if (slice_shift[mid_id] > slice_id)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  return mid_id;
}

/*----------------------------------------------------------------------------
 * Exchange periodic couple info between processors providing the data
 * and processors handling the related global numbering interval slices.
 *
 * _periodic_couples_slice_id() should have been called first.
 *
 * parameters:
 *   slice_size          <-- size of the slice handled by each processor
 *   equiv_id            <-- equivalence id for each slice element (-1 if none)
 *   equiv               <-- temporary equivalence structure for slice
 *   n_slice_couples     <-- number of couples in current slice
 *   slice_couples       <-- couple information for slice: for each couple,
 *                           {global number of local element,
 *                            global number of periodic element,
 *                            transform id}
 *   couple_slice_id     <-- local id in slice
 *   send_count          --> local number of values to send to each rank
 *   slice_count         --> local number of values to receive from each rank
 *   comm                <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_count_periodic_equiv_exchange(size_t                     slice_size,
                               const fvm_lnum_t           equiv_id[],
                               const _per_slice_equiv_t  *equiv,
                               fvm_lnum_t                 n_slice_couples,
                               const fvm_gnum_t           slice_couples[],
                               const int                  couple_slice_id[],
                               int                        send_count[],
                               int                        recv_count[],
                               MPI_Comm                   comm)
{
  int          size;
  int          rank;
  fvm_lnum_t   couple_id;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  if (equiv != NULL && equiv_id != NULL) {

    /* Compute list sizes to send to distant processors */

    for (couple_id = 0; couple_id < n_slice_couples; couple_id++) {

      int e_mult;
      fvm_gnum_t num_2 = slice_couples[couple_id*3 + 1];
      fvm_lnum_t e_id = equiv_id[couple_slice_id[couple_id]];
      int rank_2 = (num_2 - 1) / slice_size;

      if (e_id > -1)
        e_mult = equiv->shift[e_id +1] - equiv->shift[e_id];
      else
        e_mult = 1;

      send_count[rank_2] += 3 + 2*e_mult;

    }

  }
  else { /* if (equiv == NULL || equiv_id == NULL) */

    for (couple_id = 0; couple_id < n_slice_couples; couple_id++) {

      fvm_gnum_t num_2 = slice_couples[couple_id*3 + 1];
      int rank_2 = (num_2 - 1) / slice_size;

      send_count[rank_2] += 5;

    }

  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);
}

/*----------------------------------------------------------------------------
 * Exchange periodic couple info between processors providing the data
 * and processors handling the related global numbering interval slices.
 *
 * parameters:
 *   slice_size          <-- size of the slice handled by each processor
 *   slice_shift         <-- shift in received data per rank (size: n_ranks+1)
 *   slice_global_num    <-- global numbering received
 *   slice_num           <-- local numbering received
 *   equiv_id            <-- equivalence id for each slice element (-1 if none)
 *   equiv               <-- temporary equivalence structure for slice
 *   periodicity         <-- periodicity information (NULL if none)
 *   n_slice_couples     <-- number of couples in current slice
 *   slice_couples       <-- couple information for slice: for each couple,
 *                           {global number of local element,
 *                            global number of periodic element,
 *                            transform id}
 *   send_count          <-> local number of values to send to each rank
 *   recv_count          <-> local number of values to receive from each rank
 *   send_shift          --- local index of values to send to each rank
 *                           (re-used work array, size n_ranks + 1)
 *   comm                <-- associated MPI communicator
 *
 * returns:
 *   structure defining a temporary list of periodic interfaces
 *----------------------------------------------------------------------------*/

static _per_slice_period_t
_exchange_periodic_equiv(size_t                     slice_size,
                         const fvm_lnum_t           slice_shift[],
                         const fvm_gnum_t           slice_global_num[],
                         const fvm_lnum_t           slice_num[],
                         const fvm_lnum_t           equiv_id[],
                         const _per_slice_equiv_t  *equiv,
                         const fvm_periodicity_t   *periodicity,
                         fvm_lnum_t                 n_slice_couples,
                         const fvm_gnum_t           slice_couples[],
                         int                        send_count[],
                         int                        recv_count[],
                         int                        send_shift[],
                         MPI_Comm                   comm)
{
  int          tr_id;
  int          size;
  int          rank;
  fvm_lnum_t   n_slice_elements, couple_id;

  int  n_tr = 0;
  size_t  recv_size = 0;
  int *couple_slice_id = NULL;
  int *reverse_tr_id = NULL;
  int *recv_shift = NULL;
  fvm_lnum_t *order = NULL;
  fvm_gnum_t *equiv_send = NULL, *equiv_recv = NULL;
  fvm_gnum_t *slice_recv_num = NULL;

  _per_slice_period_t pe;

  /* Initialize return structure */

  pe.count = 0;
  pe.slice_id = NULL;
  pe.tr_id = NULL;
  pe.shift = NULL;
  pe.rank  = NULL;
  pe.num   = NULL;

  if (periodicity == NULL)
    return pe;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  n_slice_elements = slice_shift[size];

  /* Build ordering array for binary search */

  order = fvm_order_local(NULL, slice_global_num, n_slice_elements);

  /* Associate id in slice for periodic couples prior to sending */

  BFT_MALLOC(couple_slice_id, n_slice_couples, fvm_lnum_t);

  _periodic_couples_slice_id(n_slice_elements,
                             order,
                             slice_global_num,
                             n_slice_couples,
                             3,
                             slice_couples,
                             couple_slice_id);

  /* build count and shift arrays for MPI_Alltoallv */

  _count_periodic_equiv_exchange(slice_size,
                                 equiv_id,
                                 equiv,
                                 n_slice_couples,
                                 slice_couples,
                                 couple_slice_id,
                                 send_count,
                                 recv_count,
                                 comm);

  BFT_MALLOC(recv_shift, size + 1, int);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  /* arrays to exchange; most of the exchanged data is of type int or
     fvm_lnum_t, using only positive values. For each periodic couple,
     one value of type fvm_gnum_t is exchanged, so to group all communication
     in one MPI call, all data is exchanged as fvm_gnum_t, even if this
     means larger messages if fvm_gnum_t is larger than fvm_lnum_t.
     As the number of elements per couple is variable (depending on prior
     equivalence info), using an MPI datatype to mix int and fvm_gnum
     types rather than casting all to fvm_gnum_t is not feasible */

  BFT_MALLOC(equiv_send, send_shift[size], fvm_gnum_t);
  BFT_MALLOC(equiv_recv, recv_shift[size], fvm_gnum_t);

  /* temporary array to find reverse transforms */

  n_tr = fvm_periodicity_get_n_transforms(periodicity);

  BFT_MALLOC(reverse_tr_id, n_tr, int);

  for (tr_id = 0; tr_id < n_tr; tr_id++)
    reverse_tr_id[tr_id] = fvm_periodicity_get_reverse_id(periodicity, tr_id);

  /* Reset send count */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  if (equiv != NULL && equiv_id != NULL) {

    /* Compute list sizes to send to distant processors */

    for (couple_id = 0; couple_id < n_slice_couples; couple_id++) {

      const fvm_gnum_t num_2 = slice_couples[couple_id*3 + 1];
      const fvm_lnum_t local_id = couple_slice_id[couple_id];
      const fvm_lnum_t e_id = equiv_id[local_id];
      const int rank_2 = (num_2 - 1) / slice_size;

      size_t i = send_shift[rank_2] + send_count[rank_2];

      if (e_id > -1) {

        int j;
        int j_start = equiv->shift[e_id];
        int j_end = equiv->shift[e_id + 1];

        equiv_send[i++] = j_end - j_start;
        equiv_send[i++] = num_2;
        equiv_send[i++] = reverse_tr_id[slice_couples[couple_id*3 + 2]];

        for (j = j_start; j < j_end; j++) {
          equiv_send[i++] = equiv->rank[j];
          equiv_send[i++] = equiv->num[j];
        }

        send_count[rank_2] += 3 + 2*(j_end - j_start);
      }
      else {

        equiv_send[i++] = 1;
        equiv_send[i++] = num_2;
        equiv_send[i++] = reverse_tr_id[slice_couples[couple_id*3 + 2]];
        equiv_send[i++] = _rank_by_slice_id(size, slice_shift, local_id);
        equiv_send[i++] = slice_num[local_id];

        send_count[rank_2] += 5;
      }

    }

  }
  else { /* if (equiv == NULL || equiv_id == NULL) */

    for (couple_id = 0; couple_id < n_slice_couples; couple_id++) {

      const fvm_gnum_t num_2 = slice_couples[couple_id*3 + 1];
      const fvm_lnum_t local_id = couple_slice_id[couple_id];
      const int rank_2 = (num_2 - 1) / slice_size;

      size_t i = send_shift[rank_2] + send_count[rank_2];

      equiv_send[i++] = 1;
      equiv_send[i++] = num_2;
      equiv_send[i++] = reverse_tr_id[slice_couples[couple_id*3 + 2]];
      equiv_send[i++] = _rank_by_slice_id(size, slice_shift, local_id);
      equiv_send[i++] = slice_num[local_id];

      send_count[rank_2] += 5;

    }

  }

  BFT_FREE(couple_slice_id);
  BFT_FREE(reverse_tr_id);

  /* Parallel exchange */

  MPI_Alltoallv(equiv_send, send_count, send_shift, FVM_MPI_GNUM,
                equiv_recv, recv_count, recv_shift, FVM_MPI_GNUM, comm);

  recv_size = recv_shift[size];

  /* Free memory */

  BFT_FREE(recv_shift);
  BFT_FREE(equiv_send);

  /* Build return structure */

  {
    size_t i, j, k, l, e_mult;

    pe.count = 0;
    i = 0;
    j = 0;

    while (i < recv_size) {
      pe.count += 1;
      j += equiv_recv[i];
      i += 3 + 2*equiv_recv[i];
    }

    BFT_MALLOC(slice_recv_num, pe.count, fvm_gnum_t);

    BFT_MALLOC(pe.tr_id, pe.count, int);

    BFT_MALLOC(pe.shift, pe.count + 1, int);

    BFT_MALLOC(pe.rank, j, int);
    BFT_MALLOC(pe.num, j, int);

    pe.shift[0] = 0;

    for (i = 0, j = 0, k = 0, l = 0; i < (size_t)(pe.count); i++) {

      e_mult = equiv_recv[k++];
      slice_recv_num[i] = equiv_recv[k++];
      pe.tr_id[i] = equiv_recv[k++];

      for (l = 0; l < e_mult; l++) {
        pe.rank[j] = equiv_recv[k++];
        pe.num[j] = equiv_recv[k++];
        pe.shift[i+1] = ++j;
      }

    }

  }

  BFT_FREE(equiv_recv);

  /* Associate id in slice for received periodic equivalences */

  BFT_MALLOC(pe.slice_id, pe.count, fvm_lnum_t);

  _periodic_couples_slice_id(n_slice_elements,
                             order,
                             slice_global_num,
                             pe.count,
                             1,
                             slice_recv_num,
                             pe.slice_id);

  /* Free remaining working arrays */

  BFT_FREE(slice_recv_num);
  BFT_FREE(couple_slice_id);
  BFT_FREE(order);

  return pe;
}

/*----------------------------------------------------------------------------
 * Merge periodic equivalent interface info with slice equivalence info.
 *
 * Expands slice equivalence info, and frees temporary list of periodic
 * interfaces.
 *
 * parameters:
 *   n_slices          <-- number of slices (communicator size)
 *   slice_shift       <-- shift in received data per rank (size: n_ranks+1)
 *   slice_num         <-- local numbering received
 *   equiv_id          <-> equivalence id for each slice element (-1 if none)
 *   equiv             <-> temporary equivalence structure for slice
 *   perio_equiv       <-> temporary list of periodic interfaces
 *
 * returns:
 *   structure defining a temporary equivalence structure
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 2 /* Crash with O3 on IA64 with icc 9.1 20070320 */
#endif

static void
_merge_periodic_equiv(int                   n_slices,
                      const fvm_lnum_t      slice_shift[],
                      const fvm_lnum_t      slice_num[],
                      fvm_lnum_t            equiv_id[],
                      _per_slice_equiv_t   *equiv,
                      _per_slice_period_t  *perio_equiv)
{
  int i;
  size_t j;

  int rank;
  int old_count, new_count;
  size_t new_size;

  int *eq_mult = NULL;
  int *new_shift = NULL;

  _per_slice_period_t *pe = perio_equiv;

  assert(equiv != NULL);

  if (perio_equiv == NULL)
    return;

  old_count = equiv->count;

  /* Note that by construction, the global numbers of elements appearing
     in the original (parallel) equivalence must appear multiple times
     in the slice (the original equivalences are built using this), and
     the global numbers of elements appearing only in the periodic
     equivalence must appear only once; thus only one equiv_id[] value
     needs to be updated when appending purely periodic equivalences */

  for (i = 0, new_count = old_count; i < pe->count; i++) {
    if (equiv_id[pe->slice_id[i]] == -1)
      equiv_id[pe->slice_id[i]] = new_count++;
  }

  BFT_MALLOC(eq_mult, new_count, int);

  for (i = 0; i < old_count; i++)
    eq_mult[i] = equiv->shift[i+1] - equiv->shift[i];

  for (i = old_count; i < new_count; i++)
    eq_mult[i] = 0;

  for (i = 0; i < pe->count; i++) {
    if (eq_mult[equiv_id[pe->slice_id[i]]] == 0)
      eq_mult[equiv_id[pe->slice_id[i]]] += pe->shift[i+1] - pe->shift[i] + 1;
    else
      eq_mult[equiv_id[pe->slice_id[i]]] += pe->shift[i+1] - pe->shift[i];
  }

  /* Build new (merged) index, resetting eq_mult to use as a counter */

  BFT_MALLOC(new_shift, new_count+1, int);

  new_shift[0] = 0;

  for (i = 0; i < new_count; i++) {
    assert(eq_mult[i] > 0);
    new_shift[i+1] = new_shift[i] + eq_mult[i];
    eq_mult[i] = 0;
  }

  new_size = new_shift[new_count];

  /* Expand previous periodicity info */
  /*----------------------------------*/

  equiv->count = new_count;

  if (old_count > 0) {

    int k;
    int *new_rank = NULL, *new_num = NULL;

    BFT_MALLOC(new_rank, new_size, int);

    for (i = 0; i < old_count; i++) {
      eq_mult[i] = equiv->shift[i+1] - equiv->shift[i];
      for (k = 0; k < eq_mult[i]; k++)
        new_rank[new_shift[i] + k] = equiv->rank[equiv->shift[i] + k];
    }

    BFT_FREE(equiv->rank);
    equiv->rank = new_rank;

    BFT_MALLOC(new_num, new_size, int);

    for (i = 0; i < old_count; i++) {
      for (k = 0; k < eq_mult[i]; k++)
        new_num[new_shift[i] + k] = equiv->num[equiv->shift[i] + k];
    }
    BFT_FREE(equiv->num);
    equiv->num = new_num;

    if (equiv->tr_id != NULL) {

      int *new_tr_id = NULL;
      BFT_MALLOC(new_tr_id, new_size, int);
      for (i = 0; i < old_count; i++) {
        for (k = 0; k < eq_mult[i]; k++)
          new_tr_id[new_shift[i] + k] = equiv->tr_id[equiv->shift[i] + k];
      }
      BFT_FREE(equiv->tr_id);
      equiv->tr_id = new_tr_id;

    }

    /* All is expanded at this stage, so old index may be replaced */

    BFT_FREE(equiv->shift);
    equiv->shift = new_shift;

  }
  else {

    BFT_FREE(equiv->shift);
    equiv->shift = new_shift;
    BFT_MALLOC(equiv->rank, new_size, int);
    BFT_MALLOC(equiv->num, new_size, int);

  }

  if (equiv->tr_id == NULL) {
    BFT_MALLOC(equiv->tr_id, new_size, int);
    for (j = 0; j < new_size; j++)
      equiv->tr_id[j] = 0;
  }

  /* Now insert periodic equivalence info */

  for (rank = 0; rank < n_slices; rank++) {

    fvm_lnum_t _slice_id;

    for (_slice_id = slice_shift[rank];
         _slice_id < slice_shift[rank+1];
         _slice_id++) {

      if (equiv_id[_slice_id] >= old_count) {

        const int eq_id = equiv_id[_slice_id];
        const int l = equiv->shift[eq_id];

        assert(eq_mult[eq_id] == 0);
        equiv->rank[l] = rank;
        equiv->num[l] = slice_num[_slice_id];
        equiv->tr_id[l] = 0;
        eq_mult[eq_id] = 1;

      }
    }
  }

  for (i = 0; i < pe->count; i++) {

    int k, l;
    const fvm_lnum_t _slice_id = pe->slice_id[i];
    const int eq_id = equiv_id[_slice_id];

    for (k = pe->shift[i]; k < pe->shift[i+1]; k++) {

      l = equiv->shift[eq_id] + eq_mult[eq_id];

      equiv->rank[l] = pe->rank[k];
      equiv->num[l] = pe->num[k];

      equiv->tr_id[l] = pe->tr_id[i] + 1;

      eq_mult[eq_id] += 1;

    }

  }

  BFT_FREE(eq_mult);

  /* Free temporary periodic equivalence structure elements */

  BFT_FREE(pe->slice_id);
  BFT_FREE(pe->tr_id);
  BFT_FREE(pe->shift);
  BFT_FREE(pe->rank);
  BFT_FREE(pe->num);
}

/*----------------------------------------------------------------------------
 * Creation of a list of interfaces between entities of a same type.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   this_interface_set  <-> pointer to structure that should be updated
 *   n_ent               <-- local number of entities
 *   global_num          <-- global number (id) associated with each entity
 *   periodicity         <-- periodicity information (NULL if none)
 *   n_periodic_lists    <-- number of periodic lists (may be local)
 *   periodicity_num     <-- periodicity number (1 to n) associated with
 *                           each periodic list (primary periodicities only)
 *   n_periodic_couples  <-- number of periodic couples associated with
 *                           each periodic list
 *   periodic_couples    <-- array indicating periodic couples (using
 *                           global numberings) for each list
 *   comm                <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_add_global_equiv_periodic(fvm_interface_set_t      *this_interface_set,
                           fvm_lnum_t                n_ent,
                           fvm_gnum_t                global_num[],
                           const fvm_periodicity_t  *periodicity,
                           int                       n_periodic_lists,
                           const int                 periodicity_num[],
                           const fvm_lnum_t          n_periodic_couples[],
                           const fvm_gnum_t   *const periodic_couples[],
                           MPI_Comm                  comm)
{
  fvm_gnum_t  global_max;
  fvm_lnum_t  i, j, n_ent_recv, n_ent_send;
  size_t      slice_size;
  int         tr_index_size, size, rank;

  _per_slice_equiv_t  e;
  _per_slice_period_t pe;

  int         n_counts = 2;
  fvm_lnum_t  n_slice_couples = 0;
  int         *send_count = NULL, *recv_count = NULL;
  int         *send_shift = NULL, *recv_shift = NULL;
  fvm_lnum_t  *equiv_id = NULL, *couple_equiv_id = NULL;
  fvm_lnum_t  *equiv_send = NULL, *equiv_recv = NULL;
  fvm_lnum_t  *recv_num = NULL, *send_num = NULL;
  fvm_gnum_t  *recv_global_num = NULL;
  fvm_gnum_t  *slice_couples = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  tr_index_size = fvm_periodicity_get_n_transforms(periodicity) + 2;

  /* Get temporary maximum global number value */

  global_max = _global_num_max(n_ent,
                               global_num,
                               comm);

  /* slice_size = ceil(global_max/size) */

  slice_size = global_max / size;
  if (global_max % size > 0)
    slice_size += 1;

  assert(sizeof(fvm_gnum_t) >= sizeof(fvm_lnum_t));

  BFT_MALLOC(send_count, size*n_counts, int);
  BFT_MALLOC(recv_count, size*n_counts, int);

  BFT_MALLOC(send_shift, size + 1, int);
  BFT_MALLOC(recv_shift, size + 1, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size*n_counts; rank++)
    send_count[rank] = 0;

  for (i = 0; i < n_ent; i++)
    send_count[((global_num[i] - 1) / slice_size)] += 1;

  /* Add periodic send count */

  _count_periodic_couple_exchange(slice_size,
                                  periodicity,
                                  n_periodic_lists,
                                  n_periodic_couples,
                                  periodic_couples,
                                  send_count + size,
                                  comm);

  _exchange_counts(n_counts,
                   size,
                   send_count,
                   recv_count,
                   comm);

  /* build shift arrays for MPI_Alltoallv */

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the global_num[] array to the correct processors */

  n_ent_recv = recv_shift[size];

  BFT_MALLOC(recv_global_num, n_ent_recv, fvm_gnum_t);
  BFT_MALLOC(recv_num, n_ent_recv, fvm_lnum_t);

  MPI_Alltoallv(global_num, send_count, send_shift, FVM_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, FVM_MPI_GNUM, comm);

  /* We also need to send the corresponding local numbers */

  n_ent_send = send_shift[size];

  BFT_MALLOC(send_num, n_ent_send, fvm_lnum_t);

  for (i = 0; i < n_ent_send; i++)
    send_num[i] = i+1;

  MPI_Alltoallv(send_num, send_count, send_shift, FVM_MPI_LNUM,
                recv_num, recv_count, recv_shift, FVM_MPI_LNUM, comm);

  BFT_FREE(send_num);

  /* Exchange periodicity information */

  _exchange_periodic_couples(slice_size,
                             send_count + size,
                             recv_count + size,
                             periodicity,
                             n_periodic_lists,
                             periodicity_num,
                             n_periodic_couples,
                             periodic_couples,
                             &n_slice_couples,
                             &slice_couples,
                             comm);

  /* no grouped parallel and periodic exchanges possible from here,
     so only the first block of "size" elements is necessary */

  BFT_REALLOC(send_count, size, int);
  BFT_REALLOC(recv_count, size, int);

  /* Combine periodic couples if necessary */

  if (fvm_periodicity_get_n_levels(periodicity) > 1)
    _combine_periodic_couples(slice_size,
                              periodicity,
                              &n_slice_couples,
                              &slice_couples,
                              comm);

  /* Build purely parallel equivalence data first */

  if (n_ent_recv > 0)
    BFT_MALLOC(equiv_id, n_ent_recv, fvm_lnum_t);

  e = _slice_global_num_to_equiv(size,
                                 n_ent_recv,
                                 recv_shift,
                                 recv_global_num,
                                 recv_num,
                                 equiv_id);

  /* Now combine periodic and parallel equivalences */

  pe = _exchange_periodic_equiv(slice_size,
                                recv_shift,
                                recv_global_num,
                                recv_num,
                                equiv_id,
                                &e,
                                periodicity,
                                n_slice_couples,
                                slice_couples,
                                send_count,
                                recv_count,
                                send_shift,
                                comm);

  BFT_FREE(recv_global_num);

  _merge_periodic_equiv(size,
                        recv_shift,
                        recv_num,
                        equiv_id,
                        &e,
                        &pe);

  /* Free all arrays not needed anymore */

  BFT_FREE(recv_num);

  BFT_FREE(couple_equiv_id);
  BFT_FREE(slice_couples);
  n_slice_couples = 0;

  /* Now that equivalences are marked, count for each rank; for each
     equivalence, we will need to send the local entity number,
     the number of equivalent entities (e.multiple[...] - 1),
     and the corresponding entity numbers, ranks, and transform ids,
     for a total number of  1 + 1 + 3*(e.multiple[...] - 1)
     = 2 + 3*(e.multiple[...]) values */

  for (rank = 0; rank < size; rank++) {
    send_count[rank] = 0;
    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {
      if (equiv_id[i] > -1) {
        const int e_id = equiv_id[i];
        const int e_multiple = e.shift[e_id+1] - e.shift[e_id];
        send_count[rank] += 2 + (3*(e_multiple - 1));
      }
    }
  }

  for (rank = 0; rank < size; rank++)
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];

  /* Now prepare new send buffer */

  n_ent_send = send_shift[size];
  BFT_MALLOC(equiv_send, n_ent_send, fvm_lnum_t);

  for (rank = 0; rank < size; rank++) {

    send_count[rank] = 0; /* reset, will be re-incremented */

    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

      if (equiv_id[i] > -1) {

        fvm_lnum_t *equiv_send_p
          = equiv_send + send_shift[rank] + send_count[rank];
        fvm_lnum_t  k = 2;

        const int         e_id = equiv_id[i];
        const int         e_multiple = e.shift[e_id+1] - e.shift[e_id];
        const int        *rank_p = e.rank + e.shift[e_id];
        const int        *tr_id_p = e.tr_id + e.shift[e_id];
        const fvm_lnum_t *num_p = e.num  + e.shift[e_id];

        send_count[rank] += 2 + (3*(e_multiple - 1));

        for (j = 0; j < e_multiple; j++) {

          if (rank_p[j] == rank && tr_id_p[j] == 0) {
            equiv_send_p[0] = num_p[j];
            equiv_send_p[1] = e_multiple - 1;
          }
          else {
            equiv_send_p[k++] = num_p[j];
            equiv_send_p[k++] = tr_id_p[j];
            equiv_send_p[k++] = rank_p[j];
          }
        }

      }
    }

    assert(send_count[rank] == send_shift[rank+1] - send_shift[rank]);

  }

  /* Free temporary (slice) equivalence info */

  e.count = 0;
  BFT_FREE(e.shift);
  BFT_FREE(e.rank);
  BFT_FREE(e.tr_id);
  BFT_FREE(e.num);

  BFT_FREE(equiv_id);

  /* Send prepared slice data to destination rank */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  n_ent_recv = recv_shift[size];

  BFT_MALLOC(equiv_recv, n_ent_recv, fvm_lnum_t);

  MPI_Alltoallv(equiv_send, send_count, send_shift, FVM_MPI_LNUM,
                equiv_recv, recv_count, recv_shift, FVM_MPI_LNUM, comm);

  /* At this stage, MPI operations are finished; we may release
     the corresponding counts and indexes */

  BFT_FREE(equiv_send);

  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);

  /* Add interface */

  _interfaces_from_flat_equiv(this_interface_set,
                              tr_index_size,
                              n_ent_recv,
                              equiv_recv);

  BFT_FREE(equiv_recv);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Prepare periodic couple info in single process mode.
 *
 * Note that the array pointed to by slice_couples is allocated here,
 * and must be freed by the calling code.
 *
 * parameters:
 *   periodicity         <-- periodicity information (NULL if none)
 *   n_periodic_lists    <-- number of periodic lists (may be local)
 *   periodicity_num     <-- periodicity number (1 to n) associated with
 *                           each periodic list (primary periodicities only)
 *   n_periodic_couples  <-- number of periodic couples associated with
 *                           each periodic list
 *   periodic_couples    <-- array indicating periodic couples (using
 *                           global numberings) for each list
 *   n_couples           --> number of couples
 *   couples             --> couple information: for each couple,
 *                           {global number of local element,
 *                            global number of periodic element,
 *                            transform id}
 *----------------------------------------------------------------------------*/

static void
_define_periodic_couples_sp(const fvm_periodicity_t  *periodicity,
                            int                       n_periodic_lists,
                            const int                 periodicity_num[],
                            const fvm_lnum_t          n_periodic_couples[],
                            const fvm_gnum_t   *const periodic_couples[],
                            fvm_lnum_t               *n_couples,
                            fvm_gnum_t              **couples)
{
  int         list_id;

  fvm_lnum_t  count = 0;
  fvm_lnum_t   _n_couples = 0;
  fvm_gnum_t  *_couples = NULL;

  /* Initialization */

  *n_couples = 0;
  *couples = NULL;

  for (list_id = 0, _n_couples = 0; list_id < n_periodic_lists; list_id++)
    _n_couples += n_periodic_couples[list_id] * 2;

  BFT_MALLOC(_couples, _n_couples*3, fvm_gnum_t);

  /* Prepare lists */

  for (list_id = 0; list_id < n_periodic_lists; list_id++) {

    fvm_lnum_t couple_id;
    fvm_gnum_t num_1, num_2;

    const int external_num = periodicity_num[list_id];
    const int direct_id = fvm_periodicity_get_transform_id(periodicity,
                                                           external_num,
                                                           1);
    const int reverse_id = fvm_periodicity_get_transform_id(periodicity,
                                                            external_num,
                                                            -1);

    const fvm_lnum_t  _n_periodic_couples = n_periodic_couples[list_id];
    const fvm_gnum_t  *_periodic_couples = periodic_couples[list_id];

    assert(direct_id >= 0 && reverse_id >= 0);

    for (couple_id = 0; couple_id < _n_periodic_couples; couple_id++) {

      num_1 = _periodic_couples[couple_id*2];
      num_2 = _periodic_couples[couple_id*2 + 1];

      _couples[count] = num_1;
      _couples[count+1] = num_2;
      _couples[count+2] = direct_id;

      _couples[count+3] = num_2;
      _couples[count+4] = num_1;
      _couples[count+5] = reverse_id;

      count += 6;

    }

  }

  /* Sort periodic couples by local correspondant, remove duplicates */

  _sort_periodic_couples(&_n_couples,
                         &_couples);

  /* Set return values */

  *n_couples = _n_couples;
  *couples = _couples;
}

/*----------------------------------------------------------------------------
 * Build eventual couples belonging to combined periodicities in single
 * process mode.
 *
 * parameters:
 *   periodicity  <-- periodicity information (NULL if none)
 *   n_couples    <-> number of couples in current slice
 *   couples      <-> couple information for this rank: for each couple,
 *                    {global number of local element,
 *                     global number of periodic element,
 *                     transform id}
 *----------------------------------------------------------------------------*/

static void
_combine_periodic_couples_sp(const fvm_periodicity_t   *periodicity,
                             fvm_lnum_t                *n_couples,
                             fvm_gnum_t               **couples)
{
  int  i, n_tr, level;

  fvm_lnum_t  start_id, end_id, j, k;

  int  add_count = 0;
  int  *tr_reverse_id = NULL;
  fvm_gnum_t *_add_couples = NULL;

  fvm_lnum_t  _n_couples = *n_couples;
  fvm_gnum_t  *_couples = *couples;

  assert (periodicity != NULL);

  /* Build periodicity related arrays for quick access */

  n_tr = fvm_periodicity_get_n_transforms(periodicity);

  BFT_MALLOC(tr_reverse_id, n_tr, int);

  for (i = 0; i < n_tr; i++)
    tr_reverse_id[i] = fvm_periodicity_get_reverse_id(periodicity, i);

  /* Loop on combination levels */

  for (level = 1;
       level < fvm_periodicity_get_n_levels(periodicity);
       level++) {

    int  n_rows = 0;
    int  *tr_combine = NULL;

    _transform_combine_info(periodicity,
                            level,
                            &n_rows,
                            &tr_combine);

    /* Count values to add */

    add_count = 0;

    start_id = 0;
    end_id = 1;

    while (end_id < _n_couples) {

      if (_couples[start_id*3] == _couples[end_id*3]) {

        end_id++;
        while (end_id < _n_couples) {
          if (_couples[end_id*3] != _couples[start_id*3])
            break;
          end_id++;
        }

        for (j = start_id; j < end_id; j++) {
          for (k = j+1; k < end_id; k++) {

            int tr_1 = tr_reverse_id[_couples[j*3 + 2]];
            int tr_2 = _couples[k*3 + 2];

            if (tr_combine[tr_1 *n_rows + tr_2] > -1)
              add_count += 6;

          }
        }

      }

      start_id = end_id;
      end_id += 1;
    }

    /* Nothing to do for this combination level if add_count = 0 */

    if (add_count == 0) {
      BFT_FREE(tr_combine);
      continue;
    }

    BFT_REALLOC(_couples, _n_couples*3 + add_count, fvm_gnum_t);

    _add_couples = _couples + _n_couples*3;

    /* Now add combined couples */

    start_id = 0;
    end_id = 1;

    while (end_id < _n_couples) {

      if (_couples[start_id*3] == _couples[end_id*3]) {

        end_id++;
        while (end_id < _n_couples) {
          if (_couples[end_id*3] != _couples[start_id*3])
            break;
          end_id++;
        }

        /* Loop on couple combinations */

        for (j = start_id; j < end_id; j++) {
          for (k = j+1; k < end_id; k++) {

            fvm_gnum_t num_1 = _couples[j*3 + 1];
            fvm_gnum_t num_2 = _couples[k*3 + 1];
            int tr_1 = tr_reverse_id[_couples[j*3 + 2]];
            int tr_2 = _couples[k*3 + 2];
            int tr_c = tr_combine[tr_1 *n_rows + tr_2];

            if (tr_c > -1) {

              _add_couples[0] = num_1;
              _add_couples[1] = num_2;
              _add_couples[2] = tr_c;

              _add_couples[3] = num_2;
              _add_couples[4] = num_1;
              _add_couples[5] = tr_reverse_id[tr_c];

              _add_couples += 6;

            }

          }
        } /* End of loop on couple combinations */

      }

      start_id = end_id;
      end_id += 1;
    }

    BFT_FREE(tr_combine);

    /* Finally, merge additional couples received for slice with
       existing periodicity information */

    assert(add_count % 3 == 0);

    _n_couples += add_count / 3;

    /* Sort and remove duplicates to update periodicity info */

    _sort_periodic_couples(&_n_couples, &_couples);

    *n_couples = _n_couples;
    *couples = _couples;

  }

  BFT_FREE(tr_reverse_id);
}

/*----------------------------------------------------------------------------
 * Creation of a list of interfaces between entities of a same type.
 *
 * This simplified algorithm is intended for single-process mode.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   this_interface_set  <-> pointer to structure that should be updated
 *   global_num          <-- global number (id) associated with each entity
 *   periodicity         <-- periodicity information (NULL if none)
 *   n_periodic_lists    <-- number of periodic lists (may be local)
 *   periodicity_num     <-- periodicity number (1 to n) associated with
 *                           each periodic list (primary periodicities only)
 *   n_periodic_couples  <-- number of periodic couples associated with
 *                           each periodic list
 *   periodic_couples    <-- array indicating periodic couples (using
 *                           global numberings) for each list
 *----------------------------------------------------------------------------*/

static void
_add_global_equiv_periodic_sp(fvm_interface_set_t      *this_interface_set,
                              const fvm_periodicity_t  *periodicity,
                              int                       n_periodic_lists,
                              const int                 periodicity_num[],
                              const fvm_lnum_t          n_periodic_couples[],
                              const fvm_gnum_t   *const periodic_couples[])
{
  fvm_lnum_t  i, couple_id;
  fvm_lnum_t  n_couples = 0;
  fvm_lnum_t  *n_ent_tr = NULL;
  fvm_gnum_t  *couples = NULL;

  fvm_interface_t  *_interface = NULL;

  assert(sizeof(fvm_gnum_t) >= sizeof(fvm_lnum_t));

  /* Organize periodicity information */

  _define_periodic_couples_sp(periodicity,
                              n_periodic_lists,
                              periodicity_num,
                              n_periodic_couples,
                              periodic_couples,
                              &n_couples,
                              &couples);

  /* Combine periodic couples if necessary */

  if (fvm_periodicity_get_n_levels(periodicity) > 1)
    _combine_periodic_couples_sp(periodicity,
                                 &n_couples,
                                 &couples);

  /* Add interface to set */

  this_interface_set->size += 1;

  BFT_REALLOC(this_interface_set->interfaces,
              this_interface_set->size,
              fvm_interface_t *);

  _interface = _fvm_interface_create();

  this_interface_set->interfaces[this_interface_set->size - 1] = _interface;

  /* Build interface */

  _interface->rank = 0;
  _interface->size = n_couples;

  _interface->tr_index_size
    = fvm_periodicity_get_n_transforms(periodicity) + 2;

  BFT_MALLOC(_interface->tr_index, _interface->tr_index_size, fvm_lnum_t);
  BFT_MALLOC(_interface->local_num, _interface->size, fvm_lnum_t);
  BFT_MALLOC(_interface->distant_num, _interface->size, fvm_lnum_t);

  /* Count couples for each transform */

  BFT_MALLOC(n_ent_tr, _interface->tr_index_size - 2, fvm_lnum_t);

  for (i = 0; i < _interface->tr_index_size - 2; i++)
    n_ent_tr[i] = 0;

  for (couple_id = 0; couple_id < n_couples; couple_id++)
    n_ent_tr[couples[couple_id*3 + 2]] += 1;

  /* Build index */

  _interface->tr_index[0] = 0;
  _interface->tr_index[1] = 0;

  for (i = 2; i < _interface->tr_index_size; i++) {
    _interface->tr_index[i] = _interface->tr_index[i-1] + n_ent_tr[i-2];
    n_ent_tr[i-2] = 0;
  }

  /* Build local and distant correspondants */

  for (couple_id = 0; couple_id < n_couples; couple_id++) {

    const int tr_id = couples[couple_id*3 + 2];
    const fvm_lnum_t j = _interface->tr_index[tr_id + 1] + n_ent_tr[tr_id];

    _interface->local_num[j] = couples[couple_id*3];
    _interface->distant_num[j] = couples[couple_id*3 + 1];

    n_ent_tr[tr_id] += 1;

  }

  /* Free temporary arrays */

  BFT_FREE(n_ent_tr);
  BFT_FREE(couples);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Apply renumbering of entities referenced by an interface set.
 *
 * For any given entity i, a negative old_to_new[i] value means that that
 * entity does not appear anymore in the new numbering.
 *
 * parameters:
 *   this_interface_set <-> pointer to interface set structure
 *   old_to_new         <-- renumbering array (0 to n-1 numbering)
 *   comm               <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_set_renumber(fvm_interface_set_t  *this_interface_set,
              const fvm_lnum_t      old_to_new[],
              MPI_Comm              comm)
{
  int i;
  int local_rank = 0;
  int n_interfaces = 0;
  int n_tr = 0;
  int request_count = 0;
  int *rev_tr_id = NULL;
  fvm_lnum_t *if_index = NULL;
  fvm_lnum_t *tr_send_buf = NULL;
  MPI_Request  *request = NULL;
  MPI_Status  *status  = NULL;

  assert(this_interface_set != NULL);
  assert(old_to_new != NULL);

  n_interfaces = this_interface_set->size;

  /* Additionnal arrays for differently ordered sections and their
     arrays are needed in case of periodic transforms */

  if (this_interface_set->periodicity != NULL) {

    const fvm_periodicity_t *p = this_interface_set->periodicity;
    n_tr = fvm_periodicity_get_n_transforms(p);

    BFT_MALLOC(rev_tr_id, n_tr + 1, int);
    rev_tr_id[0] = 0;
    for (i = 0; i < n_tr; i++)
      rev_tr_id[i+1] = fvm_periodicity_get_reverse_id(p, i) + 1;

    BFT_MALLOC(if_index, n_interfaces + 1, fvm_lnum_t);
    if_index[0] = 0;
    for (i = 0; i < n_interfaces; i++) {
      const fvm_interface_t *_if = this_interface_set->interfaces[i];
      if_index[i+1] = if_index[i] + _if->size;
    }

    BFT_MALLOC(tr_send_buf, if_index[n_interfaces], fvm_lnum_t);
  }

  MPI_Comm_rank(comm, &local_rank);

  /* Update numbering */

  for (i = 0; i < n_interfaces; i++) {

    fvm_lnum_t j;
    fvm_interface_t *_if = this_interface_set->interfaces[i];

    fvm_lnum_t *loc_num = _if->local_num;
    fvm_lnum_t *dist_num = _if->distant_num;

    for (j = 0; j < _if->size; j++)
      loc_num[j] = old_to_new[loc_num[j] - 1] + 1;

    if (local_rank == _if->rank) {
      for (j = 0; j < _if->size; j++)
        dist_num[j] = old_to_new[dist_num[j] - 1] + 1;
    }
  }

  /* Exchange distant face numbers */

  BFT_MALLOC(request, n_interfaces*2, MPI_Request);
  BFT_MALLOC(status, n_interfaces*2, MPI_Status);

  for (i = 0; i < n_interfaces; i++) {

    fvm_interface_t *_if = this_interface_set->interfaces[i];
    const int distant_rank = _if->rank;

    if (distant_rank != local_rank)
      MPI_Irecv(_if->distant_num,
                _if->size,
                FVM_MPI_LNUM,
                distant_rank,
                distant_rank,
                comm,
                &(request[request_count++]));
  }

  for (i = 0; i < n_interfaces; i++) {

    fvm_interface_t *_if = this_interface_set->interfaces[i];
    fvm_lnum_t *send_buf = _if->local_num;
    const int distant_rank = _if->rank;

    if (distant_rank != local_rank) {

      if (rev_tr_id != NULL) {
        int j;
        fvm_lnum_t k, l;
        send_buf = tr_send_buf + if_index[i];
        assert(_if->tr_index_size == n_tr + 2);
        for (j = 0, l = 0; j < _if->tr_index_size - 1; j++) {
          for (k = _if->tr_index[rev_tr_id[j]];
               k < _if->tr_index[rev_tr_id[j] + 1];
               k++)
            send_buf[l++] = _if->local_num[k];
        }
      }

      MPI_Isend(send_buf,
                _if->size,
                FVM_MPI_LNUM,
                distant_rank,
                local_rank,
                comm,
                &(request[request_count++]));
    }
  }

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  if (this_interface_set->periodicity != NULL) {
    BFT_FREE(rev_tr_id);
    BFT_FREE(if_index);
    BFT_FREE(tr_send_buf);
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Apply renumbering of entities referenced by an interface set.
 *
 * This simplified algorithm is intended for single-process mode.
 *
 * For any given entity i, a negative old_to_new[i] value means that that
 * entity does not appear anymore in the new numbering.
 *
 * parameters:
 *   this_interface_set <-> pointer to interface set structure
 *   old_to_new         <-- renumbering array (0 to n-1 numbering)
 *   comm               <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_set_renumber_sp(fvm_interface_set_t  *this_interface_set,
                 const fvm_lnum_t      old_to_new[])
{
  int i;

  assert(this_interface_set != NULL);
  assert(old_to_new != NULL);

  /* Update numbering */

  for (i = 0; i < this_interface_set->size; i++) {

    fvm_lnum_t j;
    fvm_interface_t *_if = this_interface_set->interfaces[i];

    fvm_lnum_t *loc_num = _if->local_num;
    fvm_lnum_t *dist_num = _if->distant_num;

    for (j = 0; j < _if->size; j++) {
      loc_num[j] = old_to_new[loc_num[j] - 1] + 1;
      dist_num[j] = old_to_new[dist_num[j] - 1] + 1;
    }
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return process rank associated with an interface's distant entities.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   process rank associated with the interface's distant entities
 *----------------------------------------------------------------------------*/

int
fvm_interface_rank(const fvm_interface_t  *this_interface)
{
  int retval = -1;

  if (this_interface != NULL)
    retval = this_interface->rank;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return number of local and distant entities defining an interface.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   number of local and distant entities defining the interface
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_interface_size(const fvm_interface_t  *this_interface)
{
  fvm_lnum_t retval = 0;

  if (this_interface != NULL)
    retval = this_interface->size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to array of local entity numbers defining an interface.
 *
 * The size of the array may be obtained by fvm_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   pointer to array of local entity numbers defining the interface
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_interface_get_local_num(const fvm_interface_t  *this_interface)
{
  const fvm_lnum_t *retval = NULL;

  if (this_interface != NULL)
    retval = this_interface->local_num;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to array of distant entity numbers defining an interface.
 *
 * The size of the array may be obtained by fvm_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   pointer to array of local entity numbers defining the interface
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_interface_get_distant_num(const fvm_interface_t  *this_interface)
{
  const fvm_lnum_t *retval = NULL;

  if (this_interface != NULL)
    retval = this_interface->distant_num;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return size of index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * Its size should thus be equal to 1 + number of periodic transforms + 1,
 * In absence of periodicity, it may be 0, as the index is not needed.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   transform index size for the interface
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_interface_get_tr_index_size(const fvm_interface_t  *this_interface)
{
  fvm_lnum_t retval = 0;

  if (this_interface != NULL)
    retval = this_interface->tr_index_size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * In absence of periodicity, it may be NULL, as it is not needed.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   pointer to transform index for the interface
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_interface_get_tr_index(const fvm_interface_t  *this_interface)
{
  const fvm_lnum_t *tr_index = 0;

  if (this_interface != NULL)
    tr_index = this_interface->tr_index;

  return tr_index;
}

/*----------------------------------------------------------------------------
 * Creation of a list of interfaces between entities of a same type.
 *
 * These interfaces may be used to identify equivalent vertices or faces using
 * domain splitting, as well as periodic entities (on the same or on
 * distant ranks).
 *
 * Note that periodicity information will be completed and made consistent
 * based on the input, so that if a periodic couple is defined on a given rank,
 * the reverse couple wil be defined, whether it is also defined on the same
 * or a different rank.
 *
 * In addition, multiple periodicity interfaces will be built automatically
 * if the periodicity structure provides for composed periodicities, so they
 * need not be defined prior to this function.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   n_ent                <-- number of local entities considered
 *                            (size of parent_entity_number[])
 *   parent_entity_number <-- pointer to list of selected entities local
 *                            numbers (1 to n), or NULL if all first n_ent
 *                            entities are used
 *   global_number        <-- pointer to list of global (i.e. domain splitting
 *                            independent) entity numbers
 *   periodicity          <-- periodicity information (NULL if none)
 *   n_periodic_lists     <-- number of periodic lists (may be local)
 *   periodicity_num      <-- periodicity number (1 to n) associated with
 *                            each periodic list (primary periodicities only)
 *   n_periodic_couples   <-- number of periodic couples associated with
 *                            each periodic list
 *   periodic_couples     <-- array indicating periodic couples (interlaced,
 *                            using global numberings) for each list
 *
 * returns:
 *  pointer to list of interfaces (possibly NULL in serial mode)
 *----------------------------------------------------------------------------*/

fvm_interface_set_t *
fvm_interface_set_create(fvm_lnum_t                n_ent,
                         const fvm_lnum_t          parent_entity_number[],
                         const fvm_gnum_t          global_number[],
                         const fvm_periodicity_t  *periodicity,
                         int                       n_periodic_lists,
                         const int                 periodicity_num[],
                         const fvm_lnum_t          n_periodic_couples[],
                         const fvm_gnum_t   *const periodic_couples[])
{
  fvm_interface_set_t  *this_interface_set;

  /* Initial checks */

  if (   fvm_order_local_test(parent_entity_number,
                              global_number,
                              n_ent) == false
      && global_number != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Trying to build an interface with unordered elements\n"
                "not currently handled."));

  if (   (fvm_parall_get_size() < 2)
      && (periodicity == NULL || n_periodic_lists == 0))
    return NULL;

  /* Create structure */

  BFT_MALLOC(this_interface_set, 1, fvm_interface_set_t);
  this_interface_set->size = 0;
  this_interface_set->interfaces = NULL;
  this_interface_set->periodicity = periodicity;

#if defined(HAVE_MPI)

  if (fvm_parall_get_size() > 1) {

    size_t  i;
    fvm_gnum_t  *global_num = NULL;

    if (n_ent > 0) {

      BFT_MALLOC(global_num, n_ent, fvm_gnum_t);

      /* Assign initial global numbers */

      if (parent_entity_number != NULL) {
        for (i = 0 ; i < (size_t)n_ent ; i++)
          global_num[i] = global_number[parent_entity_number[i]-1];
      }
      else {
        for (i = 0 ; i < (size_t)n_ent ; i++)
          global_num[i] = global_number[i];
      }

    }

    /* Build interfaces */

    if (fvm_parall_get_size() > 1 && periodicity == NULL)

      _fvm_interface_add_global_equiv(this_interface_set,
                                      n_ent,
                                      global_num,
                                      fvm_parall_get_mpi_comm());

    else if (fvm_parall_get_size() > 1)

      _add_global_equiv_periodic(this_interface_set,
                                 n_ent,
                                 global_num,
                                 periodicity,
                                 n_periodic_lists,
                                 periodicity_num,
                                 n_periodic_couples,
                                 periodic_couples,
                                 fvm_parall_get_mpi_comm());

    if (fvm_parall_get_size() > 1)
      BFT_FREE(global_num);

 }

#endif /* defined(HAVE_MPI) */

  if (   (fvm_parall_get_size() == 1)
      && (periodicity != NULL && n_periodic_lists > 0)) {

    _add_global_equiv_periodic_sp(this_interface_set,
                                  periodicity,
                                  n_periodic_lists,
                                  periodicity_num,
                                  n_periodic_couples,
                                  periodic_couples);


  }

  return this_interface_set;
}

/*----------------------------------------------------------------------------
 * Destruction of an interface list.
 *
 * parameters:
 *   this_interface_set <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_interface_set_t *
fvm_interface_set_destroy(fvm_interface_set_t  *this_interface_set)
{
  int i;

  if (this_interface_set != NULL) {
    for (i = 0; i < this_interface_set->size; i++) {
      _fvm_interface_destroy(this_interface_set->interfaces[i]);
    }
    BFT_FREE(this_interface_set->interfaces);
    BFT_FREE(this_interface_set);
  }

  return this_interface_set;
}

/*----------------------------------------------------------------------------
 * Return number of interfaces associated with an interface set.
 *
 * parameters:
 *   this_interface_set <-- pointer to interface set structure
 *
 * returns:
 *   number of interfaces in set
 *----------------------------------------------------------------------------*/

int
fvm_interface_set_size(const fvm_interface_set_t  *this_interface_set)
{
  int retval = 0;

  if (this_interface_set != NULL)
    retval = this_interface_set->size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to a given interface in an interface set.
 *
 * parameters:
 *   this_interface_set <-- pointer to interface set structure
 *   interface_id       <-- index of interface in set (0 to n-1)
 *
 * returns:
 *   pointer to interface structure
 *----------------------------------------------------------------------------*/

const fvm_interface_t *
fvm_interface_set_get(const fvm_interface_set_t  *this_interface_set,
                      int                         interface_id)
{
  const fvm_interface_t  *retval = NULL;

  if (this_interface_set != NULL) {
    if (interface_id > -1 && interface_id < this_interface_set->size)
      retval = this_interface_set->interfaces[interface_id];
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to the periocicity structure associated of an interface set.
 *
 * parameters:
 *   this_interface_set <-- pointer to interface set structure
 *
 * returns:
 *   pointer to periodicity structure, or NULL
 *----------------------------------------------------------------------------*/

const fvm_periodicity_t *
fvm_interface_set_periodicity(const fvm_interface_set_t  *this_interface_set)
{
  const fvm_periodicity_t  *retval = NULL;

  if (this_interface_set != NULL)
    retval = this_interface_set->periodicity;

  return retval;
}

/*----------------------------------------------------------------------------
 * Apply renumbering of entities referenced by an interface set.
 *
 * For any given entity i, a negative old_to_new[i] value means that that
 * entity does not appear anymore in the new numbering.
 *
 * parameters:
 *   this_interface_set <-> pointer to interface set structure
 *   old_to_new         <-- renumbering array (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

void
fvm_interface_set_renumber(fvm_interface_set_t  *this_interface_set,
                           const fvm_lnum_t      old_to_new[])
{
  int i;
  int n_interfaces = 0;

  assert(this_interface_set != NULL);
  assert(old_to_new != NULL);

#if defined(HAVE_MPI)

  if (fvm_parall_get_size() > 1)
    _set_renumber(this_interface_set,
                  old_to_new,
                  fvm_parall_get_mpi_comm());

#endif /* defined(HAVE_MPI) */

  if (fvm_parall_get_size() == 1)
    _set_renumber_sp(this_interface_set,
                     old_to_new);

  /* Remove references to entities not appearing anymore */

  for (i = 0; i < this_interface_set->size; i++) {

    fvm_lnum_t j, k;
    fvm_interface_t *_if = this_interface_set->interfaces[i];

    fvm_lnum_t *loc_num = _if->local_num;
    fvm_lnum_t *dist_num = _if->distant_num;

    if (_if->tr_index_size == 0) {
      for (j = 0, k = 0; j < _if->size; j++) {
        if (loc_num[j] > -1 && dist_num[j] > -1) {
          loc_num[k] = loc_num[j];
          dist_num[k] = dist_num[j];
          k += 1;
        }
      }
    }
    else {
      int tr_id;
      fvm_lnum_t end_id = _if->tr_index[0];
      k = 0;
      for (tr_id = 0; tr_id < _if->tr_index_size - 1; tr_id++) {
        fvm_lnum_t start_id = end_id;
        end_id = _if->tr_index[tr_id + 1];
        for (j = start_id; j < end_id; j++) {
          if (loc_num[j] > -1 && dist_num[j] > -1) {
            loc_num[k] = loc_num[j];
            dist_num[k] = dist_num[j];
            k += 1;
          }
        }
        _if->tr_index[tr_id + 1] = k;
      }
    }

    if (k < _if->size) {
      if (k > 0) {
        _if->size = k;
        BFT_REALLOC(_if->local_num, k, fvm_lnum_t);
        BFT_REALLOC(_if->distant_num, k, fvm_lnum_t);
      }
      else {
        BFT_FREE(_if->local_num);
        BFT_FREE(_if->distant_num);
        BFT_FREE(this_interface_set->interfaces[i]);
      }
    }
  }

  for (i = 0, n_interfaces = 0; i < this_interface_set->size; i++) {
    if (this_interface_set->interfaces[i] != NULL)
      this_interface_set->interfaces[n_interfaces++]
        = this_interface_set->interfaces[i];
  }

  if (n_interfaces < this_interface_set->size) {
    BFT_REALLOC(this_interface_set->interfaces,
                n_interfaces,
                fvm_interface_t *);
    this_interface_set->size = n_interfaces;
  }
}

/*----------------------------------------------------------------------------
 * Dump printout of an interface list.
 *
 * parameters:
 *   this_interface_set <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_interface_set_dump(const fvm_interface_set_t  *this_interface_set)
{
  int i;

  if (this_interface_set == NULL) {
    bft_printf("  interface list: nil\n");
    return;
  }

  bft_printf("  interface list: %p\n"
             "  n interfaces:   %d\n",
             this_interface_set, this_interface_set->size);

  for (i = 0; i < this_interface_set->size; i++) {
    bft_printf("\n  interface %d:\n", i);
    _fvm_interface_dump(this_interface_set->interfaces[i]);
  }

  if (this_interface_set->periodicity != NULL)
    bft_printf("\n  periodicity %p:\n", this_interface_set->periodicity);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
