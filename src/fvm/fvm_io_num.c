/*============================================================================
 * Main structure for an I/O numbering scheme associated with mesh entities
 * (such as cells, faces, and vertices);
 *
 * In parallel mode, such a scheme is important so as to redistribute
 * locally numbered entities on n processes to files written by p
 * processes, with p <= n.
 *
 * Only the case where p = 1 is presently implemented, so the numbering
 * scheme is simply based on entity's global labels.
 *
 * For p > 1, it would probably be necessary to extend the numbering
 * schemes so as to account for the fact that a given entity may have
 * a main index on its main associated domain, but may be present
 * as a ghost entity with another index on neighboring domains.
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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_hilbert.h"
#include "fvm_morton.h"

#include "cs_all_to_all.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_sort_partition.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_io_num.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

/*
 * Notes:
 *
 * This structure currently only contains a global numbering array containing
 * each entity's global number (a 1 to n index). In the future, it may
 * also contain information relative to ghost zones, for I/O to file formats
 * enabling domain splitting with (with multiple groups of processes writing
 * different subsets, using ghost zones for entities on boundaries between
 * mesh parts assigned to different process groups). In such a case, the
 * main numbering would only be global as regards a process group, at least
 * for use with formats such as that of EnSight Gold.
 *
 * In some cases, a global number may appear on multiple processes (for a face
 * or vertex on a processor boundary). This means that multiple processes may
 * update the corresponding value in a gather-type operation preceding or
 * associated with I/O, in an undefined order. If the associated data is up to
 * date on each process, it should be identical, so this should not be a
 * problem. MPI-IO file writes or MPI-2 one-sided communication PUT operations
 * also authorize this, though for the latter, the MPI-2 standard indicates
 * that this is authorized if processes overlapping PUT operations should use
 * the same predefined datatype, which seems to exclude similar indexed
 * datatypes with different indexes. To avoid problems if we wish to use
 * MPI-2 one-sided communication, one relatively simple solution would be to
 * consider that for processes other than that of lowest rank in which an
 * entity appears, it appears after all entities not occuring in processes
 * of lower rank. In that case, we would have two array sizes:
 * global_num_size_tot defining the array's full size, and global_num_size
 * defining the size of the portion to use in gather type operations.
 */

struct _fvm_io_num_t {

  cs_gnum_t          global_count;    /* Global number of entities */
  cs_lnum_t          global_num_size; /* Local size of global numbering array */
  const cs_gnum_t   *global_num;      /* Global (possibly shared) entity
                                         numbers (1 to n) */
  cs_gnum_t         *_global_num;     /* Global entity numbers if owner,
                                         NULL otherwise */

};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Names of space-filling curve types */

const char  *fvm_io_num_sfc_type_name[] = {N_("Morton (in bounding box)"),
                                           N_("Morton (in bounding cube)"),
                                           N_("Hilbert (in bounding box)"),
                                           N_("Hilbert (in bounding cube)")};

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer for conversion of a double precision value in
 * range [0, 1] to the same value.
 *
 * This is a trivial function, passed through a function pointer
 * to cs_sort_partition_dest_rank_id.
 *
 * parameters:
 *   s     <-- coordinate between 0 and 1
 *   elt   -->  pointer to element
 *   input <-- pointer to optional (untyped) value or structure.
 */
/*----------------------------------------------------------------------------*/

static void
_s_to_real(double       s,
           void        *elt,
           const void  *input)
{
  CS_UNUSED(input);

  cs_real_t  *v = elt;
  *v = s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief function pointer for comparison of 2
 *
 * This function is the same type as that used by qsort_r.
 *
 * \param[in]  elt1   coordinate between 0 and 1
 * \param[in]  elt2   pointer to optional (untyped) value or structure.
 * \param[in]  input  pointer to optional (untyped) value or structure.
 *
 * \return < 0 if elt1 < elt2, 0 if elt1 == elt2, > 0 if elt1 > elt2
 */
/*----------------------------------------------------------------------------*/

static int
_s_compare(const void  *elt1,
           const void  *elt2,
           const void  *input)
{
  CS_UNUSED(input);

  int retval = 0;
  if (  *(const cs_real_t *)elt1
      < *(const cs_real_t *)elt2)
    retval = -1;
  else if (  *(const cs_real_t *)elt1
           > *(const cs_real_t *)elt2)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Use bubble sort on an expectedly short sequence of coordinates
 * to ensure lexicographical ordering.
 *
 * parameters:
 *   dim        <-- spatial dimension
 *   start_id   <-- start id in array
 *   end_id     <-- past-the-end id in array
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   order      <-> ordering array base on Morton encoding, or
 *                  lexicographical coordinate ordering for ties
 *----------------------------------------------------------------------------*/

inline static void
_reorder_coords_lexicographic(int               dim,
                              size_t            start_id,
                              size_t            end_id,
                              const cs_coord_t  coords[],
                              cs_lnum_t         order[])
{
  size_t  i;
  _Bool g_swap;

  do {

    g_swap = false;

    for (i = start_id + 1; i < end_id; i++) {

      size_t j_prev = order[i-1], j = order[i];
      _Bool l_swap = false;

      if (dim == 3) {
        if (coords[j_prev*3] < coords[j*3])
          continue;
        else if (coords[j_prev*3] > coords[j*3])
          l_swap = true;
        else if (coords[j_prev*3 + 1] < coords[j*3 + 1])
          continue;
        else if (   coords[j_prev*3 + 1] > coords[j*3 + 1]
                 || coords[j_prev*3 + 2] > coords[j*3 + 2])
          l_swap = true;
      }
      else if (dim == 2) {
        if (coords[j_prev*2] < coords[j*2 + 1])
          continue;
        else if (   coords[j_prev*2]     > coords[j*2]
                 || coords[j_prev*2 + 1] > coords[j*2 + 1])
          l_swap = true;
      }
      else { /* if (dim == 1) */
        if (coords[j_prev] > coords[j])
          l_swap = true;
      }

      if (l_swap) {
        cs_lnum_t o_save = order[i-1];
        order[i-1] = order[i];
        order[i] = o_save;
        g_swap = true;
      }
    }

  } while (g_swap);
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be different, but their order is undetermined.
 *
 * parameters:
 *   dim        <-- spatial dimension
 *   n_entities <-- number of entities considered
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   m_code     <-- Morton code associated with each entity
 *   order      <-> ordering array base on Morton encoding, or
 *                  lexicographical coordinate ordering for ties
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

static void
_check_morton_ordering(int                      dim,
                       size_t                   n_entities,
                       const cs_coord_t         coords[],
                       const fvm_morton_code_t  m_code[],
                       cs_lnum_t                order[])
{
  size_t  i_prev = 0, i = 1;

  if (n_entities == 0)
    return;

  /* Check ordering; if two entities have the same Morton codes,
     use lexicographical coordinates ordering to ensure the
     final order is deterministic. */

  for (i = 1; i < n_entities; i++) {

    size_t j_prev = order[i_prev], j = order[i];

    if (   m_code[j_prev].X[0] != m_code[j].X[0]
        || m_code[j_prev].X[1] != m_code[j].X[1]
        || m_code[j_prev].X[2] != m_code[j].X[2]) {

      /* If successive values have the same Morton code,
         order them lexicographically */
      if (i_prev < i - 1)
        _reorder_coords_lexicographic(dim, i_prev, i-1, coords, order);

      i_prev = i;
    }
  }

  if (i_prev < n_entities - 1)
    _reorder_coords_lexicographic(dim, i_prev, n_entities - 1, coords, order);
}

/*----------------------------------------------------------------------------
 * Maximum local global number associated with an I/O numbering structure.
 *
 * This function is to be used for ordered global numberings.
 *
 * parameters:
 *   this_io_num <-- pointer to partially initialized I/O numbering structure.
 *
 * returns:
 *   maximum global number associated with the I/O numbering
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_fvm_io_num_local_max(const fvm_io_num_t  *this_io_num)
{
  cs_gnum_t   local_max;

  /* Get maximum global number value */

  size_t n_ent = this_io_num->global_num_size;
  if (n_ent > 0)
    local_max = this_io_num->global_num[n_ent - 1];
  else
    local_max = 0;

  return local_max;
}

/*----------------------------------------------------------------------------
 * When sub-entities have been added to an I/O numbering, switch from
 * a numbering on the initial entities (shifted by number of sub-entities) to
 * a numbering on the final sub-entities.
 *
 * Also update whether the numbering may be shared, and discard copy if not
 * needed.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   n_sub_entities <-- optional number of sub-entities per initial entity,
 *                      or NULL if unused
 *   may_be_shared  <-- indicate if structure may be shared at this stage
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_order_finalize(fvm_io_num_t     *this_io_num,
                           const cs_lnum_t   n_sub_entities[],
                           _Bool             may_be_shared)
{
  if (n_sub_entities != NULL) {

    cs_lnum_t i, j, k;
    cs_gnum_t *_global_num;

    for (i = 0, j = 0; i < this_io_num->global_num_size; i++)
      j += n_sub_entities[i];

    BFT_MALLOC(_global_num, j, cs_gnum_t);

    for (i = 0, j = 0; i < this_io_num->global_num_size; i++) {
      for (k = 0; k < n_sub_entities[i]; j++, k++)
        _global_num[j] = this_io_num->_global_num[i] - n_sub_entities[i] + k + 1;
    }

    BFT_FREE(this_io_num->_global_num);
    this_io_num->_global_num = _global_num;

    if (this_io_num->global_num_size != (cs_lnum_t)j) {
      this_io_num->global_num_size = j;
      may_be_shared = false;
    }

    if (may_be_shared == false)
      this_io_num->global_num = this_io_num->_global_num;
  }

  /* If numbering was initially shared, check if it was changed or if it
     may remain shared (in which case the copy may be discarded) */

  if (may_be_shared == true) {
    cs_lnum_t i;
    for (i = 0; i < this_io_num->global_num_size; i++)
      if (this_io_num->_global_num[i] != this_io_num->global_num[i])
        break;
    if (i < this_io_num->global_num_size)
      this_io_num->global_num = this_io_num->_global_num;
    else
      BFT_FREE(this_io_num->_global_num);
  }
}

/*----------------------------------------------------------------------------
 * Local ordering associated with an I/O numbering structure.
 *
 * The structure should contain an initial ordering, which should
 * be sorted, but need not be contiguous. On output, the numbering
 * will be contiguous.
 *
 * As an option, a number of sub-entities per initial entity may be
 * given, in which case sub-entities of a same entity will have contiguous
 * numbers in the final ordering.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   n_sub_entities <-- optional number of sub-entities per initial entity,
 *                      or NULL if unused
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_local_order(fvm_io_num_t     *this_io_num,
                        const cs_lnum_t   n_sub_entities[])
{
  cs_gnum_t   num_prev, num_cur;

  _Bool       may_be_shared = false;

  cs_gnum_t   current_gnum = 0;

  /* If numbering is shared, we will check later it was changed or if
     it can remain shared (in which case the copy may be discarded) */

  if (this_io_num->global_num != this_io_num->_global_num)
    may_be_shared = true;
  else
    may_be_shared = false;

  size_t n_ent = this_io_num->global_num_size;

  if (n_ent > 0) {

    cs_lnum_t *b_order;
    cs_gnum_t *b_gnum = this_io_num->_global_num;
    const cs_lnum_t *b_nsub = n_sub_entities;

    BFT_MALLOC(b_order, n_ent, cs_lnum_t);

    cs_order_gnum_allocated(NULL,
                            b_gnum,
                            b_order,
                            n_ent);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    if (b_nsub != NULL) {

      current_gnum = b_nsub[b_order[0]];
      num_prev = b_gnum[b_order[0]];
      b_gnum[b_order[0]] = current_gnum;

      for (size_t i = 1; i < n_ent; i++) {
        num_cur = b_gnum[b_order[i]];
        if (num_cur > num_prev)
          current_gnum += b_nsub[b_order[i]];
        b_gnum[b_order[i]] = current_gnum;
        num_prev = num_cur;
      }

    }
    else { /* if (b_n_sub == NULL) */

      current_gnum = 1;
      num_prev = b_gnum[b_order[0]];
      b_gnum[b_order[0]] = current_gnum;

      for (size_t i = 1; i < n_ent; i++) {
        num_cur = b_gnum[b_order[i]];
        if (num_cur > num_prev)
          current_gnum += 1;
        b_gnum[b_order[i]] = current_gnum;
        num_prev = num_cur;
      }

    }

    BFT_FREE(b_order);

  }

  /* When sub-entities have been added, now switch from a numbering on
     the initial entities (shifted by number of sub-entities) to
     a numbering on the final sub-entities */

  _fvm_io_num_order_finalize(this_io_num,
                             n_sub_entities,
                             may_be_shared);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvm_io_num_local_max(this_io_num);
}

/*----------------------------------------------------------------------------
 * Copy selected shared global ordering information to private ordering
 * information for an I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to numbering structure
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_copy_on_write(fvm_io_num_t  *const this_io_num)
{
  if (this_io_num->_global_num == NULL) {
    cs_lnum_t i;
    BFT_MALLOC(this_io_num->_global_num,
               this_io_num->global_num_size,
               cs_gnum_t);
    for (i = 0; i < this_io_num->global_num_size; i++)
      this_io_num->_global_num[i] = this_io_num->global_num[i];
    this_io_num->global_num = this_io_num->_global_num;
  }
  assert(this_io_num->global_num == this_io_num->_global_num);
}

/*----------------------------------------------------------------------------
 * Copy selected shared global ordering information to private ordering
 * information for an I/O numbering structure.
 *
 * parameters:
 *   this_io_num          <-> pointer to numbering structure
 *   parent_global_number <-- pointer to shared list of global parent
 *                            entity numbers
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_try_to_set_shared(fvm_io_num_t      *const this_io_num,
                              const cs_gnum_t          parent_global_number[])
{
  if (this_io_num->_global_num != NULL && parent_global_number != NULL) {
    cs_lnum_t i;
    for (i = 0; i < this_io_num->global_num_size; i++)
      if (this_io_num->_global_num[i] != parent_global_number[i])
        break;
    if (i < this_io_num->global_num_size)
      this_io_num->global_num = this_io_num->_global_num;
    else {
      this_io_num->global_num = parent_global_number;
      BFT_FREE(this_io_num->_global_num);
    }
  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Maximum global number associated with an I/O numbering structure.
 *
 * This function is to be used for ordered global numberings.
 *
 * parameters:
 *   this_io_num <-- pointer to partially initialized I/O numbering structure.
 *   comm        <-- associated MPI communicator
 *
 * returns:
 *   maximum global number associated with the I/O numbering
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_fvm_io_num_global_max(const fvm_io_num_t  *const this_io_num,
                       const MPI_Comm             comm)
{
  cs_gnum_t  local_max = _fvm_io_num_local_max(this_io_num);
  cs_gnum_t  global_max = 0;

  MPI_Allreduce(&local_max, &global_max, 1, CS_MPI_GNUM, MPI_MAX, comm);

  return global_max;
}

/*----------------------------------------------------------------------------
 * Maximum global number associated with an I/O numbering structure.
 *
 * This function may be used for ordered global numberings.
 *
 * parameters:
 *   this_io_num <-- pointer to partially initialized I/O numbering structure.
 *   comm        <-- associated MPI communicator
 *
 * returns:
 *   maximum global number associated with the I/O numbering
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_fvm_io_num_global_max_unordered(const fvm_io_num_t  *const this_io_num,
                                 const MPI_Comm             comm)
{
  size_t     i, n_ent;
  cs_gnum_t  local_max = 0, global_max = 0;;

  /* Get maximum global number value */

  n_ent = this_io_num->global_num_size;
  for (i = 0; i < n_ent; i++) {
    cs_gnum_t local_val = this_io_num->global_num[i];
    if (local_val > local_max)
      local_max = local_val;
  }

  MPI_Allreduce(&local_max, &global_max, 1, CS_MPI_GNUM, MPI_MAX, comm);

  return global_max;
}

/*----------------------------------------------------------------------------
 * Global ordering associated with an I/O numbering structure.
 *
 * The structure should contain an initial ordering, which should
 * be sorted, but need not be contiguous. On output, the numbering
 * will be contiguous.
 *
 * As an option, a number of sub-entities per initial entity may be
 * given, in which case sub-entities of a same entity will have contiguous
 * numbers in the final ordering.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   n_sub_entities <-- optional number of sub-entities per initial entity,
 *                      or NULL if unused
 *   comm           <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_global_order(fvm_io_num_t       *this_io_num,
                         const cs_lnum_t     n_sub_entities[],
                         MPI_Comm            comm)
{
  cs_gnum_t   num_prev, num_cur;

  _Bool       may_be_shared = false;

  cs_lnum_t  *b_nsub = NULL;
  int         have_sub_loc = 0, have_sub_glob = 0;

  int         local_rank, size;

  cs_gnum_t   current_gnum = 0, gnum_shift = 0;

  /* Initialization */

  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &size);

  num_prev = 0;    /* true initialization later for block 0, */

  /* If numbering is shared, we will check later it was changed or if
     it can remain shared (in which case the copy may be discarded) */

  if (this_io_num->global_num != this_io_num->_global_num)
    may_be_shared = true;
  else
    may_be_shared = false;

  /* Get temporary maximum global number value */

  this_io_num->global_count = _fvm_io_num_global_max(this_io_num, comm);

  /* block_size = ceil(this_io_num->global_count/size) */

  cs_block_dist_info_t
    bi = cs_block_dist_compute_sizes(local_rank,
                                     size,
                                     1,
                                     0,
                                     this_io_num->global_count);

  cs_all_to_all_t
    *d = cs_all_to_all_create_from_block(this_io_num->global_num_size,
                                         0, /* flags */
                                         this_io_num->global_num,
                                         bi,
                                         comm);

  cs_gnum_t *b_gnum = cs_all_to_all_copy_array(d,
                                               CS_GNUM_TYPE,
                                               1,
                                               false, /* reverse */
                                               this_io_num->global_num,
                                               NULL);

  cs_lnum_t b_size = cs_all_to_all_n_elts_dest(d);

  /* Do we have sub-entities ? */

  if (n_sub_entities != NULL)
    have_sub_loc = 1;

  MPI_Allreduce(&have_sub_loc, &have_sub_glob, 1, MPI_INT, MPI_MAX, comm);

  if (have_sub_glob > 0)
    b_nsub = cs_all_to_all_copy_array(d,
                                      CS_LNUM_TYPE,
                                      1,
                                      false, /* reverse */
                                      n_sub_entities,
                                      NULL);

  if (b_size > 0) {

    cs_lnum_t *b_order;

    BFT_MALLOC(b_order, b_size, cs_lnum_t);

    cs_order_gnum_allocated(NULL,
                            b_gnum,
                            b_order,
                            b_size);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    if (have_sub_glob > 0) {

      current_gnum = b_nsub[b_order[0]];
      num_prev = b_gnum[b_order[0]];
      b_gnum[b_order[0]] = current_gnum;

      for (cs_lnum_t i = 1; i < b_size; i++) {
        num_cur = b_gnum[b_order[i]];
        if (num_cur > num_prev)
          current_gnum += b_nsub[b_order[i]];
        b_gnum[b_order[i]] = current_gnum;
        num_prev = num_cur;
      }

    }
    else { /* if (have_sub_glob == 0) */

      current_gnum = 1;
      num_prev = b_gnum[b_order[0]];
      b_gnum[b_order[0]] = current_gnum;

      for (cs_lnum_t i = 1; i < b_size; i++) {
        num_cur = b_gnum[b_order[i]];
        if (num_cur > num_prev)
          current_gnum += 1;
        b_gnum[b_order[i]] = current_gnum;
        num_prev = num_cur;
      }

    }

    BFT_FREE(b_order);

  }

  /* Partial clean-up */

  BFT_FREE(b_nsub);

  /* At this stage, b_gnum[] is valid for this process, and
     current_gnum indicates the total number of entities handled
     by this process; we must now shift global numberings on different
     processes by the cumulative total number of entities handled by
     each process */

  MPI_Scan(&current_gnum, &gnum_shift, 1, CS_MPI_GNUM,
           MPI_SUM, comm);
  gnum_shift -= current_gnum;

  for (cs_lnum_t i = 0; i < b_size; i++)
    b_gnum[i] += gnum_shift;

  /* Return global order to all ranks */

  cs_all_to_all_copy_array(d,
                           CS_GNUM_TYPE,
                           1,
                           true, /* reverse */
                           b_gnum,
                           this_io_num->_global_num);

  /* Free memory */

  BFT_FREE(b_gnum);

  cs_all_to_all_destroy(&d);

  /* When sub-entities have been added, now switch from a numbering on
     the initial entities (shifted by number of sub-entities) to
     a numbering on the final sub-entities */

  _fvm_io_num_order_finalize(this_io_num,
                             n_sub_entities,
                             may_be_shared);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvm_io_num_global_max(this_io_num, comm);
}

/*----------------------------------------------------------------------------
 * Global ordering associated with an I/O numbering structure.
 *
 * The structure does not need to contain an initial ordering,
 * though the array should be allocated.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   stride         <-- values per entity
 *   global_num     <-- global numbering array
 *   comm           <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_global_order_s(fvm_io_num_t       *this_io_num,
                           size_t              stride,
                           cs_gnum_t           global_num[],
                           MPI_Comm            comm)
{
  int  local_rank, size;
  cs_gnum_t current_gnum = 0, gnum_shift = 0;

  cs_gnum_t *r_gnum = NULL;

  /* Initialization */

  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &size);

  /* Get maximum global number value for first value of each series
     (does not need to be exact, simply used to define blocks) */

  {
    cs_gnum_t   local_max = 0, global_max = 0;
    size_t      n_ent = this_io_num->global_num_size;

    if (n_ent > 0)
      local_max = global_num[(n_ent-1)*stride];
    MPI_Allreduce(&local_max, &global_max, 1, CS_MPI_GNUM, MPI_MAX, comm);
    this_io_num->global_count = global_max;
  }

  /* block_info */

  cs_block_dist_info_t
    bi = cs_block_dist_compute_sizes(local_rank,
                                     size,
                                     1,
                                     0,
                                     this_io_num->global_count);

  const cs_gnum_t block_size = bi.block_size;

  int *dest_rank;
  BFT_MALLOC(dest_rank, this_io_num->global_num_size, int);

  for (cs_lnum_t i = 0; i < this_io_num->global_num_size; i++) {
    cs_gnum_t g_elt_id = global_num[stride*i] - 1;
    dest_rank[i] = g_elt_id / block_size;
  }

  cs_all_to_all_t
    *d = cs_all_to_all_create(this_io_num->global_num_size,
                              0,      /* flags */
                              NULL,  /* dest_id */
                              dest_rank,
                              comm);

  cs_all_to_all_transfer_dest_rank(d, &dest_rank);

  cs_gnum_t *b_gnum = cs_all_to_all_copy_array(d,
                                               CS_GNUM_TYPE,
                                               stride,
                                               false, /* reverse */
                                               global_num,
                                               NULL);

  cs_lnum_t b_size = cs_all_to_all_n_elts_dest(d);

  /* Order received data based on global number sets */

  if (b_size > 0) {

    cs_lnum_t *b_order = NULL;

    BFT_MALLOC(r_gnum, b_size, cs_gnum_t);
    BFT_MALLOC(b_order, b_size, cs_lnum_t);

    cs_order_gnum_allocated_s(NULL,
                              b_gnum,
                              stride,
                              b_order,
                              b_size);

    /* We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    current_gnum = 1;
    cs_lnum_t prev_id = b_order[0];
    r_gnum[b_order[0]] = current_gnum;

    const cs_lnum_t _stride = stride;

    for (cs_lnum_t i = 1; i < b_size; i++) {
      bool greater_than_prev = false;
      cs_lnum_t cur_id = b_order[i];
      for (cs_lnum_t j = 0; j < _stride; j++) {
        if (  b_gnum[cur_id*_stride + j]
            > b_gnum[prev_id*_stride + j])
          greater_than_prev = true;
      }
      if (greater_than_prev)
        current_gnum += 1;
      r_gnum[b_order[i]] = current_gnum;
      prev_id = cur_id;
    }

    BFT_FREE(b_order);

  }

  BFT_FREE(b_gnum);

  /* At this stage, r_gnum[] is valid for this process, and
     current_gnum indicates the total number of entities handled
     by this process; we must now shift global numberings on different
     processes by the cumulative total number of entities handled by
     each process */

  MPI_Scan(&current_gnum, &gnum_shift, 1, CS_MPI_GNUM,
           MPI_SUM, comm);
  gnum_shift -= current_gnum;

  for (cs_lnum_t i = 0; i < b_size; i++)
    r_gnum[i] += gnum_shift;

  /* Return global order to all ranks */

  cs_all_to_all_copy_array(d,
                           CS_GNUM_TYPE,
                           1,
                           true, /* reverse */
                           r_gnum,
                           this_io_num->_global_num);

  /* Partial clean-up */

  BFT_FREE(r_gnum);

  cs_all_to_all_destroy(&d);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvm_io_num_global_max(this_io_num, comm);
}

/*----------------------------------------------------------------------------
 * Compare two elements in an indexed list and returns true if element in
 * position i1 is strictly greater than element in position i2.
 *
 * parameters:
 *   i1        <-- position in index for the first element
 *   i2        <-- position in index for the second element
 *   index     <-- number of values to compare for each entity
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_indexed_is_greater(size_t            i1,
                    size_t            i2,
                    const cs_lnum_t   index[],
                    const cs_gnum_t   number[])
{
  int  i;

  cs_lnum_t   i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  cs_lnum_t   i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

  if (s1 > s2) {

    for (i = 0; i < s2; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return true;
  }
  else { /* s1 <= s2 */

    for (i = 0; i < s1; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return false;
  }

}

/*----------------------------------------------------------------------------
 * Global indexed ordering associated with an I/O numbering structure.
 *
 * The structure should contain an initial ordering, which should
 * be sorted, but need not be contiguous. On output, the numbering
 * will be contiguous.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   index          <-- index on entities for global_num[]
 *   global_num     <--
 *   comm           <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvm_io_num_global_order_index(fvm_io_num_t       *this_io_num,
                               cs_lnum_t           index[],
                               cs_gnum_t           global_num[],
                               MPI_Comm            comm)
{
  int  rank, local_rank, size;
  size_t  i, shift, block_size;

  cs_gnum_t   n_ent_recv = 0, n_ent_send = 0;
  cs_gnum_t   current_global_num = 0, global_num_shift = 0;
  int  *send_count = NULL, *recv_count = NULL;
  int  *send_shift = NULL, *recv_shift = NULL;
  cs_lnum_t   *recv_order = NULL, *recv_sub_index = NULL;
  cs_lnum_t   *recv_sub_count = NULL, *send_sub_count = NULL;
  cs_gnum_t   *block_global_num = NULL, *recv_global_num = NULL;

  /* Initialization */

  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &size);

  /* Get maximum global number value for first value of each series
     (does not need to be exact, simply used to define blocks) */

  {
    cs_gnum_t   local_max = 0, global_max = 0;
    size_t      n_ent = this_io_num->global_num_size;

    if (n_ent > 0)
      local_max = global_num[index[n_ent-1]];
    MPI_Allreduce(&local_max, &global_max, 1, CS_MPI_GNUM, MPI_MAX, comm);
    this_io_num->global_count = global_max;
  }

  /* block_size = ceil(this_io_num->global_count/size) */

  block_size = this_io_num->global_count / size;
  if (this_io_num->global_count % size > 0)
    block_size += 1;

  /* Build for each block, a new ordered indexed list from the received
     elements */

  assert(sizeof(cs_gnum_t) >= sizeof(cs_lnum_t));

  BFT_MALLOC(send_count, size, int);
  BFT_MALLOC(recv_count, size, int);
  BFT_MALLOC(send_shift, size + 1, int);
  BFT_MALLOC(recv_shift, size + 1, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++) {
    rank = (global_num[index[i]] - 1) / block_size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  /* Get recv_index */

  n_ent_recv = recv_shift[size];
  n_ent_send = send_shift[size];

  BFT_MALLOC(recv_sub_count, n_ent_recv, cs_lnum_t);
  BFT_MALLOC(send_sub_count, n_ent_send, cs_lnum_t);

  assert(n_ent_send == (cs_gnum_t)this_io_num->global_num_size);

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)this_io_num->global_num_size; i++) {
    rank = (global_num[index[i]] - 1) / block_size;
    shift = send_shift[rank] + send_count[rank];
    send_sub_count[shift] = index[i+1] - index[i];
    send_count[rank] +=  1;
  }

  MPI_Alltoallv(send_sub_count, send_count, send_shift, CS_MPI_LNUM,
                recv_sub_count, recv_count, recv_shift, CS_MPI_LNUM, comm);

  BFT_MALLOC(recv_sub_index, n_ent_recv + 1, cs_lnum_t);

  recv_sub_index[0] = 0;
  for (i = 0; i < n_ent_recv; i++)
      recv_sub_index[i+1] = recv_sub_index[i] + recv_sub_count[i];

  BFT_FREE(send_sub_count);
  BFT_FREE(recv_sub_count);

  /* Get recv_global_num */

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++) {
    rank = (global_num[index[i]] - 1) / block_size;
    send_count[rank] += index[i+1] - index[i];
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  BFT_MALLOC(recv_global_num, recv_sub_index[n_ent_recv], cs_gnum_t);

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the indexed list to the correct processors */

  MPI_Alltoallv(global_num, send_count, send_shift, CS_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, CS_MPI_GNUM, comm);

  if (n_ent_recv > 0) { /* Order received elements of the indexed list */

    size_t prev_id, cur_id;

    BFT_MALLOC(recv_order, n_ent_recv, cs_lnum_t);

    cs_order_gnum_allocated_i(NULL,
                              recv_global_num,
                              recv_sub_index,
                              recv_order,
                              n_ent_recv);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of elements */

    BFT_MALLOC(block_global_num, n_ent_recv, cs_gnum_t);

    current_global_num = 1;
    prev_id = recv_order[0];
    block_global_num[recv_order[0]] = current_global_num;

    for (i = 1; i < n_ent_recv; i++) {

      cur_id = recv_order[i];

      if (_indexed_is_greater(cur_id, prev_id, recv_sub_index, recv_global_num))
        current_global_num += 1;

      block_global_num[recv_order[i]] = current_global_num;
      prev_id = cur_id;

    }

  } /* End if n_ent_recv > 0 */

  /* Partial clean-up */

  BFT_FREE(recv_order);
  BFT_FREE(recv_sub_index);
  BFT_FREE(recv_global_num);

  /* At this stage, block_global_num[] is valid for this process, and
     current_global_num indicates the total number of entities handled
     by this process; we must now shift global numberings on different
     processes by the cumulative total number of entities handled by
     each process */

  MPI_Scan(&current_global_num, &global_num_shift, 1, CS_MPI_GNUM,
           MPI_SUM, comm);
  global_num_shift -= current_global_num;

  for (i = 0; i < n_ent_recv; i++)
    block_global_num[i] += global_num_shift;

  /* Return global order to all processors */

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++) {
    rank = (global_num[index[i]] - 1) / block_size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  MPI_Alltoallv(block_global_num, recv_count, recv_shift, CS_MPI_GNUM,
                this_io_num->_global_num, send_count, send_shift, CS_MPI_GNUM,
                comm);

  /* Free memory */

  BFT_FREE(block_global_num);
  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvm_io_num_global_max(this_io_num, comm);
}

/*----------------------------------------------------------------------------
 * Return the global number of sub-entities associated with an initial
 * entity whose global numbering is known, given the number of
 * sub-entities per initial entity.
 *
 * parameters:
 *   this_io_num    <-- pointer to base io numbering
 *   n_sub_entities <-- number of sub-entities per initial entity
 *   comm           <-- associated MPI communicator
 *
 * returns:
 *   global number of sub-entities
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_fvm_io_num_global_sub_size(const fvm_io_num_t  *this_io_num,
                            const cs_lnum_t      n_sub_entities[],
                            MPI_Comm             comm)
{

  cs_gnum_t   global_count, n_ent_recv, num_prev, num_cur;
  size_t      i, block_size;
  int         rank;

  cs_gnum_t   *recv_global_num = NULL;
  cs_gnum_t   *send_global_num = NULL;
  cs_lnum_t   *recv_n_sub = NULL, *recv_order = NULL;
  int         *send_count = NULL, *recv_count = NULL;
  int         *send_shift = NULL, *recv_shift = NULL;
  int         have_sub_loc = 0, have_sub_glob = 0;

  int         size;

  cs_gnum_t   current_global_num = 0;
  cs_gnum_t   retval = 0;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  num_prev = 0;    /* true initialization later for block 0, */

  /* Get temporary maximum global number value */

  global_count = _fvm_io_num_global_max(this_io_num, comm);

  /* block_size = ceil(this_io_num->global_count/size) */

  block_size = global_count / size;
  if (global_count % size > 0)
    block_size += 1;

  assert(sizeof(cs_gnum_t) >= sizeof(cs_lnum_t));

  BFT_MALLOC(send_count, size, int);
  BFT_MALLOC(recv_count, size, int);

  BFT_MALLOC(send_shift, size, int);
  BFT_MALLOC(recv_shift, size, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
    send_count[(this_io_num->global_num[i] - 1) / block_size] += 1;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 1; rank < size; rank++) {
    send_shift[rank] = send_shift[rank - 1] + send_count[rank -1];
    recv_shift[rank] = recv_shift[rank - 1] + recv_count[rank -1];
  }

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the n_sub_entities[] array to the correct processors */

  n_ent_recv = recv_shift[size - 1] + recv_count[size - 1];

  BFT_MALLOC(recv_global_num, n_ent_recv, cs_gnum_t);
  BFT_MALLOC(recv_order, n_ent_recv, cs_lnum_t);

  if (this_io_num->_global_num != NULL)
    send_global_num = this_io_num->_global_num;
  else {
    BFT_MALLOC(send_global_num,
               this_io_num->global_num_size,
               cs_gnum_t);
    memcpy(send_global_num,
           this_io_num->global_num,
           this_io_num->global_num_size * sizeof(cs_gnum_t));
  }

  MPI_Alltoallv(send_global_num, send_count, send_shift, CS_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, CS_MPI_GNUM, comm);

  if (send_global_num != this_io_num->_global_num)
    BFT_FREE(send_global_num);

  /* Do we have sub-entities ? */

  if (n_sub_entities != NULL)
    have_sub_loc = 1;

  MPI_Allreduce(&have_sub_loc, &have_sub_glob, 1, MPI_INT, MPI_MAX, comm);

  if (have_sub_glob > 0) {

    cs_lnum_t   *send_n_sub;

    BFT_MALLOC(send_n_sub, this_io_num->global_num_size, cs_lnum_t);
    BFT_MALLOC(recv_n_sub, n_ent_recv, cs_lnum_t);

    if (n_sub_entities != NULL) {
      for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
        send_n_sub[i] = n_sub_entities[i];
    }
    else {
      for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
        send_n_sub[i] = 1;
    }

    MPI_Alltoallv(send_n_sub, send_count, send_shift, CS_MPI_LNUM,
                  recv_n_sub, recv_count, recv_shift, CS_MPI_LNUM, comm);

    BFT_FREE(send_n_sub);
  }

  if (n_ent_recv > 0) {

    cs_order_gnum_allocated(NULL,
                            recv_global_num,
                            recv_order,
                            n_ent_recv);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    current_global_num = recv_n_sub[recv_order[0]];
    num_prev = recv_global_num[recv_order[0]];
    recv_global_num[recv_order[0]] = current_global_num;

    for (i = 1; i < n_ent_recv; i++) {
      num_cur = recv_global_num[recv_order[i]];
      if (num_cur > num_prev)
        current_global_num += recv_n_sub[recv_order[i]];
      num_prev = num_cur;
    }

  }

  /* Partial clean-up */

  BFT_FREE(recv_n_sub);
  BFT_FREE(recv_order);
  BFT_FREE(recv_global_num);

  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);

  /* At this stage, current_global_num indicates the total number of
     entities handled by this process; we must now shift global
     numberings on different processes by the cumulative total
     number of entities handled by each process */

  MPI_Allreduce(&current_global_num, &retval, 1, CS_MPI_GNUM, MPI_SUM, comm);

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Stretch extents in some directions.
 *
 * If the multiplication factor for a given axis direction is less than 1,
 * the extents in that direction will be defined so that the extent box's
 * size in that direction is its maximum size.
 *
 * parameters:
 *   extents     <-> box extents
 *   box_to_cube <-- if 1, transform bounding box to boundign cube
 *----------------------------------------------------------------------------*/

static void
_adjust_extents(cs_coord_t  extents[6],
                int         box_to_cube)
{
  size_t  i;
  cs_coord_t max_width = 0.;
  const double epsilon = 1e-12;

  for (i = 0; i < 3; i++) {
    double  w = fabs(extents[i+3] - extents[i]);
    max_width = CS_MAX(max_width, w);
  }

  for (i = 0; i < 3; i++) {
    double  mult = 1.0;
    double  m = (extents[i] + extents[i+3])*0.5;
    double  w = fabs(extents[i+3] - extents[i]);
    if (box_to_cube > 0) {
      if (w > epsilon)
        mult = max_width / w;
    }
    if (mult < 1.0 + epsilon)
      mult= 1.0 + epsilon;
    extents[i] = m - (w*0.5*mult);
    extents[i+3] = m + (w*0.5*mult);
  }
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * parameters:
 *   coords      <-- pointer to entity coordinates (interlaced)
 *   dim         <-- spatial dimension
 *   n_entities  <-- number of entities considered
 *   box_to_cube <-- if 1, use bounding cube instead of bounding box
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

static fvm_io_num_t *
_create_from_coords_morton(const cs_coord_t  coords[],
                           int               dim,
                           size_t            n_entities,
                           int               box_to_cube)
{
  size_t i;
  cs_coord_t extents[6];

#if defined(HAVE_MPI)
  MPI_Comm comm = cs_glob_mpi_comm;
#endif

  const int level = sizeof(fvm_morton_int_t)*8 - 1;
  const int n_ranks = cs_glob_n_ranks;

  fvm_io_num_t  *this_io_num = NULL;

  /* Create structure */

  BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

  this_io_num->global_num_size = n_entities;

  BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  /* Build Morton encoding and order it */

#if defined(HAVE_MPI)
  fvm_morton_get_coord_extents(dim, n_entities, coords, extents, comm);
#else
  fvm_morton_get_coord_extents(dim, n_entities, coords, extents);
#endif

  _adjust_extents(extents, box_to_cube);

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    int *dest_rank = NULL;
    cs_lnum_t *order = NULL;
    fvm_morton_code_t *m_code = NULL;
    int input[1] = {dim};

    BFT_MALLOC(m_code, n_entities, fvm_morton_code_t);
    BFT_MALLOC(order, n_entities, cs_lnum_t);
    BFT_MALLOC(dest_rank, n_entities, int);

    fvm_morton_encode_coords(dim, level, extents, n_entities, coords, m_code);
    fvm_morton_local_order(n_entities, m_code, order);

    cs_sort_partition_dest_rank_id(_sampling_factors[dim],
                                   sizeof(fvm_morton_code_t),
                                   n_entities,
                                   m_code,
                                   NULL, /* weight */
                                   order,
                                   dest_rank,
                                   fvm_morton_s_to_code,
                                   fvm_morton_compare_o,
                                   input,
                                   comm);

    BFT_FREE(order);
    BFT_FREE(m_code);

    cs_all_to_all_t
      *d = cs_all_to_all_create(this_io_num->global_num_size,
                                0,     /* flags */
                                NULL,  /* dest_id */
                                dest_rank,
                                comm);

    cs_all_to_all_transfer_dest_rank(d, &dest_rank);

    cs_real_t *b_coords
      = cs_all_to_all_copy_array(d,
                                 CS_REAL_TYPE,
                                 3,
                                 false, /* reverse */
                                 coords,
                                 NULL);

    size_t b_size = cs_all_to_all_n_elts_dest(d);

    /* Now re-build Morton codes on block distribution
       (exchanging and ordering Morton codes directly would be more
       efficient, but using coordinates allows breaking ties for
       possibly identical codes through lexicographical ordering). */

    BFT_MALLOC(order, b_size, cs_lnum_t);
    BFT_MALLOC(m_code, b_size, fvm_morton_code_t);

    fvm_morton_encode_coords(dim,
                             level,
                             extents,
                             b_size,
                             b_coords,
                             m_code);

    fvm_morton_local_order(b_size, m_code, order);

    /* Check ordering; if two entities have the same Morton codes,
       use lexicographical coordinates ordering to ensure the
       final order is deterministic. */

    _check_morton_ordering(dim, b_size, b_coords, m_code, order);

    BFT_FREE(m_code);
    BFT_FREE(b_coords);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    cs_gnum_t *b_gnum;
    BFT_MALLOC(b_gnum, b_size, cs_gnum_t);

    for (i = 0; i < b_size; i++)
      b_gnum[order[i]] = i+1;

    BFT_FREE(order);

    cs_gnum_t gnum_shift = 0, current_gnum = b_size;

    /* At this stage, b_gnum[] is valid for this process, and
       current_gnum indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    MPI_Scan(&current_gnum, &gnum_shift, 1, CS_MPI_GNUM,
             MPI_SUM, comm);
    gnum_shift -= current_gnum;

    for (i = 0; i < b_size; i++)
      b_gnum[i] += gnum_shift;

    /* Return global order to all ranks */

    cs_all_to_all_copy_array(d,
                             CS_GNUM_TYPE,
                             1,
                             true, /* reverse */
                             b_gnum,
                             this_io_num->_global_num);

    /* Free memory */

    BFT_FREE(b_gnum);

    cs_all_to_all_destroy(&d);

    /* Get final maximum global number value */

    this_io_num->global_count
      = _fvm_io_num_global_max_unordered(this_io_num, comm);

  }

#endif /* HAVE_MPI */

  if (n_ranks == 1) {

    cs_lnum_t *order = NULL;
    fvm_morton_code_t *m_code = NULL;

    BFT_MALLOC(m_code, n_entities, fvm_morton_code_t);
    BFT_MALLOC(order, n_entities, cs_lnum_t);

    fvm_morton_encode_coords(dim, level, extents, n_entities, coords, m_code);
    fvm_morton_local_order(n_entities, m_code, order);

    _check_morton_ordering(dim, n_entities, coords, m_code, order);

    BFT_FREE(m_code);

    for (i = 0; i < n_entities; i++)
      this_io_num->_global_num[order[i]] = i+1;

    BFT_FREE(order);

    this_io_num->global_count = n_entities;

  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Hilbert curve, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Hilbert code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * parameters:
 *   coords      <-- pointer to entity coordinates (interlaced)
 *   dim         <-- spatial dimension
 *   n_entities  <-- number of entities considered
 *   box_to_cube <-- if 1, use bounding cube instead of bounding box
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

static fvm_io_num_t *
_create_from_coords_hilbert(const cs_coord_t  coords[],
                            int               dim,
                            size_t            n_entities,
                            int               box_to_cube)
{
  size_t i;
  cs_coord_t extents[6];

#if defined(HAVE_MPI)
  MPI_Comm comm = cs_glob_mpi_comm;
#endif

  const int n_ranks = cs_glob_n_ranks;

  fvm_io_num_t  *this_io_num = NULL;

  /* Create structure */

  BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

  this_io_num->global_num_size = n_entities;

  BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  /* Build Hilbert encoding and order it */

#if defined(HAVE_MPI)
  fvm_hilbert_get_coord_extents(dim, n_entities, coords, extents, comm);
#else
  fvm_hilbert_get_coord_extents(dim, n_entities, coords, extents);
#endif

  _adjust_extents(extents, box_to_cube);

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    int *dest_rank = NULL;
    cs_lnum_t *order = NULL;
    fvm_hilbert_code_t *h_code = NULL;

    BFT_MALLOC(h_code, n_entities, fvm_hilbert_code_t);
    BFT_MALLOC(order, n_entities, cs_lnum_t);
    BFT_MALLOC(dest_rank, n_entities, int);

    fvm_hilbert_encode_coords(dim, extents, n_entities, coords, h_code);
    fvm_hilbert_local_order(n_entities, h_code, order);

    cs_sort_partition_dest_rank_id(_sampling_factors[dim],
                                   sizeof(fvm_hilbert_code_t),
                                   n_entities,
                                   h_code,
                                   NULL, /* weight */
                                   order,
                                   dest_rank,
                                   fvm_hilbert_s_to_code,
                                   fvm_hilbert_compare,
                                   NULL,
                                   comm);

    BFT_FREE(order);
    BFT_FREE(h_code);

    cs_all_to_all_t
      *d = cs_all_to_all_create(this_io_num->global_num_size,
                                0,     /* flags */
                                NULL,  /* dest_id */
                                dest_rank,
                                comm);

    cs_all_to_all_transfer_dest_rank(d, &dest_rank);

    cs_real_t *b_coords
      = cs_all_to_all_copy_array(d,
                                 CS_REAL_TYPE,
                                 3,
                                 false, /* reverse */
                                 coords,
                                 NULL);

    size_t b_size = cs_all_to_all_n_elts_dest(d);

    /* Now re-build Hilbert coords on block distribution
       (exchanging and ordering Hilbert codes directly would be more
       efficient, but using coordinates allows breaking ties for
       possibly identical codes through lexicographical ordering). */

    BFT_MALLOC(order, b_size, cs_lnum_t);

    fvm_hilbert_local_order_coords(dim,
                                   extents,
                                   b_size,
                                   b_coords,
                                   order);

    BFT_FREE(b_coords);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    cs_gnum_t *b_gnum;
    BFT_MALLOC(b_gnum, b_size, cs_gnum_t);

    for (i = 0; i < b_size; i++)
      b_gnum[order[i]] = i+1;

    BFT_FREE(order);

    cs_gnum_t gnum_shift = 0, current_gnum = b_size;

    /* At this stage, b_gnum[] is valid for this process, and
       current_gnum indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    MPI_Scan(&current_gnum, &gnum_shift, 1, CS_MPI_GNUM,
             MPI_SUM, comm);
    gnum_shift -= current_gnum;

    for (i = 0; i < b_size; i++)
      b_gnum[i] += gnum_shift;

    /* Return global order to all ranks */

    cs_all_to_all_copy_array(d,
                             CS_GNUM_TYPE,
                             1,
                             true, /* reverse */
                             b_gnum,
                             this_io_num->_global_num);

    /* Free memory */

    BFT_FREE(b_gnum);

    cs_all_to_all_destroy(&d);

    /* Get final maximum global number value */

    this_io_num->global_count
      = _fvm_io_num_global_max_unordered(this_io_num, comm);

  }

#endif /* HAVE_MPI */

  if (n_ranks == 1) {

    cs_lnum_t *order = NULL;
    BFT_MALLOC(order, n_entities, cs_lnum_t);

    fvm_hilbert_local_order_coords(dim,
                                   extents,
                                   n_entities,
                                   coords,
                                   order);

    for (i = 0; i < n_entities; i++)
      this_io_num->_global_num[order[i]] = i+1;

    BFT_FREE(order);

    this_io_num->global_count = n_entities;

  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Test if an array of global numbers is lexicographically ordered.
 *
 * parameters:
 *   parent_entity_id <-- optional list (0 to n-1 numbering) of selected
 *                        entities (or NULL if all n_entities are selected).
 *                        This list may contain element ids in any order
 *   number           <-- array of all entity numbers (number of entity i given
 *                        by number[i] or number[list[i] - 1]) if list exists
 *                        (if NULL, a default 1 to n numbering is considered)
 *   stride           <-- stride of number array (number of values to compare)
 *   n_entities       <-- number of entities considered
 *
 * returns:
 *   1 if ordered, 0 otherwise.
 *----------------------------------------------------------------------------*/

static int
_is_gnum_ordered_s(const cs_lnum_t  parent_entity_id[],
                   const cs_gnum_t  number[],
                   size_t           stride,
                   size_t           n_entities)
{
  size_t j;
  size_t i = 0;

  /* If numbering is explicit */

  if (number != NULL) {

    if (parent_entity_id != NULL) {
      for (i = 1 ; i < n_entities ; i++) {
        size_t j_prev, k;
        bool unordered = false;
        j_prev = parent_entity_id[i-1];
        j = parent_entity_id[i];
        for (k = 0; k < stride; k++) {
          if (number[j_prev*stride + k] < number[j*stride + k])
            break;
          else if (number[j_prev*stride + k] > number[j*stride + k])
            unordered = true;
        }
        if (unordered == true)
          break;
      }
    }
    else {
      for (i = 1 ; i < n_entities ; i++) {
        size_t i_prev, k;
        bool unordered = false;
        i_prev = i-1;
        for (k = 0; k < stride; k++) {
          if (number[i_prev*stride + k] < number[i*stride + k])
            break;
          else if (number[i_prev*stride + k] > number[i*stride + k])
            unordered = true;
        }
        if (unordered == true)
          break;
      }
    }

  /* If numbering is implicit */

  }
  else {

    if (parent_entity_id != NULL) {
      for (i = 1 ; i < n_entities ; i++) {
        if (parent_entity_id[i] < parent_entity_id[i-1])
          break;
      }
    }
    else
      i = n_entities;
  }

  if (i == n_entities || n_entities == 0)
    return 1;
  else
    return 0;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure.
 *
 * This function is similar to fvm_io_num_create_from_select, albeit
 * using parent entity numbers (1 to n) instead of ids (0 to n-1).
 *
 * parameters:
 *   parent_entity_number <-- pointer to list of selected entitie's parent's
 *                            numbers, or NULL if all first n_ent entities
 *                            are used
 *   parent_global_number <-- pointer to list of global (i.e. domain splitting
 *                            independent) parent entity numbers
 *   n_entities           <-- number of entities considered
 *   share_parent_global  <-- if non zero, try to share parent_global_number
 *                            instead of using a local copy
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create(const cs_lnum_t   parent_entity_number[],
                  const cs_gnum_t   parent_global_number[],
                  size_t            n_entities,
                  int               share_parent_global)
{
  cs_lnum_t *parent_entity_id = NULL;

  if (parent_entity_number != NULL) {
    BFT_MALLOC(parent_entity_id, n_entities, cs_lnum_t);
    for (size_t i = 0; i < n_entities; i++)
      parent_entity_id[i] = parent_entity_number[i] - 1;
  }

  fvm_io_num_t  *this_io_num
    = fvm_io_num_create_from_select(parent_entity_id,
                                    parent_global_number,
                                    n_entities,
                                    share_parent_global);

  BFT_FREE(parent_entity_id);

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a selection of entities.
 *
 * parameters:
 *   parent_entity_id     <-- pointer to list of selected entitie's parent's
 *                            ids, or NULL if all first n_ent entities
 *                            are used
 *   parent_global_number <-- pointer to list of global (i.e. domain splitting
 *                            independent) parent entity numbers
 *   n_entities           <-- number of entities considered
 *   share_parent_global  <-- if non zero, try to share parent_global_number
 *                            instead of using a local copy
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_select(const cs_lnum_t   parent_entity_id[],
                              const cs_gnum_t   parent_global_number[],
                              size_t            n_entities,
                              int               share_parent_global)
{
  size_t  i;
  cs_lnum_t  *order = NULL;

  fvm_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (cs_glob_n_ranks < 2 && parent_global_number == NULL)
    return NULL;

  /* Create structure */

  BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

  this_io_num->global_num_size = n_entities;

  BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  if (n_entities > 0) {

    /* Assign initial global numbers */

    if (parent_entity_id != NULL) {
      for (i = 0 ; i < n_entities ; i++)
        this_io_num->_global_num[i] = parent_global_number[parent_entity_id[i]];
    }
    else {
      for (i = 0 ; i < n_entities ; i++)
        this_io_num->_global_num[i] = parent_global_number[i];
    }

    if (cs_order_gnum_test(NULL,
                           this_io_num->_global_num,
                           n_entities) == false) {
      cs_gnum_t *tmp_num;
      order = cs_order_gnum(NULL,
                            this_io_num->_global_num,
                            n_entities);
      BFT_MALLOC(tmp_num, n_entities, cs_gnum_t);
      for (i = 0; i < n_entities; i++)
        tmp_num[i] = this_io_num->_global_num[order[i]];
      memcpy(this_io_num->_global_num, tmp_num, n_entities*sizeof(cs_gnum_t));
      BFT_FREE(tmp_num);
    }
  }

  this_io_num->global_count = n_entities;

  _fvm_io_num_copy_on_write(this_io_num);

  /* Order globally */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _fvm_io_num_global_order(this_io_num,
                             NULL,
                             cs_glob_mpi_comm);

#endif

  if (cs_glob_n_ranks == 1)
    _fvm_io_num_local_order(this_io_num,
                            NULL);

  if (order != NULL) {
    cs_gnum_t *tmp_num;
    BFT_MALLOC(tmp_num, n_entities, cs_gnum_t);
    for (i = 0; i < n_entities; i++)
      tmp_num[order[i]] = this_io_num->_global_num[i];
    memcpy(this_io_num->_global_num, tmp_num, n_entities*sizeof(cs_gnum_t));
    BFT_FREE(tmp_num);
    BFT_FREE(order);
  }

  if (share_parent_global != 0)
    _fvm_io_num_try_to_set_shared(this_io_num,
                                  parent_global_number);

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure,
 * sharing a given global numbering array.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   global_number <-- pointer to list of global (i.e. domain splitting
 *                     independent) entity numbers
 *   global_count  <-- global number of entities
 *   n_entities    <-- number of local entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_shared(const cs_gnum_t   global_number[],
                         cs_gnum_t         global_count,
                         size_t            n_entities)
{
  fvm_io_num_t  *this_io_num;

  /* Create structure */

  BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

  this_io_num->global_count = global_count;
  this_io_num->global_num_size = n_entities;

  this_io_num->global_num = global_number;
  this_io_num->_global_num = NULL;

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on an an initial
 * I/O numbering and a number of new entities per base entity.
 *
 * This is useful for example to create an I/O numbering for
 * triangles based on split polygons, whose I/O numbering is defined.
 *
 * parameters:
 *   base_io_num    <-- pointer to base I/O numbering structure
 *   n_sub_entities <-- number of new entities per base entity
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_sub(const fvm_io_num_t  *base_io_num,
                           const cs_lnum_t      n_sub_entities[])
{
  fvm_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (base_io_num == NULL)
    return NULL;

  /* Create structure */

  BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

  cs_lnum_t n_ent = base_io_num->global_num_size;

  BFT_MALLOC(this_io_num->_global_num, n_ent, cs_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  this_io_num->global_num_size = n_ent;

  /* Assign initial global numbers */

  for (cs_lnum_t i = 0 ; i < n_ent ; i++)
    this_io_num->_global_num[i] = base_io_num->global_num[i];

  this_io_num->global_count = n_ent;

  _fvm_io_num_copy_on_write(this_io_num);

  /* Order globally */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _fvm_io_num_global_order(this_io_num,
                             n_sub_entities,
                             cs_glob_mpi_comm);

#endif

  if (cs_glob_n_ranks == 1)
    _fvm_io_num_local_order(this_io_num,
                            n_sub_entities);

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a strided adjacency.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   parent_entity_id <-- pointer to list of selected entitie's parent's ids,
 *                        or NULL if all first n_ent entities are used
 *   index            <-- index on entities for adjacency
 *   adjacency        <-- entity adjacency (1 to n global numbering)
 *   n_entities       <-- number of entities considered
 *
 * returns:
 *   pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_adj_s(const cs_lnum_t   parent_entity_id[],
                             const cs_gnum_t   adjacency[],
                             size_t            n_entities,
                             size_t            stride)
{
  fvm_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (cs_glob_n_ranks < 2)
    return NULL;

  assert(_is_gnum_ordered_s(parent_entity_id,
                            adjacency,
                            stride,
                            n_entities) == true);

#if defined(HAVE_MPI)
  {
    cs_gnum_t *_adjacency = NULL;

    /* Create structure */

    BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

    this_io_num->global_num_size = n_entities;

    BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    if (n_entities > 0) {

      size_t  i, j;

      /* Assign initial global numbers */

      BFT_MALLOC(_adjacency, n_entities*stride, cs_gnum_t);

      if (parent_entity_id != NULL) {
        for (i = 0 ; i < n_entities ; i++) {
          for (j = 0; j < stride; j++)
            _adjacency[i*stride + j]
              = adjacency[parent_entity_id[i]*stride + j];
        }
      }
      else
        memcpy(_adjacency, adjacency, n_entities*stride*sizeof(cs_gnum_t));

    }

    /* Order globally */

    this_io_num->global_count = n_entities;

    _fvm_io_num_global_order_s(this_io_num,
                               stride,
                               _adjacency,
                               cs_glob_mpi_comm);

    BFT_FREE(_adjacency);
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on an indexed adjacency.
 *
 * The corresponding entities do not need to be locally ordered.
 *
 * parameters:
 *   parent_entity_id <-- pointer to list of selected entitie's parent's ids,
 *                        or NULL if all first n_ent entities are used
 *   index            <-- index on entities for adjacency
 *   adjacency        <-- entity adjacency (1 to n global numbering)
 *   n_entities       <-- number of entities considered
 *
 * returns:
 *   pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_adj_i(const cs_lnum_t   parent_entity_id[],
                             const cs_lnum_t   index[],
                             const cs_gnum_t   adjacency[],
                             cs_lnum_t         n_entities)
{
  fvm_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (cs_glob_n_ranks < 2)
    return NULL;

#if defined(HAVE_MPI)
  {
    cs_lnum_t   *_index = NULL;
    cs_gnum_t *_adjacency = NULL;

#if defined(DEBUG) && !defined(NDEBUG)
    const char no_adjacent_elt_msg[]
      = " Error during the creation of a fvm_io_num_t structure\n"
        " for an indexed adjacency.\n\n"
        " At least one entity has no adjacent element, and\n"
        " this case is not currently handled.\n\n"
        " Check if this definition is correct.";
#endif

    /* Create structure */

    BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

    this_io_num->global_num_size = n_entities;

    BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    if (n_entities > 0) {

      cs_lnum_t   i, j, k, ent_id, _shift;

      /* Assign initial global numbers */

      BFT_MALLOC(_index, n_entities + 1, cs_lnum_t);
      _index[0] = 0;

      if (parent_entity_id != NULL) {

#if defined(DEBUG) && !defined(NDEBUG)
        for (i = 0 ; i < n_entities ; i++) {
          ent_id = parent_entity_id[i];
          if ((index[ent_id+1] - index[ent_id]) == 0)
            bft_error(__FILE__, __LINE__, 0, no_adjacent_elt_msg);
        }
#endif

        /* Count reduced size */

        for (i = 0 ; i < n_entities ; i++) {
          ent_id = parent_entity_id[i];
          _index[i+1] = index[ent_id+1] - index[ent_id];
        }

        for (i = 0 ; i < n_entities ; i++)
          _index[i+1] += _index[i];

        BFT_MALLOC(_adjacency, _index[n_entities], cs_gnum_t);

        /* Define reduced index and adjacency */

        for (i = 0 ; i < n_entities ; i++) {

          ent_id = parent_entity_id[i];
          _shift = _index[i];

          for (j = index[ent_id], k = 0; j < index[ent_id+1]; j++, k++)
            _adjacency[_shift + k] = adjacency[j];

        }

      }
      else {

#if defined(DEBUG) && !defined(NDEBUG)
        for (i = 0 ; i < n_entities ; i++) {
          if ((index[i+1] - index[i]) == 0)
            bft_error(__FILE__, __LINE__, 0, no_adjacent_elt_msg);
        }
#endif

        BFT_MALLOC(_adjacency, index[n_entities], cs_gnum_t);

        memcpy(_index, index, (n_entities+1)*sizeof(cs_lnum_t));
        memcpy(_adjacency, adjacency, index[n_entities]*sizeof(cs_gnum_t));

      }

    }

    /* Order globally */

    this_io_num->global_count = n_entities;

    _fvm_io_num_global_order_index(this_io_num,
                                   _index,
                                   _adjacency,
                                   cs_glob_mpi_comm);

    if (_adjacency != NULL)
      BFT_FREE(_adjacency);
    if (_index != NULL)
      BFT_FREE(_index);
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a space-filling curve.
 *
 * It is expected that entities are unique (i.e. not duplicated on 2 or
 * more ranks). If 2 entities have a same Morton or Hilbert code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * parameters:
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   dim        <-- spatial dimension
 *   n_entities <-- number of entities considered
 *   sfc_type   <-- type of space-filling curve (Morton or Hilbert)
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_sfc(const cs_coord_t  coords[],
                           int               dim,
                           size_t            n_entities,
                           fvm_io_num_sfc_t  sfc_type)
{
  fvm_io_num_t  *this_io_num = NULL;

  switch(sfc_type) {
  case FVM_IO_NUM_SFC_MORTON_BOX:
    this_io_num = _create_from_coords_morton(coords, dim, n_entities, 0);
    break;
  case FVM_IO_NUM_SFC_MORTON_CUBE:
    this_io_num = _create_from_coords_morton(coords, dim, n_entities, 1);
    break;
  case FVM_IO_NUM_SFC_HILBERT_BOX:
    this_io_num = _create_from_coords_hilbert(coords, dim, n_entities, 0);
    break;
  case FVM_IO_NUM_SFC_HILBERT_CUBE:
    this_io_num = _create_from_coords_hilbert(coords, dim, n_entities, 1);
    break;
  default:
    assert(0);
  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on real values, assuming
 * ordering by increasing values.
 *
 * It is expected that entities are unique (i.e. not duplicated on 2 or
 * more ranks). If 2 entities have a same value, their global
 * number will be determined by their initial order.
 *
 * parameters:
 *   val        <-- pointer to real values
 *   n_entities <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_real(const cs_real_t  val[],
                            size_t           n_entities)
{
  size_t i;

#if defined(HAVE_MPI)
  MPI_Comm comm = cs_glob_mpi_comm;
#endif

  const int n_ranks = cs_glob_n_ranks;

  fvm_io_num_t  *this_io_num = NULL;

  /* Create structure */

  BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

  this_io_num->global_num_size = n_entities;

  BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  /* Scale values */

  double v_min = DBL_MAX;
  double v_max = -DBL_MAX;

  for (i = 0; i < n_entities; i++) {
    if (val[i] < v_min) v_min = val[i];
    if (val[i] > v_max) v_max = val[i];
  }

  cs_parall_min(1, CS_DOUBLE, &v_min);
  cs_parall_max(1, CS_DOUBLE, &v_max);

  if (v_max <= v_min)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: point set is empty or contains identical values."),
              __func__);

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    double scale = (1.0 - 1.e12) / (v_max - v_min);

    cs_real_t *s_val;
    BFT_MALLOC(s_val, n_entities, cs_real_t);

    for (i = 0; i < n_entities; i++)
      s_val[i] = (val[i] - v_min)*scale;

    int *dest_rank = NULL;
    cs_lnum_t *order = NULL;

    BFT_MALLOC(order, n_entities, cs_lnum_t);
    BFT_MALLOC(dest_rank, n_entities, int);

    cs_order_real_allocated(NULL, s_val, order, n_entities);

    cs_sort_partition_dest_rank_id(_sampling_factors[1],
                                   sizeof(cs_real_t),
                                   n_entities,
                                   s_val,
                                   NULL, /* weight */
                                   order,
                                   dest_rank,
                                   _s_to_real,
                                   _s_compare,
                                   NULL,
                                   comm);

    BFT_FREE(order);

    cs_all_to_all_t
      *d = cs_all_to_all_create(this_io_num->global_num_size,
                                0,     /* flags */
                                NULL,  /* dest_id */
                                dest_rank,
                                comm);

    cs_all_to_all_transfer_dest_rank(d, &dest_rank);

    cs_real_t *b_val
      = cs_all_to_all_copy_array(d,
                                 CS_REAL_TYPE,
                                 1,
                                 false, /* reverse */
                                 s_val,
                                 NULL);

    BFT_FREE(s_val);

    size_t b_size = cs_all_to_all_n_elts_dest(d);

    /* Now re-build order on block distribution. */

    BFT_MALLOC(order, b_size, cs_lnum_t);

    cs_order_real_allocated(NULL, b_val, order, b_size);

    BFT_FREE(b_val);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (blocks associated with each process are
       already sorted, but the whole "gathered" block is not).
       We build an initial global order based on the initial global numbering,
       such that for each block, the global number of an entity is equal to
       the cumulative number of sub-entities */

    cs_gnum_t *b_gnum;
    BFT_MALLOC(b_gnum, b_size, cs_gnum_t);

    for (i = 0; i < b_size; i++)
      b_gnum[order[i]] = i+1;

    BFT_FREE(order);

    cs_gnum_t gnum_shift = 0, current_gnum = b_size;

    /* At this stage, b_gnum[] is valid for this process, and
       current_gnum indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    MPI_Scan(&current_gnum, &gnum_shift, 1, CS_MPI_GNUM,
             MPI_SUM, comm);
    gnum_shift -= current_gnum;

    for (i = 0; i < b_size; i++)
      b_gnum[i] += gnum_shift;

    /* Return global order to all ranks */

    cs_all_to_all_copy_array(d,
                             CS_GNUM_TYPE,
                             1,
                             true, /* reverse */
                             b_gnum,
                             this_io_num->_global_num);

    /* Free memory */

    BFT_FREE(b_gnum);

    cs_all_to_all_destroy(&d);

    /* Get final maximum global number value */

    this_io_num->global_count
      = _fvm_io_num_global_max_unordered(this_io_num, comm);

  }

#endif /* HAVE_MPI */

  if (n_ranks == 1) {

    cs_lnum_t *order = NULL;
    BFT_MALLOC(order, n_entities, cs_lnum_t);

    cs_order_real_allocated(NULL, val, order, n_entities);

    for (i = 0; i < n_entities; i++)
      this_io_num->_global_num[order[i]] = i+1;

    BFT_FREE(order);

    this_io_num->global_count = n_entities;

  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a simple accumulation
 * (i.e. scan) of counts on successive ranks.
 *
 * parameters:
 *   n_entities <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_scan(size_t  n_entities)
{
  fvm_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (cs_glob_n_ranks < 2)
    return NULL;

#if defined(HAVE_MPI)
  {
    size_t  i;
    cs_gnum_t gnum_base = n_entities;
    cs_gnum_t gnum_sum = n_entities;
    cs_gnum_t gnum_shift = 0;

    MPI_Comm comm = cs_glob_mpi_comm;

    /* Create structure */

    BFT_MALLOC(this_io_num, 1, fvm_io_num_t);

    BFT_MALLOC(this_io_num->_global_num, n_entities, cs_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    this_io_num->global_num_size = n_entities;

    MPI_Scan(&gnum_base, &gnum_shift, 1, CS_MPI_GNUM, MPI_SUM, comm);

    gnum_base = gnum_shift - gnum_base + 1;

    for (i = 0; i < n_entities; i++)
      this_io_num->_global_num[i] = gnum_base + i;

    gnum_base = n_entities;

    MPI_Allreduce(&gnum_base, &gnum_sum, 1, CS_MPI_GNUM, MPI_SUM, comm);

    this_io_num->global_count = gnum_sum;
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Destruction of a I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_destroy(fvm_io_num_t  * this_io_num)
{
  if (this_io_num != NULL) {
    BFT_FREE(this_io_num->_global_num);
    BFT_FREE(this_io_num);
  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Transfer ownership of global numbering array from IO numbering structure.
 *
 * parameters:
 *   this_io_num <-> pointer to structure transferring array ownership.
 *
 * returns:
 *   pointer to transferred array
 *----------------------------------------------------------------------------*/

cs_gnum_t *
fvm_io_num_transfer_global_num(fvm_io_num_t  * this_io_num)
{
  cs_gnum_t *retval = NULL;

  if (this_io_num != NULL) {
    retval = this_io_num->_global_num;
    this_io_num->_global_num = NULL;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return local number of entities associated with an I/O numbering
 * structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  local number of associated entities
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_io_num_get_local_count(const fvm_io_num_t  *const this_io_num)
{
  assert(this_io_num != NULL);

  return this_io_num->global_num_size;
}

/*----------------------------------------------------------------------------
 * Return global number of entities associated with an I/O numbering
 * structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  global number of associated entities
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_io_num_get_global_count(const fvm_io_num_t  *const this_io_num)
{
  assert(this_io_num != NULL);

  return this_io_num->global_count;
}

/*----------------------------------------------------------------------------
 * Return global numbering associated with an I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  pointer to array of global numbers associated with local entities
 *  (1 to n numbering)
 *----------------------------------------------------------------------------*/

const cs_gnum_t *
fvm_io_num_get_global_num(const fvm_io_num_t  *const this_io_num)
{
  assert(this_io_num != NULL);

  return this_io_num->global_num;
}

/*----------------------------------------------------------------------------
 * Return the global number of sub-entities associated with an initial
 * entity whose global numbering is known, given the number of
 * sub-entities per initial entity.
 *
 * parameters:
 *   this_io_num    <-> pointer to base io numbering
 *   n_sub_entities <-- number of sub-entities per initial entity
 *   comm           <-- associated MPI communicator
 *
 * returns:
 *   global number of sub-entities
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_io_num_global_sub_size(const fvm_io_num_t  *this_io_num,
                           const cs_lnum_t      n_sub_entities[])
{
  cs_gnum_t   retval = 0;

  /* Initial checks */

  if (this_io_num == NULL)
    return retval;

#if defined(HAVE_MPI)
 if (cs_glob_n_ranks > 1) {
   int  have_sub_loc = 0, have_sub_glob = 0;

   /* Caution: we may have sub-entities on some ranks and not on others */

   if (n_sub_entities != NULL)
     have_sub_loc = 1;

   MPI_Allreduce(&have_sub_loc, &have_sub_glob, 1, MPI_INT, MPI_MAX,
                 cs_glob_mpi_comm);

   if (have_sub_glob > 0)
     retval = _fvm_io_num_global_sub_size(this_io_num,
                                          n_sub_entities,
                                          cs_glob_mpi_comm);
 }
#endif

 if (cs_glob_n_ranks == 1 && n_sub_entities != NULL) {
   for (size_t i = 0; i < (size_t)(this_io_num->global_num_size); i++)
     retval += n_sub_entities[i];
 }

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump printout of a I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_io_num_dump(const fvm_io_num_t  *const this_io_num)
{
  cs_lnum_t i;

  if (this_io_num == NULL) {
    bft_printf("  global numbering: nil\n");
    return;
  }

  bft_printf("  global numbering size:            %u\n",
             (unsigned)this_io_num->global_num_size);

  bft_printf("\n"
             "  pointer to shareable array:\n"
             "    global_num:                     %p\n",
             (const void *)this_io_num->global_num);

  bft_printf("\n"
             "  pointer to local array:\n"
             "    _global_num:                    %p\n",
             (const void *)this_io_num->global_num);

  if (this_io_num->global_num_size > 0) {

    bft_printf("\n  global number:\n\n");
    for (i = 0 ; i < this_io_num->global_num_size ; i++)
      bft_printf("  %10u : %10llu\n",
                 (unsigned)i + 1,
                 (unsigned long long)this_io_num->global_num[i]);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
