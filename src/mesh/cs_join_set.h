#ifndef __CS_JOIN_SET_H__
#define __CS_JOIN_SET_H__

/*============================================================================
 * Subroutines useful to manage list structures
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include  <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

typedef struct { /* Definition of a global indexed list of global elements */

  cs_lnum_t   n_elts;
  cs_gnum_t   n_g_elts;

  cs_gnum_t  *g_elts;   /* Global numbering of elements */

  cs_lnum_t  *index;    /* Index on elements from */
  cs_gnum_t  *g_list;   /* Global numbering of entities linked with g_elts */

} cs_join_gset_t;

typedef struct { /* Resizable array structure */

  cs_lnum_t    n_max_elts;
  cs_lnum_t    n_elts;
  cs_lnum_t   *array;

} cs_join_rset_t;

/* ------------------------------------------------------------------ *
 * Definition of a structure defining a set of equivalence between
 * vertices for instance
 * ------------------------------------------------------------------ */

typedef struct {

  cs_lnum_t   n_max_equiv;    /* max. number of equiv. allocated */
  cs_lnum_t   n_equiv;        /* number of equivalences */
  cs_lnum_t  *equiv_couple;   /* ids of the two equivalent entities.
                                 size = 2 * n_equiv */
} cs_join_eset_t;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Allocate a resizable array.
 *
 * parameters:
 *   max_size <-- initial number of elements to allocate
 *
 * returns:
 *   pointer to a new alloacted resizable array
 *---------------------------------------------------------------------------*/

cs_join_rset_t *
cs_join_rset_create(cs_lnum_t  max_size);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_rset_t structure.
 *
 * parameter:
 *   set <-- pointer to pointer to the cs_join_rset_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_rset_destroy(cs_join_rset_t  **set);

/*----------------------------------------------------------------------------
 * Check if we need to resize the current cs_join_rset_t structure and do
 * it if necessary.
 *
 * parameters:
 *   set       <-- pointer to pointer to the cs_join_rset_t structure to test
 *   test_size <-- target size
 *---------------------------------------------------------------------------*/

void
cs_join_rset_resize(cs_join_rset_t  **set,
                    cs_lnum_t         test_size);

/*----------------------------------------------------------------------------
 * Create a new cs_join_eset_t structure.
 *
 * parameters:
 *   init_size <-- number of initial equivalences to allocate
 *
 * returns:
 *   a pointer to a new cs_join_eset_t structure
 *---------------------------------------------------------------------------*/

cs_join_eset_t *
cs_join_eset_create(cs_lnum_t  init_size);

/*----------------------------------------------------------------------------
 * Check if the requested size if allocated in the structure.
 *
 * Reallocate cs_join_eset_t structure if necessary.
 *
 * parameters:
 *   request_size <-- necessary size
 *   equiv_set    <-> pointer to pointer to the cs_join_eset_t struct.
 *---------------------------------------------------------------------------*/

void
cs_join_eset_check_size(cs_lnum_t         request_size,
                        cs_join_eset_t  **equiv_set);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_eset_t structure.
 *
 * parameter:
 *   equiv_set <-- pointer to pointer to the structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_eset_destroy(cs_join_eset_t  **equiv_set);

/*----------------------------------------------------------------------------
 * Clean a cs_join_eset_t structure.
 *
 * If necessary, create a new cs_join_eset_t structure with no redundancy.
 *
 * parameters:
 *   eset <-- pointer to pointer to the cs_join_eset_t structure to clean
 *---------------------------------------------------------------------------*/

void
cs_join_eset_clean(cs_join_eset_t  **eset);

/*----------------------------------------------------------------------------
 * Create a cs_join_gset_t structure (indexed list on global numbering)
 *
 * parameters:
 *   n_elts <-- number of elements composing the list
 *
 * returns:
 *   a new allocated pointer to a cs_join_gset_t structure.
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_create(cs_lnum_t  n_elts);

/*----------------------------------------------------------------------------
 * Build a cs_join_gset_t structure to store all the potential groups
 * between elements.
 *
 * Values in g_elts are the tag values and values in g_list
 * are position in tag array.
 *
 * parameters:
 *   n_elts <-- number of elements in tag array
 *   tag    <-- tag array used to define a new cs_join_gset_t
 *
 * returns:
 *   a new allocated cs_join_gset_t structure
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_create_from_tag(cs_lnum_t        n_elts,
                             const cs_gnum_t  tag[]);

/*----------------------------------------------------------------------------
 * Create a new cs_join_gset_t which holds equivalences between elements of
 * g_list in cs_join_gset_t.
 *
 * For a subset of equivalences, we store their initial value in the return
 * cs_join_gset_t structure. A subset is defined if at least two elements
 * are equivalent.
 *
 * The behavior of this function is near from cs_join_gset_create_from_tag
 * but we don't store the position in init_array but its value in init_array.
 *
 * parameters:
 *   set        <-- pointer to a cs_join_gset_t structure
 *   init_array <-- initial values of set->g_list
 *
 * returns:
 *   a new allocated cs_join_gset_t structure, or NULL if it would be empty
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_create_by_equiv(const cs_join_gset_t  *set,
                             const cs_gnum_t        init_array[]);

/*----------------------------------------------------------------------------
 * Copy a cs_join_gset_t structure.
 *
 * parameters:
 *   src <-- pointer to the cs_join_gset_t structure to copy
 *
 * returns:
 *   a new allocated cs_join_gset_t structure.
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_copy(const cs_join_gset_t  *src);

/*----------------------------------------------------------------------------
 * Destroy a cs_join_gset_t structure.
 *
 * parameters:
 *   set <-- pointer to pointer to the cs_join_gset_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_gset_destroy(cs_join_gset_t  **set);

/*----------------------------------------------------------------------------
 * Sort a cs_join_gset_t structure according to the global numbering of
 * the g_elts in cs_join_gset_t structure.
 *
 * parameters:
 *   set <-> pointer to the structure to order
 *---------------------------------------------------------------------------*/

void
cs_join_gset_sort_elts(cs_join_gset_t  *set);

/*----------------------------------------------------------------------------
 * Sort each sub-list of the g_list array in a cs_join_gset_t structure.
 *
 * parameters:
 *   p_set <-> pointer to the structure to sort
 *---------------------------------------------------------------------------*/

void
cs_join_gset_sort_sublist(cs_join_gset_t  *set);

/*----------------------------------------------------------------------------
 * Invert a cs_join_gset_t structure.
 *
 * parameters:
 *   set <-- pointer to the cs_join_gset_t structure to work with
 *
 * returns:
 *   the new allocated and inverted set structure
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_invert(const cs_join_gset_t  *set);

/*----------------------------------------------------------------------------
 * Delete redudancies in a cs_join_gset_t structure.
 *
 * The output set has an ordered sub-list for each element in the set.
 *
 * parameters:
 *   set <-> pointer to the structure to clean
 *---------------------------------------------------------------------------*/

void
cs_join_gset_clean(cs_join_gset_t  *set);

/*----------------------------------------------------------------------------
 * Delete redudancies in g_list array of a cs_join_gset_t structure.
 *
 * parameters:
 *   set          <-> pointer to the structure to clean
 *   linked_array <-> array for which redundancies are scanned
 *---------------------------------------------------------------------------*/

void
cs_join_gset_clean_from_array(cs_join_gset_t  *set,
                              cs_gnum_t        linked_array[]);

/*----------------------------------------------------------------------------
 * Concatenate the two g_elts and g_list arrays.
 *
 * Order the new concatenated array and delete redundant elements.
 * We get a single ordered array.
 *
 * parameters:
 *   set       <-- pointer to the structure to work with
 *   n_elts    --> number of elements in the new set
 *   new_array --> pointer to the new created array
 *---------------------------------------------------------------------------*/

void
cs_join_gset_single_order(const cs_join_gset_t  *set,
                          cs_lnum_t             *n_elts,
                          cs_gnum_t             *new_array[]);

/*----------------------------------------------------------------------------
 * Compress a g_list such as for each element "e" in g_elts:
 *  - there is no redundancy for the linked elements of set->g_list
 *  - there is no element in set->g_list < e except if this element is not
 *    present in g_elts
 *
 * g_list and g_elts need to be ordered before calling this function.
 *
 * parameters:
 *   set <-> pointer to the structure to work with
 *---------------------------------------------------------------------------*/

void
cs_join_gset_compress(cs_join_gset_t  *set);

/*----------------------------------------------------------------------------
 * Delete redundancies in set->g_elts.
 *
 * Merge sub-arrays associated to a common set->g_elts[i].
 *
 * parameters:
 *   set       <-- pointer to the structure to work with
 *   order_tag <-- 0: set->g_elts is not ordered, 1: ordered
 *---------------------------------------------------------------------------*/

void
cs_join_gset_merge_elts(cs_join_gset_t  *set,
                        int              order_tag);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Synchronize a cs_join_gset_t structure and distribute the resulting set
 * over the rank thanks to a round-robin distribution. Elements in sync_set
 * are ordered and there is no redundancy but list may have redundancies.
 * Use cs_join_gset_clean() to remove redundancies in g_list.
 *
 * parameters:
 *   loc_set  <-> pointer to the local structure to work with
 *   comm     <-- mpi_comm on which synchro. and distribution take place
 *
 * returns:
 *   a synchronized and distributed cs_join_gset_t structure.
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_robin_sync(cs_join_gset_t   *loc_set,
                        MPI_Comm          comm);

/*----------------------------------------------------------------------------
 * Update a local cs_join_gset_t structure from a distributed and
 * synchronized cs_join_gset_t structure. Round-robin distribution is used
 * to store synchronized elements.
 *
 * parameters:
 *   sync_set <-- pointer to the structure which holds a synchronized block
 *   loc_set  <-> pointer to a local structure holding elements to update
 *   comm     <-- comm on which synchronization and distribution take place
 *---------------------------------------------------------------------------*/

void
cs_join_gset_robin_update(const cs_join_gset_t   *sync_set,
                          cs_join_gset_t         *loc_set,
                          MPI_Comm                comm);

/*----------------------------------------------------------------------------
 * Synchronize a cs_join_gset_t structure and distribute the resulting set
 * over the rank by block
 *
 * parameters:
 *   max_gnum <-- max global number in global element numbering
 *   loc_set  <-> pointer to the local structure to work with
 *   comm     <-- mpi_comm on which synchro. and distribution take place
 *
 * returns:
 *   a synchronized and distributed cs_join_gset_t structure.
 *---------------------------------------------------------------------------*/

cs_join_gset_t *
cs_join_gset_block_sync(cs_gnum_t        max_gnum,
                        cs_join_gset_t  *loc_set,
                        MPI_Comm         comm);

/*----------------------------------------------------------------------------
 * Update a local cs_join_gset_t structure from a distributed and
 * synchronized cs_join_gset_t structure.
 *
 * loc_set should not have redundant elements.
 *
 * parameters:
 *   max_gnum <-- max global number in global element numbering
 *   sync_set <-- pointer to the structure which holds a synchronized block
 *   loc_set  <-> pointer to a local structure holding elements to update
 *   comm     <-- comm on which synchronization and distribution take place
 *---------------------------------------------------------------------------*/

void
cs_join_gset_block_update(cs_gnum_t              max_gnum,
                          const cs_join_gset_t  *sync_set,
                          cs_join_gset_t        *loc_set,
                          MPI_Comm               comm);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Dump an array (int or double).
 *
 * This function is called according to the verbosity.
 *
 * parameters:
 *   f       <-- handle to output file
 *   type    <-- type of the array to display
 *   header  <-- header to display in front of the array
 *   n_elts  <-- number of elements to display
 *   array   <-- array to display
 *---------------------------------------------------------------------------*/

void
cs_join_dump_array(FILE        *f,
                   const char  *type,
                   const char  *header,
                   int          n_elts,
                   const void  *array);

/*----------------------------------------------------------------------------
 * Dump a cs_join_gset_t structure.
 *
 * parameters:
 *   f    <-- handle to output file
 *   set  <-- pointer to the cs_join_gset_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_gset_dump(FILE                  *f,
                  const cs_join_gset_t  *set);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_SET_H__ */
