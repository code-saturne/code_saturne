/*============================================================================
 * Definition of entity groups
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_group.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro and Type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

#define _GROUP_TAG      (int)('G'+'R'+'O'+'U'+'P') /* MPI tag */

#endif

/*----------------------------------------------------------------------------
 * Opaque structure describing a group class
 *----------------------------------------------------------------------------*/

struct _fvm_group_class_t {
  int              n_groups;            /* Number of groups in class */
  char           **group_name;          /* Array of group names */
};

/*----------------------------------------------------------------------------
 * Opaque structure defining a set of group classes
 *----------------------------------------------------------------------------*/

struct _fvm_group_class_set_t {

  int size;                             /* Number of group classes */

  fvm_group_class_t   *class;           /* Array of group classes */
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compare strings (qsort function).
 *
 * parameters:
 *   a <-> pointer to first string
 *   b <-> pointer to second string
 *
 * returns:
 *   result of strcmp() on strings
 *----------------------------------------------------------------------------*/

static int
_compare_names(const void  *a,
               const void  *b)
{
  return strcmp(*((const char *const *)a), *((const char *const *)b));
}

#if 0 && defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Send group class set definition to distant rank.
 *
 * parameters:
 *   class_set <-- pointer to group class set structure
 *   dest_rank <-- destination rank
 *   comm      <-- associated mpi communicator
 *----------------------------------------------------------------------------*/

static void
_group_class_set_send(const fvm_group_class_set_t  *class_set,
                      int                           dest_rank,
                      MPI_Comm                      comm)
{
  int i, j;

  int n_ints = 0;
  int n_chars = 0;

  int send_count[3];
  int *send_ints = NULL;
  char *send_chars = NULL;

  /* Counting pass */

  for (i = 0; i < class_set->size; i++) {
    const fvm_group_class_t  *gc = class_set->class + i;
    n_ints += 2 + gc->n_groups;
    for (j = 0; j < gc->n_groups; j++)
      n_chars += strlen(gc->group_name[j]) + 1;
  }

  /* Allocate and prepare buffers */

  BFT_MALLOC(send_ints, n_ints, int);
  BFT_MALLOC(send_chars, n_chars, char);

  send_count[0] = class_set->size;
  send_count[1] = n_ints;
  send_count[2] = n_chars;

  n_ints = 0;
  n_chars = 0;

  for (i = 0; i < class_set->size; i++) {
    const fvm_group_class_t  *gc = class_set->class + i;
    send_ints[n_ints++] = gc->n_groups;
    for (j = 0; j < gc->n_groups; j++) {
      strcpy(send_chars + n_chars, gc->group_name[j]);
      n_chars += strlen(gc->group_name[j]);
      send_chars[n_chars] = '\0';
      n_chars++;
    }
  }

  assert(n_ints == send_count[1]);
  assert(n_chars == send_count[2]);

  MPI_Send(send_count, 3, MPI_INT, dest_rank, _GROUP_TAG, comm);

  if (n_ints > 0)
    MPI_Send(send_ints, n_ints, MPI_INT, dest_rank, _GROUP_TAG, comm);
  if (n_chars > 0)
    MPI_Send(send_chars, n_chars, MPI_CHAR, dest_rank, _GROUP_TAG, comm);

  BFT_FREE(send_ints);
  BFT_FREE(send_chars);
}

/*----------------------------------------------------------------------------
 * Receive group class set definition from distant rank.
 *
 * parameters:
 *   class_set --> pointer to initially empty class set structure
 *   src_rank  <-- source rank
 *   comm      <-- associated mpi communicator
 *----------------------------------------------------------------------------*/

static void
_group_class_set_recv(fvm_group_class_set_t  *class_set,
                      int                     src_rank,
                      MPI_Comm                comm)
{
  int i, j;

  MPI_Status status;

  int n_ints = 0;
  int n_chars = 0;

  int recv_count[3];
  int *recv_ints = NULL;
  char *recv_chars = NULL;

  assert(class_set != NULL);
  assert(class_set->size == 0);

  MPI_Recv(recv_count, 3, MPI_INT, src_rank, _GROUP_TAG, comm, &status);

  /* Allocate and receive buffers */

  if (recv_count[1] > 0) {
    BFT_MALLOC(recv_ints, recv_count[1], int);
    MPI_Recv(recv_ints, recv_count[1], MPI_INT,
             src_rank, _GROUP_TAG, comm, &status);
  }

  if (recv_count[2] > 0) {
    BFT_MALLOC(recv_chars, recv_count[2], char);
    MPI_Recv(recv_chars, recv_count[2], MPI_CHAR,
             src_rank, _GROUP_TAG, comm, &status);
  }

  /* Decode buffers */

  class_set->size = recv_count[0];
  BFT_MALLOC(class_set->class, class_set->size, fvm_group_class_t);

  for (i = 0; i < class_set->size; i++) {
    fvm_group_class_t  *gc = class_set->class + i;
    gc->n_groups = recv_ints[n_ints++];
    if (gc->n_groups > 0)
      BFT_MALLOC(gc->group_name, gc->n_groups, char *);
    for (j = 0; j < gc->n_groups; j++) {
      size_t l = strlen(recv_chars + n_chars) + 1;
      BFT_MALLOC(gc->group_name[j], l, char);
      strcpy(recv_chars + n_chars, gc->group_name[j]);
      n_chars += l;
    }
  }

  /* Free receive buffers */

  BFT_FREE(recv_ints);
  BFT_FREE(recv_chars);
}

#endif

/*----------------------------------------------------------------------------
 * Copy a group class
 *
 * parameters:
 *   src  <-- pointer to source group class
 *   dest --> pointer to group class copy
 *----------------------------------------------------------------------------*/

static void
_group_class_copy(const fvm_group_class_t  *src,
                  fvm_group_class_t        *dest)
{
  int i;

  if (src == NULL) {
    dest->n_groups = 0;
    dest->group_name = NULL;
  }
  else {
    dest->n_groups = src->n_groups;
    BFT_MALLOC(dest->group_name, dest->n_groups, char *);
    for (i = 0; i < src->n_groups; i++) {
      BFT_MALLOC(dest->group_name[i], strlen(src->group_name[i]) + 1, char);
      strcpy(dest->group_name[i], src->group_name[i]);
    }
  }
}

/*----------------------------------------------------------------------------
 * Dump printout of a group class
 *
 * parameters:
 *   group_class         <-- pointer to group class to be dumped
 *   id                  <-- index of group class in set
 *----------------------------------------------------------------------------*/

static void
_group_class_dump(const fvm_group_class_t  *this_group_class,
                  int                       id)
{
  int i;

  if (this_group_class == NULL) {
    bft_printf("\n    _group_class[%d]: nil\n", id);
    return;
  }

  bft_printf("\n"
             "    _group_class[%3d]: %p\n"
             "    n_groups:          %d\n",
             id, (const void *)this_group_class,
             this_group_class->n_groups);

  if (this_group_class->n_groups > 0) {
    bft_printf("    group names:\n");
    for (i = 0; i < this_group_class->n_groups; i++)
      bft_printf("     \" %s\"\n", this_group_class->group_name[i]);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the number of groups of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   number of groups in group class
 *----------------------------------------------------------------------------*/

int
fvm_group_class_get_n_groups(const fvm_group_class_t  *this_group_class)
{
  int retval = 0;

  if (this_group_class != NULL) {
    retval = this_group_class->n_groups;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the array of group names of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   pointer to array of group names in group class
 *----------------------------------------------------------------------------*/

const char **
fvm_group_class_get_group_names(const fvm_group_class_t  *this_group_class)
{
  const char **retval = NULL;

  if (this_group_class != NULL) {
    retval = (const char **)(this_group_class->group_name);
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Creation of a group class set structure.
 *
 * returns:
 *   pointer to created group class set structure
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
fvm_group_class_set_create(void)
{
  fvm_group_class_set_t *class_set;

  BFT_MALLOC(class_set, 1, fvm_group_class_set_t);

  class_set->size = 0;

  class_set->class = NULL;

  return class_set;
}

/*----------------------------------------------------------------------------
 * Add a group class to a set.
 *
 * Group names are automatically sorted in the group class description.
 *
 * parameters:
 *   this_group_class_set <-> pointer to group class set structure
 *   n_groups             <-- number of groups in class
 *   group_names          <-- array of group names
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_add(fvm_group_class_set_t   *this_group_class_set,
                        int                      n_groups,
                        const char             **group_names)
{
  int i;
  fvm_group_class_set_t *class_set = this_group_class_set;
  fvm_group_class_t *_class = NULL;

  assert(class_set != NULL);

  /* Resize array of group class descriptors */

  BFT_REALLOC(class_set->class, class_set->size + 1, fvm_group_class_t);

  /* Initialize new descriptor */

  _class = class_set->class + class_set->size;

  _class->n_groups = n_groups;

  BFT_MALLOC(_class->group_name, n_groups, char *);

  /* Copy group names */

  if (n_groups > 0) {
    for (i = 0; i < n_groups; i++) {
      BFT_MALLOC(_class->group_name[i], strlen(group_names[i]) + 1, char);
      strcpy(_class->group_name[i], group_names[i]);
    }
    qsort(_class->group_name, n_groups, sizeof(char *), &_compare_names);
  }

  /* Update the number of classes in set */

  class_set->size += 1;
}

/*----------------------------------------------------------------------------
 * Destruction of a group class set structure.
 *
 * parameters:
 *   this_class_set <-- pointer to structure which should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
fvm_group_class_set_destroy(fvm_group_class_set_t  *this_group_class_set)
{
  if (this_group_class_set == NULL)
    return NULL;

  for (int i = 0; i < this_group_class_set->size; i++) {

    fvm_group_class_t *_class = this_group_class_set->class + i;

    for (int j = 0; j < _class->n_groups; j++)
      BFT_FREE(_class->group_name[j]);

    _class->n_groups = 0;

    BFT_FREE(_class->group_name);

  }

  BFT_FREE(this_group_class_set->class);
  BFT_FREE(this_group_class_set);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Return number of classes in a group class set.
 *
 * parameters:
 *   this_group_class_set <-> pointer to group class set structure
 *
 * returns:
 *   number of classes in a group class set
 *----------------------------------------------------------------------------*/

int
fvm_group_class_set_size(const fvm_group_class_set_t  *this_group_class_set)
{
  int retval = 0;

  if (this_group_class_set != NULL)
    retval = this_group_class_set->size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to a given group class in a group class set.
 *
 * parameters:
 *   this_group_class_set <-- pointer to group class set structure
 *   group_class_id       <-- index of group class in set (0 to n-1)
 *
 * returns:
 *   pointer to group class structure
 *----------------------------------------------------------------------------*/

const fvm_group_class_t *
fvm_group_class_set_get(const fvm_group_class_set_t  *this_group_class_set,
                        int                           group_class_id)
{
  const fvm_group_class_t  *retval = NULL;

  if (this_group_class_set != NULL) {
    if (group_class_id > -1 && group_class_id < this_group_class_set->size)
      retval = this_group_class_set->class + group_class_id;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Copy a group class set, or a subset thereof
 *
 * parameters:
 *   this_class_set <-- pointer to group class set to be copied
 *   n_gcs          <-- number of group classes
 *   gc_id          <-- group class list (0 to n-1), or NULL
 *
 * returns:
 *   pointer to new copy of group class set
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t  *
fvm_group_class_set_copy(const fvm_group_class_set_t  *this_group_class_set,
                         int                           n_gcs,
                         int                           gc_id[])
{
  int i;
  fvm_group_class_set_t  *class_set = NULL;

  BFT_MALLOC(class_set, 1, fvm_group_class_set_t);

  if (n_gcs == 0)
    class_set->size = this_group_class_set->size;
  else
    class_set->size = n_gcs;

  BFT_MALLOC(class_set->class, class_set->size, fvm_group_class_t);

  if (n_gcs == 0) {
    for (i = 0; i < class_set->size; i++)
      _group_class_copy(this_group_class_set->class + i,
                        class_set->class + i);
  }
  else {
    for (i = 0; i < n_gcs; i++)
      _group_class_copy(this_group_class_set->class + gc_id[i],
                        class_set->class + i);
  }

  return class_set;
}

/*----------------------------------------------------------------------------
 * Dump printout of a group class set
 *
 * parameters:
 *   this_class_set      <-- pointer to group class set to be dumped
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_dump(const fvm_group_class_set_t  *this_group_class_set)
{
  int i;
  const fvm_group_class_set_t  *class_set = this_group_class_set;

  if (this_group_class_set == NULL) {
    bft_printf("  group_class_set: nil\n");
    return;
  }

  bft_printf("  _group_class_set: %p\n"
             "  size:             %d\n",
             (const void *)class_set,
             class_set->size);

  if (class_set->size > 0) {
    bft_printf("\n  group_classes:");
    for (i = 0; i < class_set->size; i++)
      _group_class_dump(class_set->class + i, i);
  }

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
