/*============================================================================
 * Mechanism for entity selection based on groups or attributes
 *============================================================================*/

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_timer.h"

#include "fvm_defs.h"
#include "fvm_group.h"
#include "fvm_selector_postfix.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_selector.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local structure types
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Description (caching) of all previously interpreted criteria
 *----------------------------------------------------------------------------*/

typedef struct {

  int     n_operations;                /* Number of previously interpreted
                                          operations */
  int     n_max_operations;            /* Maximum size before reallocation */

  fvm_selector_postfix_t **postfix;    /* Array of postfix operations */

  size_t *n_calls;                     /* Number of calls per operation */

  int    *n_group_classes;             /* Array of group class numbers
                                          for each operation */
  int   **group_class_set;             /* Array of group class lists
                                          for each operation */
} _operation_list_t;

/*----------------------------------------------------------------------------
 * Opaque type for management of all necessary info for elements selection
 *----------------------------------------------------------------------------*/

struct _fvm_selector_t {

  int                dim;                      /* Spatial dimension */
  cs_lnum_t          n_elements;               /* Number of elements */

  const int         *group_class_id;           /* Element group class ids */
  int               *_group_class_id;          /* private group_class_id,
                                                  or NULL */
  int                group_class_id_base;      /* Starting base (usually
                                                  0 or 1) of group class ids */

  int                n_group_classes;          /* Number of group classes */

  int                n_groups;                 /* Total number of groups */
  int                n_attributes;             /* Total number of attributes */

  char             **group_name;               /* Ordered group names */
  int               *attribute;                /* Ordered attributes */

  int               *n_class_groups;           /* Number of groups per class */
  int              **group_ids;                /* Id of groups per class in
                                                  group_name */

  int               *n_class_attributes;       /* Number of attrs. per class */
  int              **attribute_ids;            /* Id of attributes per class in
                                                  attribute */

  const double      *coords;                   /* Element coordinates
                                                  (i.e. centers), interlaced */
  double            *_coords;                  /* private coords, or NULL */

  const double      *normals;                  /* Element normals, interlaced */
  double            *_normals;                 /* private normals, or NULL */

  _operation_list_t *_operations;              /* Array which caches all
                                                  previously interpreted
                                                  strings (operations) */

  cs_lnum_t         *_n_group_class_elements;  /* Number of elements per
                                                  group class */
  cs_lnum_t         **_group_class_elements;   /* Group class elements array */

  int                 n_evals;                 /* Number of evaluations */
  double              eval_wtime;              /* Wall-clock time
                                                  evaluating selections */
};

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compare groups (qsort wrapper to strcmp).
 *
 * parameters:
 *   x <-> pointer to first group name
 *   y <-> pointer to second group name
 *
 * returns:
 *   result of strcmp() on group names
 *----------------------------------------------------------------------------*/

static int _compare_groups(const void *x, const void *y)
{
  return strcmp(*(char * const *)x, *(char * const *)y);
}

/*----------------------------------------------------------------------------
 * Add group information to a selector.
 *
 * A global list of group names from all group classes is built and sorted,
 * and the id's of groups of each class are determined.
 *
 * parameters:
 *   this_selector   <-> pointer to selector
 *   group_class_set <-- group classes description
 *----------------------------------------------------------------------------*/

static void
_assign_groups(fvm_selector_t               *this_selector,
               const fvm_group_class_set_t  *group_class_set)
{
  int i, j;
  char **_set_group_names = NULL;
  const char **set_group_names = NULL;
  int n_groups_tot = 0;

  /* Build basic arrays which exist in any case */

  BFT_MALLOC(this_selector->n_class_groups,
             this_selector->n_group_classes,
             int);
  BFT_MALLOC(this_selector->group_ids, this_selector->n_group_classes, int *);

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvm_group_class_t *gc = fvm_group_class_set_get(group_class_set, i);
    const int n_groups = fvm_group_class_get_n_groups(gc);

    n_groups_tot += n_groups;

    this_selector->n_class_groups[i] = n_groups;
    if (n_groups > 0)
      BFT_MALLOC(this_selector->group_ids[i], n_groups, int);
    else
      this_selector->group_ids[i] = NULL;
  }

  if (n_groups_tot == 0)
    return;

  /* Fill arrays with unsorted group class info */

  BFT_MALLOC(_set_group_names, n_groups_tot, char *);
  set_group_names = (const char **)_set_group_names;

  n_groups_tot = 0;

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvm_group_class_t *gc = fvm_group_class_set_get(group_class_set, i);
    const int n_groups = fvm_group_class_get_n_groups(gc);
    const char **group_names = fvm_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++)
      set_group_names[n_groups_tot + j] = group_names[j];

    n_groups_tot += n_groups;
  }

  qsort(set_group_names, n_groups_tot, sizeof(char *), &_compare_groups);

  BFT_MALLOC(this_selector->group_name, n_groups_tot, char *);

  BFT_MALLOC(this_selector->group_name[0],
             strlen(set_group_names[0]) + 1,
             char);
  strcpy(this_selector->group_name[0], set_group_names[0]);
  for (i = 1, j = 1; i < n_groups_tot; i++) {
    const char *name = set_group_names[i];
    if (strcmp(name, set_group_names[i-1]) != 0) {
      BFT_MALLOC(this_selector->group_name[j], strlen(name) + 1, char);
      strcpy(this_selector->group_name[j], name);
      j++;
    }
  }

  set_group_names = NULL;
  BFT_FREE(_set_group_names);

  this_selector->n_groups = j;
  BFT_REALLOC(this_selector->group_name, this_selector->n_groups, char *);

  /* Now build group_id arrays */

  for (i = 0; i < this_selector->n_group_classes; i++) {

    int mid_id;
    const char *name;
    const fvm_group_class_t *gc = fvm_group_class_set_get(group_class_set, i);
    const int n_groups = fvm_group_class_get_n_groups(gc);
    const char **group_names = fvm_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++) {

      /* use binary search */

      int start_id = 0;
      int end_id = this_selector->n_groups - 1;
      mid_id = (end_id -start_id) / 2;
      name = group_names[j];

      while (start_id <= end_id) {
        int att_cmp = strcmp(this_selector->group_name[mid_id], name);
        if (att_cmp < 0)
          start_id = mid_id + 1;
        else if (att_cmp > 0)
          end_id = mid_id - 1;
        else
          break;
        mid_id = start_id + ((end_id -start_id) / 2);
      }

      assert(strcmp(this_selector->group_name[mid_id], name) == 0);
      this_selector->group_ids[i][j] = mid_id;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compare attributes (for qsort).
 *
 * parameters:
 *   x <-> pointer to first attribute
 *   y <-> pointer to second attribute
 *
 * returns:
 *   result of strcmp() on group names
 *----------------------------------------------------------------------------*/

static int _compare_attributes(const void *x, const void *y)
{
  return (*(const int *)x - *(const int *)y);
}

/*----------------------------------------------------------------------------
 * Check if a string defines an integer and scan its value
 *
 * parameters:
 *   str   <-- string parsed
 *   value --> integer conversion
 *
 * returns:
 *   true if the string defines an integer, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_is_int(const char  *str,
        int         *value)
{
  int _value;
  int retcode, int_len;

  *value = 0;
  retcode = (bool)(sscanf(str, "%i%n", &_value, &int_len));

  if (retcode) {
    if (int_len != (int)strlen(str))
      retcode = 0;
  }

  if (retcode)
    *value = _value;

  return retcode;
}

/*----------------------------------------------------------------------------
 * Add attribute information to a selector.
 *
 * A global list of attributes from all group classes is built and sorted,
 * and the id's of attributes of each class are determined.
 *
 * parameters:
 *   this_selector   <-> pointer to selector
 *   group_class_set <-- group classes description
 *----------------------------------------------------------------------------*/

static void
_assign_attributes(fvm_selector_t               *this_selector,
                   const fvm_group_class_set_t  *group_class_set)
{
  int i, j, group_int;
  int *set_attributes = NULL, *attributes = NULL;
  int n_attributes_tot = 0, n_attributes_max = 0;

  /* Build basic arrays which exist in any case */

  BFT_MALLOC(this_selector->n_class_attributes,
             this_selector->n_group_classes,
             int);
  BFT_MALLOC(this_selector->attribute_ids,
             this_selector->n_group_classes,
             int *);

  for (i = 0; i < this_selector->n_group_classes; i++) {

    int n_attributes = 0;
    const fvm_group_class_t *gc = fvm_group_class_set_get(group_class_set, i);
    const int n_groups = fvm_group_class_get_n_groups(gc);
    const char **group_names = fvm_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++) {
      if (_is_int(group_names[j], &group_int))
        n_attributes += 1;
    }

    n_attributes_tot += n_attributes;
    if (n_attributes > n_attributes_max)
      n_attributes_max = n_attributes;

    this_selector->n_class_attributes[i] = n_attributes;
    if (n_attributes > 0)
      BFT_MALLOC(this_selector->attribute_ids[i], n_attributes, int);
    else
      this_selector->attribute_ids[i] = NULL;
  }

  if (n_attributes_tot == 0)
    return;

  BFT_MALLOC(set_attributes, n_attributes_tot, int);

  n_attributes_tot = 0;

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvm_group_class_t *gc = fvm_group_class_set_get(group_class_set, i);
    const int n_groups = fvm_group_class_get_n_groups(gc);
    const char **group_names = fvm_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++) {
      if (_is_int(group_names[j], &group_int)) {
        set_attributes[n_attributes_tot++] = group_int;
      }
    }

  }

  qsort(set_attributes, n_attributes_tot, sizeof(int), &_compare_attributes);

  BFT_MALLOC(this_selector->attribute, n_attributes_tot, int);

  this_selector->attribute[0] = set_attributes[0];
  for (i = 1, j = 1; i < n_attributes_tot; i++) {
    if (set_attributes[i] != set_attributes[i-1])
      this_selector->attribute[j++] = set_attributes[i];
  }

  BFT_FREE(set_attributes);

  this_selector->n_attributes = j;
  BFT_REALLOC(this_selector->attribute, this_selector->n_attributes, int);

  /* Now build attribute_id arrays */

  BFT_MALLOC(attributes, n_attributes_max, int);

  for (i = 0; i < this_selector->n_group_classes; i++) {

    int mid_id;
    int n_attributes = 0;
    const fvm_group_class_t *gc = fvm_group_class_set_get(group_class_set, i);
    const int n_groups = fvm_group_class_get_n_groups(gc);
    const char **group_names = fvm_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++) {
      if (_is_int(group_names[j], &group_int)) {
        attributes[n_attributes++] = group_int;
      }
    }

    for (j = 0; j < n_attributes; j++) {

      /* use binary search */

      int start_id = 0;
      int end_id = this_selector->n_attributes - 1;
      int val = attributes[j];
      mid_id = (end_id -start_id) / 2;

      /* use binary search */

      while (start_id <= end_id) {
        int att_cmp = this_selector->attribute[mid_id];
        if (att_cmp < val)
          start_id = mid_id + 1;
        else if (att_cmp > val)
          end_id = mid_id - 1;
        else
          break;
        mid_id = start_id + ((end_id -start_id) / 2);
      }

      assert(this_selector->attribute[mid_id] == val);
      this_selector->attribute_ids[i][j] = mid_id;

    }
  }

  BFT_FREE(attributes);
}

/*----------------------------------------------------------------------------
 * Create an operations list.
 *
 * The default number of operations is 30.
 *
 * returns:
 *   pointer to the operations list
 *----------------------------------------------------------------------------*/

static _operation_list_t *
_operation_list_allocate(void)
{
  _operation_list_t *ops;
  int i;
  const int n_operations = 16;

  BFT_MALLOC(ops, 1, _operation_list_t);

  /*  Definitions */

  ops->n_operations = 0;
  ops->n_max_operations = n_operations;

  BFT_MALLOC(ops->postfix,
             ops->n_max_operations,
             fvm_selector_postfix_t *);

  BFT_MALLOC(ops->n_calls, ops->n_max_operations, size_t);

  BFT_MALLOC(ops->n_group_classes, ops->n_max_operations, int);
  BFT_MALLOC(ops->group_class_set, ops->n_max_operations, int *);

  for (i = 0; i < ops->n_max_operations; i++) {
    ops->postfix[i] = NULL;
    ops->group_class_set[i] = NULL;
    ops->n_calls[i] = 0;
    ops->n_group_classes[i] = 0;
  }

  return ops;
}

/*----------------------------------------------------------------------------
 * Increase operations list length.
 *
 * Length is multiplied by 2.
 *
 * parameters:
 *   ops <-> operations list to be updated
 *----------------------------------------------------------------------------*/

static void
_operation_list_reallocate(_operation_list_t  *ops)
{
  int old_size;

  int i = 0;

  old_size = ops->n_max_operations;
  ops->n_max_operations *= 2;

  /* Reallocation */

  BFT_REALLOC(ops->postfix,
              ops->n_max_operations,
              fvm_selector_postfix_t *);

  BFT_REALLOC(ops->n_calls, ops->n_max_operations, size_t);

  BFT_REALLOC(ops->n_group_classes, ops->n_max_operations, int);
  BFT_REALLOC(ops->group_class_set, ops->n_max_operations, int *);

  for (i = old_size; i < ops->n_max_operations; i++) {
    ops->postfix[i] = NULL;
    ops->group_class_set[i] = NULL;
    ops->n_calls[i] = 0;
    ops->n_group_classes[i] = 0;
  }
}

/*----------------------------------------------------------------------------
 * Delete operations list.
 *
 * parameters:
 *   ops <-> operations list to be deleted
 *
 * returns:
 *   NULL pointer;
 *----------------------------------------------------------------------------*/

static _operation_list_t *
_operation_list_free(_operation_list_t  *ops)
{
  int i = 0;

  if (ops != NULL) {
    BFT_FREE(ops->n_calls);
    BFT_FREE(ops->n_group_classes);
    for (i = 0; i < ops->n_max_operations; i++) {
      if (ops->group_class_set[i] != NULL)
        BFT_FREE(ops->group_class_set[i]);
      if (ops->postfix[i] != NULL)
        fvm_selector_postfix_destroy(ops->postfix + i);
    }
    BFT_FREE(ops->postfix);
    BFT_FREE(ops->group_class_set);
    BFT_FREE(ops);
  }

  return NULL;
}

/*----------------------------------------------------------------------------
 * Interpret the postfix string for the last operation of operations list to
 * build the group class list which corresponds to this operation.
 *
 * parameters:
 *   this_selector <-- selector structure
 *   operations    <-> operations list to be updated
 *----------------------------------------------------------------------------*/

static void
_create_operation_group_class_set(const fvm_selector_t  *this_selector,
                                  _operation_list_t     *operations)
{
  int gc_id;
  int *group_class_set;

  int n_group_classes = 0;

  const fvm_selector_postfix_t  *pf
    = operations->postfix[operations->n_operations -1];

  BFT_MALLOC(operations->group_class_set[operations->n_operations - 1],
             this_selector->n_group_classes,
             int);

  group_class_set
    = operations->group_class_set[operations->n_operations - 1];

  for (gc_id = 0; gc_id < this_selector->n_group_classes; gc_id++) {

    /* update group class list for current operation */

    if (fvm_selector_postfix_eval(pf,
                                  this_selector->n_class_groups[gc_id],
                                  this_selector->n_class_attributes[gc_id],
                                  this_selector->group_ids[gc_id],
                                  this_selector->attribute_ids[gc_id],
                                  NULL,
                                  NULL))
      group_class_set[n_group_classes++] = gc_id;
  }

  operations->n_group_classes[operations->n_operations-1] = n_group_classes;

  BFT_REALLOC(operations->group_class_set[operations->n_operations-1],
              n_group_classes,
              int);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("  - group_classes list: ");
  {
    int i ;
    /* fvm_group_class_set_dump(class_defs); */
    for (i = 0; i < n_group_classes; i++)
      bft_printf("%i ",
                 operations->group_class_set[operations->n_operations-1][i]);
    bft_printf("\n");
  }
#endif /* defined(DEBUG) && !defined(NDEBUG) */
}

/*----------------------------------------------------------------------------
 * Add a new operation in the operation_list :
 *
 * - infix string is parsed into postfix string
 * - interpret the postfix string
 *
 * parameters:
 *   selector     <-> selector to be updated
 *   infix_string <-- string parsed
 *----------------------------------------------------------------------------*/

static void
_add_new_operation(fvm_selector_t  *selector,
                   const char      *infix_string)
{
  fvm_selector_postfix_t *pf = NULL;

  /* reallocation  */

  if (   selector->_operations->n_operations
      >= selector->_operations->n_max_operations)
    _operation_list_reallocate(selector->_operations);

  /* Parse infix_string */

  pf = fvm_selector_postfix_create(infix_string,
                                   selector->n_groups,
                                   selector->n_attributes,
                                   (const char **)selector->group_name,
                                   selector->attribute);

  /* update n_operations */

  selector->_operations->postfix[selector->_operations->n_operations] = pf;
  selector->_operations->n_operations++;


  /* Create group class list if there are no numerical tests in postfix */

  if (   fvm_selector_postfix_coords_dep(pf) == false
      && fvm_selector_postfix_normals_dep(pf) == false)
    _create_operation_group_class_set(selector,
                                      selector->_operations);
}

/*----------------------------------------------------------------------------
 * Get the test number in the operation_list which corresponds to
 * the infix string "test_str". If this string doesn't correspond to a an
 * operation already parsed, a new operation is added.
 *
 * parameters:
 *   teststr  <-- string parsed
 *   selector <-> current selector
 *
 * returns:
 *   the number of the operation which corresponds to "test_str" in the
 *   "selector" operations list
 *----------------------------------------------------------------------------*/

static int
_get_criteria_id(fvm_selector_t  *selector,
                 const char      *teststr)
{
  int op = 0;

  assert(teststr != NULL);

  /* Search for teststr in the operations list */

  if (selector->_operations == NULL)
    selector->_operations = _operation_list_allocate();

  for (op = 0; op < selector->_operations->n_operations; op++) {
    const fvm_selector_postfix_t *pf = selector->_operations->postfix[op];
    if (!strcmp(fvm_selector_postfix_get_infix(pf), teststr))
      break;
  }

  /* if teststr is not in the list : add teststrcpy in the list */
  if (op == selector->_operations->n_operations)
    _add_new_operation(selector, teststr);

  return op;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a selector object.
 *
 * parameters:
 *   dim                  <-- spatial dimension (coordinates and normals)
 *   n_elements           <-- number of selectable elements
 *   group_class_set      <-- pointer to group class set definition
 *   group_class_id       <-- group class id associated with each element
 *                            (size: n_elements)
 *   group_class_id_base; <-- Starting group class id base (usually 0 or 1)
 *   coords               <-- coordinates (interlaced) associated with each
 *                            element, whether vertex, face or cell center, ...
 *                            (size: n_elements * dim)
 *   normals              <-- normals (interlaced) associated with each element
 *                            if applicable (such as for face normals), or NULL
 *
 * returns:
 *   pointer to new selector
 *----------------------------------------------------------------------------*/

fvm_selector_t *
fvm_selector_create(int                           dim,
                    cs_lnum_t                     n_elements,
                    const fvm_group_class_set_t  *group_class_set,
                    const int                     group_class_id[],
                    int                           group_class_id_base,
                    const double                  coords[],
                    const double                  normals[])
{
  int i;
  cs_lnum_t j;
  fvm_selector_t *selector;

  int n_group_classes = fvm_group_class_set_size(group_class_set);

  BFT_MALLOC(selector, 1, fvm_selector_t);

  selector->dim = dim;
  selector->n_elements = n_elements;

  selector->group_class_id = group_class_id;
  selector->_group_class_id = NULL;
  selector->group_class_id_base = group_class_id_base;

  selector->n_group_classes = fvm_group_class_set_size(group_class_set);

  selector->n_groups = 0;
  selector->n_attributes = 0;
  selector->group_name = NULL;
  selector->attribute = NULL;

  selector->n_class_groups = NULL;
  selector->group_ids = NULL;
  selector->n_class_attributes = NULL;
  selector->attribute_ids = NULL;

  selector->coords = coords;
  selector->_coords = NULL;
  selector->normals = normals;
  selector->_normals = NULL;

  selector->_operations = NULL;

  selector->_n_group_class_elements = NULL;
  selector->_group_class_elements = NULL;

  _assign_groups(selector, group_class_set);
  _assign_attributes(selector, group_class_set);

  if (group_class_id != NULL && n_group_classes > 0) {

    BFT_MALLOC(selector->_n_group_class_elements, n_group_classes, cs_lnum_t);
    BFT_MALLOC(selector->_group_class_elements, n_group_classes, cs_lnum_t *);

    /* Counting loop and allocation */

    for (i = 0; i < n_group_classes; i++)
      selector->_n_group_class_elements[i] = 0;

    for (j = 0; j < n_elements; j++)
      selector->_n_group_class_elements[  group_class_id[j]
                                        - group_class_id_base] += 1;

    for (i = 0; i < n_group_classes; i++)
      BFT_MALLOC(selector->_group_class_elements[i],
                 selector->_n_group_class_elements[i],
                 int);

    /* Definition loop */

    for (i = 0; i < n_group_classes; i++)
      selector->_n_group_class_elements[i] = 0;

    for (j = 0; j < n_elements; j++) {
      selector->_group_class_elements
                   [group_class_id[j]-1]
                   [selector->_n_group_class_elements
                                [  group_class_id[j]
                                 - group_class_id_base]++]
        = j;

    }
  }

  selector->n_evals = 0;
  selector->eval_wtime = 0.0;

  return selector;
}

/*----------------------------------------------------------------------------
 * Destruction of a selector structure.
 *
 * parameters:
 *   this_selector <-> selector to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_selector_t *
fvm_selector_destroy(fvm_selector_t  *this_selector)
{
  int i;

  /* Delete private arrays */

  _operation_list_free(this_selector->_operations);

  if (this_selector->_coords != NULL)
    BFT_FREE(this_selector->_coords);
  if (this_selector->_normals != NULL)
    BFT_FREE(this_selector->_normals);

  for (i = 0; i < this_selector->n_groups; i++)
    BFT_FREE(this_selector->group_name[i]);
  BFT_FREE(this_selector->group_name);

  BFT_FREE(this_selector->attribute);

  BFT_FREE(this_selector->n_class_groups);
  BFT_FREE(this_selector->n_class_attributes);

  for (i = 0; i < this_selector->n_group_classes; i++) {
    if (this_selector->group_ids[i] != NULL)
      BFT_FREE(this_selector->group_ids[i]);
    if (this_selector->attribute_ids[i] != NULL)
      BFT_FREE(this_selector->attribute_ids[i]);
  }

  BFT_FREE(this_selector->group_ids);
  BFT_FREE(this_selector->attribute_ids);

  if (this_selector->_group_class_elements != NULL) {
    for (i = 0; i < this_selector->n_group_classes; i++)
      BFT_FREE(this_selector->_group_class_elements[i]);

    BFT_FREE(this_selector->_n_group_class_elements);
    BFT_FREE(this_selector->_group_class_elements);
  }

  /* Delete selector */

  BFT_FREE(this_selector);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Define the list of the elements verifying the criteria described
 * by a character string.
 *
 * The selected_element[] array must be pre-allocated, and be of sufficient
 * size to contain all elements associated with the selector.
 *
 * parameters:
 *   this_selector       <-> pointer to selector
 *   str                 <-- string defining selection criteria
 *   elt_id_base         <-- element id base (usually 0 or 1)
 *   n_selected_elements <-- number of elements selected
 *   selected_elements   <-> selected elements list (1 to n numbering)
 *
 * returns:
 *   criteria id associated by selector with str
 *----------------------------------------------------------------------------*/

int
fvm_selector_get_list(fvm_selector_t  *this_selector,
                      const char      *str,
                      cs_lnum_t        elt_id_base,
                      cs_lnum_t       *n_selected_elements,
                      cs_lnum_t       *selected_elements)

{
  int  c_id, gc_id;
  cs_lnum_t   i;
  const fvm_selector_postfix_t *pf = NULL;
  fvm_selector_t  *ts = this_selector;
  double t0 = cs_timer_wtime();

  assert(this_selector != NULL);

  *n_selected_elements = 0;

  /* Add or find the test number in the the cached operations list */

  c_id = _get_criteria_id(ts, str);

  ts->_operations->n_calls[c_id] += 1;
  pf = ts->_operations->postfix[c_id];

  /* Case without geometrical test: get group class list without the
     interpretration of the postfix writing of the criteria string */

  if (   fvm_selector_postfix_coords_dep(pf) == false
      && fvm_selector_postfix_normals_dep(pf) == false) {

    if (ts->_operations->group_class_set[c_id] != NULL) {

      int n_criteria_group_classes
        = ts->_operations->n_group_classes[c_id];
      const int *_criteria_group_class_set
        = ts->_operations->group_class_set[c_id];

      if (ts->_n_group_class_elements != NULL) {

        for (gc_id = 0; gc_id < n_criteria_group_classes; gc_id++) {
          for (i = 0;
               i < ts->_n_group_class_elements
                         [_criteria_group_class_set[gc_id]];
               i++) {
            selected_elements[(*n_selected_elements)++]
              = ts->_group_class_elements
                      [_criteria_group_class_set[gc_id]][i] + elt_id_base;
          }
        }

      }
    }

  }

  /* Case with geometrical test:
     evaluation of the postfix expression for each element */

  else if (ts->n_elements > 0) {

    const int dim = ts->dim;

    assert(ts->group_class_id != NULL);

    if (fvm_selector_postfix_coords_dep(pf) == true && ts->coords == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Selection criteria:\n"
                  "\"%s\"\n"
                  "depends on coordinates, but the current selector\n"
                  "has no associated coordinates."),
                str);
    else if (   fvm_selector_postfix_normals_dep(pf) == true
             && ts->normals == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Selection criteria:\n"
                  "\"%s\"\n"
                  "depends on normals, but the current selector\n"
                  "has no associated normals."),
                str);
    if (dim != 3)
        bft_error(__FILE__, __LINE__, 0,
                  _("Selection criteria:\n"
                  "\"%s\"\n"
                    "is associated with %d spatial dimensions, but\n"
                    "geometric conditions are only currently implemented\n"
                    "for 3 spatial dimension."),
                  str, dim);

    /* Loop on all elements */

    for (i = 0; i < ts->n_elements; i++) {

      gc_id = ts->group_class_id[i] - ts->group_class_id_base;

      if (fvm_selector_postfix_eval(pf,
                                    ts->n_class_groups[gc_id],
                                    ts->n_class_attributes[gc_id],
                                    ts->group_ids[gc_id],
                                    ts->attribute_ids[gc_id],
                                    ts->coords + (i*dim),
                                    ts->normals + (i*dim)))
        selected_elements[(*n_selected_elements)++] = i + elt_id_base;

    }
  }

  ts->n_evals += 1;
  ts->eval_wtime += (cs_timer_wtime() - t0);

  return c_id;
}

/*----------------------------------------------------------------------------
 * Define the list of group classes verifying the criteria described
 * by a character string.
 *
 * The selected_gc[] array must be pre-allocated, and be of sufficient
 * size to contain all elements associated with the selector.
 *
 * parameters:
 *   this_selector  <-> pointer to selector
 *   str            <-- string defining selection criteria
 *   n_selected_gcs <-- number of group classes selected
 *   selected_gcs   <-> selected group class list (0 to n numbering,
 *                      as group class "zero" may exist)
 *
 * returns:
 *   criteria id associated by selector with str
 *----------------------------------------------------------------------------*/

int
fvm_selector_get_gc_list(fvm_selector_t  *this_selector,
                         const char      *str,
                         int             *n_selected_gcs,
                         int              selected_gcs[])
{
  int  c_id, gc_id;
  const fvm_selector_postfix_t *pf = NULL;
  fvm_selector_t  *ts = this_selector;
  double t0 = cs_timer_wtime();

  assert(this_selector != NULL);

  *n_selected_gcs = 0;

  /* Add or find the test number in the the cached operations list */

  c_id = _get_criteria_id(ts, str);

  ts->_operations->n_calls[c_id] += 1;
  pf = ts->_operations->postfix[c_id];

  if (   fvm_selector_postfix_coords_dep(pf) == true
      || fvm_selector_postfix_normals_dep(pf) == true)
    bft_error(__FILE__, __LINE__, 0,
              _("Selection of group classes by criteria:\n"
                "\"%s\"\n"
                "must not depend on coordinates or normals."),
              str);

  /* copy group class list */

  if (ts->_operations->group_class_set[c_id] != NULL) {

    int n_criteria_group_classes
      = ts->_operations->n_group_classes[c_id];
    const int *_criteria_group_class_set
      = ts->_operations->group_class_set[c_id];

    for (gc_id = 0; gc_id < n_criteria_group_classes; gc_id++)
      selected_gcs[gc_id] = _criteria_group_class_set[gc_id];

    *n_selected_gcs = n_criteria_group_classes;

  }

  ts->n_evals += 1;
  ts->eval_wtime += (cs_timer_wtime() - t0);

  return c_id;
}

/*----------------------------------------------------------------------------
 * Return the number of operands associated with a selection criteria
 * which are missing in the selector's associated group class set.
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   criteria_id   <-- id of criteria returned by fvm_selector_get_list()
 *
 * returns:
 *   number of missing operands
 *----------------------------------------------------------------------------*/

int
fvm_selector_n_missing(const fvm_selector_t  *this_selector,
                       int                    criteria_id)
{
  int retval = 0;

  if (this_selector != NULL && criteria_id >= 0) {
    if (   this_selector->_operations != NULL
        && this_selector->_operations->n_operations > criteria_id) {
      const fvm_selector_postfix_t *pf
        = this_selector->_operations->postfix[criteria_id];

      retval = fvm_selector_postfix_n_missing(pf);
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return a pointer to the name of an of operand associated with a selection
 * criteria which is missing in the selector's associated group class set.
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   criteria_id   <-- id of criteria returned by fvm_selector_get_list()
 *   missing_id    <-- id of missing operand for this criteria
 *
 * returns:
 *   pointer to name of missing operand
 *----------------------------------------------------------------------------*/

const char *
fvm_selector_get_missing(const fvm_selector_t  *this_selector,
                         int                    criteria_id,
                         int                    missing_id)
{
  const char *retval = NULL;

  if (this_selector != NULL && criteria_id >= 0) {
    if (   this_selector->_operations != NULL
        && this_selector->_operations->n_operations > criteria_id) {
      const fvm_selector_postfix_t *pf
        = this_selector->_operations->postfix[criteria_id];

      retval = fvm_selector_postfix_get_missing(pf, missing_id);
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Get statistics on selector usage
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   n_evals       <-> number of evaluations, or NULL
 *   eval_wtime    <-> evaluation wall-clock time, or NULL
 *----------------------------------------------------------------------------*/

void
fvm_selector_get_stats(const fvm_selector_t  *this_selector,
                       int                   *n_evals,
                       double                *eval_wtime)
{
  if (this_selector != NULL) {
    if (n_evals != NULL)
      *n_evals = this_selector->n_evals;
    if (eval_wtime != NULL)
      *eval_wtime = this_selector->eval_wtime;
  }
}

/*----------------------------------------------------------------------------
 * Dump the contents of a selector structure in human readable form
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *----------------------------------------------------------------------------*/

void
fvm_selector_dump(const fvm_selector_t  *this_selector)
{
  int i, j;
  const fvm_selector_t  *ts = this_selector;

  if (ts == NULL) {
    bft_printf("\n"
               "Null selector dump:\n");
    return;
  }

  bft_printf("\n"
             "Selector dump:\n"
             "  Dimension:                          %d\n"
             "  Number of selectable elements:      %d\n"
             "  Shared group class id's:            %p\n"
             "  Private group class id's:           %p\n"
             "  Group class id base:                %d\n"
             "  Number of associated group classes: %d\n"
             "  Number of associated groups:        %d\n"
             "  Number of associated attributes:    %d\n"
             "  Number of evaluations:              %d\n"
             "  Wall-clock time in evaluations:     %f\n",
             ts->dim, (int)ts->n_elements,
             (const void *)ts->group_class_id,
             (const void *)ts->_group_class_id,
             ts->group_class_id_base,
             ts->n_group_classes, ts->n_groups, ts->n_attributes,
             ts->n_evals, ts->eval_wtime);

  if (ts->n_groups > 0) {
    bft_printf("  Group names:\n");
    for (i = 0; i < ts->n_groups; i++)
      bft_printf("    \"%s\"\n", ts->group_name[i]);
  }
  if (ts->n_attributes > 0) {
    bft_printf("  Attributes:\n");
    for (i = 0; i < ts->n_attributes; i++)
      bft_printf("    %d\n", ts->attribute[i]);
  }

  if (ts->n_group_classes > 0) {
    bft_printf("  Group classes:\n");
    for (i = 0; i < ts->n_group_classes; i++) {
      bft_printf("    Group class %d\n", (int)i);
      if (ts->n_groups > 0) {
        bft_printf("      Number of groups: %d\n",
                   ts->n_class_groups[i]);
        for (j = 0; j < ts->n_class_groups[i]; j++)
          bft_printf("        %d\n", ts->group_ids[i][j]);
      }
      if (ts->n_attributes > 0) {
        bft_printf("      Number of attributes: %d\n",
                   ts->n_class_attributes[i]);
        for (j = 0; j < ts->n_class_attributes[i]; j++)
          bft_printf("        %d\n", ts->attribute_ids[i][j]);
      }
    }
  }

  bft_printf("  Shared coordinates:                 %p\n"
             "  Private coordinates:                %p\n"
             "  Shared normals;                     %p\n"
             "  Private normals:                    %p\n"
             "  Operations list:                    %p\n",
             (const void *)ts->coords, (const void *)ts->_coords,
             (const void *)ts->normals, (const void *)ts->_normals,
             (const void *)ts->_operations);

  if (ts->n_group_classes > 0) {
    bft_printf("  Number of elements per group class:\n");
    for (i = 0; i < ts->n_group_classes; i++) {
      bft_printf("    %d (%p)\n",
                 (int)ts->_n_group_class_elements[i],
                 (const void *)ts->_group_class_elements[i]);
    }
  }

  if (ts->_operations != NULL) {

    bft_printf("\n");

    for (i = 0; i < ts->_operations->n_operations; i++) {
      bft_printf("  Operation %d (cached, n_calls = %llu)\n",
                 i, (unsigned long long) ts->_operations->n_calls[i]);
      fvm_selector_postfix_dump(ts->_operations->postfix[i],
                                ts->n_groups, ts->n_attributes,
                                (const char **)ts->group_name,
                                ts->attribute);
    }

  }

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
