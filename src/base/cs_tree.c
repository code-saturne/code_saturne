/*============================================================================
 * Routines to handle a tree structure used to store data and settings
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_tree.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type and structure definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search for a node located at path from root
 *
 * \param[in] root   pointer to the root node where we start searching
 * \param[in] path   string describing the path access
 *
 * \return a pointer to the node
 */
/*----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_or_create_node(cs_tree_node_t   *root,
                     const char       *path)
{
  cs_tree_node_t  *nodes = root;
  cs_tree_node_t  *node = NULL;

  const size_t  path_len = strlen(path);
  char  _name[128];
  char *name = NULL;
  size_t  start = 0, level_len = 0;

  while (start < path_len) {

    const char *p = path + start;
    level_len = strcspn(p, "/");
    if (level_len == 0) { /* path begins with / or // appears */
      start += 1;
      continue;
    }
    if (level_len + 1 == path_len)
      level_len += 1;

    if (level_len > 128) {
      BFT_MALLOC(name, level_len, char);
      strncpy(name, p, level_len);
    }
    else {
      strncpy(_name, p, level_len);
      _name[level_len] ='\0';
      name = _name;
    }

    /* Search for the node with the given name */
    if (nodes->children == NULL)
      nodes = cs_tree_add_child(nodes, name);
    else
      nodes = nodes->children;

    for (node = nodes; node != NULL; node = node->next)
      if (strcmp(node->name, name) == 0)
        break;

    if (node == NULL)
      nodes = cs_tree_add_sibling(nodes, name);
    else
      nodes = node;

    if (name != _name)
      BFT_FREE(name);

    start += level_len + 1;

  } /* Manipulate the path */

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search for a node located at path from root
 *
 * \param[in] root   pointer to the root node where we start searching
 * \param[in] path   string describing the path access
 *
 * \return a pointer to the node
 */
/*----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_node(const cs_tree_node_t   *root,
             const char             *path)
{
  cs_tree_node_t  *nodes = root;
  cs_tree_node_t  *node = NULL;

  const size_t  path_len = strlen(path);
  char  _name[128];
  char *name = NULL;
  size_t  start = 0, level_len = 0;

  while (start < path_len) {

    const char *p = path + start;
    level_len = strcspn(p, "/");
    if (level_len == 0) { /* path begins with / or // appears */
      start += 1;
      continue;
    }
    if (level_len + 1 == path_len)
      level_len += 1;

    /* Search for the node with the given name */
    nodes = nodes->children;
    if (nodes == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Fail to reach the requested node located at %s\n",
                __func__, path);

    if (level_len > 128) {
      BFT_MALLOC(name, level_len, char);
      strncpy(name, p, level_len);
    }
    else {
      strncpy(_name, p, level_len);
      _name[level_len] ='\0';
      name = _name;
    }

    for (node = nodes; node != NULL; node = node->next)
      if (strcmp(node->name, name) == 0)
        break;

    nodes = node;

    if (name != _name)
      BFT_FREE(name);
    if (nodes == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Fail to reach the requested node located at %s\n",
                __func__, path);

    start += level_len + 1;

  } /* Manipulate the path */

  return node;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty node. Only the name is assigned if given
 *
 * \param[in]  name    name of the node or NULL
 *
 * \return a pointer to a new allocated cs_tree_node_t structure
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_node_create(const char  *name)
{
  cs_tree_node_t  *n = NULL;

  BFT_MALLOC(n, 1, cs_tree_node_t);
  if (name != NULL) {
    int  len = strlen(name);
    BFT_MALLOC(n->name, len + 1, char);
    strcpy(n->name, name);
  }
  else
    n->name = NULL;

  n->desc = NULL;
  n->flag = 0;
  n->value = NULL;
  n->size = 0;

  n->parent = n->children = NULL;
  n->prev = n->next = NULL;

  return n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_tree_node_t structure.
 *
 * \param[in, out] pnode  pointer to a pointer to a cs_tree_node_t to free
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_free(cs_tree_node_t  **pnode)
{
  if (pnode == NULL)
    return;

  cs_tree_node_t  *node = *pnode;

  if (node->name != NULL)
    BFT_FREE(node->name);
  if (node->desc != NULL)
    BFT_FREE(node->desc);
  if (node->value != NULL)
    BFT_FREE(node->value);

  BFT_FREE(node);
  node = NULL;
  pnode = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a string
 *
 * \param[in, out] node   pointer to a cs_tree_node_t to modify
 * \param[in]      val    value of the string
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_string_val(cs_tree_node_t  *node,
                            const char      *val)
{
  if (val == NULL)
    return;
  if (node == NULL)
    node = cs_tree_node_create(NULL);

  node->size = 1;
  BFT_MALLOC(node->value, strlen(val) + 1, char);
  sprintf((char *)node->value, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a boolean
 *
 * \param[in, out] node   pointer to a cs_tree_node_t to modify
 * \param[in]      val    boolean
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_bool(cs_tree_node_t  *node,
                      bool             val)
{
  if (node == NULL)
    node = cs_tree_node_create(NULL);

  node->size = 1;
  node->flag |= CS_TREE_NODE_BOOL;
  BFT_MALLOC(node->value, 1, bool);
  ((bool *)node->value)[0] = val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a boolean
 *
 * \param[in, out] node   pointer to a cs_tree_node_t to modify
 * \param[in]      size   number of elements in val
 * \param[in]      val    array of boolean
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_bool_val(cs_tree_node_t  *node,
                          int              size,
                          const bool      *val)
{
  if (val == NULL)
    return;
  if (node == NULL)
    node = cs_tree_node_create(NULL);

  node->size = size;
  node->flag |= CS_TREE_NODE_BOOL;
  BFT_MALLOC(node->value, size, bool);
  memcpy(node->value, val, size*sizeof(bool));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it an integer
 *
 * \param[in, out] node   pointer to a cs_tree_node_t to modify
 * \param[in]      size   number of integers in val
 * \param[in]      val    array of integers
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_int_val(cs_tree_node_t  *node,
                         int              size,
                         const int       *val)
{
  if (val == NULL)
    return;
  if (node == NULL)
    node = cs_tree_node_create(NULL);

  node->size = size;
  node->flag |= CS_TREE_NODE_INTEGER;
  BFT_MALLOC(node->value, size, int);
  memcpy(node->value, val, size*sizeof(int));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it an array of
 *         real values
 *
 * \param[in, out] node   pointer to a cs_tree_node_t to modify
 * \param[in]      size   number of elements in val
 * \param[in]      val    array of real values
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_real_val(cs_tree_node_t   *node,
                          int               size,
                          const cs_real_t  *val)
{
  if (val == NULL)
    return;
  if (node == NULL)
    node = cs_tree_node_create(NULL);

  node->size = size;
  node->flag |= CS_TREE_NODE_REAL;
  BFT_MALLOC(node->value, size, cs_real_t);
  memcpy(node->value, val, size*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_tree_node_t structure.
 *
 * \param[in] log    indicate which log file to use
 * \param[in] depth  shift to apply when printing
 * \param[in] node   pointer to a cs_tree_node_t to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_dump(cs_log_t                log,
                  int                     depth,
                  const cs_tree_node_t   *node)
{
  const int  n_element_by_line = 9;
  char  *shift = NULL;
  char  _shift[65] = "";
  if (depth > 32) {
    BFT_MALLOC(shift, depth*2+1, char);
    for (int i = 0; i < 2*depth; i++)
      shift[i] = ' ';
  }
  else {
    int i = depth;
    while (i > 0) {
      strcat(_shift, "  ");
      i--;
    }
    shift = _shift;
  }

  cs_log_printf(log, "%snode_pointer: %p\n", shift, (const void *)node);
  if (node == NULL) {
    if (shift != _shift)
      BFT_FREE(shift);
    return;
  }

  if (node->name == NULL)
    cs_log_printf(log, "%sname: NULL\n", shift);
  else
    cs_log_printf(log, "%sname: %s\n", shift, node->name);

  if (node->value != NULL) {

    switch (node->size) {

    case 0:
      bft_error(__FILE__, __LINE__, 0,
                " Incompatibility: node->value != NULL and node->size = 0.\n");
      break;

    case 1:
      if (node->flag & CS_TREE_NODE_INTEGER)
        cs_log_printf(log, "%svalue: %d\n", shift, ((int *)node->value)[0]);
      else if (node->flag & CS_TREE_NODE_REAL)
        cs_log_printf(log, "%svalue: %-6.4e\n", shift,
                      ((cs_real_t *)node->value)[0]);
      else if (node->flag & CS_TREE_NODE_BOOL)
        cs_log_printf(log, "%svalue: %s\n", shift,
                      ((bool *)node->value)[0] == true ? "true" : "false");
      else
        cs_log_printf(log, "%svalue: %s\n", shift, (char *)node->value);
      break;

    default:
      {
        const int  n_pass = node->size / n_element_by_line;
        const int  n_last = node->size - n_pass*n_element_by_line;
        cs_log_printf(log, "%svalue: >\n", shift);

        if (node->flag & CS_TREE_NODE_INTEGER) {

          int  *v = (int *)node->value;
          for (int i = 0; i < n_pass; i++) {
            cs_log_printf(log, "%s", shift);
            for (int j = 0; j < n_element_by_line; j++)
              cs_log_printf(log, "%d", v[n_element_by_line*i + j]);
            cs_log_printf(log, "\n");
          }
          if (n_last > 0) {
            cs_log_printf(log, "%s", shift);
            for (int j = 0; j < n_last; j++)
              cs_log_printf(log, "%d", v[n_element_by_line*n_pass + j]);
            cs_log_printf(log, "\n");
          }

        }
        else if (node->flag & CS_TREE_NODE_REAL) {

          cs_real_t  *v = (cs_real_t *)node->value;
          for (int i = 0; i < n_pass; i++) {
            cs_log_printf(log, "%s", shift);
            for (int j = 0; j < n_element_by_line; j++)
              cs_log_printf(log, "%-6.4e", v[n_element_by_line*i + j]);
            cs_log_printf(log, "\n");
          }
          if (n_last > 0) {
            cs_log_printf(log, "%s", shift);
            for (int j = 0; j < n_last; j++)
              cs_log_printf(log, "%-6.4e", v[n_element_by_line*n_pass + j]);
            cs_log_printf(log, "\n");
          }

        }
        else if (node->flag & CS_TREE_NODE_BOOL) {

          bool  *v = (bool *)node->value;
          for (int i = 0; i < n_pass; i++) {
            cs_log_printf(log, "%s", shift);
            for (int j = 0; j < n_element_by_line; j++)
              cs_log_printf(log, "%s", (v[n_element_by_line*i + j] == true ?
                                        "true" : "false"));
            cs_log_printf(log, "\n");
          }
          if (n_last > 0) {
            cs_log_printf(log, "%s", shift);
            for (int j = 0; j < n_last; j++)
              cs_log_printf(log, "%s",
                            (v[n_element_by_line*n_pass + j] == true) ?
                            "true" : "false");
            cs_log_printf(log, "\n");
          }
        }
        else /* value is a string */
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Array of strings is not handled\n", __func__);

      }
      break;

    } /* Switch on node->size */

  } /* Dump value */

  cs_log_printf(log, "%sflag: %d\n", shift, node->flag);
  if (node->desc != NULL)
    cs_log_printf(log, "%sdesc: |\n%s\n", shift, node->desc);

  if (shift != _shift)
    BFT_FREE(shift);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node with read/write access.
 *         This node is located at "path" from the given root node
 *         level switch is indicated by a "/" in path
 *         Exit on error if not found
 *
 * \param[in, out] root   pointer to the root node where we start searching
 * \param[in]      path   string describing the path access
 *
 * \return a pointer to the node
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_get_node_rw(cs_tree_node_t     *root,
                    const char         *path)
{
  if (path == NULL)
    return root;
  if (strlen(path) == 0)
    return root;
  if (root == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: root is NULL\n", __func__);

  return _find_or_create_node(root, path);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node in read-only access.
 *         This node is located at "path" from the given root node
 *         level switch is indicated by a "/" in path
 *         Exit on error if not found
 *
 * \param[in] root   pointer to the root node where we start searching
 * \param[in] path   string describing the path access
 *
 * \return a pointer to the node
 */
/*----------------------------------------------------------------------------*/

const cs_tree_node_t *
cs_tree_get_node(const cs_tree_node_t   *root,
                 const char             *path)
{
  if (path == NULL)
    return root;
  if (strlen(path) == 0)
    return root;
  if (root == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: root is NULL\n", __func__);

  return _find_node(root, path);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree below the given parent node
 *
 * \param[in, out] parent    pointer to the parent node to handle
 * \param[in]      name      name of the node to add
 *
 * \return a pointer to the node in the new tree structure
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_add_child(cs_tree_node_t  *parent,
                  const char      *name)
{
  /* Allocate a new node */
  cs_tree_node_t  *node = cs_tree_node_create(name);

  if (parent == NULL) { /* node is a root */
    node->parent = NULL;
    node->prev = node->next = NULL;
    return node;
  }

  /* The current node is a child of parent with no successor */
  node->parent = parent;
  node->next = NULL;

  /* Locate the child among other children */
  cs_tree_node_t  *child = parent->children; /* First child */
  if (child == NULL) { /* First child is this node */
    parent->children = node;
    node->prev = NULL;
  }
  else {
    while (child->next != NULL)
      child = child->next;
    child->next = node;
    node->prev = child;
  }

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree at the right of the given  node
 *
 * \param[in, out] sibling   pointer to the sibling node to handle
 * \param[in]      name      name of the node to add
 *
 * \return a pointer to the node in the new tree structure
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_add_sibling(cs_tree_node_t  *sibling,
                    const char      *name)
{
  /* Allocate a new node */
  cs_tree_node_t  *node = cs_tree_node_create(name);

  if (sibling == NULL) { /* node is a root */
    node->parent = NULL;
    node->prev = node->next = NULL;
    return node;
  }

  node->parent = sibling->parent;

  node->next = sibling->next;
  node->prev = sibling;
  sibling->next = node;

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a branch in a tree starting from root. If root is the real
 *         root of the tree, the whole tree is freed.
 *
 * \param[in, out] proot  pointer to a pointer to a cs_tree_node_t to free
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_free(cs_tree_node_t  **proot)
{
  if (proot == NULL)
    return;

  cs_tree_node_t  *root = *proot;
  if (root == NULL)
    return;

  if (root->children != NULL) { /* There is at least one child */
    cs_tree_node_t  *next_child = root->children->next;
    while (next_child != NULL) {
      cs_tree_node_t  *tmp = next_child->next;
      cs_tree_free(&next_child);
      next_child = tmp;
    }
    cs_tree_free(&(root->children));
  }

  cs_tree_node_free(&root);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_tree_node_t structure starting from the node "root"
 *
 * \param[in] log     indicate which log file to use
 * \param[in] depth   starting depth in the tree
 * \param[in] root    pointer to a cs_tree_node_t to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_dump(cs_log_t                log,
             int                     depth,
             const cs_tree_node_t   *root)
{
  if (depth < 0)
    depth = 0;
  cs_tree_node_dump(log, depth, root);
  if (root == NULL)
    return;

  if (root->children != NULL) { /* There is at least one child */
    cs_tree_node_t  *next_child = root->children;
    while (next_child != NULL) {
      cs_tree_dump(log, depth+1, next_child);
      next_child = next_child->next;
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
