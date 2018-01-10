#ifndef __CS_TREE_H__
#define __CS_TREE_H__

/*============================================================================
 * Tree structure used to store data and settings.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Set of bit masks to tune finely the behavior of a node.
   By default node value is assumed to be a string.
 */

#define CS_TREE_NODE_INTEGER  (1 << 0)  /* 1: value is an integer */
#define CS_TREE_NODE_REAL     (1 << 1)  /* 2: value is a cs_real_t */
#define CS_TREE_NODE_BOOL     (1 << 2)  /* 4: value is a bool */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_tree_node_t  cs_tree_node_t;

struct _cs_tree_node_t {

  char       *name;   /* name of the node */
  char       *desc;   /* NULL or short description/help about this node */
  int         flag;   /* metadata used to specified the node behavior */

  void       *value;  /* value related to this node. Cast on-the-fly */
  int         size;   /* size > 1 if it is an array */

  /* Pointers to other nodes to allow an easy navigation among the tree */

  cs_tree_node_t  *parent;    /* Pointer to the parent node or NULL if root */
  cs_tree_node_t  *children;  /* Pointer to the first child or NULL */
  cs_tree_node_t  *prev;      /* Pointer to a previous node sharing the same
                                 parent or NULL if this is the first child */
  cs_tree_node_t  *next;      /* Pointer to a next node sharing the same parent
                                 or NULL if there is no next node */

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty node.
 *
 * Only the name is assigned if given
 *
 * \param[in]  name  name of the node, or NULL
 *
 * \return  pointer to a new allocated cs_tree_node_t structure
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_node_create(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a branch in a tree starting from a node.
 *
 *  If the node is the root of the tree, the whole tree is freed.
 *
 * \param[in, out]  pnode  pointer to a pointer to a cs_tree_node_t to free
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_free(cs_tree_node_t  **pnode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Name or rename a node.
 *
 * \param[in, out]  node    pointer to a cs_tree_node_t to modify
 * \param[in]       name    name to set
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_name(cs_tree_node_t  *node,
                      const char      *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a string.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   value of the string
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_val_string(cs_tree_node_t  *node,
                            const char      *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a string.
 *
 * \deprecated Use cs_tree_node_set_val_string instead.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   value of the string
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_tree_node_set_string_val(cs_tree_node_t  *node,
                            const char      *val)
{
  cs_tree_node_set_val_string(node, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a boolean.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   boolean
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_bool(cs_tree_node_t  *node,
                      bool             val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it a boolean.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       size  number of elements in val
 * \param[in]       val   array of boolean
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_bool_val(cs_tree_node_t  *node,
                          int              size,
                          const bool      *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it an integer
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       size  number of integers in val
 * \param[in]       val   array of integers
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_int_val(cs_tree_node_t  *node,
                         int              size,
                         const int       *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the value member of a node and assign to it an array of
 *         real values
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       size  number of elements in val
 * \param[in]       val   array of real values
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_real_val(cs_tree_node_t   *node,
                          int               size,
                          const cs_real_t  *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_tree_node_t structure.
 *
 * \param[in]  log    indicate which log file to use
 * \param[in]  depth  shift to apply when printing
 * \param[in]  node   pointer to a cs_tree_node_t to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_dump(cs_log_t                log,
                  int                     depth,
                  const cs_tree_node_t   *node);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the value member of the node.
 *
 * This node can be modified.
 *
 * This node is located at "path" from the given root node
 * level switch is indicated by a "/" in path
 *
 * Exits on error if not found.
 *
 * \param[in, out]  root  pointer to the root node where we start searching
 * \param[in]       path  string describing the path access
 *
 * \return  pointer to the node
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_get_node_rw(cs_tree_node_t     *root,
                    const char         *path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node with read-only access.
 *
 * This node is located at "path" from the given root node
 * level switch is indicated by a "/" in path.
 *
 * Exits on error if not found.
 *
 * \param[in]  root  pointer to the root node where we start searching
 * \param[in]  path  string describing the path access
 *
 * \return  pointer to the node
 */
/*----------------------------------------------------------------------------*/

const cs_tree_node_t *
cs_tree_get_node(const cs_tree_node_t   *root,
                 const char             *path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a read-only pointer to the value member of the node
 *         in case of a string.
 *
 * This node is located at "path" from the given root node
 * level switch is indicated by a "/" in path.
 *
 * \param[in, out]  root  pointer to the root node where we start searching
 * \param[in]       path  string describing the path access
 *
 * \return  pointer to the string
 */
/*----------------------------------------------------------------------------*/

static inline const char *
cs_tree_get_node_string_val(const cs_tree_node_t   *root,
                            const char             *path)
{
  const cs_tree_node_t  *node = cs_tree_get_node(root, path);
  return (const char *)(node->value);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree below the given parent node.
 *
 * \param[in, out]  parent  pointer to the parent node to handle.
 * \param[in]       name    name of the node to add
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_add_child(cs_tree_node_t  *parent,
                  const char      *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree below the given parent node.
 *
 * This node has a string value set to val_str.
 *
 * \param[in, out]  parent   pointer to the parent node to handle
 * \param[in]       name     name of the node to add
 * \param[in]       val_str  node value to set
 *
 * \return a pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_child_str(cs_tree_node_t  *parent,
                      const char      *name,
                      const char      *val_str)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_string_val(child, val_str);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree below the given parent node.
 *
 * This node has a string value set to val_str.
 *
 * \deprecated Use cs_tree_add_child_str instead.
 *
 * \param[in, out]  parent   pointer to the parent node to handle
 * \param[in]       name     name of the node to add
 * \param[in]       val_str  node value to set
 *
 * \return a pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_string_child(cs_tree_node_t  *parent,
                         const char      *name,
                         const char      *val_str)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_string_val(child, val_str);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree at the right of the given  node.
 *
 * \param[in, out]  sibling  pointer to the sibling node to handle
 * \param[in]       name     name of the node to add
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_add_sibling(cs_tree_node_t  *sibling,
                    const char      *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a boolean-valued node in a tree below
 *         the given parent node.
 *
 * \param[in, out]  parent  pointer to the parent node to handle
 * \param[in]       name    name of the node to add
 * \param[in]       val     node value to assign
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_child_bool(cs_tree_node_t  *parent,
                       const char      *name,
                       bool             val)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_bool_val(child, 1, &val);
  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a boolean-valued node in a tree below
 *         the given parent node.
 *
 * \deprecated Use cs_tree_add_child_bool instead.
 *
 * \param[in, out]  parent  pointer to the parent node to handle
 * \param[in]       name    name of the node to add
 * \param[in]       val     node value to assign
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_bool_child(cs_tree_node_t  *parent,
                       const char      *name,
                       bool             val)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_bool_val(child, 1, &val);
  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add an integer-valued node in a tree
 *         below the given parent node.
 *
 * \param[in, out]  parent  pointer to the parent node to handle
 * \param[in]       name    name of the node to add
 * \param[in]       val     node value to assign
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_child_int(cs_tree_node_t  *parent,
                      const char      *name,
                      int              val)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_int_val(child, 1, &val);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add an integer-valued node in a tree
 *         below the given parent node.
 *
 * \deprecated Use cs_tree_add_child_int instead.
 *
 * \param[in, out]  parent  pointer to the parent node to handle
 * \param[in]       name    name of the node to add
 * \param[in]       val     node value to assign
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_int_child(cs_tree_node_t  *parent,
                      const char      *name,
                      int              val)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_int_val(child, 1, &val);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a real-valued node in a tree
 *         below the given parent node.
 *
 * \param[in, out]  parent  pointer to the parent node to handle
 * \param[in]       name    name of the node to add
 * \param[in]       val     node value to assign
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_child_real(cs_tree_node_t  *parent,
                       const char      *name,
                       cs_real_t        val)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_real_val(child, 1, &val);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a real-valued node in a tree
 *         below the given parent node.
 *
 * \deprecated Use cs_tree_add_child_real instead.
 *
 * \param[in, out]  parent  pointer to the parent node to handle
 * \param[in]       name    name of the node to add
 * \param[in]       val     node value to assign
 *
 * \return  pointer to the new node in the tree structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_add_real_child(cs_tree_node_t  *parent,
                       const char      *name,
                       cs_real_t        val)
{
  cs_tree_node_t  *child = cs_tree_add_child(parent, name);
  cs_tree_node_set_real_val(child, 1, &val);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_tree_node_t structure starting from the node "root".
 *
 * \param[in] log    indicate which log file to use
 * \param[in] depth  starting depth in the tree
 * \param[in] root   pointer to a cs_tree_node_t to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_dump(cs_log_t                log,
             int                     depth,
             const cs_tree_node_t   *root);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TREE_H__ */
