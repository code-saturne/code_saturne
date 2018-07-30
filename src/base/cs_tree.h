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

#define CS_TREE_NODE_CHAR     (1 << 0)  /* 1: value is a character string */
#define CS_TREE_NODE_INT      (1 << 1)  /* 2: value is an integer */
#define CS_TREE_NODE_REAL     (1 << 2)  /* 4: value is a cs_real_t */
#define CS_TREE_NODE_BOOL     (1 << 3)  /* 8: value is a bool */

#define CS_TREE_NODE_TAG      (1 << 4)  /* 16: node is a tag
                                           (metadata for XML conversion) */

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
 * \brief  Return a child node with a given name.
*
 * The child node must be located directly under the given node (i.e. it is
 * a child, not a grand-child or beyond).
 *
 * This function is similar to \ref cs_tree_get_node, but is simpler
 * (albeit more restricted in scope) and may be faster in cases where
 * one level of the tree sis searched at a time.
 *
 * In case of multiple children sharing the given name, the first such node
 * is returned.
 *
 * \param[in]  node  pointer to the given node
 * \param[in]  name  name of child node
 *
 * \return string value associated to tag if found, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_node_get_child(cs_tree_node_t  *node,
                       const char      *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the next sibling node with the same name (type)
 *         as a given node.
 *
 * The first node of a series is obtained using \ref cs_tree_get_node.
 *
 * \param[in]  node  pointer to the starting node
 *
 * \return pointer to next sibling node with same name, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_node_get_next_of_name(cs_tree_node_t  *node);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search for a child node (used as a tag) with a given name,
 *         and return its associated string value.
 *
 * The child node must be located directly under the given node (i.e. it is
 * a child, not a grand-child or beyond).
 *
 * If the child "tag" node does not exist, NULL is returned.
 *
 * The CS_TREE_NODE_TAG flag is set for child nodes accessed by this function.
 * It is currently only relevant for possible mapping to XML.
 *
 * \param[in]  node  pointer to the given node
 * \param[in]  tag   name of child node used as tag
 *
 * \return string value associated to tag if found, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_tree_node_get_tag(cs_tree_node_t  *node,
                     const char      *tag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a tag to a given node.
 *
 * A tag is simply a string-valued child node.
 *
 * The CS_TREE_NODE_TAG flag is also set for this child.
 * It is currently only relevant for possible mapping to XML.
 *
 * \param[in, out]  node     pointer to the given node
 * \param[in]       tag      name of child node used as tag
 * \param[in]       tag_str  character string to be copied to tag
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_tag(cs_tree_node_t  *node,
                     const char      *tag,
                     const char      *tag_str);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a character string value associated to a node if present.
 *
 * If the node was never accessed before and the value type was not defined,
 * it is set to CS_TREE_NODE_CHAR. If it was previously converted to
 * a different type, an error is returned.
 *
 * \param[in]  node  pointer to a cs_tree_node_t to access, or NULL
 *
 * \return  associated string, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_tree_node_get_value_str(cs_tree_node_t  *node);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a character string value to a node.
 *
 * The assigned value is copied to the node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   pointer to character string
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_value_str(cs_tree_node_t  *node,
                           const char      *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return array of boolean values associated to a node if present.
 *
 * If the value type was not defined, or defined as a string, values are
 * converted and the type flag set to CS_TREE_NODE_BOOL. If it was previously
 * accessed (and converted) using  a different type, an error is returned.
 *
 * The following strings (case-independent) are converted to "true":
 *   "true", "yes", "on", "1".
 * All other strings are converted to "false".
 *
 * \param[in]  node  pointer to a cs_tree_node_t to access, or NULL
 *
 * \return  pointer to associated values, or NULL
 */
/*----------------------------------------------------------------------------*/

const bool *
cs_tree_node_get_values_bool(cs_tree_node_t  *node);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign an array of boolean values to node.
 *
 * The assigned array is copied to the node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       n     number of elements in val
 * \param[in]       val   array of boolean
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_values_bool(cs_tree_node_t  *node,
                             int              n,
                             const bool      *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a single boolean value to a node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   boolean value to assign
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_tree_node_set_value_bool(cs_tree_node_t  *node,
                            bool             val)
{
  cs_tree_node_set_values_bool(node, 1, &val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return an array of integer values associated to a node if present.
 *
 * If the value type was not defined, or defined as a string, values are
 * converted and the type flag set to CS_TREE_NODE_INT. If it was previously
 * accessed (and converted) using  a different type, an error is returned.
 *
 * \param[in]  node  pointer to a cs_tree_node_t to access, or NULL
 *
 * \return  pointer to associated array, or NULL
 */
/*----------------------------------------------------------------------------*/

const int *
cs_tree_node_get_values_int(cs_tree_node_t  *node);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign an array of integer values to a node.
 *
 * The array values are copied to the node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       n     number of elements in val
 * \param[in]       val   array of integers
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_values_int(cs_tree_node_t  *node,
                            int              n,
                            const int       *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a single integer value to a node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   integer value to assign
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_tree_node_set_value_int(cs_tree_node_t  *node,
                           int              val)
{
  cs_tree_node_set_values_int(node, 1, &val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return an array of real values associated to a node if present.
 *
 * If the value type was not defined, or defined as a string, values are
 * converted and the type flag set to CS_TREE_NODE_REAL. If it was previously
 * accessed (and converted) using  a different type, an error is returned.
 *
 * \param[in]  node  pointer to a cs_tree_node_t to access, or NULL
 *
 * \return  pointer to associated array, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_tree_node_get_values_real(cs_tree_node_t  *node);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign an array of real values to a node.
 *
 * The array values are copied to the node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       n     number of elements in val
 * \param[in]       val   array of real values
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_values_real(cs_tree_node_t   *node,
                             int               n,
                             const cs_real_t  *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a single real value to a node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   real value to assign
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_tree_node_set_value_real(cs_tree_node_t  *node,
                            cs_real_t        val)
{
  cs_tree_node_set_values_real(node, 1, &val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string value associated to a child node if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_value_str.

 * \param[in]  node         pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated values, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_tree_node_get_child_value_str(cs_tree_node_t  *node,
                                 const char      *child_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return array of boolean values associated to a child node if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_values_bool.
 *
 * \param[in]  node         pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated values, or NULL
 */
/*----------------------------------------------------------------------------*/

const bool *
cs_tree_node_get_child_values_bool(cs_tree_node_t  *node,
                                   const char      *child_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return an array of integer values associated to a child node
 *         if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_values_int.
 *
 * \param[in]  node        pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated array, or NULL
 */
/*----------------------------------------------------------------------------*/

const int *
cs_tree_node_get_child_values_int(cs_tree_node_t  *node,
                                  const char      *child_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return an array of real values associated to a child node if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_values_real.
 *
 * \param[in]  node         pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated array, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_tree_node_get_child_values_real(cs_tree_node_t  *node,
                                   const char      *child_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node with a child having a given
 *         (character string) tag value.
 *
 * This node is searched for among siblings of a given node sharing the
 * same path (i.e. the same name).
 *
 * Using the following example tree:
 *
 * /
 *   section1
 *   section2
 *     entry
 *       label
 *         (value = a)
 *     entry
 *       label
 *         (value = b)
 *
 * Using \ref cs_tree_get_node(node, "section2/entry") will return
 * the first node with path "section2/entry" (which has a child named
 * "label" with value a).
 *
 * Using \ref cs_tree_get_sibling_with_tag(node, "label", "a") from that
 * node will return the same node, while
 * \ref cs_tree_get_sibling_with_tag(node, "label", "b") will return
 * the second "section2/entry" node.
 *
 * This function can be called from any sibling (not necessarily the
 * first).
 *
 * \param[in]  node       pointer to the starting node
 * \param[in]  tag        name of the required "tag" child
 * \param[in]  tag_value  value of the required "tag" child
 *
 * \return  pointer to the node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_node_get_sibling_with_tag(cs_tree_node_t  *node,
                                  const char      *tag,
                                  const char      *tag_value);

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
 * \brief  Add a node to a tree.
 *
 * This node is located at "path" from the given node
 * level switch is indicated by a "/" in path
 *
 * Exits on error if a node already exists on this path.
 *
 * \param[in, out]  node  pointer to the node where we start searching
 * \param[in]       path  string describing the path access
 *
 * \return  pointer to the new node
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_add_node(cs_tree_node_t  *node,
                 const char      *path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node.
 *
 * This node is located at "path" from the given node
 * level switch is indicated by a "/" in path.
 *
 * In case of multiple nodes sharing the given path, the first such node
 * is returned.
 *
 * \param[in]  node  pointer to the node where we start searching
 * \param[in]  path  string describing the path access
 *
 * \return  pointer to the node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_get_node(cs_tree_node_t   *node,
                 const char       *path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Count number of nodes sharing a given path.
 *
 * \param[in]  node  pointer to the node where we start searching
 * \param[in]  path  string describing the path access
 *
 * \return  number of nodes sharing path
 */
/*----------------------------------------------------------------------------*/

int
cs_tree_get_node_count(cs_tree_node_t  *node,
                       const char      *path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node with a child having a given
 *         (character string) tag value.
 *
 * This node is located at "path" from the given node
 * level switch is indicated by a "/" in path.
 *
 * \param[in]  node       pointer to the node where we start searching
 * \param[in]  path       string describing the path access
 * \param[in]  tag        name of the required "tag" child
 * \param[in]  tag_value  value of the required "tag" child
 *
 * \return  pointer to the node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_get_node_with_tag(cs_tree_node_t   *node,
                          const char       *path,
                          const char       *tag,
                          const char       *tag_value)
{
  cs_tree_node_t *_node = cs_tree_get_node(node, path);
  if (_node != NULL)
    _node = cs_tree_node_get_sibling_with_tag(_node, tag, tag_value);

  return _node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node, adding it if not present.
 *
 * This node is located at "path" from the given node
 * level switch is indicated by a "/" in path.
 *
 * In case of multiple nodes sharing the given path, the first such node
 * is returned.
 *
 * \param[in]  node  pointer to the node where we start searching
 * \param[in]  path  string describing the path access
 *
 * \return  pointer to the node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

static inline cs_tree_node_t *
cs_tree_get_or_add_node(cs_tree_node_t   *node,
                        const char       *path)
{
  cs_tree_node_t  *_node = cs_tree_get_node(node, path);
  if (_node == NULL)
    node = cs_tree_add_node(_node, path);

  return _node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a node in a tree below the given node.
 *
 * \param[in, out]  parent  pointer to the parent node
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
 * \brief  Create and add a string-valued node in a tree below the given node.
 *
 * This node has a string value set to val_str.
 *
 * \param[in, out]  parent   pointer to the parent node
 * \param[in]       name     name of the node to add
 * \param[in]       val_str  value to assign
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
  cs_tree_node_set_value_str(child, val_str);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a boolean-valued node in a tree below the given node.
 *
 * \param[in, out]  parent  pointer to the parent node
 * \param[in]       name    name of the node to add
 * \param[in]       val     value to assign
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
  cs_tree_node_set_values_bool(child, 1, &val);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add an integer-valued node in a tree below the given node.
 *
 * \param[in, out]  parent  pointer to the parent node
 * \param[in]       name    name of the node to add
 * \param[in]       val     value to assign
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
  cs_tree_node_set_values_int(child, 1, &val);

  return child;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add an real-valued node in a tree below the given node.
 *
 * \param[in, out]  parent  pointer to the parent node
 * \param[in]       name    name of the node to add
 * \param[in]       val     value to assign
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
  cs_tree_node_set_values_real(child, 1, &val);

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
 * \brief  Dump a cs_tree_node_t structure starting from a given node
 *
 * \param[in] log    indicate which log file to use
 * \param[in] depth  starting depth in the tree
 * \param[in] node   pointer to a cs_tree_node_t to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_dump(cs_log_t                log,
             int                     depth,
             const cs_tree_node_t   *node);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TREE_H__ */
