/*============================================================================
 * Tree structure used to store data and settings.
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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_tree.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_tree.c
        Tree structure used to store data and settings.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type and structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static const int _any_type
  = (  CS_TREE_NODE_CHAR | CS_TREE_NODE_INT
     | CS_TREE_NODE_REAL | CS_TREE_NODE_BOOL);

static const int _no_char_type
  = (CS_TREE_NODE_INT | CS_TREE_NODE_REAL | CS_TREE_NODE_BOOL);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Normalize a string, converting it to lowercase and
 *         transforming all line control/tab characters to single whitespace.
 *
 * \param[in, out]  s  string to normalize
 *
 * \return  normalized string length
 */
/*----------------------------------------------------------------------------*/

static size_t
_normalize_string(char *s)
{
  size_t l = strlen(s);
  size_t j = 0;

  for (size_t i = 0; i < l; i++)
    s[i] = tolower(s[i]);

  for (size_t i = 0; i < l; i++) {
    if (s[i] == '\t' || s[i] == '\n' || s[i] == '\r') {
      if (j > 0) {
        if (s[j] != ' ')
          s[j++] = ' ';
      }
    }
    else
      s[j++] = s[i];
  }

  if (j > 0) {
    if (s[j-1] == ' ')
      j -= 1;
  }

  s[j] = '\0';

  return j;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search for a node located at path from node
 *
 * If the node does not exist, it is created.
 *
 * \param[in] node  pointer to the node where we start searching
 * \param[in] path  string describing the path access
 *
 * \return a pointer to the node
 */
/*----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_or_create_node(cs_tree_node_t   *node,
                     const char       *path)
{
  cs_tree_node_t  *_nodes = node;
  cs_tree_node_t  *_node = NULL;

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
    if (_nodes->children == NULL)
      _nodes = cs_tree_add_child(_nodes, name);
    else
      _nodes = _nodes->children;

    for (_node = _nodes; _node != NULL; _node = _node->next)
      if (strcmp(_node->name, name) == 0)
        break;

    if (_node == NULL)
      _nodes = cs_tree_add_sibling(_nodes, name);
    else
      _nodes = _node;

    if (name != _name)
      BFT_FREE(name);

    start += level_len + 1;

  } /* Manipulate the path */

  return _node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search for a node located at path from node
 *
 * If the node does not exist, NULL is returned.
 *
 * \param[in] parent  pointer to the node where we start searching
 * \param[in] path    string describing the path access
 *
 * \return a pointer to the node, or NULL
 */
/*----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_node(cs_tree_node_t   *parent,
           const char       *path)
{
  cs_tree_node_t  *_nodes = (cs_tree_node_t *)parent;
  cs_tree_node_t  *_node = NULL;

  const char *p = path;

  while (p[0] != '\0') {

    if (p[0] == '/') { /* path begins with "/" */
      p++;
      continue;
    }

    _nodes = _nodes->children;
    if (_nodes == NULL)
      return _nodes;

    size_t level_len = 0;
    while (p[level_len] != '/' && p[level_len] != '\0')
      level_len++;

    /* Search for the node with the given name */

    for (_node = _nodes; _node != NULL; _node = _node->next)
      if (strncmp(_node->name, p, level_len) == 0) {
        if (strlen(_node->name) == level_len)
          break;
      }

    _nodes = _node;

    if (_nodes == NULL)
      return NULL;

    p += level_len;

  } /* Manipulate the path */

  return _node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node matching a given sub-path.
 *
 * This node is located at "path" from the given node or one of its
 * descendants, with the path separator indicated by a "/".
 *
 * In case of multiple nodes sharing the given path, the first such node
 * is returned, using a depth-first search.
 *
 * \param[in]  root      pointer to the root node where we start searching
 * \param[in]  sub_path  string describing the path access
 *
 * \return  pointer to the first matching node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_sub_node(cs_tree_node_t  *root,
               const char      *sub_path)
{
  cs_tree_node_t *retval = NULL;
  cs_tree_node_t *tn = root->children;

  /* Check root first */
  retval = cs_tree_get_node(root, sub_path);

  /* Recursively search descendants */
  while (retval == NULL && tn != NULL) {
    retval = _find_sub_node(tn, sub_path);
    tn = tn->next;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a with a given name.
 *
 * In case of multiple nodes sharing the given name, the first such node
 * is returned, using a depth-first search.
 *
 * \param[in]  root  pointer to the root node where we start searching
 * \param[in]  name  node name searched for
 *
 * \return  pointer to the first matching node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_sub_node_simple(cs_tree_node_t  *root,
                      const char      *name)
{
  cs_tree_node_t *retval = NULL;
  cs_tree_node_t *tn = root->children;

  /* Check root first */
  retval = cs_tree_node_get_child(root, name);

  /* Recursively search descendants */
  while (retval == NULL && tn != NULL) {
    retval = _find_sub_node_simple(tn, name);
    tn = tn->next;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_tree_node_t structure.
 *
 * \param[in, out]  pnode  pointer to a pointer to a cs_tree_node_t to free
 */
/*----------------------------------------------------------------------------*/

static void
_node_free(cs_tree_node_t  **pnode)
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

  BFT_FREE(*pnode);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
cs_tree_node_create(const char  *name)
{
  cs_tree_node_t  *n = NULL;

  BFT_MALLOC(n, 1, cs_tree_node_t);
  if (name != NULL) {
    size_t  len = strlen(name);
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
 * \brief  Free a branch in a tree starting from a node.
 *
 *  If the node is the root of the tree, the whole tree is freed.
 *
 * \param[in, out]  pnode  pointer to a pointer to a cs_tree_node_t to free
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_free(cs_tree_node_t  **pnode)
{
  if (pnode == NULL)
    return;

  cs_tree_node_t  *root = *pnode;
  if (root == NULL)
    return;

  if (root->children != NULL) { /* There is at least one child */
    cs_tree_node_t  *next_child = root->children->next;
    while (next_child != NULL) {
      cs_tree_node_t  *tmp = next_child->next;
      cs_tree_node_free(&next_child);
      next_child = tmp;
    }
    cs_tree_node_free(&(root->children));
  }

  _node_free(&root);
}

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
                      const char      *name)
{
  if (name == NULL)
    BFT_FREE(node->name);

  else {
    BFT_REALLOC(node->name, strlen(name) + 1, char);
    strcpy(node->name, name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a child node with a given name.
*
 * The child node must be located directly under the given node (i.e. it is
 * a child, not a grand-child or beyond).
 *
 * This function is similar to \ref cs_tree_get_node, but is simpler
 * (albeit more restricted in scope) and may be faster in cases where
 * one level of the tree is searched at a time.
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
                       const char      *name)
{
  cs_tree_node_t  *child = NULL;
  if (node != NULL) {

    child = node->children;
    while (child != NULL) { /* Search for the node with the given name */
      if (strcmp(child->name, name) == 0)
        break;
      else
        child = child->next;
    }

  }
  return child;
}

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
cs_tree_node_get_next_of_name(cs_tree_node_t  *node)
{
  cs_tree_node_t  *sibling = NULL;
  if (node != NULL) {

    sibling = node->next;
    while (sibling != NULL) { /* Search for the node with the given name */
      if (strcmp(sibling->name, node->name) == 0)
        break;
      else
        sibling = sibling->next;
    }

  }
  return sibling;
}

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
                     const char      *tag)
{
  const char *retval = NULL;

  if (node != NULL) {

    cs_tree_node_t  *child = node->children;
    while (child != NULL) { /* Search for the node with the given name */
      if (strcmp(child->name, tag) != 0)
        child = child->next;
      else
        break;
    }

    if (child != NULL) {
      if (child->flag & CS_TREE_NODE_CHAR)
        retval = (const char *)(child->value);

      else if (child->flag & _no_char_type)
        bft_error(__FILE__, __LINE__, 0,
                  "Tree node %s accessed as type %d (string),\n"
                  "but previously accessed as type %d.",
                  child->name, CS_TREE_NODE_CHAR, (child->flag & _no_char_type));

      else {
        retval = (const char *)(child->value);
        child->flag =   ((child->flag | _any_type) - _any_type)
                      | CS_TREE_NODE_CHAR;
      }
      if (! (child->flag &CS_TREE_NODE_TAG))
        child->flag = child->flag | CS_TREE_NODE_TAG;
    }

  }

  return retval;
}

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
                     const char      *tag_str)
{
  cs_tree_node_t  *child = cs_tree_node_get_child(node, tag);
  if (child == NULL)
    child = cs_tree_add_child(node, tag);

  cs_tree_node_set_value_str(child, tag_str);
  child->flag |= CS_TREE_NODE_TAG;
}

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
cs_tree_node_get_value_str(cs_tree_node_t  *node)
{
  const char *retval = NULL;

  if (node != NULL) {

    if (node->flag & CS_TREE_NODE_CHAR)
      retval = (const char *)(node->value);

    else if (node->flag & _no_char_type)
      bft_error(__FILE__, __LINE__, 0,
                "Tree node %s accessed as type %d (string),\n"
                "but previously accessed as type %d.",
                node->name, CS_TREE_NODE_CHAR, (node->flag & _no_char_type));

    else {
      retval = (const char *)(node->value);
      node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_CHAR;
    }

  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a character string value to a node.
 *
 * The assigned string is copied to the node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       val   pointer to character string
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_value_str(cs_tree_node_t  *node,
                           const char      *val)
{
  assert(node != NULL);

  node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_CHAR;

  if (val == NULL) {
    BFT_FREE(node->value);
    return;
  }

  node->size = 1;
  BFT_REALLOC(node->value, strlen(val) + 1, char);
  strcpy((char *)node->value, val);
}

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
 * \return  pointer to associated array, or NULL
 */
/*----------------------------------------------------------------------------*/

const bool *
cs_tree_node_get_values_bool(cs_tree_node_t  *node)
{
  const bool *retval = NULL;

  if (node != NULL) {

    if (node->flag & CS_TREE_NODE_BOOL)
      retval = (const bool *)(node->value);

    else if (node->flag & _no_char_type)
      bft_error(__FILE__, __LINE__, 0,
                "Tree node %s accessed as type %d (boolean),\n"
                "but previously accessed as type %d.",
                node->name, CS_TREE_NODE_BOOL, (node->flag & _no_char_type));

    else {
      char *s = node->value;
      size_t l = _normalize_string(s);
      bool *v = NULL;
      if (l > 0) {
        node->size = 1;
        for (size_t i = 0; i < l; i++) {
          if (s[i] == ' ')
            node->size += 1;
        }
        BFT_MALLOC(v, node->size, bool);
        int n = 0;
        size_t i = 0;
        while (i < l) {
          size_t j;
          for (j = i; j < l+1; j++) {
            if (s[j] == ' ' || s[j] == '\0') {
              s[j] = '\0';
              break;
            }
          }
          const char *_s = s + i;
          if (   strcmp(_s, "true") == 0
              || strcmp(_s, "yes") == 0
              || strcmp(_s, "on") == 0
              || strcmp(s, "1") == 0)
            v[n] = true;
          else
            v[n] = false;
          n++;
          i = j;
        }
        assert(node->size == n);
      }
      BFT_FREE(node->value);
      node->value = v;
      node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_BOOL;
      retval = (const bool *)(node->value);
    }

  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign an array of boolean values to node.
 *
 * \param[in, out]  node  pointer to a cs_tree_node_t to modify
 * \param[in]       n     number of elements in val
 * \param[in]       val   array of boolean
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_node_set_values_bool(cs_tree_node_t  *node,
                             int              n,
                             const bool      *val)
{
  assert(node != NULL);

  if (val == NULL)
    n = 0;

  node->size = n;
  node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_BOOL;
  BFT_REALLOC(node->value, node->size, bool);

  if (node->size > 0)
    memcpy(node->value, val, node->size*sizeof(bool));
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
cs_tree_node_get_values_int(cs_tree_node_t  *node)
{
  const int *retval = NULL;

  if (node != NULL) {

    if (node->flag & CS_TREE_NODE_INT)
      retval = (const int *)(node->value);

    else if (node->flag & _no_char_type)
      bft_error(__FILE__, __LINE__, 0,
                "Tree node %s accessed as type %d (integer),\n"
                "but previously accessed as type %d.",
                node->name, CS_TREE_NODE_INT, (node->flag & _no_char_type));

    else {
      char *s = node->value;
      size_t l = _normalize_string(s);
      int *v = NULL;
      if (l > 0) {
        node->size = 1;
        for (size_t i = 0; i < l; i++) {
          if (s[i] == ' ')
            node->size += 1;
        }
        BFT_MALLOC(v, node->size, int);
        int n = 0;
        size_t i = 0;
        while (i < l) {
          size_t j;
          for (j = i; j < l+1; j++) {
            if (s[j] == ' ' || s[j] == '\0') {
              s[j] = '\0';
              break;
            }
          }
          errno = 0;
          v[n] = strtol(s+i, NULL, 10);
          if (errno != 0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Error parsing \"%s\" as integer:\n\n"
                        "  %s"), s+i, strerror(errno));
          n++;
          i = j;
        }
        assert(node->size == n);
      }
      BFT_FREE(node->value);
      node->value = v;
      node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_INT;
      retval = (const int *)(node->value);
    }

  }

  return retval;
}

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
                            const int       *val)
{
  assert(node != NULL);

  if (val == NULL)
    n = 0;

  node->size = n;
  node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_INT;
  BFT_REALLOC(node->value, node->size, int);

  if (node->size > 0)
    memcpy(node->value, val, node->size*sizeof(int));
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
cs_tree_node_get_values_real(cs_tree_node_t  *node)
{
  const cs_real_t *retval = NULL;

  if (node != NULL) {

    if (node->flag & CS_TREE_NODE_REAL)
      retval = (const cs_real_t *)(node->value);

    else if (node->flag & _no_char_type)
      bft_error(__FILE__, __LINE__, 0,
                "Tree node %s accessed as type %d (real),\n"
                "but previously accessed as type %d.",
                node->name, CS_TREE_NODE_REAL, (node->flag & _no_char_type));

    else {
      char *s = node->value;
      size_t l = _normalize_string(s);
      cs_real_t *v = NULL;
      if (l > 0) {
        node->size = 1;
        for (size_t i = 0; i < l; i++) {
          if (s[i] == ' ')
            node->size += 1;
        }
        BFT_MALLOC(v, node->size, cs_real_t);
        int n = 0;
        size_t i = 0;
        while (i < l) {
          size_t j;
          for (j = i; j < l+1; j++) {
            if (s[j] == ' ' || s[j] == '\0') {
              s[j] = '\0';
              break;
            }
          }
          errno = 0;
          v[n] = strtod(s+i, NULL);
          if (errno != 0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Error parsing \"%s\" as real:\n\n"
                        "  %s"), s+i, strerror(errno));
          n++;
          i = j;
        }
        assert(node->size == n);
      }
      BFT_FREE(node->value);
      node->value = v;
      node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_REAL;
      retval = (const cs_real_t *)(node->value);
    }

  }

  return retval;
}

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
                             const cs_real_t  *val)
{
  assert(node != NULL);

  if (val == NULL)
    n = 0;

  node->size = n;
  node->flag = ((node->flag | _any_type) - _any_type) | CS_TREE_NODE_REAL;
  BFT_REALLOC(node->value, node->size, cs_real_t);

  if (node->size > 0)
    memcpy(node->value, val, node->size*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string value associated to a child node if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_value_str.

 * \param[in]  node        pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated values, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_tree_node_get_child_value_str(cs_tree_node_t  *node,
                                 const char      *child_name)
{
  const char *retval = NULL;

  cs_tree_node_t *child = cs_tree_node_get_child(node, child_name);
  if (child != NULL)
    retval = cs_tree_node_get_value_str(child);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return array of boolean values associated to a child node if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_values_bool.
 *
 * \param[in]  node        pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated values, or NULL
 */
/*----------------------------------------------------------------------------*/

const bool *
cs_tree_node_get_child_values_bool(cs_tree_node_t  *node,
                                   const char      *child_name)
{
  const bool *retval = NULL;

  cs_tree_node_t *child = cs_tree_node_get_child(node, child_name);
  if (child != NULL)
    retval = cs_tree_node_get_values_bool(child);

  return retval;
}

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
                                  const char      *child_name)
{
  const int *retval = NULL;

  cs_tree_node_t *child = cs_tree_node_get_child(node, child_name);
  if (child != NULL)
    retval = cs_tree_node_get_values_int(child);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return an array of real values associated to a child node if present.
 *
 * The behavior is similar to that of \ref cs_tree_node_get_values_real.
 *
 * \param[in]  node        pointer to a cs_tree_node_t to access, or NULL
 * \param[in]  child_name  name of child node
 *
 * \return  pointer to associated array, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_tree_node_get_child_values_real(cs_tree_node_t  *node,
                                   const char      *child_name)
{
  const cs_real_t *retval = NULL;

  cs_tree_node_t *child = cs_tree_node_get_child(node, child_name);
  if (child != NULL)
    retval = cs_tree_node_get_values_real(child);

  return retval;
}

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
 * Using \ref cs_tree_node_get_sibling_with_tag(node, "label", "a") from that
 * node will return the same node, while
 * \ref cs_tree_node_get_sibling_with_tag(node, "label", "b") will return
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
                                  const char      *tag_value)
{
  if (node != NULL) {

    cs_tree_node_t *sn = node;
    cs_tree_node_t *cn = sn;

    do {

      if (strcmp(cn->name, node->name) == 0) {
        const char *s = cs_tree_node_get_tag(cn, tag);
        if (s != NULL) {
          if (strcmp(s, tag_value) == 0)
            return cn;
        }
      }

      cn = cn->next;
      if (cn == NULL) {
        cn = sn;
        while (cn->prev != NULL)
          cn = cn->prev;
      }

    } while (cn != sn);

  }

  return NULL;
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
  char  _shift[65] = "";
  char  *shift = _shift;
  if (depth > 31)
    BFT_MALLOC(shift, (depth+1)*2+1, char);

  for (int i = 0; i < 2*depth; i++)
    shift[i] = ' ';
  shift[2*depth] = '\0';

  cs_log_printf(log, "%snode_pointer: %p\n", shift, (const void *)node);
  if (node == NULL) {
    if (shift != _shift)
      BFT_FREE(shift);
    return;
  }

  strcat(shift, "  ");

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
      if (node->flag & CS_TREE_NODE_INT)
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

        if (node->flag & CS_TREE_NODE_INT) {

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
 * \brief  Add a node to a tree.
 *
 * This node is located at "path" from the given node, with the path
 * separator indicated by a "/".
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
                 const char      *path)
{
  cs_tree_node_t *_node = cs_tree_get_node(node, path);

  if (_node != NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: node %s already exists.",
              __func__, path);

  return _find_or_create_node(_node, path);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node matching a given path.
 *
 * This node is located at "path" from the given node, with the path
 * separator indicated by a "/".
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
cs_tree_get_node(cs_tree_node_t  *node,
                 const char      *path)
{
  if (node == NULL)
    return NULL;
  if (path == NULL)
    return node;
  if (strlen(path) == 0)
    return node;

  return _find_node(node, path);
}

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
                       const char      *path)
{
  int retval = 0;

  if (node != NULL && path != NULL) {
    cs_tree_node_t  *tn = node;

    if (strlen(path) != 0)
      tn = cs_tree_get_node(node, path);

    while (tn != NULL) {
      retval += 1;
      tn = cs_tree_node_get_next_of_name(tn);
    }

  }

  return retval;
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
 * \brief  Retrieve the pointer to a node matching a given sub-path.
 *
 * This node is located at "path" from the given node or one of its
 * descendants, with the path separator indicated by a "/".
 *
 * In case of multiple nodes sharing the given path, the first such node
 * is returned, using a depth-first search.
 *
 * \param[in]  root      pointer to the root node where we start searching
 * \param[in]  sub_path  string describing the path access
 *
 * \return  pointer to the first matching node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_find_node(cs_tree_node_t  *root,
                  const char      *sub_path)
{
  if (root == NULL)
    return NULL;
  if (sub_path == NULL)
    return root;
  if (strlen(sub_path) == 0)
    return root;

  return _find_sub_node(root, sub_path);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the next node matching a given sub-path
 *         and following a given node in a depth-first order.
 *
 * This node is located at "path" from the given node or one of its
 * descendants, with the path separator indicated by a "/".
 *
 * If current is NULL, this function behaves as \ref cs_tree_find_node.
 *
 * \param[in]  root      pointer to the root node where we start searching
 * \param[in]  current   pointer to the current node
 * \param[in]  sub_path  string describing the path access
 *
 * \return  pointer to the next matching node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_find_node_next(cs_tree_node_t  *root,
                       cs_tree_node_t  *current,
                       const char      *sub_path)
{
  if (root == NULL)
    return NULL;
  if (sub_path == NULL)
    return root;
  if (strlen(sub_path) == 0)
    return root;

  cs_tree_node_t *retval = NULL;

  if (current != NULL) {

    const char *p = sub_path;
    while (*p == '/')
      p++;

    /* Search descendants first */
    if (current->children != NULL)
      retval = _find_sub_node(current->children, p);

    /* Then search siblings or parents */
    if (retval == NULL) {

      cs_tree_node_t *tn = current;

      /* Determine sub-path to check if first part of path matches */
      size_t l = 0;
      while (p[l] != '/' && p[l] != '\0')
        l++;

      while (retval == NULL && tn != root && tn != NULL) {
        if (tn->next != NULL) {
          tn = tn->next;
          if (tn != NULL) {
            /* Search on same level if first part of path matches */
            if (l > 0) {
              if (strncmp(tn->name, p, l) == 0) {
                if (strlen(tn->name) == l) {
                  if (p[l] == '\0')
                    retval = tn;
                  else
                    retval = _find_node(tn, p+l);
                }
              }
            }
            if (retval == NULL)
              retval = _find_sub_node(tn, p);
          }
        }
        else
          tn = tn->parent;
      }
    }

  }

  else
    retval = _find_sub_node(root, sub_path);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a node's descendants matching a given name.
 *
 * This function is similar to \ref cs_tree_find_node,
 * but is simpler (as it assumes a simple name instead of a more general path)
 * and should thus be faster.
 *
 * In case of multiple nodes sharing the given path, the first such node
 * is returned, using a depth-first search.
 *
 * \param[in]  root  pointer to the root node where we start searching
 * \param[in]  name      node name searched for
 *
 * \return  pointer to the first matching node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_find_node_simple(cs_tree_node_t  *root,
                         const char      *name)
{
  if (root == NULL)
    return NULL;
  if (name == NULL)
    return root;
  if (strlen(name) == 0)
    return root;

  return _find_sub_node_simple(root, name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the next node with a given name
 *         and following a given node in a depth-first order.
 *
 * This function is similar to \ref cs_tree_find_node_next,
 * but is simpler (as it assumes a simple name instead of a more general path)
 * and should thus be faster.
 *
 * If current is NULL, this function behaves as \ref cs_tree_find_node.
 *
 * \param[in]  root      pointer to the root node where we start searching
 * \param[in]  current   pointer to the current node
 * \param[in]  name      node name searched for
 *
 * \return  pointer to the next matching node, or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_tree_node_t *
cs_tree_find_node_next_simple(cs_tree_node_t  *root,
                              cs_tree_node_t  *current,
                              const char      *name)
{
  if (root == NULL)
    return NULL;
  if (name == NULL)
    return root;
  if (strlen(name) == 0)
    return root;

  cs_tree_node_t *retval = NULL;

  if (current != NULL) {

    /* Search descendants first */
    if (current->children != NULL)
      retval = _find_sub_node_simple(current->children, name);

    /* Then search siblings or parents */
    if (retval == NULL) {

      cs_tree_node_t *tn = current;

      /* Determine sub-path to check if first part of path matches */

      while (retval == NULL && tn != root && tn != NULL) {
        if (tn->next != NULL) {
          tn = tn->next;
          if (tn != NULL) {
            /* Check for match on same level first */
            if (strcmp(tn->name, name) == 0)
              retval = tn;
            else
              retval = _find_sub_node_simple(tn, name);
          }
        }
        else
          tn = tn->parent;
      }
    }

  }

  else
    retval = _find_sub_node_simple(root, name);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Count a node's descendants matching a given sub-path.
 *
 * These nodes are located at "path" from the given node or one of its
 * descendants, with the path separator indicated by a "/".
 *
 * \param[in]  root      pointer to the root node where we start searching
 * \param[in]  sub_path  string describing the path access
 *
 * \return  number of matching nodes
 */
/*----------------------------------------------------------------------------*/

int
cs_tree_get_sub_node_count(cs_tree_node_t  *root,
                           const char      *sub_path)
{
  int retval = 0;

  for (cs_tree_node_t *tn = cs_tree_find_node(root, sub_path);
       tn != NULL;
       tn = cs_tree_find_node_next(root, tn, sub_path)) {
    retval += 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Count a node's descendants with a given name
 *
 * This function is similar to \ref cs_tree_get_sub_node_count,
 * but is simpler (as it assumes a simple name instead of a more general path)
 * and should thus be faster.
 *
 * \param[in]  root      pointer to the root node where we start searching
 * \param[in]  name      node name searched for
 *
 * \return  number of matching nodes
 */
/*----------------------------------------------------------------------------*/

int
cs_tree_get_sub_node_count_simple(cs_tree_node_t  *root,
                                  const char      *name)
{
  int retval = 0;

  for (cs_tree_node_t *tn = cs_tree_find_node_simple(root, name);
       tn != NULL;
       tn = cs_tree_find_node_next_simple(root, tn, name)) {
    retval += 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_tree_node_t structure starting from the node "root".
 *
 * \param[in] log    indicate which log file to use
 * \param[in] depth  starting depth in the tree
 * \param[in] node   pointer to a cs_tree_node_t to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_dump(cs_log_t                log,
             int                     depth,
             const cs_tree_node_t   *node)
{
  if (depth < 0)
    depth = 0;
  cs_tree_node_dump(log, depth, node);
  if (node == NULL)
    return;

  if (node->children != NULL) { /* There is at least one child */
    cs_tree_node_t  *next_child = node->children;
    while (next_child != NULL) {
      cs_tree_dump(log, depth+1, next_child);
      next_child = next_child->next;
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
