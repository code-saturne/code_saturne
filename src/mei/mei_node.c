/*!
 * \file mei_node.c
 *
 * \brief Nodal structure of the interpreter
 */

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
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "mei_node.h"

/*----------------------------------------------------------------------------
 * External global variable declarations for the parser
 *----------------------------------------------------------------------------*/

#include "mei_parser_glob.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a node for a constant.
 *
 * \param [in] value real value of the constant
 *
 * \return built node
 */
/*----------------------------------------------------------------------------*/


mei_node_t*
mei_const_node(const double value)
{
  mei_node_t *node = NULL;

  BFT_MALLOC(node, 1, mei_node_t);

  BFT_MALLOC(node->type, sizeof(const_node_t), node_type_t);

  node->flag = CONSTANT;
  node->ht = NULL;
  node->type->con.value = value;

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a node for a variable.
 *
 * \param [in] variable label of the variable
 *
 * \return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t*
mei_id_node(const char *variable)
{
  mei_node_t *node = NULL;
  size_t length;

  BFT_MALLOC(node, 1, mei_node_t);

  BFT_MALLOC(node->type, sizeof(id_node_t), node_type_t);

  length = strlen(variable)+1;

  BFT_MALLOC(node->type->id.i, length, char);

  node->flag = ID;
  node->ht = NULL;
  strncpy(node->type->id.i, variable, length);
  node->type->id.c = mei_glob_column-length+1;
  node->type->id.l = mei_glob_line;

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a node for a function of a single variable.
 *
 * \param [in] function label of the function
 * \param [in] expr node that represents the variable of the function
 *
 * \return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t*
mei_func_node(const char *function,
              mei_node_t *const expr)
{
  mei_node_t *node = NULL;
  size_t length;
  size_t nodeSize;

  nodeSize = sizeof(func_node_t) + sizeof(mei_node_t);

  BFT_MALLOC(node, 1, mei_node_t);

  BFT_MALLOC(node->type, nodeSize, node_type_t);

  length = strlen(function)+1;
  BFT_MALLOC(node->type->func.name, length, char);

  node->flag = FUNC1;
  node->ht = NULL;
  strncpy(node->type->func.name, function, length);
  node->type->func.op = expr;
  node->type->func.c  = mei_glob_column-length+1;
  node->type->func.l  = mei_glob_line;

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a node for a function of a several variables.
 *
 * \param [in] function label of the function
 * \param [in] nops number of variables
 * \param [in] ...  list of nodes which represent variables of the function
 *
 * \return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t*
mei_funcx_node(const char *function,
               const int nops,
               ...)
{
  va_list ap;
  mei_node_t *node = NULL;
  size_t length;
  size_t nodeSize;
  int i;

  nodeSize = sizeof(func2_node_t) + nops*sizeof(mei_node_t);

  BFT_MALLOC(node, 1, mei_node_t);

  BFT_MALLOC(node->type, nodeSize, node_type_t);

  length = strlen(function)+1;
  BFT_MALLOC(node->type->funcx.name, length, char);
  strncpy(node->type->funcx.name, function, length);

  if (nops == 2) {
    node->flag = FUNC2;
  }
  else if (nops == 3) {
    node->flag = FUNC3;
  }
  else if (nops == 4) {
    node->flag = FUNC4;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "Error: number of arguments for the function is too long\n");

  node->ht = NULL;

  node->type->funcx.nops = nops;
  va_start(ap, nops);
  for (i = 0; i < nops; i++)
    node->type->funcx.op[i] = va_arg(ap, mei_node_t*);
  va_end(ap);

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a node for an operators and its operands.
 *
 * \param [in] oper operator
 * \param [in] nops number of operand
 * \param [in] ...  list of nodes which represent operands
 *
 * \return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t*
mei_opr_node(const int oper,
             const int nops,
             ...)
{
  va_list ap;
  mei_node_t *node = NULL;
  size_t nodeSize;
  int i;

  nodeSize = sizeof(opr_node_t) + nops*sizeof(mei_node_t);

  BFT_MALLOC(node, 1, mei_node_t);

  BFT_MALLOC(node->type, nodeSize, node_type_t);

  node->flag = OPR;
  node->ht = NULL;
  node->type->opr.oper = oper;
  node->type->opr.nops = nops;
  va_start(ap, nops);
  for (i = 0; i < nops; i++)
    node->type->opr.op[i] = va_arg(ap, mei_node_t*);
  va_end(ap);

  return node;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return label of a node.
 *
 * \param [in] n node
 *
 * \return label of a node
 */
/*----------------------------------------------------------------------------*/

char *
mei_label_node(mei_node_t *n)
{
  char *buff;

  if (n->flag == CONSTANT) {
    BFT_MALLOC(buff, 256, char);
    sprintf(buff, "%f", n->type->con.value);
    return buff;
  }
  else if (n->flag == ID) {
    return n->type->id.i;
  }
  else if (n->flag == FUNC1) {
    return n->type->func.name;
  }
  else if (n->flag == FUNC2 || n->flag == FUNC3 || n->flag == FUNC4) {
    return n->type->funcx.name;
  }
  else if (n->flag == OPR) {
    BFT_MALLOC(buff, 256, char);
    sprintf(buff, "operator number: %d", n->type->opr.nops);
    return buff;
  }

  BFT_MALLOC(buff, 256, char);
  sprintf(buff, "%s", " ");
  return buff;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory.
 *
 * \param [in] n node
 */
/*----------------------------------------------------------------------------*/

void
mei_free_node(mei_node_t *n)
{
  int i;

  if (!n) return;

  if (n->flag == ID) {
    BFT_FREE(n->type->id.i);
  }
  else if (n->flag == FUNC1) {
    BFT_FREE(n->type->func.name);
    mei_free_node(n->type->func.op);
  }
  else if (n->flag == FUNC2 || n->flag == FUNC3 || n->flag == FUNC4) {
    BFT_FREE(n->type->funcx.name);
    for (i = 0; i < n->type->funcx.nops; i++)
      mei_free_node(n->type->funcx.op[i]);
  }
  else if (n->flag == OPR) {
    for (i = 0; i < n->type->opr.nops; i++)
      mei_free_node(n->type->opr.op[i]);
  }
  BFT_FREE(n->type);
  BFT_FREE(n);
  return;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
