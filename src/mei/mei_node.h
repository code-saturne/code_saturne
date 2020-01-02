#ifndef __MEI_NODE_H__
#define __MEI_NODE_H__

/*!
 * \file mei_node.h
 *
 * \brief Nodal structure of the interpreter
 */

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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "mei_hash_table.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \brief Constants node
 */

typedef struct
{
  double value;  /*!< value of constant */
} const_node_t;

/*!
 * \brief Identifiers node
 */

typedef struct
{
  char *i;    /*!< label of the identifier */
  int l;      /*!< line number */
  int c;      /*!< column number */
} id_node_t;

/*!
 * \brief Function with single argument node
 */

typedef struct
{
  char *name;              /*!< label of the function */
  int l;                   /*!< line number */
  int c;                   /*!< column number */
  struct _mei_node_t *op;  /*!< variable of the function */
} func_node_t;

/*!
 * \brief Function with two arguments node
 */

typedef struct
{
  char *name;                 /*!< label of the operand  */
  int l;                      /*!< line number  */
  int c;                      /*!< column number  */
  int nops;                   /*!< number of operands */
  struct _mei_node_t *op[];   /*!< list of variable of the function */
} func2_node_t;

/*!
 * \brief Operators node
 */

typedef struct
{
  int oper;                   /*!< operator */
  int nops;                   /*!< number of operands */
  struct _mei_node_t *op[];   /*!< list of operands (expandable) */
} opr_node_t;

/*!
 * \brief Type of a node
 */

typedef union
{
  const_node_t con;         /* constants */
  id_node_t    id;          /* identifiers */
  func_node_t  func;        /* function with one argument  */
  func2_node_t funcx;       /* function with two, three, or four arguments */
  opr_node_t   opr;         /* operators */
} node_type_t;

/*!
 * \brief General node definition
 */

/* union must be last entry in _mei_node_t */
/* because operNodeType may dynamically increase */

struct _mei_node_t
{
  mei_flag_t      flag;       /*!< type of node */
  hash_table_t   *ht;         /*!< pointer of the hash table of the expression */
  node_type_t    *type;       /*!< type of node */
};

/*!
 * \brief General node definition
 */

typedef struct _mei_node_t mei_node_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* Build a node for a constant.
 *
 * param [in] value real value of the constant
 *
 * return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t *
mei_const_node(const double value);

/*----------------------------------------------------------------------------*/
/* Build a node for a variable.
 *
 * param [in] variable label of the variable
 *
 * return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t *
mei_id_node(const char *variable);

/*----------------------------------------------------------------------------*/
/* Build a node for a function of a single variable.
 *
 * param [in] function label of the function
 * param [in] expr node that represents the variable of the function
 *
 * return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t *
mei_func_node(const char *const,
              mei_node_t *const expr);

/*----------------------------------------------------------------------------*/
/* Build a node for a function of a several variables.
 *
 * param [in] function label of the function
 * param [in] nops number of variables
 * param [in] ...  list of nodes which represent variables of the function
 *
 * return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t *
mei_funcx_node(const char *function,
               const int nops,
               ...);

/*----------------------------------------------------------------------------*/
/*
 * Build a node for an operators and its operands.
 *
 * param [in] oper operator
 * param [in] nops number of operand
 * param [in] ...  list of nodes which represent operands
 *
 * return built node
 */
/*----------------------------------------------------------------------------*/

mei_node_t *
mei_opr_node(const int oper,
             const int nops,
             ...);

/*----------------------------------------------------------------------------*/
/*
 * Return label of a node.
 *
 * param [in] n node
 *
 * return label of a node
 */
/*----------------------------------------------------------------------------*/

char *
mei_label_node(mei_node_t *p);

/*----------------------------------------------------------------------------*/
/*
 * Free memory.
 *
 * param [in] n node
 */
/*----------------------------------------------------------------------------*/

void
mei_free_node(mei_node_t *p);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
 #endif /* __cplusplus */

#endif /* __NODE_H__ */

