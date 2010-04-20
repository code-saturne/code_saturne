#ifndef __MEI_NODE_H__
#define __MEI_NODE_H__

/*!
 * \file mei_node.h
 *
 * \brief Nodal structure of the interpreter
 */

/*
  This file is part of the "Mathematical Expression Interpreter" library.

  Copyright (C) 2008-2009  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


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
    double value;               /* value of constant */
} const_node_t;

/*!
 * \brief Identifiers node
 */

typedef struct
{
    char *i;
    int l;
    int c;
} id_node_t;

/*!
 * \brief Function with single argument node
 */

typedef struct
{
    char *name;
    int l;
    int c;
    struct _mei_node_t *op;
} func_node_t;

/*!
 * \brief Function with two arguments node
 */

typedef struct
{
    char *name;
    int l;
    int c;
    int nops;                   /* number of operands */
    struct _mei_node_t *op[];   /* operands (expandable) */
} func2_node_t;

/*!
 * \brief Operators node
 */

typedef struct
{
    int oper;                   /* operator */
    int nops;                   /* number of operands */
    struct _mei_node_t *op[];   /* operands (expandable) */
} opr_node_t;

/*!
 * \brief Type of a node
 */

typedef union
{
    const_node_t con;        /* constants */
    id_node_t    id;         /* identifiers */
    func_node_t  func;       /* function with one argument  */
    func2_node_t funcx;      /* function with two, three, or four arguments */
    opr_node_t   opr;        /* operators */
} node_type_t;

/*!
 * \brief General node definition
 */

/* union must be last entry in _mei_node_t */
/* because operNodeType may dynamically increase */

struct _mei_node_t
{
    mei_flag_t      flag;       /* type of node */
    hash_table_t   *ht;
    node_type_t    *type;
};

/*!
 * \brief General node definition
 */

typedef struct _mei_node_t mei_node_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

mei_node_t *mei_const_node(const double value);

mei_node_t *mei_id_node(char *const variable);

mei_node_t *mei_func_node(char *const name, mei_node_t *const expr);

mei_node_t *mei_funcx_node(char *const function, const int nops, ...);

mei_node_t *mei_opr_node(const int oper, const int nops, ...);

void mei_free_node(mei_node_t *p);

char* mei_label_node(mei_node_t *p);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __NODE_H__ */

