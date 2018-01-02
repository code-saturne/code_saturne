#ifndef __MEI_EVALUATE_H__
#define __MEI_EVALUATE_H__

/*!
 * \file mei_evaluate.h
 *
 * \brief Build an interpreter for a mathematical expression
 */

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "mei_hash_table.h"
#include "mei_node.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \brief Structure defining an interpreter for a mathematical expression
 */

struct _mei_tree_t {
  char          *string;  /*!< String characters of the mathematical expression */
  int            errors;  /*!< Number of errors                                 */
  int           *columns; /*!< Errors position: list of columns                 */
  int           *lines;   /*!< Errors position: list of lines                   */
  char         **labels;  /*!< Array of the error description                   */
  hash_table_t  *symbol;  /*!< Table of symbols                                 */
  mei_node_t    *node;    /*!< Root node of the interpreter                     */
};

/*!
 * Type definition for an interpreter for a mathematical expression
 */

typedef struct _mei_tree_t mei_tree_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Returns a new interpreter.
 *
 * The node member is empty, the string member contains the mathematical
 * expression and the symbol table is initialize with the standard symbols.
 *
 * parameters:
 *   expr <-- string characters of the mathematical expression
 *
 * returns:
 *   pointer to new interpreter structure.
 *----------------------------------------------------------------------------*/

mei_tree_t *
mei_tree_new(const char  *expr);

/*----------------------------------------------------------------------------
 * Returns a new interpreter. The node member is empty,
 * the string member contains the mathematical expression and
 * the symbol table is initialize with a existing table, which is shared by
 * multiple interpreter.
 *
 * parameters:
 *   expr         <-- string characters of the mathematical expression
 *   symbol_table <-- shared table of symbols
 *
 * returns:
 *   pointer to new interpreter structure.
 *----------------------------------------------------------------------------*/

mei_tree_t *
mei_tree_new_with_shared_symbols(const char    *expr,
                                 hash_table_t  *symbol_table);

/*----------------------------------------------------------------------------
 * Returns a new table of symbols. The table contains standard mathematical
 * symbols
 *----------------------------------------------------------------------------*/

hash_table_t *
mei_table_symbols_new(void);

/*----------------------------------------------------------------------------
 * Calls the yacc parser.
 * Returns 0 if the parsing is ok, or the number of errors, if errors occurs
 *
 * parameters:
 *   ev         <-- interpreter
 *----------------------------------------------------------------------------*/

int
mei_tree_builder(mei_tree_t  *ev);

/*----------------------------------------------------------------------------
 * Inserts a constant (label and value) in the table of symbols associated to
 * an interpreter.
 *
 * parameters:
 *   ev         <-- interpreter
 *   str        <-- label of the constant
 *   value      <-- value associated to the constant
 *----------------------------------------------------------------------------*/

void
mei_tree_insert(mei_tree_t   *ev,
                const char   *str,
                const double  value);

/*----------------------------------------------------------------------------
 * Inserts a constant (label and value) in a table of symbols.
 *
 * parameters:
 *   symbo_table <-- table of symbols
 *   str         <-- label of the constant
 *   value       <-- value associated to the constant
 *----------------------------------------------------------------------------*/

void
mei_symbol_table_insert(hash_table_t  *symbol_table,
                        const char    *str,
                        const double   value);

/*----------------------------------------------------------------------------
 * Checks if the symbol 'str' exists in the expression stored in 'ev'.
 *
 * parameters:
 *   ev         <-- interpreter in which we want to know if str exists
 *   str        <-- symbol to find
 *
 * returns:
 *   0 if the symbol exists in the symbol table, 1 otherwise.
 *----------------------------------------------------------------------------*/

int
mei_tree_find_symbol(mei_tree_t  *ev,
                     const char  *str);

/*----------------------------------------------------------------------------
 * Checks if the symbols in a list exist in the expression.
 *
 * parameters:
 *   ev         <-- interpreter in which we want to know if str exists
 *   str        <-- symbol to find
 *
 * returns:
 *   number of missing symbols
 *----------------------------------------------------------------------------*/

int
mei_tree_find_symbols(mei_tree_t    *ev,
                      const int      size,
                      const char   **symbol);

/*----------------------------------------------------------------------------
 * Returns a value of a symbol (variable or constant) from an interpreter
 * symbol table.
 *
 * parameters:
 *   ev         <-- interpreter
 *   str        <-- name of the symbol
 *----------------------------------------------------------------------------*/

double
mei_tree_lookup(mei_tree_t  *ev,
                const char  *str);

/*----------------------------------------------------------------------------
 * Evaluates the expression :
 *   1) computes all values for variables inside the expression
 *   2) returns a value from the intrepreted expression.
 *
 * parameters:
 *   ev         <-- interpreter
 *----------------------------------------------------------------------------*/

double
mei_evaluate(mei_tree_t  *ev);

/*----------------------------------------------------------------------------
 * Free memory and return NULL.
 *
 * parameters:
 *   ev         <-- interpreter
 *----------------------------------------------------------------------------*/

void
mei_tree_destroy(mei_tree_t  *ev);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MEI_EVALUATE_H__ */

