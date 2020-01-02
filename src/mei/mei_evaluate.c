/*!
 * \file mei_evaluate.c
 *
 * \brief Build an interpreter for a mathematical expression
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"
#include "bft_error.h"

#include "mei_node.h"
#include "mei_parser_glob.h"
#include "mei_parser.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mei_evaluate.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Local macro definitions
 *----------------------------------------------------------------------------*/

/*!
 * \brief Default modulo for hash table.
 */

#define HASHSIZE 701

/*=============================================================================
 * Specific pragmas to disable some unrelevant warnings
 *============================================================================*/

/* Globally disable warning on float-comparisons (equality) for GCC and Intel
   compilers as we do it on purpose (infinite precision algorithm). */

#if defined(__GNUC__) && !defined(__ICC)
#pragma GCC diagnostic ignored "-Wfloat-equal"
#elif defined(__ICC)
#pragma warning disable 1572
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy the symbol table pointer into the each node of the interpreter.
 *
 * \param [in] n node of the interpreter
 * \param [in] h table of symbols
 */
/*----------------------------------------------------------------------------*/

static void
_init_symbol_table(mei_node_t    *n,
                   hash_table_t  *h)
{
  int i;

  if (!n) return;

  /* Copy the symbol table pointer into the current node */

  n->ht = h;

  /* Search the next node in the interpreter */

  if (n->flag == FUNC1) {
    _init_symbol_table(n->type->func.op, h);
  }
  else if (n->flag == FUNC2 || n->flag == FUNC3 || n->flag == FUNC4) {
    for (i = 0; i < n->type->funcx.nops; i++)
      _init_symbol_table(n->type->funcx.op[i], h);
  }
  else if (n->flag == OPR) {
    for (i = 0; i < n->type->opr.nops; i++)
      _init_symbol_table(n->type->opr.op[i], h);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of errors related to the use of undefined symbol
 * in the expression.
 *
 * \param [in] p node of an interpreter
 * \return number of undefined symbol
 */
/*----------------------------------------------------------------------------*/

static int
_check_symbol(mei_node_t  *p)
{
  int iok = 0;
  int l;

  if (!p) return 0;

  switch(p->flag) {

  case CONSTANT:
    return 0;

  case ID:
    if (mei_hash_table_lookup(p->ht, p->type->id.i) == NULL) {

      /* bft_printf("Warning: identifier %s is unknown\n", p->type->id.i); */
      /* bft_printf("- line:   %i \n", p->type->id.l); */
      /* bft_printf("- column: %i \n", p->type->id.c); */

      BFT_REALLOC(mei_glob_label_list,  mei_glob_ierr_list+1, char*);
      BFT_REALLOC(mei_glob_line_list,   mei_glob_ierr_list+1, int);
      BFT_REALLOC(mei_glob_column_list, mei_glob_ierr_list+1, int);

      l = strlen("Warning, identifier ")+1;
      BFT_MALLOC(mei_glob_label_list[mei_glob_ierr_list], l, char);
      strncpy(mei_glob_label_list[mei_glob_ierr_list], "Warning: identifier ", l);

      l = l + strlen(p->type->id.i);
      BFT_REALLOC(mei_glob_label_list[mei_glob_ierr_list], l, char);
      strncat(mei_glob_label_list[mei_glob_ierr_list], p->type->id.i, l);

      l = l + strlen(" is unknown.\n");
      BFT_REALLOC(mei_glob_label_list[mei_glob_ierr_list], l, char);
      strncat(mei_glob_label_list[mei_glob_ierr_list], " is unknown.\n", l);

      mei_glob_line_list[mei_glob_ierr_list]   = p->type->id.l;
      mei_glob_column_list[mei_glob_ierr_list] = p->type->id.c;

      mei_glob_ierr_list++;

      return 1;
    }
    return 0;

  case FUNC1:
    if (mei_hash_table_lookup(p->ht, p->type->func.name) == NULL) {

      /* it is normally impossible to arrive here */
      /* because the parser has already checked function names */

      bft_error(__FILE__, __LINE__, 0, _("Error: _check_symbol\n"));

      return 1;
    }
    return _check_symbol(p->type->func.op);

  case FUNC2:
    if (mei_hash_table_lookup(p->ht, p->type->funcx.name) == NULL) {

      /* It is normally impossible to arrive here */
      /* because the parser has already checked function names */

      bft_error(__FILE__, __LINE__, 0, _("Error: _check_symbol\n"));

      return 1;
    }
    iok  = _check_symbol(p->type->funcx.op[0]);
    iok += _check_symbol(p->type->funcx.op[1]);
    return iok;

  case FUNC3:
    bft_error(__FILE__, __LINE__, 0, _("not implemented yet \n"));
    break;

  case FUNC4:
    bft_error(__FILE__, __LINE__, 0, _("not implemented yet \n"));
    break;

  case OPR:

    switch (p->type->opr.oper) {

    case IF:
      iok  = _check_symbol(p->type->opr.op[0]);
      iok += _check_symbol(p->type->opr.op[1]);
      if (p->type->opr.nops > 2)
        iok += _check_symbol(p->type->opr.op[2]);
      return iok;

    case PRINT:
      return _check_symbol(p->type->opr.op[0]);

    case '=':
      mei_hash_table_insert(p->ht,
                            p->type->opr.op[0]->type->id.i,
                            CONSTANT,
                            0,
                            NULL,
                            NULL);
      return _check_symbol(p->type->opr.op[1]);

    case UPLUS:
      return _check_symbol(p->type->opr.op[0]);

    case UMINUS:
      return _check_symbol(p->type->opr.op[0]);

    case '!':
      return _check_symbol(p->type->opr.op[0]);

    default:
      iok  = _check_symbol(p->type->opr.op[0]);
      iok += _check_symbol(p->type->opr.op[1]);
      return iok;
    }
  }

  bft_error(__FILE__, __LINE__, 0, _("Error: _check_symbol\n"));

  return iok;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the evaluation of the expression.
 *
 * \param [in] p node of an interpreter
 * \return value of evaluated expression
 */
/*----------------------------------------------------------------------------*/

static double
_evaluate(mei_node_t  *p)
{
  double t1, t2;
  func1_t f1;
  func2_t f2;

  if (!p) return 0;

  switch(p->flag) {

  case CONSTANT:
    return p->type->con.value;

  case ID:
    return (mei_hash_table_lookup(p->ht, p->type->id.i))->data->value;

  case FUNC1:
    f1 = (mei_hash_table_lookup(p->ht, p->type->func.name))->data->func;
    return f1(_evaluate(p->type->func.op));

  case FUNC2:
    f2 = (mei_hash_table_lookup(p->ht, p->type->funcx.name))->data->f2;
    return f2(_evaluate(p->type->funcx.op[0]), _evaluate(p->type->funcx.op[1]));

  case FUNC3:
    bft_error(__FILE__, __LINE__, 0, _("not implemented\n"));
    break;

  case FUNC4:
    bft_error(__FILE__, __LINE__, 0, _("not implemented\n"));
    break;

  case OPR:

    switch(p->type->opr.oper) {

    case WHILE:
      while(_evaluate(p->type->opr.op[0])) _evaluate(p->type->opr.op[1]);
      return 0;

    case IF:
      if (_evaluate(p->type->opr.op[0]))
        _evaluate(p->type->opr.op[1]);
      else if (p->type->opr.nops > 2)
        _evaluate(p->type->opr.op[2]);
      return 0;

    case PRINT:
      bft_printf("PRINT %f\n", _evaluate(p->type->opr.op[0]));
      return 0;

    case ';':
      _evaluate(p->type->opr.op[0]);
      return _evaluate(p->type->opr.op[1]);

    case '=':
      t2 = _evaluate(p->type->opr.op[1]);
      mei_hash_table_insert(p->ht,
                            p->type->opr.op[0]->type->id.i,
                            CONSTANT,
                            t2,
                            NULL,
                            NULL);
      return 0;

    case UPLUS:
      return _evaluate(p->type->opr.op[0]);

    case UMINUS:
      return -_evaluate(p->type->opr.op[0]);

    case '+':
      return _evaluate(p->type->opr.op[0]) + _evaluate(p->type->opr.op[1]);

    case '-':
      return _evaluate(p->type->opr.op[0]) - _evaluate(p->type->opr.op[1]);

    case '*':
      return _evaluate(p->type->opr.op[0]) * _evaluate(p->type->opr.op[1]);

    case '/':
      t1 = _evaluate(p->type->opr.op[0]);
      t2 = _evaluate(p->type->opr.op[1]);
      if (t2)
        return t1 / t2;
      else
        bft_error(__FILE__, __LINE__, 0, _("Error: floating point exception\n"));
      break;

    case '^':
      return pow(_evaluate(p->type->opr.op[0]), _evaluate(p->type->opr.op[1]));

    case '<':
      return _evaluate(p->type->opr.op[0]) < _evaluate(p->type->opr.op[1]);

    case '>':
      return _evaluate(p->type->opr.op[0]) > _evaluate(p->type->opr.op[1]);

    case '!':
      return ! _evaluate(p->type->opr.op[0]);

    case GE:
      return _evaluate(p->type->opr.op[0]) >= _evaluate(p->type->opr.op[1]);

    case LE:
      return _evaluate(p->type->opr.op[0]) <= _evaluate(p->type->opr.op[1]);

    case NE:
      return _evaluate(p->type->opr.op[0]) != _evaluate(p->type->opr.op[1]);

    case EQ:
      return _evaluate(p->type->opr.op[0]) == _evaluate(p->type->opr.op[1]);

    case AND:
      return _evaluate(p->type->opr.op[0]) && _evaluate(p->type->opr.op[1]);

    case OR:
      return _evaluate(p->type->opr.op[0]) || _evaluate(p->type->opr.op[1]);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Store error message.
 *
 * \param [in] ev interpreter
 */
/*----------------------------------------------------------------------------*/

static void
_manage_error(mei_tree_t  *ev)
{
  int i;
  size_t l;

  assert(ev != NULL);

  ev->errors = mei_glob_ierr_list;

  BFT_MALLOC(ev->labels,  mei_glob_ierr_list, char*);
  BFT_MALLOC(ev->lines,   mei_glob_ierr_list, int);
  BFT_MALLOC(ev->columns, mei_glob_ierr_list, int);

  for (i=0; i<ev->errors; i++) {
    ev->lines[i]   = mei_glob_line_list[i];
    ev->columns[i] = mei_glob_column_list[i];

    l = strlen(mei_glob_label_list[i]) + 1;
    BFT_MALLOC(ev->labels[i], l, char);
    strncpy(ev->labels[i], mei_glob_label_list[i], l);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the symbol \em str exists in the expression stored
 *  in \em ev. Return 0 if the symbol exists in the symbol table, 1 if not.
 *
 * \param [in] ev  interpreter in which we want to know if \em str exists
 * \param [in] str symbol to find
 * \return 0 if the symbol exists in the symbol table, 1 if not
 */
/*----------------------------------------------------------------------------*/

static int
_find_symbol(mei_tree_t  *ev,
             const char  *str)
{
  int i, iok = 0;
  size_t l;

  assert(ev != NULL);
  assert(str != NULL);

  if (!mei_hash_table_lookup(ev->symbol, str)) {

    /* the symbol is not found: this information is stored */

    iok = 1;

    i = ev->errors;
    ev->errors++;

    if (!ev->labels)
      BFT_MALLOC(ev->labels, ev->errors, char*);
    else
      BFT_REALLOC(ev->labels, ev->errors, char*);

    l = strlen(str) + 1;
    BFT_MALLOC(ev->labels[i], l, char);
    strncpy(ev->labels[i], str, l);
  }

  return iok;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a new interpreter. The node member is empty,
 * the string member contains the mathematical expression and
 * the symbol table is initialize with the standard symbols.
 *
 * \param [in] expr string characters of the mathematical expression
 * \return new empty interpreter
 */
/*----------------------------------------------------------------------------*/

mei_tree_t*
mei_tree_new(const char *expr)
{
  mei_tree_t *ev = NULL;
  size_t length;

  if (expr == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: mathematical expression string is empty."));

  BFT_MALLOC(ev, 1, mei_tree_t);

  BFT_MALLOC(ev->symbol, 1, hash_table_t);

  length = strlen(expr)+1;
  BFT_MALLOC(ev->string, length, char);

  strncpy(ev->string, expr, length);

  mei_hash_table_create(ev->symbol, HASHSIZE);
  ev->symbol->n_inter = 1;
  mei_hash_table_init(ev->symbol);

  /* Default initializations */

  ev->errors  = 0;
  ev->columns = NULL;
  ev->lines   = NULL;
  ev->labels  = NULL;
  ev->node    = NULL;

  return ev;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a new interpreter. The node member is empty,
 * the string member contains the mathematical expression and
 * the symbol table is initialize with a existing table, which is shared by
 * multiple interpreter.
 *
 * \param [in] expr         string characters of the mathematical expression
 * \param [in] symbol_table shared table of symbols
 * \return new empty interpreter
 */
/*----------------------------------------------------------------------------*/

mei_tree_t*
mei_tree_new_with_shared_symbols(const char   *const expr,
                                 hash_table_t *const symbol_table)
{
  mei_tree_t *ev = NULL;
  size_t length;

  assert(symbol_table != NULL);

  if (expr == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: mathematical expression string is empty."));

  BFT_MALLOC(ev, 1, mei_tree_t);

  length = strlen(expr)+1;
  BFT_MALLOC(ev->string, length, char);

  strncpy(ev->string, expr, length);

  /* copy the adress of the shared table of symbols and
     incremente the number of interpreters that share this table */

  ev->symbol = symbol_table;
  ev->symbol->n_inter++;

  /* Default initializations */

  ev->errors  = 0;
  ev->columns = NULL;
  ev->lines   = NULL;
  ev->labels  = NULL;
  ev->node    = NULL;

  return ev;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a new table of symbols.
 * The table contains standard mathematical symbols.
 *
 * \return table of symbols
 */
/*----------------------------------------------------------------------------*/

hash_table_t*
mei_table_symbols_new(void)
{
    hash_table_t *ht = NULL;

    BFT_MALLOC(ht, 1, hash_table_t);

    mei_hash_table_create(ht, HASHSIZE);
    ht->n_inter = 0;
    mei_hash_table_init(ht);

    return ht;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call the yacc parser.
 * Return 0 if the parsing is ok, or the number of errors, if errors occurs.
 *
 * \param [in] ev interpreter
 * \return number of parsing errors
 */
/*----------------------------------------------------------------------------*/

int
mei_tree_builder(mei_tree_t *ev)
{
  int i;

  assert(ev != NULL);

  /* Initialize the global variables of the parser */
  /*-----------------------------------------------*/

  /* Initialize the pointer on the node returned by the parser */

  mei_glob_root = NULL;

  /* Begin/end of the expression is set in these global variables */

  mei_glob_string_begin = ev->string;
  mei_glob_string_end   = ev->string + strlen(ev->string);

  /* Initialize line and column counters,
     incrementation is done in scanner.l */

  mei_glob_line   = 1;
  mei_glob_column = 1;

  /* Initialize error counter
     incrementation is done in parser.y (yyerror function) */

  mei_glob_ierr_list = 0;

  /* Parse the expression string and construct the interpreter */
  /*-----------------------------------------------------------*/

  yyparse();

  if (mei_glob_ierr_list) {

    /* Parsing error: copy informations for display */

    _manage_error(ev);

    /* Free memory. Warning: if there is more than one error,
       not all pointers are cleaning */

    mei_free_node(mei_glob_root);
  }
  else {
    /* If the parsing is ok, copy the data in the current interpreter */

    ev->node = mei_glob_root;

    /* For all nodes of the interpreter, copy the symbols table pointer */

    _init_symbol_table(ev->node, ev->symbol);

    /* Check if all new symbols are defined in the expression string */

    mei_glob_ierr_list = _check_symbol(ev->node);

    if (mei_glob_ierr_list) {

      /* Check symbol error: copy informations for display */

      /* if (mei_glob_ierr_list == 1) */
      /*   bft_printf("Warning: there is %i error.\n",  mei_glob_ierr_list); */
      /* else if (mei_glob_ierr_list > 1) */
      /*   bft_printf("Warning: there are %i errors.\n", mei_glob_ierr_list); */

      _manage_error(ev);
    }

  }

  /* Free memory of the parser global variables if necessary */

  for (i=0; i < mei_glob_ierr_list; i++)
    BFT_FREE(mei_glob_label_list[i]);

  BFT_FREE(mei_glob_label_list);
  BFT_FREE(mei_glob_line_list);
  BFT_FREE(mei_glob_column_list);

  return mei_glob_ierr_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inserts a constant (label and value) in the table of symbols
 * associated to an interpreter.
 *
 * \param [in] ev    interpreter
 * \param [in] str   label of the constant
 * \param [in] value value associated to the constant
 */
/*----------------------------------------------------------------------------*/

void
mei_tree_insert(mei_tree_t    *ev,
                const char    *str,
                const double   value)
{
  assert(ev != NULL);
  assert(str != NULL);

  mei_hash_table_insert(ev->symbol,
                        str,
                        CONSTANT,
                        value,
                        NULL,
                        NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inserts a constant (label and value) in a table of symbols.
 *
 * \param [in]  symbol_table table of symbols
 * \param [in]  str          label of the constant
 * \param [in]  value        value associated to the constant
 */
/*----------------------------------------------------------------------------*/

void
mei_symbol_table_insert(hash_table_t  *symbol_table,
                        const char    *str,
                        const double   value)
{
  assert(symbol_table!= NULL);
  assert(str != NULL);

  mei_hash_table_insert(symbol_table,
                        str,
                        CONSTANT,
                        value,
                        NULL,
                        NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the symbol \em str exists in the expression stored
 * in \em ev.
 *
 * \param [in] ev  interpreter in which we want to know if str exists
 * \param [in] str symbol to find
 * \return 0 if the symbol exists in the symbol table, 1 if not
 */
/*----------------------------------------------------------------------------*/

int
mei_tree_find_symbol(mei_tree_t  *ev,
                     const char  *str)
{
  int i;

  assert(ev != NULL);
  assert(str != NULL);

  /* restart error's counter */

  for (i=0; i < ev->errors; i++)
    BFT_FREE(ev->labels[i]);

  BFT_FREE(ev->labels);
  BFT_FREE(ev->lines);
  BFT_FREE(ev->columns);
  ev->errors = 0;

  return _find_symbol(ev, str);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the symbol \em str from a list exists in the expression.
 * The list of missing symbols is stored in \em ev->labels.
 *
 * \param [in]  ev     interpreter
 * \param [in]  size   number of symbols to check
 * \param [in]  symbol list of symbols
 * \return 0 if all symbols exist, otherwise return the number of errors
 */
/*----------------------------------------------------------------------------*/

int
mei_tree_find_symbols(mei_tree_t   *ev,
                      const int     size,
                      const char  **symbol)
{
  int i, iok = 0;

  assert(ev != NULL);

  /* restart error's counter */

  for (i=0; i < ev->errors; i++)
    BFT_FREE(ev->labels[i]);

  BFT_FREE(ev->labels);
  BFT_FREE(ev->lines);
  BFT_FREE(ev->columns);
  ev->errors = 0;

  /* check if each symbol from the list exists */

  for (i=0; i<size; i++)
    iok += _find_symbol(ev, symbol[i]);

  return iok;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns a value of the \em str symbol (variable or constant)
 *  from table of symbols of \em ev interpreter.
 *
 * \param [in]   ev  interpreter
 * \param [in]   str name of the symbol
 * \return value of a symbol
 */
/*----------------------------------------------------------------------------*/

double
mei_tree_lookup(mei_tree_t  *ev,
                const char  *str)
{
  struct item *item = NULL;

  assert(ev != NULL);
  assert(str != NULL);

  item = mei_hash_table_lookup(ev->symbol, str);

  if (!item) {
    bft_error(__FILE__, __LINE__, 0,
              _("Error in mei_tree_lookup function: "
                "%s does not exist in the symbol table\n"), str);
  }
  return item->data->value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluates the expression \em ev :
 *   1) computes all values for variables inside the expression
 *   2) returns a value from the intrepreted expression.
 *
 * \param [in] ev interpreter
 * \return value from the intrepreted expression
 */
/*----------------------------------------------------------------------------*/

double
mei_evaluate(mei_tree_t  *ev)
{
  assert(ev != NULL);
  return _evaluate(ev->node);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory and return NULL.
 *
 * \param [in] ev interpreter
 */
/*----------------------------------------------------------------------------*/

void
mei_tree_destroy(mei_tree_t  *ev)
{
  int i;

  if (ev != NULL) {
    if (ev->symbol->n_inter == 1) { /* symbol tables for a single interpreter */
      mei_hash_table_free(ev->symbol);
      BFT_FREE(ev->symbol);
    }
    else {
      ev->symbol->n_inter--; /* shared symbol tables */
    }
    BFT_FREE(ev->string);
    mei_free_node(ev->node);

    for (i=0; i < ev->errors; i++)
      BFT_FREE(ev->labels[i]);

    BFT_FREE(ev->labels);
    BFT_FREE(ev->lines);
    BFT_FREE(ev->columns);
    BFT_FREE(ev);
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
