/*============================================================================
 * Expression handling for entity selection based on groups or attributes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_selector_postfix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Allowed operand types (bit mask)
 */

#define _OPERAND_INT          (1 << 0)
#define _OPERAND_DOUBLE       (1 << 1)
#define _OPERAND_STRING       (1 << 2)
#define _OPERAND_GEOMETRIC    (1 << 3)

/* Base stack length */

#define BASE_STACK_SIZE 32

/* Geometric operation macros*/

#define _DOT_PRODUCT(v1, v2) \
  (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define _MODULE(v) \
  sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Types of arguments or operands for an operator
 *----------------------------------------------------------------------------*/

typedef enum {

  OT_L_PAREN,          /* Left parenthesis */
  OT_R_PAREN,          /* Right parenthesis */
  OT_UNARY,            /* Unary function (arg to the right) */
  OT_BINARY,           /* Binary function (args left and right) */
  OT_FUNCTION,         /* Built-in functions (args between brackets) */
  OT_COORD_CONDITION,  /* Condition on coordinates (<, <=, >, >=) */
  OT_DEFINITION,       /* Sub-expression definition (not implemented) */
  OT_MATH_FUNCTION,    /* Mathematical function (not implemented) */
  OT_NONE              /* Not an operaror */

} _operator_type_t;

/*----------------------------------------------------------------------------
 * List of operator codes and names
 *----------------------------------------------------------------------------*/

typedef enum {

  OC_L_PAREN,
  OC_R_PAREN,

  OC_NOT,
  OC_AND,
  OC_OR,
  OC_XOR,

  OC_ALL,
  OC_NO_GROUP,
  OC_RANGE,

  OC_NORMAL,
  OC_PLANE,
  OC_BOX,
  OC_CYLINDER,
  OC_SPHERE,

  OC_GT,
  OC_LT,
  OC_GE,
  OC_LE,

  OC_NONE

} _operator_code_t;

static const char *_operator_name[] = {"(", ")",
                                       "not", "and", "or", "xor",
                                       "all", "no_group", "range",
                                       "normal", "plane", "box",
                                       "cylinder", "sphere",
                                       ">", "<", ">=", "<=",
                                       "not_an_operator"};

/*----------------------------------------------------------------------------
 * Types of postfix elements
 *----------------------------------------------------------------------------*/

typedef enum {

  PF_OPCODE,
  PF_GROUP_ID,
  PF_ATTRIBUTE_ID,
  PF_INT,
  PF_FLOAT

} _postfix_type_t;

/*----------------------------------------------------------------------------
 * Definition of a tokenized expression
 *----------------------------------------------------------------------------*/

typedef struct {

  int     n_tokens;               /* Number of tokens in expression */
  int    *infix_id;               /* Starting id in infix for each token */
  int    *token_id;               /* Starting id for each token */
  bool   *protected;              /* Indicates if a token was protected
                                     by quotes or a backslash character
                                     (and thus cannot be a keyword) */

  int     size;                   /* Current string size */
  int     max_size;               /* Maximum string size */

  char   *tokens;                 /* Tokens array */

} _tokenized_t;

/*----------------------------------------------------------------------------
 * Operator definition
 *----------------------------------------------------------------------------*/

typedef struct {

  _operator_code_t   code;        /* Operator code number */
  _operator_type_t   type;        /* Operator type */

  int                priority;    /* Priority */

  char               name[16];    /* Operator name string */

} _operator_t;

/*----------------------------------------------------------------------------
 * Parser object
 *----------------------------------------------------------------------------*/

typedef struct {

  int           n_operators;       /* Number of possible operators */
  _operator_t  *operators;         /* Array of allowed operators */

  int           n_keywords;        /* Total number of keywords */
  int          *keyword_op_id;     /* Operator id for each keyword */
  char        **keyword;           /* Pointer to individual keywords */

  size_t        keywords_size;     /* Size of keyword buffer */
  char         *keywords;          /* Pointer to keyword buffer */

} _parser_t;

/*----------------------------------------------------------------------------
 * Stack of parsed operator elements
 *----------------------------------------------------------------------------*/

typedef struct {

  const _operator_t  *op;        /* Pointer to operator definition */
  int                 token_id;  /* Associated token's id */

} _stack_entry_t;

typedef struct {

  size_t  size;                                /* Current memory size */
  size_t  max_size;                            /* Maximum memory size */

  _stack_entry_t  _elements[BASE_STACK_SIZE];  /* Base working stack */
  _stack_entry_t  *elements;                   /* Pointer to _elements
                                                  or allocated */

} _stack_t;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Postfix expression object
 *----------------------------------------------------------------------------*/

struct _fvm_selector_postfix_t {

  bool    coords_dependency;      /* Does evaluation require coordinates ? */
  bool    normals_dependency;     /* Does evaluation require normals ? */

  size_t  size;                   /* Current memory size */
  size_t  max_size;               /* Maximum memory size */

  char   *infix;                  /* Copy of original infix expression */
  unsigned char  *elements;       /* Contents array */

  int     n_missing_operands;     /* Number of operands corresponding to
                                     missing group names */
  char **missing_operand;         /* Array of missing group names */

};

/*============================================================================
 * Global variables
 *============================================================================*/

/* Reference count and pointer to global parser object */

static int _n_parser_references = 0;
static _parser_t *_parser = NULL;

/* Compute sizes so as to have aligned data */

static const size_t _postfix_type_size
  = sizeof(size_t)*(  (sizeof(_postfix_type_t)/sizeof(size_t))
                    + (sizeof(_postfix_type_t)%sizeof(size_t) ? 1 : 0));
static const size_t _postfix_opcode_size
  = sizeof(size_t)*(  (sizeof(_operator_code_t)/sizeof(size_t))
                    + (sizeof(_operator_code_t)%sizeof(size_t) ? 1 : 0));
static const size_t _postfix_int_size
  = sizeof(size_t)*(  (sizeof(int)/sizeof(size_t))
                    + (sizeof(int)%sizeof(size_t) ? 1 : 0));
static const size_t _postfix_float_size
  = sizeof(size_t)*(  (sizeof(double)/sizeof(size_t))
                    + (sizeof(double)%sizeof(size_t) ? 1 : 0));

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add an operator definition to a parser.
 *
 * parameters:
 *   this_parser   <-> Parser to which operator should be added
 *   name          <-- Operator name (truncated at 15 characters)
 *   operator_code <-- Unique operator code
 *   operator_type <-- Operator type (unary, binary, function)
 *   priority      <-- Operator priority
 *   n_keywords    <-- Number of associated keywords
 *   keywords      <-- Array of associated keywords
 *----------------------------------------------------------------------------*/

static void
_add_operator(_parser_t          *this_parser,
              const char         *name,
              _operator_code_t    operator_code,
              _operator_type_t    operator_type,
              int                 priority,
              int                 n_keywords,
              const char        **keywords)
{
  int i;

  /* Initial checks */

  assert(this_parser != NULL);
  assert(n_keywords > 0);

  /* Resize structures */

  BFT_REALLOC(this_parser->operators,
              this_parser->n_operators + 1,
              _operator_t);

  if (n_keywords > 0) {

    size_t keywords_size = 0;

    for (i = 0; i < n_keywords; i++)
      keywords_size += (strlen(keywords[i]) + 1);

    BFT_REALLOC(this_parser->keyword_op_id,
                this_parser->n_keywords + n_keywords,
                int);

    BFT_REALLOC(this_parser->keyword,
                this_parser->n_keywords + n_keywords,
                char *);

    BFT_REALLOC(this_parser->keywords,
                this_parser->keywords_size + keywords_size,
                char);

  }

  /* Fill structures */

  i = this_parser->n_operators;

  this_parser->operators[i].code = operator_code;
  this_parser->operators[i].priority = priority;
  this_parser->operators[i].type = operator_type;

  strncpy(this_parser->operators[i].name, name, 15);
  this_parser->operators[i].name[15] = '\0';

  for (i = 0; i < n_keywords; i++) {

    size_t  l = strlen(keywords[i]) + 1;

    this_parser->keyword_op_id[this_parser->n_keywords]
      = this_parser->n_operators;

    memcpy(this_parser->keywords + this_parser->keywords_size,
           keywords[i],
           l);

    this_parser->n_keywords += 1;
    this_parser->keywords_size += l;

  }

  this_parser->n_operators += 1;
}

/*----------------------------------------------------------------------------
 * Creation of a parser object.
 *
 * parameters:
 *
 * returns:
 *   pointer to new parser
 *----------------------------------------------------------------------------*/

static _parser_t *
_parser_create(void)
{
  int i;

  const char *kw_l_paren[] = {"("};
  const char *kw_r_paren[] = {")"};

  const char *kw_not[] = {"not", "!", "!=", "NOT"};
  const char *kw_and[] = {"and", "&", "&&", "AND"};
  const char *kw_or[] =  {"or", "|", "||", ",", ";", "OR"};
  const char *kw_xor[] = {"xor", "^", "XOR"};
  const char *kw_all[] = {"all", "ALL"};
  const char *kw_ngr[] = {"no_group", "NO_GROUP"};
  const char *kw_rng[] = {"range", "RANGE"};

  const char *kw_nrm[] = {"normal", "NORMAL"};
  const char *kw_pln[] = {"plane", "PLANE"};
  const char *kw_box[] = {"box", "BOX"};
  const char *kw_cyl[] = {"cylinder", "CYLINDER"};
  const char *kw_sph[] = {"sphere", "SPHERE"};

  const char *kw_gt[] = {">"};
  const char *kw_lt[] = {"<"};
  const char *kw_ge[] = {">="};
  const char *kw_le[] = {"<="};

  _parser_t *p = NULL;

  /* Create base structure */

  BFT_MALLOC(p, 1, _parser_t);

  p->n_operators = 0;
  p->operators = NULL;

  p->n_keywords = 0;
  p->keyword_op_id = NULL;
  p->keyword = NULL;

  p->keywords_size = 0;
  p->keywords = NULL;

  /* Add operator definitions */

  _add_operator(p, "(", OC_L_PAREN, OT_L_PAREN, 0, 1, kw_l_paren);
  _add_operator(p, ")", OC_R_PAREN, OT_R_PAREN, 0, 1, kw_r_paren);

  _add_operator(p, "not", OC_NOT, OT_UNARY, 3, 4, kw_not);
  _add_operator(p, "and", OC_AND, OT_BINARY, 2, 4, kw_and);
  _add_operator(p, "or", OC_OR, OT_BINARY, 1, 6, kw_or);
  _add_operator(p, "xor", OC_XOR, OT_BINARY, 1, 3, kw_xor);

  _add_operator(p, "all", OC_ALL, OT_FUNCTION, 4, 2, kw_all);
  _add_operator(p, "no_group", OC_NO_GROUP, OT_FUNCTION, 4, 2, kw_ngr);
  _add_operator(p, "range", OC_RANGE, OT_FUNCTION, 4, 2, kw_rng);

  _add_operator(p, "normal", OC_NORMAL, OT_FUNCTION, 4, 2, kw_nrm);
  _add_operator(p, "plane", OC_PLANE, OT_FUNCTION, 4, 2, kw_pln);
  _add_operator(p, "box", OC_BOX, OT_FUNCTION, 4, 2, kw_box);
  _add_operator(p, "cylinder", OC_CYLINDER, OT_FUNCTION, 4, 2, kw_cyl);
  _add_operator(p, "sphere", OC_SPHERE, OT_FUNCTION, 4, 2, kw_sph);

  _add_operator(p, ">", OC_GT, OT_COORD_CONDITION, 4, 1, kw_gt);
  _add_operator(p, "<", OC_LT, OT_COORD_CONDITION, 4, 1, kw_lt);
  _add_operator(p, ">=", OC_GE, OT_COORD_CONDITION, 4, 1, kw_ge);
  _add_operator(p, "<=", OC_LE, OT_COORD_CONDITION, 4, 1, kw_le);

  /* Build keyword pointers */

  p->keyword[0] = p->keywords;

  for (i = 1; i < p->n_keywords; i++)
    p->keyword[i]
      = p->keyword[i-1] + strlen(p->keyword[i-1]) + 1;

  /* Return parser object */

  return p;
}

/*----------------------------------------------------------------------------
 * Destruction of a parser.
 *
 * parameters:
 *   this_parser <-> parser to destroy
 *----------------------------------------------------------------------------*/

static void
_parser_destroy(_parser_t  **this_parser)
{
  if (*this_parser != NULL) {

    BFT_FREE((*this_parser)->operators);
    BFT_FREE((*this_parser)->keyword_op_id);
    BFT_FREE((*this_parser)->keyword);
    BFT_FREE((*this_parser)->keywords);

    BFT_FREE(*this_parser);

  }
}

/*----------------------------------------------------------------------------
 * Print content of a parser object.
 *
 * parameters:
 *   this_parser <-- parser object
 *----------------------------------------------------------------------------*/

static void
_parser_dump(const _parser_t  *this_parser)
{
  int i;

  const char *type_name[] = {"(", ")", "unary", "binary", "function",
                             "coord condition", "definition", "math_func"};

  if (this_parser == NULL)
    return;

  /* Global indicators */
  /*-------------------*/

  bft_printf("\n"
             "Parser:\n\n"
             "Number of operators:  %d\n"
             "Number of keywords:   %d\n\n",
             this_parser->n_operators,
             this_parser->n_keywords);


  /* Operators and related part of keywords array */

  if (this_parser->n_operators > 0)
    bft_printf("Operators:\n"
               "    id  | name     | code | pri | type  \n"
               "    ------------------------------------\n");

  for (i = 0; i < this_parser->n_operators; i++) {
    const _operator_t *op = this_parser->operators + i;
    bft_printf("   %4d | %8s | %4d | %3d | %s\n",
               i, op->name, op->code, op->priority, type_name[op->type]);
  }

  if (this_parser->n_keywords > 0)
    bft_printf("\n"
               "Keywords:\n"
               "    id  | op_id | name\n"
               "    ------------------\n");

  for (i = 0; i < this_parser->n_keywords; i++)
    bft_printf("   %4d | %5d | %s\n",
               i, this_parser->keyword_op_id[i], this_parser->keyword[i]);

  bft_printf("\n");
}

/*----------------------------------------------------------------------------
 * Tokenize an infix expression string.
 *
 * parameters:
 *   infix    <-- string parsed
 *
 * returns:
 *   a tokenized_t structure corresponding to the infix string
 *----------------------------------------------------------------------------*/

static _tokenized_t
_tokenize(const char  *infix)
{
  int i, l, tok_len, start_quote_id;
  int  protected; /* 0 for unprotected, 1 for protected char, 2 for substring,
                     3 for protected char in substring */

  _tokenized_t te;

  /* Initialization */

  te.n_tokens = 0;
  te.infix_id = NULL;
  te.token_id = NULL;
  te.protected = NULL;

  te.size = 0;
  te.max_size = 0;
  te.tokens = NULL;

  if (infix == NULL)
    return te;

  /* Allocate tokenization structure */

  l = strlen(infix);

  te.max_size = l*2 + 1; /* Allows for insertion of NULL chars after
                            each character (worst case), so no size
                            increase will be necessary */

  BFT_MALLOC(te.infix_id, l, int);
  BFT_MALLOC(te.token_id, l, int);
  BFT_MALLOC(te.protected, l, bool);
  BFT_MALLOC(te.tokens, te.max_size, char);

  for (i = 0; i < l; i++)
    te.protected[i] = false;

  i = 0;                 /* String position marker */
  start_quote_id = l;    /* Start position marker for quoted strings
                            (unset if j == j, set if j < l) */

  protected = 0;
  tok_len = 0;

  while (i < l) {

    char c = infix[i];

    /* Regular case, where previous character was not protected */

    if (protected == 0) {

      /* Protection character */

      if (c == '\\') {
        protected = 1;
        te.protected[te.n_tokens] = true;
      }

      /* Fully protected string */

      else if (c == '"' || c == '\'') {
        protected = 2;
        start_quote_id = i;
        te.protected[te.n_tokens] = true;
      }

      /* Whitespace */

      else if (c == ' ' || c == '\t' || c == '\n'|| c == '\r') {
        if (tok_len > 0) { /* Finish previous token */
          te.token_id[te.n_tokens] = te.size;
          te.tokens[te.size + tok_len] = '\0';
          te.n_tokens += 1;
          te.size += tok_len+1;
          tok_len = 0;
        }
      }

      /* Other punctuation characters */

      else if (   c == '(' || c == ')'
               || c == '[' || c == ']'
               || c == ',' || c == ';'
               || c == '!' || c == '^' || c == '|' || c == '&'
               || c == '=' || c == '<' || c == '>') {

        char c2 = '\0';

        if (tok_len > 0) { /* Finish previous token */
          te.token_id[te.n_tokens] = te.size;
          te.tokens[te.size + tok_len] = '\0';
          te.n_tokens += 1;
          te.size += tok_len+1;
          tok_len = 0;
        }

        /* Add token, handling possible 2-character punctuation */

        te.token_id[te.n_tokens] = te.size;
        te.tokens[te.size] = c;

        if (i+1 < l)
          c2 = infix[i+1];

        if (   ((c == '=' || c == '<' || c == '>') && c2 == '=')
            || (c == '!' && c2 == '=')
            || (c == '|' && c2 == '|')
            || (c == '&' && c2 == '&')) {
          i += 1;
          te.tokens[te.size + 1] = '=';
          te.tokens[te.size + 2] = '\0';
          te.size += 3;
        }
        else {
          te.tokens[te.size + 1] = '\0';
          te.size += 2;
        }
        te.infix_id[te.n_tokens] = i;
        te.n_tokens += 1;

      }

      /* Regular characters (token) */

      else {
        te.tokens[te.size + tok_len] = infix[i];
        if (tok_len == 0)
          te.infix_id[te.n_tokens] = i;
        tok_len++;
      }

    }

    /* Cases where previous character was protected */

    else if (protected == 1) {

      protected = 0;
      te.tokens[te.size + tok_len] = infix[i];
      if (tok_len == 0)
        te.infix_id[te.n_tokens] = i;
      tok_len++;

    }

    else if (protected == 2) {

      /* Single protection character */

      if (c == '\\')
        protected = 3;

      /* End of string protection */

      else if (c == infix[start_quote_id]) {
        protected = 0;
        start_quote_id = l;
      }

      else {
        te.tokens[te.size + tok_len] = infix[i];
        if (tok_len == 0)
          te.infix_id[te.n_tokens] = i;
        tok_len++;
      }

    }

    else { /* if (protected == 3) */

      te.tokens[te.size + tok_len] = infix[i];
      if (tok_len == 0)
        te.infix_id[te.n_tokens] = i;
      tok_len++;
      protected = 2;

    }

    i+= 1;

  } /* End of loop in infix string characters */

  if (tok_len > 0) { /* Finish previous token */
    te.token_id[te.n_tokens] = te.size;
    te.tokens[te.size + tok_len] = '\0';
    te.n_tokens += 1;
    te.size += tok_len+1;
    tok_len = 0;
  }

  if (protected == 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Error tokenizing expression:\n"
                "%s\n"
                "Missing character after \\\n"),
              infix);
  else if (protected >= 2)
    bft_error(__FILE__, __LINE__, 0,
              _("Error tokenizing expression:\n"
                "%s\n"
                "Missing closing quote for subexpression:\n"
                "%s\n"),
              infix, infix + start_quote_id);

  /* Resize to adjusted size */

  BFT_REALLOC(te.infix_id, te.n_tokens, int);
  BFT_REALLOC(te.token_id, te.n_tokens, int);
  BFT_REALLOC(te.protected, te.n_tokens, bool);
  BFT_REALLOC(te.tokens, te.size, char);

  /* Return tokenization structure */

  return te;
}

/*----------------------------------------------------------------------------
 * Empty a tokenized expression string.
 *
 * parameters:
 *   te    <-> tokenized expression
 *----------------------------------------------------------------------------*/

static void
_empty_tokenized(_tokenized_t  *te)
{
  te->n_tokens = 0;
  te->size = 0;
  te->max_size = 0;
  BFT_FREE(te->infix_id);
  BFT_FREE(te->token_id);
  BFT_FREE(te->protected);
  BFT_FREE(te->tokens);
}

/*----------------------------------------------------------------------------
 * Dump a tokenized expression string.
 *
 * parameters:
 *   infix <-- string parsed
 *   te    <-- tokenized expression
 *----------------------------------------------------------------------------*/

static void
_dump_tokenized(const char          *infix,
                const _tokenized_t   te)
{
  int i;

  bft_printf("\n"
             "Tokenization:\n\n"
             "Infix:\n%s\n"
             "Tokens: %d\n",
             infix, te.n_tokens);

  for (i = 0; i < te.n_tokens; i++) {
    bft_printf("  %3d: %-20s", i, te.tokens + te.token_id[i]);
    bft_printf(" (%d bytes from infix start", te.infix_id[i]);
    if (te.protected[i] != 0)
      bft_printf(", protected)\n");
    else
      bft_printf(")\n");
  }
}

/*----------------------------------------------------------------------------
 * Initialize a stack
 *
 * parameters:
 *   stack <-> stack
 *----------------------------------------------------------------------------*/

static void
_stack_initialize(_stack_t *stack)
{
  stack->size = 0;
  stack->max_size = BASE_STACK_SIZE;
  stack->elements = stack->_elements;
}

/*----------------------------------------------------------------------------
 * Empty a stack
 *
 * parameters:
 *   stack <-> stack
 *----------------------------------------------------------------------------*/

static void
_stack_empty(_stack_t *stack)
{
  stack->size = 0;
  stack->max_size = BASE_STACK_SIZE;
  if (stack->elements != stack->_elements) {
    BFT_FREE(stack->elements);
    stack->elements = stack->_elements;
  }
}

/*----------------------------------------------------------------------------
 * Push operator info to the top of the stack
 *
 * parameters:
 *   stack    <-> stack
 *   op       <-- Pointer to operator definition
 *   token_id <-- associated token id
 *----------------------------------------------------------------------------*/


inline static void
_stack_push(_stack_t           *stack,
            const _operator_t  *op,
            int                 token_id)
{
  _stack_entry_t  *e;

  assert(stack != NULL);

  if (stack->size >=  stack->max_size) {
    stack->max_size *= 2;
    if (stack->max_size > BASE_STACK_SIZE) {
      if (stack->elements == stack->_elements) {
        BFT_MALLOC(stack->elements, stack->max_size, _stack_entry_t);
        memcpy(stack->elements, stack->_elements, stack->size*sizeof(_stack_entry_t));
      }
      else
        BFT_REALLOC(stack->elements, stack->max_size, _stack_entry_t);
    }
  }

  e = stack->elements + stack->size;

  e->op = op;
  e->token_id = token_id;

  stack->size++;
}

/*----------------------------------------------------------------------------
 * Pop an entry from the top of the stack
 *
 * parameters:
 *   stack <-> stack
 *   id    --> optional id associated with value, or NULL
 *
 * returns:
 *   the value on the top of the stack
 *----------------------------------------------------------------------------*/

inline static _stack_entry_t
_stack_pop(_stack_t  *stack)
{
  assert(stack != NULL);

  if (stack->size == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Trying to pop an entry from empty stack."));

  return stack->elements[--(stack->size)];
}

/*----------------------------------------------------------------------------
 * Return the number of elements in a stack
 *
 * parameters:
 *   stack <-- stack
 *
 * returns:
 *   number of elements in the stack
 *----------------------------------------------------------------------------*/

inline static size_t
_stack_size(_stack_t  *stack)
{
  return stack->size;
}

/*----------------------------------------------------------------------------
 * Return pointer to operator definition for the entry at the top of
 * a stack (without modifying the stack)
 *
 * parameters:
 *   stack <-> stack
 *
 * returns:
 *   the value on the top of the stack
 *----------------------------------------------------------------------------*/

inline static const _operator_t *
_stack_top(const _stack_t  *stack)
{
  assert(stack != NULL);

  if (stack->size == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Querying a top value from empty stack."));

  return (stack->elements + stack->size - 1)->op;
}

/*----------------------------------------------------------------------------
 * Initialize a postfix expression
 *
 * parameters:
 *   infix <-> initially infix expression
 *
 * returns:
 *   new postfix structure
 *----------------------------------------------------------------------------*/

static fvm_selector_postfix_t *
_postfix_create(const char *infix)
{
  size_t i;
  size_t size = strlen(infix);
  fvm_selector_postfix_t *pf = NULL;

  BFT_MALLOC(pf, 1, fvm_selector_postfix_t);

  pf->coords_dependency = false;
  pf->normals_dependency = false;

  pf->size = 0;
  pf->max_size = size * sizeof(size_t);

  BFT_MALLOC(pf->infix, size + 1, char);
  strcpy(pf->infix, infix);

  BFT_MALLOC(pf->elements, pf->max_size, unsigned char);
  for (i = 0; i < pf->max_size; pf->elements[i++] = '\0');

  pf->n_missing_operands = 0;
  pf->missing_operand = NULL;

  return pf;
}

/*----------------------------------------------------------------------------
 * Destroy a postfix expression
 *
 * parameters:
 *   pf <-> pointer to structure pointer that should be destroyed
 *----------------------------------------------------------------------------*/

static void
_postfix_destroy(fvm_selector_postfix_t **pf)
{
  fvm_selector_postfix_t *_pf = *pf;

  if (*pf != NULL) {

    BFT_FREE(_pf->infix);
    BFT_FREE(_pf->elements);

    if (_pf->n_missing_operands > 0) {
      int i;
      for (i = 0; i < _pf->n_missing_operands; i++)
        BFT_FREE(_pf->missing_operand[i]);
      BFT_FREE(_pf->missing_operand);
    }

    BFT_FREE(_pf);
    *pf = _pf;
  }
}

/*----------------------------------------------------------------------------
 * Grow postfix string for future additions
 *
 * parameters:
 *   pf     <-> posfix string
 *   size_t <-- new_size
 *----------------------------------------------------------------------------*/

static void
_postfix_grow(fvm_selector_postfix_t  *pf,
              size_t                   new_size)
{
  size_t i;
  size_t old_max_size = pf->max_size;
  if (old_max_size*2 > new_size)
    pf->max_size *= 2;
  else
    pf->max_size = new_size;

  BFT_REALLOC(pf->elements, pf->max_size, unsigned char);
  for (i = old_max_size; i < pf->max_size; pf->elements[i++] = '\0');
}

/*----------------------------------------------------------------------------
 * Adjust max size of postfix structure to current size
 *
 * parameters:
 *   pf     <-> posfix string
 *----------------------------------------------------------------------------*/

static void
_postfix_adjust(fvm_selector_postfix_t  *pf)
{
  if (pf->size != pf->max_size) {
    pf->max_size = pf->size;
    BFT_REALLOC(pf->elements, pf->max_size, unsigned char);
  }
}

/*----------------------------------------------------------------------------
 * Add an operator code to a postfix expression
 *
 * parameters:
 *   pf   <-> pointer to postfix structure
 *   code <-- operator code
 *----------------------------------------------------------------------------*/

static void
_postfix_add_opcode(fvm_selector_postfix_t  *pf,
                    _operator_code_t         code)
{
  size_t add_size = _postfix_type_size + _postfix_opcode_size;

  /* Update postfix */

  if (pf->size + add_size > pf->max_size)
    _postfix_grow(pf, pf->size + add_size);

  *((_postfix_type_t *)(pf->elements + pf->size)) = PF_OPCODE;
  *((_operator_code_t *)(pf->elements + pf->size + _postfix_type_size))
    = code;

  pf->size += add_size;
}

/*----------------------------------------------------------------------------
 * Add an integer to a postfix expression
 *
 * parameters:
 *   pf   <-> pointer to postfix structure
 *   val  <-- integer value
 *   type <-- integer subtype (PF_GROUP_ID, PF_ATTRIBUTE_ID, or PF_INT)
 *----------------------------------------------------------------------------*/

static void
_postfix_add_int(fvm_selector_postfix_t  *pf,
                 int                      val,
                 _postfix_type_t          type)
{
  size_t add_size = _postfix_type_size + _postfix_int_size;

  /* Update postfix */

  if (pf->size + add_size > pf->max_size)
    _postfix_grow(pf, pf->size + add_size);

  *((_postfix_type_t *)(pf->elements + pf->size)) = type;
  *((int *)(pf->elements + pf->size + _postfix_type_size)) = val;

  pf->size += add_size;
}

/*----------------------------------------------------------------------------
 * Add an float to a postfix expression
 *
 * parameters:
 *   pf  <-> pointer to postfix structure
 *   val <-- integer value
 *----------------------------------------------------------------------------*/

static void
_postfix_add_float(fvm_selector_postfix_t  *pf,
                   double                   val)
{
  size_t add_size = _postfix_type_size + _postfix_float_size;

  /* Update postfix */

  if (pf->size + add_size > pf->max_size)
    _postfix_grow(pf, pf->size + add_size);

  *((_postfix_type_t *)(pf->elements + pf->size)) = PF_FLOAT;
  *((double *)(pf->elements + pf->size + _postfix_type_size)) = val;

  pf->size += add_size;
}

/*----------------------------------------------------------------------------
 * Add missing operand to postfix string structure
 *
 * parameters:
 *   pf      <-> posfix string
 *   missing <-- name of missing operand
 *----------------------------------------------------------------------------*/

static void
_postfix_add_missing(fvm_selector_postfix_t  *pf,
                     const char              *missing)
{
  int n = pf->n_missing_operands;

  BFT_REALLOC(pf->missing_operand, n + 1, char *);
  BFT_MALLOC(pf->missing_operand[n], strlen(missing) + 1, char);
  strcpy(pf->missing_operand[pf->n_missing_operands], missing);
  pf->n_missing_operands++;
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
 * Check if a string defines a floating-point number and scan its value
 *
 * parameters:
 *   str   <-- string parsed
 *   value --> floating-point conversion
 *
 * returns:
 *   true if the string defines a floating-point number, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_is_float(const char  *str,
          double      *value)
{
  float _value;
  int retcode, flt_len;

  *value = 0.0;
  retcode = (bool)(sscanf(str, "%f%n", &_value, &flt_len));

  if (retcode) {
    if (flt_len != (int)strlen(str))
      retcode = false;
  }

  if (retcode)
    *value = _value;

  return retcode;
}

/*----------------------------------------------------------------------------
 * Handle a parsing error.
 *
 * parameters:
 *   err_str      <-- error string
 *   valid_syntax <-- optional valid syntax info, or NULL
 *   infix        <-- string parsed
 *   te           <-- tokenized expression
 *   token_id     <-- token id (or -1 if unknown)
 *   stack        <-> stack (emptied)
 *   postfix      <-> postfix expression (destroyed)
 *----------------------------------------------------------------------------*/

static void
_parse_error(const char               *err_str,
             const char               *valid_syntax,
             const char               *infix,
             const _tokenized_t       *te,
             int                       token_id,
             _stack_t                 *stack,
             fvm_selector_postfix_t  **postfix)
{
  int infix_pos = -1;

  if (token_id > -1)
    infix_pos = te->infix_id[token_id];

  _stack_empty(stack);
  _postfix_destroy(postfix);

  if (getenv("FVM_SELECTOR_DEBUG")) {
    _parser_dump(_parser);
    _dump_tokenized(infix, *te);
  }

  if (infix_pos > -1) {

    int i;
    char *infix_string_marker = NULL;

    BFT_MALLOC(infix_string_marker, infix_pos + 2, char);
    for (i = 0; i < infix_pos; i++)
      infix_string_marker[i] = ' ';
    infix_string_marker[infix_pos] = '^';
    infix_string_marker[infix_pos + 1] = '\0';

    if (valid_syntax != NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Error parsing expression:\n"
                  "%s\n"
                  "%s\n"
                  "%s\n\n"
                  "Valid (expected) syntax:\n\n"
                  "%s"),
                infix, infix_string_marker, err_str, valid_syntax);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error parsing expression:\n"
                  "%s\n"
                  "%s\n"
                  "%s"),
                infix, infix_string_marker, err_str);

    BFT_FREE(infix_string_marker);
  }

  else {

    if (valid_syntax != NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Error parsing expression:\n"
                  "%s\n"
                  "%s\n"
                  "Valid (expected) syntax:\n\n"
                  "%s"),
                infix, err_str, valid_syntax);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error parsing expression:\n"
                  "%s\n"
                  "%s"),
                infix, err_str);

  }
}

/*----------------------------------------------------------------------------
 * Check for left / right operands for unary or binary operators
 *
 * parameters:
 *   infix         <-- string parsed
 *   te            <-- tokenized expression
 *   token_id      <-- token id in tokenized expression
 *   otype         <-- type of operator
 *   has_r_operand <-- does the operator have a right operand ?
 *   os            <-- operator stack (for clean-up in case of error)
 *   postfix       <-> postfix expression (destroyed)
 *----------------------------------------------------------------------------*/

static void
_check_left_right(const char               *infix,
                  const _tokenized_t       *te,
                  int                       token_id,
                  _operator_type_t          otype,
                  bool                      has_r_operand,
                  _stack_t                 *os,
                  fvm_selector_postfix_t  **postfix)
{
  int i;

  if (otype != OT_L_PAREN && otype != OT_UNARY && otype != OT_BINARY)
    return;

  i = token_id+1;

  if (   i == te->n_tokens
      || (   te->protected[i] == false
          && te->tokens[te->token_id[i]] == ')'))
    _parse_error(_("Operator needs a right operand."),
                 NULL, infix, te, token_id, os, postfix);

  if (otype == OT_BINARY && has_r_operand == false)
    _parse_error(_("Operator needs a left operand."),
                 NULL, infix, te, token_id, os, postfix);
  else if ((otype == OT_L_PAREN || otype == OT_UNARY) && has_r_operand == true)
    _parse_error(_("Operator should not have a left operand."),
                 NULL, infix, te, token_id, os, postfix);
}

/*----------------------------------------------------------------------------
 * Find a group or attribute corresponding to a token
 *
 * parameters:
 *   n_groups     <-- number of groups
 *   n_attributes <-- number of attributes
 *   group_name   <-- array group names (sorted)
 *   attribute    <-- array of attribute numbers (sorted)
 *   token        <-- token (string) to associated with group or attribute
 *   protected    <-- is token protected (i.e. not interpretable) ?
 *   group_id     --> -1, or id of group corresponding to token
 *   attribute_id --> -1, or id of attribute corresponding to token
 *----------------------------------------------------------------------------*/

static void
_find_group_or_attribute(int          n_groups,
                         int          n_attributes,
                         const char  *group_name[],
                         const int    attribute[],
                         const char  *token,
                         bool         protected,
                         int         *group_id,
                         int         *attribute_id)
{
  int att_cmp, start_id, end_id;

  /* Initialize return values */

  *group_id = -1;
  *attribute_id = -1;

  /* Test for attributes first */

  if (protected == false) {

    int val;

    if (_is_int(token, &val) && n_attributes > 0) {

      start_id = 0;
      end_id = n_attributes - 1;

      /* use binary search */

      while (start_id <= end_id) {
        int mid_id = start_id + ((end_id - start_id) / 2);
        att_cmp = attribute[mid_id];
        if (att_cmp < val)
          start_id = mid_id + 1;
        else if (att_cmp > val)
          end_id = mid_id - 1;
        else {
          *attribute_id = mid_id;
          return;
        }

      }

    }

  }

  /* Test for groups if no attributes found */

  if (n_groups == 0)
    return;

  start_id = 0;
  end_id = n_groups - 1;

  /* use binary search */

  while (start_id <= end_id) {
    int mid_id = start_id + ((end_id - start_id) / 2);
    att_cmp = strcmp(group_name[mid_id], token);
    if (att_cmp < 0)
      start_id = mid_id + 1;
    else if (att_cmp > 0)
      end_id = mid_id - 1;
    else {
      *group_id = mid_id;
      return;
    }
  }
}

/*----------------------------------------------------------------------------
 * Determine if tokens define a group range
 *
 * parameters:
 *   n_groups     <-- number of groups
 *   group_name   <-- array of group names (sorted)
 *   token        <-> tokens (string) associated with groups
 *   group_id     --> -1, or id of group corresponding to token
 *----------------------------------------------------------------------------*/

static void
_group_range(int          n_groups,
             const char  *group_name[],
             const char  *token[2],
             int          group_id[2])
{
  int i, att_cmp, start_id, end_id, mid_id;

  /* Initialize return values */

  group_id[0] = -1;
  group_id[1] = -1;

  /* Test for groups */

  if (n_groups > 0) {

    if (strcmp(token[0], token[1]) > 0) {
      const char *_tmp = token[0];
      token[0] = token[1];
      token[1] = _tmp;
    }

    for (i = 0; i < 2; i++) {

      start_id = 0;
      end_id = n_groups - 1;
      mid_id = (end_id - start_id) / 2;

      /* use binary search */

      while (start_id < end_id) {
        att_cmp = strcmp(group_name[mid_id], token[i]);
        if (att_cmp < 0)
          start_id = mid_id + 1;
        else if (att_cmp > 0)
          end_id = mid_id - 1;
        else
          break;
        mid_id = start_id + ((end_id - start_id) / 2);
      }

      att_cmp = strcmp(group_name[mid_id], token[i]);
      group_id[i] = mid_id;
      if (i == 0 && att_cmp < 0)
        group_id[i] = mid_id + 1;
      else if (i == 1 && att_cmp > 0)
        group_id[i] = mid_id - 1;

    }

    if (   group_id[0] >= n_groups
        || group_id[1] < 0) {
      group_id[0] = -1;
      group_id[1] = -1;
    }
    else if (group_id[1] >= n_groups)
      group_id[1] -= 1;

  }

}

/*----------------------------------------------------------------------------
 * Determine if tokens define an atribute range
 *
 * parameters:
 *   n_attributes <-- number of attributes
 *   attribute    <-- array of attribute numbers (sorted)
 *   token        <-- tokens (string) associated with attributes
 *   attribute_id --> -1, or id of attribute corresponding to token
 *----------------------------------------------------------------------------*/

static void
_attribute_range(int          n_attributes,
                 const int    attribute[],
                 const char  *token[2],
                 int          attribute_id[2])
{
  int i, att_cmp, start_id, end_id, mid_id;
  int val[2];

  /* Initialize return values */

  attribute_id[0] = -1;
  attribute_id[1] = -1;

  /* Test for attributes */

  if (   n_attributes > 0
      && _is_int(token[0], val) && _is_int(token[1], val + 1)) {

    if (val[0] > val[1]) {
      int _tmp = val[0];
      val[0] = val[1];
      val[1] = _tmp;
    }

    for (i = 0; i < 2; i++) {

      start_id = 0;
      end_id = n_attributes - 1;
      mid_id = (end_id - start_id) / 2;

      /* use binary search */

      while (start_id < end_id) {
        att_cmp = attribute[mid_id];
        if (att_cmp < val[i])
          start_id = mid_id + 1;
        else if (att_cmp > val[i])
          end_id = mid_id - 1;
        else
          break;
        mid_id = start_id + ((end_id - start_id) / 2);
      }

      attribute_id[i] = mid_id;
      if (i == 0 && attribute[mid_id] < val[i])
        attribute_id[i] = mid_id + 1;
      else if (i == 1 && attribute[mid_id] > val[i])
        attribute_id[i] = mid_id - 1;

    }

    if (   attribute_id[0] >= n_attributes
        || attribute_id[1] < 0) {
      attribute_id[0] = -1;
      attribute_id[1] = -1;
    }
    else if (attribute_id[1] >= n_attributes)
      attribute_id[1] -= 1;

  }

}

/*----------------------------------------------------------------------------
 * Handle geometric functions in a tokenized expression.
 *
 * parameters:
 *   opcode       <-- operator code
 *   infix        <-- string parsed
 *   te           <-- tokenized expression
 *   token_start  <-> start id in tokenized expression
 *   token_end t  <-> past the-end id in tokenized expression
 *   os           <-> operator stack
 *   postfix      <-> pointer to postfix expression,
 *----------------------------------------------------------------------------*/

static void
_parse_geometric_args(_operator_code_t          opcode,
                      const char               *infix,
                      const _tokenized_t       *te,
                      int                       token_start,
                      int                       token_end,
                      _stack_t                 *os,
                      fvm_selector_postfix_t  **postfix)
{
  const char *tok;
  int i = token_start;

  int n_vals = 0, n_opts = 0;
  double val[13]; /* 12 values max for box(x0, dx, dy, dz),
                      1 extra slot for error checking */
  int inout = 0;   /* 0: undefined; -1: inside; 1: outside */
  double epsilon = 1.e-2, norm = 1.0;
  bool  error = false;
  bool  have_epsilon = false;

  const char *func_syntax = NULL;
  const char *normals_syntax
    = N_("  normal[<x>, <y>, <z>, <epsilon>]\n"
         "  normal[<x>, <y>, <z>, epsilon = <epsilon>]");
  const char *plane_syntax
    = N_("  For ax + by + cz + d = 0 form:\n\n"
         "    plane[<a>, <b>, <c>, <d>, <epsilon>]\n"
         "    plane[<a>, <b>, <c>, <d>, epsilon = <epsilon>]\n"
         "    plane[<a>, <b>, <c>, <d>, inside]\n"
         "    plane[<a>, <b>, <c>, <d>, outside]\n\n"
         "  For {normal, point in plane} form:\n\n"
         "    plane[<nx>, <ny>, <nz>, <x>, <y>, <z>, <epsilon>]\n"
         "    plane[<nx>, <ny>, <nz>, <x>, <y>, <z>, epsilon = <epsilon>]\n"
         "    plane[<nx>, <ny>, <nz>, <x>, <y>, <z>, inside]\n"
         "    plane[<nx>, <ny>, <nz>, <x>, <y>, <z>, outside]");
  const char *box_syntax
    = N_("  For x_min, y_min, z_min, x_max, y_max, z_max form:\n\n"
         "    box[<xmin>, <ymin>, <zmin>, <xmax>, <ymax>, <zmax>]\n\n"
         "  For xyz_0, dxyz_1, dxyz_2, dxyz_3 form:\n\n"
         "    box[<x0>, <y0>, <z0>, <dx1>, <dy1>, <dz1>\n"
         "        <dx2>, <dy2>, <dz2>, <dx3>, <dy3>, <dz3>]");
  const char *cylinder_syntax
    = N_("  cylinder[<x0>, <y0>, <z0>, <x1>, <y1>, <z1>, <radius>]");
  const char *sphere_syntax
    = N_("  sphere[<xc>, <yc>, <zc>, <radius>]");

  /* First parse for floating-point values */

  while (i < token_end && error == false) {

    tok =  te->tokens + te->token_id[i];

    if (_is_float(tok, val + n_vals)) {
      tok =  te->tokens + te->token_id[++i];
      n_vals++;
      if (n_vals > 12)
        error = true; /* Will be handled depending on opcode */
      else if (i < token_end && strcmp(te->tokens + te->token_id[i], ",")) {
        _parse_error(_("Missing or wrong argument separator."),
                     NULL, infix, te, i, os, postfix);
        error = true;
      }
      else
        i++;
    }
    else
      break;

  }

  if (i == token_end && *(te->tokens + te->token_id[i-1]) == ',')
    _parse_error(_("Missing argument after separator."),
                 NULL, infix, te, i, os, postfix);

  /* Initialize error reporting function syntax,
     check number of floating-point arguments,
     and normalize normals. */

  switch(opcode) {

  case OC_NORMAL:
    func_syntax = normals_syntax;
    if (n_vals == 4) {
      have_epsilon = true;
      epsilon = val[n_vals-1];
      n_vals--;
    }
    if (n_vals != 3)
      error = true;
    norm = sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]);
    val[0] /= norm; val[1] /= norm; val[2] /= norm;
    break;

  case OC_PLANE:
    func_syntax = plane_syntax;
    if (n_vals == 5 || n_vals == 7) {
      have_epsilon = true;
      epsilon = val[n_vals-1];
      n_vals--;
    }
    if (n_vals != 4 && n_vals != 6)
      error = true;
    norm = sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]);
    val[0] /= norm; val[1] /= norm; val[2] /= norm;
    /* Normalize d in ax + by + cz + d form, or convert to this form */
    if (n_vals == 4)
      val[3] /= norm;
    else { /* if (n_vals == 6) */
      val[3] = - (val[0]*val[3] + val[1]*val[4] + val[2]*val[5]);
      n_vals = 4;
    }
    break;

  case OC_BOX:
    func_syntax = box_syntax;
    if (n_vals != 6 && n_vals != 12)
      error = true;
    break;

  case OC_CYLINDER:
    func_syntax = cylinder_syntax;
    if (n_vals != 7)
      error = true;
    break;

  case OC_SPHERE:
    func_syntax = sphere_syntax;
    if (n_vals != 4)
      error = true;
    break;

  default:
    _parse_error(_("This geometric function is not implemented."),
                 NULL, infix, te, token_start - 2, os, postfix);
    break;
  }

  if (error == true)
    _parse_error(_("Wrong number of floating-point arguments."),
                 _(func_syntax),
                 infix, te, token_start, os, postfix);

  /* Check for one additional (key, value) arguments */

  if (i < token_end && error == false) {

    if ((opcode == OC_NORMAL || opcode == OC_PLANE) && have_epsilon == false) {

      if (strcmp(te->tokens + te->token_id[i], "epsilon") == 0) {
        if (strcmp(te->tokens + te->token_id[i+1], "="))
          error = true;
        tok =  te->tokens + te->token_id[i+2];
        if (_is_float(tok, &epsilon)) {
          have_epsilon = true;
          i += 3;
        }
        else {
          _parse_error(_("Expected syntax:\n"
                         "  epsilon = <value>\n"),
                       NULL, infix, te, i, os, postfix);
          error = true;
        }
      }

    }

    if (opcode == OC_PLANE && have_epsilon == false) {
      if (strcmp(te->tokens + te->token_id[i], "inside") == 0) {
        inout = -1;
        i++;
      }
      else if (strcmp(te->tokens + te->token_id[i], "outside") == 0) {
        inout = 1;
        i++;
      }
    }

  }

  if (i < token_end) {
    _parse_error(_("Unexpected argument(s)."),
                 _(func_syntax),
                 infix, te, i, os, postfix);
  }

  /* Update postfix */

  if (opcode == OC_NORMAL) {
    (*postfix)->normals_dependency = true;
    /* sign*square of cosine compared to square of (1-epsilon) */
    val[n_vals++] = (1 - 2*epsilon + epsilon*epsilon);
    assert(n_vals == 4);
  }
  else {
    (*postfix)->coords_dependency = true;
    if (opcode == OC_PLANE) {
      if (inout == 0)
        val[n_vals++] = epsilon;
      else
        n_opts++;
      assert(n_vals + n_opts == 5);
    }
  }

  _postfix_add_opcode(*postfix, opcode);
  _postfix_add_int(*postfix, n_vals + n_opts, PF_INT);
  for (i = 0; i < n_vals; i++)
    _postfix_add_float(*postfix, val[i]);
  if (n_opts == 1) {
    assert (opcode == OC_PLANE && inout != 0);
    _postfix_add_int(*postfix, inout, PF_INT);
  }

}

/*----------------------------------------------------------------------------
 * Check for and handle coordinate conditions in a tokenized expression.
 *
 * parameters:
 *   this_parser   <-- parser object
 *   infix         <-- string parsed
 *   te            <-- tokenized expression
 *   n_groups      <-- number of groups
 *   n_attributes  <-- number of attributes
 *   group_name    <-- array group names (sorted)
 *   attribute     <-- array of attribute numbers (sorted)
 *   token_id      <-> current token id in tokenized expression
 *   has_r_operand <-> indicates if the current token has a right operand
 *   os            <-> operator stack
 *   postfix       <-> pointer to postfix expression,
 *----------------------------------------------------------------------------*/

static void
_parse_for_function(const _parser_t          *this_parser,
                    const char               *infix,
                    const _tokenized_t       *te,
                    int                       n_groups,
                    int                       n_attributes,
                    const char               *group_name[],
                    const int                 attribute[],
                    int                      *token_id,
                    bool                     *has_r_operand,
                    _stack_t                 *os,
                    fvm_selector_postfix_t  **postfix)
{
  const char *tok;
  int i = *token_id + 1, j = 0, k = 0;
  const _operator_t *op = NULL;

  if (te->n_tokens <= i)
    return;

  tok = te->tokens + te->token_id[i];

  /* Pre-check syntax */

  if (te->protected[i] == true || strlen(tok) != 1)
    return;

  if (tok[0] != '[')
    return;

  /* If we arrived here, we have a function start */

  k = 1;
  for (j = i + 1; j < te->n_tokens; j++) {
    tok = te->tokens + te->token_id[j];
    if (te->protected[j] == false && strlen(tok) == 1) {
      if (tok[0] == '[')
        k++;
      else if (tok[0] == ']')
        k--;
    }
    if (k == 0)
      break;
  }
  if (j == te->n_tokens) {
    _parse_error(_("Missing closing ]."),
                 NULL, infix, te, i, os, postfix);
    return;
  }

  /* We have a function-type syntax, so find the corresponding function */

  tok = te->tokens + te->token_id[*token_id];

  if (te->protected[i] == false) {
    for (k = 0; k < this_parser->n_keywords; k++) {
      if (strcmp(tok, this_parser->keyword[k]) == 0) {
        op = this_parser->operators + this_parser->keyword_op_id[k];
        break;
      }
    }
  }

  if (op == NULL) {
    _parse_error(_("Function arguments used with an unknown operator."),
                 NULL, infix, te, *token_id, os, postfix);
    return;
  }
  else if (op->type != OT_FUNCTION) {
    _parse_error(_("Operator does not accept function arguments."),
                 NULL, infix, te, *token_id, os, postfix);
    return;
  }

  *token_id = j + 1;
  *has_r_operand = true;

  /* Handle general group selection operators */

  if (op->code == OC_ALL || op->code == OC_NO_GROUP) {
    if (j != i+1)
      _parse_error(_("Function requires 0 arguments."),
                   NULL, infix, te, *token_id, os, postfix);
    else
      _postfix_add_opcode(*postfix, op->code);
  }

  /* Handle range operator */

  else if (op->code == OC_RANGE) {

    const char *t[3] = {NULL, NULL, NULL};
    bool  force_group = false, force_attrib = false, error = false;

    i++;
    k = 0;
    while (i < j && k < 3) {
      t[k++] = te->tokens + te->token_id[i];
      i++;
      if (i < j) {
        if (strcmp(te->tokens + te->token_id[i], ","))
          error = true;
        else {
          i++;
          if (i == j)
            error = true;
        }
      }
    }
    if (i < j)
      error = true;

    if (k == 3) {
      if (strcmp(t[2], "group") == 0)
        force_group = true;
      else if (strcmp(t[2], "attribute") == 0)
        force_attrib = true;
      else
        error = true;
      k = 2;
    }

    if (k != 2 || error)
      _parse_error(_("range[] argument error"),
                   _("  range[<first>, <last>]\n"
                     "  range[<first>, <last>, group]\n"
                     "  range[<first>, <last>, attribute]"),
                   infix, te, i - 1, os, postfix);


    else {
      int ga_id[2] = {-1, -1};
      _postfix_type_t pf_type = PF_GROUP_ID;
      if (force_group == false)
        _attribute_range(n_attributes, attribute, t, ga_id);
      if (force_attrib == true || ga_id[0] > -1)
        pf_type = PF_ATTRIBUTE_ID;
      else
        _group_range(n_groups, group_name, t, ga_id);
      if (ga_id[0] > -1) {
        _postfix_add_opcode(*postfix, op->code);
        _postfix_add_int(*postfix, ga_id[0], pf_type);
        _postfix_add_int(*postfix, ga_id[1], pf_type);
      }
      else {
        _postfix_add_int(*postfix, -1, pf_type);
        _postfix_add_missing(*postfix, t[0]);
        _postfix_add_missing(*postfix, t[1]);
      }
    }

  }

  /* Handle geometric operators */

  else
    _parse_geometric_args(op->code, infix, te, i+1, j, os, postfix);

}

/*----------------------------------------------------------------------------
 * Check for and handle coordinate conditions in a tokenized expression.
 *
 * parameters:
 *   infix         <-- string parsed
 *   te            <-- tokenized expression
 *   token id      <-> current token id in tokenized expression
 *   has_r_operand <-> indicates if the current token has a right operand
 *   os            <-> operator stack
 *   postfix       <-> pointer to postfix expression,
 *----------------------------------------------------------------------------*/

static void
_parse_for_coord_conditions(const char               *infix,
                            const _tokenized_t       *te,
                            int                      *token_id,
                            bool                     *has_r_operand,
                            _stack_t                 *os,
                            fvm_selector_postfix_t  **postfix)
{
  const char *t1;
  size_t      t1_len;
  double      val;
  bool  has_coord_cond = false;
  int coord_id = -1;
  int i = *token_id + 1, j = 0;

  if (te->n_tokens <= i)
    return;

  t1 = te->tokens + te->token_id[i];
  t1_len = strlen(t1);

  /* Pre-check syntax */

  if (te->protected[i] == true || t1_len == 0 || t1_len > 2)
    return;

  if ((t1[0] != '<' && t1[0] != '>') || (t1_len == 2 && t1[1] != '='))
    return;

  /* If we arrived here, we have a >, <, >=, or <= operator */

  if (te->n_tokens == i+1) {
    _parse_error(_("Operator needs a right operand."),
                 NULL, infix, te, i, os, postfix);
    return;
  }

  /* Try for coord_id, operator, value or value, operator, coord_id */

  for (j = 0; j < 2; j++) {

    const char *t2 = te->tokens + te->token_id[i - 1 + j*2];

    if (strlen(t2) == 1) {
      if (*t2 == 'x' || *t2 == 'X')
        coord_id = 0;
      else if (*t2 == 'y' || *t2 == 'Y')
        coord_id = 1;
      else if (*t2 == 'z' || *t2 == 'Z')
        coord_id = 2;
    }
    if (coord_id != -1)
      break;
  }

  if (j == 0)
    has_coord_cond = _is_float(te->tokens + te->token_id[i + 1], &val);
  else if (j == 1)
    has_coord_cond = _is_float(te->tokens + te->token_id[i - 1], &val);

  /* If we have a valid syntax, add it to postfix */

  if (has_coord_cond) {

    _operator_code_t oc;

    /* Permute operator if necessery to always have a
       {coord_id, operator, value} expression */

    if (j == 0) {
      if (t1[0] == '<' )
        oc = OC_LT;
      else
        oc = OC_GT;
    }
    else {
      assert(j == 1);
      if (t1[0] == '<' )
        oc = OC_GT;
      else
        oc = OC_LT;
    }
    if (t1_len == 2) {
      if (oc == OC_LT)
        oc = OC_LE;
      else
        oc = OC_GE;
    }

    _postfix_add_opcode(*postfix, oc);
    _postfix_add_int(*postfix, coord_id, PF_INT);
    _postfix_add_float(*postfix, val);

    i += 2;
    *token_id = i;
    *has_r_operand = true;
    (*postfix)->coords_dependency = true;

  }
  else {
    _parse_error(_("Operator needs a floating point operand\n"
                   "on one side, x, y, or z on the other"),
                 NULL, infix, te, i, os, postfix);
    return;
  }

  /* If we have a valid syntax with the coord_id on the
     left, we may have a "1 < x <= 2" type syntax */

  if (has_coord_cond && j == 1 && te->n_tokens > i) {

    const char *t3 = te->tokens + te->token_id[i];
    size_t      t3_len = strlen(t3);

    if (te->protected[i] == true || t3_len == 0 || t3_len > 2)
      return;

    if (   (t3[0] != '<' && t3[0] != '>')
        || (t3_len == 2 && t3[1] != '='))
      return;

    if (t3[0] != t1[0]) {
      _parse_error(_("Inconsistant interval specification."),
                   NULL, infix, te, i, os, postfix);
      return;
    }
    else if (te->n_tokens == i+1) {
      _parse_error(_("Operator needs a right operand."),
                   NULL, infix, te, i, os, postfix);
      return;
    }

    has_coord_cond = _is_float(te->tokens + te->token_id[i + 1], &val);

    /* If we have a valid syntax, add it to postfix */

    if (has_coord_cond) {

      _operator_code_t oc;

      /* No permutation necessary here */

      if (t3[0] == '<' )
        oc = OC_LT;
      else
        oc = OC_GT;
      if (t1_len == 2) {
        if (oc == OC_LT)
          oc = OC_LE;
        else
        oc = OC_GE;
      }

      _postfix_add_opcode(*postfix, oc);
      _postfix_add_int(*postfix, coord_id, PF_INT);
      _postfix_add_float(*postfix, val);

      /* Add implicit and here */

      _postfix_add_opcode(*postfix, OC_AND);

      i += 2;
      *token_id = i;

    }
    else {
      _parse_error(_("Operator needs a floating point operand"),
                   NULL, infix, te, i, os, postfix);
      return;
    }

  }
}

/*----------------------------------------------------------------------------
 * Parse a tokenized expression string.
 *
 * parameters:
 *   this_parser  <-- parser object
 *   infix        <-- string parsed
 *   te           <-- tokenized expression
 *   n_groups     <-- number of groups
 *   n_attributes <-- number of attributes
 *   group_name   <-- array group names (sorted)
 *   attribute    <-- array of attribute numbers (sorted)
 *
 * returns:
 *   pointer to created postfix structure
 *----------------------------------------------------------------------------*/

static fvm_selector_postfix_t *
_parse_tokenized(const _parser_t     *this_parser,
                 const char          *infix,
                 const _tokenized_t  *te,
                 int                  n_groups,
                 int                  n_attributes,
                 const char          *group_name[],
                 const int            attribute[])
{
  int i, j;
  _stack_t os;
  bool  has_r_operand = false;
  fvm_selector_postfix_t *pf = NULL;

  /* Initialization */

  _stack_initialize(&os);

  pf = _postfix_create(infix);

  /* Loop on tokens */

  i = 0;

  while (i < te->n_tokens) {

    const char *tok;
    const _operator_t *op_1 = NULL;

    /* Look ahead to handle functions and "<", "<=", ">", ">=" syntax */

    _parse_for_function(this_parser,
                        infix,
                        te,
                        n_groups,
                        n_attributes,
                        group_name,
                        attribute,
                        &i,
                        &has_r_operand,
                        &os,
                        &pf);
    _parse_for_coord_conditions(infix, te, &i, &has_r_operand, &os, &pf);

    if (i == te->n_tokens)
      break;

    /* Now handle token */

    tok = te->tokens + te->token_id[i];

    if (te->protected[i] == false) {
      for (j = 0; j < this_parser->n_keywords; j++) {
        if (strcmp(tok, this_parser->keyword[j]) == 0) {
          op_1 = this_parser->operators + this_parser->keyword_op_id[j];
          break;
        }
      }
    }

    /* Basic check for left / right operands to operators */

    if (op_1 != NULL)
      _check_left_right(infix, te, i, op_1->type, has_r_operand, &os, &pf);

    /* Now add to postfix or stack */

    if (op_1 == NULL) {

      int  group_id, attribute_id;

      /* If two operands follow each other, we have a parse error */

      if (has_r_operand == true)
        _parse_error(_("Expected operator instead of operand."),
                     NULL, infix, te, i, &os, &pf);

      /* Now add entry to postfix */

      _find_group_or_attribute(n_groups,
                               n_attributes,
                               group_name,
                               attribute,
                               tok,
                               te->protected[i],
                               &group_id,
                               &attribute_id);

      if (attribute_id > -1)
        _postfix_add_int(pf, attribute_id, PF_ATTRIBUTE_ID);
      else
        _postfix_add_int(pf, group_id, PF_GROUP_ID);

      if (attribute_id == -1 && group_id == -1)
        _postfix_add_missing(pf, tok);

    }
    else {

      switch(op_1->type) {

      case OT_L_PAREN:
        _stack_push(&os, op_1, i);
        break;

      case OT_R_PAREN:
        {
          const _operator_t *op_2 = NULL;
          bool  matched = false;
          while (_stack_size(&os) > 0) {
            _stack_entry_t e = _stack_pop(&os);
            op_2 = e.op;
            if (op_2->type == OT_L_PAREN) {
              matched = true;
              break;
            }
            else {
              _postfix_add_opcode(pf, op_2->code);
            }
          }
          if (matched == false)
            _parse_error(_("Parenthesis mismatch"),
                         NULL, infix, te, i, &os, &pf);
        }
        break;

        /* At this stage, a function or "<", "<=", ">", ">=" syntax should
           have been handled, so we have a syntax error here*/

      case OT_COORD_CONDITION:
      case OT_FUNCTION:
        _parse_error(_("Syntax error, probably due to misplaced operands."),
                     NULL, infix, te, i, &os, &pf);
        break;

      default:
        {
          while (_stack_size(&os) > 0) {
            const _operator_t *op_2 = _stack_top(&os);
            if (op_1->priority < op_2->priority) {
              _stack_entry_t e = _stack_pop(&os);
              _postfix_add_opcode(pf, (e.op)->code);
            }
            else
              break;
          }
          _stack_push(&os, op_1, i);
        }
      } /* End of switch on operator type */

    }

    i++;

    has_r_operand = true;
    if (op_1 != NULL) {
      if (op_1->type != OT_R_PAREN && op_1->type != OT_FUNCTION)
        has_r_operand = false;
    }

  } /* End of loop on tokens */

  while (_stack_size(&os) > 0) {

    _stack_entry_t e = _stack_pop(&os);
    const _operator_t *op = e.op;

    if (op->type == OT_L_PAREN)
      _parse_error(_("Parenthesis mismatch"),
                   NULL, infix, te, e.token_id, &os, &pf);
    else {
      _postfix_add_opcode(pf, op->code);
    }

  }

  _postfix_adjust(pf);

  return pf;
}

/*----------------------------------------------------------------------------
 * Evaluate normal[] function conditions in a postfix expression
 *
 * parameters:
 *   pf        <-- pointer to postfix structure
 *   normal    <-- normal associated with evaluation
 *   i         <-> current position in expression being evaluated
 *
 * returns:
 *   true or false depending on evaluation
 *----------------------------------------------------------------------------*/

static inline bool
_eval_normal(const fvm_selector_postfix_t  *pf,
             const double                   normal[],
             size_t                        *i)
{
  double dotp;
  double val[4];
  int j;
  const int n_vals = 4; /* 3 values for coordinates, 1 for tolerance */
  bool  retval = false;

  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  assert(*((int *)(pf->elements + *i + _postfix_type_size)) == n_vals);

  *i += _postfix_type_size + _postfix_int_size + _postfix_type_size;

  for (j = 0; j < n_vals; j++) {
    assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
           == PF_FLOAT);
    val[j] = *((double *)(pf->elements + *i));
    *i += _postfix_float_size + _postfix_type_size;
  }
  *i -= _postfix_type_size;

  dotp = normal[0]*val[0] + normal[1]*val[1] + normal[2]*val[2];

  if (dotp > 0) {
    double cos2 = dotp*dotp / (  normal[0]*normal[0]
                               + normal[1]*normal[1]
                               + normal[2]*normal[2]);
    if (cos2 > val[3])
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Evaluate plane[] function conditions in a postfix expression
 *
 * parameters:
 *   pf        <-- pointer to postfix structure
 *   coords    <-- coordinates associated with evaluation
 *   i         <-> current position in expression being evaluated
 *
 * returns:
 *   true or false depending on evaluation
 *----------------------------------------------------------------------------*/

static inline bool
_eval_plane(const fvm_selector_postfix_t  *pf,
            const double                   coords[],
            size_t                        *i)
{
  double pfunc;
  double val[4];
  int j;
  _postfix_type_t pf_type;
  bool  retval = false;

  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  assert(*((int *)(pf->elements + *i + _postfix_type_size)) == 5);

  /* 4 first arguments are always  floating-point
     (plane equation coefficients) */

  *i += _postfix_type_size + _postfix_int_size + _postfix_type_size;
  for (j = 0; j < 4; j++) {
    assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
           == PF_FLOAT);
    val[j] = *((double *)(pf->elements + *i));
    *i += _postfix_float_size + _postfix_type_size;
  }
  *i -= _postfix_type_size;

  pfunc = val[0]*coords[0] + val[1]*coords[1] + val[2]*coords[2] + val[3];

  /* Last argument may be floating-point or interger */

  pf_type = *((_postfix_type_t *)(pf->elements + *i));
  *i += _postfix_type_size;

  if (pf_type == PF_INT) {
    int inout = *((int *)(pf->elements + *i));
    *i += _postfix_int_size;
    assert(inout == -1 || inout == 1);
    if (   (inout == -1 && pfunc <= 0)
        || (inout == 1 && pfunc >= 0))
      retval = true;
  }
  else {
    double epsilon = *((double *)(pf->elements + *i));
    assert(pf_type == PF_FLOAT);
    *i += _postfix_float_size;
    if (CS_ABS(pfunc) < epsilon)
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Evaluate box[] function conditions in a postfix expression
 *
 * parameters:
 *   pf        <-- pointer to postfix structure
 *   coords    <-- coordinates associated with evaluation
 *   i         <-> current position in expression being evaluated
 *
 * returns:
 *   true or false depending on evaluation
 *----------------------------------------------------------------------------*/

static inline bool
_eval_box(const fvm_selector_postfix_t  *pf,
          const double                   coords[],
          size_t                        *i)
{
  double val[12];
  int j;
  int n_vals; /* number of box coefficients */
  bool  retval = false;

  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  n_vals = *((int *)(pf->elements + *i + _postfix_type_size));
  assert(n_vals == 6 || n_vals == 12);

  /* all arguments are floating-point */

  *i += _postfix_type_size + _postfix_int_size + _postfix_type_size;
  for (j = 0; j < n_vals; j++) {
    assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
           == PF_FLOAT);
    val[j] = *((double *)(pf->elements + *i));
    *i += _postfix_float_size + _postfix_type_size;
  }
  *i -= _postfix_type_size;

  /* Geometric test */

  if (n_vals == 6) {
    if (   val[0] <= coords[0] && val[1] <= coords[1] && val[2] <= coords[2]
        && val[3] >= coords[0] && val[4] >= coords[1] && val[5] >= coords[2])
      retval = true;
  }
  else {
    double _coords[3], l12, l22, l32, dp1, dp2, dp3;
    _coords[0] = coords[0] - val[0];
    _coords[1] = coords[1] - val[1];
    _coords[2] = coords[2] - val[2];
    l12 = val[3]*val[3] + val[4]*val[4] + val[5]*val[5];
    l22 = val[6]*val[6] + val[7]*val[7] + val[8]*val[8];
    l32 = val[9]*val[9] + val[10]*val[10] + val[11]*val[11];
    dp1 = _coords[0]*val[3] + _coords[1]*val[4] + _coords[2]*val[5];
    dp2 = _coords[0]*val[6] + _coords[1]*val[7] + _coords[2]*val[8];
    dp3 = _coords[0]*val[9] + _coords[1]*val[10] + _coords[2]*val[11];
    if (   dp1 >= 0   && dp2 >= 0   && dp3 >= 0
        && dp1 <= l12 && dp2 <= l22 && dp3 <= l32)
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Evaluate cylinder[] function conditions in a postfix expression
 *
 * parameters:
 *   pf        <-- pointer to postfix structure
 *   coords    <-- coordinates associated with evaluation
 *   i         <-> current position in expression being evaluated
 *
 * returns:
 *   true or false depending on evaluation
 *----------------------------------------------------------------------------*/

static inline bool
_eval_cylinder(const fvm_selector_postfix_t  *pf,
               const double                   coords[],
               size_t                        *i)
{
  double val[7];
  int j;
  const int n_vals = 7; /* number of cylinder coefficients */
  bool  retval = false;

  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  assert(*((int *)(pf->elements + *i + _postfix_type_size)) == n_vals);

  /* all arguments are floating-point */

  *i += _postfix_type_size + _postfix_int_size + _postfix_type_size;
  for (j = 0; j < n_vals; j++) {
    assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
           == PF_FLOAT);
    val[j] = *((double *)(pf->elements + *i));
    *i += _postfix_float_size + _postfix_type_size;
  }
  *i -= _postfix_type_size;

  /* Geometric test */

  {
    double _coords[3], _axis[3], dotp, len2;
    _coords[0] = coords[0] - val[0];
    _coords[1] = coords[1] - val[1];
    _coords[2] = coords[2] - val[2];
    _axis[0] = val[3] - val[0];
    _axis[1] = val[4] - val[1];
    _axis[2] = val[5] - val[2];
    dotp = _coords[0]*_axis[0] + _coords[1]*_axis[1] + _coords[2]*_axis[2];
    len2 = _axis[0]*_axis[0] + _axis[1]*_axis[1] + _axis[2]*_axis[2];
    if (dotp >= 0 && dotp <= len2) {
      double _proj[3], r2;
      double mult = dotp / len2;
      _proj[0] = _coords[0] - mult*_axis[0];
      _proj[1] = _coords[1] - mult*_axis[1];
      _proj[2] = _coords[2] - mult*_axis[2];
      r2 = _proj[0]*_proj[0] + _proj[1]*_proj[1] + _proj[2]*_proj[2];
      if (r2 <= val[6]*val[6])
        retval = true;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Evaluate sphere[] function conditions in a postfix expression
 *
 * parameters:
 *   pf        <-- pointer to postfix structure
 *   coords    <-- coordinates associated with evaluation
 *   i         <-> current position in expression being evaluated
 *
 * returns:
 *   true or false depending on evaluation
 *----------------------------------------------------------------------------*/

static inline bool
_eval_sphere(const fvm_selector_postfix_t  *pf,
             const double                   coords[],
             size_t                        *i)
{
  double val[4];
  int j;
  const int n_vals = 4; /* number of sphere coefficients */
  bool  retval = false;

  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  assert(*((int *)(pf->elements + *i + _postfix_type_size)) == n_vals);

  /* all arguments are floating-point */

  *i += _postfix_type_size + _postfix_int_size + _postfix_type_size;
  for (j = 0; j < n_vals; j++) {
    assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
           == PF_FLOAT);
    val[j] = *((double *)(pf->elements + *i));
    *i += _postfix_float_size + _postfix_type_size;
  }
  *i -= _postfix_type_size;

  /* Geometric test */

  {
    double _coords[3], len2;
    _coords[0] = coords[0] - val[0];
    _coords[1] = coords[1] - val[1];
    _coords[2] = coords[2] - val[2];
    len2
      = _coords[0]*_coords[0] + _coords[1]*_coords[1] + _coords[2]*_coords[2];
    if (len2 <= val[3]*val[3])
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Evaluate coordinate conditions in a postfix expression
 *
 * parameters:
 *   pf        <-- pointer to postfix structure
 *   coords    <-- coordinates associated with evaluation
 *   i         <-> current position in expression being evaluated
 *
 * returns:
 *   true or false depending on evaluation
 *----------------------------------------------------------------------------*/

static inline bool
_eval_coord_gt(const fvm_selector_postfix_t  *pf,
               const double                   coords[],
               size_t                        *i)
{
  double cmp_val;
  int coord_id;
  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  *i += _postfix_type_size;
  coord_id = *((int *)(pf->elements + *i));
  *i += _postfix_int_size + _postfix_type_size;
  assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
         == PF_FLOAT);
  cmp_val = *((double *)(pf->elements + *i));
  *i += _postfix_float_size;
  return (coords[coord_id] > cmp_val ? true : false);
}

static inline bool
_eval_coord_lt(const fvm_selector_postfix_t  *pf,
               const double                   coords[],
               size_t                        *i)
{
  double cmp_val;
  int coord_id;
  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  *i += _postfix_type_size;
  coord_id = *((int *)(pf->elements + *i));
  *i += _postfix_int_size + _postfix_type_size;
  assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
         == PF_FLOAT);
  cmp_val = *((double *)(pf->elements + *i));
  *i += _postfix_float_size;
  return (coords[coord_id] < cmp_val ? true : false);
}

static inline bool
_eval_coord_ge(const fvm_selector_postfix_t  *pf,
               const double                   coords[],
               size_t                        *i)
{
  double cmp_val;
  int coord_id;
  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  *i += _postfix_type_size;
  coord_id = *((int *)(pf->elements + *i));
  *i += _postfix_int_size + _postfix_type_size;
  assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
         == PF_FLOAT);
  cmp_val = *((double *)(pf->elements + *i));
  *i += _postfix_float_size;
  return (coords[coord_id] >= cmp_val ? true : false);
}

static inline bool
_eval_coord_le(const fvm_selector_postfix_t  *pf,
               const double                   coords[],
               size_t                        *i)
{
  double cmp_val;
  int coord_id;
  assert(*((_postfix_type_t *)(pf->elements + *i)) == PF_INT);
  *i += _postfix_type_size;
  coord_id = *((int *)(pf->elements + *i));
  *i += _postfix_int_size + _postfix_type_size;
  assert(*((_postfix_type_t *)(pf->elements + *i - _postfix_type_size))
         == PF_FLOAT);
  cmp_val = *((double *)(pf->elements + *i));
  *i += _postfix_float_size;
  return (coords[coord_id] <= cmp_val ? true : false);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a postfix expression from an infix expression
 *
 * parameters:
 *   infix        <-- infix expression
 *   n_groups     <-- number of groups
 *   n_attributes <-- number of attributes
 *   group_name   <-- array group names (sorted)
 *   attribute    <-- array of attribute numbers (sorted)
 *
 * returns:
 *   pointer to created postfix structure
 *----------------------------------------------------------------------------*/

fvm_selector_postfix_t *
fvm_selector_postfix_create(const char  *infix,
                            int          n_groups,
                            int          n_attributes,
                            const char  *group_name[],
                            const int    attribute[])
{
  _tokenized_t te = _tokenize(infix);
  fvm_selector_postfix_t * pf = NULL;

  if (_n_parser_references == 0)
    _parser = _parser_create();
  _n_parser_references++;

  pf = _parse_tokenized(_parser,
                        infix,
                        &te,
                        n_groups,
                        n_attributes,
                        group_name,
                        attribute);

  _empty_tokenized(&te);

  /* Return postfix expression pointer */

  return pf;
}

/*----------------------------------------------------------------------------
 * Destroy a postfix expression
 *
 * parameters:
 *   pf <-> pointer to postfix structure pointer
 *----------------------------------------------------------------------------*/

void
fvm_selector_postfix_destroy(fvm_selector_postfix_t  **postfix)
{
  assert(postfix != NULL);
  assert(*postfix != NULL);

  _n_parser_references--;
  if (_n_parser_references == 0)
    _parser_destroy(&_parser);

  _postfix_destroy(postfix);
}

/*----------------------------------------------------------------------------
 * Return a pointer to the infix string associated with a postfix expression
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   pointer to original infix string
 *----------------------------------------------------------------------------*/

const char *
fvm_selector_postfix_get_infix(const fvm_selector_postfix_t  *pf)
{
  assert(pf != NULL);

  return pf->infix;
}

/*----------------------------------------------------------------------------
 * Indicate if a postfix expression depends on coordinates
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   true if expression depends on coordinates, false otherwise
 *----------------------------------------------------------------------------*/

bool
fvm_selector_postfix_coords_dep(const fvm_selector_postfix_t  *pf)
{
  assert(pf != NULL);

  return pf->coords_dependency;
}

/*----------------------------------------------------------------------------
 * Indicate if a postfix expression depends on normals
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   true if expression depends on normals, false otherwise
 *----------------------------------------------------------------------------*/

bool
fvm_selector_postfix_normals_dep(const fvm_selector_postfix_t  *pf)
{
  assert(pf != NULL);

  return pf->normals_dependency;
}

/*----------------------------------------------------------------------------
 * Return the number of operands associated with a postfix expression
 * missing in the associated group class set
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   number of missing operands
 *----------------------------------------------------------------------------*/

int
fvm_selector_postfix_n_missing(const fvm_selector_postfix_t  *pf)
{
  assert(pf != NULL);

  return pf->n_missing_operands;
}

/*----------------------------------------------------------------------------
 * Return a pointer to the name of an of operand associated with a postfix
 * expression but missing in the associated group class set
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *   id <-- id of missing operand (0 to fvm_selector_postfix_n_missing())
 *
 * returns:
 *   pointer to name of missing operand
 *----------------------------------------------------------------------------*/

const char *
fvm_selector_postfix_get_missing(const fvm_selector_postfix_t  *pf,
                                 int                            id)
{
  const char *retval = NULL;
  assert(pf != NULL);

  if (id > -1 && id < pf->n_missing_operands)
    retval = pf->missing_operand[id];

  return retval;
}

/*----------------------------------------------------------------------------
 * Evaluate a postfix expression
 *
 * parameters:
 *   pf           <-- pointer to postfix structure
 *   n_groups     <-- number of groups associated with group class
 *   n_attributes <-- number of attributes associated with group class
 *   group_id     <-- array group ids associated with group class
 *   attribute_id <-- array of attribute ids associated with group class
 *   coords       <-- coordinates associated with evaluation, or NULL
 *   normal       <-- normal associated with evaluation, or NULL
 *
 * returns:
 *   true or false base on expression evaluation
 *----------------------------------------------------------------------------*/

bool
fvm_selector_postfix_eval(const fvm_selector_postfix_t  *pf,
                          int                            n_groups,
                          int                            n_attributes,
                          const int                      group_id[],
                          const int                      attribute_id[],
                          const double                   coords[],
                          const double                   normal[])
{
  bool  retval;
  bool  _eval_stack[BASE_STACK_SIZE];
  bool  *eval_stack = _eval_stack;
  size_t i = 0, eval_size = 0, eval_max_size = BASE_STACK_SIZE;

  /* Evaluate postfix_string */

  i = 0;

  while (i < pf->size) {

    _postfix_type_t type = *((_postfix_type_t *)(pf->elements + i));

    i += _postfix_type_size;

    switch(type) {

    case PF_GROUP_ID:
      {
        int j;
        int val = *((int *)(pf->elements + i));
        i += _postfix_int_size;
        eval_stack[eval_size] = false;
        for (j = 0; j < n_groups; j++) {
          if (val == group_id[j]) {
            eval_stack[eval_size] = true;
            break;
          }
        }
        eval_size++;
      }
      break;
    case PF_ATTRIBUTE_ID:
      {
        int j;
        int val = *((int *)(pf->elements + i));
        i += _postfix_int_size;
        eval_stack[eval_size] = false;
        for (j = 0; j < n_attributes; j++) {
          if (val == attribute_id[j]) {
            eval_stack[eval_size] = true;
            break;
          }
        }
        eval_size++;
      }
      break;
    case PF_OPCODE:
      {
        size_t min_eval_size;
        _operator_code_t oc = *((_operator_code_t *)(pf->elements + i));
        i += _postfix_opcode_size;

        if (oc == OC_NOT)
          min_eval_size = 1;
        else if (oc >= OC_AND && oc <= OC_XOR)
          min_eval_size = 2;
        else
          min_eval_size = 0;

        if (eval_size < min_eval_size) {
          fvm_selector_postfix_dump(pf, 0, 0, NULL, NULL);
          bft_error(__FILE__, __LINE__, 0,
                    _("Postfix evaluation error."));
        }

        switch(oc) {

        case OC_NOT:
          if (eval_stack[eval_size-1] == false)
            eval_stack[eval_size-1] = true;
          else
            eval_stack[eval_size-1] = false;
          break;
        case OC_AND:
          eval_stack[eval_size-2] = (   eval_stack[eval_size-2]
                                     && eval_stack[eval_size-1]);
          eval_size--;
          break;
        case OC_OR:
          eval_stack[eval_size-2] = (   eval_stack[eval_size-2]
                                     || eval_stack[eval_size-1]);
          eval_size--;
          break;
        case OC_XOR:
          if (eval_stack[eval_size-2] != eval_stack[eval_size-1])
            eval_stack[eval_size-2] = true;
          else
            eval_stack[eval_size-2] = false;
          eval_size--;
          break;

        case OC_ALL:
          eval_stack[eval_size] = true;
          eval_size++;
          break;
        case OC_NO_GROUP:
          if (n_groups == 0 && n_attributes == 0)
            eval_stack[eval_size] = true;
          else
            eval_stack[eval_size] = false;
          eval_size++;
          break;

        case OC_RANGE:
          {
            _postfix_type_t type1, type2;
            int val1, val2;

            type1 = *((_postfix_type_t *)(pf->elements + i));
            i += _postfix_type_size;
            val1 = *((int *)(pf->elements + i));
            i += _postfix_int_size;
            type2 = *((_postfix_type_t *)(pf->elements + i));
            i += _postfix_type_size;
            val2 = *((int *)(pf->elements + i));
            i += _postfix_int_size;

            if (type1 == PF_GROUP_ID && type1 == type2) {
              int j;
              eval_stack[eval_size] = false;
              for (j = 0; j < n_groups; j++) {
                if (group_id[j] >= val1 && group_id[j] <= val2) {
                  eval_stack[eval_size] = true;
                  break;
                }
              }
            }
            else if (type1 == PF_ATTRIBUTE_ID && type1 == type2) {
              int j;
              eval_stack[eval_size] = false;
              for (j = 0; j < n_attributes; j++) {
                if (attribute_id[j] >= val1 && attribute_id[j] <= val2) {
                  eval_stack[eval_size] = true;
                  break;
                }
              }
            }
            else {
              fvm_selector_postfix_dump(pf, 0, 0, NULL, NULL);
              bft_error(__FILE__, __LINE__, 0,
                        _("Postfix error: "
                          "range arguments of different or incorrect type."));
            }

          }
          eval_size++;
          break;

        case OC_NORMAL:
          eval_stack[eval_size++] = _eval_normal(pf, normal, &i);
          break;
        case OC_PLANE:
          eval_stack[eval_size++] = _eval_plane(pf, coords, &i);
          break;
        case OC_BOX:
          eval_stack[eval_size++] = _eval_box(pf, coords, &i);
          break;
        case OC_CYLINDER:
          eval_stack[eval_size++] = _eval_cylinder(pf, coords, &i);
          break;
        case OC_SPHERE:
          eval_stack[eval_size++] = _eval_sphere(pf, coords, &i);
          break;

        case OC_GT:
          eval_stack[eval_size++] = _eval_coord_gt(pf, coords, &i);
          break;
        case OC_LT:
          eval_stack[eval_size++] = _eval_coord_lt(pf, coords, &i);
          break;
        case OC_GE:
          eval_stack[eval_size++] = _eval_coord_ge(pf, coords, &i);
          break;
        case OC_LE:
          eval_stack[eval_size++] = _eval_coord_le(pf, coords, &i);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _("Operator %s not currently implemented."),
                    _operator_name[oc]);

        } /* End of inside (operator) switch */

      }
      break;

    default:
      fvm_selector_postfix_dump(pf, 0, 0, NULL, NULL);
      bft_error(__FILE__, __LINE__, 0,
                _("Postfix evaluation error."));
    }

    if (eval_size == eval_max_size) {
      eval_max_size *= 2;
      if (eval_stack == _eval_stack) {
        BFT_MALLOC(eval_stack, eval_max_size, bool);
        memcpy(eval_stack, _eval_stack, BASE_STACK_SIZE*sizeof(bool));
      }
      else
        BFT_REALLOC(eval_stack, eval_max_size, bool);
    }

  } /* End of loop on postfix elements */

  if (eval_size != 1) {
    fvm_selector_postfix_dump(pf, 0, 0, NULL, NULL);
    bft_error(__FILE__, __LINE__, 0,
              _("Postfix evaluation error."));
  }

  retval = eval_stack[0];

  if (eval_stack != _eval_stack)
    BFT_FREE(eval_stack);

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump the contents of a postfix structure in human readable form
 *
 * parameters:
 *   pf           <-> pointer to postfix structure
 *   n_groups     <-- number of groups
 *   n_attributes <-- number of attributes
 *   group_name   <-- array group names (sorted)
 *   attribute    <-- array of attribute numbers (sorted)
 *----------------------------------------------------------------------------*/

void
fvm_selector_postfix_dump(const fvm_selector_postfix_t  *pf,
                          int                            n_groups,
                          int                            n_attributes,
                          const char                    *group_name[],
                          const int                      attribute[])
{
  size_t i = 0;

  bft_printf("\n"
             "Postfix expression dump:\n"
             "  Coordinates dependency:   %d\n"
             "  Normals dependency:       %d\n"
             "  Infix:\n"
             "    %s\n"
             "  Elements:\n",
             (int)pf->coords_dependency,
             (int)pf->normals_dependency,
             pf->infix);

  /* Dump postfix_string */

  i = 0;

  while (i < pf->size) {

    _postfix_type_t type = *((_postfix_type_t *)(pf->elements + i));

    i += _postfix_type_size;

    switch(type) {
    case PF_OPCODE:
      {
        _operator_code_t oc = *((_operator_code_t *)(pf->elements + i));
        bft_printf("    %s\n", _operator_name[oc]);
        i += _postfix_opcode_size;
      }
      break;
    case PF_GROUP_ID:
    case PF_ATTRIBUTE_ID:
    case PF_INT:
      {
        int val = *((int *)(pf->elements + i));
        if (type == PF_GROUP_ID) {
          if (val < 0)
            bft_printf("    %d (non-existing group id)\n", val);
          else if (n_groups > 0)
            bft_printf("    %d (group: \"%s\")\n", val, group_name[val]);
          else
            bft_printf("    %d (group id)\n", val);
        }
        else if (type == PF_ATTRIBUTE_ID) {
          if (val < 0)
            bft_printf("    %d (non-existing attribute id)\n", val);
          else if (n_attributes > 0)
            bft_printf("    %d (attribute: %d)\n", val, attribute[val]);
          else
            bft_printf("    %d (attribute id)\n", val);
        }
        else
          bft_printf("    %d\n", val);
        i += _postfix_int_size;
      }
      break;
    case PF_FLOAT:
      {
        double val = *((double *)(pf->elements + i));
        bft_printf("    %g\n", val);
        i += _postfix_float_size;
      }
      break;
    default:
      assert(0);
      break;
    }

  }

  if (pf->n_missing_operands > 0) {
    bft_printf("  Missing operands:         %d\n",
               pf->n_missing_operands);
    for (i = 0; i < (size_t)(pf->n_missing_operands); i++)
      bft_printf("    %s\n", pf->missing_operand[i]);
  }

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

