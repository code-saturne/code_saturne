#ifndef __MEI_PARSER_GLOB_H__
#define __MEI_PARSER_GLOB_H__

/*!
 * \file mei_parser_glob.h
 *
 * \brief Define global variables usefull for the mathematical expression parsing
 */

/*
  This file is part of the "Mathematical Expression Interpreter" library.

  Copyright (C) 2008-2010  EDF

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

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "mei_node.h"

/*----------------------------------------------------------------------------
 * Max length of the variable name
 *----------------------------------------------------------------------------*/

#define TOKEN_SIZE 200

/*============================================================================
 * Global variables
 *============================================================================*/

/* Line and column counters */

extern int  mei_glob_column;
extern int  mei_glob_line;

/* Error counter and markers */

extern int     mei_glob_ierr_list;
extern int    *mei_glob_column_list;
extern int    *mei_glob_line_list;
extern char  **mei_glob_label_list;

/* Start/end of the expression is set in these global variables */

extern char  *mei_glob_string_begin;
extern char  *mei_glob_string_end;

/* Root node */

extern mei_node_t  *mei_glob_root;

/*============================================================================
 * Public function prototypes defined by Yacc
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Parse the expression
 *----------------------------------------------------------------------------*/

int yyparse(void);

/*----------------------------------------------------------------------------
 * Give an error message if needed
 *----------------------------------------------------------------------------*/

void yyerror(const char *s);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MEI_PARSER_GLOB_H__ */

