#ifndef __MEI_PARSER_GLOB_H__
#define __MEI_PARSER_GLOB_H__

/*!
 * \file mei_parser_glob.h
 *
 * \brief Define global variables useful for the mathematical expression parsing
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
 * Scan the expression
 *----------------------------------------------------------------------------*/

extern int yylex(void);

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

