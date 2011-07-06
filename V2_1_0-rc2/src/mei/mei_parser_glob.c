/*!
 * \file mei_parser_glob.c
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mei_node.h"

/*----------------------------------------------------------------------------
 * External global variables declarations for the parser
 *----------------------------------------------------------------------------*/

#include "mei_parser_glob.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Global variables
 *============================================================================*/

/* Line and column counters */

int  mei_glob_column = 0;
int  mei_glob_line = 0;

/* Error counter and markers */

int    mei_glob_ierr_list = 0;
int   *mei_glob_column_list = NULL;
int   *mei_glob_line_list = NULL;
char **mei_glob_label_list = NULL;

/* Start/end of the expression is set in these global variables */

char  *mei_glob_string_begin = NULL;
char  *mei_glob_string_end = NULL;

/* Pointer to the node returned by the parser */

mei_node_t *mei_glob_root = NULL;

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
