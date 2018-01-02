/*!
 * \file mei_parser_glob.c
 *
 * \brief Define global variables useful for the mathematical expression parsing
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

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
