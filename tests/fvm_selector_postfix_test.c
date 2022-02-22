/*============================================================================
 * Mechanism for interpretation of a user syntax
 *============================================================================*/

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <cs_defs.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include <fvm_selector_postfix.h>

/*----------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  int n_groups = 6;
  int n_attributes = 6;

#if 1
  const char *group_names[] = {"entry", "entry_01", "entry_02",
                               "exchange", "outlet", "wall"};
#else
  const char *group_names[] = {"03", "04", "05", "12", "15", "16"};
#endif
  int   attributes[] = {3, 4, 5, 12, 15, 16};

  fvm_selector_postfix_t  *pf = NULL;

#if 0
  const char infix[] = "all [\t] all \\all \t '\\to' '\\\\all' "
    "\"other 'example' \" \t ";
  const char infix[] = "(((all[]))) and 3.1 >= z >= -2 or not (15 or entry)";
  const char infix[] = "range[04, 13, attribute]";
#else
  const char infix[] = "sphere[0, 0, 0, 2] and (not no_group[])";
#endif

  bft_mem_init(getenv("CS_MEM_LOG"));

  pf = fvm_selector_postfix_create(infix,
                                   n_groups,
                                   n_attributes,
                                   group_names,
                                   attributes);

  fvm_selector_postfix_dump(pf,
                            n_groups, n_attributes,
                            group_names, attributes);

  fvm_selector_postfix_destroy(&pf);

  bft_mem_end();

  exit(EXIT_SUCCESS);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
