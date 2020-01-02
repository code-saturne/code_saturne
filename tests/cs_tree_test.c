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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_tree.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  cs_tree_node_t *n;

  /* Build a small tree */
  cs_tree_node_t  *root = cs_tree_node_create("code_saturne");
  cs_tree_node_set_value_str(root, CS_APP_VERSION);
  cs_tree_add_child(root, "fields");
  cs_tree_add_child(root, "properties");
  cs_tree_node_t  *neqs = cs_tree_add_child(root, "equations");

  /* Set the first equation */
  cs_tree_node_t  *neq1 = cs_tree_add_child(neqs, "eq1");
  n = cs_tree_add_child_str(neq1, "space_scheme", "cdovb");
  n = cs_tree_add_child_bool(neq1, "do_diffusion", true);
  n = cs_tree_add_child_str(neq1, "diffusion_hodge", "dga");
  cs_tree_node_t  *m = cs_tree_add_child_str(neq1, "iterative_solver", "cg");
  n = cs_tree_add_child_real(m, "epsilon", 1e-12);
  n = cs_tree_add_child_int(m, "n_max_iterations", 2500);
  n = cs_tree_add_child_str(m, "precond", "amg");

  /* Set the second equation */
  cs_tree_node_t  *neq2 = cs_tree_add_child(neqs, "eq2");
  n = cs_tree_add_child_str(neq2, "space_scheme", "cdofb");
  n = cs_tree_add_child_bool(neq2, "do_diffusion", true);
  n = cs_tree_add_child_bool(neq2, "do_time", true);

  /* Set the third equation */
  cs_tree_node_t  *neq3 = cs_tree_add_child(neqs, "eq3");
  n = cs_tree_add_child_str(neq3, "space_scheme", "hho_p1");
  n = cs_tree_add_child_bool(neq3, "do_reaction", true);
  n = cs_tree_add_child_bool(neq3, "do_diffusion", true);

  /* Dump the tree */
  cs_tree_dump(CS_LOG_DEFAULT, 0, root);

  /* Search for nodes */

  n = cs_tree_get_or_add_node(root, "/equations");
  cs_tree_node_dump(CS_LOG_DEFAULT, 0, n);
  n = cs_tree_get_or_add_node(root, "equations/eq1/iterative_solver/precond");
  cs_tree_node_dump(CS_LOG_DEFAULT, 0, n);

  cs_tree_node_t  *n1 = cs_tree_get_node(root, "equations/eq1");
  n = cs_tree_get_or_add_node(n1, "iterative_solver/precond");
  cs_tree_node_dump(CS_LOG_DEFAULT, 0, n);
  cs_log_printf(CS_LOG_DEFAULT, "eq1 precond = %s\n",
                cs_tree_node_get_value_str(n));
  cs_tree_node_dump(CS_LOG_DEFAULT, 0, n1);

  /* Free all the tree structure */
  cs_tree_node_free(&root);

  exit(EXIT_SUCCESS);
}
