/*============================================================================
 * Unit test for fvm_group.c and fvm_selector.c;
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include <stdlib.h>

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "fvm_group.h"
#include "fvm_selector.h"

/*---------------------------------------------------------------------------*/

static void
test_1 (void)
{
  fvm_group_class_set_t *gcset = NULL;

  const char *grp_1[] = {"group_1", "g2", "g3", "12", "56"};
  const char *grp_2[] = {"group_4", "g2", "g5", "group_6"};
  const char *grp_3[] = {"g7", "g8", "57"};
  const char *grp_4[] = {"12"};

  gcset = fvm_group_class_set_create();

  fvm_group_class_set_add(gcset, 5, grp_1);
  fvm_group_class_set_add(gcset, 4, grp_2);
  fvm_group_class_set_add(gcset, 3, grp_3);
  fvm_group_class_set_add(gcset, 1, grp_4);

  fvm_group_class_set_dump(gcset);

  gcset = fvm_group_class_set_destroy(gcset);
}


static void
test_2 (void)
{
  int ii;
  int criteria_id, n_missing;

  fvm_group_class_set_t *gcset = NULL;
  fvm_selector_t *s = NULL;

  cs_lnum_t n_se = 0;
  cs_lnum_t se[12];

  const char *att_02[] = {"2"};   /* 1 */
  const char *att_06[] = {"6"};   /* 2 */
  const char *att_10[] = {"10"};  /* 3 */
  const char *att_01[] = {"1"};   /* 4 */
  const char *att_05[] = {"5"};   /* 5 */
  const char *att_03[] = {"3"};   /* 6 */
  const char *att_11[] = {"11"};  /* 7 */

  const int f_gc_id[] = {5, 1, 2, 7, 7, 7, 3, 3, 4, 4, 8, 6};

  const double coords[] = { 1.0, 0.0, 0.0,
                            2.0, 0.0, 0.0,
                            3.0, 0.0, 0.0,
                            4.0, 0.0, 0.0,
                            5.0, 0.0, 0.0,
                            6.0, 0.0, 0.0,
                            7.0, 0.0, 0.0,
                            8.0, 0.0, 0.0,
                            9.0, 0.0, 0.0,
                           10.0, 0.0, 0.0,
                           11.0, 0.0, 0.0,
                           12.0, 0.0, 0.0};

  const double norms[] = {1.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          -1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 1.0, 0.0};

  gcset = fvm_group_class_set_create();

  fvm_group_class_set_add(gcset, 1, att_02);
  fvm_group_class_set_add(gcset, 1, att_06);
  fvm_group_class_set_add(gcset, 1, att_10);
  fvm_group_class_set_add(gcset, 1, att_01);
  fvm_group_class_set_add(gcset, 1, att_05);
  fvm_group_class_set_add(gcset, 1, att_03);
  fvm_group_class_set_add(gcset, 1, att_11);
  fvm_group_class_set_add(gcset, 0, NULL);

  fvm_group_class_set_dump(gcset);

  s = fvm_selector_create(3,
                          12,
                          gcset,
                          f_gc_id,
                          1,
                          coords,
                          norms);

  fvm_selector_dump(s);

  criteria_id
    = fvm_selector_get_list(s,
                            "11 or (1, inlet; outlet and 6)",
                            1,
                            &n_se,
                            se);

  bft_printf("selection:\n");
  for (ii = 0; ii < n_se; ii++)
    bft_printf(" %d", se[ii]);
  bft_printf("\n\n");

  n_missing = fvm_selector_n_missing(s, criteria_id);
  bft_printf("missing operands (%d): \n", n_missing);
  for (ii = 0; ii < n_missing; ii++)
    bft_printf("  \"%s\"\n",
               fvm_selector_get_missing(s, criteria_id, ii));

  criteria_id
    = fvm_selector_get_list(s,
                            "x < 5",
                            1,
                            &n_se,
                            se);

  bft_printf("selection:\n");
  for (ii = 0; ii < n_se; ii++)
    bft_printf(" %d", se[ii]);
  bft_printf("\n\n");

  criteria_id
    = fvm_selector_get_list(s,
                            "range[1, 3, attribute]",
                            1,
                            &n_se,
                            se);

  bft_printf("selection:\n");
  for (ii = 0; ii < n_se; ii++)
    bft_printf(" %d", se[ii]);
  bft_printf("\n\n");

  criteria_id
    = fvm_selector_get_list(s,
                            "sphere[4.1, 0, 0, 2] and (not no_group[])",
                            1,
                            &n_se,
                            se);

  bft_printf("selection:\n");
  for (ii = 0; ii < n_se; ii++)
    bft_printf(" %d", se[ii]);
  bft_printf("\n\n");

  /* fvm_selector_dump(s); */

  s = fvm_selector_destroy(s);

  gcset = fvm_group_class_set_destroy(gcset);
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  bft_mem_init(getenv("FVM_MEM_TRACE"));

  test_1();

  test_2();

  bft_mem_end();

  exit (EXIT_SUCCESS);
}
