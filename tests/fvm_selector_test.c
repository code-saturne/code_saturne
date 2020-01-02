/*============================================================================
 * Unit test for fvm_group.c and fvm_selector.c;
 *============================================================================*/

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

static void
test_3 (void)
{
  int ii;

  fvm_group_class_set_t *gcset = NULL;
  fvm_selector_t *s = NULL;

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

  int n_groups = 201;
  const char *group_names[] = {
    "EBA1102BW",
    "EBA1103BW",
    "EBA1104BW",
    "EBA1105BW",
    "EBA1106BW",
    "EBA1201BH",
    "EBA2101BW",
    "EBA2102BW",
    "EBA2103BH",
    "EBA2104BW",
    "EBA2105BW",
    "EBA2107BW",
    "EBA2113BW",
    "EBA2116BW",
    "EBA2152BH",
    "EBA2153BH",
    "EBA2191BH",
    "EBA2191BW",
    "EBA2192BH",
    "EBA2192BW",
    "EBA2196BW",
    "EBA2201BH",
    "EBA2202BH",
    "EBA2203BH",
    "EBA2204BH",
    "EBA2205BH",
    "EBA2206BH",
    "EBA2207BH",
    "EBA2208BH",
    "EBA2210BH",
    "EBA2211BH",
    "EBA2214BH",
    "EBA2215BH",
    "EBA2221BH",
    "EBA2222BH",
    "EBA2223BH",
    "EBA2224BH",
    "EBA2225BH",
    "EBA2227BH",
    "EBA2252BW",
    "EBA2253BW",
    "EVF1110RA",
    "EVF1421BW",
    "EVF1430RA",
    "EVR1101BH_GV1",
    "EVR1101BW_GV1",
    "EVR1111BH",
    "EVR1111BW",
    "EVR1112BW",
    "EVR1113BW",
    "EVR1115BW",
    "EVR1116BW",
    "EVR1118BW",
    "EVR1119BW",
    "EVR1120BH",
    "EVR1120BW",
    "EVR1121BW",
    "EVR1122BW",
    "EVR1123BW",
    "EVR1124BW",
    "EVR1125BW",
    "EVR1126BW",
    "EVR1130BH",
    "EVR1130BW",
    "EVR1131BW",
    "EVR1132BW",
    "EVR1133BW",
    "EVR1140BH",
    "EVR1140BW",
    "EVR1150BH",
    "EVR1150BW",
    "EVR1150ZV-BH",
    "EVR1150ZV-BW",
    "EVR1151BH",
    "EVR1151BW",
    "EVR1152BH",
    "EVR1152BW",
    "EVR1160ZV-BH",
    "EVR1160ZV-BW",
    "EVR1190BW",
    "EVR1192BW",
    "EVR1193BW",
    "EVR1194BW",
    "EVR1201BH_GV4",
    "EVR1201BW",
    "EVR1201BW_GV4",
    "EVR1202BW",
    "EVR1203BW",
    "EVR1204BW",
    "EVR1205BW",
    "EVR1206BW",
    "EVR1211BH",
    "EVR1211BW",
    "EVR1212BW",
    "EVR1213BW",
    "EVR1220BW",
    "EVR1221BW",
    "EVR1222BW",
    "EVR1230BH",
    "EVR1231BW",
    "EVR1232BW",
    "EVR1240BH",
    "EVR1240BW",
    "EVR1250BH",
    "EVR1250BW",
    "EVR1250ZV-BH",
    "EVR1250ZV-BW",
    "EVR1251BH",
    "EVR1251BW",
    "EVR1254BHy",
    "EVR1254BW",
    "EVR1255BH",
    "EVR1255BW",
    "EVR1260ZV-BH",
    "EVR1260ZV-BW",
    "EVR1290BW",
    "EVR1293BW",
    "EVR1294BW",
    "EVR1295BW",
    "EVR1296BH",
    "EVR1296BW",
    "EVR2110ZV-BH",
    "EVR2120ZV-BH",
    "EVR2156TY",
    "EVR2157TY",
    "EVR2158TY",
    "EVR2159TY",
    "EVR2160TY",
    "EVR2161TY",
    "EVR2162TY",
    "EVR2163TY",
    "EVR2196BH",
    "EVR2196BW",
    "EVR2197BH",
    "EVR2197BW",
    "EVR2198BH",
    "EVR2198BW",
    "EVR2199BH",
    "EVR2199BW",
    "EVR2210ZV-BH",
    "EVR2220ZV-BH",
    "EVR2256TY",
    "EVR2257TY",
    "EVR2258TY",
    "EVR2259TY",
    "EVR2260TY",
    "EVR2261TY",
    "EVR2262TY",
    "EVR2263TY",
    "EVR3110BH",
    "EVR3110BW",
    "EVR3120BH",
    "EVR3120BW",
    "EVR3210BH",
    "EVR3210BW",
    "EVR3220BH",
    "EVR3220BW",
    "EVR3230BH",
    "EVR3230BW",
    "EVR3240BH",
    "EVR3240BW",
    "EVR3310BH",
    "EVR3310BW",
    "EVR3311BW",
    "EVR3320BH",
    "EVR3320BW",
    "EVR3410BW",
    "EVR3420BW",
    "EVR3430BH",
    "EVR3430BW",
    "EVR3440BH",
    "EVR3440BW",
    "EVR_P1BH",
    "EVR_P1BW",
    "EVR_P2BH",
    "EVR_P2BW",
    "EVR_P3BH",
    "EVR_P3BW",
    "EVR_P4BH",
    "EVR_P4BW",
    "EVR_P5BH",
    "EVR_P5BW",
    "EVR_P6BH",
    "EVR_P6BW",
    "KRC1101TY",
    "KRC1102TY",
    "KRC1103TY",
    "KRC1104TY",
    "KRC1105TY",
    "KRC1106TY",
    "KRC1201TY",
    "KRC1202TY",
    "KRC1203TY",
    "KRC1204TY",
    "KRC1205TY",
    "KRC1206TY",
    "KRC1601TY",
    "KRC1602TY",
    "fluid",
    "wall",
    "wall2"};

  gcset = fvm_group_class_set_create();

  for (ii = 0; ii < n_groups; ii++)
    fvm_group_class_set_add(gcset, 1, group_names + ii);

  fvm_group_class_set_dump(gcset);

  const char criteria[] = "not (EBA1102BW or EBA1103BW or EBA1104BW or EBA1105BW or EBA1106BW or EBA1201BH or EBA2101BW or EBA2102BW or EBA2103BH or EBA2104BW or EBA2105BW or EBA2107BW or EBA2113BW or EBA2116BW or EBA2192BH or EBA2192BW or EBA2196BW or EBA2201BH or EBA2202BH or EBA2203BH or EBA2204BH or EBA2205BH or EBA2206BH or EBA2207BH or EBA2208BH or EBA2210BH or EBA2211BH or EBA2214BH or EBA2215BH or EBA2221BH or EBA2222BH or EBA2223BH or EBA2224BH or EBA2225BH or EBA2227BH or EVF1110RA or EVF1421BW or EVF1430RA or EVR1101BH_GV1 or EVR1101BW_GV1 or EVR1111BH or EVR1111BW or EVR1112BW or EVR1113BW or EVR1115BW or EVR1116BW or EVR1118BW or EVR1119BW or EVR1120BH or EVR1120BW or EVR1121BW or EVR1122BW or EVR1123BW or EVR1124BW or EVR1125BW or EVR1126BW or EVR1131BW or EVR1132BW or EVR1133BW or EVR1150ZV-BH or EVR1150ZV-BW or EVR1160ZV-BH or EVR1160ZV-BW or EVR1190BW or EVR1192BW or EVR1193BW or EVR1194BW or EVR1201BH_GV4 or EVR1201BW or EVR1201BW_GV4 or EVR1202BW or EVR1203BW or EVR1204BW or EVR1205BW or EVR1206BW or EVR1211BH or EVR1211BW or EVR1212BW or EVR1213BW or EVR1220BW or EVR1221BW or EVR1222BW or EVR1230BH or EVR1231BW or EVR1232BW or EVR1250ZV-BH or EVR1250ZV-BW or EVR1260ZV-BH or EVR1260ZV-BW or EVR1290BW or EVR1293BW or EVR1294BW or EVR1295BW or EVR2110ZV-BH or EVR2120ZV-BH or EVR2156TY or EVR2157TY or EVR2158TY or EVR2159TY or EVR2160TY or EVR2161TY or EVR2162TY or EVR2163TY or EVR2210ZV-BH or EVR2220ZV-BH or EVR2256TY or EVR2257TY or EVR2258TY or EVR2259TY or EVR2260TY or EVR2261TY or EVR2262TY or EVR2263TY or EVR3110BH or EVR3110BW or EVR3120BH or EVR3120BW or EVR3210BH or EVR3210BW or EVR3220BH or EVR3220BW or EVR3230BH or EVR3230BW or EVR3240BH or EVR3240BW or EVR3310BH or EVR3310BW or EVR3311BW or EVR3320BH or EVR3320BW or EVR3410BW or EVR3420BW or EVR3430BH or EVR3430BW or EVR3440BH or EVR3440BW or EVR_P1BH or EVR_P1BW or EVR_P2BH or EVR_P2BW or EVR_P3BH or EVR_P3BW or EVR_P4BH or EVR_P4BW or EVR_P5BH or EVR_P5BW or EVR_P6BH or EVR_P6BW or KRC1101TY or KRC1102TY or KRC1103TY or KRC1104TY or KRC1105TY or KRC1106TY or KRC1201TY or KRC1202TY or KRC1203TY or KRC1204TY or KRC1205TY or KRC1206TY or KRC1601TY or KRC1602TY or wall or wall2)";
  s = fvm_selector_create(3,
                          12,
                          gcset,
                          f_gc_id,
                          1,
                          coords,
                          norms);

  fvm_selector_dump(s);

  int n_gc;
  int gc[205];

  int criteria_id = fvm_selector_get_gc_list(s, criteria, &n_gc, gc);

  bft_printf("groups selection:\n");
  for (ii = 0; ii < n_gc; ii++)
    bft_printf(" %d", gc[ii]);
  bft_printf("\n\n");
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  bft_mem_init(getenv("CS_MEM_LOG"));

  test_1();

  test_2();

  test_3();

  bft_mem_end();

  exit (EXIT_SUCCESS);
}
