/*============================================================================
 * Test program for mei
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem_usage.h>
#include <bft_mem.h>
#include <bft_error.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "mei_evaluate.h"

/*============================================================================
 * External function prototype
 *============================================================================*/

extern int graphik(mei_node_t*);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Handle errors
 *----------------------------------------------------------------------------*/

static void
_base_error_handler(const char  *file_name,
                    int          line_num,
                    int          sys_err_code,
                    const char  *format,
                    va_list      arg_ptr)
{
  fflush(stdout);
  fflush(stderr);

  fprintf(stderr, "\n");

  if (sys_err_code != 0)
    fprintf(stderr, "\nSystem error: %s\n", strerror(sys_err_code));

  fprintf(stderr, "\n%s:%d: Fatal error.\n\n", file_name, line_num);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stderr, "\n\n");

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 * Initialize memory handling
 *----------------------------------------------------------------------------*/

static void
_base_mem_init(void)
{
  char  *log_name = NULL;

  /* Initialization of memory usage count */

  bft_mem_usage_init();

  /* Initialization of memory management */

  log_name = getenv("CS_FIC_MEM");

  bft_mem_init(log_name);
}

/*----------------------------------------------------------------------------
 * Finalize memory handling
 *----------------------------------------------------------------------------*/

static void
_base_mem_finalize(void)
{
  int        ind_bil, itot;
  double     valreal[2];

  int        ind_val[2] = {1, 1};
  char       unit[]    = {'k', 'm', 'g', 't', 'p'};

  const char  * type_bil[] = {"Total measured memory usage:            ",
                              "Theoretical instrumented dynamic memory:"};

  /* Memory summary */

  printf("\nMemory usage summary:\n\n");

  valreal[0] = (double) bft_mem_usage_max_pr_size();
  valreal[1] = (double) bft_mem_size_max();

  /* We will ignore non-consistent measures */

  if (valreal[2] < valreal[1] || valreal[2] < valreal[3])
    ind_val[2] = 0;

  for (ind_bil = 0; ind_bil < 2; ind_bil++) {
    if (valreal[ind_bil] < 1.0)
      ind_val[ind_bil] = 0;
  }

 /* Similar handling for several instrumentation methods */

  for (ind_bil = 0; ind_bil < 2; ind_bil++) {

    /* If an instrumentation method returns an apparently consistent
       result, print it. */

    if (ind_val[ind_bil] == 1) {

      for (itot = 0;
           valreal[ind_bil] > 1024. && unit[itot] != 'p';
           itot++)
        valreal[ind_bil] /= 1024.;

      /* Impressions */

      printf ("  %s %12.3f %co\n",
              type_bil[ind_bil], valreal[ind_bil], unit[itot]);

    }
  }

  /* Finalize memory handling */

  bft_mem_end();

  /* Finalize memory usage count */

  bft_mem_usage_end();
}

/*============================================================================
 * Main program
 *============================================================================*/

int main(void)
{
  int iok;
  const char *v[] = {"X", "Y", "Z"};
  const char *w[] = {"yy", "zz"};

  mei_tree_t *e1;
  mei_tree_t *e2;
  hash_table_t *sym;

  /* Initialize memory handling */

  _base_mem_init();
  bft_error_handler_set(_base_error_handler);

  /* Check functionality */
  /*---------------------*/

  printf("\n------------------------------------------------------------------\n");

  /* return two empty interpreter */

  e1 = mei_tree_new("x=(+2.5); y=8.7; zz=y+x; yy  = zz+0.8; \n");
  e2 = mei_tree_new("y = cos(-pi) ; tt = K+y ; abs(tt);");


  /* try to build the two interpreter */

  printf("\nBuilding interpreter for: %s\n", "x=(+2.5); y=8.7; zz=y+x; yy  = zz+0.8; \n");
  if (mei_tree_builder(e1)) {
    mei_tree_destroy(e1);
    printf("failed...\n");
  } else
    printf("OK\n");

  printf("\n------------------------------------------------------------------\n");

  printf("\nBuilding interpreter for: %s\n", e2->string);
  if (mei_tree_builder(e2)) {
    mei_tree_destroy(e2);
    printf("failed...\n");
  } else
    printf("OK\n");

  printf("\n------------------------------------------------------------------\n");

  /* complete the symbol table with required */

  /* return an empty interpreter */
  e2 = mei_tree_new("y = cos(-pi) ; tt = K+y ; abs(tt);");

  mei_tree_insert(e2, "K", -5.3e2);

  /* try to re-build the interpreter --> memory problem */

  printf("\nBuilding interpreter for: %s\n", e2->string);
  if (mei_tree_builder(e2)) {
    mei_tree_destroy(e2);
    printf("failed...\n");
  } else
    printf("OK\n");

  printf("\n------------------------------------------------------------------\n");

  printf("\nFind symbols in: %s\n", e1->string);

  printf("Try to find %s\n", "toto");

  if (mei_tree_find_symbol(e1, "toto"))
    printf("not found...\n");
  else
    printf("found\n");

  printf("Try to find %s\n", "yy");

  if (mei_tree_find_symbol(e1, "yy"))
    printf("not found...\n");
  else
    printf("found\n");

  printf("\nFind symbols in: %s\n", e1->string);

  printf("Try to find X Y Z\n");

  if (mei_tree_find_symbols(e1, 3, v))
    printf("not found...\n");
  else
    printf("found\n");

  printf("Try to find yy zz\n");

  if (mei_tree_find_symbols(e1, 2, w))
    printf("not found...\n");
  else
    printf("found\n");

  printf("\n------------------------------------------------------------------\n");

  printf("\nInterpret expression: \n%s\n", e1->string);
  mei_evaluate(e1);
  printf("Evaluate: %s = %f\n", "zz", mei_tree_lookup(e1, "zz"));
  printf("Evaluate: %s = %f\n", "yy", mei_tree_lookup(e1, "yy"));

  printf("\nInterpret expression: \n%s\n", e2->string);
  printf("Evaluate: [%s] = %f\n", e2->string, mei_evaluate(e2));

  mei_tree_destroy(e1);
  mei_tree_destroy(e2);

  printf("\n------------------------------------------------------------------\n");

  printf("\nBuild a shared table of symbols: \n");

  sym = mei_table_symbols_new();
  mei_symbol_table_insert(sym, "x", 3.2);
  e1 = mei_tree_new_with_shared_symbols("y= 0.8; x+y;", sym);

  if (!mei_tree_builder(e1)) {
    /* graphik(e1->node); */
    printf("\nExpression: \n%s\n", e1->string);
    printf("Evaluate: %f\n", mei_evaluate(e1));
  }

  mei_tree_insert(e1, "x", 1.2);
  e2 = mei_tree_new_with_shared_symbols("z=0.6; k = x*2; k+z;", sym);

  if (!mei_tree_builder(e2)) {
    graphik(e2->node);
    printf("\nExpression: \n%s\n", e2->string);
    printf("Evaluate: %f\n", mei_evaluate(e2));
  }

  mei_tree_destroy(e1);
  mei_tree_destroy(e2);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("abc = -1.1; \nwhile (abc < 3) {\n  print abc;\n  abc = abc + 1;\n};");
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    printf("\nExpression: \n%s\n", e1->string);
    printf("Evaluate: %f\n", mei_evaluate(e1));
  }
  mei_tree_destroy(e1);


  printf("\n------------------------------------------------------------------\n");

  /* Check error detection */
  /*-----------------------*/

  e1 = mei_tree_new("x = A+3; A=5;");
  printf("\nExpression: \n%s\n", e1->string);
  iok = mei_tree_builder(e1);
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("cst = 4; # toto \n u = cos(cst);");
  printf("\nExpression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    mei_evaluate(e1);
    printf("Evaluate: u = %f\n", mei_tree_lookup(e1, "u"));
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("%u = cos(pi);");
  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    mei_evaluate(e1);
    printf("Evaluate: u = %f\n", mei_tree_lookup(e1, "u"));
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("v = max(pi,+3);");
  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("v = 1+2*+3;");
  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("v = +1000.9 -( 0.049245 +0.0041579*Temp.C)*Temp.C;");
  mei_tree_insert(e1, "Temp.C", 18.0);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("v = 1.5652e-3 + (-3.3003e-5 + 2.5135e-7*Temp-C)*Temp-C;");
  mei_tree_insert(e1, "Temp-C", 18.0);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  e1 = mei_tree_new("v = 1.3806e-4 + 3.1768e-7*Temp_C;");
  mei_tree_insert(e1, "Temp_C", 18.0);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
  }
  mei_tree_destroy(e1);


  printf("\n------------------------------------------------------------------\n");

  /* http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/logical.html */

  e1 = mei_tree_new("v = ! Something && Another;");
  mei_tree_insert(e1, "Something", 1);
  mei_tree_insert(e1, "Another", 0);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
    if (mei_tree_lookup(e1, "v")) printf("Test failed");
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  /* http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/logical.html */

  e1 = mei_tree_new("v = ! (Something && Another);");
  mei_tree_insert(e1, "Something", 1);
  mei_tree_insert(e1, "Another", 0);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
    if (!mei_tree_lookup(e1, "v")) printf("Test failed");
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  /* http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/logical.html */

  e1 = mei_tree_new("v = ! a || ! b && c;");
  mei_tree_insert(e1, "a", 1);
  mei_tree_insert(e1, "b", 1);
  mei_tree_insert(e1, "c", 0);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
    if (mei_tree_lookup(e1, "v")) printf("Test failed");
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  /* http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/logical.html */

  e1 = mei_tree_new("v = n^2 + 1 > 10 && ! n < 3;");
  mei_tree_insert(e1, "n", 4);

  printf("\nInterpret expression: \n%s\n", e1->string);
  if (!mei_tree_builder(e1)) {
    graphik(e1->node);
    mei_evaluate(e1);
    printf("Evaluate: v = %f\n", mei_tree_lookup(e1, "v"));
    if (!mei_tree_lookup(e1, "v")) printf("Test failed");
  }
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  printf("\n--------- WARNING: this test must be the last one ----------------\n");
  printf("\n--------- because it is corrupted the parser.     ----------------\n");

  e1 = mei_tree_new("cst = (2*(pi)/2.0); u = coos(cst); c=uu;");
  printf("\nExpression: \n%s\n", e1->string);
  iok = mei_tree_builder(e1);
  mei_tree_destroy(e1);

  printf("\n------------------------------------------------------------------\n");

  /* Finalization of memory management */
  _base_mem_finalize();

  return 0;
}

