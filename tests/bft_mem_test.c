/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include <sys/time.h>
#include <unistd.h>

#include "bft_error.h"
#include "bft_mem_usage.h"
#include "bft_mem.h"

static void
bft_mem_error_handler_test(const char     *const file_name,
                           const int             line_num,
                           const int             sys_error_code,
                           const char     *const format,
                           va_list               arg_ptr)
{
  CS_UNUSED(file_name);
  CS_UNUSED(line_num);
  CS_UNUSED(sys_error_code);
  CS_UNUSED(format);
  CS_UNUSED(arg_ptr);

  fprintf(stderr, "test memory error handler (empty).\n");
}

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  bft_error_handler_t *errhandler_save;

  void *p1, *p2, *p3, *p4;

  /* BFT initialization and environment */

  bft_mem_usage_init();

  bft_mem_init("bft_mem_log_file");

  errhandler_save = bft_mem_error_handler_get();

  bft_mem_error_handler_set(bft_mem_error_handler_test);
  printf("test memory error handler set\n");

  BFT_MALLOC(p1, 100000, long);
  printf("p1 = %p\n", p1);
  BFT_MALLOC(p2, 100000, double);
  printf("p2 = %p\n", p2);

  p3 = NULL;
  BFT_REALLOC(p3, 100000, double);
  printf("p3 = %p\n", p3);
  BFT_REALLOC(p3, 10000, double);
  printf("p3 = %p\n", p3);

  BFT_MALLOC(p4, 5000000000, double);
  printf("p4 = %p\n", p4);

  printf("default memory error handler set\n");
  bft_mem_error_handler_set(errhandler_save);

  BFT_FREE(p1);
  BFT_FREE(p2);
  BFT_MALLOC(p2, 100000, double);
  printf("p2 = %p\n", p2);
  BFT_FREE(p2);
  BFT_REALLOC(p3, 0, double);
  printf("p3 = %p\n", p3);

  if (bft_mem_have_memalign() == 1) {
    void *pa;
    BFT_MEMALIGN(pa, 128, 100, double);
    printf("pa (aligned 128) = %p\n", pa);
    BFT_FREE(pa);
  }

  bft_mem_end();

  printf("max memory usage: %lu kB\n",
         (unsigned long) bft_mem_usage_max_pr_size());

  BFT_MALLOC(p1, 10000, long);
  printf("p1 = %p\n", p1);
  BFT_FREE(p1);
  printf("p1 = %p\n", p1);
  BFT_MALLOC(p1, 1000000000, double);
  printf("p1 = %p\n", p1);

  /* End */

  exit (EXIT_SUCCESS);
}
