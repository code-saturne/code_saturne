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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

#include "bft_error.h"
#include "bft_mem_usage.h"

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  size_t count[3];

  void *p1, *p2, *p3;

  /* BFT initialization and environment */

  bft_mem_usage_init();

  p1 = malloc(1000000);

  bft_mem_usage_init(); /* 2nd call to bft_mem_usage_init() to test if safe */

  p2 = malloc(1000000);

  p3 = malloc(100000);
  if (p1 != NULL) free(p1);
  if (p2 != NULL) free(p2);
  if (p3 != NULL) free(p3);

  printf("memory usage: %lu kB\n", (unsigned long) bft_mem_usage_pr_size());

  printf("max memory usage: %lu kB\n",
         (unsigned long) bft_mem_usage_max_pr_size());

  bft_mem_usage_n_calls(count);
  if (count[0] != 0) {
    printf("%lu calls to malloc\n"
           "%lu calls to realloc\n"
           "%lu calls to free\n",
           (unsigned long) count[0],
           (unsigned long) count[1],
           (unsigned long) count[2]);
  }

  bft_mem_usage_end();

  /* End */

  exit (EXIT_SUCCESS);
}
