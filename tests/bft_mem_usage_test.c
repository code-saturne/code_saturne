/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004  EDF

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
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

#include "bft_error.h"
#include "bft_mem_usage.h"
#include "bft_sys_info.h"

int
main (int argc, char *argv[])
{
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

  bft_mem_usage_end();

  /* End */

  exit (EXIT_SUCCESS);
}
