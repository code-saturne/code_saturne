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

#include "bft_backtrace.h"

/* function prototypes */

void
level_1(void);

void
level_2(void);

/* function definitions */

static void
level_2_static(void)
{
  level_2();
}

void
level_1(void)
{
  level_2_static();
}

void
level_2(void)
{
  int i;
  bft_backtrace_t *bt = bft_backtrace_create();

  printf("backtrace:\n");
  for (i = 0; i < bft_backtrace_size(bt); i++) {
    printf("  %d: %s [%s]\n", i, bft_backtrace_function(bt, i),
           bft_backtrace_address(bt, i));
  }

  bft_backtrace_demangle(bt);

  printf("\n"
         "backtrace after demangle (should be modified from\n"
         "previous only if using C++ compiler and libs):\n");
  for (i = 0; i < bft_backtrace_size(bt); i++) {
    printf("  %d: %s [%s]\n", i, bft_backtrace_function(bt, i),
           bft_backtrace_address(bt, i));
  }

  bft_backtrace_destroy(bt);
}

int
main (int argc, char *argv[])
{
  level_1();

  exit (EXIT_SUCCESS);
}
