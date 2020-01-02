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

#include <sys/time.h>
#include <unistd.h>

#include "bft_printf.h"

static int
bft_printf_test(const char  *const format,
                va_list            arg_ptr)
{
  int ret;

  fprintf(stdout, "test printf replacement start\n");

  ret = vfprintf(stdout, format, arg_ptr);

  fprintf(stdout, "test printf replacement end\n");

  return ret;
}

static int
bft_printf_flush_test(void)
{
  fprintf(stdout, "test printf flush replacement start\n");

  fflush(NULL);

  fprintf(stdout, "test printf flush replacement end\n");

  return 0;
}

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  bft_printf_proxy_t *printf_proxy_save;
  bft_printf_flush_proxy_t *printf_flush_proxy_save;

  printf_proxy_save = bft_printf_proxy_get();
  printf_flush_proxy_save = bft_printf_flush_proxy_get();

  bft_printf_proxy_set(bft_printf_test);
  bft_printf_flush_proxy_set(bft_printf_flush_test);

  bft_printf("replacement printf(): %p\n", (void *)bft_printf_proxy_get());
  bft_printf_flush();

  bft_printf_proxy_set(printf_proxy_save);
  bft_printf_flush_proxy_set(printf_flush_proxy_save);

  bft_printf("saved printf(): %p\n", (void *)bft_printf_proxy_get());
  bft_printf_flush();

  exit (EXIT_SUCCESS);
}
