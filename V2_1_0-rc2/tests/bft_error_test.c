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

static void
bft_error_handler_test(const char  *const file_name,
                       const int          line_num,
                       const int          sys_error_code,
                       const char  *const format,
                       va_list            arg_ptr)
{
  fprintf(stdout, "test error handler start\n");

  if (sys_error_code != 0)
    fprintf(stdout, "\n\nSystem error : %s\n", strerror(sys_error_code));

  fprintf(stderr, "\n%s:%d : Fatal error\n", file_name, line_num);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stdout, "\ntest error handler end\n");
}

int
main (int argc, char *argv[])
{
  bft_error_handler_t *errhandler_save;

  errhandler_save = bft_error_handler_get();

  bft_error_handler_set(bft_error_handler_test);

  bft_error(__FILE__, __LINE__, 1, "fake error, handler = %p",
            bft_error_handler_get());

  bft_error_handler_set(errhandler_save);

  bft_error(__FILE__, __LINE__, 1, "fake error, handler = %p",
            bft_error_handler_get());

  exit (EXIT_SUCCESS);
}
