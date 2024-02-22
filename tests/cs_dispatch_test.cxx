/*============================================================================
 * Dispatch test, CUDA implementations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "stdlib.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_dispatch.h"

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local macros
 *============================================================================*/

/*=============================================================================
 * Local definitions
 *============================================================================*/

void
cs_dispatch_test_cuda(void);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Test dispatch class.
 *----------------------------------------------------------------------------*/

static void
_cs_dispatch_test(void)
{
  const cs_lnum_t n = 100;

  cs_real_t *a0, *a1;
  CS_MALLOC_HD(a0, n, cs_real_t, CS_ALLOC_HOST);

  if (cs_get_device_id() > -1) {
    CS_MALLOC_HD(a1, n, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);
  }
  else {
    CS_MALLOC_HD(a1, n, cs_real_t, CS_ALLOC_HOST);
  }

#if 0
  for (cs_lnum_t ii = 0; ii < n; ii++) {
    cout << ii << " " << a0[ii] " " << a1[ii] << endl;
  }
#endif

  CS_FREE_HD(a0);
  CS_FREE_HD(a1);
}

/*----------------------------------------------------------------------------*/

extern "C" int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  _cs_dispatch_test();
  cs_dispatch_test_cuda();

  exit(EXIT_SUCCESS);
}
