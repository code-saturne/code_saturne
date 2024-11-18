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

#include "math.h"
#include "stdlib.h"

#include <iostream>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_base_accel.h"

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

void
cs_dispatch_test(void)
{
  const cs_lnum_t n = 100;

  //cs_dispatch_context ctx(cs_device_context(), {});
  cs_dispatch_context ctx;

  cs_alloc_mode_t amode = CS_ALLOC_HOST_DEVICE_SHARED;
  cs_real_t *a0, *a1;
  CS_MALLOC_HD(a0, n, cs_real_t, amode);
  CS_MALLOC_HD(a1, n, cs_real_t, amode);
  cs_real_3_t *a2;
  CS_MALLOC_HD(a2, n/10, cs_real_3_t, amode);

  cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);
#if defined(HAVE_ACCEL)
  cs_device_context &d_ctx = static_cast<cs_device_context&>(ctx);
#endif

  for (int i = 0; i < 3; i++) {

    if (i == 1) {
      ctx.set_use_gpu(false);
      ctx.set_n_min_per_cpu_thread(20);
    }
    else if (i == 2) {
      ctx.set_use_gpu(true);
    }

    ctx.parallel_for(n, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      cs_lnum_t c_id = ii;
      // Test to show whether we are on GPU or CPU...
#if defined( __CUDA_ARCH__) || defined( __SYCL_DEVICE_ONLY__)
      a0[ii] = c_id*0.1;
#else
      a0[ii] = -c_id*0.1;
#endif
      a1[ii] = cos(a0[ii]);
    });

    ctx.wait();

    for (cs_lnum_t ii = 0; ii < n/10; ii++) {
      std::cout << ii << " " << a0[ii] << " " << a1[ii] << std::endl;
    }

    for (cs_lnum_t ii = 0; ii < n/10; ii++) {
      a2[ii][0] = 0;
      a2[ii][1] = 0;
      a2[ii][2] = 0;
    }

    ctx.parallel_for(n, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
#if defined( __CUDA_ARCH__) || defined( __SYCL_DEVICE_ONLY__)
      cs_real_t s[3] = {0, -1, -2};
#else
      cs_real_t s[3] = {0, 1, 2};
#endif

      cs_dispatch_sum<3>(a2[ii/10], s, CS_DISPATCH_SUM_ATOMIC);
    });

    ctx.wait();

    for (cs_lnum_t ii = 0; ii < n/10; ii++) {
      std::cout << ii << " " << a2[ii][0]
                      << " " << a2[ii][1]
                      << " " << a2[ii][2] << std::endl;
    }

    double s1 = 0;
    ctx.parallel_for_reduce_sum
      (n, s1, [=] CS_F_HOST_DEVICE (cs_lnum_t ii,
                                    CS_DISPATCH_SUM_DOUBLE &sum) {
#if defined( __CUDA_ARCH__) || defined( __SYCL_DEVICE_ONLY__)
      {sum += (double)ii;}
#else
      {sum += -(double)ii;}
#endif
    });

    ctx.wait();

    std::cout << "reduction (sum)" << s1 << std::endl;
  }

#ifdef __NVCC__
  std::cout << "device_id " << cs_base_cuda_get_device() << std::endl;
#endif

  CS_FREE_HD(a0);
  CS_FREE_HD(a1);
}

/*----------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  cs_dispatch_test();

  exit(EXIT_SUCCESS);
}
