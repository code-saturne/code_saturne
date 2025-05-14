/*============================================================================
 * Dispatch queue test, CUDA implementations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

#include "base/cs_mem.h"
#include "math.h"
#include "stdlib.h"

#include <climits>
#include <iostream>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_base_accel.h"
#include "base/cs_math.h"

#include "base/cs_dispatch_queue.h"

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local macros
 *============================================================================*/

/*=============================================================================
 * Local definitions
 *============================================================================*/

struct cs_data_1r_2i { // struct: class with only public members

  // Members
  cs_real_t r[1];
  cs_lnum_t i[2];
#if 0
  // Constructors

  cs_data_1r_2i(void) {}

  cs_data_1r_2i(cs_real_t a_r0,
                cs_lnum_t a_i0,
                cs_lnum_t a_i1)
    : r{a_r0}, i{a_i0, a_i1} {}
#endif
};

struct cs_reduce_sum1r1i_max1i {    // struct: class with only public members
  using T = cs_data_1r_2i;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] = 0.;
    a.i[0] = 0;
    a.i[1] = -INT_MAX;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] += b.r[0];
    a.i[0] += b.i[0];
    a.i[1] = CS_MAX(a.i[1], b.i[1]);
  }
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Test dispatch_queue class.
 *----------------------------------------------------------------------------*/

static void
_cs_dispatch_queue_test(void)
{
  const cs_lnum_t n = 100, n_sum = 100;

  cs_dispatch_queue queue;

  cs_alloc_mode_t amode = CS_ALLOC_HOST_DEVICE_SHARED;
  cs_real_t *a0, *a1;
  CS_MALLOC_HD(a0, n, cs_real_t, amode);
  CS_MALLOC_HD(a1, n, cs_real_t, amode);
  cs_real_3_t *a2;
  CS_MALLOC_HD(a2, n/10, cs_real_3_t, amode);

  // cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);
#if defined(HAVE_ACCEL)
  // cs_device_context &d_ctx = static_cast<cs_device_context&>(ctx);
#endif

  for (int test_iteration = 0; test_iteration < 3; test_iteration++) {
    if (test_iteration == 1) {
      queue.initializer_context.set_use_gpu(false);
      queue.initializer_context.set_n_min_per_cpu_thread(20);
    }
    else if (test_iteration == 2) {
      queue.initializer_context.set_use_gpu(true);
    }

    auto task_0 = queue.parallel_for(n, [=] CS_F_HOST_DEVICE(cs_lnum_t i) {
      cs_lnum_t c_id = i;
      // Test to show whether we are on GPU or CPU...
#if defined(__CUDA_ARCH__) || defined(__SYCL_DEVICE_ONLY__)
      a0[i] = c_id * 0.1;
#else
      a0[i] = -c_id*0.1;
#endif
      a1[i] = cos(a0[i]);
    });

    // ctx.wait();

    auto task_1 = queue.single_task({ task_0 }, [&]() {
      for (cs_lnum_t i = 0; i < n / 10; i++) {
        std::cout << i << " " << a0[i] << " " << a1[i] << std::endl;
      }

      for (cs_lnum_t i = 0; i < n / 10; i++) {
        a2[i][0] = 0;
        a2[i][1] = 0;
        a2[i][2] = 0;
      }
    });

    auto task_2 =
      queue.parallel_for(n, { task_1 }, [=] CS_F_HOST_DEVICE(cs_lnum_t i) {
#if defined(__CUDA_ARCH__) || defined(__SYCL_DEVICE_ONLY__)
        cs_real_t s[3] = { 0, -1, -2 };
#else
        cs_real_t s[3] = {0, 1, 2};
#endif
        cs_dispatch_sum<3>(a2[i / 10], s, CS_DISPATCH_SUM_ATOMIC);
      });

    cs_real_t pi = cs_math_pi;

    double r_sum  = 0;
    auto   task_3 = queue.single_task({ task_2 }, [&]() {
      for (cs_lnum_t i = 0; i < n / 10; i++) {
        std::cout << i << " " << a2[i][0] << " " << a2[i][1] << " " << a2[i][2]
                  << std::endl;
      }

      // reference sum
      for (cs_lnum_t i = 0; i < n_sum; i++) {
        cs_real_t x = (i % 10 - 3) * pi;
        r_sum -= (double)x;
      };
    });

    double s1     = 0;
    auto   task_4 = queue.parallel_for_reduce_sum(
      n_sum,
      { task_3 },
      s1,
      [=] CS_F_HOST_DEVICE(cs_lnum_t i,
                           CS_DISPATCH_REDUCER_TYPE(double) & sum) {
        cs_real_t x = (i % 10 - 3) * pi;
#if defined(__CUDA_ARCH__) || defined(__SYCL_DEVICE_ONLY__)
        {
          sum += (double)x;
        }
#else
        {
          sum += -(double)x;
        }
#endif
      });

    auto task_5 = queue.single_task({ task_4 }, [&] {
      std::cout << "reduction (sum) " << s1 << " (ref " << r_sum << ")"
                << std::endl;
    });

    struct cs_data_1r_2i           rd;
    struct cs_reduce_sum1r1i_max1i reducer;

    auto task_6 =
      queue.parallel_for_reduce(n_sum,
                                { task_5 },
                                rd,
                                reducer,
                                [=] CS_F_HOST_DEVICE(cs_lnum_t i,
                                                     cs_data_1r_2i & res) {
#if defined(__CUDA_ARCH__) || defined(__SYCL_DEVICE_ONLY__)
                                  cs_real_t x = (i % 10 - 3) * pi;
#else
        cs_real_t x = -(i%10 - 3)*pi;
#endif
                                  cs_lnum_t y = (i % 10 + 1);

                                  // The following is not allowed with CUDA
                                  // (or needs __host__ __device__ constructor).
                                  // res = cs_data_1r_2i(x, y, y);
                                  res.r[0] = x, res.i[0] = y, res.i[1] = y;
                                });

    task_6.wait();

    std::cout << "reduction (mixed) " << rd.r[0] << " " << rd.i[0] << " "
              << rd.i[1] << std::endl;
  }

#ifdef __NVCC__
  std::cout << "device_id " << cs_base_cuda_get_device() << std::endl;
#endif

  CS_FREE_HD(a0);
  CS_FREE_HD(a1);
}

inline void
_cs_sycl_like_example()
{
  // Inspired by:
  // https://enccs.github.io/sycl-workshop/task-graphs-synchronization/#how-to-specify-dependencies

  cs_dispatch_queue Q;
  Q.initializer_context.set_use_gpu(true);

  constexpr std::size_t N = 16 * 1024 * 1024;
  int                  *a, *b;

  CS_MALLOC_HD(a, N, int, CS_ALLOC_HOST_DEVICE_SHARED);
  CS_MALLOC_HD(b, N, int, CS_ALLOC_DEVICE);

  {
    // Task A
    auto task_a =
      Q.parallel_for(N, [=] CS_F_HOST_DEVICE(std::size_t id) { a[id] = 1; });

    // Task B
    auto task_b =
      Q.parallel_for(N, [=] CS_F_HOST_DEVICE(std::size_t id) { b[id] = 2; });

    // Task C
    auto task_c =
      Q.parallel_for(N,
                     { task_a, task_b },
                     [=] CS_F_HOST_DEVICE(std::size_t id) { a[id] += b[id]; });

    // Task D
    int  value  = 3;
    auto task_d = Q.single_task(
      { task_c },
      [=](int *data) -> void {
        for (int i = 1; i < N; i++) {
          data[0] += data[i];
        }

        data[0] /= value;
      },
      a);

    task_d.wait();
  }

  CS_FREE_HD(a);
  CS_FREE_HD(b);
}

/*----------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  cs_glob_n_threads = omp_get_max_threads();
#endif

#if defined(HAVE_CUDA)
  cs_base_cuda_select_default_device();
#endif
#if defined(HAVE_SYCL)
  cs_sycl_select_default_device();
#endif
#if defined(HAVE_OPENMP_TARGET)
  cs_omp_target_select_default_device();
#endif

  _cs_dispatch_queue_test();

  _cs_sycl_like_example();

  exit(EXIT_SUCCESS);
}

/*----------------------------------------------------------------------------*/
