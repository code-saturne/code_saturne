/*============================================================================
 * Unit tests for atomic and dispatch atomic operations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include <cmath>
#include <cstdlib>
#include <climits>
#include <cassert>
#include <iostream>
#include <algorithm>

/*----------------------------------------------------------------------------
 * Headers for current tests
 *----------------------------------------------------------------------------*/

#include "base/cs_base_accel.h"
#include "base/cs_math.h"
#include "base/cs_atomic.h"
#include "base/cs_dispatch.h"

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*============================================================================
 * Test 1: Atomic compare and swap (cs::atomic::compare_exchange)
 *============================================================================*/

static void
test_atomic_compare_exchange(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs::atomic::compare_exchange..." << std::endl;

  // 1. Sequential sanity check
  {
    cs_lnum_t val = 42;
    bool res1 = cs::atomic::compare_exchange(&val, (cs_lnum_t)10, (cs_lnum_t)99);
    assert(res1 == false);
    assert(val == 42);

    bool res2 = cs::atomic::compare_exchange(&val, (cs_lnum_t)42, (cs_lnum_t)99);
    assert(res2 == true);
    assert(val == 99);
  }

  // 2. Multi-threaded race condition check:
  // Exactly one thread out of N should succeed in updating from 0 to its id+1.
  {
    const cs_lnum_t n_threads = 500;
    cs_alloc_mode_t amode = ctx.alloc_mode();

    cs_lnum_t *target = nullptr;
    cs_lnum_t *success_count = nullptr;
    CS_MALLOC_HD(target, 1, cs_lnum_t, amode);
    CS_MALLOC_HD(success_count, 1, cs_lnum_t, amode);

    *target = 0;
    *success_count = 0;

    ctx.parallel_for(n_threads, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      bool won = cs::atomic::compare_exchange(target, (cs_lnum_t)0, i + 1);
      if (won) {
        cs::atomic::fetch_add(success_count, (cs_lnum_t)1);
      }
    });
    ctx.wait();

    assert(*success_count == 1);
    assert(*target >= 1 && *target <= n_threads);

    CS_FREE(target);
    CS_FREE(success_count);
  }
}

/*============================================================================
 * Test 2: Atomic fetch and add (cs::atomic::fetch_add)
 *============================================================================*/

static void
test_atomic_fetch_add(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs::atomic::fetch_add..." << std::endl;

  // 1. Sequential check
  {
    cs_lnum_t val = 10;
    cs_lnum_t prev = cs::atomic::fetch_add(&val, (cs_lnum_t)5);
    assert(prev == 10);
    assert(val == 15);
  }

  // 2. Multi-threaded concurrent accumulation into shared bins
  {
    const cs_lnum_t n_items = 2000;
    const cs_lnum_t n_bins = 8;
    cs_alloc_mode_t amode = ctx.alloc_mode();

    cs_lnum_t *bins = nullptr;
    CS_MALLOC_HD(bins, n_bins, cs_lnum_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++)
      bins[b] = 0;

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs::atomic::fetch_add(&bins[bin_id], (cs_lnum_t)1);
    });
    ctx.wait();

    cs_lnum_t total = 0;
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(bins[b] == n_items / n_bins);
      total += bins[b];
    }
    assert(total == n_items);

    CS_FREE(bins);
  }
}

/*============================================================================
 * Test 3: Atomic fetch and minimum (cs::atomic::fetch_min)
 *============================================================================*/

static void
test_atomic_fetch_min(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs::atomic::fetch_min..." << std::endl;

  // 1. Sequential check for integral and floating-point types
  {
    int i_val = 100;
    int i_prev = cs::atomic::fetch_min(&i_val, 150);
    assert(i_prev == 100);
    assert(i_val == 100);
    i_prev = cs::atomic::fetch_min(&i_val, 42);
    assert(i_prev == 100);
    assert(i_val == 42);

    double d_val = 100.5;
    double d_prev = cs::atomic::fetch_min(&d_val, 120.0);
    assert(std::fabs(d_prev - 100.5) < 1e-12);
    assert(std::fabs(d_val - 100.5) < 1e-12);
    d_prev = cs::atomic::fetch_min(&d_val, 42.25);
    assert(std::fabs(d_prev - 100.5) < 1e-12);
    assert(std::fabs(d_val - 42.25) < 1e-12);
    d_prev = cs::atomic::fetch_min(&d_val, -15.75);
    assert(std::fabs(d_prev - 42.25) < 1e-12);
    assert(std::fabs(d_val - (-15.75)) < 1e-12);

    float f_val = 50.0f;
    float f_prev = cs::atomic::fetch_min(&f_val, 80.0f);
    assert(std::fabs(f_prev - 50.0f) < 1e-6f);
    assert(std::fabs(f_val - 50.0f) < 1e-6f);
    f_prev = cs::atomic::fetch_min(&f_val, 12.5f);
    assert(std::fabs(f_prev - 50.0f) < 1e-6f);
    assert(std::fabs(f_val - 12.5f) < 1e-6f);
  }

  // 2. Multi-threaded race condition for cs_real_t (double)
  {
    const cs_lnum_t n_items = 2000;
    const cs_lnum_t n_bins = 10;
    cs_alloc_mode_t amode = ctx.alloc_mode();

    cs_real_t *bins = nullptr;
    CS_MALLOC_HD(bins, n_bins, cs_real_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++)
      bins[b] = 1.0e9;

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      // Deterministic non-trivial value with negative entries
      cs_real_t v = (cs_real_t)((i * 37) % 500) - 200.0;
      cs::atomic::fetch_min(&bins[bin_id], v);
    });
    ctx.wait();

    // Verify against sequential reference
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      cs_real_t ref_min = 1.0e9;
      for (cs_lnum_t i = b; i < n_items; i += n_bins) {
        cs_real_t v = (cs_real_t)((i * 37) % 500) - 200.0;
        if (v < ref_min)
          ref_min = v;
      }
      assert(std::fabs(bins[b] - ref_min) < 1e-9);
    }

    CS_FREE(bins);
  }
}


/*============================================================================
 * Test 3b: Atomic fetch and maximum (cs::atomic::fetch_max)
 *============================================================================*/

static void
test_atomic_fetch_max(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs::atomic::fetch_max..." << std::endl;

  // 1. Sequential check for integral and floating-point types
  {
    int i_val = 100;
    int i_prev = cs::atomic::fetch_max(&i_val, 50);
    assert(i_prev == 100);
    assert(i_val == 100);
    i_prev = cs::atomic::fetch_max(&i_val, 142);
    assert(i_prev == 100);
    assert(i_val == 142);

    double d_val = 100.5;
    double d_prev = cs::atomic::fetch_max(&d_val, 80.0);
    assert(std::fabs(d_prev - 100.5) < 1e-12);
    assert(std::fabs(d_val - 100.5) < 1e-12);
    d_prev = cs::atomic::fetch_max(&d_val, 142.25);
    assert(std::fabs(d_prev - 100.5) < 1e-12);
    assert(std::fabs(d_val - 142.25) < 1e-12);
    d_prev = cs::atomic::fetch_max(&d_val, 200.75);
    assert(std::fabs(d_prev - 142.25) < 1e-12);
    assert(std::fabs(d_val - 200.75) < 1e-12);

    float f_val = 50.0f;
    float f_prev = cs::atomic::fetch_max(&f_val, 20.0f);
    assert(std::fabs(f_prev - 50.0f) < 1e-6f);
    assert(std::fabs(f_val - 50.0f) < 1e-6f);
    f_prev = cs::atomic::fetch_max(&f_val, 82.5f);
    assert(std::fabs(f_prev - 50.0f) < 1e-6f);
    assert(std::fabs(f_val - 82.5f) < 1e-6f);
  }

  // 2. Multi-threaded race condition for cs_real_t (double)
  {
    const cs_lnum_t n_items = 2000;
    const cs_lnum_t n_bins = 10;
    cs_alloc_mode_t amode = ctx.alloc_mode();

    cs_real_t *bins = nullptr;
    CS_MALLOC_HD(bins, n_bins, cs_real_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++)
      bins[b] = -1.0e9;

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      // Deterministic non-trivial value with negative entries
      cs_real_t v = (cs_real_t)((i * 37) % 500) - 200.0;
      cs::atomic::fetch_max(&bins[bin_id], v);
    });
    ctx.wait();

    // Verify against sequential reference
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      cs_real_t ref_max = -1.0e9;
      for (cs_lnum_t i = b; i < n_items; i += n_bins) {
        cs_real_t v = (cs_real_t)((i * 37) % 500) - 200.0;
        if (v > ref_max)
          ref_max = v;
      }
      assert(std::fabs(bins[b] - ref_max) < 1e-9);
    }

    CS_FREE(bins);
  }
}

/*============================================================================
 * Test 4: cs_dispatch_sum (scalar and strided / vector)
 *============================================================================*/

static void
test_cs_dispatch_sum(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs_dispatch_sum..." << std::endl;

  const cs_lnum_t n_items = 2400;
  const cs_lnum_t n_bins = 12;
  cs_alloc_mode_t amode = ctx.alloc_mode();

  // 1. Scalar cs_dispatch_sum with CS_DISPATCH_SUM_ATOMIC
  {
    cs_real_t *bins = nullptr;
    CS_MALLOC_HD(bins, n_bins, cs_real_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++)
      bins[b] = 0.0;

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs_real_t inc = 1.25;
      cs_dispatch_sum(&bins[bin_id], inc, CS_DISPATCH_SUM_ATOMIC);
    });
    ctx.wait();

    cs_real_t expected = (cs_real_t)(n_items / n_bins) * 1.25;
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(bins[b] - expected) < 1e-9);
    }

    // Scalar cs_dispatch_sum with CS_DISPATCH_SUM_SIMPLE
    ctx.parallel_for(n_bins, [=] CS_F_HOST_DEVICE (cs_lnum_t b) {
      cs_dispatch_sum(&bins[b], (cs_real_t)10.0, CS_DISPATCH_SUM_SIMPLE);
    });
    ctx.wait();

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(bins[b] - (expected + 10.0)) < 1e-9);
    }

    CS_FREE(bins);
  }

  // 2. Strided cs_dispatch_sum<3> with CS_DISPATCH_SUM_ATOMIC
  {
    cs_real_3_t *vec_bins = nullptr;
    CS_MALLOC_HD(vec_bins, n_bins, cs_real_3_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      vec_bins[b][0] = 0.0;
      vec_bins[b][1] = 0.0;
      vec_bins[b][2] = 0.0;
    }

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs_real_t v[3] = {1.0, 2.0, 3.0};
      cs_dispatch_sum<3>(vec_bins[bin_id], v, CS_DISPATCH_SUM_ATOMIC);
    });
    ctx.wait();

    cs_real_t count = (cs_real_t)(n_items / n_bins);
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(vec_bins[b][0] - count * 1.0) < 1e-9);
      assert(std::fabs(vec_bins[b][1] - count * 2.0) < 1e-9);
      assert(std::fabs(vec_bins[b][2] - count * 3.0) < 1e-9);
    }

    // Strided cs_dispatch_sum<3> with CS_DISPATCH_SUM_SIMPLE
    ctx.parallel_for(n_bins, [=] CS_F_HOST_DEVICE (cs_lnum_t b) {
      cs_real_t add[3] = {10.0, 20.0, 30.0};
      cs_dispatch_sum<3>(vec_bins[b], add, CS_DISPATCH_SUM_SIMPLE);
    });
    ctx.wait();

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(vec_bins[b][0] - (count * 1.0 + 10.0)) < 1e-9);
      assert(std::fabs(vec_bins[b][1] - (count * 2.0 + 20.0)) < 1e-9);
      assert(std::fabs(vec_bins[b][2] - (count * 3.0 + 30.0)) < 1e-9);
    }

    CS_FREE(vec_bins);
  }
}

/*============================================================================
 * Test 5: cs_dispatch_min (scalar and strided / vector)
 *============================================================================*/

static void
test_cs_dispatch_min(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs_dispatch_min..." << std::endl;

  const cs_lnum_t n_items = 2400;
  const cs_lnum_t n_bins = 12;
  cs_alloc_mode_t amode = ctx.alloc_mode();

  // 1. Scalar cs_dispatch_min with CS_DISPATCH_MIN_ATOMIC
  {
    cs_real_t *bins = nullptr;
    CS_MALLOC_HD(bins, n_bins, cs_real_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++)
      bins[b] = 1.0e9;

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs_real_t val = (cs_real_t)((i * 41 + 17) % 733) - 350.0;
      cs_dispatch_min(&bins[bin_id], val, CS_DISPATCH_MIN_ATOMIC);
    });
    ctx.wait();

    // Verify against sequential reference
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      cs_real_t ref_min = 1.0e9;
      for (cs_lnum_t i = b; i < n_items; i += n_bins) {
        cs_real_t val = (cs_real_t)((i * 41 + 17) % 733) - 350.0;
        if (val < ref_min)
          ref_min = val;
      }
      assert(std::fabs(bins[b] - ref_min) < 1e-9);
    }

    // Scalar cs_dispatch_min with CS_DISPATCH_MIN_SIMPLE
    ctx.parallel_for(n_bins, [=] CS_F_HOST_DEVICE (cs_lnum_t b) {
      cs_dispatch_min(&bins[b], (cs_real_t)-999.0, CS_DISPATCH_MIN_SIMPLE);
    });
    ctx.wait();

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(bins[b] - (-999.0)) < 1e-9);
    }

    CS_FREE(bins);
  }

  // 2. Strided cs_dispatch_min<3> with CS_DISPATCH_MIN_ATOMIC
  {
    cs_real_3_t *vec_bins = nullptr;
    CS_MALLOC_HD(vec_bins, n_bins, cs_real_3_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      vec_bins[b][0] = 1.0e9;
      vec_bins[b][1] = 1.0e9;
      vec_bins[b][2] = 1.0e9;
    }

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs_real_t v[3] = {
        (cs_real_t)((i * 13 + 5) % 311) - 100.0,
        (cs_real_t)((i * 29 + 11) % 433) - 200.0,
        (cs_real_t)((i * 47 + 23) % 521) - 300.0
      };
      cs_dispatch_min<3>(vec_bins[bin_id], v, CS_DISPATCH_MIN_ATOMIC);
    });
    ctx.wait();

    // Verify against sequential reference
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      cs_real_t ref_min[3] = {1.0e9, 1.0e9, 1.0e9};
      for (cs_lnum_t i = b; i < n_items; i += n_bins) {
        cs_real_t v[3] = {
          (cs_real_t)((i * 13 + 5) % 311) - 100.0,
          (cs_real_t)((i * 29 + 11) % 433) - 200.0,
          (cs_real_t)((i * 47 + 23) % 521) - 300.0
        };
        if (v[0] < ref_min[0]) ref_min[0] = v[0];
        if (v[1] < ref_min[1]) ref_min[1] = v[1];
        if (v[2] < ref_min[2]) ref_min[2] = v[2];
      }
      assert(std::fabs(vec_bins[b][0] - ref_min[0]) < 1e-9);
      assert(std::fabs(vec_bins[b][1] - ref_min[1]) < 1e-9);
      assert(std::fabs(vec_bins[b][2] - ref_min[2]) < 1e-9);
    }

    // Strided cs_dispatch_min<3> with CS_DISPATCH_MIN_SIMPLE
    ctx.parallel_for(n_bins, [=] CS_F_HOST_DEVICE (cs_lnum_t b) {
      cs_real_t v[3] = {-1000.0, -2000.0, -3000.0};
      cs_dispatch_min<3>(vec_bins[b], v, CS_DISPATCH_MIN_SIMPLE);
    });
    ctx.wait();

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(vec_bins[b][0] - (-1000.0)) < 1e-9);
      assert(std::fabs(vec_bins[b][1] - (-2000.0)) < 1e-9);
      assert(std::fabs(vec_bins[b][2] - (-3000.0)) < 1e-9);
    }

    CS_FREE(vec_bins);
  }
}


/*============================================================================
 * Test 6: cs_dispatch_max (scalar and strided / vector)
 *============================================================================*/

static void
test_cs_dispatch_max(cs_dispatch_context &ctx)
{
  std::cout << "  Testing cs_dispatch_max..." << std::endl;

  const cs_lnum_t n_items = 2400;
  const cs_lnum_t n_bins = 12;
  cs_alloc_mode_t amode = ctx.alloc_mode();

  // 1. Scalar cs_dispatch_max with CS_DISPATCH_MAX_ATOMIC
  {
    cs_real_t *bins = nullptr;
    CS_MALLOC_HD(bins, n_bins, cs_real_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++)
      bins[b] = -1.0e9;

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs_real_t val = (cs_real_t)((i * 41 + 17) % 733) - 350.0;
      cs_dispatch_max(&bins[bin_id], val, CS_DISPATCH_MAX_ATOMIC);
    });
    ctx.wait();

    // Verify against sequential reference
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      cs_real_t ref_max = -1.0e9;
      for (cs_lnum_t i = b; i < n_items; i += n_bins) {
        cs_real_t val = (cs_real_t)((i * 41 + 17) % 733) - 350.0;
        if (val > ref_max)
          ref_max = val;
      }
      assert(std::fabs(bins[b] - ref_max) < 1e-9);
    }

    // Scalar cs_dispatch_max with CS_DISPATCH_MAX_SIMPLE
    ctx.parallel_for(n_bins, [=] CS_F_HOST_DEVICE (cs_lnum_t b) {
      cs_dispatch_max(&bins[b], (cs_real_t)999.0, CS_DISPATCH_MAX_SIMPLE);
    });
    ctx.wait();

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(bins[b] - 999.0) < 1e-9);
    }

    CS_FREE(bins);
  }

  // 2. Strided cs_dispatch_max<3> with CS_DISPATCH_MAX_ATOMIC
  {
    cs_real_3_t *vec_bins = nullptr;
    CS_MALLOC_HD(vec_bins, n_bins, cs_real_3_t, amode);

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      vec_bins[b][0] = -1.0e9;
      vec_bins[b][1] = -1.0e9;
      vec_bins[b][2] = -1.0e9;
    }

    ctx.parallel_for(n_items, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      cs_lnum_t bin_id = i % n_bins;
      cs_real_t v[3] = {
        (cs_real_t)((i * 13 + 5) % 311) - 100.0,
        (cs_real_t)((i * 29 + 11) % 433) - 200.0,
        (cs_real_t)((i * 47 + 23) % 521) - 300.0
      };
      cs_dispatch_max<3>(vec_bins[bin_id], v, CS_DISPATCH_MAX_ATOMIC);
    });
    ctx.wait();

    // Verify against sequential reference
    for (cs_lnum_t b = 0; b < n_bins; b++) {
      cs_real_t ref_max[3] = {-1.0e9, -1.0e9, -1.0e9};
      for (cs_lnum_t i = b; i < n_items; i += n_bins) {
        cs_real_t v[3] = {
          (cs_real_t)((i * 13 + 5) % 311) - 100.0,
          (cs_real_t)((i * 29 + 11) % 433) - 200.0,
          (cs_real_t)((i * 47 + 23) % 521) - 300.0
        };
        if (v[0] > ref_max[0]) ref_max[0] = v[0];
        if (v[1] > ref_max[1]) ref_max[1] = v[1];
        if (v[2] > ref_max[2]) ref_max[2] = v[2];
      }
      assert(std::fabs(vec_bins[b][0] - ref_max[0]) < 1e-9);
      assert(std::fabs(vec_bins[b][1] - ref_max[1]) < 1e-9);
      assert(std::fabs(vec_bins[b][2] - ref_max[2]) < 1e-9);
    }

    // Strided cs_dispatch_max<3> with CS_DISPATCH_MAX_SIMPLE
    ctx.parallel_for(n_bins, [=] CS_F_HOST_DEVICE (cs_lnum_t b) {
      cs_real_t v[3] = {1000.0, 2000.0, 3000.0};
      cs_dispatch_max<3>(vec_bins[b], v, CS_DISPATCH_MAX_SIMPLE);
    });
    ctx.wait();

    for (cs_lnum_t b = 0; b < n_bins; b++) {
      assert(std::fabs(vec_bins[b][0] - 1000.0) < 1e-9);
      assert(std::fabs(vec_bins[b][1] - 2000.0) < 1e-9);
      assert(std::fabs(vec_bins[b][2] - 3000.0) < 1e-9);
    }

    CS_FREE(vec_bins);
  }
}

/*============================================================================
 * Main runner
 *============================================================================*/

int
main(int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

#if defined(HAVE_OPENMP)
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

  std::cout << "Starting cs_dispatch_atomic_test..." << std::endl;

  cs_dispatch_context ctx;

  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      ctx.set_use_gpu(false);
      ctx.set_n_min_per_cpu_thread(5);
      std::cout << "Running tests in Host / CPU mode:" << std::endl;
    }
    else {
      std::cout << "Running tests in default context mode:" << std::endl;
    }

    test_atomic_compare_exchange(ctx);
    test_atomic_fetch_add(ctx);
    test_atomic_fetch_min(ctx);
    test_atomic_fetch_max(ctx);
    test_cs_dispatch_sum(ctx);
    test_cs_dispatch_min(ctx);
    test_cs_dispatch_max(ctx);
  }

  std::cout << "All atomic and dispatch atomic tests PASSED successfully!" << std::endl;

  return EXIT_SUCCESS;
}
