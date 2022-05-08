/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

#include <iostream>

#include <cooperative_groups.h>
#if (CUDART_VERSION >= 11000)
#include <cooperative_groups/reduce.h>
#endif

namespace cg = cooperative_groups;

/*---------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

#include "cs_blas_cuda.h"
#include "cs_math.h"
#include "cs_timer.h"

#include "cs_base_accel.h"
#include "cs_base_cuda.h"
#include "cs_cuda_contrib.h"

#include "cs_system_info.h"

// Number of threads per block. Should be high enough threads to keep block
// be busy, but not too high, otherwise the likelihood of getting more blocks
// scheduled for the same SM will be decreased.
// Must be a multiple of 32.

#define BLOCKSIZE 256

/*----------------------------------------------------------------------------*/

/* Size of arrays for tests */
static const int _n_sizes = 6;
static const cs_lnum_t _n_elts[]
  = {1000000, 750000, 500000, 300000, 200000, 100000};

static double _pi = 3.14159265358979323846;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for sum reduction within a block.
 *
 * \param[in, out]  n        number of values to reduce
 * \param[in, out]  g_idata  input array (size n)
 * \param[in, out]  g_odata  onput array (size 1)
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_reduce_single_block(size_t   n,
                     T       *g_idata,
                     T       *g_odata)
{
  __shared__ T sdata[blockSize];

  size_t tid = threadIdx.x;
  T r_s = 0;

  for (int i = threadIdx.x; i < n; i+= blockSize)
    r_s += g_idata[i];

  sdata[tid] = r_s;
  __syncthreads();

  for (int j = blockSize/2; j > CS_CUDA_WARP_SIZE; j /= 2) {
    if (tid < j) {
      sdata[tid] += sdata[tid + j];
    }
    __syncthreads();
  }

  if (tid < 32) cs_blas_cuda_warp_reduce_sum<blockSize>(sdata, tid);
  if (tid == 0) *g_odata = sdata[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.y, summing result over all threads of a block.
 *        The sum over blocks may be done by a separate kernel.
 *
 * blockSize must be a power of 2.
 *
 * \param[int]  n      array size
 * \param[in]   x      x vector
 * \param[in]   y      y vector
 * \param[out]  b_res  result of s = x.x for a block
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_xy_stage_1_of_2(cs_lnum_t    n,
                     const T     *x,
                     const T     *y,
                     double      *b_res)
{
  __shared__ double stmp[blockSize];

  int tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  stmp[tid] = 0;
  for (int i = blockIdx.x*(blockDim.x) + tid;
       i < n;
       i += grid_size) {
    stmp[tid] += static_cast<double>(x[i] * y[i]);
  }

  __syncthreads();

  // Loop might be unrolled since blockSize is a template parameter.

  for (int j = blockSize/2; j > CS_CUDA_WARP_SIZE; j /= 2) {
    if (tid < j) {
      stmp[tid] += stmp[tid + j];
    }
    __syncthreads();
  }

  if (tid < 32)
    cs_blas_cuda_warp_reduce_sum<blockSize>(stmp, tid);

  // Output: b_res for this block

  if (tid == 0) b_res[blockIdx.x] = stmp[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.y, summing result over all threads of a block.
 *        The sum over blocks may be done by a separate kernel.
 *
 * \param[in]   n    array size
 * \param[in]   x    x vector
 * \param[in]   y    y vector
 * \param[out]  res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_xy_stage_1_of_2_wr(cs_lnum_t    n,
                        const T     *x,
                        const T     *y,
                        double      *b_res)
{
  __shared__ double stmp[blockSize/CS_CUDA_WARP_SIZE + 1];

  cs_lnum_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;
  unsigned int mask_length = (blockSize & 31);  // 31 = warp_size-1
  mask_length = (mask_length > 0) ? (32 - mask_length) : mask_length;
  const unsigned int mask = (0xffffffff) >> mask_length;

  cs_real_t t_sum = 0;
  for (int i = blockIdx.x*blockDim.x + tid;
       i < n;
       i += grid_size) {
    t_sum += static_cast<double>(x[i] * y[i]);
  }

  // Reduce within warp using shuffle
  // or reduce_add if T==int & CUDA_ARCH == SM 8.0
  t_sum = warpReduceSum<double>(mask, t_sum);

  // Each thread puts its local sum into shared memory
  if ((tid % CS_CUDA_WARP_SIZE) == 0) {
    stmp[tid / CS_CUDA_WARP_SIZE] = t_sum;
  }

  __syncthreads();

  const unsigned int shmem_extent
    = (blockSize / CS_CUDA_WARP_SIZE) > 0 ? (blockSize / CS_CUDA_WARP_SIZE) : 1;
  const unsigned int ballot_result = __ballot_sync(mask, tid < shmem_extent);
  if (tid < shmem_extent) {
    t_sum = stmp[tid];
    // Reduce final warp using shuffle
    // or reduce_add if T==int & CUDA_ARCH == SM 8.0
    t_sum = warpReduceSum<cs_real_t>(ballot_result, t_sum);
  }

  // Write result for this block to global mem
  if (tid == 0) {
    b_res[blockIdx.x] = t_sum;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize result for atomic add-based algorithm.
 *
 * Can be run on a single block.
 *
 * \param[out]  res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_dot_xy_block_aa_pre(double  *res)
{
  if (blockIdx.x + threadIdx.x < 1)
    res[0] = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.y on device using atomic addition algorithm.
 *
 * \param[in]   n    array size
 * \param[in]   x    x vector
 * \param[in]   x    y vector
 * \param[out]  res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_xy_block_aa(cs_lnum_t    n,
                 T           *x,
                 T           *y,
                 double      *res)
{
  /* We use 1 shared-memory element per warp. */
  __shared__ double stmp[blockSize / CS_CUDA_WARP_SIZE];

  /*  We expect a 1D block. */
  cs_lnum_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  double t_xy = 0;
  for (int i = blockIdx.x*blockDim.x + tid;
       i < n;
       i += grid_size) {
    t_xy += static_cast<double>(x[i] * y[i]);
  }

  /* Reduce across the warp */

  const unsigned int full_mask = 0xffffffff;
  for (int offset = CS_CUDA_WARP_SIZE/2; offset > 0; offset /= 2)
    t_xy += __shfl_down_sync(full_mask, t_xy, offset);

  const unsigned n_warps = blockSize / CS_CUDA_WARP_SIZE;
  const unsigned lane_id = threadIdx.x % CS_CUDA_WARP_SIZE;
  const unsigned warp_id = threadIdx.x / CS_CUDA_WARP_SIZE;

  // First lane of each warp records its contribution to shared memory.
  if (lane_id == 0) {
    stmp[warp_id] = t_xy;
  }

  __syncthreads();

  // We only need the first warp from now on.
  if (warp_id)
    return;

  // Each thread in the warp picks one element from shared memory.
  t_xy = 0.0;
  if (lane_id < n_warps) {
    t_xy = stmp[lane_id];
  }

  // Make a reduction across the warp, but now we only need to cover n_warps.
  for (cs_lnum_t kk = 1; kk < n_warps; kk *= 2) {
    t_xy += __shfl_down_sync(0xffffffff, t_xy, kk);
  }

  // Write results atomically to memory.
  if (!lane_id) {
#if (__CUDA_ARCH__ < 600)
    atomicAddDouble(res, t_xy);
#else
    atomicAdd(res, t_xy);
#endif
  }
}

#if (CUDART_VERSION >= 11000)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Parallel sum reduction using shared memory.
 *
 * This version adds multiple elements per thread sequentially. This reduces
 * the overall cost of the algorithm while keeping the work complexity O(n) and
 * the step complexity O(log n). (Brent's Theorem optimization)
 *
 * takes log(n) steps for n input elements
 * uses n/2 threads
 * only works for power-of-2 arrays
 *
 * Note, this kernel needs a minimum of:
 * sizeof(double) * ((blockSize/32) + 1); bytes of shared memory.
 *
 * \param[in]       cta    cooperative thread block
 * \param[in, out]  sdata  shared data
 */
/*----------------------------------------------------------------------------*/

__device__ static void
_reduce_block_d(const cg::thread_block  &cta,
                double                  *sdata)
{
  const unsigned int tid = cta.thread_rank();
  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);

  sdata[tid] = cg::reduce(tile32, sdata[tid], cg::plus<double>());
  cg::sync(cta);

  double w = 0.0;
  if (cta.thread_rank() == 0) {
    w = 0;
    for (int i = 0; i < blockDim.x; i += tile32.size()) {
      w += sdata[i];
    }
    sdata[0] = w;
  }
  cg::sync(cta);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Kernel for dot product x.y on device using cooperative groups in
 * a single kernel invocation.
 *
 * For more details on the reduction algorithm (notably the multi-pass
 * approach), see the "reduction" sample in the CUDA SDK.
 *
 * This version adds multiple elements per thread sequentially.  This reduces
 * the overall cost of the algorithm while keeping the work complexity O(n) and
 * the step complexity O(log n). (Brent's Theorem optimization).
 *
 * Note, this kernel needs a minimum of:
 * sizeof(double) * ((blockSize/32) + 1); bytes of shared memory.
 *
 * \param[in]   n    array size
 * \param[in]   x    x vector
 * \param[in]   x    y vector
 * \param[out]  res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <class T>
__global__ static void
_dot_xy_cg_single_pass(cs_lnum_t   n,
                       const T    *x,
                       const T    *y,
                       double     *res)
{
  cg::thread_block block = cg::this_thread_block();
  cg::grid_group grid = cg::this_grid();

  extern double __shared__ sdata[];

  // Stride over grid and add the values to a shared memory buffer
  sdata[block.thread_rank()] = 0;

  for (int i = grid.thread_rank(); i < n; i += grid.size()) {
    sdata[block.thread_rank()] += static_cast<double>(x[i] * y[i]);
  }

  cg::sync(block);

  // Reduce each block (called once per block)
  _reduce_block_d(block, sdata);

  // Write out the result to global memory
  if (block.thread_rank() == 0) {
    res[blockIdx.x] = sdata[0];
  }
  cg::sync(grid);

  if (grid.thread_rank() == 0) {
    for (int block = 1; block < gridDim.x; block++) {
      res[0] += res[block];
    }
  }
}

#endif /* (CUDART_VERSION >= 11000) */

/*----------------------------------------------------------------------------
 * Count number of operations.
 *
 * parameters:
 *   n_runs       <-- Local number of runs
 *   n_ops        <-- Local number of operations
 *   n_ops_single <-- Single-processor equivalent number of operations
 *                    (without ghosts); ignored if 0
 *   wt           <-- wall-clock time
 *----------------------------------------------------------------------------*/

static void
_print_stats(long    n_runs,
             long    n_ops,
             long    n_ops_single,
             double  wt)
{
  double fm = 1.0 * n_runs / (1.e9 * CS_MAX(wt, 1e-12));

  if (cs_glob_n_ranks == 1)
    bft_printf("  N ops:       %12ld\n"
               "  Wall clock:  %12.5e\n"
               "  GFLOPS:      %12.5e\n",
               n_ops, wt/n_runs, n_ops*fm);

#if defined(HAVE_MPI)

  else {

    long n_ops_min, n_ops_max, n_ops_tot;
    double loc_count[2], glob_sum[2], glob_min[2], glob_max[2], fmg;

    loc_count[0] = wt;
    loc_count[1] = n_ops*fm;

    MPI_Allreduce(&n_ops, &n_ops_min, 1, MPI_LONG, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_max, 1, MPI_LONG, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_tot, 1, MPI_LONG, MPI_SUM,
                  cs_glob_mpi_comm);

    MPI_Allreduce(loc_count, glob_min, 2, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_max, 2, MPI_DOUBLE, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_sum, 2, MPI_DOUBLE, MPI_SUM,
                  cs_glob_mpi_comm);

    /* global flops multiplier */
    fmg = n_runs / (1.e9 * CS_MAX(glob_max[0], 1));

    glob_sum[0] /= n_runs;
    glob_min[0] /= n_runs;
    glob_max[0] /= n_runs;

    if (n_ops_single == 0)
      bft_printf
        ("               Mean         Min          Max          Total\n"
         "  N ops:       %12ld %12ld %12ld %12ld\n"
         "  Wall clock:  %12.5e %12.5e %12.5e\n"
         "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e\n",
         n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
         glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
         glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

    else
      bft_printf
        ("               Mean         Min          Max          Total"
         "        Single\n"
         "  N ops:       %12ld %12ld %12ld %12ld %12ld\n"
         "  Wall clock:  %12.5e %12.5e %12.5e\n"
         "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e %12.5e\n",
         n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
         n_ops_single,
         glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
         glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1],
         n_ops_tot*fmg, n_ops_single*fmg);
  }

#endif

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Simple dot product.
 *
 * parameters:
 *   t_measure <-- minimum time for each measure (< 0 for single pass)
 *   global    <-- 0 for local use, 1 for MPI sum
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_dot_products_1(double   t_measure,
                int      global)
{
  long n_ops;

  double test_sum = 0.0;

  int _global = global;

  cs_real_t *x = NULL;
  cs_real_t *y = NULL;
  cs_real_t *r = NULL, *r_grid = NULL;

  if (cs_glob_n_ranks == 1)
    _global = 0;

  /* Local dot product */

  const char *dot_name[] = {
    "2-stage",
    "2-stage, forced unroll",
    "2-stage, warpReduceSum",
    "1-stage, atomic block addition",
    "cublas",
    "cooperative groups"};

  for (int sub_id = 0; sub_id < _n_sizes; sub_id++) {

    cs_lnum_t n = _n_elts[sub_id];
    n_ops = n*2 - 1;

    double ref_s = 165 * (n/10) * _pi;

    const unsigned int block_size = 256;
    unsigned int grid_size = (n % block_size) ?
      n/block_size : n/block_size + 1;

    std::cout << std::endl << "block_size = " << block_size
              << ", grid_size = " << grid_size << std::endl;

    CS_MALLOC_HD(x, n, double, CS_ALLOC_HOST_DEVICE_PINNED);
    CS_MALLOC_HD(y, n, double, CS_ALLOC_HOST_DEVICE_PINNED);
    CS_MALLOC_HD(r, 1, double, CS_ALLOC_HOST_DEVICE_SHARED);

    CS_MALLOC_HD(r_grid, grid_size, double, CS_ALLOC_DEVICE);

    cudaDeviceSynchronize();

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 - 3)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);

    const cs_real_t *x_d = (const cs_real_t *)cs_get_device_ptr_const(x);
    const cs_real_t *y_d = (const cs_real_t *)cs_get_device_ptr_const(y);
    cs_real_t *r_d = (cs_real_t *)cs_get_device_ptr(r);

#if (CUDART_VERSION >= 11000)
    cudaDeviceProp prop = {0};
    cudaGetDeviceProperties(&prop, cs_glob_cuda_device_id);
    int n_max_blocks
      = prop.multiProcessorCount * (  prop.maxThreadsPerMultiProcessor
                                    / prop.maxThreadsPerBlock);

    int n_threads_cg = 1, n_blocks_cg = 1;
    if (n > 1) {
      cudaOccupancyMaxPotentialBlockSize(&n_blocks_cg, &n_threads_cg,
                                         _dot_xy_cg_single_pass<double>);
      n_blocks_cg = min(n_blocks_cg, n_max_blocks);
    }

    int smem_size_cg = n_threads_cg * sizeof(double);
#endif

    /* Events for timing */
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    for (int j = 0; j < 6; j++) {

      bool v_enabled = true;
      int run_id = 0;
      int n_runs = (t_measure > 0) ? 8 : 1;
      double s1 = 0;

      double wt0 = cs_timer_wtime(), wt1 = wt0;

      while (run_id < n_runs && v_enabled) {
        double test_sum_mult = 1.0/n_runs;
        while (run_id < n_runs) {
          switch(j) {
          case 0:
            _dot_xy_stage_1_of_2<block_size><<<grid_size, block_size, 0>>>
              (n, x_d, y_d, r_grid);
            _reduce_single_block<block_size><<<1, block_size, 0>>>
              (grid_size, r_grid, r_d);
            cudaStreamSynchronize(0);
            break;
          case 1:
            cs_blas_cuda_dot(n, x_d, y_d);
            break;
          case 2:
            _dot_xy_stage_1_of_2_wr<block_size><<<grid_size, block_size, 0>>>
              (n, x_d, y_d, r_grid);
            _reduce_single_block<block_size><<<1, block_size, 0>>>
              (grid_size, r_grid, r_d);
            cudaStreamSynchronize(0);
            break;
          case 3:
            _dot_xy_block_aa_pre<<<1, block_size>>>(r_d);
            cudaStreamSynchronize(0);
            _dot_xy_block_aa<block_size><<<grid_size, block_size>>>
              (n, x_d, y_d, r_d);
            cudaStreamSynchronize(0);
            break;
          case 4:
#if defined(HAVE_CUBLAS)
            cs_blas_cublas_dot(n, x_d, y_d);
#else
            run_id = n_runs;
            v_enabled = false;
#endif
            break;
          case 5:
#if (CUDART_VERSION >= 11000)
            {
              dim3 dim_grid(n_blocks_cg, 1, 1);
              dim3 dim_block(n_threads_cg, 1, 1);

              void *kernel_args[] = {
                (void *)&n,  (void *)&x_d, (void *)&y_d, (void *)&r_d,
              };

              cudaEventRecord(start, 0);
              cudaLaunchCooperativeKernel((void *)_dot_xy_cg_single_pass<double>,
                                          dim_grid, dim_block, kernel_args,
                                          smem_size_cg, NULL);
              cudaEventRecord(stop, 0);
              cudaStreamSynchronize(0);

              float k_time = 0;
              cudaEventElapsedTime(&k_time, start, stop);
            }
#else
            run_id = n_runs;
            v_enabled = false;
#endif
            break;
          }

          if (run_id == 0) {
            cs_sync_d2h(r);
            s1 = r[0];
          }
          else
            s1 = r[0];

#if defined(HAVE_MPI)
          if (_global && v_enabled) {
            double s1_glob = 0.0;
            MPI_Allreduce(&s1, &s1_glob, 1, MPI_DOUBLE, MPI_SUM,
                          cs_glob_mpi_comm);
            test_sum += test_sum_mult*s1_glob;
          }
#endif
          test_sum += test_sum_mult*s1;
          run_id++;
        }
        wt1 = cs_timer_wtime();
        if (wt1 - wt0 < t_measure)
          n_runs *= 2;
      }

      if (v_enabled == false)
        continue;

      if (_global == 0)
        bft_printf("\n"
                   "CUDA dot product X.Y %s (%d elts., %g err)\n"
                   "----------------\n",
                   dot_name[j], (int)n, fabs(s1 -ref_s));
      else
        bft_printf("\n"
                   "CUDA + MPI dot product X.Y %s (%d elts/rank)\n"
                   "----------------------\n",
                   dot_name[j], (int)n);

      bft_printf("  (calls: %d)\n", n_runs);

      _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    }

    cudaDeviceSynchronize();

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cs_blas_cuda_finalize();

    CS_FREE_HD(y);
    CS_FREE_HD(x);
    CS_FREE_HD(r);

    CS_FREE_HD(r_grid);

  }

  return test_sum;
}

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/

void
main_cuda(void)
{
  double t_measure = 10.0;
  double test_sum = 0.0;

  /* Initialization and environment */

  _pi = 4 * atan(1);

  bft_printf("\n");

  cs_base_cuda_select_default_device();

  /* Performance tests */
  /*-------------------*/

  if (cs_glob_n_ranks > 1)
    test_sum += _dot_products_1(t_measure, 1);

  test_sum += _dot_products_1(t_measure, 0);

  if (isnan(test_sum))
    bft_printf("test_sum is NaN\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
