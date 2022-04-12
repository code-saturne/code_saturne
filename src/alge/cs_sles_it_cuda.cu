/*============================================================================
 * Low-level functions and global variables definition for CUDA.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_log.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>
#include <stdio.h>
#include <algorithm>
#include <list>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_benchmark_cuda.h"
#include "cs_matrix.h"
#include "cs_numbering.h"
#include "cs_matrix_priv.h"
#include <cublas_v2.h>
#include "cs_base_accel.h"
#include <cusparse.h>
#include <cuda_runtime_api.h>
#include "cs_sles_it.h"
#include "cs_sles_pc.h"
#include "cs_sles_it_priv.h"
#include "cs_sles_it_cuda.h"
#include "cs_halo.h"
#include "cs_parall.h"
#include "mpi-ext.h"
#if (CUDART_VERSION > 10000)
#include <cooperative_groups/reduce.h>
#endif
#include <cooperative_groups.h>
namespace cg = cooperative_groups;
/*----------------------------------------------------------------------------*/


/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

 /* SIMD unit size to ensure SIMD alignement (2 to 8 required on most
  * current architectures, so 16 should be enough on most architectures) */

#define CS_SIMD_SIZE(s) (((s-1)/32+1)*32)
#define BLOCKSIZE 256

# if defined(HAVE_MPI)
  static MPI_Comm _comm;
# endif

/* in the file cs_sles_it.c GPUdirect is used for choosing the variant*/
int GPUdirect = getenv("VARIANT_CUDA")==NULL?0:atoi(getenv ("VARIANT_CUDA"));

/* in the file cs_sles_default.c Methode is used for choosing  the solver iterative*/
int Methode = getenv("METHODE_SOLVER")==NULL?0:atoi(getenv ("METHODE_SOLVER"));


template<typename T>
__device__ static __forceinline__ T ldg(const T* ptr) {
#if __CUDA_ARCH__ >= 350
    return __ldg(ptr);
#else
    return *ptr;
#endif
}


/* fonction for reduction within a warp*/
 template <size_t blockSize, typename T>
__device__ static void __forceinline__ _warpReduce(volatile T *sdata, size_t tid)
{
    if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] += sdata[tid +  8];
    if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
    if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
    if (blockSize >=  2) sdata[tid] += sdata[tid +  1];
}

/* fonction for reduction within a block, reduce0 in report*/
template <size_t blockSize, typename T>
__global__ static void _reduceCUDA(T* g_idata, T* g_odata, size_t n)
{
    __shared__ T sdata[blockSize];

    size_t tid = threadIdx.x;
    size_t i = blockIdx.x*(blockSize) + tid;
    T mySum = 0;
    while (i < n) { mySum += g_idata[i]; i += blockSize; }

    sdata[tid] = mySum;
    __syncthreads();

    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }

    if (tid < 32) _warpReduce<blockSize>(sdata, tid);
    if (tid == 0) *g_odata = sdata[0];
}



#if (CUDART_VERSION > 10000)
/* do the reduction within a block*/
__device__ void _reduceBlock(double *sdata, const cg::thread_block &cta)
{
    const unsigned int tid = cta.thread_rank();
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);

    sdata[tid] = cg::reduce(tile32, sdata[tid], cg::plus<double>());
    cg::sync(cta);

    double beta = 0.0;
    if (cta.thread_rank() == 0) {
        beta  = 0;
        for (int i = 0; i < blockDim.x; i += tile32.size()) {
            beta  += sdata[i];
        }
        sdata[0] = beta;
    }
    cg::sync(cta);
}

/* fonction for reduction within a block, reduce2 in report*/
template <size_t blockSize, typename T>
__global__ static void _reduceCUDA2(T* g_idata, T* g_odata, size_t n)
{
    cg::thread_block block = cg::this_thread_block();
    __shared__ T sdata[blockSize];

    size_t tid = threadIdx.x;
    size_t i = blockIdx.x*(blockSize) + tid;
    sdata[tid] = 0;

    while (i < n) { sdata[tid] += g_idata[i]; i += blockSize; }
    __syncthreads();

    _reduceBlock(sdata, block);

    if (tid == 0)
      *g_odata = sdata[0];
}
#endif

/* fonction for reduction within a block, reduce1 in report*/
template <size_t blockSize, typename T>
__global__ static void _reduceCUDA1(T* g_idata, T* g_odata, size_t n)
{
    __shared__ T sdata[blockSize];

    size_t tid = threadIdx.x;
    size_t i = blockIdx.x*(blockSize) + tid;
    T mySum = 0;
    while (i < n) { mySum += g_idata[i]; i += blockSize; }
    sdata[tid] = mySum;
    __syncthreads();

    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; mySum = sdata[tid]; } __syncthreads(); }

    cg::thread_block cta = cg::this_thread_block();
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);
    if (cta.thread_rank() < 32) {
      // Fetch final intermediate sum from 2nd warp
      if (blockSize >= 64) mySum += sdata[tid + 32];
        // Reduce final warp using shuffle
      for (int offset = tile32.size() / 2; offset > 0; offset /= 2) {
        mySum += tile32.shfl_down(mySum, offset);
      }
    }

    if (tid == 0) *g_odata = mySum;
}

/*----------------------------------------------------------------------------
 * Local residue compute
 *
 * parameters:
 *   rhs          <--  right hand side
 *   vx           <--> system solution
 *   rk           <--  old system solution
 *   ad_inv       <--  inverse of diagonal
 *   ad           <--  diagonal
 *   res2         -->  residue
 *----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_compute_residue( const cs_real_t           *rhs,
                  cs_real_t                 *restrict vx,
                  cs_real_t                 *restrict rk,
                  const cs_real_t           *restrict ad_inv,
                  const cs_real_t           *restrict ad,
                  T                         *sum_block,
                  cs_lnum_t                  n_rows)
{
    __shared__ T sdata[blockSize];
    cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
    size_t tid = threadIdx.x;
    if (ii < n_rows){
      vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
      double r = ad[ii] * (vx[ii]-rk[ii]);
      sdata[tid] = r * r;
      rk[ii] = vx[ii];
    }else
      sdata[tid] = 0.0;
    __syncthreads();
    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    if (tid < 32) _warpReduce<blockSize>(sdata, tid);
    if (tid == 0) sum_block[blockIdx.x] = sdata[0];
}

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute forward and backward to solve an LU 3*3 system.
 *
 * parameters:
 *   mat   <-- 3*3*dim matrix
 *   x     --> solution
 *   b     --> 1st part of RHS (c - b)
 *   c     --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

__device__ static void  __forceinline__ _fw_and_bw_lu33_cuda(const cs_real_t       mat[],
                                                             cs_real_t        x[restrict],
                                                             cs_real_t        b[restrict],
                                                             const cs_real_t  c[restrict])
{
  cs_real_t  aux[3];

  aux[0] = (c[0] - b[0]);
  aux[1] = (c[1] - b[1]) - aux[0]*mat[3];
  aux[2] = (c[2] - b[2]) - aux[0]*mat[6] - aux[1]*mat[7];

  x[2] = aux[2]/mat[8];
  x[1] = (aux[1] - mat[5]*x[2])/mat[4];
  x[0] = (aux[0] - mat[1]*x[1] - mat[2]*x[2])/mat[0];
}


/*----------------------------------------------------------------------------
 * Local residue compute uesd by Block Jacobi
 *----------------------------------------------------------------------------*/
__global__ static void
_block_3_compute_residue( const cs_real_t           *rhs,
                          cs_real_t                 *restrict vx,
                          cs_real_t                 *restrict rk,
                          cs_real_t                 *restrict vxx,
                          const cs_real_t           *restrict ad_inv,
                          const cs_real_t           *restrict ad,
                          cs_real_t                 *sum_block,
                          cs_lnum_t                  n_blocks)
{
    __shared__ cs_real_t sdata[BLOCKSIZE];
    cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
    size_t tid = threadIdx.x;
    sdata[tid] = 0.0;
    if (ii < n_blocks){
      _fw_and_bw_lu33_cuda(ad_inv + 9*ii, vx + 3*ii, vxx + 3*ii, rhs + 3*ii);
      #pragma unroll
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        double r = 0.0;
        #pragma unroll
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          r +=  ad[ii*9 + jj*3 + kk] * (vx[ii*3 + kk] - rk[ii*3 + kk]);

        sdata[tid] += (r*r);
      }
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        rk[ii*3 + kk] = vx[ii*3 + kk];
    }
    __syncthreads();
    if (BLOCKSIZE >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (BLOCKSIZE >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (BLOCKSIZE >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (BLOCKSIZE >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    if (tid < 32) _warpReduce<BLOCKSIZE>(sdata, tid);
    if (tid == 0) sum_block[blockIdx.x] = sdata[0];
}


/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/
 __global__ static void _mat_vect_p_l_msr(bool                         exclude_diag,
                                         const cs_lnum_t       *restrict col_id_gpu,
                                         const cs_lnum_t    *restrict row_index_gpu,
                                         const cs_real_t        *restrict x_val_gpu,
                                         const cs_real_t        *restrict d_val_gpu,
                                         const cs_real_t                *restrict x,
                                         cs_real_t                      *restrict y,
                                         cs_lnum_t                           n_rows)
{
  cs_lnum_t ii;
  ii = blockIdx.x*blockDim.x + threadIdx.x;
  if( ii < n_rows){
    const cs_lnum_t *restrict col_id = col_id_gpu + row_index_gpu[ii];
    const cs_real_t *restrict m_row = x_val_gpu + row_index_gpu[ii];
    cs_lnum_t n_cols = row_index_gpu[ii+1] - row_index_gpu[ii];
    cs_real_t sii = 0.0;
    #pragma unroll
    for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*ldg(&x[col_id[jj]]));
    if(!exclude_diag && d_val_gpu != NULL)
      y[ii] = sii + d_val_gpu[ii]*x[ii];
    else
      y[ii] = sii;
  }
}


/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, 3x3 blocked version.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/
__global__ static void _b_3_3_mat_vect_p_l_msr(bool                         exclude_diag,
                                              const cs_lnum_t       *restrict col_id_gpu,
                                              const cs_lnum_t    *restrict row_index_gpu,
                                              const cs_real_t        *restrict x_val_gpu,
                                              const cs_real_t        *restrict d_val_gpu,
                                              const cs_real_t                *restrict x,
                                              cs_real_t                      *restrict y,
                                              cs_lnum_t                           n_rows)
{
  cs_lnum_t ii;
  ii = blockIdx.x*blockDim.x + threadIdx.x;
  if( ii < n_rows){
    const cs_lnum_t *restrict col_id = col_id_gpu + row_index_gpu[ii];
    const cs_real_t *restrict m_row = x_val_gpu + row_index_gpu[ii];
    cs_lnum_t n_cols = row_index_gpu[ii+1] - row_index_gpu[ii];
    cs_real_t sii[3];

    if (!exclude_diag && d_val_gpu != NULL){
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < 3; kk++)
         sii[kk]  =   d_val_gpu[ii*9+kk*3]*x[ii*3]
                    + d_val_gpu[ii*9+kk*3+1]*x[ii*3+1]
                    + d_val_gpu[ii*9+kk*3+2]*x[ii*3+2];
    }else{
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < 3; kk++)
         sii[kk] = 0.;
    }
    #pragma unroll(2)
    for (cs_lnum_t jj = 0; jj < n_cols; jj++){
      for (cs_lnum_t kk = 0; kk < 3; kk++)
         sii[kk] += (m_row[jj]*ldg(&x[col_id[jj]*3 + kk]));
    }

    y[ii*3] = sii[0];
    y[ii*3+1] = sii[1];
    y[ii*3+2] = sii[2];
  }
}

/* fonction for MPI */
#if defined(HAVE_MPI)
static void
CUDART_CB _residue_mpi(void*          res)
{

  double *res1 = (double *)res;
  double res2 = *res1;
  double _sum;
  MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, _comm);
  *res1 = _sum;
}
#endif

 /*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

#if (CUDART_VERSION > 9020)
cs_sles_convergence_state_t
jacobi_graph( cs_sles_it_t              *c,
             const cs_matrix_t          *a,
             cs_lnum_t                  diag_block_size,
             cs_sles_it_convergence_t  *convergence,
             const cs_real_t           *rhs,
             cs_real_t                 *restrict vx,
             size_t                     aux_size,
             void                      *aux_vectors)
{
    int DeviceId;
    cudaGetDevice(&DeviceId);
    cudaStream_t stream1;
    cudaStreamCreate(&stream1);
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    size_t size = n_cols*sizeof(cs_real_t);
    cudaMemPrefetchAsync(vx, size, DeviceId, stream);
    cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);
    cs_sles_convergence_state_t cvg;
    double  *res;
    double res2;
    double residue;
    const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
    const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
    const cs_lnum_t  *restrict col_id
      = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
    const cs_lnum_t  *restrict row_index
      = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
    const cs_real_t  *restrict val
      = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
    const cs_real_t  *restrict d_val
      = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));
    CS_MALLOC_HD(res, 1, double, CS_ALLOC_HOST_DEVICE_SHARED);
    cs_real_t *_aux_vectors;
    cs_real_t *restrict rk;
    unsigned n_iter = 0;

    /* Allocate or map work arrays */
    /*-----------------------------*/

    assert(c->setup_data != NULL);
    const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv_gpu;
    const cs_lnum_t n_rows = c->setup_data->n_rows;
    bool destroy_aux_vectors = false;

    {
      const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
      const size_t n_wa = 1;
      const size_t wa_size = CS_SIMD_SIZE(n_cols);

      if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
        destroy_aux_vectors = true;
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
        cudaMalloc(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#else
        cudaMallocManaged(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#endif
      }
      else
        _aux_vectors = (cs_real_t *)aux_vectors;

      rk = _aux_vectors;
    }


    cudaMemcpyAsync(rk, vx, n_rows* sizeof(cs_real_t), cudaMemcpyDeviceToDevice, stream);


    const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);
    cvg = CS_SLES_ITERATING;

    const unsigned int blocksize = 256;
    unsigned int gridsize = (unsigned int)ceil((double)n_rows/blocksize);
    unsigned int gridsize2 = 1;
    double  *sum_block;
    CS_MALLOC_HD(sum_block, gridsize, double, CS_ALLOC_DEVICE);

    /* Build graph */
    /*-------------------*/
    cudaGraph_t graph;
    cudaGraphExec_t execGraph = NULL;
    cudaGraphNode_t errorNode;
    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
    _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(true, col_id, row_index, val, d_val, rk, vx, n_rows);
    _compute_residue<blocksize><<<gridsize, blocksize, 0, stream>>>(rhs, vx, rk, ad_inv, ad, sum_block, n_rows);
    _reduceCUDA<blocksize><<<gridsize2, blocksize, 0, stream>>>(sum_block, res, gridsize);
    cudaStreamEndCapture(stream, &graph);
    cudaError_t status = cudaGraphInstantiate(&execGraph, graph, &errorNode, nullptr, 0);
    assert(status == cudaSuccess);

    /* Current iteration */
    /*-------------------*/
    while (cvg == CS_SLES_ITERATING) {
      n_iter += 1;
      if (cs_glob_n_ranks == 1 && a->halo != NULL)
        local_halo_cuda(a->halo, CS_HALO_STANDARD, rk, 1, stream);
      else if(a->halo != NULL)
        cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, rk, 1);
      cudaGraphLaunch(execGraph, stream);
      cudaStreamSynchronize(stream);
      res2 = *res;
      #if defined(HAVE_MPI)
          if (c->comm != MPI_COMM_NULL) {
            double _sum;
            MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
            res2 = _sum;
          }
      #endif /* defined(HAVE_MPI) */

      residue = sqrt(res2); /* Actually, residue of previous iteration */
      /* Convergence test */
      if (n_iter == 1)
            c->setup_data->initial_residue = residue;

      cvg = convergence_test_cuda(c, n_iter, residue, convergence);
    }
    cudaGraphDestroy(graph);
    cudaGraphExecDestroy(execGraph);
    if (destroy_aux_vectors)
      cudaFree(_aux_vectors);

    CS_FREE_HD(res);
    CS_FREE_HD(sum_block);
    cudaStreamDestroy(stream);
    cudaStreamDestroy(stream1);
    return cvg;

}
#endif

cs_sles_convergence_state_t
jacobi_cuda( cs_sles_it_t              *c,
             const cs_matrix_t         *a,
             cs_lnum_t                  diag_block_size,
             cs_sles_it_convergence_t  *convergence,
             const cs_real_t           *rhs,
             cs_real_t                 *restrict vx,
             size_t                     aux_size,
             void                      *aux_vectors)
{
    int DeviceId;
    cudaGetDevice(&DeviceId);
    cudaStream_t stream1;
    cudaStreamCreate(&stream1);
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    size_t size = n_cols*sizeof(cs_real_t);
    cudaMemPrefetchAsync(vx, size, DeviceId, stream);
    cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);
    cs_sles_convergence_state_t cvg;
    double  residue;
    double  res2;
    double  *res;
    const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
    const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
    const cs_lnum_t  *restrict col_id
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
    const cs_lnum_t  *restrict row_index
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
    const cs_real_t  *restrict val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
    const cs_real_t  *restrict d_val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));
    CS_MALLOC_HD(res, 1, double, CS_ALLOC_HOST_DEVICE_SHARED);
    cs_real_t *_aux_vectors;
    cs_real_t *restrict rk;
    unsigned n_iter = 0;

    /* Allocate or map work arrays */
    /*-----------------------------*/

    assert(c->setup_data != NULL);

    const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv_gpu;
    const cs_lnum_t n_rows = c->setup_data->n_rows;
    bool destroy_aux_vectors = false;


    {
      const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
      const size_t n_wa = 1;
      const size_t wa_size = CS_SIMD_SIZE(n_cols);

      if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
        destroy_aux_vectors = true;
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
        cudaMalloc(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#else
        cudaMallocManaged(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#endif
      }
      else
        _aux_vectors = (cs_real_t *)aux_vectors;

      rk = _aux_vectors;
    }


    cudaMemcpyAsync(rk, vx, n_rows* sizeof(cs_real_t), cudaMemcpyDeviceToDevice, stream);


    const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);

    cvg = CS_SLES_ITERATING;

    const unsigned int blocksize = 256;
    unsigned int gridsize = (unsigned int)ceil((double)n_rows/blocksize);
    unsigned int gridsize2 = 1;
    double  *sum_block;
    CS_MALLOC_HD(sum_block, gridsize, double, CS_ALLOC_DEVICE);
    /* Current iteration */
    /*-------------------*/

    while (cvg == CS_SLES_ITERATING) {

      n_iter += 1;

      /* Compute Vx <- Vx - (A-diag).Rk and residue. */
      if (cs_glob_n_ranks == 1 && a->halo != NULL)
        local_halo_cuda(a->halo, CS_HALO_STANDARD, rk, 1, stream);
      else if(a->halo != NULL)
        cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, rk, 1);

      _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(true, col_id, row_index, val, d_val, rk, vx, n_rows);

      _compute_residue<blocksize><<<gridsize, blocksize, 0, stream>>>(rhs, vx, rk, ad_inv, ad, sum_block, n_rows);

      _reduceCUDA<blocksize><<<gridsize2, blocksize, 0, stream>>>(sum_block, res, gridsize);

      cudaStreamSynchronize(stream);
      res2 = *res;

  #if defined(HAVE_MPI)

      if (c->comm != MPI_COMM_NULL) {
        double _sum;
        MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
        res2 = _sum;
      }

  #endif /* defined(HAVE_MPI) */

      residue = sqrt(res2); /* Actually, residue of previous iteration */

      /* Convergence test */

      if (n_iter == 1)
        c->setup_data->initial_residue = residue;

      cvg = convergence_test_cuda(c, n_iter, residue, convergence);
    }

    if (destroy_aux_vectors)
      cudaFree(_aux_vectors);

    CS_FREE_HD(res);
    CS_FREE_HD(sum_block);
    cudaStreamDestroy(stream);
    cudaStreamDestroy(stream1);
    return cvg;

}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- block size of diagonal elements
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              --> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/
 #if (CUDART_VERSION > 9020)
 cs_sles_convergence_state_t
 block_3_jacobi_graph(cs_sles_it_t             *c,
                     const cs_matrix_t         *a,
                     cs_lnum_t                  diag_block_size,
                     cs_sles_it_convergence_t  *convergence,
                     const cs_real_t           *rhs,
                     cs_real_t                 *restrict vx,
                     size_t                     aux_size,
                     void                      *aux_vectors)
{
  int DeviceId;
  cudaGetDevice(&DeviceId);
  cudaStream_t stream1;
  cudaStreamCreate(&stream1);
  cudaStream_t stream;
  cudaStreamCreate(&stream);
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
  size_t size = n_cols*sizeof(cs_real_t);
  cudaMemPrefetchAsync(vx, size, DeviceId, stream);
  cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);

  assert(diag_block_size == 3);
  const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
  const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
  const cs_lnum_t  *restrict col_id
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
  const cs_lnum_t  *restrict row_index
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
  const cs_real_t  *restrict val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
  const cs_real_t  *restrict d_val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));
  cs_sles_convergence_state_t cvg;
  double  res2, residue;
  double  *res;
  CS_MALLOC_HD(res, 1, double, CS_ALLOC_HOST_DEVICE_SHARED);
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vxx;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv_gpu;

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = c->setup_data->n_rows / 3;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
       cudaMalloc(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#else
       cudaMallocManaged(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#endif
    }else
      _aux_vectors = (cs_real_t *)aux_vectors;

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  cudaMemcpyAsync(rk, vx, n_rows* sizeof(cs_real_t), cudaMemcpyDeviceToDevice, stream);
  const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);

  cvg = CS_SLES_ITERATING;

  const unsigned int blocksize = BLOCKSIZE;
  unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
  unsigned int gridsize2 = 1;
  double  *sum_block;
  CS_MALLOC_HD(sum_block, gridsize, double, CS_ALLOC_DEVICE);
  /* Build graph */
  /*-------------------*/
  cudaGraph_t graph;
  cudaGraphExec_t execGraph = NULL;
  cudaGraphNode_t errorNode;
  cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
  _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(true, col_id, row_index, val, d_val, rk, vxx, n_blocks);
  _block_3_compute_residue<<<gridsize, blocksize, 0, stream>>>(rhs, vx, rk, vxx, ad_inv, ad, sum_block, n_blocks);
  _reduceCUDA<blocksize><<<gridsize2, blocksize, 0, stream>>>(sum_block, res, gridsize);
  //cudaLaunchHostFunc(stream, _residue_plus_mpi, c, res, n_iter, convergence, cvg);
  cudaStreamEndCapture(stream, &graph);
  cudaError_t status = cudaGraphInstantiate(&execGraph, graph, &errorNode, nullptr, 0);
  assert(status == cudaSuccess);
  while (cvg == CS_SLES_ITERATING) {
    n_iter += 1;
    if (cs_glob_n_ranks == 1 && a->halo != NULL)
       local_halo_cuda(a->halo, CS_HALO_STANDARD, rk, 3, stream);
    else if(a->halo != NULL)
       cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, rk, 3);
    cudaGraphLaunch(execGraph, stream);
    cudaStreamSynchronize(stream);
    res2 = *res;
    #if defined(HAVE_MPI)
        if (c->comm != MPI_COMM_NULL) {
          double _sum;
          MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
          res2 = _sum;
        }
    #endif /* defined(HAVE_MPI) */

    residue = sqrt(res2); /* Actually, residue of previous iteration */
    /* Convergence test */
    if (n_iter == 1)
          c->setup_data->initial_residue = residue;

    cvg = convergence_test_cuda(c, n_iter, residue, convergence);
  }
  cudaGraphDestroy(graph);
  cudaGraphExecDestroy(execGraph);
  if (_aux_vectors != aux_vectors)
    cudaFree(_aux_vectors);

  CS_FREE_HD(res);
  CS_FREE_HD(sum_block);
  cudaStreamDestroy(stream);
  cudaStreamDestroy(stream1);
  return cvg;
}
#endif

 cs_sles_convergence_state_t
 block_3_jacobi_cuda(cs_sles_it_t              *c,
                     const cs_matrix_t         *a,
                     cs_lnum_t                  diag_block_size,
                     cs_sles_it_convergence_t  *convergence,
                     const cs_real_t           *rhs,
                     cs_real_t                 *restrict vx,
                     size_t                     aux_size,
                     void                      *aux_vectors)
{
  int DeviceId;
  cudaGetDevice(&DeviceId);
  cudaStream_t stream1;
  cudaStreamCreate(&stream1);
  cudaStream_t stream;
  cudaStreamCreate(&stream);
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
  size_t size = n_cols*sizeof(cs_real_t);
  cudaMemPrefetchAsync(vx, size, DeviceId, stream);
  cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);

  assert(diag_block_size == 3);
  cs_sles_convergence_state_t cvg;
  const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
  const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
  const cs_lnum_t  *restrict col_id
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
  const cs_lnum_t  *restrict row_index
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
  const cs_real_t  *restrict val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
  const cs_real_t  *restrict d_val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));
  double  res2, residue;
  double  *res;
  CS_MALLOC_HD(res, 1, double, CS_ALLOC_HOST_DEVICE_SHARED);

  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vxx;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv_gpu;

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = c->setup_data->n_rows / 3;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      cudaMalloc(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#else
      cudaMallocManaged(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#endif
    }else
      _aux_vectors = (cs_real_t *)aux_vectors;

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  cudaMemcpyAsync(rk, vx, n_rows* sizeof(cs_real_t), cudaMemcpyDeviceToDevice, stream);
  const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);

  cvg = CS_SLES_ITERATING;

  const unsigned int blocksize = BLOCKSIZE;
  unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
  unsigned int gridsize2 = 1;
  double  *sum_block;
  CS_MALLOC_HD(sum_block, gridsize, double, CS_ALLOC_DEVICE);
  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* Compute vxx <- vx - (a-diag).rk and residue. */
    if (cs_glob_n_ranks == 1 && a->halo != NULL)
      local_halo_cuda(a->halo, CS_HALO_STANDARD, rk, 3, stream);
    else if(a->halo != NULL)
      cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, rk, 3);

    _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(true, col_id, row_index, val, d_val, rk, vxx, n_blocks);
    /* Compute vx <- diag^-1 . (vxx - rhs) and residue. */
    _block_3_compute_residue<<<gridsize, blocksize, 0, stream>>>(rhs, vx, rk, vxx, ad_inv, ad, sum_block, n_blocks);

    _reduceCUDA<blocksize><<<gridsize2, blocksize, 0, stream>>>(sum_block, res, gridsize);

    cudaStreamSynchronize(stream);

    res2 = *res;


#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residue = sqrt(res2); /* Actually, residue of previous iteration */

    if (n_iter == 1)
      c->setup_data->initial_residue = residue;

    /* Convergence test */

    cvg = convergence_test_cuda(c, n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    cudaFree(_aux_vectors);

  CS_FREE_HD(res);
  CS_FREE_HD(sum_block);
  cudaStreamDestroy(stream);
  cudaStreamDestroy(stream1);
  return cvg;
}

/*----------------------------------------------------------------------------
 * Block Gauss-Seidel utilities.
 * Compute forward and backward to solve an LU P*P system.
 *
 * parameters:
 *   mat     <-- P*P*dim matrix
 *   db_size <-- matrix size
 *   x       --> solution
 *   b       <-> RHS in, work array
 *----------------------------------------------------------------------------*/

__device__ void __forceinline__ _fw_and_bw_lu_gs_cuda(const cs_real_t mat[], int db_size,
                                      cs_real_t x[], const cs_real_t b[])
{

  /* forward */
  #pragma unroll
  for (int ii = 0; ii < db_size; ii++) {
    x[ii] = b[ii];
    for (int jj = 0; jj < ii; jj++)
      x[ii] -= x[jj] * mat[ii * db_size + jj];
  }

  /* backward */
  for (int ii = db_size - 1; ii >= 0; ii--) {
    for (int jj = db_size - 1; jj > ii; jj--)
      x[ii] -= x[jj] * mat[ii * db_size + jj];
    x[ii] /= mat[ii * (db_size + 1)];
  }
}

/*----------------------------------------------------------------------------
 * Gauss-Seidel utilities.
 * Compute Vx <- Vx - (A-diag).Rk and residue: forward step .
 *----------------------------------------------------------------------------*/

__global__ static void
_mv_seidel_forward_cuda(const cs_lnum_t   *a_col_id,
                        const cs_real_t   *a_x_val,
                        const cs_lnum_t   *a_row_index,
                        const cs_real_t   *rhs,
                        cs_real_t         *vx,
                        const cs_real_t   *ad_inv,
                        const cs_lnum_t    n_rows)
{
  cs_lnum_t ii;
  ii = blockIdx.x*blockDim.x + threadIdx.x;
  if( ii < n_rows){
    const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
    const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
    const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
    cs_real_t vx0 = rhs[ii];
    #pragma unroll
    for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*ldg(&vx[col_id[jj]]));
    vx[ii] = vx0 * ad_inv[ii];
  }
}

/*----------------------------------------------------------------------------
 * Gauss-Seidel utilities.
 * Compute Vx <- Vx - (A-diag).Rk and residue: backward step .
 *----------------------------------------------------------------------------*/

__global__ static void
_mv_seidel_backward_cuda(const cs_lnum_t  *a_col_id,
                        const cs_real_t   *a_x_val,
                        const cs_lnum_t   *a_row_index,
                        const cs_real_t   *rhs,
                        const cs_real_t   *ad,
                        cs_real_t         *vx,
                        const cs_real_t   *ad_inv,
                        const cs_lnum_t    n_rows,
                        cs_real_t         *sum_block)
{
  cs_lnum_t ii;
  ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  __shared__ cs_real_t sdata[BLOCKSIZE];
  sdata[tid] = 0.0;
  if( ii < n_rows){

    const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
    const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
    const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

    cs_real_t vxm1 = vx[ii];
    cs_real_t vx0 = rhs[ii];

    #pragma unroll
    for (cs_lnum_t jj = 0; jj < n_cols; jj++)
      vx0 -= (m_row[jj]*ldg(&vx[col_id[jj]]));

    vx0 *= ad_inv[ii];
    double r = ad[ii] * (vx0-vxm1);
    sdata[tid] = (r*r);

    vx[ii] = vx0;
  }
  __syncthreads();
  if (BLOCKSIZE >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (BLOCKSIZE >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (BLOCKSIZE >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (BLOCKSIZE >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
  if (tid < 32) _warpReduce<BLOCKSIZE>(sdata, tid);
  if (tid == 0) sum_block[blockIdx.x] = sdata[0];
}

/*----------------------------------------------------------------------------
 * Block Gauss-Seidel utilities.
 * Compute Vx <- Vx - (A-diag).Rk and residue: forward step .
 *----------------------------------------------------------------------------*/
__global__ static void
_mv_seidel_block_forward_cuda(const cs_lnum_t   *a_col_id,
                              const cs_real_t   *a_x_val,
                              const cs_lnum_t   *a_row_index,
                              const cs_real_t   *rhs,
                              cs_real_t         *vx,
                              const cs_real_t   *ad_inv,
                              const cs_lnum_t    n_rows,
                              cs_lnum_t          diag_block_size,
                              cs_lnum_t          *db_size)
{

  cs_lnum_t ii;
  ii = blockIdx.x*blockDim.x + threadIdx.x;
  if( ii < n_rows){
    const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
    const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
    const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

    cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];
    #pragma unroll
    for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
        vx0[kk] = rhs[ii*db_size[1] + kk];
    #pragma unroll(2)
    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
       for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
           vx0[kk] -= (m_row[jj]*ldg(&vx[col_id[jj]*db_size[1] + kk]));
    }

    _fw_and_bw_lu_gs_cuda(ad_inv + db_size[3]*ii,
                         db_size[0],
                         _vx,
                         vx0);
    #pragma unroll
    for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
       vx[ii*db_size[1] + kk] = _vx[kk];
   }

}

/*----------------------------------------------------------------------------
 * Block Gauss-Seidel utilities.
 * Compute Vx <- Vx - (A-diag).Rk and residue: backward step .
 *----------------------------------------------------------------------------*/
__global__ static void
_mv_seidel_block_backward_cuda(const cs_lnum_t  *a_col_id,
                              const cs_real_t   *a_x_val,
                              const cs_lnum_t   *a_row_index,
                              const cs_real_t   *rhs,
                              cs_real_t         *vx,
                              const cs_real_t   *ad_inv,
                              const cs_real_t   *ad,
                              const cs_lnum_t    n_rows,
                              cs_lnum_t          diag_block_size,
                              cs_lnum_t          *db_size,
                              cs_real_t          *sum_block)
{
  cs_lnum_t ii;
  ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  __shared__ cs_real_t sdata[BLOCKSIZE];
  sdata[tid] = 0.0;
  if( ii < n_rows){
    const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
    const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
    const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

    cs_real_t vx0[DB_SIZE_MAX], vxm1[DB_SIZE_MAX], _vx[DB_SIZE_MAX];
    #pragma unroll
    for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
      vxm1[kk] = vx[ii*db_size[1] + kk];
      vx0[kk] = rhs[ii*db_size[1] + kk];
    }
    #pragma unroll(2)
    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      for (cs_lnum_t kk = 0; kk < db_size[0]; kk++)
        vx0[kk] -= (m_row[jj]*ldg(&vx[col_id[jj]*db_size[1] + kk]));
    }

    _fw_and_bw_lu_gs_cuda(ad_inv + db_size[3]*ii,
                         db_size[0],
                         _vx,
                         vx0);

    double rr = 0;
    #pragma unroll
    for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
       double r = ad[ii*db_size[1] + kk] * (_vx[kk]-vxm1[kk]);
       rr += (r*r);
       vx[ii*db_size[1] + kk] = _vx[kk];
    }
    sdata[tid] = rr;
  }
  __syncthreads();
  if (BLOCKSIZE >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (BLOCKSIZE >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (BLOCKSIZE >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (BLOCKSIZE >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
  if (tid < 32) _warpReduce<BLOCKSIZE>(sdata, tid);
  if (tid == 0) sum_block[blockIdx.x] = sdata[0];

}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local symmetric Gauss-Seidel.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (unused here)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

 cs_sles_convergence_state_t
 p_sym_gauss_seidel_msr_cuda(cs_sles_it_t              *c,
                             const cs_matrix_t         *a,
                             cs_lnum_t                  diag_block_size,
                             cs_sles_it_convergence_t  *convergence,
                             const cs_real_t           *rhs,
                             cs_real_t                 *restrict vx,
                             size_t                     aux_size,
                             void                      *aux_vectors)
 {
   CS_UNUSED(aux_size);
   CS_UNUSED(aux_vectors);
   const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
   int DeviceId;
   cudaGetDevice(&DeviceId);
   cudaStream_t stream;
   cudaStreamCreate(&stream);
   size_t sz = n_rows*sizeof(cs_real_t);
   cudaMemPrefetchAsync(vx, sz, DeviceId, stream);
   cudaMemPrefetchAsync(rhs, sz, DeviceId, stream);
   cs_sles_convergence_state_t cvg;
   const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
   const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
   double  res2, residue;

   double  *res;
   CS_MALLOC_HD(res, 1, double, CS_ALLOC_HOST_DEVICE_SHARED);

   /* Check matrix storage type */

   if (cs_matrix_get_type(a) != CS_MATRIX_MSR)
     bft_error
       (__FILE__, __LINE__, 0,
        _("Symmetric Gauss-Seidel Jacobi hybrid solver only supported with a\n"
          "matrix using %s (%s) storage."),
        cs_matrix_type_name[CS_MATRIX_MSR],
        _(cs_matrix_type_fullname[CS_MATRIX_MSR]));

   unsigned n_iter = 0;

   const cs_halo_t *halo = cs_matrix_get_halo(a);

   const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv_gpu;

   const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);

   const cs_lnum_t  *a_col_id
    = (const cs_lnum_t  *)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
   const cs_lnum_t  *a_row_index
    = (const cs_lnum_t  *)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
   const cs_real_t  *a_x_val
    = (const cs_real_t  *)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
   const cs_real_t  *a_d_val
    = (const cs_real_t  *)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));

   const cs_lnum_t *db_size = cs_matrix_get_diag_block_size(a);
   cs_lnum_t *size = NULL;
   CS_MALLOC_HD(size, 4, cs_lnum_t, CS_ALLOC_DEVICE);
   cudaMemcpy(size, db_size, 4*sizeof(cs_lnum_t), cudaMemcpyHostToDevice);

   cvg = CS_SLES_ITERATING;

   const unsigned int blocksize = BLOCKSIZE;
   unsigned int gridsize = (unsigned int)ceil((double)n_rows/blocksize);
   unsigned int gridsize2 = 1;
   double  *sum_block;
   CS_MALLOC_HD(sum_block, gridsize, double, CS_ALLOC_DEVICE);

   /* Current iteration */
   /*-------------------*/

   while (cvg == CS_SLES_ITERATING) {

     n_iter += 1;

     /* Synchronize ghost cells first */

     if (halo != NULL){
      if (a->db_size[3] == 1)
         cs_halo_sync_var_cuda_aware(halo, CS_HALO_STANDARD, vx, 1);
      /* Blocked version */
      else
         cs_halo_sync_var_cuda_aware(halo, CS_HALO_STANDARD, vx,
                                             a->db_size[1]);
     }

     /* Compute Vx <- Vx - (A-diag).Rk and residue: forward step */

     if (diag_block_size == 1)
       _mv_seidel_forward_cuda<<<gridsize, blocksize, 0, stream>>>(a_col_id, a_x_val, a_row_index,
                                rhs, vx, ad_inv, n_rows);
     else
       _mv_seidel_block_forward_cuda<<<gridsize, blocksize, 0, stream>>>(a_col_id, a_x_val, a_row_index,
      rhs, vx, ad_inv, n_rows, diag_block_size, size);

     /* Synchronize ghost cells again */

     if (halo != NULL){
       if (a->db_size[3] == 1)
          cs_halo_sync_var_cuda_aware(halo, CS_HALO_STANDARD, vx, 1);
       else
          cs_halo_sync_var_cuda_aware(halo, CS_HALO_STANDARD, vx,
                                             a->db_size[1]);
     }


     /* Compute Vx <- Vx - (A-diag).Rk and residue: backward step */

     if (diag_block_size == 1)
        _mv_seidel_backward_cuda<<<gridsize, blocksize, 0, stream>>>(a_col_id, a_x_val, a_row_index,
        rhs, ad, vx, ad_inv, n_rows, sum_block);
     else
        _mv_seidel_block_backward_cuda<<<gridsize, blocksize, 0, stream>>>(a_col_id, a_x_val, a_row_index,
         rhs, vx, ad_inv, ad, n_rows, diag_block_size, size, sum_block);

     _reduceCUDA<blocksize><<<gridsize2, blocksize, 0, stream>>>(sum_block, res, gridsize);
     cudaStreamSynchronize(stream);
     res2 = *res;
     if (convergence->precision > 0. || c->plot != NULL) {

 #if defined(HAVE_MPI)

       if (c->comm != MPI_COMM_NULL) {
         double _sum;
         MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         res2 = _sum;
       }

 #endif /* defined(HAVE_MPI) */

       residue = sqrt(res2); /* Actually, residue of previous iteration */

       /* Convergence test */

       if (n_iter == 1)
         c->setup_data->initial_residue = residue;

       cvg = convergence_test_cuda(c, n_iter, residue, convergence);

     }
     else if (n_iter >= convergence->n_iterations_max) {
       convergence->n_iterations = n_iter;
       cvg = CS_SLES_MAX_ITERATION;
     }

   }

   CS_FREE_HD(res);
   CS_FREE_HD(sum_block);
   CS_FREE_HD(size);
   cudaStreamDestroy(stream);
   return cvg;
 }
/*----------------------------------------------------------------------------
 * Function for application of a Jacobi preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

 __global__ static void
 _sles_pc_poly_apply_jacobi(const cs_lnum_t      n_rows,
                            const cs_real_t      *ad_inv,
                            cs_real_t            *x_in,
                            cs_real_t            *x_out)
 {
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
   size_t gridsize = blockDim.x*gridDim.x;
   while(ii < n_rows){
    x_out[ii] = x_in[ii] * ad_inv[ii];
    ii += gridsize;
  }
 }

 /*----------------------------------------------------------------------------
  * Compute dot product x.x, summing result over all threads of a block.
  *
  * parameters:
  *   x      <-- vector in s = x.x
  *
  * returns:
  *   result of s = x.x
  *----------------------------------------------------------------------------*/
  __global__ static void dot_xx(cs_lnum_t  size,
                                cs_real_t  *x,
                                cs_real_t  *sum)
 {
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
   cs_lnum_t tid = threadIdx.x;
   __shared__ cs_real_t sdata[BLOCKSIZE];
   size_t gridsize = blockDim.x*gridDim.x;
   cs_real_t mySum = 0;
   while(ii < size){
     mySum += x[ii]*x[ii];
     ii += gridsize;
   }
   sdata[threadIdx.x] = mySum;
   __syncthreads();
   if (BLOCKSIZE >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
   if (BLOCKSIZE >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
   if (BLOCKSIZE >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
   if (BLOCKSIZE >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
   if (tid < 32) _warpReduce<BLOCKSIZE>(sdata, tid);
   if (tid == 0) sum[blockIdx.x] = sdata[0];
 }

  /*----------------------------------------------------------------------------
  * Compute dot product x.y, summing result over all threads of a block.
  *
  * parameters:
  *   x      <-- vector in s = x.y
  *
  * returns:
  *   result of s = x.y
  *----------------------------------------------------------------------------*/
  __global__ static void dot_xy(cs_lnum_t  size,
                         cs_real_t  *x,
                         cs_real_t  *y,
                         cs_real_t  *sum)
 {
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
   cs_lnum_t tid = threadIdx.x;
   __shared__ cs_real_t sdata[BLOCKSIZE];
   size_t gridsize = blockDim.x*gridDim.x;
   cs_real_t mySum = 0;
   while(ii < size){
     mySum += x[ii]*y[ii];
     ii += gridsize;
   }
   sdata[threadIdx.x] = mySum;
   __syncthreads();
   if (BLOCKSIZE >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
   if (BLOCKSIZE >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
   if (BLOCKSIZE >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
   if (BLOCKSIZE >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
   if (tid < 32) _warpReduce<BLOCKSIZE>(sdata, tid);
   if (tid == 0) sum[blockIdx.x] = sdata[0];
 }

 /*----------------------------------------------------------------------------
  * Compute x = x - y
  *----------------------------------------------------------------------------*/
  __global__ static void _x_minus_y (cs_lnum_t     size,
                                    cs_real_t        *x,
                                    const cs_real_t  *y)
 {
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
   size_t gridsize = blockDim.x*gridDim.x;
   while(ii < size){
      x[ii] -= y[ii];
      ii += gridsize;
   }
 }

  /*----------------------------------------------------------------------------
  * Compute x = x - alpha * y
  *----------------------------------------------------------------------------*/

  __global__ static void _x_minus_ay (cs_lnum_t  size,
                                      cs_real_t  *x,
                                      cs_real_t  *y,
                                      cs_real_t  *a)
 {
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
   double alpha = *a;
   size_t gridsize = blockDim.x*gridDim.x;
   while(ii < size){
     x[ii] -= alpha * y[ii];
     ii += gridsize;
   }
 }

  /*----------------------------------------------------------------------------
  * Compute x = x / alpha
  *----------------------------------------------------------------------------*/
 __global__ static void _x_divise_a (cs_lnum_t  size,
                                     cs_real_t  *x,
                                     cs_real_t  *a)
 {
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
   size_t gridsize = blockDim.x*gridDim.x;
   cs_real_t alpha = sqrt(*a);
   while(ii < size){
    x[ii] /= alpha;
    ii += gridsize;
  }
 }


  /*----------------------------------------------------------------------------
  * Compute solution of gcr
  *----------------------------------------------------------------------------*/
__global__ static void _compute_vx(cs_lnum_t     n_rows,
                                   cs_lnum_t     n_iter,
                                   cs_real_t     *vx,
                                   cs_real_t     *alpha,
                                   cs_real_t     *zk,
                                   size_t         wa_size,
                                   cs_real_t     *gkj_inv)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if(ii < n_rows){
    double sii = 0.0;
    #pragma unroll(2)
    for(cs_lnum_t kk = 0; kk < n_iter; kk++){
      for(cs_lnum_t jj = 0; jj <= kk; jj++){
        const cs_real_t *zk_j = zk + jj*wa_size;
        sii += alpha[kk] * zk_j[ii] * gkj_inv[(kk + 1) * kk / 2 + jj];
      }
    }
    vx[ii] -= sii;
  }
}


__global__ static void _compute_inv(cs_real_t *gkj_inv, cs_real_t *gkj)
{
  cs_lnum_t kk = threadIdx.x;

  cs_real_t alpha = sqrt(gkj[(kk + 1) * kk / 2 + kk]);

  for(cs_lnum_t ii = 0; ii < kk; ii++) {
    for (cs_lnum_t jj = 0; jj < kk; jj++)
      gkj_inv[(kk + 1) * kk / 2 + ii]
        +=   ((ii <= jj) ? gkj_inv[(jj + 1) * jj / 2 + ii] : 0.0)
           * gkj[(kk + 1) * kk / 2  + jj];
  }

  for (cs_lnum_t jj = 0; jj < kk; jj++)
    gkj_inv[(kk + 1) * kk / 2 + jj] /= - gkj[(kk + 1) * kk / 2 + kk];

  gkj_inv[(kk + 1) * kk / 2 + kk] = 1.0 / alpha;

}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using optimised preconditioned GCR.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size (unused here)
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/
 #if (CUDART_VERSION > 9020)
 cs_sles_convergence_state_t
 gcr_graph(cs_sles_it_t             *c,
          const cs_matrix_t         *a,
          cs_lnum_t                  diag_block_size,
          cs_sles_it_convergence_t  *convergence,
          const cs_real_t           *rhs,
          cs_real_t                 *restrict vx,
          size_t                     aux_size,
          void                      *aux_vectors)
{
  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = n_rows / 3;
  int DeviceId;
  cudaGetDevice(&DeviceId);
  cudaStream_t stream, stream1;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream);
  size_t size = n_rows*sizeof(cs_real_t);
  cudaMemPrefetchAsync(vx, size, DeviceId, stream);
  cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);
  const cs_real_t  *ad_inv = cs_sles_pc_get_ad_inv(c->setup_data->pc_context);
  cudaMemPrefetchAsync(ad_inv, size, DeviceId, stream1);
  const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
  const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
  const cs_lnum_t  *restrict col_id
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
  const cs_lnum_t  *restrict row_index
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
  const cs_real_t  *restrict val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
  const cs_real_t  *restrict d_val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));

  cs_sles_convergence_state_t cvg;
  double  residue;
  double  *res, *_recv_buf;
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  cudaMalloc(&res, sizeof(cs_real_t));
  cudaMalloc(&_recv_buf, sizeof(cs_real_t));
#else
  cudaMallocManaged(&res, sizeof(cs_real_t));
  cudaMallocManaged(&_recv_buf, sizeof(cs_real_t));
#endif

  cs_real_t *_aux_vectors, *alpha;
  cs_real_t *restrict rk, *restrict zk, *restrict ck;
  cs_real_t *restrict gkj, *restrict gkj_inv;

  /* In case of the standard GCR, n_k_per_restart --> Inf,
   * or stops until convergence*/
  unsigned n_iter = 0;
  const unsigned n_k_per_restart = c->restart_interval;
  size_t wa_size;
  unsigned n_restart = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1 + n_k_per_restart * 2;
    wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
      CS_MALLOC_HD(_aux_vectors, wa_size, cs_real_t, CS_ALLOC_DEVICE);
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      CS_MALLOC_HD(zk, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_DEVICE);
#else
      CS_MALLOC_HD(zk, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);
#endif
      CS_MALLOC_HD(ck, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_DEVICE);
    }
    else{
      _aux_vectors = (cs_real_t *)aux_vectors;                            /* store residuals  */
      zk = _aux_vectors + wa_size;                     /* store inv(M)*r   */
      ck = _aux_vectors + wa_size * (1 + n_k_per_restart);   /* store A*zk */
    }
    rk = _aux_vectors;
  }

  /* gkj stores the upper triangle matrix Gamma of crossed inner-products
   * Not necessary to allocate the full matrix size
   * gkj_inv stores the inverse of gkj */
  cudaMallocManaged(&alpha, n_k_per_restart*sizeof(cs_real_t));
  cudaMallocManaged(&gkj, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));

  cudaMalloc(&gkj_inv, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));

  cvg = CS_SLES_ITERATING;

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);
  cs_real_t *zk_n;
  cs_real_t *ck_n;

  unsigned int gridsize_dot;
  /* 80 * 2048 is the number of threads used for product vector */
  if(n_rows > 80*2048){
    gridsize_dot =  (unsigned int)ceil((double)80*2048/BLOCKSIZE);
  }else
    gridsize_dot = (unsigned int)ceil((double)n_rows/BLOCKSIZE);
  const unsigned int blocksize = BLOCKSIZE;
  const unsigned int blocksize2 = 640;
  unsigned int gridsize = (unsigned int)ceil((double)n_rows/blocksize);
  unsigned int gridsize2 = 1;
  double  *sum_block;
  cudaMalloc(&sum_block, gridsize*sizeof(double));
  bool commMPI = false;
#if defined(HAVE_MPI)
  if (c->comm != MPI_COMM_NULL){
       _comm = c->comm;
       commMPI = true;
  }
#endif /* defined(HAVE_MPI) */

  cudaGraph_t graph;
  cudaGraphExec_t graphExec = NULL;
  cudaGraphNode_t errorNode;
  cudaGraphExecUpdateResult updateResult_out;


  /* Current Restart */
  while (cvg == CS_SLES_ITERATING) {

    n_iter = 0;

    /* Initialize iterative calculation */
    /*----------------------------------*/
    if (a->db_size[3] == 1){
      cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, vx, 1);
      _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, vx, rk, n_rows);
    }else{
      cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, vx, a->db_size[1]);
      unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
      _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, vx, rk, n_blocks);
    }
    _x_minus_y<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, rhs);
    dot_xx<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, sum_block);
    _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, res, gridsize_dot);
    cudaMemcpyAsync(&residue, res, sizeof(double), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

#if defined(HAVE_MPI)
    if (c->comm != MPI_COMM_NULL){
      MPI_Allreduce(res, _recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      cudaMemcpy(&residue, _recv_buf, sizeof(double), cudaMemcpyDeviceToHost);
    }
#endif

    residue = sqrt(residue);
    //printf("residue,: %f\n", residue);
    if (n_restart == 0)
      c->setup_data->initial_residue = residue;

    /* Current Iteration on k */
    /* ---------------------- */

    while (cvg == CS_SLES_ITERATING && n_iter < n_k_per_restart) {

      /* Preconditionning */
      zk_n = zk + n_iter * wa_size;
      ck_n = ck + n_iter * wa_size;
      _sles_pc_poly_apply_jacobi<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ad_inv, rk, zk_n);
      cudaStreamSynchronize(stream);
      /* Halo and Matrix.vector product y = A.x with MSR matrix*/
      if (a->db_size[3] == 1){
         cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, zk_n, 1);
         _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, zk_n, ck_n, n_rows);
      }else{
         cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, zk_n, a->db_size[1]);
         unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
         _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, zk_n, ck_n, n_blocks);
      }

      /* compute the direction ck_n*/
      for(cs_lnum_t jj = 0; jj < n_iter; jj++) {
        cs_real_t *ck_j = ck + jj * wa_size;
        dot_xy<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_j, ck_n, sum_block);
        _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, &gkj[(n_iter + 1) * n_iter / 2 + jj], gridsize_dot);
        if(commMPI)
          cudaLaunchHostFunc(stream, _residue_mpi, (void *)&gkj[(n_iter+1) * n_iter / 2 + jj]);
        _x_minus_ay<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, ck_j, &gkj[(n_iter + 1) * n_iter / 2 + jj]);
      }

     /* begin to capture graph*/
      CS_CUDA_CHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal));
      dot_xx<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, sum_block);
      _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, &gkj[(n_iter+1) * n_iter / 2 + n_iter], gridsize_dot);
      if (commMPI)
        cudaLaunchHostFunc(stream, _residue_mpi, (void *)&gkj[(n_iter+1) * n_iter / 2 + n_iter]);
      _x_divise_a<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, &gkj[(n_iter+1) * n_iter / 2 + n_iter]);
      dot_xy<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, rk, sum_block);
      _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, &alpha[n_iter], gridsize_dot);
      if (commMPI)
        cudaLaunchHostFunc(stream, _residue_mpi, (void *)&alpha[n_iter]);
      _x_minus_ay<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, ck_n, &alpha[n_iter]);
      dot_xx<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, sum_block);
      _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, res, gridsize_dot);
      cudaMemcpyAsync(&residue, res, sizeof(double), cudaMemcpyDeviceToHost, stream);
      cudaStreamEndCapture(stream, &graph);
      /* end to capture graph*/

      /* update and launch graph*/
      if(graphExec == NULL){
        CS_CUDA_CHECK(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));
      }
      else{
        CS_CUDA_CHECK(cudaGraphExecUpdate(graphExec, graph, NULL, &updateResult_out));
        if (updateResult_out != cudaGraphExecUpdateSuccess){
          if (graphExec != NULL)
              CS_CUDA_CHECK(cudaGraphExecDestroy(graphExec));
          printf("k = %d graph update failed with error - %d\n", n_iter, updateResult_out);
          CS_CUDA_CHECK(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));
        }
      }
      CS_CUDA_CHECK(cudaGraphLaunch(graphExec, stream));
      CS_CUDA_CHECK(cudaStreamSynchronize(stream));
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         MPI_Allreduce(res, _recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&residue, _recv_buf, sizeof(double), cudaMemcpyDeviceToHost);
      }
#endif

      residue = sqrt(residue);
      //printf("%d niter : residue : %f\n", n_iter, residue);
      n_iter += 1;

      /* Convergence test of current iteration */
      cvg = convergence_test_cuda(c, (n_restart * n_k_per_restart) + n_iter,
                              residue, convergence);
      if (cvg != CS_SLES_ITERATING)
        break;

    } /* Needs iterating or k < n_restart */

    /* Inversion of Gamma */
    cudaMemsetAsync(gkj_inv, 0, ((n_k_per_restart + 1) * n_k_per_restart / 2) * sizeof(cs_real_t), stream);
    _compute_inv<<<1, n_iter, 0, stream>>>(gkj_inv, gkj);
    /* Compute the solution */
    _compute_vx<<<gridsize, blocksize, 0, stream>>>(n_rows, n_iter, vx, alpha, zk, wa_size, gkj_inv);

    n_restart += 1;

  } /* Needs iterating */

  if (_aux_vectors != aux_vectors){
    CS_FREE_HD(_aux_vectors);
    CS_FREE_HD(zk);
    CS_FREE_HD(ck);
  }
  cudaGraphDestroy(graph);
  cudaGraphExecDestroy(graphExec);
  cudaFree(alpha);
  cudaFree(gkj);
  cudaFree(_recv_buf);
  cudaFree(res);
  cudaFree(gkj_inv);
  cudaStreamDestroy(stream);
  cudaStreamDestroy(stream1);
  return cvg;

}
#endif


cs_sles_convergence_state_t
gcr_cuda_blas(cs_sles_it_t             *c,
              const cs_matrix_t         *a,
              cs_lnum_t                  diag_block_size,
              cs_sles_it_convergence_t  *convergence,
              const cs_real_t           *rhs,
              cs_real_t                 *restrict vx,
              size_t                     aux_size,
              void                      *aux_vectors)
{
  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = n_rows / 3;
  int DeviceId;
  cudaGetDevice(&DeviceId);
  cudaStream_t stream, stream1;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream);
  size_t size = n_rows*sizeof(cs_real_t);
  cudaMemPrefetchAsync(vx, size, DeviceId, stream);
  cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);
  const cs_real_t  *ad_inv = cs_sles_pc_get_ad_inv(c->setup_data->pc_context);
  cudaMemPrefetchAsync(ad_inv, size, DeviceId, stream1);
  const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
  const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
  const cs_lnum_t  *restrict col_id
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
  const cs_lnum_t  *restrict row_index
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
  const cs_real_t  *restrict val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
  const cs_real_t  *restrict d_val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));

  cs_sles_convergence_state_t cvg;
  double  residue;
  double  *res, *recv_buf;
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  cudaMalloc(&recv_buf, sizeof(cs_real_t));
  cudaMalloc(&res, sizeof(cs_real_t));
#else
  cudaMallocManaged(&recv_buf, sizeof(cs_real_t));
  cudaMallocManaged(&res, sizeof(cs_real_t));
#endif
  cs_real_t *_aux_vectors, *alpha;
  cs_real_t *restrict rk, *restrict zk, *restrict ck;
  cs_real_t *restrict gkj, *restrict gkj_inv;

  /* In case of the standard GCR, n_k_per_restart --> Inf,
   * or stops until convergence*/
  unsigned n_iter = 0;
  const unsigned n_k_per_restart = c->restart_interval;
  size_t wa_size;
  unsigned n_restart = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1 + n_k_per_restart * 2;
    wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
      CS_MALLOC_HD(_aux_vectors, wa_size, cs_real_t, CS_ALLOC_DEVICE);
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      CS_MALLOC_HD(zk, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_DEVICE);
#else
      CS_MALLOC_HD(zk, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);
#endif
      CS_MALLOC_HD(ck, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_DEVICE);
    }
    else{
      _aux_vectors = (cs_real_t *)aux_vectors;                            /* store residuals  */
      zk = _aux_vectors + wa_size;                     /* store inv(M)*r   */
      ck = _aux_vectors + wa_size * (1 + n_k_per_restart);   /* store A*zk */
    }
    rk = _aux_vectors;
  }

#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  cudaMalloc(&alpha, n_k_per_restart*sizeof(cs_real_t));
  cudaMalloc(&gkj, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));
#else
  cudaMallocManaged(&alpha, n_k_per_restart*sizeof(cs_real_t));
  cudaMallocManaged(&gkj, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));
#endif
  cudaMalloc(&gkj_inv, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));

  cvg = CS_SLES_ITERATING;
  const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);
  cs_real_t *zk_n;
  cs_real_t *ck_n;
  unsigned int gridsize_dot;
  if(n_rows > 80*2048){
     gridsize_dot =  (unsigned int)ceil((double)80*2048/BLOCKSIZE);
  }else
     gridsize_dot = (unsigned int)ceil((double)n_rows/BLOCKSIZE);
  const unsigned int blocksize = BLOCKSIZE;
  unsigned int gridsize = (unsigned int)ceil((double)n_rows/blocksize);

  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetStream(handle, stream);

  /* Current Restart */
  while (cvg == CS_SLES_ITERATING) {

    n_iter = 0;

    /* Initialize iterative calculation */
    /*----------------------------------*/
    if (a->db_size[3] == 1){
      cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, vx, 1);
      _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, vx, rk, n_rows);
    }else{
      cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, vx, a->db_size[1]);
      unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
      _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, vx, rk, n_blocks);
    }
    _x_minus_y<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, rhs);
    cublasDdot(handle, n_rows, rk, 1, rk, 1, res);
    cudaMemcpyAsync(&residue, res, sizeof(double), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
#if defined(HAVE_MPI)
    if (c->comm != MPI_COMM_NULL){
      MPI_Allreduce(res, recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      cudaMemcpy(&residue, recv_buf, sizeof(double), cudaMemcpyDeviceToHost);
    }
#endif /* defined(HAVE_MPI) */
    residue = sqrt(residue);
    if (n_restart == 0)
      c->setup_data->initial_residue = residue;

    /* Current Iteration on k */
    /* ---------------------- */

    while (cvg == CS_SLES_ITERATING && n_iter < n_k_per_restart) {

      /* Preconditionning */
      zk_n = zk + n_iter * wa_size;
      ck_n = ck + n_iter * wa_size;
      _sles_pc_poly_apply_jacobi<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ad_inv, rk, zk_n);
      cudaStreamSynchronize(stream);

      /* Halo and Matrix.vector product y = A.x with MSR matrix*/
      if (a->db_size[3] == 1){
         cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, zk_n, 1);
         _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, zk_n, ck_n, n_rows);
      }else{
         cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, zk_n, a->db_size[1]);
         unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
         _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, zk_n, ck_n, n_blocks);
      }

      /* compute the direction ck_n*/
      for(cs_lnum_t jj = 0; jj < n_iter; jj++) {
        cs_real_t *ck_j = ck + jj * wa_size;
        cublasDdot(handle, n_rows, ck_j, 1, ck_n, 1, &gkj[(n_iter + 1) * n_iter / 2 + jj]);
#if defined(HAVE_MPI)
        if (c->comm != MPI_COMM_NULL){
           cudaStreamSynchronize(stream);
           MPI_Allreduce(&gkj[(n_iter + 1) * n_iter / 2 + jj], recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
           cudaMemcpy(&gkj[(n_iter + 1) * n_iter / 2 + jj], recv_buf, sizeof(double),cudaMemcpyDeviceToDevice);
        }
 #endif
        _x_minus_ay<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, ck_j, &gkj[(n_iter + 1) * n_iter / 2 + jj]);
      }
      cublasDdot(handle, n_rows, ck_n, 1, ck_n, 1, &gkj[(n_iter + 1) * n_iter / 2 + n_iter]);
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         cudaStreamSynchronize(stream);
         MPI_Allreduce(&gkj[(n_iter + 1) * n_iter / 2 + n_iter], recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&gkj[(n_iter + 1) * n_iter / 2 + n_iter], recv_buf, sizeof(double),cudaMemcpyDeviceToDevice);
      }
#endif
      _x_divise_a<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, &gkj[(n_iter+1) * n_iter / 2 + n_iter]);
      cublasDdot(handle, n_rows, ck_n, 1, rk, 1, &alpha[n_iter]);
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         cudaStreamSynchronize(stream);
         MPI_Allreduce(&alpha[n_iter], recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&alpha[n_iter], recv_buf, sizeof(double), cudaMemcpyDeviceToDevice);
      }
#endif
      _x_minus_ay<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, ck_n, &alpha[n_iter]);
      cublasDdot(handle, n_rows, rk, 1, rk, 1, res);
      cudaMemcpyAsync(&residue, res, sizeof(double), cudaMemcpyDeviceToHost, stream);
      cudaStreamSynchronize(stream);
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         MPI_Allreduce(res, recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&residue, recv_buf, sizeof(double), cudaMemcpyDeviceToHost);
      }
#endif
      residue = sqrt(residue);
      n_iter += 1;

      /* Convergence test of current iteration */
      cvg = convergence_test_cuda(c, (n_restart * n_k_per_restart) + n_iter,
                              residue, convergence);
      if (cvg != CS_SLES_ITERATING)
        break;

    } /* Needs iterating or k < n_restart */

    /* Inversion of Gamma */
    cudaMemsetAsync(gkj_inv, 0, ((n_k_per_restart + 1) * n_k_per_restart / 2) * sizeof(cs_real_t), stream);
    _compute_inv<<<1, n_iter, 0, stream>>>(gkj_inv, gkj);
    /* Compute the solution */
    _compute_vx<<<gridsize, blocksize, 0, stream>>>(n_rows, n_iter, vx, alpha, zk, wa_size, gkj_inv);

    n_restart += 1;

  } /* Needs iterating */

  if (_aux_vectors != aux_vectors){
    CS_FREE_HD(_aux_vectors);
    CS_FREE_HD(zk);
    CS_FREE_HD(ck);
  }
  cudaFree(alpha);
  cudaFree(gkj);
  cudaFree(recv_buf);
  cudaFree(res);
  cudaFree(gkj_inv);
  cudaStreamDestroy(stream);
  cudaStreamDestroy(stream1);
  return cvg;

}


cs_sles_convergence_state_t
gcr_cuda(cs_sles_it_t              *c,
         const cs_matrix_t         *a,
         cs_lnum_t                  diag_block_size,
         cs_sles_it_convergence_t  *convergence,
         const cs_real_t           *rhs,
         cs_real_t                 *restrict vx,
         size_t                     aux_size,
         void                      *aux_vectors)
{
  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = n_rows / 3;
  int DeviceId;
  cudaGetDevice(&DeviceId);
  cudaStream_t stream, stream1;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream);
  size_t size = n_rows*sizeof(cs_real_t);
  cudaMemPrefetchAsync(vx, size, DeviceId, stream);
  cudaMemPrefetchAsync(rhs, size, DeviceId, stream1);
  const cs_real_t  *ad_inv = cs_sles_pc_get_ad_inv(c->setup_data->pc_context);
  cudaMemPrefetchAsync(ad_inv, size, DeviceId, stream1);
  const cs_matrix_struct_csr_t  *ms = (const cs_matrix_struct_csr_t  *)a->structure;
  const cs_matrix_coeff_msr_t  *mc = (const cs_matrix_coeff_msr_t  *)a->coeffs;
  const cs_lnum_t  *restrict col_id
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->col_id));
  const cs_lnum_t  *restrict row_index
     = (const cs_lnum_t  *restrict)cs_get_device_ptr(const_cast<cs_lnum_t *>(ms->row_index));
  const cs_real_t  *restrict val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->x_val));
  const cs_real_t  *restrict d_val
     = (const cs_real_t  *restrict)cs_get_device_ptr(const_cast<cs_real_t *>(mc->d_val));

  cs_sles_convergence_state_t cvg;
  double  residue;
  double  *res, *recv_buf;
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  cudaMalloc(&recv_buf, sizeof(cs_real_t));
  cudaMalloc(&res, sizeof(cs_real_t));
#else
  cudaMallocManaged(&recv_buf, sizeof(cs_real_t));
  cudaMallocManaged(&res, sizeof(cs_real_t));
#endif
  cs_real_t *_aux_vectors, *alpha;
  cs_real_t *restrict rk, *restrict zk, *restrict ck;
  cs_real_t *restrict gkj, *restrict gkj_inv;

  /* In case of the standard GCR, n_k_per_restart --> Inf,
   * or stops until convergence*/
  unsigned n_iter = 0;
  const unsigned n_k_per_restart = c->restart_interval;
  size_t wa_size;
  unsigned n_restart = 0;


  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1 + n_k_per_restart * 2;
    wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)){
      CS_MALLOC_HD(_aux_vectors, wa_size, cs_real_t, CS_ALLOC_DEVICE);
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      CS_MALLOC_HD(zk, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_DEVICE);
#else
      CS_MALLOC_HD(zk, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);
#endif
      CS_MALLOC_HD(ck, wa_size * n_k_per_restart, cs_real_t, CS_ALLOC_DEVICE);
    }
    else{
      _aux_vectors = (cs_real_t *)aux_vectors;
      zk = _aux_vectors + wa_size;                     /* store inv(M)*r   */
      ck = _aux_vectors + wa_size * (1 + n_k_per_restart);   /* store A*zk */
    }
    rk = _aux_vectors;                               /* store residuals  */
  }


#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  cudaMalloc(&alpha, n_k_per_restart*sizeof(cs_real_t));
#else
  cudaMallocManaged(&alpha, n_k_per_restart*sizeof(cs_real_t));
#endif
  /* gkj stores the upper triangle matrix Gamma of crossed inner-products
   * Not necessary to allocate the full matrix size
   * gkj_inv stores the inverse of gkj */
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
   cudaMalloc(&gkj, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));
#else
   cudaMallocManaged(&gkj, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));
#endif

  cudaMalloc(&gkj_inv, ((n_k_per_restart + 1) * n_k_per_restart / 2)*sizeof(cs_real_t));

  cvg = CS_SLES_ITERATING;
  const cs_real_t  *restrict ad = cs_matrix_get_diagonal_gpu(a);
  cs_real_t *zk_n;
  cs_real_t *ck_n;
  unsigned int gridsize_dot;
  if(n_rows > 80*2048){
     gridsize_dot =  (unsigned int)ceil((double)80*2048/BLOCKSIZE);
  }else
     gridsize_dot = (unsigned int)ceil((double)n_rows/BLOCKSIZE);
  const unsigned int blocksize = BLOCKSIZE;
  const unsigned int blocksize2 = 640;
  unsigned int gridsize = (unsigned int)ceil((double)n_rows/blocksize);
  unsigned int gridsize2 = 1;
  //gridsize_dot = gridsize;
  double  *sum_block;
  cudaMalloc(&sum_block, gridsize*sizeof(double));

  /* Current Restart */
  while (cvg == CS_SLES_ITERATING) {

    n_iter = 0;
    /* Initialize iterative calculation */
    /*----------------------------------*/
    if (a->db_size[3] == 1){
      if (cs_glob_n_ranks == 1 && a->halo != NULL)
        local_halo_cuda(a->halo, CS_HALO_STANDARD, vx, 1, stream);
      else if (a->halo != NULL)
        cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, vx, 1);
      _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, vx, rk, n_rows);
    }else{
      if (cs_glob_n_ranks == 1 && a->halo != NULL)
        local_halo_cuda(a->halo, CS_HALO_STANDARD, vx, a->db_size[1], stream);
      else if (a->halo != NULL)
        cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, vx, a->db_size[1]);
      unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
      _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, vx, rk, n_blocks);
    }
    _x_minus_y<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, rhs);
    dot_xx<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, sum_block);
    _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, res, gridsize_dot);
    cudaMemcpyAsync(&residue, res, sizeof(double), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
#if defined(HAVE_MPI)
    if (c->comm != MPI_COMM_NULL){
      MPI_Allreduce(res, recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      cudaMemcpy(&residue, recv_buf, sizeof(double), cudaMemcpyDeviceToHost);
    }
#endif /* defined(HAVE_MPI) */
    residue = sqrt(residue);
    if (n_restart == 0)
      c->setup_data->initial_residue = residue;

    /* Current Iteration on k */
    /* ---------------------- */

    while (cvg == CS_SLES_ITERATING && n_iter < n_k_per_restart) {

      /* Preconditionning */
      zk_n = zk + n_iter * wa_size;
      ck_n = ck + n_iter * wa_size;

      /* Preconditionning */
      _sles_pc_poly_apply_jacobi<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ad_inv, rk, zk_n);
      cudaStreamSynchronize(stream);

      /* Halo and Matrix.vector product y = A.x with MSR matrix*/
      if (a->db_size[3] == 1){
        if (cs_glob_n_ranks == 1 && a->halo != NULL)
          local_halo_cuda(a->halo, CS_HALO_STANDARD, zk_n, 1, stream);
        else if (a->halo != NULL)
          cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, zk_n, 1);
        _mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, zk_n, ck_n, n_rows);
      }else{
        if (cs_glob_n_ranks == 1 && a->halo != NULL)
          local_halo_cuda(a->halo, CS_HALO_STANDARD, zk_n, a->db_size[1], stream);
        else if (a->halo != NULL)
          cs_halo_sync_var_cuda_aware(a->halo, CS_HALO_STANDARD, zk_n, a->db_size[1]);
        unsigned int gridsize = (unsigned int)ceil((double)n_blocks/blocksize);
        _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, stream>>>(false, col_id, row_index, val, d_val, zk_n, ck_n, n_blocks);
      }

      /* compute the direction ck_n*/
      for(cs_lnum_t jj = 0; jj < n_iter; jj++) {
        cs_real_t *ck_j = ck + jj * wa_size;
        dot_xy<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_j, ck_n, sum_block);
        _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, &gkj[(n_iter + 1) * n_iter / 2 + jj], gridsize_dot);
#if defined(HAVE_MPI)
       if (c->comm != MPI_COMM_NULL){
          cudaStreamSynchronize(stream);
          MPI_Allreduce(&gkj[(n_iter + 1) * n_iter / 2 + jj], recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
          cudaMemcpy(&gkj[(n_iter + 1) * n_iter / 2 + jj], recv_buf, sizeof(double),cudaMemcpyDeviceToDevice);
       }
#endif
        _x_minus_ay<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, ck_j, &gkj[(n_iter + 1) * n_iter / 2 + jj]);
      }
       dot_xx<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, sum_block);
      _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, &gkj[(n_iter+1) * n_iter / 2 + n_iter], gridsize_dot);
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         cudaStreamSynchronize(stream);
         MPI_Allreduce(&gkj[(n_iter + 1) * n_iter / 2 + n_iter], recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&gkj[(n_iter + 1) * n_iter / 2 + n_iter], recv_buf, sizeof(double),cudaMemcpyDeviceToDevice);
      }
#endif
      _x_divise_a<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, &gkj[(n_iter+1) * n_iter / 2 + n_iter]);
      dot_xy<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, ck_n, rk, sum_block);
      _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, &alpha[n_iter], gridsize_dot);
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         cudaStreamSynchronize(stream);
         MPI_Allreduce(&alpha[n_iter], recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&alpha[n_iter], recv_buf, sizeof(double), cudaMemcpyDeviceToDevice);
      }
#endif
      _x_minus_ay<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, ck_n, &alpha[n_iter]);
      dot_xx<<<gridsize_dot, blocksize, 0, stream>>>(n_rows, rk, sum_block);
      _reduceCUDA<blocksize2><<<gridsize2, blocksize2, 0, stream>>>(sum_block, res, gridsize_dot);
      cudaMemcpyAsync(&residue, res, sizeof(double), cudaMemcpyDeviceToHost, stream);
      cudaStreamSynchronize(stream);
#if defined(HAVE_MPI)
      if (c->comm != MPI_COMM_NULL){
         MPI_Allreduce(res, recv_buf, 1, MPI_DOUBLE, MPI_SUM, c->comm);
         cudaMemcpy(&residue, recv_buf, sizeof(double), cudaMemcpyDeviceToHost);
      }
#endif
      residue = sqrt(residue);
      n_iter += 1;

      /* Convergence test of current iteration */
      cvg = convergence_test_cuda(c, (n_restart * n_k_per_restart) + n_iter,
                              residue, convergence);

      if (cvg != CS_SLES_ITERATING)
        break;

    } /* Needs iterating or k < n_restart */

    /* Inversion of Gamma */

    cudaMemsetAsync(gkj_inv, 0, ((n_k_per_restart + 1) * n_k_per_restart / 2) * sizeof(cs_real_t), stream);
    _compute_inv<<<1, n_iter, 0, stream>>>(gkj_inv, gkj);
    /* Compute the solution */
    _compute_vx<<<gridsize, blocksize, 0, stream>>>(n_rows, n_iter, vx, alpha, zk, wa_size, gkj_inv);

    n_restart += 1;

  } /* Needs iterating */

  if (_aux_vectors != aux_vectors){
    CS_FREE_HD(_aux_vectors);
    CS_FREE_HD(ck);
    CS_FREE_HD(zk);
  }

  cudaFree(alpha);
  cudaFree(gkj);
  cudaFree(recv_buf);
  cudaFree(res);
  cudaFree(gkj_inv);
  cudaStreamDestroy(stream);
  cudaStreamDestroy(stream1);
  return cvg;

}

/*----------------------------------------------------------------------------
 * build buffer for MPI.
 *
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   var       <-> pointer to variable value array
 *----------------------------------------------------------------------------*/

__global__ static void _build_buffer_cuda(cs_lnum_t     start,
                                         cs_lnum_t      length,
                                         cs_lnum_t      *send_list,
                                         cs_real_t      *var,
                                         cs_real_t      *build_buffer)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if(ii < length)
     build_buffer[start + ii] = var[send_list[start + ii]];
}

__global__ static void _build_buffer_cuda_stride(cs_lnum_t     start,
                                                 cs_lnum_t      length,
                                                 cs_lnum_t      *send_list,
                                                 cs_real_t      *var,
                                                 cs_real_t      *build_buffer,
                                                 int             stride)
{
  cs_lnum_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i < length){
    for (cs_lnum_t j = 0; j < stride; j++)
      build_buffer[(start + i)*stride + j] = var[send_list[start + i]*stride + j];
  }

}

void
build_buffer_cuda( const cs_halo_t  *halo,
                   cs_halo_type_t    sync_mode,
                   cs_real_t         var[],
                   cs_real_t         *build_buffer)
{
  cs_lnum_t start, length;
  const cs_lnum_t end_shift = (sync_mode == CS_HALO_STANDARD) ? 1 : 2;
  cs_lnum_t *send_list = halo->send_list_gpu;
  cudaStream_t nstream[halo->n_c_domains];
  for (cs_lnum_t ii = 0; ii < halo->n_c_domains; ii++)
      cudaStreamCreate(&nstream[ii]);

  for (cs_lnum_t rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
      start = halo->send_index[2*rank_id];
      length =   halo->send_index[2*rank_id + end_shift]
               - halo->send_index[2*rank_id];
      if(length > 0){
        unsigned int blocksize = (length < 256)? length:256;
        unsigned int gridsize = (unsigned int)ceil((double)length/blocksize);
        _build_buffer_cuda<<<gridsize, blocksize, 0, nstream[rank_id]>>>(start, length, send_list, var, build_buffer);
      }
  }

  for (cs_lnum_t ii = 0; ii < halo->n_c_domains; ii++){
    CS_CUDA_CHECK(cudaStreamSynchronize(nstream[ii]));
    cudaStreamDestroy(nstream[ii]);
  }
}



void
build_buffer_cuda_stride( const cs_halo_t  *halo,
                          cs_halo_type_t    sync_mode,
                          cs_real_t         var[],
                          cs_real_t         *build_buffer,
                          int               stride,
                          cs_lnum_t         end_shift)
{
  cs_lnum_t start, length;
  cs_lnum_t *send_list = halo->send_list_gpu;
  cudaStream_t nstream[halo->n_c_domains];

  for (cs_lnum_t ii = 0; ii < halo->n_c_domains; ii++)
    cudaStreamCreate(&nstream[ii]);

  for (cs_lnum_t rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start = halo->send_index[2*rank_id];
      length =   halo->send_index[2*rank_id + end_shift]
               - halo->send_index[2*rank_id];
      if(length > 0){
        unsigned int blocksize = (length < 256)? length:256;
        unsigned int gridsize = (unsigned int)ceil((double)length/blocksize);
        _build_buffer_cuda_stride<<<gridsize, blocksize, 0, nstream[rank_id]>>>(start, length, send_list, var, build_buffer, stride);
      }
  }

  for (cs_lnum_t ii = 0; ii < halo->n_c_domains; ii++){
    CS_CUDA_CHECK(cudaStreamSynchronize(nstream[ii]));
    cudaStreamDestroy(nstream[ii]);
  }
}

/*----------------------------------------------------------------------------
 * Update array of variable (floating-point) halo values in case of periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   var       <-> pointer to variable value array
 *----------------------------------------------------------------------------*/

__global__ static void _local_halo_cuda(cs_lnum_t      start,
                                        cs_lnum_t      length,
                                        cs_lnum_t      *send_list,
                                        cs_real_t      *var,
                                        cs_real_t      *recv_var)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if(ii < length)
     recv_var[ii] = var[send_list[start + ii]];
}

__global__ static void _local_halo_stride_cuda(cs_lnum_t      start,
                                               cs_lnum_t      length,
                                               cs_lnum_t      *send_list,
                                               cs_real_t      *var,
                                               cs_real_t      *recv_var,
                                               int             stride)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if(ii < length){
     for (cs_lnum_t j = 0; j < stride; j++)
        recv_var[ii] = var[send_list[start + ii]*stride + j];
  }
}


void
local_halo_cuda( const cs_halo_t  *halo,
                 cs_halo_type_t    sync_mode,
                 cs_real_t         var[],
                 int               stride,
                 cudaStream_t      stream)
{
  cs_lnum_t start, length, end_shift;
  if (sync_mode == CS_HALO_STANDARD)
    end_shift = 1;

  else if (sync_mode == CS_HALO_EXTENDED)
    end_shift = 2;

  cs_real_t *recv_var
  = var + (halo->n_local_elts + halo->index[0])*stride;

  start = halo->send_index[0];
  length =  halo->send_index[end_shift]
         - halo->send_index[0];

  unsigned int blocksize = (length < 256)? length:256;
  unsigned int gridsize = (unsigned int)ceil((double)length/blocksize);
  cs_lnum_t *send_list = halo->send_list;
  if(stride == 1)
     _local_halo_cuda<<<gridsize, blocksize, 0, stream>>>(start, length, send_list, var, recv_var);
  else
     _local_halo_stride_cuda<<<gridsize, blocksize, 0, stream>>>(start, length, send_list, var, recv_var, stride);
  CS_CUDA_CHECK(cudaGetLastError());

}

END_C_DECLS