#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <chrono>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <cuda_runtime_api.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base_accel.h"
#include "cs_base_cuda.h"
#include "cs_blas.h"
#include "cs_cell_to_vertex.h"
#include "cs_ext_neighborhood.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"

#include "cs_convection_diffusion.h"
#include "cs_convection_diffusion_priv.h"

#include "cs_convection_diffusion_cuda_scatter.cuh"
#include "cs_convection_diffusion_cuda_gather.cuh"



template <typename T>
static void
_sync_or_copy_real_h2d_c(const  T   *val_h,
                       cs_lnum_t           n_vals,
                       int                 device_id,
                       cudaStream_t        stream,
                       const T   **val_d,
                       void              **buf_d)
{
  const T  *_val_d = NULL;
  void             *_buf_d = NULL;

  cs_alloc_mode_t alloc_mode = cs_check_device_ptr(val_h);
  size_t size = n_vals * sizeof(T);

  if (alloc_mode == CS_ALLOC_HOST) {
    CS_CUDA_CHECK(cudaMalloc(&_buf_d, size));
    cs_cuda_copy_h2d(_buf_d, val_h, size);
    _val_d = (const T *)_buf_d;
  }
  else {
    _val_d = (const T *)cs_get_device_ptr((void *)val_h);

    if (alloc_mode == CS_ALLOC_HOST_DEVICE_SHARED)
      cudaMemPrefetchAsync(val_h, size, device_id, stream);
    else
      cs_sync_h2d(val_h);
  }

  *val_d = _val_d;
  *buf_d = _buf_d;
}

/*----------------------------------------------------------------------------
 * _gradient_vector the gradient of a vector using a given gradient of
 * this vector (typically lsq).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   r_grad         --> gradient used for reconstruction
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/
extern "C" void
cs_convection_diffusion_vector_cuda(const cs_mesh_t             *mesh,
                                    const cs_mesh_adjacencies_t *madj,
                                    const cs_mesh_quantities_t  *fvq,
                                    const cs_real_3_t  *restrict pvar,
                                    const cs_real_t              i_massflux[],
                                    const cs_real_33_t          *grad,
                                    cs_real_33_t                *grdpa,
                                    const cs_real_3_t  *restrict coefav,
                                    const cs_real_33_t *restrict coefbv,
                                    const int                    inc,
                                    const bool                   flag1,
                                    const bool                   flag2,
                                    const bool                   perf)
{
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_b_cells = mesh->n_b_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  int device_id;
  cudaGetDevice(&device_id);

  cudaStream_t stream;
  cudaStreamCreate(&stream);
  
  cudaEvent_t start, mem_h2d, init, f_i, f_b, f_f, stop;
  float msec = 0.0f;
  CS_CUDA_CHECK(cudaEventCreate(&start));
  CS_CUDA_CHECK(cudaEventCreate(&mem_h2d));
  CS_CUDA_CHECK(cudaEventCreate(&init));
  CS_CUDA_CHECK(cudaEventCreate(&f_i));
  CS_CUDA_CHECK(cudaEventCreate(&f_b));
  CS_CUDA_CHECK(cudaEventCreate(&f_f));
  CS_CUDA_CHECK(cudaEventCreate(&stop));


  // Record the start event
  CS_CUDA_CHECK(cudaEventRecord(start, stream));

  unsigned int blocksize = 256;

  cs_real_33_t *grad_d = NULL;
  CS_CUDA_CHECK(cudaMalloc(&grad_d, n_cells_ext * sizeof(cs_real_33_t)));

  cs_real_33_t *grdpa_d;
  CS_CUDA_CHECK(cudaMalloc(&grdpa_d, n_cells_ext * sizeof(cs_real_33_t)));

  cs_gnum_t n_upwind;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)cs_get_device_ptr_const_pf(mesh->i_face_cells);

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(mesh->b_face_cells);

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_cen);

  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->diipb);

  const cs_real_3_t *restrict b_f_face_normal
  = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_f_face_normal);

  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_cells_idx);

  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_cells);

  const cs_lnum_t *restrict b_cells
    = (cs_lnum_t *restrict)cs_get_device_ptr_const_pf(mesh->b_cells);

  const cs_lnum_t *restrict cell_b_faces
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_b_faces);

  const cs_lnum_t *restrict cell_b_faces_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_b_faces_idx);

  cs_real_t *restrict i_massflux_d;
  CS_CUDA_CHECK(cudaMalloc(&i_massflux_d, sizeof(cs_real_t)*n_i_faces));
  cs_cuda_copy_h2d(i_massflux_d, (void *)i_massflux, sizeof(cs_real_t)*n_i_faces);

  cs_real_3_t *restrict i_face_cog;
  CS_CUDA_CHECK(cudaMalloc(&i_face_cog, sizeof(cs_real_3_t)*n_i_faces));
  cs_cuda_copy_h2d(i_face_cog, (void *)fvq->i_face_cog, sizeof(cs_real_3_t)*n_i_faces);

  cs_real_3_t *restrict i_f_face_normal;
  CS_CUDA_CHECK(cudaMalloc(&i_f_face_normal, sizeof(cs_real_3_t)*n_i_faces));
  cs_cuda_copy_h2d(i_f_face_normal, (void *)fvq->i_f_face_normal, sizeof(cs_real_3_t)*n_i_faces);
  
  cs_real_t *restrict cell_vol;
  CS_CUDA_CHECK(cudaMalloc(&cell_vol, sizeof(cs_real_t)*n_cells));
  cs_cuda_copy_h2d(cell_vol, (void *)fvq->cell_vol, sizeof(cs_real_t)*n_cells);

  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t n_cells_i_face = (madj->cell_cells_idx[n_cells]);

  cs_lnum_t *restrict cell_i_faces;
  CS_CUDA_CHECK(cudaMalloc(&cell_i_faces, sizeof(cs_lnum_t)*n_cells_i_face));
  cs_cuda_copy_h2d(cell_i_faces, madj->cell_i_faces, sizeof(cs_lnum_t)*n_cells_i_face);

  short int *restrict cell_i_faces_sgn;
  CS_CUDA_CHECK(cudaMalloc(&cell_i_faces_sgn, sizeof(short int)*n_cells_i_face));
  cs_cuda_copy_h2d(cell_i_faces_sgn, madj->cell_i_faces_sgn, sizeof(short int)*n_cells_i_face);


  void *_coefb_d, *_coefa_d, *_pvar_d;

  const cs_real_3_t * coefa_d = NULL;
  const cs_real_3_t * pvar_d = NULL;
  const cs_real_33_t * coefb_d = NULL;

  /* Initialization */

  _sync_or_copy_real_h2d_c(pvar, n_cells_ext, device_id, stream,
    &pvar_d, &_pvar_d);
  _sync_or_copy_real_h2d_c(coefav, n_b_faces, device_id, stream,
        &coefa_d, &_coefa_d);
  _sync_or_copy_real_h2d_c(coefbv, n_b_faces, device_id, stream,
        &coefb_d, &_coefb_d);

  if(flag1){
    cs_cuda_copy_h2d(grad_d, grad, sizeof(cs_real_33_t)*n_cells_ext);
  }
  else{
    cudaMemset(grad_d, 0, n_cells_ext * sizeof(cs_real_33_t));
  }

  CS_CUDA_CHECK(cudaEventRecord(mem_h2d, stream));

  cudaMemset(grdpa_d, 0, n_cells_ext * sizeof(cs_real_33_t));

  CS_CUDA_CHECK(cudaEventRecord(init, stream));

  if (flag2) {
    // cs_slope_test_gradient_vector_cuda_i<<<(unsigned int)ceil((double)n_i_faces / blocksize), blocksize, 0, stream>>>
    //                               (n_i_faces,
    //                                i_face_cells,
    //                                i_face_cog,
    //                                cell_cen,
    //                                pvar_d,
    //                                i_massflux_d,
    //                                i_f_face_normal,
    //                                grad_d,
    //                                grdpa_d);
    

    cs_slope_test_gradient_vector_cuda_i_gather<<<(unsigned int)ceil((double)n_cells / blocksize), blocksize, 0, stream>>>
                                  (n_cells,
                                  i_face_cog,
                                  cell_cen,
                                  pvar_d,
                                  i_massflux_d,
                                  i_f_face_normal,
                                  cell_cells_idx,
                                  cell_cells,
                                  cell_i_faces,
                                  cell_i_faces_sgn,
                                  grad_d,
                                  grdpa_d);

    CS_CUDA_CHECK(cudaEventRecord(f_i, stream));

    // cs_slope_test_gradient_vector_cuda_b<<<(unsigned int)ceil((double)n_b_faces / blocksize), blocksize, 0, stream>>>
    //                                 (n_b_faces,
    //                                  pvar_d,
    //                                  b_face_cells,
    //                                  diipb,
    //                                  inc,
    //                                  coefa_d,
    //                                  coefb_d,
    //                                  b_f_face_normal,
    //                                  grad_d,
    //                                  grdpa_d);


    cs_slope_test_gradient_vector_cuda_b_gather<<<(unsigned int)ceil((double)n_b_cells / blocksize), blocksize, 0, stream>>>
                                    (n_b_cells,
                                    pvar_d,
                                    diipb,
                                    inc,
                                    coefa_d,
                                    coefb_d,
                                    b_f_face_normal,
                                    b_cells,
                                    cell_b_faces,
                                    cell_b_faces_idx,
                                    grad_d,
                                    grdpa_d);

    CS_CUDA_CHECK(cudaEventRecord(f_b, stream));

    cs_slope_test_gradient_vector_cuda_f<<<(unsigned int)ceil((double)n_cells * 3 * 3 / blocksize), blocksize, 0, stream>>>
                                    (n_cells * 3 * 3,
                                     cell_vol,
                                     grdpa_d);

    CS_CUDA_CHECK(cudaEventRecord(f_f, stream));
    
  }
 
  n_upwind = 0;
  
  /* Sync to host */
  if (grdpa_d != NULL) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * 3 * 3;
    cs_cuda_copy_d2h(grdpa, grdpa_d, size);
  }
  else
    cs_sync_d2h(grdpa);
  
  CS_CUDA_CHECK(cudaEventRecord(stop, stream));
  CS_CUDA_CHECK(cudaEventSynchronize(stop));

  cudaStreamSynchronize(stream);
  cudaStreamDestroy(stream);

  if(perf){
    printf("convection_diffusion Kernels:\n");
    printf("Execution time in us: \t");

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, mem_h2d, init));
    printf("Init = %f\t", msec*1000.f);
    
    if (flag2) {
      msec = 0.0f;
      CS_CUDA_CHECK(cudaEventElapsedTime(&msec, init, f_i));
      printf("f_i = %f\t", msec*1000.f);

      msec = 0.0f;
      CS_CUDA_CHECK(cudaEventElapsedTime(&msec, f_i, f_b));
      printf("f_b = %f\t", msec*1000.f);

      msec = 0.0f;
      CS_CUDA_CHECK(cudaEventElapsedTime(&msec, f_b, f_f));
      printf("f_f = %f\t", msec*1000.f);

    }

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, start, stop));
    printf("Total = %f\n", msec*1000.f);
  }

  if (!flag1){
    CS_CUDA_CHECK(cudaFree(grad_d));
  }

  if (_pvar_d != NULL)
    CS_CUDA_CHECK(cudaFree(_pvar_d));
  if (_coefa_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefa_d));
  if (_coefb_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefb_d));
    
  CS_CUDA_CHECK(cudaFree(grdpa_d));
  CS_CUDA_CHECK(cudaFree(i_massflux_d));
  CS_CUDA_CHECK(cudaFree(i_f_face_normal));
  CS_CUDA_CHECK(cudaFree(cell_vol));
  CS_CUDA_CHECK(cudaFree(cell_i_faces));
  CS_CUDA_CHECK(cudaFree(cell_i_faces_sgn));
}
