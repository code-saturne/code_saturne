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

#include "cs_convection_diffusion_cuda.cuh"


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
                                    const cs_mesh_quantities_t  *fvq,
                                    const cs_real_3_t  *restrict pvar,
                                    const cs_real_t              i_massflux[],
                                    const cs_real_33_t          *grad,
                                    cs_real_33_t                *grdpa,
                                    cs_real_3_t       *restrict  rhs,
                                    const cs_real_3_t  *restrict coefav,
                                    const cs_real_33_t *restrict coefbv,
                                    const int                    inc,
                                    const bool                   flag1,
                                    const bool                   flag2,
                                    const bool                   perf)
{
  const cs_lnum_t n_cells = mesh->n_cells;
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

  cs_real_33_t *grdpa_d;
  CS_CUDA_CHECK(cudaMalloc(&grdpa_d, n_cells_ext * sizeof(cs_real_33_t)));

  cs_real_3_t *rhs_d;
  CS_CUDA_CHECK(cudaMalloc(&rhs_d, n_cells_ext * sizeof(cs_real_3_t)));
  cs_cuda_copy_h2d(rhs_d, rhs, sizeof(cs_real_3_t)*n_cells_ext);

  cs_gnum_t n_upwind;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)cs_get_device_ptr_const_pf(mesh->i_face_cells);

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(mesh->b_face_cells);

  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->i_face_cog);

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_cen);

  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->diipb);

  const cs_real_3_t *restrict b_f_face_normal
  = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_f_face_normal);

  cs_real_3_t *restrict i_f_face_normal;
  CS_CUDA_CHECK(cudaMalloc(&i_f_face_normal, sizeof(cs_real_3_t)*n_i_faces));
  cs_cuda_copy_h2d(i_f_face_normal, (void *)fvq->i_f_face_normal, sizeof(cs_real_3_t)*n_i_faces);
  
  cs_real_3_t *restrict coefa_d;
  CS_CUDA_CHECK(cudaMalloc(&coefa_d, sizeof(cs_real_3_t)*n_b_faces));
  cs_cuda_copy_h2d(coefa_d, (void *)coefav, sizeof(cs_real_3_t)*n_b_faces);

  cs_real_33_t *restrict coefb_d;
  CS_CUDA_CHECK(cudaMalloc(&coefb_d, sizeof(cs_real_33_t)*n_b_faces));
  cs_cuda_copy_h2d(coefb_d, (void *)coefbv, sizeof(cs_real_33_t)*n_b_faces);

  cs_real_t *restrict cell_vol;
  CS_CUDA_CHECK(cudaMalloc(&cell_vol, n_cells * sizeof(cs_real_t)));
  cs_cuda_copy_h2d(cell_vol, (void *)fvq->cell_vol, sizeof(cs_real_t)*n_cells);

  cs_real_3_t *restrict pvar_d;
  CS_CUDA_CHECK(cudaMalloc(&pvar_d, sizeof(cs_real_3_t)*n_cells));
  cs_cuda_copy_h2d(pvar_d, (void *)pvar, sizeof(cs_real_3_t)*n_cells);

  /* Initialization */

  CS_CUDA_CHECK(cudaEventRecord(mem_h2d, stream));

  if(flag1){
    cs_cuda_copy_h2d(grad_d, grad, sizeof(cs_real_33_t)*n_cells_ext);
  }else{
    CS_CUDA_CHECK(cudaMalloc(&grad_d, n_cells_ext * sizeof(cs_real_33_t)));
    cudaMemset(grad_d, 0, n_cells_ext * sizeof(cs_real_33_t));
  }

  cudaMemset(grdpa_d, 0, n_cells_ext * sizeof(cs_real_33_t));

  CS_CUDA_CHECK(cudaEventRecord(init, stream));

  if (flag2) {
    printf(" On passe dans cs_slope_test");
    cs_slope_test_gradient_vector_cuda_i<<<n_i_faces / blocksize, blocksize, 0, stream>>>
                                  (n_i_faces,
                                   i_face_cells,
                                   i_face_cog,
                                   cell_cen,
                                   pvar_d,
                                   i_massflux,
                                   i_f_face_normal,
                                   grad_d,
                                   grdpa_d);

    CS_CUDA_CHECK(cudaEventRecord(f_i, stream));

    cs_slope_test_gradient_vector_cuda_b<<<n_b_faces / blocksize, blocksize, 0, stream>>>
                                    (n_b_faces,
                                     pvar_d,
                                     b_face_cells,
                                     diipb,
                                     inc,
                                     coefa_d,
                                     coefb_d,
                                     b_f_face_normal,
                                     grad_d,
                                     grdpa_d);

    CS_CUDA_CHECK(cudaEventRecord(f_b, stream));

    // cs_slope_test_gradient_vector_cuda_f<<<(n_cells / blocksize) * 3 * 3, blocksize, 0, stream>>>
    //                                 (n_cells,
    //                                  cell_vol,
    //                                  grdpa_d);

    CS_CUDA_CHECK(cudaEventRecord(f_f, stream));
    
  }
 
  n_upwind = 0;

  if(n_cells_ext > n_cells){
    cudaMemset(rhs_d[n_cells], 0, (n_cells_ext-n_cells) * sizeof(cs_real_3_t));
  }

  
  /* Sync to host */
  if (grdpa_d != NULL) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * 3 * 3;
    cs_cuda_copy_d2h(grdpa, grdpa_d, size);
  }
  else
    cs_sync_d2h(grdpa_d);

  
  /* Sync to host */
  if (rhs_d != NULL) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * 3;
    cs_cuda_copy_d2h(rhs, rhs_d, size);
  }
  else
    cs_sync_d2h(rhs_d);
  
  CS_CUDA_CHECK(cudaEventRecord(stop, stream));
  CS_CUDA_CHECK(cudaEventSynchronize(stop));

  cudaStreamSynchronize(stream);
  cudaStreamDestroy(stream);

  if(perf){
    printf("convection_diffusion Kernels times:\t");

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, mem_h2d, init));
    printf("Kernels execution time in us: \t");
    printf("Init = %f\t", msec*1000.f);
    
    if (flag2) {
      msec = 0.0f;
      CS_CUDA_CHECK(cudaEventElapsedTime(&msec, init, f_i));
      printf("Kernels execution time in us: \t");
      printf("Init = %f\t", msec*1000.f);

      msec = 0.0f;
      CS_CUDA_CHECK(cudaEventElapsedTime(&msec, f_i, f_b));
      printf("Kernels execution time in us: \t");
      printf("Init = %f\t", msec*1000.f);

      msec = 0.0f;
      CS_CUDA_CHECK(cudaEventElapsedTime(&msec, f_b, f_f));
      printf("Kernels execution time in us: \t");
      printf("Init = %f\t", msec*1000.f);

    }

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, start, stop));
    printf("Total = %f\t", msec*1000.f);

    printf("\n");
  }

  if (!flag1){
    CS_CUDA_CHECK(cudaFree(grad_d));
  }

  CS_CUDA_CHECK(cudaFree(grdpa_d));
  CS_CUDA_CHECK(cudaFree(rhs_d));
  CS_CUDA_CHECK(cudaFree(i_f_face_normal));
  CS_CUDA_CHECK(cudaFree(coefa_d));
  CS_CUDA_CHECK(cudaFree(coefb_d));
  CS_CUDA_CHECK(cudaFree(cell_vol));
  CS_CUDA_CHECK(cudaFree(pvar_d));
}
