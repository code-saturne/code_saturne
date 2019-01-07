/*============================================================================
 * CUDA offloading support
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) IBM Corp. 2017, 2018

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

#ifndef __CS_CUDA_H__
#define __CS_CUDA_H__

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

//
// Utility to execute function controlled by a conditional.
//
#define CS_CUDA_GPU_THRESHOLD GPU_THRESHOLD

//
// If CUDA offload is not defined, the CUDA related entry points
// do not produce any action
//

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifdef HAVE_CUDA_OFFLOAD

void cs_cuda_initialize(void);
void cs_cuda_finalize(void);
void cs_cuda_map_alloc(const void *Pointer, size_t Size);
void cs_cuda_map_to(const void *Pointer, size_t Size);
void cs_cuda_map_from(const void *Pointer, size_t Size);
void cs_cuda_map_from_sync(const void *Pointer, size_t Size);
void cs_cuda_map_release(const void *Pointer, size_t Size);

int cs_cuda_mat_vec_p_l_csr(bool exclude_diag,
                            const cs_lnum_t *restrict row_index,
                            const cs_lnum_t *restrict col_id,
                            const cs_real_t *restrict val,
                            const cs_real_t *restrict x, cs_real_t *restrict y,
                            cs_lnum_t n_rows, cs_lnum_t n_cols);

int cs_cuda_mat_vec_p_l_msr(bool exclude_diag,
                            const cs_lnum_t *restrict row_index,
                            const cs_lnum_t *restrict col_id,
                            const cs_real_t *restrict x_val,
                            const cs_real_t *restrict d_val,
                            const cs_real_t *restrict x, cs_real_t *restrict y,
                            cs_lnum_t n_rows, cs_lnum_t n_cols);

int cs_cuda_seidel_forward(const cs_lnum_t *restrict row_index,
                           const cs_lnum_t *restrict col_id,
                           const cs_real_t *restrict val,
                           const cs_real_t *restrict ad_inv,
                           const cs_real_t *restrict rhs,
                           cs_real_t *restrict vx, cs_lnum_t diag_block_size,
                           const int *db_size, cs_lnum_t n_rows,
                           cs_lnum_t n_cols);

int cs_cuda_seidel_backward(
    const cs_lnum_t *restrict row_index, const cs_lnum_t *restrict col_id,
    const cs_real_t *restrict val, const cs_real_t *restrict ad_inv,
    const cs_real_t *restrict ad, const cs_real_t *restrict rhs,
    cs_real_t *restrict vx, cs_real_t *restrict red, cs_lnum_t diag_block_size,
    const int *db_size, cs_lnum_t n_rows, cs_lnum_t n_cols);

int cs_cuda_dot_product_xx(cs_real_t *restrict xx, const cs_real_t *restrict x,
                           cs_lnum_t n_rows);
int cs_cuda_dot_product_xx_xy(cs_real_t *restrict xx, cs_real_t *restrict xy,
                              const cs_real_t *restrict x,
                              const cs_real_t *restrict y, cs_lnum_t n_rows);
int cs_cuda_dot_product_xy_yz(cs_real_t *restrict xy, cs_real_t *restrict yz,
                              const cs_real_t *restrict x,
                              const cs_real_t *restrict y,
                              const cs_real_t *restrict z, cs_lnum_t n_rows);
int cs_cuda_vector_vc_equal_zero_if_exists(cs_lnum_t n_rows, cs_lnum_t n_elems,
                                           cs_real_t *restrict vc);
int cs_cuda_vector_vc_equal_va(cs_lnum_t n_rows, cs_real_t *restrict vc,
                               const cs_real_t *restrict va);
int cs_cuda_vector_vc_sub_equal_va(cs_lnum_t n_rows, cs_real_t *restrict vc,
                                   const cs_real_t *restrict va);
int cs_cuda_vector_vc_mul_equal_va(cs_lnum_t n_rows, cs_real_t *restrict vc,
                                   const cs_real_t *restrict va);
int cs_cuda_vector_vc_equal_va_mul_vb(cs_lnum_t n_rows, cs_real_t *restrict vc,
                                      const cs_real_t *restrict va,
                                      const cs_real_t *restrict vb);
int cs_cuda_vector_vc_add_equal_s_mul_vb(cs_lnum_t n_rows, cs_real_t s,
                                         cs_real_t *restrict vc1,
                                         const cs_real_t *restrict vb1,
                                         cs_real_t *restrict vc2,
                                         const cs_real_t *restrict vb2);
int cs_cuda_vector_vc_equal_va_add_s_mul_vb(cs_lnum_t n_rows, cs_real_t s,
                                            cs_real_t *restrict vc,
                                            const cs_real_t *restrict va,
                                            const cs_real_t *restrict vb);
int cs_cuda_move_to_device_if_exists(cs_lnum_t n_rows, cs_lnum_t n_elems,
                                     cs_real_t *restrict vc);
int cs_cuda_move_from_device_if_exists(cs_lnum_t n_rows, cs_lnum_t n_elems,
                                       cs_real_t *restrict vc);
void *cs_cuda_host_alloc(size_t size_val);
int cs_cuda_host_free(void *ptr);
#endif // HAVE_CUDA_OFFLOAD

// Functions that must be defined regardless of the CUDA support
// being activated or not. They are used by the FORTRAN modules.
void cs_cuda_attempt_host_alloc(cs_real_t **ptr, int elems);
void cs_cuda_attempt_host_free(cs_real_t *ptr);

#ifdef __cplusplus
}
#endif //__cplusplus
#endif // __CS_CUDA_H__
