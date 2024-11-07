#pragma once

/*============================================================================
 * Private memory handling wrappersfor CUDA
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

#if defined(HAVE_CUDA)

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the common
 * host-device memory management functions from cs_mem.cpp, and
 * not directly by the user.
 *============================================================================*/

void
cs_mem_cuda_set_prefetch_stream(cudaStream_t  stream);

void *
cs_mem_cuda_malloc_device(size_t        n,
                          const char   *var_name,
                          const char   *file_name,
                          int           line_num);

void *
cs_mem_cuda_malloc_host(size_t        n,
                        const char   *var_name,
                        const char   *file_name,
                        int           line_num);

void *
cs_mem_cuda_malloc_managed(size_t        n,
                           const char   *var_name,
                           const char   *file_name,
                           int           line_num);

void
cs_mem_cuda_free(void         *p,
                 const char   *var_name,
                 const char   *file_name,
                 int           line_num);

void
cs_mem_cuda_free_host(void         *p,
                      const char   *var_name,
                      const char   *file_name,
                      int           line_num);

void
cs_mem_cuda_copy_h2d(void         *dst,
                     const void   *src,
                     size_t        size);

void
cs_mem_cuda_copy_h2d_async(void        *dst,
                           const void  *src,
                           size_t       size);

void
cs_mem_cuda_copy_d2h(void        *dst,
                     const void  *src,
                     size_t       size);

void
cs_mem_cuda_copy_d2h_async(void        *dst,
                           const void  *src,
                           size_t       size);

void
cs_mem_cuda_prefetch_h2d(const void  *dst,
                         size_t       size);

void
cs_mem_cuda_prefetch_d2h(const void  *dst,
                         size_t       size);

void
cs_mem_cuda_copy_d2d(void        *dst,
                     const void  *src,
                     size_t       size);

void *
cs_mem_cuda_get_host_ptr(const void  *ptr);

bool
cs_mem_cuda_is_device_ptr(const void  *ptr);

void
cs_mem_cuda_set_advise_read_mostly(const void  *ptr,
                                   size_t       size);

void
cs_mem_cuda_unset_advise_read_mostly(const void  *ptr,
                                     size_t       size);

#endif  /* CS_HAVE_CUDA */

/*----------------------------------------------------------------------------*/
