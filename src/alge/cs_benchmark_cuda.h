#ifndef __CS_BENCHMARK_CUDA_H__
#define __CS_BENCHMARK_CUDA_H__

/*============================================================================
 * Low-level operator benchmarking with CUDA.
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * SpMV extradiagonal terms using native to face-based array and scatter
 * approach, handling conflicts through atomic add.
 *
 * Symmetric matrix case.
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   y               <-> vector
 *----------------------------------------------------------------------------*/

void
cs_mat_vec_exdiag_native_sym_cuda(cs_lnum_t           n_faces,
                                  const cs_lnum_2_t  *face_cell,
                                  const cs_real_t    *xa,
                                  const cs_real_t    *x,
                                  cs_real_t          *y);

/*----------------------------------------------------------------------------*/

END_C_DECLS


#endif /* __CS_BENCHMARK_CUDA_H__ */
