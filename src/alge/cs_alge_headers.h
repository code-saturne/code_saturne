#ifndef __CS_ALGE_HEADERS_H__
#define __CS_ALGE_HEADERS_H__

/*============================================================================
 * Global code_saturne headers file for easier include
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_bad_cells_regularisation.h"
#include "alge/cs_balance_by_zone.h"
#include "alge/cs_balance.h"
#include "alge/cs_benchmark.h"
#include "alge/cs_blas.h"
#include "alge/cs_bw_time_diff.h"
#ifdef __cplusplus
#include "alge/cs_cell_to_vertex.h"
#endif
#include "alge/cs_convection_diffusion.h"
#include "alge/cs_divergence.h"
#include "alge/cs_face_viscosity.h"
#include "alge/cs_gradient.h"
#include "alge/cs_gradient_boundary.h"
#include "alge/cs_grid.h"
#include "alge/cs_matrix_assembler.h"
#include "alge/cs_matrix_building.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_tuning.h"
#include "alge/cs_matrix_util.h"
#include "alge/cs_multigrid.h"
#include "alge/cs_multigrid_smoother.h"
#include "alge/cs_param_amg.h"
#include "alge/cs_param_mumps.h"
#include "alge/cs_param_saddle.h"
#include "alge/cs_param_sles.h"
#include "alge/cs_param_sles_setup.h"
#include "alge/cs_saddle_solver.h"
#include "alge/cs_saddle_solver_setup.h"
#include "alge/cs_sles_default.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_it.h"
#include "alge/cs_sles_pc.h"
#ifdef __cplusplus
#include "alge/cs_vertex_to_cell.h"
#endif

#if defined(HAVE_HYPRE)
#include "alge/cs_sles_hypre.h"
#endif

#if defined(HAVE_MUMPS)
#include "alge/cs_sles_mumps.h"
#endif

#if defined(HAVE_PETSC)
#include "alge/cs_sles_petsc.h"
#endif

#if defined(HAVE_AMGX)
#include "alge/cs_sles_amgx.h"
#endif

#if defined(HAVE_CUDSS)
#include "alge/cs_sles_cudss.h"
#endif

/*----------------------------------------------------------------------------*/

#endif /* __CS_ALGE_HEADERS_H__ */
