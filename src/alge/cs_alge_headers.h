#ifndef __CS_ALGE_HEADERS_H__
#define __CS_ALGE_HEADERS_H__

/*============================================================================
 * Global code_saturne headers file for easier include
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_bad_cells_regularisation.h"
#include "cs_balance_by_zone.h"
#include "cs_balance.h"
#include "cs_benchmark.h"
#include "cs_blas.h"
#include "cs_convection_diffusion.h"
#include "cs_divergence.h"
#include "cs_face_viscosity.h"
#include "cs_gradient.h"
#include "cs_grid.h"
#include "cs_matrix_assembler.h"
#include "cs_matrix_building.h"
#include "cs_matrix_default.h"
#include "cs_matrix.h"
#include "cs_matrix_tuning.h"
#include "cs_matrix_util.h"
#include "cs_multigrid.h"
#include "cs_multigrid_smoother.h"
#include "cs_param_sles.h"
#include "cs_sles_default.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_sles_pc.h"

#if defined(HAVE_HYPRE)
#include "cs_sles_hypre.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

#if defined(HAVE_AMGX)
#include "cs_sles_amgx.h"
#endif

/*----------------------------------------------------------------------------*/

#endif /* __CS_ALGE_HEADERS_H__ */
