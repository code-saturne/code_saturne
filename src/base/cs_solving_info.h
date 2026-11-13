#ifndef CS_SOLVING_INFO_H
#define CS_SOLVING_INFO_H

/*============================================================================
 * Variable solution info (statistics).
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

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Maximum dimension for solving information clipping arrays */
#define CS_SOLVING_INFO_MAX_DIM 9

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure of the solving info
 *----------------------------------------------------------------------------*/

class cs_solving_info_t {

public:

  /* Members */

  int        n_it;
  double     rhs_norm;
  double     res_norm;
  double     derive;
  double     l2residual;

  /* Clipping information (global MPI values after reduction) */

  cs_gnum_t  n_clip_min;
  cs_gnum_t  n_clip_max;
  cs_gnum_t  n_clip_min_comp[CS_SOLVING_INFO_MAX_DIM];
  cs_gnum_t  n_clip_max_comp[CS_SOLVING_INFO_MAX_DIM];
  /*! Pre-clip minimum (valid only if n_clip_min > 0) */
  cs_real_t  min_pre_clip[CS_SOLVING_INFO_MAX_DIM];
  /*! Pre-clip maximum (valid only if n_clip_max > 0) */
  cs_real_t  max_pre_clip[CS_SOLVING_INFO_MAX_DIM];

  /* Constructors */

  cs_solving_info_t();

  cs_solving_info_t
  (
    double  rhs_norm,
    double  res_norm,
    double  derive = 0.,
    double  l2residual = 0.
  );

  /* Public functions */

  void
  reset();

  /* Custom new and delete operators*/

  void *
  operator new(std::size_t size);

  void
  operator delete(void *ptr);

};

/*----------------------------------------------------------------------------*/

#endif /* CS_SOLVING_INFO_H */
