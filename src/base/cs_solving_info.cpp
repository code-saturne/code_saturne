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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_mem.h"
#include "base/cs_solving_info.h"

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default cs_solving_info_t constructor
 */
/*----------------------------------------------------------------------------*/

cs_solving_info_t::cs_solving_info_t()
{
  n_it = 0;
  rhs_norm = 0.;
  res_norm = 0.;
  derive = 0.;
  l2residual = 0.;

  n_clip_min = 0;
  n_clip_max = 0;

  for (int i = 0; i < CS_SOLVING_INFO_MAX_DIM; i++) {
    n_clip_min_comp[i] = 0;
    n_clip_max_comp[i] = 0;
    min_pre_clip[i] = 0.;
    max_pre_clip[i] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_solving_info_t constructor with main convergence values
 *
 * \param[in]  rhs_norm_    right hand side norm
 * \param[in]  res_norm_    normed residual
 * \param[in]  derive_      norm of the time derivative
 * \param[in]  l2residual_  L2 time residual
 */
/*----------------------------------------------------------------------------*/

cs_solving_info_t::cs_solving_info_t
(
  double  rhs_norm_,
  double  res_norm_,
  double  derive_,
  double  l2residual_
)
{
  *this = cs_solving_info_t();

  rhs_norm = rhs_norm_;
  res_norm = res_norm_;
  derive = derive_;
  l2residual = l2residual_;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset all solving info values to their defaults
 */
/*----------------------------------------------------------------------------*/

void
cs_solving_info_t::reset()
{
  *this = cs_solving_info_t();
}

/*----------------------------------------------------------------------------
 * cs_solving_info_t custom new and delete operators
 *----------------------------------------------------------------------------*/

void *
cs_solving_info_t::operator new
(
  std::size_t size   /*!<[in] Size to allocate */
)
{
  return cs_mem_malloc(1, size, "cs_solving_info::new",
                       __FILE__, __LINE__);
}

void
cs_solving_info_t::operator delete
(
  void *ptr   /*!<[in] Size to allocate */
)
{
  cs_mem_free(ptr, "cs_solving_info::delete",
              __FILE__, __LINE__);
}

/*----------------------------------------------------------------------------*/
