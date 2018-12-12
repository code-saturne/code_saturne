#ifndef __CS_HGN_SOURCE_TERMS_STEP_H__
#define __CS_HGN_SOURCE_TERMS_STEP_H__

/*============================================================================
 * Source terms computation for compressible homogeneous two-phase flow model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the source terms for the two-phase flow model
 *
 * Update of these three fractions to account for the return to
 * equilibrium source terms. This update is deduced by solving the ODE
 * associated to the source term for each fraction.
 *
 * First we compute the equilibrium fractions and second the fractions are
 * relaxed towards those equilibrium values. Relaxation timescale \f$ \tau \f$
 * is to be provided by the user (equal to 1e-30 by default).
 *
 * Note that the energy, velocity and density remain constant throuh this
 * fractional step but the pressure and the temperature depend on the
 * fractions and thus evolve. They are updated at the end of the step
 * using the thermodynamic relation defined in \ref cs_hg_thermo.c.
 *
 * \param[in]     m   pointer to mesh
 */
/*----------------------------------------------------------------------------*/

void cs_hgn_source_terms_step(const cs_mesh_t *m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HGN_SOURCE_TERMS_STEP_H__ */
