/*============================================================================
 * Source terms computation for compressible homogeneous two-phase flow model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_config.h"
#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_hgn_thermo.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hgn_source_terms_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_hgn_source_terms_step.c
 *
 *  \brief Return to equilibrium source terms computation for volume, mass,
 *         energy fractions in compressible homogeneous two-phase model.
 *
 * A fractional step method used to solve the whole homogeneous two-phase. The
 * last two steps solve the equations on the volume, mass, energy fractions:
 *
 *    1 ) a convection step perfomed by solving a pure convection equation on
 *        the volume, mass, energy fractions.
 *
 *    2 ) an update of these three fractions to account for the return to
 *        equilibrium source terms. This update is deduced by solving the ODE
 *        associated to the source term for each fraction.
 *
 *  Note that the energy, velocity and density remain constant throuh this
 *  fractional step but the pressure and the temperature depend on the
 *  fractions and thus evolve. They are updated at the end of the step
 *  using the thermodynamic relation defined in \ref cs_hg_thermo.c.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_hgn_source_terms_step(void);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_hgn_source_terms_step(void)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  cs_hgn_source_terms_step(m);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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

void
cs_hgn_source_terms_step(const cs_mesh_t *m)
{
  const cs_lnum_t n_cells = m->n_cells;

  cs_real_t *dt = CS_F_(dt)->val;
  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t *cvar_pr = CS_F_(p)->val;
  cs_real_t *cvar_energ = CS_F_(e_tot)->val;
  cs_real_t *cvar_fracv = CS_F_(volume_f)->val;
  cs_real_t *cvar_fracm = CS_F_(mass_f)->val;
  cs_real_t *cvar_frace = CS_F_(energy_f)->val;
  cs_real_t *cvar_tempk = CS_F_(t_kelvin)->val;

  cs_real_t *ei, *v;
  BFT_MALLOC(ei, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(v, m->n_cells_with_ghosts, cs_real_t);

  cs_real_t *alpha_eq, *y_eq, *z_eq;
  BFT_MALLOC(alpha_eq, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(y_eq, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(z_eq, m->n_cells_with_ghosts, cs_real_t);

  cs_real_t *relax_tau;
  BFT_MALLOC(relax_tau, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t u2 = cs_math_3_norm(vel[cell_id]);
    ei[cell_id] = cvar_energ[cell_id] - 0.5*u2;
    v[cell_id] = 1. / crom[cell_id];

    /* initialize relaxation time scale: instantaneous return to equilibrium */
    relax_tau[cell_id] = 1.e-30;

    /* Compute the equilibrium fractions alpha_eq, y_eq, z_eq for a given
       specific volume v and a given specific energy e */
    cs_hgn_thermo_eq(ei[cell_id],
                     v[cell_id],
                     alpha_eq + cell_id,
                     y_eq + cell_id,
                     z_eq + cell_id);
  }

  cs_user_hgn_thermo_relax_time(m,
                                alpha_eq,
                                y_eq,
                                z_eq,
                                ei,
                                v,
                                relax_tau);

  /* Update the volume fraction, mass fraction and energy fraction using the
   * equilibrium fractions alpha_eq, y_eq, z_eq computed above. */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t exp_tau = exp(-dt[cell_id]/relax_tau[cell_id]);

    cvar_fracv[cell_id] =  exp_tau*cvar_fracv[cell_id]
                         - (exp_tau-1.)*alpha_eq[cell_id];

    cvar_fracm[cell_id] =   exp_tau*cvar_fracm[cell_id]
                         - (exp_tau-1.)*y_eq[cell_id];

    cvar_frace[cell_id] = exp_tau*cvar_frace[cell_id]
                         - (exp_tau-1.)*z_eq[cell_id];
  }

  /* Update pressure and temperature */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_hgn_thermo_pt(cvar_fracv[cell_id],
                     cvar_fracm[cell_id],
                     cvar_frace[cell_id],
                     ei[cell_id],
                     v[cell_id],
                     &cvar_tempk[cell_id],
                     &cvar_pr[cell_id]);
  }

  BFT_FREE(ei);
  BFT_FREE(v);
  BFT_FREE(alpha_eq);
  BFT_FREE(y_eq);
  BFT_FREE(z_eq);
  BFT_FREE(relax_tau);

  /* synchronize */
  const cs_halo_t  *halo = cs_glob_mesh->halo;
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cvar_fracv);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cvar_fracm);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cvar_frace);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cvar_tempk);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cvar_pr);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
