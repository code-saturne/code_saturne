#ifndef CS_TURBULENCE_RIJ_H
#define CS_TURBULENCE_RIJ_H

/*============================================================================
 * Rij-epsilon turbulence model.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Solve the \f$ R_{ij} - \epsilon \f$ for incompressible flows or
 *         slightly compressible flows for one time step.
 *
 * Please refer to the
 * <a href="../../theory.pdf#rijeps"><b>\f$ R_{ij} - \epsilon \f$ model</b></a>
 * section of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#turrij"><b>turrij</b></a> section.
 *
 * \param[in]     phase_id     turbulent phase id (-1 for single phase flow)
 !*/
/*-----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij(int phase_id);

/*----------------------------------------------------------------------------*/
/*! \brief Solve the equation on alpha in the framework of the Rij-EBRSM model.
 *
 * Also called for alpha of scalars for EB-DFM.
 *
 * \param[in]  f_id          field id of alpha variable
 * \param[in]  phase_id      turbulent phase id (-1 for single phase flow)
 * \param[in]  c_durbin_l    constant for the Durbin length
 !*/
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_solve_alpha(int        f_id,
                              int        phase_id,
                              cs_real_t  c_durbin_l);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize Rij-epsilon variables based on reference quantities.
 *
 * If uref is not provided (0 or negative), values are set at a large
 * negative value (-cs_math_big_r) to allow for later checks.
 *
 * \param[in]  uref    characteristic flow velocity
 * \param[in]  almax   characteristic macroscopic length of the domain
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_init_by_ref_quantities(cs_real_t  uref,
                                         cs_real_t  almax);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (coupled components version).
 *
 * \param[in]  phase_id   turbulent phase id (-1 for single phase flow)
 * \param[in]  n_cells    number of cells
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_clip(int        phase_id,
                       cs_lnum_t  n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the turbulent viscosity for the Reynolds Stress model.
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_mu_t(int  phase_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute anisotropic turbulent viscosity for RSM models
 *
 * \param[in]  mq        mesh quantities
 * \param[in]  n_cells   number of cells
 * \param[in]  phase_id  turbulent phase id (-1 for single phase flow)
 * \param[in]  idfm      use DFM model ?
 * \param[in]  iggafm    use AFM model ?
 * \param[in]  iebdfm    use EBDFM model ?
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_anisotropic_mu_t
(
 const cs_mesh_quantities_t  *mq,
 cs_lnum_t                    n_cells,
 int                          phase_id,
 bool                         idfm,
 bool                         iggafm,
 bool                         iebdfm
);

/*----------------------------------------------------------------------------*/
/*! \brief Compute Rusanov equivalent diffusivity of the model.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_compute_rusanov(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute, once per time step, the exact Riemann interface state
 *        for the coupled {u, R} system on every interior and boundary
 *        face, and store it in the shared fields "i_velocity",
 *        "i_reynolds_stress", "b_velocity", "b_reynolds_stress"
 *        (these fields are created in cs_setup.cpp when
 *        rij_discretization_scheme == CS_RIJ_SCHEME_GODUNOV).
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_godunov_interface_states(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute, once per time step, the exact Riemann interface state
 *        for the scalar system on every interior and boundary face.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_godunov_interface_states_scalar(cs_field_t *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the exact Riemann divergence of R flux.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rij_godunov_div_rij_flux(const cs_real_t  crom[],
                                       const cs_real_t  brom[],
                                       cs_real_3_t     *tflmas,
                                       cs_real_3_t     *tflmab);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENCE_RIJ_H */
