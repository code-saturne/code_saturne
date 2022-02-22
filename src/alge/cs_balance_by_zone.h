#ifndef __CS_BALANCE_BY_ZONE_H__
#define __CS_BALANCE_BY_ZONE_H__

/*============================================================================
 * Scalar balance on zones.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Balance terms */
/* -------------- */

typedef enum {

  CS_BALANCE_VOLUME,              /*!< volume contribution of unsteady terms */
  CS_BALANCE_DIV,                 /*!< volume contribution due to div(rho.u) */
  CS_BALANCE_UNSTEADY,            /*!< contribution of unsteady terms
                                    (volume + div(rho.u)) */
  CS_BALANCE_MASS,                /*!< in and out-flowing mass */
  CS_BALANCE_MASS_IN,             /*!< inflowing mass */
  CS_BALANCE_MASS_OUT,            /*!< outflowing mass */
  CS_BALANCE_INTERIOR_IN,         /*!< inflow through interior faces */
  CS_BALANCE_INTERIOR_OUT,        /*!< outflow through interior faces */
  CS_BALANCE_BOUNDARY_IN,         /*!< inflow through boundary faces */
  CS_BALANCE_BOUNDARY_OUT,        /*!< outflow through boundary faces */
  CS_BALANCE_BOUNDARY_SYM,        /*!< symmetry faces */
  CS_BALANCE_BOUNDARY_WALL,       /*!< wall faces (rough + smooth) */
  CS_BALANCE_BOUNDARY_WALL_S,     /*!< smooth wall faces */
  CS_BALANCE_BOUNDARY_WALL_R,     /*!< rough wall faces */
  CS_BALANCE_BOUNDARY_COUPLED,    /*!< coupled boundary faces (all) */
  CS_BALANCE_BOUNDARY_COUPLED_E,  /*!< externally coupled boundary faces */
  CS_BALANCE_BOUNDARY_COUPLED_I,  /*!< internally coupled boundary faces */
  CS_BALANCE_BOUNDARY_OTHER,      /*!< other boundary faces */
  CS_BALANCE_TOTAL,               /*!< total balance */
  CS_BALANCE_TOTAL_NORMALIZED,    /*!< total balance normalized by square root
                                    of total volume */

  /* End of balance terms */

  CS_BALANCE_N_TERMS

} cs_balance_term_t;

/*! Pressure drop balance terms */
/* ---------------------------- */

typedef enum {

  CS_BALANCE_P_IN,                /*!< inlet pressure: p u */
  CS_BALANCE_P_OUT,               /*!< outlet pressure: p u  */
  CS_BALANCE_P_U2_IN,             /*!< inlet contribution: u^2/2 rho u */
  CS_BALANCE_P_U2_OUT,            /*!< outlet contribution: u^2/2 rho u */
  CS_BALANCE_P_RHOGX_IN,          /*!< inlet contribution: -rho(g.x) u */
  CS_BALANCE_P_RHOGX_OUT,         /*!< outlet contribution: -rho(g.x) u */
  CS_BALANCE_P_U_IN,              /*!< inlet contribution:  u */
  CS_BALANCE_P_U_OUT,             /*!< outlet contribution:  u */
  CS_BALANCE_P_RHOU_IN,           /*!< inlet mass flow:  rho u */
  CS_BALANCE_P_RHOU_OUT,          /*!< outlet mass flow:  rho u */

  /* End of balance terms */

  CS_BALANCE_P_N_TERMS

} cs_balance_p_term_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the different terms of the balance of a given scalar,
 *        on a volume zone defined by selected cell ids/
 *
 * This function computes the balance relative to a given scalar
 * on a selected zone of the mesh.
 * We assume that we want to compute balances (convective and diffusive)
 * at the boundaries of the calculation domain represented below
 * (with different boundary types).
 *
 * In the case of the temperature, the energy balance in Joules will be
 * computed by multiplying by the specific heat.
 *
 * \param[in]     scalar_name         scalar name
 * \param[in]     n_cells_sel         number of selected cells
 * \param[in]     cell_sel_ids        ids of selected cells
 * \param[out]    balance             array of computed balance terms
 *                                    (see \ref cs_balance_term_t)
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone_compute(const char      *scalar_name,
                           cs_lnum_t        n_cells_sel,
                           const cs_lnum_t  cell_sel_ids[],
                           cs_real_t        balance[CS_BALANCE_N_TERMS]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and log the different terms of the balance of a given scalar,
 *        on a volumic zone defined by selection criteria.
 *        The different contributions to the balance are printed in the
 *        run_solver.log.
 *
 * This function computes the balance relative to a given scalar
 * on a selected zone of the mesh.
 * We assume that we want to compute balances (convective and diffusive)
 * at the boundaries of the calculation domain represented below
 * (with different boundary types).
 *
 * The scalar and the zone are selected at the top of the routine
 * by the user.
 * In the case of the temperature, the energy balance in Joules will be
 * computed by multiplying by the specific heat.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone(const char  *selection_crit,
                   const char  *scalar_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 *        on a volume zone defined by selected cell ids/
 *
 * \param[in]     n_cells_sel         number of selected cells
 * \param[in]     cell_sel_ids        ids of selected cells
 * \param[out]    balance             array of computed balance terms
 *                                    (see \ref cs_balance_p_term_t)
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone_compute(cs_lnum_t        n_cells_sel,
                                 const cs_lnum_t  cell_sel_ids[],
                                 cs_real_t        balance[CS_BALANCE_P_N_TERMS]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 * volumic zone defined by the criterion also given as argument.
 * The different contributions are printed in the run_solver.log.
 *
 * \param[in]     selection_crit      zone selection criterion
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone(const char  *selection_crit);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface balance of a given scalar.
 *
 * For interior faces, the flux is counted negatively relative to the given
 * normal (as neighboring interior faces may have differently-aligned normals).
 *
 * For boundary faces, the flux is counted negatively in the outwards-facing
 * direction.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 * \param[in]     normal              outwards normal direction
 */
/*----------------------------------------------------------------------------*/

void
cs_surface_balance(const char       *selection_crit,
                   const char       *scalar_name,
                   const cs_real_t   normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the face by face surface flux of a given scalar, through a
 *        surface area defined by the given face ids.
 *
 * For interior faces, the flux is counted negatively relative to the given
 * normal (as neighboring interior faces may have differently-aligned normals).
 *
 * For boundary faces, the flux is counted negatively in the outwards-facing
 * direction.
 *
 * \param[in]   scalar_name       scalar name
 * \param[in]   normal            outwards normal direction
 * \param[in]   n_b_faces_sel     number of selected boundary faces
 * \param[in]   n_i_faces_sel     number of selected internal faces
 * \param[in]   b_face_sel_ids    ids of selected boundary faces
 * \param[in]   i_face_sel_ids    ids of selected internal faces
 * \param[out]  balance           optional array of computed balance terms
 *                                (see \ref cs_balance_term_t), of
 *                                size CS_BALANCE_N_TERMS, or NULL
 * \param[out]  flux_b_faces      optional surface flux through selected
 *                                boundary faces, or NULL
 * \param[out]  flux_i_faces      optional surface flux through selected
 *                                interior faces, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_flux_through_surface(const char         *scalar_name,
                        const cs_real_t     normal[3],
                        cs_lnum_t           n_b_faces_sel,
                        cs_lnum_t           n_i_faces_sel,
                        const cs_lnum_t     b_face_sel_ids[],
                        const cs_lnum_t     i_face_sel_ids[],
                        cs_real_t          *balance,
                        cs_real_t          *flux_b_faces,
                        cs_real_t          *flux_i_faces);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BALANCE_BY_ZONE_H__ */
