#ifndef __CS_ALE_H__
#define __CS_ALE_H__

/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_domain.h"
#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * ALE type
 *----------------------------------------------------------------------------*/

/*! ALE computation type */

typedef enum {

  CS_ALE_NONE = 0,    /*!< no ALE */
  CS_ALE_LEGACY = 1,  /*!< Mesh deformation computed with legacy operators */
  CS_ALE_CDO = 2      /*!< Mesh deformation computed with CDO */

} cs_ale_type_t;

/*! ALE data */

typedef struct {

  int          *impale;    /*!< 1st component of low-level BC */
  int          *bc_type;   /*!< ALE BC type code */

} cs_ale_data_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

extern cs_ale_type_t  cs_glob_ale;

extern cs_ale_data_t  *cs_glob_ale_data;

/* defined in albase.f90 (bind(C, name='') :: ...) */

extern int cs_glob_ale_n_ini_f; /*!< Number of sub-iterations for fluid
                                  flow initialization */

extern int cs_glob_ale_need_init; /*!< Indicate whether an iteration to
                                    initialize ALE is required */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocation of ialtyb and impale for the ALE structure.
 */
/*----------------------------------------------------------------------------*/

void cs_ale_allocate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell and face centers of gravity, cell volumes
 *         and update bad cells.
 *
 * \param[out]       min_vol        Minimum cell volume
 * \param[out]       max_vol        Maximum cell volume
 * \param[out]       tot_vol        Total cell volume
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_mesh_quantities(cs_real_t  *min_vol,
                              cs_real_t  *max_vol,
                              cs_real_t  *tot_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project the displacement on mesh vertices (solved on cell center).
 *
 * \param[in]       ale_bc_type   Type of boundary for ALE
 * \param[in]       meshv         Mesh velocity
 * \param[in]       gradm         Mesh velocity gradient
 *                                (du_i/dx_j : gradv[][i][j])
 * \param[in]       claale        Boundary conditions A
 * \param[in]       clbale        Boundary conditions B
 * \param[in]       dt            Time step
 * \param[out]      disp_proj     Displacement projected on vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_project_displacement(const int           ale_bc_type[],
                            const cs_real_3_t  *meshv,
                            const cs_real_33_t  gradm[],
                            const cs_real_3_t  *claale,
                            const cs_real_33_t *clbale,
                            const cs_real_t    *dt,
                            cs_real_3_t        *disp_proj);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update mesh in the ALE framework.
 *
 * \param[in]       itrale        number of the current ALE iteration
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_mesh(int  itrale);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update ALE BCs for required for the fluid
 *
 * \param[out]      ale_bc_type   type of ALE bcs
 * \param[out]      b_fluid_vel   Fluid velocity at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_bcs(int         *ale_bc_type,
                  cs_real_3_t *b_fluid_vel);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a Poisson equation on the mesh velocity in ALE framework.
 *
 * It also updates the mesh displacement
 * so that it can be used to update mass fluxes (due to mesh displacement).
 *
 * \param[in]  iterns  Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_solve_mesh_velocity(int  iterns);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the mesh velocity solving with CDO
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_activate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if mesh velocity solving with CDO is activated
 *
 * \return true ifmesh velocity solving with CDO is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_ale_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup the equations related to mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_init_setup(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup the equations solving the mesh velocity
 *
 * \param[in]   domain       pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_setup_boundaries(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for the equation of the mesh velocity
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_finalize_setup(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the ALE mesh velocity solving
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read ALE data from restart file.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_restart_read(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write ALE data from restart file.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_restart_write(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ALE_H__ */
