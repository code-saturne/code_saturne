#ifndef __CS_1D_WALL_THERMAL_H__
#define __CS_1D_WALL_THERMAL_H__

/*============================================================================
 * Modelling the thermal wall with 1D approach
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * 1D Wall Thermal local model type (local to a coupled boundary face)
 *----------------------------------------------------------------------------*/

typedef struct {

  int nppt1d;       /* Number of discretisation cells in the 1D wall
                     * at a given boundary face which is coupled
                     * with a 1D wall thermal module.*/

  int iclt1d;       /* Boundary condition type at the external (pseudo)
                     * wall:
                     - 1: Dirichlet,
                     - 2: Flux condition. */

  cs_real_t eppt1d; /* Thickness of the 1D wall for the nfpt1d boundary faces
                     * which are coupled with a 1D wall thermal module. */

  cs_real_t rgpt1d; /* Geometry of the pseudo wall mesh (refined as a fluid
                     * if rgt1d is smaller than 1.*/

  cs_real_t tept1d; /* External temperature of the pseudo wall in the Dirichlet
                     * case. */

  cs_real_t hept1d; /* External coefficient of transfer in the pseudo wall
                     * under Dirichlet conditions, in W.m^-2.K. */

  cs_real_t fept1d; /* External heat flux in the pseudo wall under the flux
                     * conditions (in W.m^-2, negative value for energy
                     * entering the wall). */

  cs_real_t xlmbt1; /* Thermal diffusivity. */

  cs_real_t rcpt1d; /* Volumetric heat capacity rho C_p of the wall uniform
                     * in thickness (in J.m^-3.K^-1). */

  cs_real_t dtpt1d; /* Physical time step associated with the solved 1D equation
                     * of the pseudo wall (which can be different from the time
                     * step of the calculation). */

  cs_real_t *z;     /* Discretization points coordinates. */

  cs_real_t *t;     /* Temperature at each point of discretization. */

} cs_1d_wall_thermal_local_model_t;

/*----------------------------------------------------------------------------
 * 1D Wall Thermal module type
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_lnum_t nfpt1d;  /* Number of boundary faces which are coupled
                      * with a 1D wall thermal module. */

  cs_gnum_t nfpt1t;  /* Global number of boundary faces which are coupled with
                      * a 1D wall thermal module, i.e. sum of nfpt1d over all
                      * ranks */

  int nmxt1d;        /* Maximum number of discretization cells in 1d wall */

  cs_lnum_t *izft1d; /* Zones of t1d, dimensioned with nfabor */

  cs_lnum_t *ifpt1d; /* Array allowing to mark out the numbers of
                      * the nfpt1d boundary faces which are coupled with
                      * a 1D wall. */

  cs_real_t *tppt1d; /* Initialisation temperature of the wall (uniform in
                      * thickness). During the calculation, the array stores
                      * the temperature of the solid at the fluid/solid
                      * interface. */

  cs_1d_wall_thermal_local_model_t *local_models; /* Array of structures */

} cs_1d_wall_thermal_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to the 1D wall thermal structure */

extern const cs_1d_wall_thermal_t  *cs_glob_1d_wall_thermal;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

 /*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the cs_glob_1d_wall_thermal structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_local_models_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the discretization points coordinates array and
          the temperature at each point of discretization.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_local_models_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the 1D mesh for each face and initialize the temperature.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_mesh_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the 1D equation for a given face
 *
 * parameters:
 * \param[in]   ii   face number
 * \param[in]   tf   fluid temperature at the boundarys
 * \param[in]   hf   exchange coefficient for the fluid
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_solve(cs_lnum_t ii,
                         cs_real_t tf,
                         cs_real_t hf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read the restart file of the 1D-wall thermal module.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_read(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write the restart file of the 1D-wall thermal module.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_write(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the array of structures local_models.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy the global 1d wall thermal structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_1d_wall_thermal.
 */
/*----------------------------------------------------------------------------*/

cs_1d_wall_thermal_t *
cs_get_glob_1d_wall_thermal(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_1D_WALL_THERMAL_H__ */
