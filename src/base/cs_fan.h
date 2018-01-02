#ifndef __CS_FAN_H__
#define __CS_FAN_H__

/*============================================================================
 * Fan modeling through velocity source terms.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Structure definition
 *============================================================================*/

typedef struct _cs_fan_t cs_fan_t;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the flows through the fans
 *
 * Fortran interface:
 *
 * subroutine debvtl
 * *****************
 *
 * double precision flumas(*)      : <-- : interior faces mass flux
 * double precision flumab(*)      : <-- : boundary faces mass flux
 * double precision rhofac(*)      : <-- : density at cells
 * double precision rhofab(*)      : <-- : density at boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (debvtl, DEBVTL)
(
 cs_real_t  flumas[],
 cs_real_t  flumab[],
 cs_real_t  rhofac[],
 cs_real_t  rhofab[]
);

/*----------------------------------------------------------------------------
 * Compute the force induced by the fans (needs a previous calculation
 * of the flows through each fan).
 *
 * The induced force is added to the array crvxep (which can have other
 * contributions).
 *
 * Fortran interface:
 *
 * subroutine tsvvtl
 * *****************
 *
 * parameters:
 *  crvexp         <-> Explicit source term (velocity)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tsvvtl, TSVVTL)
(
 cs_real_3_t  crvexp[]
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fan definition (added to the ones previously defined)
 *
 * Fans are handled as explicit momentum source terms at the given location,
 * based on the fan's axis and diameter.
 * The fan's pressure characteristic curve is defined by 3 coefficients,
 * such that:
 *   delta P = C_0 + C_1.flow + C_2.flow^2
 * An axial torque may also be defined for the 3D model.
 *
 * parameters:
 *   fan_dim             <-- fan dimension:
 *                           2: pseudo-2D (extruded mesh)
 *                           3: 3D (standard)
 *   inlet_axis_coords   <-- intersection coords. of axis and inlet face
 *   outlet_axis_coords  <-- intersection coords. od axis and outlet face
 *   fan_radius          <-- fan radius
 *   blades_radius       <-- blades radius
 *   hub_radius          <-- hub radius
 *   curve_coeffs        <-- coefficients of degre 0, 1 and 2 of
 *                           the pressure drop/flow rate
 *                           characteristic curve
 *   axial_torque        <-- fan axial torque
 *----------------------------------------------------------------------------*/

void
cs_fan_define(int              fan_dim,
              const cs_real_t  inlet_axis_coords[3],
              const cs_real_t  outlet_axis_coords[3],
              cs_real_t        fan_radius,
              cs_real_t        blades_radius,
              cs_real_t        hub_radius,
              const cs_real_t  curve_coeffs[3],
              cs_real_t        axial_torque);

/*----------------------------------------------------------------------------
 * Destroy the structures associated with fans.
 *----------------------------------------------------------------------------*/

void
cs_fan_destroy_all(void);

/*----------------------------------------------------------------------------
 * Return number of fans.
 *
 * returns:
 *   number of defined fans
 *----------------------------------------------------------------------------*/

int
cs_fan_n_fans(void);

/*----------------------------------------------------------------------------
 * Log fans definition setup information.
 *----------------------------------------------------------------------------*/

void
cs_fan_log_setup(void);

/*----------------------------------------------------------------------------
 * Log fan information for a given iteration.
 *----------------------------------------------------------------------------*/

void
cs_fan_log_iteration(void);

/*----------------------------------------------------------------------------
 * Define the cells belonging to the different fans.
 *
 * parameters:
 *   mesh            <-- associated mesh structure
 *   mesh_quantities <-- mesh quantities
 *----------------------------------------------------------------------------*/

void
cs_fan_build_all(const cs_mesh_t             *mesh,
                 const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute the flows through the fans.
 *
 * parameters:
 *   mesh            <-- mesh structure
 *   mesh_quantities <-- mesh quantities
 *   i_mass_flux     <-- interior faces mass flux
 *   b_mass_flux     <-- boundary faces mass flux
 *   c_rho           <-- density at cells
 *   b_rho           <-- density at boundary faces
 *----------------------------------------------------------------------------*/

void
cs_fan_compute_flows(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                     const cs_real_t              i_mass_flux[],
                     const cs_real_t              b_mass_flux[],
                     const cs_real_t              c_rho[],
                     const cs_real_t              b_rho[]);

/*----------------------------------------------------------------------------
 * Compute the force induced by the fans (needs a previous calculation
 * of the flows through each fan).
 *
 * The induced force is added to the array CRVXEP (which can have other
 * other contributions).
 *
 * parameters:
 *   mesh_quantities <-- mesh quantities
 *   source_t        <-> explicit source term for the velocity
 *----------------------------------------------------------------------------*/

void
cs_fan_compute_force(const cs_mesh_quantities_t  *mesh_quantities,
                     cs_real_3_t                  source_t[]);

/*----------------------------------------------------------------------------
 * Flag the cells belonging to the different fans
 * (by the fan id, -1 otherwise)
 *
 * parameters:
 *   mesh        <-- associated mesh structure
 *   cell_fan_id --> indicator by cell
 *----------------------------------------------------------------------------*/

void
cs_fan_flag_cells(const cs_mesh_t  *mesh,
                  int               cell_fan_id[]);

/*----------------------------------------------------------------------------
 * Selection function for cells belonging to fans.
 *
 * This function may be used for the definition of postprocessing meshes.
 *
 * param
 * \param[in, out]  input    pointer to input (unused here)
 * \param[out]      n_cells  number of selected cells
 * \param[out]      cell_ids array of selected cell ids (0 to n-1 numbering)
 */
/*----------------------------------------------------------------------------*/

void
cs_fan_cells_select(void         *input,
                    cs_lnum_t    *n_cells,
                    cs_lnum_t   **cell_ids);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FAN_H__ */
