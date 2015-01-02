#ifndef __CS_FAN_H__
#define __CS_FAN_H__

/*============================================================================
 * Fan modeling through velocity source terms.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
 * Get the number of fans.
 *
 * Fortran interface:
 *
 * subroutine tstvtl
 * *****************
 *
 * integer          nbrvtl         : --> : number of fans
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstvtl, TSTVTL)
(
 cs_int_t  *const nbrvtl
);

/*----------------------------------------------------------------------------
 * Adds a fan.
 *
 * Fortran interface:
 *
 * subroutine defvtl
 * *****************
 *
 * integer          dimmod     : <-- : fan model dimension:
 *                             :     : 1: constant_f; 2: force_profile;
 *                             :     : 3: force_profile + tangential couple
 * integer          dimvtl     : <-- : fan dimension:
 *                             :     : 2: pseudo-2D (extruded mesh)
 *                             :     : 3: 3D (standard)
 * double precision xyzvt1(3)  : <-- : coo. of the axis point in inlet face
 * double precision xyzvt2(3)  : <-- : coo. of the axis point in outlet face
 * double precision rvvt       : <-- : fan radius
 * double precision rpvt       : <-- : blades radius
 * double precision rmvt       : <-- : hub radius
 * double precision ccarac(3)  : <-- : coefficients of degre 0, 1 and 2
 *                             :     : of the characteristic curve
 * double precision tauvt      : <-- : fan axial couple
 *----------------------------------------------------------------------------*/

void CS_PROCF (defvtl, DEFVTL)
(
 const cs_int_t   *dimmod,
 const cs_int_t   *dimvtl,
 const cs_real_t   xyzvt1[3],
 const cs_real_t   xyzvt2[3],
 const cs_real_t  *rvvt,
 const cs_real_t  *rpvt,
 const cs_real_t  *rmvt,
 const cs_real_t   ccarac[3],
 const cs_real_t  *tauvt
);

/*----------------------------------------------------------------------------
 * Build the list of cells associated to the fans
 *
 * Fotrtran interface:
 *
 * subroutine inivtl
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (inivtl, INIVTL)
(
 void
);

/*----------------------------------------------------------------------------
 * Mark the fans and associate the fan number to the cells belonging to
 * thus fan, 0 otherwise.
 *
 * Fortran interface:
 *
 * SUBROUTINE NUMVTL (INDIC)
 * *****************
 *
 * INTEGER INDIC(NCELET)       : --> : Fan number (0 if outside the fan)
 *----------------------------------------------------------------------------*/

void CS_PROCF (numvtl, NUMVTL)
(
 cs_int_t  indic[]
);

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
 * double precision debent(nbrvtl) : --> : inlet flow through the fan
 * double precision debsor(nbrvtl) : --> : Outlet flow through the fan
 *----------------------------------------------------------------------------*/

void CS_PROCF (debvtl, DEBVTL)
(
 cs_real_t  flumas[],
 cs_real_t  flumab[],
 cs_real_t  rhofac[],
 cs_real_t  rhofab[],
 cs_real_t  debent[],
 cs_real_t  debsor[]
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
 *  idimts         <-- Dimension associated to the source
 *                     term of velocity (1: X; 2: Y; 3: Z)
 *  crvexp         <-> Explicit source term (velocity)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tsvvtl, TSVVTL)
(
 cs_int_t  *idimts,
 cs_real_t  crvexp[]
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fan definition (added to the ones previously defined)
 *
 * parameters:
 *   model_dim           <-- fan model dimension:
 *                           1: constant_f
 *                           2: force_profile
 *                           3: force_profile + tangential couple
 *   fan_dim             <-- fan dimension:
 *                           2: pseudo-2D (extruded mesh)
 *                           3: 3D (standard)
 *   inlet_axis_coords   <-- intersection coords. of axis and inlet face
 *   outlet_axis_coords  <-- intersection coords. od axis and outlet face
 *   fan_radius          <-- fan radius
 *   blades_radius       <-- blades radius
 *   hub_radius          <-- hub radius
 *   curve_coeffs        <-- coefficients of degre 0, 1 and 2 of
                             the characteristic curve
 *   axial_torque        <-- fan axial torque
 *----------------------------------------------------------------------------*/

void
cs_fan_define(int              model_dim,
              int              fan_dim,
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
 *   source_coo_id   <-- coordinate associated to the source term of velocity
 *                        (0: X; 0: Y; 0: Z)
 *   source_t        <-> explicit source term for the velocity
 *----------------------------------------------------------------------------*/

void
cs_fan_compute_force(const cs_mesh_quantities_t  *mesh_quantities,
                     int                          source_coo_id,
                     cs_real_t                    source_t[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FAN_H__ */
