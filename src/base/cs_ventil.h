#ifndef __CS_VENTIL_H__
#define __CS_VENTIL_H__

/*============================================================================
 * Management of fans
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

typedef struct _cs_ventil_t cs_ventil_t;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get the number of fans.
 *
 * Fortran interface:
 *
 * SUBROUTINE TSTVTL
 * *****************
 *
 * INTEGER          NBRVTL         : --> : number of fans
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
 * SUBROUTINE DEFVTL
 * *****************
 *
 * INTEGER          DIMMOD     : <-- : Fan model dimension:
 *                             :     : 1: constant_f; 2: force_profile;
 *                             :     : 3: force_profile + tangential couple
 * INTEGER          DIMVTL     : <-- : Fan dimension:
 *                             :     : 2: pseudo-2D (extruded mesh)
 *                             :     : 3: 3D (standard)
 * DOUBLE PRECISION XYZVT1(3)  : <-- : Coo. of the axis point in upstream face
 * DOUBLE PRECISION XYZVT2(3)  : <-- : Coo. of the axis point in downstream face
 * DOUBLE PRECISION RVVT       : <-- : Fan radius
 * DOUBLE PRECISION RPVT       : <-- : Blades radius
 * DOUBLE PRECISION RMVT       : <-- : Hub radius
 * DOUBLE PRECISION CCARAC(3)  : <-- : Coefficients of degre 0, 1 and 2
 *                             :     : of the characteristic curve
 * DOUBLE PRECISION TAUVT      : <-- : Fan axial couple
 *----------------------------------------------------------------------------*/

void CS_PROCF (defvtl, DEFVTL)
(
 const cs_int_t  *const  dimmod,
 const cs_int_t  *const  dimvtl,
 const cs_real_t         xyzvt1[3],
 const cs_real_t         xyzvt2[3],
 const cs_real_t  *const rvvt,
 const cs_real_t  *const rpvt,
 const cs_real_t  *const rmvt,
 const cs_real_t         ccarac[3],
 const cs_real_t  *const tauvt
);

/*----------------------------------------------------------------------------
 * Build the list of cells associated to the fans
 *
 * Fotrtran interface:
 *
 * SUBROUTINE INIVTL
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
 * Calculate the flows through the fans
 *
 * Fortran interface:
 *
 * SUBROUTINE DEBVTL
 * *****************
 *
 * DOUBLE PRECISION FLUMAS(*)      : <-- : Interior faces mass flux
 * DOUBLE PRECISION FLUMAB(*)      : <-- : Boundary faces mass flux
 * DOUBLE PRECISION RHOFAC(*)      : <-- : Density at cells
 * DOUBLE PRECISION RHOFAB(*)      : <-- : Density at boundary faces
 * DOUBLE PRECISION DEBENT(NBRVTL) : --> : Inlet flow through the fan
 * DOUBLE PRECISION DEBSOR(NBRVTL) : --> : Outlet flow through the fan
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
 * Calculate the force induced by the fans (needs a previous calculation
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
 *   dim_modele    <-- Fan model dimension:
 *                     1: constant_f
 *                     2: force_profile
 *                     3: force_profile + tangential couple
 *   dim_ventil    <-- Fan dimension:
 *                     2: pseudo-2D (extruded mesh)
 *                     3: 3D (standard)
 *   coo_axe_amont <-- Coo. of the axis point in upstream face
 *   coo_axe_aval  <-- Coo. of the axis point in downstream face
 *   ray_ventil    <-- Fan radius
 *   ray_pales     <-- Blades radius
 *   ray_moyeu     <-- Hub radius
 *   coeff_carac   <-- Coefficients of degre 0, 1 and 2 of
                       the characteristic curve
 *   couple_axial  <-- Fan axial couple
 *----------------------------------------------------------------------------*/

void
cs_ventil_definit(const cs_int_t   dim_modele,
                  const cs_int_t   dim_ventil,
                  const cs_real_t  coo_axe_amont[3],
                  const cs_real_t  coo_axe_aval[3],
                  const cs_real_t  ray_ventil,
                  const cs_real_t  ray_pales,
                  const cs_real_t  ray_moyeu,
                  const cs_real_t  coeff_carac[3],
                  const cs_real_t  couple_axial);

/*----------------------------------------------------------------------------
 * Destroy the structures associated to fans
 *----------------------------------------------------------------------------*/

void
cs_ventil_detruit_tous(void);

/*----------------------------------------------------------------------------
 * Looks for the cells belonging to the different fans.
 *
 * parameters:
 *   mesh            <-- associated mesh structure
 *   mesh_quantities <-- mesh quantities
 *----------------------------------------------------------------------------*/

void
cs_ventil_cree_listes(const cs_mesh_t             *mesh,
                      const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Calculate the flows through the fans
 *
 * parameters:
 *   mesh           <-- mesh structure
 *   mesh_qantities <-- mesh quantities
 *   flux_masse_fac <-- interior faces mass flux
 *   flux_masse_fbr <-- boundary faces mass flux
 *   densite_cel    <-- density at cells
 *   densite_fbr    <-- density at boundary faces
 *----------------------------------------------------------------------------*/

void
cs_ventil_calcul_debits(const cs_mesh_t             *mesh,
                        const cs_mesh_quantities_t  *mesh_quantities,
                        const cs_real_t             flux_masse_fac[],
                        const cs_real_t             flux_masse_fbr[],
                        const cs_real_t             densite_cel[],
                        const cs_real_t             densite_fbr[]);

/*----------------------------------------------------------------------------
 * Calculate the force induced by the fans (needs a previous calculation
 * of the flows through each fan).
 *
 * The induced force is added to the array CRVXEP (which can have other
 * other contributions).
 *
 * parameters:
 *   mesh_quantities <-- mesh quantities
 *   idim_source     <-- Dimension associated to the source term of velocity
 *                       (1: X; 2: Y; 3: Z)
 *   t_source        <-> Explicit source term for the velocity
 *----------------------------------------------------------------------------*/

void
cs_ventil_calcul_force(const cs_mesh_quantities_t  *mesh_quantities,
                       const cs_int_t               idim_source,
                       cs_real_t                    t_source[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VENTIL_H__ */
